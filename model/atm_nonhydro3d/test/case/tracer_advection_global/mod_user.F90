!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort  
  use scale_const, only: &
    PI  => CONST_PI,        &
    P00 => CONST_PRE00,     &
    Grav => CONST_GRAV,     &
    Rdry => CONST_Rdry,     &
    RPlanet => CONST_RADIUS

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D   
  use scale_cubedsphere_cnv, only: &
    CubedSphereCnv_LonLat2CSVec
  
  use mod_exp, only: experiment    
  use mod_atmos_component, only: &
    AtmosComponent
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_mkinit
  public :: USER_setup
  public :: USER_calc_tendency
  public :: USER_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  type, private, extends(experiment) :: Exp_tracer_advection_global
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_tracer_advection
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_tracer_advection_global), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  subroutine USER_mkinit ( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    call exp_manager%Init('tracer_advection_global')

    call exp_manager%SetInitCond( atm%mesh,                &
      atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager, &
      atm%vars%QTRCVARS_manager                            )
    
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

  subroutine USER_setup( atm )
    use scale_tracer, only: &
       TRACER_regist    
    implicit none
    
    class(AtmosComponent), intent(inout) :: atm

    namelist / PARAM_USER / &
       USER_do

    integer :: ierr    
    integer :: iq        
    !------------------------------------------


    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    !-
    call TRACER_REGIST( iq,                   & ! [OUT]
                        1,                    & ! [IN]
                        (/'PTracer'/),        & ! [IN]
                        (/'Passive tracer'/), & ! [IN]
                        (/'1'/)               ) ! [IN]
        
    return
  end subroutine USER_setup

  subroutine USER_calc_tendency( atm )
    use scale_time_manager, only:  TIME_NOWSTEP
    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars   
    use scale_prc 
    implicit none

    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    type(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke2D
    real(RP) :: time

    real(RP), allocatable :: svec(:,:,:)
    real(RP), allocatable :: W(:)
    real(RP), allocatable :: lon3D(:)
    real(RP), allocatable :: lat3D(:)  
    !------------------------------------------

    time = atm%time_manager%dtsec * real( TIME_NOWSTEP - 1, kind=RP )
    !if (PRC_myrank==0) write(*,*) "time=", time
    
    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, lcmesh                               )
      elem3D => lcmesh%refElem3D

      allocate( svec(elem3D%Np,lcmesh%Ne,2), W(elem3D%Np) )
      allocate( lon3D(elem3D%Np), lat3D(elem3D%Np) )

      !$omp parallel do private( ke2D, W, lon3D, lat3D )
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)
        lon3D(:) = lcmesh%lon2D(elem3D%IndexH2Dto3D(:),ke2D)
        lat3D(:) = lcmesh%lat2D(elem3D%IndexH2Dto3D(:),ke2D)
  
        call DCMIP2012_deformation_flow( &
          svec(:,ke,1), svec(:,ke,2), W(:),        &
          lon3D(:), lat3D(:), lcmesh%zlev(:,ke),   &
          time, elem3D%Np                          )  
          
        svec    (:,ke,1) = DENS_hyd%val(:,ke) * svec(:,ke,1) / cos(lat3D(:))
        svec    (:,ke,2) = DENS_hyd%val(:,ke) * svec(:,ke,2)
        MOMZ%val(:,ke  ) = DENS_hyd%val(:,ke) * W(:)
      end do

      call CubedSphereCnv_LonLat2CSVec( &
        lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),          &
        lcmesh%Ne * elem3D%Np, RPlanet, svec(:,:,1), svec(:,:,2),            &
        MOMX%val(:,lcmesh%NeS:lcmesh%NeE), MOMY%val(:,lcmesh%NeS:lcmesh%NeE) )  
        
      deallocate( svec )
      deallocate( lon3D, lat3D )
    end do

    return
  end subroutine USER_calc_tendency

  subroutine USER_update( atm )
    implicit none
    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_update

  !------

  subroutine exp_SetInitCond_tracer_advection( this,                       &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,        &
      GRAV => CONST_GRAV,    &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      CVdry => CONST_CVdry,  &
      PRES00 => CONST_PRE00, &
      Pstd   => CONST_Pstd
    use scale_tracer, only: &
      TRACER_inq_id

    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constT
    use mod_exp, only: &
      TracerLocalMeshField_ptr
  
    implicit none

    class(Exp_tracer_advection_global), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
    type(TracerLocalMeshField_ptr), intent(inout) :: tracer_field_list(:)
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_ymin, dom_ymax
    real(RP), intent(in) :: dom_zmin, dom_zmax
    
    ! Surface state
     real(RP)               :: SFC_PRES           ! surface pressure [Pa]
    ! Environment state
    real(RP)               :: ENV_TEMP = 300.0_RP ! potential temperature of environment [K]
   
    namelist /PARAM_EXP/ &
      SFC_PRES, &
      ENV_TEMP     

    integer :: ke
    integer :: ke2D
    integer :: ierr
    integer :: iq

    real(RP) :: svec(elem%Np,lcmesh%Ne,2)
    real(RP) :: W(elem%Np)
    real(RP) :: lon3D(elem%Np), lat3D(elem%Np)
    !-----------------------------------------------------------------------------

    SFC_PRES = Pstd

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("TRACER_ADVECTION_global_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("TRACER_ADVECTION_global_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    
    ENV_TEMP = 300.0_RP      
    call hydrostatic_calc_basicstate_constT( DENS_hyd, PRES_hyd,                            &
      ENV_TEMP, SFC_PRES, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3), &
      lcmesh, elem )

    !---

    call TRACER_inq_id( "PTracer", iq )
    !$omp parallel do private(ke2D, lon3D, lat3D, W)
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)
      lon3D(:) = lcmesh%lon2D(elem%IndexH2Dto3D(:),ke2D)
      lat3D(:) = lcmesh%lat2D(elem%IndexH2Dto3D(:),ke2D)

      call DCMIP2012_tracer( &
        tracer_field_list(iq)%ptr%val(:,ke),         &
        'q1', lon3D(:), lat3D(:), lcmesh%zlev(:,ke), &
        elem%Np                                      )

      call DCMIP2012_deformation_flow( &
        svec(:,ke,1), svec(:,ke,2), W(:),        &
        lon3D(:), lat3D(:), lcmesh%zlev(:,ke),   &
        0.0_RP, elem%Np                          )
      
      svec(:,ke,1) = DENS_hyd(:,ke) * svec(:,ke,1) / cos(lat3D(:))
      svec(:,ke,2) = DENS_hyd(:,ke) * svec(:,ke,2)
      MOMZ(:,ke) = DENS_hyd(:,ke) * W(:)
    end do
        
    call CubedSphereCnv_LonLat2CSVec( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),  &
      lcmesh%Ne * elem%Np, RPlanet, svec(:,:,1), svec(:,:,2),      &
      MOMX(:,lcmesh%NeS:lcmesh%NeE), MOMY(:,lcmesh%NeS:lcmesh%NeE) )

    return
  end subroutine exp_SetInitCond_tracer_advection

  subroutine exp_geostrophic_balance_correction( this,  &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, &
    lcmesh, elem )
    
    implicit none

    class(Exp_tracer_advection_global), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(inout) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DRHOT(elem%Np,lcmesh%NeA)

    !---------------------------------------------------
    return
  end subroutine exp_geostrophic_balance_correction 

!- private --

!OCL SERIAL
  subroutine DCMIP2012_tracer( q, &
    qtrcname, lon, lat, z, Np )

    implicit none

    integer, intent(in) :: Np
    real(RP), intent(out) :: q(Np)
    character(*), intent(in) :: qtrcname
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: z(Np)
    
    integer :: i
    real(RP) :: d(Np,2)
    real(RP) :: r(Np)

    real(RP) :: Rt                           !< Horizontal half-width of tracers
    real(RP), parameter :: Zt = 1000.0_RP    !< Vertical half-width of tracers
    real(RP) :: lon_c(2)                     !< Initial longitude of first and second tracers
    real(RP), parameter :: lat_c = 0.0_RP    !< Initial latitude of tracers
    real(RP), parameter :: z_c   = 5000.0_RP !< Initial altitude of tracers
    !----------------------------------------------

    Rt = 0.5_RP * RPlanet 
    lon_c(:) = (/ 5.0_RP, 7.0_RP /) * PI / 6.0_RP

    do i=1, 2
      r(:) = RPlanet * acos(sin(lat_c)*sin(lat(:)) + cos(lat_c)*cos(lat(:))*cos(lon(:)-lon_c(i)))
      d(:,i) = min( 1.0_RP, (r(:) / Rt)**2 + ((z(:) - z_c) / Zt)**2 )
    end do

    select case(qtrcname)
    case('q1')
      q(:) = 1.0_RP + 0.5_RP * ( cos(PI * d(:,1)) + cos(PI * d(:,2)) )
    case default
      LOG_ERROR('DCMIP2012_tracer',*) trim(qtrcname)//' is not supported. Check!'
      call PRC_abort
    end select
    return
  end subroutine DCMIP2012_tracer

!OCL SERIAL
  subroutine DCMIP2012_deformation_flow( U, V, W, &
    lon, lat, z, time, Np )

    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: U(Np)
    real(RP), intent(out) :: V(Np)
    real(RP), intent(out) :: W(Np)
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: z(Np)
    real(RP), intent(in) :: time

    real(RP) :: lon2(Np)
    real(RP) :: omg(Np)
    real(RP) :: pres(Np)
    real(RP) :: ud(Np), ua(Np)
    real(RP) :: vd(Np), va(Np)

    real(RP), parameter :: tau = 1036800.0_RP  !> Period of motion [sec]
    real(RP) :: OMG0                           !> Maximum of the vertical pressure velocity in units Pa/s
    real(RP), parameter :: b    = 0.2_RP       !> Normalized pressure depth of the divergent layer
    real(RP), parameter :: ptop = 254.944E2_RP
    real(RP), parameter :: T0   = 300.0_RP     !> Isothermal atmospheric temperature [K]
    !----------------------------------------------

    lon2(:) = lon(:) - 2.0_RP * PI * time / tau
    OMG0 = 23000.0_RP * PI / tau

    pres(:) = P00 * exp(- Grav * z(:) / ( Rdry * T0 ) )
    omg(:) = OMG0 * sin(lon2(:)) * cos(lat(:)) * cos(2.0_RP * PI * time / tau)             &
      * ( 1.0_RP + exp( (ptop - P00) / (b * ptop) ) - exp( (pres(:) - P00) / (b * ptop) )  &
                                                    - exp( (ptop - pres(:)) / (b * ptop) ) )

    ua(:) = 10.0_RP * RPlanet / tau * sin(lon2(:))**2 * sin(2.0_RP*lat(:)) * cos(PI * time / tau) &
          + 2.0_RP * PI * RPlanet / tau *  cos(lat(:))
    va(:) = 10.0_RP * RPlanet / tau * sin(2.0_RP * lon2(:)) * cos(lat(:)) * cos(PI * time / tau)

    ud(:) = OMG0 * RPlanet / (b * ptop) * cos(lon2(:)) * cos(lat(:))**2 * cos(2.0_RP * PI * time / tau) &
            * ( - exp( (pres(:) - P00) / (b * ptop) ) + exp( (ptop - pres(:)) / (b * ptop) )  )
    vd(:) = 0.0_RP

    U(:) = ua(:) + ud(:)
    V(:) = va(:) + vd(:)
    W(:) = - omg(:) / (Grav * pres(:) / ( Rdry * T0 ) )

    return
  end subroutine DCMIP2012_deformation_flow

end module mod_user
