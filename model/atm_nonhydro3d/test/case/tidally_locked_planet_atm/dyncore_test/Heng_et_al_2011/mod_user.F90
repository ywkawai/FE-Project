!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a benchmark test of dynamical core for tidally locked Earth atmosphere proposed by Heng et al. (2011, Monthly Notices of the Royal Astronomical Society).
!!          
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
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  & 
    CVdry => CONST_CVdry,  &      
    PRES00 => CONST_PRE00, &
    Grav => CONST_GRAV,    &
    PI => CONST_PI
  
  use mod_exp, only: experiment

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_line, only: LineElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: MeshField3D

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

  type, private, extends(experiment) :: Exp_Heng2011_TidallyLockedEarth
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_Heng2011_TidallyLockedEarth
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_Heng2011_TidallyLockedEarth), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  real(RP), private, parameter :: kf          = 1.0_RP / ( 86400.0_RP * 1.0_RP  )
  real(RP), private, parameter :: ka          = 1.0_RP / ( 86400.0_RP * 40.0_RP )
  real(RP), private, parameter :: ks          = 1.0_RP / ( 86400.0_RP * 4.0_RP  )
  real(RP), private, parameter :: TEMP_strato = 200.0_RP
  real(RP), private, parameter :: SFCTEMP_eq  = 315.0_RP
  real(RP), private, parameter :: DelT_y      = 60.0_RP
  real(RP), private, parameter :: DelPT_z     = 10.0_RP
  real(RP), private, parameter :: sigb        = 0.7_RP

  !-----------------------------------------------------------------------------

contains

!OCL SERIAL
  subroutine USER_mkinit ( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    call exp_manager%Init('TidallyLockedEarth')

    call exp_manager%SetInitCond( atm%mesh,                &
      atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager, &
      atm%vars%QTRCVARS_manager                            )
    
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL  
  subroutine USER_setup( atm )
    implicit none
    
    class(AtmosComponent), intent(inout) :: atm

    namelist / PARAM_USER / &
       USER_do

    integer :: ierr    
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

    return
  end subroutine USER_setup

!OCL SERIAL  
  subroutine USER_calc_tendency( atm )  
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_calc_tendency

!OCL SERIAL
  subroutine USER_update( atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars, &
      MOMX_p  => ATMOS_PHYTEND_MOMX_ID, &
      MOMY_p  => ATMOS_PHYTEND_MOMY_ID, &
      MOMZ_p  => ATMOS_PHYTEND_MOMZ_ID, &
      RHOH_p  => ATMOS_PHYTEND_RHOH_ID
    use scale_localmeshfield_base, only: LocalMeshFieldBase

    implicit none

    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    type(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke2D

    real(RP), allocatable :: DENS(:), T(:), Teq(:), sig(:), PRES_sfc(:)
    real(RP), allocatable :: rtauT(:), rtauV(:)
    real(RP), allocatable :: lon(:), lat(:)

    real(DP) :: dt
    real(RP) :: Gamm
    !----------------------------------------------------------

    dt = atm%time_manager%dtsec
    gamm = CpDry / CvDry 

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      

      call AtmosVars_GetLocalMeshPhyAuxVars( n, atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                         )
      
      elem3D => lcmesh%refElem3D

      allocate( DENS(elem3D%Np), T(elem3D%Np), Teq(elem3D%Np), sig(elem3D%Np) )
      allocate( PRES_sfc(elem3D%Nnode_h1D**2) )
      allocate( rtauT(elem3D%Np), rtauV(elem3D%Np) )
      allocate( lon(elem3D%Np), lat(elem3D%Np) )

      !$omp parallel do private(DENS, T, Teq, PRES_sfc, sig, lon, lat, rtauT, rtauV, ke2D)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)

        PRES_sfc(:) = PRES%val(elem3D%Hslice(:,1),ke2D)
        sig(:) = PRES%val(:,ke) / PRES_sfc(elem3D%IndexH2Dto3D)
        lon(:) = lcmesh%lon2D(elem3D%IndexH2Dto3D,ke2D)
        lat(:) = lcmesh%lat2D(elem3D%IndexH2Dto3D,ke2D)

        rtauT(:) = ka + (ks - ka) * max( 0.0_RP, (sig(:) - sigb)/(1.0_RP - sigb) ) * cos(lat(:))**4
        rtauV(:) = kf * max( 0.0_RP, (sig(:) - sigb)/(1.0_RP - sigb) )

        Teq(:) = max( TEMP_strato, &
          ( SFCTEMP_eq + DelT_y * cos(PI - lon(:)) * cos(lat(:)) - DelPT_z * log(PRES%val(:,ke)/PRES00) * cos(lat(:))**2 ) &
          * (PRES%val(:,ke)/PRES00)**(Rdry/CPDry)                                                                          )
        
        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)
        T(:) =  PRES%val(:,ke) / ( Rdry * DENS(:) )

        MOMX%val(:,ke) = MOMX%val(:,ke) / ( 1.0_RP + dt * rtauV )
        MOMY%val(:,ke) = MOMY%val(:,ke) / ( 1.0_RP + dt * rtauV )

        !-- For the case of d DRHOT / dt = dens * Cv * ( Teq - T ) / tauT <- ( Ullrich and Jablonowski (2012, JCP) )
        ! DRHOT%val(:,ke) = DRHOT%val(:,ke) &
        !                 - dt * rtauT(:) / gamm * ( 1.0_RP - Teq(:) / T(:) ) * DENS(:) * PT%val(:,ke)     &
        !                 / ( 1.0_RP + dt * rtauT(:) / gamm * ( 1.0_RP + (gamm - 1.0_RP) * Teq(:) / T(:) ) )

        !- For the case of d DRHOT /dt = dens * Cp * ( Teq - T ) / tauT     
        !  <- It is based on the forcing form in Held and Surez in which the temperature evolution equation is assumed to be dT/dt = R/Cp * T/p * dp/dt + (Teq - T) / tauT
        DRHOT%val(:,ke) = DRHOT%val(:,ke) &
                        - dt * rtauT(:) * ( 1.0_RP - Teq(:) / T(:) ) * DENS(:) * PT%val(:,ke)     &
                        / ( 1.0_RP + dt * rtauT(:) * ( 1.0_RP + (gamm - 1.0_RP) * Teq(:) / T(:) ) )  
      end do

      deallocate( DENS, T, Teq, sig, PRES_sfc )
      deallocate( rtauT, rtauV )
      deallocate( lon, lat )
    end do
    
    return
  end subroutine USER_update

  !------

!OCL SERIAL
  subroutine exp_SetInitCond_Heng2011_TidallyLockedEarth( this,                            &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use mod_exp, only: &
      TracerLocalMeshField_ptr

    use scale_const, only: &
      PRES00 => CONST_PRE00, &
      RPlanet => CONST_Radius, &
      CPdry => CONST_CPdry, &     
      CVdry => CONST_CVdry

    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constT
  
    implicit none

    class(Exp_Heng2011_TidallyLockedEarth), intent(inout) :: this
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
    
    real(RP) :: TEMP0    = 264.0_RP
    namelist /PARAM_EXP/ &
      TEMP0
  
    integer :: ierr
    !-----------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("Heng2011_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("Heng2011_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    call hydrostatic_calc_basicstate_constT( DENS_hyd, PRES_hyd, & ! (out)
      TEMP0, PRES00,                                             & ! (in)
      x, y, z, lcmesh, elem                                      ) ! (in)

    return
  end subroutine exp_SetInitCond_Heng2011_TidallyLockedEarth

!OCL SERIAL
  subroutine exp_geostrophic_balance_correction( this,  &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, &
    lcmesh, elem )
    
    implicit none

    class(Exp_Heng2011_TidallyLockedEarth), intent(inout) :: this
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

end module mod_user
