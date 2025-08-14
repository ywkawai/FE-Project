!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of equatorial wave (Test Case 4 in Tomita and Satoh (2004))
!!
!! @author Yuta Kawai, Team SCALE
!!  - Tomita, H and M. Satoh, 2004:
!!    A New Dynamical Framework of Nonhydrostatic Global Model Using the Icosahedral Grid. 
!!    Fluid Dyn. Res., 34, 357 
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

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: MeshField3D

  use mod_user_base, only: UserBase
  use mod_experiment, only: Experiment

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(UserBase) :: User
  contains
    procedure :: mkinit_ => USER_mkinit
    generic :: mkinit => mkinit_
    procedure :: setup_ => USER_setup
    generic :: setup => setup_
    procedure :: calc_tendency => USER_calc_tendency
  end type User

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

  logical, private :: USER_do = .false. !< do user step?

  type(MeshField3D), private :: q_heat
  logical :: is_Qheat_calculated

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init( 'equatorial_wave_global' )
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_equatorial_wave )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL
  subroutine USER_setup( this, atm )
    implicit none
    
    class(User), intent(inout) :: this
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

    call this%UserBase%Setup( atm, USER_do )

    !-
    call q_heat%Init( 'Qheat', 'J/kg.s-1', atm%mesh%ptr_mesh )
    is_Qheat_calculated = .false.

    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    use scale_const, only: &
      Rdry => CONST_Rdry,  &
      CpDry => CONST_CPdry
    use scale_time_manager, only:  TIME_NOWSTEP      
    use scale_localmeshfield_base, only: LocalMeshFieldBase

    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      MOMX_p  => PHYTEND_MOMX_ID, &
      MOMY_p  => PHYTEND_MOMY_ID, &
      MOMZ_p  => PHYTEND_MOMZ_ID, &
      RHOH_p  => PHYTEND_RHOH_ID

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars

    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: PRES, PT

    real(RP) :: tsec    
    real(RP), parameter :: rtau = 1.0_RP / ( 10.0_RP * 86400.0_RP ) ! (10 day)^-1
    real(RP), parameter :: rtau2 = 1.0_RP / ( 20.0_RP * 86400.0_RP ) ! (20 day)^-1
    real(RP) :: fac
    integer :: n
    integer :: ke

    real(RP), allocatable :: DENS(:)

    type(Experiment) :: exp_manager
    !------------------------------------------

    if ( .not. is_Qheat_calculated ) then
      call exp_manager%Init( 'equatorial_wave_global' )
      call exp_manager%Regist_SetInitCond( exp_SetInitCond_equatorial_wave )
      call exp_manager%SetInitCond( atm%mesh, &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager, atm%vars%QTRCVARS_manager )
      call exp_manager%Final()

      is_Qheat_calculated = .true.
    end if

    call FILE_HISTORY_meshfield_in( q_heat, "heating source" )

    tsec = atm%time_manager%dtsec * real( TIME_NOWSTEP - 1, kind=RP )
    fac = 1.0_RP!exp(-tsec * rtau2)

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      

      call AtmosVars_GetLocalMeshPhyAuxVars( n, atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                         )
      
      allocate( DENS(lcmesh%refElem3D%Np) )
      !$omp parallel do private(DENS)
      do ke=lcmesh%NeS, lcmesh%NeE
        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)

        atm%vars%PHY_TEND(MOMX_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMX_p)%local(n)%val(:,ke)   &
          - rtau * MOMX%val(:,ke)
        atm%vars%PHY_TEND(MOMY_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMY_p)%local(n)%val(:,ke)   &
          - rtau * MOMY%val(:,ke)
        atm%vars%PHY_TEND(MOMZ_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMZ_p)%local(n)%val(:,ke)   &
          - rtau * MOMZ%val(:,ke)

        atm%vars%PHY_TEND(RHOH_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(RHOH_p)%local(n)%val(:,ke)   &
          + DENS(:) * ( q_heat%local(n)%val(:,ke) * fac                                                                         &
                      - rtau * CPtot%val(:,ke) * ( PRES%val(:,ke) / DENS(:) - PRES_hyd%val(:,ke) / DENS_hyd%val(:,ke) ) / Rdry  )
      end do

      deallocate( DENS )
    end do

    return
  end subroutine USER_calc_tendency

  !------
!OCL SERIAL
  subroutine exp_SetInitCond_equatorial_wave( this,                      &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,        &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      PRES00 => CONST_PRE00, &
      RPlanet => CONST_RADIUS
    
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constBVFreq
    use mod_mkinit_util, only: &
      mkinitutil_GalerkinProjection_global
    use mod_experiment, only: &
      TracerLocalMeshField_ptr
    
    implicit none

    class(Experiment), intent(inout) :: this
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
    
    real(RP) :: THETA0           = 300.0_RP
    real(RP) :: BruntVaisalaFreq = 0.01_RP  
    real(RP) :: Q0   = 5.0_RP / 86400.0_RP !< the maximum heating ratio (5 K / day)
    real(RP) :: DLon                       !< The halfwidths of the heating region in the zonal direction
    real(RP) :: DLat                       !< The halfwidths of the heating region in the meridional direction
    integer  :: nv     = 1
    real(RP) :: Zt
    namelist /PARAM_EXP/ &
      Q0, DLon, DLat,    &
      nv,                &
      BruntVaisalaFreq
    integer, parameter :: IntrpPolyOrder_h = 8
    integer, parameter :: IntrpPolyOrder_v = 8

    integer :: ke
    integer :: ierr
    !-----------------------------------------------------------------------------

    DLon = 30.0_RP / 180.0_RP * PI
    DLat = 10.0_RP / 180.0_RP * PI

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("EQUATORIAL_WAVE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("EQUATORIAL_WAVE_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    Zt = dom_zmax - dom_zmin

    call hydrostatic_calc_basicstate_constBVFreq( DENS_hyd, PRES_hyd, & ! (out)
      BruntVaisalaFreq, THETA0, PRES00,                               & ! (in)
      x, y, z, lcmesh, elem                                           ) ! (in)

    if ( USER_do ) then
      call mkinitutil_GalerkinProjection_global( q_heat%local(lcmesh%lcdomID)%val(:,:), & ! (out)
        func_qheat, IntrpPolyOrder_h, IntrpPolyOrder_v, lcmesh, elem, RPlanet           ) ! (in)
    end if

    return
  contains
!OCL SERIAL
    subroutine func_qheat( q_intrp,           &
      lon, lat, zlev, elem_intrp, rplanet_    )
      implicit none
      class(ElementBase3D), intent(in) :: elem_intrp
      real(RP), intent(out) :: q_intrp(elem_intrp%Np)
      real(RP), intent(in) :: lon(elem_intrp%Np)
      real(RP), intent(in) :: lat(elem_intrp%Np)
      real(RP), intent(in) :: zlev(elem_intrp%Np)
      real(RP), intent(in) :: rplanet_

      real(RP) :: lon_(elem_intrp%Np)
      !------------------------------------------

      lon_(:) = PI - lon(:)
      where ( abs(lon_(:)) < DLon .and. abs(lat(:)) < DLat )
        q_intrp(:) = ( CpDry - Rdry ) * Q0 * cos(0.5_RP * PI * lon_(:) / DLon)**2  &
                                      * cos(0.5_RP * PI * lat (:) / DLat)**2       &
                                      * sin( real(nv,kind=RP) * PI * zlev(:) / Zt )
      elsewhere
        q_intrp(:) = 0.0_RP
      end where

      return
    end subroutine func_qheat
  end subroutine exp_SetInitCond_equatorial_wave

end module mod_user
