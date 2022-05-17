!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of rising warm bubble
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
  use mod_exp, only: experiment

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D  
  
  
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

  type, private, extends(experiment) :: Exp_rising_therm_bubble_dry
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_rising_therm_bubble
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_rising_therm_bubble_dry), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    call exp_manager%Init('rising_therm_bubble_dry')

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

    !-

    return
  end subroutine USER_setup

  subroutine USER_calc_tendency( atm )
    implicit none
    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_calc_tendency

  subroutine USER_update( atm )
    implicit none
    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_update

  !------

!OCL SERIAL  
  subroutine exp_SetInitCond_rising_therm_bubble( this,                    &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,       &
      GRAV => CONST_GRAV,   &
      Rdry => CONST_Rdry,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constPT
    
    use mod_mkinit_util, only: &
      mkinitutil_calc_cosinebell
    use mod_exp, only: &
      TracerLocalMeshField_ptr
    
    implicit none

    class(Exp_rising_therm_bubble_dry), intent(inout) :: this
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
    
    real(RP) :: THETA0
    real(RP) :: DTHETA
    real(RP) :: x_c, y_c, z_c
    real(RP) :: r_x, r_y, r_z

    namelist /PARAM_EXP/ &
      THETA0, DTHETA,            &
      x_c, y_c, z_c,            &
      r_x, r_y, r_z

    integer, parameter :: IntrpPolyOrder_h = 6
    integer, parameter :: IntrpPolyOrder_v = 8
    real(RP), allocatable :: THETA_purtub(:,:)
    
    real(RP) :: RovCp
    real(RP) :: PT  (elem%Np)
    real(RP) :: DENS(elem%Np)

    integer :: ke
    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 500.0_RP; y_c = 500.0_RP; z_c = 350.0_RP
    r_x = 250.0_RP; r_y = 250.0_RP; r_z = 250.0_RP;
    THETA0    = 300.0_RP
    DTHETA    = 0.5_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("RISING_THERMAL_BUBBL_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("RISING_THERMAL_BUBBL_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    allocate( THETA_purtub(elem%Np,lcmesh%NeA) )
    call mkinitutil_calc_cosinebell( &
      THETA_purtub,                          &
      DTHETA, r_x, r_y, r_z, x_c, y_c, z_c,  &
      x, y, z, lcmesh, elem,                 &
      IntrpPolyOrder_h, IntrpPolyOrder_v     )  
    
    call hydrostatic_calc_basicstate_constPT( DENS_hyd, PRES_hyd,                       &
      THETA0, PRES00, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3), &
      lcmesh, elem )
    
    !---
    RovCp = Rdry / CpDry

    !$omp parallel do private(PT, DENS)
    do ke=lcmesh%NeS, lcmesh%NeE
      PT(:) = THETA0 + THETA_purtub(:,ke)
      DENS(:) = PRES_hyd(:,ke) / ( Rdry * PT(:) * (PRES_hyd(:,ke)/PRES00)**(RovCp) )
      DDENS(:,ke) = DENS(:) - DENS_hyd(:,ke)
      DRHOT(:,ke) = DENS(:) * PT(:) - DENS_hyd(:,ke)*THETA0
    end do

    return
  end subroutine exp_SetInitCond_rising_therm_bubble

!OCL SERIAL
  subroutine exp_geostrophic_balance_correction( this,                              &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    lcmesh, elem )
    
    implicit none

    class(Exp_rising_therm_bubble_dry), intent(inout) :: this
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
