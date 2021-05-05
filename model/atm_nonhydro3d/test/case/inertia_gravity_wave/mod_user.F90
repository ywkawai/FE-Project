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

  type, private, extends(experiment) :: Exp_inertia_gravity_wave
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_inertia_gravity_wave
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_inertia_gravity_wave), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  subroutine USER_mkinit( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm

    !------------------------------------------

    call exp_manager%Init('inertia_gravity_wave')
    call exp_manager%SetInitCond( &
      atm%mesh, atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager )
    call exp_manager%Final()
        
    return
  end subroutine USER_mkinit

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

  subroutine USER_calc_tendency
    implicit none
    !------------------------------------------

    return
  end subroutine USER_calc_tendency

  subroutine USER_update
    implicit none
    !------------------------------------------

    return
  end subroutine USER_update

  !------
  subroutine exp_SetInitCond_inertia_gravity_wave( this,                      &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,       &
      GRAV => CONST_GRAV,   &
      Rdry => CONST_Rdry,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constBVFreq
    use mod_mkinit_util, only: &
      mkinitutil_gen_GPMat

    implicit none

    class(Exp_inertia_gravity_wave), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_ymin, dom_ymax
    real(RP), intent(in) :: dom_zmin, dom_zmax
    
    real(RP) :: THETA0
    real(RP) :: DTHETA
    real(RP) :: BruntVaisalaFreq = 1.0E-2_RP ! [s-1]

    real(RP) :: x_c, y_c
    real(RP) :: r_x, r_y

    namelist /PARAM_EXP/ &
      BruntVaisalaFreq,         &
      THETA0, DTHETA,           &
      x_c, y_c,                 &
      r_x, r_y

    integer, parameter :: IntrpPolyOrder_h = 6
    integer, parameter :: IntrpPolyOrder_v = 6

    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: IntrpMat(:,:)
    real(RP), allocatable :: x_intrp(:), y_intrp(:), z_intrp(:)
    real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)    
    
    real(RP) :: RovCp
    real(RP) :: PT  (elem%Np)
    real(RP) :: DENS(elem%Np)

    integer :: ke
    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 0.0_RP; y_c = 0.0_RP
    r_x = 5.0E3_RP; r_y = 5.0E3_RP
    THETA0  = 300.0_RP
    DTHETA  = 0.01_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("INERTIA_GRAVITY_WAVE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("INERTIA_GRAVITY_WAVE_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, .false. )
    
    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    call mkinitutil_gen_GPMat( IntrpMat, elem_intrp, elem )

    allocate( x_intrp(elem_intrp%Np), y_intrp(elem_intrp%Np), z_intrp(elem_intrp%Np) )
    
    call hydrostatic_calc_basicstate_constBVFreq( DENS_hyd, PRES_hyd,     &
      BruntVaisalaFreq, THETA0, PRES00,                                   &
      lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3),   &
      lcmesh, elem )
    
    !---
    RovCp = Rdry / CpDry

    !$omp parallel do private( PT, DENS,        &
    !$omp vx, vy, vz, x_intrp, y_intrp, z_intrp )
    do ke=lcmesh%NeS, lcmesh%NeE
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
      vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)
      x_intrp(:) = vx(1) + 0.5_RP*(elem_intrp%x1(:) + 1.0_RP)*(vx(2) - vx(1))
      y_intrp(:) = vy(1) + 0.5_RP*(elem_intrp%x2(:) + 1.0_RP)*(vy(4) - vy(1))
      z_intrp(:) = vz(1) + 0.5_RP*(elem_intrp%x3(:) + 1.0_RP)*(vz(5) - vz(1))

      PT(:) = matmul( IntrpMat, &
                THETA0 * exp( BruntVaisalaFreq**2 / Grav * z_intrp(:) )                  &
              + DTHETA * sin( PI * z_intrp(:) / (dom_zmax - dom_zmin) )                  &
                / ( 1.0_RP + ((x_intrp(:) - x_c)/r_x)**2 + ((y_intrp(:) - y_c)/r_y)**2 ) )
      
      DENS(:) = PRES_hyd(:,ke) / ( Rdry * PT(:) * (PRES_hyd(:,ke)/PRES00)**(RovCp) )
      DDENS(:,ke) = DENS(:) - DENS_hyd(:,ke)
    end do

    call elem_intrp%Final()

    return
  end subroutine exp_SetInitCond_inertia_gravity_wave

  subroutine exp_geostrophic_balance_correction( this,                              &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    lcmesh, elem )
    
    implicit none

    class(Exp_inertia_gravity_wave), intent(inout) :: this
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
