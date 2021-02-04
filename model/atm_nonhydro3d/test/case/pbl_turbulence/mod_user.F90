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

  type, private, extends(experiment) :: Exp_pbl_turblence
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_pbl_turblence
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_pbl_turblence), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  subroutine USER_mkinit( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    call exp_manager%Init('idealized_pbl_turbulence')
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
  subroutine exp_SetInitCond_pbl_turblence( this,                      &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,        &
      GRAV => CONST_GRAV,    &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      CVdry => CONST_CVdry,  &
      PRES00 => CONST_PRE00, &
      Pstd   => CONST_Pstd
    use scale_random, only: &
      RANDOM_uniform
    implicit none

    class(Exp_pbl_turblence), intent(inout) :: this
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

    real(RP) :: ENV_PRES_SFC    
    real(RP) :: ENV_U          = 5.0_RP
    real(RP) :: ENV_V          = 0.0_RP
    real(RP) :: ENV_THETA_SFC  = 298.0_RP 
    real(RP) :: ENV_THETA_LAPS = 4.0E-3_RP
    real(RP) :: RANDOM_THETA   = 1.0_RP
    real(RP) :: RANDOM_U       = 0.0_RP
    real(RP) :: RANDOM_V       = 0.0_RP        
    logical :: InitCond_GalerkinProjFlag = .false.

    namelist /PARAM_EXP/ &
      ENV_U,            &
      ENV_THETA_SFC,    &
      ENV_THETA_LAPS,   &
      ENV_PRES_SFC,     &
      RANDOM_THETA,     &
      InitCond_GalerkinProjFlag


    integer :: ke
    real(RP) :: EXNER_sfc
    real(RP) :: EXNER(elem%Np)
    real(RP) :: THETA(elem%Np), THETA0(elem%Np)
    real(RP) :: DENS(elem%Np)
    real(RP) :: rndm(elem%Np)  
    integer :: ierr
    !-----------------------------------------------------------------------------

    ENV_PRES_SFC = Pstd

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("exp_SetInitCond_densitycurrent",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("exp_SetInitCond_densitycurrent",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    EXNER_sfc = (Pstd/PRES00)**(Rdry/Cpdry)

    !$omp parallel do private(             &
    !$omp EXNER, THETA, THETA0, DENS, rndm )
    do ke=1, lcmesh%Ne
      THETA0(:) = ENV_THETA_SFC + ENV_THETA_LAPS * z(:,ke)       
      EXNER(:) = EXNER_sfc  &
               - Grav / (CpDry * ENV_THETA_LAPS ) * log(1.0_RP + ENV_THETA_LAPS / ENV_THETA_SFC * z(:,ke))
      PRES_hyd(:,ke) = PRES00 * EXNER(:)**(CpDry/Rdry)
      DENS_hyd(:,ke) = PRES_hyd(:,ke) / ( Rdry * EXNER(:) * THETA0(:) )


      call RANDOM_uniform( rndm )
      THETA(:) = THETA0(:) + (rndm(:) * 2.0_RP - 1.0_RP ) * RANDOM_THETA     
      
      DENS(:) = PRES_hyd(:,ke) / ( Rdry * EXNER(:) * THETA(:) )
      DDENS(:,ke) = DENS(:) - DENS_hyd(:,ke)
      DRHOT(:,ke) = DENS(:) * THETA(:) - DENS_hyd(:,ke) * THETA0(:)

      MOMX(:,ke) = DENS(:) * (ENV_U + (rndm(:) * 2.0_RP - 1.0_RP ) * RANDOM_U )
      MOMY(:,ke) = DENS(:) * (ENV_V + (rndm(:) * 2.0_RP - 1.0_RP ) * RANDOM_V )    
      MOMZ(:,ke) = 0.0_RP
    end do
    
    return
  end subroutine exp_SetInitCond_pbl_turblence

  subroutine exp_geostrophic_balance_correction( this,                              &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    lcmesh, elem )
    
    implicit none

    class(Exp_pbl_turblence), intent(inout) :: this
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
