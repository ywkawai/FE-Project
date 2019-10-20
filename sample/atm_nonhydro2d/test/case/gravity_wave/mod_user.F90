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

  type, private, extends(experiment) :: Exp_density_current
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_gravcurrent
  end type
  type(Exp_density_current), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  subroutine USER_mkinit
    implicit none
    !------------------------------------------

    return
  end subroutine USER_mkinit

  subroutine USER_setup
    implicit none

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
    call exp_manager%Init('gravity_wave')
    call exp_manager%SetInitCond()
    call exp_manager%Final()

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
  subroutine exp_SetInitCond_gravcurrent( this, &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMZ, DRHOT, &
    x, z, dom_xmin, dom_xmax, dom_zmin, dom_zmax, lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,       &
      GRAV => CONST_GRAV,   &
      Rdry => CONST_Rdry,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    
    use scale_element_base, only: ElementBase2D
    use scale_element_quadrilateral, only: QuadrilateralElement
    use scale_localmesh_2d, only: LocalMesh2D   
    implicit none

    class(Exp_density_current), intent(inout) :: this
    type(LocalMesh2D), intent(in) :: lcmesh
    class(QuadrilateralElement), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_zmin, dom_zmax
    
    real(RP) :: T0
    real(RP) :: DTemp
    real(RP) :: x_c, r_d
    real(RP) :: U0

    namelist /PARAM_EXP/ &
      T0, DTemp,        &
      x_c, r_d,         &
      U0

    integer :: k
    real(RP) :: DENS00
    real(RP) :: H0, Htop
    real(RP) :: Theta_hyd(elem%Np)
    real(RP) :: DENSb(elem%Np), Tb(elem%Np), DT(elem%Np), DENS(elem%Np)
    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 0.5_RP * (dom_xmin + dom_xmax)
    r_d = (dom_xmax - dom_xmin)/60.0_RP
    T0    = 250.0_RP
    DTemp = 0.01_RP
    U0 = 0.0_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("exp_SetInitCond_gravwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("exp_SetInitCond_gravwave",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    H0 = Rdry*T0/Grav
    Htop = dom_zmax - dom_zmin
    DENS00 = PRES00/(H0*Grav)
    do k=1, lcmesh%Ne
      
      DENS_hyd(:,k) = DENS00*exp(-z(:,k)/H0)
      PRES_hyd(:,k) = PRES00*exp(-z(:,k)/H0)
      Theta_hyd(:) = T0*exp(z(:,k)/H0)**(Rdry/CPdry)

      Tb(:) = DTemp * exp(- ((x(:,k) - x_c)/r_d)**2) * sin(PI*z(:,k)/Htop)
      DT(:) = Tb(:) * exp(0.5_RP*z(:,k)/H0)

      DENSb(:) = - PRES00*Tb(:)/(Rdry*T0**2)
      DDENS(:,k) = DENSb(:) * exp(- 0.5_RP*z(:,k)/H0)
      DENS(:) = DENS_hyd(:,k) + DDENS(:,k)

      MOMX(:,k) = DENS(:)*U0
      MOMZ(:,k) = 0.0_RP
      DRHOT(:,k) =   DENS(:) * (T0 + DT(:)) * (PRES00/PRES_hyd(:,k))**(Rdry/CPdry) &
                   - DENS_hyd(:,k)*Theta_hyd(:)
    end do

    return
  end subroutine exp_SetInitCond_gravcurrent

end module mod_user
