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
    procedure :: setInitCond_lc => exp_SetInitCond_density_current
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
    call exp_manager%Init('density_current')
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
  subroutine exp_SetInitCond_density_current( this,              &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMZ, DRHOT,                &
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
    
    real(RP) :: THETA0
    real(RP) :: DTHETA
    real(RP) :: x_c, z_c, r_x, r_z
    logical :: InitCond_GalerkinProjFlag

    namelist /PARAM_EXP/ &
      THETA0, DTHETA,            &
      x_c, z_c, r_x, r_z,        &
      InitCond_GalerkinProjFlag


    integer :: k
    real(RP) :: THETA(elem%Np), DENS(elem%Np), dens_zfunc(elem%Np), RHOT(elem%Np)
    real(RP) :: r(elem%Np)

    integer, parameter :: IntrpPolyOrder = 12
    type(QuadrilateralElement) :: elem_intrp
    real(RP), allocatable :: x_intrp(:), z_intrp(:)
    real(RP) :: vx(elem%Nv), vz(elem%Nv)
    real(RP), allocatable :: IntrpMat(:,:), InvV_intrp(:,:)
    integer :: p1, p2, p_, p_intrp

    real(RP), allocatable :: r_intrp(:)
    real(RP), allocatable :: THETA_intrp(:)
  
    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 0.0_RP; z_c = 3.0E3_RP
    r_x = 4.0E3_RP; r_z = 2.0E3_RP
    THETA0    = 300.0_RP
    DTHETA    = -15.0_RP
    InitCond_GalerkinProjFlag = .false.

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
    call elem_intrp%Init( IntrpPolyOrder, .false. )
    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    allocate( InvV_intrp(elem%Np,elem_intrp%Np) )
    allocate( x_intrp(elem_intrp%Np), z_intrp(elem_intrp%Np) )
  
    allocate( r_intrp(elem_intrp%Np) )
    allocate( THETA_intrp(elem_intrp%Np) )


    InvV_intrp(:,:) = 0.0_RP
    do p2=1, elem%PolyOrder+1
    do p1=1, elem%PolyOrder+1
      p_ = p1 + (p2-1)*(elem%PolyOrder + 1)
      p_intrp = p1 + (p2-1)*(elem_intrp%PolyOrder + 1)
      InvV_intrp(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do    
    IntrpMat(:,:) = matmul(elem%V, InvV_intrp)

    !----
    do k=1, lcmesh%Ne
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),1)
      vz(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),2)
      x_intrp(:) = vx(1) + 0.5_RP*(elem_intrp%x1(:) + 1.0_RP)*(vx(2) - vx(1))
      z_intrp(:) = vz(1) + 0.5_RP*(elem_intrp%x2(:) + 1.0_RP)*(vz(3) - vz(1))

      dens_zfunc(:) = (1.0_RP - Grav*z(:,k)/(CpDry*THETA0))**(CVdry/Rdry)
      DENS_hyd(:,k) = PRES00/(THETA0*Rdry) * dens_zfunc(:)
      PRES_hyd(:,k) = PRES00 * (Rdry*DENS_hyd(:,k)*THETA0/PRES00)**(CPdry/Cvdry)

      r(:) = min(1.0_RP, sqrt(((x(:,k) - x_c)/r_x)**2 + ((z(:,k) - z_c)/r_z)**2))
      r_intrp(:) = min(1.0_RP, sqrt(((x_intrp(:) - x_c)/r_x)**2 + ((z_intrp(:) - z_c)/r_z)**2))
      
      THETA_intrp(:) = THETA0                          &
         + DTHETA*0.5_RP*(1.0_RP + cos(PI*r_intrp(:))) &
           / (1.0_RP - Grav*z_intrp(:)/(CpDry*THETA0)) 
      !THETA(:) = THETA0 + DTHETA*0.5_RP*(1.0_RP + cos(PI*r(:)))
      THETA(:) = matmul(IntrpMat, THETA_intrp)

      DENS(:) = PRES00/(THETA(:)*Rdry) * dens_zfunc(:)
      DDENS(:,k) = DENS(:) - DENS_hyd(:,k)

      DRHOT(:,k) = DENS(:)*THETA(:) - DENS_hyd(:,k)*THETA0


      MOMX(:,k) = 0.0_RP
      MOMZ(:,k) = 0.0_RP
    end do
    
    call elem_intrp%Final()

    return
  end subroutine exp_SetInitCond_density_current
end module mod_user
