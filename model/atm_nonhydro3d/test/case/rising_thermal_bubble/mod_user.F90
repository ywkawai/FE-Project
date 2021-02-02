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

  type, private, extends(experiment) :: Exp_rising_therm_bubble
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_rising_therm_bubble
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_rising_therm_bubble), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  subroutine USER_mkinit( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    call exp_manager%Init('rising_therm_bubble')
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
  subroutine exp_SetInitCond_rising_therm_bubble( this,                  &
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
    
    implicit none

    class(Exp_rising_therm_bubble), intent(inout) :: this
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
    real(RP) :: x_c, y_c, z_c
    real(RP) :: r_x, r_y, r_z
    logical :: InitCond_GalerkinProjFlag

    namelist /PARAM_EXP/ &
      THETA0, DTHETA,            &
      x_c, y_c, z_c,             &
      r_x, r_y, r_z,             &
      InitCond_GalerkinProjFlag


    integer :: ke
    real(RP) :: THETA(elem%Np), DENS(elem%Np), dens_zfunc(elem%Np), RHOT(elem%Np)

    integer, parameter :: IntrpPolyOrder = 10
    real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)
    real(RP), allocatable :: r(:)

    real(RP), allocatable :: IntrpMat(:,:), InvV_intrp(:,:)
    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: x_intrp(:), y_intrp(:), z_intrp(:)
    real(RP), allocatable :: r_intrp(:)
    real(RP), allocatable :: THETA_intrp(:)
  
    integer :: p1, p2, p3
    integer :: p_, p_intrp

    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 500.0_RP; y_c = 500.0_RP; z_c = 350.0_RP
    r_x = 250.0_RP; r_y = 250.0_RP; r_z = 250.0_RP;
    THETA0    = 300.0_RP
    DTHETA    = 0.5_RP
    InitCond_GalerkinProjFlag = .false.

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("exp_SetInitCond_risingwarmbubble",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("exp_SetInitCond_risingwarmbubble",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    call elem_intrp%Init( IntrpPolyOrder, IntrpPolyOrder, .false. )
    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    allocate( InvV_intrp(elem%Np,elem_intrp%Np) )
    allocate( x_intrp(elem_intrp%Np), y_intrp(elem_intrp%Np), z_intrp(elem_intrp%Np) )
  
    allocate( r(elem%Np) )
    allocate( r_intrp(elem_intrp%Np) )
    allocate( THETA_intrp(elem_intrp%Np) )


    InvV_intrp(:,:) = 0.0_RP
    do p3=1, elem%Nnode_v
    do p2=1, elem%Nnode_h1D
    do p1=1, elem%Nnode_h1D
      p_ = p1 + (p2-1)*elem%PolyOrder_h + (p3-1)*elem%PolyOrder_h**2
      p_intrp = p1 + (p2-1)*elem_intrp%PolyOrder_h + (p3-1)*elem_intrp%PolyOrder_h**2
      InvV_intrp(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do    
    end do    
    IntrpMat(:,:) = matmul(elem%V, InvV_intrp)

    !----
    do ke=1, lcmesh%Ne
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
      vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)
      x_intrp(:) = vx(1) + 0.5_RP*(elem_intrp%x1(:) + 1.0_RP)*(vx(2) - vx(1))
      y_intrp(:) = vy(1) + 0.5_RP*(elem_intrp%x2(:) + 1.0_RP)*(vy(3) - vy(1))
      z_intrp(:) = vz(1) + 0.5_RP*(elem_intrp%x3(:) + 1.0_RP)*(vz(5) - vz(1))

      dens_zfunc(:) = (1.0_RP - Grav*z(:,ke)/(CpDry*THETA0))**(CVdry/Rdry)
      DENS_hyd(:,ke) = PRES00/(THETA0*Rdry) * dens_zfunc(:)
      PRES_hyd(:,ke) = PRES00 * (Rdry*DENS_hyd(:,ke)*THETA0/PRES00)**(CPdry/Cvdry)

      if (InitCond_GalerkinProjFlag) then
        r_intrp(:) = min( sqrt(  ((x_intrp(:) - x_c)/r_x)**2           &
                               + ((y_intrp(:) - y_c)/r_y)**2           &
                               + ((z_intrp(:) - z_c)/r_z)**2 ), 1.0_RP )
        THETA_intrp(:) = THETA0 + DTHETA * 0.5_RP*( 1.0_RP + cos(PI*r_intrp(:)) )
        THETA(:) = matmul(IntrpMat, THETA_intrp)
      else
        r(:) = min( sqrt(  ((x(:,ke) - x_c)/r_x)**2           &
                         + ((y(:,ke) - y_c)/r_y)**2           &
                         + ((z(:,ke) - z_c)/r_z)**2 ), 1.0_RP )
        THETA(:) = THETA0 + DTHETA * 0.5_RP*( 1.0_RP + cos(PI * r(:)) ) 
      end if

      DENS(:) = PRES00/(THETA(:)*Rdry) * dens_zfunc(:)
      DDENS(:,ke) = DENS(:) - DENS_hyd(:,ke)

      DRHOT(:,ke) = DENS(:)*THETA(:) - DENS_hyd(:,ke)*THETA0

      MOMX(:,ke) = 0.0_RP
      MOMY(:,ke) = 0.0_RP
      MOMZ(:,ke) = 0.0_RP
    end do
    
    call elem_intrp%Final()

    return
  end subroutine exp_SetInitCond_rising_therm_bubble

  subroutine exp_geostrophic_balance_correction( this,                              &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    lcmesh, elem )
    
    implicit none

    class(Exp_rising_therm_bubble), intent(inout) :: this
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
