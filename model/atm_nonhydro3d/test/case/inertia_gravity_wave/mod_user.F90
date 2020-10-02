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
  subroutine USER_mkinit
    implicit none
    !------------------------------------------

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
    call exp_manager%Init('inertia_gravity_wave')
    call exp_manager%SetInitCond( &
      atm%mesh, atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager )
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
    
    real(RP) :: THETA0 = 300.0_RP
    real(RP) :: BruntVaisalaFreq = 1.0E-2_RP ! [s-1]
    real(RP) :: DTHETA = 0.01_RP
    real(RP) :: H0
    real(RP) :: x_c, y_c
    real(RP) :: r_x, r_y
    logical :: InitCond_GalerkinProjFlag = .false.

    namelist /PARAM_EXP/ &
      THETA0, DTHETA,           &
      x_c, y_c,                 &
      r_x, r_y,                 &
      InitCond_GalerkinProjFlag


    integer :: k
    real(RP) :: THETA(elem%Np), DENS(elem%Np), dens_zfunc(elem%Np), RHOT(elem%Np)
    real(RP) :: r(elem%Np)

    integer, parameter :: IntrpPolyOrder_h = 6
    integer, parameter :: IntrpPolyOrder_v = 8
    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: x_intrp(:), y_intrp(:), z_intrp(:)
    real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)
    real(RP), allocatable :: IntrpMat(:,:), InvV_intrp(:,:)
    real(RP), allocatable :: IntrpVm1Mat(:,:), InvV_intrpVm1(:,:)
    integer :: p1, p2, p3, p_, p_intrp

    real(RP), allocatable :: THETA_hyd_intrp(:)
    real(RP), allocatable :: DENS_intrp(:)
    real(RP), allocatable :: DTHETA_intrp(:)
    real(RP), allocatable :: RHOT_intrp(:)    
    real(RP), allocatable :: EXNER_hyd_intrp(:) 
    real(RP), allocatable :: DRHOT_intrp(:)
    real(RP), allocatable :: r_intrp(:)   
  
    integer :: ierr
    !-----------------------------------------------------------------------------


    InitCond_GalerkinProjFlag = .false.

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("exp_SetInitCond_inertia_gravity_wave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("exp_SetInitCond_inertia_gravity_wave",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, .false. )
    allocate( IntrpMat(elem%Np,elem_intrp%Np), InvV_intrp(elem%Np,elem_intrp%Np) )
    allocate( IntrpVm1Mat(elem%Np,elem_intrp%Np), InvV_intrpVm1(elem%Np,elem_intrp%Np) )

    allocate( x_intrp(elem_intrp%Np), y_intrp(elem_intrp%Np), z_intrp(elem_intrp%Np) )

    allocate( THETA_hyd_intrp(elem_intrp%Np) )
    allocate( EXNER_hyd_intrp(elem_intrp%Np) )

    allocate( DTHETA_intrp(elem_intrp%Np) )
    allocate( DENS_intrp(elem_intrp%Np) )
    allocate( RHOT_intrp(elem_intrp%Np) )

    allocate( DRHOT_intrp(elem_intrp%Np) )    

    InvV_intrp(:,:) = 0.0_RP
    do p3=1, elem%PolyOrder_v+1
    do p2=1, elem%PolyOrder_h+1
    do p1=1, elem%PolyOrder_h+1
      p_ = p1 + (p2-1)*(elem%PolyOrder_h + 1) + (p3-1)*(elem%PolyOrder_h + 1)**2
      p_intrp = p1 + (p2-1)*(elem_intrp%PolyOrder_h + 1) + (p3-1)*(elem_intrp%PolyOrder_h + 1)**2
      InvV_intrp(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do
    end do
    IntrpMat(:,:) = matmul(elem%V, InvV_intrp)

    InvV_intrpVm1(:,:) = 0.0_RP
    do p3=1, elem%PolyOrder_v
    do p2=1, elem%PolyOrder_h+1
    do p1=1, elem%PolyOrder_h+1
      p_ = p1 + (p2-1)*(elem%PolyOrder_h + 1) + (p3-1)*(elem%PolyOrder_h + 1)**2
      p_intrp = p1 + (p2-1)*(elem_intrp%PolyOrder_h + 1) + (p3-1)*(elem_intrp%PolyOrder_h + 1)**2
      InvV_intrpVm1(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do
    end do
    IntrpVm1Mat(:,:) = matmul(elem%V, InvV_intrpVm1)

    !----
    H0 = Rdry*THETA0/Grav

    do k=1, lcmesh%Ne
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),2)
      vz(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),3)
      x_intrp(:) = vx(1) + 0.5_RP*(elem_intrp%x1(:) + 1.0_RP)*(vx(2) - vx(1))
      y_intrp(:) = vy(1) + 0.5_RP*(elem_intrp%x2(:) + 1.0_RP)*(vy(4) - vy(1))
      z_intrp(:) = vz(1) + 0.5_RP*(elem_intrp%x3(:) + 1.0_RP)*(vz(5) - vz(1))

      THETA_hyd_intrp(:) = THETA0 * exp(BruntVaisalaFreq**2*z_intrp(:)/Grav)
      EXNER_hyd_intrp(:) = 1.0_RP + Grav**2/(CpDry*BruntVaisalaFreq**2) * (1.0_RP/THETA_hyd_intrp(:) - 1.0_RP/THETA0)

      PRES_hyd(:,k) = matmul( IntrpMat, PRES00 * (EXNER_hyd_intrp(:))**(CpDry/Rdry) )
      DENS_hyd(:,k) = matmul( IntrpVm1Mat, PRES00 / (Rdry * THETA_hyd_intrp(:)) * (EXNER_hyd_intrp(:))**(CvDry/Rdry) )

      DTHETA_intrp(:) = DTHETA * sin( PI*z_intrp(:)/10.0E3_RP) &
                        / ( 1.0_RP + ((x_intrp(:) - x_c)/r_x)**2 + ((y_intrp(:) - y_c)/r_y)**2 )

      DDENS(:,k) =  matmul(IntrpVm1Mat, PRES00 / (Rdry * (THETA_hyd_intrp(:) + DTHETA_intrp(:))) * (EXNER_hyd_intrp(:))**(CvDry/Rdry))  &
                  - DENS_hyd(:,k)

      MOMX(:,k) = 0.0_RP
      MOMY(:,k) = 0.0_RP      
      MOMZ(:,k) = 0.0_RP
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
