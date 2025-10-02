!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of density current based on Straka et al. (1993)
!!
!! @author Yuta Kawai, Team SCALE
!!
!! @par Reference
!!  - Straka, J.M., Wilhelmson, R.B., Wicker, L.J., Anderson, J.R., Droegemeier, K.K. 1993: 
!!    Numerical solutions of a non-linear density current: a benchmark solution and comparison. 
!!    Int. J. Numer. Methods Fluids., 17, 1–22. 
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

  type(MeshField3D), private :: PT_diff

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init( 'density_current' )
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_density_current )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL
  subroutine USER_setup( this, atm )
    implicit none
    
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    logical :: USER_do                   = .false. !< do user step?
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
    call this%UserBase%Setup( atm, USER_do )

    !-
    if ( USER_do ) call PT_diff%Init( 'PT_diff', 'K', atm%mesh%ptr_mesh )

    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    implicit none

    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    if ( this%USER_do ) then
      call atm%vars%Calc_diagVar( 'PT_diff', PT_diff )
      call FILE_HISTORY_meshfield_in( PT_diff, "perturbation of PT" )
    end if

    return
  end subroutine USER_calc_tendency
  
  !------
!OCL SERIAL
  subroutine exp_SetInitCond_density_current( this,                      &
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
    
    real(RP) :: THETA0
    real(RP) :: DTHETA
    real(RP) :: x_c, y_c, z_c
    real(RP) :: r_x, r_y, r_z

    integer :: IntrpPolyOrder_h
    integer :: IntrpPolyOrder_v

    namelist /PARAM_EXP/ &
      THETA0, DTHETA,    &
      x_c, y_c, z_c,     &
      r_x, r_y, r_z,     &
      IntrpPolyOrder_h,  &
      IntrpPolyOrder_v

    real(RP), allocatable :: THETA_purtub(:,:)
    
    real(RP) :: RovCp
    real(RP) :: PT  (elem%Np)
    real(RP) :: DENS(elem%Np)

    integer :: ke
    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 0.0_RP; y_c = 0.0_RP; z_c = 3.0E3_RP
    r_x = 4.0E3_RP; r_y = 4.0E3_RP; r_z = 2.0E3_RP
    THETA0    = 300.0_RP
    DTHETA    = -15.0_RP

    IntrpPolyOrder_h = elem%PolyOrder_h
    IntrpPolyOrder_v = elem%PolyOrder_v

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("DENSITY_CURRENT_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("DENSITY_CURRENT_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
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
    do ke=1, lcmesh%Ne
      PT(:) = THETA0 + THETA_purtub(:,ke)
      DENS(:) = PRES_hyd(:,ke) / ( Rdry * PT(:) * (PRES_hyd(:,ke)/PRES00)**(RovCp) )
      DDENS(:,ke) = DENS(:) - DENS_hyd(:,ke)
      DRHOT(:,ke) = DENS(:) * PT(:) - DENS_hyd(:,ke)*THETA0
    end do

    return
  end subroutine exp_SetInitCond_density_current

end module mod_user
