!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of Rayleigh–Bénard convection
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

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D   

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

  integer :: IntrpPolyOrder_h = 8
  integer :: IntrpPolyOrder_v = 8

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('Rayleigh--Bénard_convection')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_RB_convection )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL
  subroutine USER_setup( this, atm )
    implicit none
    
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    logical :: USER_do = .false. !< do user step?
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

    return
  end subroutine USER_setup


  !------
!OCL SERIAL  
  subroutine exp_SetInitCond_RB_convection( this,                   &
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
      hydrostatic_calc_basicstate_constBVFreq
    
    use mod_mkinit_util, only: &
      mkinitutil_gen_GPMat
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
    real(RP) :: BruntVaisalaFreq = 0.0E-2_RP ! [s-1]

    real(RP) :: x_c, y_c, z_c
    real(RP) :: r_x, r_y, r_z

    namelist /PARAM_EXP/ &
      BruntVaisalaFreq,         &
      THETA0, DTHETA,           &
      x_c, y_c, z_c,            &
      r_x, r_y, r_z,            &
      IntrpPolyOrder_h,         &
      IntrpPolyOrder_v

    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: IntrpMat(:,:)
    real(RP), allocatable :: x_intrp(:), y_intrp(:), z_intrp(:)
    real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)    
    
    real(RP) :: RovCp
    real(RP) :: PT  (elem%Np)
    real(RP) :: DENS(elem%Np)
    real(RP) :: exner_sfc
    real(RP) :: EXNER(elem%Np)

    integer :: ke
    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 0.0_RP; y_c = 0.0_RP; z_c = 5.E3_RP
    r_x = 5.0E3_RP; r_y = 5.0E3_RP; r_z = 1E2_RP

    THETA0  = 300.0_RP
    DTHETA  = 0.01_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("RB_convection_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("RB_convection_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, .false. )
    
    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    call mkinitutil_gen_GPMat( IntrpMat, elem_intrp, elem )

    allocate( x_intrp(elem_intrp%Np), y_intrp(elem_intrp%Np), z_intrp(elem_intrp%Np) )
    
    call hydrostatic_calc_basicstate_constBVFreq( DENS_hyd, PRES_hyd, & ! (out)
      BruntVaisalaFreq, THETA0, PRES00, x, y, z, lcmesh, elem         ) ! (in)
    
    RovCp = Rdry / CPdry
    
    exner_sfc = (PRES00 / PRES00)**RovCP

    !$omp parallel do private(PT, exner)
    do ke=lcmesh%NeS, lcmesh%NeE
      ! d exner / dz = - g / ( Cp * PT0 ) * exp (- N2/g * z)
      ! exner = exner(zs) - g^2 / (Cp * N^2) [ 1/PT (z) - 1/PT(zs) ] 
      PT(:) = THETA0
      exner(:) = exner_sfc - Grav / ( CpDry * THETA0 ) * lcmesh%zlev(:,ke)

      PRES_hyd(:,ke) = PRES00 * exner(:)**(CPdry/Rdry)
      DENS_hyd(:,ke) =  PRES_hyd(:,ke) / ( Rdry * exner(:) * PT(:) )
    end do
    
    !---

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
              + DTHETA                   &
                / ( 1.0_RP + ((x_intrp(:) - x_c)/r_x)**2 + ((y_intrp(:) - y_c)/r_y)**2   + ((z_intrp(:) - z_c)/r_z)**2 ) )
      
      DENS(:) = PRES_hyd(:,ke) / ( Rdry * PT(:) * (PRES_hyd(:,ke)/PRES00)**(RovCp) )
      DDENS(:,ke) = DENS(:) - DENS_hyd(:,ke)
    end do

    call elem_intrp%Final()

    return
  end subroutine exp_SetInitCond_RB_convection

end module mod_user
