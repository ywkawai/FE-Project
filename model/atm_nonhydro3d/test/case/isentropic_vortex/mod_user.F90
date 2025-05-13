!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of isentropic vortex
!!
!! @author Yuta Kawai, Team SCALE
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


  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( this, atm )
    implicit none

    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('isentropic_vortex')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_isentropic_vortex )
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

    call this%UserBase%Setup( atm, USER_do )
    !-

    return
  end subroutine USER_setup

  !------
!OCL SERIAL
  subroutine exp_SetInitCond_isentropic_vortex( this,                      &
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
    
    real(RP) :: BETA
    real(RP) :: x_c, y_c
    real(RP) :: U0, V0

    namelist /PARAM_EXP/ &
      BETA, &
      x_c, y_c, &
      U0, V0

    integer, parameter :: IntrpPolyOrder_h = 6
    integer, parameter :: IntrpPolyOrder_v = 8
    real(RP), allocatable :: THETA_purtub(:,:)
    
    real(RP) :: CpOvCv, CvOvCp
    real(RP) :: rP0
    real(RP) :: DENS(elem%Np)
    real(RP) :: R2(elem%Np)

    integer :: ke
    integer :: ierr
    !-----------------------------------------------------------------------------

    BETA = 5.0_RP
    U0   = 5.0_RP
    V0   = 0.0_RP
    x_c  = 0.5_RP * ( dom_xmax - dom_xmin  )
    y_c  = 0.5_RP * ( dom_ymax - dom_ymin  )

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("isentropic_vortex_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("isentropic_vortex_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    CpOvCv = CPdry / CVdry
    CvOvCp = CVdry / CPdry
    rP0    = 1.0_RP / PRES00

    !$omp parallel do private(DENS, R2)
    do ke=1, lcmesh%Ne

      DENS_hyd(:,ke) = 1.0_RP
      PRES_hyd(:,ke) = DENS_hyd(:,ke)**CpOvCv

      R2(:) = (x(:,ke) - x_c)**2 + (y(:,ke) - y_c)**2

      DENS(:) = ( 1.0_RP - BETA**2 * ( CpOvCv - 1.0_RP ) / ( 16.0_RP * CpOvCv * PI**2 ) &
                           * exp( 2.0_RP * ( 1.0_RP - R2(:) ) )                         &
                )**(CvDry/Rdry)
      
      DDENS(:,ke) = DENS(:) - DENS_hyd(:,ke)
      MOMX(:,ke) = DENS(:) * ( U0 - BETA * ( y(:,ke) - y_c ) / ( 2.0_RP * PI ) * exp( 1.0_RP - R2(:) ) )
      MOMY(:,ke) = DENS(:) * ( V0 + BETA * ( x(:,ke) - x_c ) / ( 2.0_RP * PI ) * exp( 1.0_RP - R2(:) ) )

      DRHOT(:,ke) = PRES00 / Rdry * ( ( DENS(:)**CpOvCv * rP0 )**CvOvCp &
                                    - ( PRES_hyd(:,ke)  * rP0 )**CvOvCp )
    end do

    return
  end subroutine exp_SetInitCond_isentropic_vortex

end module mod_user
