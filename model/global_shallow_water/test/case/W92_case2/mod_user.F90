!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!          Test case 2 of Williamson et al. (1992)
!!          Steady-state geostrophic flow
!!
!! @par Reference
!!  - Williamson et al. 1992: 
!!    A Standard Test Set for Numerical Approximations to the Shallow Water Equations in Spherical Geometry.
!!    Journal of Computational Physics, 102, 211-224.
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

  use mod_globalsw_component, only: &
    GlobalSWComponent

  use scale_element_base, only: ElementBase2D
  use scale_localmesh_2d, only: LocalMesh2D   

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

  type, private, extends(experiment) :: Exp_W92_case2
  contains 
   procedure :: setInitCond_lc => exp_SetInitCond_W92_case2
  end type
  type(Exp_W92_case2), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  subroutine USER_mkinit ( swmodel )
    implicit none

    class(GlobalSWComponent), intent(inout) :: swmodel
    !------------------------------------------

    call exp_manager%Init('W92_case2')
    call exp_manager%SetInitCond( &
      swmodel%mesh, swmodel%vars%PROGVARS_manager, swmodel%vars%AUXVARS_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

  subroutine USER_setup( atm )
    implicit none
    
    class(GlobalSWComponent), intent(inout) :: atm

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
    class(GlobalSWComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_calc_tendency

  subroutine USER_update( atm )
    implicit none
    class(GlobalSWComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_update

  !------
!OCL SERIAL
  subroutine exp_SetInitCond_W92_case2( this,  &
    h, U, V, hs, u1, u2,                       &
    x, y, lcmesh, elem                         )
    
    use scale_const, only: &
      PI => CONST_PI,          &
      GRAV => CONST_GRAV,      &
      RPlanet => CONST_RADIUS, &
      OMG => CONST_OHM 
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSVec
    implicit none

    class(Exp_W92_case2), intent(inout) :: this
    type(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: h(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: U(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: V(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: hs(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: u1(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: u2(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)  

    real(RP), parameter :: PHI0 = 2.94E4_RP
    real(RP) :: U0
    real(RP) :: alpha
    namelist /PARAM_EXP/ &
      alpha

    integer, parameter :: IntrpPolyOrder_h = 6
    integer, parameter :: IntrpPolyOrder_v = 6
    real(RP), allocatable :: PRES_purtub(:,:)
  
    real(RP) :: gam(elem%Np,lcmesh%Ne)
    real(RP) :: VelLon(elem%Np,lcmesh%Ne)
    real(RP) :: VelLat(elem%Np,lcmesh%Ne)
    integer :: ke

    integer :: ierr
    !-----------------------------------------------------------------------------

    U0 = 2.0_RP * PI * RPlanet / ( 12.0_RP * 2.0_RP * PI / OMG )

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("W92Case2_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("W92Case2_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    !$omp parallel do
    do ke=lcmesh%NeS, lcmesh%NeE
      h(:,ke) = ( PHI0                                                                                              &
          - 0.5_RP * U0  * ( 2.0_RP * RPlanet * OMG + U0 )                                                          &
          * ( sin(lcmesh%lat(:,ke)) * cos(alpha) - cos(lcmesh%lon(:,ke)) * cos(lcmesh%lat(:,ke)) * sin(alpha) )**2  &
        ) / Grav
      hs(:,ke) = 0.0_RP
      VelLon(:,ke) =   U0 * ( cos(alpha) * cos(lcmesh%lat(:,ke)) + sin(alpha) * cos(lcmesh%lon(:,ke)) * sin(lcmesh%lat(:,ke)) )
      VelLat(:,ke) = - U0 * sin(alpha) * sin(lcmesh%lon(:,ke))
    end do

    gam(:,:) = 1.0_RP
    call CubedSphereCoordCnv_LonLat2CSVec( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), & ! (in)
      gam(:,:), elem%Np * lcmesh%Ne,                              & ! (in)
      VelLon(:,:), VelLat(:,:),                                   & ! (in)
      U(:,lcmesh%NeS:lcmesh%NeE), V(:,lcmesh%NeS:lcmesh%NeE)      ) ! (out)

    !$omp parallel do
    do ke=lcmesh%NeS, lcmesh%NeE
      u1(:,ke) = lcmesh%G_ij(:,ke,1,1) * U(:,ke) + lcmesh%G_ij(:,ke,1,2) * V(:,ke)
      u2(:,ke) = lcmesh%G_ij(:,ke,2,1) * U(:,ke) + lcmesh%G_ij(:,ke,2,2) * V(:,ke)
    end do

    return
  end subroutine exp_SetInitCond_W92_case2
end module mod_user
