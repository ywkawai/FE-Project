!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of tracer advection
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

  !-----------------------------------------------------------------------------
contains
  subroutine USER_mkinit ( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('tracer_advection')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_tracer_advection )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

  subroutine USER_setup( this, atm )
    use scale_tracer, only: &
       TRACER_regist    
    class(User), intent(inout) :: this    
    class(AtmosComponent), intent(inout) :: atm

    logical :: USER_do        = .false. !< do user step?
    namelist / PARAM_USER / &
       USER_do

    integer :: ierr    
    integer :: iq        
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
    call TRACER_REGIST( iq,                   & ! [OUT]
                        1,                    & ! [IN]
                        (/'PTracer'/),        & ! [IN]
                        (/'Passive tracer'/), & ! [IN]
                        (/'1'/)               ) ! [IN]
        
    return
  end subroutine USER_setup

  !------
  subroutine exp_SetInitCond_tracer_advection( this,                       &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,        &
      GRAV => CONST_GRAV,    &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      CVdry => CONST_CVdry,  &
      PRES00 => CONST_PRE00, &
      Pstd   => CONST_Pstd
    use scale_tracer, only: &
      TRACER_inq_id

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
    
    ! Surface state
    real(RP)               :: SFC_THETA = 300.0_RP ! surface potential temperature [K]
    real(RP)               :: SFC_PRES             ! surface pressure [Pa]
    ! Environment state
    real(RP)               :: ENV_THETA = 300.0_RP ! potential temperature of environment [K]
    real(RP)               :: ENV_U     =   0.0_RP ! velocity u of environment [m/s]
    real(RP)               :: ENV_V     =   0.0_RP ! velocity v of environment [m/s]
    ! Bubble
    character(len=H_SHORT) :: SHAPE_PTracer = 'BUBBLE' ! BUBBLE or RECT
    real(RP)               :: BBL_PTracer   = 1.0_RP   ! extremum of passive tracer in bubble [kg/kg]
    real(RP)               :: x_c, y_c, z_c
    real(RP)               :: r_x, r_y, r_z

    namelist /PARAM_EXP/ &
      SFC_THETA,     &
      SFC_PRES,      &
      ENV_THETA,     &
      ENV_U,         &
      ENV_V,         &
      SHAPE_PTracer, &
      BBL_PTracer,   &
      x_c, y_c, z_c, &
      r_x, r_y, r_z          

    integer, parameter :: IntrpPolyOrder_h = 8
    integer, parameter :: IntrpPolyOrder_v = 8

    integer :: ke
    integer :: ierr
    integer :: iq
    real(RP) :: TEMP0
    !-----------------------------------------------------------------------------

    x_c = 0.5_RP * (dom_xmax + dom_xmin)
    y_c = 0.5_RP * (dom_ymax + dom_ymin)
    z_c = 0.5_RP * (dom_zmax + dom_zmin)

    r_x = 0.1_RP * (dom_xmax - dom_xmin)
    r_y = 0.1_RP * (dom_ymax - dom_ymin)
    r_z = 0.1_RP * (dom_zmax - dom_zmin)

    SFC_PRES = Pstd

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("TRACER_ADVECTION_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("TRACER_ADVECTION_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    
    TEMP0 = SFC_THETA * ( SFC_PRES / PRES00 )**(Rdry/CpDry)      
    call hydrostatic_calc_basicstate_constPT( DENS_hyd, PRES_hyd,                      &
      TEMP0, SFC_PRES, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3), &
      lcmesh, elem )

    !---

    !$omp parallel do
    do ke=lcmesh%NeS, lcmesh%NeE
      MOMX(:,ke) = ENV_U * DENS_hyd(:,ke)
      MOMY(:,ke) = ENV_V * DENS_hyd(:,ke)
    end do

    !---

    call TRACER_inq_id( "PTracer", iq )
        
    call mkinitutil_calc_cosinebell( &
      tracer_field_list(iq)%ptr%val(:,:),    & ! (out)
      BBL_PTracer,                           & ! (in)
      r_x, r_y, r_z, x_c, y_c, z_c,          & ! (in)
      x, y, z, lcmesh, elem,                 & ! (in)
      IntrpPolyOrder_h, IntrpPolyOrder_v     ) ! (in)  
        
    return
  end subroutine exp_SetInitCond_tracer_advection

end module mod_user
