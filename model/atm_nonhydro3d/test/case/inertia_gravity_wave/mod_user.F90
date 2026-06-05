!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of inertia gravity wave
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
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

  type(MeshField3D), private :: PRES_diff

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('inertia_gravity_wave')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_inertia_gravity_wave )
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
    if ( USER_do ) call PRES_diff%Init( 'PRES_diff', 'Pa', atm%mesh%ptr_mesh )

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
      call atm%vars%Calc_diagVar( 'PRES_diff', PRES_diff )
      call FILE_HISTORY_meshfield_in( PRES_diff, "perturbation of PRES" )
    end if

    return
  end subroutine USER_calc_tendency

  !------

!OCL SERIAL  
  subroutine exp_SetInitCond_inertia_gravity_wave( this,                   &
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
      mkinitutil_GalerkinProjection
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
    real(RP) :: BruntVaisalaFreq = 1.0E-2_RP ! [s-1]

    real(RP) :: x_c, y_c
    real(RP) :: r_x, r_y

    integer :: IntrpPolyOrder_h
    integer :: IntrpPolyOrder_v

    namelist /PARAM_EXP/ &
      BruntVaisalaFreq,         &
      IntrpPolyOrder_h,         &
      IntrpPolyOrder_v,         &
      THETA0, DTHETA,           &
      x_c, y_c,                 &
      r_x, r_y

    real(RP) :: PT(elem%Np,lcmesh%NeA)
    
    real(RP) :: RovCp
    real(RP) :: DENS

    integer :: ke, p
    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 0.0_RP; y_c = 0.0_RP
    r_x = 5.0E3_RP; r_y = 5.0E3_RP
    THETA0  = 300.0_RP
    DTHETA  = 0.01_RP

    IntrpPolyOrder_h = elem%PolyOrder_h
    IntrpPolyOrder_v = elem%PolyOrder_v

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("INERTIA_GRAVITY_WAVE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("INERTIA_GRAVITY_WAVE_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    !$acc data create( PT )

    call hydrostatic_calc_basicstate_constBVFreq( DENS_hyd, PRES_hyd, & ! (out)
      BruntVaisalaFreq, THETA0, PRES00, x, y, z, lcmesh, elem         ) ! (in)

    call mkinitutil_GalerkinProjection( PT,                     & ! (out)
      calc_PT, IntrpPolyOrder_h, IntrpPolyOrder_v, lcmesh, elem ) ! (in)
    
    !---
    RovCp = Rdry / CpDry

    !$omp parallel do private( DENS )
    !$acc parallel loop gang present(DENS_hyd, PRES_hyd, PT, DDENS, x, y, z, lcmesh, elem)
    do ke=lcmesh%NeS, lcmesh%NeE
     !$acc loop vector
      do p=1, elem%Np        
        DENS = PRES_hyd(p,ke) / ( Rdry * PT(p,ke) * (PRES_hyd(p,ke)/PRES00)**(RovCp) )
        DDENS(p,ke) = DENS - DENS_hyd(p,ke)
      end do
    end do

    !$acc end data
    return
  contains
!OCL SERIAL
    subroutine calc_PT( PT_intrp, &
        x_, y_, z_, lcmesh3D, elem_intrp   )
      implicit none
      class(LocalMesh3D), intent(in) :: lcmesh3D
      class(ElementBase3D), intent(in) :: elem_intrp
      real(RP), intent(out) :: PT_intrp(elem_intrp%Np,lcmesh3D%Ne)
      real(RP), intent(in) :: x_(elem_intrp%Np,lcmesh3D%Ne)
      real(RP), intent(in) :: y_(elem_intrp%Np,lcmesh3D%Ne)
      real(RP), intent(in) :: z_(elem_intrp%Np,lcmesh3D%Ne)

      integer :: ke, p
      !-------------------------------------------
      !$omp parallel do
      !$acc parallel loop gang present(x_, y_, z_, lcmesh3D, elem_intrp)
      do ke=lcmesh3D%NeS, lcmesh3D%NeE
        !$acc loop vector
        do p=1, elem_intrp%Np
          PT_intrp(p,ke) = THETA0 * exp( BruntVaisalaFreq**2 / Grav * z_(p,ke) )              &
                         + DTHETA * sin( PI * z_(p,ke) / (dom_zmax - dom_zmin) )              &
                           / ( 1.0_RP + ((x_(p,ke) - x_c)/r_x)**2 + ((y_(p,ke) - y_c)/r_y)**2 )
        end do
      end do
    end subroutine calc_PT
  end subroutine exp_SetInitCond_inertia_gravity_wave

end module mod_user
