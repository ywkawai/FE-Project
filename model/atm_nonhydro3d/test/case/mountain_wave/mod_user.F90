!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for mountain wave experiments(linear hydrostatic case, linear nonhydrostatic case, and Schaer case)
!!
!! @author Yuta Kawai, Team SCALE
!!
!! @par Reference
!!  - Schär, C., Leuenberger, D., Fuhrer, O., Lüthi, D., and Girard, C. 2002: 
!!    A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models. 
!!    Mon. Wea. Rev., 130, 2459–2480
!!  - Giraldo, F. X., and Restelli, M. 2008:
!!    A study of spectral element and discontinuous Galerkin methods for the Navier–Stokes equations in nonhydrostatic mesoscale atmospheric modeling: Equation sets and test cases. 
!!    J. Comput. Phys., 227(8), 3849-3877.
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
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_mesh_base, only: MeshBase
  use scale_meshfield_base, only: MeshField3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D  

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
    procedure :: update => USER_update
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
  real(RP), private :: U0            
  real(RP), private :: zTop
  real(RP), private :: SPONGE_HEIGHT
  real(RP), private :: SPONGE_LATERAL_WIDTH
  real(RP), private :: SPONGE_EFOLD_SEC      = 600.0_RP
  real(RP), private :: Lx

  type(MeshField3D), private :: PT_diff

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit ( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init( 'mountain_wave' )
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_mountain_wave )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL  
  subroutine USER_setup( this, atm )
    implicit none
    class(User), intent(inout) :: this    
    class(AtmosComponent), intent(inout) :: atm

    logical :: USER_do        = .false. !< do user step?
    namelist / PARAM_USER / &
       USER_do,              &
       U0,                   &
       zTop,                 &
       SPONGE_HEIGHT,        &
       SPONGE_LATERAL_WIDTH, &
       SPONGE_EFOLD_SEC

    class(MeshBase), pointer :: ptr_mesh
    class(MeshCubeDom3D), pointer :: mesh

    integer :: ierr    
    !------------------------------------------


    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    U0                   = 10.0_RP
    zTop                 = 30.E3_RP
    SPONGE_HEIGHT        = 15.E3_RP
    SPONGE_EFOLD_SEC     = 900.0_RP
    SPONGE_LATERAL_WIDTH = 0.0_RP

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
    if ( this%USER_do ) call PT_diff%Init( 'PT_diff', 'K', atm%mesh%ptr_mesh )

    call atm%mesh_rm%GetModelMesh( ptr_mesh )
    select type( ptr_mesh )
    type is (MeshCubeDom3D)
      mesh => ptr_mesh
    end select

    Lx = mesh%xmax_gl

    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_update( this, atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in

    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CpDry => CONST_CPdry,    &
      RPlanet => CONST_RADIUS, &
      PI => CONST_PI

    use scale_localmeshfield_base, only: LocalMeshFieldBase
  
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      MOMX_p  => PHYTEND_MOMX_ID, &
      MOMY_p  => PHYTEND_MOMY_ID, &
      MOMZ_p  => PHYTEND_MOMZ_ID, &
      RHOH_p  => PHYTEND_RHOH_ID

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars

    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(ElementBase3D), pointer :: elem
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot

    real(RP) :: rtau

    integer :: n
    integer :: ke

    real(RP), allocatable :: DENS(:)
    real(RP), allocatable :: sfac(:)
    real(RP), allocatable :: rsfac(:)

    real(DP) :: dt
    real(RP) :: sponge_lateral_x00
    real(RP) :: sponge_lateral_x0
    !------------------------------------------
    
    rtau = 1.0_RP / SPONGE_EFOLD_SEC
    dt = atm%time_manager%dtsec
    sponge_lateral_x00 = Lx - SPONGE_LATERAL_WIDTH
    sponge_lateral_x0 = sponge_lateral_x00 + 0.5_RP * SPONGE_LATERAL_WIDTH
    
    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      
      
      call AtmosVars_GetLocalMeshPhyAuxVars( n,  atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                          )
      
      elem => lcmesh%refElem3D
      allocate( DENS(elem%Np), sfac(elem%Np), rsfac(elem%Np) )

      !$omp parallel do private(DENS, sfac, rsfac)
      do ke=lcmesh%NeS, lcmesh%NeE
        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)

        sfac(:) = 0.0_RP 
        where ( lcmesh%pos_en(:,ke,1) > sponge_lateral_x00 )
          sfac(:) = dt * rtau * 0.5_RP * ( 1.0_RP - cos( 2.0_RP * PI * ( ( lcmesh%pos_en(:,ke,1) - sponge_lateral_x00 ) / SPONGE_LATERAL_WIDTH ) ) )
        end where        
        where ( lcmesh%pos_en(:,ke,3) > SPONGE_HEIGHT )
          sfac(:) = sfac(:) + dt * rtau * 0.5_RP * ( 1.0_RP - cos( PI * ( lcmesh%pos_en(:,ke,3) - SPONGE_HEIGHT ) / ( zTop - SPONGE_HEIGHT ) ) )
        end where

        rsfac(:) = 1.0_RP / ( 1.0_RP + sfac(:) )

        MOMX%val(:,ke) = ( MOMX%val(:,ke) + sfac(:) * DENS(:) * U0 ) * rsfac(:)
        MOMY%val(:,ke) = MOMY%val(:,ke) * rsfac(:)
        MOMZ%val(:,ke) = MOMZ%val(:,ke) * rsfac(:)
        DRHOT%val(:,ke) = DRHOT%val(:,ke) * rsfac(:)
      end do
      deallocate( DENS, sfac, rsfac )
    end do

    return
  end subroutine USER_update

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
      call FILE_HISTORY_meshfield_in( PT_diff, "perturbation of potential temperature" )
    end if
    
    return
  end subroutine USER_calc_tendency

  !------

!OCL SERIAL  
  subroutine exp_SetInitCond_mountain_wave( this,                 &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,        &
      GRAV => CONST_GRAV,    &
      OHM => CONST_OHM,      &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      PRES00 => CONST_PRE00, &
      RPlanet => CONST_RADIUS
    
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constBVFreq
    
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
    
    real(RP) :: THETA0 = 280.0_RP
    real(RP) :: BruntVaisalaFreq = 0.01_RP  
    namelist /PARAM_EXP/ &
      U0,                &
      BruntVaisalaFreq,  &
      THETA0
    integer, parameter :: IntrpPolyOrder_h = 8
    integer, parameter :: IntrpPolyOrder_v = 8

    integer :: ke
    integer :: ierr

    !-----------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MOUNTAIN_WAVE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MOUNTAIN_WAVE_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)
    !---

    call hydrostatic_calc_basicstate_constBVFreq( DENS_hyd, PRES_hyd, & ! (out)
      BruntVaisalaFreq, THETA0, PRES00,                               & ! (in)
      x, y, lcmesh%zlev, lcmesh, elem                                 ) ! (in)
    !---

    !$omp parallel do
    do ke=lcmesh%NeS, lcmesh%NeE
      MOMX(:,ke) = DENS_hyd(:,ke) * U0
      MOMY(:,ke) = 0.0_RP
    end do

    return
  end subroutine exp_SetInitCond_mountain_wave

end module mod_user
