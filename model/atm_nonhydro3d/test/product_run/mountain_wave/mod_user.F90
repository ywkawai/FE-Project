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
  use scale_prc, only: &
    PRC_abort, PRC_myrank
  use scale_const, only: &
    PI => CONST_PI,        &
    GRAV => CONST_GRAV,    &
    OHM => CONST_OHM,      &
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  &
    CVdry => CONST_CVdry,  &
    PRES00 => CONST_PRE00, &
    RPlanet => CONST_RADIUS
  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D, ElementBase2D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_cubedom3d, only: MeshCubeDom3D  

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
    procedure :: update => USER_update
    procedure :: update_pre => USER_update_pre
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

  logical :: sponge_layer_flag = .false.
  integer, private :: SLFUNC_TYPEID
  integer, private, parameter :: SLFUNC_COSBELL_TYPEID = 1
  integer, private, parameter :: SLFUNC_TANH_TYPEID    = 2
  real(RP), private :: zTop
  real(RP), private :: SPONGE_HEIGHT
  real(RP), private :: SPONGE_LATERAL_WIDTH
  real(RP), private :: SPONGE_EFOLD_SEC = 600.0_RP
  logical :: lateral_sponge_layer_flag = .false.
  real(RP), private :: LATERAL_SPONGE_EFOLD_SEC = 600.0_RP
  real(RP), private :: SL_TANH_NONDIM_WIDTH     = 0.16_RP
  logical, private :: SL_APPLY_DENS = .false.

  logical :: ini_bg_force_flag = .false.
  real(RP) :: ini_bg_force_tscale = 10.0_RP
  real(RP) :: ini_bg_force_turnoff_tstart = 50.0_RP
  real(RP) :: ini_bg_force_turnoff_tscale = 100.0_RP

  type(MeshField3D), private :: PT_diff
  type(MeshField3D), private :: GsqrtDENS
  type(MeshField3D), private :: GsqrtDPRES
  type(MeshField3D), private :: GsqrtG13DPRES
  type(MeshField3D), private :: DENS_dt_1
  type(MeshField3D), private :: DENS_dt_2
  type(MeshField3D), private :: DENS_dt_3
  type(MeshField3D), private :: DENS_dt_4
  type(MeshField3D), private :: DENS_dt_5
  type(MeshField3D), private :: DENS_dt_6
  type(MeshField3D), private :: GsqrtMOMX
  type(MeshField3D), private :: GsqrtMOMZ
  type(MeshField3D), private :: GsqrtMOMW
  type(MeshField3D), private :: sfac_save

  real(RP), private :: Lx
  real(RP), allocatable :: sfac(:,:)
  real(RP), allocatable :: sfac_h(:,:)  
  logical :: is_PREShyd_ref_set

  real(RP) :: TEMP0 = 280.0_RP
  logical :: U_tappering_flag
  real(RP) :: U_TAPPERING_HEIGHT
  real(RP) :: U_TAPPERING_WIDTH    

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit ( this, atm )
    use scale_file_base_meshfield, only: &
      FILE_base_meshfield
    use scale_mesh_base3d, only: &
      MeshBase3D_DIMTYPEID_XY, &
      MeshBase3D_DIMTYPEID_XYZ
    use scale_mesh_base3d, only: &
      MeshBase3D
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    type(FILE_base_meshfield) :: writer
    logical :: fileexisted

    class(MeshBase3D), pointer :: mesh3D
    type(MeshField3D) :: Gsqrt, G13
    type(LocalMesh3D), pointer :: lcmesh3D
    integer :: n, ke

    integer :: date(6)
    !------------------------------------------


    call exp_manager%Init( 'mountain_wave' )
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_mountain_wave )
    call exp_manager%Regist_geostrophic_balance_correction( exp_geostrophic_balance_correction_lc )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    mesh3D => atm%mesh%ptr_mesh
    call Gsqrt%Init("Gsqrt", "1", mesh3D)
    call G13%Init("G13", "1", mesh3D)
    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(n)
      do ke=lcmesh3D%NeS, lcmesh3D%NeE
        Gsqrt%local(n)%val(:,ke) = lcmesh3D%Gsqrt(:,ke)
        G13%local(n)%val(:,ke) = lcmesh3D%GI3(:,ke,1)
      end do
    end do
    select type(mesh3D)
    class is (MeshCubeDom3D)
      call writer%Init(2, mesh3D=mesh3D)
    end select
    date(:) = -1
    call writer%Create( "metricInfo", "Metric Information", "REAL8", &
      fileexisted, myrank=PRC_myrank, tunits="s" )
    call writer%Put_GlobalAttribute_time( date, 0.0_RP )
    call writer%Def_Var( Gsqrt, "Gsqrt", 1, MeshBase3D_DIMTYPEID_XYZ, "XYZ" )
    call writer%Def_Var( G13, "G13", 2, MeshBase3D_DIMTYPEID_XYZ, "XYZ" )
    call writer%End_def()
    call writer%Write_var3D(1, Gsqrt, 0.0_RP, 1.0_RP)
    call writer%Write_var3D(2, G13, 0.0_RP, 1.0_RP)
    call writer%Close()
    call writer%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL  
  subroutine USER_setup( this, atm )
    implicit none
    class(User), intent(inout) :: this    
    class(AtmosComponent), intent(inout) :: atm

    logical :: USER_do        = .false. !< do user step?
    character(len=H_SHORT) :: SPONGE_LAYER_FUNC_NAME

    namelist / PARAM_USER /  &
       USER_do,              &
       sponge_layer_flag,    &
       U0,                   &
       zTop,                 &
       SPONGE_HEIGHT,        &
       SPONGE_LATERAL_WIDTH, &
       SPONGE_EFOLD_SEC,     &
       lateral_sponge_layer_flag, &
       LATERAL_SPONGE_EFOLD_SEC,  &
       SPONGE_LAYER_FUNC_NAME,    &
       SL_TANH_NONDIM_WIDTH,      &
       ini_bg_force_flag,         &
       ini_bg_force_tscale,       &
       ini_bg_force_turnoff_tstart, &
       ini_bg_force_turnoff_tscale, &
       SL_APPLY_DENS

    class(MeshBase), pointer :: ptr_mesh
    class(MeshCubeDom3D), pointer :: mesh

    integer :: ierr    
    !------------------------------------------


    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    U0             = 10.0_RP
    zTop           = 30.E3_RP
    SPONGE_HEIGHT  = 15.E3_RP
    SPONGE_EFOLD_SEC = 900.0_RP
    SPONGE_LATERAL_WIDTH = 0.0_RP

    SPONGE_LAYER_FUNC_NAME = "COSBELL"

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
    
    select case( trim(SPONGE_LAYER_FUNC_NAME) )
    case ( "COSBELL")
      SLFUNC_TYPEID = SLFUNC_COSBELL_TYPEID
    case ( "TANH")
      SLFUNC_TYPEID = SLFUNC_TANH_TYPEID
    case default
      LOG_ERROR("USER_setup",*) 'Not supported function for SPONGE_LAYER_FUNC. Check!', SPONGE_LAYER_FUNC_NAME
      call PRC_abort
    end select

    is_PREShyd_ref_set = .false.
    
    call sfac_save%Init( "sfac", "s-1", atm%mesh%ptr_mesh )
    call setup_sfac( this, atm )
    
    call GsqrtDENS%Init( "GsqrtDDENS", "kg.m-3", atm%mesh%ptr_mesh )
    call GsqrtDPRES%Init( "GsqrtDPRES", "kg.m-3", atm%mesh%ptr_mesh )
    call GsqrtG13DPRES%Init( "GsqrtG13DPRES", "kg.m-3", atm%mesh%ptr_mesh )
    call GsqrtMOMX%Init( "GsqrtMOMX", "kg.m-2.s-1", atm%mesh%ptr_mesh )
    call GsqrtMOMZ%Init( "GsqrtMOMZ", "kg.m-2.s-1", atm%mesh%ptr_mesh )
    call GsqrtMOMW%Init( "GsqrtMOMW", "kg.m-2.s-1", atm%mesh%ptr_mesh )
    call DENS_dt_1%Init( "DENS_dt_1", "kg.m-2.s-1", atm%mesh%ptr_mesh )
    call DENS_dt_2%Init( "DENS_dt_2", "kg.m-2.s-1", atm%mesh%ptr_mesh )
    call DENS_dt_3%Init( "DENS_dt_3", "kg.m-2.s-1", atm%mesh%ptr_mesh )
    call DENS_dt_4%Init( "DENS_dt_4", "kg.m-2.s-1", atm%mesh%ptr_mesh )
    call DENS_dt_5%Init( "DENS_dt_5", "kg.m-2.s-1", atm%mesh%ptr_mesh )
    call DENS_dt_6%Init( "DENS_dt_6", "kg.m-2.s-1", atm%mesh%ptr_mesh )

    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine setup_sfac( this, atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in

    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CpDry => CONST_CPdry,    &
      RPlanet => CONST_RADIUS, &
      PI => CONST_PI

    use scale_localmeshfield_base, only: LocalMeshFieldBase

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
    real(RP) :: rtau_sponge
    real(RP) :: rtau_lateral_sponge

    integer :: n
    integer :: ke

    real(RP) :: sponge_lateral_x00
    real(RP) :: sponge_lateral_x0
    
    !------------------------------------------
    
    rtau_sponge = 0.0_RP
    if ( sponge_layer_flag ) then
      rtau_sponge = 1.0_RP / SPONGE_EFOLD_SEC
    end if

    rtau_lateral_sponge = 0.0_RP
    if ( lateral_sponge_layer_flag ) then
      rtau_lateral_sponge = 1.0_RP / LATERAL_SPONGE_EFOLD_SEC
    end if

    sponge_lateral_x00 = Lx - SPONGE_LATERAL_WIDTH
    sponge_lateral_x0 = sponge_lateral_x00 + 0.5_RP * SPONGE_LATERAL_WIDTH
    
    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM      
      lcmesh => atm%mesh%ptr_mesh%lcmesh_list(n) 
      elem => lcmesh%refElem3D
      allocate( sfac(elem%Np,lcmesh%Ne), sfac_h(elem%Np,lcmesh%Ne) )

      !$omp parallel do
      do ke=lcmesh%NeS, lcmesh%NeE
        select case( SLFUNC_TYPEID )
        case ( SLFUNC_COSBELL_TYPEID )
    
          sfac(:,ke) = 0.0_RP 
          sfac_h(:,ke) = 0.0_RP 
          where ( lcmesh%pos_en(:,ke,1) > sponge_lateral_x00 )
            sfac_h(:,ke) = sfac_h(:,ke) + rtau_lateral_sponge * 0.5_RP * ( 1.0_RP - cos( 2.0_RP * PI * ( ( lcmesh%pos_en(:,ke,1) - sponge_lateral_x00 ) / SPONGE_LATERAL_WIDTH ) ) )
          end where        
          where ( lcmesh%zlev(:,ke) > SPONGE_HEIGHT )
            sfac(:,ke) = sfac(:,ke) + rtau_sponge * 0.5_RP * ( 1.0_RP - cos( PI * ( lcmesh%zlev(:,ke) - SPONGE_HEIGHT ) / ( zTop - SPONGE_HEIGHT ) ) )
          end where
        case ( SLFUNC_TANH_TYPEID )
          sfac(:,ke) = &
            + rtau_sponge * 0.5_RP * ( 1.0_RP + tanh( ( lcmesh%zlev(:,ke) - 0.5_RP * ( zTop + SPONGE_HEIGHT ) ) / ( SL_TANH_NONDIM_WIDTH * ( zTop - SPONGE_HEIGHT ) ) ) )
          sfac_h(:,ke) = &
              rtau_lateral_sponge * 0.5_RP * ( 1.0_RP - tanh( ( lcmesh%pos_en(:,ke,1) - 0.5_RP * SPONGE_LATERAL_WIDTH ) / ( SL_TANH_NONDIM_WIDTH * SPONGE_LATERAL_WIDTH ) ) ) &
            + rtau_lateral_sponge * 0.5_RP * ( 1.0_RP + tanh( ( lcmesh%pos_en(:,ke,1) -             sponge_lateral_x0 ) / ( SL_TANH_NONDIM_WIDTH * SPONGE_LATERAL_WIDTH ) ) )
        end select

        sfac_save%local(1)%val(:,ke) = sfac(:,ke)
      end do
    end do

    return
  end subroutine setup_sfac

!OCL SERIAL
  subroutine USER_update( this, atm )
    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CpDry => CONST_CPdry,    &
      RPlanet => CONST_RADIUS, &
      PI => CONST_PI

    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use scale_time_manager, only:  TIME_NOWSTEP

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

    integer :: n
    integer :: ke

    real(RP), allocatable :: DENS(:), rsfac(:)
    real(DP) :: dt_

    real(RP) :: tsec
    real(RP) :: U_bg
    real(RP) :: ini_bg_sfac, ini_bg_off_tsec
    real(RP) :: sw
    !------------------------------------------
        
    tsec = atm%time_manager%dtsec * ( real( TIME_NOWSTEP, kind=RP ) - 1.5_RP )

    dt_ = 0.5_RP * atm%time_manager%dtsec
    ini_bg_off_tsec = ini_bg_force_turnoff_tstart + ini_bg_force_turnoff_tscale
    if ( ini_bg_force_flag .and. tsec < ini_bg_off_tsec ) then
      U_bg = U0 * ( 1.0_RP - exp(-tsec/ini_bg_force_tscale) )
      ini_bg_sfac = 1.0_RP / ini_bg_force_tscale
      if ( tsec > ini_bg_force_turnoff_tstart ) then
        sw = 0.5_RP * ( 1.0_RP - cos( PI * ( ( tsec - ini_bg_force_turnoff_tstart ) / ini_bg_force_turnoff_tscale - 1.0_RP ) ) )
      else
        sw = 1.0_RP
      end if
    else
      U_bg = U0
      ini_bg_sfac = 0.0_RP
      sw = 0.0_RP
    end if
    ini_bg_sfac = 0.0_RP
    LOG_INFO("USER_up",*) "time=", tsec, U_bg, sw

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      
      
      call AtmosVars_GetLocalMeshPhyAuxVars( n,  atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                          )
      
      elem => lcmesh%refElem3D
      allocate( DENS(elem%Np), rsfac(elem%Np) )

      !$omp parallel do private(DENS, rsfac)
      do ke=lcmesh%NeS, lcmesh%NeE
        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)
        rsfac(:) = 1.0_RP / ( 1.0_RP + dt_ * ( sw * ini_bg_sfac + ( 1.0_RP - sw ) * ( sfac_h(:,ke) + sfac(:,ke) ) ) )

        ! MOMX%val(:,ke) = ( MOMX%val(:,ke) +  dt_ * ( sw * ini_bg_sfac + (1.0_RP - sw ) * ( sfac_h(:,ke) + sfac(:,ke) ) ) * DENS(:) * U0 ) * rsfac(:) 
        ! if ( sw > 0.0_RP .or. SL_APPLY_DENS ) then
        !   DDENS%val(:,ke) = DDENS%val(:,ke) * rsfac(:)          
        ! end if
        ! MOMY%val(:,ke) = MOMY%val(:,ke) * rsfac(:)
        ! MOMZ%val(:,ke) = MOMZ%val(:,ke) * rsfac(:)
        ! DRHOT%val(:,ke) = DRHOT%val(:,ke) * rsfac(:)
      end do
      deallocate( DENS, rsfac )
    end do

    return
  end subroutine USER_update

!OCL SERIAL
  subroutine USER_update_pre( this, atm )
    use scale_time_manager, only:  TIME_NOWSTEP    
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in

    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CpDry => CONST_CPdry,    &
      RPlanet => CONST_RADIUS, &
      PI => CONST_PI

    use scale_localmeshfield_base, only: LocalMeshFieldBase
  
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      MOMX_p  => PHYTEND_MOMX_ID, &
      MOMY_p  => PHYTEND_MOMY_ID, &
      MOMZ_p  => PHYTEND_MOMZ_ID, &
      RHOH_p  => PHYTEND_RHOH_ID, &
      PRESHYD_REF_VID => AUXVAR_PRESHYDRO_REF_ID, &
      PRESHYD_VID => AUXVAR_PRESHYDRO_ID

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars

    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(ElementBase3D), pointer :: elem
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd_
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot


    integer :: n
    integer :: ke

    class(MeshField3D), pointer :: PRES_hyd
    class(MeshField3D), pointer :: PRES_hyd_ref

    real(RP), allocatable :: DENS(:), rsfac(:)

    real(DP) :: dt, dt_
    real(RP) :: tsec
    real(RP) :: U_bg
    real(RP) :: ini_bg_sfac, ini_bg_off_tsec
    real(RP) :: sw
    !------------------------------------------

    tsec = atm%time_manager%dtsec * ( real( TIME_NOWSTEP, kind=RP ) - 2.0_RP )

    dt = atm%time_manager%dtsec
    dt_ = 0.5_RP * atm%time_manager%dtsec
    
    ini_bg_off_tsec = ini_bg_force_turnoff_tstart + ini_bg_force_turnoff_tscale
    if ( ini_bg_force_flag .and. tsec < ini_bg_off_tsec ) then
      U_bg = U0 * ( 1.0_RP - exp(-tsec/ini_bg_force_tscale) )
      ini_bg_sfac = 1.0_RP / ini_bg_force_tscale
      if ( tsec > ini_bg_force_turnoff_tstart ) then
        sw = 0.5_RP * ( 1.0_RP - cos( PI * ( ( tsec - ini_bg_force_turnoff_tstart ) / ini_bg_force_turnoff_tscale - 1.0_RP ) ) )
      else
        sw = 1.0_RP
      end if
    else
      U_bg = U0
      ini_bg_sfac = 0.0_RP
      sw = 0.0_RP
    end if
    ini_bg_sfac = 0.0_RP
    LOG_INFO("USER_up",*) "time=", tsec, U_bg, sw

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd_, Rtot, CVtot, CPtot, lcmesh          )      
      
      call AtmosVars_GetLocalMeshPhyAuxVars( n,  atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                          )
      
      elem => lcmesh%refElem3D
      allocate( DENS(elem%Np), rsfac(elem%Np) )

      !$omp parallel do private(DENS, rsfac)
      do ke=lcmesh%NeS, lcmesh%NeE
        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)
        rsfac(:) = 1.0_RP / ( 1.0_RP + dt_ * ( sw * ini_bg_sfac + ( 1.0_RP - sw ) * ( sfac_h(:,ke) + sfac(:,ke) ) ) )

        ! MOMX%val(:,ke) = ( MOMX%val(:,ke) +  dt_ * ( sw * ini_bg_sfac + (1.0_RP - sw ) * ( sfac_h(:,ke) + sfac(:,ke) ) ) * DENS(:) * U0 ) * rsfac(:) 
        ! if ( sw > 0.0_RP .or. SL_APPLY_DENS ) then
        !   DDENS%val(:,ke) = DDENS%val(:,ke) * rsfac(:)          
        ! end if
        ! MOMY%val(:,ke) = MOMY%val(:,ke) * rsfac(:)
        ! MOMZ%val(:,ke) = MOMZ%val(:,ke) * rsfac(:)
        ! DRHOT%val(:,ke) = DRHOT%val(:,ke) * rsfac(:)
      end do
      deallocate( DENS, rsfac )
    end do

    return
  end subroutine USER_update_pre

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in

    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CpDry => CONST_CPdry,    &
      RPlanet => CONST_RADIUS, &
      PI => CONST_PI

    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use scale_sparsemat, only: sparsemat, sparsemat_matmul
  
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      MOMX_p  => PHYTEND_MOMX_ID, &
      MOMY_p  => PHYTEND_MOMY_ID, &
      MOMZ_p  => PHYTEND_MOMZ_ID, &
      RHOH_p  => PHYTEND_RHOH_ID, &
      PRESHYD_REF_VID => AUXVAR_PRESHYDRO_REF_ID, &
      PRESHYD_VID => AUXVAR_PRESHYDRO_ID, &
      DDENS_VID => PRGVAR_DDENS_ID, &
      MOMX_VID => PRGVAR_MOMX_ID,   &
      MOMZ_VID => PRGVAR_MOMZ_ID

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars

    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      atm_dyn_dgm_nonhydro3d_common_DRHOT2PRES
    implicit none

    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(ElementBase3D), pointer :: elem

    integer :: n
    integer :: ke

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT, Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd_
    class(MeshField3D), pointer :: PRES_hyd, PRES_hyd_ref


    real(RP), allocatable :: DENS(:), rsfac(:), var_new(:)
    real(RP), allocatable :: Fx(:), Fz(:), LiftDelFlx(:)
    real(RP), allocatable :: del_flux(:,:,:)
    real(RP), allocatable :: u_(:), wt_(:)
    real(RP), allocatable :: PRES(:,:), DPRES(:,:)
    real(DP) :: dt, dt_
    !------------------------------------------
    
    if ( this%USER_do ) then
      call atm%vars%Calc_diagVar( 'PT_diff', PT_diff )
      call FILE_HISTORY_meshfield_in( PT_diff, "perturbation of potential temperature" )
    end if

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd_, Rtot, CVtot, CPtot, lcmesh          )      
      elem => lcmesh%refElem3D

      allocate( Fx(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np), del_flux(elem%NfpTot,lcmesh%Ne,6) )
      allocate( u_(elem%Np), wt_(elem%Np) )

      ! call cal_bndflux( del_flux, &
      !   DDENS%val, MOMX%val, MOMY%val, MOMZ%val, &
      !   DENS_hyd%val, PRES_hyd_%val, & ! (in)
      !   lcmesh%Gsqrt, lcmesh%GI3(:,:,1), lcmesh%GI3(:,:,2),                        & ! (in)    
      !   lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), & ! (in)
      !   lcmesh%vmapM, lcmesh%vmapP,                                               & ! (in)
      !   lcmesh, elem, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D )  
      !                                           ! (in)

      ! allocate( PRES(elem%Np,lcmesh%NeA), DPRES(elem%Np,lcmesh%NeA) )
      ! call atm_dyn_dgm_nonhydro3d_common_DRHOT2PRES( PRES, DPRES, &
      !   DRHOT%val, PRES_hyd_%val, Rtot%val, CVtot%val, CPtot%val, lcmesh, elem )
      
      !$omp parallel do private(ke, Fx, Fz, LiftDelFlx, u_, wt_)
      do ke=lcmesh%NeS, lcmesh%NeE
        GsqrtDENS%local(n)%val(:,ke) = lcmesh%Gsqrt(:,ke) * ( DENS_hyd%val(:,ke) + DDENS%val(:,ke) )
        GsqrtMOMX%local(n)%val(:,ke) = lcmesh%Gsqrt(:,ke) *  MOMX%val(:,ke)
        GsqrtMOMZ%local(n)%val(:,ke) = lcmesh%Gsqrt(:,ke) *  MOMZ%val(:,ke)
        GsqrtMOMW%local(n)%val(:,ke) = MOMZ%val(:,ke) +  lcmesh%Gsqrt(:,ke) * lcmesh%GI3(:,ke,1) * MOMX%val(:,ke)
        ! GsqrtDPRES%local(n)%val(:,ke) = lcmesh%Gsqrt(:,ke) * DPRES(:,ke)
        ! GsqrtG13DPRES%local(n)%val(:,ke) = lcmesh%Gsqrt(:,ke) * lcmesh%GI3(:,ke,1) * DPRES(:,ke)

    !     u_(:) = MOMX%val(:,ke) / DENS_hyd%val(:,ke) 
    !     wt_(:) = ( MOMZ%val(:,ke) / lcmesh%Gsqrt(:,ke) + lcmesh%GI3(:,ke,1) * MOMX%val(:,ke) ) / DENS_hyd%val(:,ke)

        ! call sparsemat_matmul( atm%mesh%DOptrMat(1), GsqrtMOMX%local(n)%val(:,ke), Fx(:) )
        ! call sparsemat_matmul( atm%mesh%LiftOptrMat, lcmesh%Fscale(:,ke) * del_flux(:,ke,1), LiftDelFlx )
        ! DENS_dt_1%local(n)%val(:,ke) = - ( lcmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) ) !/ lcmesh%Gsqrt(:,ke)

        ! call sparsemat_matmul( atm%mesh%DOptrMat(3), GsqrtMOMW%local(n)%val(:,ke), Fz(:) )
        ! call sparsemat_matmul( atm%mesh%LiftOptrMat, lcmesh%Fscale(:,ke) * del_flux(:,ke,2), LiftDelFlx )
        ! DENS_dt_2%local(n)%val(:,ke) = - ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) ) !/ lcmesh%Gsqrt(:,ke)

        ! call sparsemat_matmul( atm%mesh%DOptrMat(3), MOMZ%val(:,ke), Fz(:) )
        ! call sparsemat_matmul( atm%mesh%LiftOptrMat, lcmesh%Fscale(:,ke) * del_flux(:,ke,3), LiftDelFlx )
        ! DENS_dt_3%local(n)%val(:,ke) = - ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) ) !/ lcmesh%Gsqrt(:,ke)

        ! call sparsemat_matmul( atm%mesh%DOptrMat(3), GsqrtMOMW%local(n)%val(:,ke) - MOMZ%val(:,ke), Fz(:) )
        ! call sparsemat_matmul( atm%mesh%LiftOptrMat, lcmesh%Fscale(:,ke) * del_flux(:,ke,4), LiftDelFlx )
        ! DENS_dt_4%local(n)%val(:,ke) = - ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) ) !/ lcmesh%Gsqrt(:,ke)

        ! call sparsemat_matmul( atm%mesh%LiftOptrMat, lcmesh%Fscale(:,ke) * del_flux(:,ke,5), LiftDelFlx )
        ! DENS_dt_5%local(n)%val(:,ke) = - ( LiftDelFlx(:) ) !/ lcmesh%Gsqrt(:,ke)

        ! call sparsemat_matmul( atm%mesh%LiftOptrMat, lcmesh%Fscale(:,ke) * del_flux(:,ke,6), LiftDelFlx )
        ! DENS_dt_6%local(n)%val(:,ke) = - ( LiftDelFlx(:) ) !/ lcmesh%Gsqrt(:,ke)
      end do
    end do
    call FILE_HISTORY_meshfield_in( GsqrtDENS, "dens" )
    call FILE_HISTORY_meshfield_in( GsqrtMOMX, "momx" )
    call FILE_HISTORY_meshfield_in( GsqrtMOMZ, "momz" )
    call FILE_HISTORY_meshfield_in( GsqrtMOMW, "momw" )
    ! call FILE_HISTORY_meshfield_in( GsqrtDPRES, "GsqrtDPRES" )
    ! call FILE_HISTORY_meshfield_in( GsqrtG13DPRES, "GsqrtG13DPRES" )
    ! call FILE_HISTORY_meshfield_in( DENS_dt_1, "DENS_dt_1" )
    ! call FILE_HISTORY_meshfield_in( DENS_dt_2, "DENS_dt_2" )
    ! call FILE_HISTORY_meshfield_in( DENS_dt_3, "DENS_dt_3" )
    ! call FILE_HISTORY_meshfield_in( DENS_dt_4, "DENS_dt_4" )
    ! call FILE_HISTORY_meshfield_in( DENS_dt_5, "DENS_dt_5" )
    ! call FILE_HISTORY_meshfield_in( DENS_dt_6, "DENS_dt_6" )


    ! Set reference hydrostatic pressure
    if ( .not. is_PREShyd_ref_set ) then
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_VID, PRES_hyd )
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_REF_VID, PRES_hyd_ref )

      do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
        lcmesh => atm%mesh%ptr_mesh%lcmesh_list(n)
        do ke=lcmesh%NeS, lcmesh%NeE
!          PRES_hyd_ref%local(n)%val(:,ke) = PRES_hyd%local(n)%val(:,ke)
!          PRES_hyd_ref%local(n)%val(:,ke) = PRES00 * exp( - lcmesh%zlev(:,ke) / ( Rdry * 300.0_RP / Grav ) )
        end do
      end do

      is_PREShyd_ref_set = .true.
    end if

    call FILE_HISTORY_meshfield_in( sfac_save, "sfac" )

    return
  end subroutine USER_calc_tendency

  subroutine cal_bndflux( del_flux, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DENS_hyd, PRES_hyd, & ! (in)
    Gsqrt, G13, G23, nx, ny, nz,                                     & ! (in)
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D                       ) ! (in)
    use scale_element_base, only: ElementBase2D
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne,6)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)

    integer :: ke, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: GsqrtDensM(elem%NfpTot), GsqrtDensP(elem%NfpTot)
    real(RP) :: GsqrtDDENS_P(elem%NfpTot), GsqrtDDENS_M(elem%NfpTot)
    real(RP) :: GsqrtMOMX_P(elem%NfpTot), GsqrtMOMX_M(elem%NfpTot)
    real(RP) :: GsqrtMOMY_P(elem%NfpTot), GsqrtMOMY_M(elem%NfpTot)
    real(RP) :: GsqrtMOMZ_P(elem%NfpTot), GsqrtMOMZ_M(elem%NfpTot)
    real(RP) :: VelP_v(elem%NfpTot), VelM_v(elem%NfpTot)
    real(RP) :: VelP_h(elem%NfpTot), VelM_h(elem%NfpTot)
    real(RP) :: alpha(elem%NfpTot)
    real(RP) :: Gsqrt_P(elem%NfpTot), Gsqrt_M(elem%NfpTot)
    real(RP) :: GsqrtV_P(elem%NfpTot), GsqrtV_M(elem%NfpTot)
    real(RP) :: G13_M(elem%NfpTot), G13_P(elem%NfpTot)
    real(RP) :: G23_M(elem%NfpTot), G23_P(elem%NfpTot)

    !---------
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      Gsqrt_M(:) = Gsqrt(iM)
      Gsqrt_P(:) = Gsqrt(iP)
      GsqrtV_M(:) = Gsqrt_M(:)
      GsqrtV_P(:) = Gsqrt_P(:)

      G13_M(:) = G13(iM)
      G13_P(:) = G13(iP)
      G23_M(:) = G23(iM)
      G23_P(:) = G23(iP)
      GsqrtDDENS_M(:) = Gsqrt_M(:) * DDENS_(iM)
      GsqrtDDENS_P(:) = Gsqrt_P(:) * DDENS_(iP)
      GsqrtMOMX_M (:) = Gsqrt_M(:) * MOMX_ (iM)
      GsqrtMOMX_P (:) = Gsqrt_P(:) * MOMX_ (iP)
      GsqrtMOMY_M (:) = Gsqrt_M(:) * MOMY_ (iM)
      GsqrtMOMY_P (:) = Gsqrt_P(:) * MOMY_ (iP)
      GsqrtMOMZ_M (:) = Gsqrt_M(:) * MOMZ_ (iM)
      GsqrtMOMZ_P (:) = Gsqrt_P(:) * MOMZ_ (iP)

      GsqrtDensM(:) = GsqrtDDENS_M(:) + Gsqrt_M(:) * DENS_hyd(iM)
      GsqrtDensP(:) = GsqrtDDENS_P(:) + Gsqrt_P(:) * DENS_hyd(iP)

      VelM_h(:) = ( GsqrtMOMX_M(:) * nx(:,ke) + GsqrtMOMY_M(:) * ny(:,ke) ) / GsqrtDensM(:)
      VelP_h(:) = ( GsqrtMOMX_P(:) * nx(:,ke) + GsqrtMOMY_P(:) * ny(:,ke) ) / GsqrtDensP(:)
      
      VelM_v(:) = ( ( ( GsqrtMOMZ_M(:) / GsqrtV_M(:)                                       &
                    + G13_M(:) * GsqrtMOMX_M(:) + G23_M(:) * GsqrtMOMY_M(:) ) * nz(:,ke) ) &
                  ) / GsqrtDensM(:)
      VelP_v(:) = ( ( ( GsqrtMOMZ_P(:) / GsqrtV_P(:)                                       &
                    + G13_P(:) * GsqrtMOMX_P(:) + G23_P(:) * GsqrtMOMY_P(:) ) * nz(:,ke) ) &
                  ) / GsqrtDensP(:)

      alpha(:) = max( sqrt(1.41_RP * PRES_hyd(iM) / DENS_hyd(iM)) + abs(VelM_h+VelM_v), &
                      sqrt(1.41_RP * PRES_hyd(iP) / DENS_hyd(iP)) + abs(VelP_h+VelP_v) )

      del_flux(:,ke,1) = 0.5_RP * ( &
                    ( GsqrtDensP(:) * VelP_h(:) - GsqrtDensM(:) * VelM_h(:) )  &
                    - alpha(:) * (1.0_RP-nz(:,ke)**2) * ( GsqrtDDENS_P(:) - GsqrtDDENS_M(:) )     )
      del_flux(:,ke,2) = 0.5_RP * ( &
                    ( GsqrtDensP(:) * VelP_v(:) - GsqrtDensM(:) * VelM_v(:) )  &
                    - alpha(:) * nz(:,ke)**2 * ( GsqrtDDENS_P(:) - GsqrtDDENS_M(:) )     )
      where ( nz(:,ke) < -0.9_RP .and. iP(:) > elem%Np * lmesh%Ne )
        del_flux(:,ke,2) = - GsqrtDensM(:) * VelM_v(:)
      end where


      VelM_v(:) = ( ( ( GsqrtMOMZ_M(:) / GsqrtV_M(:)  ) * nz(:,ke) ) &
                  ) / GsqrtDensM(:)
      VelP_v(:) = ( ( ( GsqrtMOMZ_P(:) / GsqrtV_P(:) ) * nz(:,ke) ) &
                  ) / GsqrtDensP(:)
      del_flux(:,ke,3) = 0.5_RP * ( &
                    ( GsqrtDensP(:) * VelP_v(:) - GsqrtDensM(:) * VelM_v(:) )  &
                    - 0.0_RP*alpha(:) * nz(:,ke)**2 * ( GsqrtDDENS_P(:) - GsqrtDDENS_M(:) )     )
      where ( nz(:,ke) < -0.9_RP .and. iP(:) > elem%Np * lmesh%Ne )
        del_flux(:,ke,3) = - GsqrtDensM(:) * VelM_v(:)
      end where


      VelM_v(:) = ( ( ( G13_M(:) * GsqrtMOMX_M(:) + G23_M(:) * GsqrtMOMY_M(:) ) * nz(:,ke) ) &
                  ) / GsqrtDensM(:)
      VelP_v(:) = ( ( ( G13_P(:) * GsqrtMOMX_P(:) + G23_P(:) * GsqrtMOMY_P(:) ) * nz(:,ke) ) &
                  ) / GsqrtDensP(:)
      del_flux(:,ke,4) = 0.5_RP * ( &
                    ( GsqrtDensP(:) * VelP_v(:) - GsqrtDensM(:) * VelM_v(:) )  &
                    - 0.0_RP*alpha(:) * nz(:,ke)**2 * ( GsqrtDDENS_P(:) - GsqrtDDENS_M(:) )     )
      where ( nz(:,ke) < -0.9_RP .and. iP(:) > elem%Np * lmesh%Ne )
        del_flux(:,ke,4) = - GsqrtDensM(:) * VelM_v(:)
      end where
      
      del_flux(:,ke,5) = 0.5_RP * ( &
                    - alpha(:) * nz(:,ke)**2 * ( GsqrtDDENS_P(:) - GsqrtDDENS_M(:) )     )
      del_flux(:,ke,6) = 0.5_RP * ( &
                    - alpha(:) * (1.0_RP-nz(:,ke)**2) * ( GsqrtDDENS_P(:) - GsqrtDDENS_M(:) )     )
    end do

    return
  end subroutine cal_bndflux

  !------

!OCL SERIAL  
  subroutine exp_SetInitCond_mountain_wave( this,                 &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constT, &
      hydrostatic_calc_basicstate_constPT, &
      hydrostatic_calc_basicstate_constBVFreq

    use mod_experiment, only: &
      TracerLocalMeshField_ptr
    use mod_mkinit_util, only: &
      mkinitutil_gen_GPMat

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

    integer :: IntrpPolyOrder_h
    integer :: IntrpPolyOrder_v   
    namelist /PARAM_EXP/ &
      U0,                &
      TEMP0,             &
      U_tappering_flag,  &
      U_TAPPERING_HEIGHT, U_TAPPERING_WIDTH, &
      IntrpPolyOrder_h, IntrpPolyOrder_v

    integer :: ke
    integer :: ierr

    real(RP) :: Usfac(elem%Np)

    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: IntrpMat(:,:)

    real(RP), allocatable :: Gsqrt_ip(:,:), DENS_hyd_ip(:,:), PRES_hyd_ip(:,:), MOMX_ip(:,:), RHOT_ip(:,:), DDENS_ip(:,:)
    real(RP), allocatable :: topo_ip(:,:)
    real(RP) :: DENS_gp(elem%Np)
    real(RP) :: PRES_gp(elem%Np)
    real(RP) :: RHOT_gp(elem%Np)
    real(RP) :: Gsqrt_gp(elem%Np)    
    !-----------------------------------------------------------------------------

    U_tappering_flag = .false.
    U_TAPPERING_HEIGHT = 6E3_RP
    U_TAPPERING_WIDTH = .5E3_RP

    IntrpPolyOrder_h = elem%PolyOrder_h
    IntrpPolyOrder_v = elem%PolyOrder_v

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

    U_TAPPERING_HEIGHT = 3E3_RP
    U_TAPPERING_WIDTH  = 1E3_RP

    call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, elem%IsLumpedMatrix() )

    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    call mkinitutil_gen_GPMat( IntrpMat, elem_intrp, elem )

    !--
    call hydrostatic_calc_basicstate_constT( DENS_hyd, PRES_hyd, & ! (out)
      TEMP0, PRES00,                                             & ! (in)
      x, y, lcmesh%zlev, lcmesh, elem                            ) ! (in)
    ! call hydrostatic_calc_basicstate_constPT( DENS_hyd, PRES_hyd, & ! (out)
    !   TEMP0, PRES00,                                             & ! (in)
    !   x, y, lcmesh%zlev, lcmesh, elem                            ) ! (in)
    ! call hydrostatic_calc_basicstate_constBVFreq( DENS_hyd, PRES_hyd, & ! (out)
    !   1.8E-4_RP, TEMP0, PRES00,                                             & ! (in)
    !   x, y, lcmesh%zlev, lcmesh, elem                            ) ! (in)

    ! PRES_hyd(:,:) = PRES00
    ! DENS_hyd(:,:) = PRES00 / ( Rdry * TEMP0 )
!    return
    !---
    allocate( Gsqrt_ip(elem_intrp%Np,lcmesh%NeA), PRES_hyd_ip(elem_intrp%Np,lcmesh%NeA), DENS_hyd_ip(elem_intrp%Np,lcmesh%NeA) )
    allocate( DDENS_ip(elem_intrp%Np,lcmesh%NeA) )
    allocate( MOMX_ip(elem_intrp%Np,lcmesh%NeA), RHOT_ip(elem_intrp%Np,lcmesh%NeA) )

    call calc_init_profile_GP( Gsqrt_ip, DENS_hyd_ip, PRES_hyd_ip, DDENS_ip, MOMX_ip, RHOT_ip, &
      lcmesh, elem_intrp%x1, elem_intrp%x2, elem_intrp%x3, elem_intrp%Np, elem_intrp%Nv )

    !----
    !$omp parallel do private(Usfac, DENS_gp, RHOT_gp, PRES_gp, Gsqrt_gp)
    do ke=lcmesh%NeS, lcmesh%NeE
      if ( U_tappering_flag ) then
        Usfac(:) = 0.5_RP * ( 1.0_RP + tanh( ( lcmesh%zlev(:,ke) - U_TAPPERING_HEIGHT ) / U_TAPPERING_WIDTH ) )
      else
        Usfac(:) = 1.0_RP
      end if
      MOMX(:,ke) = DENS_hyd(:,ke) * U0! * Usfac(:)

      ! MOMX(:,ke) = matmul( IntrpMat, Gsqrt_ip(:,ke) * MOMX_ip(:,ke) ) / lcmesh%Gsqrt(:,ke)
      DENS_gp(:) = matmul( IntrpMat, Gsqrt_ip(:,ke) * DENS_hyd_ip(:,ke) ) / lcmesh%Gsqrt(:,ke)
!      DENS_gp(:) = matmul( IntrpMat, DENS_hyd_ip(:,ke) )
      ! MOMY(:,ke) = 0.0_RP
      
      RHOT_gp(:) = matmul( IntrpMat, Gsqrt_ip(:,ke) * RHOT_ip(:,ke) ) / lcmesh%Gsqrt(:,ke)
!      RHOT_gp(:) = matmul( IntrpMat, RHOT_ip(:,ke) )

!      DDENS(:,ke) = DENS_gp(:) - DENS_hyd(:,ke) 
      ! DDENS(:,ke) = DDENS(:,ke) + matmul( IntrpMat, Gsqrt_ip(:,ke) * DDENS_ip(:,ke) ) / lcmesh%Gsqrt(:,ke)
!      DRHOT(:,ke) = RHOT_gp(:) - PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CvDry/CpDry)
!       PRES_hyd(:,ke) = PRES00 * ( Rdry/PRES00 * RHOT_gp(:) )**(CpDry/CvDry)
!      PRES_hyd(:,ke) = matmul( IntrpMat, PRES_gp(:,ke) )

      ! DENS_hyd(:,ke) = matmul( IntrpMat, DENS_hyd_ip(:,ke) )
      ! PRES_hyd(:,ke) = matmul( IntrpMat, PRES_hyd_ip(:,ke) ) !/ lcmesh%Gsqrt(:,ke)
      DRHOT(:,ke) = RHOT_gp(:) - PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CvDry/CpDry)
      DDENS(:,ke) = DENS_gp(:) - DENS_hyd(:,ke)

    !  DDENS(:,ke) = &
    !     5E-3_RP * cos(4.0_RP*PI*x(:,ke)/240.0E3_RP) * cos(4.0_RP*PI*lcmesh%zlev(:,ke)/30.0E3_RP)
      !  DDENS(:,ke) = &
      !     1E-4_RP * DENS_hyd(:,ke)
      !DDENS(:,ke) = 5.E-4_RP * RHOT_gp(:)
    end do

    ! call calc_hydrostatic_pressure( PRES_hyd, &
    !   DENS_hyd, lcmesh, elem, lcmesh%lcmesh2D%refElem2D )

    call elem_intrp%Final()

    return
  end subroutine exp_SetInitCond_mountain_wave

  subroutine exp_geostrophic_balance_correction_lc( this,                &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    lcmesh, elem )
    implicit none
    class(Experiment), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(inout) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DRHOT(elem%Np,lcmesh%NeA)

    integer :: ke
    real(RP) :: DENS(elem%Np), Fz(elem%Np)
    real(RP) :: LiftDelFlux(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lcmesh%Ne)
    !--------------------------

    call cal_bndflux_ini( del_flux, &
      PRES_hyd, lcmesh%normal_fn(:,:,3), lcmesh%VMapM, lcmesh%VMapP, &
      lcmesh, elem, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D )
    
    !$omp parallel do private(Fz, LiftDelFlux, DENS)
    do ke=lcmesh%NeS, lcmesh%NeE
      Fz(:) = matmul( elem%Dx3, PRES_hyd(:,ke) )
      LiftDelFlux(:) = matmul( elem%Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke) )

      ! DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      ! DENS_hyd(:,ke) = - ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlux(:) ) / ( lcmesh%Gsqrt(:,ke ) * GRAV )
      ! DDENS(:,ke) = DENS(:) - DENS_hyd(:,ke)
    end do

    return
  end subroutine exp_geostrophic_balance_correction_lc

  subroutine cal_bndflux_ini( del_flux, &
    PRES_hyd, & ! (in)
    nz,                                     & ! (in)
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D                       ) ! (in)
    use scale_element_base, only: ElementBase2D
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)

    integer :: ke, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D

    !---------
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      del_flux(:,ke) = 0.5_RP * ( &
                    ( PRES_hyd(iP(:)) - PRES_hyd(iM(:))  ) * nz(:,ke) )
    end do

    return
  end subroutine cal_bndflux_ini

  subroutine calc_init_profile_GP( Gsqrt, DENS_hyd, PRES_hyd, DDENS, MOMX, RHOT, &
    lcmesh, elem_x1, elem_x2, elem_x3, Np, Nv ) ! (in)
    implicit none
    type(LocalMesh3D), intent(in) :: lcmesh
    integer, intent(in) :: Np, Nv
    real(RP), intent(out) :: Gsqrt(Np,lcmesh%NeA)
    real(RP), intent(out) :: DENS_hyd(Np,lcmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(Np,lcmesh%NeA)
    real(RP), intent(out) :: DDENS(Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMX(Np,lcmesh%NeA)
    real(RP), intent(out) :: RHOT(Np,lcmesh%NeA)
    real(RP), intent(in) :: elem_x1(Np)
    real(RP), intent(in) :: elem_x2(Np)
    real(RP), intent(in) :: elem_x3(Np)

    real(RP) :: vx(Nv), vy(Nv), vz(Nv)
    real(RP) :: x_(Np,lcmesh%Ne), y_(Np,lcmesh%Ne), zeta_(Np,lcmesh%Ne)
    real(RP) :: zlev(Np,lcmesh%Ne), topo(Np,lcmesh%Ne)

    integer :: ke
    !---------------------------
    
    !$omp parallel do private(vx, vy, vz)
    do ke=lcmesh%NeS, lcmesh%NeE
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
      vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)
      x_(:,ke) = vx(1) + 0.5_RP * ( elem_x1(:) + 1.0_RP ) * ( vx(2) - vx(1) ) 
      y_(:,ke) = vy(1) + 0.5_RP * ( elem_x2(:) + 1.0_RP ) * ( vy(4) - vy(1) )
      zeta_(:,ke) = vz(1) + 0.5_RP * ( elem_x3(:) + 1.0_RP ) * ( vz(5) - vz(1) )

      topo(:,ke) = 25.0_RP * exp( - ( x_(:,ke) - 120.E3_RP )**2 / 5.0E3_RP**2 ) &
!        * ( cos( PI * ( x_(:,ke) - 120.0E3_RP ) / 4.0E10_RP) )**2
        * ( cos( PI * ( x_(:,ke) - 120.0E3_RP ) / 4.0E3_RP) )**2
      zlev(:,ke) = ( 1.0E0_RP - topo(:,ke) / 30.0E3_RP ) * zeta_(:,ke) + topo(:,ke)
    end do

    !$omp parallel do
    do ke=lcmesh%NeS, lcmesh%NeE
      Gsqrt(:,ke) =  1.0_RP - topo(:,ke) / 30.0E3_RP
      PRES_hyd(:,ke) = PRES00 * exp( - zlev(:,ke) / ( Rdry * 300.0_RP / Grav ) )
      DENS_hyd(:,ke) = PRES_hyd(:,ke) / ( Rdry * 300.0_RP )
      MOMX(:,ke) = DENS_hyd(:,ke) &
        * U0 !* 0.5_RP * ( 1.0_RP + tanh( ( zlev(:,ke) - U_TAPPERING_HEIGHT ) / U_TAPPERING_WIDTH ) )
      RHOT(:,ke) =  DENS_hyd(:,ke) &
        * 300.0_RP * ( PRES00 / PRES_hyd(:,ke) )**(Rdry/CPdry)

!      DDENS(:,ke) = 0.01E0_RP * exp( - ( (x_(:,ke) - 110E3_RP)**2 + (zlev(:,ke) - 7E3_RP)**2 ) / 4E3_RP**2 )
    end do

    return
  end subroutine calc_init_profile_GP

  subroutine calc_hydrostatic_pressure( PRES_hyd, &
    DENS_hyd, lmesh, elem, elem2D )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lmesh%NeA)

    integer :: ke, ke_xy, ke_z, ke2D
    integer :: kl, ku, nz_1D

    real(RP) :: LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
    real(RP) :: DPDz_bc(elem%Np)

    real(RP) :: GsqrtV_z(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%NeX*lmesh%NeY)

    real(RP), allocatable :: PmatBnd(:,:,:)
    integer :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer :: vmapP(elem%NfpTot,lmesh%NeZ)
    integer :: ij
    integer :: ColMask(elem%Nnode_v)
    real(RP) :: b1D(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2,lmesh%NeX*lmesh%NeY)
    integer :: info
    integer :: ipiv(elem%Nnode_v*1*lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: SfcPRES(elem%Np,lmesh%NeA)
    !---------------------------

    nz_1D = elem%Nnode_v * 1 * lmesh%NeZ
    kl = elem%Nnode_v
    ku = kl
    allocate( PmatBnd   (2*kl+ku+1,nz_1D,elem%Nnode_h1D**2) )    

    do ke=lmesh%NeS, lmesh%NeE
      SfcPRES(:,ke) = DENS_hyd(:,ke) * Rdry * 300.0_RP
    end do
    call cal_del_flux( del_flux, &
      SfcPRES, lmesh%Gsqrt, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2), lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),                                     & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem, lmesh%lcmesh2D, elem2D                       ) ! (in)

    do ke_xy=1, lmesh%NeX*lmesh%NeY
    do ke_z=1, lmesh%NeZ
      ke = ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
      ke2D = lmesh%EMap3Dto2D(ke)
      nz(:,ke_z,ke_xy) = lmesh%normal_fn(:,ke,3)
      GsqrtV_z(:,ke_z,ke_xy) = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2D)

      LiftDelFlx(:) = matmul(elem%Lift, lmesh%Fscale(:,ke) * del_flux(:,ke))
      DPDz_bc(:) = LiftDelFlx(:) / lmesh%Gsqrt(:,ke)

      do ij=1, elem%Nnode_h1D**2          
        ColMask(:) = elem%Colmask(:,ij)
        b1D(:,ke_z,ij,ke_xy) = - DPDz_bc(ColMask(:)) - Grav * DENS_hyd(ColMask(:),ke)
      end do
    end do
    end do
      
    call set_vmapZ1D( vmapM, vmapP, &
      elem, lmesh )
    
    do ke_xy=1, lmesh%NeX * lmesh%NeY
      call atm_dyn_dgm_nonhydro3d_forcing_construct_matbnd( &
        PmatBnd,                                & ! (out)
        kl, ku, nz_1D,                          & ! (in)
        GsqrtV_z(:,:,ke_xy),                    & ! (in)
        lmesh, elem,                            & ! (in)
        nz(:,:,ke_xy), vmapM, vmapP, ke_xy, 1 )   ! (in)
      do ij=1, elem%Nnode_h1D**2
        call dgbsv( nz_1D, kl, ku, 1, PmatBnd(:,:,ij), 2*kl+ku+1, ipiv(:,ij), b1D(:,:,ij,ke_xy), nz_1D, info)
        ColMask(:) = elem%Colmask(:,ij)
        do ke_z=1, lmesh%NeZ
          ke = ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
          PRES_hyd(ColMask(:),ke) = b1D(:,ke_z,ij,ke_xy)
        end do
      end do ! for ij
    end do
    
    return
  end subroutine calc_hydrostatic_pressure

!OCL SERIAL
  subroutine set_vmapZ1D( vmapM, vmapP, &
    elem, lmesh )

    implicit none
    type(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(out) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(out) :: vmapP(elem%NfpTot,lmesh%NeZ)    

    integer :: ke_z
    integer :: f
    integer :: vs, ve
    !------------------------------

    do ke_z=1, lmesh%NeZ
      do f=1, elem%Nfaces_h
        vs = 1 + (f-1)*elem%Nfp_h
        ve = vs + elem%Nfp_h - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_h(:,f) + (ke_z-1)*elem%Np
      end do
      do f=1, elem%Nfaces_v
        vs = elem%Nfp_h*elem%Nfaces_h + 1 + (f-1)*elem%Nfp_v
        ve = vs + elem%Nfp_v - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_v(:,f) + (ke_z-1)*elem%Np
      end do
      vmapP(:,ke_z) = vmapM(:,ke_z)
    end do

    do ke_z=1, lmesh%NeZ
      vs = elem%Nfp_h*elem%Nfaces_h + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z > 1) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,2) + (ke_z-2)*elem%Np

      vs = elem%Nfp_h*elem%Nfaces_h + elem%Nfp_v + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z < lmesh%NeZ) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1) + ke_z*elem%Np
    end do

    return
  end subroutine set_vmapZ1D

  subroutine cal_del_flux( del_flux,                       & ! (in)
    SfcPRES, Gsqrt, G13, G23, nx, ny, nz,                                     & ! (in)
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D                       ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  SfcPRES(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)

    integer :: ke, ke_z, ke_xy, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: Gsqrt_P(elem%NfpTot), Gsqrt_M(elem%NfpTot)
    real(RP) :: GsqrtV_P(elem%NfpTot), GsqrtV_M(elem%NfpTot)
    real(RP) :: G13_M(elem%NfpTot), G13_P(elem%NfpTot)
    real(RP) :: G23_M(elem%NfpTot), G23_P(elem%NfpTot)
    real(RP) :: Gnn_M(elem%NfpTot), Gnn_P(elem%NfpTot)
    real(RP) :: SfcPRES_(elem%NfpTot)

    real(RP) :: gamm
    integer :: fp
    !------------------------------------------------------------------------

    gamm  = CPDry / CvDry

    do ke_z=1, lmesh%NeZ
    do ke_xy=1, lmesh%NeX*lmesh%NeY
      ke = ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      Gsqrt_M(:) = Gsqrt(iM)
      Gsqrt_P(:) = Gsqrt(iP)
      GsqrtV_M(:) = Gsqrt_M(:)
      GsqrtV_P(:) = Gsqrt_P(:)

      G13_M(:) = G13(iM)
      G13_P(:) = G13(iP)
      G23_M(:) = G23(iM)
      G23_P(:) = G23(iP)
      SfcPRES_(:) = SfcPRES(iM(:))
      
      where (nz(:,ke) < -0.5_RP .and. iP(:) > lmesh%Ne * elem%Np )
        del_flux(:,ke) = 0.5_RP * ( &
                      ( SfcPRES_(:) - 0.0_RP ) * nz(:,ke)  )
      elsewhere
        del_flux(:,ke) =  0.0_RP
      end where
    end do
    end do

    return
  end subroutine cal_del_flux

!OCL SERIAL  
  subroutine atm_dyn_dgm_nonhydro3d_forcing_construct_matbnd( &
    PmatBnd,                                & ! (out)
    kl, ku, nz_1D,                          & ! (in)
    GsqrtV,                                 & ! (in)
    lmesh, elem,                            & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y )            ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(in) :: kl, ku, nz_1D
    real(RP), intent(out) :: PmatBnd(2*kl+ku+1,elem%Nnode_v,1,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP), intent(in) ::  GsqrtV(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y

    integer :: ke_z, ke_z2
    integer :: v, ke, p, f1, f2, fp, fp2, FmV
    real(RP) :: fac_dz_p(elem%Nnode_v)
    real(RP) :: PmatD(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: PmatL(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: PmatU(elem%Nnode_v,elem%Nnode_v)

    integer :: Colmask(elem%Nnode_v)
    real(RP) :: tmp1
    real(RP) :: fac
    integer :: ij, v1, v2, pv1, pv2, pv11, pv22, g_kj, g_kjp1, g_kjm1, pb1
    logical :: bc_flag
    !--------------------------------------------------------

    !$omp parallel private(Colmask)
    !$omp workshare
    PmatD(:,:) = 0.0_RP
    PmatL(:,:) = 0.0_RP
    PmatU(:,:) = 0.0_RP  
    !$omp end workshare
    !$omp do
    do ij=1, elem%Nnode_h1D**2
      PmatBnd   (:,:,:,:,ij) = 0.0_RP
    end do
    !$omp end parallel

    !$omp parallel private(ke_z, ke, ColMask, p, fp, fp2, v, f1, f2, ke_z2, fac_dz_p,    &
    !$omp fac, tmp1, FmV,  pv11, pv22,                                                             &
    !$omp ij, v1, v2, pv1, pv2, pb1, g_kj, g_kjp1, g_kjm1, bc_flag)                      &
    !$omp firstprivate(PmatD, PmatL, PmatU)

    !$omp do collapse(2)
    do ij=1, elem%Nnode_h1D**2
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
      Colmask(:) = elem%Colmask(:,ij)

      !-----  
      do p=1, elem%Nnode_v
        fac_dz_p(:) = lmesh%Escale(Colmask(:),ke,3,3) / GsqrtV(Colmask(:),ke_z) &
                    * elem%Dx3(Colmask(:),Colmask(p))
        ! DDENS
        PmatD(:,p) = fac_dz_p(:)
      end do

      do f1=1, 2
        if (f1==1) then
          ke_z2 = max(ke_z-1,1)
          pv1 = 1; pv2 = elem%Nnode_v
          f2 = 2
        else
          ke_z2 = min(ke_z+1,lmesh%NeZ)
          pv1 = elem%Nnode_v; pv2 = 1
          f2 = 1
        end if
        fac  = 0.5_RP / GsqrtV(Colmask(pv1),ke_z)
        if ( (ke_z == 1 .and. f1==1) .or. (ke_z == lmesh%NeZ .and. f1==elem%Nfaces_v) ) then
          bc_flag = .true.
          pv2 = pv1; f2 = f1
        else 
          bc_flag = .false.    
        end if

        FmV  = elem%Fmask_v(ij,f1)
        fp  = elem%Nfp_h * elem%Nfaces_h + (f1-1)*elem%Nfp_v + ij
        fp2 = elem%Nfp_h * elem%Nfaces_h + (f2-1)*elem%Nfp_v + ij

        !--
        ! tmp1 = fac * elem%Lift(FmV,fp) * lmesh%Fscale(fp,ke) &
        !      * max( alph(fp,ke_z), alph(fp2,ke_z2) )
        ! if (bc_flag) then
        ! else         
        !   PmatD(pv1,pv1) = PmatD(pv1,pv1) + tmp1            
        !   if (f1 == 1) then
        !     PmatL(pv1,pv2) = - tmp1                                
        !   else
        !     PmatU(pv1,pv2) = - tmp1
        !   end if
        ! end if 

        !--
        ! tmp1 = fac * elem%Lift(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)
        ! if (bc_flag) then
        !   if ( ke_z == 1 ) then
        !     PmatD(pv1,pv1) = PmatD(pv1,pv1) - tmp1
        !   end if
        ! else 
        !   PmatD(pv1,pv1) = PmatD(pv1,pv1) - tmp1
        !   if (f1 == 1) then
        !     PmatL(pv1,pv2) = + tmp1
        !   else
        !     PmatU(pv1,pv2) = + tmp1
        !   end if
        ! end if

        do pv11=1, elem%Nnode_v
          tmp1 = fac * elem%Lift(Colmask(pv11),fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)
          if (bc_flag) then
            if ( ke_z == 1 ) then
              PmatD(pv11,pv1) = PmatD(pv11,pv1) - tmp1
            end if
          else 
            PmatD(pv11,pv1) = PmatD(pv11,pv1) - tmp1
            if (f1 == 1) then
              PmatL(pv11,pv2) = + tmp1
            else
              PmatU(pv11,pv2) = + tmp1
            end if
          end if          
        end do
      end do

      do pv2=1, elem%Nnode_v
        g_kj   = pv2 + (ke_z-1)*elem%Nnode_v
        g_kjm1 = pv2 + (ke_z-2)*elem%Nnode_v
        g_kjp1 = pv2 + (ke_z  )*elem%Nnode_v

        do pv1=1, elem%Nnode_v            
          pb1 = pv1 + (ke_z-1)*elem%Nnode_v
          if (ke_z > 1 .and. pv2 == elem%Nnode_v ) then
            PmatBnd(kl+ku+1+pb1-g_kjm1,pv2,1,ke_z-1, ij) = PmatL(pv1,pv2)
          end if
          PmatBnd(kl+ku+1+pb1-g_kj,pv2,1,ke_z, ij) = PmatD(pv1,pv2)
          if (ke_z < lmesh%NeZ .and. pv2 == 1 ) then
            PmatBnd(kl+ku+1+pb1-g_kjp1,pv2,1,ke_z+1, ij) = PmatU(pv1,pv2)
          end if
        end do
      end do            

    end do  
    end do  
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine atm_dyn_dgm_nonhydro3d_forcing_construct_matbnd

end module mod_user
