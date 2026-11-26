!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for an idealized regional moist convection experiment
!!          (horizontally doubly-periodic domain, f = 0), which is used for investigating
!!          physics-dynamics coupling in SCALE-DG.
!!
!!          The setup is based on "preset A" (moist Held-Suarez-like tropical environment):
!!           - Horizontally uniform basic state given by smooth vertical profiles of
!!             potential temperature and relative humidity
!!           - Newtonian relaxation of temperature toward a horizontally uniform reference profile T_ref(z)
!!             (TREF_RELAX_TYPE = 'THETA' (default): rho*theta is relaxed toward rho*theta_eq with
!!              theta_eq = T_ref (P00/PRE_ref)^(Rtot/CPtot), the form used in the experiments on Fugaku;
!!              TREF_RELAX_TYPE = 'HS': T is relaxed toward T_ref at fixed density, as in the moist Held-Suarez test)
!!           - Surface sensible heat and moisture fluxes given by the bulk formula (ATMOS_PHY_SF, SIMPLE)
!!             with a fixed ground temperature
!!           - Simple boundary-layer vertical mixing (mod_user_sub_BLmixing)
!!           - Rayleigh damping of momentum in the boundary layer and a sponge layer near the model top
!!
!!          This module was reconstructed in Nov. 2026 because the original source code (used on Fugaku)
!!          was lost. See model/atm_nonhydro3d/test/product_run/moist_conv/NOTES.md for details.
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_user

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
    PRC_abort, PRC_myrank, PRC_LOCAL_COMM_WORLD
  use scale_const, only: &
    Rdry => CONST_Rdry,    &
    Rvap => CONST_Rvap,    &
    CPdry => CONST_CPdry,  &
    CVdry => CONST_CVdry,  &
    PRES00 => CONST_PRE00, &
    Grav => CONST_GRAV,    &
    PI => CONST_PI
  use scale_tracer, only: &
    TRACER_inq_id

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_meshfield_base, only: &
    MeshField3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use mod_user_base, only: UserBase
  use mod_experiment, only: Experiment

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    PHYTEND_NUM1 => PHYTEND_NUM, &
    MOMX_p  => PHYTEND_MOMX_ID,  &
    MOMY_p  => PHYTEND_MOMY_ID,  &
    RHOT_p  => PHYTEND_RHOT_ID
  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    AUXVAR_THERMHYDRO_ID

  use mod_atmos_vars, only: &
    AtmosVarsContainer,               &
    AtmosVars_GetLocalMeshPrgVars,    &
    AtmosVars_GetLocalMeshPhyAuxVars, &
    ATM_VARS_CONTAINER_PRIMARY_ID

  use mod_user_sub_BLmixing, only: &
    USER_sub_BLmixing_Init, USER_sub_BLmixing_calc_tendency
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
    procedure :: setup_vars => USER_setup_vars
    procedure :: calc_tendency => USER_calc_tendency
    procedure :: update_pre => USER_update_pre
    procedure :: update => USER_update
  end type User

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  !-- Parameters of the basic state (namelist PARAM_EXP) ------------------
  !   PARAM_EXP is read both in init.conf (making initial data) and in run.conf.
  !   NOTE: Following the original run template (run_case1_template.conf),
  !         SFC_THETA has different meanings in the two phases:
  !          - init phase: potential temperature of the basic state at the surface (theta_s; 300 K in preset A)
  !          - run  phase: ground temperature used by the surface flux scheme (305 K in the original run template)

  real(RP) :: SFC_THETA       = 300.0_RP   !< [init] Surface potential temperature of the basic state [K]
                                           !< [run ] Ground temperature for the surface flux scheme [K]
  real(RP) :: SFC_PRES        = 101325.0_RP !< Surface pressure [Pa]

  real(RP) :: ENV_PT_LAPS_BL  = 0.3E-3_RP  !< Vertical gradient of potential temperature in the boundary layer [K/m]
  real(RP) :: ENV_PT_LAPS_FT  = 3.0E-3_RP  !< Vertical gradient of potential temperature in the free troposphere [K/m]
  real(RP) :: ENV_PT_LAPS_STR = 6.0E-3_RP  !< Vertical gradient of potential temperature in the upper layer [K/m] (restoration candidate of preset A; the original A used 1.5E-3)
  real(RP) :: ENV_PBL_Z       = 1200.0_RP  !< Height of the transition at the PBL top [m]
  real(RP) :: ENV_PBL_W       =  400.0_RP  !< Width of the transition at the PBL top [m]
  real(RP) :: ENV_TRO_Z       = 11000.0_RP !< Height of the transition at the tropopause [m]
  real(RP) :: ENV_TRO_W       = 1000.0_RP  !< Width of the transition at the tropopause [m]
  real(RP) :: ENV_MID_Z       = 5000.0_RP  !< Height of the transition of RH in the middle troposphere [m]
  real(RP) :: ENV_MID_W       = 2000.0_RP  !< Width of the transition of RH in the middle troposphere [m]
  real(RP) :: ENV_RH_SFC      = 0.85_RP    !< Relative humidity in the lower layer (ratio, not percent)
  real(RP) :: ENV_RH_MID      = 0.55_RP    !< Relative humidity in the middle layer (ratio, not percent)
  real(RP) :: ENV_RH_UP       = 0.20_RP    !< Relative humidity in the upper layer (ratio, not percent)
  real(RP) :: ENV_U           = 0.0_RP     !< Background wind (x) [m/s]
  real(RP) :: ENV_V           = 0.0_RP     !< Background wind (y) [m/s]
  real(RP) :: RANDOM_THETA      = 1.0_RP   !< Amplitude of random perturbation of potential temperature [K]
  real(RP) :: RANDOM_THETA_ZTOP = 1200.0_RP !< Perturbation is added below this height [m]
  integer  :: RANDOM_ISEED      = 1        !< Seed of random numbers (the MPI rank is added)

  namelist / PARAM_EXP / &
    SFC_THETA, SFC_PRES,                                &
    ENV_PT_LAPS_BL, ENV_PT_LAPS_FT, ENV_PT_LAPS_STR,    &
    ENV_PBL_Z, ENV_PBL_W, ENV_TRO_Z, ENV_TRO_W,         &
    ENV_MID_Z, ENV_MID_W,                               &
    ENV_RH_SFC, ENV_RH_MID, ENV_RH_UP,                  &
    ENV_U, ENV_V,                                       &
    RANDOM_THETA, RANDOM_THETA_ZTOP, RANDOM_ISEED

  !-- Parameters of the forcing (namelist PARAM_USER) ---------------------
  integer  :: BL_VARS_CONTAINER_ID                   !< ID of variable container used by BL mixing
  character(len=H_LONG) :: in_filebase_Tref  = ''    !< Basename of history file from which T_ref is computed.
                                                     !< If empty, T_ref is computed from the initial state of this run.
  integer  :: Tref_step     = 1                      !< Time step in the history file used for T_ref (1: initial time)
  real(RP) :: TAU_RAD_BL    = 6.0_RP  * 86400.0_RP   !< Relaxation time of temperature in the boundary layer [s]
  real(RP) :: TAU_RAD_FT    = 20.0_RP * 86400.0_RP   !< Relaxation time of temperature in the free troposphere [s]
  real(RP) :: TAU_RAD_STR   = 30.0_RP * 86400.0_RP   !< Relaxation time of temperature in the upper layer [s]
  real(RP) :: PBL_MOM_TAU   = 86400.0_RP             !< Time scale of Rayleigh damping of momentum in the PBL [s] (<=0: off)
  character(len=H_SHORT) :: PBL_MOM_TYPE = 'HEIGHT'  !< Vertical profile of the PBL Rayleigh damping
                                                     !<  'HEIGHT': (1 - S(z; ENV_PBL_Z, ENV_PBL_W)) / PBL_MOM_TAU, applied to MOMX, MOMY, and MOMZ
                                                     !<            (the form used in the experiments on Fugaku)
                                                     !<  'SIGMA' : max(0, (sig - PBL_MOM_SIGB)/(1 - PBL_MOM_SIGB)) / PBL_MOM_TAU, applied to MOMX and MOMY
                                                     !<            (as in Held and Suarez 1994)
  real(RP) :: PBL_MOM_SIGB  = 0.7_RP                 !< Top of the Rayleigh damping layer in sigma (only for PBL_MOM_TYPE = 'SIGMA')
  real(RP) :: SPONGE_HEIGHT = 15.0E3_RP              !< Bottom height of the sponge layer [m]
  real(RP) :: TAU_SPONGE    = 43200.0_RP             !< E-folding time of the sponge layer at the model top [s] (<=0: off)
  logical  :: SPONGE_DAMP_TEMP = .true.              !< If true, temperature is also relaxed toward T_ref in the sponge layer
                                                     !< (the form used in the experiments on Fugaku)
  character(len=H_SHORT) :: TREF_RELAX_TYPE = 'THETA' !< Form of the Newtonian relaxation toward T_ref
                                                     !<  'THETA': d(rho theta)/dt = -(rho theta - rho theta_eq)/tau at fixed density,
                                                     !<           theta_eq = T_ref (P00/PRE_ref)^(Rtot/CPtot), PRE_ref = horizontal mean pressure
                                                     !<           at the same time as T_ref (default; the form used in the experiments on Fugaku).
                                                     !<           With 'HS', the model top tended to become unstable owing to strong horizontal winds.
                                                     !<  'HS'   : dT/dt = -(T - T_ref)/tau at fixed density (moist Held-Suarez form)

  !-- Reference temperature profile T_ref(z) -----------------------------
  logical :: is_Tref_ready
  real(RP), allocatable :: Tref_z(:,:)               !< T_ref at the vertical nodes (Nnode_v, NeZ)
  real(RP), allocatable :: PREref_z(:,:)             !< PRE_ref (horizontal mean pressure) at the vertical nodes (Nnode_v, NeZ)
  logical :: relax_theta                             !< .true. if TREF_RELAX_TYPE = 'THETA'

  logical :: is_setup_bl_vars_container

contains

!OCL SERIAL
  subroutine USER_mkinit ( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init( 'IdealizedMoistConv' )
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_idealized_moist_conv )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL
  subroutine USER_setup( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout), target :: atm

    logical :: USER_do = .false. !< do user step?

    namelist / PARAM_USER / &
       USER_do,              &
       BL_VARS_CONTAINER_ID, &
       in_filebase_Tref,     &
       Tref_step,            &
       TAU_RAD_BL, TAU_RAD_FT, TAU_RAD_STR, &
       PBL_MOM_TAU, PBL_MOM_TYPE, PBL_MOM_SIGB, &
       SPONGE_HEIGHT, TAU_SPONGE, SPONGE_DAMP_TEMP, &
       TREF_RELAX_TYPE
    integer :: ierr
    !------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    BL_VARS_CONTAINER_ID = ATM_VARS_CONTAINER_PRIMARY_ID

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

    select case( trim(PBL_MOM_TYPE) )
    case( 'HEIGHT', 'SIGMA' )
    case default
      LOG_ERROR("USER_setup",*) 'Unsupported PBL_MOM_TYPE: '//trim(PBL_MOM_TYPE)//'. Check!'
      call PRC_abort
    end select

    select case( trim(TREF_RELAX_TYPE) )
    case( 'HS' )
      relax_theta = .false.
    case( 'THETA' )
      relax_theta = .true.
    case default
      LOG_ERROR("USER_setup",*) 'Unsupported TREF_RELAX_TYPE: '//trim(TREF_RELAX_TYPE)//'. Check!'
      call PRC_abort
    end select

    !- In the run phase, SFC_THETA means the ground temperature (see the comment on PARAM_EXP).
    call read_param_exp( "USER_setup" )
    LOG_INFO("USER_setup",*) 'Ground temperature for the surface flux (SFC_THETA) [K]: ', SFC_THETA

    call this%UserBase%Setup( atm, USER_do )

    !- Setup BL mixing
    call USER_sub_BLmixing_Init( atm%mesh%ptr_mesh )
    is_setup_bl_vars_container = .false.

    !- T_ref is determined just before the first time step (USER_update_pre),
    !  because restart data has not been read yet at this point.
    is_Tref_ready = .false.

    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_setup_vars( this, atm )
    use scale_meshfield_base, only: MeshFieldBase
    use mod_atmos_phy_sfc_vars, only: ATMOS_PHY_SF_SVAR_TEMP_ID
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    integer :: n, ke2D
    class(MeshBase2D), pointer :: mesh2D
    class(LocalMesh2D), pointer :: lcmesh2D

    class(MeshFieldBase), pointer :: SfcTemp
    class(LocalMeshFieldBase), pointer :: SfcTemp_lc
    !---------------------------------------------------------------------------

    if ( .not. atm%phy_sfc_proc%IsActivated() ) return

    call atm%mesh%ptr_mesh%GetMesh2D( mesh2D )
    call atm%phy_sfc_proc%vars%SFCVARS_manager%Get( ATMOS_PHY_SF_SVAR_TEMP_ID, SfcTemp )

    do n=1, mesh2D%LOCAL_MESH_NUM
      lcmesh2D => mesh2D%lcmesh_list(n)
      call SfcTemp%GetLocalMeshField( n, SfcTemp_lc )
      !$omp parallel do
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
        SfcTemp_lc%val(:,ke2D) = SFC_THETA
      end do
    end do

    return
  end subroutine USER_setup_vars

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(AtmosVarsContainer), pointer :: vars_container
    class(AtmosVarsContainer), pointer :: vars_container_BL
    !----------------------------

    if ( .not. this%USER_do ) return

    call atm%vars%Get_container( ATM_VARS_CONTAINER_PRIMARY_ID, vars_container )

    if ( .not. is_setup_bl_vars_container )  then
        call atm%vars%Setup_container( BL_VARS_CONTAINER_ID, atm%mesh )
        is_setup_bl_vars_container = .true.

        call atm%vars%container_list(BL_VARS_CONTAINER_ID)%Preproc_operation_for_phys( atm%vars%container, atm%dyn_proc%dyncore_driver )
    end if
    call atm%vars%Get_container( BL_VARS_CONTAINER_ID, vars_container_BL )

    !-- Boundary layer mixing
    call USER_sub_BLmixing_calc_tendency( &
      vars_container%PHY_TEND(MOMX_p), vars_container%PHY_TEND(MOMY_p), &
      vars_container%PHY_TEND(RHOT_p), vars_container%PHY_TEND(PHYTEND_NUM1+1), &
      vars_container, vars_container_BL, &
      atm%mesh%DOptrMat(3), atm%mesh%LiftOptrMat, atm%mesh%ptr_mesh )

    return
  end subroutine USER_calc_tendency

  !> Determine the reference temperature profile T_ref(z) before the first time step
!OCL SERIAL
  subroutine USER_update_pre( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm
    !----------------------------------------------------------

    if ( .not. this%USER_do ) return
    if ( is_Tref_ready ) return

    call setup_Tref( atm )
    is_Tref_ready = .true.

    return
  end subroutine USER_update_pre

  !> Apply Newtonian relaxation of temperature, Rayleigh damping in the PBL, and the sponge layer
  !! These are treated with the implicit (backward Euler) method after the dynamical step.
!OCL SERIAL
  subroutine USER_update( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: THERM_hyd
    class(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke2D, ke_z
    integer :: p, p_v, Nh2D

    real(RP), allocatable :: DENS(:), T(:), Tref(:), PREref(:), PTeq(:), sig(:), PRES_sfc(:)
    real(RP), allocatable :: rtauT(:), rtauV(:), sfac(:), zlev(:)

    real(RP) :: dom_zmax
    real(DP) :: dt
    real(RP) :: gamm
    !----------------------------------------------------------

    if ( .not. this%USER_do ) return

    dt = atm%time_manager%dtsec
    gamm = CpDry / CvDry

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )

      call AtmosVars_GetLocalMeshPhyAuxVars( n, atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                         )

      !- Hydrostatic part of RHOT (RHOT = THERM_hyd + DRHOT)
      THERM_hyd => atm%vars%container%AUX_VARS(AUXVAR_THERMHYDRO_ID)%local(n)

      elem3D => lcmesh%refElem3D
      Nh2D = elem3D%Nnode_h1D**2
      dom_zmax = lcmesh%zmax

      allocate( DENS(elem3D%Np), T(elem3D%Np), Tref(elem3D%Np), PREref(elem3D%Np), PTeq(elem3D%Np), sig(elem3D%Np) )
      allocate( PRES_sfc(Nh2D) )
      allocate( rtauT(elem3D%Np), rtauV(elem3D%Np), sfac(elem3D%Np), zlev(elem3D%Np) )

      !$omp parallel do private(ke2D, ke_z, p, p_v, DENS, T, Tref, PREref, PTeq, PRES_sfc, sig, rtauT, rtauV, sfac, zlev)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)
        ke_z = (ke - 1) / lcmesh%Ne2D + 1
        zlev(:) = lcmesh%zlev(:,ke)

        do p=1, elem3D%Np
          p_v = (p - 1) / Nh2D + 1
          Tref(p) = Tref_z(p_v,ke_z)
          PREref(p) = PREref_z(p_v,ke_z)
        end do

        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)
        T(:) =  PRES%val(:,ke) / ( Rtot%val(:,ke) * DENS(:) )

        !- Damping rate of the sponge layer near the model top: rtauS (stored in sfac)
        sfac(:) = 0.0_RP
        if ( TAU_SPONGE > 0.0_RP ) then
          where ( zlev(:) > SPONGE_HEIGHT )
            sfac(:) = 1.0_RP / TAU_SPONGE &
                    * 0.5_RP * ( 1.0_RP - cos( PI * ( zlev(:) - SPONGE_HEIGHT ) / ( dom_zmax - SPONGE_HEIGHT ) ) )
          end where
        end if

        !- Damping rate of momentum in the PBL: rtauV
        rtauV(:) = 0.0_RP
        if ( PBL_MOM_TAU > 0.0_RP ) then
          select case( trim(PBL_MOM_TYPE) )
          case( 'HEIGHT' )
            rtauV(:) = ( 1.0_RP - trans_func( zlev(:), ENV_PBL_Z, ENV_PBL_W ) ) / PBL_MOM_TAU
          case( 'SIGMA' )
            PRES_sfc(:) = PRES%val(elem3D%Hslice(:,1),ke2D)
            sig(:) = PRES%val(:,ke) / PRES_sfc(elem3D%IndexH2Dto3D)
            rtauV(:) = max( 0.0_RP, (sig(:) - PBL_MOM_SIGB)/(1.0_RP - PBL_MOM_SIGB) ) / PBL_MOM_TAU
          end select
        end if

        !- Newtonian relaxation of temperature
        !  As in the experiments on Fugaku, the sponge rate is added to the relaxation rate (SPONGE_DAMP_TEMP).
        rtauT(:) = 1.0_RP / calc_tau_rad( zlev(:), elem3D%Np )
        if ( SPONGE_DAMP_TEMP ) rtauT(:) = rtauT(:) + sfac(:)

        if ( relax_theta ) then
          !- 'THETA' (Fugaku form): rho*theta is relaxed toward rho*theta_eq at fixed density, where
          !  theta_eq = T_ref (P00/PRE_ref)^(Rtot/CPtot) is a fixed function of height.
          !  Backward Euler: RHOT^{n+1} = ( RHOT^n + dt/tau * DENS * theta_eq ) / ( 1 + dt/tau ), RHOT = THERM_hyd + DRHOT
          PTeq(:) = Tref(:) * ( PRES00 / PREref(:) )**( Rtot%val(:,ke) / CPtot%val(:,ke) )
          DRHOT%val(:,ke) = ( DRHOT%val(:,ke) + dt * rtauT(:) * ( DENS(:) * PTeq(:) - THERM_hyd%val(:,ke) ) ) &
                          / ( 1.0_RP + dt * rtauT(:) )
        else
          !- 'HS' (moist Held-Suarez form): T is relaxed toward T_ref at fixed density.
          DRHOT%val(:,ke) = DRHOT%val(:,ke) &
                          - dt * rtauT(:) * ( 1.0_RP - Tref(:) / T(:) ) * DENS(:) * PT%val(:,ke)     &
                          / ( 1.0_RP + dt * rtauT(:) * ( 1.0_RP + (gamm - 1.0_RP) * Tref(:) / T(:) ) )
        end if

        !- Rayleigh damping of momentum (PBL + sponge), treated with a single implicit step as on Fugaku
        !  MOMZ is damped by the PBL term only for PBL_MOM_TYPE = 'HEIGHT' (the Fugaku form).
        MOMX%val(:,ke) = ( MOMX%val(:,ke) + dt * ( rtauV(:) + sfac(:) ) * DENS(:) * ENV_U ) / ( 1.0_RP + dt * ( rtauV(:) + sfac(:) ) )
        MOMY%val(:,ke) = ( MOMY%val(:,ke) + dt * ( rtauV(:) + sfac(:) ) * DENS(:) * ENV_V ) / ( 1.0_RP + dt * ( rtauV(:) + sfac(:) ) )
        if ( trim(PBL_MOM_TYPE) == 'HEIGHT' ) then
          MOMZ%val(:,ke) = MOMZ%val(:,ke) / ( 1.0_RP + dt * ( rtauV(:) + sfac(:) ) )
        else
          MOMZ%val(:,ke) = MOMZ%val(:,ke) / ( 1.0_RP + dt * sfac(:) )
        end if
      end do

      deallocate( DENS, T, Tref, PREref, PTeq, sig, PRES_sfc )
      deallocate( rtauT, rtauV, sfac, zlev )
    end do
    call atm%vars%Calc_diagnostics()

    return
  end subroutine USER_update

  !------

  !> Read the namelist PARAM_EXP
!OCL SERIAL
  subroutine read_param_exp( caller )
    implicit none
    character(*), intent(in) :: caller
    integer :: ierr
    !------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO(caller,*) 'Not found namelist PARAM_EXP. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR(caller,*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    return
  end subroutine read_param_exp

  !> Smooth transition function S(z; z0, w) = 0.5 * [ 1 + tanh( (z - z0) / w ) ]
!OCL SERIAL
  elemental function trans_func( z, z0, w ) result( s )
    implicit none
    real(RP), intent(in) :: z, z0, w
    real(RP) :: s
    !------------------------------------------
    s = 0.5_RP * ( 1.0_RP + tanh( ( z - z0 ) / w ) )
    return
  end function trans_func

  !> Integral of S from 0 to z: I(z; z0, w) = z/2 + (w/2) log[ cosh((z-z0)/w) / cosh(z0/w) ]
  !! The form of log(cosh(a)) - log(cosh(b)) is evaluated with a numerically stable formula.
!OCL SERIAL
  elemental function trans_func_integ( z, z0, w ) result( val )
    implicit none
    real(RP), intent(in) :: z, z0, w
    real(RP) :: val
    !------------------------------------------
    val = 0.5_RP * z + 0.5_RP * w * ( log_cosh( (z - z0)/w ) - log_cosh( z0/w ) )
    return
  end function trans_func_integ

!OCL SERIAL
  elemental function log_cosh( x ) result( val )
    implicit none
    real(RP), intent(in) :: x
    real(RP) :: val
    !------------------------------------------
    val = abs(x) + log( 0.5_RP * ( 1.0_RP + exp( - 2.0_RP * abs(x) ) ) )
    return
  end function log_cosh

  !> Potential temperature of the basic state theta_ref(z)
!OCL SERIAL
  elemental function calc_theta_ref( z ) result( pt )
    implicit none
    real(RP), intent(in) :: z
    real(RP) :: pt
    !------------------------------------------
    pt = SFC_THETA + ENV_PT_LAPS_BL * z                                                &
       + ( ENV_PT_LAPS_FT  - ENV_PT_LAPS_BL ) * trans_func_integ( z, ENV_PBL_Z, ENV_PBL_W ) &
       + ( ENV_PT_LAPS_STR - ENV_PT_LAPS_FT ) * trans_func_integ( z, ENV_TRO_Z, ENV_TRO_W )
    return
  end function calc_theta_ref

  !> Relative humidity of the basic state RH_ref(z) (ratio)
!OCL SERIAL
  elemental function calc_rh_ref( z ) result( rh )
    implicit none
    real(RP), intent(in) :: z
    real(RP) :: rh
    !------------------------------------------
    rh = ENV_RH_UP                                                              &
       + ( ENV_RH_MID - ENV_RH_UP  ) * ( 1.0_RP - trans_func( z, ENV_MID_Z, ENV_MID_W ) ) &
       + ( ENV_RH_SFC - ENV_RH_MID ) * ( 1.0_RP - trans_func( z, ENV_PBL_Z, ENV_PBL_W ) )
    return
  end function calc_rh_ref

  !> Relaxation time of temperature tau_rad(z)
!OCL SERIAL
  function calc_tau_rad( z, n ) result( tau )
    implicit none
    integer, intent(in) :: n
    real(RP), intent(in) :: z(n)
    real(RP) :: tau(n)
    !------------------------------------------
    tau(:) = TAU_RAD_BL                                                     &
           + ( TAU_RAD_FT  - TAU_RAD_BL ) * trans_func( z(:), ENV_PBL_Z, ENV_PBL_W ) &
           + ( TAU_RAD_STR - TAU_RAD_FT ) * trans_func( z(:), ENV_TRO_Z, ENV_TRO_W )
    return
  end function calc_tau_rad

  !> Set the reference temperature profile T_ref(z) and the reference pressure profile PRE_ref(z)
  !!
  !! T_ref(z) (PRE_ref(z)) is the horizontal mean of temperature (pressure) at the initial time:
  !!  - If in_filebase_Tref is specified, 'T' (and 'PRES' for TREF_RELAX_TYPE = 'THETA') at the time step Tref_step
  !!    is read from the history file.
  !!    The history file must have been written with the same mesh and the same process decomposition.
  !!  - Otherwise, the current model state (the initial state of this run) is used.
  !! Since the mesh is horizontally uniform with a flat bottom, all columns share the same vertical node positions.
!OCL SERIAL
  subroutine setup_Tref( atm )
    use scale_file_base_meshfield, only: FILE_base_meshfield
    use scale_mesh_cubedom3d, only: MeshCubeDom3D
    use scale_mesh_base3d, only: &
      MF3D_XYZT => MeshBase3D_DIMTYPEID_XYZT
    implicit none
    class(AtmosComponent), intent(inout) :: atm

    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lcmesh
    class(ElementBase3D), pointer :: elem3D
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot

    type(FILE_base_meshfield) :: file
    type(MeshField3D) :: T_in
    type(MeshField3D) :: P_in
    logical :: use_file

    integer :: n, ke, ke_z, p, p_v, Nh2D, Nnode_v, NeZ
    real(DP), allocatable :: sum_wT(:,:), sum_wP(:,:), sum_w(:,:)
    real(DP), allocatable :: sum_wT_g(:,:), sum_wP_g(:,:), sum_w_g(:,:)
    real(RP), allocatable :: T(:), Pr(:)
    real(RP) :: w
    integer :: ierr
    !------------------------------------------

    mesh3D => atm%mesh%ptr_mesh
    lcmesh => mesh3D%lcmesh_list(1)
    elem3D => lcmesh%refElem3D
    Nh2D    = elem3D%Nnode_h1D**2
    Nnode_v = elem3D%Nnode_v
    NeZ     = lcmesh%NeZ

    use_file = ( trim(in_filebase_Tref) /= '' )

    !- T_in is always allocated. NOTE: Referencing an unallocated field even in a branch which is not executed
    !  can cause a segmentation fault with aggressive optimization (e.g., -Kpreex of the Fujitsu compiler).
    call T_in%Init( "T_in", "K", mesh3D )
    call P_in%Init( "P_in", "Pa", mesh3D )
    do n=1, mesh3D%LOCAL_MESH_NUM
      T_in%local(n)%val(:,:) = 0.0_RP
      P_in%local(n)%val(:,:) = 0.0_RP
    end do

    if ( use_file ) then
      LOG_INFO("USER_setup_Tref",*) 'T_ref is computed from the history file: ', trim(in_filebase_Tref), ', step=', Tref_step
      select type( mesh3D )
      type is ( MeshCubeDom3D )
        call file%Init( 1, mesh3D=mesh3D )
      class default
        LOG_ERROR("USER_setup_Tref",*) 'This experiment supports only the regional (CubeDom3D) mesh. Check!'
        call PRC_abort
      end select
      call file%Open( in_filebase_Tref, myrank=PRC_myrank )
      call file%Read_Var( MF3D_XYZT, 'T', T_in, step=Tref_step )
      if ( relax_theta ) then
        call file%Read_Var( MF3D_XYZT, 'PRES', P_in, step=Tref_step )
      end if
      call file%Close()
      call file%Final()
    else
      LOG_INFO("USER_setup_Tref",*) 'T_ref is computed from the initial state of this run.'
      do n=1, mesh3D%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshPrgVars( n, mesh3D,                 &
          atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,         &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                              &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh               )
        call AtmosVars_GetLocalMeshPhyAuxVars( n, mesh3D,              &
          atm%vars%AUXVARS_manager, PRES, PT                           )
        do ke=lcmesh%NeS, lcmesh%NeE
          T_in%local(n)%val(:,ke) = PRES%val(:,ke) / ( Rtot%val(:,ke) * ( DENS_hyd%val(:,ke) + DDENS%val(:,ke) ) )
          P_in%local(n)%val(:,ke) = PRES%val(:,ke)
        end do
      end do
    end if

    allocate( sum_wT(Nnode_v,NeZ), sum_wP(Nnode_v,NeZ), sum_w(Nnode_v,NeZ) )
    allocate( sum_wT_g(Nnode_v,NeZ), sum_wP_g(Nnode_v,NeZ), sum_w_g(Nnode_v,NeZ) )
    sum_wT(:,:) = 0.0_DP
    sum_wP(:,:) = 0.0_DP
    sum_w (:,:) = 0.0_DP

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      elem3D => lcmesh%refElem3D

      allocate( T(elem3D%Np), Pr(elem3D%Np) )
      do ke=lcmesh%NeS, lcmesh%NeE
        ke_z = (ke - 1) / lcmesh%Ne2D + 1
        T(:) = T_in%local(n)%val(:,ke)
        Pr(:) = P_in%local(n)%val(:,ke)
        do p=1, elem3D%Np
          p_v = (p - 1) / Nh2D + 1
          w = elem3D%IntWeight_lgl(p) * lcmesh%J(p,ke) * lcmesh%Gsqrt(p,ke)
          sum_wT(p_v,ke_z) = sum_wT(p_v,ke_z) + real( w * T(p), kind=DP )
          sum_wP(p_v,ke_z) = sum_wP(p_v,ke_z) + real( w * Pr(p), kind=DP )
          sum_w (p_v,ke_z) = sum_w (p_v,ke_z) + real( w, kind=DP )
        end do
      end do
      deallocate( T, Pr )
    end do

    call MPI_Allreduce( sum_wT, sum_wT_g, Nnode_v*NeZ, MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_Allreduce( sum_wP, sum_wP_g, Nnode_v*NeZ, MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_Allreduce( sum_w,  sum_w_g,  Nnode_v*NeZ, MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr )

    if ( allocated(Tref_z) ) deallocate( Tref_z )
    allocate( Tref_z(Nnode_v,NeZ) )
    Tref_z(:,:) = real( sum_wT_g(:,:) / sum_w_g(:,:), kind=RP )

    !- PRE_ref is always allocated (used only for TREF_RELAX_TYPE = 'THETA').
    !  For 'HS' with in_filebase_Tref, PRES is not read and PRE_ref is set to PRES00 (dummy, not used).
    if ( allocated(PREref_z) ) deallocate( PREref_z )
    allocate( PREref_z(Nnode_v,NeZ) )
    PREref_z(:,:) = real( sum_wP_g(:,:) / sum_w_g(:,:), kind=RP )
    if ( use_file .and. ( .not. relax_theta ) ) PREref_z(:,:) = PRES00

    LOG_INFO("USER_setup_Tref",*) 'TREF_RELAX_TYPE = ', trim(TREF_RELAX_TYPE)
    LOG_INFO("USER_setup_Tref",*) 'T_ref at the vertical nodes (ke_z, p_v, T_ref [K], PRE_ref [Pa], theta_eq(dry) [K]):'
    do ke_z=1, NeZ
    do p_v=1, Nnode_v
      LOG_INFO_CONT('(2I4,F10.3,F14.3,F10.3)') ke_z, p_v, Tref_z(p_v,ke_z), PREref_z(p_v,ke_z), &
        Tref_z(p_v,ke_z) * ( PRES00 / PREref_z(p_v,ke_z) )**( Rdry / CPdry )
    end do
    end do

    call T_in%Final()
    call P_in%Final()
    deallocate( sum_wT, sum_wP, sum_w, sum_wT_g, sum_wP_g, sum_w_g )

    return
  end subroutine setup_Tref

  !------

  !> Set the initial condition
  !!
  !! 1. Build a dry hydrostatic state with theta_ref(z).
  !! 2. Diagnose QV from RH_ref(z) with the temperature and pressure of the current estimate,
  !!    and rebuild a moist hydrostatic state (iterated NITER_QV times).
  !! 3. Add random perturbation of potential temperature near the surface.
!OCL SERIAL
  subroutine exp_SetInitCond_idealized_moist_conv( this,                   &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )

    use scale_const, only: &
      EPSvap => CONST_EPSvap
    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_psat_all
    use scale_atmos_hydrometeor, only: &
      CV_VAPOR, CP_VAPOR
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constPT, &
      hydrostatic_build_rho_XYZ
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

    integer, parameter :: NITER_QV = 3

    real(RP) :: PT_tmp(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: Rtot  (elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: CPtot (elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: CVtot (elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: CPtot_ov_CVtot(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: QV    (elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: TEMP  (elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: PRES_ (elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)

    real(RP) :: DENS(elem%Np)
    real(RP) :: PT(elem%Np)
    real(RP) :: rh, psat, e
    real(RP) :: rnd(elem%Np)

    real(RP), allocatable :: bnd_SFC_PRES(:,:)
    real(RP) :: sfc_rhot(elem%Nfp_v)

    integer :: ke, ke2D, ke_x, ke_y, ke_z
    integer :: p
    integer :: iq_QV
    integer :: itr
    integer :: nseed
    integer, allocatable :: seed(:)
    !-----------------------------------------------------------------------------

    call read_param_exp( "IDEALIZED_MOIST_CONV_setup" )

    !-- Step 1: dry hydrostatic state ------------------------------------

    call hydrostatic_calc_basicstate_constPT( DENS_hyd, PRES_hyd,                            &
      SFC_THETA, SFC_PRES, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3), &
      lcmesh, elem )

    allocate( bnd_SFC_PRES(elem%Nfp_v,lcmesh%Ne2DA) )

    !$omp parallel do collapse(3) private( ke, ke2D, ke_x, ke_y, ke_z, sfc_rhot )
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
    do ke_z=1, lcmesh%NeZ
      ke2D = ke_x + (ke_y-1)*lcmesh%NeX
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D

      PT_tmp(:,ke_z,ke_x,ke_y) = calc_theta_ref( lcmesh%zlev(:,ke) )

      if ( ke_z == 1 ) then
        sfc_rhot(:) = DENS_hyd(elem%Hslice(:,1),ke2D) * PT_tmp(elem%Hslice(:,1),ke_z,ke_x,ke_y)
        bnd_SFC_PRES(:,ke2D) = PRES00 * ( Rdry * sfc_rhot(:) / PRES00 )**( CpDry / CVdry )
      end if
    end do
    end do
    end do

    call hydrostatic_build_rho_XYZ( DDENS, & ! (out)
      DENS_hyd, PRES_hyd, PT_tmp,          & ! (in)
      x, y, z, lcmesh, elem, bnd_SFC_PRES  ) ! (in)

    !$omp parallel do collapse(3) private( ke, ke_x, ke_y, ke_z, DENS )
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
    do ke_z=1, lcmesh%NeZ
      ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1)*lcmesh%Ne2D

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      PRES_hyd(:,ke) = PRES00 * ( Rdry / PRES00 * DENS(:) * PT_tmp(:,ke_z,ke_x,ke_y) )**( CpDry / CvDry )
      DENS_hyd(:,ke) = DENS(:)
      DDENS(:,ke) = 0.0_RP

      TEMP (:,ke_z,ke_x,ke_y) = PRES_hyd(:,ke) / ( Rdry * DENS_hyd(:,ke) )
      PRES_(:,ke_z,ke_x,ke_y) = PRES_hyd(:,ke)
    end do
    end do
    end do

    !-- Step 2: moist hydrostatic state ----------------------------------

    call TRACER_inq_id( "QV", iq_QV )

    do itr=1, NITER_QV

      !$omp parallel do collapse(3) private( ke, ke_x, ke_y, ke_z, p, rh, psat, e )
      do ke_y=1, lcmesh%NeY
      do ke_x=1, lcmesh%NeX
      do ke_z=1, lcmesh%NeZ
        ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1)*lcmesh%Ne2D

        do p=1, elem%Np
          rh = calc_rh_ref( lcmesh%zlev(p,ke) )
          call ATMOS_SATURATION_psat_all( TEMP(p,ke_z,ke_x,ke_y), psat )
          e = rh * psat
          ! specific humidity of air without condensates
          QV(p,ke_z,ke_x,ke_y) = EPSvap * e / ( PRES_(p,ke_z,ke_x,ke_y) - ( 1.0_RP - EPSvap ) * e )
        end do

        Rtot (:,ke_z,ke_x,ke_y) = Rdry  * ( 1.0_RP - QV(:,ke_z,ke_x,ke_y) ) + Rvap     * QV(:,ke_z,ke_x,ke_y)
        CPtot(:,ke_z,ke_x,ke_y) = CPdry * ( 1.0_RP - QV(:,ke_z,ke_x,ke_y) ) + CP_VAPOR * QV(:,ke_z,ke_x,ke_y)
        CVtot(:,ke_z,ke_x,ke_y) = CVdry * ( 1.0_RP - QV(:,ke_z,ke_x,ke_y) ) + CV_VAPOR * QV(:,ke_z,ke_x,ke_y)
        CPtot_ov_CVtot(:,ke_z,ke_x,ke_y) = CPtot(:,ke_z,ke_x,ke_y) / CVtot(:,ke_z,ke_x,ke_y)
      end do
      end do
      end do

      call hydrostatic_build_rho_XYZ( DDENS,  & ! (out)
        DENS_hyd, PRES_hyd, PT_tmp,           & ! (in)
        Rtot, CPtot_ov_CVtot,                 & ! (in)
        x, y, z, lcmesh, elem, bnd_SFC_PRES   ) ! (in)

      !- Update temperature and pressure for the next iteration
      !$omp parallel do collapse(3) private( ke, ke_x, ke_y, ke_z, DENS )
      do ke_y=1, lcmesh%NeY
      do ke_x=1, lcmesh%NeX
      do ke_z=1, lcmesh%NeZ
        ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1)*lcmesh%Ne2D

        DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
        PRES_(:,ke_z,ke_x,ke_y) = PRES00 &
          * ( Rtot(:,ke_z,ke_x,ke_y) / PRES00 * DENS(:) * PT_tmp(:,ke_z,ke_x,ke_y) )**CPtot_ov_CVtot(:,ke_z,ke_x,ke_y)
        TEMP (:,ke_z,ke_x,ke_y) = PRES_(:,ke_z,ke_x,ke_y) / ( Rtot(:,ke_z,ke_x,ke_y) * DENS(:) )
      end do
      end do
      end do
    end do

    !-- Step 3: set prognostic variables with random perturbation ---------

    call random_seed( size=nseed )
    allocate( seed(nseed) )
    do p=1, nseed
      seed(p) = RANDOM_ISEED + 7919 * PRC_myrank + 104729 * p
    end do
    call random_seed( put=seed )
    deallocate( seed )

    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
    do ke_z=1, lcmesh%NeZ
      ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1)*lcmesh%Ne2D

      call random_number( rnd )

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      PT(:) = PT_tmp(:,ke_z,ke_x,ke_y)
      where ( lcmesh%zlev(:,ke) < RANDOM_THETA_ZTOP )
        PT(:) = PT(:) + RANDOM_THETA * ( 2.0_RP * rnd(:) - 1.0_RP )
      end where

      DRHOT(:,ke) = DENS(:) * PT(:) &
                  - PRES00 / Rdry * ( PRES_hyd(:,ke) / PRES00 )**( CVdry / CPdry )
      MOMX(:,ke) = DENS(:) * ENV_U
      MOMY(:,ke) = DENS(:) * ENV_V
      MOMZ(:,ke) = 0.0_RP

      tracer_field_list(iq_QV)%ptr%val(:,ke) = QV(:,ke_z,ke_x,ke_y)
    end do
    end do
    end do

    deallocate( bnd_SFC_PRES )

    return
  end subroutine exp_SetInitCond_idealized_moist_conv

end module mod_user
