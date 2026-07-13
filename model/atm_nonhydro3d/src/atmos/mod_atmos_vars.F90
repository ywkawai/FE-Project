!-------------------------------------------------------------------------------
!> module Atmosphere / Variables
!!
!! @par Description
!!          Module to manage variables with atmospheric component
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_debug
  use scale_tracer, only: QA

  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: &
    MeshBase2D
  use scale_mesh_base3d, only: &
    MeshBase3D,                              &
    DIMTYPE_XYZ  => MeshBase3D_DIMTYPEID_XYZ

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField2D, MeshField3D, MeshField3DList
  
  use scale_file_restart_meshfield, only: &
    FILE_restart_meshfield_component
  
  use scale_meshfieldcomm_base, only: MeshFieldContainer

  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    PRGVAR_NUM, AUXVAR_NUM, PHYTEND_NUM1 => PHYTEND_NUM,                                        &
    PRGVAR_DDENS_ID, PRGVAR_THERM_ID, PRGVAR_MOMZ_ID, PRGVAR_MOMX_ID, PRGVAR_MOMY_ID,           &
    AUXVAR_DENSHYDRO_ID, AUXVAR_PRESHYDRO_ID, AUXVAR_THERMHYDRO_ID, AUXVAR_PRESHYDRO_REF_ID,    &
    AUXVAR_Rtot_ID, AUXVAR_CPtot_ID, AUXVAR_CVtot_ID,                                           &
    AUXVAR_PRES_ID, AUXVAR_PT_ID, AUXVAR_Qdry_ID,                                               &
    PHYTEND_DENS_ID, PHYTEND_MOMX_ID, PHYTEND_MOMY_ID, PHYTEND_MOMZ_ID, PHYTEND_RHOT_ID,        &
    PHYTEND_RHOH_ID    

  use mod_atmos_mesh, only: AtmosMesh

  use mod_atmos_vars_container, only: &
    AtmosVarsContainer, &
    AtmosVars_GetLocalMeshPrgVar, AtmosVars_GetLocalMeshPrgVars,       &
    AtmosVars_GetLocalMeshSfcVar,                                      &
    AtmosVars_GetLocalMeshQTRCVar, AtmosVars_GetLocalMeshQTRCVarList,  &
    AtmosVars_GetLocalMeshQTRC_Qv,                                     &
    AtmosVars_GetLocalMeshPhyAuxVars,                                  &
    AtmosVars_GetLocalMeshPhyTends, AtmosVars_GetLocalMeshQTRCPhyTend, &
    ATMOS_AUXVARS2D_NUM,                                               &
    ATM_VARS_CONTAINER_PRIMARY_ID
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  integer, parameter :: ATM_VARS_CONTAINER_LIST_MAX = 16

  !> Derived type to manage variables with atmospheric component
  type, public :: AtmosVars
    class(AtmosMesh), pointer :: mesh                !< Pointer to an object to manage a primary mesh for atmospheric component

    type(AtmosVarsContainer) :: container_list(ATM_VARS_CONTAINER_LIST_MAX) !< Array of containers to manage variables. The first one is the primary container.
    logical :: init_containers_item_flag(ATM_VARS_CONTAINER_LIST_MAX)       !< Flag whether the item of container_list is initialized or not. 
    
    type(AtmosVarsContainer), pointer :: container              !< Pointer to a primary container to manage variables
    type(ModelVarManager), pointer :: PROGVARS_manager          !< Pointer to manage prognostic variables with primary container
    type(ModelVarManager), pointer :: QTRCVARS_manager          !< Pointer to manage tracer variables with primary container
    type(ModelVarManager), pointer :: AUXVARS_manager           !< Pointer to manage auxiliary variables with primary container
    type(ModelVarManager), pointer :: AUXVARS2D_manager         !< Pointer to manage 2D auxiliary variables with primary container
    type(ModelVarManager), pointer :: PHYTENDS_manager          !< Pointer to manage tendency variables with physics
   
    !-
    type(ModelVarManager), pointer :: ptr_MP_AUXVARS2D_manager

    logical :: moist
    type(MeshField3D), pointer :: QV
    type(MeshField3D) :: zero

    !--
    integer, allocatable :: DIAGVARS2D_HISTID(:)
    integer, allocatable :: DIAGVARS3D_HISTID(:)

    type(FILE_restart_meshfield_component) :: restart_file !< Object to manage restart file for atmospheric component
    
    character(len=H_MID) :: phy_preproc_file_basename !< Basename of configuration file for preprocesses before physics

    logical :: check_range
    logical :: check_total

  contains
    procedure :: Init => AtmosVars_Init
    procedure :: Final => AtmosVars_Final
    procedure :: Setup_container => AtmosVars_Setup_container   
    procedure :: Get_container => AtmosVars_Get_container 
    procedure :: Calc_diagnostics => AtmosVars_CalculateDiagnostics
    procedure :: Calc_diagVar => AtmosVars_CalcDiagvar
    procedure :: Calc_diagVar2D => AtmosVars_CalcDiagvar2D
    procedure :: PreprocOperationForPhys => AtmosVars_PreprocOperationForPhys
    procedure :: History => AtmosVars_History
    procedure :: Check   => AtmosVars_Check
    procedure :: Monitor => AtmosVars_Monitor
    procedure :: Read_restart_file => AtmosVar_Read_restart_file
    procedure :: Write_restart_file_prep => AtmosVar_Write_restart_file_prep
    procedure :: Write_restart_file => AtmosVar_Write_restart_file
    procedure :: Regist_physvar_manager => AtmosVars_Regist_physvar_manager
  end type AtmosVars

  ! Cascade
  public :: AtmosVarsContainer
  public :: ATM_VARS_CONTAINER_PRIMARY_ID
  public :: AtmosVars_GetLocalMeshPrgVar
  public :: AtmosVars_GetLocalMeshPrgVars
  public :: AtmosVars_GetLocalMeshSfcVar
  public :: AtmosVars_GetLocalMeshQTRCVar  
  public :: AtmosVars_GetLocalMeshQTRCVarList
  public :: AtmosVars_GetLocalMeshQTRC_Qv
  public :: AtmosVars_GetLocalMeshPhyAuxVars
  public :: AtmosVars_GetLocalMeshPhyTends
  public :: AtmosVars_GetLocalMeshQTRCPhyTend

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), parameter :: PROGVARS_check_min(PRGVAR_NUM) = (/ -1.0_RP, -100.0_RP, -200.0_RP, -200.0_RP, -200.0_RP /)
  real(RP), parameter :: PROGVARS_check_max(PRGVAR_NUM) = (/  1.0_RP,   100.0_RP, 200.0_RP,  200.0_RP,  200.0_RP /)
                      

  ! Diagnostic variables

  integer, public, parameter :: ATMOS_DIAGVARS_DENS_ID   = 1
  integer, public, parameter :: ATMOS_DIAGVARS_U_ID      = 2
  integer, public, parameter :: ATMOS_DIAGVARS_V_ID      = 3
  integer, public, parameter :: ATMOS_DIAGVARS_W_ID      = 4
  integer, public, parameter :: ATMOS_DIAGVARS_T_ID      = 5  
  integer, public, parameter :: ATMOS_DIAGVARS_Umet_ID   = 6
  integer, public, parameter :: ATMOS_DIAGVARS_Vmet_ID   = 7
  integer, public, parameter :: ATMOS_DIAGVARS_QDRY_ID   = 8
  integer, public, parameter :: ATMOS_DIAGVARS_RH_ID     = 9
  integer, public, parameter :: ATMOS_DIAGVARS_ENGT_ID   = 10
  integer, public, parameter :: ATMOS_DIAGVARS_ENGP_ID   = 11
  integer, public, parameter :: ATMOS_DIAGVARS_ENGK_ID   = 12
  integer, public, parameter :: ATMOS_DIAGVARS_ENGI_ID   = 13
  integer, public, parameter :: ATMOS_DIAGVARS3D_NUM     = 13

  type(VariableInfo), public :: ATMOS_DIAGVARS3D_VINFO(ATMOS_DIAGVARS3D_NUM)
  DATA ATMOS_DIAGVARS3D_VINFO / &
    VariableInfo( ATMOS_DIAGVARS_DENS_ID, 'DENS',  'density'              , 'kg/m3', 3, 'XYZ', 'air_density'          ), &
    VariableInfo( ATMOS_DIAGVARS_U_ID   , 'U'   ,  'velocity u'           , 'm/s'  , 3, 'XYZ', 'x_wind'               ), &
    VariableInfo( ATMOS_DIAGVARS_V_ID   , 'V'   ,  'velocity v'           , 'm/s'  , 3, 'XYZ', 'y_wind'               ), &  
    VariableInfo( ATMOS_DIAGVARS_W_ID   , 'W'   ,  'velocity w'           , 'm/s'  , 3, 'XYZ', 'upward_air_velocity'  ), &
    VariableInfo( ATMOS_DIAGVARS_T_ID   , 'T'   ,  'temperature'          , 'K'    , 3, 'XYZ', 'air_temperature'      ), &
    VariableInfo( ATMOS_DIAGVARS_Umet_ID, 'Umet',  'eastward velocity'    , 'm/s'  , 3, 'XYZ', 'x_wind'               ), &
    VariableInfo( ATMOS_DIAGVARS_Vmet_ID, 'Vmet',  'northward velocity'   , 'm/s'  , 3, 'XYZ', 'y_wind'               ), &  
    Variableinfo( ATMOS_DIAGVARS_QDRY_ID, 'QDRY',  'dry air'              , 'kg/kg', 3, 'XYZ', ''                     ), &
    Variableinfo( ATMOS_DIAGVARS_RH_ID  ,   'RH', 'relative humidity(liq)',     '%', 3, 'XYZ', 'relative_humidity'    ), &
    Variableinfo( ATMOS_DIAGVARS_ENGT_ID, 'ENGT',  'total energy'         , 'J/m3' , 3, 'XYZ', ''                     ), &
    Variableinfo( ATMOS_DIAGVARS_ENGP_ID, 'ENGP',  'potential energy'     , 'J/m3' , 3, 'XYZ', ''                     ), &
    Variableinfo( ATMOS_DIAGVARS_ENGK_ID, 'ENGK',  'kinetic energy'       , 'J/m3' , 3, 'XYZ', ''                     ), &
    Variableinfo( ATMOS_DIAGVARS_ENGI_ID, 'ENGI',   'internal energy'     , 'J/m3' , 3, 'XYZ', ''                     )  /


  integer, public, parameter :: ATMOS_DIAGVARS_RAIN_ID   = 1
  integer, public, parameter :: ATMOS_DIAGVARS_SNOW_ID   = 2
  integer, public, parameter :: ATMOS_DIAGVARS2D_NUM     = 2
  
  type(VariableInfo), public :: ATMOS_DIAGVARS2D_VINFO(ATMOS_DIAGVARS2D_NUM)
  DATA ATMOS_DIAGVARS2D_VINFO / &
    VariableInfo( ATMOS_DIAGVARS_RAIN_ID, 'RAIN', 'surface rain flux'              , 'kg/m2/s', 2, 'XY', 'rainfall_flux'       ), &
    VariableInfo( ATMOS_DIAGVARS_SNOW_ID, 'SNOW', 'surface snow flux'              , 'kg/m2/s', 2, 'XY', 'snowfall_flux'       )  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  
  private :: vars_calc_diagnoseVar2D_lc

  ! for monitor

  integer, private, parameter   :: IM_QDRY         =  1
  integer, private, parameter   :: IM_QTOT         =  2
  integer, private, parameter   :: IM_ENGT         =  3
  integer, private, parameter   :: IM_ENGP         =  4
  integer, private, parameter   :: IM_ENGK         =  5
  integer, private, parameter   :: IM_ENGI         =  6
  integer, private, parameter   :: DVM_nmax        =  6
  integer, private              :: DV_MONIT_id(DVM_nmax)

contains

!> Setup an object to manage variables with atmospheric component
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!!
!OCL SERIAL
  subroutine AtmosVars_Init( this, atm_mesh )
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRGVAR_SCALAR_NUM, PRGVAR_HVEC_NUM    
    use scale_file_monitor_meshfield, only:    &
      MONITOR_reg => FILE_monitor_meshfield_reg
    implicit none

    class(AtmosVars), target, intent(inout) :: this
    class(AtmosMesh), target, intent(inout) :: atm_mesh

    integer :: iv
    logical :: reg_file_hist

    type(ModelVarManager) :: diagvar_manager               ! dummy
    type(MeshField2D) :: diag_vars2D(ATMOS_DIAGVARS2D_NUM) ! dummy
    type(MeshField3D) :: diag_vars3D(ATMOS_DIAGVARS3D_NUM) ! dummy

    logical :: CHECK_RANGE    = .false.  !< Flag whether the range of values is checked
    logical :: CHECK_TOTAL    = .false.

    character(len=H_MID) :: PHY_PREPROC_FILE_BASENAME = 'phy_preoperation' !< Basename of configuration file for preprocesses before physics

    namelist / PARAM_ATMOS_VARS / &
      PHY_PREPROC_FILE_BASENAME, &
      CHECK_RANGE, &
      CHECK_TOTAL

    character(len=H_LONG) :: IN_BASENAME           = ''        !< Basename of the input  file
    logical :: IN_POSTFIX_TIMELABEL                = .false.   !< Add timelabel to the basename of input  file?
    character(len=H_LONG) :: OUT_BASENAME          = ''        !< Basename of the output file
    logical :: OUT_POSTFIX_TIMELABEL               = .true.    !< Add timelabel to the basename of output file?
    character(len=H_MID) :: OUT_TITLE              = ''        !< Title    of the output file
    character(len=H_SHORT) :: OUT_DTYPE            = 'DEFAULT' !< REAL4 or REAL8  

    namelist / PARAM_ATMOS_VARS_RESTART / &
      IN_BASENAME,           &
      IN_POSTFIX_TIMELABEL,  &
      OUT_BASENAME,          &
      OUT_POSTFIX_TIMELABEL, &
      OUT_TITLE,             &
      OUT_DTYPE    
    
    integer :: ierr
    logical :: is_specified

    class(MeshBase3D), pointer :: mesh3D
    class(MeshBase2D), pointer :: mesh2D
    !--------------------------------------------------

    LOG_INFO('AtmosVars_Init',*)

    !- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_vars_setup",*) 'Invalid names in namelist PARAM_ATMOS_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_VARS)

    this%phy_preproc_file_basename = PHY_PREPROC_FILE_BASENAME

    !- Set the pointer of mesh
    this%mesh => atm_mesh
    mesh3D => atm_mesh%ptr_mesh
    call mesh3D%GetMesh2D( mesh2D )

    !- Initialize variables associated with dynamical core 
    !  (prognostic variables, tracer variables, 3D auxiliary variables, and tendencies of physical processes)

    this%init_containers_item_flag(:) = .false.
    call this%Setup_container( 1, atm_mesh )  

    !- Initialize diagnostic variables for output

    ! 3D
    call diagvar_manager%Init()
    allocate( this%DIAGVARS3D_HISTID(ATMOS_DIAGVARS3D_NUM) )

    reg_file_hist = .true.
    do iv = 1, ATMOS_DIAGVARS3D_NUM
      call diagvar_manager%Regist( &
        ATMOS_DIAGVARS3D_VINFO(iv), atm_mesh%ptr_mesh, & ! (in) 
        diag_vars3D(iv), reg_file_hist                   ) ! (out)
      
      this%DIAGVARS3D_HISTID(iv) = diag_vars3D(iv)%hist_id
    end do
    call diagvar_manager%Final()

    ! 2D
    call diagvar_manager%Init()
    allocate( this%DIAGVARS2D_HISTID(ATMOS_DIAGVARS2D_NUM) )

    reg_file_hist = .true.
    do iv = 1, ATMOS_DIAGVARS2D_NUM
      call diagvar_manager%Regist( &
        ATMOS_DIAGVARS2D_VINFO(iv), mesh2D,  & ! (in) 
        diag_vars2D(iv), reg_file_hist       ) ! (out)
      
      this%DIAGVARS2D_HISTID(iv) = diag_vars2D(iv)%hist_id
    end do
    call diagvar_manager%Final()

    !-- Setup information for input/output restart files. 

    is_specified = .true.
    !- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_VARS_RESTART,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("AtmosVars_Init",*) 'Not found namelist PARAM_ATMOS_VARS_RESTART. Default used.'
       is_specified = .false.
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("AtmosVars_Init",*) 'Not appropriate names in namelist PARAM_ATMOS_VARS_RESTART. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_VARS_RESTART)

    if (is_specified) then
      call atm_mesh%Setup_restartfile( this%restart_file,           &
        IN_BASENAME, IN_POSTFIX_TIMELABEL,                          &
        OUT_BASENAME, OUT_POSTFIX_TIMELABEL, OUT_DTYPE, OUT_TITLE,  &
        PRGVAR_NUM + AUXVAR_NUM, ""                                 )
    else
      call atm_mesh%Setup_restartfile( this%restart_file, &
        PRGVAR_NUM + AUXVAR_NUM                           )
    end if

    !-----< monitor output setup >-----
    
    call MONITOR_reg( 'QTOT',         'water mass',            'kg', & ! (in)
                      DV_MONIT_id(IM_QTOT),                          & ! (out)
                      dim_type='ATM3D', is_tendency=.false.          ) ! (in)
    
    call MONITOR_reg( 'ENGT',         'total     energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGT),                          & ! (out)
                      dim_type='ATM3D', is_tendency=.false.          ) ! (in)
    call MONITOR_reg( 'ENGP',         'potential energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGP),                          & ! (out)
                      dim_type='ATM3D', is_tendency=.false.          ) ! (in)
    call MONITOR_reg( 'ENGK',         'kinetic   energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGK),                          & ! (out)
                      dim_type='ATM3D', is_tendency=.false.          ) ! (in)
    call MONITOR_reg( 'ENGI',         'internal  energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGI),                          & ! (out)
                      dim_type='ATM3D', is_tendency=.false.          ) ! (in)


    !-----< check the range of values >-----

    this%check_range = CHECK_RANGE
    this%check_total = CHECK_TOTAL
    LOG_INFO("ATMOS_vars_setup",*) 'Check value range of variables?     : ', CHECK_RANGE
    LOG_INFO("ATMOS_vars_setup",*) 'Check total value of variables?     : ', CHECK_TOTAL

    return
  end subroutine AtmosVars_Init

  !> Setup a container to manage variables with atmospheric component
  !!
!OCL SERIAL
  subroutine AtmosVars_Setup_container( this, container_type, atm_mesh )
    implicit none
    class(AtmosVars), intent(inout), target :: this
    integer, intent(in) :: container_type                   !< Type of container to setup. 1: primary container, >=2: other containers
    class(AtmosMesh), target, intent(inout) :: atm_mesh
    !---------------------------------------------
    
    if ( container_type < ATM_VARS_CONTAINER_PRIMARY_ID .or. container_type > ATM_VARS_CONTAINER_LIST_MAX ) then
      LOG_ERROR("ATMOS_vars_setup_container",*) 'container_type is out of range. Check!', container_type
      call PRC_abort
    end if

    if ( .not. this%init_containers_item_flag(container_type) ) then
      LOG_INFO("ATMOS_vars_setup_container",*) 'container_type: ', container_type

      this%init_containers_item_flag(container_type) = .true.
      call this%container_list(container_type)%Init( container_type, this%phy_preproc_file_basename, atm_mesh )

      if ( container_type == ATM_VARS_CONTAINER_PRIMARY_ID ) then
        this%container => this%container_list(container_type)
        this%PROGVARS_manager => this%container%PROGVARS_manager
        this%QTRCVARS_manager => this%container%QTRCVARS_manager
        this%AUXVARS_manager  => this%container%AUXVARS_manager
        this%AUXVARS2D_manager  => this%container%AUXVARS2D_manager
        this%PHYTENDS_manager  => this%container%PHYTENDS_manager
      end if
    end if
    return    
  end subroutine AtmosVars_Setup_container

  !> Get a container to manage variables with atmospheric component
  !!
!OCL SERIAL
  subroutine AtmosVars_Get_container( this, container_type, &
    container )
    implicit none
    class(AtmosVars), intent(inout), target :: this
    integer, intent(in) :: container_type            !< Type of container to get. 1: primary container, >=2: other containers
    class(AtmosVarsContainer), pointer :: container  !< Pointer to a container to manage variables with atmospheric component
    !---------------------------------------------

    if ( container_type < ATM_VARS_CONTAINER_PRIMARY_ID .or. container_type > ATM_VARS_CONTAINER_LIST_MAX ) then
      LOG_ERROR("ATMOS_vars_get_container",*) 'container_type is out of range. Check!', container_type
      call PRC_abort
    end if

    if ( this%init_containers_item_flag(container_type) ) then
      container => this%container_list(container_type)
    else
      LOG_ERROR("ATMOS_vars_get_container",*) 'Not initialized container. Check!'
      call PRC_abort
    end if
    return
  end subroutine AtmosVars_Get_container

  !> Finalize an object to manage variables with atmospheric component
  !!
!OCL SERIAL
  subroutine AtmosVars_Final( this )
    implicit none
    class(AtmosVars), intent(inout) :: this

    integer :: i    
    !--------------------------------------------------

    LOG_INFO('AtmosVars_Final',*)

    call this%restart_file%Final()

    do i=1, ATM_VARS_CONTAINER_LIST_MAX
      if ( this%init_containers_item_flag(i) ) call this%container_list(i)%Final()
    end do

    deallocate( this%DIAGVARS3D_HISTID )

    return
  end subroutine AtmosVars_Final

!OCL SERIAL
  subroutine AtmosVars_Regist_physvar_manager( this, &
    mp_AUXVARS2D_manager )
    implicit none

    class(AtmosVars), target, intent(inout) :: this
    type(ModelVarManager), intent(in), target :: mp_AUXVARS2D_manager
    !----------------------------------------------

    this%ptr_MP_AUXVARS2D_manager => mp_AUXVARS2D_manager

    return
  end subroutine AtmosVars_Regist_physvar_manager

  !> Put data with atmospheric variables to history file
  !!
!OCL SERIAL
  subroutine AtmosVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosVars), intent(inout), target :: this
  
    integer :: v
    integer :: hst_id

    type(MeshField3D) :: tmp_field3D
    class(MeshBase3D), pointer :: mesh3D

    type(MeshField2D) :: tmp_field2D
    class(MeshBase2D), pointer :: mesh2D
    !-------------------------------------------------------------------------

    mesh3D => this%container%PROG_VARS(1)%mesh
    call mesh3D%GetMesh2D(mesh2D)

    !-
    call this%container%Calc_diagnostics()

    do v = 1, PRGVAR_NUM
      hst_id = this%container%PROG_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%container%PROG_VARS(v) )
    end do

    do v = 1, QA
      hst_id = this%container%QTRC_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%container%QTRC_VARS(v) )
    end do

    do v = 1, AUXVAR_NUM
      hst_id = this%container%AUX_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%container%AUX_VARS(v) )
    end do
    do v = 1, ATMOS_AUXVARS2D_NUM
      hst_id = this%container%AUX_VARS2D(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%container%AUX_VARS2D(v) )
    end do

    do v = 1, this%container%PHYTEND_NUM_TOT
      hst_id = this%container%PHY_TEND(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%container%PHY_TEND(v) )
    end do

    !- Output diagnostic variables

    ! 3D
    call tmp_field3D%Init( "tmp_field", "", mesh3D)
    do v = 1, ATMOS_DIAGVARS3D_NUM
      hst_id = this%DIAGVARS3D_HISTID(v)
      if ( hst_id > 0 ) then
        call AtmosVars_CalcDiagvar( this, ATMOS_DIAGVARS3D_VINFO(v)%NAME, tmp_field3D )
        call FILE_HISTORY_meshfield_put( hst_id, tmp_field3D )
      end if
    end do
    call tmp_field3D%Final()

    ! 2D
    call tmp_field2D%Init( "tmp_field", "", mesh2D)
    do v = 1, ATMOS_DIAGVARS2D_NUM
      hst_id = this%DIAGVARS2D_HISTID(v)
      if ( hst_id > 0 ) then
        call AtmosVars_CalcDiagvar2D( this, ATMOS_DIAGVARS2D_VINFO(v)%NAME, tmp_field2D )
        call FILE_HISTORY_meshfield_put( hst_id, tmp_field2D )
      end if
    end do
    call tmp_field2D%Final()

    return
  end subroutine AtmosVars_history

!> Read data with atmospheric variables from restart file
!!
!OCL SERIAL
  subroutine AtmosVar_Read_restart_file( this, atmos_mesh, dyncore  )

    use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
    use scale_meshfieldcomm_base, only: MeshFieldContainer
    use scale_atm_dyn_dgm_driver_nonhydro3d, only: AtmDynDGMDriver_nonhydro3d

    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      AUXVAR_DENSHYDRO_ID    
    implicit none
    
    class(AtmosVars), intent(inout), target :: this
    class(AtmosMesh), intent(in) :: atmos_mesh
    class(AtmDynDGMDriver_nonhydro3d), intent(inout) :: dyncore

    integer :: iv

    integer :: domid
    integer :: ke
    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lcmesh3D
    class(MeshField3D), pointer :: Phyd_ref
    !---------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOSVar_read_restart_file",*) 'Open restart file (ATMOS) '
        
    !- Open restart file
    call this%restart_file%Open()

    !- Read restart file
    
    do iv=1, PRGVAR_NUM
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%container%PROG_VARS(iv)%varname, &
        this%container%PROG_VARS(iv)                                                      )
    end do
    do iv=1, AUXVAR_DENSHYDRO_ID
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%container%AUX_VARS(iv)%varname, &
        this%container%AUX_VARS(iv)                                                      )
    end do
    do iv=1, QA
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%container%QTRC_VARS(iv)%varname, &
        this%container%QTRC_VARS(iv)                                                      )
    end do

    !- Close restart file
    LOG_INFO("ATMOSVar_read_restart_file",*) 'Close restart file (ATMOS) '
    call this%restart_file%Close()

    !-- Prepare diagnostic variables

    ! Calculate specific heat
    call this%container%Calc_SpecificHeat()

    ! Set a basic state of thermodynamics variable
    call dyncore%update_therm_hyd( this%container%AUXVARS_manager )

    ! Calculate pressure
    call dyncore%calc_pressure( this%container%AUX_VARS(AUXVAR_PRES_ID), &
      this%container%PROGVARS_manager, this%container%AUXVARS_manager    )

    ! Set reference value of hydrostatic pressure
    Phyd_ref => this%container%AUX_VARS(AUXVAR_PRESHYDRO_REF_ID)
    mesh3D => Phyd_ref%mesh
    do domid=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(domid)
      !$omp parallel do
      do ke=lcmesh3D%NeS, lcmesh3D%NeE
        Phyd_ref%local(domid)%val(:,ke) = 0.0_RP
      end do
      !$acc update device(Phyd_ref%local(domid)%val)
    end do

    !-- Check read data
    call this%Check( force = .true. )

    !-- Calculate diagnostic variables
    call this%Calc_diagnostics()   

    !-- Communicate halo data of hydrostatic & diagnostic variables
    call this%container%AUXVARS_manager%MeshFieldComm_Exchange()

    !-- Set horizontal gradient of hydrostatic pressure
    call dyncore%update_phyd_hgrad( this%container%AUX_VARS(AUXVAR_PRESHYDRO_ID), Phyd_ref, &
      mesh3D, atmos_mesh%element3D_operation )

    return
  end subroutine AtmosVar_Read_restart_file

!> Write data with atmospheric variables to restart file
!!
!OCL SERIAL
  subroutine AtmosVar_Write_restart_file_prep( this )
    use scale_tracer, only: &
      TRACER_DESC    
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      atm_dyn_dgm_nonhydro3d_common_get_varinfo
    
    implicit none
    class(AtmosVars), intent(inout) :: this

    type(VariableInfo) :: prgvar_info(PRGVAR_NUM)
    type(VariableInfo) :: auxvar_info(AUXVAR_NUM)

    integer :: iv, rf_vid 
    !---------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("ATMOSVar_Write_restart_file",*) 'Create restart file (ATMOS) '

    !- Check data which will be written to restart file
    call this%Check( force = .true. )

    !- Create restart file
    call this%restart_file%Create()
    call PRC_mpibarrier()

    !- Define variables

    call atm_dyn_dgm_nonhydro3d_common_get_varinfo( prgvar_info, auxvar_info )

    do iv=1, PRGVAR_NUM
      rf_vid = iv
      call this%restart_file%Def_var( this%container%PROG_VARS(iv),  &
        prgvar_info(iv)%DESC, rf_vid, DIMTYPE_XYZ                    )
    end do
    do iv=1, AUXVAR_DENSHYDRO_ID
      rf_vid = PRGVAR_NUM + iv
      call this%restart_file%Def_var( this%container%AUX_VARS(iv),   &
        auxvar_info(iv)%DESC, rf_vid, DIMTYPE_XYZ                    )
    end do
    do iv=1, QA
      rf_vid = rf_vid + 1
      call this%restart_file%Def_var( this%container%QTRC_VARS(iv), &
        TRACER_DESC(iv), rf_vid, DIMTYPE_XYZ                        )    
    end do
    return
  end subroutine AtmosVar_Write_restart_file_prep

!> Write data with atmospheric variables to restart file
!!
!OCL SERIAL
  subroutine AtmosVar_Write_restart_file( this )
    implicit none
    class(AtmosVars), intent(inout) :: this

    integer :: iv, rf_vid 
    !---------------------------------------

    call this%restart_file%End_def()

    !- Write restart file
    do iv=1, PRGVAR_NUM
      rf_vid = iv
      call this%restart_file%Write_var(rf_vid, this%container%PROG_VARS(iv) )
    end do
    do iv=1, AUXVAR_DENSHYDRO_ID
      rf_vid = PRGVAR_NUM + iv
      call this%restart_file%Write_var(rf_vid, this%container%AUX_VARS(iv) )
    end do
    do iv=1, QA
      rf_vid = rf_vid + 1
      call this%restart_file%Write_var(rf_vid, this%container%QTRC_VARS(iv) )
    end do

    !- Close restart file
    LOG_INFO("ATMOSVar_Write_restart_file",*) 'Close restart file (ATMOS) '
    call this%restart_file%Close()

    return
  end subroutine AtmosVar_Write_restart_file

!> Check the range of values with atmospheric variables
!!
!OCL SERIAL
  subroutine AtmosVars_Check( this, force )

    use scale_meshfield_statistics, only: &
      MeshField_statistics_total,         &
      MeshField_statistics_detail
    
    implicit none
    class(AtmosVars), intent(inout) :: this
    logical, intent(in), optional :: force

    integer :: iv
    integer :: iv_diag
    integer :: n
    logical  :: check

    class(MeshBase3D), pointer :: mesh3D
    class(LocalMeshBase), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: lcfield
    type(ElementBase), pointer :: elem
    character(len=H_MID) :: varname

    type(MeshField3D), pointer :: var    
    type(MeshField3D) :: vel_fields(3)
    type(MeshField3D) :: work

    integer :: ldomID
    !--------------------------------------------------------------------------

    if ( present(force) ) then
      check = force
    else
      check = this%check_range
    end if

    if (check) then
      mesh3D => this%container%PROG_VARS(1)%mesh

      do iv=1, PRGVAR_NUM
        if ( iv == PRGVAR_THERM_ID ) cycle
        
        var => this%container%PROG_VARS(iv)        
        do n=1, mesh3D%LOCAL_MESH_NUM
          lcmesh => mesh3D%lcmesh_list(n)
          elem => lcmesh%refElem

          call var%GetLocalMeshField(n, lcfield)

          write(varname,'(a,i3.3,a)') var%varname//'(domID=', n, ')' 

          ! Note : *acc update host* is called in VALCHECK of SCALE library. 
          call VALCHECK( elem%Np, 1, elem%Np, lcmesh%NeA, lcmesh%NeS, lcmesh%NeE, lcfield%val, &
            PROGVARS_check_min(iv), PROGVARS_check_max(iv), trim(varname), __FILE__, __LINE__  )
        end do
      end do

      do iv=1, 3
        iv_diag = ATMOS_DIAGVARS_U_ID + iv - 1
        call vel_fields(iv)%Init( ATMOS_DIAGVARS3D_VINFO(iv_diag)%NAME, "",  mesh3D )
        call AtmosVars_CalcDiagvar( this, vel_fields(iv)%varname, vel_fields(iv) )
#ifdef _OPENACC        
        do ldomID=1, mesh3D%LOCAL_MESH_NUM
          !$acc update host(vel_fields(iv)%local(ldomID)%val)
        end do
#endif        
      end do
      call MeshField_statistics_detail( vel_fields )
      do iv=1, 3
        call vel_fields(iv)%Final()
      end do

      ! Check total
      call work%Init("tmp", "", mesh3D)
      call work%Final()
    end if

    return
  end subroutine AtmosVars_Check

!> Put the stastics with atmospheric variables
!!
!OCL SERIAL
  subroutine AtmosVars_Monitor( this )
    implicit none
    class(AtmosVars), intent(inout) :: this
    !--------------------------------------------------------------------------
    call AtmosVars_Monitor_core( this%container )
    return
  end subroutine AtmosVars_Monitor
!OCL SERIAL
  subroutine AtmosVars_Monitor_core( this )  
    use scale_file_monitor_meshfield, only: &
      FILE_monitor_meshfield_put
    implicit none
    class(AtmosVarsContainer), intent(inout) :: this

    integer :: iv
    class(MeshBase3D), pointer :: mesh3D
    type(MeshField3D) :: work

    integer :: n
    integer :: ke
    class(LocalMesh3D), pointer :: lcmesh
    !--------------------------------------------------------------------------

    mesh3D => this%PROG_VARS(1)%mesh
    call work%Init("tmp", "", mesh3D)

    do iv=1, PRGVAR_NUM
      call FILE_monitor_meshfield_put( this%PROG_VARS(iv)%monitor_id, this%PROG_VARS(iv) )
    end do
  
    do iv=1, QA
      if ( this%QTRC_VARS(iv)%monitor_id > 0 ) then
        do n=1, mesh3D%LOCAL_MESH_NUM
          lcmesh => mesh3D%lcmesh_list(n)
          !$omp parallel do
          do ke=lcmesh%NeS, lcmesh%NeE
            work%local(n)%val(:,ke) = ( this%AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val(:,ke) &
                                      + this%PROG_VARS(PRGVAR_DDENS_ID  )%local(n)%val(:,ke) &
                                      ) * this%QTRC_VARS(iv)%local(n)%val(:,ke)
          end do
        end do
        call FILE_monitor_meshfield_put( this%QTRC_VARS(iv)%monitor_id, work )
      end if
    end do
    if ( DV_MONIT_id(IM_QTOT) > 0 ) then
      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh => mesh3D%lcmesh_list(n)
        !$omp parallel do
        do ke=lcmesh%NeS, lcmesh%NeE
          work%local(n)%val(:,ke) = ( this%AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val(:,ke) &
                                    + this%PROG_VARS(PRGVAR_DDENS_ID   )%local(n)%val(:,ke) &
                                    ) * ( 1.0_RP - this%AUX_VARS(AUXVAR_Qdry_ID)%local(n)%val(:,ke) )
        end do
      end do
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_QTOT), work )
    end if

    !##### Energy Budget #####

    if ( DV_MONIT_id(IM_ENGT) > 0 ) then
      call this%Calc_diagVar( 'ENGT', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGT), work )
    end if
    if ( DV_MONIT_id(IM_ENGP) > 0 ) then
      call this%Calc_diagVar( 'ENGP', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGP), work )
    end if
    if ( DV_MONIT_id(IM_ENGK) > 0 ) then
      call this%Calc_diagVar( 'ENGK', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGK), work )
    end if
    if ( DV_MONIT_id(IM_ENGI) > 0 ) then
      call this%Calc_diagVar( 'ENGI', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGI), work )
    end if

    call work%Final()

    return
  end subroutine AtmosVars_Monitor_core

  !> Preprocess operation for physical processes
  !!
!OCL SERIAL
  subroutine AtmosVars_PreprocOperationForPhys( this, dyncore )
    use scale_atm_dyn_dgm_driver_nonhydro3d, only: AtmDynDGMDriver_nonhydro3d
    implicit none
    class(AtmosVars), intent(inout) :: this
    class(AtmDynDGMDriver_nonhydro3d), intent(in) :: dyncore

    integer :: ic
    !--------------------------------------------------------------------------

    do ic=ATM_VARS_CONTAINER_PRIMARY_ID+1, ATM_VARS_CONTAINER_LIST_MAX
      if ( this%init_containers_item_flag(ic) ) then
        call this%container_list(ic)%Preproc_operation_for_phys( this%container, dyncore )
      end if
    end do
    return
  end subroutine AtmosVars_PreprocOperationForPhys

  !-----------------------------------------------------------------------------
  !> Calculate diagnostic variables
!OCL SERIAL  
  subroutine AtmosVars_CalculateDiagnostics( this )
    implicit none
    class(AtmosVars), intent(inout), target :: this
    !-------------------------------------------------------
    call this%container%Calc_diagnostics()
    return
  end subroutine AtmosVars_CalculateDiagnostics

!OCL SERIAL
  subroutine AtmosVars_CalcDiagvar( this, field_name, field_work ) 
    implicit none
    class(AtmosVars), intent(inout) :: this
    character(*), intent(in) :: field_name
    type(MeshField3D), intent(inout) :: field_work
    !--------------------------------------------------
    call this%container%Calc_diagVar( field_name, field_work )
    return
  end subroutine AtmosVars_CalcDiagvar

!OCL SERIAL
  subroutine AtmosVars_CalcDiagvar2D( this, field_name, field_work ) 
    implicit none
    class(AtmosVars), intent(inout) :: this
    character(*), intent(in) :: field_name
    type(MeshField2D), intent(inout), target :: field_work

    class(LocalMesh2D), pointer :: lcmesh2D

    integer :: n
    !--------------------------------------------------

    field_work%varname = field_name

    do n=1, field_work%mesh%LOCAL_MESH_NUM
      lcmesh2D => field_work%mesh%lcmesh_list(n)
      call vars_calc_diagnoseVar2D_lc( field_name, field_work%local(n)%val,  &
        this%ptr_MP_AUXVARS2D_manager,                                       &
        field_work%mesh, lcmesh2D, lcmesh2D%refElem2D                        )
    end do
    !$acc wait(1)
    return
  end subroutine AtmosVars_CalcDiagvar2D

!--- private -----

!OCL SERIAL
  subroutine vars_calc_diagnoseVar2D_lc( field_name, & ! (in)
    var_out,                                         & ! (out)
    MP_auxvars2D, mesh2D, lcmesh, elem )               ! (in)

    use mod_atmos_phy_mp_vars, only: &
      AtmosPhyMpVars_GetLocalMeshFields_sfcflx

    implicit none
    class(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem
    character(*), intent(in) :: field_name
    real(RP), intent(out) :: var_out(elem%Np,lcmesh%NeA)
    class(ModelVarManager), intent(inout) :: MP_auxvars2D
    class(MeshBase2D), intent(in) :: mesh2D

    integer :: ke, p

    class(LocalMeshFieldBase), pointer :: SFLX_rain_MP, SFLX_snow_MP, SFLX_ENGI_MP
    !-------------------------------------------------------------------------

    select case(trim(field_name))
    case('RAIN', 'SNOW')
      call AtmosPhyMpVars_GetLocalMeshFields_sfcflx( &
        lcmesh%lcdomID, mesh2D, MP_auxvars2D,        &
        SFLX_rain_MP, SFLX_snow_MP, SFLX_ENGI_MP     )
    end select

    select case(trim(field_name))
    case('RAIN')
      !$omp parallel do
      !$acc parallel loop collapse(2) present(SFLX_rain_MP%val, var_out) async(1)
      do ke=lcmesh%NeS, lcmesh%NeE
      do p=1, elem%Np
        var_out(p,ke) = SFLX_rain_MP%val(p,ke)
      end do
      end do
    case('SNOW')
      !$omp parallel do
      !$acc parallel loop collapse(2) present(var_out, SFLX_snow_MP) async(1)
      do ke=lcmesh%NeS, lcmesh%NeE
      do p=1, elem%Np
        var_out(p,ke) = SFLX_snow_MP%val(p,ke)
      end do
      end do
    case default
      LOG_ERROR("AtmosVars_calc_diagnoseVar2D_lc",*) 'The name of diagnostic variable is not suported. Check!', field_name
      call PRC_abort
    end select

    return
  end subroutine vars_calc_diagnoseVar2D_lc  
end module mod_atmos_vars