!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Variables
!!
!! @par Description
!!          Container for atmospheric variables
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
  use scale_tracer, only: &
    QA, TRACER_NAME, TRACER_DESC, TRACER_UNIT
  use scale_atmos_hydrometeor, only: &
    ATMOS_HYDROMETEOR_dry

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
    AUXVAR_DENSHYDRO_ID, AUXVAR_PRESHYDRO_ID, AUXVAR_PRESHYDRO_REF_ID,                          &
    AUXVAR_Rtot_ID, AUXVAR_CPtot_ID, AUXVAR_CVtot_ID,                                           &
    AUXVAR_PRES_ID, AUXVAR_PT_ID, AUXVAR_Qdry_ID,                                               &
    PHYTEND_DENS_ID, PHYTEND_MOMX_ID, PHYTEND_MOMY_ID, PHYTEND_MOMZ_ID, PHYTEND_RHOT_ID,        &
    PHYTEND_RHOH_ID    

  use mod_atmos_mesh, only: AtmosMesh
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: AtmosVars
    class(AtmosMesh), pointer :: mesh

    !- prognostic variables
    type(MeshField3D), allocatable :: PROG_VARS(:)
    type(ModelVarManager) :: PROGVARS_manager
    integer :: PROG_VARS_commID

    !- tracer variables
    type(MeshField3D), allocatable :: QTRC_VARS(:)
    type(ModelVarManager) :: QTRCVARS_manager
    integer :: QTRC_VARS_commID 

    !- auxiliary variables    
    type(MeshField3D), allocatable :: AUX_VARS(:)
    type(ModelVarManager) :: AUXVARS_manager 
    integer :: AUX_VARS_commID

    !- auxiliary variables (2D)
    type(MeshField2D), allocatable :: AUX_VARS2D(:)
    type(ModelVarManager) :: AUXVARS2D_manager 
    
    !-
    type(ModelVarManager), pointer :: ptr_MP_AUXVARS2D_manager

    logical :: moist
    type(MeshField3D), pointer :: QV
    type(MeshField3D) :: zero

    !- Tendency with physics
    type(MeshField3D), allocatable :: PHY_TEND(:)
    type(ModelVarManager) :: PHYTENDS_manager 
    integer :: PHYTENDS_commID
    integer :: PHYTEND_NUM_TOT

    !--
    integer, allocatable :: DIAGVARS2D_HISTID(:)
    integer, allocatable :: DIAGVARS3D_HISTID(:)

    type(FILE_restart_meshfield_component) :: restart_file
    
    logical :: check_range
    logical :: check_total

  contains
    procedure :: Init => AtmosVars_Init
    procedure :: Final => AtmosVars_Final
    procedure :: Calc_diagnostics => AtmosVars_CalculateDiagnostics
    procedure :: Calc_diagVar => AtmosVars_CalcDiagvar
    procedure :: Calc_diagVar2D => AtmosVars_CalcDiagvar2D
    procedure :: History => AtmosVars_History
    procedure :: Check   => AtmosVars_Check
    procedure :: Monitor => AtmosVars_Monitor
    procedure :: Read_restart_file => AtmosVar_Read_restart_file
    procedure :: Write_restart_file => AtmosVar_Write_restart_file
    procedure :: Regist_physvar_manager => AtmosVars_Regist_physvar_manager
  end type AtmosVars

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
              

  ! Surface variables
  
  integer, public, parameter :: ATMOS_AUXVARS2D_PREC_ID      = 1
  integer, public, parameter :: ATMOS_AUXVARS2D_PREC_ENGI_ID = 2
  integer, public, parameter :: ATMOS_AUXVARS2D_NUM          = 2

  type(VariableInfo), public :: ATMOS_AUXVARS2D_VINFO(ATMOS_AUXVARS2D_NUM)
  DATA ATMOS_AUXVARS2D_VINFO / &
    VariableInfo( ATMOS_AUXVARS2D_PREC_ID     ,      'PREC', 'surface precipitaion flux'        , 'kg/m2/s', 2, 'XY', 'precipitation_flux'  ), &
    VariableInfo( ATMOS_AUXVARS2D_PREC_ENGI_ID, 'PREC_ENGI', 'internal energy of precipitation' ,    'J/m2', 2, 'XY', ''  )                    /
  
                    

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
  !-------------------

  private :: vars_calc_diagnoseVar_lc
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

!OCL SERIAL
  subroutine AtmosVars_Init( this, atm_mesh )
    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_dry
    use scale_file_monitor_meshfield, only:    &
      MONITOR_reg => FILE_monitor_meshfield_reg

    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRGVAR_SCALAR_NUM, PRGVAR_HVEC_NUM,          &
      atm_dyn_dgm_nonhydro3d_common_setup_variables
    implicit none

    class(AtmosVars), target, intent(inout) :: this
    class(AtmosMesh), target, intent(inout) :: atm_mesh

    integer :: n
    integer :: iv
    integer :: iq
    logical :: reg_file_hist

    type(ModelVarManager) :: diagvar_manager               ! dummy
    type(MeshField2D) :: diag_vars2D(ATMOS_DIAGVARS2D_NUM) ! dummy
    type(MeshField3D) :: diag_vars3D(ATMOS_DIAGVARS3D_NUM) ! dummy

    logical :: CHECK_RANGE    = .false.
    logical :: CHECK_TOTAL    = .false.

    namelist / PARAM_ATMOS_VARS / &
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

    type(VariableInfo) :: prgvar_info(PRGVAR_NUM)
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

    !- Set the pointer of mesh
    this%mesh => atm_mesh
    mesh3D => atm_mesh%ptr_mesh
    call mesh3D%GetMesh2D( mesh2D )

    !- Initialize variables associated with dynamical core 
    !  (prognostic variables, tracer variables, 3D auxiliary variables, and tendencies of physical processes)

    call this%PROGVARS_manager%Init()
    call this%QTRCVARS_manager%Init()
    call this%AUXVARS_manager%Init()
    call this%PHYTENDS_manager%Init()

    allocate( this%PROG_VARS(PRGVAR_NUM) )
    allocate( this%QTRC_VARS(0:QA) )
    allocate( this%AUX_VARS(AUXVAR_NUM) )

    this%PHYTEND_NUM_TOT = PHYTEND_NUM1 + max(1,QA)
    allocate( this%PHY_TEND(this%PHYTEND_NUM_TOT) )

    call atm_dyn_dgm_nonhydro3d_common_setup_variables( &
      this%PROG_VARS, this%QTRC_VARS, this%AUX_VARS, this%PHY_TEND,                              & ! (inout)
      this%PROGVARS_manager, this%QTRCVARS_manager, this%AUXVARS_manager, this%PHYTENDS_manager, & ! (inout)
      this%PHYTEND_NUM_TOT, mesh3D,                                                              & ! (in)
      prgvar_info ) ! (out)
 
    ! Setup communicator
    
    call atm_mesh%Create_communicator( &
      PRGVAR_SCALAR_NUM, PRGVAR_HVEC_NUM, 0,              & ! (in)
      this%PROGVARS_manager,                              & ! (inout)
      this%PROG_VARS(:),                                  & ! (in)
      this%PROG_VARS_commID                               ) ! (out)
    
    if ( QA > 0 ) then
      call atm_mesh%Create_communicator( &
        QA, 0, 0,                        & ! (in)
        this%QTRCVARS_manager,           & ! (inout)
        this%QTRC_VARS(1:QA),            & ! (in)
        this%QTRC_VARS_commID            ) ! (out)
    end if

    call atm_mesh%Create_communicator( &
      AUXVAR_NUM, 0, 0,                & ! (in)
      this%AUXVARS_manager,            & ! (inout)
      this%AUX_VARS(:),                & ! (in)
      this%AUX_VARS_commID             ) ! (out)

    ! Output list of prognostic variables

    LOG_NEWLINE
    LOG_INFO("ATMOS_vars_setup",*) 'List of prognostic variables (ATMOS) '
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, PRGVAR_NUM
      LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
      'NO.',iv,'|',prgvar_info(iv)%NAME,'|', prgvar_info(iv)%DESC,'[', prgvar_info(iv)%UNIT,']'
    end do
    do iv = 1, QA
      LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
      'NO.',PRGVAR_NUM+iv,'|',TRACER_NAME(iv),'|', TRACER_DESC(iv),'[', TRACER_UNIT(iv),']'
    end do
    LOG_NEWLINE


    !- Initialize 2D auxiliary variables
    call this%AUXVARS2D_manager%Init()
    allocate( this%AUX_VARS2D(ATMOS_AUXVARS2D_NUM) )
    
    reg_file_hist = .true.
    do iv = 1, ATMOS_AUXVARS2D_NUM
      call this%AUXVARS2D_manager%Regist(    &
        ATMOS_AUXVARS2D_VINFO(iv), mesh2D,   & ! (in) 
        this%AUX_VARS2D(iv),                 & ! (inout)
        reg_file_hist, fill_zero=.true.      ) ! (in)
    end do
  

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
        PRGVAR_NUM + AUXVAR_NUM                                     )
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

!OCL SERIAL
  subroutine AtmosVars_Final( this )
    implicit none
    class(AtmosVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosVars_Final',*)

    call this%restart_file%Final()

    call this%PROGVARS_manager%Final()
    deallocate( this%PROG_VARS )

    call this%QTRCVARS_manager%Final()
    deallocate( this%QTRC_VARS )

    call this%AUXVARS_manager%Final()
    deallocate( this%AUX_VARS )

    call this%AUXVARS2D_manager%Final()
    deallocate( this%AUX_VARS2D )

    call this%PHYTENDS_manager%Final()
    deallocate( this%PHY_TEND )

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

!OCL SERIAL
  subroutine AtmosVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosVars), intent(inout) :: this
  
    integer :: v
    integer :: hst_id

    type(MeshField3D) :: tmp_field3D
    class(MeshBase3D), pointer :: mesh3D

    type(MeshField2D) :: tmp_field2D
    class(MeshBase2D), pointer :: mesh2D
    !-------------------------------------------------------------------------

    mesh3D => this%PROG_VARS(1)%mesh
    call mesh3D%GetMesh2D(mesh2D)

    do v = 1, PRGVAR_NUM
      hst_id = this%PROG_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%PROG_VARS(v) )
    end do

    do v = 1, QA
      hst_id = this%QTRC_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%QTRC_VARS(v) )
    end do

    call this%Calc_diagnostics()
    
    do v = 1, AUXVAR_NUM
      hst_id = this%AUX_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%AUX_VARS(v) )
    end do
    do v = 1, ATMOS_AUXVARS2D_NUM
      hst_id = this%AUX_VARS2D(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%AUX_VARS2D(v) )
    end do

    do v = 1, this%PHYTEND_NUM_TOT
      hst_id = this%PHY_TEND(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%PHY_TEND(v) )
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
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%PROG_VARS(iv)%varname, &
        this%PROG_VARS(iv)                                                      )
    end do
    do iv=1, AUXVAR_DENSHYDRO_ID
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%AUX_VARS(iv)%varname, &
        this%AUX_VARS(iv)                                                      )
    end do
    do iv=1, QA
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%QTRC_VARS(iv)%varname, &
        this%QTRC_VARS(iv)                                                      )
    end do

    !- Close restart file
    LOG_INFO("ATMOSVar_read_restart_file",*) 'Close restart file (ATMOS) '
    call this%restart_file%Close()

    !-- Diagnostic pressure
    call vars_calc_specific_heat( this )
    
    call dyncore%calc_pressure( this%AUX_VARS(AUXVAR_PRES_ID), &
      this%PROGVARS_manager, this%AUXVARS_manager              )

    ! Set reference value of hydrostatic pressure
    Phyd_ref => this%AUX_VARS(AUXVAR_PRESHYDRO_REF_ID)
    mesh3D => Phyd_ref%mesh
    do domid=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(domid)
      !$omp parallel do
      do ke=lcmesh3D%NeS, lcmesh3D%NeE
        Phyd_ref%local(domid)%val(:,ke) = 0.0_RP
      end do
    end do

    !-- Check read data
    call this%Check( force = .true. )

    !-- Calculate diagnostic variables
    call this%Calc_diagnostics()   

    !-- Communicate halo data of hydrostatic & diagnostic variables
    call this%AUXVARS_manager%MeshFieldComm_Exchange()

    !-- Set horizontal gradient of hydrostatic pressure
    call dyncore%update_phyd_hgrad( this%AUX_VARS(AUXVAR_PRESHYDRO_ID), Phyd_ref, &
      mesh3D, atmos_mesh%element3D_operation )

    return
  end subroutine AtmosVar_Read_restart_file

!OCL SERIAL
  subroutine AtmosVar_write_restart_file( this )
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      atm_dyn_dgm_nonhydro3d_common_get_varinfo
    
    implicit none
    class(AtmosVars), intent(inout) :: this

    type(VariableInfo) :: prgvar_info(PRGVAR_NUM)
    type(VariableInfo) :: auxvar_info(AUXVAR_NUM)

    integer :: iv, rf_vid 
    !---------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("ATMOSVar_write_restart_file",*) 'Create restart file (ATMOS) '

    !- Check data which will be written to restart file
    call this%Check( force = .true. )

    !- Create restart file
    call this%restart_file%Create()
    call PRC_mpibarrier()

    !- Define variables

    call atm_dyn_dgm_nonhydro3d_common_get_varinfo( prgvar_info, auxvar_info )

    do iv=1, PRGVAR_NUM
      rf_vid = iv
      call this%restart_file%Def_var( this%PROG_VARS(iv),  &
        prgvar_info(iv)%DESC, rf_vid, DIMTYPE_XYZ          )
    end do
    do iv=1, AUXVAR_DENSHYDRO_ID
      rf_vid = PRGVAR_NUM + iv
      call this%restart_file%Def_var( this%AUX_VARS(iv),   &
        auxvar_info(iv)%DESC, rf_vid, DIMTYPE_XYZ          )
    end do
    do iv=1, QA
      rf_vid = rf_vid + 1
      call this%restart_file%Def_var( this%QTRC_VARS(iv), &
        TRACER_DESC(iv), rf_vid, DIMTYPE_XYZ              )    
    end do

    call this%restart_file%End_def()

    !- Write restart file
    do iv=1, PRGVAR_NUM
      rf_vid = iv
      call this%restart_file%Write_var(rf_vid, this%PROG_VARS(iv) )
    end do
    do iv=1, AUXVAR_DENSHYDRO_ID
      rf_vid = PRGVAR_NUM + iv
      call this%restart_file%Write_var(rf_vid, this%AUX_VARS(iv) )
    end do
    do iv=1, QA
      rf_vid = rf_vid + 1
      call this%restart_file%Write_var(rf_vid, this%QTRC_VARS(iv) )
    end do

    !- Close restart file
    LOG_INFO("ATMOSVar_write_restart_file",*) 'Close restart file (ATMOS) '
    call this%restart_file%Close()

    return
  end subroutine AtmosVar_write_restart_file

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

    type(MeshField3D) :: vel_fields(3)
    type(MeshField3D) :: work
    !--------------------------------------------------------------------------

    if ( present(force) ) then
      check = force
    else
      check = this%check_range
    end if

    if (check) then
      do iv=1, PRGVAR_NUM
        if ( iv == PRGVAR_THERM_ID ) cycle
        
        mesh3D => this%PROG_VARS(iv)%mesh
        do n=1, mesh3D%LOCAL_MESH_NUM
          lcmesh => mesh3D%lcmesh_list(n)
          elem => lcmesh%refElem
          call this%PROG_VARS(iv)%GetLocalMeshField(n, lcfield)
          write(varname,'(a,i3.3,a)') this%PROG_VARS(iv)%varname//'(domID=', n, ')' 
          call VALCHECK( elem%Np, 1, elem%Np, lcmesh%NeA, lcmesh%NeS, lcmesh%NeE, lcfield%val(:,:), &
            PROGVARS_check_min(iv), PROGVARS_check_max(iv), trim(varname), __FILE__, __LINE__       )
        end do
      end do

      mesh3D => this%PROG_VARS(1)%mesh
      do iv=1, 3
        iv_diag = ATMOS_DIAGVARS_U_ID + iv - 1
        call vel_fields(iv)%Init( ATMOS_DIAGVARS3D_VINFO(iv_diag)%NAME, "",  mesh3D )
        call AtmosVars_CalcDiagvar( this, vel_fields(iv)%varname, vel_fields(iv) )
      end do
      call MeshField_statistics_detail( vel_fields(:) )
      do iv=1, 3
        call vel_fields(iv)%Final()
      end do
    end if

    if ( present(force) ) then
      check = force
    else
      check = this%check_total
    end if
    if (check) then
      mesh3D => this%PROG_VARS(1)%mesh
      call work%Init("tmp", "", mesh3D)
      call work%Final()
    end if

    return
  end subroutine AtmosVars_Check

!OCL SERIAL
  subroutine AtmosVars_Monitor( this )
    use scale_file_monitor_meshfield, only: &
      FILE_monitor_meshfield_put

    implicit none
    class(AtmosVars), intent(inout) :: this

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
      call AtmosVars_CalcDiagvar( this, 'ENGT', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGT), work )
    end if
    if ( DV_MONIT_id(IM_ENGP) > 0 ) then
      call AtmosVars_CalcDiagvar( this, 'ENGP', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGP), work )
    end if
    if ( DV_MONIT_id(IM_ENGK) > 0 ) then
      call AtmosVars_CalcDiagvar( this, 'ENGK', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGK), work )
    end if
    if ( DV_MONIT_id(IM_ENGI) > 0 ) then
      call AtmosVars_CalcDiagvar( this, 'ENGI', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGI), work )
    end if

    call work%Final()

    return
  end subroutine AtmosVars_Monitor

  !----  Getter ---------------------------------------------------------------------------

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshPrgVar( domID, mesh, prgvars_list, auxvars_list, &
     varid,                                                                         &
     var, DENS_hyd, PRES_hyd, lcmesh3D                                              )
   
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    integer, intent(in) :: varid
    class(LocalMeshFieldBase), pointer, intent(out) :: var
    class(LocalMeshFieldBase), pointer, intent(out), optional :: DENS_hyd, PRES_hyd
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call prgvars_list%Get(varid, field)
    call field%GetLocalMeshField(domID, var)

    if (present(DENS_hyd)) then
      call auxvars_list%Get(AUXVAR_DENSHYDRO_ID, field)
      call field%GetLocalMeshField(domID, DENS_hyd)
    end if
    if (present(PRES_hyd)) then
      call auxvars_list%Get(AUXVAR_PRESHYDRO_ID, field)
      call field%GetLocalMeshField(domID, PRES_hyd)
    end if

    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshPrgVar

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshPrgVars( domID, mesh, prgvars_list, auxvars_list, &
    DDENS, MOMX, MOMY, MOMZ, THERM,                                                  &
    DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot,                                          &
    lcmesh3D                                                                         )
    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: DDENS, MOMX, MOMY, MOMZ, THERM
    class(LocalMeshFieldBase), pointer, intent(out) :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer, intent(out) :: Rtot, CVtot, CPtot
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call prgvars_list%Get(PRGVAR_DDENS_ID, field)
    call field%GetLocalMeshField(domID, DDENS)

    call prgvars_list%Get(PRGVAR_MOMX_ID, field)
    call field%GetLocalMeshField(domID, MOMX)
    
    call prgvars_list%Get(PRGVAR_MOMY_ID, field)
    call field%GetLocalMeshField(domID, MOMY)

    call prgvars_list%Get(PRGVAR_MOMZ_ID, field)
    call field%GetLocalMeshField(domID, MOMZ)

    call prgvars_list%Get(PRGVAR_THERM_ID, field)
    call field%GetLocalMeshField(domID, THERM)
  
    !--
    call auxvars_list%Get(AUXVAR_DENSHYDRO_ID, field)
    call field%GetLocalMeshField(domID, DENS_hyd)

    call auxvars_list%Get(AUXVAR_PRESHYDRO_ID, field)
    call field%GetLocalMeshField(domID, PRES_hyd)

    call auxvars_list%Get(AUXVAR_Rtot_ID, field)
    call field%GetLocalMeshField(domID, Rtot)

    call auxvars_list%Get(AUXVAR_CVtot_ID, field)
    call field%GetLocalMeshField(domID, CVtot)

    call auxvars_list%Get(AUXVAR_CPtot_ID, field)
    call field%GetLocalMeshField(domID, CPtot)

    !---
    
    if ( present(lcmesh3D) ) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshPrgVars

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshSfcVar( domID, mesh, auxvars2D_list, &
    PREC, PREC_ENGI, lcmesh2D                                           )

   implicit none
   integer, intent(in) :: domID
   class(MeshBase), intent(in) :: mesh
   class(ModelVarManager), intent(inout) :: auxvars2D_list
   class(LocalMeshFieldBase), pointer, intent(out) :: PREC, PREC_ENGI
   class(LocalMesh2D), pointer, intent(out), optional :: lcmesh2D

   class(MeshFieldBase), pointer :: field
   class(LocalMeshBase), pointer :: lcmesh
   !-------------------------------------------------------

   !--
   call auxvars2D_list%Get(ATMOS_AUXVARS2D_PREC_ID, field)
   call field%GetLocalMeshField(domID, PREC)

   call auxvars2D_list%Get(ATMOS_AUXVARS2D_PREC_ENGI_ID, field)
   call field%GetLocalMeshField(domID, PREC_ENGI)

   if (present(lcmesh2D)) then
     call mesh%GetLocalMesh( domID, lcmesh )
     nullify( lcmesh2D )

     select type(lcmesh)
     type is (LocalMesh2D)
       if (present(lcmesh2D)) lcmesh2D => lcmesh
     end select
   end if

   return
 end subroutine AtmosVars_GetLocalMeshSfcVar

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshQTRCVar( domID, mesh, trcvars_list,  &
    varid,                                                              &
    var, lcmesh3D                                                       )

   implicit none
   integer, intent(in) :: domID
   class(MeshBase), intent(in) :: mesh
   class(ModelVarManager), intent(inout) :: trcvars_list
   integer, intent(in) :: varid
   class(LocalMeshFieldBase), pointer, intent(out) :: var
   class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

   class(MeshFieldBase), pointer :: field
   class(LocalMeshBase), pointer :: lcmesh
   !-------------------------------------------------------

   !--
   call trcvars_list%Get(varid, field)
   call field%GetLocalMeshField(domID, var)

   if (present(lcmesh3D)) then
     call mesh%GetLocalMesh( domID, lcmesh )
     nullify( lcmesh3D )

     select type(lcmesh)
     type is (LocalMesh3D)
       if (present(lcmesh3D)) lcmesh3D => lcmesh
     end select
   end if

   return
  end subroutine AtmosVars_GetLocalMeshQTRCVar

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshQTRC_Qv( domID, mesh, trcvars_list, forcing_list, &
    var, var_tp, lcmesh3D                                               )

    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_dry, &
      I_QV    
   implicit none
   integer, intent(in) :: domID
   class(MeshBase), intent(in) :: mesh
   class(ModelVarManager), intent(inout) :: trcvars_list
   class(ModelVarManager), intent(inout) :: forcing_list
   class(LocalMeshFieldBase), pointer, intent(out) :: var
   class(LocalMeshFieldBase), pointer, intent(out) :: var_tp
   class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

   class(MeshFieldBase), pointer :: field
   class(LocalMeshBase), pointer :: lcmesh

   integer :: iq, tend_iq
   !-------------------------------------------------------

   !--
   if ( ATMOS_HYDROMETEOR_dry ) then
     iq =0; tend_iq = PHYTEND_NUM1+1
   else
     iq = I_QV; tend_iq = PHYTEND_NUM1 + I_QV
   end if

   call trcvars_list%Get(iq, field)
   call field%GetLocalMeshField(domID, var)

   call forcing_list%Get(tend_iq, field)
   call field%GetLocalMeshField(domID, var_tp)

   if (present(lcmesh3D)) then
     call mesh%GetLocalMesh( domID, lcmesh )
     nullify( lcmesh3D )

     select type(lcmesh)
     type is (LocalMesh3D)
       if (present(lcmesh3D)) lcmesh3D => lcmesh
     end select
   end if

   return
  end subroutine AtmosVars_GetLocalMeshQTRC_Qv

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshQTRCVarList( domID, mesh, trcvars_list,  &
    varid_s,                                                                &
    var_list, lcmesh3D                                                      )

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: trcvars_list
    integer, intent(in) :: varid_s
    type(LocalMeshFieldBaseList), intent(out) :: var_list(:)
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh

    integer :: iq
    !-------------------------------------------------------

    !--
    do iq = varid_s, varid_s + size(var_list) - 1
      call trcvars_list%Get(iq, field)
      call field%GetLocalMeshField(domID, var_list(iq-varid_s+1)%ptr)
    end do
    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshQTRCVarList

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshPhyAuxVars( domID, mesh, phyauxvars_list, &
    PRES, PT,                                                                &
    lcmesh3D                                                                 )

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: phyauxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: PRES, PT
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--    
    call phyauxvars_list%Get(AUXVAR_PRES_ID, field)
    call field%GetLocalMeshField(domID, PRES)

    call phyauxvars_list%Get(AUXVAR_PT_ID, field)
    call field%GetLocalMeshField(domID, PT)

    !---
    
    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshPhyAuxVars

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshPhyTends( domID, mesh, phytends_list,  &
    DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p,                  &
    RHOQ_tp,                                                              &
    lcmesh3D                                                              )

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: phytends_list
    class(LocalMeshFieldBase), pointer, intent(out) :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp
    class(LocalMeshFieldBase), pointer, intent(out) :: RHOH_p
    type(LocalMeshFieldBaseList), intent(inout), optional :: RHOQ_tp(QA)
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh

    integer :: iq
    !-------------------------------------------------------

    !--
    call phytends_list%Get(PHYTEND_DENS_ID, field)
    call field%GetLocalMeshField(domID, DENS_tp)

    call phytends_list%Get(PHYTEND_MOMX_ID, field)
    call field%GetLocalMeshField(domID, MOMX_tp)
    
    call phytends_list%Get(PHYTEND_MOMY_ID, field)
    call field%GetLocalMeshField(domID, MOMY_tp)

    call phytends_list%Get(PHYTEND_MOMZ_ID, field)
    call field%GetLocalMeshField(domID, MOMZ_tp)

    call phytends_list%Get(PHYTEND_RHOT_ID, field)
    call field%GetLocalMeshField(domID, RHOT_tp)

    call phytends_list%Get(PHYTEND_RHOH_ID, field)
    call field%GetLocalMeshField(domID, RHOH_p)

    if ( present(RHOQ_tp) ) then
      do iq = 1, QA
        call phytends_list%Get(PHYTEND_NUM1+iq, field)
        call field%GetLocalMeshField(domID, RHOQ_tp(iq)%ptr)  
      end do
    end if

    !---
    if ( present(lcmesh3D) ) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if ( present(lcmesh3D) ) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshPhyTends

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshQTRCPhyTend( domID, mesh, phytends_list,  &
    qtrcid,                                                                  &
    RHOQ_tp                                                                  )

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: phytends_list
    integer, intent(in) :: qtrcid
    class(LocalMeshFieldBase), pointer, intent(out) :: RHOQ_tp

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    call phytends_list%Get(PHYTEND_NUM1 + qtrcid, field)
    call field%GetLocalMeshField(domID, RHOQ_tp)

    return
  end subroutine AtmosVars_GetLocalMeshQTRCPhyTend  

  !-----------------------------------------------------------------------------
  !> Calculate diagnostic variables
!OCL SERIAL  
  subroutine AtmosVars_CalculateDiagnostics( this )
    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CPdry => CONST_CPdry,    &
      CVdry => CONST_CVdry,    &
      PRES00 => CONST_PRE00
    use scale_tracer, only: &
      TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP    
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat
    implicit none
    class(AtmosVars), intent(inout), target :: this

    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    integer :: varid
    integer :: ke

    class(MeshField3D), pointer :: field
    class(ElementBase3D), pointer :: elem3D

    type(LocalMeshFieldBaseList) :: QTRC(QA)
    !-------------------------------------------------------

    ! Calculate specific heat
    call vars_calc_specific_heat( this )
    
    ! Calculate diagnostic variables
    do varid=AUXVAR_DENSHYDRO_ID+1, AUXVAR_PT_ID
      field => this%AUX_VARS(varid)
      do n=1, field%mesh%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshQTRCVarList( n, &
          field%mesh, this%QTRCVARS_manager,       &
          1, QTRC, lcmesh3D )

        elem3D => lcmesh3D%refElem3D

        call vars_calc_diagnoseVar_lc( &
          field%varname, field%local(n)%val,                       &
          this%PROG_VARS(PRGVAR_DDENS_ID)%local(n)%val,            &
          this%PROG_VARS(PRGVAR_MOMX_ID)%local(n)%val,             &
          this%PROG_VARS(PRGVAR_MOMY_ID)%local(n)%val,             &
          this%PROG_VARS(PRGVAR_MOMZ_ID)%local(n)%val,             &
          this%AUX_VARS(AUXVAR_PRES_ID)%local(n)%val,              &
          this%AUX_VARS(AUXVAR_QDRY_ID)%local(n)%val,              &
          QTRC,                                                    &
          this%AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val,         & 
          this%AUX_VARS(AUXVAR_PRESHYDRO_ID)%local(n)%val,         &
          this%AUX_VARS(AUXVAR_Rtot_ID )%local(n)%val,             & 
          this%AUX_VARS(AUXVAR_CVtot_ID)%local(n)%val,             & 
          this%AUX_VARS(AUXVAR_CPtot_ID)%local(n)%val,             & 
          lcmesh3D, lcmesh3D%refElem3D )
      end do
    end do

    return
  end subroutine AtmosVars_CalculateDiagnostics

!OCL SERIAL
  subroutine AtmosVars_CalcDiagvar( this, field_name, field_work ) 
    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CPdry => CONST_CPdry,    &
      CVdry => CONST_CVdry,    &
      PRES00 => CONST_PRE00
    use scale_tracer, only: &
      TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP    
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat

    implicit none
    class(AtmosVars), intent(inout) :: this
    character(*), intent(in) :: field_name
    type(MeshField3D), intent(inout) :: field_work

    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D
    integer :: n
    integer :: ke
    integer :: iq

    type(MeshField3D) :: field_work_UVmet(2)
    logical :: is_UVmet
    integer :: UVmet_i

    type(LocalMeshFieldBaseList) :: QTRC(QA)
    !--------------------------------------------------

    is_UVmet = .false.
    if ( field_name == 'Umet' ) then
      is_UVmet = .true.; UVmet_i = 1
    else if ( field_name == 'Vmet' ) then      
      is_UVmet = .true.; UVmet_i = 2
    end if

    field_work%varname = field_name

    do n=1, field_work%mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshQTRCVarList( n, &
        field_work%mesh, this%QTRCVARS_manager,  &
        1, QTRC, lcmesh3D )
      
      if ( .not. is_UVmet ) then
        call vars_calc_diagnoseVar_lc( field_name, field_work%local(n)%val,  &
          this%PROG_VARS(PRGVAR_DDENS_ID)%local(n)%val,               &
          this%PROG_VARS(PRGVAR_MOMX_ID)%local(n)%val,                &
          this%PROG_VARS(PRGVAR_MOMY_ID)%local(n)%val,                &
          this%PROG_VARS(PRGVAR_MOMZ_ID)%local(n)%val,                &
          this%AUX_VARS(AUXVAR_PRES_ID)%local(n)%val,                 &
          this%AUX_VARS(AUXVAR_QDRY_ID)%local(n)%val,                 &
          QTRC,                                                       &
          this%AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val,            & 
          this%AUX_VARS(AUXVAR_PRESHYDRO_ID)%local(n)%val,            &
          this%AUX_VARS(AUXVAR_Rtot_ID )%local(n)%val,                & 
          this%AUX_VARS(AUXVAR_CVtot_ID)%local(n)%val,                & 
          this%AUX_VARS(AUXVAR_CPtot_ID)%local(n)%val,                & 
          lcmesh3D, lcmesh3D%refElem3D )
      else
        call field_work_UVmet(1)%Init( 'Umet', '', field_work%mesh )
        call field_work_UVmet(2)%Init( 'Vmet', '', field_work%mesh )
        call this%mesh%Calc_UVmet( &
          this%PROG_VARS(PRGVAR_MOMX_ID), this%PROG_VARS(PRGVAR_MOMY_ID), & ! (in)
          field_work_UVmet(1), field_work_UVmet(2)                        ) ! (inout)
        !$omp parallel do
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          field_work%local(n)%val(:,ke) = field_work_UVmet(UVmet_i)%local(n)%val(:,ke) &
            / ( this%AUX_VARS (AUXVAR_DENSHYDRO_ID)%local(n)%val(:,ke)                 &
              + this%PROG_VARS(PRGVAR_DDENS_ID   )%local(n)%val(:,ke)                  )
        end do
        call field_work_UVmet(1)%Final()
        call field_work_UVmet(2)%Final()
      end if
    end do

    return
  end subroutine AtmosVars_CalcDiagvar

!OCL SERIAL
  subroutine AtmosVars_CalcDiagvar2D( this, field_name, field_work ) 
    implicit none
    class(AtmosVars), intent(inout) :: this
    character(*), intent(in) :: field_name
    type(MeshField2D), intent(inout), target :: field_work

    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase2D), pointer :: elem2D

    integer :: n
    !--------------------------------------------------

    field_work%varname = field_name

    do n=1, field_work%mesh%LOCAL_MESH_NUM
      lcmesh2D => field_work%mesh%lcmesh_list(n)
      call vars_calc_diagnoseVar2D_lc( field_name, field_work%local(n)%val,  &
        this%ptr_MP_AUXVARS2D_manager,                                       &
        field_work%mesh, lcmesh2D, lcmesh2D%refElem2D                        )
    end do

    return
  end subroutine AtmosVars_CalcDiagvar2D

!-- private -----------------------------------------------------------------------
    
!OCL SERIAL
  subroutine vars_calc_diagnoseVar_lc( field_name, var_out,  &
    DDENS_, MOMX_, MOMY_, MOMZ_, PRES_, QDRY_, QTRC,         &
    DENS_hyd, PRES_hyd, Rtot, CVtot, CPTot,                  &
    lcmesh, elem )

    use scale_const, only: &
      GRAV => CONST_GRAV,   &
      Rdry => CONST_Rdry,   &
      Rvap => CONST_Rvap,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    use scale_tracer, only: &
      TRACER_inq_id,          &
      TRACER_CV, TRACER_ENGI0
    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_psat_liq
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    character(*), intent(in) :: field_name
    real(RP), intent(out) :: var_out(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: QDRY_(elem%Np,lcmesh%NeA)
    type(LocalMeshFieldBaseList), intent(in) :: QTRC(QA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot (elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem%Np,lcmesh%NeA)

    integer :: ke, ke2D
    integer :: iq
    real(RP) :: DENS(elem%Np), TEMP(elem%Np)
    real(RP) :: mom_u1(elem%Np), mom_u2(elem%Np), G_11(elem%Np), G_12(elem%Np), G_22(elem%Np)
    real(RP) :: PSAT(elem%Np)

    integer :: iq_QV
    !-------------------------------------------------------------------------

    select case(trim(field_name))
    case('DENS')
      !$omp parallel do
      do ke=1, lcmesh%Ne
        var_out(:,ke) = DDENS_(:,ke) + DENS_hyd(:,ke)
      end do
    
    case('U')
      !$omp parallel do private (DENS)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = MOMX_(:,ke) / DENS(:)
      end do
    
    case('V')
      !$omp parallel do private (DENS)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = MOMY_(:,ke) / DENS(:)
      end do  
          
    case('W')
      !$omp parallel do private (DENS)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = MOMZ_(:,ke) / DENS(:)
      end do
    
    case ( 'PRES' )
    case('PRES_diff')  
      !$omp parallel do
      do ke=1, lcmesh%Ne
        var_out(:,ke) = PRES_(:,ke) - PRES_hyd(:,ke)
      end do
    
    case('T')
      !$omp parallel do
      do ke=1, lcmesh%Ne
        var_out(:,ke) = PRES_(:,ke) / (Rtot(:,ke) * (DDENS_(:,ke) + DENS_hyd(:,ke)) )
      end do
    
    case('T_diff')
      !$omp parallel do
      do ke=1, lcmesh%Ne
        var_out(:,ke) = PRES_(:,ke) / ( Rtot(:,ke) * (DDENS_(:,ke) + DENS_hyd(:,ke)) ) &
                      - PRES_hyd(:,ke) / ( Rdry * DENS_hyd(:,ke) )
      end do
    
    case('PT')
      !$omp parallel do private( DENS )
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = PRES_(:,ke) / (Rtot(:,ke) * DENS(:) ) * ( PRES00 / PRES_(:,ke) )**( Rtot(:,ke) / CPtot(:,ke) )
      end do 
    
    case('PT_diff')
      !$omp parallel do private( DENS )
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = PRES_(:,ke) / (Rtot(:,ke) * DENS(:) ) * ( PRES00 / PRES_(:,ke) )**( Rtot(:,ke) / CPtot(:,ke) ) &
                      - PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) / DENS_hyd(:,ke)
      end do 
    
    case( 'RH', 'RHL' )
      if ( ATMOS_HYDROMETEOR_dry ) then
        var_out(:,ke) = 0.0_RP
      else
        call TRACER_inq_id( "QV", iq_QV )

        !$omp parallel do private (TEMP, PSAT)
        do ke=1, lcmesh%Ne
          TEMP(:) = PRES_(:,ke) / (Rtot(:,ke) * (DDENS_(:,ke) + DENS_hyd(:,ke)) )

          call ATMOS_SATURATION_psat_liq( &
            elem%Np, 1, elem%Np, TEMP(:),     & ! (in)
            PSAT(:)                           ) ! (out)

          var_out(:,ke) = ( DDENS_(:,ke) + DENS_hyd(:,ke) ) * QTRC(iq_QV)%ptr%val(:,ke) &
                        / PSAT(:) * Rvap * TEMP(:) * 100.0_RP
        end do 
      end if
    
    case('ENGK')
      !$omp parallel do private (ke2D, DENS, mom_u1, mom_u2, G_11, G_12, G_22)
      do ke=1, lcmesh%Ne
        ke2D = lcmesh%EMap3Dto2D(ke)

        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        G_11(:) = lcmesh%G_ij(elem%IndexH2Dto3D,ke2D,1,1)
        G_12(:) = lcmesh%G_ij(elem%IndexH2Dto3D,ke2D,1,2)
        G_22(:) = lcmesh%G_ij(elem%IndexH2Dto3D,ke2D,2,2)

        mom_u1(:) = G_11(:) * MOMX_(:,ke) + G_12(:) * MOMY_(:,ke)
        mom_u2(:) = G_12(:) * MOMX_(:,ke) + G_22(:) * MOMY_(:,ke)

        var_out(:,ke) = 0.5_RP * ( MOMX_(:,ke) * mom_u1(:) + MOMY_(:,ke) * mom_u2(:) + MOMZ_(:,ke)**2 ) / DENS(:)
      end do
    
    case('ENGP')
      !$omp parallel do private (DENS)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = DENS(:) * Grav * lcmesh%zlev(:,ke)
      end do
    
    case('ENGI')
      !$omp parallel do private (DENS, iq)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = QDRY_(:,ke) * PRES_(:,ke) / Rtot(:,ke) * CVdry
        do iq = 1, QA
          var_out(:,ke) = var_out(:,ke) &
            + QTRC(iq)%ptr%val(:,ke) * ( PRES_(:,ke) / Rtot(:,ke) * TRACER_CV(iq) + DENS(:) * TRACER_ENGI0(iq) )
        end do
      end do
    
    case('ENGT')
      !$omp parallel do private (ke2D, DENS, mom_u1, mom_u2, iq, G_11, G_12, G_22)
      do ke=1, lcmesh%Ne
        ke2D = lcmesh%EMap3Dto2D(ke)

        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)

        G_11(:) = lcmesh%G_ij(elem%IndexH2Dto3D,ke2D,1,1)
        G_12(:) = lcmesh%G_ij(elem%IndexH2Dto3D,ke2D,1,2)
        G_22(:) = lcmesh%G_ij(elem%IndexH2Dto3D,ke2D,2,2)
        mom_u1(:) = G_11(:) * MOMX_(:,ke) + G_12(:) * MOMY_(:,ke)
        mom_u2(:) = G_12(:) * MOMX_(:,ke) + G_22(:) * MOMY_(:,ke)

        ! ENGI
        var_out(:,ke) = QDRY_(:,ke) * PRES_(:,ke) / Rtot(:,ke) * CVdry
        do iq = 1, QA
          var_out(:,ke) = var_out(:,ke) &
            + QTRC(iq)%ptr%val(:,ke) * ( PRES_(:,ke) / Rtot(:,ke) * TRACER_CV(iq) + DENS(:) * TRACER_ENGI0(iq) )
        end do
        ! ENGT
        var_out(:,ke) = &
            0.5_RP * ( MOMX_(:,ke) * mom_u1(:) + MOMY_(:,ke) * mom_u2(:) + MOMZ_(:,ke)**2 ) / DENS(:) & ! ENGK       
          + var_out(:,ke)                                                                             & ! ENGI
          + DENS(:) * Grav * lcmesh%pos_en(:,ke,3)                                                      ! ENGP
      end do
    
    case default
      LOG_ERROR("AtmosVars_calc_diagnoseVar_lc",*) 'The name of diagnostic variable is not suported. Check!', field_name
      call PRC_abort
    
    end select

    return
  end subroutine vars_calc_diagnoseVar_lc

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

    integer :: ke

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
      do ke=lcmesh%NeS, lcmesh%NeE
        var_out(:,ke) = SFLX_rain_MP%val(:,ke)
      end do
    case('SNOW')
      !$omp parallel do
      do ke=lcmesh%NeS, lcmesh%NeE
        var_out(:,ke) = SFLX_snow_MP%val(:,ke)
      end do
    case default
      LOG_ERROR("AtmosVars_calc_diagnoseVar2D_lc",*) 'The name of diagnostic variable is not suported. Check!', field_name
      call PRC_abort
    end select

    return
  end subroutine vars_calc_diagnoseVar2D_lc

!OCL SERIAL  
  subroutine vars_calc_specific_heat( this )
    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CPdry => CONST_CPdry,    &
      CVdry => CONST_CVdry,    &
      PRES00 => CONST_PRE00
    use scale_tracer, only: &
      TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP    
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat
    implicit none
    class(AtmosVars), intent(inout), target :: this

    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    integer :: varid
    integer :: ke
    integer :: iq

    class(ElementBase3D), pointer :: elem3D

    real(RP), allocatable :: q_tmp(:,:)
    !-------------------------------------------------------

    ! Calculate specific heat
    do n=1, this%AUX_VARS(1)%mesh%LOCAL_MESH_NUM
      lcmesh3D => this%AUX_VARS(1)%mesh%lcmesh_list(n)
      elem3D => lcmesh3D%refElem3D
      allocate( q_tmp(elem3D%Np,QA) )

      !$omp parallel do private(ke, iq, q_tmp)
      do ke = lcmesh3D%NeS, lcmesh3D%NeE
        do iq = 1, QA
          q_tmp(:,iq) = this%QTRC_VARS(iq)%local(n)%val(:,ke)
        end do
        call ATMOS_THERMODYN_specific_heat( &
          elem3D%Np, 1, elem3D%Np, QA,                                         & ! (in)
          q_tmp(:,:), TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! (in)
          this%AUX_VARS(AUXVAR_QDRY_ID )%local(n)%val(:,ke),            & ! (out)
          this%AUX_VARS(AUXVAR_Rtot_ID )%local(n)%val(:,ke),            & ! (out)
          this%AUX_VARS(AUXVAR_CVtot_ID)%local(n)%val(:,ke),            & ! (out)
          this%AUX_VARS(AUXVAR_CPtot_ID)%local(n)%val(:,ke)             ) ! (out)
      end do
      deallocate(q_tmp)
    end do
  end subroutine vars_calc_specific_heat

end module mod_atmos_vars