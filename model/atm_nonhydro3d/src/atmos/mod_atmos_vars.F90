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

  use scale_element_base, only: &
    ElementBase, ElementBase3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_meshfield_base, only: MeshFieldBase, MeshField3D
  
  use scale_file_restart_meshfield, only: &
    FILE_restart_meshfield_component
  use scale_file_common_meshfield, only: &
    DIMTYPE_XYZ  => FILE_COMMON_MESHFILED3D_DIMTYPEID_XYZ
  
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
  use scale_meshfieldcomm_base, only: MeshFieldContainer

  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo

  use mod_atmos_mesh, only: AtmosMesh
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: AtmosVars
    type(MeshField3D), allocatable :: PROG_VARS(:)
    type(ModelVarManager) :: PROGVARS_manager
    type(MeshFieldCommCubeDom3D) :: PROGVARS_comm
    
    type(MeshField3D), allocatable :: AUX_VARS(:)
    type(ModelVarManager) :: AUXVARS_manager 
    
    type(MeshField3D), allocatable :: PHY_TEND(:)
    type(ModelVarManager) :: PHYTENDS_manager 

    type(ModelVarManager) :: DIAGVARS_manager    
    integer, allocatable :: DIAGVARS_HISTID(:)

    type(FILE_restart_meshfield_component) :: restart_file
    
    logical :: check_range
    logical :: check_total
  contains
    procedure :: Init => AtmosVars_Init
    procedure :: Final => AtmosVars_Final
    procedure :: History => AtmosVars_History
    procedure :: Check   => AtmosVars_Check
    procedure :: Monitor => AtmosVars_Monitor
    procedure :: Read_restart_file => AtmosVar_Read_restart_file
    procedure :: Write_restart_file => AtmosVar_Write_restart_file
  end type AtmosVars

  public :: AtmosVars_GetLocalMeshPrgVar
  public :: AtmosVars_GetLocalMeshPrgVars
  public :: AtmosVars_GetLocalMeshPhyTends

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  ! Prognostic variables in dynamical process  

  integer, public, parameter :: ATMOS_PROGVARS_DDENS_ID = 1
  integer, public, parameter :: ATMOS_PROGVARS_MOMX_ID  = 2
  integer, public, parameter :: ATMOS_PROGVARS_MOMY_ID  = 3
  integer, public, parameter :: ATMOS_PROGVARS_MOMZ_ID  = 4
  integer, public, parameter :: ATMOS_PROGVARS_DRHOT_ID = 5
  integer, public, parameter :: ATMOS_PROGVARS_NUM      = 5

  type(VariableInfo), public :: ATMOS_PROGVARS_VINFO(ATMOS_PROGVARS_NUM)

  DATA ATMOS_PROGVARS_VINFO / &
    VariableInfo( ATMOS_PROGVARS_DDENS_ID, 'DDENS', 'deviation of density',       &
                  'kg/m3',  3, 'XYZ',  'air_density'                           ), &
    VariableInfo( ATMOS_PROGVARS_MOMX_ID , 'MOMX', 'momentum x',                  &
                  'kg/m2/s', 3, 'XYZ', 'upward_mass_flux_of_air'               ), &
    VariableInfo( ATMOS_PROGVARS_MOMY_ID , 'MOMY', 'momentum y',                  &
                  'kg/m2/s', 3, 'XYZ', 'eastward_mass_flux_of_air'             ), &
    VariableInfo( ATMOS_PROGVARS_MOMZ_ID , 'MOMZ', 'momentum z',                  &
                  'kg/m2/s', 3, 'XYZ', 'northward_mass_flux_of_air'            ), &
    VariableInfo( ATMOS_PROGVARS_DRHOT_ID, 'DRHOT', 'deviation of rho * theta',   &
                  'kg/m3*K', 3, 'XYZ',  ''                                    )   /

  real(RP), parameter :: PROGVARS_check_min(ATMOS_PROGVARS_NUM) = (/ -1.0_RP, -200.0_RP, -200.0_RP, -200.0_RP, -100.0_RP /)
  real(RP), parameter :: PROGVARS_check_max(ATMOS_PROGVARS_NUM) = (/  1.0_RP,  200.0_RP,  200.0_RP,  200.0_RP,  100.0_RP /)
              
  ! Reference state
  
  integer, public, parameter :: ATMOS_AUXVARS_PRESHYDRO_ID = 1
  integer, public, parameter :: ATMOS_AUXVARS_DENSHYDRO_ID = 2
  integer, public, parameter :: ATMOS_AUXVARS_NUM          = 2

  type(VariableInfo), public :: ATMOS_AUXVARS_VINFO(ATMOS_AUXVARS_NUM)
  DATA ATMOS_AUXVARS_VINFO / &
    VariableInfo( ATMOS_AUXVARS_PRESHYDRO_ID, 'PRES_hyd', 'hydrostatic part of pressure',  &
                  'Pa',  3, 'XYZ',  ''                                                  ), &
    VariableInfo( ATMOS_AUXVARS_DENSHYDRO_ID , 'DENS_hyd', 'hydrostatic part of density',  &
                  'kg/m3', 3, 'XYZ', ''                                                 )  /
  
  
  ! Tendency by physical processes
  
  integer, public, parameter :: ATMOS_PHYTEND_DENS_ID     = 1
  integer, public, parameter :: ATMOS_PHYTEND_MOMX_ID     = 2
  integer, public, parameter :: ATMOS_PHYTEND_MOMY_ID     = 3
  integer, public, parameter :: ATMOS_PHYTEND_MOMZ_ID     = 4
  integer, public, parameter :: ATMOS_PHYTEND_RHOH_ID     = 5
  integer, public, parameter :: ATMOS_PHYTEND_NUM         = 5

  type(VariableInfo), public :: ATMOS_PHYTEND_VINFO(ATMOS_PHYTEND_NUM)
  DATA ATMOS_PHYTEND_VINFO / &
    VariableInfo( ATMOS_PHYTEND_DENS_ID, 'DENS_tp', 'DENS_tp',                        &
                  'kg/m3/s',  3, 'XYZ',  'tendency of physical process for DENS' ),   &
    VariableInfo( ATMOS_PHYTEND_MOMX_ID, 'MOMX_tp', 'MOMX_tp',                        &
                  'kg/m2/s',  3, 'XYZ',  'tendency of physical process for MOMX' ),   &
    VariableInfo( ATMOS_PHYTEND_MOMY_ID, 'MOMY_tp', 'MOMY_tp',                        &
                  'kg/m2/s',  3, 'XYZ',  'tendency of physical process for MOMY' ),   &
    VariableInfo( ATMOS_PHYTEND_MOMZ_ID, 'MOMZ_tp', 'MOMZ_tp',                        &
                  'kg/m2/s',  3, 'XYZ',  'tendency of physical process for MOMZ' ),   &
    VariableInfo( ATMOS_PHYTEND_RHOH_ID, 'RHOT_tp', 'RHOT_tp',                        &
                  'kg/m3.K/s',  3, 'XYZ',  'tendency of physical process for RHOT' )  /

  ! Diagnostic variables

  integer, public, parameter :: ATMOS_DIAGVARS_U_ID      = 1
  integer, public, parameter :: ATMOS_DIAGVARS_V_ID      = 2
  integer, public, parameter :: ATMOS_DIAGVARS_W_ID      = 3
  integer, public, parameter :: ATMOS_DIAGVARS_PRES_ID   = 4
  integer, public, parameter :: ATMOS_DIAGVARS_T_ID      = 5  
  integer, public, parameter :: ATMOS_DIAGVARS_THETA_ID  = 6
  integer, public, parameter :: ATMOS_DIAGVARS_QDRY      = 7
  integer, public, parameter :: ATMOS_DIAGVARS_ENGT      = 8
  integer, public, parameter :: ATMOS_DIAGVARS_ENGP      = 9
  integer, public, parameter :: ATMOS_DIAGVARS_ENGK      = 10
  integer, public, parameter :: ATMOS_DIAGVARS_ENGI      = 11 
  integer, public, parameter :: ATMOS_DIAGVARS_NUM       = 11

  type(VariableInfo), public :: ATMOS_DIAGVARS_VINFO(ATMOS_DIAGVARS_NUM)
  DATA ATMOS_DIAGVARS_VINFO / &
    VariableInfo( ATMOS_DIAGVARS_U_ID    , 'U'   , 'velocity u'           , 'm/s'  , 3, 'XYZ', 'x_wind'               ), &
    VariableInfo( ATMOS_DIAGVARS_V_ID    , 'V'   , 'velocity v'           , 'm/s'  , 3, 'XYZ', 'y_wind'               ), &  
    VariableInfo( ATMOS_DIAGVARS_W_ID    , 'W'   , 'velocity w'           , 'm/s'  , 3, 'XYZ', 'upward_air_velocity'  ), &
    VariableInfo( ATMOS_DIAGVARS_PRES_ID , 'PRES', 'pressure'             , 'Pa'   , 3, 'XYZ', 'air_pressure'         ), &
    VariableInfo( ATMOS_DIAGVARS_T_ID    , 'T'   , 'temperature'          , 'K'    , 3, 'XYZ', 'air_temperature'      ), &
    VariableInfo( ATMOS_DIAGVARS_THETA_ID, 'PT'  , 'potential temperature', 'K'    , 3, 'XYZ', 'potential_temperature'), &
    Variableinfo( ATMOS_DIAGVARS_QDRY    , 'QDRY', 'dry air'              , 'kg/kg', 3, 'XYZ', ''                     ), &
    Variableinfo( ATMOS_DIAGVARS_ENGT    , 'ENGT', 'total energy'         , 'J/m3' , 3, 'XYZ', ''                     ), &
    Variableinfo( ATMOS_DIAGVARS_ENGP    , 'ENGP', 'potential energy'     , 'J/m3' , 3, 'XYZ', ''                     ), &
    Variableinfo( ATMOS_DIAGVARS_ENGK    , 'ENGK', 'kinetic energy'       , 'J/m3' , 3, 'XYZ', ''                     ), &
    Variableinfo( ATMOS_DIAGVARS_ENGI    , 'ENGI',  'internal energy'     , 'J/m3' , 3, 'XYZ', ''                    )  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  ! for monitor

  integer, private, parameter   :: IM_QDRY         =  1
  integer, private, parameter   :: IM_ENGT         =  2
  integer, private, parameter   :: IM_ENGP         =  3
  integer, private, parameter   :: IM_ENGK         =  4
  integer, private, parameter   :: IM_ENGI         =  5
  integer, private, parameter   :: DVM_nmax        =  5
  integer, private              :: DV_MONIT_id(DVM_nmax)

contains
  subroutine AtmosVars_Init( this, atm_mesh )

    use scale_file_monitor_meshfield, only:     &
      MONITOR_reg => FILE_monitor_meshfield_reg
    implicit none

    class(AtmosVars), target, intent(inout) :: this
    class(AtmosMesh), intent(in) :: atm_mesh

    integer :: iv
    integer :: n
    logical :: reg_file_hist

    type(MeshField3D) :: diag_vars(ATMOS_DIAGVARS_NUM)

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

    integer :: DV_id
    !--------------------------------------------------

    LOG_INFO('AtmosVars_Init',*)

    !- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_VARS)

    !- Initialize prognostic variables

    call this%PROGVARS_manager%Init()
    allocate( this%PROG_VARS(ATMOS_PROGVARS_NUM) )

    reg_file_hist = .true.    
    do iv = 1, ATMOS_PROGVARS_NUM

      call this%PROGVARS_manager%Regist(         &
        ATMOS_PROGVARS_VINFO(iv), atm_mesh%mesh, & ! (in) 
        this%PROG_VARS(iv),                      & ! (inout)
        reg_file_hist,  monitor_flag=.true.      ) ! (out)

      do n = 1, atm_mesh%mesh%LOCAL_MESH_NUM
        this%PROG_VARS(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    call this%PROGVARS_comm%Init(ATMOS_PROGVARS_NUM, 0, atm_mesh%mesh)
    call this%PROGVARS_manager%MeshFieldComm_Prepair( this%PROGVARS_comm, this%PROG_VARS(:) )

    LOG_NEWLINE
    LOG_INFO("ATMOS_vars_setup",*) 'List of prognostic variables (ATMOS) '
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, ATMOS_PROGVARS_NUM
      LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
      'NO.',iv,'|',ATMOS_PROGVARS_VINFO(iv)%NAME,'|', ATMOS_PROGVARS_VINFO(iv)%DESC,'[', ATMOS_PROGVARS_VINFO(iv)%UNIT,']'
    end do

    !- Initialize auxiliary variables

    call this%AUXVARS_manager%Init()
    allocate( this%AUX_VARS(ATMOS_AUXVARS_NUM) )
    
    reg_file_hist = .true.
    do iv = 1, ATMOS_AUXVARS_NUM
      call this%AUXVARS_manager%Regist(        &
        ATMOS_AUXVARS_VINFO(iv), atm_mesh%mesh, & ! (in) 
        this%AUX_VARS(iv), reg_file_hist        ) ! (out)
      do n = 1, atm_mesh%mesh%LOCAL_MESH_NUM
        this%AUX_VARS(iv)%local(n)%val(:,:) = 1.0_RP
      end do             
    end do

    !- Initialize the tendency of physical processes

    call this%PHYTENDS_manager%Init()
    allocate( this%PHY_TEND(ATMOS_PHYTEND_NUM) )
    
    reg_file_hist = .true.
    do iv = 1, ATMOS_PHYTEND_NUM
      call this%PHYTENDS_manager%Regist(       &
        ATMOS_PHYTEND_VINFO(iv), atm_mesh%mesh, & ! (in) 
        this%PHY_TEND(iv), reg_file_hist        ) ! (out)
      do n = 1, atm_mesh%mesh%LOCAL_MESH_NUM
        this%PHY_TEND(iv)%local(n)%val(:,:) = 0.0_RP
      end do             
    end do

    !- Initialize diagnostic variables for output
    call this%DIAGVARS_manager%Init()
    allocate( this%DIAGVARS_HISTID(ATMOS_DIAGVARS_NUM) )

    reg_file_hist = .true.
    do iv = 1, ATMOS_DIAGVARS_NUM
      call this%DIAGVARS_manager%Regist(        &
        ATMOS_DIAGVARS_VINFO(iv), atm_mesh%mesh, & ! (in) 
        diag_vars(iv), reg_file_hist             ) ! (out)
      
      this%DIAGVARS_HISTID(iv) = diag_vars(iv)%hist_id
    end do
    call this%DIAGVARS_manager%Final()

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
      call this%restart_file%Init('ATMOS',                          &
        IN_BASENAME, IN_POSTFIX_TIMELABEL,                          &
        OUT_BASENAME, OUT_POSTFIX_TIMELABEL, OUT_DTYPE, OUT_TITLE,  &
        ATMOS_PROGVARS_NUM + ATMOS_AUXVARS_NUM,                     &
        mesh3D=atm_mesh%mesh                    )
    else
      call this%restart_file%Init('ATMOS',      &
        ATMOS_PROGVARS_NUM + ATMOS_AUXVARS_NUM, &
        mesh3D=atm_mesh%mesh                    )
    end if

    !-----< monitor output setup >-----
    
    call MONITOR_reg( 'QDRY',         'dry air mass',          'kg', & ! (in)
                      DV_MONIT_id(IM_QDRY),                          & ! (out)
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

  subroutine AtmosVars_Final( this )
    implicit none
    class(AtmosVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosVars_Final',*)

    call this%restart_file%Final()

    call this%PROGVARS_comm%Final()
    call this%PROGVARS_manager%Final()
    call this%PHYTENDS_manager%Final()
    call this%AUXVARS_manager%Final()


    deallocate( this%DIAGVARS_HISTID )

    return
  end subroutine AtmosVars_Final

  subroutine AtmosVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosVars), intent(in) :: this
  
    integer :: n
    integer :: v
    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lcmesh
    integer :: hst_id

    type(MeshField3D) :: tmp_field
    !-------------------------------------------------------------------------

    mesh3D => this%PROG_VARS(1)%mesh
    do v = 1, ATMOS_PROGVARS_NUM
      hst_id = this%PROG_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%PROG_VARS(v) )
    end do

    do v = 1, ATMOS_AUXVARS_NUM
      hst_id = this%AUX_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%AUX_VARS(v) )
    end do

    call tmp_field%Init( "tmp_field", "", mesh3D)
    do v = 1, ATMOS_DIAGVARS_NUM
      hst_id = this%DIAGVARS_HISTID(v)
      if ( hst_id > 0 ) then
        call vars_calc_diagvar( this, ATMOS_DIAGVARS_VINFO(v)%NAME, tmp_field )
        call FILE_HISTORY_meshfield_put( hst_id, tmp_field )
      end if
    end do
    call tmp_field%Final()

    return
  end subroutine AtmosVars_history

  subroutine AtmosVar_Read_restart_file( this, atmos_mesh )

    use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
    use scale_meshfieldcomm_base, only: MeshFieldContainer
    implicit none
    
    class(AtmosVars), intent(inout), target :: this
    class(AtmosMesh), intent(in) :: atmos_mesh

    integer :: v

    type(MeshFieldCommCubeDom3D) :: auxvars_comm
    type(MeshFieldContainer) :: auxvars_comm_list(ATMOS_AUXVARS_NUM)
    !---------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOSVar_read_restart_file",*) 'Open restart file (ATMOS) '

    !- Open restart file
    call this%restart_file%Open()

    !- Read restart file
    do v=1, ATMOS_PROGVARS_NUM
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%PROG_VARS(v)%varname, &
        this%PROG_VARS(v)                                                      )
    end do
    do v=1, ATMOS_AUXVARS_NUM
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%AUX_VARS(v)%varname, &
        this%AUX_VARS(v)                                                      )
    end do

    !- Close restart file
    LOG_INFO("ATMOSVar_read_restart_file",*) 'Close restart file (ATMOS) '
    call this%restart_file%Close()

    !-- Check read data
    call this%Check( force = .true. )

    !- Communicate halo data of hydrostatic variables

    call auxvars_comm%Init( ATMOS_AUXVARS_NUM, 0, atmos_mesh%mesh )
    do v=1, ATMOS_AUXVARS_NUM
      auxvars_comm_list(v)%field3d => this%AUX_VARS(v)
    end do
    call auxvars_comm%Put( auxvars_comm_list, 1 )
    call auxvars_comm%Exchange()
    call auxvars_comm%Get( auxvars_comm_list, 1 )
    call auxvars_comm%Final()

    return
  end subroutine AtmosVar_Read_restart_file

  subroutine AtmosVar_write_restart_file( this )

    implicit none
    class(AtmosVars), intent(inout) :: this

    integer :: v, rf_vid
    !---------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("ATMOSVar_write_restart_file",*) 'Create restart file (ATMOS) '

    !- Check data which will be written to restart file
    call this%Check( force = .true. )

    !- Create restart file
    call this%restart_file%Create()

    !- Define variables
    do v=1, ATMOS_PROGVARS_NUM
      rf_vid = v
      call this%restart_file%Def_var( this%PROG_VARS(v),  &
        ATMOS_PROGVARS_VINFO(v)%DESC, rf_vid, DIMTYPE_XYZ )
    end do
    do v=1, ATMOS_AUXVARS_NUM
      rf_vid = ATMOS_PROGVARS_NUM + v
      call this%restart_file%Def_var( this%AUX_VARS(v),   &
        ATMOS_AUXVARS_VINFO(v)%DESC, rf_vid, DIMTYPE_XYZ  )
    end do
    call this%restart_file%End_def()

    !- Write restart file
    do v=1, ATMOS_PROGVARS_NUM
      rf_vid = v
      call this%restart_file%Write_var(rf_vid, this%PROG_VARS(v) )
    end do
    do v=1, ATMOS_AUXVARS_NUM
      rf_vid = ATMOS_PROGVARS_NUM + v
     call this%restart_file%Write_var(rf_vid, this%AUX_VARS(v) )
    end do

    !- Close restart file
    LOG_INFO("ATMOSVar_write_restart_file",*) 'Close restart file (ATMOS) '
    call this%restart_file%Close()

    return
  end subroutine AtmosVar_write_restart_file

  subroutine AtmosVars_Check( this, force )

    use scale_meshfield_statistics, only: &
      MeshField_statistics_total,         &
      MeshField_statistics_detail
    
    implicit none
    class(AtmosVars), intent(in) :: this
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
      do iv=1, ATMOS_PROGVARS_NUM
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
        call vel_fields(iv)%Init( ATMOS_DIAGVARS_VINFO(iv_diag)%NAME, "",  mesh3D )
        call vars_calc_diagvar( this, vel_fields(iv)%varname, vel_fields(iv) )
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

  subroutine AtmosVars_Monitor( this )
    use scale_file_monitor_meshfield, only: &
      FILE_monitor_meshfield_put
    
    implicit none
    class(AtmosVars), intent(in) :: this

    integer :: iv
    class(MeshBase3D), pointer :: mesh3D
    type(MeshField3D) :: work
    !--------------------------------------------------------------------------

    do iv=1, ATMOS_PROGVARS_NUM
      call FILE_monitor_meshfield_put( this%PROG_VARS(iv)%monitor_id, this%PROG_VARS(iv) )
    end do
  
    !##### Energy Budget #####
    mesh3D => this%PROG_VARS(1)%mesh
    call work%Init("tmp", "", mesh3D)

    if ( DV_MONIT_id(IM_ENGT) > 0 ) then
      call vars_calc_diagvar( this, 'ENGT', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGT), work )
    end if
    if ( DV_MONIT_id(IM_ENGP) > 0 ) then
      call vars_calc_diagvar( this, 'ENGP', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGP), work )
    end if
    if ( DV_MONIT_id(IM_ENGK) > 0 ) then
      call vars_calc_diagvar( this, 'ENGK', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGK), work )
    end if
    if ( DV_MONIT_id(IM_ENGI) > 0 ) then
      call vars_calc_diagvar( this, 'ENGI', work )
      call FILE_monitor_meshfield_put( DV_MONIT_id(IM_ENGI), work )
    end if

    call work%Final()

    return
  end subroutine AtmosVars_Monitor

  !----  Getter ---------------------------------------------------------------------------

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
      call auxvars_list%Get(ATMOS_AUXVARS_DENSHYDRO_ID, field)
      call field%GetLocalMeshField(domID, DENS_hyd)
    end if
    if (present(PRES_hyd)) then
      call auxvars_list%Get(ATMOS_AUXVARS_PRESHYDRO_ID, field)
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

  subroutine AtmosVars_GetLocalMeshPrgVars( domID, mesh, prgvars_list, auxvars_list, &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,                                                  &
    DENS_hyd, PRES_hyd, lcmesh3D                                                     &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    use scale_localmesh_base, only: LocalMeshBase
    use scale_localmesh_3d, only: LocalMesh3D

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer, intent(out) :: DENS_hyd, PRES_hyd
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call prgvars_list%Get(ATMOS_PROGVARS_DDENS_ID, field)
    call field%GetLocalMeshField(domID, DDENS)

    call prgvars_list%Get(ATMOS_PROGVARS_MOMX_ID, field)
    call field%GetLocalMeshField(domID, MOMX)
    
    call prgvars_list%Get(ATMOS_PROGVARS_MOMY_ID, field)
    call field%GetLocalMeshField(domID, MOMY)

    call prgvars_list%Get(ATMOS_PROGVARS_MOMZ_ID, field)
    call field%GetLocalMeshField(domID, MOMZ)

    call prgvars_list%Get(ATMOS_PROGVARS_DRHOT_ID, field)
    call field%GetLocalMeshField(domID, DRHOT)
  
    !--
    call auxvars_list%Get(ATMOS_AUXVARS_DENSHYDRO_ID, field)
    call field%GetLocalMeshField(domID, DENS_hyd)

    call auxvars_list%Get(ATMOS_AUXVARS_PRESHYDRO_ID, field)
    call field%GetLocalMeshField(domID, PRES_hyd)

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
  end subroutine AtmosVars_GetLocalMeshPrgVars

  subroutine AtmosVars_GetLocalMeshPhyTends( domID, mesh, phytends_list,  &
    DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOH_p,                           &
    lcmesh3D                                                              )

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: phytends_list
    class(LocalMeshFieldBase), pointer, intent(out) :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOH_p
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call phytends_list%Get(ATMOS_PHYTEND_DENS_ID, field)
    call field%GetLocalMeshField(domID, DENS_tp)

    call phytends_list%Get(ATMOS_PHYTEND_MOMX_ID, field)
    call field%GetLocalMeshField(domID, MOMX_tp)
    
    call phytends_list%Get(ATMOS_PHYTEND_MOMY_ID, field)
    call field%GetLocalMeshField(domID, MOMY_tp)

    call phytends_list%Get(ATMOS_PHYTEND_MOMZ_ID, field)
    call field%GetLocalMeshField(domID, MOMZ_tp)

    call phytends_list%Get(ATMOS_PHYTEND_RHOH_ID, field)
    call field%GetLocalMeshField(domID, RHOH_p)

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
  end subroutine AtmosVars_GetLocalMeshPhyTends

!-- private -----------------------------------------------------------------------
    
  subroutine vars_calc_diagvar( this, field_name, field_work ) 
    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CPdry => CONST_CPdry,    &
      CVdry => CONST_CVdry,    &
      PRES00 => CONST_PRE00

    implicit none
    class(AtmosVars), intent(in) :: this
    character(*), intent(in) :: field_name
    type(MeshField3D), intent(inout) :: field_work

    type(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    !--------------------------------------------------

    field_work%varname = field_name

    do n=1, field_work%mesh%LOCAL_MESH_NUM
      lcmesh3D => field_work%mesh%lcmesh_list(n)
      call vars_calc_diagnoseVar( field_name, field_work%local(n)%val,  &
        this%PROG_VARS(ATMOS_PROGVARS_DDENS_ID)%local(n)%val,           &
        this%PROG_VARS(ATMOS_PROGVARS_MOMX_ID)%local(n)%val,            &
        this%PROG_VARS(ATMOS_PROGVARS_MOMY_ID)%local(n)%val,            &
        this%PROG_VARS(ATMOS_PROGVARS_MOMZ_ID)%local(n)%val,            &
        this%PROG_VARS(ATMOS_PROGVARS_DRHOT_ID)%local(n)%val,           &
        this%AUX_VARS(ATMOS_AUXVARS_DENSHYDRO_ID)%local(n)%val,         & 
        this%AUX_VARS(ATMOS_AUXVARS_PRESHYDRO_ID)%local(n)%val,         &
        lcmesh3D, lcmesh3D%refElem3D )
    end do

    return
  end subroutine vars_calc_diagvar

  subroutine vars_calc_diagnoseVar( field_name, var_out,        &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,    &
    lcmesh, elem )

    use scale_const, only: &
      GRAV => CONST_GRAV,  &
      Rdry => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    implicit none

    type(LocalMesh3D), intent(in) :: lcmesh
    type(ElementBase3D), intent(in) :: elem
    character(*), intent(in) :: field_name
    real(RP), intent(out) :: var_out(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)

    integer :: ke
    real(RP) :: RHOT(elem%Np), DENS(elem%Np), PRES(elem%Np), THETA(elem%Np), TEMP(elem%Np)

    !-------------------------------------------------------------------------

    select case(trim(field_name))
    case('U')
      !$omp parallel do private (DENS)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = MOMX_(:,ke)/DENS(:)
      end do
    case('V')
      !$omp parallel do private (DENS)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = MOMY_(:,ke)/DENS(:)
      end do        
    case('W')
      !$omp parallel do private (DENS)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = MOMZ_(:,ke)/DENS(:)
      end do
    case('PRES')  
      !$omp parallel do private (RHOT)
      do ke=1, lcmesh%Ne
        RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
        var_out(:,ke) = PRES00 * (Rdry*RHOT(:)/PRES00)**(CPdry/Cvdry)
      end do
    case('T')
      !$omp parallel do private (RHOT, PRES)
      do ke=1, lcmesh%Ne
        RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
        PRES(:) = PRES00 * (Rdry*RHOT(:)/PRES00)**(CPdry/Cvdry)
        var_out(:,ke) = PRES(:) / (Rdry * (DDENS_(:,ke) + DENS_hyd(:,ke)) )
      end do
    case('PT')
      !$omp parallel do private (RHOT)
      do ke=1, lcmesh%Ne
        RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
        var_out(:,ke) = RHOT(:) / (DDENS_(:,ke) + DENS_hyd(:,ke))
      end do 
    case('ENGK')
      !$omp parallel do private (DENS)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = 0.5_RP * ( MOMX_(:,ke)**2 + MOMY_(:,ke)**2 + MOMZ_(:,ke)**2 ) / DENS(:)
      end do
    case('ENGP')
      !$omp parallel do private (DENS)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = DENS(:) * Grav * lcmesh%pos_en(:,ke,3)
      end do
    case('ENGI')
      !$omp parallel do private (RHOT, PRES)
      do ke=1, lcmesh%Ne
        RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
        PRES(:) = PRES00 * (Rdry*RHOT(:)/PRES00)**(CPdry/Cvdry)
        var_out(:,ke) = PRES(:) * CvDry / Rdry 
      end do
    case('ENGT')
      !$omp parallel do private (DENS, RHOT, PRES)
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
        PRES(:) = PRES00 * (Rdry*RHOT(:)/PRES00)**(CPdry/Cvdry)
        var_out(:,ke) = &
            0.5_RP * ( MOMX_(:,ke)**2 + MOMY_(:,ke)**2 + MOMZ_(:,ke)**2 ) / DENS(:)  & ! ENGK       
          + PRES(:) * CvDry / Rdry                                                   & ! ENGI
          + DENS(:) * Grav * lcmesh%pos_en(:,ke,3)                                     ! ENGP
      end do
    end select

    return
  end subroutine vars_calc_diagnoseVar  
end module mod_atmos_vars