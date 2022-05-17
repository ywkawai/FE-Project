!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_sw_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_debug
  use scale_const, only: &
    Grav => CONST_GRAV,     &
    RPlanet => CONST_RADIUS

  use scale_element_base, only: &
    ElementBase, ElementBase2D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: &
    MeshBase2D,                             &
    DIMTYPE_XY  => MeshBase2D_DIMTYPEID_XYT

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField2D  
  use scale_file_restart_meshfield, only: &
    FILE_restart_meshfield_component
  use scale_meshfieldcomm_cubedspheredom2d, only: MeshFieldCommCubedSphereDom2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer

  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo

  use mod_sw_mesh, only: SWMesh
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: SWVars
    type(MeshField2D), allocatable :: PROG_VARS(:)
    type(ModelVarManager) :: PROGVARS_manager
    type(MeshFieldCommCubedSphereDom2D) :: PROGVARS_comm
    
    type(ModelVarManager) :: TRCVARS_manager ! Dummy

    type(MeshField2D), allocatable :: AUX_VARS(:)
    type(ModelVarManager) :: AUXVARS_manager 
    type(MeshFieldCommCubedSphereDom2D) :: AUXVARS_comm
    
    type(MeshField2D), allocatable :: PHY_TEND(:)
    type(ModelVarManager) :: PHYTENDS_manager 

    type(ModelVarManager) :: DIAGVARS_manager    
    integer, allocatable :: DIAGVARS_HISTID(:)

    type(FILE_restart_meshfield_component) :: restart_file
    
    logical :: check_range
    logical :: check_total
  contains
    procedure :: Init => SWVars_Init
    procedure :: Final => SWVars_Final
    procedure :: Calc_diagnostics => SWVars_CalculateDiagnostics
    procedure :: Calc_Vorticity => SWVars_CalculateVor
    procedure :: History => SWVars_History
    procedure :: Check   => SWVars_Check
    procedure :: Monitor => SWVars_Monitor
    procedure :: Read_restart_file => SWVar_Read_restart_file
    procedure :: Write_restart_file => SWVar_Write_restart_file
  end type SWVars

  public :: SWVars_GetLocalMeshPrgVar
  public :: SWVars_GetLocalMeshPrgVars

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  ! Prognostic variables in dynamical process  

  integer, public, parameter :: SW_PROGVARS_h_ID  = 1
  integer, public, parameter :: SW_PROGVARS_U_ID  = 2
  integer, public, parameter :: SW_PROGVARS_V_ID  = 3
  integer, public, parameter :: SW_PROGVARS_NUM   = 3

  type(VariableInfo), public :: SW_PROGVARS_VINFO(SW_PROGVARS_NUM)

  DATA SW_PROGVARS_VINFO / &
    VariableInfo( SW_PROGVARS_h_ID, 'h', 'fluid thickness',                        &
                  'm',  2, 'XY',  'depth of fluid'                              ), &
    VariableInfo( SW_PROGVARS_U_ID , 'U', 'velocity x1 (contravariant vector)',    &
                  's-1', 2, 'XY', 'x1_component_contravariant_velocity_vector'  ), &
    VariableInfo( SW_PROGVARS_V_ID , 'V', 'velocity x2 (contravariant vector)',    &
                  's-1', 2, 'XY', 'x1_component_contravariant_velocity_vector'  )  / 

  real(RP) :: PROGVARS_check_min(SW_PROGVARS_NUM) 
  real(RP) :: PROGVARS_check_max(SW_PROGVARS_NUM)
              
  ! Reference state
  
  integer, public, parameter :: SW_AUXVARS_u1_ID   = 1
  integer, public, parameter :: SW_AUXVARS_u2_ID   = 2
  integer, public, parameter :: SW_AUXVARS_hs_ID   = 3
  integer, public, parameter :: SW_AUXVARS_VOR_ID  = 4
  integer, public, parameter :: SW_AUXVARS_NUM     = 4

  type(VariableInfo), public :: SW_AUXVARS_VINFO(SW_AUXVARS_NUM)
  DATA SW_AUXVARS_VINFO / &
    VariableInfo(   SW_AUXVARS_u1_ID, 'u1', 'velocity x1 (covariant vector)',    &
                  'm2/s', 2, 'XY', 'x1_component_covariant_velocity_vector'  ),  &
    VariableInfo(   SW_AUXVARS_u2_ID, 'u2', 'velocity x2 (covariant vector)',    &
                  'm2/s', 2, 'XY', 'x1_component_covariant_velocity_vector'  ),  & 
    VariableInfo(   SW_AUXVARS_hs_ID, 'hs', 'topography',                        &
                     'm', 2, 'XY', 'topography'                              ),  & 
    VariableInfo( SW_AUXVARS_VOR_ID, 'VOR', 'relative vorticity',                &
                    's-1', 2, 'XY', 'relative_vorticity'                     )   /
  
  ! Tendency by physical processes
  
  integer, public, parameter :: SW_PHYTEND_h_ID     = 1
  integer, public, parameter :: SW_PHYTEND_U_ID     = 2
  integer, public, parameter :: SW_PHYTEND_V_ID     = 3
  integer, public, parameter :: SW_PHYTEND_NUM      = 3

  type(VariableInfo), public :: SW_PHYTEND_VINFO(SW_PHYTEND_NUM)
  DATA SW_PHYTEND_VINFO / &
    VariableInfo( SW_PHYTEND_h_ID, 'h_tp', 'h_tp',                            &
                  'm/s',  2, 'XY',  'tendency of physical process for h'   ), &
    VariableInfo( SW_PHYTEND_U_ID, 'V_tp', 'U_tp',                            &
                  's-2',  2, 'XY',  'tendency of physical process for U'   ), &
    VariableInfo( SW_PHYTEND_V_ID, 'U_tp', 'V_tp',                            &
                  's-2',  2, 'XY',  'tendency of physical process for V' )    /

  ! Diagnostic variables

  integer, public, parameter :: SW_DIAGVARS_Vellon_ID   = 1
  integer, public, parameter :: SW_DIAGVARS_Vellat_ID   = 2
  integer, public, parameter :: SW_DIAGVARS_Height_ID   = 3
  integer, public, parameter :: SW_DIAGVARS_ENGT        = 4
  integer, public, parameter :: SW_DIAGVARS_ENGP        = 5
  integer, public, parameter :: SW_DIAGVARS_ENGK        = 6
  integer, public, parameter :: SW_DIAGVARS_NUM         = 6

  type(VariableInfo), public :: SW_DIAGVARS_VINFO(SW_DIAGVARS_NUM)
  DATA SW_DIAGVARS_VINFO / &
    VariableInfo( SW_DIAGVARS_Vellon_ID, 'Vel_lon', 'velocity u'    , 'm/s'  , 2, 'XY', 'lon_wind' ), &
    VariableInfo( SW_DIAGVARS_Vellat_ID, 'Vel_lat', 'velocity v'    , 'm/s'  , 2, 'XY', 'lat_wind' ), &  
    VariableInfo( SW_DIAGVARS_Height_ID,  'Height', 'Height'        , 'm'    , 2, 'XY', 'height'   ), &  
    Variableinfo( SW_DIAGVARS_ENGT   , 'ENGT', 'total energy'          , 'J/m3' , 2, 'XY', ''      ), &
    Variableinfo( SW_DIAGVARS_ENGP   , 'ENGP', 'potential energy'      , 'J/m3' , 2, 'XY', ''      ), &
    Variableinfo( SW_DIAGVARS_ENGK   , 'ENGK', 'kinetic energy'        , 'J/m3' , 2, 'XY', ''      )  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  ! for monitor

  integer, private, parameter   :: IM_ENGT         =  1
  integer, private, parameter   :: IM_ENGP         =  2
  integer, private, parameter   :: IM_ENGK         =  3
  integer, private, parameter   :: DVM_nmax        =  3
  integer, private              :: DV_MONIT_id(DVM_nmax)

contains
  subroutine SWVars_Init( this, sw_mesh )

    use scale_file_monitor_meshfield, only:     &
      MONITOR_reg => FILE_monitor_meshfield_reg
    implicit none

    class(SWVars), target, intent(inout) :: this
    class(SWMesh), intent(in) :: sw_mesh

    integer :: iv
    integer :: n
    logical :: reg_file_hist

    type(MeshField2D) :: diag_vars(SW_DIAGVARS_NUM)

    logical :: CHECK_RANGE    = .false.
    logical :: CHECK_TOTAL    = .false.

    namelist / PARAM_SW_VARS / &
      CHECK_RANGE, &
      CHECK_TOTAL

    character(len=H_LONG) :: IN_BASENAME           = ''        !< Basename of the input  file
    logical :: IN_POSTFIX_TIMELABEL                = .false.   !< Add timelabel to the basename of input  file?
    character(len=H_LONG) :: OUT_BASENAME          = ''        !< Basename of the output file
    logical :: OUT_POSTFIX_TIMELABEL               = .true.    !< Add timelabel to the basename of output file?
    character(len=H_MID) :: OUT_TITLE              = ''        !< Title    of the output file
    character(len=H_SHORT) :: OUT_DTYPE            = 'DEFAULT' !< REAL4 or REAL8  

    namelist / PARAM_SW_VARS_RESTART / &
      IN_BASENAME,           &
      IN_POSTFIX_TIMELABEL,  &
      OUT_BASENAME,          &
      OUT_POSTFIX_TIMELABEL, &
      OUT_TITLE,             &
      OUT_DTYPE    
    
    integer :: ierr
    logical :: is_specified

    logical :: monitor_flag(SW_PROGVARS_NUM)
    integer :: DV_id
    !--------------------------------------------------

    LOG_INFO('SWVars_Init',*)

    !- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SW_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SW_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SW_vars_setup",*) 'Not appropriate names in namelist PARAM_SW_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_SW_VARS)

    !- Initialize prognostic variables

    call this%PROGVARS_manager%Init()
    allocate( this%PROG_VARS(SW_PROGVARS_NUM) )

    reg_file_hist = .true.  
    monitor_flag(:) = .false.  
    monitor_flag(SW_PROGVARS_h_ID) = .true.

    do iv = 1, SW_PROGVARS_NUM

      call this%PROGVARS_manager%Regist( &
        SW_PROGVARS_VINFO(iv), sw_mesh%mesh, & ! (in) 
        this%PROG_VARS(iv),                  & ! (inout)
        reg_file_hist, monitor_flag(iv)      ) ! (in)

      do n = 1, sw_mesh%mesh%LOCAL_MESH_NUM
        this%PROG_VARS(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    call this%PROGVARS_comm%Init( 1, 1, sw_mesh%mesh )
    call this%PROGVARS_manager%MeshFieldComm_Prepair( this%PROGVARS_comm, this%PROG_VARS(:) )

    LOG_NEWLINE
    LOG_INFO("SW_vars_setup",*) 'List of prognostic variables (SW) '
    LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
               '      |', 'VARNAME                 ','|', &
               'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
    do iv = 1, SW_PROGVARS_NUM
      LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
      'NO.',iv,'|',SW_PROGVARS_VINFO(iv)%NAME,'|', SW_PROGVARS_VINFO(iv)%DESC,'[', SW_PROGVARS_VINFO(iv)%UNIT,']'
    end do

    !- Initialize auxiliary variables

    call this%AUXVARS_manager%Init()
    allocate( this%AUX_VARS(SW_AUXVARS_NUM) )
    
    reg_file_hist = .true.
    do iv = 1, SW_AUXVARS_NUM
      call this%AUXVARS_manager%Regist(         &
        SW_AUXVARS_VINFO(iv), sw_mesh%mesh,     & ! (in) 
        this%AUX_VARS(iv), reg_file_hist        ) ! (out)
      do n = 1, sw_mesh%mesh%LOCAL_MESH_NUM
        this%AUX_VARS(iv)%local(n)%val(:,:) = 1.0_RP
      end do             
    end do

    call this%AUXVARS_comm%Init(SW_AUXVARS_NUM, 0, sw_mesh%mesh)
    call this%AUXVARS_manager%MeshFieldComm_Prepair( this%AUXVARS_comm, this%AUX_VARS(:) )

    call this%PROGVARS_comm%SetCovariantVec( 1, &
      this%AUX_VARS(SW_AUXVARS_u1_ID),          &
      this%AUX_VARS(SW_AUXVARS_u2_ID)           )
    
    !- Initialize the tendency of physical processes

    call this%PHYTENDS_manager%Init()
    allocate( this%PHY_TEND(SW_PHYTEND_NUM) )
    
    reg_file_hist = .true.
    do iv = 1, SW_PHYTEND_NUM
      call this%PHYTENDS_manager%Regist(       &
        SW_PHYTEND_VINFO(iv), sw_mesh%mesh,    & ! (in) 
        this%PHY_TEND(iv), reg_file_hist       ) ! (out)
      do n = 1, sw_mesh%mesh%LOCAL_MESH_NUM
        this%PHY_TEND(iv)%local(n)%val(:,:) = 0.0_RP
      end do             
    end do

    !- Initialize diagnostic variables for output
    call this%DIAGVARS_manager%Init()
    allocate( this%DIAGVARS_HISTID(SW_DIAGVARS_NUM) )

    reg_file_hist = .true.
    do iv = 1, SW_DIAGVARS_NUM
      call this%DIAGVARS_manager%Regist(        &
        SW_DIAGVARS_VINFO(iv), sw_mesh%mesh,    & ! (in) 
        diag_vars(iv), reg_file_hist            ) ! (out)
      
      this%DIAGVARS_HISTID(iv) = diag_vars(iv)%hist_id
    end do
    call this%DIAGVARS_manager%Final()

    !-- Setup information for input/output restart files. 

    is_specified = .true.
    !- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SW_VARS_RESTART,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SWVars_Init",*) 'Not found namelist PARAM_SW_VARS_RESTART. Default used.'
       is_specified = .false.
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SWVars_Init",*) 'Not appropriate names in namelist PARAM_SW_VARS_RESTART. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_SW_VARS_RESTART)

    if (is_specified) then
      call this%restart_file%Init('SW',                             &
        IN_BASENAME, IN_POSTFIX_TIMELABEL,                          &
        OUT_BASENAME, OUT_POSTFIX_TIMELABEL, OUT_DTYPE, OUT_TITLE,  &
        SW_PROGVARS_NUM + SW_AUXVARS_NUM,                           &
        meshCubedSphere2D=sw_mesh%mesh                              )
    else
      call this%restart_file%Init('SW',      &
        SW_PROGVARS_NUM + SW_AUXVARS_NUM,    &
        meshCubedSphere2D=sw_mesh%mesh       )
    end if

    !-----< monitor output setup >-----
    
    call MONITOR_reg( 'ENGT',         'total     energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGT),                          & ! (out)
                      dim_type='ATM2D', is_tendency=.false.          ) ! (in)
    call MONITOR_reg( 'ENGP',         'potential energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGP),                          & ! (out)
                      dim_type='ATM2D', is_tendency=.false.          ) ! (in)
    call MONITOR_reg( 'ENGK',         'kinetic   energy',       'J', & ! (in)
                      DV_MONIT_id(IM_ENGK),                          & ! (out)
                      dim_type='ATM2D', is_tendency=.false.          ) ! (in)

    !-----< check the range of values >-----

    this%check_range = CHECK_RANGE
    this%check_total = CHECK_TOTAL
    LOG_INFO("SW_vars_setup",*) 'Check value range of variables?     : ', CHECK_RANGE
    LOG_INFO("SW_vars_setup",*) 'Check total value of variables?     : ', CHECK_TOTAL
      
    PROGVARS_check_min(:) = (/ -10000.0_RP, -1.0_RP, -1.0_RP  /)
    PROGVARS_check_max(:) = (/  15000.0_RP,  1.0_RP,  1.0_RP  /)

    return
  end subroutine SWVars_Init

  subroutine SWVars_Final( this )
    implicit none
    class(SWVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('SWVars_Final',*)

    call this%restart_file%Final()

    call this%PROGVARS_comm%Final()
    call this%AUXVARS_comm%Final()

    call this%PROGVARS_manager%Final()
    call this%AUXVARS_manager%Final()
    call this%PHYTENDS_manager%Final()

    deallocate( this%DIAGVARS_HISTID )

    return
  end subroutine SWVars_Final

  subroutine SWVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(SWVars), intent(inout) :: this
  
    integer :: n
    integer :: v
    class(MeshBase2D), pointer :: mesh2D
    class(LocalMesh2D), pointer :: lcmesh
    integer :: hst_id

    type(MeshField2D) :: tmp_field
    !-------------------------------------------------------------------------

    mesh2D => this%PROG_VARS(1)%mesh
    do v = 1, SW_PROGVARS_NUM
      hst_id = this%PROG_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%PROG_VARS(v) )
    end do

    do v = 1, SW_AUXVARS_NUM
      hst_id = this%AUX_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%AUX_VARS(v) )
    end do

    do v = 1, SW_PHYTEND_NUM
      hst_id = this%PHY_TEND(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%PHY_TEND(v) )
    end do

    call tmp_field%Init( "tmp_field", "", mesh2D)
    do v = 1, SW_DIAGVARS_NUM
      hst_id = this%DIAGVARS_HISTID(v)
      if ( hst_id > 0 ) then
        call vars_calc_diagvar( this, SW_DIAGVARS_VINFO(v)%NAME, tmp_field )
        call FILE_HISTORY_meshfield_put( hst_id, tmp_field )
      end if
    end do
    call tmp_field%Final()

    return
  end subroutine SWVars_history

  subroutine SWVar_Read_restart_file( this, sw_mesh )

    use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
    use scale_meshfieldcomm_base, only: MeshFieldContainer
    implicit none
    
    class(SWVars), intent(inout), target :: this
    class(SWMesh), intent(in) :: sw_mesh

    integer :: v
    !---------------------------------------

    LOG_NEWLINE
    LOG_INFO("SWVar_read_restart_file",*) 'Open restart file (SW) '

    !- Open restart file
    call this%restart_file%Open()

    !- Read restart file
    do v=1, SW_PROGVARS_NUM
      call this%restart_file%Read_var( DIMTYPE_XY, this%PROG_VARS(v)%varname, &
        this%PROG_VARS(v)                                                     )
    end do
    do v=1, SW_AUXVARS_hs_ID
      call this%restart_file%Read_var( DIMTYPE_XY, this%AUX_VARS(v)%varname, &
        this%AUX_VARS(v)                                                     )
    end do

    !- Close restart file
    LOG_INFO("SWVar_read_restart_file",*) 'Close restart file (SW) '
    call this%restart_file%Close()

    !-- Check read data
    call this%Check( force = .true. )

    !-- Calculate diagnostic variables
    call this%Calc_diagnostics( sw_mesh )   

    !-- Communicate halo data of hydrostatic & diagnostic variables
    call this%AUXVARS_manager%MeshFieldComm_Exchange()
    
    return
  end subroutine SWVar_Read_restart_file

  subroutine SWVar_write_restart_file( this )

    implicit none
    class(SWVars), intent(inout) :: this

    integer :: v, rf_vid
    !---------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("SWVar_write_restart_file",*) 'Create restart file (SW) '

    !- Check data which will be written to restart file
    call this%Check( force = .true. )

    !- Create restart file
    call this%restart_file%Create()

    !- Define variables
    do v=1, SW_PROGVARS_NUM
      rf_vid = v
      call this%restart_file%Def_var( this%PROG_VARS(v), &
        SW_PROGVARS_VINFO(v)%DESC, rf_vid, DIMTYPE_XY    )
    end do
    do v=1, SW_AUXVARS_hs_ID
      rf_vid = SW_PROGVARS_NUM + v
      call this%restart_file%Def_var( this%AUX_VARS(v), &
        SW_AUXVARS_VINFO(v)%DESC, rf_vid, DIMTYPE_XY    )
    end do
    call this%restart_file%End_def()

    !- Write restart file
    do v=1, SW_PROGVARS_NUM
      rf_vid = v
      call this%restart_file%Write_var( rf_vid, this%PROG_VARS(v) )
    end do
    do v=1, SW_AUXVARS_hs_ID
      rf_vid = SW_PROGVARS_NUM + v
     call this%restart_file%Write_var( rf_vid, this%AUX_VARS(v) )
    end do

    !- Close restart file
    LOG_INFO("SWVar_write_restart_file",*) 'Close restart file (SW) '
    call this%restart_file%Close()

    return
  end subroutine SWVar_write_restart_file

  subroutine SWVars_Check( this, force )

    use scale_meshfield_statistics, only: &
      MeshField_statistics_total,         &
      MeshField_statistics_detail
    
    implicit none
    class(SWVars), intent(in) :: this
    logical, intent(in), optional :: force

    integer :: iv
    integer :: iv_diag
    integer :: n
    logical  :: check

    class(MeshBase2D), pointer :: mesh2D
    class(LocalMeshBase), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: lcfield
    type(ElementBase), pointer :: elem
    character(len=H_MID) :: varname

    type(MeshField2D) :: vel_fields(2)
    type(MeshField2D) :: work
    !--------------------------------------------------------------------------

    if ( present(force) ) then
      check = force
    else
      check = this%check_range
    end if

    if (check) then
      do iv=1, SW_PROGVARS_NUM
        mesh2D => this%PROG_VARS(iv)%mesh
        do n=1, mesh2D%LOCAL_MESH_NUM
          lcmesh => mesh2D%lcmesh_list(n)
          elem => lcmesh%refElem
          call this%PROG_VARS(iv)%GetLocalMeshField(n, lcfield)
          write(varname,'(a,i3.3,a)') this%PROG_VARS(iv)%varname//'(domID=', n, ')' 
          call VALCHECK( elem%Np, 1, elem%Np, lcmesh%NeA, lcmesh%NeS, lcmesh%NeE, lcfield%val(:,:), &
            PROGVARS_check_min(iv), PROGVARS_check_max(iv), trim(varname), __FILE__, __LINE__       )
        end do
      end do

      mesh2D => this%PROG_VARS(1)%mesh
      do iv=1, 2
        iv_diag = SW_DIAGVARS_Vellon_ID + iv - 1
        call vel_fields(iv)%Init( SW_DIAGVARS_VINFO(iv_diag)%NAME, "",  mesh2D )
        call vars_calc_diagvar( this, vel_fields(iv)%varname, vel_fields(iv) )
      end do
      call MeshField_statistics_detail( vel_fields(:) )
      do iv=1, 2
        call vel_fields(iv)%Final()
      end do
    end if

    if ( present(force) ) then
      check = force
    else
      check = this%check_total
    end if
    if (check) then
      mesh2D => this%PROG_VARS(1)%mesh
      call work%Init("tmp", "", mesh2D)
      call work%Final()
    end if

    return
  end subroutine SWVars_Check

  subroutine SWVars_Monitor( this )
    use scale_file_monitor_meshfield, only: &
      FILE_monitor_meshfield_put
    
    implicit none
    class(SWVars), intent(in) :: this

    integer :: iv
    class(MeshBase2D), pointer :: mesh2D
    type(MeshField2D) :: work
    !--------------------------------------------------------------------------

    do iv=1, SW_PROGVARS_h_ID
      call FILE_monitor_meshfield_put( this%PROG_VARS(iv)%monitor_id, this%PROG_VARS(iv) )
    end do
  
    !##### Energy Budget #####
    mesh2D => this%PROG_VARS(1)%mesh
    call work%Init("tmp", "", mesh2D)

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

    call work%Final()

    return
  end subroutine SWVars_Monitor

  !----  Getter ---------------------------------------------------------------------------

  subroutine SWVars_GetLocalMeshPrgVar( domID, mesh, prgvars_list, auxvars_list,    &
     varid,                                                                         &
     var, hs, lcmesh2D                                                              )

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    integer, intent(in) :: varid
    class(LocalMeshFieldBase), pointer, intent(out) :: var
    class(LocalMeshFieldBase), pointer, intent(out), optional :: hs
    class(LocalMesh2D), pointer, intent(out), optional :: lcmesh2D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call prgvars_list%Get(varid, field)
    call field%GetLocalMeshField(domID, var)

    if (present(hs)) then
      call auxvars_list%Get(SW_AUXVARS_hs_ID, field)
      call field%GetLocalMeshField(domID, hs)
    end if

    if (present(lcmesh2D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh2D )

      select type(lcmesh)
      type is (LocalMesh2D)
        if (present(lcmesh2D)) lcmesh2D => lcmesh
      end select
    end if

    return
  end subroutine SWVars_GetLocalMeshPrgVar

  subroutine SWVars_GetLocalMeshPrgVars( domID, mesh, prgvars_list, auxvars_list, &
    h, U, V,                                                                      &
    hs, u1, u2, lcmesh2D                                                          &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    use scale_localmesh_base, only: LocalMeshBase
    use scale_localmesh_2d, only: LocalMesh2D

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: h, U, V
    class(LocalMeshFieldBase), pointer, intent(out) :: hs, u1, u2
    class(LocalMesh2D), pointer, intent(out), optional :: lcmesh2D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call prgvars_list%Get(SW_PROGVARS_h_ID, field)
    call field%GetLocalMeshField(domID, h)

    call prgvars_list%Get(SW_PROGVARS_U_ID, field)
    call field%GetLocalMeshField(domID, U)
    
    call prgvars_list%Get(SW_PROGVARS_V_ID, field)
    call field%GetLocalMeshField(domID, V)
  
    !--
    call auxvars_list%Get(SW_AUXVARS_hs_ID, field)
    call field%GetLocalMeshField(domID, hs)

    call auxvars_list%Get(SW_AUXVARS_u1_ID, field)
    call field%GetLocalMeshField(domID, u1)

    call auxvars_list%Get(SW_AUXVARS_u2_ID, field)
    call field%GetLocalMeshField(domID, u2)

    !---
    
    if (present(lcmesh2D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh2D )

      select type(lcmesh)
      type is (LocalMesh2D)
        if (present(lcmesh2D)) lcmesh2D => lcmesh
      end select
    end if

    return
  end subroutine SWVars_GetLocalMeshPrgVars

  !-----------------------------------------------------------------------------
  !> Calculate diagnostic variables  
  subroutine SWVars_CalculateDiagnostics( this, model_mesh )
    implicit none
    class(SWVars), intent(inout), target :: this
    class(SWMesh), intent(in) :: model_mesh

    type(LocalMesh2D), pointer :: lcmesh2D
    integer :: n
    integer :: varid

    class(MeshField2D), pointer :: field
    !-------------------------------------------------------

    call this%PROGVARS_manager%MeshFieldComm_Exchange()

    do varid=SW_AUXVARS_VOR_ID+1, SW_AUXVARS_NUM
      field => this%AUX_VARS(varid)
      do n=1, field%mesh%LOCAL_MESH_NUM
        lcmesh2D => field%mesh%lcmesh_list(n)
        call vars_calc_diagnoseVar( &
          field%varname, field%local(n)%val,                  &
          this%PROG_VARS(SW_PROGVARS_h_ID )%local(n)%val,     &
          this%PROG_VARS(SW_PROGVARS_U_ID)%local(n)%val,      &
          this%PROG_VARS(SW_PROGVARS_V_ID)%local(n)%val,      &
          this%AUX_VARS(SW_AUXVARS_hs_ID)%local(n)%val,       & 
          this%AUX_VARS(SW_AUXVARS_u1_ID )%local(n)%val,      &
          this%AUX_VARS(SW_AUXVARS_u2_ID )%local(n)%val,      &
          lcmesh2D, lcmesh2D%refElem2D )
      end do
    end do

    ! VOR
    call this%Calc_Vorticity( this%AUX_VARS(SW_AUXVARS_VOR_ID), model_mesh )

    return
  end subroutine SWVars_CalculateDiagnostics

  subroutine SWVars_CalculateVor( this, vor, model_mesh )
    implicit none
    class(SWVars), intent(in) :: this
    class(SWMesh), intent(in) :: model_mesh    
    class(MeshField2D), intent(inout) :: vor

    type(LocalMesh2D), pointer :: lcmesh2D    
    integer :: n
    !-------------------------------------

    do n=1, vor%mesh%LOCAL_MESH_NUM
      lcmesh2D => vor%mesh%lcmesh_list(n)

      call vars_eval_vor_lc( vor%local(n)%val, &
        this%PROG_VARS(SW_PROGVARS_U_ID)%local(n)%val, this%PROG_VARS(SW_PROGVARS_V_ID)%local(n)%val, &
        this%AUX_VARS(SW_AUXVARS_u1_ID )%local(n)%val, this%AUX_VARS(SW_AUXVARS_u2_ID )%local(n)%val, &
        model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%LiftOptrMat,                       &
        lcmesh2D, lcmesh2D%refElem2D )
    end do

    return
  end subroutine SWVars_CalculateVor

!-- private -----------------------------------------------------------------------
    
!OCL SERIAL
  subroutine vars_calc_diagvar( this, field_name, field_work ) 
    implicit none
    
    class(SWVars), intent(in) :: this
    character(*), intent(in) :: field_name
    type(MeshField2D), intent(inout) :: field_work

    type(LocalMesh2D), pointer :: lcmesh2D
    integer :: n
    !--------------------------------------------------

    field_work%varname = field_name

    do n=1, field_work%mesh%LOCAL_MESH_NUM
      lcmesh2D => field_work%mesh%lcmesh_list(n)
      call vars_calc_diagnoseVar( &
        field_name, field_work%local(n)%val,                &
        this%PROG_VARS(SW_PROGVARS_h_ID )%local(n)%val,     &
        this%PROG_VARS(SW_PROGVARS_U_ID)%local(n)%val,      &
        this%PROG_VARS(SW_PROGVARS_V_ID)%local(n)%val,      &
        this%AUX_VARS(SW_AUXVARS_hs_ID)%local(n)%val,       & 
        this%AUX_VARS(SW_AUXVARS_u1_ID )%local(n)%val,      &
        this%AUX_VARS(SW_AUXVARS_u2_ID )%local(n)%val,      &
        lcmesh2D, lcmesh2D%refElem2D )
    end do

    return
  end subroutine vars_calc_diagvar

!OCL SERIAL
  subroutine vars_eval_vor_lc( vor,           &
    U, V, u1, u2, Dx, Dy, Lift, lcmesh, elem  )
      
    use scale_sparsemat, only: sparsemat, &
      sparsemat_matmul
    implicit none

    type(LocalMesh2D), intent(in) :: lcmesh
    type(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: vor(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: U(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: V(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: u1(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: u2(elem%Np,lcmesh%NeA)
    type(SparseMat), intent(in) :: Dx, Dy
    type(SparseMat), intent(in) :: Lift

    integer :: ke
    real(RP) :: del_flux(elem%NfpTot,lcmesh%Ne,1)
    real(RP) :: Fx(elem%Np), Fy(elem%Np), LiftDelFlx(elem%Np)
    !----------------------------------------------------

    call eval_vor_del_flux( del_flux, &
      U, V, u1, u2, lcmesh%Gsqrt, lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), &
      lcmesh%VMapM, lcmesh%VMapP, lcmesh, lcmesh%refElem2D )
    
    !$omp parallel do private(Fx, Fy, LiftDelFlx)
    do ke=lcmesh%NeS, lcmesh%NeE
      call sparsemat_matmul(Dx, u2(:,ke), Fx)
      call sparsemat_matmul(Dy, u1(:,ke), Fy)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke,1), LiftDelFlx)

      vor(:,ke) = ( lcmesh%Escale(:,ke,1,1) * Fx(:) &
                  - lcmesh%Escale(:,ke,2,2) * Fy(:) &
                  + LiftDelFlx(:) ) / lcmesh%Gsqrt(:,ke)
    end do

    return
  end subroutine vars_eval_vor_lc

!OCL SERIAL
  subroutine eval_vor_del_flux( del_flux_aux,  &
    U_, V_, u1_, u2_,                          &
    Gsqrt_, nx, ny, vmapM, vmapP, lmesh, elem  )

    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux_aux(elem%NfpTot,lmesh%Ne,1)
    real(RP), intent(in) ::  U_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  V_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  u1_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  u2_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt_(elem%Np*lmesh%Ne) 
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    
    integer :: ke, iP(elem%NfpTot), iM(elem%NfpTot)
    !------------------------------------------------------------------------

    !$omp parallel do private( iM, iP )
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)

      del_flux_aux(:,ke,1) = 0.5_RP * ( &
                     ( u2_(iP) - u2_(iM) ) * nx(:,ke) &
                   - ( u1_(iP) - u1_(iM) ) * ny(:,ke) )
    end do

    return
  end subroutine eval_vor_del_flux

!OCL SERIAL
  subroutine vars_calc_diagnoseVar( field_name, var_out,  &
    h_, U_, V_, hs_, u1_, u2_,                            &
    lcmesh, elem )

    use scale_const, only: &
      GRAV => CONST_GRAV,  &
      Rdry => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    use scale_cubedsphere_cnv, only: &
      CubedSphereCnv_CS2LonLatVec
    implicit none

    type(LocalMesh2D), intent(in) :: lcmesh
    type(ElementBase2D), intent(in) :: elem
    character(*), intent(in) :: field_name
    real(RP), intent(out) :: var_out(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: h_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: U_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: V_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: hs_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: u1_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: u2_(elem%Np,lcmesh%NeA)

    integer :: ke
    real(RP) :: dummy(elem%Np,lcmesh%NeA)
    !-------------------------------------------------------------------------

    select case(trim(field_name))
    case('Vel_lon')
      call CubedSphereCnv_CS2LonLatVec( &
        lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), elem%Np * lcmesh%Ne, RPlanet, & ! (in)
        U_(:,lcmesh%NeS:lcmesh%NeE), V_(:,lcmesh%NeS:lcmesh%NeE),                                 & ! (in)
        var_out(:,lcmesh%NeS:lcmesh%NeE), dummy(:,lcmesh%NeS:lcmesh%NeE)                          ) ! (out)
       !$omp parallel do
        do ke=1, lcmesh%Ne
          var_out(:,ke) = var_out(:,ke) * cos(lcmesh%lat(:,ke))
        end do        
    case('Vel_lat')
      call CubedSphereCnv_CS2LonLatVec( &
        lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), elem%Np * lcmesh%Ne, RPlanet, & ! (in)
        U_(:,lcmesh%NeS:lcmesh%NeE), V_(:,lcmesh%NeS:lcmesh%NeE),                                 & ! (in)
        dummy(:,lcmesh%NeS:lcmesh%NeE), var_out(:,lcmesh%NeS:lcmesh%NeE)                          ) ! (out)
    case('Height')
      !$omp parallel do
      do ke=1, lcmesh%Ne
        var_out(:,ke) = h_(:,ke) + hs_(:,ke)
      end do
    case('ENGK')
      !$omp parallel do
      do ke=1, lcmesh%Ne
        var_out(:,ke) = 0.5_RP * h_(:,ke) * ( u1_(:,ke) * U_(:,ke) + u2_(:,ke) * V_(:,ke) )
      end do
    case('ENGP')
      !$omp parallel do
      do ke=1, lcmesh%Ne
        var_out(:,ke) = Grav * ( h_(:,ke) + hs_(:,ke) )
      end do
    case('ENGT')
      !$omp parallel do
      do ke=1, lcmesh%Ne
        var_out(:,ke) = &
            0.5_RP * h_(:,ke) * ( u1_(:,ke) * U_(:,ke) + u2_(:,ke) * V_(:,ke) )  & ! ENGK       
          + Grav * ( h_(:,ke) + hs_(:,ke) )                                        ! ENGP
      end do
    end select

    return
  end subroutine vars_calc_diagnoseVar

end module mod_sw_vars