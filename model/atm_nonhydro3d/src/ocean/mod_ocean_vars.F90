!-------------------------------------------------------------------------------
!> module Ocean / Variables
!!
!! @par Description
!!          Module to manage variables with ocean component
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_ocean_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_debug

  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: &
    MeshBase2D
  use scale_mesh_base3d, only: &
    MeshBase3D,                              &
    DIMTYPE_XY  => MeshBase3D_DIMTYPEID_XY,  &
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

  use mod_ocean_mesh, only: OceanMesh
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  !> Derived type to manage variables with ocean component
  type, public :: OceanVars
    class(OceanMesh), pointer :: mesh                !< Pointer to an object to manage a primary mesh for ocean component

    !- prognostic variables
    type(MeshField3D), allocatable :: PROG_VARS(:) !< Array of 3D prognostic variables
    type(ModelVarManager) :: PROGVARS_manager      !< An object to manage prognostic variables

    !- auxiliary variables (2D)
    type(MeshField2D), allocatable :: AUX_VARS2D(:) !< Array of 2D auxiliary variables
    type(ModelVarManager) :: AUXVARS2D_manager      !< An object to manage 2D auxiliary variables

    !- atmospheric variables (2D)
    type(MeshField2D), allocatable :: ATM_VARS2D(:) !< Array of 2D atmospheric variables

    !- history file
    integer :: hist_comp_id = -1 !< Component ID for history file

    !- restart file
    type(FILE_restart_meshfield_component) :: restart_file !< Object to manage restart file for ocean component
  contains
    procedure :: Init => OceanVars_Init
    procedure :: Final => OceanVars_Final
    procedure :: History => OceanVars_History
    procedure :: Read_restart_file => OceanVars_Read_restart_file
    procedure :: Write_restart_file_prep => OceanVars_Write_restart_file_prep
    procedure :: Write_restart_file => OceanVars_Write_restart_file
    procedure :: Write_restart_file_post => OceanVars_Write_restart_file_post
  end type OceanVars

  public :: OceanVars_GetLocalMeshPrgVars

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: PRGVAR_THERM_ID   = 1 ! Variable associated with energy equation (DRHOT or ETOT)
  integer, public, parameter :: PRGVAR_DRHOT_ID   = 1
  integer, public, parameter :: PRGVAR_SALT_ID    = 2
  integer, public, parameter :: PRGVAR_W_ID       = 3
  integer, public, parameter :: PRGVAR_U_ID       = 4
  integer, public, parameter :: PRGVAR_V_ID       = 5
  integer, public, parameter :: PRGVAR_SCALAR_NUM = 3
  integer, public, parameter :: PRGVAR_HVEC_NUM   = 1
  integer, public, parameter :: PRGVAR_NUM        = 5

  integer, public, parameter :: PHYTEND_U_ID        = 1
  integer, public, parameter :: PHYTEND_V_ID        = 2
  integer, public, parameter :: PHYTEND_W_ID        = 3
  integer, public, parameter :: PHYTEND_THERM_ID    = 4
  integer, public, parameter :: PHYTEND_SALT_ID     = 5
  integer, public, parameter :: PHYTEND_NUM         = 5

  integer, public, parameter :: AUXVAR2D_SFC_TEMP_ID        = 1
  integer, public, parameter :: AUXVAR2D_SFC_ALB_IR_dir_ID  = 2
  integer, public, parameter :: AUXVAR2D_SFC_ALB_IR_dif_ID  = 3
  integer, public, parameter :: AUXVAR2D_SFC_ALB_NIR_dir_ID = 4
  integer, public, parameter :: AUXVAR2D_SFC_ALB_NIR_dif_ID = 5
  integer, public, parameter :: AUXVAR2D_SFC_ALB_VIS_dir_ID = 6
  integer, public, parameter :: AUXVAR2D_SFC_ALB_VIS_dif_ID = 7
  integer, public, parameter :: AUXVAR2D_NUM                = 7

  ! Lower atmosphere variables received from CPL buffer
  integer, public, parameter :: ATMVAR2D_SFLX_RD_SW_DIR_ID = 1
  integer, public, parameter :: ATMVAR2D_SFLX_RD_LW_DIF_ID = 2
  integer, public, parameter :: ATMVAR2D_NUM               = 2

  ! Diagnostic variables

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !

  type(VariableInfo) :: PRGVAR_VARINFO(PRGVAR_NUM)
  DATA PRGVAR_VARINFO / &
    VariableInfo( PRGVAR_THERM_ID , 'OCEAN_THERM', 'ocean temperature',     &
                  'K', 1, 'XYZ_O', 'sea_water_temperature'               ), &
    VariableInfo( PRGVAR_SALT_ID , 'OCEAN_SALT', 'ocean salinity',          &
                  'PSU', 1, 'XYZ_O', 'sea_water_salinity'                ), &
    VariableInfo( PRGVAR_W_ID , 'OCEAN_WVEL', 'ocean w-velocity',           &
                  'm/s', 1, 'XYZ_O', 'vertical_sea_water_velocity'       ), &
    VariableInfo( PRGVAR_U_ID , 'OCEAN_UVEL', 'ocean u-velocity',           &
                  'm/s', 1, 'XYZ_O', 'eastward_sea_water_velocity'       ), &
    VariableInfo( PRGVAR_V_ID , 'OCEAN_VVEL', 'ocean v-velocity',           &
                  'm/s', 1, 'XYZ_O', 'northward_sea_water_velocity'      )  /

  type(VariableInfo) :: AUXVAR2D_VARINFO(AUXVAR2D_NUM)
  DATA AUXVAR2D_VARINFO / &
    VariableInfo( AUXVAR2D_SFC_TEMP_ID , 'OCEAN_SFC_TEMP', 'ocean surface temperature', &
                  'K', 1, 'XY_O', 'sea_surface_temperature'                       ), &
    VariableInfo( AUXVAR2D_SFC_ALB_IR_dir_ID , 'OCEAN_SFC_ALB_IR_DIR', 'ocean surface albedo IR (direct)',     &
                  '1', 1, 'XY_O', 'sea_surface_albedo_ir_dir'                     ), &
    VariableInfo( AUXVAR2D_SFC_ALB_IR_dif_ID , 'OCEAN_SFC_ALB_IR_DIF', 'ocean surface albedo IR (diffuse)',    &
                  '1', 1, 'XY_O', 'sea_surface_albedo_ir_dif'                     ), &
    VariableInfo( AUXVAR2D_SFC_ALB_NIR_dir_ID , 'OCEAN_SFC_ALB_NIR_DIR', 'ocean surface albedo NIR (direct)',  &
                  '1', 1, 'XY_O', 'sea_surface_albedo_nir_dir'                     ), &
    VariableInfo( AUXVAR2D_SFC_ALB_NIR_dif_ID , 'OCEAN_SFC_ALB_NIR_DIF', 'ocean surface albedo NIR (diffuse)', &
                  '1', 1, 'XY_O', 'sea_surface_albedo_nir_dif'                     ), &
    VariableInfo( AUXVAR2D_SFC_ALB_VIS_dir_ID , 'OCEAN_SFC_ALB_VIS_DIR', 'ocean surface albedo VIS (direct)',  &
                  '1', 1, 'XY_O', 'sea_surface_albedo_vis_dir'                     ), &
    VariableInfo( AUXVAR2D_SFC_ALB_VIS_dif_ID , 'OCEAN_SFC_ALB_VIS_DIF', 'ocean surface albedo VIS (diffuse)', &
                  '1', 1, 'XY_O', 'sea_surface_albedo_vis_dif'                     )  /

contains

  !> Setup an object to manage variables with ocean component
  !!
  !! @param model_mesh Object to manage computational mesh of ocean model 
  !!
!OCL SERIAL
  subroutine OceanVars_Init( this, ocean_mesh )
    use scale_file_monitor_meshfield, only:    &
      MONITOR_reg => FILE_monitor_meshfield_reg
    use scale_mesh_base3d, only: MeshBase3D
    use scale_meshfield_base, only: MeshField3D
    implicit none

    class(OceanVars), target, intent(inout) :: this
    class(OceanMesh), target, intent(inout) :: ocean_mesh

    integer :: iv
    logical :: reg_file_hist

    character(len=H_LONG) :: IN_BASENAME           = ''        !< Basename of the input  file
    logical :: IN_POSTFIX_TIMELABEL                = .false.   !< Add timelabel to the basename of input  file?
    character(len=H_LONG) :: OUT_BASENAME          = ''        !< Basename of the output file
    logical :: OUT_POSTFIX_TIMELABEL               = .true.    !< Add timelabel to the basename of output file?
    character(len=H_MID) :: OUT_TITLE              = ''        !< Title    of the output file
    character(len=H_SHORT) :: OUT_DTYPE            = 'DEFAULT' !< REAL4 or REAL8  

    namelist / PARAM_OCEAN_VARS_RESTART / &
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

    !---------------------------------------------------------

    this%mesh => ocean_mesh
    mesh3D => ocean_mesh%ptr_mesh
    call mesh3D%GetMesh2D( mesh2D )

    !-
    call this%PROGVARS_manager%Init()
    call this%AUXVARS2D_manager%Init()

    allocate( this%PROG_VARS(PRGVAR_NUM) )
    allocate( this%AUX_VARS2D(AUXVAR2D_NUM) )
    allocate( this%ATM_VARS2D(ATMVAR2D_NUM) )
    !$acc enter data create( this%PROG_VARS, this%AUX_VARS2D, this%ATM_VARS2D )

    !- Initialize prognostic variables

    reg_file_hist = .true.
    do iv = 1, PRGVAR_NUM
      call this%PROGVARS_manager%Regist(  &
        PRGVAR_VARINFO(iv), mesh3D,                           & ! (in) 
        this%PROG_VARS(iv),                                   & ! (inout)
        reg_file_hist,  monitor_flag=.true., fill_zero=.true. ) ! (out)
      !$acc update device( this%PROG_VARS(iv) )
    end do

    !- Initialize auxiliary variables (2D)

    reg_file_hist = .true.
    do iv = 1, AUXVAR2D_NUM
      call this%AUXVARS2D_manager%Regist( &
        AUXVAR2D_VARINFO(iv), mesh2D,         & ! (in) 
        this%AUX_VARS2D(iv),                  & ! (inout)
        reg_file_hist, fill_zero=.true.       ) ! (in)
      !$acc update device( this%AUX_VARS2D(iv) )
    end do

    !- Initialize atmospheric variables (2D)

    do iv =1, ATMVAR2D_NUM
      call this%ATM_VARS2D(iv)%Init( "", "", mesh2D )
      !$acc update device( this%ATM_VARS2D(iv) )
    end do

    !-- Setup information for input/output restart files. 

    is_specified = .true.
    !- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_VARS_RESTART,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OceanVars_Init",*) 'Not found namelist PARAM_OCEAN_VARS_RESTART. Default used.'
       is_specified = .false.
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OceanVars_Init",*) 'Not appropriate names in namelist PARAM_OCEAN_VARS_RESTART. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_VARS_RESTART)

    if (is_specified) then
      call ocean_mesh%Setup_restartfile( this%restart_file,         &
        IN_BASENAME, IN_POSTFIX_TIMELABEL,                          &
        OUT_BASENAME, OUT_POSTFIX_TIMELABEL, OUT_DTYPE, OUT_TITLE,  &
        PRGVAR_NUM + AUXVAR2D_NUM, "_O" )
    else
      call ocean_mesh%Setup_restartfile( this%restart_file, &
        PRGVAR_NUM + AUXVAR2D_NUM, "_O" )
    end if

    return
  end subroutine OceanVars_Init

  !> Finalize an object to manage variables with ocean component
  !!
!OCL SERIAL
  subroutine OceanVars_Final( this )
    implicit none
    class(OceanVars), intent(inout) :: this

    integer :: iv
    !--------------------------------------------------

    LOG_INFO('OceanVars_Final',*)

    !$acc exit data delete( this%PROG_VARS, this%AUX_VARS2D, this%ATM_VARS2D )

    call this%PROGVARS_manager%Final()
    deallocate( this%PROG_VARS )

    call this%AUXVARS2D_manager%Final()
    deallocate( this%AUX_VARS2D )

    do iv=1, ATMVAR2D_NUM
      call this%ATM_VARS2D(iv)%Final()
    end do
    deallocate( this%ATM_VARS2D )

    call this%restart_file%Final()
    return
  end subroutine OceanVars_Final

  !> Put data with oceanic variables to history file
  !!
!OCL SERIAL
  subroutine OceanVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(OceanVars), intent(inout), target :: this
  
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
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%PROG_VARS(v), this%hist_comp_id )
    end do
    do v = 1, AUXVAR2D_NUM
      hst_id = this%AUX_VARS2D(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%AUX_VARS2D(v), this%hist_comp_id )
    end do
    
    return
  end subroutine OceanVars_history

!> Read data with oceanic variables from restart file
!!
!OCL SERIAL
  subroutine OceanVars_Read_restart_file( this, ocean_mesh )
    use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
    use scale_meshfieldcomm_base, only: MeshFieldContainer
    implicit none
    
    class(OceanVars), intent(inout), target :: this
    class(OceanMesh), intent(in) :: ocean_mesh

    integer :: iv
    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lcmesh3D
    !---------------------------------------

    LOG_NEWLINE
    LOG_INFO("OceanVar_Read_restart_file",*) 'Open restart file (OCEAN) '
        
    !- Open restart file
    call this%restart_file%Open()

    !- Read restart file
    
    do iv=1, PRGVAR_NUM
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%PROG_VARS(iv)%varname, &
        this%PROG_VARS(iv)                                                      )
    end do
    do iv=1, AUXVAR2D_NUM
      call this%restart_file%Read_var( DIMTYPE_XY, this%AUX_VARS2D(iv)%varname, &
        this%AUX_VARS2D(iv)                                                      )
    end do

    !- Close restart file
    LOG_INFO("OceanVar_Read_restart_file",*) 'Close restart file (OCEAN) '
    call this%restart_file%Close()

    !-- Prepare diagnostic variables

    !-- Check read data
    ! call this%Check( force = .true. )

    !-- Calculate diagnostic variables
    ! call this%Calc_diagnostics()   

    !-- Communicate halo data of diagnostic variables
    ! call this%AUXVARS_manager%MeshFieldComm_Exchange()

    return
  end subroutine OceanVars_Read_restart_file

  !> Write data with oceanic variables to restart file
  !!
!OCL SERIAL
  subroutine OceanVars_Write_restart_file_prep( this )
    ! use scale_ocn_dyn_dgm_nonhydro3d_common, only: &
    !   ocn_dyn_dgm_nonhydro3d_common_get_varinfo
    
    implicit none
    class(OceanVars), intent(inout) :: this

    ! type(VariableInfo) :: prgvar_info(PRGVAR_NUM)
    ! type(VariableInfo) :: auxvar2D_info(AUXVAR2D_NUM)

    integer :: iv, rf_vid 
    !---------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("OceanVars_Write_restart_file",*) 'Create restart file (OCEAN) '

    !- Check data which will be written to restart file
    ! call this%Check( force = .true. )

    !- Create restart file
    call this%restart_file%Create()
    call PRC_mpibarrier()

    !- Define variables

    ! call ocn_dyn_dgm_nonhydro3d_common_get_varinfo( prgvar_info, auxvar2D_info )

    do iv=1, PRGVAR_NUM
      rf_vid = iv
      call this%restart_file%Def_var( this%PROG_VARS(iv),  &
        PRGVAR_VARINFO(iv)%DESC, rf_vid, DIMTYPE_XYZ       )
    end do
    do iv=1, AUXVAR2D_NUM
      rf_vid = PRGVAR_NUM + iv
      call this%restart_file%Def_var( this%AUX_VARS2D(iv),   &
        AUXVAR2D_VARINFO(iv)%DESC, rf_vid, DIMTYPE_XY        )
    end do

    return
  end subroutine OceanVars_Write_restart_file_prep

  !> Write data with oceanic variables to restart file
  !!
!OCL SERIAL
  subroutine OceanVars_Write_restart_file( this )
    ! use scale_ocn_dyn_dgm_nonhydro3d_common, only: &
    !   ocn_dyn_dgm_nonhydro3d_common_get_varinfo
    
    implicit none
    class(OceanVars), intent(inout) :: this

    ! type(VariableInfo) :: prgvar_info(PRGVAR_NUM)
    ! type(VariableInfo) :: auxvar2D_info(AUXVAR2D_NUM)

    integer :: iv, rf_vid 
    !---------------------------------------

    call this%restart_file%End_def()

    !- Write restart file
    do iv=1, PRGVAR_NUM
      rf_vid = iv
      call this%restart_file%Write_var( rf_vid, this%PROG_VARS(iv) )
    end do
    do iv=1, AUXVAR2D_NUM
      rf_vid = PRGVAR_NUM + iv
      call this%restart_file%Write_var( rf_vid, this%AUX_VARS2D(iv) )
    end do
    return
  end subroutine OceanVars_Write_restart_file

  !> Close restart file for oceanic variables
!OCL SERIAL
  subroutine OceanVars_Write_restart_file_post( this )
    implicit none
    class(OceanVars), intent(inout) :: this
    !---------------------------------------
    LOG_INFO("OceanVars_Write_restart_file",*) 'Close restart file (OCEAN) '
    call this%restart_file%Close()
    return
  end subroutine OceanVars_Write_restart_file_post

!OCL SERIAL
!OCL SERIAL
  subroutine OceanVars_GetLocalMeshPrgVars( domID, mesh, prgvars_list, auxvars2D_list, &
    U, V, W, THERM, SALT,                                                            &
    SFC_TEMP, &
    SFC_ALB_IR_dir, SFC_ALB_IR_dif, SFC_ALB_NIR_dir, SFC_ALB_NIR_dif,                &
    SFC_ALB_VIS_dir, SFC_ALB_VIS_dif,                                                &
    lcmesh3D                                                                         )
    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars2D_list
    class(LocalMeshFieldBase), pointer, intent(out) :: U, V, W, THERM, SALT
    class(LocalMeshFieldBase), pointer, intent(out) :: SFC_TEMP
    class(LocalMeshFieldBase), pointer, intent(out), optional :: SFC_ALB_IR_dir, SFC_ALB_IR_dif
    class(LocalMeshFieldBase), pointer, intent(out), optional :: SFC_ALB_NIR_dir, SFC_ALB_NIR_dif
    class(LocalMeshFieldBase), pointer, intent(out), optional :: SFC_ALB_VIS_dir, SFC_ALB_VIS_dif
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call prgvars_list%Get(PRGVAR_U_ID, field)
    call field%GetLocalMeshField(domID, U)

    call prgvars_list%Get(PRGVAR_V_ID, field)
    call field%GetLocalMeshField(domID, V)

    call prgvars_list%Get(PRGVAR_W_ID, field)
    call field%GetLocalMeshField(domID, W)

    call prgvars_list%Get(PRGVAR_THERM_ID, field)
    call field%GetLocalMeshField(domID, THERM)
  
    call prgvars_list%Get(PRGVAR_SALT_ID, field)
    call field%GetLocalMeshField(domID, SALT)

    call auxvars2D_list%Get(AUXVAR2D_SFC_TEMP_ID, field)
    call field%GetLocalMeshField(domID, SFC_TEMP)

    if ( present(SFC_ALB_IR_dir) ) then
      call auxvars2D_list%Get(AUXVAR2D_SFC_ALB_IR_dir_ID, field)
      call field%GetLocalMeshField(domID, SFC_ALB_IR_dir)
    end if
    if ( present(SFC_ALB_IR_dif) ) then
      call auxvars2D_list%Get(AUXVAR2D_SFC_ALB_IR_dif_ID, field)
      call field%GetLocalMeshField(domID, SFC_ALB_IR_dif)
    end if
    if ( present(SFC_ALB_NIR_dir) ) then
      call auxvars2D_list%Get(AUXVAR2D_SFC_ALB_NIR_dir_ID, field)
      call field%GetLocalMeshField(domID, SFC_ALB_NIR_dir)
    end if
    if ( present(SFC_ALB_NIR_dif) ) then
      call auxvars2D_list%Get(AUXVAR2D_SFC_ALB_NIR_dif_ID, field)
      call field%GetLocalMeshField(domID, SFC_ALB_NIR_dif)
    end if
    if ( present(SFC_ALB_VIS_dir) ) then
      call auxvars2D_list%Get(AUXVAR2D_SFC_ALB_VIS_dir_ID, field)
      call field%GetLocalMeshField(domID, SFC_ALB_VIS_dir)
    end if
    if ( present(SFC_ALB_VIS_dif) ) then
      call auxvars2D_list%Get(AUXVAR2D_SFC_ALB_VIS_dif_ID, field)
      call field%GetLocalMeshField(domID, SFC_ALB_VIS_dif)
    end if  

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
  end subroutine OceanVars_GetLocalMeshPrgVars
end module mod_ocean_vars