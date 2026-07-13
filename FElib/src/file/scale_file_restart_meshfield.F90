!-------------------------------------------------------------------------------
!> module FElib / File / Restart
!!
!! @par Description
!!           A module for outputting data to restart simulations
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_file_restart_meshfield
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: &
    PRC_myrank, &
    PRC_abort
  use scale_file_h, only: &
     FILE_FILE_MAX
  use scale_file_common_meshfield, only: &
    FILE_common_meshfield_diminfo
  
  use scale_element_base, only: ElementBase1D, ElementBase2D, ElementBase3D
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D

  use scale_meshfield_base, only: MeshFieldBase, MeshField1D, MeshField2D, MeshField3D
  use scale_variableinfo, only: VariableInfo

  use scale_file_base_meshfield, only: &
    FILE_base_meshfield

  use scale_mesh_base1d, only: &
    MF1D_DIMTYPE_X => MeshBase1D_DIMTYPEID_X, &
    MF1D_DTYPE_NUM => MeshBase1D_DIMTYPE_NUM 
  use scale_mesh_base2d, only: &
    MF2D_DIMTYPE_X => MeshBase2D_DIMTYPEID_X, &
    MF2D_DIMTYPE_Y => MeshBase2D_DIMTYPEID_Y, &
    MF2D_DTYPE_NUM => MeshBase2D_DIMTYPE_NUM 
  use scale_mesh_base3d, only: &
    MF3D_DIMTYPE_X => MeshBase3D_DIMTYPEID_X, &
    MF3D_DIMTYPE_Y => MeshBase3D_DIMTYPEID_Y, &
    MF3D_DIMTYPE_Z => MeshBase3D_DIMTYPEID_Z, &
    MF3D_DTYPE_NUM => MeshBase3D_DIMTYPE_NUM    
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  !> Base derived type to manage restart file with each component
  type :: FILE_restart_meshfield
    type(FILE_base_meshfield), pointer :: file_ptr => null() !< Pointer to the file object for managing restart file with each component
    logical :: is_file_ptr_main                              !< Flag whether the file_ptr is the main file object

    logical :: flag_output                 !< Flag whether output restart file
    character(len=H_LONG) :: in_basename   !< Basename of the input file
    logical :: in_postfix_timelabel        !< Flag whether timelabel is added to the basename of input file

    character(len=H_LONG) :: out_basename  !< Basename of the output file
    logical :: out_postfix_timelabel       !< Flag whether timelabel is added to the basename of output file
    character(len=H_MID) :: out_title      !< Title of the output file
    character(len=H_SHORT) :: out_dtype    !< Output data type (DEFAULT, REAL4 or REAL8)
  contains
    procedure, private :: Set_information => FILE_restart_meshfield_component_set_information
  end type FILE_restart_meshfield
  type(FILE_restart_meshfield), public :: restart_file    !< An object to manage main restart file

  !> Derived type to manage restart file with each component
  type, extends(FILE_restart_meshfield), public :: FILE_restart_meshfield_component
    character(len=H_SHORT) :: comp_name  !< Name of the component
    integer :: registered_comp_id        !< ID of the component registered for restart file
  contains
    procedure :: Init1 => FILE_restart_meshfield_component_Init1
    procedure :: Init2 => FILE_restart_meshfield_component_Init2
    generic :: Init => Init1, Init2
    procedure :: Open => FILE_restart_meshfield_component_open    
    procedure :: FILE_restart_meshfield_component_create
    generic :: Create => FILE_restart_meshfield_component_create
    procedure :: FILE_restart_meshfield_component_def_var
    generic :: Def_var => FILE_restart_meshfield_component_def_var
    procedure :: End_def => FILE_restart_meshfield_component_enddef
    procedure :: FILE_restart_meshfield_component_write_var2d
    procedure :: FILE_restart_meshfield_component_write_var3d
    generic :: Write_var => &
      FILE_restart_meshfield_component_write_var2d, &
      FILE_restart_meshfield_component_write_var3d
    procedure :: Close => FILE_restart_meshfield_component_close
    procedure :: FILE_restart_meshfield_component_read_var2d
    procedure :: FILE_restart_meshfield_component_read_var3d
    generic :: Read_Var => &
      FILE_restart_meshfield_component_read_var2d, &
      FILE_restart_meshfield_component_read_var3d
    procedure :: Final => FILE_restart_meshfield_component_Final
  end type

  ! Public procedures for main restart file
  public :: FILE_restart_meshfield_setup
  public :: FILE_restart_meshfield_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  type(FILE_base_meshfield), target :: file_main  !< A file object to manage main restart file

  logical :: is_file_main_opened          !< Flag whether the main restart file is opened
  logical :: is_file_main_created         !< Flag whether the main restart file is created
  logical :: is_file_main_enddef          !< Flag whether the main restart file is in enddef mode
  logical :: is_file_main_closed          !< Flag whether the main restart file is closed

contains

!----------------
  !> Setup the main restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_setup()
    implicit none

    logical :: OUTPUT_FLAG                         = .false.   !< Flag whether output restart file
    character(len=H_LONG) :: IN_BASENAME           = ''        !< Basename of the input  file
    logical :: IN_POSTFIX_TIMELABEL                = .false.   !< Flag whether timelabel is added to the basename of input  file
    character(len=H_LONG) :: OUT_BASENAME          = ''        !< Basename of the output file
    logical :: OUT_POSTFIX_TIMELABEL               = .true.    !< Flag whether timelabel is added to the basename of output file
    character(len=H_MID) :: OUT_TITLE              = ''        !< Title of the output file
    character(len=H_SHORT) :: OUT_DTYPE            = 'DEFAULT' !< Output data type (DEFAULT, REAL4 or REAL8)

    namelist / PARAM_RESTART / &
      OUTPUT_FLAG,           &
      IN_BASENAME,           &
      IN_POSTFIX_TIMELABEL,  &
      OUT_BASENAME,          &
      OUT_POSTFIX_TIMELABEL, &
      OUT_TITLE,             &
      OUT_DTYPE

    integer :: ierr  
    !----------------------------------------

    LOG_NEWLINE

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_RESTART,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("FILE_restart_meshfield_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("FILE_restart_meshfield_setup",*) 'Not appropriate names in namelist PARAM_RESTART. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_RESTART)

    ! Setup the information of the main restart file

    call restart_file%Set_information( IN_BASENAME, IN_POSTFIX_TIMELABEL, &
      OUT_BASENAME, OUT_POSTFIX_TIMELABEL, OUT_DTYPE, OUT_TITLE           )

    restart_file%flag_output = OUTPUT_FLAG

    restart_file%file_ptr => null() ! Initialize the pointer to null
    is_file_main_opened  = .false.
    is_file_main_created = .false.
    is_file_main_enddef  = .false.
    is_file_main_closed  = .false.

    return
  end subroutine FILE_restart_meshfield_setup

  !> Finalize the main restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_finalize()
    implicit none
    !--------------------------------------------------
    if ( associated(restart_file%file_ptr) ) then
      call restart_file%file_ptr%Final()
      nullify( restart_file%file_ptr )
    end if
    return
  end subroutine FILE_restart_meshfield_finalize

  !> Initialize an object to manage restart file with each component, and register the component to the main restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_Init1( this, &
    comp_name, var_num, dim_name_postfix,                  &
    mesh1D, mesh2D, meshCubedSphere2D,                     &
    mesh3D, meshCubedSphere3D                              )

    use scale_file_common_meshfield, only: &
      File_common_meshfield_get_dims  
    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this
    character(*), intent(in) :: comp_name
    integer, intent(in) :: var_num
    character(len=H_SHORT), optional, intent(in) :: dim_name_postfix
    class(MeshBase1D), target, optional, intent(in) :: mesh1D
    class(MeshRectDom2D), target, optional, intent(in) :: mesh2D
    class(MeshCubedSphereDom2D), target, optional, intent(in) :: meshCubedSphere2D
    class(MeshCubeDom3D), target, optional, intent(in) :: mesh3D
    class(MeshCubedSphereDom3D), target, optional, intent(in) :: meshCubedSphere3D

    character(len=H_SHORT) :: dim_name_postfix_
    !--------------------------------------------------

    ! Setup an object to manage restart file with each component, and register the component to the main restart file

    if (present(dim_name_postfix)) then
      dim_name_postfix_ = trim(dim_name_postfix)
    else
      dim_name_postfix_ = ""
    end if
    
    if ( .not. associated(restart_file%file_ptr) ) then
      ! Initialize the main restart file
      call file_main%Init( var_num, comp_name, dim_name_postfix_,    &
        mesh1D, mesh2D, meshCubedSphere2D, mesh3D, meshCubedSphere3D )
        
      restart_file%file_ptr => file_main
      ! The first registered component is the main component, and its ID is 1
      this%registered_comp_id = 1
    else
      ! Register the subsequent component to the main restart file
      call file_main%Register_comp( this%registered_comp_id,         & ! (out)
        comp_name, var_num, dim_name_postfix_,                       & ! (in)
        mesh1D, mesh2D, meshCubedSphere2D, mesh3D, meshCubedSphere3D ) ! (in)
    end if

    this%file_ptr => file_main
    this%is_file_ptr_main = .true.
    
    ! Setup the information of the component restart file based on the main restart file
    call this%Set_information( &
      restart_file%in_basename, restart_file%in_postfix_timelabel,   &
      restart_file%out_basename, restart_file%out_postfix_timelabel, &
      restart_file%out_dtype, restart_file%out_title                 )

    return
  end subroutine FILE_restart_meshfield_component_Init1

  !> Initialize an object to manage restart file with each component. 
  !! Restart data for this component are individually outputted to a separate file, not to the main restart file.
  !!
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_Init2( this,  &
      comp_name,                                            &
      in_basename, in_postfix_timelabel,                    &
      out_basename, out_postfix_timelabel,                  &
      out_dtype, out_title,                                 &
      var_num, dim_name_postfix,                            &
      mesh1D,                                               &
      mesh2D, meshCubedSphere2D,                            &
      mesh3D, meshCubedSphere3D                             )

    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this
    character(*), intent(in) :: comp_name
    character(*), intent(in) :: in_basename
    logical, intent(in) :: in_postfix_timelabel
    character(*), intent(in) :: out_basename
    logical, intent(in) :: out_postfix_timelabel
    character(*), intent(in) :: out_title
    character(*), intent(in) :: out_dtype
    integer, intent(in) :: var_num
    character(len=H_SHORT), intent(in) :: dim_name_postfix
    class(MeshBase1D), target, optional, intent(in) :: mesh1D
    class(MeshRectDom2D), target, optional, intent(in) :: mesh2D
    class(MeshCubedSphereDom2D), target, optional, intent(in) :: meshCubedSphere2D
    class(MeshCubeDom3D), target, optional, intent(in) :: mesh3D
    class(MeshCubedSphereDom3D), target, optional, intent(in) :: meshCubedSphere3D
    !--------------------------------------------------

    call this%Set_information( in_basename, in_postfix_timelabel, &
      out_basename, out_postfix_timelabel, out_dtype, out_title  )
    
    if ( .not. associated(this%file_ptr) ) then
      ! In this case, restart data for this component is individually outputted to a separate file, not to the main restart file. 
      ! Therefore, a new file object is created for this component.
      allocate( this%file_ptr )
      call this%file_ptr%Init( var_num, comp_name, dim_name_postfix, &
        mesh1D, mesh2D, meshCubedSphere2D, mesh3D, meshCubedSphere3D )
      
      this%registered_comp_id = 1
      this%is_file_ptr_main = .false.
    else
      LOG_INFO("FILE_restart_meshfield_component_Init2",*) 'File object for this component is already associated, which is unexpected. Check!'
      call PRC_abort
    end if

    return
  end subroutine FILE_restart_meshfield_component_Init2

  !> Open a restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_open( &
    this )

    use scale_time, only: &
      TIME_gettimelabel    
    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename  
    !--------------------------------------------------------------

    if ( this%in_basename == '' ) then
      LOG_INFO(trim(this%comp_name)//"_vars_restart_open",*) 'restart file is not specified. Check!'
      call PRC_abort
      return
    end if

    if ( this%is_file_ptr_main ) then
      if ( is_file_main_opened ) then
        LOG_INFO(trim(this%comp_name)//"_vars_restart_open",*) 'Main restart file is already opened. Skipping opening.'
      else
        is_file_main_opened = .true.
      end if
    end if

    !--------------------------------

    if ( this%in_postfix_timelabel ) then
      call TIME_gettimelabel( timelabel )
      basename = trim(this%in_basename)//'_'//trim(timelabel)
    else
      basename = trim(this%in_basename)
    endif

    !--------------------------------

    LOG_NEWLINE
    LOG_INFO(trim(this%comp_name)//"_vars_restart_open",*) 'Open restart file'
    call this%file_ptr%open( basename, myrank=PRC_myrank )

    return
  end subroutine FILE_restart_meshfield_component_open

  !> Create a restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_create( &
      this )

    use scale_time, only: &
      TIME_gettimelabel,         &
      NOWDATE   => TIME_NOWDATE, &
      NOWSUBSEC => TIME_NOWSUBSEC
    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this

    character(len=19)     :: timelabel
    character(len=H_LONG) :: basename  
    character(len=34) :: tunits
    character(len=H_SHORT) :: calendar
    logical :: fileexisted
    !--------------------------------------------------------------

    if ( this%out_basename == '' ) return

    if ( this%is_file_ptr_main ) then
      if ( is_file_main_created ) then
        LOG_INFO(trim(this%comp_name)//"_vars_restart_create",*) 'Main restart file is already created. Skipping creation.'
      else
        is_file_main_created = .true.
      end if
    end if

    !--------------------------------

    LOG_NEWLINE
    LOG_INFO(trim(this%comp_name)//"_vars_restart_create",*) 'Create restart file'

    if ( this%out_postfix_timelabel ) then
        call TIME_gettimelabel( timelabel )
        basename = trim(this%out_basename)//'_'//trim(timelabel)
    else
        basename = trim(this%out_basename)
    endif

    LOG_INFO(trim(this%comp_name)//"_vars_restart_create",*) 'basename: ', trim(basename)

    
    call this%file_ptr%Create( basename, this%out_title, this%out_dtype,       & ! (in)
                           fileexisted,                                        & ! (out)
                           myrank=PRC_myrank, tunits=tunits, calendar=calendar ) ! (in)
    
    if ( .not. fileexisted ) then
      call this%file_ptr%Put_GlobalAttribute_time( NOWDATE, NOWSUBSEC )
    end if

    return
  end subroutine FILE_restart_meshfield_component_create

  !> Define a variable in the restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_def_var( this,      &
    field, desc, vid, dim_type_id                                 )
    implicit none
  
    class(FILE_restart_meshfield_component), intent(inout) :: this
    class(MeshFieldBase), intent(in) :: field
    character(len=*), intent(in) :: desc
    integer, intent(in) :: vid
    integer, intent(in) :: dim_type_id
    !------------------------------------------------------------------

    call this%file_ptr%Def_var( &
      field, desc, vid, dim_type_id, this%out_dtype, &
      comp_id=this%registered_comp_id                )   
     
    return 
  end subroutine FILE_restart_meshfield_component_def_var

  !> Finalize the definition of variables in the restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_enddef( this )
    implicit none
    class(FILE_restart_meshfield_component), intent(inout) :: this
    !--------------------------------------------------------------

    if ( this%is_file_ptr_main ) then
      if ( is_file_main_enddef ) then
        LOG_INFO(trim(this%comp_name)//"_vars_restart_close",*) 'Main restart file is already in enddef mode. Skipping enddef.'
      else
        is_file_main_enddef = .true.
      end if
    end if
      
    call this%file_ptr%End_def()
    return
  end subroutine FILE_restart_meshfield_component_enddef

  !> Write a 2D variable to the restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_write_var2d( this, &
    vid, field2d )
    
    use scale_time, only: TIME_NOWDAYSEC
    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this
    integer, intent(in) :: vid
    class(MeshField2D), intent(in) :: field2d
    !--------------------------------------------------

    call this%file_ptr%Write_var2D( vid, field2d, TIME_NOWDAYSEC, TIME_NOWDAYSEC, &
      comp_id=this%registered_comp_id )
    
    return
  end subroutine FILE_restart_meshfield_component_write_var2d

  !> Write a 3D variable to the restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_write_var3d( this, &
    vid, field3d )
    
    use scale_time, only: TIME_NOWDAYSEC
    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this
    integer, intent(in) :: vid
    class(MeshField3D), intent(in) :: field3d
    !--------------------------------------------------

    call this%file_ptr%Write_var3D( vid, field3d, TIME_NOWDAYSEC, TIME_NOWDAYSEC, &
      comp_id=this%registered_comp_id )

    return
  end subroutine FILE_restart_meshfield_component_write_var3d

  !> Read a 2D variable from the restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_read_var2d( this,  &
    dim_typeid, varname, field2d, step, allow_missing  )
  
  
    implicit none
  
    class(FILE_restart_meshfield_component), intent(inout) :: this
    integer, intent(in) :: dim_typeid
    character(*), intent(in) :: varname
    class(MeshField2D), intent(inout) :: field2d
    integer, intent(in), optional :: step
    logical, intent(in), optional :: allow_missing
    !------------------------------------------------------

    call this%file_ptr%Read_Var( &
          dim_typeid, varname, field2d, step, allow_missing, &
          comp_id=this%registered_comp_id                    )
    
    return
  end subroutine FILE_restart_meshfield_component_read_var2d

  !> Read a 3D variable from the restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_read_var3d( this,  &
    dim_typeid, varname, field3d, step, allow_missing  )
  
  
    implicit none
  
    class(FILE_restart_meshfield_component), intent(inout) :: this
    integer, intent(in) :: dim_typeid
    character(*), intent(in) :: varname
    class(MeshField3D), intent(inout) :: field3d
    integer, intent(in), optional :: step
    logical, intent(in), optional :: allow_missing
    !------------------------------------------------------

    call this%file_ptr%Read_Var( &
          dim_typeid, varname, field3d, step, allow_missing, &
          comp_id=this%registered_comp_id                    )
    
    return
  end subroutine FILE_restart_meshfield_component_read_var3d

  !> Close a restart file
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_close( this )
    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this
    !--------------------------------------------------

    if ( this%is_file_ptr_main ) then
      if ( is_file_main_closed ) then
        LOG_INFO(trim(this%comp_name)//"_vars_restart_close",*) 'Main restart file is already closed. Skipping closure.'
      else
        is_file_main_closed = .true.
      end if
    end if
        
    if ( this%file_ptr%fid /= -1 ) then
      LOG_NEWLINE
      LOG_INFO(trim(this%comp_name)//"_vars_restart_close",*) 'Close restart file'
      call this%file_ptr%Close()
    end if

    return
  end subroutine FILE_restart_meshfield_component_close

  !> Finalize an object to manage restart file with each component
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_Final( this )
    implicit none
    class(FILE_restart_meshfield_component), intent(inout) :: this
    !--------------------------------------------------

    ! Release the file object for managing restart file with each component

    if ( .not. this%is_file_ptr_main ) then
      call this%file_ptr%Final()
      deallocate( this%file_ptr )
    end if

    nullify( this%file_ptr )
    return
  end subroutine FILE_restart_meshfield_component_Final

  !- Private -----------

  !> Set the information of an object to manage restart file with each component
!OCL SERIAL
  subroutine FILE_restart_meshfield_component_set_information( this,  &
      in_basename, in_postfix_timelabel,                    &
      out_basename, out_postfix_timelabel,                  &
      out_dtype, out_title                                  )
    implicit none

    class(FILE_restart_meshfield), intent(inout) :: this
    character(*), intent(in) :: in_basename
    logical, intent(in) :: in_postfix_timelabel
    character(*), intent(in) :: out_basename
    logical, intent(in) :: out_postfix_timelabel
    character(*), intent(in) :: out_title
    character(*), intent(in) :: out_dtype
    !--------------------------------------------------

    this%in_basename = in_basename
    this%in_postfix_timelabel = in_postfix_timelabel

    this%out_basename = out_basename
    this%out_postfix_timelabel = out_postfix_timelabel
    this%out_title = out_title
    this%out_dtype = out_dtype
    return
  end subroutine FILE_restart_meshfield_component_set_information

end module scale_file_restart_meshfield
