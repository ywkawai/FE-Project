!> module file_restart_meshfield
!!
!! @par Description
!!           module for outputing data to restart simulations
!!
!! @author Team SCALE
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
  
  use scale_element_base, only: elementbase1D, elementbase2D, elementbase3D
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
  type :: FILE_restart_meshfield
    type(FILE_base_meshfield) :: base

    logical :: flag_output
    character(len=H_LONG) :: in_basename
    logical :: in_postfix_timelabel

    character(len=H_LONG) :: out_basename
    logical :: out_postfix_timelabel
    character(len=H_MID) :: out_title
    character(len=H_SHORT) :: out_dtype
  end type FILE_restart_meshfield
  type(FILE_restart_meshfield), public :: restart_file

  public :: FILE_restart_meshfield_setup


  type, extends(FILE_restart_meshfield), public :: FILE_restart_meshfield_component
    character(len=H_SHORT) :: comp_name
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

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains

!----------------

  subroutine FILE_restart_meshfield_setup()
    implicit none

    logical :: OUTPUT_FLAG                         = .false.   !< Output restart file?
    character(len=H_LONG) :: IN_BASENAME           = ''        !< Basename of the input  file
    logical :: IN_POSTFIX_TIMELABEL                = .false.   !< Add timelabel to the basename of input  file?
    character(len=H_LONG) :: OUT_BASENAME          = ''        !< Basename of the output file
    logical :: OUT_POSTFIX_TIMELABEL               = .true.    !< Add timelabel to the basename of output file?
    character(len=H_MID) :: OUT_TITLE              = ''        !< Title    of the output file
    character(len=H_SHORT) :: OUT_DTYPE            = 'DEFAULT' !< REAL4 or REAL8

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

    restart_file%flag_output = OUTPUT_FLAG

    restart_file%in_basename = IN_BASENAME
    restart_file%in_postfix_timelabel = IN_POSTFIX_TIMELABEL
    
    restart_file%out_basename = OUT_BASENAME
    restart_file%out_postfix_timelabel = OUT_POSTFIX_TIMELABEL
    restart_file%out_title = OUT_TITLE
    restart_file%out_dtype = OUT_DTYPE

    restart_file%base%fid = -1

    return
  end subroutine FILE_restart_meshfield_setup

  subroutine FILE_restart_meshfield_component_Init1( this,  &
    comp_name,                                              &
    var_num, mesh1D, mesh2D, meshCubedSphere2D,             &
    mesh3D, meshCubedSphere3D                               )

  use scale_file_common_meshfield, only: &
    File_common_meshfield_get_dims  

    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this
    character(*), intent(in) :: comp_name
    integer, intent(in) :: var_num
    class(MeshBase1D), target, optional, intent(in) :: mesh1D
    class(MeshRectDom2D), target, optional, intent(in) :: mesh2D
    class(MeshCubedSphereDom2D), target, optional, intent(in) :: meshCubedSphere2D
    class(MeshCubeDom3D), target, optional, intent(in) :: mesh3D
    class(MeshCubedSphereDom3D), target, optional, intent(in) :: meshCubedSphere3D

    !--------------------------------------------------

    call this%Init2( &
      comp_name,                                                            &
      restart_file%in_basename, restart_file%in_postfix_timelabel,          &
      restart_file%out_basename, restart_file%out_postfix_timelabel,        &
      restart_file%out_dtype, restart_file%out_title,                       &
      var_num, mesh1D, mesh2D, meshCubedSphere2D, mesh3D, meshCubedSphere3D )

    return
  end subroutine FILE_restart_meshfield_component_Init1

  subroutine FILE_restart_meshfield_component_Init2( this,  &
      comp_name,                                            &
      in_basename, in_postfix_timelabel,                    &
      out_basename, out_postfix_timelabel,                  &
      out_dtype, out_title,                                 &
      var_num,                                              &
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
    class(MeshBase1D), target, optional, intent(in) :: mesh1D
    class(MeshRectDom2D), target, optional, intent(in) :: mesh2D
    class(MeshCubedSphereDom2D), target, optional, intent(in) :: meshCubedSphere2D
    class(MeshCubeDom3D), target, optional, intent(in) :: mesh3D
    class(MeshCubedSphereDom3D), target, optional, intent(in) :: meshCubedSphere3D
    !--------------------------------------------------
    this%comp_name = comp_name

    this%in_basename = in_basename
    this%in_postfix_timelabel = in_postfix_timelabel

    this%out_basename = out_basename
    this%out_postfix_timelabel = out_postfix_timelabel
    this%out_title = out_title
    this%out_dtype = out_dtype

    !-
    call this%base%Init( var_num, mesh1D, mesh2D, meshCubedSphere2D, mesh3D, meshCubedSphere3D )

    return
  end subroutine FILE_restart_meshfield_component_Init2

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

    if ( this%in_postfix_timelabel ) then
      call TIME_gettimelabel( timelabel )
      basename = trim(this%in_basename)//'_'//trim(timelabel)
    else
      basename = trim(this%in_basename)
    endif

    !--------------------------------

    LOG_NEWLINE
    LOG_INFO(trim(this%comp_name)//"_vars_restart_open",*) 'Open restart file'
    call this%base%open( basename, myrank=PRC_myrank )

    return
  end subroutine FILE_restart_meshfield_component_open

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

    
    call this%base%Create( basename, this%out_title, this%out_dtype,           & ! (in)
                           fileexisted,                                        & ! (out)
                           myrank=PRC_myrank, tunits=tunits, calendar=calendar ) ! (in)
    
    if ( .not. fileexisted ) then
      call this%base%Put_GlobalAttribute_time( NOWDATE, NOWSUBSEC )
    end if

    return
  end subroutine FILE_restart_meshfield_component_create

  subroutine FILE_restart_meshfield_component_def_var( this,      &
    field, desc, vid, dim_type_id                                 )
  
    implicit none
  
    class(FILE_restart_meshfield_component), intent(inout) :: this
    class(MeshFieldBase), intent(in) :: field
    character(len=*), intent(in) :: desc
    integer, intent(in) :: vid
    integer, intent(in) :: dim_type_id
    !------------------------------------------------------------------

    call this%base%Def_var( &
      field, desc, vid, dim_type_id, this%out_dtype )   
     
    return 
  end subroutine FILE_restart_meshfield_component_def_var

  subroutine FILE_restart_meshfield_component_enddef( this )

    implicit none
    class(FILE_restart_meshfield_component), intent(inout) :: this
    !--------------------------------------------------------------
  
    call this%base%End_def()
    return
  end subroutine FILE_restart_meshfield_component_enddef

  subroutine FILE_restart_meshfield_component_write_var2d( this, &
    vid, field2d )
    
    use scale_time, only: TIME_NOWDAYSEC
    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this
    integer, intent(in) :: vid
    class(MeshField2D), intent(in) :: field2d
    !--------------------------------------------------

    call this%base%Write_var2D( vid, field2d, TIME_NOWDAYSEC, TIME_NOWDAYSEC )

    return
  end subroutine FILE_restart_meshfield_component_write_var2d

  subroutine FILE_restart_meshfield_component_write_var3d( this, &
    vid, field3d )
    
    use scale_time, only: TIME_NOWDAYSEC
    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this
    integer, intent(in) :: vid
    class(MeshField3D), intent(in) :: field3d
    !--------------------------------------------------

    call this%base%Write_var3D( vid, field3d, TIME_NOWDAYSEC, TIME_NOWDAYSEC )

    return
  end subroutine FILE_restart_meshfield_component_write_var3d

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

    call this%base%Read_Var( &
          dim_typeid, varname, field2d, step, allow_missing  )
    
    return
  end subroutine FILE_restart_meshfield_component_read_var2d

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

    call this%base%Read_Var( &
          dim_typeid, varname, field3d, step, allow_missing  )
    
    return
  end subroutine FILE_restart_meshfield_component_read_var3d

  subroutine FILE_restart_meshfield_component_close( this )
    implicit none

    class(FILE_restart_meshfield_component), intent(inout) :: this
    !--------------------------------------------------

    if ( this%base%fid /= -1 ) then
      LOG_NEWLINE
      LOG_INFO(trim(this%comp_name)//"_vars_restart_close",*) 'Close restart file'
      call this%base%Close()
    end if

    return
  end subroutine FILE_restart_meshfield_component_close


  subroutine FILE_restart_meshfield_component_Final( this )
    implicit none
    class(FILE_restart_meshfield_component), intent(inout) :: this
    !--------------------------------------------------

    call this%base%Final()

    return
  end subroutine FILE_restart_meshfield_component_Final

  !------------

end module scale_file_restart_meshfield
