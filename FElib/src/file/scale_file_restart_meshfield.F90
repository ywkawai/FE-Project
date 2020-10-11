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
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D

  use scale_meshfield_base, only: MeshFieldBase, MeshField1D, MeshField2D, MeshField3D
  use scale_variableinfo, only: VariableInfo

  use scale_file_common_meshfield, only: &
    MF1D_DTYPE_NUM => FILE_COMMON_MESHFILED1D_DIMTYPE_NUM, &
    MF2D_DTYPE_NUM => FILE_COMMON_MESHFILED2D_DIMTYPE_NUM, &
    MF3D_DTYPE_NUM => FILE_COMMON_MESHFILED3D_DIMTYPE_NUM, &
    MF3D_DIMTYPE_X => FILE_COMMON_MESHFILED3D_DIMTYPEID_X, &
    MF3D_DIMTYPE_Y => FILE_COMMON_MESHFILED3D_DIMTYPEID_Y, &
    MF3D_DIMTYPE_Z => FILE_COMMON_MESHFILED3D_DIMTYPEID_Z
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type :: FILE_restart_meshfield
    logical :: flag_output
    character(len=H_LONG) :: in_basename
    logical :: in_postfix_timelabel

    character(len=H_LONG) :: out_basename
    logical :: out_postfix_timelabel
    character(len=H_MID) :: out_title
    character(len=H_SHORT) :: out_dtype

    integer :: fid
    integer :: write_buf_amount
    logical :: File_axes_written

    class(MeshBase1D), pointer :: mesh1D
    class(MeshRectDom2D), pointer :: mesh2D
    class(MeshCubeDom3D), pointer :: mesh3D  
    type(FILE_common_meshfield_diminfo), allocatable :: dimsinfo(:)
  end type
  public :: FILE_restart_meshfield_setup


  type, extends(FILE_restart_meshfield), public :: FILE_restart_meshfield_component
    character(len=H_SHORT) :: comp_name
    integer, allocatable :: vars_ncid(:)
  contains
    procedure :: Init1 => FILE_restart_meshfield_component_Init1
    procedure :: Init2 => FILE_restart_meshfield_component_Init2
    generic :: Init => Init1, Init2
    procedure :: Open => FILE_restart_meshfield_component_open
   procedure :: Read_var3D => FILE_restart_meshfield_component_read_var3d
   generic :: Read_var => Read_var3D
    procedure :: Create => FILE_restart_meshfield_component_create
    procedure :: Def_var => FILE_restart_meshfield_component_def_var
    procedure :: End_def => FILE_restart_meshfield_component_enddef
    procedure :: Write_var3D => FILE_restart_meshfield_component_write_var3d
    generic :: Write_var => Write_var3D
    procedure :: Close => FILE_restart_meshfield_component_close
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

  type(FILE_restart_meshfield) :: restart_file

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

  restart_file%fid = -1

  return
end subroutine FILE_restart_meshfield_setup

subroutine FILE_restart_meshfield_component_Init1( this,  &
  comp_name,                                              &
  var_num, mesh1D, mesh2D, mesh3D )

use scale_file_common_meshfield, only: &
  File_common_meshfield_get_dims  

  implicit none

  class(FILE_restart_meshfield_component), intent(inout) :: this
  character(*), intent(in) :: comp_name
  integer, intent(in) :: var_num
  class(MeshBase1D), target, optional, intent(in) :: mesh1D
  class(MeshRectDom2D), target, optional, intent(in) :: mesh2D
  class(MeshCubeDom3D), target, optional, intent(in) :: mesh3D

  !--------------------------------------------------

  call this%Init( &
    comp_name,                                                     &
    restart_file%in_basename, restart_file%in_postfix_timelabel,   &
    restart_file%out_basename, restart_file%out_postfix_timelabel, &
    restart_file%out_dtype, restart_file%out_title,                &
    var_num, mesh1D, mesh2D, mesh3D )

  return
end subroutine FILE_restart_meshfield_component_Init1

subroutine FILE_restart_meshfield_component_Init2( this,  &
    comp_name,                                            &
    in_basename, in_postfix_timelabel,                    &
    out_basename, out_postfix_timelabel,                  &
    out_dtype, out_title,                                 &
    var_num, mesh1D, mesh2D, mesh3D )

  use scale_file_common_meshfield, only: &
    File_common_meshfield_get_dims
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
  class(MeshCubeDom3D), target, optional, intent(in) :: mesh3D

  logical :: check_specify_mesh

  !--------------------------------------------------

  this%comp_name = comp_name

  this%in_basename = in_basename
  this%in_postfix_timelabel = in_postfix_timelabel

  this%out_basename = out_basename
  this%out_postfix_timelabel = out_postfix_timelabel
  this%out_title = out_title
  this%out_dtype = out_dtype

  this%fid = -1

  allocate( this%vars_ncid(var_num) )
  this%vars_ncid(:) = -1

  !-
  check_specify_mesh = .false.
  nullify( this%mesh1D, this%mesh2D, this%mesh3D )

  if (present(mesh1D)) then
    this%mesh1D => mesh1D
    check_specify_mesh = .true.

    allocate( this%dimsinfo(MF1D_DTYPE_NUM) )
    call File_common_meshfield_get_dims( mesh1D, this%dimsinfo(:) )
  end if
  if (present(mesh2D)) then
    this%mesh2D => mesh2D
    check_specify_mesh = .true.

    allocate( this%dimsinfo(MF2D_DTYPE_NUM) )
    call File_common_meshfield_get_dims( mesh2D, this%dimsinfo(:) )
  end if
  if (present(mesh3D)) then
    this%mesh3D => mesh3D
    check_specify_mesh = .true.

    allocate( this%dimsinfo(MF3D_DTYPE_NUM) )
    call File_common_meshfield_get_dims( mesh3D, this%dimsinfo(:) )
  end if

  if (.not. check_specify_mesh) then
    LOG_ERROR("FILE_restart_meshfield_component_Init",*) 'Specify a mesh among mesh1D, 2D, and 3D. Check!'
    call PRC_abort
  end if

  !-
  this%File_axes_written = .false.
  this%write_buf_amount = 0

  return
end subroutine FILE_restart_meshfield_component_Init2

subroutine FILE_restart_meshfield_component_open( &
  this )

  use scale_file, only: &
    FILE_open
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

  LOG_INFO(trim(this%comp_name)//"_vars_restart_create",*) 'basename: ', trim(basename)

  !--------------------------------

  LOG_NEWLINE
  LOG_INFO(trim(this%comp_name)//"_vars_restart_open",*) 'Open restart file'
  call FILE_open( basename, this%fid )

  return
end subroutine FILE_restart_meshfield_component_open


subroutine FILE_restart_meshfield_component_create( &
    this )
  
  use scale_file, only: &
    FILE_Create
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

  call get_tunits_and_calendarname( NOWDATE, &
    tunits, calendar )
  
  call FILE_Create( basename,                & ! [IN]
                    this%out_title,          & ! [IN]
                    H_SOURCE,                & ! [IN]
                    H_INSTITUTE,             & ! [IN]
                    this%fid,                & ! [OUT]
                    fileexisted,             & ! [OUT]
                    rankid     = PRC_myrank, & ! [IN]
                    time_units = tunits,     & ! [IN]
                    calendar   = calendar    ) ! [IN]
  
  if ( .not. fileexisted ) then
    call put_global_attribute( this%fid, NOWSUBSEC, tunits, calendar )
    call def_axes( this )
  end if

  return
end subroutine FILE_restart_meshfield_component_create

subroutine FILE_restart_meshfield_component_def_var( this, &
  field, desc, vid, dim_type_id, datatype, standard_name, timeinv, nsteps )

  use scale_file, only: &
    FILE_opened, &
    FILE_Def_Variable, &
    FILE_Set_Attribute  
  implicit none

  class(FILE_restart_meshfield_component), intent(inout) :: this
  class(MeshFieldBase), intent(in) :: field
  character(len=*), intent(in) :: desc
  integer, intent(in) :: dim_type_id
  integer, intent(in) :: vid
  character(len=*), optional, intent(in) :: datatype
  character(len=*), optional, intent(in) :: standard_name
  real(DP), optional, intent(in) :: timeinv
  integer, optional, intent(in) :: nsteps

  integer :: dtype
  integer :: v
  integer :: ndim
  character(len=H_MID)   :: standard_name_  
  !--------------------------------------------------------------

  if ( present(datatype) ) then
    dtype = get_dtype(datatype)
  else
    dtype = get_dtype(this%out_dtype)
  end if

  if ( present(nsteps) ) then
    this%write_buf_amount = this%write_buf_amount + this%dimsinfo(dim_type_id)%size * nsteps
  else
    this%write_buf_amount = this%write_buf_amount + this%dimsinfo(dim_type_id)%size
  end if
  if ( present(standard_name) ) then
    standard_name_ = standard_name
  else
    standard_name_ = ""
  end if

  ndim = this%dimsinfo(dim_type_id)%ndim
  if ( present(timeinv) ) then
    call FILE_Def_Variable( this%fid, field%varname, desc, field%unit, standard_name_, &
      ndim, this%dimsinfo(dim_type_id)%dims(1:ndim), dtype, this%vars_ncid(vid),       &
      time_int=timeinv )
  else
    call FILE_Def_Variable( this%fid, field%varname, desc, field%unit, standard_name_, &
      ndim, this%dimsinfo(dim_type_id)%dims(1:ndim), dtype, this%vars_ncid(vid)        )    
  end if

  return
end subroutine FILE_restart_meshfield_component_def_var

subroutine FILE_restart_meshfield_component_enddef( this )

  use scale_file, only: &
    FILE_opened,    &
    FILE_EndDef

  implicit none
  class(FILE_restart_meshfield_component), intent(inout) :: this

  integer :: start(3)
  !--------------------------------------------------------------

  if (this%fid == -1) return

  call FILE_EndDef( this%fid )
  
  if ( .not. this%File_axes_written ) then
    start(:) = 1
    call write_axes( this, start(:) )
    this%File_axes_written = .true.
  end if

  return
end subroutine FILE_restart_meshfield_component_enddef

subroutine FILE_restart_meshfield_component_write_var3d( this, &
    vid, field3d )

  use scale_file, only: &
      FILE_opened, &
      FILE_Write
  use scale_time, only: &
    NOWDAYSEC  => TIME_NOWDAYSEC  
  use scale_file_common_meshfield, only: &
    File_common_meshfield_put_field3D_cartesbuf

  class(FILE_restart_meshfield_component), intent(inout) :: this
  integer, intent(in) :: vid
  class(MeshField3D), intent(in) :: field3d

  real(RP), allocatable :: buf(:,:,:)
  integer :: dims(3)
  integer :: start(3)
  !-------------------------------------------------

  if ( this%fid /= -1 ) then
    start(:) = 1
    dims(1) = this%dimsinfo(MF3D_DIMTYPE_X)%size
    dims(2) = this%dimsinfo(MF3D_DIMTYPE_Y)%size
    dims(3) = this%dimsinfo(MF3D_DIMTYPE_Z)%size
    allocate( buf(dims(1),dims(2),dims(3)) )
    call File_common_meshfield_put_field3D_cartesbuf( this%mesh3D, field3d, buf(:,:,:) )

    call FILE_Write( this%vars_ncid(vid), buf(:,:,:),     & ! (in)
      NOWDAYSEC, NOWDAYSEC, start  ) ! (in)
  end if

  return
end subroutine FILE_restart_meshfield_component_write_var3d

subroutine FILE_restart_meshfield_component_read_var3d( this, &
  dim_typeid, varname, field3d, step, allow_missing )

  use scale_file, only: &
    FILE_Read
  use scale_time, only: &
    NOWDAYSEC  => TIME_NOWDAYSEC  
  use scale_file_common_meshfield, only: &
    File_common_meshfield_set_cartesbuf_field3D

  implicit none

  class(FILE_restart_meshfield_component), intent(inout) :: this
  integer, intent(in) :: dim_typeid
  character(*), intent(in) :: varname
  class(MeshField3D), intent(inout) :: field3d
  integer, intent(in), optional :: step
  logical, intent(in), optional :: allow_missing

  real(RP), allocatable :: buf(:,:,:)
  integer :: dims(3)
  integer :: start(3)   ! start offset of globale variable
  !-------------------------------------------------

  if ( this%fid /= -1 ) then
    start(:) = 1
    dims(1) = this%dimsinfo(MF3D_DIMTYPE_X)%size
    dims(2) = this%dimsinfo(MF3D_DIMTYPE_Y)%size
    dims(3) = this%dimsinfo(MF3D_DIMTYPE_Z)%size
    allocate( buf(dims(1),dims(2),dims(3)) )

    call FILE_Read( this%fid, varname,                       & ! (in)
      buf(:,:,:),                                            & ! (out)
      step=step, allow_missing=allow_missing                 ) ! (in)

    call File_common_meshfield_set_cartesbuf_field3D( this%mesh3D, buf(:,:,:), &
      field3d )
  end if

  return
end subroutine FILE_restart_meshfield_component_read_var3d

subroutine FILE_restart_meshfield_component_close( this )
  use scale_file, only: FILE_Close
  
  implicit none

  class(FILE_restart_meshfield_component), intent(inout) :: this
  !--------------------------------------------------

  if ( this%fid /= -1 ) then
    LOG_NEWLINE
    LOG_INFO(trim(this%comp_name)//"_vars_restart_close",*) 'Close restart file'

    call FILE_Close( this%fid ) ! [IN]
    this%fid = -1
  end if

  return
end subroutine FILE_restart_meshfield_component_close


subroutine FILE_restart_meshfield_component_Final( this )
  implicit none
  class(FILE_restart_meshfield_component), intent(inout) :: this
  !--------------------------------------------------

  if ( allocated(this%vars_ncid) ) deallocate( this%vars_ncid )
  if ( allocated(this%dimsinfo) ) deallocate( this%dimsinfo )
  nullify( this%mesh1D, this%mesh2D, this%mesh3D )

  return
end subroutine FILE_restart_meshfield_component_Final

!- private -----------------------------------------

subroutine file_restart_setup( this )
  use scale_const, only: &
    EPS => CONST_EPS
  use scale_file, only: &
    FILE_setup    
  implicit none

  class(FILE_restart_meshfield_component), intent(inout) :: this
  
  integer :: ierr
  !--------------------------------------------------
  
  call FILE_setup( PRC_myrank )

  return
end subroutine file_restart_setup

subroutine file_restart_open( &
  basename, &
  fid,      &
  aggregate )
  use scale_file_h, only: &
    FILE_FREAD
  use scale_file, only: &
    FILE_AGGREGATE, &
    FILE_Open
  implicit none

  character(len=*), intent(in)  :: basename !< basename of the file
  integer,          intent(out) :: fid      !< file ID
  logical,          intent(in), optional :: aggregate
  !---------------------------------------------------------------------------

  call FILE_Open( basename,          & ! [IN]
                fid,                 & ! [OUT]
                aggregate=aggregate, & ! [IN]
                rankid=PRC_myrank    ) ! [IN]

  return
end subroutine file_restart_open

subroutine def_axes( this )
  use scale_const, only: &
    UNDEF => CONST_UNDEF
  use scale_file, only: &
    FILE_Def_Axis,                 &
    FILE_Set_Attribute,            &
    FILE_Def_AssociatedCoordinate, &
    FILE_Add_AssociatedVariable

  implicit none

  class(FILE_restart_meshfield_component), intent(in) :: this
  integer :: d, n, ndim
  integer :: dtype
  !------------

  dtype = get_dtype( this%out_dtype )

  if ( associated(this%mesh1D) ) then
    do d=1, 1
      call FILE_Def_Axis( this%fid, &
        this%dimsinfo(d)%name, this%dimsinfo(d)%desc, this%dimsinfo(d)%unit, &
        this%dimsinfo(d)%name, dtype, this%dimsinfo(d)%size                  )
    end do
  end if

  if ( associated(this%mesh2D) ) then
    do d=1, 2
      call FILE_Def_Axis( this%fid, &
        this%dimsinfo(d)%name, this%dimsinfo(d)%desc, this%dimsinfo(d)%unit, &
        this%dimsinfo(d)%name, dtype, this%dimsinfo(d)%size                  )
    end do
  end if

  if ( associated(this%mesh3D) ) then
    do d=1, 3
      call FILE_Def_Axis( this%fid, &
        this%dimsinfo(d)%name, this%dimsinfo(d)%desc, this%dimsinfo(d)%unit, &
        this%dimsinfo(d)%name, dtype, this%dimsinfo(d)%size                  )
    end do
  end if

  return
end subroutine def_axes

subroutine write_axes( this, start )
  use scale_const, only: &
    UNDEF => CONST_UNDEF
  use scale_file, only: &
    FILE_Write_Axis,    &
    FILE_Write_AssociatedCoordinate
  use scale_file_common_meshfield, only: &
    File_common_meshfield_get_axis
  implicit none

  class(FILE_restart_meshfield_component), intent(in) :: this
  integer, intent(in) :: start(3)
  real(RP), allocatable :: x(:), y(:), z(:)
  !------------

  if ( associated(this%mesh1D) ) then
    allocate( x(this%dimsinfo(1)%size) )
    call File_common_meshfield_get_axis(this%mesh1D, this%dimsinfo, x(:))

    call FILE_Write_Axis( this%fid, this%dimsinfo(1)%name, x(:), start(1:1) )
  end if

  if ( associated(this%mesh2D) ) then
    allocate( x(this%dimsinfo(1)%size), y(this%dimsinfo(2)%size) )
    call File_common_meshfield_get_axis(this%mesh2D, this%dimsinfo, x(:), y(:))

    call FILE_Write_Axis( this%fid, this%dimsinfo(1)%name, x(:), start(1:1) )
    call FILE_Write_Axis( this%fid, this%dimsinfo(2)%name, y(:), start(2:2) )
  end if

  if ( associated(this%mesh3D) ) then
    allocate( x(this%dimsinfo(1)%size), y(this%dimsinfo(2)%size), z(this%dimsinfo(3)%size) )
    call File_common_meshfield_get_axis(this%mesh3D, this%dimsinfo, x(:), y(:), z(:))

    call FILE_Write_Axis( this%fid, this%dimsinfo(1)%name, x(:), start(1:1) )
    call FILE_Write_Axis( this%fid, this%dimsinfo(2)%name, y(:), start(2:2) )
    call FILE_Write_Axis( this%fid, this%dimsinfo(3)%name, z(:), start(3:3) )
  end if  

  return
end subroutine write_axes

subroutine put_global_attribute( &
  fid, time, tunits, calendar  )

  use scale_file, only: &
   FILE_Set_Attribute
  implicit none

  integer, intent(in) :: fid
  real(DP), intent(in)     :: time
  character(*), intent(in) :: tunits
  character(*), intent(in) :: calendar

  !------------------------------------

  call FILE_Set_Attribute( fid, "global", "Conventions", "CF-1.6" ) ! [IN]
  call FILE_Set_Attribute( fid, "global", "grid_name", "hoge" ) ! [IN]

  if ( calendar /= "" ) call FILE_Set_Attribute( fid, "global", "calendar", calendar )
  call FILE_Set_Attribute( fid, "global", "time_units", tunits )
  call FILE_Set_Attribute( fid, "global", "time_start", (/time/) )

  return
end subroutine put_global_attribute

!------------

function get_dtype( datatype ) result( dtype )

  use scale_file_h, only: &
    FILE_REAL8, FILE_REAL4
  implicit none

  character(*), intent(in) :: datatype
  integer :: dtype
  !--------------------------

  ! dtype is used to define the data type of axis variables in file
  if    ( datatype == 'REAL8' ) then
    dtype = FILE_REAL8
  elseif( datatype == 'REAL4' ) then
      dtype = FILE_REAL4
  else
    if    ( RP == 8 ) then
      dtype = FILE_REAL8
    elseif( RP == 4 ) then
      dtype = FILE_REAL4
    else
      LOG_ERROR("file_restart_meshfield_get_dtype",*) 'unsupported data type. Check!', trim(datatype)
      call PRC_abort
    endif
  endif

  return
end function get_dtype

subroutine get_tunits_and_calendarname( date, &
    tunits, calendar_name )

  use scale_file, only: &
    FILE_get_CFtunits
  use scale_calendar, only: &
    CALENDAR_get_name    
  implicit none

  integer, intent(in) :: date(6)
  character(len=34), intent(out) :: tunits
  character(len=H_SHORT), intent(out) :: calendar_name
  !--------------------------------------------------

  if ( date(1) > 0 ) then
    call FILE_get_CFtunits( date(:), tunits )
    call CALENDAR_get_name( calendar_name )
  else
    tunits        = 'seconds'
    calendar_name = ''
  endif

  return
end subroutine get_tunits_and_calendarname

end module scale_file_restart_meshfield
