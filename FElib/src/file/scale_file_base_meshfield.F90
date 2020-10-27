!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_file_base_meshfield
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_file_h
  use scale_prc, only: &
    PRC_abort
  use scale_file_h, only: &
    FILE_FILE_MAX
  
  use scale_file_common_meshfield, only: &
    FILE_common_meshfield_diminfo,                         &
    MF1D_DTYPE_NUM => FILE_COMMON_MESHFILED1D_DIMTYPE_NUM, &
    MF2D_DTYPE_NUM => FILE_COMMON_MESHFILED2D_DIMTYPE_NUM, &
    MF3D_DTYPE_NUM => FILE_COMMON_MESHFILED3D_DIMTYPE_NUM, &
    MF3D_DIMTYPE_X => FILE_COMMON_MESHFILED3D_DIMTYPEID_X, &
    MF3D_DIMTYPE_Y => FILE_COMMON_MESHFILED3D_DIMTYPEID_Y, &
    MF3D_DIMTYPE_Z => FILE_COMMON_MESHFILED3D_DIMTYPEID_Z, &
    get_dtype => File_common_meshfield_get_dtype

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
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
    
  type, public :: FILE_base_meshfield
    integer :: fid
    integer, allocatable :: vars_ncid(:)

    integer :: write_buf_amount
    logical :: File_axes_written

    class(MeshBase1D), pointer :: mesh1D
    class(MeshRectDom2D), pointer :: mesh2D
    class(MeshCubeDom3D), pointer :: mesh3D  
    type(FILE_common_meshfield_diminfo), allocatable :: dimsinfo(:)
  contains
    procedure :: Init => FILE_base_meshfield_Init
    procedure :: Open => FILE_base_meshfield_open
    procedure :: Create => FILE_base_meshfield_create
    !-
    procedure :: FILE_base_meshfield_def_var1
    procedure :: FILE_base_meshfield_def_var2
    generic :: Def_Var => FILE_base_meshfield_def_var1, FILE_base_meshfield_def_var2
    procedure :: End_def => FILE_base_meshfield_enddef
    !-
    procedure :: FILE_base_meshfield_write_var3d
    generic :: Write_var3D => FILE_base_meshfield_write_var3d
    !-
    procedure :: FILE_base_meshfield_read_var3d
    procedure :: FILE_base_meshfield_read_var3d_local
    generic :: Read_Var => FILE_base_meshfield_read_var3d, FILE_base_meshfield_read_var3d_local
    !-
    procedure :: Get_commonInfo => FILE_base_meshfield_get_commonInfo
    procedure :: Get_dataInfo => FILE_base_meshfield_get_dataInfo
    procedure :: Get_VarStepSize => FILE_base_meshfield_get_VarStepSize
    !-
    procedure :: Close => FILE_base_meshfield_close
    procedure :: Final => FILE_base_meshfield_Final
  end type FILE_base_meshfield

contains

  subroutine FILE_base_meshfield_Init( this,  &
    var_num, mesh1D, mesh2D, mesh3D )

    use scale_file_common_meshfield, only: &
      File_common_meshfield_get_dims  

    implicit none

    class(FILE_base_meshfield), intent(inout) :: this
    integer, intent(in) :: var_num
    class(MeshBase1D), target, optional, intent(in) :: mesh1D
    class(MeshRectDom2D), target, optional, intent(in) :: mesh2D
    class(MeshCubeDom3D), target, optional, intent(in) :: mesh3D

    logical :: check_specify_mesh

    !--------------------------------------------------

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
      LOG_ERROR("FILE_base_meshfield_Init",*) 'Specify a mesh among mesh1D, 2D, and 3D. Check!'
      call PRC_abort
    end if

    !-
    this%File_axes_written = .false.
    this%write_buf_amount = 0
  
    return
  end subroutine FILE_base_meshfield_Init

  subroutine FILE_base_meshfield_open( &
    this, basename, myrank )
  
    use scale_file, only: &
      FILE_open 
    implicit none
  
    class(FILE_base_meshfield), intent(inout) :: this
    character(*), intent(in) :: basename
    integer, intent(in), optional :: myrank
    !--------------------------------------------------------------
    
    call FILE_open( basename,         & ! [in]
                    this%fid,         & ! [out]
                    rankid=myrank     ) ! [in]
  
    return
  end subroutine FILE_base_meshfield_open

  subroutine FILE_base_meshfield_create( &
    this, basename, title, dtype,        &
    fileexisted,                         &
    myrank, tunits, calendar )
  
    use scale_file, only: &
      FILE_Create
    implicit none

    class(FILE_base_meshfield), intent(inout) :: this
    character(*), intent(in) :: basename  
    character(*), intent(in) :: title    
    character(*), intent(in) :: dtype
    logical, intent(out) :: fileexisted
    integer, intent(in), optional :: myrank
    character(*), intent(in), optional :: calendar
    character(*), intent(in), optional :: tunits
    !--------------------------------------------------------------
    
    call FILE_Create( basename,                & ! [IN]
                      title,                   & ! [IN]
                      H_SOURCE,                & ! [IN]
                      H_INSTITUTE,             & ! [IN]
                      this%fid,                & ! [OUT]
                      fileexisted,             & ! [OUT]
                      rankid     = myrank,     & ! [IN]
                      time_units = tunits,     & ! [IN]
                      calendar   = calendar    ) ! [IN]
    
    if ( .not. fileexisted ) then
      call def_axes( this, dtype )
    end if

    return
  end subroutine FILE_base_meshfield_create

  subroutine FILE_base_meshfield_def_var1( this, &
    field, desc, vid, dim_type_id, datatype, standard_name, timeinv, nsteps )
    implicit none
  
    class(FILE_base_meshfield), intent(inout) :: this
    class(MeshFieldBase), intent(in) :: field
    character(len=*), intent(in) :: desc
    integer, intent(in) :: dim_type_id
    integer, intent(in) :: vid
    character(len=*), intent(in) :: datatype
    character(len=*), optional, intent(in) :: standard_name
    real(DP), optional, intent(in) :: timeinv
    integer, optional, intent(in) :: nsteps
    !--------------------------------------------------------------

    call this%Def_Var( field%varname, field%unit, &
      desc, vid, dim_type_id, datatype, standard_name, timeinv, nsteps )
  
    return
  end subroutine FILE_base_meshfield_def_var1  

  subroutine FILE_base_meshfield_def_var2( this, &
    varname, units, desc, vid, dim_type_id, datatype, standard_name, timeinv, nsteps )
  
    use scale_file, only: &
      FILE_opened, &
      FILE_Def_Variable, &
      FILE_Set_Attribute  
    implicit none
  
    class(FILE_base_meshfield), intent(inout) :: this    
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: desc
    integer, intent(in) :: dim_type_id
    integer, intent(in) :: vid
    character(len=*), intent(in) :: datatype
    character(len=*), optional, intent(in) :: standard_name
    real(DP), optional, intent(in) :: timeinv
    integer, optional, intent(in) :: nsteps
  
    integer :: i_dtype
    integer :: ndim
    character(len=H_MID)   :: standard_name_  
    !--------------------------------------------------------------
  
    i_dtype = get_dtype(datatype)
  
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
      call FILE_Def_Variable( this%fid, varname, desc, units, standard_name_,            &
        ndim, this%dimsinfo(dim_type_id)%dims(1:ndim), i_dtype, this%vars_ncid(vid),     &
        time_int=timeinv )
    else
      call FILE_Def_Variable( this%fid, varname, desc, units, standard_name_,            &
        ndim, this%dimsinfo(dim_type_id)%dims(1:ndim), i_dtype, this%vars_ncid(vid)      )    
    end if
  
    return
  end subroutine FILE_base_meshfield_def_var2  

  subroutine FILE_base_meshfield_enddef( this )

    use scale_file, only: &
      FILE_opened,    &
      FILE_EndDef
    implicit none
    class(FILE_base_meshfield), intent(inout) :: this
  
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
  end subroutine FILE_base_meshfield_enddef
  
  subroutine FILE_base_meshfield_write_var3d( this, &
      vid, field3d, sec_str, sec_end )
  
    use scale_file, only: &
        FILE_opened, &
        FILE_Write 
    use scale_file_common_meshfield, only: &
      File_common_meshfield_put_field3D_cartesbuf
  
    class(FILE_base_meshfield), intent(inout) :: this
    integer, intent(in) :: vid
    class(MeshField3D), intent(in) :: field3d
    real(DP), intent(in) :: sec_str
    real(DP), intent(in) :: sec_end
    
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
        sec_str, sec_end, start                             ) ! (in)
    end if
  
    return
  end subroutine FILE_base_meshfield_write_var3d

  subroutine FILE_base_meshfield_get_commonInfo( this, &
    title, source, institution )
    use scale_file, only: &
      FILE_get_attribute
    implicit none
    
    class(FILE_base_meshfield), intent(in) :: this
    character(len=FILE_HMID), intent(out), optional :: title
    character(len=FILE_HMID), intent(out), optional :: source
    character(len=FILE_HMID), intent(out), optional :: institution
    !-------------------------------------------------

    if ( present(title) ) call FILE_get_attribute( this%fid, 'global', 'title', title )
    if ( present(source) ) call FILE_get_attribute( this%fid, 'global', 'source', source )
    if ( present(institution) ) call FILE_get_attribute( this%fid, 'global', 'institution', institution )

    return
  end subroutine FILE_base_meshfield_get_commonInfo

  subroutine FILE_base_meshfield_get_VarStepSize( this, varname, &
    len )
    use scale_file, only: &
      FILE_get_stepSize    
    implicit none

    class(FILE_base_meshfield), intent(in) :: this
    character(*), intent(in) :: varname
    integer, intent(out) :: len
    !-------------------------------------------------

    call FILE_get_stepSize( this%fid, varname, & ! (in)
      len )                                      ! (out)
    
    return    
  end subroutine FILE_base_meshfield_get_VarStepSize

  subroutine FILE_base_meshfield_get_dataInfo( this, varname, istep, &
    description, units, standard_name,                               &
    time_start, time_end, time_units, calendar )
    use scale_file, only: &
      FILE_get_dataInfo
    implicit none
    
    class(FILE_base_meshfield), intent(in) :: this
    character(*), intent(in) :: varname
    integer, intent(in), optional :: istep
    character(len=FILE_HMID), intent(out), optional :: description
    character(len=FILE_HSHORT), intent(out), optional :: units
    character(len=FILE_HMID), intent(out), optional :: standard_name
    real(DP), intent(out), optional :: time_start
    real(DP), intent(out), optional :: time_end
    character(len=FILE_HMID), intent(out), optional :: time_units
    character(len=FILE_HSHORT), intent(out), optional :: calendar
    !-------------------------------------------------

    call FILE_get_dataInfo( this%fid, varname, istep=istep,                               & ! (in)
      description=description, units=units, standard_name=standard_name,                  & ! (out)
      time_start=time_start, time_end=time_end, time_units=time_units, calendar=calendar  ) ! (out)

    return
  end subroutine FILE_base_meshfield_get_dataInfo

  subroutine FILE_base_meshfield_read_var3d( this,            &
    dim_typeid, varname, field3d, step, allow_missing )
  
    use scale_file, only: &
      FILE_Read
    use scale_file_common_meshfield, only: &
      File_common_meshfield_set_cartesbuf_field3D
  
    implicit none
  
    class(FILE_base_meshfield), intent(inout) :: this
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
  end subroutine FILE_base_meshfield_read_var3d

  subroutine FILE_base_meshfield_read_var3d_local( this, &
    dim_typeid, varname, lcmesh, i0_s, j0_s, k0_s, val,  &
    step, allow_missing  )
  
    use scale_file, only: &
      FILE_Read
    use scale_file_common_meshfield, only: &
      File_common_meshfield_set_cartesbuf_field3D_local
  
    implicit none
  
    class(FILE_base_meshfield), intent(inout) :: this
    integer, intent(in) :: dim_typeid
    character(*), intent(in) :: varname
    class(LocalMesh3D), intent(in) :: lcmesh
    integer, intent(in) :: i0_s, j0_s, k0_s
    real(RP), intent(out) :: val(lcmesh%refElem3D%Np,lcmesh%NeA)
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
  
      call File_common_meshfield_set_cartesbuf_field3D_local( &
        lcmesh, buf(:,:,:), i0_s, j0_s, k0_s,                 &
        val(:,:) )
    end if
  
    return
  end subroutine FILE_base_meshfield_read_var3d_local 

  subroutine FILE_base_meshfield_close( this )
    use scale_file, only: FILE_Close
    
    implicit none
    class(FILE_base_meshfield), intent(inout) :: this
    !--------------------------------------------------
  
    if ( this%fid /= -1 ) then
      call FILE_Close( this%fid ) ! [IN]
      this%fid = -1
    end if
  
    return
  end subroutine FILE_base_meshfield_close

  
  subroutine FILE_base_meshfield_Final( this )
    implicit none
    class(FILE_base_meshfield), intent(inout) :: this
    !--------------------------------------------------
  
    if ( allocated(this%vars_ncid) ) deallocate( this%vars_ncid )
    if ( allocated(this%dimsinfo) ) deallocate( this%dimsinfo )
    nullify( this%mesh1D, this%mesh2D, this%mesh3D )
  
    return
  end subroutine FILE_base_meshfield_Final

  !- private -----------------------------------------

  subroutine def_axes( this, dtype )
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    use scale_file, only: &
      FILE_Def_Axis,                 &
      FILE_Set_Attribute,            &
      FILE_Def_AssociatedCoordinate, &
      FILE_Add_AssociatedVariable

    implicit none

    class(FILE_base_meshfield), intent(in) :: this
    character(*), intent(in) :: dtype
    integer :: d
    integer :: i_dtype
    !------------

    i_dtype = get_dtype( dtype )

    if ( associated(this%mesh1D) ) then
      do d=1, 1
        call FILE_Def_Axis( this%fid, &
          this%dimsinfo(d)%name, this%dimsinfo(d)%desc, this%dimsinfo(d)%unit, &
          this%dimsinfo(d)%name, i_dtype, this%dimsinfo(d)%size                )
      end do
    end if

    if ( associated(this%mesh2D) ) then
      do d=1, 2
        call FILE_Def_Axis( this%fid, &
          this%dimsinfo(d)%name, this%dimsinfo(d)%desc, this%dimsinfo(d)%unit, &
          this%dimsinfo(d)%name, i_dtype, this%dimsinfo(d)%size                )
      end do
    end if

    if ( associated(this%mesh3D) ) then
      do d=1, 3
        call FILE_Def_Axis( this%fid, &
          this%dimsinfo(d)%name, this%dimsinfo(d)%desc, this%dimsinfo(d)%unit, &
          this%dimsinfo(d)%name, i_dtype, this%dimsinfo(d)%size                )
      end do
    end if

    return
  end subroutine def_axes

  subroutine write_axes( this, start )
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    use scale_file, only: &
      FILE_Write_Axis
    use scale_file_common_meshfield, only: &
      File_common_meshfield_get_axis
    implicit none

    class(FILE_base_meshfield), intent(in) :: this
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

  
end module scale_file_base_meshfield
