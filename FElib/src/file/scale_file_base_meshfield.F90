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
    FILE_common_meshfield_diminfo,               &
    get_dtype => File_common_meshfield_get_dtype
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
    class(MeshCubedSphereDom2D), pointer :: meshCS2D
    class(MeshCubeDom3D), pointer :: mesh3D  
    class(MeshCubedSphereDom3D), pointer :: meshCS3D
    type(FILE_common_meshfield_diminfo), allocatable :: dimsinfo(:)

    logical :: force_uniform_grid
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
    procedure :: FILE_base_meshfield_write_var1d
    generic :: Write_var1D => FILE_base_meshfield_write_var1d
    procedure :: FILE_base_meshfield_write_var2d
    generic :: Write_var2D => FILE_base_meshfield_write_var2d
    procedure :: FILE_base_meshfield_write_var3d
    generic :: Write_var3D => FILE_base_meshfield_write_var3d
    !-
    procedure :: Put_GlobalAttribute_time => FILE_base_meshfield_put_global_attribute_time
    !-
    procedure :: FILE_base_meshfield_read_var1d
    procedure :: FILE_base_meshfield_read_var1d_local
    procedure :: FILE_base_meshfield_read_var2d
    procedure :: FILE_base_meshfield_read_var2d_local
    procedure :: FILE_base_meshfield_read_var3d
    procedure :: FILE_base_meshfield_read_var3d_local
    generic :: Read_Var => &
      FILE_base_meshfield_read_var1d, FILE_base_meshfield_read_var1d_local, &
      FILE_base_meshfield_read_var2d, FILE_base_meshfield_read_var2d_local, &
      FILE_base_meshfield_read_var3d, FILE_base_meshfield_read_var3d_local
    
    !-  
    procedure :: Get_commonInfo => FILE_base_meshfield_get_commonInfo
    procedure :: Get_dataInfo => FILE_base_meshfield_get_dataInfo
    procedure :: Get_VarStepSize => FILE_base_meshfield_get_VarStepSize
    !-
    procedure :: Close => FILE_base_meshfield_close
    procedure :: Final => FILE_base_meshfield_Final
  end type FILE_base_meshfield

contains

  subroutine FILE_base_meshfield_Init( this, & ! (inout)
    var_num,                                 & ! (in)
    mesh1D,                                  & ! (in)
    mesh2D, meshCubedSphere2D,               & ! (in)
    mesh3D, meshCubedSphere3D,               & ! (in)
    force_uniform_grid )                       ! (in)

    use scale_file_common_meshfield, only: &
      File_common_meshfield_get_dims  

    implicit none

    class(FILE_base_meshfield), intent(inout) :: this
    integer, intent(in) :: var_num
    class(MeshBase1D), target, optional, intent(in) :: mesh1D
    class(MeshRectDom2D), target, optional, intent(in) :: mesh2D
    class(MeshCubedSphereDom2D), target, optional, intent(in) :: meshCubedSphere2D    
    class(MeshCubeDom3D), target, optional, intent(in) :: mesh3D
    class(MeshCubedSphereDom3D), target, optional, intent(in) :: meshCubedSphere3D
    logical, intent(in), optional :: force_uniform_grid

    logical :: check_specify_mesh
    !--------------------------------------------------

    this%fid = -1

    allocate( this%vars_ncid(var_num) )
    this%vars_ncid(:) = -1
  
    !-
    check_specify_mesh = .false.
    nullify( this%mesh1D, this%mesh2D, this%mesh3D )
    nullify( this%meshCS2D, this%meshCS3D )
  
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
    if (present(meshCubedSphere2D)) then
      this%meshCS2D => meshCubedSphere2D
      check_specify_mesh = .true.
  
      allocate( this%dimsinfo(MF2D_DTYPE_NUM) )
      call File_common_meshfield_get_dims( meshCubedSphere2D, this%dimsinfo(:) )
    end if
    if (present(mesh3D)) then
      this%mesh3D => mesh3D
      check_specify_mesh = .true.
  
      allocate( this%dimsinfo(MF3D_DTYPE_NUM) )
      call File_common_meshfield_get_dims( mesh3D, this%dimsinfo(:) )
    end if
    if (present(meshCubedSphere3D)) then
      this%meshCS3D => meshCubedSphere3D
      check_specify_mesh = .true.
  
      allocate( this%dimsinfo(MF3D_DTYPE_NUM) )
      call File_common_meshfield_get_dims( meshCubedSphere3D, this%dimsinfo(:) )
    end if

    if ( present(force_uniform_grid) ) then
      this%force_uniform_grid = force_uniform_grid
    else
      this%force_uniform_grid = .false.
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

  subroutine FILE_base_meshfield_open( this, & ! (inout)
    basename, myrank                         ) ! (in)
  
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
      this%File_axes_written = .false.      
    end if

    return
  end subroutine FILE_base_meshfield_create

  subroutine FILE_base_meshfield_def_var1( this, & ! (inout)
    field, desc, vid, dim_type_id, datatype,     & ! (in)
    standard_name, timeinv, nsteps               ) ! (in)
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

  subroutine FILE_base_meshfield_def_var2( this,      & ! (inout)
    varname, units, desc, vid, dim_type_id, datatype, & ! (in)
    standard_name, timeinv, nsteps                    ) ! (in)
  
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

  subroutine FILE_base_meshfield_enddef( this ) ! (inout)

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
  
!OCL_SERIAL
  subroutine FILE_base_meshfield_write_var1d( this, & ! (inout)
    vid, field1d, sec_str, sec_end                  ) ! (in)

    use scale_file, only: &
        FILE_opened, &
        FILE_Write 
    use scale_file_common_meshfield, only: &
      File_common_meshfield_put_field1D_cartesbuf
    implicit none

    class(FILE_base_meshfield), intent(inout) :: this
    integer, intent(in) :: vid
    class(MeshField1D), intent(in) :: field1d
    real(DP), intent(in) :: sec_str
    real(DP), intent(in) :: sec_end
    
    real(RP), allocatable :: buf(:)
    integer :: dims(1)
    integer :: start(1)
    !-------------------------------------------------

    if ( this%fid /= -1 ) then
      start(:) = 1
      dims(1) = this%dimsinfo(MF1D_DIMTYPE_X)%size
      allocate( buf(dims(1)) )
      call File_common_meshfield_put_field1D_cartesbuf( this%mesh1D, field1d, buf(:), &
        this%force_uniform_grid )

      call FILE_Write( this%vars_ncid(vid), buf(:),     & ! (in)
        sec_str, sec_end, start=start                   ) ! (in)
    end if

    return
  end subroutine FILE_base_meshfield_write_var1d

!OCL_SERIAL
  subroutine FILE_base_meshfield_write_var2d( this, & ! (inout)
    vid, field2d, sec_str, sec_end                  ) ! (in)

    use scale_file, only: &
        FILE_opened, &
        FILE_Write 
    use scale_file_common_meshfield, only: &
      File_common_meshfield_put_field2D_cartesbuf,             &
      File_common_meshfield_put_field2D_cubedsphere_cartesbuf
    implicit none

    class(FILE_base_meshfield), intent(inout) :: this
    integer, intent(in) :: vid
    class(MeshField2D), intent(in) :: field2d
    real(DP), intent(in) :: sec_str
    real(DP), intent(in) :: sec_end
    
    real(RP), allocatable :: buf(:,:)
    integer :: dims(2)
    integer :: start(2)
    !-------------------------------------------------

    if ( this%fid /= -1 ) then
      start(:) = 1
      dims(1) = this%dimsinfo(MF2D_DIMTYPE_X)%size
      dims(2) = this%dimsinfo(MF2D_DIMTYPE_Y)%size
      allocate( buf(dims(1),dims(2)) )
      if ( associated(this%mesh2D) ) then
        call File_common_meshfield_put_field2D_cartesbuf( this%mesh2D, field2d, buf(:,:), &
          this%force_uniform_grid )
      else if ( associated(this%meshCS2D) ) then
        call File_common_meshfield_put_field2D_cubedsphere_cartesbuf( &
          this%meshCS2D, field2d, buf(:,:)                            )
      end if

      call FILE_Write( this%vars_ncid(vid), buf(:,:),   & ! (in)
        sec_str, sec_end, start=start                   ) ! (in)
    end if

    return
  end subroutine FILE_base_meshfield_write_var2d

!OCL_SERIAL
  subroutine FILE_base_meshfield_write_var3d( this, & ! (inout)
      vid, field3d, sec_str, sec_end                ) ! (in)
  
    use scale_file, only: &
        FILE_opened, &
        FILE_Write 
    use scale_file_common_meshfield, only: &
      File_common_meshfield_put_field3D_cartesbuf,            &
      File_common_meshfield_put_field3D_cubedsphere_cartesbuf
    use scale_prof
    implicit none

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

      if ( associated(this%mesh3D) ) then
        call File_common_meshfield_put_field3D_cartesbuf( this%mesh3D, field3d, buf(:,:,:), &
          this%force_uniform_grid )
      else if ( associated(this%meshCS3D) ) then
        call File_common_meshfield_put_field3D_cubedsphere_cartesbuf( &
          this%meshCS3D, field3d, buf(:,:,:)                          )
      end if

      call FILE_Write( this%vars_ncid(vid), buf(:,:,:),     & ! (in)
        sec_str, sec_end, start                             ) ! (in)
    end if
  
    return
  end subroutine FILE_base_meshfield_write_var3d

  subroutine FILE_base_meshfield_get_commonInfo( this, & ! (in)
    title, source, institution                         ) ! (out)
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

  subroutine FILE_base_meshfield_get_VarStepSize( this, varname, & ! (in)
    len                                                          ) ! (out)
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

  subroutine FILE_base_meshfield_get_dataInfo( this, varname, istep, & ! (in)
    description, units, standard_name,                               & ! (out)
    time_start, time_end, time_units, calendar                       ) ! (out)
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

!OCL_SERIAL
  subroutine FILE_base_meshfield_read_var1d( this,    & ! (inout)
    dim_typeid, varname,                              & ! (in)
    field1d,                                          & ! (inout)
    step, allow_missing )                               ! (in)
  
    use scale_file, only: &
      FILE_Read
    use scale_file_common_meshfield, only: &
      File_common_meshfield_set_cartesbuf_field1D

    implicit none
  
    class(FILE_base_meshfield), intent(inout) :: this
    integer, intent(in) :: dim_typeid
    character(*), intent(in) :: varname
    class(MeshField1D), intent(inout) :: field1d
    integer, intent(in), optional :: step
    logical, intent(in), optional :: allow_missing
  
    real(RP), allocatable :: buf(:)
    integer :: dims(1)
    integer :: start(1)   ! start offset of globale variable
    !-------------------------------------------------
  
    if ( this%fid /= -1 ) then
      start(:) = 1
      dims(1) = this%dimsinfo(dim_typeid)%size
      allocate( buf(dims(1)) )
  
      call FILE_Read( this%fid, varname,                       & ! (in)
        buf(:),                                                & ! (out)
        step=step, allow_missing=allow_missing                 ) ! (in)
  
      call File_common_meshfield_set_cartesbuf_field1D( this%mesh1D, buf(:), &
        field1d )
    end if
  
    return
  end subroutine FILE_base_meshfield_read_var1d

!OCL_SERIAL
  subroutine FILE_base_meshfield_read_var1d_local( this, & ! (inout)
    dim_typeid, varname, lcmesh, i0_s,                   & ! (in)
    val,                                                 & ! (out)
    step, allow_missing  )                                 ! (in)
  
    use scale_file, only: &
      FILE_Read
    use scale_file_common_meshfield, only: &
      File_common_meshfield_set_cartesbuf_field1D_local
  
    implicit none
  
    class(FILE_base_meshfield), intent(inout) :: this
    integer, intent(in) :: dim_typeid
    character(*), intent(in) :: varname
    class(LocalMesh1D), intent(in) :: lcmesh
    integer, intent(in) :: i0_s
    real(RP), intent(out) :: val(lcmesh%refElem1D%Np,lcmesh%NeA)
    integer, intent(in), optional :: step
    logical, intent(in), optional :: allow_missing
  
    real(RP), allocatable :: buf(:)
    integer :: dims(1)
    integer :: start(1)   ! start offset of globale variable
    !-------------------------------------------------
  
    if ( this%fid /= -1 ) then
      start(:) = 1
      dims(1) = this%dimsinfo(dim_typeid)%size
      allocate( buf(dims(1)) )
  
      call FILE_Read( this%fid, varname,                       & ! (in)
        buf(:),                                                & ! (out)
        step=step, allow_missing=allow_missing                 ) ! (in)
  
      call File_common_meshfield_set_cartesbuf_field1D_local( &
        lcmesh, buf(:), i0_s,                                 &
        val(:,:) )
    end if
  
    return
  end subroutine FILE_base_meshfield_read_var1d_local 

!OCL_SERIAL
  subroutine FILE_base_meshfield_read_var2d( this, & ! (inout)
    dim_typeid, varname,                           & ! (in)
    field2d,                                       & ! (inout)
    step, allow_missing                            ) ! (in)
  
    use scale_file, only: &
      FILE_Read
    use scale_file_common_meshfield, only: &
      File_common_meshfield_set_cartesbuf_field2D,           &
      File_common_meshfield_set_cartesbuf_field2D_cubedsphere
  
    implicit none
  
    class(FILE_base_meshfield), intent(inout) :: this
    integer, intent(in) :: dim_typeid
    character(*), intent(in) :: varname
    class(MeshField2D), intent(inout) :: field2d
    integer, intent(in), optional :: step
    logical, intent(in), optional :: allow_missing
  
    real(RP), allocatable :: buf(:,:)
    integer :: dims(2)
    integer :: start(2)   ! start offset of globale variable
    !-------------------------------------------------
  
    if ( this%fid /= -1 ) then
      start(:) = 1
      dims(1) = this%dimsinfo(MF2D_DIMTYPE_X)%size
      dims(2) = this%dimsinfo(MF2D_DIMTYPE_Y)%size
      allocate( buf(dims(1),dims(2)) )
  
      call FILE_Read( this%fid, varname,                       & ! (in)
        buf(:,:),                                              & ! (out)
        step=step, allow_missing=allow_missing                 ) ! (in)
  
      if ( associated( this%meshCS2D) ) then
        call File_common_meshfield_set_cartesbuf_field2D_cubedsphere( &
          this%meshCS2D, buf(:,:),                                    &
          field2d )
      else if ( associated( this%mesh2D) ) then
        call File_common_meshfield_set_cartesbuf_field2D( this%mesh2D, buf(:,:),   &
          field2d )
      end if
    end if
  
    return
  end subroutine FILE_base_meshfield_read_var2d

!OCL_SERIAL
  subroutine FILE_base_meshfield_read_var2d_local( this, & ! (inout)
    dim_typeid, varname, lcmesh, i0_s, j0_s,             & ! (in)
    val,                                                 & ! (out)
    step, allow_missing  )                                 ! (in)
  
    use scale_file, only: &
      FILE_Read
    use scale_file_common_meshfield, only: &
      File_common_meshfield_set_cartesbuf_field2D_local
  
    implicit none
  
    class(FILE_base_meshfield), intent(inout) :: this
    integer, intent(in) :: dim_typeid
    character(*), intent(in) :: varname
    class(LocalMesh2D), intent(in) :: lcmesh
    integer, intent(in) :: i0_s, j0_s
    real(RP), intent(out) :: val(lcmesh%refElem2D%Np,lcmesh%NeA)
    integer, intent(in), optional :: step
    logical, intent(in), optional :: allow_missing
  
    real(RP), allocatable :: buf(:,:)
    integer :: dims(2)
    integer :: start(2)   ! start offset of globale variable
    !-------------------------------------------------
  
    if ( this%fid /= -1 ) then
      start(:) = 1
      dims(1) = this%dimsinfo(MF2D_DIMTYPE_X)%size
      dims(2) = this%dimsinfo(MF2D_DIMTYPE_Y)%size
      allocate( buf(dims(1),dims(2)) )
  
      call FILE_Read( this%fid, varname,                       & ! (in)
        buf(:,:),                                              & ! (out)
        step=step, allow_missing=allow_missing                 ) ! (in)
  
      call File_common_meshfield_set_cartesbuf_field2D_local( &
        lcmesh, buf(:,:), i0_s, j0_s,                         &
        val(:,:) )
    end if
  
    return
  end subroutine FILE_base_meshfield_read_var2d_local

!OCL_SERIAL
  subroutine FILE_base_meshfield_read_var3d( this, & ! (inout)
    dim_typeid, varname,                           & ! (in)
    field3d,                                       & ! (inout)
    step, allow_missing                            ) ! (in)
  
    use scale_file, only: &
      FILE_Read
    use scale_file_common_meshfield, only: &
      File_common_meshfield_set_cartesbuf_field3D,            &
      File_common_meshfield_set_cartesbuf_field3D_cubedsphere

  
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
  
      if ( associated(this%meshCS3D) ) then
        call File_common_meshfield_set_cartesbuf_field3D_cubedsphere( &
          this%meshCS3D, buf(:,:,:),                                  &
          field3d )
      else if ( associated(this%mesh3D) ) then
        call File_common_meshfield_set_cartesbuf_field3D( this%mesh3D, buf(:,:,:), &
          field3d )
      end if
    end if
  
    return
  end subroutine FILE_base_meshfield_read_var3d

!OCL_SERIAL
  subroutine FILE_base_meshfield_read_var3d_local( this, & ! (inout)
    dim_typeid, varname, lcmesh, i0_s, j0_s, k0_s,       & ! (in)
    val,                                                 & ! (out)
    step, allow_missing  )                                 ! (in)
  
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

  subroutine FILE_base_meshfield_close( this ) ! (inout)
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

  
  subroutine FILE_base_meshfield_Final( this ) ! (inout)
    implicit none
    class(FILE_base_meshfield), intent(inout) :: this
    !--------------------------------------------------
  
    if ( allocated(this%vars_ncid) ) deallocate( this%vars_ncid )
    if ( allocated(this%dimsinfo) ) deallocate( this%dimsinfo )
    nullify( this%mesh1D, this%mesh2D, this%mesh3D )
  
    return
  end subroutine FILE_base_meshfield_Final

  subroutine FILE_base_meshfield_put_global_attribute_time( &
    this, date, subsec  )

    use scale_file, only: &
      FILE_Set_Attribute, &
      FILE_get_CFtunits
    use scale_calendar, only: &
      CALENDAR_get_name
       
    implicit none

    class(FILE_base_meshfield), intent(inout) :: this 
    integer, intent(in) :: date(6)
    real(DP), intent(in) :: subsec

    character(34) :: tunits
    character(len=H_SHORT) :: calendar_name
    !------------------------------------

    call FILE_Set_Attribute( this%fid, "global", "Conventions", "CF-1.6" ) ! [IN]
    call FILE_Set_Attribute( this%fid, "global", "grid_name", "hoge"     ) ! [IN]

    !- time

    if ( date(1) > 0 ) then
      call FILE_get_CFtunits( date(:), tunits )
      call CALENDAR_get_name( calendar_name )
    else
      tunits        = 'seconds'
      calendar_name = ''
    endif
        
    if ( calendar_name /= "" ) &
      call FILE_Set_Attribute( this%fid, "global", "calendar", calendar_name )
    call FILE_Set_Attribute( this%fid, "global", "time_units", tunits )
    call FILE_Set_Attribute( this%fid, "global", "time_start", (/ subsec /) )

    return
  end subroutine FILE_base_meshfield_put_global_attribute_time

  !- private -----------------------------------------

  subroutine def_axes( this, & ! (in)
    dtype                    ) ! (in)
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

    if (      associated(this%mesh2D)   &
         .or. associated(this%meshCS2D) ) then
      do d=1, 2
        call FILE_Def_Axis( this%fid, &
          this%dimsinfo(d)%name, this%dimsinfo(d)%desc, this%dimsinfo(d)%unit, &
          this%dimsinfo(d)%name, i_dtype, this%dimsinfo(d)%size                )
      end do
    end if

    if (      associated(this%mesh3D)   &
         .or. associated(this%meshCS3D) ) then
      do d=1, 3
        call FILE_Def_Axis( this%fid, &
          this%dimsinfo(d)%name, this%dimsinfo(d)%desc, this%dimsinfo(d)%unit, &
          this%dimsinfo(d)%name, i_dtype, this%dimsinfo(d)%size                )
      end do
    end if

    return
  end subroutine def_axes

  subroutine write_axes( this, & ! (in)
    start                      ) ! (in)
    use scale_const, only: &
      UNDEF => CONST_UNDEF
    use scale_file, only: &
      FILE_Write_Axis, &
      FILE_Set_Attribute
    use scale_file_common_meshfield, only: &
      File_common_meshfield_get_axis
    implicit none

    class(FILE_base_meshfield), intent(in) :: this
    integer, intent(in) :: start(3)

    real(RP), allocatable :: x(:)
    real(RP), allocatable :: y(:)
    real(RP), allocatable :: z(:)
    !------------

    if ( associated(this%mesh1D) ) then
      allocate( x(this%dimsinfo(1)%size) )
      call File_common_meshfield_get_axis( this%mesh1D, this%dimsinfo, x(:), this%force_uniform_grid )

      call FILE_Write_Axis( this%fid, this%dimsinfo(1)%name, x(:), start(1:1) )
    end if

    if ( associated(this%mesh2D)  ) then
      allocate( x(this%dimsinfo(1)%size), y(this%dimsinfo(2)%size) )
      call File_common_meshfield_get_axis( this%mesh2D, this%dimsinfo, x(:), y(:), this%force_uniform_grid )

      call FILE_Write_Axis( this%fid, this%dimsinfo(1)%name, x(:), start(1:1) )
      call FILE_Write_Axis( this%fid, this%dimsinfo(2)%name, y(:), start(2:2) )
    end if

    if ( associated(this%meshCS2D)  ) then
      allocate( x(this%dimsinfo(1)%size), y(this%dimsinfo(2)%size) )
      call File_common_meshfield_get_axis( this%meshCS2D, this%dimsinfo, x(:), y(:) )

      call FILE_Write_Axis( this%fid, this%dimsinfo(1)%name, x(:), start(1:1) )
      call FILE_Write_Axis( this%fid, this%dimsinfo(2)%name, y(:), start(2:2) )
    end if

    if ( associated(this%mesh3D) ) then
      allocate( x(this%dimsinfo(1)%size), y(this%dimsinfo(2)%size), z(this%dimsinfo(3)%size) )
      call File_common_meshfield_get_axis( this%mesh3D, this%dimsinfo, x(:), y(:), z(:), this%force_uniform_grid )

      call FILE_Write_Axis( this%fid, this%dimsinfo(1)%name, x(:), start(1:1) )
      call FILE_Write_Axis( this%fid, this%dimsinfo(2)%name, y(:), start(2:2) )
      call FILE_Write_Axis( this%fid, this%dimsinfo(3)%name, z(:), start(3:3) )
      if ( this%dimsinfo(3)%positive_down(1) ) &
        call FILE_Set_Attribute( this%fid, this%dimsinfo(3)%name, "positive", "down" )
    end if  

    if ( associated(this%meshCS3D) ) then
      allocate( x(this%dimsinfo(1)%size), y(this%dimsinfo(2)%size), z(this%dimsinfo(3)%size) )
      call File_common_meshfield_get_axis( this%meshCS3D, this%dimsinfo, x(:), y(:), z(:) )

      call FILE_Write_Axis( this%fid, this%dimsinfo(1)%name, x(:), start(1:1) )
      call FILE_Write_Axis( this%fid, this%dimsinfo(2)%name, y(:), start(2:2) )
      call FILE_Write_Axis( this%fid, this%dimsinfo(3)%name, z(:), start(3:3) )
      if ( this%dimsinfo(3)%positive_down(1) ) &
        call FILE_Set_Attribute( this%fid, this%dimsinfo(3)%name, "positive", "down" )
    end if  

    return
  end subroutine write_axes

end module scale_file_base_meshfield
