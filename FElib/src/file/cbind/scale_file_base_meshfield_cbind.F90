#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_file_base_meshfield_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_io
  use scale_prc

  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f, &
    logical_f2c => cbinding_util_logical_f2c, &
    string_c2f => cbinding_util_string_c2f
  use scale_file_base_meshfield, only: FILE_base_meshfield

  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CFILE_base_meshfield
    type(FILE_base_meshfield) :: obj
  end type CFILE_base_meshfield

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CFILE_base_meshfield_Init( var_num, mesh1D, mesh2D, meshCubedSphere2D, mesh3D, meshCubedSphere3D, &
    force_uniform_grid ) result(ptr) bind(C, name="CFILE_base_meshfield_Init")
    use scale_mesh_base1d_cbind, only: CMeshBase1DPtr
    use scale_mesh_rectdom2d_cbind, only: CMeshRectDom2D
    use scale_mesh_cubedom3d_cbind, only: CMeshCubeDom3D
    implicit none
    integer(c_int), value :: var_num
    type(c_ptr), value :: mesh1D
    type(c_ptr), value :: mesh2D
    type(c_ptr), value :: meshCubedSphere2D
    type(c_ptr), value :: mesh3D
    type(c_ptr), value :: meshCubedSphere3D
    logical(c_bool), intent(in), value :: force_uniform_grid
    type(c_ptr) :: ptr

    type(CFILE_base_meshfield), pointer :: handle

    type(CMeshBase1DPtr), pointer :: mesh1D_fptr
    type(CMeshRectDom2D), pointer :: mesh2D_fptr
    type(CMeshCubeDom3D), pointer :: mesh3D_fptr
    character(:), allocatable :: varname_f, units_f
    !------------------------------------
    
    call create_handle(handle, ptr)

    if ( c_associated(mesh1D) ) then
      call c_f_pointer( mesh1D, mesh1D_fptr )
      call handle%obj%Init( var_num, mesh1D=mesh1D_fptr%obj, force_uniform_grid=logical_c2f(force_uniform_grid) )
    else if (c_associated(mesh2D) ) then
      call c_f_pointer( mesh2D, mesh2D_fptr )
      call handle%obj%Init( var_num, mesh2D=mesh2D_fptr%obj, force_uniform_grid=logical_c2f(force_uniform_grid) )
    else if (c_associated(mesh3D) ) then
      call c_f_pointer( mesh3D, mesh3D_fptr )
      call handle%obj%Init( var_num, mesh3D=mesh3D_fptr%obj, force_uniform_grid=logical_c2f(force_uniform_grid) )
    else
      LOG_INFO("CFILE_base_meshfield_Init",*) "Unsupported mesh type is specified. Check!"
      call PRC_abort
    end if
    return
  end function CFILE_base_meshfield_Init

  subroutine CFILE_base_meshfield_Final( ptr ) bind(C, name="CFILE_base_meshfield_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle(ptr)
    return
  end subroutine CFILE_base_meshfield_Final

  subroutine CFILE_base_meshfield_Open( ptr, basename, myrank ) bind(C, name="CFILE_base_meshfield_Open")
    implicit none
    type(c_ptr), value :: ptr
    type(c_ptr), value :: basename
    integer(c_int), value :: myrank

    type(CFILE_base_meshfield), pointer :: handle

    character(:), allocatable :: basename_f
    integer :: myrank_
    !------------------------------------

    call c_f_pointer(ptr, handle)

    if ( myrank < 0 ) then
      myrank_ = PRC_myrank
    else
      myrank_ = myrank
    end if
    basename_f = string_c2f(basename)

    call handle%obj%Open( basename_f, myrank_ )
    return
  end subroutine CFILE_base_meshfield_Open

  subroutine CFILE_base_meshfield_Create( ptr, basename, title, dtype, &
    fileexisted,             &
    myrank, tunits, calendar ) bind(C, name="CFILE_base_meshfield_Create")
    implicit none
    type(c_ptr), value :: ptr
    type(c_ptr), value :: basename
    type(c_ptr), value :: title
    type(c_ptr), value :: dtype
    logical(c_bool) :: fileexisted
    integer(c_int), value :: myrank
    type(c_ptr), value :: tunits
    type(c_ptr), value :: calendar

    type(CFILE_base_meshfield), pointer :: handle

    character(:), allocatable :: basename_f, title_f, dtype_f
    character(:), allocatable :: tunits_f, calendar_f
    integer :: myrank_
    logical :: fileexisted_f
    !------------------------------------

    call c_f_pointer(ptr, handle)

    if ( myrank < 0 ) then
      myrank_ = PRC_myrank
    else
      myrank_ = myrank
    end if
    basename_f = string_c2f(basename)
    title_f = string_c2f(title)
    dtype_f = string_c2f(dtype)
    tunits_f = string_c2f(tunits)
    calendar_f = string_c2f(calendar)

    call handle%obj%Create( basename_f, title_f, dtype_f, & ! (in)
        fileexisted_f,                & ! (out)
        myrank_, tunits_f, calendar_f ) ! (in)
        
    fileexisted = logical_f2c(fileexisted_f)
    return
  end subroutine CFILE_base_meshfield_Create

  subroutine CFILE_base_meshfield_def_var1( ptr, field, desc, vid, dim_type_id, datatype, &
    standard_name, timeinv, nsteps ) bind(C, name="CFILE_base_meshfield_Def_var1")
    use scale_meshfield_base_cbind, only: CMeshFieldBasePtr
    implicit none
    type(c_ptr), value :: ptr
    type(c_ptr), value :: field
    type(c_ptr), value :: desc
    integer(c_int), value :: vid
    integer(c_int), value :: dim_type_id
    type(c_ptr), value :: datatype
    type(c_ptr), value :: standard_name
    real(c_double), value :: timeinv
    integer(c_int), value :: nsteps

    type(CFILE_base_meshfield), pointer :: handle
    type(CMeshFieldBasePtr), pointer :: field_fptr
    character(:), allocatable :: desc_f, datatype_f
    character(:), allocatable :: standard_name_f
    logical :: fileexisted_f
    !------------------------------------

    call c_f_pointer(ptr, handle)

    call c_f_pointer(field, field_fptr)
    desc_f = string_c2f(desc)
    datatype_f = string_c2f(datatype)
    standard_name_f = string_c2f(standard_name)

    call handle%obj%Def_Var( field_fptr%obj, desc_f, vid, dim_type_id, datatype_f, &
        standard_name=standard_name_f, timeinv=timeinv, nsteps=nsteps)
    return
  end subroutine CFILE_base_meshfield_def_var1

  subroutine CFILE_base_meshfield_write_var1d( ptr, vid, field1d, sec_str, sec_end &
    ) bind(C, name="CFILE_base_meshfield_Write_var1d")
    use scale_meshfield_base_cbind, only: CMeshField1D
    implicit none
    type(c_ptr), value :: ptr
    integer(c_int), value :: vid
    type(c_ptr), value :: field1d
    real(c_double), value :: sec_str
    real(c_double), value :: sec_end

    type(CFILE_base_meshfield), pointer :: handle
    type(CMeshField1D), pointer :: field1d_fptr
    !------------------------------------

    call c_f_pointer(ptr, handle)

    call c_f_pointer(field1d, field1d_fptr)
    call handle%obj%Write_var1D(vid, field1d_fptr%obj, sec_str, sec_end)
    return
  end subroutine CFILE_base_meshfield_write_var1d

  subroutine CFILE_base_meshfield_write_var2d( ptr, vid, field2d, sec_str, sec_end &
    ) bind(C, name="CFILE_base_meshfield_Write_var2d")
    use scale_meshfield_base_cbind, only: CMeshField2D
    implicit none
    type(c_ptr), value :: ptr
    integer(c_int), value :: vid
    type(c_ptr), value :: field2d
    real(c_double), value :: sec_str
    real(c_double), value :: sec_end

    type(CFILE_base_meshfield), pointer :: handle
    type(CMeshField2D), pointer :: field2d_fptr
    !------------------------------------

    call c_f_pointer(ptr, handle)

    call c_f_pointer(field2d, field2d_fptr)
    call handle%obj%Write_var2D(vid, field2d_fptr%obj, sec_str, sec_end)
    return
  end subroutine CFILE_base_meshfield_write_var2d

  subroutine CFILE_base_meshfield_write_var3d( ptr, vid, field3d, sec_str, sec_end &
    ) bind(C, name="CFILE_base_meshfield_Write_var3d")
    use scale_meshfield_base_cbind, only: CMeshField3D
    implicit none
    type(c_ptr), value :: ptr
    integer(c_int), value :: vid
    type(c_ptr), value :: field3d
    real(c_double), value :: sec_str
    real(c_double), value :: sec_end

    type(CFILE_base_meshfield), pointer :: handle
    type(CMeshField3D), pointer :: field3d_fptr
    !------------------------------------

    call c_f_pointer(ptr, handle)

    call c_f_pointer(field3d, field3d_fptr)
    call handle%obj%Write_var3D(vid, field3d_fptr%obj, sec_str, sec_end)
    return
  end subroutine CFILE_base_meshfield_write_var3d

!*****
DEF_C_BIND(CFILE_base_meshfield,"CFILE_base_meshfield",)
DEF_C_BIND_SUB(CFILE_base_meshfield,End_def,"CFILE_base_meshfield_End_def")
DEF_C_BIND_SUB(CFILE_base_meshfield,Close,"CFILE_base_meshfield_Close")
!***** Getter
!******

end module scale_file_base_meshfield_cbind