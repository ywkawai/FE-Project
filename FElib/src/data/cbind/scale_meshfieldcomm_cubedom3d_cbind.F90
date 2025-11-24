#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_meshfieldcomm_cubedom3d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f, &
    string_c2f => cbinding_util_string_c2f
  use scale_meshfieldcomm_cubedom3d, only: &
    MeshFieldCommCubeDom3D
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_base_cbind, only: CMeshFieldContainer

  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CMeshFieldCommCubeDom3D
    type(MeshFieldCommCubeDom3D) :: obj
  end type CMeshFieldCommCubeDom3D
  type, public, extends(CBindingBase) :: CMeshFieldCommCubeDom3DPtr
    type(MeshFieldCommCubeDom3D), pointer :: obj
  end type CMeshFieldCommCubeDom3DPtr  

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CMeshFieldCommCubeDom3D_Init( sfield_num, hvfield_num, htensorfield_num, mesh3d_ptr) result(ptr) bind(C, name="CMeshFieldCommCubeDom3D_Init")
    use scale_mesh_cubedom3d_cbind, only: CMeshCubeDom3D
    implicit none
    integer(c_int), value :: sfield_num
    integer(c_int), value :: hvfield_num
    integer(c_int), value :: htensorfield_num
    type(c_ptr), value :: mesh3D_ptr
    type(c_ptr) :: ptr

    type(CMeshFieldCommCubeDom3D), pointer :: handle
    type(CMeshCubeDom3D), pointer :: mesh3D_fptr
    !------------------------------------
    
    call create_handle(handle, ptr)
    call c_f_pointer( mesh3D_ptr, mesh3D_fptr )

    call handle%obj%Init(sfield_num, hvfield_num, htensorfield_num, mesh3D_fptr%obj)
    return
  end function CMeshFieldCommCubeDom3D_Init

  subroutine CMeshFieldCommCubeDom3D_Final( ptr ) bind(C, name="CMeshFieldCommCubeDom3D_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle(ptr)
    return
  end subroutine CMeshFieldCommCubeDom3D_Final

  subroutine CMeshFieldCommCubeDom3D_put( ptr, field_list_ptr, varnum, varid_s ) bind(C, name="CMeshFieldCommCubeDom3D_Put")
    use scale_meshfieldcomm_base, only: MeshFieldContainer
    implicit none
    type(c_ptr), value :: ptr
    type(c_ptr), value :: field_list_ptr
    integer(c_int), value :: varnum
    integer(c_int), value :: varid_s

    type(CMeshFieldCommCubeDom3D), pointer :: handle

    type(c_ptr), pointer :: field_list_tmp(:)
    type(CMeshFieldContainer), pointer :: item_fptr 
    type(MeshFieldContainer) :: field_list(varnum)
    integer :: i
    !------------------------------------

    call c_f_pointer(ptr, handle)

    call c_f_pointer( field_list_ptr, field_list_tmp, [varnum] )
    do i=1, varnum
      call c_f_pointer( field_list_tmp(i), item_fptr )
      field_list(i)%field1d => item_fptr%obj_1D%obj
    end do
    call handle%obj%Put( field_list, varid_s )
    return
  end subroutine CMeshFieldCommCubeDom3D_put

  subroutine CMeshFieldCommCubeDom3D_get( ptr, field_list_ptr, varnum, varid_s ) bind(C, name="CMeshFieldCommCubeDom3D_Get")
    implicit none
    type(c_ptr), value :: ptr
    type(c_ptr), value :: field_list_ptr
    integer(c_int), value :: varnum
    integer(c_int), value :: varid_s

    type(CMeshFieldCommCubeDom3D), pointer :: handle

    type(c_ptr), pointer :: field_list_tmp(:)
    type(CMeshFieldContainer), pointer :: item_fptr 
    type(MeshFieldContainer) :: field_list(varnum)
    integer :: i
    !------------------------------------

    call c_f_pointer(ptr, handle)

    call c_f_pointer( field_list_ptr, field_list_tmp, [varnum] )
    do i=1, varnum
      call c_f_pointer( field_list_tmp(i), item_fptr )
      field_list(i)%field1d => item_fptr%obj_1D%obj
    end do
    call handle%obj%Get( field_list, varid_s )
    return
  end subroutine CMeshFieldCommCubeDom3D_get

  subroutine CMeshFieldCommCubeDom3D_exchange( ptr, do_wait ) bind(C, name="CMeshFieldCommCubeDom3D_Exchange")
    implicit none
    type(c_ptr), value :: ptr
    logical(c_bool), value :: do_wait

    type(CMeshFieldCommCubeDom3D), pointer :: handle
    !------------------------------------

    call c_f_pointer(ptr, handle)
    call handle%obj%Exchange( logical_c2f(do_wait) )
    return
  end subroutine CMeshFieldCommCubeDom3D_exchange

!-- private -

!*****
DEF_C_BIND(CMeshFieldCommCubeDom3D,"CMeshFieldCommCubeDom3D",)
!***** Getter
!******

end module scale_meshfieldcomm_cubedom3d_cbind