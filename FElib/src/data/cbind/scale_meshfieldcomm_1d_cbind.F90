#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_meshfieldcomm_1d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f, &
    string_c2f => cbinding_util_string_c2f
  use scale_meshfieldcomm_1d, only: &
    MeshFieldComm1D
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
  type, public, extends(CBindingBase) :: CMeshFieldComm1D
    type(MeshFieldComm1D) :: obj
  end type CMeshFieldComm1D
  type, public, extends(CBindingBase) :: CMeshFieldComm1DPtr
    type(MeshFieldComm1D), pointer :: obj
  end type CMeshFieldComm1DPtr  

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CMeshFieldComm1D_Init( sfield_num, hvfield_num, mesh1D_ptr) result(ptr) bind(C, name="CMeshFieldComm1D_Init")
    use scale_mesh_linedom1d_cbind, only: CMeshLineDom1D
    implicit none
    integer(c_int), value :: sfield_num
    integer(c_int), value :: hvfield_num
    type(c_ptr), value :: mesh1D_ptr
    type(c_ptr) :: ptr

    type(CMeshFieldComm1D), pointer :: handle
    type(CMeshLineDom1D), pointer :: mesh1D_fptr
    !------------------------------------
    
    call create_handle(handle, ptr)
    call c_f_pointer( mesh1D_ptr, mesh1D_fptr )

    call handle%obj%Init(sfield_num, hvfield_num, mesh1D_fptr%obj)
    write(*,*) "0 LocalMesh_NUM", handle%obj%mesh%LOCAL_MESH_NUM
    return
  end function CMeshFieldComm1D_Init

  subroutine CMeshFieldComm1D_Final( ptr ) bind(C, name="CMeshFieldComm1D_Final")
    implicit none
    type(c_ptr), value :: ptr

    type(CMeshFieldComm1D), pointer :: handle
    !------------------------------------
    call c_f_pointer(ptr, handle)
    write(*,*) "1 LocalMesh_NUM", handle%obj%mesh%LOCAL_MESH_NUM
    write(*,*) "1 nfaces_comm", handle%obj%nfaces_comm
    call handle%obj%Final()
    deallocate( handle )
    ! call destroy_handle(ptr)
    return
  end subroutine CMeshFieldComm1D_Final

  subroutine CMeshFieldComm1D_put( ptr, field_list_ptr, varnum, varid_s ) bind(C, name="CMeshFieldComm1D_Put")
    use scale_meshfieldcomm_base, only: MeshFieldContainer
    implicit none
    type(c_ptr), value :: ptr
    type(c_ptr), value :: field_list_ptr
    integer(c_int), value :: varnum
    integer(c_int), value :: varid_s

    type(CMeshFieldComm1D), pointer :: handle

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
  end subroutine CMeshFieldComm1D_put

  subroutine CMeshFieldComm1D_get( ptr, field_list_ptr, varnum, varid_s ) bind(C, name="CMeshFieldComm1D_Get")
    implicit none
    type(c_ptr), value :: ptr
    type(c_ptr), value :: field_list_ptr
    integer(c_int), value :: varnum
    integer(c_int), value :: varid_s

    type(CMeshFieldComm1D), pointer :: handle

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
  end subroutine CMeshFieldComm1D_get

  subroutine CMeshFieldComm1D_exchange( ptr, do_wait ) bind(C, name="CMeshFieldComm1D_Exchange")
    implicit none
    type(c_ptr), value :: ptr
    logical(c_bool), value :: do_wait

    type(CMeshFieldComm1D), pointer :: handle
    !------------------------------------

    call c_f_pointer(ptr, handle)
    call handle%obj%Exchange( logical_c2f(do_wait) )
    return
  end subroutine CMeshFieldComm1D_exchange

!-- private --

!*****
DEF_C_BIND(CMeshFieldComm1D,"CMeshFieldComm1D",)
!***** Getter
!******

end module scale_meshfieldcomm_1d_cbind