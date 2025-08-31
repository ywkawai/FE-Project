#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_localmeshfield_base_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_localmeshfield_base, only: &
    LocalMeshField1D, LocalMeshField2D, LocalMeshField3D
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CLocalMeshField1D
    type(LocalMeshField1D) :: obj
  end type CLocalMeshField1D
  type, public, extends(CBindingBase) :: CLocalMeshField1DPtr
    type(LocalMeshField1D), pointer :: obj
  end type CLocalMeshField1DPtr  

  type, public, extends(CBindingBase) :: CLocalMeshField2D
    type(LocalMeshField2D) :: obj
  end type CLocalMeshField2D
  type, public, extends(CBindingBase) :: CLocalMeshField2DPtr
    type(LocalMeshField2D), pointer :: obj
  end type CLocalMeshField2DPtr  

  type, public, extends(CBindingBase) :: CLocalMeshField3D
    type(LocalMeshField3D) :: obj
  end type CLocalMeshField3D
  type, public, extends(CBindingBase) :: CLocalMeshField3DPtr
    type(LocalMeshField3D), pointer :: obj
  end type CLocalMeshField3DPtr  

  type, bind(C) :: f2dview
    type(c_ptr)       :: data
    integer(c_size_t) :: n1, n2
    integer(c_size_t) :: s0, s1
  end type f2dview

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  subroutine CLocalMeshField1D_get_val_view(ptr, view) bind(C, name="CLocalMeshField1D_get_val_view")
    implicit none
    type(c_ptr), value :: ptr
    type(f2dview),  intent(out) :: view
    type(CLocalMeshField1DPtr), pointer :: ptr_f
    !-----------------------------------------
    call c_f_pointer(ptr, ptr_f)
    call get_view( ptr_f%obj%val, view)
    return
  end subroutine CLocalMeshField1D_get_val_view   

  subroutine CLocalMeshField2D_get_val_view(ptr, view) bind(C, name="CLocalMeshField2D_get_val_view")
    implicit none
    type(c_ptr), value :: ptr
    type(f2dview),  intent(out) :: view
    type(CLocalMeshField2DPtr), pointer :: ptr_f
    !-----------------------------------------
    call c_f_pointer(ptr, ptr_f)
    call get_view( ptr_f%obj%val, view)
    return
  end subroutine CLocalMeshField2D_get_val_view   

  subroutine CLocalMeshField3D_get_val_view(ptr, view) bind(C, name="CLocalMeshField3D_get_val_view")
    implicit none
    type(c_ptr), value :: ptr
    type(f2dview),  intent(out) :: view
    type(CLocalMeshField3DPtr), pointer :: ptr_f
    !-----------------------------------------
    call c_f_pointer(ptr, ptr_f)
    call get_view( ptr_f%obj%val, view)
    return
  end subroutine CLocalMeshField3D_get_val_view   

  subroutine get_view(array2d, view)
    implicit none
    real(RP), intent(in), target :: array2d(:,:)

    real(RP), pointer :: alias(:,:)
    type(f2dview), intent(out) :: view
    !-----------------------------------------
    alias => array2d

    view%data  = c_loc(alias(1,1))
    view%n1 = size(alias,1, kind=c_size_t)
    view%n2 = size(alias,2, kind=c_size_t)
    view%s0 = 1_c_size_t
    view%s1 = size(alias,1, kind=c_size_t)
    return
  end subroutine get_view

!*****
DEF_C_BIND_RELEASE_HANDLE(CLocalMeshField1DPtr,"CLocalMeshField1D_release_handle",_1D)
DEF_C_BIND_RELEASE_HANDLE(CLocalMeshField2DPtr,"CLocalMeshField2D_release_handle",_2D)
DEF_C_BIND_RELEASE_HANDLE(CLocalMeshField3DPtr,"CLocalMeshField3D_release_handle",_3D)
!***** Getter
!******

end module scale_localmeshfield_base_cbind