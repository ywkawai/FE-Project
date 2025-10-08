#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_meshfieldcomm_base_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_io
  use scale_cbinding_util, only: CBindingBase
  use scale_meshfield_base_cbind, only: &
    CMeshField1D, CMeshField2D, CMeshField3D, &
    CMeshField1DPtr, CMeshField2DPtr, CMeshField3DPtr
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CMeshFieldContainer
    type(CMeshField1D), pointer :: obj_1D
    type(CMeshField2D), pointer :: obj_2D
    type(CMeshField3D), pointer :: obj_3D
    character(len=H_SHORT) :: tag
  end type CMeshFieldContainer

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CMeshFieldContainer_Init() result(ptr) bind(C, name="CMeshFieldContainer_Init")
    implicit none
    type(c_ptr) :: ptr

    type(CMeshFieldContainer), pointer :: handle
    !------------------------------------
    allocate(handle)
    ptr = c_loc(handle)
    handle%tag = "test"
    return
  end function CMeshFieldContainer_Init

  subroutine CMeshFieldContainer_SetField1D( ptr, field1D_ptr ) bind(C, name="CMeshFieldContainer_SetField1D")
    implicit none
    type(c_ptr), value :: ptr
    type(c_ptr), value :: field1D_ptr

    type(CMeshFieldContainer), pointer :: handle
    type(CMeshField1D), pointer :: field1D_fptr
    !------------------------------------
    call c_f_pointer( ptr, handle )
    call c_f_pointer( field1D_ptr, field1D_fptr )
    handle%obj_1D => field1D_fptr
    return
  end subroutine CMeshFieldContainer_SetField1D

  subroutine CMeshFieldContainer_SetField2D( ptr, field2D_ptr ) bind(C, name="CMeshFieldContainer_SetField2D")
    implicit none
    type(c_ptr), value :: ptr
    type(c_ptr), value :: field2D_ptr

    type(CMeshFieldContainer), pointer :: handle
    type(CMeshField2D), pointer :: mesh2D_fptr
    !------------------------------------
    call c_f_pointer( ptr, handle )
    call c_f_pointer( field2D_ptr, mesh2D_fptr )
    handle%obj_2D => mesh2D_fptr
    return
  end subroutine CMeshFieldContainer_SetField2D

  subroutine CMeshFieldContainer_SetField3D( ptr, field3D_ptr ) bind(C, name="CMeshFieldContainer_SetField3D")
    implicit none
    type(c_ptr), value :: ptr
    type(c_ptr), value :: field3D_ptr

    type(CMeshFieldContainer), pointer :: handle
    type(CMeshField3D), pointer :: mesh3D_fptr
    !------------------------------------
    call c_f_pointer( ptr, handle )
    call c_f_pointer( field3D_ptr, mesh3D_fptr )
    handle%obj_3D => mesh3D_fptr
    return
  end subroutine CMeshFieldContainer_SetField3D

!*****
DEF_C_BIND_RELEASE_HANDLE(CMeshFieldContainer,"CMeshFieldContainer_release_handle",)
!***** Getter
!******

end module scale_meshfieldcomm_base_cbind