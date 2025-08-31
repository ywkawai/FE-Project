#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_mesh_base3d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase3D
  use scale_mesh_base3d, only: MeshBase3D
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CMeshBase3DPtr
    class(MeshBase3D), pointer :: obj
  end type CMeshBase3DPtr
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
!******
DEF_C_BIND_RELEASE_HANDLE(CMeshBase3DPtr,"CMeshBase3D_release_handle",)
!***** Getter
!******
!******
end module scale_mesh_base3d_cbind