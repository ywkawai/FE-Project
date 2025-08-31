#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_mesh_base1d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase1D
  use scale_mesh_base1d, only: MeshBase1D
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CMeshBase1DPtr
    class(MeshBase1D), pointer :: obj
  end type CMeshBase1DPtr
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
!******
DEF_C_BIND_RELEASE_HANDLE(CMeshBase1DPtr,"CMeshBase1D_release_handle",)
!***** Getter
DEF_C_BIND_GETTER_I(CMeshBase1DPtr,get_NeG,"CMeshBase1D_get_NeG",NeG)

!******
!******
end module scale_mesh_base1d_cbind