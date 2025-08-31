#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_element_base3d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase3D
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CElementBase3D
    type(ElementBase3D) :: obj
  end type CElementBase3D
  type, public, extends(CBindingBase) :: CElementBase3DPtr
    type(ElementBase3D), pointer :: obj
  end type CElementBase3DPtr  
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
!*****
DEF_C_BIND_RELEASE_HANDLE(CElementBase3DPtr,"CElementBase3D_release_handle",)
!***** Getter
DEF_C_BIND_GETTER_I(CElementBase3DPtr,get_Np,"CElementBase3D_get_Np",Np)
DEF_C_BIND_GETTER_I(CElementBase3DPtr,get_Nnode_h1D,"CElementBase3D_get_Nnode_h1D",Nnode_h1D)
DEF_C_BIND_GETTER_I(CElementBase3DPtr,get_Nnode_v,"CElementBase3D_get_Nnode_v",Nnode_v)
DEF_C_BIND_GETTER_I(CElementBase3DPtr,get_PolyOrder_h,"CElementBase3D_get_PolyOrder_h",PolyOrder_h)
DEF_C_BIND_GETTER_I(CElementBase3DPtr,get_PolyOrder_v,"CElementBase3D_get_PolyOrder_v",PolyOrder_v)
!******

end module scale_element_base3d_cbind