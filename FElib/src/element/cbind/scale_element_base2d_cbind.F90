#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_element_base2d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase2D
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CElementBase2D
    type(ElementBase2D) :: obj
  end type CElementBase2D
  type, public, extends(CBindingBase) :: CElementBase2DPtr
    type(ElementBase2D), pointer :: obj
  end type CElementBase2DPtr  
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
!*****
DEF_C_BIND_RELEASE_HANDLE(CElementBase2DPtr,"CElementBase2D_release_handle")
!***** Getter
DEF_C_BIND_GETTER_I(CElementBase2DPtr,get_Np,"CElementBase2D_get_Np",Np)
DEF_C_BIND_GETTER_I(CElementBase2DPtr,get_PolyOrder,"CElementBase2D_get_PolyOrder",PolyOrder)
!******

end module scale_element_base2d_cbind