#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_element_base1d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase1D
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CElementBase1D
    type(ElementBase1D) :: obj
  end type CElementBase1D
  type, public, extends(CBindingBase) :: CElementBase1DPtr
    type(ElementBase1D), pointer :: obj
  end type CElementBase1DPtr  
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
!*****
DEF_C_BIND_RELEASE_HANDLE(CElementBase1DPtr,"CElementBase1D_release_handle")
!***** Getter
DEF_C_BIND_GETTER_I(CElementBase1DPtr,get_Np,"CElementBase1D_get_Np",Np)
DEF_C_BIND_GETTER_I(CElementBase1DPtr,get_PolyOrder,"CElementBase1D_get_PolyOrder",PolyOrder)
!******

end module scale_element_base1d_cbind