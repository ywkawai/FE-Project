#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_element_line_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_line, only: LineElement
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CLineElement
    type(LineElement) :: obj
  end type CLineElement
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CLineElement_Init(elemOrder, LumpedMassMatFlag ) result(ptr) bind(C, name="CLineElement_Init")
    implicit none
    type(c_ptr) :: ptr
    integer(c_int), intent(in), value :: elemOrder
    logical(c_bool), intent(in), value :: LumpedMassMatFlag
    type(CLineElement), pointer :: handle
    !------------------------------------
    call create_handle(handle, ptr)
    call handle%obj%Init( elemOrder, logical_c2f(LumpedMassMatFlag) )
    return
  end function CLineElement_Init

  subroutine CLineElement_Final( ptr ) bind(C, name="CLineElement_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle(ptr)
    return
  end subroutine CLineElement_Final

!******
DEF_C_BIND(CLineElement,"CLineElement",)
!***** Getter
DEF_C_BIND_GETTER_I(CLineElement,get_Np,"CLineElement_get_Np",Np)
DEF_C_BIND_GETTER_I(CLineElement,get_Nfaces,"CLineElement_get_Nfaces",Nfaces)
DEF_C_BIND_GETTER_I(CLineElement,get_NfpTot,"CLineElement_get_NfpTot",NfpTot)
DEF_C_BIND_GETTER_I(CLineElement,get_Nv,"CLineElement_get_Nv",Nv)
DEF_C_BIND_GETTER_I(CLineElement,get_PolyOrder,"CLineElement_get_PolyOrder",PolyOrder)
DEF_C_BIND_GETTER_Array1D_RP(CLineElement,get_x1,"CLineElement_get_x1",x1)
DEF_C_BIND_GETTER_Array1D_RP(CLineElement,get_IntWeight_lgl,"CLineElement_get_IntWeight_lgl",IntWeight_lgl)
DEF_C_BIND_GETTER_Array2D_RP(CLineElement,get_V,"CLineElement_get_V",V)
DEF_C_BIND_GETTER_Array2D_RP(CLineElement,get_invV,"CLineElement_get_invV",invV)
DEF_C_BIND_GETTER_Array2D_RP(CLineElement,get_M,"CLineElement_get_M",M)
DEF_C_BIND_GETTER_Array2D_RP(CLineElement,get_invM,"CLineElement_get_invM",invM)
DEF_C_BIND_GETTER_Array2D_RP(CLineElement,get_Lift,"CLineElement_get_Lift",Lift)
DEF_C_BIND_GETTER_Array2D_RP(CLineElement,get_Dx1,"CLineElement_get_Dx1",Dx1)
!******

end module scale_element_line_cbind