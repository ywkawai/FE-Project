#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_element_quadrilateral_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_quadrilateral, only: QuadrilateralElement
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CQuadrilateralElement
    type(QuadrilateralElement) :: obj
  end type CQuadrilateralElement
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CQuadrilateralElement_Init(elemOrder, LumpedMassMatFlag ) result(ptr) bind(C, name="CQuadrilateralElement_Init")
    implicit none
    type(c_ptr) :: ptr
    integer(c_int), intent(in), value :: elemOrder
    logical(c_bool), intent(in), value :: LumpedMassMatFlag
    type(CQuadrilateralElement), pointer :: handle
    !------------------------------------
    call create_handle(handle, ptr)
    call handle%obj%Init( elemOrder, logical_c2f(LumpedMassMatFlag) )
    return
  end function CQuadrilateralElement_Init

  subroutine CQuadrilateralElement_Final( ptr ) bind(C, name="CQuadrilateralElement_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle(ptr)
    return
  end subroutine CQuadrilateralElement_Final

!******
DEF_C_BIND(CQuadrilateralElement,"CQuadrilateralElement",)
!***** Getter
DEF_C_BIND_GETTER_I(CQuadrilateralElement,get_Np,"CQuadrilateralElement_get_Np",Np)
DEF_C_BIND_GETTER_I(CQuadrilateralElement,get_Nfaces,"CQuadrilateralElement_get_Nfaces",Nfaces)
DEF_C_BIND_GETTER_I(CQuadrilateralElement,get_Nfp,"CQuadrilateralElement_get_Nfp",Nfp)
DEF_C_BIND_GETTER_I(CQuadrilateralElement,get_NfpTot,"CQuadrilateralElement_get_NfpTot",NfpTot)
DEF_C_BIND_GETTER_I(CQuadrilateralElement,get_Nv,"CQuadrilateralElement_get_Nv",Nv)
DEF_C_BIND_GETTER_I(CQuadrilateralElement,get_PolyOrder,"CQuadrilateralElement_get_PolyOrder",PolyOrder)
DEF_C_BIND_GETTER_Array1D_RP(CQuadrilateralElement,get_x1,"CQuadrilateralElement_get_x1",x1)
DEF_C_BIND_GETTER_Array1D_RP(CQuadrilateralElement,get_x2,"CQuadrilateralElement_get_x2",x2)
DEF_C_BIND_GETTER_Array1D_RP(CQuadrilateralElement,get_IntWeight_lgl,"CQuadrilateralElement_get_IntWeight_lgl",IntWeight_lgl)
DEF_C_BIND_GETTER_Array2D_RP(CQuadrilateralElement,get_V,"CQuadrilateralElement_get_V",V)
DEF_C_BIND_GETTER_Array2D_RP(CQuadrilateralElement,get_invV,"CQuadrilateralElement_get_invV",invV)
DEF_C_BIND_GETTER_Array2D_RP(CQuadrilateralElement,get_M,"CQuadrilateralElement_get_M",M)
DEF_C_BIND_GETTER_Array2D_RP(CQuadrilateralElement,get_invM,"CQuadrilateralElement_get_invM",invM)
DEF_C_BIND_GETTER_Array2D_RP(CQuadrilateralElement,get_Lift,"CQuadrilateralElement_get_Lift",Lift)
DEF_C_BIND_GETTER_Array2D_RP(CQuadrilateralElement,get_Dx1,"CQuadrilateralElement_get_Dx1",Dx1)
DEF_C_BIND_GETTER_Array2D_RP(CQuadrilateralElement,get_Dx2,"CQuadrilateralElement_get_Dx2",Dx2)
!******

end module scale_element_quadrilateral_cbind