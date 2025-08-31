#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_element_hexahedral_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_hexahedral, only: HexahedralElement
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CHexahedralElement
    type(HexahedralElement) :: obj
  end type CHexahedralElement
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CHexahedralElement_Init(elemOrder_h, elemOrder_v, LumpedMassMatFlag ) result(ptr) bind(C, name="CHexahedralElement_Init")
    implicit none
    type(c_ptr) :: ptr
    integer(c_int), intent(in), value :: elemOrder_h
    integer(c_int), intent(in), value :: elemOrder_v
    logical(c_bool), intent(in), value :: LumpedMassMatFlag
    type(CHexahedralElement), pointer :: handle
    !------------------------------------
    call create_handle(handle, ptr)
    call handle%obj%Init( elemOrder_h, elemOrder_v, logical_c2f(LumpedMassMatFlag) )
    return
  end function CHexahedralElement_Init

  subroutine CHexahedralElement_Final( ptr ) bind(C, name="CHexahedralElement_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle(ptr)
    return
  end subroutine CHexahedralElement_Final

!******
DEF_C_BIND(CHexahedralElement,"CHexahedralElement",)
!***** Getter
DEF_C_BIND_GETTER_I(CHexahedralElement,get_Np,"CHexahedralElement_get_Np",Np)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_Nnode_h1D,"CHexahedralElement_get_Nnode_h1D",Nnode_h1D)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_Nnode_v,"CHexahedralElement_get_Nnode_v",Nnode_v)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_Nfaces_h,"CHexahedralElement_get_Nfaces_h",Nfaces_h)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_Nfaces_v,"CHexahedralElement_get_Nfaces_v",Nfaces_v)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_Nfaces,"CHexahedralElement_get_Nfaces",Nfaces)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_Nfp_h,"CHexahedralElement_get_Nfp_h",Nfp_h)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_Nfp_v,"CHexahedralElement_get_Nfp_v",Nfp_v)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_NfpTot,"CHexahedralElement_get_NfpTot",NfpTot)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_Nv,"CHexahedralElement_get_Nv",Nv)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_PolyOrder_h,"CHexahedralElement_get_PolyOrder_h",PolyOrder_h)
DEF_C_BIND_GETTER_I(CHexahedralElement,get_PolyOrder_v,"CHexahedralElement_get_PolyOrder_v",PolyOrder_v)
DEF_C_BIND_GETTER_Array1D_RP(CHexahedralElement,get_x1,"CHexahedralElement_get_x1",x1)
DEF_C_BIND_GETTER_Array1D_RP(CHexahedralElement,get_x2,"CHexahedralElement_get_x2",x2)
DEF_C_BIND_GETTER_Array1D_RP(CHexahedralElement,get_x3,"CHexahedralElement_get_x3",x3)
DEF_C_BIND_GETTER_Array1D_RP(CHexahedralElement,get_IntWeight_lgl,"CHexahedralElement_get_IntWeight_lgl",IntWeight_lgl)
DEF_C_BIND_GETTER_Array2D_RP(CHexahedralElement,get_V,"CHexahedralElement_get_V",V)
DEF_C_BIND_GETTER_Array2D_RP(CHexahedralElement,get_invV,"CHexahedralElement_get_invV",invV)
DEF_C_BIND_GETTER_Array2D_RP(CHexahedralElement,get_M,"CHexahedralElement_get_M",M)
DEF_C_BIND_GETTER_Array2D_RP(CHexahedralElement,get_invM,"CHexahedralElement_get_invM",invM)
DEF_C_BIND_GETTER_Array2D_RP(CHexahedralElement,get_Lift,"CHexahedralElement_get_Lift",Lift)
DEF_C_BIND_GETTER_Array2D_RP(CHexahedralElement,get_Dx1,"CHexahedralElement_get_Dx1",Dx1)
DEF_C_BIND_GETTER_Array2D_RP(CHexahedralElement,get_Dx2,"CHexahedralElement_get_Dx2",Dx2)
DEF_C_BIND_GETTER_Array2D_RP(CHexahedralElement,get_Dx3,"CHexahedralElement_get_Dx3",Dx3)
!******

end module scale_element_hexahedral_cbind