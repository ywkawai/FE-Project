#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_localmesh1d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase1D
  use scale_localmesh_1d, only: LocalMesh1D
  use iso_c_binding
  use scale_element_base1d_cbind, only: CElementBase1DPtr
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CLocalMesh1D
    type(LocalMesh1D) :: obj
  end type CLocalMesh1D
  type, public, extends(CBindingBase) :: CLocalMesh1DPtr
    type(LocalMesh1D), pointer :: obj
  end type CLocalMesh1DPtr
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
!******
DEF_C_BIND_RELEASE_HANDLE(CLocalMesh1DPtr,"CLocalMesh1D_release_handle")
!***** Getter
DEF_C_BIND_GETTER_I(CLocalMesh1DPtr,get_Ne,"CLocalMesh1D_get_Ne",Ne)
DEF_C_BIND_GETTER_I(CLocalMesh1DPtr,get_NeS,"CLocalMesh1D_get_NeS",NeS)
DEF_C_BIND_GETTER_I(CLocalMesh1DPtr,get_NeE,"CLocalMesh1D_get_NeE",NeE)
DEF_C_BIND_GETTER_I(CLocalMesh1DPtr,get_NeA,"CLocalMesh1D_get_NeA",NeA)
DEF_C_BIND_GETTER_Array3D_RP(CLocalMesh1DPtr,get_pos_en,"CLocalMesh1D_get_pos_en",pos_en)
!******
DEF_C_BIND_MAKE_FOBJ_HANDLE(CLocalMesh1DPtr,get_refElem1D,"CLocalMesh1D_get_refElem1D",ElementBase1D,CElementBase1DPtr,refElem1D)
!******
end module scale_localmesh1d_cbind