#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_localmesh2d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase2D
  use scale_localmesh_2d, only: LocalMesh2D
  use iso_c_binding
  use scale_element_base2d_cbind, only: CElementBase2DPtr
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CLocalMesh2D
    type(LocalMesh2D) :: obj
  end type CLocalMesh2D
  type, public, extends(CBindingBase) :: CLocalMesh2DPtr
    type(LocalMesh2D), pointer :: obj
  end type CLocalMesh2DPtr
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
!******
DEF_C_BIND_RELEASE_HANDLE(CLocalMesh2DPtr,"CLocalMesh2D_release_handle")
!***** Getter
DEF_C_BIND_GETTER_I(CLocalMesh2DPtr,get_Ne,"CLocalMesh2D_get_Ne",Ne)
DEF_C_BIND_GETTER_I(CLocalMesh2DPtr,get_NeS,"CLocalMesh2D_get_NeS",NeS)
DEF_C_BIND_GETTER_I(CLocalMesh2DPtr,get_NeE,"CLocalMesh2D_get_NeE",NeE)
DEF_C_BIND_GETTER_I(CLocalMesh2DPtr,get_NeA,"CLocalMesh2D_get_NeA",NeA)
DEF_C_BIND_GETTER_Array3D_RP(CLocalMesh2DPtr,get_pos_en,"CLocalMesh2D_get_pos_en",pos_en)
!******
DEF_C_BIND_MAKE_FOBJ_HANDLE(CLocalMesh2DPtr,get_refElem2D,"CLocalMesh2D_get_refElem2D",ElementBase2D,CElementBase2DPtr,refElem2D)
!******
end module scale_localmesh2d_cbind