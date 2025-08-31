#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_localmesh3d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase3D
  use scale_localmesh_3d, only: LocalMesh3D
  use iso_c_binding
  use scale_element_base3d_cbind, only: CElementBase3DPtr
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CLocalMesh3D
    type(LocalMesh3D) :: obj
  end type CLocalMesh3D
  type, public, extends(CBindingBase) :: CLocalMesh3DPtr
    type(LocalMesh3D), pointer :: obj
  end type CLocalMesh3DPtr
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
!******
DEF_C_BIND_RELEASE_HANDLE(CLocalMesh3DPtr,"CLocalMesh3D_release_handle",)
!***** Getter
DEF_C_BIND_GETTER_I(CLocalMesh3DPtr,get_Ne,"CLocalMesh3D_get_Ne",Ne)
DEF_C_BIND_GETTER_I(CLocalMesh3DPtr,get_NeS,"CLocalMesh3D_get_NeS",NeS)
DEF_C_BIND_GETTER_I(CLocalMesh3DPtr,get_NeE,"CLocalMesh3D_get_NeE",NeE)
DEF_C_BIND_GETTER_I(CLocalMesh3DPtr,get_NeA,"CLocalMesh3D_get_NeA",NeA)
DEF_C_BIND_GETTER_Array3D_RP(CLocalMesh3DPtr,get_pos_en,"CLocalMesh3D_get_pos_en",pos_en)
!******
DEF_C_BIND_MAKE_FOBJ_HANDLE(CLocalMesh3DPtr,get_refElem3D,"CLocalMesh3D_get_refElem3D",CElementBase3DPtr,refElem3D)
!******
end module scale_localmesh3d_cbind