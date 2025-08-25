#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_mesh_rectdom2d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase1D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CMeshRectDom2D
    type(MeshRectDom2D) :: obj
  end type CMeshRectDom2D
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CMeshRectDom2D_Init( &
    NeGX, NeGY,                             & ! (in)
    dom_xmin, dom_xmax, dom_ymin, dom_ymax, & ! (in)
    is_PeriodicX, is_PeriodicY,             & ! (in)
    refElem_ptr, NLocalMeshPerPrc,          & ! (in)
    NprcX, NprcY, nproc, myrank )           & ! (in)
    result(ptr) bind(C, name="CMeshRectDom2D_Init")
    implicit none
    type(c_ptr) :: ptr
    integer(c_int), value :: NeGX
    integer(c_int), value :: NeGY
    real(c_double), value :: dom_xmin
    real(c_double), value :: dom_xmax   
    real(c_double), value :: dom_ymin
    real(c_double), value :: dom_ymax
    logical(c_bool), value :: is_PeriodicX
    logical(c_bool), value :: is_PeriodicY
    type(c_ptr), value :: refElem_ptr
    integer(c_int), value :: NLocalMeshPerPrc
    integer(c_int), value :: NprcX
    integer(c_int), value :: NprcY
    integer(c_int), value :: nproc
    integer(c_int), value :: myrank

    type(CMeshRectDom2D), pointer :: handle
    type(QuadrilateralElement), pointer :: refElem_fptr
    !------------------------------------

    call create_handle(handle, ptr)

    call c_f_pointer( refElem_ptr, refElem_fptr )
    call handle%obj%Init( NeGX, NeGY, dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
      logical_c2f(is_PeriodicX), logical_c2f(is_PeriodicY),                   &
      refElem_fptr, NLocalMeshPerPrc, NprcX, NprcY,                           &
      nproc, myrank )
    return
  end function CMeshRectDom2D_Init

  subroutine CMeshRectDom2D_Final( ptr ) bind(C, name="CMeshRectDom2D_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle(ptr)
    return
  end subroutine CMeshRectDom2D_Final

  function CMeshRectDom2D_get_lcmesh2D( ptr, lcmeshID ) result(lcmesh) bind(C, name="CMeshRectDom2D_get_lcmesh2D")
    use scale_localmesh2d_cbind, only: CLocalMesh2DPtr
    implicit none
    type(c_ptr), value :: ptr
    integer(c_int), value :: lcmeshID
    type(c_ptr) :: lcmesh

    type(CMeshRectDom2D), pointer :: handle
    type(CLocalMesh2DPtr), pointer :: lmesh_ptr
    !------------------------------------
    call c_f_pointer(ptr, handle)
    allocate(lmesh_ptr)
    lmesh_ptr%obj => handle%obj%lcmesh_list(lcmeshID)
    lcmesh = c_loc(lmesh_ptr)
    return
  end function CMeshRectDom2D_get_lcmesh2D

!******
DEF_C_BIND(CMeshRectDom2D,"CMeshRectDom2D")
DEF_C_BIND_SUB(CMeshRectDom2D,Generate,"CMeshRectDom2D_Generate")
!***** Getter
DEF_C_BIND_GETTER_I(CMeshRectDom2D,get_NeGX,"CMeshRectDom2D_get_NeGX",NeGX)
DEF_C_BIND_GETTER_I(CMeshRectDom2D,get_NeGY,"CMeshRectDom2D_get_NeGY",NeGY)
DEF_C_BIND_GETTER_I(CMeshRectDom2D,get_NprcX,"CMeshRectDom2D_get_NprcX",NprcX)
DEF_C_BIND_GETTER_I(CMeshRectDom2D,get_NprcY,"CMeshRectDom2D_get_NprcY",NprcY)
DEF_C_BIND_GETTER_D(CMeshRectDom2D,get_xmin_gl,"CMeshRectDom2D_get_xmin_gl",xmin_gl)
DEF_C_BIND_GETTER_D(CMeshRectDom2D,get_xmax_gl,"CMeshRectDom2D_get_xmax_gl",xmax_gl)
!******

end module scale_mesh_rectdom2d_cbind