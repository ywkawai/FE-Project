#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_mesh_cubedom3d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase1D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_mesh_cubedom3d, only: MeshCubeDom3D

  use scale_mesh_base3d_cbind, only: CMeshBase3DPtr
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CMeshCubeDom3D
    type(MeshCubeDom3D) :: obj
  end type CMeshCubeDom3D
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CMeshCubeDom3D_Init( &
    NeGX, NeGY, NeGZ,                                           & ! (in)
    dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, & ! (in)
    is_PeriodicX, is_PeriodicY, is_PeriodicZ,                   & ! (in)
    refElem_ptr, NLocalMeshPerPrc,                              & ! (in)
    NprcX, NprcY, nproc, myrank, FZ_ptr )                       & ! (in)
    result(ptr) bind(C, name="CMeshCubeDom3D_Init")
    implicit none
    type(c_ptr) :: ptr
    integer(c_int), value :: NeGX
    integer(c_int), value :: NeGY
    integer(c_int), value :: NeGZ
    real(c_double), value :: dom_xmin
    real(c_double), value :: dom_xmax   
    real(c_double), value :: dom_ymin
    real(c_double), value :: dom_ymax
    real(c_double), value :: dom_zmin
    real(c_double), value :: dom_zmax
    logical(c_bool), value :: is_PeriodicX
    logical(c_bool), value :: is_PeriodicY
    logical(c_bool), value :: is_PeriodicZ
    type(c_ptr), value :: refElem_ptr
    integer(c_int), value :: NLocalMeshPerPrc
    integer(c_int), value :: NprcX
    integer(c_int), value :: NprcY
    integer(c_int), value :: nproc
    integer(c_int), value :: myrank
    type(c_ptr), value :: FZ_ptr

    type(CMeshCubeDom3D), pointer :: handle
    type(HexahedralElement), pointer :: refElem_fptr

    real(RP), pointer :: FZ_fptr(:)
    !------------------------------------

    call create_handle(handle, ptr)

    call c_f_pointer( refElem_ptr, refElem_fptr )

    if ( c_associated(FZ_ptr) ) then
      call handle%obj%Init( NeGX, NeGY, NeGZ, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
        logical_c2f(is_PeriodicX), logical_c2f(is_PeriodicY), logical_c2f(is_PeriodicZ),                  &
        refElem_fptr, NLocalMeshPerPrc, NprcX, NprcY,                           &
        nproc, myrank )
    else
      call c_f_pointer( FZ_ptr, FZ_fptr, [NeGZ+1] )
      call handle%obj%Init( NeGX, NeGY, NeGZ, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
        logical_c2f(is_PeriodicX), logical_c2f(is_PeriodicY), logical_c2f(is_PeriodicZ),                  &
        refElem_fptr, NLocalMeshPerPrc, NprcX, NprcY,                           &
        nproc, myrank, FZ=FZ_fptr )
    end if
    return
  end function CMeshCubeDom3D_Init

  subroutine CMeshCubeDom3D_Final( ptr ) bind(C, name="CMeshCubeDom3D_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle(ptr)
    return
  end subroutine CMeshCubeDom3D_Final

  function CMeshCubeDom3D_get_lcmesh3D( ptr, lcmeshID ) result(lcmesh) bind(C, name="CMeshCubeDom3D_get_lcmesh3D")
    use scale_localmesh3d_cbind, only: CLocalMesh3DPtr
    implicit none
    type(c_ptr), value :: ptr
    integer(c_int), value :: lcmeshID
    type(c_ptr) :: lcmesh

    type(CMeshCubeDom3D), pointer :: handle
    type(CLocalMesh3DPtr), pointer :: lmesh_ptr
    !------------------------------------
    call c_f_pointer(ptr, handle)
    allocate(lmesh_ptr)
    lmesh_ptr%obj => handle%obj%lcmesh_list(lcmeshID)
    lcmesh = c_loc(lmesh_ptr)
    return
  end function CMeshCubeDom3D_get_lcmesh3D

!******
DEF_C_BIND(CMeshCubeDom3D,"CMeshCubeDom3D",)
DEF_C_BIND_SUB(CMeshCubeDom3D,Generate,"CMeshCubeDom3D_Generate")
!***** Getter
DEF_C_BIND_GETTER_I(CMeshCubeDom3D,get_NeGX,"CMeshCubeDom3D_get_NeGX",NeGX)
DEF_C_BIND_GETTER_I(CMeshCubeDom3D,get_NeGY,"CMeshCubeDom3D_get_NeGY",NeGY)
DEF_C_BIND_GETTER_I(CMeshCubeDom3D,get_NeGZ,"CMeshCubeDom3D_get_NeGZ",NeGZ)
DEF_C_BIND_GETTER_I(CMeshCubeDom3D,get_NprcX,"CMeshCubeDom3D_get_NprcX",NprcX)
DEF_C_BIND_GETTER_I(CMeshCubeDom3D,get_NprcY,"CMeshCubeDom3D_get_NprcY",NprcY)
DEF_C_BIND_GETTER_D(CMeshCubeDom3D,get_xmin_gl,"CMeshCubeDom3D_get_xmin_gl",xmin_gl)
DEF_C_BIND_GETTER_D(CMeshCubeDom3D,get_xmax_gl,"CMeshCubeDom3D_get_xmax_gl",xmax_gl)
DEF_C_BIND_GETTER_D(CMeshCubeDom3D,get_ymin_gl,"CMeshCubeDom3D_get_ymin_gl",ymin_gl)
DEF_C_BIND_GETTER_D(CMeshCubeDom3D,get_ymax_gl,"CMeshCubeDom3D_get_ymax_gl",ymax_gl)
DEF_C_BIND_GETTER_D(CMeshCubeDom3D,get_zmin_gl,"CMeshCubeDom3D_get_zmin_gl",zmin_gl)
DEF_C_BIND_GETTER_D(CMeshCubeDom3D,get_zmax_gl,"CMeshCubeDom3D_get_zmax_gl",zmax_gl)
!******
DEF_C_BIND_MAKE_FOBJ_PARENT_HANDLE(CMeshCubeDom3D,get_MeshBase3D,"CMeshCubeDom3D_get_MeshBase3D",CMeshBase3DPtr)
!******

end module scale_mesh_cubedom3d_cbind