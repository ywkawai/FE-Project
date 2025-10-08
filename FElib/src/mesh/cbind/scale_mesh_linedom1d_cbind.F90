#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_mesh_linedom1d_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f
  use scale_element_base, only: ElementBase1D
  use scale_element_line, only: LineElement
  use scale_mesh_linedom1d, only: MeshLineDom1D

  use scale_mesh_base1d_cbind, only: CMeshBase1DPtr
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CMeshLineDom1D
    type(MeshLineDom1D) :: obj
  end type CMeshLineDom1D
  type, public, extends(CBindingBase) :: CMeshLineDom1DPtr
    type(MeshLineDom1D), pointer :: obj
  end type CMeshLineDom1DPtr

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CMeshLineDom1D_Init( &
    NeG,                                    & ! (in)
    dom_xmin, dom_xmax,                     & ! (in)
    refElem_ptr, NLocalMeshPerPrc,          & ! (in)
    nproc, myrank, FX_ptr )                 & ! (in)
    result(ptr) bind(C, name="CMeshLineDom1D_Init")
    implicit none
    type(c_ptr) :: ptr
    integer(c_int), value :: NeG
    real(c_double), value :: dom_xmin
    real(c_double), value :: dom_xmax   
    type(c_ptr), value :: refElem_ptr
    integer(c_int), value :: NLocalMeshPerPrc
    integer(c_int), value :: nproc
    integer(c_int), value :: myrank
    type(c_ptr), value :: FX_ptr

    type(CMeshLineDom1D), pointer :: handle

    type(LineElement), pointer :: refElem_fptr
    real(c_double), pointer :: FX_fptr(:)
    !------------------------------------

    call create_handle(handle, ptr)

    call c_f_pointer( refElem_ptr, refElem_fptr )
    if ( c_associated(FX_ptr) ) then
      call c_f_pointer( FX_ptr, FX_fptr, [NeG+1] )
      call handle%obj%Init( NeG, dom_xmin, dom_xmax, refElem_fptr, NLocalMeshPerPrc, nproc, myrank, FX_fptr )
    else
      call handle%obj%Init( NeG, dom_xmin, dom_xmax, refElem_fptr, NLocalMeshPerPrc, nproc, myrank )
    end if
    return
  end function CMeshLineDom1D_Init

  subroutine CMeshLineDom1D_Final( ptr ) bind(C, name="CMeshLineDom1D_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle(ptr)
    return
  end subroutine CMeshLineDom1D_Final

  function CMeshLineDom1D_get_lcmesh1D( ptr, lcmeshID ) result(lcmesh) bind(C, name="CMeshLineDom1D_get_lcmesh1D")
    use scale_localmesh1d_cbind, only: CLocalMesh1DPtr
    implicit none
    type(c_ptr), value :: ptr
    integer(c_int), value :: lcmeshID
    type(c_ptr) :: lcmesh

    type(CMeshLineDom1D), pointer :: handle
    type(CLocalMesh1DPtr), pointer :: lmesh_ptr
    !------------------------------------
    call c_f_pointer(ptr, handle)
    allocate(lmesh_ptr)
    lmesh_ptr%obj => handle%obj%lcmesh_list(lcmeshID)
    lcmesh = c_loc(lmesh_ptr)
    return
  end function CMeshLineDom1D_get_lcmesh1D

!******
DEF_C_BIND(CMeshLineDom1D,"CMeshLineDom1D",)
DEF_C_BIND_SUB(CMeshLineDom1D,Generate,"CMeshLineDom1D_Generate")
!***** Getter
DEF_C_BIND_GETTER_I(CMeshLineDom1D,get_NeG,"CMeshLineDom1D_get_NeG",NeG)
DEF_C_BIND_GETTER_I(CMeshLineDom1D,get_Nprc,"CMeshLineDom1D_get_Nprc",Nprc)
DEF_C_BIND_GETTER_D(CMeshLineDom1D,get_xmin_gl,"CMeshLineDom1D_get_xmin_gl",xmin_gl)
DEF_C_BIND_GETTER_D(CMeshLineDom1D,get_xmax_gl,"CMeshLineDom1D_get_xmax_gl",xmax_gl)
!******
DEF_C_BIND_MAKE_FOBJ_PARENT_HANDLE(CMeshLineDom1D,get_MeshBase1D,"CMeshLineDom1D_get_MeshBase1D",CMeshBase1DPtr)
!******

end module scale_mesh_linedom1d_cbind