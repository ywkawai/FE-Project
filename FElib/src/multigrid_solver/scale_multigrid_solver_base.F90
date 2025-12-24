!-------------------------------------------------------------------------------
!> module FElib / Multigrid / Solver base
!!
!! @par Description
!!      Manage a multigrid solver
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_multigrid_solver_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  use scale_mesh_hierarchy_base, only: &
    MeshHierarchyBase, &
    pMG_FINEST_LEVEL => MESH_HIERARCHY_pMG_FINEST_LEVEL, &
    hMG_FINEST_LEVEL => MESH_HIERARCHY_hMG_FINEST_LEVEL, &
    MESH_HIERARCHY_TYPE_pMG,                             &
    MESH_HIERARCHY_TYPE_hMG

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, public :: MultiGridSolverBase
    integer, allocatable :: p_itr_num_list(:)
    integer, allocatable :: h_itr_num_list(:)

    integer :: current_p_lev
    integer :: current_h_lev
  end type MultiGridSolverBase

  public :: MultiGridSolverBase_Init
  public :: MultiGridSolverBase_Final

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
contains
!OCL SERIAL
  subroutine MultiGridSolverBase_Init( this, &
    mesh_hierarchy,                 &
    p_itr_num_list, h_itr_num_list  )
    implicit none
    class(MultiGridSolverBase), intent(inout) :: this
    class(MeshHierarchyBase), intent(in) :: mesh_hierarchy
    integer, intent(in) :: p_itr_num_list(mesh_hierarchy%NUM_pMG_LEVEL)
    integer, intent(in) :: h_itr_num_list(mesh_hierarchy%NUM_hMG_LEVEL)
    !-------------------------------------------------------------
    
    allocate( this%p_itr_num_list( mesh_hierarchy%NUM_pMG_LEVEL ) )
    this%p_itr_num_list(:) = p_itr_num_list(:)

    allocate( this%h_itr_num_list( mesh_hierarchy%NUM_hMG_LEVEL ) )
    this%h_itr_num_list(:) = h_itr_num_list(:)
    return
  end subroutine MultiGridSolverBase_Init

!OCL SERIAL
  subroutine MultiGridSolverBase_Final(this)
    implicit none
    class(MultiGridSolverBase), intent(inout) :: this
    !-------------------------------------------------------------

    deallocate( this%p_itr_num_list )
    deallocate( this%h_itr_num_list )
    return
  end subroutine MultiGridSolverBase_Final
end module scale_multigrid_solver_base
