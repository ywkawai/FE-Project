!-------------------------------------------------------------------------------
!> module FElib / Multigrid / Solver base
!!
!! @par Description
!!      A base module to provide a multigrid solver
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

  !> Base type for multigrid solver
  type, public :: MultiGridSolverBase
    integer :: current_p_lev  !< Current p-MG level
    integer :: current_h_lev  !< Current h-MG level

    integer :: vcyc_num_max
    real(RP) :: history_residual_l2_initial  !< Initial residual L2 norm
    real(RP) :: history_residual_max_initial !< Initial residual max norm

    real(RP) :: threshold_ratio_residual_l2  !< Threshold of decreasing ratio of residual for smoothing (L2 norm)
    real(RP) :: threshold_residual_l2  !< Residual threshold for smoothing (L2 norm)
    real(RP) :: threshold_residual_max !< Residual threshold for smoothing (max norm)
  contains
    procedure :: Set_vcycle_parameter => MultiGridSolverBase_set_vcyc_param
    procedure :: Is_converged => MultiGridSolverBase_is_converged
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
  !> Initialize a base object for multigrid solver
!OCL SERIAL
  subroutine MultiGridSolverBase_Init( this, &
    mesh_hierarchy )
    implicit none
    class(MultiGridSolverBase), intent(inout) :: this
    class(MeshHierarchyBase), intent(in) :: mesh_hierarchy               !< Object to manage mesh hierarchy
    !-------------------------------------------------------------
    
    !-
    this%vcyc_num_max = 10
    return
  end subroutine MultiGridSolverBase_Init

  !> Finalize a base object for multigrid solver
!OCL SERIAL
  subroutine MultiGridSolverBase_Final(this)
    implicit none
    class(MultiGridSolverBase), intent(inout) :: this
    !-------------------------------------------------------------
    return
  end subroutine MultiGridSolverBase_Final

  !> Set parameters for V-cycle
!OCL SERIAL
  subroutine MultiGridSolverBase_set_vcyc_param( this, &
    vcyc_num_max, threshold_ratio_residual_l2, threshold_residual_l2, threshold_residual_max )
    implicit none
    class(MultiGridSolverBase), intent(inout) :: this
    integer, intent(in) :: vcyc_num_max
    real(RP), intent(in) :: threshold_ratio_residual_l2
    real(RP), intent(in) :: threshold_residual_l2
    real(RP), intent(in) :: threshold_residual_max
    !-------------------------------------------------------------
    this%vcyc_num_max = vcyc_num_max
    this%threshold_ratio_residual_l2 = threshold_ratio_residual_l2
    this%threshold_residual_l2 = threshold_residual_l2
    this%threshold_residual_max = threshold_residual_max
    return
  end subroutine MultiGridSolverBase_set_vcyc_param

  !> Check convergence of multigrid solver
!OCL SERIAL
  function MultiGridSolverBase_is_converged( this, smoother ) result( conv_flag )
    use scale_multigrid_smoother_base, only: MGSmootherBase
    implicit none
    class(MultiGridSolverBase), intent(in) :: this
    class(MGSmootherBase), intent(in) :: smoother
    logical :: conv_flag

    real(RP) :: res_l2, res_max
    real(RP) :: res_decrease_ratio_l2
    !---------------------------------------------------------------------

    call smoother%Get_current_residual_statistics( res_l2, res_max )
    
    res_decrease_ratio_l2 = res_l2 / this%history_residual_l2_initial
    if ( ( res_decrease_ratio_l2 < this%threshold_ratio_residual_l2 ) &
          .or. (      res_l2 < this%threshold_residual_l2             &
                .and. res_max < this%threshold_residual_max )         ) then
      conv_flag = .true.
    else
      conv_flag = .false.
    end if
    return
  end function MultiGridSolverBase_is_converged  

end module scale_multigrid_solver_base
