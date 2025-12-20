!-------------------------------------------------------------------------------
!> module FElib / Mesh / 2D domain
!!
!! @par Description
!!      Manage mesh hierarchy of 2D domain for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_multigrid_smoother_base
   !-----------------------------------------------------------------------------
   !
   !++ used modules
   !
   !
   use scale_precision
   use scale_io
   use scale_prc, only: PRC_abort

   use scale_sparsemat, only: SparseMat

   use scale_element_base, only: ElementBase2D
   use scale_element_quadrilateral, only: QuadrilateralElement

   use scale_localmesh_2d, only: LocalMesh2D
   use scale_mesh_base2d, only: MeshBase2D
   use scale_mesh_rectdom2d, only: MeshRectDom2D

   use scale_meshfield_base, only: MeshField2D
   use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

   !-----------------------------------------------------------------------------
   implicit none
   private

   !-----------------------------------------------------------------------------
   !
   !++ Public type & procedure
   ! 
   type, abstract, public :: MGSmoother2DBase
   contains
     procedure(MGSmoother2DBase_Advance_itr_1step), deferred :: Advance_itr_1step
   end type MGSmoother2DBase
   public :: MGSmoother2DBase_Init
   public :: MGSmoother2DBase_Final

  abstract interface
    subroutine MGSmoother2DBase_Advance_itr_1step( this, q, res, &
      f, aux_var, aux_hvec_comp, itr, var_comm, aux_comm, Dx, Dy, mesh2D, &
      cal_res_flag, zero_initial_guess, mg_p_level, mg_h_level            )
      import MeshBase2D
      import :: MeshField2D
      import MeshFieldCommRectDom2D
      import SparseMat
      import :: MGSmoother2DBase
      class(MGSmoother2DBase), intent(inout) :: this
      type(MeshField2D), intent(inout), target :: q
      type(MeshField2D), intent(inout) :: res
      type(MeshField2D), intent(in) :: f
      type(MeshField2D), intent(inout), target :: aux_var(:)
      type(MeshField2D), intent(inout), target :: aux_hvec_comp(:)
      integer, intent(in) :: itr
      class(MeshFieldCommRectDom2D), intent(inout) :: var_comm
      class(MeshFieldCommRectDom2D), intent(inout) :: aux_comm
      type(SparseMat), intent(in) :: Dx
      type(SparseMat), intent(in) :: Dy
      class(MeshBase2D), intent(in), target :: mesh2D
      logical, intent(in) :: cal_res_flag
      logical, intent(in) :: zero_initial_guess
      integer, intent(in) :: mg_p_level
      integer, intent(in) :: mg_h_level
    end subroutine MGSmoother2DBase_Advance_itr_1step
  end interface  

contains
!OCL SERIAL
  subroutine MGSmoother2DBase_Init(this)
    implicit none
    class(MGSmoother2DBase), intent(inout) :: this
    !---------------------------------
  end subroutine MGSmoother2DBase_Init

!OCL SERIAL
  subroutine MGSmoother2DBase_Final(this)
    implicit none
    class(MGSmoother2DBase), intent(inout) :: this
    !---------------------------------
  end subroutine MGSmoother2DBase_Final
end module scale_multigrid_smoother_base
