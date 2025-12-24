!-------------------------------------------------------------------------------
!> module FElib / Multigrid / Smoother base
!!
!! @par Description
!!      Manage multigrid smoother
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

   use scale_element_quadrilateral, only: QuadrilateralElement
   use scale_element_hexahedral, only: HexahedralElement

   use scale_mesh_base2d, only: MeshBase2D
   use scale_mesh_base3d, only: MeshBase3D
   use scale_mesh_rectdom2d, only: MeshRectDom2D
   use scale_mesh_cubedom3d, only: MeshCubeDom3D

   use scale_meshfield_base, only: &
     MeshField2D, MeshField3D
   use scale_meshfieldcomm_base, only: MeshFieldCommBase
   use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D
   use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
   !-----------------------------------------------------------------------------
   implicit none
   private

   !-----------------------------------------------------------------------------
   !
   !++ Public type & procedure
   ! 
   integer, public, parameter :: MGSmoother_PRE_ID  = 1
   integer, public,parameter :: MGSmoother_POST_ID = 2

   type, abstract, public :: MGSmoother2DBase
   contains
     procedure(MGSmoother2DBase_Advance_itr_1step), deferred :: Advance_itr_1step
   end type MGSmoother2DBase
   public :: MGSmoother2DBase_Init
   public :: MGSmoother2DBase_Final

   abstract interface
    subroutine MGSmoother2DBase_Advance_itr_1step( this, q, res, &
      f, aux_var, itr, var_comm, aux_comm, Dx, Dy, mesh2D,       &
      cal_res_flag, zero_initial_guess, mg_p_level, mg_h_level,  &
      pre_or_post_smooth )
      import MeshBase2D
      import :: MeshField2D
      import MeshFieldCommBase
      import SparseMat
      import :: MGSmoother2DBase
      class(MGSmoother2DBase), intent(inout) :: this
      type(MeshField2D), intent(inout), target :: q
      type(MeshField2D), intent(inout) :: res
      type(MeshField2D), intent(in) :: f
      type(MeshField2D), intent(inout), target :: aux_var(:)
      integer, intent(in) :: itr
      class(MeshFieldCommBase), intent(inout) :: var_comm
      class(MeshFieldCommBase), intent(inout) :: aux_comm
      type(SparseMat), intent(in) :: Dx
      type(SparseMat), intent(in) :: Dy
      class(MeshBase2D), intent(in), target :: mesh2D
      logical, intent(in) :: cal_res_flag
      logical, intent(in) :: zero_initial_guess
      integer, intent(in) :: mg_p_level
      integer, intent(in) :: mg_h_level
      integer, intent(in) :: pre_or_post_smooth
    end subroutine MGSmoother2DBase_Advance_itr_1step
   end interface  

   type, abstract, public :: MGSmoother3DBase
   contains
     procedure(MGSmoother3DBase_Advance_itr_1step), deferred :: Advance_itr_1step
   end type MGSmoother3DBase
   public :: MGSmoother3DBase_Init
   public :: MGSmoother3DBase_Final

  abstract interface
    subroutine MGSmoother3DBase_Advance_itr_1step( this, q, res, &
      f, aux_var, itr, var_comm, aux_comm, Dx, Dy, Dz, mesh3D,   &
      cal_res_flag, zero_initial_guess, mg_p_level, mg_h_level,  &
      pre_or_post_smooth )
      import MeshBase3D
      import :: MeshField3D
      import MeshFieldCommBase
      import SparseMat
      import :: MGSmoother3DBase
      class(MGSmoother3DBase), intent(inout) :: this
      type(MeshField3D), intent(inout), target :: q
      type(MeshField3D), intent(inout) :: res
      type(MeshField3D), intent(in) :: f
      type(MeshField3D), intent(inout), target :: aux_var(:)
      integer, intent(in) :: itr
      class(MeshFieldCommBase), intent(inout) :: var_comm
      class(MeshFieldCommBase), intent(inout) :: aux_comm
      type(SparseMat), intent(in) :: Dx
      type(SparseMat), intent(in) :: Dy
      type(SparseMat), intent(in) :: Dz
      class(MeshBase3D), intent(in), target :: mesh3D
      logical, intent(in) :: cal_res_flag
      logical, intent(in) :: zero_initial_guess
      integer, intent(in) :: mg_p_level
      integer, intent(in) :: mg_h_level
      integer, intent(in) :: pre_or_post_smooth
    end subroutine MGSmoother3DBase_Advance_itr_1step
  end interface  

contains
!-- 2D
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

!-- 3D
!OCL SERIAL
  subroutine MGSmoother3DBase_Init(this)
    implicit none
    class(MGSmoother3DBase), intent(inout) :: this
    !---------------------------------
  end subroutine MGSmoother3DBase_Init
!OCL SERIAL
  subroutine MGSmoother3DBase_Final(this)
    implicit none
    class(MGSmoother3DBase), intent(inout) :: this
    !---------------------------------
  end subroutine MGSmoother3DBase_Final
end module scale_multigrid_smoother_base
