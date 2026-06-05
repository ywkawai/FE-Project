!-------------------------------------------------------------------------------
!> module FElib / Mesh / Hierarchy base
!!
!! @par Description
!!      A base module to provide a mesh hierarchy
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_mesh_hierarchy_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !
  use scale_precision
  use scale_io

  use scale_localmesh_base, only: LocalMeshBase

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !> Base type to manage local mesh data for multigrid
  type, public :: MeshHierarchyLocalMGDataBase
    real(RP), allocatable :: Ic2f(:,:,:)    !< Matrix for correction operation in h-MG
    integer, allocatable :: Ic2f_emap(:)    !< Mapping of elements used for correction operation in h-MG
    integer, allocatable :: If2c_emap(:,:)  !< Mapping of elements used for restriction operation in h-MG

    integer, allocatable :: CoarseLocalMesh_tileIDlist(:)
    integer, allocatable :: CoarseLocalMesh_lcdomIDlist(:)
  end type MeshHierarchyLocalMGDataBase
  public :: MeshHierarchyLocalMGDataBase_Init
  public :: MeshHierarchyLocalMGDataBase_Final

  type, public :: MeshHierarchyLevelBase
    integer :: hierarchy_type
    integer :: level_id   = 0
  end type MeshHierarchyLevelBase

  !> Base type for mesh hierarchy
  type, public :: MeshHierarchyBase
    integer :: NUM_TOT_LEVEL   !< Total number of levels
    integer :: NUM_hMG_LEVEL   !< Number of h-MG levels
    integer :: NUM_pMG_LEVEL   !< Number of p-MG levels
  contains
  end type MeshHierarchyBase
  public :: MeshHierarchyBase_Init
  public :: MeshHierarchyBase_Final
  
  public :: MeshHierarchy_construct_pMG_mat1D
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: MESH_HIERARCHY_pMG_FINEST_LEVEL = 1  !< Finest level index in p-MG
  integer, public, parameter :: MESH_HIERARCHY_hMG_FINEST_LEVEL = 1  !< Finest level index in h-MG

  integer, public, parameter :: MESH_HIERARCHY_TYPE_pMG         = 1  !< Type ID of mesh hierarchy: p-MG
  integer, public, parameter :: MESH_HIERARCHY_TYPE_hMG         = 2  !< Type ID of mesh hierarchy: h-MG

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
contains
  !> Initialize a base object for mesh hierarchy
!OCL SERIAL
  subroutine MeshHierarchyBase_Init( this, &
    p_LEVEL_NUM, &
    h_LEVEL_NUM  )
    implicit none
    class(MeshHierarchyBase), intent(inout) :: this
    integer, intent(in) :: p_LEVEL_NUM   !< Number of p-mesh levels
    integer, intent(in) :: h_LEVEL_NUM   !< Number of h-mesh levels
    !-------------------------------------------------------------

    this%NUM_TOT_LEVEL = 0
    this%NUM_hMG_LEVEL = h_LEVEL_NUM
    this%NUM_pMG_LEVEL = p_LEVEL_NUM

    this%NUM_TOT_LEVEL = this%NUM_hMG_LEVEL + this%NUM_pMG_LEVEL
    return
  end subroutine MeshHierarchyBase_Init

  !> Finalize a base object for mesh hierarchy
!OCL SERIAL
  subroutine MeshHierarchyBase_Final(this)
    implicit none
    class(MeshHierarchyBase), intent(inout) :: this
    !-------------------------------------------------------------
    return
  end subroutine MeshHierarchyBase_Final

  !> Construct a 1D p-multigrid matrix used for transfer between different p-levels
!OCL SERIAL
  subroutine MeshHierarchy_construct_pMG_mat1D( mat1D, np_i, np_o )
    use scale_polynomial, only: &
      Polynomial_GenLagrangePoly, &
      Polynomial_GenGaussLegendrePt,        &
      Polynomial_GenGaussLegendrePtIntWeight
    use scale_element_line, only: LineElement
    implicit none
    integer, intent(in) :: np_i               !< Number of points with input DOF
    integer, intent(in) :: np_o               !< Number of points with output DOF
    real(RP), intent(out) :: mat1D(np_o,np_i) !< Transfer matrix with p-MG

    type(LineElement) :: elem1D_i
    type(LineElement) :: elem1D_o

    real(RP), allocatable :: lag_i(:,:)
    real(RP), allocatable :: lag_o(:,:)

    integer :: np_int
    real(RP), allocatable :: int_w(:)
    real(RP), allocatable :: int_pts(:)

    integer :: pi, po
    !---------------------------------------------------------------------

    call elem1D_i%Init( np_i-1, .false. )
    call elem1D_o%Init( np_o-1, .false. )

    !-
    np_int = max(np_i,np_o)
    
    allocate( int_pts(np_int) )
    int_pts(:) = Polynomial_GenGaussLegendrePt(np_int)

    allocate( int_w(np_int) )     
    int_w(:) = Polynomial_GenGaussLegendrePtIntWeight(np_int)

    !-
    allocate( lag_i(np_int,np_i) )
    lag_i(:,:) = Polynomial_GenLagrangePoly(elem1D_i%PolyOrder, elem1D_i%x1, int_pts)

    allocate( lag_o(np_int,np_o) )
    lag_o(:,:) = Polynomial_GenLagrangePoly(elem1D_o%PolyOrder, elem1D_o%x1, int_pts)

    do pi=1, np_i
    do po=1, np_o
      mat1D(po,pi) = sum( int_w(:) * lag_i(:,pi) * lag_o(:,po) )
    end do
    end do
    mat1D(:,:) = matmul( elem1D_o%invM, mat1D )

    !-
    call elem1D_i%Final()
    call elem1D_o%Final()

    return
  end subroutine MeshHierarchy_construct_pMG_mat1D

!-
  !> Initialize a base object to manage local multigrid data for mesh hierarchy
!OCL SERIAL
  subroutine MeshHierarchyLocalMGDataBase_Init( this, &
    lcmesh )
    implicit none
    class(MeshHierarchyLocalMGDataBase), intent(inout) :: this
    class(LocalMeshBase), intent(in) :: lcmesh
    !-------------------------------------------------------------

    allocate( this%Ic2f(lcmesh%refElem%Np,lcmesh%refElem%Nv,lcmesh%NeA) )
    allocate( this%Ic2f_emap(lcmesh%Ne) )

    allocate( this%CoarseLocalMesh_tileIDlist(1) )
    allocate( this%CoarseLocalMesh_lcdomIDlist(1) )
    return
  end subroutine MeshHierarchyLocalMGDataBase_Init

  !> Finalize a base object to manage local multigrid data for mesh hierarchy
!OCL SERIAL
  subroutine MeshHierarchyLocalMGDataBase_Final( this )
    implicit none
    class(MeshHierarchyLocalMGDataBase), intent(inout) :: this
    !------------------------------------------------------------- 
    deallocate( this%Ic2f )
    deallocate( this%Ic2f_emap )
    deallocate( this%CoarseLocalMesh_tileIDlist )
    deallocate( this%CoarseLocalMesh_lcdomIDlist )
    return
  end subroutine MeshHierarchyLocalMGDataBase_Final
end module scale_mesh_hierarchy_base
