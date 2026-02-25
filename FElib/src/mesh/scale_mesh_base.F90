!-------------------------------------------------------------------------------
!> module FElib / Mesh / Base
!!
!! @par Description
!!      Base module to manage meshes for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
#include "scaleFElib.h"
module scale_mesh_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision 
  use scale_io 
  
  use scale_element_base, only: ElementBase
  use scale_localmesh_base, only: LocalMeshBase

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !> Derived type to manage an information of a mesh dimension
  type, public :: MeshDimInfo
    character(len=H_SHORT) :: name !< Name of the dimension
    character(len=H_MID) :: desc   !< Description of the dimension
    character(len=H_SHORT) :: unit !< Unit of the dimension
    logical :: positive_down       !< Flag whether the positive direction is downward (e.g., for vertical dimension)
  end type MeshDimInfo

  !> Base type to manage a computational mesh
  type, abstract, public :: MeshBase
    integer :: LOCAL_MESH_NUM         !< Number of local meshes in each MPI process
    integer :: PRC_NUM                !< Number of MPI processes
    integer :: LOCAL_MESH_NUM_global  !< Total number of local meshes across all MPI processes
    
    integer, allocatable :: tileID_globalMap(:,:)
    integer, allocatable :: tileFaceID_globalMap(:,:)
    integer, allocatable :: tilePanelID_globalMap(:,:)
    integer, allocatable :: tileID_global2localMap(:)
    integer, allocatable :: PRCrank_globalMap(:)

    class(ElementBase), pointer :: refElem       !< Pointer to an object with a reference element
    type(MeshDimInfo), allocatable :: dimInfo(:) !< Array of information for each dimension

    real(RP) :: dom_vol        !< Total volume of the computational domain

    logical :: isGenerated     !< Flag whether the mesh is generated
  contains
    procedure(MeshBase_get_localmesh), deferred :: GetLocalMesh
    procedure :: SetDimInfo => MeshBase_SetDimInfo
  end type MeshBase

  interface 
    subroutine MeshBase_get_localmesh( this, id, ptr_lcmesh )
      import MeshBase
      import LocalMeshBase
      class(MeshBase), target, intent(in) :: this
      integer, intent(in) :: id
      class(LocalMeshBase), pointer, intent(out) :: ptr_lcmesh
    end subroutine MeshBase_get_localmesh
  end interface

  public :: MeshBase_Init
  public :: MeshBase_Final
  public :: MeshBase_setGeometricInfo

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
  !> Initialize an object to manage a computational mesh
!OCL SERIAL
  subroutine MeshBase_Init( this, &
     ndimtype, refElem, NLocalMeshPerPrc, NsideTile, &
     nprocs                                          )
    
    use scale_prc, only: &
      PRC_nprocs
    implicit none

    class(MeshBase), intent(inout) :: this
    integer, intent(in) :: ndimtype                      !< Number of DIMTYPE
    class(ElementBase), intent(in), target :: refElem    !< An object with a reference element
    integer, intent(in) :: NLocalMeshPerPrc              !< Number of local meshes in each MPI process
    integer, intent(in) :: NsideTile                     !< Number of side tiles
    integer, intent(in), optional :: nprocs              !< MPI processes (if not provided, it will be set to the value from PRC_nprocs)

    integer :: n
    !-----------------------------------------------------------------------------
    
    if ( present(nprocs) ) then
      this%PRC_NUM = nprocs
    else
      this%PRC_NUM = PRC_nprocs
    end if

    this%LOCAL_MESH_NUM        = NLocalMeshPerPrc
    this%LOCAL_MESH_NUM_global = this%PRC_NUM * this%LOCAL_MESH_NUM

    this%refElem => refElem
        
    allocate( this%tileID_globalMap(NsideTile, this%LOCAL_MESH_NUM_global) )
    allocate( this%tileFaceID_globalMap(NsideTile, this%LOCAL_MESH_NUM_global) )
    allocate( this%tilePanelID_globalMap(NsideTile, this%LOCAL_MESH_NUM_global) )
    allocate( this%tileID_global2localMap(this%LOCAL_MESH_NUM_global) )
    allocate( this%PRCRank_globalMap(this%LOCAL_MESH_NUM_global) )
    allocate( this%dimInfo(ndimtype) )
    !$acc enter data create(this%tileID_globalMap, this%tileFaceID_globalMap, this%tilePanelID_globalMap, &
    !$acc   this%tileID_global2localMap, this%PRCRank_globalMap, this%dimInfo )
    
    this%isGenerated = .false.

    return
  end subroutine MeshBase_Init

  !> Finalize an object to manage a computational mesh
!OCL SERIAL
  subroutine MeshBase_Final( this )
    implicit none
    class(MeshBase), intent(inout) :: this
    !-----------------------------------------------------------------------------
  
    if ( allocated(this%tileID_globalMap) ) then
      !$acc exit data delete(this%tileID_globalMap, this%tileFaceID_globalMap, this%tilePanelID_globalMap, &
      !$acc   this%tileID_global2localMap, this%PRCRank_globalMap, this%dimInfo)
      deallocate( this%tileID_globalMap )
      deallocate( this%tileFaceID_globalMap )
      deallocate( this%tilePanelID_globalMap )
      deallocate( this%tileID_global2localMap )
      deallocate( this%PRCRank_globalMap )
      deallocate( this%dimInfo )
    end if
    return
  end subroutine MeshBase_Final

!OCL SERIAL
  subroutine MeshBase_setGeometricInfo( mesh, ndim )
    use scale_element_base, only: ElementBase3D
    use scale_localmesh_3d, only: LocalMesh3D
    implicit none
    
    class(LocalMeshBase), intent(inout) :: mesh
    integer, intent(in) :: ndim

    class(ElementBase), pointer :: refElem
    !-----------------------------------------------------------------------------

    refElem => mesh%refElem

    allocate( mesh%pos_en(refElem%Np,mesh%Ne,ndim) )
    !allocate( mesh%fx(refElem%Nfaces*refElem%Nfp,mesh%Ne) )
    !allocate( mesh%fy(refElem%Nfaces*refElem%Nfp,mesh%Ne) )
    allocate( mesh%normal_fn(refElem%NfpTot,mesh%Ne,ndim) )
    allocate( mesh%sJ(refElem%NfpTot,mesh%Ne) )
    allocate( mesh%J(refElem%Np,mesh%Ne) )
    allocate( mesh%Fscale(refElem%NfpTot,mesh%Ne) )
    allocate( mesh%Escale(refElem%Np,mesh%Ne,ndim,ndim) )
    allocate( mesh%Gsqrt(refElem%Np,mesh%NeA) )
    !$acc enter data create(mesh%pos_en, mesh%normal_fn, mesh%sJ, mesh%J, &
    !$acc   mesh%Fscale, mesh%Escale, mesh%Gsqrt)

    return
  end subroutine MeshBase_setGeometricInfo

!OCL SERIAL
  subroutine MeshBase_SetDimInfo( this, &
      dimID, name, unit, desc, positive_down )
    implicit none
    class(MeshBase), intent(inout) :: this
    integer, intent(in) :: dimID
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    character(len=*), intent(in) :: desc
    logical, intent(in), optional :: positive_down

    !-----------------------------------------------

    this%dimInfo(dimID)%name = name
    this%dimInfo(dimID)%unit = unit
    this%dimInfo(dimID)%desc = desc
    if ( present(positive_down) ) then
      this%dimInfo(dimID)%positive_down = positive_down      
    else
      this%dimInfo(dimID)%positive_down = .false.      
    end if

    return
  end subroutine MeshBase_SetDimInfo

end module scale_mesh_base