#include "scalelib.h"
module scale_mesh_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision  
  use scale_element_base, only: ElementBase
  use scale_localmesh_base, only: LocalMeshBase

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, abstract, public :: Meshbase
    integer :: LOCAL_MESH_NUM
    integer :: PRC_NUM
    integer :: LOCAL_MESH_NUM_global
    
    integer, allocatable :: tileID_globalMap(:,:)
    integer, allocatable :: tileFaceID_globalMap(:,:)
    integer, allocatable :: tilePanelID_globalMap(:,:)
    integer, allocatable :: tileID_global2localMap(:)
    integer, allocatable :: PRCrank_globalMap(:)

    class(ElementBase), pointer :: refElem
    logical :: isGenerated
  end type MeshBase

  public :: MeshBase_Init, MeshBase_Final
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
  subroutine MeshBase_Init(this, &
    & refElem, NLocalMeshPerPrc, NsideTile )
    
    use scale_prc, only: PRC_nprocs, PRC_myrank
    implicit none

    class(MeshBase), intent(inout) :: this
    class(ElementBase), intent(in), target :: refElem
    integer, intent(in) :: NLocalMeshPerPrc
    integer, intent(in) :: NsideTile

    integer :: n

    !-----------------------------------------------------------------------------
    
    this%PRC_NUM               = PRC_nprocs
    this%LOCAL_MESH_NUM        = NLocalMeshPerPrc
    this%LOCAL_MESH_NUM_global = PRC_nprocs * this%LOCAL_MESH_NUM
        
    allocate( this%tileID_globalMap(NsideTile, this%LOCAL_MESH_NUM_global) )
    allocate( this%tileFaceID_globalMap(NsideTile, this%LOCAL_MESH_NUM_global) )
    allocate( this%tilePanelID_globalMap(NsideTile, this%LOCAL_MESH_NUM_global) )
    allocate( this%tileID_global2localMap(this%LOCAL_MESH_NUM_global) )
    allocate( this%PRCRank_globalMap(this%LOCAL_MESH_NUM_global) )
    
    this%isGenerated = .false.
  end subroutine MeshBase_Init

  subroutine MeshBase_Final( this )
    
    class(MeshBase), intent(inout) :: this

    integer :: n

    !-----------------------------------------------------------------------------
  
    deallocate( this%tileID_globalMap )
    deallocate( this%tileFaceID_globalMap )
    deallocate( this%tilePanelID_globalMap )
    deallocate( this%tileID_global2localMap )
    deallocate( this%PRCRank_globalMap )

  end subroutine MeshBase_Final

  subroutine MeshBase_setGeometricInfo( mesh, ndim )

    implicit none
    
    class(LocalMeshBase), intent(inout) :: mesh
    integer, intent(in) :: ndim

    class(elementbase), pointer :: refElem

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
    allocate( mesh%Gsqrt(refElem%Np,mesh%Ne) )
    allocate( mesh%G_ij(refElem%Np,mesh%Ne, ndim,ndim) )
    allocate( mesh%GIJ (refElem%Np,mesh%Ne, ndim,ndim) )

  end subroutine MeshBase_setGeometricInfo

end module scale_mesh_base