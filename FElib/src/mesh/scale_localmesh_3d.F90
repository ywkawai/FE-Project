#include "scaleFElib.h"
module scale_localmesh_3d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_localmesh_base, only: &
    LocalMeshBase, LocalMeshBase_Init, LocalMeshBase_Final
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_element_base, only: elementbase, elementbase3D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(LocalMeshBase), public :: LocalMesh3D

    type(ElementBase3D), pointer :: refElem3D
    real(DP) :: xmin, xmax
    real(DP) :: ymin, ymax
    real(DP) :: zmin, zmax
    integer :: NeX, NeY, Ne2D, NeZ

    real(DP), allocatable :: Sz(:,:)
    real(DP), allocatable :: zS(:,:)

    class(LocalMesh2D), pointer :: lcmesh2D
    real(DP), allocatable :: lon2D(:,:)     
    real(DP), allocatable :: lat2D(:,:)

    integer, allocatable :: EMap3Dto2D(:)  
  contains
    procedure :: SetLocalMesh2D => LocalMesh3D_setLocalMesh2D  
  end type LocalMesh3D

  public :: LocalMesh3D_Init, LocalMesh3D_Final

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
  subroutine LocalMesh3D_Init( this, &
    refElem, PRC_myrank )
    
    implicit none

    class(LocalMesh3D), intent(inout) :: this
    class(ElementBase3D), intent(in), target :: refElem
    integer, intent(in) :: PRC_myrank
  
    this%PRC_myrank   = PRC_myrank
    this%refElem3D    => refElem
    nullify( this%lcmesh2D )

    call LocalMeshBase_Init(this, refElem, 3)

  end subroutine LocalMesh3D_Init

  subroutine LocalMesh3D_Final( this )
    type(LocalMesh3D), intent(inout) :: this

    call LocalMeshBase_Final(this)

  end subroutine LocalMesh3D_Final
  
  subroutine LocalMesh3D_setLocalMesh2D( this, lcmesh2D )
    implicit none
    class(LocalMesh3D), intent(inout) :: this
    class(LocalMesh2D), intent(in), target :: lcmesh2D
    !----------------

    this%lcmesh2D => lcmesh2D

    return
  end subroutine LocalMesh3D_setLocalMesh2D

end module scale_localmesh_3d