#include "scalelib.h"
module scale_localmesh_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_element_base, only: elementbase

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  type, public :: LocalMeshBase
    integer :: Ne
    integer :: NeS
    integer :: NeE
    integer :: NeA
    integer :: Nv

    class(ElementBase), pointer :: refElem
    
    real(RP), allocatable :: pos_ev(:,:)
    real(RP), allocatable :: pos_en(:,:,:)
    real(RP), allocatable :: normal_fn(:,:,:)

    real(RP), allocatable :: sJ(:,:)
    real(RP), allocatable :: J(:,:)
    
    real(RP), allocatable :: Escale(:,:,:,:)
    real(RP), allocatable :: Fscale(:,:)

    integer, allocatable :: EToV(:,:)
    integer, allocatable :: EToE(:,:)
    integer, allocatable :: EToF(:,:)
    integer, allocatable :: VMapM(:,:)
    integer, allocatable :: VMapP(:,:)
    integer, allocatable :: MapM(:,:)
    integer, allocatable :: MapP(:,:)

    integer, allocatable :: BCType(:,:)
    integer, allocatable :: MapB(:)
    integer, allocatable :: VMapB(:)
    integer, allocatable :: MapD(:)
    integer, allocatable :: VMapD(:)
    integer, allocatable :: MapN(:)
    integer, allocatable :: VMapN(:)

    integer :: tileID
    integer :: panelID
    integer :: PRC_myrank
  
    real(DP), allocatable :: G_ij(:,:,:,:)
    real(DP), allocatable :: GIJ(:,:,:,:)
    real(DP), allocatable :: Gsqrt(:,:)    
  end type LocalMeshBase

  public :: LocalMeshBase_Init
  public :: LocalMeshBase_Final

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: BCTYPE_INTERIOR  = 0
  integer, public, parameter :: BCTYPE_UNDEFBC   = 1
  integer, public, parameter :: BCTYPE_PERIODIC  = 2
  integer, public, parameter :: BCTYPE_DIRCHLET  = 3
  integer, public, parameter :: BCTYPE_NEUMAN    = 4
  integer, public, parameter :: BCTYPE_SHORELINE = 5

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
  subroutine LocalMeshBase_Init( this, refElem, dims )

    class(LocalMeshBase), intent(inout) :: this
    class(ElementBase), intent(in), target :: refElem
    integer, intent(in) :: dims

    this%refElem => refElem
    
  end subroutine LocalMeshBase_Init

  subroutine LocalMeshBase_Final( this )

    class(LocalMeshBase), intent(inout) :: this

  end subroutine LocalMeshBase_Final

end module scale_localmesh_base