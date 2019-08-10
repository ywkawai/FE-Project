#include "scalelib.h"
module scale_localmesh_1d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_localmesh_base, only: &
    LocalMeshBase, LocalMeshBase_Init, LocalMeshBase_Final
  use scale_element_base, only: elementbase, elementbase1D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(LocalMeshBase), public :: LocalMesh1D

    type(elementbase1D), pointer :: refElem1D
    real(RP) :: xmin, xmax 
  end type LocalMesh1D

  public :: LocalMesh1D_Init, LocalMesh1D_Final

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
  subroutine LocalMesh1D_Init( this, &
    refElem, PRC_myrank )
    
    implicit none

    type(LocalMesh1D), intent(inout) :: this
    class(elementbase1D), intent(in), target :: refElem
    integer, intent(in) :: PRC_myrank
  
    this%PRC_myrank   = PRC_myrank
    this%refElem1D    => refElem
    
    call LocalMeshBase_Init(this, refElem, 1)

  end subroutine LocalMesh1D_Init

  subroutine LocalMesh1D_Final( this )
    type(LocalMesh1D), intent(inout) :: this

    call LocalMeshBase_Final(this)

  end subroutine LocalMesh1D_Final
  
end module scale_localmesh_1d