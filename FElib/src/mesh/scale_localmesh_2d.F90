#include "scaleFElib.h"
module scale_localmesh_2d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_localmesh_base, only: &
    LocalMeshBase, LocalMeshBase_Init, LocalMeshBase_Final
  use scale_element_base, only: elementbase, elementbase2D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(LocalMeshBase), public :: LocalMesh2D

    type(ElementBase2D), pointer :: refElem2D
    real(DP) :: xmin, xmax
    real(DP) :: ymin, ymax
    integer :: NeX, NeY

    real(DP), allocatable :: lon(:,:)     
    real(DP), allocatable :: lat(:,:)     
  end type LocalMesh2D

  public :: LocalMesh2D_Init, LocalMesh2D_Final

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
  subroutine LocalMesh2D_Init( this, &
    refElem, myrank )
    
    implicit none

    class(LocalMesh2D), intent(inout) :: this
    class(ElementBase2D), intent(in), target :: refElem
    integer, intent(in), optional :: myrank
    !-------------------------------------------------

    this%refElem2D  => refElem
    call LocalMeshBase_Init(this, refElem, 2, myrank)

    return
  end subroutine LocalMesh2D_Init

  subroutine LocalMesh2D_Final( this, is_generated )
    implicit none
    type(LocalMesh2D), intent(inout) :: this
    logical, intent(in) :: is_generated
    !-------------------------------------------------

    if (is_generated) then
      deallocate( this%lon, this%lat )
    end if
    call LocalMeshBase_Final( this, is_generated )

    return
  end subroutine LocalMesh2D_Final
  
end module scale_localmesh_2d