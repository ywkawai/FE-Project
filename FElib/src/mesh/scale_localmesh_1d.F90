!-------------------------------------------------------------------------------
!> module FElib / Mesh / Local 1D
!!
!! @par Description
!!      Module to manage 1D local mesh for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
#include "scaleFElib.h"
module scale_localmesh_1d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_localmesh_base, only: &
    LocalMeshBase, LocalMeshBase_Init, LocalMeshBase_Final
  use scale_element_base, only: ElementBase, ElementBase1D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(LocalMeshBase), public :: LocalMesh1D

    type(ElementBase1D), pointer :: refElem1D
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
!OCL SERIAL
  subroutine LocalMesh1D_Init( this, &
    lcdomID, refElem, myrank )
    
    implicit none

    type(LocalMesh1D), intent(inout) :: this
    integer, intent(in) :: lcdomID
    class(ElementBase1D), intent(in), target :: refElem
    integer, intent(in), optional :: myrank
    !-------------------------------------------------

    this%refElem1D => refElem
    call LocalMeshBase_Init(this, lcdomID, refElem, 1, myrank)

    return
  end subroutine LocalMesh1D_Init

!OCL SERIAL
  subroutine LocalMesh1D_Final( this, is_generated )
    implicit none

    type(LocalMesh1D), intent(inout) :: this
    logical, intent(in) :: is_generated
    !-------------------------------------------------

    call LocalMeshBase_Final( this, is_generated )

    return
  end subroutine LocalMesh1D_Final
  
end module scale_localmesh_1d