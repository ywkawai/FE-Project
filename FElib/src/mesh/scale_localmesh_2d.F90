!-------------------------------------------------------------------------------
!> module FElib / Mesh / Local 2D
!!
!! @par Description
!!      Module to manage 2D local mesh for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
#include "scaleFElib.h"
module scale_localmesh_2d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_localmesh_base, only: &
    LocalMeshBase, LocalMeshBase_Init, LocalMeshBase_Final
  use scale_element_base, only: ElementBase, ElementBase2D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  !> Derived type representing a local mesh for 2D domain
  type, extends(LocalMeshBase), public :: LocalMesh2D
    type(ElementBase2D), pointer :: refElem2D !< Pointer to the reference element for 2D
    
    real(RP) :: xmin !< Minimum x-coordinate of the local mesh
    real(RP) :: xmax !< Maximum x-coordinate of the local mesh
    real(RP) :: ymin !< Minimum y-coordinate of the local mesh
    real(RP) :: ymax !< Maximum y-coordinate of the local mesh

    integer :: NeX !< Number of elements in the x-direction
    integer :: NeY !< Number of elements in the y-direction

    real(RP), allocatable :: lon(:,:)  !< Array to save longitude coordinates of the local mesh nodes
    real(RP), allocatable :: lat(:,:)  !< Array to save latitude coordinates of the local mesh nodes
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

  !> Initialize an object to manage a local mesh for 2D domain
!OCL SERIAL
  subroutine LocalMesh2D_Init( this, &
    lcdomID, refElem, myrank )
    
    implicit none

    class(LocalMesh2D), intent(inout) :: this
    integer, intent(in) :: lcdomID
    class(ElementBase2D), intent(in), target :: refElem
    integer, intent(in), optional :: myrank
    !-------------------------------------------------

    this%refElem2D  => refElem
    call LocalMeshBase_Init(this, lcdomID, refElem, 2, myrank)
    !$acc enter data copyin(this)
    return
  end subroutine LocalMesh2D_Init

  !> Finalize an object to manage a local mesh for 2D domain
!OCL SERIAL
  subroutine LocalMesh2D_Final( this, is_generated )
    implicit none
    type(LocalMesh2D), intent(inout) :: this
    logical, intent(in) :: is_generated
    !-------------------------------------------------

    if (is_generated) then
      !$acc exit data delete(this%G_ij, this%GIJ)
      deallocate( this%G_ij, this%GIJ )
      !$acc exit data delete(this%lon, this%lat)
      deallocate( this%lon, this%lat )
    end if
    call LocalMeshBase_Final( this, is_generated )
    !$acc exit data delete(this)
    return
  end subroutine LocalMesh2D_Final
  
end module scale_localmesh_2d