#include "scaleFElib.h"
module scale_localmeshfield_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_localmesh_1d, only: &
    LocalMesh1D, LocalMesh1D_Init, LocalMesh1D_Final

  use scale_localmesh_2d, only: &
    LocalMesh2D, LocalMesh2D_Init, LocalMesh2D_Final

  use scale_localmesh_3d, only: &
    LocalMesh3D, LocalMesh3D_Init, LocalMesh3D_Final

  use scale_element_base, only: elementbase2D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  type, public :: LocalMeshFieldBase
    real(RP), allocatable :: val(:,:)
  end type LocalMeshFieldBase

  type, extends(LocalMeshFieldBase), public :: LocalMeshField1D
    type(LocalMesh1D), pointer :: mesh => null()
  contains
    procedure :: Init => LocalMeshField1D_Init
    procedure :: Final => LocalMeshField1D_Final
  end type LocalMeshField1D

  type, extends(LocalMeshFieldBase), public :: LocalMeshField2D
    type(LocalMesh2D), pointer :: mesh => null()
  contains
    procedure :: Init => LocalMeshField2D_Init
    procedure :: Final => LocalMeshField2D_Final
  end type LocalMeshField2D

  type, extends(LocalMeshFieldBase), public :: LocalMeshField3D
    type(LocalMesh3D), pointer :: mesh => null()
  contains
    procedure :: Init => LocalMeshField3D_Init
    procedure :: Final => LocalMeshField3D_Final
  end type LocalMeshField3D

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
  subroutine LocalMeshField1D_Init(this, mesh)
    implicit none
    class(LocalMeshField1D), intent(inout) :: this
    class(LocalMesh1D), target, intent(in) :: mesh

    !-----------------------------------------------------------------------------

    this%mesh => mesh
    allocate( this%val(mesh%refElem%Np, mesh%NeA) )

    return
  end subroutine LocalMeshField1D_Init

  subroutine LocalMeshField1D_Final( this )
    implicit none
    class(LocalMeshField1D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    if (allocated(this%val)) deallocate( this%val )

    return
  end subroutine LocalMeshField1D_Final

  !**********
  
  subroutine LocalMeshField2D_Init(this, mesh)
    implicit none
    class(LocalMeshField2D), intent(inout) :: this
    class(LocalMesh2D), target, intent(in) :: mesh

    !-----------------------------------------------------------------------------

    this%mesh => mesh
    allocate( this%val(mesh%refElem%Np, mesh%NeA) )

    return
  end subroutine LocalMeshField2D_Init

  subroutine LocalMeshField2D_Final( this )
    implicit none
    class(LocalMeshField2D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    if (allocated(this%val)) deallocate( this%val )

    return
  end subroutine LocalMeshField2D_Final

  !**********
  
  subroutine LocalMeshField3D_Init(this, mesh)
    implicit none
    class(LocalMeshField3D), intent(inout) :: this
    class(LocalMesh3D), target, intent(in) :: mesh
    !-----------------------------------------------------------------------------

    this%mesh => mesh
    allocate( this%val(mesh%refElem%Np, mesh%NeA) )

    return
  end subroutine LocalMeshField3D_Init

  subroutine LocalMeshField3D_Final( this )
    implicit none
    class(LocalMeshField3D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    if (allocated(this%val)) deallocate( this%val )

    return
  end subroutine LocalMeshField3D_Final

end module scale_localmeshfield_base