#include "scaleFElib.h"
module scale_model_mesh_manager
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: &
    PRC_abort

  use scale_mesh_base, only: MeshBase
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_sparsemat, only: sparsemat

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  type, abstract, public :: ModelMeshBase
    type(SparseMat), allocatable :: DOptrMat(:)
    type(SparseMat), allocatable :: SOptrMat(:)
    type(SparseMat) :: LiftOptrMat

    integer :: communicator_num
  contains
    procedure :: ModelMeshBase_Init
    procedure :: ModelMeshBase_Final
    procedure :: Get_communicatorID => ModelMeshBase_get_communicatorID
    procedure(ModelMeshBase_get_modelmesh), public, deferred :: GetModelMesh
  end type ModelMeshBase
  
  interface
    subroutine ModelMeshBase_get_modelmesh( this, ptr_mesh )
      import ModelMeshBase
      import MeshBase
      class(ModelMeshBase), target, intent(in) :: this
      class(MeshBase), pointer, intent(out) :: ptr_mesh
    end subroutine ModelMeshBase_get_modelmesh
  end interface

  type, extends(ModelMeshBase), public :: ModelMesh1D
    class(MeshBase1D), pointer :: ptr_mesh
  contains
    procedure, public :: ModelMesh1D_Init
    procedure, public :: ModelMesh1D_Final
    procedure, public :: GetModelMesh => ModelMesh1D_get_modelmesh
  end type ModelMesh1D

  type, extends(ModelMeshBase), public :: ModelMesh2D
    class(MeshBase2D), pointer :: ptr_mesh
  contains
    procedure, public :: ModelMesh2D_Init
    procedure, public :: ModelMesh2D_Final
    procedure, public :: GetModelMesh => ModelMesh2D_get_modelmesh
  end type ModelMesh2D

  type, extends(ModelMeshBase), public :: ModelMesh3D
    class(MeshBase3D), pointer :: ptr_mesh
  contains
    procedure, public :: ModelMesh3D_Init
    procedure, public :: ModelMesh3D_Final
    procedure, public :: GetModelMesh => ModelMesh3D_get_modelmesh
  end type ModelMesh3D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !------------------

contains
  subroutine ModelMeshBase_Init( this, nDim )
    implicit none
    class(ModelMeshBase), intent(inout) :: this
    integer, intent(in) :: nDim

    integer :: d
    !--------------------------------------------

    this%communicator_num = 0    
    allocate( this%SOptrMat(nDim), this%DOptrMat(nDim) )

    return
  end subroutine ModelMeshBase_Init

  function ModelMeshBase_get_communicatorID( this, max_communicator_num ) result(commid)
    implicit none
    class(ModelMeshBase), intent(inout) :: this
    integer, intent(in) :: max_communicator_num
    integer :: commid
    !--------------------------------------------

    this%communicator_num = this%communicator_num + 1
    commid = this%communicator_num

    if ( commid > max_communicator_num ) then
      LOG_ERROR('ModelMeshBase_get_communicatorID',*) 'The number of communicator exceeds expectation. Check!' 
      call PRC_abort
    end if

    return
  end function ModelMeshBase_get_communicatorID

  subroutine ModelMeshBase_Final( this )
    implicit none
    class(ModelMeshBase), intent(inout) :: this

    integer :: d
    !--------------------------------------------

    do d = 1, size(this%DOptrMat)
      call this%DOptrMat(d)%Final()
      call this%SOptrMat(d)%Final()
    end do
    deallocate( this%SOptrMat, this%DOptrMat )

    call this%LiftOptrMat%Final()

    return
  end subroutine ModelMeshBase_Final

  !* 1D *************************************************************

  subroutine ModelMesh1D_Init( this, mesh )
    implicit none
    class(ModelMesh1D), target, intent(inout) :: this
    class(MeshBase1D), target, intent(in) :: mesh
    !-----------------------------------------------------

    this%ptr_mesh => mesh
    call this%ModelMeshBase_Init(1)

    return
  end subroutine ModelMesh1D_Init

  subroutine ModelMesh1D_Final( this )
    implicit none
    class(ModelMesh1D), target, intent(inout) :: this

    integer :: d
    !-----------------------------------------------------

    nullify( this%ptr_mesh )
    call this%ModelMeshBase_Final()

    return
  end subroutine ModelMesh1D_Final

  subroutine ModelMesh1D_get_modelmesh( this, ptr_mesh )
    implicit none
    class(ModelMesh1D), target, intent(in) :: this
    class(MeshBase), pointer, intent(out) :: ptr_mesh
    !-----------------------------------------------------

    ptr_mesh => this%ptr_mesh

    return
  end subroutine ModelMesh1D_get_modelmesh

  !* 2D *************************************************************

  subroutine ModelMesh2D_Init( this, mesh )
    implicit none
    class(ModelMesh2D), target, intent(inout) :: this
    class(MeshBase2D), target, intent(in) :: mesh
    !-----------------------------------------------------

    this%ptr_mesh => mesh
    call this%ModelMeshBase_Init(2)

    return
  end subroutine ModelMesh2D_Init

  subroutine ModelMesh2D_Final( this )
    implicit none
    class(ModelMesh2D), target, intent(inout) :: this  
    !-----------------------------------------------------

    nullify( this%ptr_mesh )
    call this%ModelMeshBase_Final()

    return
  end subroutine ModelMesh2D_Final

  subroutine ModelMesh2D_get_modelmesh( this, ptr_mesh )
    implicit none
    class(ModelMesh2D), target, intent(in) :: this
    class(MeshBase), pointer, intent(out) :: ptr_mesh
    !-----------------------------------------------------

    ptr_mesh => this%ptr_mesh
    return
  end subroutine ModelMesh2D_get_modelmesh

  !* 3D *************************************************************

  subroutine ModelMesh3D_Init( this, mesh )
    implicit none
    class(ModelMesh3D), target, intent(inout) :: this
    class(MeshBase3D), target, intent(in) :: mesh
    !-----------------------------------------------------

    this%ptr_mesh => mesh
    call this%ModelMeshBase_Init(3)

    return
  end subroutine ModelMesh3D_Init  

  subroutine ModelMesh3D_Final( this )
    implicit none
    class(ModelMesh3D), target, intent(inout) :: this
    !-----------------------------------------------------

    nullify( this%ptr_mesh )
    call this%ModelMeshBase_Final()

    return
  end subroutine ModelMesh3D_Final

  subroutine ModelMesh3D_get_modelmesh( this, ptr_mesh )
    implicit none
    class(ModelMesh3D), target, intent(in) :: this
    class(MeshBase), pointer, intent(out) :: ptr_mesh
    !-----------------------------------------------------

    ptr_mesh => this%ptr_mesh
    return
  end subroutine ModelMesh3D_get_modelmesh

end module scale_model_mesh_manager