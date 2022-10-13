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
  use scale_meshfield_base, only: MeshField3D

  use scale_sparsemat, only: sparsemat

  use scale_model_meshbase_manager, only: &
    ModelMeshBase, &
    ModelMeshBase1D, ModelMeshBase2D, ModelMeshBase3D

  use scale_model_var_manager, only: ModelVarManager
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  
  type, extends(ModelMeshBase1D), public :: ModelMesh1D
  contains
    procedure, public :: ModelMesh1D_Init
    procedure, public :: ModelMesh1D_Final
  end type ModelMesh1D

  type, extends(ModelMeshBase2D), public :: ModelMesh2D
  contains
    procedure, public :: ModelMesh2D_Init
    procedure, public :: ModelMesh2D_Final
  end type ModelMesh2D

  type, abstract, extends(ModelMeshBase3D), public :: ModelMesh3D
  contains
    procedure, public :: ModelMesh3D_Init
    procedure, public :: ModelMesh3D_Final
    procedure(ModelMesh3D_create_communicator), public, deferred :: Create_communicator
  end type ModelMesh3D

  interface
    subroutine ModelMesh3D_create_communicator( this, sfield_num, hvfield_num, var_manager, field_list, commid )
      import ModelMesh3D
      import MeshBase3D
      import ModelVarManager
      import MeshField3D
      class(ModelMesh3D), target, intent(inout) :: this
      integer, intent(in) :: sfield_num
      integer, intent(in) :: hvfield_num
      class(ModelVarManager), intent(inout) :: var_manager
      class(MeshField3D), intent(in) :: field_list(:)
      integer, intent(out) :: commid
    end subroutine ModelMesh3D_create_communicator
  end interface

  ! Cascade
  public :: ModelMeshBase

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

  !* 1D *************************************************************

  subroutine ModelMesh1D_Init( this, mesh )
    implicit none
    class(ModelMesh1D), target, intent(inout) :: this
    class(MeshBase1D), target, intent(in) :: mesh
    !-----------------------------------------------------

    call this%ModelMeshBase1D_Init(mesh)

    return
  end subroutine ModelMesh1D_Init

  subroutine ModelMesh1D_Final( this )
    implicit none
    class(ModelMesh1D), target, intent(inout) :: this

    !-----------------------------------------------------

    call this%ModelMeshBase1D_Final()

    return
  end subroutine ModelMesh1D_Final

  !* 2D *************************************************************

  subroutine ModelMesh2D_Init( this, mesh )
    implicit none
    class(ModelMesh2D), target, intent(inout) :: this
    class(MeshBase2D), target, intent(in) :: mesh
    !-----------------------------------------------------

    call this%ModelMeshBase2D_Init(mesh)

    return
  end subroutine ModelMesh2D_Init

  subroutine ModelMesh2D_Final( this )
    implicit none
    class(ModelMesh2D), target, intent(inout) :: this  
    !-----------------------------------------------------

    call this%ModelMeshBase2D_Final()

    return
  end subroutine ModelMesh2D_Final

  !* 3D *************************************************************

  subroutine ModelMesh3D_Init( this, mesh )
    implicit none
    class(ModelMesh3D), target, intent(inout) :: this
    class(MeshBase3D), target, intent(in) :: mesh
    !-----------------------------------------------------

    this%ptr_mesh => mesh
    call this%ModelMeshBase3D_Init(mesh)

    return
  end subroutine ModelMesh3D_Init  

  subroutine ModelMesh3D_Final( this )
    implicit none
    class(ModelMesh3D), target, intent(inout) :: this
    !-----------------------------------------------------

    nullify( this%ptr_mesh )
    call this%ModelMeshBase3D_Final()

    return
  end subroutine ModelMesh3D_Final

end module scale_model_mesh_manager