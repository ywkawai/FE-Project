!-------------------------------------------------------------------------------
!> FElib / model framework / mesh manager (base)
!!
!! @par Description
!!          A module for managing mesh used in models
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_model_meshbase_manager
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: &
    PRC_abort

  use scale_element_base, only: ElementBase3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_element_operation_base, only: ElementOperationBase3D
  use scale_element_operation_general, only: ElementOperationGeneral
  use scale_element_operation_tensorprod3D, only: ElementOperationTensorProd3D

  use scale_sparsemat, only: sparsemat  

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  !> Derived type to manage model mesh and spatial operators (base type)
  type, abstract, public :: ModelMeshBase
    type(SparseMat), allocatable :: DOptrMat(:)
    type(SparseMat), allocatable :: SOptrMat(:)
    type(SparseMat) :: LiftOptrMat
    class(ElementOperationBase3D), pointer :: element3D_operation

    integer :: communicator_num
  contains
    procedure :: ModelMeshBase_Init
    procedure :: ModelMeshBase_Final
    procedure :: Get_communicatorID => ModelMeshBase_get_communicatorID
    procedure(ModelMeshBase_get_modelmesh), public, deferred :: GetModelMesh
  end type ModelMeshBase
  
  abstract interface
    subroutine ModelMeshBase_get_modelmesh( this, ptr_mesh )
      import ModelMeshBase
      import MeshBase
      class(ModelMeshBase), target, intent(in) :: this
      class(MeshBase), pointer, intent(out) :: ptr_mesh
    end subroutine ModelMeshBase_get_modelmesh
  end interface

  !> Derived type to manage 1D mesh and spatial operators
  type, extends(ModelMeshBase), public :: ModelMeshBase1D
    class(MeshBase1D), pointer :: ptr_mesh
  contains
    procedure, public :: ModelMeshBase1D_Init
    procedure, public :: ModelMeshBase1D_Final
    procedure, public :: GetModelMesh => ModelMeshBase1D_get_modelmesh
  end type ModelMeshBase1D

  !> Derived type to manage 2D mesh and spatial operators
  type, extends(ModelMeshBase), public :: ModelMeshBase2D
    class(MeshBase2D), pointer :: ptr_mesh
  contains
    procedure, public :: ModelMeshBase2D_Init
    procedure, public :: ModelMeshBase2D_Final
    procedure, public :: GetModelMesh => ModelMeshBase2D_get_modelmesh
  end type ModelMeshBase2D

  !> Derived type to manage 3D mesh and spatial operators
  type, extends(ModelMeshBase), abstract, public :: ModelMeshBase3D
    class(MeshBase3D), pointer :: ptr_mesh
    type(ElementOperationGeneral) :: element_operation_general
    class(ElementOperationTensorProd3D), allocatable :: element_operation_tensorprod
    logical :: initialized_element_operation
  contains
    procedure, public :: ModelMeshBase3D_Init
    procedure, public :: ModelMeshBase3D_Final
    procedure, public :: PrepairElementOperation => ModelMeshBase3D_prepair_ElementOperation
    procedure, public :: GetModelMesh => ModelMeshBase3D_get_modelmesh
  end type ModelMeshBase3D

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
!OCL SERIAL
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

!OCL SERIAL
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

!OCL SERIAL
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

  subroutine ModelMeshBase1D_Init( this, mesh )
    implicit none
    class(ModelMeshBase1D), target, intent(inout) :: this
    class(MeshBase1D), target, intent(in) :: mesh
    !-----------------------------------------------------

    this%ptr_mesh => mesh
    call this%ModelMeshBase_Init(1)

    return
  end subroutine ModelMeshBase1D_Init

!OCL SERIAL
  subroutine ModelMeshBase1D_Final( this )
    implicit none
    class(ModelMeshBase1D), target, intent(inout) :: this

    integer :: d
    !-----------------------------------------------------

    nullify( this%ptr_mesh )
    call this%ModelMeshBase_Final()

    return
  end subroutine ModelMeshBase1D_Final

!OCL SERIAL
  subroutine ModelMeshBase1D_get_modelmesh( this, ptr_mesh )
    implicit none
    class(ModelMeshBase1D), target, intent(in) :: this
    class(MeshBase), pointer, intent(out) :: ptr_mesh
    !-----------------------------------------------------

    ptr_mesh => this%ptr_mesh

    return
  end subroutine ModelMeshBase1D_get_modelmesh

  !* 2D *************************************************************

!OCL SERIAL
  subroutine ModelMeshBase2D_Init( this, mesh )
    implicit none
    class(ModelMeshBase2D), target, intent(inout) :: this
    class(MeshBase2D), target, intent(in) :: mesh
    !-----------------------------------------------------

    this%ptr_mesh => mesh
    call this%ModelMeshBase_Init(2)

    return
  end subroutine ModelMeshBase2D_Init

!OCL SERIAL
  subroutine ModelMeshBase2D_Final( this )
    implicit none
    class(ModelMeshBase2D), target, intent(inout) :: this  
    !-----------------------------------------------------

    nullify( this%ptr_mesh )
    call this%ModelMeshBase_Final()

    return
  end subroutine ModelMeshBase2D_Final

!OCL SERIAL
  subroutine ModelMeshBase2D_get_modelmesh( this, ptr_mesh )
    implicit none
    class(ModelMeshBase2D), target, intent(in) :: this
    class(MeshBase), pointer, intent(out) :: ptr_mesh
    !-----------------------------------------------------

    ptr_mesh => this%ptr_mesh
    return
  end subroutine ModelMeshBase2D_get_modelmesh

  !* 3D *************************************************************

!OCL SERIAL
  subroutine ModelMeshBase3D_Init( this, mesh )
    implicit none
    class(ModelMeshBase3D), target, intent(inout) :: this
    class(MeshBase3D), target, intent(in) :: mesh
    !-----------------------------------------------------

    this%ptr_mesh => mesh
    call this%ModelMeshBase_Init(3)

    this%initialized_element_operation = .false.

    return
  end subroutine ModelMeshBase3D_Init  

!OCL SERIAL
  subroutine ModelMeshBase3D_Final( this )
    implicit none
    class(ModelMeshBase3D), target, intent(inout) :: this
    !-----------------------------------------------------

    if ( this%initialized_element_operation ) then
      call this%element3D_operation%Final()
    end if

    nullify( this%ptr_mesh )
    call this%ModelMeshBase_Final()

    return
  end subroutine ModelMeshBase3D_Final

!OCL SERIAL
  subroutine ModelMeshBase3D_prepair_ElementOperation( this, element_operation_type, &
    SpMV_storage_format_ )
    use scale_prc, only: PRC_abort
    use scale_element_operation_tensorprod3D, only: ElementOperationTensorprod3D_create
    implicit none
    class(ModelMeshBase3D), intent(inout), target :: this
    character(len=*), intent(in), optional :: element_operation_type
    character(len=*), intent(in), optional :: SpMV_storage_format_

    character(len=H_SHORT) :: element_operation_type_
    character(len=H_SHORT) :: SpMV_storage_format

    class(ElementBase3D), pointer :: elem3D
    !-----------------------------------------------------

    if ( present(element_operation_type) ) then
      element_operation_type_ = element_operation_type
    else
      element_operation_type_ = "General"
    end if
    SpMV_storage_format = "ELL"
    elem3D => this%ptr_mesh%refElem3D
    call this%DOptrMat(1)%Init( elem3D%Dx1, storage_format=SpMV_storage_format )
    call this%DOptrMat(2)%Init( elem3D%Dx2, storage_format=SpMV_storage_format )
    call this%DOptrMat(3)%Init( elem3D%Dx3, storage_format=SpMV_storage_format )

    call this%SOptrMat(1)%Init( elem3D%Sx1, storage_format=SpMV_storage_format )
    call this%SOptrMat(2)%Init( elem3D%Sx2, storage_format=SpMV_storage_format )
    call this%SOptrMat(3)%Init( elem3D%Sx3, storage_format=SpMV_storage_format )  

    call this%LiftOptrMat%Init( elem3D%Lift, storage_format=SpMV_storage_format )

    select case(element_operation_type_)
    case ("General")
      call this%element_operation_general%Init( elem3D, this%DOptrMat(1), this%DOptrMat(2), this%DOptrMat(3), this%LiftOptrMat )
      this%element3D_operation => this%element_operation_general
    case ("TensorProd3D")
      call ElementOperationTensorprod3D_create( elem3D, &
        this%element_operation_tensorprod ) ! (out)
        this%element3D_operation => this%element_operation_tensorprod
    case default
      LOG_INFO("ModelMeshBase3D_prepair_ElementOperation",*) "The specified element_operation_type is not supported. Check!", trim(element_operation_type_)
      call PRC_abort
    end select

    this%initialized_element_operation = .true.

    return
  end subroutine ModelMeshBase3D_prepair_ElementOperation

!OCL SERIAL
  subroutine ModelMeshBase3D_get_modelmesh( this, ptr_mesh )
    implicit none
    class(ModelMeshBase3D), target, intent(in) :: this
    class(MeshBase), pointer, intent(out) :: ptr_mesh
    !-----------------------------------------------------

    ptr_mesh => this%ptr_mesh
    return
  end subroutine ModelMeshBase3D_get_modelmesh

end module scale_model_meshbase_manager