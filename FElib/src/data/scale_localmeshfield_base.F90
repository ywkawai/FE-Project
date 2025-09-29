!> module FElib / Data / base
!!
!! @par Description
!!           A module for managing field data with subdomains included in each MPI process 
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_localmeshfield_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_localmesh_base, only: &
    LocalMeshBase
  
  use scale_localmesh_1d, only: &
    LocalMesh1D, LocalMesh1D_Init, LocalMesh1D_Final

  use scale_localmesh_2d, only: &
    LocalMesh2D, LocalMesh2D_Init, LocalMesh2D_Final

  use scale_localmesh_3d, only: &
    LocalMesh3D, LocalMesh3D_Init, LocalMesh3D_Final

  use scale_element_base, only: ElementBase2D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !> Derived type representing a field with local mesh (base type)
  type, public :: LocalMeshFieldBase
    real(RP), allocatable :: val(:,:)       !< Buffer to save field data at all nodes
    real(RP), allocatable :: face_val(:,:)  !< Buffer to save field data at face nodes
  end type LocalMeshFieldBase

  type, public :: LocalMeshFieldBaseList
    class(LocalMeshFieldBase), pointer :: ptr
  end type LocalMeshFieldBaseList

  !> Derived type representing a field with 1D local mesh
  type, extends(LocalMeshFieldBase), public :: LocalMeshField1D
    type(LocalMesh1D), pointer :: mesh => null()
  contains
    procedure :: Init => LocalMeshField1D_Init
    procedure :: Final => LocalMeshField1D_Final
  end type LocalMeshField1D

  !> Derived type representing a field with 2D local mesh
  type, extends(LocalMeshFieldBase), public :: LocalMeshField2D
    type(LocalMesh2D), pointer :: mesh => null()
  contains
    procedure :: Init => LocalMeshField2D_Init
    procedure :: Final => LocalMeshField2D_Final
  end type LocalMeshField2D

  !> Derived type representing a field with 3D local mesh
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
  integer, public, parameter :: LOCAL_MESHFIELD_TYPE_NODES_VAL     = 1  !< ID of data type with a field on all nodes in local mesh 
  integer, public, parameter :: LOCAL_MESHFIELD_TYPE_NODES_FACEVAL = 2  !< ID of data type with a field on face nodes in local mesh 

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  private :: LocalMeshFieldBase_Init
  private :: LocalMeshFieldBase_Final

contains

!> Initialize an object to manage a field data (base type)
!OCL SERIAL
  subroutine LocalMeshFieldBase_Init( this, lcmesh, data_type )
    use scale_prc, only: PRC_abort
    implicit none
    class(LocalMeshFieldBase), intent(inout) :: this
    class(LocalMeshBase), intent(in) :: lcmesh   !< Object to manage a field data on a local mesh
    integer, intent(in), optional :: data_type   !< ID of data type (all nodes or face nodes)

    integer :: data_type_
    !-----------------------------------------------------------------------------

    if ( present(data_type) ) then
      data_type_ = data_type
    else
      data_type_ = LOCAL_MESHFIELD_TYPE_NODES_VAL
    end if

    select case( data_type_ )
    case (LOCAL_MESHFIELD_TYPE_NODES_VAL)
      allocate( this%val(lcmesh%refElem%Np,lcmesh%NeA) )
    case (LOCAL_MESHFIELD_TYPE_NODES_FACEVAL)
      allocate( this%face_val(lcmesh%refElem%NfpTot,lcmesh%Ne) )
    case default
      LOG_ERROR("LocalMeshFieldBase_Init",*) "Unexcepted data_type", data_type_
      call PRC_abort        
    end select

    return
  end subroutine LocalMeshFieldBase_Init

!> Finalize an object to manage a field data (base type)
!OCL SERIAL
  subroutine LocalMeshFieldBase_Final( this )
    implicit none
    class(LocalMeshFieldBase), intent(inout) :: this
    !-----------------------------------------------------------------------------

    if ( allocated(this%val) ) deallocate( this%val )
    if ( allocated(this%face_val) ) deallocate( this%face_val )

    return
  end subroutine LocalMeshFieldBase_Final

  !* 1D *********

!> Setup an object to manage a field data on 1D local computational mesh
!OCL SERIAL
  subroutine LocalMeshField1D_Init( this, mesh, data_type )
    implicit none
    class(LocalMeshField1D), intent(inout) :: this
    class(LocalMesh1D), target, intent(in) :: mesh   !< Object to manage a field data on a 1D local mesh
    integer, intent(in), optional :: data_type       !< ID of data type (all nodes or face nodes)
    !-----------------------------------------------------------------------------

    this%mesh => mesh
    call LocalMeshFieldBase_Init( this, mesh, data_type )

    return
  end subroutine LocalMeshField1D_Init

!> Finalize an object to manage a field data on 1D local computational mesh
!OCL SERIAL
  subroutine LocalMeshField1D_Final( this )
    implicit none
    class(LocalMeshField1D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call LocalMeshFieldBase_Final( this )

    return
  end subroutine LocalMeshField1D_Final

  !* 2D *********
  
!> Setup an object to manage a field data on 2D local computational mesh
!OCL SERIAL
  subroutine LocalMeshField2D_Init( this, mesh, data_type )
    implicit none
    class(LocalMeshField2D), intent(inout) :: this
    class(LocalMesh2D), target, intent(in) :: mesh   !< Object to manage a field data on a 2D local mesh
    integer, intent(in), optional :: data_type       !< ID of data type (all nodes or face nodes)
    !-----------------------------------------------------------------------------

    this%mesh => mesh
    call LocalMeshFieldBase_Init( this, mesh, data_type )

    return
  end subroutine LocalMeshField2D_Init

!> Finalize an object to manage a field data on 2D local computational mesh
!OCL SERIAL
  subroutine LocalMeshField2D_Final( this )
    implicit none
    class(LocalMeshField2D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call LocalMeshFieldBase_Final( this )
    return
  end subroutine LocalMeshField2D_Final

  !* 3D *********
  
!> Setup an object to manage a field data on 3D local computational mesh
!OCL SERIAL
  subroutine LocalMeshField3D_Init( this, mesh, data_type )
    implicit none
    class(LocalMeshField3D), intent(inout) :: this
    class(LocalMesh3D), target, intent(in) :: mesh   !< Object to manage a field data on a 3D local mesh
    integer, optional, intent(in) :: data_type
    !-----------------------------------------------------------------------------

    this%mesh => mesh
    call LocalMeshFieldBase_Init( this, mesh, data_type )

    return
  end subroutine LocalMeshField3D_Init

!> Finalize an object to manage a field data on 3D local computational mesh
!OCL SERIAL
  subroutine LocalMeshField3D_Final( this )
    implicit none
    class(LocalMeshField3D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call LocalMeshFieldBase_Final( this )
    return
  end subroutine LocalMeshField3D_Final

end module scale_localmeshfield_base