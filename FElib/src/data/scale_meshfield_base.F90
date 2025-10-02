!> module FElib / Data / base
!!
!! @par Description
!!           A module for managing field data with FEM 
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_meshfield_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_localmesh_1d, only: LocalMesh1D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D

  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase,                                  &
    LocalMeshField1D, LocalMeshField2D, LocalMeshField3D
    
  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D


  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !> Derived type representing a field (base type)
  type, abstract, public :: MeshFieldBase
    character(len=H_SHORT) :: varname  !< Variable name
    character(len=H_SHORT) :: unit     !< Unit of variable
    integer :: hist_id                 !< ID for outputting history data
    integer :: monitor_id              !< ID for monitoring variable
  contains
    procedure(MeshFieldBase_get_LocalMeshField), deferred, public :: GetLocalMeshField
  end type MeshFieldBase

  interface
    subroutine MeshFieldBase_get_LocalMeshField( this, domID, ptr_lcmeshField )
      import MeshFieldBase
      import LocalMeshFieldBase
      class(MeshFieldBase), target, intent(in) :: this
      integer, intent(in) :: domID
      class(LocalMeshFieldBase), pointer, intent(out) :: ptr_lcmeshfield
    end subroutine MeshFieldBase_get_LocalMeshField
  end interface

  !------
  
  !> Derived type representing a field with 1D mesh
  type, extends(MeshFieldBase), public :: MeshField1D
    type(LocalMeshField1D), allocatable :: local(:) !< Array of objects with 1D field data
    class(MeshBase1D), pointer :: mesh              !< Pointer to an object with a 1D computational mesh
  contains
    procedure :: Init => MeshField1D_Init
    procedure :: Final => MeshField1D_Final
    procedure :: GetLocalMeshField =>  MeshField1D_get_LocalMeshField
  end type MeshField1D

  type, public :: MeshField1DList
    class(MeshField1D), pointer :: ptr
  end type MeshField1DList

  !> Derived type representing a field with 2D mesh
  type, extends(MeshFieldBase), public :: MeshField2D
    type(LocalMeshField2D), allocatable :: local(:) !< Array of objects with 2D field data
    class(MeshBase2D), pointer :: mesh              !< Pointer to an object with a 2D computational mesh  
  contains
    procedure :: Init => MeshField2D_Init
    procedure :: Final => MeshField2D_Final
    procedure :: GetLocalMeshField =>  MeshField2D_get_LocalMeshField
  end type MeshField2D

  type, public :: MeshField2DList
    class(MeshField2D), pointer :: ptr
  end type MeshField2DList

  !> Derived type representing a field with 3D mesh
  type, extends(MeshFieldBase), public :: MeshField3D
    type(LocalMeshField3D), allocatable :: local(:) !< Array of objects with 3D field data
    class(MeshBase3D), pointer :: mesh              !< Pointer to an object with a 3D computational mesh   
  contains
    procedure :: Init => MeshField3D_Init
    procedure :: Final => MeshField3D_Final  
    procedure :: GetLocalMeshField =>  MeshField3D_get_LocalMeshField
  end type MeshField3D

  type, public :: MeshField3DList
    class(MeshField3D), pointer :: ptr
  end type MeshField3DList

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
  !**** 1D **********************************************************************

!> Setup an object to manage a field data with a 1D computational mesh
!OCL SERIAL
  subroutine MeshField1D_Init( this, varname, units, mesh, data_type )
    implicit none
    class(MeshField1D), intent(inout) :: this
    character(len=*), intent(in) :: varname       !< Variable name
    character(len=*), intent(in) :: units         !< Unit of variable
    class(MeshBase1D), target, intent(in) :: mesh !< Pointer to an object for a 1D computational mesh
    integer, intent(in), optional :: data_type    !< ID of data type (all nodes or face nodes)

    integer :: n
    !-----------------------------------------------------------------------------
    
    this%varname = varname
    this%unit    = units
    this%mesh => mesh
    this%hist_id    = -1
    this%monitor_id = -1
    
    allocate( this%local(mesh%LOCAL_MESH_NUM) )  
    do n=1, mesh%LOCAL_MESH_NUM
      call this%local(n)%Init( mesh%lcmesh_list(n), data_type )
    end do

    return
  end subroutine MeshField1D_Init

!> Finalize an object to manage a field data with a 1D computational mesh
!OCL SERIAL
  subroutine MeshField1D_Final(this)
    implicit none
    class(MeshField1D), intent(inout) :: this
    
    integer :: n
    !-----------------------------------------------------------------------------

    do n=1, size(this%local)
      call this%local(n)%Final()
    end do
    deallocate( this%local )

    return
  end subroutine MeshField1D_Final

!OCL SERIAL
  subroutine MeshField1D_get_LocalMeshField( this, domID, ptr_lcmeshField )
    implicit none

    class(MeshField1D), target, intent(in) :: this
    integer, intent(in) :: domID
    class(LocalMeshFieldBase), pointer, intent(out) :: ptr_lcmeshfield
    !-----------------------------------------------------------------------------

    ptr_lcmeshfield => this%local(domID)
    return
  end subroutine MeshField1D_get_LocalMeshField

  !**** 2D **********************************************************************

!> Setup an object to manage a field data with a 2D computational mesh
!OCL SERIAL
  subroutine MeshField2D_Init( this, varname, units, mesh, data_type )
    implicit none

    class(MeshField2D), intent(inout) :: this
    character(len=*), intent(in) :: varname       !< Variable name
    character(len=*), intent(in) :: units         !< Unit of variable
    class(MeshBase2D), target, intent(in) :: mesh !< Pointer to an object for a 2D computational mesh
    integer, intent(in), optional :: data_type

    integer :: n
    !-----------------------------------------------------------------------------
    
    this%varname = varname
    this%unit = units
    this%mesh => mesh
    this%hist_id = -1

    allocate( this%local(mesh%LOCAL_MESH_NUM) )  
    do n=1, mesh%LOCAL_MESH_NUM
      call this%local(n)%Init( mesh%lcmesh_list(n), data_type )
    end do

    return
  end subroutine MeshField2D_Init
  
!> Finalize an object to manage a field data with a 2D computational mesh
!OCL SERIAL
  subroutine MeshField2D_Final(this)
    implicit none

    class(MeshField2D), intent(inout) :: this
    
    integer :: n
    !-----------------------------------------------------------------------------

    do n=1, size(this%local)
      call this%local(n)%Final()
    end do
    deallocate( this%local )

    return
  end subroutine MeshField2D_Final

!OCL SERIAL
  subroutine MeshField2D_get_LocalMeshField( this, domID, ptr_lcmeshField )
    implicit none

    class(MeshField2D), target, intent(in) :: this
    integer, intent(in) :: domID
    class(LocalMeshFieldBase), pointer, intent(out) :: ptr_lcmeshfield
    !-----------------------------------------------------------------------------

    ptr_lcmeshfield => this%local(domID)
    return
  end subroutine MeshField2D_get_LocalMeshField

  !**** 3D **********************************************************************

!> Setup an object to manage a field data with a 3D computational mesh
!OCL SERIAL
  subroutine MeshField3D_Init( this, varname, units, mesh, data_type )
    implicit none

    class(MeshField3D), intent(inout) :: this
    character(len=*), intent(in) :: varname       !< Variable name
    character(len=*), intent(in) :: units         !< Unit of variable
    class(MeshBase3D), target, intent(in) :: mesh !< Pointer to an object for a 3D computational mesh
    integer, intent(in), optional :: data_type

    integer :: n
    !-----------------------------------------------------------------------------
    
    this%varname = varname
    this%unit = units
    this%mesh => mesh
    this%hist_id = -1
    
    allocate( this%local(mesh%LOCAL_MESH_NUM) )  
    do n=1, mesh%LOCAL_MESH_NUM
      call this%local(n)%Init( mesh%lcmesh_list(n), data_type )
    end do

    return
  end subroutine MeshField3D_Init
  
!> Finalize an object to manage a field data with a 3D computational mesh
!OCL SERIAL
  subroutine MeshField3D_Final(this)
    implicit none

    class(MeshField3D), intent(inout) :: this
    
    integer :: n
    !-----------------------------------------------------------------------------

    do n=1, size(this%local)
      call this%local(n)%Final()
    end do
    deallocate( this%local )

    return
  end subroutine MeshField3D_Final

  subroutine MeshField3D_get_LocalMeshField( this, domID, ptr_lcmeshField )
    implicit none

    class(MeshField3D), target, intent(in) :: this
    integer, intent(in) :: domID
    class(LocalMeshFieldBase), pointer, intent(out) :: ptr_lcmeshfield
    !-----------------------------------------------------------------------------

    ptr_lcmeshfield => this%local(domID)
    return
  end subroutine MeshField3D_get_LocalMeshField

end module scale_meshfield_base