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
    elementbase1D, elementbase2D, elementbase3D


  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  type, abstract, public :: MeshFieldBase
    character(len=H_SHORT) :: varname
    character(len=H_SHORT) :: unit
    integer :: hist_id
    integer :: monitor_id
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
  
  type, extends(MeshFieldBase), public :: MeshField1D
    type(LocalMeshField1D), allocatable :: local(:)
    class(MeshBase1D), pointer :: mesh 
  contains
    procedure :: Init => MeshField1D_Init
    procedure :: Final => MeshField1D_Final
    procedure :: GetLocalMeshField =>  MeshField1D_get_LocalMeshField
  end type MeshField1D

  type, extends(MeshFieldBase), public :: MeshField2D
    type(LocalMeshField2D), allocatable :: local(:)
    class(MeshBase2D), pointer :: mesh   
  contains
    procedure :: Init => MeshField2D_Init
    procedure :: Final => MeshField2D_Final
    procedure :: GetLocalMeshField =>  MeshField2D_get_LocalMeshField
  end type MeshField2D

  type, extends(MeshFieldBase), public :: MeshField3D
    type(LocalMeshField3D), allocatable :: local(:)
    class(MeshBase3D), pointer :: mesh   
  contains
    procedure :: Init => MeshField3D_Init
    procedure :: Final => MeshField3D_Final  
    procedure :: GetLocalMeshField =>  MeshField3D_get_LocalMeshField
  end type MeshField3D

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

!OCL SERIAL
  subroutine MeshField1D_Init( this, varname, units, mesh, data_type )
    implicit none
    class(MeshField1D), intent(inout) :: this
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: units    
    class(MeshBase1D), target, intent(in) :: mesh
    integer, intent(in), optional :: data_type

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

!OCL SERIAL
  subroutine MeshField2D_Init( this, varname, units, mesh, data_type )
    implicit none

    class(MeshField2D), intent(inout) :: this
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: units    
    class(MeshBase2D), target, intent(in) :: mesh
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

!OCL SERIAL
  subroutine MeshField3D_Init( this, varname, units, mesh, data_type )
    implicit none

    class(MeshField3D), intent(inout) :: this
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: units    
    class(MeshBase3D), target, intent(in) :: mesh
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