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

  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  
  use scale_localmeshfield_base, only: &
    LocalMeshField1D, LocalMeshField2D
    
  use scale_element_base, only: &
    elementbase1D, elementbase2D


  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  type, public :: MeshFieldBase
    character(len=H_SHORT) :: varname
    character(len=H_SHORT) :: unit
  end type MeshFieldBase
  
  type, extends(MeshFieldBase), public :: MeshField1D
    type(LocalMeshField1D), allocatable :: local(:)
    class(MeshBase1D), pointer :: mesh 
  contains
    procedure :: Init => MeshField1D_Init
    procedure :: Final => MeshField1D_Final  
  end type MeshField1D

  type, extends(MeshFieldBase), public :: MeshField2D
    type(LocalMeshField2D), allocatable :: local(:)
    class(MeshBase2D), pointer :: mesh   
  contains
    procedure :: Init => MeshField2D_Init
    procedure :: Final => MeshField2D_Final  
  end type MeshField2D

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
  subroutine MeshField1D_Init(this, varname, units, mesh)
    class(MeshField1D), intent(inout) :: this
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: units    
    class(MeshBase1D), target, intent(in) :: mesh

    integer :: n
    !-----------------------------------------------------------------------------
    
    this%varname = varname
    this%unit = units
    this%mesh => mesh
    
    allocate( this%local(mesh%LOCAL_MESH_NUM) )  
    do n=1, mesh%LOCAL_MESH_NUM
      call this%local(n)%Init( mesh%lcmesh_list(n) )
    end do

  end subroutine MeshField1D_Init

  subroutine MeshField1D_Final(this)
    class(MeshField1D), intent(inout) :: this
    
    integer :: n
    !-----------------------------------------------------------------------------

    do n=1, size(this%local)
      call this%local(n)%Final()
    end do
    deallocate( this%local )

  end subroutine MeshField1D_Final

  !*******

  subroutine MeshField2D_Init(this, varname, units, mesh)
    class(MeshField2D), intent(inout) :: this
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: units    
    class(MeshBase2D), target, intent(in) :: mesh

    integer :: n
    !-----------------------------------------------------------------------------
    
    this%varname = varname
    this%unit = units
    this%mesh => mesh

    allocate( this%local(mesh%LOCAL_MESH_NUM) )  
    do n=1, mesh%LOCAL_MESH_NUM
      call this%local(n)%Init( mesh%lcmesh_list(n) )
    end do

  end subroutine MeshField2D_Init
  
  subroutine MeshField2D_Final(this)
    class(MeshField2D), intent(inout) :: this
    
    integer :: n
    !-----------------------------------------------------------------------------

    do n=1, size(this%local)
      call this%local(n)%Final()
    end do
    deallocate( this%local )

  end subroutine MeshField2D_Final

end module scale_meshfield_base