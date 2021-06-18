#include "scaleFElib.h"
module scale_mesh_topography

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_element_base, only: ElementBase3D
  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase, MeshFieldContainer

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, abstract, public :: MeshTopography
    type(MeshField2D) :: topo
  contains
    procedure :: Init => MeshTopography_Init
    procedure :: Final => MeshTopography_Final
    procedure :: SetVCoordinate => MeshTopography_set_vcoordinate
  end type MeshTopography

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
  subroutine MeshTopography_Init( this, mesh )
    implicit none

    class(MeshTopography), intent(inout) :: this
    class(MeshBase2D), intent(in), target :: mesh
    !-----------------------------------------------------------------------------

    call this%topo%Init( "TOPO", "m", mesh )

    return
  end subroutine MeshTopography_Init

  subroutine MeshTopography_Final( this )
    implicit none

    class(MeshTopography), intent(inout) :: this
    !-----------------------------------------------------------------------------
          
    call this%topo%Final()

    return
  end subroutine MeshTopography_Final
  
  subroutine MeshTopography_set_vcoordinate( this, mesh3D, &
      vcoord_name, zTop, comm3D, comm2D                    )
    
    use scale_sparsemat, only: SparseMat
    use scale_meshutil_vcoord, only: MeshUtil_VCoord_GetMetric
    implicit none

    class(MeshTopography), target, intent(in) :: this
    class(MeshBase3D), intent(inout), target :: mesh3D
    character(len=*), intent(in) :: vcoord_name
    real(RP), intent(in) :: zTop
    class(MeshFieldCommBase) :: comm3D
    class(MeshFieldCommBase) :: comm2D

    type(SparseMat) :: Dx2D, Dy2D, Lift2D
    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMesh2D), pointer :: lcmesh2D

    integer :: n
    integer :: ke

    type(MeshFieldContainer) :: comm2d_varlist(1)
    type(MeshFieldContainer) :: comm3d_varlist(3)

    type(MeshField3D), target :: Gsqrt, G13, G23
    !-------------------------------------------------------------

    lcmesh => mesh3D%lcmesh_list(1)
    lcmesh2D => lcmesh%lcmesh2D
    call Dx2D  %Init( lcmesh2D%refElem2D%Dx1  )
    call Dy2D  %Init( lcmesh2D%refElem2D%Dx2  )
    call Lift2D%Init( lcmesh2D%refElem2D%Lift )

    !--
    call Gsqrt%Init( "Gsqrt", "", mesh3D )
    call G13%Init( "G13", "", mesh3D )
    call G23%Init( "G23", "", mesh3D )

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      lcmesh2D => lcmesh%lcmesh2D
  
      Gsqrt%local(n)%val(:,:) = lcmesh%Gsqrt(:,:)
      call MeshUtil_VCoord_GetMetric( &
        G13%local(n)%val, G23%local(n)%val, lcmesh%zlev(:,:),   & ! (out)
        Gsqrt%local(n)%val,                                     & ! (inout)
        this%topo%local(n)%val(:,:), zTop, vcoord_name,         & ! (in)
        lcmesh, lcmesh%refElem3D, lcmesh2D, lcmesh2D%refElem2D, & ! (in)
        Dx2D, Dy2D, Lift2D                                      ) ! (in)
    end do

    !-- Communicate halo data

    ! 2D

    comm2d_varlist(1)%field2d => this%topo
    call comm2D%Put(comm2d_varlist, 1)
    call comm2D%Exchange()
    call comm2D%Get(comm2d_varlist, 1)

    ! 3D 
    comm3d_varlist(1)%field3d => Gsqrt
    comm3d_varlist(2)%field3d => G13
    comm3d_varlist(3)%field3d => G23

    call comm3D%Put(comm3d_varlist, 1)
    call comm3D%Exchange()
    call comm3D%Get(comm3d_varlist, 1)

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      !$omp parallel do
      do ke=lcmesh%NeS, lcmesh%NeA
        lcmesh%Gsqrt(:,ke) = Gsqrt%local(n)%val(:,ke)
        lcmesh%GI3(:,ke,1) = G13%local(n)%val(:,ke)
        lcmesh%GI3(:,ke,2) = G23%local(n)%val(:,ke)
      end do
    end do

    !---
    call Dx2D%Final()
    call Dy2D%Final()
    call Lift2D%Final()

    return
  end subroutine MeshTopography_set_vcoordinate

end module scale_mesh_topography