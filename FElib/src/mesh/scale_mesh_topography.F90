!-------------------------------------------------------------------------------
!> module FElib / Mesh / Topography
!!
!! @par Description
!!          A module to manage vertical coordinate with topography 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
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
  type, public :: MeshTopography
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
!OCL SERIAL
  subroutine MeshTopography_Init( this, varname, mesh )
    implicit none

    class(MeshTopography), intent(inout) :: this
    character(len=*), intent(in) :: varname
    class(MeshBase2D), intent(in), target :: mesh

    integer :: n
    !-----------------------------------------------------------------------------

    call this%topo%Init( varname, "m", mesh )

    do n=1, mesh%LOCAL_MESH_NUM
!$omp parallel workshare
      this%topo%local(n)%val(:,:) = 0.0_RP
!$omp end parallel workshare
    end do

    return
  end subroutine MeshTopography_Init

!OCL SERIAL
  subroutine MeshTopography_Final( this )
    implicit none

    class(MeshTopography), intent(inout) :: this
    !-----------------------------------------------------------------------------
          
    call this%topo%Final()

    return
  end subroutine MeshTopography_Final
  
!OCL SERIAL
  subroutine MeshTopography_set_vcoordinate( this, mesh3D, &
      vcoord_id, zTop, comm3D, comm2D                      )
    
    use scale_sparsemat, only: SparseMat
    use scale_meshutil_vcoord, only: MeshUtil_VCoord_GetMetric
    use scale_meshfieldcomm_cubedspheredom3d, only: MeshFieldCommCubedSphereDom3D
    implicit none

    class(MeshTopography), target, intent(inout) :: this
    class(MeshBase3D), intent(inout), target :: mesh3D
    integer, intent(in) :: vcoord_id
    real(RP), intent(in) :: zTop
    class(MeshFieldCommBase), intent(inout) :: comm3D
    class(MeshFieldCommBase), intent(inout) :: comm2D

    type(SparseMat) :: Dx2D, Dy2D, Lift2D
    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke2D

    type(MeshFieldContainer) :: comm2d_varlist(1)
    type(MeshFieldContainer) :: comm3d_varlist(3)

    type(MeshField3D), target :: Gsqrt, G13, G23
    type(MeshField3D), target :: tmp_G13, tmp_G23

    logical :: flag_covariantvec
    real(RP), allocatable :: G_ij(:,:,:,:)
    !-------------------------------------------------------------

    lcmesh => mesh3D%lcmesh_list(1)
    lcmesh2D => lcmesh%lcmesh2D
    call Dx2D  %Init( lcmesh2D%refElem2D%Dx1  )
    call Dy2D  %Init( lcmesh2D%refElem2D%Dx2  )
    call Lift2D%Init( lcmesh2D%refElem2D%Lift )

    ! Exchange topography data to fill halo

    comm2d_varlist(1)%field2d => this%topo
    call comm2D%Put(comm2d_varlist, 1)
    call comm2D%Exchange()
    call comm2D%Get(comm2d_varlist, 1)

    ! Calculate metric factors associated with general vertical coordinates

    call Gsqrt%Init( "Gsqrt", "", mesh3D )
    call G13%Init( "G13", "", mesh3D )
    call G23%Init( "G23", "", mesh3D )

    call tmp_G13%Init( "tmp_G13", "", mesh3D )
    call tmp_G23%Init( "tmp_G23", "", mesh3D )

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      lcmesh2D => lcmesh%lcmesh2D
      elem3D => lcmesh%refElem3D

      Gsqrt%local(n)%val(:,:) = lcmesh%Gsqrt(:,:)
      call MeshUtil_VCoord_GetMetric( &
        G13%local(n)%val, G23%local(n)%val, lcmesh%zlev(:,:),   & ! (out)
        Gsqrt%local(n)%val,                                     & ! (inout)
        this%topo%local(n)%val(:,:), zTop, vcoord_id,           & ! (in)
        lcmesh, lcmesh%refElem3D, lcmesh2D, lcmesh2D%refElem2D, & ! (in)
        Dx2D, Dy2D, Lift2D                                      ) ! (in)

      !$omp parallel do private(ke2D)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)
        tmp_G13%local(n)%val(:,ke) = &
          lcmesh%GIJ(elem3D%IndexH2Dto3D(:),ke2D,1,1) * G13%local(n)%val(:,ke) &
        + lcmesh%GIJ(elem3D%IndexH2Dto3D(:),ke2D,1,2) * G23%local(n)%val(:,ke)
        tmp_G23%local(n)%val(:,ke) = &
          lcmesh%GIJ(elem3D%IndexH2Dto3D(:),ke2D,2,1) * G13%local(n)%val(:,ke) &
        + lcmesh%GIJ(elem3D%IndexH2Dto3D(:),ke2D,2,2) * G23%local(n)%val(:,ke)
      end do        
    end do

    ! Exchange metric data to fill halo

    comm3d_varlist(1)%field3d => Gsqrt
    comm3d_varlist(2)%field3d => tmp_G13
    comm3d_varlist(3)%field3d => tmp_G23

    flag_covariantvec = .false.
    select type(comm3D)
    type is (MeshFieldCommCubedSphereDom3D)
      call comm3D%SetCovariantVec( 1, G13, G23 )
      flag_covariantvec = .true.
    end select

    call comm3D%Put(comm3d_varlist, 1)
    call comm3D%Exchange()
    call comm3D%Get(comm3d_varlist, 1)

    ! Update metric data managed by local mesh

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      elem3D => lcmesh%refElem3D

      !$omp parallel do
      do ke=lcmesh%NeS, lcmesh%NeA
        lcmesh%Gsqrt(:,ke) = Gsqrt%local(n)%val(:,ke)

        if ( flag_covariantvec ) then
          lcmesh%GI3(:,ke,1) = G13%local(n)%val(:,ke)
          lcmesh%GI3(:,ke,2) = G23%local(n)%val(:,ke)
        else
          lcmesh%GI3(:,ke,1) = tmp_G13%local(n)%val(:,ke)
          lcmesh%GI3(:,ke,2) = tmp_G23%local(n)%val(:,ke)
        end if
      end do
    end do

    !---
    call Gsqrt%Final()
    call G13%Final()
    call G23%Final()
    call tmp_G13%Final()
    call tmp_G23%Final()

    call Dx2D%Final()
    call Dy2D%Final()
    call Lift2D%Final()

    return
  end subroutine MeshTopography_set_vcoordinate

end module scale_mesh_topography