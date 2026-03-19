!-------------------------------------------------------------------------------
!> module FElib / Mesh / Topography
!!
!! @par Description
!!          A module to manage vertical coordinate with topography 
!!
!! @author Yuta Kawai, Team SCALE
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

  use scale_sparsemat, only: SparseMat

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  !> Derived type to manage topography datat and setup vertical coordinate metric
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
      !$acc update device(this%topo%local(n)%val)
    end do

    return
  end subroutine MeshTopography_Init

  !> Finalize an object of MeshTopography
!OCL SERIAL
  subroutine MeshTopography_Final( this )
    implicit none

    class(MeshTopography), intent(inout) :: this
    !-----------------------------------------------------------------------------
          
    call this%topo%Final()
    return
  end subroutine MeshTopography_Final
  
  !> Setup vertical coordinate metric with topography data
!OCL SERIAL
  subroutine MeshTopography_set_vcoordinate( this, mesh3D, &
      vcoord_id, zTop, comm3D, comm2D                      )
    use scale_meshfieldcomm_cubedspheredom3d, only: MeshFieldCommCubedSphereDom3D
    implicit none

    class(MeshTopography), target, intent(inout) :: this   !< Object of MeshTopography to manage topography data and vertical coordinate metric
    class(MeshBase3D), intent(inout), target :: mesh3D     !< 3D mesh object
    integer, intent(in) :: vcoord_id                       !< Vertical coordinate ID
    real(RP), intent(in) :: zTop                           !< Height of the model domain
    class(MeshFieldCommBase), intent(inout) :: comm3D      !< Object for 3D data communication
    class(MeshFieldCommBase), intent(inout) :: comm2D      !< Object for 2D data communication

    type(SparseMat) :: Dx2D, Dy2D, Lift2D
    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke2D

    type(MeshFieldContainer) :: comm2d_varlist(1)
    type(MeshFieldContainer) :: comm3d_varlist(4)

    type(MeshField3D), target :: zlev, GsqrtV, G13, G23
    type(MeshField3D), target :: tmp_G13, tmp_G23

    logical :: flag_covariantvec
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

    call zlev%Init( "zlev", "", mesh3D )
    call GsqrtV%Init( "GsqrtV", "", mesh3D )
    call G13%Init( "G13", "", mesh3D )
    call G23%Init( "G23", "", mesh3D )

    call tmp_G13%Init( "tmp_G13", "", mesh3D )
    call tmp_G23%Init( "tmp_G23", "", mesh3D )

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      lcmesh2D => lcmesh%lcmesh2D
      
      call calc_vcoordinate_metrics_lc( &
        G13%local(n)%val, G23%local(n)%val, tmp_G13%local(n)%val, tmp_G23%local(n)%val,     & ! (out)
        zlev%local(n)%val, GsqrtV%local(n)%val,                                             & ! (out)
        this%topo%local(n)%val, zTop, vcoord_id,                                            & ! (in)
        lcmesh, lcmesh%refElem3D, lcmesh2D, lcmesh2D%refElem2D,                             & ! (in)
        lcmesh%GIJ(:,:,1,1), lcmesh%GIJ(:,:,1,2), lcmesh%GIJ(:,:,2,1), lcmesh%GIJ(:,:,2,2), & ! (in)
        Dx2D, Dy2D, Lift2D )                                                                  ! (in)
    end do

    ! Exchange metric data to fill halo

    comm3d_varlist(1)%field3d => zlev
    comm3d_varlist(2)%field3d => GsqrtV
    comm3d_varlist(3)%field3d => tmp_G13
    comm3d_varlist(4)%field3d => tmp_G23

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

      if ( flag_covariantvec ) then
        call mesh3D%Set_geometric_with_vcoord( n, &
          GsqrtV%local(n)%val, zlev%local(n)%val, G13%local(n)%val, G23%local(n)%val )
      else
        call mesh3D%Set_geometric_with_vcoord( n, &
          GsqrtV%local(n)%val, zlev%local(n)%val, tmp_G13%local(n)%val, tmp_G23%local(n)%val )
      end if
    end do

    !---
    call zlev%Final()
    call GsqrtV%Final()
    call G13%Final()
    call G23%Final()
    call tmp_G13%Final()
    call tmp_G23%Final()

    call Dx2D%Final()
    call Dy2D%Final()
    call Lift2D%Final()

    return
  end subroutine MeshTopography_set_vcoordinate

!-- prviate ----

!OCL SERIAL
  subroutine calc_vcoordinate_metrics_lc( G13, G23, tmp_G13, tmp_G23, zlev, GsqrtV, &
    topo, zTop, vcoord_id, lcmesh, elem3D, lcmesh2D, elem2D,      &
    G11, G12, G21, G22, Dx2D, Dy2D, Lift2D )
    use scale_element_base, only: ElementBase2D
    use scale_meshutil_vcoord, only: MeshUtil_VCoord_GetMetric
    implicit none
    type(LocalMesh3D), intent(in) :: lcmesh
    type(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase3D), intent(in) :: elem3D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: G13(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: G23(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: tmp_G13(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: tmp_G23(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: zlev(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: GsqrtV(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: topo(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: zTop
    integer, intent(in) :: vcoord_id
    real(RP), intent(in) :: G11(elem2D%Np,lcmesh2D%Ne)
    real(RP), intent(in) :: G12(elem2D%Np,lcmesh2D%Ne)
    real(RP), intent(in) :: G21(elem2D%Np,lcmesh2D%Ne)
    real(RP), intent(in) :: G22(elem2D%Np,lcmesh2D%Ne)
    type(SparseMat), intent(in) :: Dx2D, Dy2D, Lift2D

    integer :: IndexH2Dto3D(elem3D%Np)
    integer :: ke, ke2D, p
    !------------------------------------------------------------------------------

    call MeshUtil_VCoord_GetMetric( &
      G13, G23, zlev, GsqrtV,                                 & ! (out)
      topo, zTop, vcoord_id,                                  & ! (in)
      lcmesh, lcmesh%refElem3D, lcmesh2D, lcmesh2D%refElem2D, & ! (in)
      Dx2D, Dy2D, Lift2D                                      ) ! (in)

    IndexH2Dto3D(:) = elem3D%IndexH2Dto3D(:)
    !$omp parallel do private(ke2D)
    !$acc parallel loop private(ke2D) collapse(2) &
    !$acc   present(G11, G12, G21, G22, G13, G23, tmp_G13, tmp_G23, lcmesh, elem3D) copyin(IndexH2Dto3D)
    do ke=lcmesh%NeS, lcmesh%NeE
    do p=1, elem3D%Np
      ke2D = lcmesh%EMap3Dto2D(ke)
      tmp_G13(p,ke) = G11(IndexH2Dto3D(p),ke2D) * G13(p,ke) + G12(IndexH2Dto3D(p),ke2D) * G23(p,ke)
      tmp_G23(p,ke) = G21(IndexH2Dto3D(p),ke2D) * G13(p,ke) + G22(IndexH2Dto3D(p),ke2D) * G23(p,ke)
    end do
    end do
    return
  end subroutine calc_vcoordinate_metrics_lc

end module scale_mesh_topography