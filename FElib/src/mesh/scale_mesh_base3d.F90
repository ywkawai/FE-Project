#include "scaleFElib.h"
module scale_mesh_base3d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_localmesh_3d, only: &
    LocalMesh3D, LocalMesh3D_Init, LocalMesh3D_Final

  use scale_mesh_base, only: &
    MeshBase, MeshBase_Init, MeshBase_Final
  use scale_mesh_base2d, only: &
    MeshBase2D
  use scale_element_base, only: elementbase3D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, abstract, public, extends(MeshBase) :: MeshBase3D
    type(LocalMesh3D), allocatable :: lcmesh_list(:)
    type(elementbase3D), pointer :: refElem3D
  contains
    procedure(MeshBase3D_generate), deferred :: Generate 
    procedure(MeshBase3D_getMesh2D), deferred  :: GetMesh2D
    procedure :: GetLocalMesh => MeshBase3D_get_localmesh
  end type MeshBase3D

  interface 
    subroutine MeshBase3D_generate(this)
      import MeshBase3D
      class(MeshBase3D), intent(inout), target :: this
    end subroutine MeshBase3D_generate

    subroutine MeshBase3D_getMesh2D(this, ptr_mesh2D)
      import MeshBase3D
      import MeshBase2D
      class(MeshBase3D), intent(in), target :: this
      class(MeshBase2D), pointer, intent(out) :: ptr_mesh2D
    end subroutine MeshBase3D_getMesh2D
  end interface

  public :: MeshBase3D_Init, MeshBase3D_Final
  public :: MeshBase3D_setGeometricInfo

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
  subroutine MeshBase3D_Init(this,        &
    refElem, NLocalMeshPerPrc, NsideTile, &
    nproc, myrank )
    
    implicit none

    class(MeshBase3D), intent(inout) :: this
    class(elementbase3D), intent(in), target :: refElem
    integer, intent(in) :: NLocalMeshPerPrc
    integer, intent(in) :: NsideTile
    integer, intent(in), optional :: nproc
    integer, intent(in), optional :: myrank

    integer :: n
    !-----------------------------------------------------------------------------
    
    this%refElem3D => refElem
    call MeshBase_Init( this, refElem, NLocalMeshPerPrc, NsideTile, nproc )

    allocate( this%lcmesh_list(this%LOCAL_MESH_NUM) )
    do n=1, this%LOCAL_MESH_NUM
      call LocalMesh3D_Init( this%lcmesh_list(n), refElem, myrank )
    end do

    return
  end subroutine MeshBase3D_Init

  subroutine MeshBase3D_Final( this )
    
    implicit none

    class(MeshBase3D), intent(inout) :: this

    integer :: n
    !-----------------------------------------------------------------------------
  
    do n=1, this%LOCAL_MESH_NUM
      call LocalMesh3D_Final( this%lcmesh_list(n), this%isGenerated )
    end do
    deallocate( this%lcmesh_list )
    
    call MeshBase_Final(this)

    return
  end subroutine MeshBase3D_Final
  
  subroutine MeshBase3D_get_localmesh( this, id, ptr_lcmesh )
    use scale_localmesh_base, only: LocalMeshBase
    implicit none

    class(MeshBase3D), target, intent(in) :: this
    integer, intent(in) :: id
    class(LocalMeshBase), pointer, intent(out) :: ptr_lcmesh
    !-------------------------------------------------------------

    ptr_lcmesh => this%lcmesh_list(id)
    return
  end subroutine MeshBase3D_get_localmesh

  subroutine MeshBase3D_setGeometricInfo( lcmesh, coord_conv, calc_normal )
    implicit none
    
    type(LocalMesh3D), intent(inout) :: lcmesh
    interface
      subroutine coord_conv( x, y, z, xX, xY, xZ, yX, yY, yZ, zX, zY, zZ, &
        vx, vy, vz, elem )
        import elementbase3D
        import RP
        type(elementbase3D), intent(in) :: elem
        real(RP), intent(out) :: x(elem%Np), y(elem%Np), z(elem%Np)
        real(RP), intent(out) :: xX(elem%Np), xY(elem%Np), xZ(elem%Np)
        real(RP), intent(out) :: yX(elem%Np), yY(elem%Np), yZ(elem%Np)
        real(RP), intent(out) :: zX(elem%Np), zY(elem%Np), zZ(elem%Np)
        real(RP), intent(in) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)        
      end subroutine coord_conv
      subroutine calc_normal( normal_fn, &
        Escale_f, fid_h, fid_v, elem )
        import elementbase3D
        import RP
        type(elementbase3D), intent(in) :: elem
        real(RP), intent(out) :: normal_fn(elem%NfpTot,3)
        integer, intent(in) :: fid_h(elem%Nfp_h,elem%Nfaces_h)
        integer, intent(in) :: fid_v(elem%Nfp_v,elem%Nfaces_v)        
        real(RP), intent(in) :: Escale_f(elem%NfpTot,3,3)
      end subroutine calc_normal
    end interface

    class(ElementBase3D), pointer :: refElem
    integer :: n
    integer :: f
    integer :: i, j
    integer :: d
    integer :: fmask(lcmesh%refElem%NfpTot)
    integer :: fid_h(lcmesh%refElem3D%Nfp_h,lcmesh%refElem3D%Nfaces_h)
    integer :: fid_v(lcmesh%refElem3D%Nfp_v,lcmesh%refElem3D%Nfaces_v)
    real(DP) :: Escale_f(lcmesh%refElem%NfpTot,3,3)

    integer :: node_ids(lcmesh%refElem%Nv)
    real(RP) :: vx(lcmesh%refElem%Nv), vy(lcmesh%refElem%Nv), vz(lcmesh%refElem%Nv)    
    real(RP) :: xX(lcmesh%refElem%Np), xY(lcmesh%refElem%Np), xZ(lcmesh%refElem%Np)
    real(RP) :: yX(lcmesh%refElem%Np), yY(lcmesh%refElem%Np), yZ(lcmesh%refElem%Np)
    real(RP) :: zX(lcmesh%refElem%Np), zY(lcmesh%refElem%Np), zZ(lcmesh%refElem%Np)

    !-----------------------------------------------------------------------------

    refElem => lcmesh%refElem3D

    allocate( lcmesh%pos_en(refElem%Np,lcmesh%Ne,3) )
    !allocate( mesh%fx(refElem%Nfaces*refElem%Nfp,mesh%Ne) )
    !allocate( mesh%fy(refElem%Nfaces*refElem%Nfp,mesh%Ne) )
    allocate( lcmesh%normal_fn(refElem%NfpTot,lcmesh%Ne,3) )
    allocate( lcmesh%sJ(refElem%NfpTot,lcmesh%Ne) )
    allocate( lcmesh%J(refElem%Np,lcmesh%Ne) )
    allocate( lcmesh%Fscale(refElem%NfpTot,lcmesh%Ne) )
    allocate( lcmesh%Escale(refElem%Np,lcmesh%Ne,3,3) )
    allocate( lcmesh%zS(refElem%Np,lcmesh%Ne) )
    allocate( lcmesh%Sz(refElem%Np,lcmesh%Ne) )
    allocate( lcmesh%Gsqrt(refElem%Np,lcmesh%Ne) )
    allocate( lcmesh%G_ij(refElem%Np,lcmesh%Ne2D,2,2) )
    allocate( lcmesh%GIJ (refElem%Np,lcmesh%Ne2D,2,2) )
    allocate( lcmesh%lon2D(refElem%Np,lcmesh%Ne2D) )
    allocate( lcmesh%lat2D(refElem%Np,lcmesh%Ne2D) )
    
    do f=1, refElem%Nfaces_h
    do i=1, refElem%Nfp_h
      fid_h(i,f) = i + (f-1)*refElem%Nfp_h
      fmask(fid_h(i,f)) = refElem%Fmask_h(i,f)
    end do
    end do
    do f=1, refElem%Nfaces_v
    do i=1, refElem%Nfp_v
      fid_v(i,f) = i + refElem%Nfaces_h*refElem%Nfp_h + (f-1)*refElem%Nfp_v
      fmask(fid_v(i,f)) = refElem%Fmask_v(i,f)
    end do
    end do

    do n=1, lcmesh%Ne
      node_ids(:) = lcmesh%EToV(n,:)
      vx(:) = lcmesh%pos_ev(node_ids(:),1)
      vy(:) = lcmesh%pos_ev(node_ids(:),2)
      vz(:) = lcmesh%pos_ev(node_ids(:),3)
      call coord_conv( &
        lcmesh%pos_en(:,n,1), lcmesh%pos_en(:,n,2), lcmesh%pos_en(:,n,3), & ! (in)
        xX, xY, xZ, yX, yY, yZ, zX, zY, zZ,                               & ! (out)
        vx, vy, vz, refElem )                                               ! (in)

      lcmesh%J(:,n) =   xX(:)*(yY(:)*zZ(:) - zY(:)*yZ) &
                      - yX(:)*(xY(:)*zZ(:) - zY(:)*xZ) &
                      + zX(:)*(xY(:)*yZ(:) - yY(:)*xZ)

      lcmesh%Escale(:,n,1,1) =   (yY(:)*zZ(:) - zY(:)*yZ(:))/lcmesh%J(:,n)
      lcmesh%Escale(:,n,1,2) = - (xY(:)*zZ(:) - zY(:)*xZ(:))/lcmesh%J(:,n)
      lcmesh%Escale(:,n,1,3) =   (xY(:)*yZ(:) - yY(:)*xZ(:))/lcmesh%J(:,n)
      
      lcmesh%Escale(:,n,2,1) = - (yX(:)*zZ(:) - zX(:)*YZ(:))/lcmesh%J(:,n)
      lcmesh%Escale(:,n,2,2) =   (xX(:)*zZ(:) - zX(:)*xZ(:))/lcmesh%J(:,n)
      lcmesh%Escale(:,n,2,3) = - (xX(:)*yZ(:) - yX(:)*xZ(:))/lcmesh%J(:,n)
      
      lcmesh%Escale(:,n,3,1) =   (yX(:)*zY(:) - zX(:)*yY(:))/lcmesh%J(:,n)
      lcmesh%Escale(:,n,3,2) = - (xX(:)*zY(:) - zX(:)*xY(:))/lcmesh%J(:,n)
      lcmesh%Escale(:,n,3,3) =   (xX(:)*yY(:) - yX(:)*xY(:))/lcmesh%J(:,n)
      
      !* Face

      !
      !mesh%fx(:,n) = mesh%x(fmask(:),n)
      !mesh%fy(:,n) = mesh%y(fmask(:),n)

      ! Calculate normal vectors
      do j=1, 3
      do i=1, 3
        Escale_f(:,i,j) = lcmesh%Escale(fmask(:),n,i,j)
      end do
      end do
      call calc_normal( lcmesh%normal_fn(:,n,:), & ! (out)
        Escale_f, fid_h, fid_v, refElem )          ! (in)

      lcmesh%sJ(:,n) = sqrt( &
        lcmesh%normal_fn(:,n,1)**2 + lcmesh%normal_fn(:,n,2)**2 + lcmesh%normal_fn(:,n,3)**2 )
      do d=1, 3
        lcmesh%normal_fn(:,n,d) = lcmesh%normal_fn(:,n,d)/lcmesh%sJ(:,n)
      end do
      lcmesh%sJ(:,n) = lcmesh%sJ(:,n)*lcmesh%J(fmask(:),n)

      lcmesh%Fscale(:,n) = lcmesh%sJ(:,n)/lcmesh%J(fmask(:),n)
      lcmesh%Gsqrt(:,n) = 1.0_RP
    end do

    return
  end subroutine MeshBase3D_setGeometricInfo
  
end module scale_mesh_base3d