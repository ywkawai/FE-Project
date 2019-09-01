#include "scaleFElib.h"
module scale_mesh_base2d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_localmesh_2d, only: &
    LocalMesh2D, LocalMesh2D_Init, LocalMesh2D_Final

  use scale_mesh_base, only: &
    MeshBase, MeshBase_Init, MeshBase_Final
    
  use scale_element_base, only: elementbase2D
  use scale_element_quadrial, only: QuadrialElement

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, abstract, public, extends(MeshBase) :: MeshBase2D
    type(LocalMesh2D), allocatable :: lcmesh_list(:)
    type(elementbase2D), pointer :: refElem2D
  contains
    procedure(MeshBase2D_generate), deferred :: Generate 
  end type MeshBase2D

  interface 
    subroutine MeshBase2D_generate(this)
      import MeshBase2D
      class(MeshBase2D), intent(inout), target :: this
    end subroutine MeshBase2D_generate
  end interface

  public :: MeshBase2D_Init, MeshBase2D_Final
  public :: MeshBase2D_setGeometricInfo

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
  subroutine MeshBase2D_Init(this, &
    & refElem, NLocalMeshPerPrc )
    
    use scale_prc, only: PRC_nprocs, PRC_myrank
    implicit none

    class(MeshBase2D), intent(inout) :: this
    class(elementbase2D), intent(in), target :: refElem
    integer, intent(in) :: NLocalMeshPerPrc

    integer :: n

    !-----------------------------------------------------------------------------
    
    this%refElem2D => refElem
    call MeshBase_Init(this, refElem, NLocalMeshPerPrc, 4)

    allocate( this%lcmesh_list(this%LOCAL_MESH_NUM) )
    do n=1, this%LOCAL_MESH_NUM
      call LocalMesh2D_Init( this%lcmesh_list(n), refElem, PRC_myrank )
    end do
  end subroutine MeshBase2D_Init

  subroutine MeshBase2D_Final( this )
    
    class(MeshBase2D), intent(inout) :: this

    integer :: n

    !-----------------------------------------------------------------------------
  
    do n=1, this%LOCAL_MESH_NUM
      call LocalMesh2D_Final( this%lcmesh_list(n) )
    end do
    deallocate( this%lcmesh_list )
    
    call MeshBase_Final(this)

  end subroutine MeshBase2D_Final
  
  subroutine MeshBase2D_setGeometricInfo( mesh, coord_conv )

    implicit none
    
    type(LocalMesh2D), intent(inout) :: mesh
    interface
      subroutine coord_conv( x, y, xr, xs, yr, ys, &
        vx, vy, elem )
        import elementbase2D
        import RP
        type(elementbase2D), intent(in) :: elem
        real(RP), intent(out) :: x(elem%Np), y(elem%Np)
        real(RP), intent(out) :: xr(elem%Np), xs(elem%Np), yr(elem%Np), ys(elem%Np)
        real(RP), intent(in) :: vx(elem%Nv), vy(elem%Nv)        
      end subroutine coord_conv
    end interface

    class(ElementBase2D), pointer :: refElem
    integer :: n
    integer :: f
    integer :: i
    real(RP) :: vx(mesh%refElem%Nv), vy(mesh%refElem%Nv)
    real(RP) :: xr(mesh%refElem%Np), xs(mesh%refElem%Np)
    real(RP) :: yr(mesh%refElem%Np), ys(mesh%refElem%Np)
    real(DP) :: Escale(2,2,mesh%refElem%Np)
    integer :: fmask(mesh%refElem%NfpTot)
    integer :: fid(mesh%refElem2D%Nfp,mesh%refElem2D%Nfaces)
    real(DP) :: fxr(mesh%refElem%NfpTot), fxs(mesh%refElem%NfpTot)
    real(DP) :: fyr(mesh%refElem%NfpTot), fys(mesh%refElem%NfpTot)

  !-----------------------------------------------------------------------------

    refElem => mesh%refElem2D

    allocate( mesh%pos_en(refElem%Np,mesh%Ne,2) )
    !allocate( mesh%fx(refElem%Nfaces*refElem%Nfp,mesh%Ne) )
    !allocate( mesh%fy(refElem%Nfaces*refElem%Nfp,mesh%Ne) )
    allocate( mesh%normal_fn(refElem%NfpTot,mesh%Ne,2) )
    allocate( mesh%sJ(refElem%NfpTot,mesh%Ne) )
    allocate( mesh%J(refElem%Np,mesh%Ne) )
    allocate( mesh%Fscale(refElem%NfpTot,mesh%Ne) )
    allocate( mesh%Escale(refElem%Np,mesh%Ne,2,2) )
    allocate( mesh%Gsqrt(refElem%Np,mesh%Ne) )
    allocate( mesh%G_ij(refElem%Np,mesh%Ne, 2, 2) )
    allocate( mesh%GIJ (refElem%Np,mesh%Ne, 2, 2) )
    allocate( mesh%lon(refElem%Np,mesh%Ne) )
    allocate( mesh%lat(refElem%Np,mesh%Ne) )
    
    fmask(:) = reshape(refElem%Fmask, shape(fmask))
    do f=1, refElem%Nfaces
    do i=1, refElem%Nfp
       fid(i,f) = i + (f-1)*refElem%Nfp
    end do
    end do

    do n=1, mesh%Ne
       vx(:) = mesh%pos_ev(mesh%EToV(n,:),1)
       vy(:) = mesh%pos_ev(mesh%EToV(n,:),2)
       call coord_conv( &
        mesh%pos_en(:,n,1), mesh%pos_en(:,n,2), xr, xs, yr, ys, & ! (out)
        vx, vy, refElem )                                         ! (in)

       mesh%J(:,n) = - xs*yr + xr*ys
       mesh%Escale(:,n,1,1) =   ys/mesh%J(:,n)
       mesh%Escale(:,n,1,2) = - xs/mesh%J(:,n)
       mesh%Escale(:,n,2,1) = - yr/mesh%J(:,n)
       mesh%Escale(:,n,2,2) =   xr/mesh%J(:,n)

       !* Face

       !
       !mesh%fx(:,n) = mesh%x(fmask(:),n)
       !mesh%fy(:,n) = mesh%y(fmask(:),n)

       ! Calculate normal vectors
       fxr(:) = xr(fmask(:))
       fxs(:) = xs(fmask(:))
       fyr(:) = yr(fmask(:))
       fys(:) = ys(fmask(:))
       mesh%normal_fn(fid(:,1),n,1) = +fyr(fid(:,1))
       mesh%normal_fn(fid(:,1),n,2) = -fxr(fid(:,1))
       mesh%normal_fn(fid(:,2),n,1) = +fys(fid(:,2))
       mesh%normal_fn(fid(:,2),n,2) = +fxs(fid(:,2))
       mesh%normal_fn(fid(:,3),n,1) = +fyr(fid(:,3))
       mesh%normal_fn(fid(:,3),n,2) = +fxr(fid(:,3))
       mesh%normal_fn(fid(:,4),n,1) = -fys(fid(:,4))
       mesh%normal_fn(fid(:,4),n,2) = +fxs(fid(:,4))

       mesh%sJ(:,n) = sqrt(mesh%normal_fn(:,n,1)**2 + mesh%normal_fn(:,n,2)**2)
       mesh%normal_fn(:,n,1) = mesh%normal_fn(:,n,1)/mesh%sJ(:,n)
       mesh%normal_fn(:,n,2) = mesh%normal_fn(:,n,2)/mesh%sJ(:,n)
       mesh%Fscale(:,n) = mesh%sJ(:,n)/mesh%J(fmask(:),n)
       mesh%Gsqrt(:,n) = 1d0       
    end do

  end subroutine MeshBase2D_setGeometricInfo
  
end module scale_mesh_base2d