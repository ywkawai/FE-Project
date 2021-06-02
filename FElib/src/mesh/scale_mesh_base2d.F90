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
    procedure :: GetLocalMesh => MeshBase2D_get_localmesh
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
  integer, public :: MeshBase2D_DIMTYPE_NUM   = 4
  integer, public :: MeshBase2D_DIMTYPEID_X   = 1
  integer, public :: MeshBase2D_DIMTYPEID_Y   = 2
  integer, public :: MeshBase2D_DIMTYPEID_XY  = 3
  integer, public :: MeshBase2D_DIMTYPEID_XYT = 4

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
    refElem, NLocalMeshPerPrc,     &
    nprocs, myrank                 )
    
    implicit none

    class(MeshBase2D), intent(inout) :: this
    class(elementbase2D), intent(in), target :: refElem
    integer, intent(in) :: NLocalMeshPerPrc
    integer, intent(in), optional :: nprocs
    integer, intent(in), optional :: myrank

    integer :: n
    !-----------------------------------------------------------------------------
    
    this%refElem2D => refElem
    call MeshBase_Init( this, &
      MeshBase2D_DIMTYPE_NUM, refElem, &
      NLocalMeshPerPrc, 4,             &
      nprocs                           )

    allocate( this%lcmesh_list(this%LOCAL_MESH_NUM) )
    do n=1, this%LOCAL_MESH_NUM
      call LocalMesh2D_Init( this%lcmesh_list(n), n, refElem, myrank )
    end do

    call this%SetDimInfo( MeshBase2D_DIMTYPEID_X, "x", "m", "X-coordinate" )
    call this%SetDimInfo( MeshBase2D_DIMTYPEID_Y, "y", "m", "Y-coordinate" )
    call this%SetDimInfo( MeshBase2D_DIMTYPEID_XY, "xy", "m", "XY-coordinate" )
    call this%SetDimInfo( MeshBase2D_DIMTYPEID_XYT, "xyt", "m", "XY-coordinate" )

    return
  end subroutine MeshBase2D_Init

  subroutine MeshBase2D_Final( this )
    implicit none
    
    class(MeshBase2D), intent(inout) :: this
    integer :: n
    !-----------------------------------------------------------------------------
  
    do n=1, this%LOCAL_MESH_NUM
      call LocalMesh2D_Final( this%lcmesh_list(n), this%isGenerated )
    end do
    deallocate( this%lcmesh_list )
    
    call MeshBase_Final(this)

    return
  end subroutine MeshBase2D_Final
  
  subroutine MeshBase2D_get_localmesh( this, id, ptr_lcmesh )
    use scale_localmesh_base, only: LocalMeshBase
    implicit none

    class(MeshBase2D), target, intent(in) :: this
    integer, intent(in) :: id
    class(LocalMeshBase), pointer, intent(out) :: ptr_lcmesh
    !-------------------------------------------------------------

    ptr_lcmesh => this%lcmesh_list(id)
    return
  end subroutine MeshBase2D_get_localmesh

  subroutine MeshBase2D_setGeometricInfo( lcmesh, coord_conv, calc_normal )

    implicit none
    
    type(LocalMesh2D), intent(inout) :: lcmesh
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
      subroutine calc_normal( normal_fn, &
        Escale_f, fid, elem )
        import elementbase2D
        import RP
        type(elementbase2D), intent(in) :: elem
        real(RP), intent(out) :: normal_fn(elem%NfpTot,2)
        integer, intent(in) :: fid(elem%Nfp,elem%Nfaces)
        real(RP), intent(in) :: Escale_f(elem%NfpTot,2,2)
      end subroutine calc_normal      
    end interface

    class(ElementBase2D), pointer :: refElem
    integer :: n
    integer :: f
    integer :: i, j
    integer :: d
    integer :: node_ids(lcmesh%refElem%Nv)
    real(RP) :: vx(lcmesh%refElem%Nv), vy(lcmesh%refElem%Nv)
    real(RP) :: xr(lcmesh%refElem%Np), xs(lcmesh%refElem%Np)
    real(RP) :: yr(lcmesh%refElem%Np), ys(lcmesh%refElem%Np)
    integer :: fmask(lcmesh%refElem%NfpTot)
    integer :: fid(lcmesh%refElem2D%Nfp,lcmesh%refElem2D%Nfaces)
    real(RP) :: Escale_f(lcmesh%refElem%NfpTot,2,2)

    !-----------------------------------------------------------------------------

    refElem => lcmesh%refElem2D

    allocate( lcmesh%pos_en(refElem%Np,lcmesh%Ne,2) )
    !allocate( mesh%fx(refElem%Nfaces*refElem%Nfp,mesh%Ne) )
    !allocate( mesh%fy(refElem%Nfaces*refElem%Nfp,mesh%Ne) )
    allocate( lcmesh%normal_fn(refElem%NfpTot,lcmesh%Ne,2) )
    allocate( lcmesh%sJ(refElem%NfpTot,lcmesh%Ne) )
    allocate( lcmesh%J(refElem%Np,lcmesh%Ne) )
    allocate( lcmesh%Fscale(refElem%NfpTot,lcmesh%Ne) )
    allocate( lcmesh%Escale(refElem%Np,lcmesh%Ne,2,2) )
    allocate( lcmesh%Gsqrt(refElem%Np,lcmesh%Ne) )
    allocate( lcmesh%G_ij(refElem%Np,lcmesh%Ne, 2, 2) )
    allocate( lcmesh%GIJ (refElem%Np,lcmesh%Ne, 2, 2) )
    allocate( lcmesh%lon(refElem%Np,lcmesh%Ne) )
    allocate( lcmesh%lat(refElem%Np,lcmesh%Ne) )
    
    fmask(:) = reshape(refElem%Fmask, shape(fmask))
    do f=1, refElem%Nfaces
    do i=1, refElem%Nfp
       fid(i,f) = i + (f-1)*refElem%Nfp
    end do
    end do

    do n=1, lcmesh%Ne
      node_ids(:) = lcmesh%EToV(n,:)
      vx(:) = lcmesh%pos_ev(node_ids(:),1)
      vy(:) = lcmesh%pos_ev(node_ids(:),2)
      call coord_conv( &
        lcmesh%pos_en(:,n,1), lcmesh%pos_en(:,n,2), xr, xs, yr, ys, & ! (out)
        vx, vy, refElem )                                             ! (in)
      
      lcmesh%J(:,n) = - xs*yr + xr*ys

      lcmesh%Escale(:,n,1,1) =   ys/lcmesh%J(:,n)
      lcmesh%Escale(:,n,1,2) = - xs/lcmesh%J(:,n)
      lcmesh%Escale(:,n,2,1) = - yr/lcmesh%J(:,n)
      lcmesh%Escale(:,n,2,2) =   xr/lcmesh%J(:,n)
      
      !* Face
      
      !
      !mesh%fx(:,n) = mesh%x(fmask(:),n)
      !mesh%fy(:,n) = mesh%y(fmask(:),n)

      ! Calculate normal vectors
      do j=1, 2
      do i=1, 2
        Escale_f(:,i,j) = lcmesh%Escale(fmask(:),n,i,j)
      end do
      end do
      
      call calc_normal( lcmesh%normal_fn(:,n,:), & ! (out)
        Escale_f, fid, refElem )                   ! (in)
      
      lcmesh%sJ(:,n) = sqrt( lcmesh%normal_fn(:,n,1)**2 + lcmesh%normal_fn(:,n,2)**2 )
      do d=1, 2
        lcmesh%normal_fn(:,n,d) = lcmesh%normal_fn(:,n,d)/lcmesh%sJ(:,n)
      end do
      lcmesh%sJ(:,n) = lcmesh%sJ(:,n)*lcmesh%J(fmask(:),n)

      lcmesh%Fscale(:,n) = lcmesh%sJ(:,n)/lcmesh%J(fmask(:),n)       
      lcmesh%Gsqrt(:,n) = 1.0_RP     
    end do

    return
  end subroutine MeshBase2D_setGeometricInfo
  
end module scale_mesh_base2d