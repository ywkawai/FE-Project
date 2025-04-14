!-------------------------------------------------------------------------------
!> module FElib / Mesh / Local 3D
!!
!! @par Description
!!      Module to manage 3D local mesh for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
#include "scaleFElib.h"
module scale_localmesh_3d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_localmesh_base, only: &
    LocalMeshBase, LocalMeshBase_Init, LocalMeshBase_Final
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_element_base, only: ElementBase, ElementBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(LocalMeshBase), public :: LocalMesh3D

    type(ElementBase3D), pointer :: refElem3D
    real(RP) :: xmin, xmax
    real(RP) :: ymin, ymax
    real(RP) :: zmin, zmax
    integer :: NeX
    integer :: NeY
    integer :: NeZ
    integer :: Ne2D
    integer :: Ne2DA

    real(RP), allocatable :: Sz(:,:)
    real(RP), allocatable :: zS(:,:)
    real(RP), allocatable :: GI3(:,:,:)    !< The contravariant component of metric tensor with vertical general coordinate 
    real(RP), allocatable :: GsqrtH(:,:)   !< The Jacobian of horizontal transformation in the computational coordinate
    real(RP), allocatable :: zlev(:,:)
    real(RP), allocatable :: gam(:,:)      !< Factor for approximation with spherical shell domain (= r/a)

    class(LocalMesh2D), pointer :: lcmesh2D
    real(RP), allocatable :: lon2D(:,:)     
    real(RP), allocatable :: lat2D(:,:)

    integer, allocatable :: EMap3Dto2D(:)  
  contains
    procedure :: SetLocalMesh2D => LocalMesh3D_setLocalMesh2D 
    procedure :: GetVmapZ1D => LocalMesh3D_getVmapZ1D 
    procedure :: GetVmapZ3D => LocalMesh3D_getVmapZ3D
  end type LocalMesh3D

  public :: LocalMesh3D_Init, LocalMesh3D_Final

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
  subroutine LocalMesh3D_Init( this, &
    lcdomID, refElem, myrank )
    implicit none

    class(LocalMesh3D), intent(inout) :: this
    integer, intent(in) :: lcdomID
    class(ElementBase3D), intent(in), target :: refElem
    integer, intent(in), optional :: myrank
    !-------------------------------------------------

    this%refElem3D    => refElem
    nullify( this%lcmesh2D )

    call LocalMeshBase_Init( this, lcdomID, refElem, 3, myrank )

    return
  end subroutine LocalMesh3D_Init

!OCL SERIAL
  subroutine LocalMesh3D_Final( this, is_generated )
    implicit none
    type(LocalMesh3D), intent(inout) :: this
    logical, intent(in) :: is_generated
    !-------------------------------------------------

    call LocalMeshBase_Final( this, is_generated )
    if (is_generated) then
      deallocate( this%zS, this%Sz )
      deallocate( this%GI3, this%GsqrtH )
      deallocate( this%zlev )
      deallocate( this%gam )
      deallocate( this%lon2D, this%lat2D )
      deallocate( this%EMap3Dto2D )
    end if

    return
  end subroutine LocalMesh3D_Final
  
!OCL SERIAL
  subroutine LocalMesh3D_setLocalMesh2D( this, lcmesh2D )
    implicit none
    class(LocalMesh3D), intent(inout) :: this
    class(LocalMesh2D), intent(in), target :: lcmesh2D
    !-------------------------------------------------

    this%lcmesh2D => lcmesh2D

    return
  end subroutine LocalMesh3D_setLocalMesh2D

!OCL SERIAL
  subroutine LocalMesh3D_getVmapZ1D( this, vmapM, vmapP )

    implicit none
    class(LocalMesh3D), intent(in), target :: this
    integer, intent(out) :: vmapM(this%refElem3D%NfpTot,this%NeZ)
    integer, intent(out) :: vmapP(this%refElem3D%NfpTot,this%NeZ)    

    integer :: ke_z
    integer :: f
    integer :: vs, ve

    class(ElementBase3D), pointer :: elem
    !------------------------------

    elem => this%refElem3D

    do ke_z=1, this%NeZ
      do f=1, elem%Nfaces_h
        vs = 1 + (f-1)*elem%Nfp_h
        ve = vs + elem%Nfp_h - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_h(:,f) + (ke_z-1)*elem%Np
      end do
      do f=1, elem%Nfaces_v
        vs = elem%Nfp_h*elem%Nfaces_h + 1 + (f-1)*elem%Nfp_v
        ve = vs + elem%Nfp_v - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_v(:,f) + (ke_z-1)*elem%Np
      end do
      vmapP(:,ke_z) = vmapM(:,ke_z)
    end do

    do ke_z=1, this%NeZ
      vs = elem%Nfp_h*elem%Nfaces_h + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z > 1) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,2) + (ke_z-2)*elem%Np

      vs = elem%Nfp_h*elem%Nfaces_h + elem%Nfp_v + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z < this%NeZ) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1) + ke_z*elem%Np
    end do

    return
  end subroutine LocalMesh3D_getVmapZ1D

!OCL SERIAL
  subroutine LocalMesh3D_getVmapZ3D( this, vmapM, vmapP )

    implicit none
    class(LocalMesh3D), intent(in), target :: this
    integer, intent(out) :: vmapM(this%refElem3D%NfpTot,this%Ne)
    integer, intent(out) :: vmapP(this%refElem3D%NfpTot,this%Ne)    

    integer :: ke, ke_xy, ke_z
    integer :: ke_neigh
    integer :: f
    integer :: vs, ve

    class(ElementBase3D), pointer :: elem
    !------------------------------

    elem => this%refElem3D

    !$omp parallel do collapse(2) private(ke,ke_xy,ke_z,ke_neigh,f,vs,ve)
    do ke_z=1, this%NeZ
    do ke_xy=1, this%NeX * this%NeY
      ke = ke_xy + (ke_z-1) * this%NeX * this%NeY

      do f=1, elem%Nfaces_h
        vs = 1 + (f-1)*elem%Nfp_h
        ve = vs + elem%Nfp_h - 1
        vmapM(vs:ve,ke) = elem%Fmask_h(:,f) + (ke-1)*elem%Np
      end do
      do f=1, elem%Nfaces_v
        vs = elem%Nfp_h*elem%Nfaces_h + 1 + (f-1)*elem%Nfp_v
        ve = vs + elem%Nfp_v - 1
        vmapM(vs:ve,ke) = elem%Fmask_v(:,f) + (ke-1)*elem%Np
      end do
      vmapP(:,ke) = vmapM(:,ke)

      !--
      vs = elem%Nfp_h*elem%Nfaces_h + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z > 1) then
        ke_neigh = ke_xy + (ke_z-2) * this%NeX * this%NeY
        vmapP(vs:ve,ke) = elem%Fmask_v(:,2) + (ke_neigh-1)*elem%Np
      end if
      vs = elem%Nfp_h*elem%Nfaces_h + elem%Nfp_v + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z < this%NeZ) then
        ke_neigh = ke_xy + ke_z * this%NeX * this%NeY
        vmapP(vs:ve,ke) = elem%Fmask_v(:,1) + (ke_neigh-1)*elem%Np
      end if
    end do
    end do

    return
  end subroutine LocalMesh3D_getVmapZ3D

end module scale_localmesh_3d