!-------------------------------------------------------------------------------
!> module FElib / Mesh / utility for general vertical coordinate
!!
!! @par Description
!!          A module useful for general vertical coordinate
!!
!! @author Yuta Kawai, Team SCALE
!<
#include "scaleFElib.h"
module scale_meshutil_vcoord
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_const, only: &
    PI => CONST_PI,      &
    EPS => CONST_EPS
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prc

  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_sparsemat, only:&
    SparseMat, sparsemat_matmul
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  public :: MeshUtil_get_VCoord_TypeID
  public :: MeshUtil_VCoord_GetMetric

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(*), public, parameter :: MESH_VCOORD_TERRAIN_FOLLOWING_NAME = "TERRAIN_FOLLOWING"
  integer, public, parameter :: MESH_VCOORD_TERRAIN_FOLLOWING_ID        = 1

  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !  
  !-----------------------------------------------------------------------------

contains
!OCL SERIAL
  function MeshUtil_get_VCoord_TypeID( vcoord_type ) result( vcoord_id )
    implicit none

    character(len=*), intent(in) :: vcoord_type
    integer :: vcoord_id

    select case( vcoord_type )
    case( MESH_VCOORD_TERRAIN_FOLLOWING_NAME )
      vcoord_id = MESH_VCOORD_TERRAIN_FOLLOWING_ID
    case default
      LOG_ERROR("MeshUtil_VCoord_TypeID",*) "vcoord_type is inappropriate. Check!", vcoord_type
      call PRC_abort
    end select

    return
  end function MeshUtil_get_VCoord_TypeID

!OCL SERIAL
  subroutine MeshUtil_VCoord_GetMetric( G13, G23, zlev, GsqrtV, &
    topo, zTop, vcoord_id, lcmesh, elem, lcmesh2D, elem2D,      &
    Dx2D, Dy2D, Lift2D )
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem   
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: G13(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: G23(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: zlev(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: GsqrtV(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: topo(elem2D%Np,lcmesh2D%NeA)
    integer, intent(in) :: vcoord_id
    real(RP), intent(in) :: zTop
    type(SparseMat), intent(in) :: Dx2D
    type(SparseMat), intent(in) :: Dy2D
    type(SparseMat), intent(in) :: Lift2D

    integer :: ke, ke2D
    real(RP) :: del_flux(elem2D%NfpTot,lcmesh2D%Ne,2)
    real(RP) :: Fx2D(elem2D%Np), Fy2D(elem2D%Np), LiftDelFlux2D(elem2D%Np)
    real(RP) :: GradZs(elem2D%Np,lcmesh2D%Ne,2)
    real(RP) :: coef3D(elem%NP)
    !------------------------------------------------

    if ( vcoord_id == MESH_VCOORD_TERRAIN_FOLLOWING_ID ) then
      
      ! * z    = topo + (1 - topo / zTop ) * zeta
      ! * zeta = zTop * (z - topo)/(zTop - topo)
      ! * Gi3  = (dzeta(x1,x2,z)/dxi)_z = d (zeta,z) / d (xi,z) = - d(z,zeta)/ d(xi,zeta) *  d(xi,zeta)/d(xi,z)
      !                                 = - (dz/dxi)_zeta * dzeta/dz
      !                                 = (GsqrtV)^-1 * [ - 1 + zeta / zTop ] * d topo /dxi     (i=1, 2)

      call cal_del_flux( del_flux, &
        topo, lcmesh2D%normal_fn(:,:,1), lcmesh2D%normal_fn(:,:,2), &
        lcmesh2D%VMapM, lcmesh2D%VMapP, lcmesh2D, elem2D            )
      
      !$omp parallel private(ke2D, ke,        &
      !$omp Fx2D, Fy2D, LiftDelFlux2D, coef3D )

      !$omp do
      do ke2D=1, lcmesh2D%Ne
        call sparsemat_matmul( Dx2D, topo(:,ke2D), Fx2D )
        call sparsemat_matmul( Lift2D, lcmesh2D%Fscale(:,ke2D) * del_flux(:,ke2D,1), LiftDelFlux2D)
        GradZs(:,ke2D,1) = lcmesh2D%Escale(:,ke2D,1,1) * Fx2D(:) + LiftDelFlux2D(:)

        call sparsemat_matmul( Dy2D, topo(:,ke2D), Fy2D )
        call sparsemat_matmul( Lift2D, lcmesh2D%Fscale(:,ke2D) * del_flux(:,ke2D,2), LiftDelFlux2D)
        GradZs(:,ke2D,2) = lcmesh2D%Escale(:,ke2D,2,2) * Fy2D(:) + LiftDelFlux2D(:)
      end do
      !$omp end do

      !$omp do
      do ke=1, lcmesh%Ne
        ke2D = lcmesh%EMap3Dto2D(ke)
        coef3D(:) = 1.0_RP - lcmesh%pos_en(:,ke,3) / zTop
        zlev(:,ke) = lcmesh%pos_en(:,ke,3)                  &
                   + coef3D(:) * topo(elem%IndexH2Dto3D,ke2D)

        GsqrtV(:,ke) = 1.0_RP - topo(elem%IndexH2Dto3D,ke2D) / zTop ! dz/dzeta
        coef3D(:) = - coef3D(:) / GsqrtV(:,ke)
        G13(:,ke) = coef3D(:) * GradZs(elem%IndexH2Dto3D(:),ke2D,1)
        G23(:,ke) = coef3D(:) * GradZs(elem%IndexH2Dto3D(:),ke2D,2)
      end do
      !$omp end do
      !$omp end parallel
    else
      LOG_ERROR("Mesh_VCoord_GetMetric",*) "vcoord_id is inappropriate. Check!", vcoord_id
      call PRC_abort        
    end if

    return
  end subroutine MeshUtil_VCoord_GetMetric

!OCL SERIAL
  subroutine cal_del_flux( del_flux, &
    topo, nx, ny, vmapM, vmapP, lmesh, elem )
    implicit none
    type(LocalMesh2D), intent(in) :: lmesh
    type(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: del_flux(elem%NfpTot*lmesh%Ne,2)
    real(RP), intent(in) :: topo(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)

    integer :: i
    integer :: iP, iM
    real(RP) :: dtopo
    !-------------------------------------

    !$omp parallel do private( i, iM, iP, dtopo )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      dtopo = 0.5_RP * ( topo(iP) - topo(iM) )

      del_flux(i,1) = dtopo * nx(i)
      del_flux(i,2) = dtopo * ny(i)
    end do

    return
  end subroutine cal_del_flux
end module scale_meshutil_vcoord