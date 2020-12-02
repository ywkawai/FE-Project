!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics HEVI 
!!
!! @par Description
!!      Modal filter for Atmospheric dynamical process. 
!!      The modal filter surpresses the numerical instability due to the aliasing errors. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_modalfilter
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_modalfilter_Init
  public :: atm_dyn_modalfilter_Final
  public :: atm_dyn_modalfilter_apply

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  real(RP), private, allocatable :: FilterMat(:,:)

contains

  subroutine atm_dyn_modalfilter_Init(  &
    elem,                               &
    etac_h, alpha_h, ord_h,             &
    etac_v, alpha_v, ord_v  )

    use scale_element_line, only: LineElement

    implicit none
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(in) :: etac_h
    real(RP), intent(in) :: alpha_h
    integer, intent(in) :: ord_h
    real(RP), intent(in) :: etac_v
    real(RP), intent(in) :: alpha_v
    integer, intent(in) :: ord_v

    real(RP) :: filter1D_h(elem%Nnode_h1D)
    real(RP) :: filter1D_v(elem%Nnode_v)
    type(LineElement) :: elem1D

    real(RP) :: eta
    integer :: p1, p2, p3
    integer :: l

    !----------------------------------------------------

    filter1D_h(:) = 1.0_RP
    do p1=1, elem%Nnode_h1D
      eta = dble(p1-1)/dble(elem%PolyOrder_h)
      if ( eta >  etac_h .and. p1 /= 1) then
        filter1D_h(p1) = exp( - alpha_h * ( ((eta - etac_h)/(1.0_RP - etac_h))**ord_h ))
      end if
    end do

    filter1D_v(:) = 1.0_RP
    do p3=1, elem%Nnode_v
      eta = dble(p3-1)/dble(elem%PolyOrder_v)
      if ( eta >  etac_v .and. p3 /= 1) then
        filter1D_v(p3) = exp( -  alpha_v * ( ((eta - etac_v)/(1.0_RP - etac_v))**ord_v ))
      end if
    end do

    allocate( FilterMat(elem%Np,elem%Np) )
    FilterMat(:,:) = 0.0_RP
    do p3=1, elem%Nnode_v
    do p2=1, elem%Nnode_h1D
    do p1=1, elem%Nnode_h1D
      l = p1 + (p2-1)*elem%Nnode_h1D + (p3-1)*elem%Nnode_h1D**2
      FilterMat(l,l) = filter1D_h(p1) * filter1D_h(p2) * filter1D_v(p3)
    end do  
    end do
    end do
    FilterMat(:,:) = matmul(FilterMat, elem%invV)
    FilterMat(:,:) = matmul(elem%V, FilterMat)
    
    !--

    return
  end subroutine atm_dyn_modalfilter_Init


  subroutine atm_dyn_modalfilter_Final()
    implicit none
    !--------------------------------------------

    if( allocated(FilterMat) ) deallocate( FilterMat )
    
    return
  end subroutine atm_dyn_modalfilter_Final  

  subroutine atm_dyn_modalfilter_apply( &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, lmesh, elem  )
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem    
    real(RP), intent(inout)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: DRHOT_(elem%Np,lmesh%NeA)
    
    integer :: k
    !------------------------------------

    !$omp parallel do
    do k=1, lmesh%Ne
      DDENS_(:,k) = matmul(FilterMat,DDENS_(:,k))
      MOMX_ (:,k) = matmul(FilterMat,MOMX_(:,k))
      MOMY_ (:,k) = matmul(FilterMat,MOMY_(:,k))
      MOMZ_ (:,k) = matmul(FilterMat,MOMZ_(:,k))
      DRHOT_(:,k) = matmul(FilterMat,DRHOT_(:,k))
    end do
    
    return
  end subroutine atm_dyn_modalfilter_apply

end module scale_atm_dyn_modalfilter
