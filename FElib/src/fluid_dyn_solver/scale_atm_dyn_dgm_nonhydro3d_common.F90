!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!      HEVE DGM scheme for Atmospheric dynamical process. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,  &
    Rdry => CONST_Rdry,  &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00

  use scale_sparsemat  
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base2d, only: MeshBase2D
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
  public :: atm_dyn_dgm_nonhydro3d_common_Init
  public :: atm_dyn_dgm_nonhydro3d_common_Final

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  integer, public, parameter :: DENS_VID       = 1
  integer, public, parameter :: MOMX_VID       = 2
  integer, public, parameter :: MOMY_VID       = 3
  integer, public, parameter :: MOMZ_VID       = 4
  integer, public, parameter :: RHOT_VID       = 5
  integer, public, parameter :: PROG_VARS_NUM  = 5

  real(RP), public, allocatable :: IntrpMat_VPOrdM1(:,:)
  integer, public, allocatable :: iM2Dto3D(:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------
  
contains
  subroutine atm_dyn_dgm_nonhydro3d_common_Init( mesh )

    implicit none
    class(MeshBase3D), intent(in) :: mesh

    integer :: p1, p2, p_
    type(ElementBase3D), pointer :: elem
    real(RP) :: invV_VPOrdM1(mesh%refElem3D%Np,mesh%refElem3D%Np)

    integer :: f_h, f_v
    integer :: fp, fp_h1, fp_h2, fp_v
    type(ElementBase2D), pointer :: elem2D
    class(MeshBase2D), pointer :: mesh2D    
    !--------------------------------------------

    elem => mesh%refElem3D
    allocate( IntrpMat_VPOrdM1(elem%Np,elem%Np) )
    
    InvV_VPOrdM1(:,:) = elem%invV
    do p2=1, elem%Nnode_h1D
    do p1=1, elem%Nnode_h1D
      p_ = p1 + (p2-1)*elem%Nnode_h1D + (elem%Nnode_v-1)*elem%Nnode_h1D**2
      InvV_VPOrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem%V, invV_VPOrdM1)

    !--

    allocate( iM2Dto3D(elem%NfpTot) )
    call mesh%GetMesh2D( mesh2D ) 
    elem2D => mesh2D%refElem2D

    do f_h=1, 4
      do fp_v=1, elem%Nnode_v
      do fp_h1=1, elem%Nnode_h1D
        fp = fp_h1 + (fp_v-1)*elem%Nnode_h1D + (f_h-1)*elem%Nfp_h
        iM2Dto3D(fp) = elem2D%Fmask(fp_h1,f_h)
      end do  
      end do
    end do
    do f_v=1, 2
      do fp_h2=1, elem%Nnode_h1D
      do fp_h1=1, elem%Nnode_h1D
        fp = fp_h1 + (fp_h2-1)*elem%Nnode_h1D    &
           + (f_v-1) * elem%Nfp_v                &
           + 4 * elem%Nnode_h1D * elem%Nnode_v
        iM2Dto3D(fp) = fp_h1 + (fp_h2-1)*elem%Nnode_h1D
      end do  
      end do
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_common_Init


  subroutine atm_dyn_dgm_nonhydro3d_common_Final()
    implicit none
    !--------------------------------------------
    
    deallocate( IntrpMat_VPOrdM1 )
    deallocate( iM2Dto3D )
    
    return
  end subroutine atm_dyn_dgm_nonhydro3d_common_Final  

end module scale_atm_dyn_dgm_nonhydro3d_common
