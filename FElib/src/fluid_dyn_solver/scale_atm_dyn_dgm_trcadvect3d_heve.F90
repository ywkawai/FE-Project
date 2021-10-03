!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics HEVE
!!
!! @par Description
!!      HEVE DGM scheme for tracer advection. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_trcadvect3d_heve
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

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
  public :: atm_dyn_dgm_trcadvect3d_heve_Init
  public :: atm_dyn_dgm_trcadvect3d_heve_Final
  public :: atm_dyn_dgm_trcadvect3d_heve_cal_tend
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  private :: atm_dyn_dgm_trcadvect3d_heve_numflux_get_generalhvc

contains

!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_heve_Init(mesh)
    implicit none
    class(MeshBase3D), intent(in) :: mesh
    !--------------------------------------------

    return
  end subroutine atm_dyn_dgm_trcadvect3d_heve_Init
  
!OCL SERIAL  
  subroutine atm_dyn_dgm_trcadvect3d_heve_Final()
    implicit none
    !--------------------------------------------

    return    
  end subroutine atm_dyn_dgm_trcadvect3d_heve_Final

!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_heve_cal_tend( &
    QTRC_dt,                                                   & ! (out)
    QTRC_, MOMX_, MOMY_, MOMZ_,                                & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift   
    
    real(RP), intent(out) :: QTRC_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: QTRC_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lmesh%NeA)    
    real(RP), intent(in) :: MOMZ_(elem%Np,lmesh%NeA)    

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)    
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
    real(RP) :: momwt_(elem%Np)    
    
    integer :: ke, ke2d    
    !---------------------------------------------------------------------------

    call PROF_rapstart('cal_trcadv_tend_bndflux', 3)
    call atm_dyn_dgm_trcadvect3d_heve_numflux_get_generalhvc( &
      del_flux,                                                                & ! (out)
      QTRC_, MOMX_, MOMY_, MOMZ_,                                              & ! (in)
      lmesh%Gsqrt, lmesh%GsqrtH, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),           & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%vmapM, lmesh%vmapP, elem%IndexH2Dto3D_bnd,                         & ! (in)
      lmesh, elem, lmesh2D, elem2D                                             ) ! (in)
    call PROF_rapend('cal_trcadv_tend_bndflux', 3)

    call PROF_rapstart('cal_trcadv_tend_interior', 3)
    !$omp parallel do private( &
    !$omp momwt_, Fx, Fy, Fz, LiftDelFlx,  &    
    !$omp ke, ke2d                         )  
    do ke=lmesh%NeS, lmesh%NeE
      ke2d = lmesh%EMap3Dto2D(ke)
      
      momwt_(:) = MOMZ_(:,ke) * lmesh%GsqrtH(elem%IndexH2Dto3D,ke2d) / lmesh%Gsqrt(:,ke) &
                + lmesh%GI3(:,ke,1) * MOMX_(:,ke) + lmesh%GI3(:,ke,2) * MOMY_(:,ke)

      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * MOMX_(:,ke) * QTRC_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * MOMY_(:,ke) * QTRC_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * momwt_(:) * QTRC_(:,ke)  , Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx)

      QTRC_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)     &
          + lmesh%Escale(:,ke,2,2) * Fy(:)     &
          + lmesh%Escale(:,ke,3,3) * Fz(:)     &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
    end do
  
    call PROF_rapend('cal_trcadv_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_trcadvect3d_heve_cal_tend

!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_heve_numflux_get_generalhvc( &
    del_flux,                                                    & ! (out)
    QTRC_, MOMX_, MOMY_, MOMZ_,                                  & ! (in)
    Gsqrt, GsqrtH, G13, G23, nx, ny, nz,                         & ! (in)
    vmapM, vmapP, iM2Dto3D, lmesh, elem, lmesh2D, elem2D         ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  QTRC_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GsqrtH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: iM2Dto3D(elem%NfpTot)
    
    integer :: ke, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: MomFlxP(elem%NfpTot), MomFlxM(elem%NfpTot), alpha(elem%NfpTot)
    real(RP) :: QTRC_P(elem%NfpTot), QTRC_M(elem%NfpTot)
    real(RP) :: GsqrtMOMX_P(elem%NfpTot), GsqrtMOMX_M(elem%NfpTot)
    real(RP) :: GsqrtMOMY_P(elem%NfpTot), GsqrtMOMY_M(elem%NfpTot)
    real(RP) :: GsqrtMOMZ_P(elem%NfpTot), GsqrtMOMZ_M(elem%NfpTot)
    real(RP) :: Gsqrt_P(elem%NfpTot), Gsqrt_M(elem%NfpTot)
    real(RP) :: GsqrtV_P(elem%NfpTot), GsqrtV_M(elem%NfpTot)
    real(RP) :: G13_M(elem%NfpTot), G13_P(elem%NfpTot)
    real(RP) :: G23_M(elem%NfpTot), G23_P(elem%NfpTot)   
    !------------------------------------------------------------------------

    !$omp parallel do private( &
    !$omp ke, iM, iP, ke2D,                                                                 &
    !$omp alpha, MomFlxM, MomFlxP,                                                          &
    !$omp GsqrtMOMX_M, GsqrtMOMX_P, GsqrtMOMY_M, GsqrtMOMY_P, GsqrtMOMZ_M, GsqrtMOMZ_P,     &
    !$omp QTRC_M, QTRC_P, Gsqrt_P, Gsqrt_M, GsqrtV_P, GsqrtV_M, G13_P, G13_M, G23_P, G23_M  )
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      Gsqrt_M(:) = Gsqrt(iM)
      Gsqrt_P(:) = Gsqrt(iP)
      GsqrtV_M(:) = Gsqrt_M(:) / GsqrtH(iM2Dto3D(:),ke2D)
      GsqrtV_P(:) = Gsqrt_P(:) / GsqrtH(iM2Dto3D(:),ke2D)

      G13_M(:) = G13(iM)
      G13_P(:) = G13(iP)
      G23_M(:) = G23(iM)
      G23_P(:) = G23(iP)

      QTRC_M(:) = QTRC_(iM)
      QTRC_P(:) = QTRC_(iP)
      GsqrtMOMX_M (:) = Gsqrt_M(:) * MOMX_ (iM)
      GsqrtMOMX_P (:) = Gsqrt_P(:) * MOMX_ (iP)
      GsqrtMOMY_M (:) = Gsqrt_M(:) * MOMY_ (iM)
      GsqrtMOMY_P (:) = Gsqrt_P(:) * MOMY_ (iP)
      GsqrtMOMZ_M (:) = Gsqrt_M(:) * MOMZ_ (iM)
      GsqrtMOMZ_P (:) = Gsqrt_P(:) * MOMZ_ (iP)

      MomFlxM(:) = ( GsqrtMOMX_M(:) * nx(:,ke) + GsqrtMOMY_M(:) * ny(:,ke)                  &
                 + ( ( GsqrtMOMZ_M(:) / GsqrtV_M(:)                                         &
                     + G13_M(:) * GsqrtMOMX_M(:) + G23_M(:) * GsqrtMOMY_M(:) ) * nz(:,ke) ) &
                 ) 
      MomFlxP(:) = ( GsqrtMOMX_P(:) * nx(:,ke) + GsqrtMOMY_P(:) * ny(:,ke)                  &
                 + ( ( GsqrtMOMZ_P(:) / GsqrtV_P(:)                                         &
                     + G13_P(:) * GsqrtMOMX_P(:) + G23_P(:) * GsqrtMOMY_P(:) ) * nz(:,ke) ) &
                 ) 

      alpha(:) = max( abs(MomFlxM(:)), abs(MomFlxP(:)) )
      
      del_flux(:,ke) = 0.5_RP * ( &
                    ( QTRC_P(:) * MomFlxP(:) - QTRC_M(:) * MomFlxM(:) )  &
                    - alpha(:) * ( QTRC_P(:) - QTRC_M(:) )               )
    end do

    return
  end subroutine atm_dyn_dgm_trcadvect3d_heve_numflux_get_generalhvc

end module scale_atm_dyn_dgm_trcadvect3d_heve