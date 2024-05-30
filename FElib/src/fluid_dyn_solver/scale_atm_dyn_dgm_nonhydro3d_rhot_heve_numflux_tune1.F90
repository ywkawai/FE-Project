!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Nonhydrostatic model / HEVE / Numflux
!!
!! @par Description
!!      HEVE DGM scheme for Atmospheric dynamical process. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux
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
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    DENS_VID => PRGVAR_DDENS_ID, RHOT_VID => PRGVAR_DRHOT_ID, &
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID,                               &
    PRGVAR_NUM
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

contains
 
!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc( &
    del_flux, del_flux_hyd,                                          & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, & ! (in)
    Rtot, CVtot, CPtot,                                              & ! (in)
    Gsqrt, G13, G23, nx, ny, nz,                                     & ! (in)
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D                       ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP), intent(out) ::  del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DPRES_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Rtot (elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CVtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CPtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    
    integer :: ke, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: VelP(elem%NfpTot), VelM(elem%NfpTot), alpha(elem%NfpTot)
    real(RP) :: dpresP(elem%NfpTot), dpresM(elem%NfpTot)
    real(RP) :: GsqrtDensM(elem%NfpTot), GsqrtDensP(elem%NfpTot)
    real(RP) :: GsqrtRhotM(elem%NfpTot), GsqrtRhotP(elem%NfpTot)
    real(RP) :: GsqrtDDENS_P(elem%NfpTot), GsqrtDDENS_M(elem%NfpTot)
    real(RP) :: GsqrtMOMX_P(elem%NfpTot), GsqrtMOMX_M(elem%NfpTot)
    real(RP) :: GsqrtMOMY_P(elem%NfpTot), GsqrtMOMY_M(elem%NfpTot)
    real(RP) :: GsqrtMOMZ_P(elem%NfpTot), GsqrtMOMZ_M(elem%NfpTot)
    real(RP) :: GsqrtDRHOT_P(elem%NfpTot), GsqrtDRHOT_M(elem%NfpTot)
    real(RP) :: Phyd_P(elem%NfpTot), Phyd_M(elem%NfpTot)
    real(RP) :: Gsqrt_P(elem%NfpTot), Gsqrt_M(elem%NfpTot)
    real(RP) :: GsqrtV_P(elem%NfpTot), GsqrtV_M(elem%NfpTot)
    real(RP) :: G13_M(elem%NfpTot), G13_P(elem%NfpTot)
    real(RP) :: G23_M(elem%NfpTot), G23_P(elem%NfpTot)
    real(RP) :: Gnn_M(elem%NfpTot), Gnn_P(elem%NfpTot)

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR     
    !------------------------------------------------------------------------

    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( &
    !$omp ke, iM, iP, ke2D,                                                             &
    !$omp alpha, VelM, VelP,                                                            &
    !$omp dpresM, dpresP, GsqrtDensM, GsqrtDensP, GsqrtRhotM, GsqrtRhotP,               &
    !$omp GsqrtMOMX_M, GsqrtMOMX_P, GsqrtMOMY_M, GsqrtMOMY_P, GsqrtMOMZ_M, GsqrtMOMZ_P, &
    !$omp GsqrtDDENS_M, GsqrtDDENS_P, GsqrtDRHOT_M, GsqrtDRHOT_P,                       &
    !$omp Phyd_M, Phyd_P,                                                               &
    !$omp Gsqrt_P, Gsqrt_M, GsqrtV_P, GsqrtV_M, G13_P, G13_M, G23_P, G23_M,             &
    !$omp Gnn_P, Gnn_M                                                                  )
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      Gsqrt_M(:) = Gsqrt(iM)
      Gsqrt_P(:) = Gsqrt(iP)
      GsqrtV_M(:) = Gsqrt_M(:)
      GsqrtV_P(:) = Gsqrt_P(:)

      G13_M(:) = G13(iM)
      G13_P(:) = G13(iP)
      G23_M(:) = G23(iM)
      G23_P(:) = G23(iP)

      GsqrtDDENS_M(:) = Gsqrt_M(:) * DDENS_(iM)
      GsqrtDDENS_P(:) = Gsqrt_P(:) * DDENS_(iP)
      GsqrtMOMX_M (:) = Gsqrt_M(:) * MOMX_ (iM)
      GsqrtMOMX_P (:) = Gsqrt_P(:) * MOMX_ (iP)
      GsqrtMOMY_M (:) = Gsqrt_M(:) * MOMY_ (iM)
      GsqrtMOMY_P (:) = Gsqrt_P(:) * MOMY_ (iP)
      GsqrtMOMZ_M (:) = Gsqrt_M(:) * MOMZ_ (iM)
      GsqrtMOMZ_P (:) = Gsqrt_P(:) * MOMZ_ (iP)
      GsqrtDRHOT_M(:) = Gsqrt_M(:) * DRHOT_(iM)
      GsqrtDRHOT_P(:) = Gsqrt_P(:) * DRHOT_(iP)
      Phyd_M(:) = PRES_hyd(iM)
      Phyd_P(:) = PRES_hyd(iP)

      Gnn_M(:) = abs( nx(:,ke) ) + abs( ny(:,ke) ) &
               + ( 1.0_RP / GsqrtV_M(:)**2 + G13_M(:)**2 + G23_M(:)**2 ) * abs( nz(:,ke) )
      Gnn_P(:) = abs( nx(:,ke) ) + abs( ny(:,ke) ) &
               + ( 1.0_RP / GsqrtV_P(:)**2 + G13_P(:)**2 + G23_P(:)**2 ) * abs( nz(:,ke) )

      GsqrtDensM(:) = GsqrtDDENS_M(:) + Gsqrt_M(:) * DENS_hyd(iM)
      GsqrtDensP(:) = GsqrtDDENS_P(:) + Gsqrt_P(:) * DENS_hyd(iP)

      GsqrtRhotM(:) = Gsqrt_M(:) * P0ovR * (Phyd_M(:) * rP0)**rgamm + GsqrtDRHOT_M(:)
      GsqrtRhotP(:) = Gsqrt_P(:) * P0ovR * (Phyd_P(:) * rP0)**rgamm + GsqrtDRHOT_P(:)

      VelM(:) = ( GsqrtMOMX_M(:) * nx(:,ke) + GsqrtMOMY_M(:) * ny(:,ke)                    &
                + ( ( GsqrtMOMZ_M(:) / GsqrtV_M(:)                                         &
                    + G13_M(:) * GsqrtMOMX_M(:) + G23_M(:) * GsqrtMOMY_M(:) ) * nz(:,ke) ) &
                ) / GsqrtDensM(:)
      VelP(:) = ( GsqrtMOMX_P(:) * nx(:,ke) + GsqrtMOMY_P(:) * ny(:,ke)                    &
                + ( ( GsqrtMOMZ_P(:) / GsqrtV_P(:)                                         &
                    + G13_P(:) * GsqrtMOMX_P(:) + G23_P(:) * GsqrtMOMY_P(:) ) * nz(:,ke) ) &
                ) / GsqrtDensP(:)
        

      ! dpresM(:) = PRES00 * ( Rtot(iM) * rP0 * GsqrtRhotM(:) / Gsqrt_M(:) )**( CPtot(iM) / CVtot(iM) ) &
      !           - Phyd_M(:)
      ! dpresP(:) = PRES00 * ( Rtot(iP) * rP0 * GsqrtRhotP(:) / Gsqrt_P(:) )**( CPtot(iP) / CVtot(iP) ) &
      !           - Phyd_P(:)
      dpresM(:) = DPRES_(iM)
      dpresP(:) = DPRES_(iP)

      alpha(:) = max( sqrt( Gnn_M(:) * gamm * ( Phyd_M(:) + dpresM(:) ) * Gsqrt_M(:) / GsqrtDensM(:) ) + abs(VelM(:)), &
                      sqrt( Gnn_P(:) * gamm * ( Phyd_P(:) + dpresP(:) ) * Gsqrt_P(:) / GsqrtDensP(:) ) + abs(VelP(:))  )
      
      del_flux(:,ke,DENS_VID) = 0.5_RP * ( &
                    ( GsqrtDensP(:) * VelP(:) - GsqrtDensM(:) * VelM(:) )  &
                    - alpha(:) * ( GsqrtDDENS_P(:) - GsqrtDDENS_M(:) )     )

      del_flux(:,ke,MOMX_VID ) = 0.5_RP * ( &
                    ( GsqrtMOMX_P(:) * VelP(:) - GsqrtMOMX_M(:) * VelM(:) )           &
                    + (  Gsqrt_P(:) * ( nx(:,ke) + G13_P(:) * nz(:,ke)) * dpresP(:)   &
                       - Gsqrt_M(:) * ( nx(:,ke) + G13_M(:) * nz(:,ke)) * dpresM(:) ) &
                    - alpha(:) * ( GsqrtMOMX_P(:) - GsqrtMOMX_M(:) )                  )

      del_flux(:,ke,MOMY_VID ) = 0.5_RP * ( &
                    ( GsqrtMOMY_P(:) * VelP(:) - GsqrtMOMY_M(:) * VelM(:) ) &
                    + (  Gsqrt_P(:) * ( ny(:,ke) + G23_P(:) * nz(:,ke)) * dpresP(:)   &
                       - Gsqrt_M(:) * ( ny(:,ke) + G23_M(:) * nz(:,ke)) * dpresM(:) ) &
                    - alpha(:) * ( GsqrtMOMY_P(:) - GsqrtMOMY_M(:) )        )

      del_flux(:,ke,MOMZ_VID ) = 0.5_RP * ( &
                    ( GsqrtMOMZ_P(:) * VelP(:) - GsqrtMOMZ_M(:) * VelM(:) ) &
                    + (  Gsqrt_P(:) * dpresP(:) / GsqrtV_P(:)               &
                       - Gsqrt_M(:) * dpresM(:) / GsqrtV_M(:) ) * nz(:,ke)  &
                    - alpha(:) * ( GsqrtMOMZ_P(:) - GsqrtMOMZ_M(:) )        )
                    
      del_flux(:,ke,RHOT_VID) = 0.5_RP * ( &
                    ( GsqrtRhotP(:) * VelP(:) - GsqrtRhotM(:) * VelM(:) )   &
                    - alpha(:) * ( GsqrtDRHOT_P(:) - GsqrtDRHOT_M(:) )      )

      del_flux_hyd(:,ke,1) = 0.5_RP * ( &
          GsqrtV_P(:) * ( nx(:,ke) + G13_P(:) * nz(:,ke) ) * Phyd_P(:) &
        - GsqrtV_M(:) * ( nx(:,ke) + G13_M(:) * nz(:,ke) ) * Phyd_M(:) )
      
      del_flux_hyd(:,ke,2) = 0.5_RP * ( &
          GsqrtV_P(:) * ( ny(:,ke) + G23_P(:) * nz(:,ke) ) * Phyd_P(:) &
        - GsqrtV_M(:) * ( ny(:,ke) + G23_M(:) * nz(:,ke) ) * Phyd_M(:) )
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc( &
    del_flux, del_flux_hyd,                                           & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd,  & ! (in)
    Rtot, CVtot, CPtot,                                               & ! (in)
    Gsqrt, G11, G12, G22, GsqrtH, gam, G13, G23, nx, ny, nz,          & ! (in)
    vmapM, vmapP, iM2Dto3D, lmesh, elem, lmesh2D, elem2D              ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP), intent(out) ::  del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DPRES_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Rtot (elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CVtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CPtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G11(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G12(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G22(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  GsqrtH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: gam(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: iM2Dto3D(elem%NfpTot)
    
    integer :: ke, kp, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: VelP(elem%NfpTot), VelM(elem%NfpTot), alpha(elem%NfpTot)
    real(RP) :: dpresP(elem%NfpTot), dpresM(elem%NfpTot)
    real(RP) :: GsqrtDensM(elem%NfpTot), GsqrtDensP(elem%NfpTot)
    real(RP) :: GsqrtRhotM(elem%NfpTot), GsqrtRhotP(elem%NfpTot)
    real(RP) :: GsqrtDDENS_P(elem%NfpTot), GsqrtDDENS_M(elem%NfpTot)
    real(RP) :: GsqrtMOMX_P(elem%NfpTot), GsqrtMOMX_M(elem%NfpTot)
    real(RP) :: GsqrtMOMY_P(elem%NfpTot), GsqrtMOMY_M(elem%NfpTot)
    real(RP) :: GsqrtMOMZ_P(elem%NfpTot), GsqrtMOMZ_M(elem%NfpTot)
    real(RP) :: GsqrtDRHOT_P(elem%NfpTot), GsqrtDRHOT_M(elem%NfpTot)
    real(RP) :: Phyd_P(elem%NfpTot), Phyd_M(elem%NfpTot)
    real(RP) :: Gsqrt_P(elem%NfpTot), Gsqrt_M(elem%NfpTot)
    real(RP) :: GsqrtV_P(elem%NfpTot), GsqrtV_M(elem%NfpTot)
    real(RP) :: G13_M(elem%NfpTot), G13_P(elem%NfpTot)
    real(RP) :: G23_M(elem%NfpTot), G23_P(elem%NfpTot)
    real(RP) :: G1n_M(elem%NfpTot), G2n_M(elem%NfpTot)
    real(RP) :: Gnn_M(elem%NfpTot), Gnn_P(elem%NfpTot)
    real(RP) :: Gxz_M(elem%NfpTot), Gxz_P(elem%NfpTot)
    real(RP) :: Gyz_M(elem%NfpTot), Gyz_P(elem%NfpTot)
    real(RP) :: rgam2_M(elem%NfpTot), rgam2_P(elem%NfpTot)

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR

    !------------------------------------------------------------------------

    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !if(PRC_GLOBAL_myrank == 0) then
    !  write(*,*) "lmesh.NeA", lmesh%NeA
    !  write(*,*) "lmesh.NeS", lmesh%NeS
    !  write(*,*) "lmesh.NeE", lmesh%NeE
    !  write(*,*) "lmesh2D.Ne", lmesh2D%Ne
    !  write(*,*) "elem.NfpTot", elem%NfpTot
    !end if

    !$omp parallel do private( &
    !$omp ke, kp, iM, iP, ke2D,                                                             &
    !$omp alpha, VelM, VelP,                                                            &
    !$omp dpresM, dpresP, GsqrtDensM, GsqrtDensP, GsqrtRhotM, GsqrtRhotP,               &
    !$omp GsqrtMOMX_M, GsqrtMOMX_P, GsqrtMOMY_M, GsqrtMOMY_P, GsqrtMOMZ_M, GsqrtMOMZ_P, &
    !$omp GsqrtDDENS_M, GsqrtDDENS_P, GsqrtDRHOT_M, GsqrtDRHOT_P,                       &
    !$omp Phyd_M, Phyd_P,                                                               &
    !$omp Gsqrt_P, Gsqrt_M, GsqrtV_P, GsqrtV_M, G13_P, G13_M, G23_P, G23_M,             &
    !$omp rgam2_M, rgam2_P, Gxz_P, Gxz_M, Gyz_P, Gyz_M, G1n_M, G2n_M, Gnn_P, Gnn_M      )
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      do kp=1, elem%NfpTot
        G13_M(kp) = G13(iM(kp))
        G13_P(kp) = G13(iP(kp))
      end do
      do kp=1, elem%NfpTot
        G23_M(kp) = G23(iM(kp))
        G23_P(kp) = G23(iP(kp))
      end do

      do kp=1, elem%NfpTot
        Gsqrt_M(kp) = Gsqrt(iM(kp))
        Gsqrt_P(kp) = Gsqrt(iP(kp))
      end do

      do kp=1, elem%NfpTot
        rgam2_M(kp) = 1.0_RP / gam(iM(kp))**2
        rgam2_P(kp) = 1.0_RP / gam(iP(kp))**2
      end do

      do kp=1, elem%NfpTot
        GsqrtV_M(kp) = Gsqrt_M(kp) * rgam2_M(kp) / GsqrtH(iM2Dto3D(kp),ke2D)
        GsqrtV_P(kp) = Gsqrt_P(kp) * rgam2_P(kp) / GsqrtH(iM2Dto3D(kp),ke2D)
      end do

      do kp=1, elem%NfpTot
        Gxz_M(kp) = rgam2_M(kp) * ( G11(iM2Dto3D(kp),ke2D) * G13_M(kp) + G12(iM2Dto3D(kp),ke2D) * G23_M(kp) )
        Gxz_P(kp) = rgam2_P(kp) * ( G11(iM2Dto3D(kp),ke2D) * G13_P(kp) + G12(iM2Dto3D(kp),ke2D) * G23_P(kp) )
      end do

      do kp=1, elem%NfpTot
        Gyz_M(kp) = rgam2_M(kp) * ( G12(iM2Dto3D(kp),ke2D) * G13_M(kp) + G22(iM2Dto3D(kp),ke2D) * G23_M(kp) )
        Gyz_P(kp) = rgam2_P(kp) * ( G12(iM2Dto3D(kp),ke2D) * G13_P(kp) + G22(iM2Dto3D(kp),ke2D) * G23_P(kp) )
      end do

      do kp=1, elem%NfpTot
        G1n_M(kp) = rgam2_M(kp) * ( G11(iM2Dto3D(kp),ke2D) * nx(kp,ke) + G12(iM2Dto3D(kp),ke2D) * ny(kp,ke) )
        G2n_M(kp) = rgam2_P(kp) * ( G12(iM2Dto3D(kp),ke2D) * nx(kp,ke) + G22(iM2Dto3D(kp),ke2D) * ny(kp,ke) )
      end do

      do kp=1, elem%NfpTot
        Gnn_M(kp)  = rgam2_M(kp) * ( G11(iM2Dto3D(kp),ke2D) * abs( nx(kp,ke) ) + G22(iM2Dto3D(kp),ke2D) * abs( ny(kp,ke) ) ) &
                  + ( 1.0_RP / GsqrtV_M(kp)**2 + G13_M(kp) * Gxz_M(kp) + G23_M(kp) * Gyz_M(kp) ) * abs( nz(kp,ke) )
        Gnn_P(kp)  = rgam2_P(kp) * ( G11(iM2Dto3D(kp),ke2D) * abs( nx(kp,ke) ) + G22(iM2Dto3D(kp),ke2D) * abs( ny(kp,ke) ) ) &
                  + ( 1.0_RP / GsqrtV_P(kp)**2 + G13_P(kp) * Gxz_P(kp) + G23_P(kp) * Gyz_P(kp) ) * abs( nz(kp,ke) )
      end do

      do kp=1, elem%NfpTot
        GsqrtDDENS_M(kp) = Gsqrt_M(kp) * DDENS_(iM(kp))
        GsqrtMOMX_M (kp) = Gsqrt_M(kp) * MOMX_ (iM(kp))
        GsqrtMOMY_M (kp) = Gsqrt_M(kp) * MOMY_ (iM(kp))
        GsqrtMOMZ_M (kp) = Gsqrt_M(kp) * MOMZ_ (iM(kp))
        GsqrtDRHOT_M(kp) = Gsqrt_M(kp) * DRHOT_(iM(kp))
      end do

      do kp=1, elem%NfpTot
        GsqrtDDENS_P(kp) = Gsqrt_P(kp) * DDENS_(iP(kp))
        GsqrtMOMX_P (kp) = Gsqrt_P(kp) * MOMX_ (iP(kp))
        GsqrtMOMY_P (kp) = Gsqrt_P(kp) * MOMY_ (iP(kp))
        GsqrtMOMZ_P (kp) = Gsqrt_P(kp) * MOMZ_ (iP(kp))
        GsqrtDRHOT_P(kp) = Gsqrt_P(kp) * DRHOT_(iP(kp))
      end do

      do kp=1, elem%NfpTot
        Phyd_M(kp) = PRES_hyd(iM(kp))
        Phyd_P(kp) = PRES_hyd(iP(kp))
      end do

      do kp=1, elem%NfpTot
        GsqrtDensM(kp) = GsqrtDDENS_M(kp) + Gsqrt_M(kp) * DENS_hyd(iM(kp))
        GsqrtDensP(kp) = GsqrtDDENS_P(kp) + Gsqrt_P(kp) * DENS_hyd(iP(kp))
      end do


      do kp=1, elem%NfpTot
        GsqrtRhotM(kp) = Gsqrt_M(kp) * P0ovR * (Phyd_M(kp) * rP0)**rgamm + GsqrtDRHOT_M(kp)
      end do

      do kp=1, elem%NfpTot
        GsqrtRhotP(kp) = Gsqrt_P(kp) * P0ovR * (Phyd_P(kp) * rP0)**rgamm + GsqrtDRHOT_P(kp)
      end do

      do kp=1, elem%NfpTot
        VelM(kp) = ( GsqrtMOMX_M(kp) * nx(kp,ke) + GsqrtMOMY_M(kp) * ny(kp,ke)                    &
                  + ( ( GsqrtMOMZ_M(kp) / GsqrtV_M(kp)                                         &
                      + G13_M(kp) * GsqrtMOMX_M(kp) + G23_M(kp) * GsqrtMOMY_M(kp) ) * nz(kp,ke) ) &
                  ) / GsqrtDensM(kp)
      end do

      do kp=1, elem%NfpTot
        VelP(kp) = ( GsqrtMOMX_P(kp) * nx(kp,ke) + GsqrtMOMY_P(kp) * ny(kp,ke)                    &
                  + ( ( GsqrtMOMZ_P(kp) / GsqrtV_P(kp)                                         &
                      + G13_P(kp) * GsqrtMOMX_P(kp) + G23_P(kp) * GsqrtMOMY_P(kp) ) * nz(kp,ke) ) &
                  ) / GsqrtDensP(kp)
      end do

      ! dpresM(:) = PRES00 * ( Rtot(iM) * rP0 * GsqrtRhotM(:) / Gsqrt_M(:) )**( CPtot(iM) / CVtot(iM) ) &
      !           - Phyd_M(:)
      ! dpresP(:) = PRES00 * ( Rtot(iP) * rP0 * GsqrtRhotP(:) / Gsqrt_P(:) )**( CPtot(iP) / CVtot(iP) ) &
      !           - Phyd_P(:)
      do kp=1, elem%NfpTot
        dpresM(kp) = DPRES_(iM(kp))
        dpresP(kp) = DPRES_(iP(kp))
      end do

      do kp=1, elem%NfpTot
        alpha(kp) = max( sqrt( Gnn_M(kp) * gamm * ( Phyd_M(kp) + dpresM(kp) ) * Gsqrt_M(kp) / GsqrtDensM(kp) ) + abs(VelM(kp)), &
                        sqrt( Gnn_P(kp) * gamm * ( Phyd_P(kp) + dpresP(kp) ) * Gsqrt_P(kp) / GsqrtDensP(kp) ) + abs(VelP(kp))  )
      end do

      do kp=1, elem%NfpTot
        del_flux(kp,ke,DENS_VID) = 0.5_RP * ( &
                      ( GsqrtDensP(kp) * VelP(kp) - GsqrtDensM(kp) * VelM(kp) )  &
                      - alpha(kp) * ( GsqrtDDENS_P(kp) - GsqrtDDENS_M(kp) )     )
      end do

      do kp=1, elem%NfpTot
        del_flux(kp,ke,MOMX_VID ) = 0.5_RP * ( &
                      ( GsqrtMOMX_P(kp) * VelP(kp) - GsqrtMOMX_M(kp) * VelM(kp) )           &
                      + (  Gsqrt_P(kp) * ( G1n_M(kp) + Gxz_P(kp) * nz(kp,ke)) * dpresP(kp)   &
                        - Gsqrt_M(kp) * ( G1n_M(kp) + Gxz_M(kp) * nz(kp,ke)) * dpresM(kp) ) &
                      - alpha(kp) * ( GsqrtMOMX_P(kp) - GsqrtMOMX_M(kp) )                  )
      end do

      do kp=1, elem%NfpTot
        del_flux(kp,ke,MOMY_VID ) = 0.5_RP * ( &
                      ( GsqrtMOMY_P(kp) * VelP(kp) - GsqrtMOMY_M(kp) * VelM(kp) )            &
                      + (  Gsqrt_P(kp) * ( G2n_M(kp) + Gyz_P(kp) * nz(kp,ke) ) * dpresP(kp)   &
                        - Gsqrt_M(kp) * ( G2n_M(kp) + Gyz_M(kp) * nz(kp,ke) ) * dpresM(kp) ) &
                      - alpha(kp) * ( GsqrtMOMY_P(kp) - GsqrtMOMY_M(kp) )        )
      end do
      
      do kp=1, elem%NfpTot
          del_flux(kp,ke,MOMZ_VID ) = 0.5_RP * ( &
                        ( GsqrtMOMZ_P(kp) * VelP(kp) - GsqrtMOMZ_M(kp) * VelM(kp) ) &
                        + (  Gsqrt_P(kp) * dpresP(kp) / GsqrtV_P(kp)               &
                          - Gsqrt_M(kp) * dpresM(kp) / GsqrtV_M(kp) ) * nz(kp,ke)  &
                        - alpha(kp) * ( GsqrtMOMZ_P(kp) - GsqrtMOMZ_M(kp) )        )
      end do
                    
      do kp=1, elem%NfpTot
        del_flux(kp,ke,RHOT_VID) = 0.5_RP * ( &
                      ( GsqrtRhotP(kp) * VelP(kp) - GsqrtRhotM(kp) * VelM(kp) )   &
                      - alpha(kp) * ( GsqrtDRHOT_P(kp) - GsqrtDRHOT_M(kp) )      )
      end do

      do kp=1, elem%NfpTot
        del_flux_hyd(kp,ke,1) = 0.5_RP * ( &
            GsqrtV_P(kp) * ( nx(kp,ke) + G13_P(kp) * nz(kp,ke) ) * Phyd_P(kp) &
          - GsqrtV_M(kp) * ( nx(kp,ke) + G13_M(kp) * nz(kp,ke) ) * Phyd_M(kp) )
      end do
      
      do kp=1, elem%NfpTot
        del_flux_hyd(kp,ke,2) = 0.5_RP * ( &
            GsqrtV_P(kp) * ( ny(kp,ke) + G23_P(kp) * nz(kp,ke) ) * Phyd_P(kp) &
          - GsqrtV_M(kp) * ( ny(kp,ke) + G23_M(kp) * nz(kp,ke) ) * Phyd_M(kp) )
      end do
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc

end module scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux
