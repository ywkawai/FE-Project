!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Nonhydrostatic model / HEVE / Numflux
!!
!! @par Description
!!      HEVE DGM scheme for Atmospheric dynamical process. 
!!
!! @author Yuta Kawai, Team SCALE
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
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_asis
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc
!  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_LMARS
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc_asis
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc
  public ::  atm_dyn_dgm_nonhydro3d_rhot_heve_add_bnd_contrib_generalvc
  ! public ::  atm_dyn_dgm_nonhydro3d_rhot_heve_add_bnd_contrib2_generalvc

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
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_asis( &
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
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_asis

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_add_bnd_contrib_generalvc( &
    PRGVAR_dt,                                                       & ! (inout)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES,                      & ! (in)
    DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot,               & ! (in)
    Gsqrt, G13, G23, nx, ny, nz,                                     & ! (in)
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D, elem3D_optr          ) ! (in)
    use scale_element_operation_base, only: ElementOperationBase3D
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(inout) :: PRGVAR_dt(elem%Np,lmesh%NeA,5)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)
    real(RP), intent(inout) ::  DPRES(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  THERM_hyd(elem%Np*lmesh%NeA)
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
    class(ElementOperationBase3D), intent(in) :: elem3D_optr

    integer :: ke, fp, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: Vel(elem%NfpTot,2), alpha(elem%NfpTot)
    real(RP) :: DPRES_(elem%NfpTot,2)
    real(RP) :: GsqrtDens(elem%NfpTot,2)
    real(RP) :: GsqrtRhot(elem%NfpTot,2)
    real(RP) :: GsqrtDDENS(elem%NfpTot,2)
    real(RP) :: GsqrtMOMX(elem%NfpTot,2)
    real(RP) :: GsqrtMOMY(elem%NfpTot,2)
    real(RP) :: GsqrtMOMZ(elem%NfpTot,2)
    real(RP) :: GsqrtDRHOT(elem%NfpTot,2)
    real(RP) :: Cs2_(elem%NfpTot,2)
    real(RP) :: Gsqrt_(elem%NfpTot,2)
    real(RP) :: GsqrtV_(elem%NfpTot,2)
    real(RP) :: RGsqrtV(elem%NfpTot,2)
    real(RP) :: G13_(elem%NfpTot,2)
    real(RP) :: G23_(elem%NfpTot,2)
    real(RP) :: Gnn_M, Gnn_P

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR
    
    integer, parameter :: IN = 1
    integer, parameter :: EX = 2

    real(RP) :: tmp1, tmp2, tmp3, tmp4, tmp01, tmp02, tmp03
    real(RP) :: del_flux_tmp_mom(elem%NfpTot,3)

    integer :: iv
    integer :: p
    real(RP) :: del_flux(elem%NfpTot,5)
    real(RP) :: Lift(elem%Np,5)
    real(RP) :: RGsqrt(elem%Np)

    integer :: iim, iip, fpp
    integer, parameter :: Nfp_s = 0
    !------------------------------------------------------------------------

    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry
    call PROF_rapstart('cal_dyn_tend_bndflux1', 3)

    !$omp parallel private( &
    !$omp ke, iM, iP, ke2D, fp, i, iim, iip, fpp,                                                    &
    !$omp alpha, Vel,                                                                   &
    !$omp dpres_, GsqrtDens, GsqrtRhot,                                                 &
    !$omp GsqrtDDENS, GsqrtMOMX, GsqrtMOMY, GsqrtMOMZ, GsqrtDRHOT,                      &
    !$omp Cs2_,                                                                         &
    !$omp Gsqrt_, GsqrtV_, RGsqrtV, G13_, G23_,                                         &
    !$omp Gnn_P, Gnn_M, tmp1, tmp2, tmp3, tmp4, tmp01, tmp02, tmp03, del_flux_tmp_mom,  &
    !$omp del_flux, p, iv, Lift, RGsqrt )

    !$omp do
    do i=elem%Np*lmesh%NeE+1, elem%Np*lmesh%NeE+size(lmesh%VMapB)
      DPRES(i) = PRES00 * ( Rtot(i) * rP0 * ( THERM_hyd(i) + DRHOT_(i) ) )**( CPtot(i) / CVtot(i) ) &
               - PRES_hyd(i)
    end do
    !$omp end do
!OCL PREFETCH
    !$omp do
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      do fp=1, elem%NfpTot
        iim = iM(fp)
        Gsqrt_(fp,IN) = Gsqrt(iim)
        GsqrtV_(fp,IN) = Gsqrt_(fp,IN)
        RGsqrtV(fp,IN) = 1.0_RP / GsqrtV_(fp,IN)  

        GsqrtDDENS(fp,IN) = Gsqrt_(fp,IN) * DDENS_(iim)
        GsqrtDens(fp,IN) = GsqrtDDENS(fp,IN) + Gsqrt_(fp,IN) * DENS_hyd(iim)

        GsqrtDRHOT(fp,IN) = Gsqrt_(fp,IN) * DRHOT_(iim)
        GsqrtRhot(fp,IN) = GsqrtDRHOT(fp,IN) + Gsqrt_(fp,IN) * THERM_hyd(iim) 
      end do
      do fp=1, elem%NfpTot
        iim = iM(fp)

        DPRES_(fp,IN) = DPRES(iim)
        Cs2_(fp,IN) = gamm * ( PRES_hyd(iim) + DPRES_(fp,IN) ) * Gsqrt_(fp,IN) / GsqrtDens(fp,IN)
      end do
      do fp=1, elem%NfpTot
        iim = iM(fp)
        GsqrtMOMX (fp,IN) = Gsqrt_(fp,IN) * MOMX_ (iim)
        GsqrtMOMY (fp,IN) = Gsqrt_(fp,IN) * MOMY_ (iim)
        GsqrtMOMZ (fp,IN) = Gsqrt_(fp,IN) * MOMZ_ (iim)
        
        G13_(fp,IN) = G13(iim)
        G23_(fp,IN) = G23(iim)
      end do
      do fp=1, elem%NfpTot
        fpp = Nfp_s + fp
        Vel(fp,IN) = ( GsqrtMOMX(fp,IN) * nx(fpp,ke) + GsqrtMOMY(fp,IN) * ny(fpp,ke)                     &
                  + ( ( GsqrtMOMZ(fp,IN) * RGsqrtV(fp,IN)                                                &
                      + G13_(fp,IN) * GsqrtMOMX(fp,IN) + G23_(fp,IN) * GsqrtMOMY(fp,IN) ) * nz(fpp,ke) ) &
                  ) / GsqrtDens(fp,IN)
      end do

      !-
      do fp=1, elem%NfpTot
        iip = iP(fp)

        Gsqrt_(fp,EX) = Gsqrt(iip)
        GsqrtV_(fp,EX) = Gsqrt_(fp,EX)
        RGsqrtV(fp,EX) = 1.0_RP / GsqrtV_(fp,EX)  


        GsqrtDDENS(fp,EX) = Gsqrt_(fp,EX) * DDENS_(iip)
        GsqrtDens(fp,EX) = GsqrtDDENS(fp,EX) + Gsqrt_(fp,EX) * DENS_hyd(iip)

        GsqrtDRHOT(fp,EX) = Gsqrt_(fp,EX) * DRHOT_(iip)
        GsqrtRhot(fp,EX) = GsqrtDRHOT(fp,EX) + Gsqrt_(fp,EX) * THERM_hyd(iip)
      end do
      do fp=1, elem%NfpTot
        iip = iP(fp)
        DPRES_(fp,EX) = DPRES(iip)
        Cs2_(fp,EX) = gamm * ( PRES_hyd(iip) + DPRES_(fp,EX) ) * Gsqrt_(fp,EX) / GsqrtDens(fp,EX)
      end do
      do fp=1, elem%NfpTot
        iip = iP(fp)
        GsqrtMOMX (fp,EX) = Gsqrt_(fp,EX) * MOMX_ (iip)
        GsqrtMOMY (fp,EX) = Gsqrt_(fp,EX) * MOMY_ (iip)
        GsqrtMOMZ (fp,EX) = Gsqrt_(fp,EX) * MOMZ_ (iip)

        G13_(fp,EX) = G13(iip)
        G23_(fp,EX) = G23(iip)
      end do
      do fp=1, elem%NfpTot
        fpp = Nfp_s + fp
        Vel(fp,EX) = ( GsqrtMOMX(fp,EX) * nx(fpp,ke) + GsqrtMOMY(fp,EX) * ny(fpp,ke)                     &
                  + ( ( GsqrtMOMZ(fp,EX) * RGsqrtV(fp,EX)                                                &
                      + G13_(fp,EX) * GsqrtMOMX(fp,EX) + G23_(fp,EX) * GsqrtMOMY(fp,EX) ) * nz(fpp,ke) ) &
                  ) / GsqrtDens(fp,EX)
      end do

      do fp=1, elem%NfpTot
        tmp1 = abs( nx(fp,ke) ) + abs( ny(fp,ke) )
        Gnn_M = tmp1 &
                + ( 1.0_RP * RGsqrtV(fp,IN)**2 + G13_(fp,IN)**2 + G23_(fp,IN)**2 ) * abs( nz(fp,ke) )
        Gnn_P = tmp1 &
                + ( 1.0_RP * RGsqrtV(fp,EX)**2 + G13_(fp,EX)**2 + G23_(fp,EX)**2 ) * abs( nz(fp,ke) )

        alpha(fp) = max( sqrt( Gnn_M * Cs2_(fp,IN) ) + abs(Vel(fp,IN)), &
                         sqrt( Gnn_P * Cs2_(fp,EX) ) + abs(Vel(fp,EX))  )
      end do
      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP

        tmp2 = - alpha(fp) * ( GsqrtDDENS(fp,EX) - GsqrtDDENS(fp,IN) )
        del_flux(fp,DENS_VID) = tmp1 * ( &
                        GsqrtDens(fp,EX) * Vel(fp,EX) - GsqrtDens(fp,IN) * Vel(fp,IN)  &
                      + tmp2  )

        tmp2 = - alpha(fp) * ( GsqrtDRHOT(fp,EX) - GsqrtDRHOT(fp,IN) )
        del_flux(fp,RHOT_VID) = tmp1 * ( &
                        GsqrtRhot(fp,EX) * Vel(fp,EX) - GsqrtRhot(fp,IN) * Vel(fp,IN) &
                      + tmp2 )
      end do

      do fp=1, elem%NfpTot
        tmp3 = Gsqrt_(fp,EX) * dpres_(fp,EX)
        tmp4 = Gsqrt_(fp,IN) * dpres_(fp,IN)

        del_flux_tmp_mom(fp,1) = &
                (  tmp3 * RGsqrtV(fp,EX)               &
                 - tmp4 * RGsqrtV(fp,IN) ) * nz(fp,ke)
        
        del_flux_tmp_mom(fp,2) = &
              ( nx(fp,ke) + G13_(fp,EX) * nz(fp,ke) ) * tmp3   &
            - ( nx(fp,ke) + G13_(fp,IN) * nz(fp,ke) ) * tmp4
        
        del_flux_tmp_mom(fp,3) = &
              ( ny(fp,ke) + G23_(fp,EX) * nz(fp,ke) ) * tmp3   &
            - ( ny(fp,ke) + G23_(fp,IN) * nz(fp,ke) ) * tmp4
      end do

      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP

        tmp2 = - alpha(fp) * ( GsqrtMOMZ(fp,EX) - GsqrtMOMZ(fp,IN) )
        del_flux(fp,MOMZ_VID) = tmp1 * ( &
                        GsqrtMOMZ(fp,EX) * Vel(fp,EX) - GsqrtMOMZ(fp,IN) * Vel(fp,IN) &
                      + del_flux_tmp_mom(fp,1) &
                      + tmp2 )
      end do

      !-
      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP

        tmp2 = - alpha(fp) * ( GsqrtMOMX(fp,EX) - GsqrtMOMX(fp,IN) )
        del_flux(fp,MOMX_VID) = tmp1 * ( &
                        GsqrtMOMX(fp,EX) * Vel(fp,EX) - GsqrtMOMX(fp,IN) * Vel(fp,IN) &
                      + del_flux_tmp_mom(fp,2) &
                      + tmp2 )

        tmp2 = - alpha(fp) * ( GsqrtMOMY(fp,EX) - GsqrtMOMY(fp,IN) )
        del_flux(fp,MOMY_VID) = tmp1 * ( &
                        GsqrtMOMY(fp,EX) * Vel(fp,EX) - GsqrtMOMY(fp,IN) * Vel(fp,IN) &
                      + del_flux_tmp_mom(fp,3) &
                      + tmp2 )
      end do

      call elem3D_optr%Lift_var5( del_flux, Lift )
      RGsqrt(:) = 1.0_RP / lmesh%Gsqrt(:,ke)
      do iv=1, 5
        PRGVAR_dt(:,ke,iv) = PRGVAR_dt(:,ke,iv) - RGsqrt(:) * Lift(:,iv)
      end do
    end do
    !$omp end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_bndflux1', 3)
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_add_bnd_contrib_generalvc


! !OCL SERIAL
!   subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_add_bnd_contrib2_generalvc( &
!     PRGVAR_dt,                                                       & ! (inout)
!     DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES,                      & ! (in)
!     DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot,               & ! (in)
!     Gsqrt, G13, G23, nx, ny, nz,                                     & ! (in)
!     vmapM, vmapP, lmesh, elem, lmesh2D, elem2D, elem3D_optr          ) ! (in)
!     use scale_element_operation_base, only: ElementOperationBase3D
!     implicit none

!     class(LocalMesh3D), intent(in) :: lmesh
!     class(ElementBase3D), intent(in) :: elem  
!     class(LocalMesh2D), intent(in) :: lmesh2D
!     class(ElementBase2D), intent(in) :: elem2D
!     real(RP), intent(inout) :: PRGVAR_dt(elem%Np,lmesh%NeA,5)
!     real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
!     real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
!     real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
!     real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)
!     real(RP), intent(inout) ::  DPRES(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  THERM_hyd(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  Rtot (elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  CVtot(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  CPtot(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
!     real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
!     real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
!     real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
!     integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
!     integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
!     class(ElementOperationBase3D), intent(in) :: elem3D_optr

!     integer :: ke, fp, i, iP(elem%NfpTot), iM(elem%NfpTot)

!     integer :: iv
!     real(RP) :: del_flux(elem%NfpTot,5)
!     real(RP) :: del_flux_P(elem%NfpTot,5)
!     real(RP) :: del_flux_save_h(elem%Nfp_h,5,elem%Nfaces_h,lmesh%Ne)
!     real(RP) :: del_flux_save_v(elem%Nfp_v,5,elem%Nfaces_v,lmesh%Ne)
!     real(RP) :: Lift(elem%Np,5)
!     real(RP) :: RGsqrt(elem%Np)
!     real(RP) :: tmp1

!     integer :: color(lmesh%Ne)
!     integer :: j, k

!     integer :: ke2, f, f2, fpp
!     integer :: Nfp_s

!     real(RP) :: rP0
!     !------------------------------------------------------------------------

!     do k=1, lmesh%NeZ
!     do j=1, lmesh%NeY
!     do i=1, lmesh%NeX
!       ke = i + (j-1)*lmesh%NeX + (k-1)*lmesh%Ne2D
!       if (mod(i+j+k,2) == 0) then
!         color(ke) = 0
!       else
!         color(ke) = 1
!       end if
!     end do
!     end do
!     end do
!     call PROF_rapstart('cal_dyn_tend_bndflux1', 3)

!     rP0 = 1.0_RP / PRES00

!     !$omp parallel private( &
!     !$omp ke, iM, iP, fp, i, f,fpp,f2,ke2, Nfp_s, &
!     !$omp del_flux, del_flux_P, iv, Lift, RGsqrt, tmp1    )

!     !$omp do
!     do i=elem%Np*lmesh%NeE+1, elem%Np*lmesh%NeE+size(lmesh%VMapB)
!       DPRES(i) = PRES00 * ( Rtot(i) * rP0 * ( THERM_hyd(i) + DRHOT_(i) ) )**( CPtot(i) / CVtot(i) ) &
!                - PRES_hyd(i)
!     end do
!     !$omp end do

!     !$omp do
!     do ke=lmesh%NeS, lmesh%NeE
!       if (color(ke) == 1) cycle

!       iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
!       Nfp_s = 0
!       call calc_bnd_flux_core( del_flux, del_flux_P,         & ! (out)
!         DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES,          & ! (in)
!         DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot,   & ! (in)
!         Gsqrt, G13, G23,                                     & ! (in)
!         nx(:,ke), ny(:,ke), nz(:,ke), lmesh%Fscale(:,ke),    & ! (in)
!         iM, iP, elem%NfpTot, Nfp_s, elem%NfpTot, lmesh, elem ) ! (in)

!       !--
!       do f=1, elem%Nfaces_h
!         ke2 = lmesh%EToE(ke,f); f2 = lmesh%EToF(ke,f)
!         Nfp_s = (f-1)*elem%Nfp_h
!         do fpp=1, elem%Nfp_h                    
!           fp = Nfp_s + fpp
!           del_flux_save_h(fpp,DENS_VID,f2,ke2) = del_flux_P(fp,DENS_VID)
!           del_flux_save_h(fpp,RHOT_VID,f2,ke2) = del_flux_P(fp,RHOT_VID)
!           del_flux_save_h(fpp,MOMZ_VID,f2,ke2) = del_flux_P(fp,MOMZ_VID)
!           del_flux_save_h(fpp,MOMX_VID,f2,ke2) = del_flux_P(fp,MOMX_VID)
!           del_flux_save_h(fpp,MOMY_VID,f2,ke2) = del_flux_P(fp,MOMY_VID)
!         end do

!         ! call calc_bnd_flux_core( del_flux, del_flux_save_h(:,:,f2,ke2), & ! (out)
!         !   DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES,                   & ! (in)
!         !   DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot,            & ! (in)
!         !   Gsqrt, G13, G23,                                              & ! (in)
!         !   nx(:,ke), ny(:,ke), nz(:,ke), lmesh%Fscale(:,ke),             & ! (in)
!         !   iM(Nfp_s+1:Nfp_s+elem%Nfp_h), iP(Nfp_s+1:Nfp_s+elem%Nfp_h),   & ! (in)
!         !   elem%Nfp_h, Nfp_s, elem%NfpTot, lmesh, elem                   ) ! (in)
!       end do
!       do f=1, elem%Nfaces_v
!         ke2 = lmesh%EToE(ke,f+elem%Nfaces_h); f2 = lmesh%EToF(ke,f+elem%Nfaces_h) - elem%Nfaces_h
!         Nfp_s = elem%Nfaces_h*elem%Nfp_h + (f-1)*elem%Nfp_v
!         do fpp=1, elem%Nfp_v                    
!           fp = Nfp_s + fpp
!           del_flux_save_v(fpp,DENS_VID,f2,ke2) = del_flux_P(fp,DENS_VID)
!           del_flux_save_v(fpp,RHOT_VID,f2,ke2) = del_flux_P(fp,RHOT_VID)
!           del_flux_save_v(fpp,MOMZ_VID,f2,ke2) = del_flux_P(fp,MOMZ_VID)
!           del_flux_save_v(fpp,MOMX_VID,f2,ke2) = del_flux_P(fp,MOMX_VID)
!           del_flux_save_v(fpp,MOMY_VID,f2,ke2) = del_flux_P(fp,MOMY_VID)
!         end do

!         ! call calc_bnd_flux_core( del_flux, del_flux_save_v(:,:,f2,ke2), & ! (out)
!         !   DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES,                   & ! (in)
!         !   DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot,            & ! (in)
!         !   Gsqrt, G13, G23,                                              & ! (in)
!         !   nx(:,ke), ny(:,ke), nz(:,ke), lmesh%Fscale(:,ke),             & ! (in)
!         !   iM(Nfp_s+1:Nfp_s+elem%Nfp_v), iP(Nfp_s+1:Nfp_s+elem%Nfp_v),   & ! (in)
!         !   elem%Nfp_v, Nfp_s, elem%NfpTot, lmesh, elem                   ) ! (in)
!       end do
      
!       !--
!       call elem3D_optr%Lift_var5( del_flux, Lift )
!       RGsqrt(:) = 1.0_RP / lmesh%Gsqrt(:,ke)
!       do iv=1, 5
!         PRGVAR_dt(:,ke,iv) = PRGVAR_dt(:,ke,iv) - RGsqrt(:) * Lift(:,iv)
!       end do
!     end do
!     !$omp end do

!     !$omp do
!     do ke=lmesh%NeS, lmesh%NeE
!       if (color(ke) == 0) cycle

!       iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
!       !--
!       do f=1, elem%Nfaces_h
!         ke2 = lmesh%EToE(ke,f); f2 = lmesh%EToF(ke,f)
!        if (ke2 == ke) then
!           Nfp_s = (f-1)*elem%Nfp_h
!           call calc_bnd_flux_core( del_flux, del_flux_save_h(:,:,f2,ke2), & ! (out)
!             DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES,                   & ! (in)
!             DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot,            & ! (in)
!             Gsqrt, G13, G23,                                              & ! (in)
!             nx(:,ke), ny(:,ke), nz(:,ke), lmesh%Fscale(:,ke),             & ! (in)
!             iM(Nfp_s+1:Nfp_s+elem%Nfp_h), iP(Nfp_s+1:Nfp_s+elem%Nfp_h),   & ! (in)
!             elem%Nfp_h, Nfp_s, elem%NfpTot, lmesh, elem                   ) ! (in)
!         else
!           do fpp=1, elem%Nfp_h                    
!             fp = fpp + (f-1)*elem%Nfp_h
!             tmp1 = 0.5_RP * lmesh%Fscale(fp,ke)
!             del_flux(fp,DENS_VID) = tmp1 * del_flux_save_h(fpp,DENS_VID,f,ke)
!             del_flux(fp,RHOT_VID) = tmp1 * del_flux_save_h(fpp,RHOT_VID,f,ke)
!             del_flux(fp,MOMZ_VID) = tmp1 * del_flux_save_h(fpp,MOMZ_VID,f,ke)
!             del_flux(fp,MOMX_VID) = tmp1 * del_flux_save_h(fpp,MOMX_VID,f,ke)
!             del_flux(fp,MOMY_VID) = tmp1 * del_flux_save_h(fpp,MOMY_VID,f,ke)          
!           end do
!         end if
!       end do
!       do f=1, elem%Nfaces_v
!         ke2 = lmesh%EToE(ke,f+elem%Nfaces_h); f2 = lmesh%EToF(ke,f+elem%Nfaces_h) - elem%Nfaces_h
!         if (ke2 == ke) then
!           Nfp_s = elem%Nfaces_h*elem%Nfp_h + (f-1)*elem%Nfp_v
!           call calc_bnd_flux_core( del_flux, del_flux_save_v(:,:,f2,ke2), & ! (out)
!             DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES,                   & ! (in)
!             DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot,            & ! (in)
!             Gsqrt, G13, G23,                                              & ! (in)
!             nx(:,ke), ny(:,ke), nz(:,ke), lmesh%Fscale(:,ke),             & ! (in)
!             iM(Nfp_s+1:Nfp_s+elem%Nfp_v), iP(Nfp_s+1:Nfp_s+elem%Nfp_v),   & ! (in)
!             elem%Nfp_v, Nfp_s, elem%NfpTot, lmesh, elem                   ) ! (in)
!         else
!           do fpp=1, elem%Nfp_v
!             fp = fpp + (f-1)*elem%Nfp_v + elem%Nfaces_h*elem%Nfp_h
!             tmp1 = 0.5_RP * lmesh%Fscale(fp,ke)
!             del_flux(fp,DENS_VID) = tmp1 * del_flux_save_v(fpp,DENS_VID,f,ke)
!             del_flux(fp,RHOT_VID) = tmp1 * del_flux_save_v(fpp,RHOT_VID,f,ke)
!             del_flux(fp,MOMZ_VID) = tmp1 * del_flux_save_v(fpp,MOMZ_VID,f,ke)
!             del_flux(fp,MOMX_VID) = tmp1 * del_flux_save_v(fpp,MOMX_VID,f,ke)
!             del_flux(fp,MOMY_VID) = tmp1 * del_flux_save_v(fpp,MOMY_VID,f,ke)
!           end do
!         end if
!       end do

!       call elem3D_optr%Lift_var5( del_flux, Lift )
!       RGsqrt(:) = 1.0_RP / lmesh%Gsqrt(:,ke)
!       do iv=1, 5
!         PRGVAR_dt(:,ke,iv) = PRGVAR_dt(:,ke,iv) - RGsqrt(:) * Lift(:,iv)
!       end do
!     end do
!     !$omp end do
!     !$omp end parallel

!     call PROF_rapend('cal_dyn_tend_bndflux1', 3)
!     return
!   end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_add_bnd_contrib2_generalvc

! !OCL SERIAL
!   subroutine calc_bnd_flux_core( del_flux, del_flux_save, &
!     DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES,        & ! (in)
!     DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot, & ! (in)
!     Gsqrt, G13, G23,                                   & ! (in)
!     nx, ny, nz, Fscale, iM, iP, Nfp, Nfp_s, NfpTot, lmesh, elem )
!     integer, intent(in) :: Nfp
!     integer, intent(in) :: NfpTot
!     class(LocalMesh3D), intent(in) :: lmesh
!     class(ElementBase3D), intent(in) :: elem
!     real(RP), intent(out) :: del_flux(NfpTot,5)
!     real(RP), intent(out) :: del_flux_save(Nfp,5)
!     real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
!     real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
!     real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
!     real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)
!     real(RP), intent(inout) ::  DPRES(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  THERM_hyd(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  Rtot (elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  CVtot(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  CPtot(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
!     real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
!     real(RP), intent(in) :: nx(NfpTot)
!     real(RP), intent(in) :: ny(NfpTot)
!     real(RP), intent(in) :: nz(NfpTot)
!     real(RP), intent(in) :: Fscale(NfpTot)
!     integer, intent(in) :: iM(Nfp)
!     integer, intent(in) :: iP(Nfp)
!     integer, intent(in) :: Nfp_s
!     !-
!     real(RP) :: Vel(Nfp,2), alpha(Nfp)
!     real(RP) :: DPRES_(Nfp,2)
!     real(RP) :: GsqrtDens(Nfp,2)
!     real(RP) :: GsqrtRhot(Nfp,2)
!     real(RP) :: GsqrtDDENS(Nfp,2)
!     real(RP) :: GsqrtMOMX(Nfp,2)
!     real(RP) :: GsqrtMOMY(Nfp,2)
!     real(RP) :: GsqrtMOMZ(Nfp,2)
!     real(RP) :: GsqrtDRHOT(Nfp,2)
!     real(RP) :: Cs2_(Nfp,2)
!     real(RP) :: Gsqrt_(Nfp,2)
!     real(RP) :: GsqrtV_(Nfp,2)
!     real(RP) :: RGsqrtV(Nfp,2)
!     real(RP) :: G13_(Nfp,2)
!     real(RP) :: G23_(Nfp,2)    
!     implicit none

!     real(RP) :: Gnn_M, Gnn_P

!     real(RP) :: gamm

!     integer, parameter :: IN = 1
!     integer, parameter :: EX = 2

!     real(RP) :: tmp1, tmp2, tmp3, tmp4, tmp01, tmp02, tmp03
!     real(RP) :: del_flux_tmp_mom(Nfp,3)

!     integer :: fp, fpp
!     integer :: iv

!     integer :: iim, iip
!     !------------------------

!     gamm  = CPDry / CvDry

!     !-
! !OCL PREFETCH
!     do fp=1, Nfp
!       iim = iM(fp)
!       Gsqrt_(fp,IN) = Gsqrt(iim)
!       GsqrtV_(fp,IN) = Gsqrt_(fp,IN)
!       RGsqrtV(fp,IN) = 1.0_RP / GsqrtV_(fp,IN)  

!       GsqrtDDENS(fp,IN) = Gsqrt_(fp,IN) * DDENS_(iim)
!       GsqrtDens(fp,IN) = GsqrtDDENS(fp,IN) + Gsqrt_(fp,IN) * DENS_hyd(iim)

!       GsqrtDRHOT(fp,IN) = Gsqrt_(fp,IN) * DRHOT_(iim)
!       GsqrtRhot(fp,IN) = GsqrtDRHOT(fp,IN) + Gsqrt_(fp,IN) * THERM_hyd(iim) 
!     end do
! !OCL PREFETCH
!     do fp=1, Nfp
!       iim = iM(fp)

!       DPRES_(fp,IN) = DPRES(iim)
!       Cs2_(fp,IN) = gamm * ( PRES_hyd(iim) + DPRES_(fp,IN) ) * Gsqrt_(fp,IN) / GsqrtDens(fp,IN)
!     end do
! !OCL PREFETCH
!     do fp=1, Nfp
!       iim = iM(fp)
!       GsqrtMOMX (fp,IN) = Gsqrt_(fp,IN) * MOMX_ (iim)
!       GsqrtMOMY (fp,IN) = Gsqrt_(fp,IN) * MOMY_ (iim)
!       GsqrtMOMZ (fp,IN) = Gsqrt_(fp,IN) * MOMZ_ (iim)
      
!       G13_(fp,IN) = G13(iim)
!       G23_(fp,IN) = G23(iim)
!     end do
!     do fp=1, Nfp
!       fpp = Nfp_s + fp
!       Vel(fp,IN) = ( GsqrtMOMX(fp,IN) * nx(fpp) + GsqrtMOMY(fp,IN) * ny(fpp)                        &
!                 + ( ( GsqrtMOMZ(fp,IN) * RGsqrtV(fp,IN)                                             &
!                     + G13_(fp,IN) * GsqrtMOMX(fp,IN) + G23_(fp,IN) * GsqrtMOMY(fp,IN) ) * nz(fpp) ) &
!                 ) / GsqrtDens(fp,IN)
!     end do

!     !-
! !OCL PREFETCH
!     do fp=1, Nfp
!       iip = iP(fp)

!       Gsqrt_(fp,EX) = Gsqrt(iip)
!       GsqrtV_(fp,EX) = Gsqrt_(fp,EX)
!       RGsqrtV(fp,EX) = 1.0_RP / GsqrtV_(fp,EX)  


!       GsqrtDDENS(fp,EX) = Gsqrt_(fp,EX) * DDENS_(iip)
!       GsqrtDens(fp,EX) = GsqrtDDENS(fp,EX) + Gsqrt_(fp,EX) * DENS_hyd(iip)

!       GsqrtDRHOT(fp,EX) = Gsqrt_(fp,EX) * DRHOT_(iip)
!       GsqrtRhot(fp,EX) = GsqrtDRHOT(fp,EX) + Gsqrt_(fp,EX) * THERM_hyd(iip)
!     end do
! !OCL PREFETCH
!     do fp=1, Nfp
!       iip = iP(fp)
!       DPRES_(fp,EX) = DPRES(iip)
!       Cs2_(fp,EX) = gamm * ( PRES_hyd(iip) + DPRES_(fp,EX) ) * Gsqrt_(fp,EX) / GsqrtDens(fp,EX)
!     end do
! !OCL PREFETCH
!     do fp=1, Nfp
!       iip = iP(fp)
!       GsqrtMOMX (fp,EX) = Gsqrt_(fp,EX) * MOMX_ (iip)
!       GsqrtMOMY (fp,EX) = Gsqrt_(fp,EX) * MOMY_ (iip)
!       GsqrtMOMZ (fp,EX) = Gsqrt_(fp,EX) * MOMZ_ (iip)

!       G13_(fp,EX) = G13(iip)
!       G23_(fp,EX) = G23(iip)
!     end do
!     do fp=1, Nfp
!       fpp = Nfp_s + fp
!       Vel(fp,EX) = ( GsqrtMOMX(fp,EX) * nx(fpp) + GsqrtMOMY(fp,EX) * ny(fpp)                        &
!                 + ( ( GsqrtMOMZ(fp,EX) * RGsqrtV(fp,EX)                                             &
!                     + G13_(fp,EX) * GsqrtMOMX(fp,EX) + G23_(fp,EX) * GsqrtMOMY(fp,EX) ) * nz(fpp) ) &
!                 ) / GsqrtDens(fp,EX)
!     end do
        
!     !-
!     do fp=1, Nfp
!       fpp = Nfp_s + fp
!       tmp1 = abs( nx(fpp) ) + abs( ny(fpp) )
!       Gnn_M = tmp1 &
!               + ( 1.0_RP * RGsqrtV(fp,IN)**2 + G13_(fp,IN)**2 + G23_(fp,IN)**2 ) * abs( nz(fpp) )
!       Gnn_P = tmp1 &
!               + ( 1.0_RP * RGsqrtV(fp,EX)**2 + G13_(fp,EX)**2 + G23_(fp,EX)**2 ) * abs( nz(fpp) )

!       alpha(fp) = max( sqrt( Gnn_M * Cs2_(fp,IN) ) + abs(Vel(fp,IN)), &
!                        sqrt( Gnn_P * Cs2_(fp,EX) ) + abs(Vel(fp,EX))  )
!     end do
!     do fp=1, Nfp
!       fpp = Nfp_s + fp
!       tmp1 = Fscale(fpp) * 0.5_RP

!       !-
!       tmp2 = GsqrtDens(fp,EX) * Vel(fp,EX) - GsqrtDens(fp,IN) * Vel(fp,IN)
!       tmp3 = - alpha(fp) * ( GsqrtDDENS(fp,EX) - GsqrtDDENS(fp,IN) )

!       del_flux(fpp,DENS_VID) = tmp1 * ( tmp2 + tmp3 )
!       del_flux_save(fp,DENS_VID) = tmp2 - tmp3

!       !-
!       tmp2 = GsqrtRhot(fp,EX) * Vel(fp,EX) - GsqrtRhot(fp,IN) * Vel(fp,IN)
!       tmp3 = - alpha(fp) * ( GsqrtDRHOT(fp,EX) - GsqrtDRHOT(fp,IN) )

!       del_flux(fpp,RHOT_VID) = tmp1 * ( tmp2 + tmp3 )
!       del_flux_save(fp,RHOT_VID) = tmp2 - tmp3
!     end do

!     do fp=1, Nfp
!       fpp = Nfp_s + fp
!       tmp3 = Gsqrt_(fp,EX) * dpres_(fp,EX)
!       tmp4 = Gsqrt_(fp,IN) * dpres_(fp,IN)

!       del_flux_tmp_mom(fp,1) = &
!               (  tmp3 * RGsqrtV(fp,EX)              &
!                 - tmp4 * RGsqrtV(fp,IN) ) * nz(fpp)
      
!       del_flux_tmp_mom(fp,2) = &
!             ( nx(fpp) + G13_(fp,EX) * nz(fpp) ) * tmp3 &
!           - ( nx(fpp) + G13_(fp,IN) * nz(fpp) ) * tmp4
      
!       del_flux_tmp_mom(fp,3) = &
!             ( ny(fpp) + G23_(fp,EX) * nz(fpp) ) * tmp3 &
!           - ( ny(fpp) + G23_(fp,IN) * nz(fpp) ) * tmp4
!     end do

!     do fp=1, Nfp
!       fpp = Nfp_s + fp
!       tmp1 = Fscale(fpp) * 0.5_RP

!       tmp2 = GsqrtMOMZ(fp,EX) * Vel(fp,EX) - GsqrtMOMZ(fp,IN) * Vel(fp,IN)
!       tmp3 = - alpha(fp) * ( GsqrtMOMZ(fp,EX) - GsqrtMOMZ(fp,IN) )

!       del_flux(fpp,MOMZ_VID) = tmp1 * ( tmp2 + del_flux_tmp_mom(fp,1) + tmp3 )
!       del_flux_save(fp,MOMZ_VID) = tmp2 + del_flux_tmp_mom(fp,1) - tmp3
!     end do

!     !-
!     do fp=1, Nfp
!       fpp = Nfp_s + fp
!       tmp1 = Fscale(fpp) * 0.5_RP

!       !-
!       tmp2 = GsqrtMOMX(fp,EX) * Vel(fp,EX) - GsqrtMOMX(fp,IN) * Vel(fp,IN)
!       tmp3 = - alpha(fp) * ( GsqrtMOMX(fp,EX) - GsqrtMOMX(fp,IN) )

!       del_flux(fpp,MOMX_VID) = tmp1 * ( tmp2 + del_flux_tmp_mom(fp,2) + tmp3 )
!       del_flux_save(fp,MOMX_VID) = tmp2 + del_flux_tmp_mom(fp,2) - tmp3

!       !-
!       tmp2 = GsqrtMOMY(fp,EX) * Vel(fp,EX) - GsqrtMOMY(fp,IN) * Vel(fp,IN)
!       tmp3 = - alpha(fp) * ( GsqrtMOMY(fp,EX) - GsqrtMOMY(fp,IN) )

!       del_flux(fpp,MOMY_VID) = tmp1 * ( tmp2 + del_flux_tmp_mom(fp,3) + tmp3 )
!       del_flux_save(fp,MOMY_VID) = tmp2 + del_flux_tmp_mom(fp,3) - tmp3
!     end do

!     return
!   end subroutine calc_bnd_flux_core

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc( &
    del_flux,                                                        & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES,                      & ! (in)
    DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot,               & ! (in)
    Gsqrt, G13, G23, nx, ny, nz,                                     & ! (in)
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D                       ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DPRES(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  THERM_hyd(elem%Np*lmesh%NeA)
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
    
    integer :: ke, fp, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: Vel(elem%NfpTot,2), alpha(elem%NfpTot)
    real(RP) :: DPRES_(elem%NfpTot,2)
    real(RP) :: GsqrtDens(elem%NfpTot,2)
    real(RP) :: GsqrtRhot(elem%NfpTot,2)
    real(RP) :: GsqrtDDENS(elem%NfpTot,2)
    real(RP) :: GsqrtMOMX(elem%NfpTot,2)
    real(RP) :: GsqrtMOMY(elem%NfpTot,2)
    real(RP) :: GsqrtMOMZ(elem%NfpTot,2)
    real(RP) :: GsqrtDRHOT(elem%NfpTot,2)
    real(RP) :: Phyd_(elem%NfpTot,2)
    real(RP) :: Gsqrt_(elem%NfpTot,2)
    real(RP) :: GsqrtV_(elem%NfpTot,2)
    real(RP) :: RGsqrtV(elem%NfpTot,2)
    real(RP) :: G13_(elem%NfpTot,2)
    real(RP) :: G23_(elem%NfpTot,2)
    real(RP) :: Gnn_M, Gnn_P

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR
    
    integer, parameter :: IN = 1
    integer, parameter :: EX = 2

    real(RP) :: tmp1, tmp2, tmp3, tmp4, tmp01, tmp02, tmp03
    real(RP) :: del_flux_tmp_mom(elem%NfpTot,3)
    !------------------------------------------------------------------------

    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( &
    !$omp ke, iM, iP, ke2D, fp,                                                         &
    !$omp alpha, Vel,                                                                   &
    !$omp dpres_, GsqrtDens, GsqrtRhot,                                                 &
    !$omp GsqrtDDENS, GsqrtMOMX, GsqrtMOMY, GsqrtMOMZ, GsqrtDRHOT,                      &
    !$omp Phyd_,                                                                        &
    !$omp Gsqrt_, GsqrtV_, RGsqrtV, G13_, G23_,                                         &
    !$omp Gnn_P, Gnn_M, tmp1, tmp2, tmp3, tmp4, tmp01, tmp02, tmp03, del_flux_tmp_mom                                          )
!OCL PREFETCH
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      Gsqrt_(:,IN) = Gsqrt(iM)
      Gsqrt_(:,EX) = Gsqrt(iP)
      GsqrtV_(:,:) = Gsqrt_(:,:)
      RGsqrtV(:,:) = 1.0_RP / GsqrtV_(:,:)      

      G13_(:,IN) = G13(iM)
      G13_(:,EX) = G13(iP)
      G23_(:,IN) = G23(iM)
      G23_(:,EX) = G23(iP)

      GsqrtDDENS(:,IN) = Gsqrt_(:,IN) * DDENS_(iM)
      GsqrtDDENS(:,EX) = Gsqrt_(:,EX) * DDENS_(iP)
      GsqrtMOMX (:,IN) = Gsqrt_(:,IN) * MOMX_ (iM)
      GsqrtMOMX (:,EX) = Gsqrt_(:,EX) * MOMX_ (iP)
      GsqrtMOMY (:,IN) = Gsqrt_(:,IN) * MOMY_ (iM)
      GsqrtMOMY (:,EX) = Gsqrt_(:,EX) * MOMY_ (iP)
      GsqrtMOMZ (:,IN) = Gsqrt_(:,IN) * MOMZ_ (iM)
      GsqrtMOMZ (:,EX) = Gsqrt_(:,EX) * MOMZ_ (iP)
      GsqrtDRHOT(:,IN) = Gsqrt_(:,IN) * DRHOT_(iM)
      GsqrtDRHOT(:,EX) = Gsqrt_(:,EX) * DRHOT_(iP)

      Phyd_(:,IN) = PRES_hyd(iM)
      Phyd_(:,EX) = PRES_hyd(iP)
      DPRES_(:,IN) = DPRES(iM)
      DPRES_(:,EX) = DPRES(iP)

      GsqrtDens(:,IN) = GsqrtDDENS(:,IN) + Gsqrt_(:,IN) * DENS_hyd(iM)
      GsqrtDens(:,EX) = GsqrtDDENS(:,EX) + Gsqrt_(:,EX) * DENS_hyd(iP)

      GsqrtRhot(:,IN) = Gsqrt_(:,IN) * THERM_hyd(iM) + GsqrtDRHOT(:,IN)
      GsqrtRhot(:,EX) = Gsqrt_(:,EX) * THERM_hyd(iP) + GsqrtDRHOT(:,EX)

      Vel(:,IN) = ( GsqrtMOMX(:,IN) * nx(:,ke) + GsqrtMOMY(:,IN) * ny(:,ke)                  &
                + ( ( GsqrtMOMZ(:,IN) * RGsqrtV(:,IN)                                        &
                    + G13_(:,IN) * GsqrtMOMX(:,IN) + G23_(:,IN) * GsqrtMOMY(:,IN) ) * nz(:,ke) ) &
                ) / GsqrtDens(:,IN)
      Vel(:,EX) = ( GsqrtMOMX(:,EX) * nx(:,ke) + GsqrtMOMY(:,EX) * ny(:,ke)                  &
                + ( ( GsqrtMOMZ(:,EX) * RGsqrtV(:,EX)                                        &
                    + G13_(:,EX) * GsqrtMOMX(:,EX) + G23_(:,EX) * GsqrtMOMY(:,EX) ) * nz(:,ke) ) &
                ) / GsqrtDens(:,EX)
            
      do fp=1, elem%NfpTot
        tmp1 = abs( nx(fp,ke) ) + abs( ny(fp,ke) )
        Gnn_M = tmp1 &
                + ( 1.0_RP * RGsqrtV(fp,IN)**2 + G13_(fp,IN)**2 + G23_(fp,IN)**2 ) * abs( nz(fp,ke) )
        Gnn_P = tmp1 &
                + ( 1.0_RP * RGsqrtV(fp,EX)**2 + G13_(fp,EX)**2 + G23_(fp,EX)**2 ) * abs( nz(fp,ke) )

        alpha(fp) = max( sqrt( Gnn_M * gamm * ( Phyd_(fp,IN) + dpres_(fp,IN) ) * Gsqrt_(fp,IN) / GsqrtDens(fp,IN) ) + abs(Vel(fp,IN)), &
                         sqrt( Gnn_P * gamm * ( Phyd_(fp,EX) + dpres_(fp,EX) ) * Gsqrt_(fp,EX) / GsqrtDens(fp,EX) ) + abs(Vel(fp,EX))  )
      end do
      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP

        tmp2 = - alpha(fp) * ( GsqrtDDENS(fp,EX) - GsqrtDDENS(fp,IN) )
        del_flux(fp,DENS_VID,ke) = tmp1 * ( &
                        GsqrtDens(fp,EX) * Vel(fp,EX) - GsqrtDens(fp,IN) * Vel(fp,IN)  &
                      + tmp2  )

        tmp2 = - alpha(fp) * ( GsqrtDRHOT(fp,EX) - GsqrtDRHOT(fp,IN) )
        del_flux(fp,RHOT_VID,ke) = tmp1 * ( &
                        GsqrtRhot(fp,EX) * Vel(fp,EX) - GsqrtRhot(fp,IN) * Vel(fp,IN) &
                      + tmp2 )
      end do

      do fp=1, elem%NfpTot
        tmp3 = Gsqrt_(fp,EX) * dpres_(fp,EX)
        tmp4 = Gsqrt_(fp,IN) * dpres_(fp,IN)

        del_flux_tmp_mom(fp,1) = &
                (  tmp3 * RGsqrtV(fp,EX)               &
                 - tmp4 * RGsqrtV(fp,IN) ) * nz(fp,ke)
        
        del_flux_tmp_mom(fp,2) = &
              ( nx(fp,ke) + G13_(fp,EX) * nz(fp,ke) ) * tmp3   &
            - ( nx(fp,ke) + G13_(fp,IN) * nz(fp,ke) ) * tmp4
        
        del_flux_tmp_mom(fp,3) = &
              ( ny(fp,ke) + G23_(fp,EX) * nz(fp,ke) ) * tmp3   &
            - ( ny(fp,ke) + G23_(fp,IN) * nz(fp,ke) ) * tmp4
      end do
      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP

        tmp2 = - alpha(fp) * ( GsqrtMOMZ(fp,EX) - GsqrtMOMZ(fp,IN) )
        del_flux(fp,MOMZ_VID,ke) = tmp1 * ( &
                        GsqrtMOMZ(fp,EX) * Vel(fp,EX) - GsqrtMOMZ(fp,IN) * Vel(fp,IN) &
                      + del_flux_tmp_mom(fp,1) &
                      + tmp2 )
      end do
      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP

        tmp2 = - alpha(fp) * ( GsqrtMOMX(fp,EX) - GsqrtMOMX(fp,IN) )
        del_flux(fp,MOMX_VID,ke) = tmp1 * ( &
                        GsqrtMOMX(fp,EX) * Vel(fp,EX) - GsqrtMOMX(fp,IN) * Vel(fp,IN) &
                      + del_flux_tmp_mom(fp,2) &
                      + tmp2 )

        tmp2 = - alpha(fp) * ( GsqrtMOMY(fp,EX) - GsqrtMOMY(fp,IN) )
        del_flux(fp,MOMY_VID,ke) = tmp1 * ( &
                        GsqrtMOMY(fp,EX) * Vel(fp,EX) - GsqrtMOMY(fp,IN) * Vel(fp,IN) &
                      + del_flux_tmp_mom(fp,3) &
                      + tmp2 )
      end do
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_LMARS( &
    del_flux,                                                        & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES, DENS_hyd, PRES_hyd,  & ! (in)
    Rtot, CVtot, CPtot,                                              & ! (in)
    Gsqrt, G13, G23, nx, ny, nz,                                     & ! (in)
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D                       ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DPRES(elem%Np*lmesh%NeA)
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
    
    integer :: ke, fp, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: Vel(elem%NfpTot,2), alpha(elem%NfpTot), alpha2(elem%NfpTot)
    real(RP) :: DPRES_(elem%NfpTot,2)
    real(RP) :: GsqrtDens(elem%NfpTot,2)
    real(RP) :: GsqrtRhot(elem%NfpTot,2)
    real(RP) :: GsqrtDDENS(elem%NfpTot,2)
    real(RP) :: GsqrtMOMX(elem%NfpTot,2)
    real(RP) :: GsqrtMOMY(elem%NfpTot,2)
    real(RP) :: GsqrtMOMZ(elem%NfpTot,2)
    real(RP) :: GsqrtDRHOT(elem%NfpTot,2)
    real(RP) :: Phyd_(elem%NfpTot,2)
    real(RP) :: Gsqrt_(elem%NfpTot,2)
    real(RP) :: GsqrtV_(elem%NfpTot,2)
    real(RP) :: RGsqrtV(elem%NfpTot,2)
    real(RP) :: G13_(elem%NfpTot,2)
    real(RP) :: G23_(elem%NfpTot,2)
    real(RP) :: Gnn_M, Gnn_P

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR
    
    integer, parameter :: IN = 1
    integer, parameter :: EX = 2

    real(RP) :: tmp1, tmp2, tmp3, tmp4, tmp01, tmp02, tmp03
    real(RP) :: del_flux_tmp_mom(elem%NfpTot,3)

    real(RP) :: Vel_n(elem%NfpTot)
    real(RP) :: Cs(elem%NfpTot,2)
    !------------------------------------------------------------------------

    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( &
    !$omp ke, iM, iP, ke2D, fp,                                                         &
    !$omp alpha, alpha2, Vel, Vel_n, Cs,                                                &
    !$omp dpres_, GsqrtDens, GsqrtRhot,                                                 &
    !$omp GsqrtDDENS, GsqrtMOMX, GsqrtMOMY, GsqrtMOMZ, GsqrtDRHOT,                      &
    !$omp Phyd_,                                                                        &
    !$omp Gsqrt_, GsqrtV_, RGsqrtV, G13_, G23_,                                         &
    !$omp Gnn_P, Gnn_M, tmp1, tmp2, tmp3, tmp4, tmp01, tmp02, tmp03, del_flux_tmp_mom                                          )
!OCL PREFETCH
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      Gsqrt_(:,IN) = Gsqrt(iM)
      Gsqrt_(:,EX) = Gsqrt(iP)
      GsqrtV_(:,:) = Gsqrt_(:,:)
      RGsqrtV(:,:) = 1.0_RP / GsqrtV_(:,:)      

      G13_(:,IN) = G13(iM)
      G13_(:,EX) = G13(iP)
      G23_(:,IN) = G23(iM)
      G23_(:,EX) = G23(iP)

      GsqrtDDENS(:,IN) = Gsqrt_(:,IN) * DDENS_(iM)
      GsqrtDDENS(:,EX) = Gsqrt_(:,EX) * DDENS_(iP)
      GsqrtMOMX (:,IN) = Gsqrt_(:,IN) * MOMX_ (iM)
      GsqrtMOMX (:,EX) = Gsqrt_(:,EX) * MOMX_ (iP)
      GsqrtMOMY (:,IN) = Gsqrt_(:,IN) * MOMY_ (iM)
      GsqrtMOMY (:,EX) = Gsqrt_(:,EX) * MOMY_ (iP)
      GsqrtMOMZ (:,IN) = Gsqrt_(:,IN) * MOMZ_ (iM)
      GsqrtMOMZ (:,EX) = Gsqrt_(:,EX) * MOMZ_ (iP)
      GsqrtDRHOT(:,IN) = Gsqrt_(:,IN) * DRHOT_(iM)
      GsqrtDRHOT(:,EX) = Gsqrt_(:,EX) * DRHOT_(iP)

      Phyd_(:,IN) = PRES_hyd(iM)
      Phyd_(:,EX) = PRES_hyd(iP)
      DPRES_(:,IN) = DPRES(iM)
      DPRES_(:,EX) = DPRES(iP)

      GsqrtDens(:,IN) = GsqrtDDENS(:,IN) + Gsqrt_(:,IN) * DENS_hyd(iM)
      GsqrtDens(:,EX) = GsqrtDDENS(:,EX) + Gsqrt_(:,EX) * DENS_hyd(iP)

      GsqrtRhot(:,:) = Gsqrt_(:,:) * P0ovR * (Phyd_(:,:) * rP0)**rgamm + GsqrtDRHOT(:,:)

      Vel(:,IN) = ( GsqrtMOMX(:,IN) * nx(:,ke) + GsqrtMOMY(:,IN) * ny(:,ke)                  &
                + ( ( GsqrtMOMZ(:,IN) * RGsqrtV(:,IN)                                        &
                    + G13_(:,IN) * GsqrtMOMX(:,IN) + G23_(:,IN) * GsqrtMOMY(:,IN) ) * nz(:,ke) ) &
                ) / GsqrtDens(:,IN)
      Vel(:,EX) = ( GsqrtMOMX(:,EX) * nx(:,ke) + GsqrtMOMY(:,EX) * ny(:,ke)                  &
                + ( ( GsqrtMOMZ(:,EX) * RGsqrtV(:,EX)                                        &
                    + G13_(:,EX) * GsqrtMOMX(:,EX) + G23_(:,EX) * GsqrtMOMY(:,EX) ) * nz(:,ke) ) &
                ) / GsqrtDens(:,EX)
            
      do fp=1, elem%NfpTot
        tmp1 = abs( nx(fp,ke) ) + abs( ny(fp,ke) )
        Gnn_M = tmp1 &
                + ( 1.0_RP * RGsqrtV(fp,IN)**2 + G13_(fp,IN)**2 + G23_(fp,IN)**2 ) * abs( nz(fp,ke) )
        Gnn_P = tmp1 &
                + ( 1.0_RP * RGsqrtV(fp,EX)**2 + G13_(fp,EX)**2 + G23_(fp,EX)**2 ) * abs( nz(fp,ke) )

        Cs(fp,IN) = sqrt( Gnn_M * gamm * ( Phyd_(fp,IN) + dpres_(fp,IN) ) * Gsqrt_(fp,IN) / GsqrtDens(fp,IN) )
        Cs(fp,EX) = sqrt( Gnn_P * gamm * ( Phyd_(fp,EX) + dpres_(fp,EX) ) * Gsqrt_(fp,EX) / GsqrtDens(fp,EX) )
        alpha(fp) = max( Cs(fp,IN) + abs(Vel(fp,IN)), Cs(fp,EX) + abs(Vel(fp,EX))  )
      end do
      Vel_n(:) = 0.5_RP * ( &
                 ( Vel(:,IN) + Vel(:,EX) ) &
               - 2.0_RP * ( DPRES_(:,EX) - DPRES_(:,IN) ) /  ( 0.5_RP * ( GsqrtDens(:,IN) / Gsqrt_(:,IN) * Cs(:,IN) + GsqrtDens(:,EX) / Gsqrt_(:,EX)  * Cs(:,EX) ) ) &
               )


      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP
        tmp3 = Vel_n(fp) - 2.0_RP * Vel(fp,IN)

        tmp2 = - abs( Vel_n(fp) ) * ( GsqrtDDENS(fp,EX) - GsqrtDDENS(fp,IN) )
        del_flux(fp,DENS_VID,ke) = tmp1 * ( &
                        GsqrtDens(fp,EX) * Vel_n(fp) + GsqrtDens(fp,IN) * tmp3 &
                      + tmp2  )

        tmp2 = - abs( Vel_n(fp) ) * ( GsqrtDRHOT(fp,EX) - GsqrtDRHOT(fp,IN) )
        del_flux(fp,RHOT_VID,ke) = tmp1 * ( &
                        GsqrtRhot(fp,EX) * Vel_n(fp) + GsqrtRhot(fp,IN) * tmp3 &
                      + tmp2 )
      end do

      do fp=1, elem%NfpTot
        tmp1 = 0.5_RP * ( Cs(fp,IN) + Cs(fp,EX) )
        tmp3 = Gsqrt_(fp,EX) * dpres_(fp,EX)
        tmp4 = Gsqrt_(fp,IN) * dpres_(fp,IN)

        del_flux_tmp_mom(fp,1) = &
                (  tmp3 * RGsqrtV(fp,EX)               &
                 - tmp4 * RGsqrtV(fp,IN) ) * nz(fp,ke) &
              - tmp1 * ( GsqrtMOMZ(fp,EX) - GsqrtMOMZ(fp,IN) ) !* abs(nz(fp,ke))
        
        del_flux_tmp_mom(fp,2) = &
              ( nx(fp,ke) + G13_(fp,EX) * nz(fp,ke) ) * tmp3   &
            - ( nx(fp,ke) + G13_(fp,IN) * nz(fp,ke) ) * tmp4   &
            - tmp1 * ( GsqrtMOMX(fp,EX) - GsqrtMOMX(fp,IN) ) !* abs(nx(fp,ke))
        
        del_flux_tmp_mom(fp,3) = &
              ( ny(fp,ke) + G23_(fp,EX) * nz(fp,ke) ) * tmp3   &
            - ( ny(fp,ke) + G23_(fp,IN) * nz(fp,ke) ) * tmp4   &
            - tmp1 * ( GsqrtMOMY(fp,EX) - GsqrtMOMY(fp,IN) ) !* abs(ny(fp,ke))

      end do
      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP
        tmp3 = Vel_n(fp) - 2.0_RP * Vel(fp,IN)

        tmp2 = - abs( Vel_n(fp) ) * ( GsqrtMOMZ(fp,EX) - GsqrtMOMZ(fp,IN) )
        del_flux(fp,MOMZ_VID,ke) = tmp1 * ( &
                        GsqrtMOMZ(fp,EX) * Vel_n(fp) + GsqrtMOMZ(fp,IN) * tmp3 &
                      + del_flux_tmp_mom(fp,1) &
                      + tmp2 )
      end do
      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP
        tmp3 = Vel_n(fp) - 2.0_RP * Vel(fp,IN)

        tmp2 = - abs( Vel_n(fp) ) * ( GsqrtMOMX(fp,EX) - GsqrtMOMX(fp,IN) )
        del_flux(fp,MOMX_VID,ke) = tmp1 * ( &
                        GsqrtMOMX(fp,EX) * Vel_n(fp) + GsqrtMOMX(fp,IN) * tmp3 &
                      + del_flux_tmp_mom(fp,2) &
                      + tmp2 )

        tmp2 = - abs( Vel_n(fp) ) * ( GsqrtMOMY(fp,EX) - GsqrtMOMY(fp,IN) )
        del_flux(fp,MOMY_VID,ke) = tmp1 * ( &
                        GsqrtMOMY(fp,EX) * Vel_n(fp) + GsqrtMOMY(fp,IN) * tmp3 &
                      + del_flux_tmp_mom(fp,3) &
                      + tmp2 )
      end do
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_LMARS  
!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc_asis( &
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
    
    integer :: ke, iP(elem%NfpTot), iM(elem%NfpTot)
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

    !$omp parallel do private( &
    !$omp ke, iM, iP, ke2d,                                                             &
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

      Gsqrt_M(:) = Gsqrt(iM)
      Gsqrt_P(:) = Gsqrt(iP)

      rgam2_M(:) = 1.0_RP / gam(iM)**2
      rgam2_P(:) = 1.0_RP / gam(iP)**2
      GsqrtV_M(:) = Gsqrt_M(:) * rgam2_M(:) / GsqrtH(iM2Dto3D(:),ke2D)
      GsqrtV_P(:) = Gsqrt_P(:) * rgam2_P(:) / GsqrtH(iM2Dto3D(:),ke2D)

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

      Gxz_M(:) = rgam2_M(:) * ( G11(iM2Dto3D(:),ke2D) * G13_M(:) + G12(iM2Dto3D(:),ke2D) * G23_M(:) )
      Gxz_P(:) = rgam2_P(:) * ( G11(iM2Dto3D(:),ke2D) * G13_P(:) + G12(iM2Dto3D(:),ke2D) * G23_P(:) )

      Gyz_M(:) = rgam2_M(:) * ( G12(iM2Dto3D(:),ke2D) * G13_M(:) + G22(iM2Dto3D(:),ke2D) * G23_M(:) )
      Gyz_P(:) = rgam2_P(:) * ( G12(iM2Dto3D(:),ke2D) * G13_P(:) + G22(iM2Dto3D(:),ke2D) * G23_P(:) )
               
      G1n_M(:) = rgam2_M(:) * ( G11(iM2Dto3D(:),ke2D) * nx(:,ke) + G12(iM2Dto3D(:),ke2D) * ny(:,ke) )
      G2n_M(:) = rgam2_P(:) * ( G12(iM2Dto3D(:),ke2D) * nx(:,ke) + G22(iM2Dto3D(:),ke2D) * ny(:,ke) )

      Gnn_M(:)  = rgam2_M(:) * ( G11(iM2Dto3D(:),ke2D) * abs( nx(:,ke) ) + G22(iM2Dto3D(:),ke2D) * abs( ny(:,ke) ) ) &
                + ( 1.0_RP / GsqrtV_M(:)**2 + G13_M(:) * Gxz_M(:) + G23_M(:) * Gyz_M(:) ) * abs( nz(:,ke) )
      Gnn_P(:)  = rgam2_P(:) * ( G11(iM2Dto3D(:),ke2D) * abs( nx(:,ke) ) + G22(iM2Dto3D(:),ke2D) * abs( ny(:,ke) ) ) &
                + ( 1.0_RP / GsqrtV_P(:)**2 + G13_P(:) * Gxz_P(:) + G23_P(:) * Gyz_P(:) ) * abs( nz(:,ke) )

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
                    + (  Gsqrt_P(:) * ( G1n_M(:) + Gxz_P(:) * nz(:,ke)) * dpresP(:)   &
                       - Gsqrt_M(:) * ( G1n_M(:) + Gxz_M(:) * nz(:,ke)) * dpresM(:) ) &
                    - alpha(:) * ( GsqrtMOMX_P(:) - GsqrtMOMX_M(:) )                  )

      del_flux(:,ke,MOMY_VID ) = 0.5_RP * ( &
                    ( GsqrtMOMY_P(:) * VelP(:) - GsqrtMOMY_M(:) * VelM(:) )            &
                    + (  Gsqrt_P(:) * ( G2n_M(:) + Gyz_P(:) * nz(:,ke) ) * dpresP(:)   &
                       - Gsqrt_M(:) * ( G2n_M(:) + Gyz_M(:) * nz(:,ke) ) * dpresM(:) ) &
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
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc_asis

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc( &
    del_flux,                                                         & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES,                       & ! (in)
    DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot,                & ! (in)
    Gsqrt, G11, G12, G22, GsqrtH, gam, G13, G23, nx, ny, nz,          & ! (in)
    vmapM, vmapP, iM2Dto3D, lmesh, elem, lmesh2D, elem2D              ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DPRES(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  THERM_hyd(elem%Np*lmesh%NeA)
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
    
    integer :: ke, fp, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: Vel(elem%NfpTot,2), alpha(elem%NfpTot)
    real(RP) :: DPRES_(elem%NfpTot,2)
    real(RP) :: GsqrtDens(elem%NfpTot,2)
    real(RP) :: GsqrtRhot(elem%NfpTot,2)
    real(RP) :: GsqrtDDENS(elem%NfpTot,2)
    real(RP) :: GsqrtMOMX(elem%NfpTot,2)
    real(RP) :: GsqrtMOMY(elem%NfpTot,2)
    real(RP) :: GsqrtMOMZ(elem%NfpTot,2)
    real(RP) :: GsqrtDRHOT(elem%NfpTot,2)
    real(RP) :: Phyd_(elem%NfpTot,2)
    real(RP) :: Gsqrt_(elem%NfpTot,2)
    real(RP) :: GsqrtV_(elem%NfpTot,2)
    real(RP) :: RGsqrtV(elem%NfpTot,2)
    real(RP) :: G13_(elem%NfpTot,2)
    real(RP) :: G23_(elem%NfpTot,2)
    real(RP) :: Gxz_(elem%NfpTot,2)
    real(RP) :: Gyz_(elem%NfpTot,2)
    real(RP) :: G1n_(elem%NfpTot,2)
    real(RP) :: G2n_(elem%NfpTot,2)

    real(RP) :: Gnn_M, Gnn_P
    real(RP) :: rgam2(elem%NfpTot,2)

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR
    
    integer, parameter :: IN = 1
    integer, parameter :: EX = 2

    real(RP) :: tmp1, tmp2, tmp3, tmp4, tmp01, tmp02, tmp03
    real(RP) :: del_flux_tmp_mom(elem%NfpTot,3)

    real(RP) :: G11_, G12_, G22_
    !------------------------------------------------------------------------

    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( &
    !$omp ke, iM, iP, ke2D, fp,                                                         &
    !$omp alpha, Vel,                                                                   &
    !$omp dpres_, GsqrtDens, GsqrtRhot,                                                 &
    !$omp GsqrtDDENS, GsqrtMOMX, GsqrtMOMY, GsqrtMOMZ, GsqrtDRHOT,                      &
    !$omp Phyd_,                                                                        &
    !$omp Gsqrt_, GsqrtV_, RGsqrtV, G11_, G12_, G22_, G13_, G23_, rgam2,                &
    !$omp Gxz_, Gyz_, G1n_, G2n_,                                                       &
    !$omp Gnn_P, Gnn_M, tmp1, tmp2, tmp3, tmp4, tmp01, tmp02, tmp03, del_flux_tmp_mom   )
!OCL PREFETCH
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      Gsqrt_(:,IN) = Gsqrt(iM)
      Gsqrt_(:,EX) = Gsqrt(iP)

      rgam2(:,IN) = 1.0_RP / gam(iM)**2
      rgam2(:,EX) = 1.0_RP / gam(iP)**2
      GsqrtV_(:,IN) = Gsqrt_(:,IN) * rgam2(:,IN) / GsqrtH(iM2Dto3D(:),ke2D)
      GsqrtV_(:,EX) = Gsqrt_(:,EX) * rgam2(:,EX) / GsqrtH(iM2Dto3D(:),ke2D)
      RGsqrtV(:,:) = 1.0_RP / GsqrtV_(:,:)      

      G13_(:,IN) = G13(iM)
      G13_(:,EX) = G13(iP)
      G23_(:,IN) = G23(iM)
      G23_(:,EX) = G23(iP)

      GsqrtDDENS(:,IN) = Gsqrt_(:,IN) * DDENS_(iM)
      GsqrtDDENS(:,EX) = Gsqrt_(:,EX) * DDENS_(iP)
      GsqrtMOMX (:,IN) = Gsqrt_(:,IN) * MOMX_ (iM)
      GsqrtMOMX (:,EX) = Gsqrt_(:,EX) * MOMX_ (iP)
      GsqrtMOMY (:,IN) = Gsqrt_(:,IN) * MOMY_ (iM)
      GsqrtMOMY (:,EX) = Gsqrt_(:,EX) * MOMY_ (iP)
      GsqrtMOMZ (:,IN) = Gsqrt_(:,IN) * MOMZ_ (iM)
      GsqrtMOMZ (:,EX) = Gsqrt_(:,EX) * MOMZ_ (iP)
      GsqrtDRHOT(:,IN) = Gsqrt_(:,IN) * DRHOT_(iM)
      GsqrtDRHOT(:,EX) = Gsqrt_(:,EX) * DRHOT_(iP)

      Phyd_(:,IN) = PRES_hyd(iM)
      Phyd_(:,EX) = PRES_hyd(iP)
      DPRES_(:,IN) = DPRES(iM)
      DPRES_(:,EX) = DPRES(iP)

      GsqrtDens(:,IN) = GsqrtDDENS(:,IN) + Gsqrt_(:,IN) * DENS_hyd(iM)
      GsqrtDens(:,EX) = GsqrtDDENS(:,EX) + Gsqrt_(:,EX) * DENS_hyd(iP)

      GsqrtRhot(:,IN) = Gsqrt_(:,IN) * THERM_hyd(iM) + GsqrtDRHOT(:,IN)
      GsqrtRhot(:,EX) = Gsqrt_(:,EX) * THERM_hyd(iP) + GsqrtDRHOT(:,EX)

      Vel(:,IN) = ( GsqrtMOMX(:,IN) * nx(:,ke) + GsqrtMOMY(:,IN) * ny(:,ke)                  &
                + ( ( GsqrtMOMZ(:,IN) * RGsqrtV(:,IN)                                        &
                    + G13_(:,IN) * GsqrtMOMX(:,IN) + G23_(:,IN) * GsqrtMOMY(:,IN) ) * nz(:,ke) ) &
                ) / GsqrtDens(:,IN)
      Vel(:,EX) = ( GsqrtMOMX(:,EX) * nx(:,ke) + GsqrtMOMY(:,EX) * ny(:,ke)                  &
                + ( ( GsqrtMOMZ(:,EX) * RGsqrtV(:,EX)                                        &
                    + G13_(:,EX) * GsqrtMOMX(:,EX) + G23_(:,EX) * GsqrtMOMY(:,EX) ) * nz(:,ke) ) &
                ) / GsqrtDens(:,EX)
            
      do fp=1, elem%NfpTot
        G11_ = G11(iM2Dto3D(fp),ke2D);  G12_ = G12(iM2Dto3D(fp),ke2D);  G22_ = G22(iM2Dto3D(fp),ke2D)

        Gxz_(fp,IN) = rgam2(fp,IN) * ( G11_ * G13_(fp,IN) + G12_ * G23_(fp,IN) )
        Gxz_(fp,EX) = rgam2(fp,EX) * ( G11_ * G13_(fp,EX) + G12_ * G23_(fp,EX) )

        Gyz_(fp,IN) = rgam2(fp,IN) * ( G12_ * G13_(fp,IN) + G22_ * G23_(fp,IN) )
        Gyz_(fp,EX) = rgam2(fp,EX) * ( G12_ * G13_(fp,EX) + G22_ * G23_(fp,EX) )

        G1n_(fp,IN) = rgam2(fp,IN) * ( G11_ * nx(fp,ke) + G12_ * ny(fp,ke) )
        G1n_(fp,EX) = rgam2(fp,EX) * ( G11_ * nx(fp,ke) + G12_ * ny(fp,ke) )

        G2n_(fp,IN) = rgam2(fp,IN) * ( G12_ * nx(fp,ke) + G22_ * ny(fp,ke) )
        G2n_(fp,EX) = rgam2(fp,EX) * ( G12_ * nx(fp,ke) + G22_ * ny(fp,ke) )
      end do
      do fp=1, elem%NfpTot
        G11_ = G11(iM2Dto3D(fp),ke2D);  G22_ = G22(iM2Dto3D(fp),ke2D)
        tmp1 = abs( G11_ * nx(fp,ke) ) + abs( G22_ * ny(fp,ke) )

        Gnn_M = rgam2(fp,IN) * tmp1 &
                + ( 1.0_RP * RGsqrtV(fp,IN)**2 + G13_(fp,IN) * Gxz_(fp,IN) + G23_(fp,IN) * Gyz_(fp,IN) ) * abs( nz(fp,ke) )
        Gnn_P = rgam2(fp,EX) * tmp1 &
                + ( 1.0_RP * RGsqrtV(fp,EX)**2 + G13_(fp,EX) * Gxz_(fp,EX) + G23_(fp,EX) * Gyz_(fp,EX) ) * abs( nz(fp,ke) )

        alpha(fp) = max( sqrt( Gnn_M * gamm * ( Phyd_(fp,IN) + dpres_(fp,IN) ) * Gsqrt_(fp,IN) / GsqrtDens(fp,IN) ) + abs(Vel(fp,IN)), &
                         sqrt( Gnn_P * gamm * ( Phyd_(fp,EX) + dpres_(fp,EX) ) * Gsqrt_(fp,EX) / GsqrtDens(fp,EX) ) + abs(Vel(fp,EX))  )
      end do

      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP

        tmp2 = - alpha(fp) * ( GsqrtDDENS(fp,EX) - GsqrtDDENS(fp,IN) )
        del_flux(fp,DENS_VID,ke) = tmp1 * ( &
                        GsqrtDens(fp,EX) * Vel(fp,EX) - GsqrtDens(fp,IN) * Vel(fp,IN)  &
                      + tmp2  )

        tmp2 = - alpha(fp) * ( GsqrtDRHOT(fp,EX) - GsqrtDRHOT(fp,IN) )
        del_flux(fp,RHOT_VID,ke) = tmp1 * ( &
                        GsqrtRhot(fp,EX) * Vel(fp,EX) - GsqrtRhot(fp,IN) * Vel(fp,IN) &
                      + tmp2 )
      end do

      do fp=1, elem%NfpTot
        tmp3 = Gsqrt_(fp,EX) * dpres_(fp,EX)
        tmp4 = Gsqrt_(fp,IN) * dpres_(fp,IN)

        del_flux_tmp_mom(fp,1) = &
                (  tmp3 * RGsqrtV(fp,EX)               &
                 - tmp4 * RGsqrtV(fp,IN) ) * nz(fp,ke)
        
        del_flux_tmp_mom(fp,2) = &
              ( G1n_(fp,EX) + Gxz_(fp,EX) * nz(fp,ke) ) * tmp3 &
            - ( G1n_(fp,IN) + Gxz_(fp,IN) * nz(fp,ke) ) * tmp4
        
        del_flux_tmp_mom(fp,3) = &
              ( G2n_(fp,EX) + Gyz_(fp,EX) * nz(fp,ke) ) * tmp3 &
            - ( G2n_(fp,IN) + Gyz_(fp,IN) * nz(fp,ke) ) * tmp4
      end do
      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP

        tmp2 = - alpha(fp) * ( GsqrtMOMZ(fp,EX) - GsqrtMOMZ(fp,IN) )
        del_flux(fp,MOMZ_VID,ke) = tmp1 * ( &
                        GsqrtMOMZ(fp,EX) * Vel(fp,EX) - GsqrtMOMZ(fp,IN) * Vel(fp,IN) &
                      + del_flux_tmp_mom(fp,1) &
                      + tmp2 )
      end do
      do fp=1, elem%NfpTot
        tmp1 = lmesh%Fscale(fp,ke) * 0.5_RP

        tmp2 = - alpha(fp) * ( GsqrtMOMX(fp,EX) - GsqrtMOMX(fp,IN) )
        del_flux(fp,MOMX_VID,ke) = tmp1 * ( &
                        GsqrtMOMX(fp,EX) * Vel(fp,EX) - GsqrtMOMX(fp,IN) * Vel(fp,IN) &
                      + del_flux_tmp_mom(fp,2) &
                      + tmp2 )

        tmp2 = - alpha(fp) * ( GsqrtMOMY(fp,EX) - GsqrtMOMY(fp,IN) )
        del_flux(fp,MOMY_VID,ke) = tmp1 * ( &
                        GsqrtMOMY(fp,EX) * Vel(fp,EX) - GsqrtMOMY(fp,IN) * Vel(fp,IN) &
                      + del_flux_tmp_mom(fp,3) &
                      + tmp2 )
      end do
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc
  
end module scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux
