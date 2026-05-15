!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Nonhydrostatic model / HEVE / Numflux
!!
!! @par Description
!!      An module for numerical fluxes for HEVE atmospheric dynamical process which runs on GPU.
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_gpu
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
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_gpu

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
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_gpu( &
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
    
    integer :: ke, fp, i, iP, iM
    integer :: ke2D
    real(RP) :: Vel_M, Vel_P, alpha
    real(RP) :: DPRES_M, DPRES_P
    real(RP) :: GsqrtDens_M, GsqrtDens_P
    real(RP) :: GsqrtRhot_M, GsqrtRhot_P
    real(RP) :: GsqrtDDENS_M, GsqrtDDENS_P
    real(RP) :: GsqrtMOMX_M, GsqrtMOMX_P
    real(RP) :: GsqrtMOMY_M, GsqrtMOMY_P
    real(RP) :: GsqrtMOMZ_M, GsqrtMOMZ_P
    real(RP) :: GsqrtDRHOT_M, GsqrtDRHOT_P
    real(RP) :: Gsqrt_M, Gsqrt_P
    real(RP) :: GsqrtV_M, GsqrtV_P
    real(RP) :: RGsqrtV_M, RGsqrtV_P
    real(RP) :: G13_M, G13_P
    real(RP) :: G23_M, G23_P
    real(RP) :: Gnn_M, Gnn_P

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR
    
    real(RP) :: tmp1
    real(RP) :: del_flux_tmp_momz, del_flux_tmp_momx, del_flux_tmp_momy
    real(RP) :: GsqrtDPres_M, GsqrtDPres_P
    real(RP) :: nx_, ny_, nz_, fscale_

    integer :: NeS, NeE, NfpTot
    !------------------------------------------------------------------------

    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    NeS = lmesh%NeS; NeE = lmesh%NeE
    NfpTot = elem%NfpTot

    !$acc parallel present( DDENS_,MOMX_,MOMY_,MOMZ_,DRHOT_,DPRES, &
    !$acc                             DENS_hyd,PRES_hyd,THERM_hyd,           &
    !$acc                             Gsqrt,G13,G23,nx,ny,nz,vmapM,vmapP,    &
    !$acc                             del_flux, lmesh, elem )
    !$acc loop gang 
    do ke=NeS, NeE
      !$acc loop vector
      do fp=1, NfpTot

        iM = vmapM(fp,ke); iP = vmapP(fp,ke)
        ke2D = lmesh%EMap3Dto2D(ke)

        nx_ = nx(fp,ke); ny_ = ny(fp,ke); nz_ = nz(fp,ke)

        !-
        Gsqrt_M = Gsqrt(iM)
        Gsqrt_P = Gsqrt(iP)
        GsqrtV_M = Gsqrt_M
        GsqrtV_P = Gsqrt_P

        RGsqrtV_M = 1.0_RP / GsqrtV_M
        RGsqrtV_P = 1.0_RP / GsqrtV_P

        G13_M = G13(iM)
        G13_P = G13(iP)
        G23_M = G23(iM)
        G23_P = G23(iP)

        GsqrtDDENS_M = Gsqrt_M * DDENS_(iM)
        GsqrtDDENS_P = Gsqrt_P * DDENS_(iP)
        GsqrtMOMX_M = Gsqrt_M * MOMX_(iM)
        GsqrtMOMX_P = Gsqrt_P * MOMX_(iP)
        GsqrtMOMY_M = Gsqrt_M * MOMY_(iM)
        GsqrtMOMY_P = Gsqrt_P * MOMY_(iP)
        GsqrtMOMZ_M = Gsqrt_M * MOMZ_(iM)
        GsqrtMOMZ_P = Gsqrt_P * MOMZ_(iP)
        GsqrtDRHOT_M = Gsqrt_M * DRHOT_(iM)
        GsqrtDRHOT_P = Gsqrt_P * DRHOT_(iP)

        DPRES_M = DPRES(iM)
        DPRES_P = DPRES(iP)

        GsqrtDens_P = GsqrtDDENS_P + Gsqrt_P * DENS_hyd(iP)
        GsqrtDens_M = GsqrtDDENS_M + Gsqrt_M * DENS_hyd(iM)

        GsqrtRhot_P = Gsqrt_P * THERM_hyd(iP) + GsqrtDRHOT_P
        GsqrtRhot_M = Gsqrt_M * THERM_hyd(iM) + GsqrtDRHOT_M

        Vel_M = ( GsqrtMOMX_M * nx_ + GsqrtMOMY_M * ny_                     &
                  + ( ( GsqrtMOMZ_M * RGsqrtV_M                             &
                      + G13_M * GsqrtMOMX_M + G23_M * GsqrtMOMY_M ) * nz_ ) &
                  ) / GsqrtDens_M

        Vel_P = ( GsqrtMOMX_P * nx_ + GsqrtMOMY_P * ny_                     &
                  + ( ( GsqrtMOMZ_P * RGsqrtV_P                             &
                      + G13_P * GsqrtMOMX_P + G23_P * GsqrtMOMY_P ) * nz_ ) &
                  ) / GsqrtDens_P
 

        !--

        tmp1 = abs( nx_ ) + abs( ny_ )
        Gnn_M = tmp1 &
              + ( 1.0_RP * RGsqrtV_M**2 + G13_M**2 + G23_M**2 ) * abs( nz_ )        

        Gnn_P = tmp1 &
              + ( 1.0_RP * RGsqrtV_P**2 + G13_P**2 + G23_P**2 ) * abs( nz_ )

        alpha = max( sqrt( Gnn_M * gamm * ( PRES_hyd(iM) + DPRES_M ) * Gsqrt_M / GsqrtDens_M ) + abs(Vel_M), &
                     sqrt( Gnn_P * gamm * ( PRES_hyd(iP) + DPRES_P ) * Gsqrt_P / GsqrtDens_P ) + abs(Vel_P)  )


        !-

        !---------------------------------------
        fscale_ = 0.5_RP * lmesh%Fscale(fp,ke)

        !- density and potential temperature
      
        del_flux(fp,DENS_VID,ke) = fscale_ * ( &
                        GsqrtDens_P * Vel_P - GsqrtDens_M * Vel_M &
                      - alpha * ( GsqrtDDENS_P - GsqrtDDENS_M )   )

        del_flux(fp,RHOT_VID,ke) = fscale_ * ( &
                        GsqrtRhot_P * Vel_P - GsqrtRhot_M * Vel_M &
                      - alpha * ( GsqrtDRHOT_P - GsqrtDRHOT_M )   )

        !-
        GsqrtDPres_M = Gsqrt_M * DPRES_M
        GsqrtDPres_P = Gsqrt_P * DPRES_P

        del_flux(fp,MOMZ_VID,ke) = fscale_ * ( &
                        GsqrtMOMZ_P * Vel_P - GsqrtMOMZ_M * Vel_M &
                      + (  GsqrtDPres_P * RGsqrtV_P          &
                         - GsqrtDPres_M * RGsqrtV_M ) * nz_  &
                      - alpha * ( GsqrtMOMZ_P - GsqrtMOMZ_M ) )

        del_flux(fp,MOMX_VID,ke) = fscale_ * ( &
                        GsqrtMOMX_P * Vel_P - GsqrtMOMX_M * Vel_M &
                      + ( nx_ + G13_P * nz_ ) * GsqrtDPres_P   &
                      - ( nx_ + G13_M * nz_ ) * GsqrtDPres_M   &
                      - alpha * ( GsqrtMOMX_P - GsqrtMOMX_M )  )

        del_flux(fp,MOMY_VID,ke) = fscale_ * ( &
                        GsqrtMOMY_P * Vel_P - GsqrtMOMY_M * Vel_M &
                      + ( ny_ + G23_P * nz_ ) * GsqrtDPres_P   &
                      - ( ny_ + G23_M * nz_ ) * GsqrtDPres_M   &
                      - alpha * ( GsqrtMOMY_P - GsqrtMOMY_M ) )
      end do
    end do
    !$acc end parallel

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_gpu  
  
end module scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_gpu
