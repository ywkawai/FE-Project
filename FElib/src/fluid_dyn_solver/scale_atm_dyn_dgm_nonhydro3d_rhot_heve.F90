!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Regional nonhydrostatic model / HEVE
!!
!! @par Description
!!      HEVE DGM scheme for Atmospheric dynamical process. 
!!      The governing equations is a fully compressible nonhydrostatic equations, 
!!      which consist of mass, momentum, and thermodynamics (density * potential temperature conservation) equations. 
!!
!! @author Yuta Kawai, Xuanzhengbo Ren, and Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_rhot_heve
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
  use scale_element_operation_base, only: ElementOperationBase3D
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    atm_dyn_dgm_nonhydro3d_common_Init,                       &
    atm_dyn_dgm_nonhydro3d_common_Final,                      &
    DENS_VID => PRGVAR_DDENS_ID, RHOT_VID => PRGVAR_DRHOT_ID, &
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID,                               &
    PRGVAR_NUM, IntrpMat_VPOrdM1


  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_Init
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_Final
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_asis
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_cco

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
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_Init( mesh )
    implicit none
    class(MeshBase3D), intent(in) :: mesh
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Init( mesh )
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_Init


  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_Final()
    implicit none
    !--------------------------------------------
    
    call atm_dyn_dgm_nonhydro3d_common_Final()
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_Final  

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_asis( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,         & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_,         & ! (in)
    DENS_hyd, PRES_hyd, PRES_hyd_ref, THERM_hyd,         & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,      & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift,   & ! (in)
    lmesh, elem, lmesh2D, elem2D )                         ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_asis

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    class(ElementOperationBase3D), intent(in) :: element3D_operation
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift
    real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DPRES_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd_ref(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: THERM_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DPRES_hyd(elem%Np), GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: RHOT_(elem%Np)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%Np)
    real(RP) :: Cori(elem%Np)
    real(RP) :: drho(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np)

    integer :: ke, ke2d

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR       
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux, del_flux_hyd,                                                 & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd,        & ! (in)
      Rtot, CVtot, CPtot,                                                     & ! (in)
      lmesh%Gsqrt, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),                        & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, lmesh2D, elem2D )                                            ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( ke, ke2d, Cori,     &
    !$omp RHOT_, rdens_, u_, v_, w_, wt_,          &
    !$omp drho, DPRES_hyd, GradPhyd_x, GradPhyd_y, &
    !$omp GsqrtV, RGsqrtV,                         &
    !$omp Fx, Fy, Fz, LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      Cori(:) = CORIOLIS(elem%IndexH2Dto3D(:),ke2d)

      GsqrtV(:)  = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2d)
      RGsqrtV(:) = 1.0_RP / GsqrtV(:)

      !--
      RHOT_(:) = P0ovR * (PRES_hyd(:,ke) * rP0)**rgamm + DRHOT_(:,ke)
      ! DPRES_(:) = PRES00 * ( Rtot(:,ke) * rP0 * RHOT_(:) )**( CPtot(:,ke) / CVtot(:,ke) ) &
      !           - PRES_hyd(:,ke)

      rdens_(:) = 1.0_RP / (DDENS_(:,ke) + DENS_hyd(:,ke))
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)
      wt_(:) = w_(:) * RGsqrtV(:) + lmesh%GI3(:,ke,1) * u_(:) + lmesh%GI3(:,ke,2) * v_(:) 

      drho(:) = matmul(IntrpMat_VPOrdM1, DDENS_(:,ke))

      !-- Gradient hydrostatic pressure
      
      DPRES_hyd(:) = PRES_hyd(:,ke) - PRES_hyd_ref(:,ke)

      call sparsemat_matmul(Dx, GsqrtV(:) * DPRES_hyd(:), Fx)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,1) * DPRES_hyd(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,1), LiftDelFlx)
      GradPhyd_x(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      call sparsemat_matmul(Dy, GsqrtV(:) * DPRES_hyd(:), Fy)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,2) * DPRES_hyd(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,2), LiftDelFlx)
      GradPhyd_y(:) = lmesh%Escale(:,ke,2,2) * Fy(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      !-- DENS
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * MOMY_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( DDENS_(:,ke) + DENS_hyd(:,ke) ) * wt_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DENS_VID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)    &
          + lmesh%Escale(:,ke,2,2) * Fy(:)    &
          + lmesh%Escale(:,ke,3,3) * Fz(:)    &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      !-- MOMX
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * (  u_(:) * MOMX_(:,ke) + DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *    v_(:) * MOMX_(:,ke)                 , Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMX_(:,ke)                     &
                                                    + lmesh%GI3(:,ke,1) * DPRES_(:,ke)    ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) &
          - GradPhyd_x(:) * RGsqrtV(:)            &
          + Cori(:) * MOMY_(:,ke)

      !-- MOMY
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *    u_(:) * MOMY_(:,ke)                 , Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * (  v_(:) * MOMY_(:,ke) + DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMY_(:,ke)                     &
                                                    + lmesh%GI3(:,ke,2) * DPRES_(:,ke)    ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) &
          - GradPhyd_y(:) * RGsqrtV(:)            &
          - Cori(:) * MOMX_(:,ke)

      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *   u_(:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *   v_(:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMZ_(:,ke)           &
                                                    + RGsqrtV(:) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)       &
            + lmesh%Escale(:,ke,2,2) * Fy(:)       &
            + lmesh%Escale(:,ke,3,3) * Fz(:)       &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)  &
          - Grav * drho(:)

      !-- RHOT
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * u_ (:) * RHOT_(:), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * v_ (:) * RHOT_(:), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * wt_(:) * RHOT_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,RHOT_VID), LiftDelFlx)
      
      RHOT_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 
    end do

    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_asis


!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,         & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_,         & ! (in)
    DENS_hyd, PRES_hyd, PRES_hyd_ref, THERM_hyd,         & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,      & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift,   & ! (in)
    lmesh, elem, lmesh2D, elem2D )                         ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    class(ElementOperationBase3D), intent(in) :: element3D_operation
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift
    real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DPRES_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd_ref(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: THERM_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    real(RP) :: Flux(elem%Np,3,5), DFlux(elem%Np,4,5)
    real(RP) :: del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP) :: u_, v_, w_, pt_
    real(RP) :: cor
    real(RP) :: drho(elem%Np)
    real(RP) :: RDENS_(elem%Np), GsqrtV(elem%Np), RGsqrtV(elem%Np), RGsqrt(elem%Np)
    real(RP) :: Gsqrt_, GsqrtDPRES_, E11, E22, E33

    integer :: ke, ke2d
    integer :: p

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR    
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux,                                                                   & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, THERM_hyd, & ! (in)
      Rtot, CVtot, CPtot,                                                         & ! (in)
      lmesh%Gsqrt, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),                            & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),     & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                                   & ! (in)
      lmesh, elem, lmesh2D, elem2D )                                                ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)

    call PROF_rapstart('cal_dyn_tend_interior', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( ke, ke2d, p, cor,                 &
    !$omp u_, v_, w_, pt_,                                       &
    !$omp drho,                                                  &
    !$omp RDENS_, RGsqrt, GsqrtV, RGsqrtV, Gsqrt_, GsqrtDPRES_,  &
    !$omp E11, E22, E33, DFlux, Flux    )
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      
      do p=1, elem%Np
        GsqrtV(p)  = lmesh%Gsqrt(p,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D(p),ke2d)
        RGsqrtV(p) = 1.0_RP / GsqrtV(p)
        RGsqrt(p) = 1.0_RP / lmesh%Gsqrt(p,ke)
        RDENS_(p) = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))
      end do

      !--

      do p=1, elem%Np
        Gsqrt_ = lmesh%Gsqrt(p,ke)
        Flux(p,1,DENS_VID) = Gsqrt_ * MOMX_(p,ke)
        Flux(p,2,DENS_VID) = Gsqrt_ * MOMY_(p,ke)
        Flux(p,3,DENS_VID) = Gsqrt_ * ( &
            MOMZ_(p,ke) * RGsqrtV(p)        &
          + lmesh%GI3(p,ke,1) * MOMX_(p,ke) &
          + lmesh%GI3(p,ke,2) * MOMY_(p,ke) )
      end do
      do p=1, elem%Np
        pt_ = ( THERM_hyd(p,ke) + DRHOT_(p,ke) ) * RDENS_(p)

        Flux(p,1,RHOT_VID) = Flux(p,1,DENS_VID) * pt_
        Flux(p,2,RHOT_VID) = Flux(p,2,DENS_VID) * pt_
        Flux(p,3,RHOT_VID) = Flux(p,3,DENS_VID) * pt_

        w_ = MOMZ_(p,ke) * RDENS_(p)
        Flux(p,1,MOMZ_VID) = Flux(p,1,DENS_VID) * w_
        Flux(p,2,MOMZ_VID) = Flux(p,2,DENS_VID) * w_
        Flux(p,3,MOMZ_VID) = Flux(p,3,DENS_VID) * w_ + lmesh%Gsqrt(p,ke) * DPRES_(p,ke) * RGsqrtV(p)
      end do
      do p=1, elem%Np
        GsqrtDPRES_ = lmesh%Gsqrt(p,ke) * DPRES_(p,ke)

        u_ = MOMX_(p,ke) * RDENS_(p)
        Flux(p,1,MOMX_VID) = Flux(p,1,DENS_VID) * u_ + GsqrtDPRES_
        Flux(p,2,MOMX_VID) = Flux(p,2,DENS_VID) * u_ 
        Flux(p,3,MOMX_VID) = Flux(p,3,DENS_VID) * u_ + GsqrtDPRES_ * lmesh%GI3(p,ke,1)

        v_ = MOMY_(p,ke) * RDENS_(p)
        Flux(p,1,MOMY_VID) = Flux(p,1,DENS_VID) * v_
        Flux(p,2,MOMY_VID) = Flux(p,2,DENS_VID) * v_ + GsqrtDPRES_ 
        Flux(p,3,MOMY_VID) = Flux(p,3,DENS_VID) * v_ + GsqrtDPRES_ * lmesh%GI3(p,ke,2) 
      end do

      call element3D_operation%Div_var5( &
        Flux, del_flux(:,:,ke), & ! (in)
        DFlux                   ) ! (out)

      do p=1, elem%Np
        E11 = lmesh%Escale(p,ke,1,1)
        E22 = lmesh%Escale(p,ke,2,2)
        E33 = lmesh%Escale(p,ke,3,3)

        DENS_dt(p,ke) = - ( &
              E11 * DFlux(p,1,DENS_VID)    &
            + E22 * DFlux(p,2,DENS_VID)    &
            + E33 * DFlux(p,3,DENS_VID)    &
            + DFlux(p,4,DENS_VID) ) * RGsqrt(p)
          
        RHOT_dt(p,ke) = - ( &
              E11 * DFlux(p,1,RHOT_VID)    &
            + E22 * DFlux(p,2,RHOT_VID)    &
            + E33 * DFlux(p,3,RHOT_VID)    &
            + DFlux(p,4,RHOT_VID) ) * RGsqrt(p)
      end do

      call element3D_operation%VFilterPM1( DDENS_(:,ke), & ! (in)
        drho ) ! (out)      

      do p=1, elem%Np
        E11 = lmesh%Escale(p,ke,1,1)
        E22 = lmesh%Escale(p,ke,2,2)
        E33 = lmesh%Escale(p,ke,3,3)
            
        MOMZ_dt(p,ke) = - ( &
              E11 * DFlux(p,1,MOMZ_VID)     &
            + E22 * DFlux(p,2,MOMZ_VID)     &
            + E33 * DFlux(p,3,MOMZ_VID)     &
            + DFlux(p,4,MOMZ_VID) ) * RGsqrt(p) &
            - Grav * drho(p)            
      end do

      !--
      do p=1, elem%Np
        cor = CORIOLIS(elem%IndexH2Dto3D(p),ke2d)
        MOMX_dt(p,ke) = - DPhydDx(p,ke) &
                        + cor * MOMY_(p,ke)
        MOMY_dt(p,ke) = - DPhydDy(p,ke) &
                        - cor * MOMX_(p,ke)
      end do

      do p=1, elem%Np
        E11 = lmesh%Escale(p,ke,1,1)
        E22 = lmesh%Escale(p,ke,2,2)
        E33 = lmesh%Escale(p,ke,3,3)
  
        MOMX_dt(p,ke) = MOMX_dt(p,ke) - ( &
              E11 * DFlux(p,1,MOMX_VID)    &
            + E22 * DFlux(p,2,MOMX_VID)    &
            + E33 * DFlux(p,3,MOMX_VID)    &
            + DFlux(p,4,MOMX_VID) ) * RGsqrt(p)
        
        MOMY_dt(p,ke) = MOMY_dt(p,ke) - ( &
              E11 * DFlux(p,1,MOMY_VID)    &
            + E22 * DFlux(p,2,MOMY_VID)    &
            + E33 * DFlux(p,3,MOMY_VID)    &
            + DFlux(p,4,MOMY_VID) ) * RGsqrt(p)
      end do
    end do

    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend


!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_cco( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,         & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_,         & ! (in)
    DENS_hyd, PRES_hyd, PRES_hyd_ref, THERM_hyd,         & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,      & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift,   & ! (in)
    lmesh, elem, lmesh2D, elem2D )                         ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    class(ElementOperationBase3D), intent(in) :: element3D_operation
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift
    real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DPRES_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd_ref(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: THERM_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    real(RP) :: Flux(elem%Np,3,5), DFlux(elem%Np,3,5)
    real(RP) :: u_, v_, w_, pt_
    real(RP) :: cor
    real(RP) :: drho(elem%Np)
    real(RP) :: RDENS_(elem%Np), GsqrtV(elem%Np), RGsqrtV(elem%Np), RGsqrt(elem%Np)
    real(RP) :: Gsqrt_, GsqrtDPRES_, E11, E22, E33

    integer :: ke, ke2d
    integer :: p

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR    
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_interior', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( ke, ke2d, p, cor,                 &
    !$omp u_, v_, w_, pt_,                                       &
    !$omp drho,                                                  &
    !$omp RDENS_, RGsqrt, GsqrtV, RGsqrtV, Gsqrt_, GsqrtDPRES_,  &
    !$omp E11, E22, E33, DFlux, Flux    )
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      
      do p=1, elem%Np
        GsqrtV(p)  = lmesh%Gsqrt(p,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D(p),ke2d)
        RGsqrtV(p) = 1.0_RP / GsqrtV(p)
        RGsqrt(p) = 1.0_RP / lmesh%Gsqrt(p,ke)
        RDENS_(p) = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))
      end do

      !--

      do p=1, elem%Np
        Gsqrt_ = lmesh%Gsqrt(p,ke)
        Flux(p,1,DENS_VID) = Gsqrt_ * MOMX_(p,ke)
        Flux(p,2,DENS_VID) = Gsqrt_ * MOMY_(p,ke)
        Flux(p,3,DENS_VID) = Gsqrt_ * ( &
            MOMZ_(p,ke) * RGsqrtV(p)        &
          + lmesh%GI3(p,ke,1) * MOMX_(p,ke) &
          + lmesh%GI3(p,ke,2) * MOMY_(p,ke) )
      end do
      do p=1, elem%Np
        pt_ = ( THERM_hyd(p,ke) + DRHOT_(p,ke) ) * RDENS_(p)

        Flux(p,1,RHOT_VID) = Flux(p,1,DENS_VID) * pt_
        Flux(p,2,RHOT_VID) = Flux(p,2,DENS_VID) * pt_
        Flux(p,3,RHOT_VID) = Flux(p,3,DENS_VID) * pt_

        w_ = MOMZ_(p,ke) * RDENS_(p)
        Flux(p,1,MOMZ_VID) = Flux(p,1,DENS_VID) * w_
        Flux(p,2,MOMZ_VID) = Flux(p,2,DENS_VID) * w_
        Flux(p,3,MOMZ_VID) = Flux(p,3,DENS_VID) * w_ + lmesh%Gsqrt(p,ke) * DPRES_(p,ke) * RGsqrtV(p)
      end do
      do p=1, elem%Np
        GsqrtDPRES_ = lmesh%Gsqrt(p,ke) * DPRES_(p,ke)

        u_ = MOMX_(p,ke) * RDENS_(p)
        Flux(p,1,MOMX_VID) = Flux(p,1,DENS_VID) * u_ + GsqrtDPRES_
        Flux(p,2,MOMX_VID) = Flux(p,2,DENS_VID) * u_ 
        Flux(p,3,MOMX_VID) = Flux(p,3,DENS_VID) * u_ + GsqrtDPRES_ * lmesh%GI3(p,ke,1)

        v_ = MOMY_(p,ke) * RDENS_(p)
        Flux(p,1,MOMY_VID) = Flux(p,1,DENS_VID) * v_
        Flux(p,2,MOMY_VID) = Flux(p,2,DENS_VID) * v_ + GsqrtDPRES_ 
        Flux(p,3,MOMY_VID) = Flux(p,3,DENS_VID) * v_ + GsqrtDPRES_ * lmesh%GI3(p,ke,2) 
      end do

      call element3D_operation%Div_var5_2( &
        Flux,                   & ! (in)
        DFlux                   ) ! (out)

      do p=1, elem%Np
        E11 = lmesh%Escale(p,ke,1,1)
        E22 = lmesh%Escale(p,ke,2,2)
        E33 = lmesh%Escale(p,ke,3,3)

        DENS_dt(p,ke) = - ( &
              E11 * DFlux(p,1,DENS_VID)    &
            + E22 * DFlux(p,2,DENS_VID)    &
            + E33 * DFlux(p,3,DENS_VID)    &
            ) * RGsqrt(p)
          
        RHOT_dt(p,ke) = - ( &
              E11 * DFlux(p,1,RHOT_VID)    &
            + E22 * DFlux(p,2,RHOT_VID)    &
            + E33 * DFlux(p,3,RHOT_VID)    &
            ) * RGsqrt(p)
      end do

      call element3D_operation%VFilterPM1( DDENS_(:,ke), & ! (in)
        drho ) ! (out)      

      do p=1, elem%Np
        E11 = lmesh%Escale(p,ke,1,1)
        E22 = lmesh%Escale(p,ke,2,2)
        E33 = lmesh%Escale(p,ke,3,3)
            
        MOMZ_dt(p,ke) = - ( &
              E11 * DFlux(p,1,MOMZ_VID)     &
            + E22 * DFlux(p,2,MOMZ_VID)     &
            + E33 * DFlux(p,3,MOMZ_VID)     &
            ) * RGsqrt(p) &
            - Grav * drho(p)            
      end do

      !--
      do p=1, elem%Np
        cor = CORIOLIS(elem%IndexH2Dto3D(p),ke2d)
        MOMX_dt(p,ke) = - DPhydDx(p,ke) &
                        + cor * MOMY_(p,ke)
        MOMY_dt(p,ke) = - DPhydDy(p,ke) &
                        - cor * MOMX_(p,ke)
      end do

      do p=1, elem%Np
        E11 = lmesh%Escale(p,ke,1,1)
        E22 = lmesh%Escale(p,ke,2,2)
        E33 = lmesh%Escale(p,ke,3,3)
  
        MOMX_dt(p,ke) = MOMX_dt(p,ke) - ( &
              E11 * DFlux(p,1,MOMX_VID)    &
            + E22 * DFlux(p,2,MOMX_VID)    &
            + E33 * DFlux(p,3,MOMX_VID)    &
            ) * RGsqrt(p)
        
        MOMY_dt(p,ke) = MOMY_dt(p,ke) - ( &
              E11 * DFlux(p,1,MOMY_VID)    &
            + E22 * DFlux(p,2,MOMY_VID)    &
            + E33 * DFlux(p,3,MOMY_VID)    &
            ) * RGsqrt(p)
      end do
    end do

    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_cco
  
end module scale_atm_dyn_dgm_nonhydro3d_rhot_heve
