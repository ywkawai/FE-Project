!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Global nonhydrostatic model / HEVE
!!
!! @par Description
!!      HEVE DGM scheme for Global Atmospheric Dynamical process. 
!!      The governing equations is a fully compressible nonhydrostatic equations, 
!!      which consist of mass, momentum, and thermodynamics (density * potential temperature conservation) equations. 
!!
!! @author Yuta Kawai, Xuanzhengbo Ren, and Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_globalnonhydro3d_rhot_heve
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,    &
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  &
    CVdry => CONST_CVdry,  &
    PRES00 => CONST_PRE00, &
    RPlanet => CONST_RADIUS

  use scale_sparsemat  
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_element_operation_base, only: ElementOperationBase3D
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base2d, only: MeshBase2D
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
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_Init
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_Final
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_shallow_atm
  public :: atm_dyn_dgm_globalnh3d_rhot_heve_cal_tend_shallow_atm_asis
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_deep_atm

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
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_Init( mesh )
    implicit none
    class(MeshBase3D), intent(in) :: mesh
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Init( mesh )

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_Init

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_Final()
    implicit none
    !--------------------------------------------
    
    call atm_dyn_dgm_nonhydro3d_common_Final()    
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_Final  

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnh3d_rhot_heve_cal_tend_shallow_atm_asis( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                   & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, PRES_hyd_ref, & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,                                & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift,                             & ! (in)
    lmesh, elem, lmesh2D, elem2D )                                                   ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc_asis
    use scale_const, only: &
      OHM => CONST_OHM
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
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot (elem%Np,lmesh%NeA)    
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DPRES_hyd(elem%Np), GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: RHOT_(elem%Np)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%np), drho(elem%Np)

    real(RP) :: G11(elem%Np), G12(elem%Np), G22(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np)
    real(RP) :: X2D(elem2D%Np,lmesh2D%Ne), Y2D(elem2D%Np,lmesh2D%Ne)
    real(RP) :: X(elem%Np), Y(elem%Np), twoOVdel2(elem%Np)
    real(RP) :: CORI(elem%Np,2)
    logical :: is_panel1to4
    real(RP) :: s

    integer :: ke, ke2d
    integer :: p, p12, p3

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux, del_flux_hyd,                                                  & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd,         & ! (in)
      Rtot, CVtot, CPtot,                                                      & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%GsqrtH, lmesh%gam, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),             & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%vmapM, lmesh%vmapP, elem%IndexH2Dto3D_bnd,                         & ! (in)
      lmesh, elem, lmesh2D, elem2D                                             ) ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    s = 1.0_RP 
    is_panel1to4 = .true.
    if ( lmesh%panelID == 5 ) then
      is_panel1to4 = .false.
    else if ( lmesh%panelID == 6 ) then
      is_panel1to4 = .false.
      s = - 1.0_RP
    end if

    !$omp parallel private(                        &
    !$omp RHOT_, rdens_, u_, v_, w_, wt_,          &
    !$omp Fx, Fy, Fz, LiftDelFlx,                  &
    !$omp drho, DPRES_hyd, GradPhyd_x, GradPhyd_y, &
    !$omp G11, G12, G22, GsqrtV, RGsqrtV,          &
    !$omp X, Y, twoOVdel2,                         &
    !$omp CORI, ke, ke2D                           )

    !$omp do
    do ke2D = lmesh2D%NeS, lmesh2D%NeE
      X2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,1))
      Y2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,2))
    end do

    !$omp do
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      G11(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,1)
      G12(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,2)
      G22(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,2,2)
      GsqrtV(:)  = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2d)
      RGsqrtV(:) = 1.0_RP / GsqrtV(:)

      !--
      RHOT_(:) = P0ovR * ( PRES_hyd(:,ke) * rP0 )**rgamm + DRHOT_(:,ke)
      ! DPRES_(:) = PRES00 * ( Rtot(:,ke) * rP0 * RHOT_(:) )**( CPtot(:,ke) / CVtot(:,ke) ) &
      !           - PRES_hyd(:,ke)

      rdens_(:) = 1.0_RP / ( DDENS_(:,ke) + DENS_hyd(:,ke) )
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)
      wt_(:) = w_(:) * RGsqrtV(:) + lmesh%GI3(:,ke,1) * u_(:) + lmesh%GI3(:,ke,2) * v_(:) 

      X(:) = X2D(elem%IndexH2Dto3D,ke2d)
      Y(:) = Y2D(elem%IndexH2Dto3D,ke2d)
      twoOVdel2(:) = 2.0_RP / ( 1.0_RP + X(:)**2 + Y(:)**2 )

      CORI(:,1) = s * OHM * twoOVdel2(:) * ( - X(:) * Y(:)          * MOMX_(:,ke) + ( 1.0_RP + Y(:)**2 ) * MOMY_(:,ke) )
      CORI(:,2) = s * OHM * twoOVdel2(:) * ( - ( 1.0_RP + X(:)**2 ) * MOMX_(:,ke) +  X(:) * Y(:)         * MOMY_(:,ke) )
      if ( is_panel1to4 ) then
        CORI(:,1) = s * Y(:) * CORI(:,1)
        CORI(:,2) = s * Y(:) * CORI(:,2)
      end if

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
            lmesh%Escale(:,ke,1,1) * Fx(:)     &
          + lmesh%Escale(:,ke,2,2) * Fy(:)     &
          + lmesh%Escale(:,ke,3,3) * Fz(:)     &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      !-- MOMX
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMX_(:,ke) + G11(:) * DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMX_(:,ke) + G12(:) * DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMX_(:,ke)                              &
                                                    + ( lmesh%GI3(:,ke,1) * G11(:) + lmesh%GI3(:,ke,2) * G12(:) ) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                         &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                         &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                         &
            + LiftDelFlx(:)                   ) / lmesh%Gsqrt(:,ke)                  &
          - twoOVdel2(:) * Y(:) *                                                    &
            ( X(:) * Y(:) * u_(:) - (1.0_RP + Y(:)**2) * v_(:) ) * MOMX_(:,ke)       &
          - ( G11(:) * GradPhyd_x(:) + G12(:) * GradPhyd_y(:) ) * RGsqrtV(:)         &
          + CORI(:,1)              

      !-- MOMY
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMY_(:,ke) + G12(:) * DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMY_(:,ke) + G22(:) * DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMY_(:,ke)                              &
                                                    + ( lmesh%GI3(:,ke,1) * G12(:) + lmesh%GI3(:,ke,2) * G22(:) ) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = &
            - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                         &
              + lmesh%Escale(:,ke,2,2) * Fy(:)                                         &
              + lmesh%Escale(:,ke,3,3) * Fz(:)                                         &
              + LiftDelFlx(:)                  ) / lmesh%Gsqrt(:,ke)                   &
            - twoOVdel2(:) * X(:) *                                                    &
              ( - (1.0_RP + X(:)**2) * u_(:) + X(:) * Y(:) * v_(:) ) * MOMY_(:,ke)     &
            - ( G12(:) * GradPhyd_x(:) + G22(:) * GradPhyd_y(:) ) * RGsqrtV(:)         &
            + CORI(:,2) 

      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *   u_ (:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *   v_ (:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMZ_(:,ke) + RGsqrtV(:) * DPRES_(:,ke) ), Fz)
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
    !$omp end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_globalnh3d_rhot_heve_cal_tend_shallow_atm_asis

  !OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_shallow_atm( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                   & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, PRES_hyd_ref, & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,                                & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift,                             & ! (in)
    lmesh, elem, lmesh2D, elem2D )                                                   ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc
    use scale_const, only: &
      OHM => CONST_OHM
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
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot (elem%Np,lmesh%NeA)    
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    real(RP) :: Flux(elem%Np,3,5), DFlux(elem%Np,4,5)
    real(RP) :: del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP) :: rdens_, u_, v_, w_, pt_
    real(RP) :: drho(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np), RGsqrt(elem%Np)
    real(RP) :: Gsqrt_, GsqrtDPRES_, E11, E22, E33

    real(RP) :: G11(elem%Np), G12(elem%Np), G22(elem%Np)
    real(RP) :: X2D(elem2D%Np,lmesh2D%Ne), Y2D(elem2D%Np,lmesh2D%Ne)
    real(RP) :: X(elem%Np), Y(elem%Np), twoOVdel2(elem%Np)
    real(RP) :: CORI(elem%Np,2)
    logical :: is_panel1to4
    real(RP) :: s

    integer :: ke, ke2d
    integer :: p

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux,                                                                & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd,         & ! (in)
      Rtot, CVtot, CPtot,                                                      & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%GsqrtH, lmesh%gam, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),             & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%vmapM, lmesh%vmapP, elem%IndexH2Dto3D_bnd,                         & ! (in)
      lmesh, elem, lmesh2D, elem2D                                             ) ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    s = 1.0_RP 
    is_panel1to4 = .true.
    if ( lmesh%panelID == 5 ) then
      is_panel1to4 = .false.
    else if ( lmesh%panelID == 6 ) then
      is_panel1to4 = .false.
      s = - 1.0_RP
    end if

    !$omp parallel private( ke, ke2d, p,                  &
    !$omp rdens_, u_, v_, w_, pt_,                        &
    !$omp drho, CORI,                                     &
    !$omp RGsqrt, GsqrtV, RGsqrtV, Gsqrt_, GsqrtDPRES_,   &
    !$omp G11, G12, G22, X, Y, twoOVdel2,                 &
    !$omp E11, E22, E33, DFlux, Flux    )

    !$omp do
    do ke2D = lmesh2D%NeS, lmesh2D%NeE
      X2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,1))
      Y2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,2))
    end do

    !$omp do
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)

      do p=1, elem%Np
        G11(p) = lmesh%GIJ(elem%IndexH2Dto3D(p),ke2d,1,1)
        G12(p) = lmesh%GIJ(elem%IndexH2Dto3D(p),ke2d,1,2)
        G22(p) = lmesh%GIJ(elem%IndexH2Dto3D(p),ke2d,2,2)
      end do
      do p=1, elem%Np
        GsqrtV(p)  = lmesh%Gsqrt(p,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D(p),ke2d)
        RGsqrtV(p) = 1.0_RP / GsqrtV(p)
        RGsqrt(p) = 1.0_RP / lmesh%Gsqrt(p,ke)
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
        pt_ = ( P0ovR * (PRES_hyd(p,ke) * rP0)**rgamm + DRHOT_(p,ke) ) &
            / ( DDENS_(p,ke) + DENS_hyd(p,ke) )

        Flux(p,1,RHOT_VID) = Flux(p,1,DENS_VID) * pt_
        Flux(p,2,RHOT_VID) = Flux(p,2,DENS_VID) * pt_
        Flux(p,3,RHOT_VID) = Flux(p,3,DENS_VID) * pt_
      end do

      do p=1, elem%Np
        w_ = MOMZ_(p,ke) / (DDENS_(p,ke) + DENS_hyd(p,ke))

        Flux(p,1,MOMZ_VID) = Flux(p,1,DENS_VID) * w_
        Flux(p,2,MOMZ_VID) = Flux(p,2,DENS_VID) * w_
        Flux(p,3,MOMZ_VID) = Flux(p,3,DENS_VID) * w_ + lmesh%Gsqrt(p,ke) * RGsqrtV(p) * DPRES_(p,ke)
      end do

      do p=1, elem%Np
        GsqrtDPRES_ = lmesh%Gsqrt(p,ke) * DPRES_(p,ke)

        rdens_ = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))
        u_ = MOMX_(p,ke) * rdens_
        v_ = MOMY_(p,ke) * rdens_

        Flux(p,1,MOMX_VID) = Flux(p,1,DENS_VID) * u_ + G11(p) * GsqrtDPRES_
        Flux(p,2,MOMX_VID) = Flux(p,2,DENS_VID) * u_ + G12(p) * GsqrtDPRES_
        Flux(p,3,MOMX_VID) = Flux(p,3,DENS_VID) * u_ + GsqrtDPRES_ * ( G11(p) * lmesh%GI3(p,ke,1) + G12(p) * lmesh%GI3(p,ke,2) ) 

        Flux(p,1,MOMY_VID) = Flux(p,1,DENS_VID) * v_ + G12(p) * GsqrtDPRES_
        Flux(p,2,MOMY_VID) = Flux(p,2,DENS_VID) * v_ + G22(p) * GsqrtDPRES_ 
        Flux(p,3,MOMY_VID) = Flux(p,3,DENS_VID) * v_ + GsqrtDPRES_ * ( G12(p) * lmesh%GI3(p,ke,1) + G22(p) * lmesh%GI3(p,ke,2) ) 
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
        X(p) = X2D(elem%IndexH2Dto3D(p),ke2d)
        Y(p) = Y2D(elem%IndexH2Dto3D(p),ke2d)
        twoOVdel2(p) = 2.0_RP / ( 1.0_RP + X(p)**2 + Y(p)**2 )
      end do

      do p=1, elem%Np
        CORI(p,1) = s * OHM * twoOVdel2(p) * ( - X(p) * Y(p)          * MOMX_(p,ke) + ( 1.0_RP + Y(p)**2 ) * MOMY_(p,ke) )
        CORI(p,2) = s * OHM * twoOVdel2(p) * ( - ( 1.0_RP + X(p)**2 ) * MOMX_(p,ke) +  X(p) * Y(p)         * MOMY_(p,ke) )
      end do
      if ( is_panel1to4 ) then
        do p=1, elem%Np
          CORI(p,1) = s * Y(p) * CORI(p,1)
          CORI(p,2) = s * Y(p) * CORI(p,2)
        end do
      end if

      do p=1, elem%Np
        rdens_ = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))
        u_ = MOMX_(p,ke) * rdens_
        v_ = MOMY_(p,ke) * rdens_

        MOMX_dt(p,ke) = - ( G11(p) * DPhydDx(p,ke) + G12(p) * DPhydDy(p,ke) ) &
                        - twoOVdel2(p) * Y(p) *                                            &
                          ( X(p) * Y(p) * u_ - (1.0_RP + Y(p)**2) * v_ ) * MOMX_(p,ke)     &    
                        + CORI(p,1)
        
        MOMY_dt(p,ke) = - ( G12(p) * DPhydDx(p,ke) + G22(p) * DPhydDy(p,ke) ) &
                        - twoOVdel2(p) * X(p) *                                            &
                          ( - (1.0_RP + X(p)**2) * u_ + X(p) * Y(p) * v_ ) * MOMY_(p,ke)   &
                        + CORI(p,2)
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
    !$omp end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_shallow_atm

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_deep_atm( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                    & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, PRES_hyd_ref,  & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,                                 & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D ) ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc_asis
    use scale_const, only: &
      OHM => CONST_OHM
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
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot (elem%Np,lmesh%NeA)    
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DPRES_hyd(elem%Np), GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: RHOT_(elem%Np)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%np), drho(elem%Np)

    real(RP) :: G11(elem%Np), G12(elem%Np), G22(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np), Rgam2(elem%Np)
    real(RP) :: X2D(elem2D%Np,lmesh2D%Ne), Y2D(elem2D%Np,lmesh2D%Ne)
    real(RP) :: X(elem%Np), Y(elem%Np), twoOVdel2(elem%Np)
    real(RP) :: OM1(elem%Np), OM2(elem%Np), OM3(elem%Np), DEL(elem%Np), R(elem%Np)
    logical :: is_panel1to4
    real(RP) :: s

    integer :: ke, ke2d
    integer :: p
    
    real(RP) :: rgamm    
    real(RP) :: rP0
    real(RP) :: P0ovR
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux, del_flux_hyd,                                                  & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd,         & ! (in)
      Rtot, CVtot, CPtot,                                                      & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%GsqrtH, lmesh%gam, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),             & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%vmapM, lmesh%vmapP, elem%IndexH2Dto3D_bnd,                         & ! (in)
      lmesh, elem, lmesh2D, elem2D                                             ) ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    P0ovR = PRES00 / Rdry

    s = 2.0_RP * OHM
    is_panel1to4 = .true.
    if ( lmesh%panelID == 5 ) then
      is_panel1to4 = .false.
    else if ( lmesh%panelID == 6 ) then
      is_panel1to4 = .false.
      s = - s
    end if

    !$omp parallel private(                           &
    !$omp RHOT_, rdens_, u_, v_, w_, wt_,             &
    !$omp Fx, Fy, Fz, LiftDelFlx,                     &
    !$omp drho, DPRES_hyd, GradPhyd_x, GradPhyd_y,    &
    !$omp G11, G12, G22, Rgam2, GsqrtV, RGsqrtV,      &
    !$omp X, Y, twoOVdel2,                            &
    !$omp OM1, OM2, OM3, DEL, R, ke, ke2D             )

    !$omp do
    do ke2D = lmesh2D%NeS, lmesh2D%NeE
      X2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,1))
      Y2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,2))
    end do

    !$omp do
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      Rgam2(:) = 1.0_RP / lmesh%gam(:,ke)**2
      G11(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,1) * Rgam2(:)
      G12(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,2) * Rgam2(:)
      G22(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,2,2) * Rgam2(:)
      GsqrtV(:)  = lmesh%Gsqrt(:,ke) * Rgam2(:) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2d)
      RGsqrtV(:) = 1.0_RP / GsqrtV(:)

      !--
      RHOT_(:) = P0ovR * ( PRES_hyd(:,ke) * rP0 )**rgamm + DRHOT_(:,ke)
      ! DPRES_(:) = PRES00 * ( Rtot(:,ke) * rP0 * RHOT_(:) )**( CPtot(:,ke) / CVtot(:,ke) ) &
      !           - PRES_hyd(:,ke)

      rdens_(:) = 1.0_RP / ( DDENS_(:,ke) + DENS_hyd(:,ke) )
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)
      wt_(:) = w_(:) * RGsqrtV(:) + lmesh%GI3(:,ke,1) * u_(:) + lmesh%GI3(:,ke,2) * v_(:) 

      X(:) = X2D(elem%IndexH2Dto3D,ke2d)
      Y(:) = Y2D(elem%IndexH2Dto3D,ke2d)
      DEL(:) = sqrt( 1.0_RP + X(:)**2 + Y(:)**2 )
      twoOVdel2(:) = 2.0_RP / ( 1.0_RP + X(:)**2 + Y(:)**2 )

      R(:) = RPlanet * lmesh%gam(:,ke)

      !-  pnl=1~4: OM1:                     0, OM2: s del / (r (1+Y^2)) ,   s OM3 : Y /del
      !-  pnl=5,6: OM1: - s X del/(r (1+X^2)), OM2: -  s Y del/(r(1+Y^2)), OM3 : s/del
      if ( is_panel1to4 ) then
        OM1(:) = 0.0_RP
        OM2(:) = s * DEL(:) / ( R(:) * ( 1.0_RP + Y(:)**2 ) )
        OM3(:) = s * Y(:) / DEL(:)
      else
        OM1(:) = - s * X(:) * DEL(:) / ( R(:) * ( 1.0_RP + X(:)**2 ) )
        OM2(:) = - s * Y(:) * DEL(:) / ( R(:) * ( 1.0_RP + Y(:)**2 ) )
        OM3(:) =   s / DEL(:)
      end if

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
            lmesh%Escale(:,ke,1,1) * Fx(:)     &
          + lmesh%Escale(:,ke,2,2) * Fy(:)     &
          + lmesh%Escale(:,ke,3,3) * Fz(:)     &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      !-- MOMX
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMX_(:,ke) + G11(:) * DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMX_(:,ke) + G12(:) * DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMX_(:,ke)                              &
                                                    + ( lmesh%GI3(:,ke,1) * G11(:) + lmesh%GI3(:,ke,2) * G12(:) ) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                                &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                                &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                                &
            + LiftDelFlx(:)                   ) / lmesh%Gsqrt(:,ke)                         &
          - twoOVdel2(:) * Y(:) *                                                           & !-> metric terms
            ( X(:) * Y(:) * u_(:) - ( 1.0_RP + Y(:)**2 ) * v_(:) ) * MOMX_(:,ke)            & !
          - 2.0_RP * u_(:) * MOMZ_(:,ke) / R(:)                                             & !<-
          - ( G11(:) * GradPhyd_x(:) + G12(:) * GradPhyd_y(:) ) * RGsqrtV(:)                & !-> gradient hydrostatic pressure
          - lmesh%Gsqrt(:,ke) * (  G11(:) * ( OM2(:) * MOMZ_(:,ke) - OM3(:) * MOMY_(:,ke) ) & !-> Coriolis term
                                 - G12(:) * ( OM1(:) * MOMZ_(:,ke) - OM3(:) * MOMX_(:,ke) ) )

      !-- MOMY
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMY_(:,ke) + G12(:) * DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMY_(:,ke) + G22(:) * DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMY_(:,ke)                              &
                                                    + ( lmesh%GI3(:,ke,1) * G12(:) + lmesh%GI3(:,ke,2) * G22(:) ) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = &
            - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                                &
              + lmesh%Escale(:,ke,2,2) * Fy(:)                                                &
              + lmesh%Escale(:,ke,3,3) * Fz(:)                                                &
              + LiftDelFlx(:)                  ) / lmesh%Gsqrt(:,ke)                          &
            - twoOVdel2(:) * X(:) *                                                           & !-> metric terms
              ( - (1.0_RP + X(:)**2) * u_(:) + X(:) * Y(:) * v_(:) ) * MOMY_(:,ke)            & !
            - 2.0_RP * v_(:) * MOMZ_(:,ke) / R(:)                                             & !<-
            - ( G12(:) * GradPhyd_x(:) + G22(:) * GradPhyd_y(:) ) * RGsqrtV(:)                & !-> gradient hydrostatic pressure
            - lmesh%Gsqrt(:,ke) * (  G12(:) * ( OM2(:) * MOMZ_(:,ke) - OM3(:) * MOMY_(:,ke) ) & !-> Coriolis term
                                   - G22(:) * ( OM1(:) * MOMZ_(:,ke) - OM3(:) * MOMX_(:,ke) ) )

      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *   u_ (:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *   v_ (:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMZ_(:,ke) + RGsqrtV(:) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                                                    &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                                                    &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                                                    &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)                                                               &
          - 0.25_RP * R(:) * twoOVdel2(:)**2 * ( 1.0_RP * X(:)**2 ) * ( 1.0_RP * Y(:)**2 )                      & !-> metric terms
            * ( - ( 1.0_RP + X(:)**2 ) * MOMX_(:,ke) * u_(:)                                                    & !
                + 2.0_RP * X(:) * Y(:) * MOMX_(:,ke) * v_(:)                                                    & !
                - ( 1.0_RP + Y(:)**2 ) * MOMY_(:,ke) * v_(:)  )                                                 & !<-
          + 2.0_RP * DPRES_(:,ke) / R(:)                                                                        & !-> metric term with gradient of pressure deviaition
          - lmesh%Gsqrt(:,ke) * ( OM1(:) * MOMY_(:,ke) - OM2(:) * MOMX_(:,ke) )                                 & !-> Coriolis term
          - Grav * Rgam2(:) * drho(:)                                                                             !-> buoyancy term

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
    !$omp end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_deep_atm

end module scale_atm_dyn_dgm_globalnonhydro3d_rhot_heve
