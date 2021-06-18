!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics HEVE
!!
!! @par Description
!!      HEVE DGM scheme for Global Atmospheric Dynamical process. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_globalnonhydro3d_heve
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
  public :: atm_dyn_dgm_globalnonhydro3d_heve_Init
  public :: atm_dyn_dgm_globalnonhydro3d_heve_Final
  public :: atm_dyn_dgm_globalnonhydro3d_heve_cal_tend

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  integer, private, parameter :: VARS_DDENS_ID  = 1
  integer, private, parameter :: VARS_MOMX_ID   = 2
  integer, private, parameter :: VARS_MOMY_ID   = 3
  integer, private, parameter :: VARS_MOMZ_ID   = 4
  integer, private, parameter :: VARS_DRHOT_ID  = 5
  integer, private, parameter :: PROG_VARS_NUM  = 5
  
  real(RP), private, allocatable :: IntrpMat_VPOrdM1(:,:)
  integer, private, allocatable :: iM2Dto3D(:)


  private :: cal_del_flux_dyn

contains
  subroutine atm_dyn_dgm_globalnonhydro3d_heve_Init( mesh )

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
  end subroutine atm_dyn_dgm_globalnonhydro3d_heve_Init


  subroutine atm_dyn_dgm_globalnonhydro3d_heve_Final()
    implicit none
    !--------------------------------------------
    
    deallocate( IntrpMat_VPOrdM1 )
    
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_heve_Final  

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_heve_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, CORIOLIS,          & ! (in)
    SL_flag, wdamp_tau, wdamp_height, hveldamp_flag,                            & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )

    use scale_atm_dyn_dgm_spongelayer, only: &
      atm_dyn_dgm_spongelayer_add_tend
    use scale_const, only: &
      OHM => CONST_OHM
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
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
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    logical, intent(in) :: SL_flag
    real(RP), intent(in) :: wdamp_tau
    real(RP), intent(in) :: wdamp_height
    logical, intent(in) :: hveldamp_flag

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: DPRES_(elem%Np)
    real(RP) :: RHOT_(elem%Np)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%np), drho(elem%Np)

    real(RP) :: G11(elem%Np), G12(elem%Np), G22(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np)
    real(RP) :: X2D(elem%Np,lmesh2D%Ne), Y2D(elem%Np,lmesh2D%Ne)
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
    call cal_del_flux_dyn( del_flux, del_flux_hyd,                             & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                 & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%GsqrtH, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),                        & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem, lmesh2D, elem2D                   ) ! (in)
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
    !$omp RHOT_, DPRES_, rdens_, u_, v_, w_, wt_,  &
    !$omp Fx, Fy, Fz, LiftDelFlx,                  &
    !$omp drho, GradPhyd_x, GradPhyd_y,            &
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
      DPRES_(:) = PRES00 * ( RovP0 * RHOT_(:) )**gamm &
                - PRES_hyd(:,ke)

      rdens_(:) = 1.0_RP / ( DDENS_(:,ke) + DENS_hyd(:,ke) )
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)
      wt_(:) = w_(:) * RGsqrtV(:) + lmesh%GI3(:,ke,1) * u_(:) + lmesh%GI3(:,ke,2) * v_(:) 

      X(:) = X2D(elem%IndexH2Dto3D,ke2d)
      Y(:) = Y2D(elem%IndexH2Dto3D,ke2d)
      twoOVdel2(:) = 2.0_RP / ( 1.0_RP + X(:)**2 + Y(:)**2 )

      CORI(:,1) = s * OHM * twoOVdel2(:) * ( - X(:) * Y(:)        * MOMX_(:,ke) + (1.0_RP + Y(:)**2) * MOMY_(:,ke) )
      CORI(:,2) = s * OHM * twoOVdel2(:) * ( - (1.0_RP + X(:)**2) * MOMX_(:,ke) +  X(:) * Y(:)       * MOMY_(:,ke) )
      if ( is_panel1to4 ) then
        CORI(:,1) = s * Y(:) * CORI(:,1)
        CORI(:,2) = s * Y(:) * CORI(:,2)
      end if

      drho(:) = matmul(IntrpMat_VPOrdM1, DDENS_(:,ke))

      !-- Gradient hydrostatic pressure
      
      call sparsemat_matmul(Dx, GsqrtV(:) * PRES_hyd(:,ke), Fx)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,1) * PRES_hyd(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,1), LiftDelFlx)
      GradPhyd_x(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      call sparsemat_matmul(Dy, GsqrtV(:) * PRES_hyd(:,ke), Fy)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,2) * PRES_hyd(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,2), LiftDelFlx)
      GradPhyd_y(:) = lmesh%Escale(:,ke,2,2) * Fy(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)
      
      !-- DENS
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * MOMY_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( DDENS_(:,ke) + DENS_hyd(:,ke) ) * wt_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_DDENS_ID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)     &
          + lmesh%Escale(:,ke,2,2) * Fy(:)     &
          + lmesh%Escale(:,ke,3,3) * Fz(:)     &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      !-- MOMX
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMX_(:,ke) + G11(:) * DPRES_(:) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMX_(:,ke) + G12(:) * DPRES_(:) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMX_(:,ke) +                             &
                                                     ( lmesh%GI3(:,ke,1) * G11(:) + lmesh%GI3(:,ke,2) * G12(:) ) * DPRES_(:) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_MOMX_ID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                         &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                         &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                         &
            + LiftDelFlx(:)                   ) / lmesh%Gsqrt(:,ke)                  &
          + twoOVdel2(:) * Y(:) *                                                    &
            ( - X(:) * Y(:) * u_(:) + (1.0_RP + Y(:)**2) * v_(:) ) * MOMX_(:,ke)     &
          - ( G11(:) * GradPhyd_x(:) + G12(:) * GradPhyd_y(:) ) * RGsqrtV(:)         &
          + CORI(:,1)              

      !-- MOMY
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMY_(:,ke) + G12(:) * DPRES_(:) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMY_(:,ke) + G22(:) * DPRES_(:) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMY_(:,ke) +                             &
                                                     ( lmesh%GI3(:,ke,1) * G12(:) + lmesh%GI3(:,ke,2) * G22(:) ) * DPRES_(:) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_MOMY_ID), LiftDelFlx)

      MOMY_dt(:,ke) = &
            - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                         &
              + lmesh%Escale(:,ke,2,2) * Fy(:)                                         &
              + lmesh%Escale(:,ke,3,3) * Fz(:)                                         &
              + LiftDelFlx(:)                  ) / lmesh%Gsqrt(:,ke)                   &
            + twoOVdel2(:) * X(:) *                                                    &
              ( (1.0_RP + X(:)**2) * u_(:) - X(:) * Y(:) * v_(:) ) * MOMY_(:,ke)       &
            - ( G12(:) * GradPhyd_x(:) + G22(:) * GradPhyd_y(:) ) * RGsqrtV(:)         &
            + CORI(:,2) 

      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *   u_ (:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *   v_ (:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMZ_(:,ke) + RGsqrtV(:) * DPRES_(:) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_MOMZ_ID), LiftDelFlx)
      
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
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_DRHOT_ID), LiftDelFlx)
      
      RHOT_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 

    end do
    !$omp end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_interior', 3)

    !- Sponge layer
    if (SL_flag) then
      call PROF_rapstart('cal_dyn_tend_sponge', 3)
      call atm_dyn_dgm_spongelayer_add_tend( &
        MOMX_dt, MOMY_dt, MOMZ_dt,                    & ! (out)
        MOMX_, MOMY_, MOMZ_, wdamp_tau, wdamp_height, & ! (in)
        hveldamp_flag, lmesh, elem                    ) ! (in)
      call PROF_rapend('cal_dyn_tend_sponge', 3)
    end if

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_heve_cal_tend

  !------

!OCL SERIAL
  subroutine cal_del_flux_dyn( del_flux, del_flux_hyd,             &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,       &
    Gsqrt, G11, G12, G22, GsqrtH, G13, G23, nx, ny, nz,            &
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D                     )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP), intent(out) ::  del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G11(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G12(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G22(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  GsqrtH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    
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
    !$omp Phyd_M, Phyd_P,                                                       &
    !$omp Gsqrt_P, Gsqrt_M, GsqrtV_P, GsqrtV_M, G13_P, G13_M, G23_P, G23_M,             &
    !$omp Gxz_P, Gxz_M, Gyz_P, Gyz_M, G1n_M, G2n_M, Gnn_P, Gnn_M                         )
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

      Gxz_M(:) = G11(iM2Dto3D(:),ke2D) * G13_M(:) + G12(iM2Dto3D(:),ke2D) * G23_M(:)
      Gxz_P(:) = G11(iM2Dto3D(:),ke2D) * G13_P(:) + G12(iM2Dto3D(:),ke2D) * G23_P(:)

      Gyz_M(:) = G12(iM2Dto3D(:),ke2D) * G13_M(:) + G22(iM2Dto3D(:),ke2D) * G23_M(:)
      Gyz_P(:) = G12(iM2Dto3D(:),ke2D) * G13_P(:) + G22(iM2Dto3D(:),ke2D) * G23_P(:)
               
      G1n_M(:)  = G11(iM2Dto3D(:),ke2D) * nx(:,ke) + G12(iM2Dto3D(:),ke2D) * ny(:,ke)
      G2n_M(:)  = G12(iM2Dto3D(:),ke2D) * nx(:,ke) + G22(iM2Dto3D(:),ke2D) * ny(:,ke)

      Gnn_M(:)  = G11(iM2Dto3D(:),ke2D) * abs( nx(:,ke) ) + G22(iM2Dto3D(:),ke2D) * abs( ny(:,ke) )    &
                + ( 1.0_RP / GsqrtV_M(:) + G13_M(:) * Gxz_M(:) + G23_M(:) * Gyz_M(:) ) * abs( nz(:,ke) )
      Gnn_P(:)  = G11(iM2Dto3D(:),ke2D) * abs( nx(:,ke) ) + G22(iM2Dto3D(:),ke2D) * abs( ny(:,ke) )    &
                + ( 1.0_RP / GsqrtV_P(:) + G13_P(:) * Gxz_P(:) + G23_P(:) * Gyz_P(:) ) * abs( nz(:,ke) )

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

      dpresM(:) = PRES00 * ( RovP0 * GsqrtRhotM(:) / Gsqrt_M(:) )**gamm &
                - Phyd_M(:)
      dpresP(:) = PRES00 * ( RovP0 * GsqrtRhotP(:) / Gsqrt_P(:) )**gamm &
                - Phyd_P(:)

      alpha(:) = max( sqrt( Gnn_M(:) * gamm * ( Phyd_M(:) + dpresM(:) ) * Gsqrt_M(:) / GsqrtDensM(:) ) + abs(VelM(:)), &
                      sqrt( Gnn_P(:) * gamm * ( Phyd_P(:) + dpresP(:) ) * Gsqrt_P(:) / GsqrtDensP(:) ) + abs(VelP(:))  )
      
      del_flux(:,ke,VARS_DDENS_ID) = 0.5_RP * ( &
                    ( GsqrtDensP(:) * VelP(:) - GsqrtDensM(:) * VelM(:) )  &
                    - alpha(:) * ( GsqrtDDENS_P(:) - GsqrtDDENS_M(:) )     )

      del_flux(:,ke,VARS_MOMX_ID ) = 0.5_RP * ( &
                    ( GsqrtMOMX_P(:) * VelP(:) - GsqrtMOMX_M(:) * VelM(:) ) &
                    + (  ( G1n_M(:) + Gxz_P(:) * nz(:,ke)) * dpresP(:)      &
                       - ( G1n_M(:) + Gxz_M(:) * nz(:,ke)) * dpresM(:) )    &
                    - alpha(:) * ( GsqrtMOMX_P(:) - GsqrtMOMX_M(:) )        )

      del_flux(:,ke,VARS_MOMY_ID ) = 0.5_RP * ( &
                    ( GsqrtMOMY_P(:) * VelP(:) - GsqrtMOMY_M(:) * VelM(:) ) &
                    + (  ( G2n_M(:) + Gyz_P(:) * nz(:,ke) ) * dpresP(:)     &
                       - ( G2n_M(:) + Gyz_M(:) * nz(:,ke) ) * dpresM(:) )   &
                    - alpha(:) * ( GsqrtMOMY_P(:) - GsqrtMOMY_M(:) )        )

      del_flux(:,ke,VARS_MOMZ_ID ) = 0.5_RP * ( &
                    ( GsqrtMOMZ_P(:) * VelP(:) - GsqrtMOMZ_M(:) * VelM(:) ) &
                    + (  dpresP(:) / GsqrtV_P(:)                            &
                       - dpresM(:) / GsqrtV_M(:) ) * nz(:,ke)               &
                    - alpha(:) * ( GsqrtMOMZ_P(:) - GsqrtMOMZ_M(:) )        )
                    
      del_flux(:,ke,VARS_DRHOT_ID) = 0.5_RP * ( &
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
  end subroutine cal_del_flux_dyn

end module scale_atm_dyn_dgm_globalnonhydro3d_heve
