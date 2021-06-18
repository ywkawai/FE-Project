!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics HEVI
!!
!! @par Description
!!      HEVI DGM scheme for Global Atmospheric Dynamical process. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_globalnonhydro3d_hevi
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
  use scale_element_modalfilter, only: ModalFilter
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use, intrinsic :: ieee_arithmetic
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_globalnonhydro3d_hevi_Init
  public :: atm_dyn_dgm_globalnonhydro3d_hevi_Final
  public :: atm_dyn_dgm_globalnonhydro3d_hevi_cal_tend
  public :: atm_dyn_dgm_globalnonhydro3d_hevi_cal_vi

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  integer, private, parameter :: DDENS_VID  = 1
  integer, private, parameter :: MOMX_VID   = 2
  integer, private, parameter :: MOMY_VID   = 3
  integer, private, parameter :: MOMZ_VID   = 4
  integer, private, parameter :: DRHOT_VID  = 5
  integer, private, parameter :: PROG_VARS_NUM  = 5
  
  real(RP), private, allocatable :: IntrpMat_VPOrdM1(:,:)
  integer, private, allocatable :: iM2Dto3D(:)


  private :: cal_del_flux_dyn

contains
  subroutine atm_dyn_dgm_globalnonhydro3d_hevi_Init( mesh )

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
  end subroutine atm_dyn_dgm_globalnonhydro3d_hevi_Init


  subroutine atm_dyn_dgm_globalnonhydro3d_hevi_Final()
    implicit none
    !--------------------------------------------
    
    deallocate( IntrpMat_VPOrdM1 )
    
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_hevi_Final  

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_hevi_cal_tend( &
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
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%Np)

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
    !$omp GradPhyd_x, GradPhyd_y,                  &
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
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( DDENS_(:,ke) + DENS_hyd(:,ke) )  &
                                                  * ( wt_(:) - w_(:) * RGsqrtV(:) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DDENS_VID), LiftDelFlx)

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
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                     &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                     &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                     &
            + LiftDelFlx(:)                   ) / lmesh%Gsqrt(:,ke)              &
          + twoOVdel2(:) * Y(:) *                                                &
            ( - X(:) * Y(:) * u_(:) + (1.0_RP + Y(:)**2) * v_(:) ) * MOMX_(:,ke) &
          - ( G11(:) * GradPhyd_x(:) + G12(:) * GradPhyd_y(:) ) * RGsqrtV(:)     &
          + CORI(:,1)              
      
      !-- MOMY
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMY_(:,ke) + G12(:) * DPRES_(:) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMY_(:,ke) + G22(:) * DPRES_(:) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMY_(:,ke) +                             &
                                                     ( lmesh%GI3(:,ke,1) * G12(:) + lmesh%GI3(:,ke,2) * G22(:) ) * DPRES_(:) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = &
            - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                   &
              + lmesh%Escale(:,ke,2,2) * Fy(:)                                   &
              + lmesh%Escale(:,ke,3,3) * Fz(:)                                   &
              + LiftDelFlx(:)                  ) / lmesh%Gsqrt(:,ke)             &
            + twoOVdel2(:) * X(:) *                                              &
              ( (1.0_RP + X(:)**2) * u_(:) - X(:) * Y(:) * v_(:) ) * MOMY_(:,ke) &
            - ( G12(:) * GradPhyd_x(:) + G22(:) * GradPhyd_y(:) ) * RGsqrtV(:)   &
            + CORI(:,2) 
      
      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * u_ (:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * v_ (:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * wt_(:) * MOMZ_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)       &
            + lmesh%Escale(:,ke,2,2) * Fy(:)       &
            + lmesh%Escale(:,ke,3,3) * Fz(:)       &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)  

      !-- RHOT
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * u_(:) * RHOT_(:), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * v_(:) * RHOT_(:), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) - w_(:) * RGsqrtV(:) ) * RHOT_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DRHOT_VID), LiftDelFlx)
      
      RHOT_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
    end do
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
  end subroutine atm_dyn_dgm_globalnonhydro3d_hevi_cal_tend

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
    real(RP) :: VelP(elem%NfpTot), VelM(elem%NfpTot)
    real(RP) :: VelhP(elem%NfpTot), VelhM(elem%NfpTot)
    real(RP) :: alpha(elem%NfpTot)
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
    real(RP) :: swV(elem%NfpTot)

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
    !$omp alpha, VelM, VelP, VelhM, VelhP,                                              &
    !$omp dpresM, dpresP, GsqrtDensM, GsqrtDensP, GsqrtRhotM, GsqrtRhotP,               &
    !$omp GsqrtMOMX_M, GsqrtMOMX_P, GsqrtMOMY_M, GsqrtMOMY_P, GsqrtMOMZ_M, GsqrtMOMZ_P, &
    !$omp GsqrtDDENS_M, GsqrtDDENS_P, GsqrtDRHOT_M, GsqrtDRHOT_P,                       &
    !$omp Phyd_M, Phyd_P,                                                               &
    !$omp Gsqrt_P, Gsqrt_M, GsqrtV_P, GsqrtV_M, G13_P, G13_M, G23_P, G23_M,             &
    !$omp Gxz_P, Gxz_M, Gyz_P, Gyz_M, G1n_M, G2n_M, Gnn_P, Gnn_M,                       &
    !$omp swV )
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
      swV(:) = 1.0_RP - nz(:,ke)**2

      Gxz_M(:) = G11(iM2Dto3D(:),ke2D) * G13_M(:) + G12(iM2Dto3D(:),ke2D) * G23_M(:)
      Gxz_P(:) = G11(iM2Dto3D(:),ke2D) * G13_P(:) + G12(iM2Dto3D(:),ke2D) * G23_P(:)

      Gyz_M(:) = G12(iM2Dto3D(:),ke2D) * G13_M(:) + G22(iM2Dto3D(:),ke2D) * G23_M(:)
      Gyz_P(:) = G12(iM2Dto3D(:),ke2D) * G13_P(:) + G22(iM2Dto3D(:),ke2D) * G23_P(:)
               
      G1n_M(:)  = G11(iM2Dto3D(:),ke2D) * nx(:,ke) + G12(iM2Dto3D(:),ke2D) * ny(:,ke)
      G2n_M(:)  = G12(iM2Dto3D(:),ke2D) * nx(:,ke) + G22(iM2Dto3D(:),ke2D) * ny(:,ke)

      Gnn_M(:)  = G11(iM2Dto3D(:),ke2D) * abs( nx(:,ke) ) + G22(iM2Dto3D(:),ke2D) * abs( ny(:,ke) )
      Gnn_P(:)  = G11(iM2Dto3D(:),ke2D) * abs( nx(:,ke) ) + G22(iM2Dto3D(:),ke2D) * abs( ny(:,ke) )

      GsqrtDensM(:) = GsqrtDDENS_M(:) + Gsqrt_M(:) * DENS_hyd(iM)
      GsqrtDensP(:) = GsqrtDDENS_P(:) + Gsqrt_P(:) * DENS_hyd(iP)

      GsqrtRhotM(:) = Gsqrt_M(:) * P0ovR * (Phyd_M(:) * rP0)**rgamm + GsqrtDRHOT_M(:)
      GsqrtRhotP(:) = Gsqrt_P(:) * P0ovR * (Phyd_P(:) * rP0)**rgamm + GsqrtDRHOT_P(:)

      VelhM(:) = ( GsqrtMOMX_M(:) * nx(:,ke) + GsqrtMOMY_M(:) * ny(:,ke)                    &
                 + ( G13_M(:) * GsqrtMOMX_M(:) + G23_M(:) * GsqrtMOMY_M(:) ) * nz(:,ke)     &
                 ) / GsqrtDensM(:)
      
      VelhP(:) = ( GsqrtMOMX_M(:) * nx(:,ke) + GsqrtMOMY_M(:) * ny(:,ke)                    &
                 + ( G13_M(:) * GsqrtMOMX_M(:) + G23_M(:) * GsqrtMOMY_M(:) ) * nz(:,ke)     &
                 ) / GsqrtDensM(:)

      VelM(:) = VelhM(:) + GsqrtMOMZ_M(:) / ( GsqrtV_M(:) * GsqrtDensM(:) ) * nz(:,ke)
      VelP(:) = VelhP(:) + GsqrtMOMZ_P(:) / ( GsqrtV_P(:) * GsqrtDensP(:) ) * nz(:,ke)

      dpresM(:) = PRES00 * ( RovP0 * GsqrtRhotM(:) / Gsqrt_M(:) )**gamm &
                - Phyd_M(:)
      dpresP(:) = PRES00 * ( RovP0 * GsqrtRhotP(:) / Gsqrt_P(:) )**gamm &
                - Phyd_P(:)

      alpha(:) = swV(:) * max( sqrt( Gnn_M(:) * gamm * ( Phyd_M(:) + dpresM(:) ) * Gsqrt_M(:) / GsqrtDensM(:) ) + abs(VelM(:)), &
                               sqrt( Gnn_P(:) * gamm * ( Phyd_P(:) + dpresP(:) ) * Gsqrt_P(:) / GsqrtDensP(:) ) + abs(VelP(:))  )
      
      del_flux(:,ke,DDENS_VID) = 0.5_RP * ( &
                    ( GsqrtDensP(:) * VelhP(:) - GsqrtDensM(:) * VelhM(:) ) &
                    - alpha(:) * ( GsqrtDDENS_P(:) - GsqrtDDENS_M(:) )      )

      del_flux(:,ke,MOMX_VID ) = 0.5_RP * ( &
                    ( GsqrtMOMX_P(:) * VelP(:) - GsqrtMOMX_M(:) * VelM(:) ) &
                    + (  ( G1n_M(:) + Gxz_P(:) * nz(:,ke)) * dpresP(:)      &
                       - ( G1n_M(:) + Gxz_M(:) * nz(:,ke)) * dpresM(:) )    &
                    - alpha(:) * ( GsqrtMOMX_P(:) - GsqrtMOMX_M(:) )        )

      del_flux(:,ke,MOMY_VID ) = 0.5_RP * ( &
                    ( GsqrtMOMY_P(:) * VelP(:) - GsqrtMOMY_M(:) * VelM(:) ) &
                    + (  ( G2n_M(:) + Gyz_P(:) * nz(:,ke) ) * dpresP(:)     &
                       - ( G2n_M(:) + Gyz_M(:) * nz(:,ke) ) * dpresM(:) )   &
                    - alpha(:) * ( GsqrtMOMY_P(:) - GsqrtMOMY_M(:) )        )

      del_flux(:,ke,MOMZ_VID ) = 0.5_RP * ( &
                    ( GsqrtMOMZ_P(:) * VelP(:) - GsqrtMOMZ_M(:) * VelM(:) ) &
                    - alpha(:) * ( GsqrtMOMZ_P(:) - GsqrtMOMZ_M(:) )        )
                    
      del_flux(:,ke,DRHOT_VID) = 0.5_RP * ( &
                    ( GsqrtRhotP(:) * VelhP(:) - GsqrtRhotM(:) * VelhM(:) ) &
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


!OCL SERIAL  
  subroutine atm_dyn_dgm_globalnonhydro3d_hevi_cal_vi( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,             & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, & ! (in)
    Dz, Lift,                                                & ! (in)
    modalFilterFlag, VModalFilter,                           & ! (in)
    impl_fac, dt,                                            & ! (in)
    lmesh, elem, lmesh2D, elem2D                             ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
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
    class(SparseMat), intent(in) :: Dz, Lift
    logical, intent(in) :: modalFilterFlag
    class(ModalFilter), intent(in) :: VModalFilter
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt

    real(RP) :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: PROG_VARS00(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: b(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: b1D(elem%Nnode_v,PROG_VARS_NUM,lmesh%NeZ,elem%Nnode_h1D**2)
    integer :: ipiv(elem%Nnode_v*PROG_VARS_NUM*lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: alph(elem%NfpTot,lmesh%NeZ)
    real(RP) :: tend(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: DENS_hyd_z(elem%Np,lmesh%NeZ)
    real(RP) :: PRES_hyd_z(elem%Np,lmesh%NeZ)
    real(RP) :: GnnM_z(elem%Np,lmesh%NeZ)
    real(RP) :: G13_z(elem%Np,lmesh%NeZ)
    real(RP) :: G23_z(elem%Np,lmesh%NeZ)
    real(RP) :: GsqrtV_z(elem%Np,lmesh%NeZ)
    real(RP) :: nz(elem%NfpTot,lmesh%NeZ)
    integer :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer :: ke_x, ke_y, ke_z, ke, ke2D, p, v
    integer :: itr_lin, itr_nlin
    integer :: f, vs, ve, kl, ku, nz_1D
    integer :: ij, info
    logical :: is_converged


    real(RP), allocatable :: PmatBnd(:,:,:)
    !------------------------------------------------------------------------

    
    call PROF_rapstart( 'hevi_cal_vi_prep', 3)

    nz_1D = elem%Nnode_v * PROG_VARS_NUM * lmesh%NeZ
    kl = 2 * elem%Nnode_v * PROG_VARS_NUM - 1
    ku = kl
    allocate( PmatBnd(2*kl+ku+1,nz_1D,elem%Nnode_h1D**2) )

    !$omp parallel private(f, vs, ve)
    !$omp do
    do ke_z=1, lmesh%NeZ
      do f=1, elem%Nfaces_h
        vs = 1 + (f-1)*elem%Nfp_h
        ve = vs + elem%Nfp_h - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_h(:,f) + (ke_z-1)*elem%Np
      end do
      do f=1, elem%Nfaces_v
        vs = elem%Nfp_h*elem%Nfaces_h + 1 + (f-1)*elem%Nfp_v
        ve = vs + elem%Nfp_v - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_v(:,f) + (ke_z-1)*elem%Np
      end do
      vmapP(:,ke_z) = vmapM(:,ke_z)
    end do
    !$omp do
    do ke_z=1, lmesh%NeZ
      vs = elem%Nfp_h*elem%Nfaces_h + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z > 1) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,2) + (ke_z-2)*elem%Np

      vs = elem%Nfp_h*elem%Nfaces_h + elem%Nfp_v + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z < lmesh%NeZ) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1) + ke_z*elem%Np
    end do
    !$omp end parallel
    call PROF_rapend( 'hevi_cal_vi_prep', 3)

    do ke_y=1, lmesh%NeY
    do ke_x=1, lmesh%NeX

      call PROF_rapstart( 'hevi_cal_vi_get_var', 3)

      !$omp parallel do private( ke, ke2D )
      do ke_z=1, lmesh%NeZ
        ke = ke_x + (ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
        ke2D = lmesh%EMap3Dto2D(ke)

        PROG_VARS(:,DDENS_VID,ke_z) = DDENS_(:,ke)
        PROG_VARS(:,MOMX_VID,ke_z) = MOMX_(:,ke)
        PROG_VARS(:,MOMY_VID,ke_z) = MOMY_(:,ke)
        PROG_VARS(:,MOMZ_VID,ke_z) = MOMZ_(:,ke)
        PROG_VARS(:,DRHOT_VID,ke_z) = DRHOT_(:,ke)
        DENS_hyd_z(:,ke_z) = DENS_hyd(:,ke)
        PRES_hyd_z(:,ke_z) = PRES_hyd(:,ke)

        PROG_VARS0(:,:,ke_z) = PROG_VARS(:,:,ke_z)
        PROG_VARS00(:,:,ke_z) = PROG_VARS(:,:,ke_z)

        nz(:,ke_z) = lmesh%normal_fn(:,ke,3)
        G13_z(:,ke_z) = lmesh%GI3(:,ke,1)
        G23_z(:,ke_z) = lmesh%GI3(:,ke,2)
        GsqrtV_z(:,ke_z) = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2D)

        GnnM_z(:,ke_z) = ( &
            1.0_RP / lmesh%Gsqrt(:,ke)                                                &
          + G13_z(:,ke_z) * ( lmesh%GIJ(elem%IndexH2Dto3D,ke2D,1,1) * G13_z(:,ke_z)   &
                            + lmesh%GIJ(elem%IndexH2Dto3D,ke2D,1,2) * G23_z(:,ke_z) ) &
          + G23_z(:,ke_z) * ( lmesh%GIJ(elem%IndexH2Dto3D,ke2D,1,2) * G13_z(:,ke_z)   &
                            + lmesh%GIJ(elem%IndexH2Dto3D,ke2D,2,2) * G23_z(:,ke_z) ) )
      end do
      call PROF_rapend( 'hevi_cal_vi_get_var', 3)
      
      if ( abs(impl_fac) > 0.0_RP ) then
        call PROF_rapstart( 'hevi_cal_vi_itr', 3)

        ! G = (q^n+1 - q^n*) + impl_fac * A(q^n+1) = 0
        ! dG/dq^n+1 del[q] = - G(q^n*)
        do itr_nlin = 1, 1

          call vi_eval_Ax( Ax(:,:,:), alph,                  & ! (out)
            PROG_VARS, PROG_VARS0, DENS_hyd_z, PRES_hyd_z,   & ! (in)
            Dz, Lift,                                        & ! (in)
            GnnM_z, G13_z, G23_z, GsqrtV_z,                  & ! (in)
            modalFilterFlag, VModalFilter%FilterMat,         & ! (in)
            impl_fac, dt,                                    & ! (in) 
            lmesh, elem,                                     & ! (in)
            nz, vmapM, vmapP, ke_x, ke_y, .false.            ) ! (in)

          do ke_z=1, lmesh%NeZ
            b(:,:,ke_z) = - Ax(:,:,ke_z) + PROG_VARS00(:,:,ke_z)
          end do

          call PROF_rapstart( 'hevi_cal_vi_matbnd', 3)

          call vi_construct_matbnd( PmatBnd,                & ! (out)
            kl, ku, nz_1D,                                  & ! (in)
            PROG_VARS0, DENS_hyd_z, PRES_hyd_z,             & ! (in)
            alph,                                           & ! (in)
            Dz, Lift,                                       & ! (in)
            modalFilterFlag, VModalFilter%FilterMat,        & ! (in)
            impl_fac, dt,                                   & ! (in)
            lmesh, elem,                                    & ! (in)
            nz, vmapM, vmapP, ke_x, ke_y                    ) ! (in)

          call PROF_rapend( 'hevi_cal_vi_matbnd', 3)
          
          call PROF_rapstart( 'hevi_cal_vi_lin', 3)
          !$omp parallel private(ij, v, ke_z, info)
          !$omp do
          do ij=1, elem%Nnode_h1D**2
            do ke_z=1, lmesh%NeZ
            do v=1, PROG_VARS_NUM
              b1D(:,v,ke_z,ij) = b(elem%Colmask(:,ij),v,ke_z)
            end do
            end do

            ! if ( lmesh%PRC_myrank==11 .and. ke_x==3 .and. ke_y==3 .and. ij == 1) then
            !   LOG_INFO("check_check_vi",*) "V1", PROG_VARS(elem%Colmask(:,ij),1,:)
            !   LOG_INFO("check_check_vi",*) "V2", PROG_VARS(elem%Colmask(:,ij),2,:)
            !   LOG_INFO("check_check_vi",*) "V3", PROG_VARS(elem%Colmask(:,ij),3,:)
            !   LOG_INFO("check_check_vi",*) "V4", PROG_VARS(elem%Colmask(:,ij),4,:)
            !   LOG_INFO("check_check_vi",*) "V5", PROG_VARS(elem%Colmask(:,ij),5,:)
            !   LOG_INFO("check_check_vi",*) "b1D1", b1D(:,1,:,ij)
            !   LOG_INFO("check_check_vi",*) "b1D2", b1D(:,2,:,ij)
            !   LOG_INFO("check_check_vi",*) "b1D3", b1D(:,3,:,ij)
            !   LOG_INFO("check_check_vi",*) "b1D4", b1D(:,4,:,ij)
            !   LOG_INFO("check_check_vi",*) "b1D5", b1D(:,5,:,ij)
            !   LOG_INFO("check_check_vi",*) "PmatBand", PmatBnd(:,:,ij)
            ! end if

            call dgbsv( nz_1D, kl, ku, 1, PmatBnd(:,:,ij), 2*kl+ku+1, ipiv(:,ij), b1D(:,:,:,ij), nz_1D, info)

            do ke_z=1, lmesh%NeZ
            do v=1, PROG_VARS_NUM
              PROG_VARS(elem%Colmask(:,ij),v,ke_z) = PROG_VARS(elem%Colmask(:,ij),v,ke_z) + b1D(:,v,ke_z,ij)
            end do
            end do
          end do
          !$omp do 
          do ke_z=1, lmesh%NeZ
            PROG_VARS0(:,:,ke_z) = PROG_VARS(:,:,ke_z)
          end do
          !$omp end parallel
          call PROF_rapend( 'hevi_cal_vi_lin', 3)
        end do ! itr nlin

        call PROF_rapend( 'hevi_cal_vi_itr', 3)
      end if

      call PROF_rapstart( 'hevi_cal_vi_retrun_var', 3)
      if ( abs(impl_fac) > 0.0_RP) then
        !$omp parallel do 
        do ke_z=1, lmesh%NeZ
          tend(:,:,ke_z) = (- PROG_VARS(:,:,ke_z) + PROG_VARS00(:,:,ke_z))/impl_fac
        end do
      else
        call vi_eval_Ax( tend(:,:,:), alph,                & ! (out)
          PROG_VARS, PROG_VARS, DENS_hyd_z, PRES_hyd_z,    & ! (in)
          Dz, Lift,                                        & ! (in)
          GnnM_z, G13_z, G23_z, GsqrtV_z,                  & ! (in)
          modalFilterFlag, VModalFilter%FilterMat,         & ! (in)
          impl_fac, dt,                                    & ! (in) 
          lmesh, elem,                                     & ! (in)
          nz, vmapM, vmapP, ke_x, ke_y, .true. )             ! (in)
      end if

!      if (lmesh%PRC_myrank==11) then
      ! do ke_z=1, lmesh%NeZ
      !   ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
      !   do p=1, elem%Np
      !       if ( IEEE_IS_NAN(tend(p,DDENS_VID,ke_z)) .and. ke_x==3 .and. ke_y==3 .and. p==1 .and. ke_z == 1) then
      !         LOG_INFO("check_nan_vi",*) impl_fac, lmesh%PRC_myrank, "Nan:", lmesh%pos_en(p,ke,:), ":", MOMX_(p,ke), MOMY_(p,ke), MOMZ_(p,ke), DRHOT_(p,ke)
      !         LOG_INFO("check_nan_vi",*) p, ke_x, ke_y, " Nan_tend:", tend(p,:,ke_z)
      !       end if
      !     end do
      !   end do
!      end if

      !$omp parallel do private(ke)
      do ke_z=1, lmesh%NeZ
        ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
        DENS_dt(:,ke) = - tend(:,DDENS_VID,ke_z)
        MOMX_dt(:,ke) = - tend(:,MOMX_VID ,ke_z)
        MOMY_dt(:,ke) = - tend(:,MOMY_VID ,ke_z)
        MOMZ_dt(:,ke) = - tend(:,MOMZ_VID ,ke_z)
        RHOT_dt(:,ke) = - tend(:,DRHOT_VID,ke_z)
      end do
      call PROF_rapend( 'hevi_cal_vi_retrun_var', 3)
    end do
    end do

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_hevi_cal_vi

  !------------------------------------------------

!OCL SERIAL  
  subroutine vi_eval_Ax( Ax, alph,              & ! (out)
    PROG_VARS, PROG_VARS0, DENS_hyd, PRES_hyd,  & ! (in)
    Dz, Lift,                                   & ! (in)
    GnnM_z, G13_z, G23_z, GsqrtV_z,             & ! (in)
    modalFilterFlag, VModalFilter,              & ! (in)
    impl_fac, dt,                               & ! (in)
    lmesh, elem,                                & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y, cal_tend_flag ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(out) :: alph(elem%NfpTot,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: GnnM_z(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: G13_z(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: G23_z(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: GsqrtV_z(elem%Np,lmesh%NeZ)
    logical, intent(in) :: modalFilterFlag
    real(RP), intent(in) :: VModalFilter(elem%Nnode_v,elem%Nnode_v)    
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y
    logical, intent(in) :: cal_tend_flag

    real(RP) :: RGsqrtV(elem%Np)    
    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%NeZ,PROG_VARS_NUM)
    real(RP) :: RHOT_hyd(elem%Np), POT(elem%Np)
    real(RP) :: DPRES(elem%Np)
    integer :: ke_z
    integer :: ke, ke2d
    integer :: v
    integer :: ij
    real(RP) :: gamm, rgamm
    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    call vi_cal_del_flux_dyn( del_flux, alph,              & ! (out)
      PROG_VARS(:,DDENS_VID,:), PROG_VARS(:,MOMX_VID,:),   & ! (in)
      PROG_VARS(:,MOMY_VID ,:), PROG_VARS(:,MOMZ_VID,:),   & ! (in)
      PROG_VARS(:,DRHOT_VID,:),                            & ! (in)
      PROG_VARS0(:,DDENS_VID,:), PROG_VARS0(:,MOMX_VID,:), & ! (in)
      PROG_VARS0(:,MOMY_VID ,:), PROG_VARS0(:,MOMZ_VID,:), & ! (in)
      PROG_VARS0(:,DRHOT_VID,:),                           & ! (in)
      DENS_hyd, PRES_hyd,                                  & ! (in)
      GnnM_z, G13_z, G23_z, GsqrtV_z, nz, vmapM, vmapP,    & ! (in)
      lmesh, elem )                                          ! (in)

    !$omp parallel do private( &
    !$omp ke, ke2d, RHOT_hyd, DPRES, POT, Fz, LiftDelFlx, &
    !$omp v, ij, RGsqrtV                                  &                                      
    !$omp )
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
      ke2d = lmesh%EMap3Dto2D(ke)

      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(:,ke_z)/PRES00)**rgamm

      DPRES(:) = PRES_hyd(:,ke_z) * ( (1.0_RP + PROG_VARS(:,DRHOT_VID,ke_z)/RHOT_hyd(:))**gamm - 1.0_RP )
      POT(:) = ( RHOT_hyd(:) + PROG_VARS(:,DRHOT_VID,ke_z) ) / ( DENS_hyd(:,ke_z) + PROG_VARS(:,DDENS_VID,ke_z) )

      RGsqrtV(:)  = 1.0_RP / GsqrtV_z(:,ke_z)

      !- DENS
      call sparsemat_matmul(Dz, PROG_VARS(:,MOMZ_VID,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,DDENS_VID), LiftDelFlx)
      Ax(:,DDENS_VID,ke_z) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) ) * RGsqrtV(:)

      !- MOMX
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,MOMX_VID), LiftDelFlx)
      Ax(:,MOMX_VID,ke_z) = LiftDelFlx(:) * RGsqrtV(:)

      !-MOMY
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,MOMY_VID), LiftDelFlx)
      Ax(:,MOMY_VID,ke_z) = LiftDelFlx(:) * RGsqrtV(:)

      !-MOMZ
      call sparsemat_matmul(Dz, DPRES(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,MOMZ_VID), LiftDelFlx)
      Ax(:,MOMZ_VID,ke_z) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) ) * RGsqrtV(:) &
                          + Grav * matmul(IntrpMat_VPOrdM1, PROG_VARS(:,DDENS_VID,ke_z)) 

      !-RHOT
      call sparsemat_matmul(Dz, POT(:)*PROG_VARS(:,MOMZ_VID,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,DRHOT_VID), LiftDelFlx)
      Ax(:,DRHOT_VID,ke_z) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) ) * RGsqrtV(:)

      !-- Modal filtering in the vertical direction      
      if ( modalFilterFlag ) then
        do v=1, PROG_VARS_NUM
          do ij=1, elem%Nnode_h1D**2
            Ax(elem%Colmask(:,ij),v,ke_z) = Ax(elem%Colmask(:,ij),v,ke_z)      &
              - matmul(VModalFilter, PROG_VARS(elem%Colmask(:,ij),v,ke_z) ) / dt
          end do
        end do
      end if

      !--
      if ( .not. cal_tend_flag ) then
        Ax(:,:,ke_z) = PROG_VARS(:,:,ke_z) + impl_fac * Ax(:,:,ke_z)
      end if 

    end do    

    return
  end subroutine vi_eval_Ax

!OCL SERIAL  
  subroutine vi_cal_del_flux_dyn( del_flux, alph,             & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,                      & ! (in)
    DDENS0_, MOMX0_, MOMY0_, MOMZ0_, DRHOT0_,                 & ! (in)
    DENS_hyd, PRES_hyd,                                       & ! (in)
    Gnn_, G13_, G23_, GsqrtV_, nz, vmapM, vmapP, lmesh, elem  ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%NeZ,PROG_VARS_NUM)
    real(RP), intent(out) :: alph(elem%NfpTot*lmesh%NeZ)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  DDENS0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMX0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMY0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMZ0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  DRHOT0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  Gnn_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  G13_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  G23_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  GsqrtV_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%NeZ)
    
    integer :: i, p, ke_z, iP, iM, iM2D
    real(RP) :: alpha0
    real(RP) :: MOMZ_P
    real(RP) :: wt0M, wt0P
    real(RP) :: rhot_hyd_M, rhot_hyd_P
    real(RP) :: dpresM, dpresP, densM, densP,  pottM, pottP
    real(RP) :: pres0M, pres0P, dens0M, dens0P
    real(RP) :: gamm, rgamm
    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry
    
    !$omp parallel do private( p, i, iM, iP,                  &
    !$omp rhot_hyd_M, rhot_hyd_P, densM, densP, pottM, pottP, &
    !$omp dpresM, dpresP, MOMZ_P, wt0M, wt0P,                 &
    !$omp dens0M, dens0P, pres0M, pres0P                      )    
    do ke_z=1, lmesh%NeZ
    do p=1, elem%NfpTot
      i = p + (ke_z-1)*elem%NfpTot
      iM = vmapM(i); iP = vmapP(i)
      iM2D = iM2Dto3D(p)

      !-
      rhot_hyd_M = PRES00/Rdry * (PRES_hyd(iM)/PRES00)**rgamm
      rhot_hyd_P = PRES00/Rdry * (PRES_hyd(iP)/PRES00)**rgamm

      densM = DENS_hyd(iM) + DDENS_(iM)
      densP = DENS_hyd(iP) + DDENS_(iP)
      
      pottM = (rhot_hyd_M + DRHOT_(iM)) / densM
      pottP = (rhot_hyd_P + DRHOT_(iP)) / densP
   
      dpresM = PRES_hyd(iM) * ((1.0_RP + DRHOT_(iM)/rhot_hyd_M)**gamm - 1.0_RP) 
      dpresP = PRES_hyd(iP) * ((1.0_RP + DRHOT_(iP)/rhot_hyd_P)**gamm - 1.0_RP) 

      !-
      dens0M = DENS_hyd(iM) + DDENS0_(iM)
      dens0P = DENS_hyd(iP) + DDENS0_(iP)

      pres0M = PRES_hyd(iM) * (1.0_RP + DRHOT0_(iM)/rhot_hyd_M)**gamm
      pres0P = PRES_hyd(iP) * (1.0_RP + DRHOT0_(iP)/rhot_hyd_P)**gamm

      wt0M = ( MOMZ0_(iM) / GsqrtV_(iM) + G13_(iM) * MOMX0_(iM) + G23_(iM) * MOMY0_(iM) ) / dens0M      
      wt0P = ( MOMZ0_(iP) / GsqrtV_(iP) + G13_(iP) * MOMX0_(iP) + G23_(iP) * MOMY0_(iP) ) / dens0M 

      alph(i) = nz(i)**2 * max( abs( wt0M ) + sqrt( Gnn_(iM) * gamm * pres0M / dens0M ), &
                                abs( wt0P ) + sqrt( Gnn_(iP) * gamm * pres0M / dens0M )  )
      
      !----                                
      if (iM==iP .and. (ke_z == 1 .or. ke_z == lmesh%NeZ)) then
        MOMZ_P = - MOMZ_(iM)
        !alph(i) = 0.0_RP
      else
        MOMZ_P = MOMZ_(iP)
      end if

      del_flux(i,DDENS_VID) = 0.5_RP * (                    &
                    + ( MOMZ_P - MOMZ_(iM) ) * nz(i)        &
                    - alph(i) * ( DDENS_(iP) - DDENS_(iM) ) )
      
      del_flux(i,MOMX_VID) = 0.5_RP * (                     &
                    - alph(i) * (  MOMX_(iP) - MOMX_(iM) )  )
      
      del_flux(i,MOMY_VID) = 0.5_RP * (                     &  
                    - alph(i) * ( MOMY_(iP) - MOMY_(iM) )   )               
      
      del_flux(i,MOMZ_VID) = 0.5_RP * (                     &
                    + ( dpresP - dpresM ) * nz(i)           &                    
                    - alph(i) * ( MOMZ_P - MOMZ_(iM) )      )
      
      del_flux(i,DRHOT_VID) = 0.5_RP * (                               &
                    + ( pottP * MOMZ_P  - pottM * MOMZ_(iM) ) * nz(i)  &
                    - alph(i) * ( DRHOT_(iP) - DRHOT_(iM) )            )
    end do
    end do

    return
  end subroutine vi_cal_del_flux_dyn

!OCL SERIAL  
  subroutine vi_construct_matbnd( PmatBnd,  & ! (out)
    kl, ku, nz_1D,                          & ! (in)
    PROG_VARS0, DENS_hyd, PRES_hyd,         & ! (in)
    alph,                                   & ! (in)
    Dz, Lift,                               & ! (in)
    modalFilterFlag, VModalFilter,          & ! (in)
    impl_fac, dt,                           & ! (in)
    lmesh, elem,                            & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y )            ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(in) :: kl, ku, nz_1D
    real(RP), intent(out) :: PmatBnd(2*kl+ku+1,elem%Nnode_v,PROG_VARS_NUM,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: alph(elem%NfpTot,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    logical, intent(in) :: modalFilterFlag
    real(RP), intent(in) :: VModalFilter(elem%Nnode_v,elem%Nnode_v)
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y

    real(RP) :: RHOT_hyd(elem%Nnode_v)
    real(RP) :: POT0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: W0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: DENS0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: DPDRHOT0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    integer :: ke_z, ke_z2
    integer :: v, ke, p, f1, f2, fp, fp2, FmV
    real(RP) :: gamm, rgamm
    real(RP) :: fac_dz_p(elem%Nnode_v)
    real(RP) :: PmatD(elem%Nnode_v,elem%Nnode_v,PROG_VARS_NUM,PROG_VARS_NUM)
    real(RP) :: PmatL(elem%Nnode_v,elem%Nnode_v,PROG_VARS_NUM,PROG_VARS_NUM)
    real(RP) :: PmatU(elem%Nnode_v,elem%Nnode_v,PROG_VARS_NUM,PROG_VARS_NUM)
    integer :: Colmask(elem%Nnode_v)
    real(RP) :: Id(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: Dd(elem%Nnode_v)
    real(RP) :: tmp1
    real(RP) :: fac

    integer :: ij, v1, v2, pv1, pv2,  g_kj, g_kjp1, g_kjm1, pb, pb1
    logical :: bc_flag
    logical :: eval_flag(PROG_VARS_NUM,PROG_VARS_NUM)
    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    eval_flag(:,:) = .false.
    do v=1, PROG_VARS_NUM
      eval_flag(v,v) = .true.
    end do
    eval_flag(DDENS_VID,MOMZ_VID) = .true.
    eval_flag(MOMZ_VID,DDENS_VID) = .true.
    eval_flag(MOMZ_VID,DRHOT_VID) = .true.
    eval_flag(DRHOT_VID,MOMZ_VID) = .true.
    eval_flag(DRHOT_VID,DDENS_VID) = .true.

    Id(:,:) = 0.0_RP
    do p=1, elem%Nnode_v
      Id(p,p) = 1.0_RP
    end do

    !$omp parallel private(RHOT_hyd, Colmask)
    !$omp do
    do v=1, PROG_VARS_NUM
      PmatD(:,:,:,v) = 0.0_RP
      PmatL(:,:,:,v) = 0.0_RP
      PmatU(:,:,:,v) = 0.0_RP  
    end do
    !$omp do
    do ij=1, elem%Nnode_h1D**2
      PmatBnd(:,:,:,:,ij) = 0.0_RP
    end do
    !$omp do collapse(2)
    do ij=1, elem%Nnode_h1D**2
    do ke_z=1, lmesh%NeZ
      Colmask(:) = elem%Colmask(:,ij)
      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(Colmask(:),ke_z)/PRES00)**rgamm

      DPDRHOT0(:,ke_z,ij) = gamm * PRES_hyd(Colmask(:),ke_z) / RHOT_hyd(:)                      &
                  * ( 1.0_RP + PROG_VARS0(Colmask(:),DRHOT_VID,ke_z) / RHOT_hyd(:) )**(gamm-1) 

      DENS0(:,ke_z,ij) = DENS_hyd(Colmask(:),ke_z) + PROG_VARS0(Colmask(:),DDENS_VID,ke_z)
      POT0(:,ke_z,ij) = ( RHOT_hyd(:) + PROG_VARS0(Colmask(:),DRHOT_VID,ke_z) ) / DENS0(:,ke_z,ij)
      W0(:,ke_z,ij) = PROG_VARS0(Colmask(:),MOMZ_VID,ke_z) / DENS0(:,ke_z,ij)
    end do
    end do
    !$omp end parallel

    !$omp parallel do private(ke_z, ke, ColMask, p, fp, fp2, v, f1, f2, ke_z2, fac_dz_p, &
    !$omp fac, tmp1, FmV,                                                                &
    !$omp ij, v1, v2, pv1, pv2, pb1, g_kj, g_kjp1, g_kjm1, bc_flag,                      &
    !$omp Dd                                                                           ) &
    !$omp firstprivate(PmatD, PmatL, PmatU)
    do ij=1, elem%Nnode_h1D**2
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
      Colmask(:) = elem%Colmask(:,ij)

      !-----  
      do p=1, elem%Nnode_v
        fac_dz_p(:) = impl_fac * lmesh%Escale(Colmask(:),ke,3,3) * elem%Dx3(Colmask(:),Colmask(p))
        if (modalFilterFlag) then
          Dd(:) = Id(:,p) - VModalFilter(:,p) * impl_fac / dt
        else
          Dd(:) = Id(:,p)
        end if

        ! DDENS
        PmatD(:,p,DDENS_VID,DDENS_VID) = Dd(:)
        PmatD(:,p,DDENS_VID,MOMZ_VID) = fac_dz_p(:) 

        ! MOMX
        PmatD(:,p,MOMX_VID,MOMX_VID) = Dd(:)

        ! MOMY
        PmatD(:,p,MOMY_VID,MOMY_VID) = Dd(:)

        ! MOMZ
        PmatD(:,p,MOMZ_VID,MOMZ_VID) = Dd(:) 
        PmatD(:,p,MOMZ_VID,DDENS_VID) = impl_fac * Grav * IntrpMat_VPOrdM1(Colmask(:),Colmask(p))
        PmatD(:,p,MOMZ_VID,DRHOT_VID) = fac_dz_p(:) * DPDRHOT0(p,ke_z,ij)

        !DRHOT
        PmatD(:,p,DRHOT_VID,DDENS_VID) = - fac_dz_p(:) * POT0(p,ke_z,ij) * W0(p,ke_z,ij)
        PmatD(:,p,DRHOT_VID,MOMZ_VID ) =   fac_dz_p(:) * POT0(p,ke_z,ij)
        PmatD(:,p,DRHOT_VID,DRHOT_VID) = Dd(:) + fac_dz_p(:) * W0(p,ke_z,ij)
      end do

      do f1=1, 2
        if (f1==1) then
          ke_z2 = max(ke_z-1,1)
          pv1 = 1; pv2 = elem%Nnode_v
          f2 = 2
        else
          ke_z2 = min(ke_z+1,lmesh%NeZ)
          pv1 = elem%Nnode_v; pv2 = 1
          f2 = 1
        end if
        fac  = 0.5_RP * impl_fac
        if ( (ke_z == 1 .and. f1==1) .or. (ke_z == lmesh%NeZ .and. f1==elem%Nfaces_v) ) then
          bc_flag = .true.
          pv2 = pv1; f2 = f1
        else 
          bc_flag = .false.    
        end if

        FmV  = elem%Fmask_v(ij,f1)
        fp  = elem%Nfp_h * elem%Nfaces_h + (f1-1)*elem%Nfp_v + ij
        fp2 = elem%Nfp_h * elem%Nfaces_h + (f2-1)*elem%Nfp_v + ij

        !--
        tmp1 = fac * elem%Lift(FmV,fp) * lmesh%Fscale(fp,ke) &
             * max( alph(fp,ke_z), alph(fp2,ke_z2) )
        if (bc_flag) then
          PmatD(pv1,pv1,MOMZ_VID,MOMZ_VID) = PmatD(pv1,pv1,MOMZ_VID,MOMZ_VID) + 2.0_RP * tmp1
        else 
          do v=1, PROG_VARS_NUM
            PmatD(pv1,pv1,v,v) = PmatD(pv1,pv1,v,v) + tmp1            
            if (f1 == 1) then
              PmatL(pv1,pv2,v,v) = - tmp1                                
            else
              PmatU(pv1,pv2,v,v) = - tmp1
            end if
          end do
        end if 

        !--
        tmp1 = fac * elem%Lift(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)

        if (bc_flag) then
          PmatD(pv1,pv1,DDENS_VID,MOMZ_VID ) = PmatD(pv1,pv1,DDENS_VID,MOMZ_VID ) - 2.0_RP * tmp1
          PmatD(pv1,pv1,DRHOT_VID,MOMZ_VID ) = PmatD(pv1,pv1,DRHOT_VID,MOMZ_VID ) - 2.0_RP * tmp1 * POT0(pv1,ke_z,ij)
          PmatD(pv1,pv1,DRHOT_VID,DDENS_VID) = PmatD(pv1,pv1,DRHOT_VID,DDENS_VID) + 2.0_RP * tmp1 * POT0(pv1,ke_z,ij) * W0(pv1,ke_z,ij)
          PmatD(pv1,pv1,DRHOT_VID,DRHOT_VID) = PmatD(pv1,pv1,DRHOT_VID,DRHOT_VID) - 2.0_RP * tmp1 * W0(pv1,ke_z,ij)
        else 
          PmatD(pv1,pv1,DDENS_VID,MOMZ_VID ) = PmatD(pv1,pv1,DDENS_VID,MOMZ_VID ) - tmp1
          PmatD(pv1,pv1,MOMZ_VID ,DRHOT_VID) = PmatD(pv1,pv1,MOMZ_VID ,DRHOT_VID) - tmp1 * DPDRHOT0(pv1,ke_z,ij)
          PmatD(pv1,pv1,DRHOT_VID,MOMZ_VID ) = PmatD(pv1,pv1,DRHOT_VID,MOMZ_VID ) - tmp1 * POT0(pv1,ke_z,ij)
          PmatD(pv1,pv1,DRHOT_VID,DDENS_VID) = PmatD(pv1,pv1,DRHOT_VID,DDENS_VID) + tmp1 * POT0(pv1,ke_z,ij) * W0(pv1,ke_z,ij)
          PmatD(pv1,pv1,DRHOT_VID,DRHOT_VID) = PmatD(pv1,pv1,DRHOT_VID,DRHOT_VID) - tmp1 * W0(pv1,ke_z,ij)

          if (f1 == 1) then
            PmatL(pv1,pv2,DDENS_VID,MOMZ_VID ) = + tmp1
            PmatL(pv1,pv2,MOMZ_VID,DRHOT_VID ) = + tmp1 * DPDRHOT0(pv2,ke_z2,ij) 
            PmatL(pv1,pv2,DRHOT_VID,MOMZ_VID ) = + tmp1 * POT0(pv2,ke_z2,ij)
            PmatL(pv1,pv2,DRHOT_VID,DDENS_VID) = - tmp1 * POT0(pv2,ke_z2,ij) * W0(pv2,ke_z2,ij)
            PmatL(pv1,pv2,DRHOT_VID,DRHOT_VID) = PmatL(pv1,pv2,DRHOT_VID,DRHOT_VID) &
                                                 + tmp1 * W0(pv2,ke_z2,ij)
          else
            PmatU(pv1,pv2,DDENS_VID,MOMZ_VID ) = + tmp1
            PmatU(pv1,pv2,MOMZ_VID,DRHOT_VID ) = + tmp1 * DPDRHOT0(pv2,ke_z2,ij)     
            PmatU(pv1,pv2,DRHOT_VID,MOMZ_VID ) = + tmp1 * POT0(pv2,ke_z2,ij)
            PmatU(pv1,pv2,DRHOT_VID,DDENS_VID) = - tmp1 * POT0(pv2,ke_z2,ij) * W0(pv2,ke_z2,ij)
            PmatU(pv1,pv2,DRHOT_VID,DRHOT_VID) = PmatU(pv1,pv2,DRHOT_VID,DRHOT_VID) &
                                                 + tmp1 * W0(pv2,ke_z2,ij)
          end if
        end if
      end do

      do v2=1, PROG_VARS_NUM
      do v1=1, PROG_VARS_NUM
        if ( eval_flag(v1,v2) ) then
          do pv2=1, elem%Nnode_v
            g_kj   = pv2 + (v2-1)*elem%Nnode_v + (ke_z-1)*elem%Nnode_v*PROG_VARS_NUM
            g_kjm1 = pv2 + (v2-1)*elem%Nnode_v + (ke_z-2)*elem%Nnode_v*PROG_VARS_NUM
            g_kjp1 = pv2 + (v2-1)*elem%Nnode_v + (ke_z  )*elem%Nnode_v*PROG_VARS_NUM

            do pv1=1, elem%Nnode_v            
              pb1 = pv1 + (v1-1)*elem%Nnode_v + (ke_z-1)*elem%Nnode_v*PROG_VARS_NUM
              if (ke_z > 1) then
                PmatBnd(kl+ku+1+pb1-g_kjm1, pv2,v2,ke_z-1, ij) = PmatL(pv1,pv2,v1,v2)
              end if
              PmatBnd(kl+ku+1+pb1-g_kj, pv2,v2,ke_z, ij) = PmatD(pv1,pv2,v1,v2)
              if (ke_z < lmesh%NeZ) then
                PmatBnd(kl+ku+1+pb1-g_kjp1, pv2,v2,ke_z+1, ij) = PmatU(pv1,pv2,v1,v2)
              end if
            end do
          end do            
        end if
      end do
      end do
      end do  
    end do  

    return
  end subroutine vi_construct_matbnd

end module scale_atm_dyn_dgm_globalnonhydro3d_hevi
