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
    SL_flag, wdamp_tau, wdamp_height,                                           & ! (in)
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

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: DPRES_(elem%Np)
    real(RP) :: RHOT_(elem%Np)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), drho(elem%Np)

    real(RP) :: GIJ(elem%Np,2,2)
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
    !$omp RHOT_, DPRES_, rdens_, u_, v_, w_,       &
    !$omp Fx, Fy, Fz, LiftDelFlx,                  &
    !$omp drho, GradPhyd_x, GradPhyd_y,            &
    !$omp GIJ, X, Y, twoOVdel2,                    &
    !$omp CORI, ke, ke2D                           )

    !$omp do
    do ke2D = lmesh2D%NeS, lmesh2D%NeE
      X2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,1))
      Y2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,2))
    end do

    !$omp do
    do ke = lmesh%NeS, lmesh%NeE
      !--
      RHOT_(:) = P0ovR * ( PRES_hyd(:,ke) * rP0 )**rgamm + DRHOT_(:,ke)
      DPRES_(:) = PRES00 * ( RovP0 * RHOT_(:) )**gamm &
                - PRES_hyd(:,ke)

      rdens_(:) = 1.0_RP / ( DDENS_(:,ke) + DENS_hyd(:,ke) )
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)

      ke2d = lmesh%EMap3Dto2D(ke)
      GIJ(:,1,1) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,1)
      GIJ(:,2,1) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,2,1)
      GIJ(:,1,2) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,2)
      GIJ(:,2,2) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,2,2)

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
      call sparsemat_matmul(Dx, PRES_hyd(:,ke), Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,1), LiftDelFlx)
      GradPhyd_x(:) = lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dy, PRES_hyd(:,ke), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,2), LiftDelFlx)
      GradPhyd_y(:) = lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:)

      !-- DENS
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * MOMY_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * MOMZ_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_DDENS_ID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)     &
          + lmesh%Escale(:,ke,2,2) * Fy(:)     &
          + lmesh%Escale(:,ke,3,3) * Fz(:)     &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      !-- MOMX
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_(:) * MOMX_(:,ke) + GIJ(:,1,1) * DPRES_(:) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_(:) * MOMX_(:,ke) + GIJ(:,1,2) * DPRES_(:) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) *   w_(:) * MOMX_(:,ke)                           , Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_MOMX_ID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                     &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                     &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                     &
            + LiftDelFlx(:)                   ) / lmesh%Gsqrt(:,ke)              &
          + twoOVdel2(:) * Y(:) *                                                &
            ( - X(:) * Y(:) * u_(:) + (1.0_RP + Y(:)**2) * v_(:) ) * MOMX_(:,ke) &
          - ( GIJ(:,1,1) * GradPhyd_x(:) + GIJ(:,1,2) * GradPhyd_y(:) )          &
          + CORI(:,1)              

      !-- MOMY
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_(:) * MOMY_(:,ke) + GIJ(:,2,1) * DPRES_(:) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_(:) * MOMY_(:,ke) + GIJ(:,2,2) * DPRES_(:) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) *   w_(:) * MOMY_(:,ke)                           , Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_MOMY_ID), LiftDelFlx)

      MOMY_dt(:,ke) = &
            - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                   &
              + lmesh%Escale(:,ke,2,2) * Fy(:)                                   &
              + lmesh%Escale(:,ke,3,3) * Fz(:)                                   &
              + LiftDelFlx(:)                  ) / lmesh%Gsqrt(:,ke)             &
            + twoOVdel2(:) * X(:) *                                              &
              ( (1.0_RP + X(:)**2) * u_(:) - X(:) * Y(:) * v_(:) ) * MOMY_(:,ke) &
            - ( GIJ(:,2,1) * GradPhyd_x(:) + GIJ(:,2,2) * GradPhyd_y(:) )        &
            + CORI(:,2) 

      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *   u_(:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *   v_(:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( w_(:) * MOMZ_(:,ke) + DPRES_(:) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_MOMZ_ID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)       &
            + lmesh%Escale(:,ke,2,2) * Fy(:)       &
            + lmesh%Escale(:,ke,3,3) * Fz(:)       &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)  &
          - Grav * drho(:)

      !-- RHOT
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * u_(:) * RHOT_(:), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * v_(:) * RHOT_(:), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * w_(:) * RHOT_(:), Fz)
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
      call atm_dyn_dgm_spongelayer_add_tend( MOMZ_dt, &
        MOMZ_, wdamp_tau, wdamp_height, lmesh, elem   )
      call PROF_rapend('cal_dyn_tend_sponge', 3)
    end if

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_heve_cal_tend

  !------

!OCL SERIAL
  subroutine cal_del_flux_dyn( del_flux, del_flux_hyd,             &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,       &
    Gsqrt, G11, G12, G22, nx, ny, nz,                              &
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
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%Ne)
    real(RP), intent(in) ::  G11(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G12(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G22(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    
    integer :: ke, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: VelP(elem%NfpTot), VelM(elem%NfpTot), alpha(elem%NfpTot)
    real(RP) :: dpres(elem%NfpTot), presM(elem%NfpTot), presP(elem%NfpTot)
    real(RP) :: densM(elem%NfpTot), densP(elem%NfpTot)
    real(RP) :: rhotM(elem%NfpTot), rhotP(elem%NfpTot)
    real(RP) :: DDENS_P(elem%NfpTot), DDENS_M(elem%NfpTot)
    real(RP) :: MOMX_P(elem%NfpTot), MOMX_M(elem%NfpTot)
    real(RP) :: MOMY_P(elem%NfpTot), MOMY_M(elem%NfpTot)
    real(RP) :: MOMZ_P(elem%NfpTot), MOMZ_M(elem%NfpTot)
    real(RP) :: DRHOT_P(elem%NfpTot), DRHOT_M(elem%NfpTot)
    real(RP) :: PRES_hyd_P(elem%NfpTot), PRES_hyd_M(elem%NfpTot)
    real(RP) :: Gsqrt_M(elem%NfpTot)
    real(RP) :: G1n_M(elem%NfpTot), G2n_M(elem%NfpTot) , Gnn_M(elem%NfpTot)
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
    !$omp ke, iM, iP, ke2d,                                    &
    !$omp alpha, VelM, VelP,                                   &
    !$omp dpres, presM, presP, densM, densP, rhotM, rhotP,     &
    !$omp MOMX_M, MOMX_P, MOMY_M, MOMY_P, MOMZ_M, MOMZ_P,      &
    !$omp DDENS_M, DDENS_P, DRHOT_M, DRHOT_P,                  &
    !$omp PRES_hyd_M, PRES_hyd_P,                              &
    !$omp Gsqrt_M, G1n_M, G2n_M, Gnn_M                         )
    do ke=lmesh%NeS, lmesh%NeE

      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)

      DDENS_M(:) = DDENS_(iM)
      DDENS_P(:) = DDENS_(iP)
      MOMX_M(:) = MOMX_(iM)
      MOMX_P(:) = MOMX_(iP)
      MOMY_M(:) = MOMY_(iM)
      MOMY_P(:) = MOMY_(iP)
      MOMZ_M(:) = MOMZ_(iM)
      MOMZ_P(:) = MOMZ_(iP)
      DRHOT_M(:) = DRHOT_(iM)
      DRHOT_P(:) = DRHOT_(iP)
      PRES_hyd_M(:) = PRES_hyd(iM)
      PRES_hyd_P(:) = PRES_hyd(iP)

      Gsqrt_M(:) = Gsqrt(iM)

      ke2D = lmesh%EMap3Dto2D(ke)
      G1n_M(:)  = G11(iM2Dto3D(:),ke2D) * nx(:,ke) + G12(iM2Dto3D(:),ke2D) * ny(:,ke)
      G2n_M(:)  = G12(iM2Dto3D(:),ke2D) * nx(:,ke) + G22(iM2Dto3D(:),ke2D) * ny(:,ke)
      Gnn_M(:)  = G11(iM2Dto3D(:),ke2D) * abs( nx(:,ke) ) + G22(iM2Dto3D(:),ke2D) * abs( ny(:,ke) ) &
                + abs( nz(:,ke) )

      densM(:) = DDENS_M(:) + DENS_hyd(iM)
      densP(:) = DDENS_P(:) + DENS_hyd(iP)

      VelM(:) = ( MOMX_M(:) * nx(:,ke) + MOMY_M(:) * ny(:,ke) + MOMZ_M(:) * nz(:,ke) ) / densM(:)
      VelP(:) = ( MOMX_P(:) * nx(:,ke) + MOMY_P(:) * ny(:,ke) + MOMZ_P(:) * nz(:,ke) ) / densP(:)

      rhotM(:) = P0ovR * (PRES_hyd_M(:) * rP0)**rgamm + DRHOT_M(:)
      rhotP(:) = P0ovR * (PRES_hyd_P(:) * rP0)**rgamm + DRHOT_P(:)

      presM(:) = PRES00 * (RovP0 * rhotM(:))**gamm
      presP(:) = PRES00 * (RovP0 * rhotP(:))**gamm

      dpres(:)  =  presP(:) - presM(:)                           &
             - ( PRES_hyd_P(:) - PRES_hyd_M(:) ) * abs( nz(:,ke) )

      alpha(:) = max( sqrt( Gnn_M(:) * gamm * presM(:) / densM(:) ) + abs(VelM(:)), &
                      sqrt( Gnn_M(:) * gamm * presP(:) / densP(:) ) + abs(VelP(:))  )
      
      del_flux(:,ke,VARS_DDENS_ID) = 0.5_RP * Gsqrt_M(:) * (     &
                    ( densP(:) * VelP(:) - densM(:) * VelM(:) )  &
                    - alpha(:) * ( DDENS_P(:) - DDENS_M(:) )     )

      del_flux(:,ke,VARS_MOMX_ID) = 0.5_RP * Gsqrt_M(:) * (        &
                    ( MOMX_P(:) * VelP(:) - MOMX_M(:) * VelM(:) )  &
                    + G1n_M(:) * dpres(:)                          &
                    - alpha(:) * ( MOMX_P(:) - MOMX_M(:) )         )
      
      del_flux(:,ke,VARS_MOMY_ID) = 0.5_RP * Gsqrt_M(:) * (        &
                    ( MOMY_P(:) * VelP(:) - MOMY_M(:) * VelM(:) )  &
                    + G2n_M(:) * dpres(:)                          &
                    - alpha(:) * ( MOMY_P(:) - MOMY_M(:) )         )               
      
      del_flux(:,ke,VARS_MOMZ_ID) = 0.5_RP * Gsqrt_M(:) * (        &
                    ( MOMZ_P(:) * VelP(:) - MOMZ_M(:) * VelM(:) )  &
                    + dpres(:) * nz(:,ke)                          &                   
                    - alpha(:) * ( MOMZ_P(:) - MOMZ_M(:) )         )
      
      del_flux(:,ke,VARS_DRHOT_ID) = 0.5_RP * Gsqrt_M(:) * (       &
                    ( rhotP(:) * VelP(:) - rhotM(:) * VelM(:) )    &
                    - alpha(:) * ( DRHOT_P(:) - DRHOT_M(:) )       )

      del_flux_hyd(:,ke,1) = 0.5_RP * ( PRES_hyd_P(:) - PRES_hyd_M(:) ) * nx(:,ke)
      del_flux_hyd(:,ke,2) = 0.5_RP * ( PRES_hyd_P(:) - PRES_hyd_M(:) ) * ny(:,ke)              
    end do

    return
  end subroutine cal_del_flux_dyn

end module scale_atm_dyn_dgm_globalnonhydro3d_heve
