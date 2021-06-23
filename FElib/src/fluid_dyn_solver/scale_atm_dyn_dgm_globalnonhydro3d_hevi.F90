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

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    atm_dyn_dgm_nonhydro3d_common_Init,               &
    atm_dyn_dgm_nonhydro3d_common_Final,              &
    DENS_VID, MOMX_VID, MOMY_VID, MOMZ_VID, RHOT_VID, &
    PROG_VARS_NUM,                                    &
    IntrpMat_VPOrdM1, iM2Dto3D

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

contains
  subroutine atm_dyn_dgm_globalnonhydro3d_hevi_Init( mesh )
    implicit none
    class(MeshBase3D), intent(in) :: mesh
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Init( mesh )
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_hevi_Init


  subroutine atm_dyn_dgm_globalnonhydro3d_hevi_Final()
    implicit none
    !--------------------------------------------
    
    call atm_dyn_dgm_nonhydro3d_common_Final()    
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_hevi_Final  

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_hevi_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, CORIOLIS,          & ! (in)
    SL_flag, wdamp_tau, wdamp_height, hveldamp_flag,                            & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )

    use scale_atm_dyn_dgm_nonhydro3d_hevi_numflux, only: &
      atm_dyn_dgm_nonhydro3d_hevi_numflux_get_generalhvc
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
    call atm_dyn_dgm_nonhydro3d_hevi_numflux_get_generalhvc( &
      del_flux, del_flux_hyd,                                                  & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                 & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%GsqrtH, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),                        & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%vmapM, lmesh%vmapP, iM2Dto3D, lmesh, elem, lmesh2D, elem2D         ) ! (in)
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
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DENS_VID), LiftDelFlx)

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
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,RHOT_VID), LiftDelFlx)
      
      RHOT_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)

      ! DENS_dt(:,ke) = 0.0_RP
      ! MOMX_dt(:,ke) = 0.0_RP
      ! MOMY_dt(:,ke) = 0.0_RP
      ! MOMZ_dt(:,ke) = 0.0_RP
      ! RHOT_dt(:,ke) = 0.0_RP
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


!OCL SERIAL  
  subroutine atm_dyn_dgm_globalnonhydro3d_hevi_cal_vi( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,             & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, & ! (in)
    Dz, Lift,                                                & ! (in)
    modalFilterFlag, VModalFilter,                           & ! (in)
    impl_fac, dt,                                            & ! (in)
    lmesh, elem, lmesh2D, elem2D                             ) ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_hevi_common, only: &
      vi_gen_vmap => atm_dyn_dgm_nonhydro3d_hevi_common_gen_vmap,                &
      vi_eval_Ax => atm_dyn_dgm_nonhydro3d_hevi_common_eval_Ax,                  &
      vi_construct_matbnd => atm_dyn_dgm_nonhydro3d_hevi_common_construct_matbnd
  
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
    integer :: kl, ku, nz_1D
    integer :: ij, info
    logical :: is_converged

    real(RP), allocatable :: PmatBnd(:,:,:)
    !------------------------------------------------------------------------

    call PROF_rapstart( 'hevi_cal_vi_prep', 3)

    nz_1D = elem%Nnode_v * PROG_VARS_NUM * lmesh%NeZ
    kl = 2 * elem%Nnode_v * PROG_VARS_NUM - 1
    ku = kl
    allocate( PmatBnd(2*kl+ku+1,nz_1D,elem%Nnode_h1D**2) )

    call vi_gen_vmap( vmapM, vmapP, & ! (out)
      lmesh, elem                   ) ! (in)

    call PROF_rapend( 'hevi_cal_vi_prep', 3)

    do ke_y=1, lmesh%NeY
    do ke_x=1, lmesh%NeX

      call PROF_rapstart( 'hevi_cal_vi_get_var', 3)

      !$omp parallel do private( ke, ke2D )
      do ke_z=1, lmesh%NeZ
        ke = ke_x + (ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
        ke2D = lmesh%EMap3Dto2D(ke)

        PROG_VARS(:,DENS_VID,ke_z) = DDENS_(:,ke)
        PROG_VARS(:,MOMX_VID,ke_z) = MOMX_(:,ke)
        PROG_VARS(:,MOMY_VID,ke_z) = MOMY_(:,ke)
        PROG_VARS(:,MOMZ_VID,ke_z) = MOMZ_(:,ke)
        PROG_VARS(:,RHOT_VID,ke_z) = DRHOT_(:,ke)
        DENS_hyd_z(:,ke_z) = DENS_hyd(:,ke)
        PRES_hyd_z(:,ke_z) = PRES_hyd(:,ke)

        PROG_VARS0 (:,:,ke_z) = PROG_VARS(:,:,ke_z)
        PROG_VARS00(:,:,ke_z) = PROG_VARS(:,:,ke_z)

        nz(:,ke_z) = lmesh%normal_fn(:,ke,3)
        G13_z(:,ke_z) = lmesh%GI3(:,ke,1)
        G23_z(:,ke_z) = lmesh%GI3(:,ke,2)
        GsqrtV_z(:,ke_z) = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2D)

        GnnM_z(:,ke_z) = ( &
            1.0_RP / GsqrtV_z(:,ke_z)**2                                              &
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
            Dz, Lift, IntrpMat_VPOrdM1,                      & ! (in)
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
            G13_z, G23_z, GsqrtV_z, alph,                   & ! (in)
            Dz, Lift, IntrpMat_VPOrdM1,                     & ! (in)
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
          tend(:,:,ke_z) = ( - PROG_VARS(:,:,ke_z) + PROG_VARS00(:,:,ke_z) ) / impl_fac
        end do
      else
        call vi_eval_Ax( tend(:,:,:), alph,                & ! (out)
          PROG_VARS, PROG_VARS, DENS_hyd_z, PRES_hyd_z,    & ! (in)
          Dz, Lift, IntrpMat_VPOrdM1,                      & ! (in)
          GnnM_z, G13_z, G23_z, GsqrtV_z,                  & ! (in)
          modalFilterFlag, VModalFilter%FilterMat,         & ! (in)
          impl_fac, dt,                                    & ! (in) 
          lmesh, elem,                                     & ! (in)
          nz, vmapM, vmapP, ke_x, ke_y, .true. )             ! (in)
      end if

      !$omp parallel do private(ke)
      do ke_z=1, lmesh%NeZ
        ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
        DENS_dt(:,ke) = - tend(:,DENS_VID,ke_z)
        MOMX_dt(:,ke) = - tend(:,MOMX_VID,ke_z)
        MOMY_dt(:,ke) = - tend(:,MOMY_VID,ke_z)
        MOMZ_dt(:,ke) = - tend(:,MOMZ_VID,ke_z)
        RHOT_dt(:,ke) = - tend(:,RHOT_VID,ke_z)
      end do
      call PROF_rapend( 'hevi_cal_vi_retrun_var', 3)
    end do
    end do

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_hevi_cal_vi

end module scale_atm_dyn_dgm_globalnonhydro3d_hevi
