!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Regional nonhydrostatic model / HEVE
!!
!! @par Description
!!      HEVE DGM scheme for Atmospheric dynamical process. 
!!      The governing equations is a fully compressibile nonhydrostic equations, 
!!      which consist of mass, momentum, and thermodynamics (density * potential temperature conservation) equations. 
!!
!! @author Team SCALE
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
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend
#ifdef SCALE_PRODUCT_RUN_RM_MOUNTAIN_WAVE
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_set_dampcoef
#endif
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

#ifdef SCALE_PRODUCT_RUN_RM_MOUNTAIN_WAVE
  type(MeshField3D), public :: forcing_U0
  type(MeshField3D), public :: forcing_W0
  real(RP) :: U0
  real(RP) :: Grav_mod
  real(RP) :: ini_bg_force_tscale
  real(RP) :: ini_bg_force_turnoff_tstart
  real(RP) :: ini_bg_force_turnoff_tscale
  real(RP) :: ini_bg_sfac

  real(RP), allocatable :: sfac(:,:)
  real(RP), allocatable :: sfac_btm(:,:)
  real(RP) :: sw
#endif

contains
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_Init( mesh )
    implicit none
    class(MeshBase3D), target, intent(in) :: mesh

#ifdef SCALE_PRODUCT_RUN_RM_MOUNTAIN_WAVE
    integer :: n, ke
    class(ElementBase3D), pointer :: elem
    class(LocalMesh3D), pointer :: lmesh3D

    real(RP) :: zTop
    real(RP) :: SPONGE_HEIGHT
    real(RP) :: SPONGE_LATERAL_WIDTH
    real(RP) :: SPONGE_EFOLD_SEC
    real(RP) :: LATERAL_SPONGE_EFOLD_SEC
    real(RP) :: SL_TANH_NONDIM_WIDTH 

    real(RP) :: rtau_sponge
    real(RP) :: rtau_lateral_sponge
    real(RP) :: sponge_lateral_x00
    real(RP) :: sponge_lateral_x0    
#endif
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Init( mesh )
#ifdef SCALE_PRODUCT_RUN_RM_MOUNTAIN_WAVE
    call forcing_U0%Init( "forcing_U0", "m/s", mesh )
    call forcing_W0%Init( "forcing_W0", "m/s", mesh )
    do n=1, mesh%LOCAL_MESH_NUM
      lmesh3D => mesh%lcmesh_list(n)
      do ke=lmesh3D%NeS, lmesh3D%NeE
        forcing_U0%local(n)%val(:,ke) = 20.0_RP
        forcing_W0%local(n)%val(:,ke) = - lmesh3D%Gsqrt(:,ke) * lmesh3D%GI3(:,ke,1) * forcing_U0%local(n)%val(:,ke) &
          * exp(-lmesh3D%pos_en(:,ke,3)/2000.0_RP)
      end do
    end do

    !---
    U0 = 20.0_RP
    ini_bg_force_tscale = 60.0_RP
    ini_bg_force_turnoff_tstart = 120.0_RP
    ini_bg_force_turnoff_tscale = 1800.0_RP  
    
    zTop = 30E3_RP
    SPONGE_HEIGHT = 15E3_RP
    SPONGE_EFOLD_SEC = 100E0_RP
    SPONGE_LATERAL_WIDTH = 120E3_RP
    LATERAL_SPONGE_EFOLD_SEC = 100E0_RP
    SL_TANH_NONDIM_WIDTH = 0.16E0_RP

    lmesh3D => mesh%lcmesh_list(1)
    elem => lmesh3D%refElem3D
    allocate( sfac(elem%Np,lmesh3D%Ne) )
    allocate( sfac_btm(elem%Np,lmesh3D%Ne) )

    rtau_sponge = 1.0_RP / SPONGE_EFOLD_SEC
    rtau_lateral_sponge = 1.0_RP / LATERAL_SPONGE_EFOLD_SEC
    sponge_lateral_x00 = 240E3_RP - SPONGE_LATERAL_WIDTH
    sponge_lateral_x0 = sponge_lateral_x00 + 0.5_RP * SPONGE_LATERAL_WIDTH

    !$omp parallel do
    do ke=lmesh3D%NeS, lmesh3D%NeE
        sfac(:,ke) = &
            rtau_sponge * 0.5_RP * ( 1.0_RP + tanh( ( lmesh3D%zlev(:,ke) - 0.5_RP * ( zTop + SPONGE_HEIGHT ) ) / ( SL_TANH_NONDIM_WIDTH * ( zTop - SPONGE_HEIGHT ) ) ) ) &
          + rtau_lateral_sponge * 0.5_RP * ( 1.0_RP - tanh( ( lmesh3D%pos_en(:,ke,1) - 0.5_RP * SPONGE_LATERAL_WIDTH ) / ( SL_TANH_NONDIM_WIDTH * SPONGE_LATERAL_WIDTH ) ) ) &
          + rtau_lateral_sponge * 0.5_RP * ( 1.0_RP + tanh( ( lmesh3D%pos_en(:,ke,1) -             sponge_lateral_x0 ) / ( SL_TANH_NONDIM_WIDTH * SPONGE_LATERAL_WIDTH ) ) )

        sfac_btm(:,ke) = rtau_sponge * 0.5_RP * ( 1.0_RP + tanh( ( lmesh3D%zlev(:,ke) - 0.5_RP * ( zTop + SPONGE_HEIGHT ) ) / ( SL_TANH_NONDIM_WIDTH * ( zTop - SPONGE_HEIGHT ) ) ) )
    end do    
#endif
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_Init


  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_Final()
    implicit none
    !--------------------------------------------
    
    call atm_dyn_dgm_nonhydro3d_common_Final()
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_Final  

#ifdef SCALE_PRODUCT_RUN_RM_MOUNTAIN_WAVE
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_set_dampcoef( tsec )
    use scale_const, only: PI => CONST_PI
    implicit none
    real(RP), intent(in) :: tsec

    real(RP) :: ini_bg_off_tsec
    !-----------------------------------

    ini_bg_off_tsec = ini_bg_force_turnoff_tstart + ini_bg_force_turnoff_tscale
    if ( tsec < ini_bg_off_tsec ) then
      if ( tsec > ini_bg_force_turnoff_tstart ) then
        sw = 0.5_RP * ( 1.0_RP - cos( PI * ( ( tsec - ini_bg_force_turnoff_tstart ) / ini_bg_force_turnoff_tscale - 1.0_RP ) ) )
      else
        sw = 1.0_RP
      end if
    else
      sw = 0.0_RP
    end if
    ini_bg_sfac = sw * 1.0_RP / ini_bg_force_tscale
    
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_set_dampcoef
#endif

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                   & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, PRES_hyd_ref, & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot,                                                  & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )                     ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
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
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)

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
    
    real(RP) :: damp_coef(elem%Np)
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
    !$omp Fx, Fy, Fz, LiftDelFlx, damp_coef )
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

#ifdef SCALE_PRODUCT_RUN_RM_MOUNTAIN_WAVE
      damp_coef(:) =  ( 1.0_RP - sw ) * sfac(:,ke )
      DENS_dt(:,ke) = DENS_dt(:,ke) &
        - ini_bg_sfac * DDENS_(:,ke)
      MOMX_dt(:,ke) = MOMX_dt(:,ke) &
        -  ( ini_bg_sfac + damp_coef(:) ) * ( MOMX_(:,ke) - ( DENS_hyd(:,ke) + DDENS_(:,ke) ) * forcing_U0%local(lmesh%lcdomID)%val(:,ke) ) &
        + ( 1.0_RP + DDENS_(:,ke)/DENS_hyd(:,ke) ) * GradPhyd_x(:) * RGsqrtV(:)
      MOMY_dt(:,ke) = MOMY_dt(:,ke) &
        + ( 1.0_RP + DDENS_(:,ke)/DENS_hyd(:,ke) ) * GradPhyd_y(:) * RGsqrtV(:)
      MOMZ_dt(:,ke) = MOMZ_dt(:,ke) &
        - ini_bg_sfac * ( MOMZ_(:,ke) - ( DENS_hyd(:,ke) + DDENS_(:,ke) ) * forcing_W0%local(lmesh%lcdomID)%val(:,ke) ) &
        - damp_coef(:) * MOMZ_(:,ke)
      RHOT_dt(:,ke) = RHOT_dt(:,ke) &
        - ( ini_bg_sfac + damp_coef(:) ) * DRHOT_(:,ke)
      ! damp_coef(:) = sfac(:,ke )
      ! MOMX_dt(:,ke) = MOMX_dt(:,ke) &
      !   -  damp_coef(:) * ( MOMX_(:,ke) - ( DENS_hyd(:,ke) + DDENS_(:,ke) ) * forcing_U0%local(lmesh%lcdomID)%val(:,ke) ) &
      !   + ( 1.0_RP + DDENS_(:,ke)/DENS_hyd(:,ke) ) * GradPhyd_x(:) * RGsqrtV(:)
      ! MOMY_dt(:,ke) = MOMY_dt(:,ke) &
      !   - damp_coef(:) * MOMY_(:,ke) &
      !   + ( 1.0_RP + DDENS_(:,ke)/DENS_hyd(:,ke) ) * GradPhyd_y(:) * RGsqrtV(:)
      ! MOMZ_dt(:,ke) = MOMZ_dt(:,ke) &
      !   - damp_coef(:) * MOMZ_(:,ke)
      ! RHOT_dt(:,ke) = RHOT_dt(:,ke) &
      !   - damp_coef(:) * DRHOT_(:,ke)
#endif
    end do

    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend

end module scale_atm_dyn_dgm_nonhydro3d_rhot_heve
