!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Global nonhydrostatic model / HEVI
!!
!! @par Description
!!      HEVI DGM scheme for Global Atmospheric Dynamical process. 
!!      The governing equations is a fully compressibile nonhydrostic equations, 
!!      which consist of mass, momentum, and thermodynamics (total energy conservation) equations. 
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_globalnonhydro3d_etot_hevi
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
    DENS_VID => PRGVAR_DDENS_ID, ETOT_VID => PRGVAR_ETOT_ID,  &
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
  public :: atm_dyn_dgm_globalnonhydro3d_etot_hevi_Init
  public :: atm_dyn_dgm_globalnonhydro3d_etot_hevi_Final
  public :: atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_tend
  public :: atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_vi

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
  subroutine atm_dyn_dgm_globalnonhydro3d_etot_hevi_Init( mesh )
    implicit none
    class(MeshBase3D), intent(in) :: mesh
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Init( mesh )
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_etot_hevi_Init

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_etot_hevi_Final()
    implicit none
    !--------------------------------------------
    
    call atm_dyn_dgm_nonhydro3d_common_Final()    
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_etot_hevi_Final  

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, EnTot_dt,                                   & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, ETOT_, DPRES_, DENS_hyd, PRES_hyd, PRES_hyd_ref,   & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot,                                                   & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )   ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_etot_hevi_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_etot_hevi_numflux_get_generalhvc
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
    real(RP), intent(out) :: EnTot_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: ETOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DPRES_(elem%Np,lmesh%NeA)    
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd_ref(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot (elem%Np,lmesh%NeA)    
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DPRES_hyd(elem%Np), GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%np)

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

    real(RP) :: entot_(elem%Np), enthalpy_(elem%Np)
    real(RP) :: u1_(elem%Np), u2_(elem%Np)
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux, del_flux_hyd,                                                  & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, ETOT_, DPRES_, DENS_hyd, PRES_hyd,          & ! (in)
      Rtot, CVtot, CPtot,                                                      & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%G_ij(:,:,1,1), lmesh%G_ij(:,:,1,2), lmesh%G_ij(:,:,2,2),           & ! (in)
      lmesh%GsqrtH, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2), lmesh%zlev(:,:),       & ! (in)
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
    !$omp rdens_, u_, v_, w_, wt_,                 &
    !$omp entot_, enthalpy_, u1_, u2_,             &
    !$omp Fx, Fy, Fz, LiftDelFlx,                  &
    !$omp DPRES_hyd, GradPhyd_x, GradPhyd_y,       &
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
      rdens_(:) = 1.0_RP / ( DDENS_(:,ke) + DENS_hyd(:,ke) )
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)
      wt_(:) = w_(:) * RGsqrtV(:) + lmesh%GI3(:,ke,1) * u_(:) + lmesh%GI3(:,ke,2) * v_(:) 
      u1_(:) = lmesh%G_ij(elem%IndexH2Dto3D(:),ke2d,1,1) * u_(:) + lmesh%G_ij(elem%IndexH2Dto3D(:),ke2d,2,1) * v_(:)
      u2_(:) = lmesh%G_ij(elem%IndexH2Dto3D(:),ke2d,2,1) * u_(:) + lmesh%G_ij(elem%IndexH2Dto3D(:),ke2d,2,2) * v_(:)

      ! DPRES_(:) = ( CPtot(:,ke) / CVtot(:,ke) - 1.0_RP ) &
      !           * (   ETOT_(:,ke) - Grav * ( DDENS_(:,ke) + DENS_hyd(:,ke) ) * lmesh%zlev(:,ke)        &
      !               - 0.5_RP * ( MOMX_(:,ke) * u1_(:) + MOMY_(:,ke) * u2_(:) + MOMZ_(:,ke) * w_(:) ) ) &
      !           - PRES_hyd(:,ke)
      
      enthalpy_(:) = ETOT_(:,ke) + PRES_hyd(:,ke) + DPRES_(:,ke)

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
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DENS_VID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)     &
          + lmesh%Escale(:,ke,2,2) * Fy(:)     &
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
          + twoOVdel2(:) * Y(:) *                                                    &
            ( - X(:) * Y(:) * u_(:) + (1.0_RP + Y(:)**2) * v_(:) ) * MOMX_(:,ke)     &
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
            + twoOVdel2(:) * X(:) *                                                    &
              ( (1.0_RP + X(:)**2) * u_(:) - X(:) * Y(:) * v_(:) ) * MOMY_(:,ke)       &
            - ( G12(:) * GradPhyd_x(:) + G22(:) * GradPhyd_y(:) ) * RGsqrtV(:)         &
            + CORI(:,2) 

      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *   u_ (:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *   v_ (:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) *   wt_(:) * MOMZ_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)       &
            + lmesh%Escale(:,ke,2,2) * Fy(:)       &
            + lmesh%Escale(:,ke,3,3) * Fz(:)       &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)

      !-- EnTot
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * u_ (:) * enthalpy_(:), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * v_ (:) * enthalpy_(:), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,ETOT_VID), LiftDelFlx)
      
      EnTot_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 

    end do
    !$omp end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_tend

!OCL SERIAL  
  subroutine atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_vi( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, ETOT_dt,             & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, ETOT_, DENS_hyd, PRES_hyd,  & ! (in)
    DDENS0_, MOMX0_, MOMY0_, MOMZ0_, ETOT0_,                 & ! (in)
    Rtot, CVtot, CPtot,                                      & ! (in)
    element3D_operation, Dz, Lift,                           & ! (in)
    modalFilterFlag, VModalFilter,                           & ! (in)
    impl_fac, dt,                                            & ! (in)
    lmesh, elem, lmesh2D, elem2D                             ) ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_etot_hevi_common, only: &
      vi_gen_vmap => atm_dyn_dgm_nonhydro3d_etot_hevi_common_gen_vmap,                  &
      vi_eval_Ax => atm_dyn_dgm_nonhydro3d_etot_hevi_common_eval_Ax,                    &
      vi_eval_Ax_uv => atm_dyn_dgm_nonhydro3d_etot_hevi_common_eval_Ax_uv,              &
      vi_construct_matbnd => atm_dyn_dgm_nonhydro3d_etot_hevi_common_construct_matbnd,  &
      vi_construct_matbnd_uv => atm_dyn_dgm_nonhydro3d_etot_hevi_common_construct_matbnd_uv
      
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: ETOT_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: ETOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS0_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX0_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY0_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ0_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: ETOT0_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    class(ElementOperationBase3D), intent(in) :: element3D_operation
    class(SparseMat), intent(in) :: Dz, Lift
    logical, intent(in) :: modalFilterFlag
    class(ModalFilter), intent(in) :: VModalFilter
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt

    real(RP) :: PROG_VARS (elem%Np,lmesh%NeZ,PRGVAR_NUM,lmesh%NeX*lmesh%NeY)
    real(RP) :: DPRES     (elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: PROG_VARS0(elem%Np,lmesh%NeZ,PRGVAR_NUM,lmesh%NeX*lmesh%NeY)
    real(RP) :: DPRES0    (elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: b1D(3,elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2,lmesh%NeX*lmesh%NeY)
    real(RP) :: GeoPot    (elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: KinHovDENS(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    integer :: ipiv(elem%Nnode_v*3*lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: b1D_uv(elem%Nnode_v,lmesh%NeZ,2,elem%Nnode_h1D**2,lmesh%NeX*lmesh%NeY)
    integer :: ipiv_uv(elem%Nnode_v*1*lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: alph(elem%NfpTot,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: Rtot_z(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: CPtot_ov_CVtot(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: DENS_hyd_z(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: PRES_hyd_z(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: GnnM_z(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: G13_z(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: G23_z(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: GsqrtV_z(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    integer :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer :: vmapP(elem%NfpTot,lmesh%NeZ)
    integer :: ColMask(elem%Nnode_v)
    integer :: ke_xy, ke_z, ke, ke2D, v
    integer :: itr_nlin
    integer :: kl, ku, nz_1D
    integer :: kl_uv, ku_uv, nz_1D_uv
    integer :: ij
    logical :: is_converged

    real(RP), allocatable :: PmatBnd(:,:,:)
    real(RP), allocatable :: PmatBnd_uv(:,:,:)
    integer :: info, info_uv

    real(RP) :: DENS_(elem%Np)
    !------------------------------------------------------------------------

    call PROF_rapstart( 'hevi_cal_vi_prep', 3)

    nz_1D = elem%Nnode_v * 3 * lmesh%NeZ
    kl = ( elem%Nnode_v + 1 ) * 3 - 1
    ku = kl
    nz_1D_uv = elem%Nnode_v * 1 * lmesh%NeZ
    kl_uv = elem%Nnode_v
    ku_uv = kl_uv
    allocate( PmatBnd   (2*kl+ku+1,nz_1D,elem%Nnode_h1D**2) )
    allocate( PmatBnd_uv(2*kl_uv+ku_uv+1,nz_1D_uv,elem%Nnode_h1D**2) )

    call vi_gen_vmap( vmapM, vmapP, & ! (out)
      lmesh, elem                   ) ! (in)

    !-
    
    !$omp parallel private( ke_xy, ke_z, ke, ke2D, DENS_ )
    !$omp do collapse(2)
    do ke_xy=1, lmesh%NeX*lmesh%NeY
    do ke_z=1, lmesh%NeZ
      ke = ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
      ke2D = lmesh%EMap3Dto2D(ke)

      PROG_VARS(:,ke_z,DENS_VID,ke_xy) = DDENS0_(:,ke)
      PROG_VARS(:,ke_z,MOMX_VID,ke_xy) = MOMX0_ (:,ke)
      PROG_VARS(:,ke_z,MOMY_VID,ke_xy) = MOMY0_ (:,ke)
      PROG_VARS(:,ke_z,MOMZ_VID,ke_xy) = MOMZ0_ (:,ke)
      PROG_VARS(:,ke_z,ETOT_VID,ke_xy) = ETOT0_ (:,ke)

      DENS_hyd_z(:,ke_z,ke_xy) = DENS_hyd(:,ke)
      PRES_hyd_z(:,ke_z,ke_xy) = PRES_hyd(:,ke)
      GeoPot    (:,ke_z,ke_xy) = GRAV * lmesh%zlev(:,ke)

      DENS_(:) = DENS_hyd(:,ke) + DDENS0_(:,ke)
      KinHovDENS(:,ke_z,ke_xy) = 0.5_RP * ( &
          MOMX0_(:,ke) * ( lmesh%G_ij(elem%IndexH2Dto3D,ke2D,1,1) * MOMX0_(:,ke) + lmesh%G_ij(elem%IndexH2Dto3D,ke2D,2,1) * MOMY0_(:,ke) ) &
        + MOMY0_(:,ke) * ( lmesh%G_ij(elem%IndexH2Dto3D,ke2D,2,1) * MOMX0_(:,ke) + lmesh%G_ij(elem%IndexH2Dto3D,ke2D,2,2) * MOMY0_(:,ke) ) &
        ) / DENS_(:)**2
        
      Rtot_z(:,ke_z,ke_xy) = Rtot(:,ke)
      CPtot_ov_CVtot(:,ke_z,ke_xy) = CPtot(:,ke) / CVtot(:,ke)
      
      DPRES(:,ke_z,ke_xy) = &
        ( CPtot_ov_CVtot(:,ke_z,ke_xy) - 1.0_RP ) &
         * ( ETOT0_(:,ke) - ( DENS_(:) * ( KinHovDENS(:,ke_z,ke_xy) + GeoPot(:,ke_z,ke_xy) ) + 0.5_RP * MOMZ0_(:,ke)**2 / DENS_(:) ) ) &
        - PRES_hyd(:,ke)

      nz(:,ke_z,ke_xy) = lmesh%normal_fn(:,ke,3)
      G13_z   (:,ke_z,ke_xy) = lmesh%GI3(:,ke,1)
      G23_z   (:,ke_z,ke_xy) = lmesh%GI3(:,ke,2)
      GsqrtV_z(:,ke_z,ke_xy) = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2D)

      GnnM_z(:,ke_z,ke_xy) = ( &
          1.0_RP / GsqrtV_z(:,ke_z,ke_xy)**2                                                    &
        + G13_z(:,ke_z,ke_xy) * ( lmesh%GIJ(elem%IndexH2Dto3D,ke2D,1,1) * G13_z(:,ke_z,ke_xy)   &
                                + lmesh%GIJ(elem%IndexH2Dto3D,ke2D,1,2) * G23_z(:,ke_z,ke_xy) ) &
        + G23_z(:,ke_z,ke_xy) * ( lmesh%GIJ(elem%IndexH2Dto3D,ke2D,1,2) * G13_z(:,ke_z,ke_xy)   &
                                + lmesh%GIJ(elem%IndexH2Dto3D,ke2D,2,2) * G23_z(:,ke_z,ke_xy) ) )
    end do
    end do
    !$omp end do
    !$omp workshare
    PROG_VARS0 (:,:,:,:) = PROG_VARS(:,:,:,:)
    DPRES0(:,:,:) = DPRES(:,:,:)
    !$omp end workshare
    !$omp end parallel
      
    call PROF_rapend( 'hevi_cal_vi_prep', 3)

    !--

    if ( abs(impl_fac) > 0.0_RP ) then
      call PROF_rapstart( 'hevi_cal_vi_itr', 3)

      ! G = (q^n+1 - q^n*) + impl_fac * A(q^n+1) = 0
      ! dG/dq^n+1 del[q] = - G(q^n*)
      do itr_nlin = 1, 1
        call PROF_rapstart( 'hevi_cal_vi_ax', 3)

        call vi_eval_Ax_uv( &
          MOMX_dt(:,:), MOMY_dt(:,:), alph(:,:,:),                              & ! (out, dummy)
          PROG_VARS, DPRES, PROG_VARS0, DPRES0,                                 & ! (in)
          DDENS_, MOMX_, MOMY_, MOMZ_, ETOT_,                                   & ! (in)
          DENS_hyd_z, PRES_hyd_z,                                               & ! (in)
          Rtot_z, CPtot_ov_CVtot,                                               & ! (in)
          Dz, Lift, IntrpMat_VPOrdM1,                                           & ! (in)
          GnnM_z, G13_z, G23_z, GsqrtV_z,                                       & ! (in)
          modalFilterFlag, VModalFilter%FilterMat,                              & ! (in)
          impl_fac, dt,                                                         & ! (in) 
          lmesh, elem, nz, vmapM, vmapP,                                        & ! (in)
          b1D_uv(:,:,:,:,:)                                                     ) ! (out)

        call PROF_rapend( 'hevi_cal_vi_ax', 3)

        do ke_xy=1, lmesh%NeX * lmesh%NeY
          call PROF_rapstart( 'hevi_cal_vi_matbnd', 3)
          call vi_construct_matbnd_uv( PmatBnd_uv(:,:,:),              & ! (out)
            kl_uv, ku_uv, nz_1D_uv,                                    & ! (in)
            PROG_VARS(:,:,:,ke_xy), KinHovDENS(:,:,ke_xy),             & ! (in)
            DENS_hyd_z(:,:,ke_xy), PRES_hyd_z(:,:,ke_xy),              & ! (in)
            G13_z(:,:,ke_xy), G23_z(:,:,ke_xy), GsqrtV_z(:,:,ke_xy),   & ! (in)
            alph(:,:,ke_xy),                                           & ! (in)
            Rtot_z(:,:,ke_xy), CPtot_ov_CVtot(:,:,ke_xy),              & ! (in)
            GeoPot(:,:,ke_xy),                                         & ! (in)
            Dz, Lift, IntrpMat_VPOrdM1,                                & ! (in)
            modalFilterFlag, VModalFilter%FilterMat,                   & ! (in)
            impl_fac, dt,                                              & ! (in)
            lmesh, elem, nz(:,:,ke_xy), vmapM, vmapP, ke_xy, 1         ) ! (in)
          call PROF_rapend( 'hevi_cal_vi_matbnd', 3)

          call PROF_rapstart( 'hevi_cal_vi_lin', 3)
          !$omp parallel private(ij, v, ke_z, info, ColMask)
          !$omp do
          do ij=1, elem%Nnode_h1D**2
            call dgbsv( nz_1D_uv, kl_uv, ku_uv, 2, PmatBnd_uv(:,:,ij), 2*kl_uv+ku_uv+1, ipiv_uv(:,ij), b1D_uv(:,:,:,ij,ke_xy), nz_1D_uv, info)

            ColMask(:) = elem%Colmask(:,ij)
            do ke_z=1, lmesh%NeZ
              PROG_VARS(ColMask(:),ke_z,MOMX_VID,ke_xy) = PROG_VARS(Colmask(:),ke_z,MOMX_VID,ke_xy) + b1D_uv(:,ke_z,1,ij,ke_xy)
              PROG_VARS(ColMask(:),ke_z,MOMY_VID,ke_xy) = PROG_VARS(Colmask(:),ke_z,MOMY_VID,ke_xy) + b1D_uv(:,ke_z,2,ij,ke_xy)
            end do
          end do ! for ij
          !$omp end do
          !$omp end parallel
          call PROF_rapend( 'hevi_cal_vi_lin', 3)

        end do ! for ke_xy
          
        call PROF_rapstart( 'hevi_cal_vi_ax', 3)
        call vi_eval_Ax( &
          DENS_dt(:,:), MOMZ_dt(:,:), ETOT_dt(:,:),                             & ! (out, dummy) 
          alph(:,:,:),                                                          & ! (in)
          PROG_VARS, DPRES, PROG_VARS0, DPRES0,                                 & ! (in)
          DDENS_, MOMX_, MOMY_, MOMZ_, ETOT_,                                   & ! (in)
          DENS_hyd_z, PRES_hyd_z,                                               & ! (in)
          Rtot_z, CPtot_ov_CVtot,                                               & ! (in)
          Dz, Lift, IntrpMat_VPOrdM1,                                           & ! (in)
          GnnM_z, G13_z, G23_z, GsqrtV_z,                                       & ! (in)
          modalFilterFlag, VModalFilter%FilterMat,                              & ! (in)
          impl_fac, dt,                                                         & ! (in) 
          lmesh, elem, nz, vmapM, vmapP,                                        & ! (in)
          b1D(:,:,:,:,:)                                                        ) ! (out)
        call PROF_rapend( 'hevi_cal_vi_ax', 3)

        do ke_xy=1, lmesh%NeX * lmesh%NeY
          call PROF_rapstart( 'hevi_cal_vi_matbnd', 3)
          call vi_construct_matbnd( PmatBnd(:,:,:),                    & ! (out)
            kl, ku, nz_1D,                                             & ! (in)
            PROG_VARS(:,:,:,ke_xy), KinHovDENS(:,:,ke_xy),             & ! (in)
            DENS_hyd_z(:,:,ke_xy), PRES_hyd_z(:,:,ke_xy),              & ! (in)
            G13_z(:,:,ke_xy), G23_z(:,:,ke_xy), GsqrtV_z(:,:,ke_xy),   & ! (in)
            alph(:,:,ke_xy),                                           & ! (in)
            Rtot_z(:,:,ke_xy), CPtot_ov_CVtot(:,:,ke_xy),              & ! (in)
            GeoPot(:,:,ke_xy),                                         & ! (in)
            Dz, Lift, IntrpMat_VPOrdM1,                                & ! (in)
            modalFilterFlag, VModalFilter%FilterMat,                   & ! (in)
            impl_fac, dt,                                              & ! (in)
            lmesh, elem, nz(:,:,ke_xy), vmapM, vmapP, ke_xy, 1         ) ! (in)
          call PROF_rapend( 'hevi_cal_vi_matbnd', 3)

          call PROF_rapstart( 'hevi_cal_vi_lin', 3)
          !$omp parallel private(ij, v, ke_z, info, ColMask)
          !$omp do
          do ij=1, elem%Nnode_h1D**2
            call dgbsv( nz_1D, kl, ku, 1, PmatBnd(:,:,ij), 2*kl+ku+1, ipiv(:,ij), b1D(:,:,:,ij,ke_xy), nz_1D, info)

            ColMask(:) = elem%Colmask(:,ij)
            do ke_z=1, lmesh%NeZ
              PROG_VARS(ColMask(:),ke_z,DENS_VID,ke_xy) = PROG_VARS(Colmask(:),ke_z,DENS_VID,ke_xy) + b1D(1,:,ke_z,ij,ke_xy)
              PROG_VARS(ColMask(:),ke_z,MOMZ_VID,ke_xy) = PROG_VARS(Colmask(:),ke_z,MOMZ_VID,ke_xy) + b1D(2,:,ke_z,ij,ke_xy)
              PROG_VARS(ColMask(:),ke_z,ETOT_VID,ke_xy) = PROG_VARS(Colmask(:),ke_z,ETOT_VID,ke_xy) + b1D(3,:,ke_z,ij,ke_xy)
            end do
          end do ! for ij
          !$omp end do
          !$omp do
          do ke_z=1, lmesh%NeZ
            DENS_(:) = DENS_hyd_z(:,ke_z,ke_xy) + PROG_VARS(:,ke_z,DENS_VID,ke_xy)
            DPRES(:,ke_z,ke_xy) = &
              ( CPtot_ov_CVtot(:,ke_z,ke_xy) - 1.0_RP ) &
              * ( PROG_VARS(:,ke_z,ETOT_VID,ke_xy)                                                                                               &
                - ( DENS_(:) * ( KinHovDENS(:,ke_z,ke_xy) + GeoPot(:,ke_z,ke_xy) ) + 0.5_RP * PROG_VARS(:,ke_z,MOMZ_VID,ke_xy)**2 / DENS_(:) ) ) &
            - PRES_hyd_z(:,ke_z,ke_xy)
          end do
          !$omp end parallel
          call PROF_rapend( 'hevi_cal_vi_lin', 3)

        end do ! for ke_xy
      end do ! itr nlin

      call PROF_rapend( 'hevi_cal_vi_itr', 3)
    end if

    call PROF_rapstart( 'hevi_cal_vi_retrun_var', 3)
    if ( abs(impl_fac) > 0.0_RP) then
      !$omp parallel do collapse(2) private(ke_xy, ke_z, ke)
      do ke_xy=1, lmesh%NeX * lmesh%NeY
      do ke_z=1, lmesh%NeZ
        ke = ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
        DENS_dt(:,ke) = ( PROG_VARS(:,ke_z,DENS_VID,ke_xy) - DDENS_(:,ke) ) / impl_fac
        MOMX_dt(:,ke) = ( PROG_VARS(:,ke_z,MOMX_VID,ke_xy) - MOMX_ (:,ke) ) / impl_fac
        MOMY_dt(:,ke) = ( PROG_VARS(:,ke_z,MOMY_VID,ke_xy) - MOMY_ (:,ke) ) / impl_fac
        MOMZ_dt(:,ke) = ( PROG_VARS(:,ke_z,MOMZ_VID,ke_xy) - MOMZ_ (:,ke) ) / impl_fac
        ETOT_dt(:,ke) = ( PROG_VARS(:,ke_z,ETOT_VID,ke_xy) - ETOT_ (:,ke) ) / impl_fac
      end do
      end do
    else
      call vi_eval_Ax_uv( &
        MOMX_dt(:,:), MOMY_dt(:,:), alph(:,:,:),                              & ! (out)
        PROG_VARS, DPRES, PROG_VARS0, DPRES0,                                 & ! (in)
        DDENS_, MOMX_, MOMY_, MOMZ_, ETOT_,                                   & ! (in)
        DENS_hyd_z, PRES_hyd_z,                                               & ! (in)
        Rtot_z, CPtot_ov_CVtot,                                               & ! (in)
        Dz, Lift, IntrpMat_VPOrdM1,                                           & ! (in)
        GnnM_z, G13_z, G23_z, GsqrtV_z,                                       & ! (in)
        modalFilterFlag, VModalFilter%FilterMat,                              & ! (in)
        impl_fac, dt,                                                         & ! (in) 
        lmesh, elem, nz, vmapM, vmapP                                         ) ! (in)
          
      call vi_eval_Ax( &
        DENS_dt(:,:), MOMZ_dt(:,:), ETOT_dt(:,:),                             & ! (out) 
        alph(:,:,:),                                                          & ! (in)
        PROG_VARS, DPRES, PROG_VARS0, DPRES0,                                 & ! (in)
        DDENS_, MOMX_, MOMY_, MOMZ_, ETOT_,                                   & ! (in)
        DENS_hyd_z, PRES_hyd_z,                                               & ! (in)
        Rtot_z, CPtot_ov_CVtot,                                               & ! (in)
        Dz, Lift, IntrpMat_VPOrdM1,                                           & ! (in)
        GnnM_z, G13_z, G23_z, GsqrtV_z,                                       & ! (in)
        modalFilterFlag, VModalFilter%FilterMat,                              & ! (in)
        impl_fac, dt,                                                         & ! (in) 
        lmesh, elem, nz, vmapM, vmapP                                         ) ! (in)
    end if
    call PROF_rapend( 'hevi_cal_vi_retrun_var', 3)


    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_vi

end module scale_atm_dyn_dgm_globalnonhydro3d_etot_hevi
