!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics HEVI 
!!
!! @par Description
!!      HEVI DGM scheme for Atmospheric dynamical process. 
!!      To improve the numerical instability due to the aliasing errors, 
!!      the split form based on Gassner et al. (2016, JCP) is used for advection terms. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_splitform_hevi
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
  use scale_element_modalfilter, only: ModalFilter
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
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
  public :: atm_dyn_dgm_nonhydro3d_hevi_splitform_Init
  public :: atm_dyn_dgm_nonhydro3d_hevi_splitform_Final
  public :: atm_dyn_dgm_nonhydro3d_hevi_splitform_cal_tend
  public :: atm_dyn_dgm_nonhydro3d_hevi_splitform_cal_vi

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  real(RP), private, allocatable :: DxT1D_(:,:)
  real(RP), private, allocatable :: DyT1D_(:,:)
  real(RP), private, allocatable :: DzT1D_(:,:)

  private :: dx_ab, dy_ab, dz_ab
  private :: dx_abc, dy_abc, dz_abc

contains
  subroutine atm_dyn_dgm_nonhydro3d_hevi_splitform_Init( mesh )

    implicit none
    class(MeshBase3D), intent(in) :: mesh

    integer :: p1, p2, p3, p_
    type(ElementBase3D), pointer :: elem
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Init( mesh )

    elem => mesh%refElem3D

    allocate( DxT1D_(elem%Nnode_h1D,elem%Nnode_h1D) )
    allocate( DyT1D_(elem%Nnode_h1D,elem%Nnode_h1D) )
    allocate( DzT1D_(elem%Nnode_v,elem%Nnode_v) )

    do p1=1, elem%Nnode_h1D
      DxT1D_(:,p1) = elem%Dx1(p1,1:elem%Nnode_h1D)
    end do

    do p2=1, elem%Nnode_h1D
    do p_=1, elem%Nnode_h1D
      DyT1D_(p_,p2) = elem%Dx2(1+(p2-1)*elem%Nnode_h1D,1+(p_-1)*elem%Nnode_h1D)
    end do
    end do

    do p3=1, elem%Nnode_v
      DzT1D_(:,p3) = elem%Dx3(elem%Colmask(p3,1),elem%Colmask(:,1))
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_hevi_splitform_Init


  subroutine atm_dyn_dgm_nonhydro3d_hevi_splitform_Final()
    implicit none
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Final()

    deallocate( DxT1D_, DyT1D_, DzT1D_ )

    return
  end subroutine atm_dyn_dgm_nonhydro3d_hevi_splitform_Final  

  !-------------------------------

  subroutine atm_dyn_dgm_nonhydro3d_hevi_splitform_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, CORIOLIS,          & ! (in)
    SL_flag, wdamp_tau, wdamp_height, hveldamp_flag,                            & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )

    use scale_atm_dyn_dgm_nonhydro3d_hevi_numflux, only: &
      atm_dyn_dgm_nonhydro3d_hevi_numflux_get_generalvc

    use scale_atm_dyn_dgm_spongelayer, only: &
      atm_dyn_dgm_spongelayer_add_tend
    
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
    real(RP) :: Fx_sp(elem%Np), Fy_sp(elem%Np), Fz_sp(elem%Np)  
    real(RP) :: GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: GsqrtDens_(elem%Np), rdens_(elem%Np), RHOT_hyd(elem%Np), RHOT_(elem%Np)
    real(RP) :: dpres_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%Np), pot_(elem%Np)
    real(RP) :: Cori(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np)

    integer :: ke, ke2d

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR       
    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 3)
    call atm_dyn_dgm_nonhydro3d_hevi_numflux_get_generalvc( &
      del_flux, del_flux_hyd,                                                 & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                & ! (in)
      lmesh%Gsqrt, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),                        & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, lmesh2D, elem2D )                                            ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( ke2d, Cori,            &
    !$omp RHOT_, DPRES_, GsqrtDens_, rdens_, u_, v_, w_, wt_,     &
    !$omp GradPhyd_x, GradPhyd_y,                     &
    !$omp GsqrtV, RGsqrtV,                            &
    !$omp Fx, Fy, Fz, Fx_sp, Fy_sp, Fz_sp, LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      Cori(:) = CORIOLIS(elem%IndexH2Dto3D(:),ke2d)

      GsqrtV(:)  = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2d)
      RGsqrtV(:) = 1.0_RP / GsqrtV(:)

      !--
      RHOT_(:) = P0ovR * (PRES_hyd(:,ke) * rP0)**rgamm + DRHOT_(:,ke)
      DPRES_(:) = PRES00 * ( RovP0 * RHOT_(:) )**gamm &
                - PRES_hyd(:,ke)

      GsqrtDens_(:) = lmesh%Gsqrt(:,ke) * ( DDENS_(:,ke) + DENS_hyd(:,ke) )
      rdens_(:) = 1.0_RP / GsqrtDens_(:)
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)
      wt_(:) = w_(:) * RGsqrtV(:) + lmesh%GI3(:,ke,1) * u_(:) + lmesh%GI3(:,ke,2) * v_(:) 
      pot_(:) = RHOT_(:) * rdens_(:)

      ke2d = lmesh%EMap3Dto2D(ke)
      Cori(:) = CORIOLIS(elem%IndexH2Dto3D(:),ke2d)

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
      call dx_ab( DxT1D_, GsqrtDens_(:), u_, elem%Nnode_h1D, elem%Nnode_v, Fx_sp )
      call dy_ab( DyT1D_, GsqrtDens_(:), v_, elem%Nnode_h1D, elem%Nnode_v, Fy_sp )
      call dz_ab( DzT1D_, GsqrtDens_(:), wt_(:) - w_(:) * RGsqrtV(:), elem%Nnode_h1D, elem%Nnode_v, Fz_sp )
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DENS_VID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx_sp(:) &
          + lmesh%Escale(:,ke,2,2) * Fy_sp(:) &
          + lmesh%Escale(:,ke,3,3) * Fz_sp(:) &
          + LiftDelFlx(:) )
      
      !-- MOMX
      call dx_abc( DxT1D_, GsqrtDens_, u_,  u_, elem%Nnode_h1D, elem%Nnode_v, Fx_sp )
      call dy_abc( DyT1D_, GsqrtDens_, u_,  v_, elem%Nnode_h1D, elem%Nnode_v, Fy_sp )
      call dz_abc( DzT1D_, GsqrtDens_, u_, wt_, elem%Nnode_h1D, elem%Nnode_v, Fz_sp )
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * dpres_(:)               , Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * ( Fx_sp(:) + Fx(:) )  &
            + lmesh%Escale(:,ke,2,2) * Fy_sp(:)              &
            + lmesh%Escale(:,ke,3,3) * Fz_sp(:)              &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)            &
          - GradPhyd_x(:) * RGsqrtV(:)                       &
          + Cori(:) * MOMY_(:,ke)

      !-- MOMY
      call dx_abc( DxT1D_, GsqrtDens_, v_,  u_, elem%Nnode_h1D, elem%Nnode_v, Fx_sp )
      call dy_abc( DyT1D_, GsqrtDens_, v_,  v_, elem%Nnode_h1D, elem%Nnode_v, Fy_sp )  
      call dz_abc( DzT1D_, GsqrtDens_, v_, wt_, elem%Nnode_h1D, elem%Nnode_v, Fz_sp )
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * dpres_(:)               , Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx_sp(:)             &
            + lmesh%Escale(:,ke,2,2) * ( Fy_sp(:) + Fy(:) ) &
            + lmesh%Escale(:,ke,3,3) * Fz_sp(:)             &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)           &
          - GradPhyd_y(:) * RGsqrtV(:)                      &
          - Cori(:) * MOMX_(:,ke)

      !-- MOMZ
      call dx_abc( DxT1D_, GsqrtDens_, w_,  u_, elem%Nnode_h1D, elem%Nnode_v, Fx_sp )
      call dy_abc( DyT1D_, GsqrtDens_, w_,  v_, elem%Nnode_h1D, elem%Nnode_v, Fy_sp )
      call dz_abc( DzT1D_, GsqrtDens_, w_, wt_, elem%Nnode_h1D, elem%Nnode_v, Fz_sp )
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx_sp(:)   &
          + lmesh%Escale(:,ke,2,2) * Fy_sp(:)   &
          + lmesh%Escale(:,ke,3,3) * Fz_sp(:)   &
          + LiftDelFlx(:)                       )   

      !-- RHOT
      call dx_abc( DxT1D_, GsqrtDens_, pot_,  u_, elem%Nnode_h1D, elem%Nnode_v, Fx_sp )
      call dy_abc( DyT1D_, GsqrtDens_, pot_,  v_, elem%Nnode_h1D, elem%Nnode_v, Fy_sp )        
      call dz_abc( DzT1D_, GsqrtDens_, pot_,  wt_(:) - w_(:) * RGsqrtV(:), elem%Nnode_h1D, elem%Nnode_v, Fz_sp )        
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,RHOT_VID), LiftDelFlx)
      
      RHOT_dt(:,ke) =  - (  &
            lmesh%Escale(:,ke,1,1) * Fx_sp(:) &
          + lmesh%Escale(:,ke,2,2) * Fy_sp(:) &
          + lmesh%Escale(:,ke,3,3) * Fz_sp(:) &
          + LiftDelFlx(:)                     ) 
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 3)

    !- Sponge layer
    if (SL_flag) then
      call PROF_rapstart( 'cal_dyn_tend_sponge', 3)
      call atm_dyn_dgm_spongelayer_add_tend( &
        MOMX_dt, MOMY_dt, MOMZ_dt,                    & ! (out)
        MOMX_, MOMY_, MOMZ_, wdamp_tau, wdamp_height, & ! (in)
        hveldamp_flag, lmesh, elem                    ) ! (in)
      call PROF_rapend( 'cal_dyn_tend_sponge', 3)
    end if

    return
  end subroutine atm_dyn_dgm_nonhydro3d_hevi_splitform_cal_tend

  !------

  subroutine atm_dyn_dgm_nonhydro3d_hevi_splitform_cal_vi( &
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

        PROG_VARS0(:,:,ke_z) = PROG_VARS(:,:,ke_z)
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
        MOMX_dt(:,ke) = - tend(:,MOMX_VID ,ke_z)
        MOMY_dt(:,ke) = - tend(:,MOMY_VID ,ke_z)
        MOMZ_dt(:,ke) = - tend(:,MOMZ_VID ,ke_z)
        RHOT_dt(:,ke) = - tend(:,RHOT_VID,ke_z)
      end do
      call PROF_rapend( 'hevi_cal_vi_retrun_var', 3)
    end do
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_hevi_splitform_cal_vi

  !-----------------------------------------

  subroutine dx_ab(DxT1D, a, b, Nnode_h1D, Nnode_v, fx)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DxT1D(Nnode_h1D,Nnode_h1D)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fx(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fx(i,j,k) = 0.5_RP * sum( DxT1D(:,i) * (a(i,j,k) + a(:,j,k)) * (b(i,j,k) + b(:,j,k)) )
    end do
    end do
    end do

    return
  end subroutine dx_ab

  subroutine dy_ab(DyT1D, a, b, Nnode_h1D, Nnode_v, fy)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DyT1D(Nnode_h1D,Nnode_h1D)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fy(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fy(i,j,k) = 0.5_RP * sum( DyT1D(:,j) * (a(i,j,k) + a(i,:,k)) * (b(i,j,k) + b(i,:,k)) )
    end do
    end do
    end do

    return
  end subroutine dy_ab
  
  subroutine dz_ab(DzT1D, a, b, Nnode_h1D, Nnode_v, fz)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DzT1D(Nnode_v,Nnode_v)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fz(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fz(i,j,k) = 0.5_RP * sum( DzT1D(:,k) * (a(i,j,k) + a(i,j,:)) * (b(i,j,k) + b(i,j,:)) )
    end do
    end do
    end do

    return
  end subroutine dz_ab

  subroutine dx_abc(DxT1D, a, b, c, Nnode_h1D, Nnode_v, fx)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DxT1D(Nnode_h1D,Nnode_h1D)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: c(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fx(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fx(i,j,k) = 0.25_RP * sum( DxT1D(:,i) * (a(i,j,k) + a(:,j,k)) * (b(i,j,k) + b(:,j,k)) * (c(i,j,k) + c(:,j,k)) )
    end do
    end do
    end do

    return
  end subroutine dx_abc

  subroutine dy_abc(DyT1D, a, b, c, Nnode_h1D, Nnode_v, fy)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DyT1D(Nnode_h1D,Nnode_h1D)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: c(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fy(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fy(i,j,k) = 0.25_RP * sum( DyT1D(:,j) * (a(i,j,k) + a(i,:,k)) * (b(i,j,k) + b(i,:,k)) * (c(i,j,k) + c(i,:,k)) )
    end do
    end do
    end do

    return
  end subroutine dy_abc

  subroutine dz_abc(DzT1D, a, b, c, Nnode_h1D, Nnode_v, fz)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DzT1D(Nnode_v,Nnode_v)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: c(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fz(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fz(i,j,k) = 0.25_RP * sum( DzT1D(:,k) * (a(i,j,k) + a(i,j,:)) * (b(i,j,k) + b(i,j,:)) * (c(i,j,k) + c(i,j,:)) )
    end do
    end do
    end do

    return
  end subroutine dz_abc

end module scale_atm_dyn_dgm_nonhydro3d_splitform_hevi
