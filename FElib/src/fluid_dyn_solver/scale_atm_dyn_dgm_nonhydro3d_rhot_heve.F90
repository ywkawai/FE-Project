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

  use scale_element_line, only: LineElement

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
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_new

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  real(RP), allocatable :: D1D(:,:), D1Dt(:,:)
  real(RP), allocatable :: IntrpMat_VPOrdM1_tr(:,:)

contains
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_Init( mesh )
    implicit none
    class(MeshBase3D), intent(in) :: mesh

    type(LineElement) :: elem1D

    integer :: ph, pz, p0

    real(RP), allocatable :: invV_VPOrdM1(:,:)
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Init( mesh )

    call elem1D%Init( mesh%refElem3D%PolyOrder_h, .false. )
    allocate( D1D(elem1D%Np,elem1D%Np) )
    allocate( D1Dt(elem1D%Np,elem1D%Np) )
    D1D(:,:) = elem1D%Dx1(:,:)
    D1Dt(:,:) = transpose( elem1D%Dx1(:,:) )

    allocate( IntrpMat_VPOrdM1_tr(elem1D%Np,elem1D%Np) )
    allocate( invV_VPOrdM1(elem1D%Np,elem1D%Np) )
    invV_VPOrdM1(:,:) = elem1D%invV(:,:)
    invV_VPOrdM1(elem1D%Np,:) = 0.0_RP
    IntrpMat_VPOrdM1_tr(:,:) = transpose( matmul(elem1D%V,invV_VPOrdM1) )

    call elem1D%Final()

    ! p0 = 1
    ! do pz=1, elem1D%Np
    !   LOG_INFO("CHECK p0=1:",'(i3,a,64E12.3)') pz, ":", IntrpMat_VPOrdM1(p0,mesh%refElem3D%Hslice(:,pz))
    ! end do
    ! p0 = 3
    ! do pz=1, elem1D%Np
    !   LOG_INFO("CHECK p0=3:",'(i3,a,64E12.3)') pz, ":", IntrpMat_VPOrdM1(p0,mesh%refElem3D%Hslice(:,pz))
    ! end do
    ! p0 = 12
    ! do pz=1, elem1D%Np
    !   LOG_INFO("CHECK p0=10:",'(i3,a,64E12.3)') pz, ":", IntrpMat_VPOrdM1(p0,mesh%refElem3D%Hslice(:,pz))
    ! end do

    ! p0 = 1
    ! do pz=1, elem1D%Np
    !   LOG_INFO("CHECK00 p0=10:",'(i3,a,64E12.3)') pz, ":", IntrpMat_VPOrdM1_tr(:,pz)
    ! end do

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
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_new( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                   & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, PRES_hyd_ref, & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot,                                                  & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )                     ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_new

    use scale_matrix_vec_kernel
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
    real(RP) :: Fx0(elem%Np), Fy0(elem%Np), Fz0(elem%Np)
    real(RP) :: DPRES_hyd, GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP) :: del_flux_hyd(elem%NfpTot,2,lmesh%Ne)
    real(RP) :: GsqrtRHOT_
    real(RP) :: rdens_, u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%Np)
    real(RP) :: Cori(elem%Np)
    real(RP) :: drho(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np)

    integer :: ke, ke2d

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR       

    integer :: p, p1, p2
    real(RP) :: tmp5(elem%Np,5)
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
    call PROF_rapstart('cal_dyn_tend_interior_test', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel private( p, p1, p2, ke, ke2d, Cori,  &
    !$omp GsqrtRHOT_, rdens_, u_, v_, w_, wt_,          &
    !$omp drho, DPRES_hyd, GradPhyd_x, GradPhyd_y,      &
    !$omp GsqrtV, RGsqrtV,                              &
    !$omp Fx, Fy, Fz, LiftDelFlx, Fx0, Fy0, Fz0, tmp5)
    !$omp do
    do ke=lmesh%NeS, lmesh%NeE
!!OCL PREFETCH_WRITE(DENS_dt(:,ke+1),level=2)
!!OCL PREFETCH_WRITE(MOMX_dt(:,ke+1),level=2)
!!OCL PREFETCH_WRITE(MOMY_dt(:,ke+1),level=2)
!!OCL PREFETCH_WRITE(MOMZ_dt(:,ke+1),level=2)
!!OCL PREFETCH_WRITE(RHOT_dt(:,ke+1),level=2)
!!OCL PREFETCH_READ(PRES_hyd(:,ke+1),level=2)
!!OCL PREFETCH_READ(DENS_hyd(:,ke+1),level=2)
!!OCL PREFETCH_READ(DDENS_(:,ke+1),level=2)
!!OCL PREFETCH_READ(MOMX_(:,ke+1),level=2)
!!OCL PREFETCH_READ(MOMY_(:,ke+1),level=2)
!!OCL PREFETCH_READ(MOMZ_(:,ke+1),level=2)
!!OCL PREFETCH_READ(DRHOT_(:,ke+1),level=2)

      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      Cori(:) = CORIOLIS(elem%IndexH2Dto3D(:),ke2d)
      GsqrtV(:)  = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D(:),ke2d)

      do p=1, elem%Np
        rdens_ = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))
        u_ (p) = MOMX_(p,ke) * rdens_
        v_ (p) = MOMY_(p,ke) * rdens_
        w_ (p) = MOMZ_(p,ke) * rdens_
        RGsqrtV(p) = 1.0_RP / GsqrtV(p)
        wt_(p) = w_(p) * RGsqrtV(p) + lmesh%GI3(p,ke,1) * u_(p) + lmesh%GI3(p,ke,2) * v_(p) 
      end do

      !-- Gradient hydrostatic pressure
      do p=1, elem%Np
        DPRES_hyd = PRES_hyd(p,ke) - PRES_hyd_ref(p,ke)
        Fx0(p) = GsqrtV(p) * DPRES_hyd
        Fz0(p) = GsqrtV(p) * lmesh%GI3(p,ke,1) * DPRES_hyd
      end do
      ! call matrix_vec_kernel_dx_8(D1D , Fx0, Fx(:))
      ! call matrix_vec_kernel_dz_8(D1Dt, Fz0, Fz(:))
      ! call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux_hyd(:,1,ke), LiftDelFlx(:))

      Fx(:) = 1.0_RP* dble(ke); Fy(:) = 4.0_RP * dble(ke); Fz(:) = 2.0_RP* dble(ke)
      GradPhyd_x(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      do p=1, elem%Np
        DPRES_hyd = PRES_hyd(p,ke) - PRES_hyd_ref(p,ke)
        Fy0(p) = GsqrtV(p) * DPRES_hyd
        Fz0(p) = GsqrtV(p) * lmesh%GI3(p,ke,2) * DPRES_hyd
      end do
      ! call matrix_vec_kernel_dy_8(D1Dt, Fy0, Fy(:))
      ! call matrix_vec_kernel_dz_8(D1Dt, Fz0, Fz(:))
      ! call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux_hyd(:,2,ke), LiftDelFlx(:))

      Fx(:) = 1.0_RP* dble(ke); Fy(:) = 4.0_RP * dble(ke); Fz(:) = 2.0_RP* dble(ke)
      GradPhyd_y(:) = lmesh%Escale(:,ke,2,2) * Fy(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      !-- DENS
      do p=1, elem%Np
        Fx0(p) = lmesh%Gsqrt(p,ke) * MOMX_(p,ke)
        Fy0(p) = lmesh%Gsqrt(p,ke) * MOMY_(p,ke)
        Fz0(p) = lmesh%Gsqrt(p,ke) * ( DDENS_(p,ke) + DENS_hyd(p,ke) ) * wt_(p)
      end do
      ! call matrix_vec_kernel_DxDyDz_8( D1D, D1Dt, Fx0, Fy0, Fz0, &
      !   Fx, Fy, Fz )
      ! call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux(:,DENS_VID,ke), LiftDelFlx(:))

      Fx(:) = 1.0_RP* dble(ke); Fy(:) = 4.0_RP * dble(ke); Fz(:) = 2.0_RP* dble(ke)
      DENS_dt(:,ke) = - ( &
!      tmp5(:,1) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)    &
          + lmesh%Escale(:,ke,2,2) * Fy(:)    &
          + lmesh%Escale(:,ke,3,3) * Fz(:)    &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 
      
      !-- MOMX
      do p=1, elem%Np
        Fx0(p) = lmesh%Gsqrt(p,ke) * (  u_(p) * MOMX_(p,ke) + DPRES_(p,ke) )
        Fy0(p) = lmesh%Gsqrt(p,ke) *    v_(p) * MOMX_(p,ke)
        Fz0(p) = lmesh%Gsqrt(p,ke) * ( wt_(p) * MOMX_(p,ke)                &
                                     + lmesh%GI3(p,ke,1) * DPRES_(p,ke)    )
      end do
      ! call matrix_vec_kernel_DxDyDz_8( D1D, D1Dt, Fx0, Fy0, Fz0, &
      !   Fx, Fy, Fz )
      ! call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux(:,MOMX_VID,ke), LiftDelFlx(:))

      Fx(:) = 1.0_RP* dble(ke); Fy(:) = 4.0_RP * dble(ke); Fz(:) = 2.0_RP* dble(ke)
      MOMX_dt(:,ke) = &
!      tmp5(:,2) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) &
          - GradPhyd_x(:) * RGsqrtV(:)            &
          + Cori(:) * MOMY_(:,ke)

      !-- MOMY
      do p=1, elem%Np
        Fx0(p) = lmesh%Gsqrt(p,ke) *    u_(p) * MOMY_(p,ke)
        Fy0(p) = lmesh%Gsqrt(p,ke) * (  v_(p) * MOMY_(p,ke) + DPRES_(p,ke) )
        Fz0(p) = lmesh%Gsqrt(p,ke) * ( wt_(p) * MOMY_(p,ke)                &
                                     + lmesh%GI3(p,ke,2) * DPRES_(p,ke)    )
      end do
    !   call matrix_vec_kernel_DxDyDz_8( D1D, D1Dt, Fx0, Fy0, Fz0, &
    !     Fx, Fy, Fz )
    !  call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux(:,MOMY_VID,ke), LiftDelFlx(:))

      Fx(:) = 1.0_RP* dble(ke); Fy(:) = 4.0_RP * dble(ke); Fz(:) = 2.0_RP* dble(ke)
      MOMY_dt(:,ke) = &
!      tmp5(:,3) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) &
          - GradPhyd_y(:) * RGsqrtV(:)            &
          - Cori(:) * MOMX_(:,ke)

      !-- MOMZ
      do p=1, elem%Np
        Fx0(p) = lmesh%Gsqrt(p,ke) *    u_(p) * MOMZ_(p,ke)
        Fy0(p) = lmesh%Gsqrt(p,ke) *    v_(p) * MOMZ_(p,ke)
        Fz0(p) = lmesh%Gsqrt(p,ke) * ( wt_(p) * MOMZ_(p,ke)      &
                                    + RGsqrtV(p) * DPRES_(p,ke) )
      end do
      ! call matrix_vec_kernel_DxDyDz_8( D1D, D1Dt, Fx0, Fy0, Fz0, &
      !   Fx, Fy, Fz )
      ! call matrix_vec_kernel_Lift_8( Lift%val, lmesh%Fscale(:,ke) * del_flux(:,MOMZ_VID,ke), LiftDelFlx(:))
      ! call matrix_vec_kernel_VPM1_8( IntrpMat_VPOrdM1_tr, DDENS_(:,ke), drho )

      Fx(:) = 1.0_RP* dble(ke); Fy(:) = 4.0_RP * dble(ke); Fz(:) = 2.0_RP* dble(ke)      
      MOMZ_dt(:,ke) = &
!      tmp5(:,4) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) &
          - Grav * drho(:)

      !-- RHOT          
      do p=1, elem%Np
        GsqrtRHOT_ = lmesh%Gsqrt(p,ke) * ( P0ovR * (PRES_hyd(p,ke) * rP0)**rgamm + DRHOT_(p,ke) )
        Fx0(p) = GsqrtRHOT_ *  u_(p)
        Fy0(p) = GsqrtRHOT_ *  v_(p)
        Fz0(p) = GsqrtRHOT_ * wt_(p)
      end do
      ! call matrix_vec_kernel_DxDyDz_8( D1D, D1Dt, Fx0, Fy0, Fz0, &
      !   Fx, Fy, Fz )
      ! call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux(:,RHOT_VID,ke), LiftDelFlx(:))

      Fx(:) = 1.0_RP* dble(ke); Fy(:) = 4.0_RP * dble(ke); Fz(:) = 2.0_RP* dble(ke)
      RHOT_dt(:,ke) = &
!      tmp5(:,5) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 

      !--
    end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_interior_test', 3)

    call PROF_rapstart('cal_dyn_tend_interior', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel private( p, p1, p2, ke, ke2d, Cori,  &
    !$omp GsqrtRHOT_, rdens_, u_, v_, w_, wt_,          &
    !$omp drho, DPRES_hyd, GradPhyd_x, GradPhyd_y,      &
    !$omp GsqrtV, RGsqrtV,                              &
    !$omp Fx, Fy, Fz, LiftDelFlx, Fx0, Fy0, Fz0)
    !$omp do
    do ke=lmesh%NeS, lmesh%NeE
!!OCL PREFETCH_WRITE(DENS_dt(:,ke+1),level=2)
!!OCL PREFETCH_WRITE(MOMX_dt(:,ke+1),level=2)
!!OCL PREFETCH_WRITE(MOMY_dt(:,ke+1),level=2)
!!OCL PREFETCH_WRITE(MOMZ_dt(:,ke+1),level=2)
!!OCL PREFETCH_WRITE(RHOT_dt(:,ke+1),level=2)
!!OCL PREFETCH_READ(PRES_hyd(:,ke+1),level=2)
!!OCL PREFETCH_READ(DENS_hyd(:,ke+1),level=2)
!!OCL PREFETCH_READ(DDENS_(:,ke+1),level=2)
!!OCL PREFETCH_READ(MOMX_(:,ke+1),level=2)
!!OCL PREFETCH_READ(MOMY_(:,ke+1),level=2)
!!OCL PREFETCH_READ(MOMZ_(:,ke+1),level=2)
!!OCL PREFETCH_READ(DRHOT_(:,ke+1),level=2)

      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      Cori(:) = CORIOLIS(elem%IndexH2Dto3D(:),ke2d)
      GsqrtV(:)  = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D(:),ke2d)

      do p=1, elem%Np
        rdens_ = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))
        u_ (p) = MOMX_(p,ke) * rdens_
        v_ (p) = MOMY_(p,ke) * rdens_
        w_ (p) = MOMZ_(p,ke) * rdens_
        RGsqrtV(p) = 1.0_RP / GsqrtV(p)
        wt_(p) = w_(p) * RGsqrtV(p) + lmesh%GI3(p,ke,1) * u_(p) + lmesh%GI3(p,ke,2) * v_(p) 
      end do
      
      !-- Gradient hydrostatic pressure
      do p=1, elem%Np
        DPRES_hyd = PRES_hyd(p,ke) - PRES_hyd_ref(p,ke)
        Fx0(p) = GsqrtV(p) * DPRES_hyd
        Fz0(p) = GsqrtV(p) * lmesh%GI3(p,ke,1) * DPRES_hyd
      end do
      call matrix_vec_kernel_dx_8(D1D , Fx0, Fx(:))
      call matrix_vec_kernel_dz_8(D1Dt, Fz0, Fz(:))
      call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux_hyd(:,1,ke), LiftDelFlx(:))

      GradPhyd_x(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      do p=1, elem%Np
        DPRES_hyd = PRES_hyd(p,ke) - PRES_hyd_ref(p,ke)
        Fy0(p) = GsqrtV(p) * DPRES_hyd
        Fz0(p) = GsqrtV(p) * lmesh%GI3(p,ke,2) * DPRES_hyd
      end do
      call matrix_vec_kernel_dy_8(D1Dt, Fy0, Fy(:))
      call matrix_vec_kernel_dz_8(D1Dt, Fz0, Fz(:))
      call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux_hyd(:,2,ke), LiftDelFlx(:))

      GradPhyd_y(:) = lmesh%Escale(:,ke,2,2) * Fy(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      !-- DENS
      do p=1, elem%Np
        Fx0(p) = lmesh%Gsqrt(p,ke) * MOMX_(p,ke)
        Fy0(p) = lmesh%Gsqrt(p,ke) * MOMY_(p,ke)
        Fz0(p) = lmesh%Gsqrt(p,ke) * ( DDENS_(p,ke) + DENS_hyd(p,ke) ) * wt_(p)
      end do
      call matrix_vec_kernel_DxDyDz_8( D1D, D1Dt, Fx0, Fy0, Fz0, &
        Fx, Fy, Fz )
      call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux(:,DENS_VID,ke), LiftDelFlx(:))

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)    &
          + lmesh%Escale(:,ke,2,2) * Fy(:)    &
          + lmesh%Escale(:,ke,3,3) * Fz(:)    &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 
      
      !-- MOMX
      do p=1, elem%Np
        Fx0(p) = lmesh%Gsqrt(p,ke) * (  u_(p) * MOMX_(p,ke) + DPRES_(p,ke) )
        Fy0(p) = lmesh%Gsqrt(p,ke) *    v_(p) * MOMX_(p,ke)
        Fz0(p) = lmesh%Gsqrt(p,ke) * ( wt_(p) * MOMX_(p,ke)                &
                                     + lmesh%GI3(p,ke,1) * DPRES_(p,ke)    )
      end do
      call matrix_vec_kernel_DxDyDz_8( D1D, D1Dt, Fx0, Fy0, Fz0, &
        Fx, Fy, Fz )
      call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux(:,MOMX_VID,ke), LiftDelFlx(:))

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) &
          - GradPhyd_x(:) * RGsqrtV(:)            &
          + Cori(:) * MOMY_(:,ke)

      !-- MOMY
      do p=1, elem%Np
        Fx0(p) = lmesh%Gsqrt(p,ke) *    u_(p) * MOMY_(p,ke)
        Fy0(p) = lmesh%Gsqrt(p,ke) * (  v_(p) * MOMY_(p,ke) + DPRES_(p,ke) )
        Fz0(p) = lmesh%Gsqrt(p,ke) * ( wt_(p) * MOMY_(p,ke)                &
                                     + lmesh%GI3(p,ke,2) * DPRES_(p,ke)    )
      end do
      call matrix_vec_kernel_DxDyDz_8( D1D, D1Dt, Fx0, Fy0, Fz0, &
        Fx, Fy, Fz )
     call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux(:,MOMY_VID,ke), LiftDelFlx(:))

      MOMY_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) &
          - GradPhyd_y(:) * RGsqrtV(:)            &
          - Cori(:) * MOMX_(:,ke)

      !-- MOMZ
      do p=1, elem%Np
        Fx0(p) = lmesh%Gsqrt(p,ke) *    u_(p) * MOMZ_(p,ke)
        Fy0(p) = lmesh%Gsqrt(p,ke) *    v_(p) * MOMZ_(p,ke)
        Fz0(p) = lmesh%Gsqrt(p,ke) * ( wt_(p) * MOMZ_(p,ke)      &
                                    + RGsqrtV(p) * DPRES_(p,ke) )
      end do
      call matrix_vec_kernel_DxDyDz_8( D1D, D1Dt, Fx0, Fy0, Fz0, &
        Fx, Fy, Fz )
      call matrix_vec_kernel_Lift_8( Lift%val, lmesh%Fscale(:,ke) * del_flux(:,MOMZ_VID,ke), LiftDelFlx(:))
      call matrix_vec_kernel_VPM1_8( IntrpMat_VPOrdM1_tr, DDENS_(:,ke), drho )
      
      MOMZ_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) &
          - Grav * drho(:)

      !-- RHOT
      do p=1, elem%Np
        GsqrtRHOT_ = lmesh%Gsqrt(p,ke) * ( P0ovR * (PRES_hyd(p,ke) * rP0)**rgamm + DRHOT_(p,ke) )
        Fx0(p) = GsqrtRHOT_ *  u_(p)
        Fy0(p) = GsqrtRHOT_ *  v_(p)
        Fz0(p) = GsqrtRHOT_ * wt_(p)
      end do
      call matrix_vec_kernel_DxDyDz_8( D1D, D1Dt, Fx0, Fy0, Fz0, &
        Fx, Fy, Fz )
      call matrix_vec_kernel_Lift_8(Lift%val, lmesh%Fscale(:,ke) * del_flux(:,RHOT_VID,ke), LiftDelFlx(:))

      RHOT_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 

      !--
    end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_new

end module scale_atm_dyn_dgm_nonhydro3d_rhot_heve
