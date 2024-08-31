!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Tracer advection
!!
!! @par Description
!!      HEVE DGM scheme for tracer advection. 
!!
!!      To preserve nonnegativity, a limiter proposed by Light and Durran (2016, MWR) is used: 
!!      we apply FCT for the lowest mode and the truncation and mass aware rescaling (TMAR).
!!
!! @par Reference
!!  - Light and Durran 2016: 
!!    Preserving Nonnegativity in Discontinuous Galerkin Approximations to Scalar Transport via Truncation and Mass Aware Rescaling (TMAR).
!!    Monthly Weather Review, 144(12), 4771–4786.
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_trcadvect3d_heve
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_sparsemat, only: &
    SparseMat, sparsemat_matmul
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
  public :: atm_dyn_dgm_trcadvect3d_heve_Init
  public :: atm_dyn_dgm_trcadvect3d_heve_Final
  public :: atm_dyn_dgm_trcadvect3d_heve_calc_fct_coef  
  public :: atm_dyn_dgm_trcadvect3d_heve_cal_tend
  public :: atm_dyn_dgm_trcadvect3d_TMAR
  public :: atm_dyn_dgm_trcadvect3d_save_massflux
  public :: atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_advtest

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  private :: atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_dyn
  private :: atm_dyn_dgm_trcadvect3d_heve_get_netOutwardFlux_generalhvc
  private :: atm_dyn_dgm_trcadvect3d_heve_get_delflux_generalhvc

contains

!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_heve_Init( mesh, FaceIntMat )
    use scale_polynominal, only: Polynominal_GenGaussLobattoPtIntWeight
    implicit none
    class(MeshBase3D), intent(in), target :: mesh
    type(SparseMat), intent(inout) :: FaceIntMat

    class(LocalMesh3D), pointer :: lcmesh
    class(ElementBase3D), pointer :: elem
    real(RP), allocatable :: intWeight_lgl1DPts_h(:)
    real(RP), allocatable :: intWeight_lgl1DPts_v(:)   
    real(RP), allocatable :: intWeight_h(:) 
    real(RP), allocatable :: intWeight_v(:)  
    
    integer :: f
    integer :: i, j, k, l
    integer :: is, ie

    real(RP), allocatable :: IntWeight(:,:)
    !-------------------------------------------------------------

    lcmesh => mesh%lcmesh_list(1)
    elem => lcmesh%refElem3D
    allocate( IntWeight(elem%Nfaces,elem%NfpTot) )
    IntWeight(:,:) = 0.0_RP

    allocate( intWeight_lgl1DPts_h(elem%Nnode_h1D) )
    allocate( intWeight_lgl1DPts_v(elem%Nnode_v) )
    allocate( intWeight_h(elem%Nnode_h1D*elem%Nnode_v) )
    allocate( intWeight_v(elem%Nnode_h1D**2) )

    intWeight_lgl1DPts_h(:) = Polynominal_GenGaussLobattoPtIntWeight(elem%PolyOrder_h)
    intWeight_lgl1DPts_v(:) = Polynominal_GenGaussLobattoPtIntWeight(elem%PolyOrder_v)

    do f=1, elem%Nfaces_h
      do k=1, elem%Nnode_v
      do i=1, elem%Nnode_h1D
        l = i + (k-1)*elem%Nnode_h1D
        intWeight_h(l) = intWeight_lgl1DPts_h(i) * intWeight_lgl1DPts_v(k)
      end do
      end do

      is = (f-1)*elem%Nfp_h + 1
      ie = is + elem%Nfp_h - 1
      IntWeight(f,is:ie) = intWeight_h(:)
    end do

    do f=1, elem%Nfaces_v
      do j=1, elem%Nnode_h1D
      do i=1, elem%Nnode_h1D
        l = i + (j-1)*elem%Nnode_h1D
        intWeight_v(l) = intWeight_lgl1DPts_h(i) * intWeight_lgl1DPts_h(j)
      end do
      end do

      is = elem%Nfaces_h*elem%Nfp_h + (f-1)*elem%Nfp_v + 1
      ie = is + elem%Nfp_v - 1
      IntWeight(elem%Nfaces_h+f,is:ie) = intWeight_v(:)
    end do

    call FaceIntMat%Init( IntWeight )

    return
  end subroutine atm_dyn_dgm_trcadvect3d_heve_Init
  
!OCL SERIAL  
  subroutine atm_dyn_dgm_trcadvect3d_heve_Final()
    implicit none
    !--------------------------------------------

    return    
  end subroutine atm_dyn_dgm_trcadvect3d_heve_Final


!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_heve_cal_tend( &
    QTRC_dt,                                        & ! (out)
    QTRC_, MOMX_, MOMY_, MOMZ_,                     & ! (in)
    alphDENS_M, alphDENS_P, fct_coef,               & ! (in)
    RHOQ_tp,                                        & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, FaceIntMat,       & ! (in) 
    lmesh, elem, lmesh2D, elem2D                    ) ! (in)

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift, FaceIntMat
    
    real(RP), intent(out) :: QTRC_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: QTRC_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lmesh%NeA)    
    real(RP), intent(in) :: MOMZ_(elem%Np,lmesh%NeA)    
    real(RP), intent(in) :: alphDENS_M(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: alphDENS_P(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: fct_coef(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: RHOQ_tp(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)    
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
    real(RP) :: momwt_(elem%Np)    
    
    integer :: ke, ke2d 

    real(RP) :: Q0, Q1, vol
    integer :: p
    !---------------------------------------------------------------------------
    
    call PROF_rapstart('cal_trcadv_tend_bndflux', 3)
    call atm_dyn_dgm_trcadvect3d_heve_get_delflux_generalhvc( &
      del_flux,                                                                & ! (out)
      QTRC_, MOMX_, MOMY_, MOMZ_, alphDens_M, alphDens_P, fct_coef,            & ! (in)
      lmesh%Gsqrt, lmesh%GsqrtH, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),           & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%J(:,:), lmesh%Fscale(:,:),                                         & ! (in) 
      lmesh%vmapM, lmesh%vmapP, elem%IndexH2Dto3D_bnd, FaceIntMat,             & ! (in)
      lmesh, elem, lmesh2D, elem2D                                             ) ! (in)
    call PROF_rapend('cal_trcadv_tend_bndflux', 3)

    call PROF_rapstart('cal_trcadv_tend_interior', 3)
    !$omp parallel do private( &
    !$omp momwt_, Fx, Fy, Fz, LiftDelFlx,  &    
    !$omp ke, ke2d                         )  
    do ke=lmesh%NeS, lmesh%NeE
      ke2d = lmesh%EMap3Dto2D(ke)
      
      momwt_(:) = MOMZ_(:,ke) * lmesh%GsqrtH(elem%IndexH2Dto3D,ke2d) / lmesh%Gsqrt(:,ke) &
                + lmesh%GI3(:,ke,1) * MOMX_(:,ke) + lmesh%GI3(:,ke,2) * MOMY_(:,ke)

      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * MOMX_(:,ke) * QTRC_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * MOMY_(:,ke) * QTRC_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * momwt_(:)   * QTRC_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx)

      QTRC_dt(:,ke) = ( &
          - lmesh%Escale(:,ke,1,1) * Fx(:)     &
          - lmesh%Escale(:,ke,2,2) * Fy(:)     &
          - lmesh%Escale(:,ke,3,3) * Fz(:)     &
          - LiftDelFlx(:)                   )  &
          / lmesh%Gsqrt(:,ke)                  &
          + RHOQ_tp(:,ke)
    end do
    call PROF_rapend('cal_trcadv_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_trcadvect3d_heve_cal_tend

!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_heve_calc_fct_coef( &
    fct_coef,                                                     & ! (out)
    QTRC_, MOMX_, MOMY_, MOMZ_, RHOQ_tp_, AlphDens_M, AlphDens_P, & ! (in)
    DENS_hyd, DDENS_, DDENS0_, rk_c_ssm1, dt,                     & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, FaceIntMat,                     & ! (in)
    lmesh, elem, lmesh2D, elem2D,                                 & ! (in)
    disable_limiter                                               ) ! (in)

    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D    
    real(RP), intent(out) :: fct_coef(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: QTRC_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lmesh%NeA)    
    real(RP), intent(in) :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: RHOQ_tp_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: alphDENS_M(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: alphDENS_P(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lmesh%NeA)     
    real(RP), intent(in) :: DDENS_ (elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DDENS0_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: rk_c_ssm1
    real(RP), intent(in) :: dt
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift, FaceIntMat
    logical, intent(in), optional :: disable_limiter

    real(RP) :: netOutwardFlux(lmesh%Ne)
    real(RP) :: momwt_(elem%Np)    
    
    integer :: ke
    real(RP) :: Q
    real(RP) :: dens_ssm1(elem%Np)
    !---------------------------------------------------------------------------

    if ( present(disable_limiter) ) then
      if ( disable_limiter ) then
        !$omp parallel do
        do ke=lmesh%NeS, lmesh%NeE
          fct_coef(:,ke) = 1.0_RP
        end do
        return
      end if 
    end if

    call PROF_rapstart('cal_trcadv_fct_coef_bndflux', 3)
    call atm_dyn_dgm_trcadvect3d_heve_get_netOutwardFlux_generalhvc( &
      netOutwardFlux,                                                                   & ! (out)
      QTRC_, MOMX_, MOMY_, MOMZ_, AlphDens_M, AlphDens_P,                               & ! (in)
      lmesh%Gsqrt, lmesh%GsqrtH, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),                    & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),           & ! (in)
      lmesh%J(:,:), lmesh%Fscale(:,:), lmesh%vmapM, lmesh%vmapP, elem%IndexH2Dto3D_bnd, & ! (in)
      FaceIntMat,                                                                       & ! (in)
      lmesh, elem, lmesh2D, elem2D                                                      ) ! (in)
    call PROF_rapend('cal_trcadv_fct_coef_bndflux', 3)

    call PROF_rapstart('cal_trcadv_fct_coef', 3)

    !$omp parallel do private( ke, Q, dens_ssm1 )  
    do ke=lmesh%NeS, lmesh%NeE

      dens_ssm1(:) = DENS_hyd(:,ke) &
                    + ( 1.0_RP - rk_c_ssm1 ) * DDENS0_(:,ke) + rk_c_ssm1 * DDENS_(:,ke)
      Q = sum( lmesh%Gsqrt(:,ke) * lmesh%J(:,ke) * elem%IntWeight_lgl(:) * ( dens_ssm1(:) *  QTRC_(:,ke) / dt + RHOQ_tp_(:,ke) ) )      
      
      fct_coef(:,ke) = max( 0.0_RP, min( 1.0_RP, Q / ( netOutwardFlux(ke) + 1.0E-10_RP ) ) )
    end do

    call PROF_rapend('cal_trcadv_fct_coef', 3)

    return
  end subroutine atm_dyn_dgm_trcadvect3d_heve_calc_fct_coef


!OCL SERIAL
!> Second Step of limiter in which nonlinear truncation and mass aware rescaling (TMAR)
  subroutine atm_dyn_dgm_trcadvect3d_TMAR( QTRC_,  & ! (inout)
    DENS_hyd, DDENS_, lmesh, elem, lmesh2D, elem2D ) ! (in)
    
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(inout) :: QTRC_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lmesh%NeA)     
    real(RP), intent(in) :: DDENS_  (elem%Np,lmesh%NeA)

    integer :: ke
    real(RP) :: Q0, Q1
    real(RP) :: Q(elem%Np)
    real(RP) :: dens(elem%Np)
    !--------------------------------------------------------

    !$omp parallel do private( dens, Q, Q0, Q1 )
    do ke=lmesh%NeS, lmesh%NeE
      dens(:) = DENS_hyd(:,ke) + DDENS_(:,ke)
      Q(:) = max( 0.0_RP, QTRC_(:,ke) )

      Q0 = sum( lmesh%Gsqrt(:,ke) * lmesh%J(:,ke) * elem%IntWeight_lgl(:) * dens(:) * QTRC_(:,ke) )
      Q1 = sum( lmesh%Gsqrt(:,ke) * lmesh%J(:,ke) * elem%IntWeight_lgl(:) * dens(:) * Q(:)        )
      QTRC_(:,ke) = Q0 / ( Q1 + 1.0E-32_RP ) * Q(:)
    end do

    return
  end subroutine atm_dyn_dgm_trcadvect3d_TMAR

!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_save_massflux( &
    MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg, alph_dens_M, alph_dens_P, &
    DDENS, MOMX, MOMY, MOMZ, DPRES, DENS_hyd, PRES_hyd,              &
    Rtot, CVtot, CPtot,                                              &
    lmesh, elem, rkstage, tavg_weight_h, tavg_weight_v               ) 
   
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) ::  MFLX_x_tavg(elem%Np,lmesh%NeA)
    real(RP), intent(inout) ::  MFLX_y_tavg(elem%Np,lmesh%NeA)
    real(RP), intent(inout) ::  MFLX_z_tavg(elem%Np,lmesh%NeA)
    real(RP), intent(inout) ::  alph_dens_M(elem%NfpTot,lmesh%Ne)
    real(RP), intent(inout) ::  alph_dens_P(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  DDENS(elem%Np,lmesh%NeA)
    real(RP), intent(in) ::  MOMX(elem%Np,lmesh%NeA)  
    real(RP), intent(in) ::  MOMY(elem%Np,lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ(elem%Np,lmesh%NeA)  
    real(RP), intent(in) ::  DPRES(elem%Np,lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) ::  Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) ::  CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) ::  CPtot(elem%Np,lmesh%NeA)
    integer, intent(in) :: rkstage
    real(RP), intent(in) :: tavg_weight_h
    real(RP), intent(in) :: tavg_weight_v

    integer :: ke
    !--------------------------------------------------------------

    !$omp parallel do private(ke)
    do ke=lmesh%NeS, lmesh%NeE
      if (rkstage == 1) then
        MFLX_x_tavg(:,ke) = tavg_weight_h * MOMX(:,ke)
        MFLX_y_tavg(:,ke) = tavg_weight_h * MOMY(:,ke)
        MFLX_z_tavg(:,ke) = tavg_weight_v * MOMZ(:,ke)
        alph_dens_M(:,ke) = 0.0_RP
        alph_dens_P(:,ke) = 0.0_RP
      else
        MFLX_x_tavg(:,ke) = MFLX_x_tavg(:,ke) + tavg_weight_h * MOMX(:,ke)
        MFLX_y_tavg(:,ke) = MFLX_y_tavg(:,ke) + tavg_weight_h * MOMY(:,ke)
        MFLX_z_tavg(:,ke) = MFLX_z_tavg(:,ke) + tavg_weight_v * MOMZ(:,ke)
      end if  
    end do

    call atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_dyn( &
      alph_dens_M, alph_dens_P,                                               & ! (out)
      DDENS, MOMX, MOMY, MOMZ, DPRES, DENS_hyd, PRES_hyd,                     & ! (in)
      Rtot, CVtot, CPtot,                                                     & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, lmesh%refElem3D,                       & ! (in)
      tavg_weight_h, tavg_weight_v                                            ) ! (in)
    
    return
  end subroutine atm_dyn_dgm_trcadvect3d_save_massflux

!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_advtest( alph_dens_M, alph_dens_P, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DENS_hyd,                                                &
    Gsqrt, nx, ny, nz, vmapM, vmapP, lmesh, elem                                          )
 
    use scale_const, only: &
     GRAV => CONST_GRAV,  &
     Rdry => CONST_Rdry,  &
     CPdry => CONST_CPdry, &
     CVdry => CONST_CVdry, &
     PRES00 => CONST_PRE00
   
    implicit none
 
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    real(RP), intent(inout) ::  alph_dens_M(elem%NfpTot*lmesh%Ne)
    real(RP), intent(inout) ::  alph_dens_P(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha, densM, densP
    !------------------------------------------------------------------------
 
 
    !$omp parallel do private( &
    !$omp iM, iP, VelP, VelM, alpha, densM, densP )
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)
  
      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iP)
 
      VelM = ( MOMX_(iM) * nx(i) + MOMY_(iM) * ny(i) + MOMZ_(iM) * nz(i) ) / densM
      VelP = ( MOMX_(iP) * nx(i) + MOMY_(iP) * ny(i) + MOMZ_(iP) * nz(i) ) / densP
 
      alpha = max( abs(VelM), abs(VelP)  )
      alph_dens_M(i) = alpha * densM * Gsqrt(iM)
      alph_dens_P(i) = alpha * densP * Gsqrt(iP)
    end do
 
    return
  end subroutine atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_advtest

!-- private ---

!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_dyn( alph_dens_M, alph_dens_P, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DPRES_, DENS_hyd, PRES_hyd,    &
    Rtot, CVtot, CPtot,                                         &
    nx, ny, nz, vmapM, vmapP, lmesh, elem,                      &
    tavg_weight_h, tavg_weight_v )
 
    use scale_const, only: &
     GRAV => CONST_GRAV,  &
     Rdry => CONST_Rdry,  &
     CPdry => CONST_CPdry, &
     CVdry => CONST_CVdry, &
     PRES00 => CONST_PRE00
   
    implicit none
 
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    real(RP), intent(inout) ::  alph_dens_M(elem%NfpTot*lmesh%Ne)
    real(RP), intent(inout) ::  alph_dens_P(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DPRES_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Rtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CVtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CPtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: tavg_weight_h
    real(RP), intent(in) :: tavg_weight_v
    
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    real(RP) :: densM, densP
    real(RP) :: gamm, rgamm
    real(RP) :: tavg_weight 
    !------------------------------------------------------------------------
 
    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry
 
    !$omp parallel do private( &
    !$omp iM, iP, VelP, VelM, alpha, densM, densP, tavg_weight )
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iP)
 
      VelM = ( MOMX_(iM)*nx(i) + MOMY_(iM)*ny(i) + MOMZ_(iM)*nz(i) ) / densM
      VelP = ( MOMX_(iP)*nx(i) + MOMY_(iP)*ny(i) + MOMZ_(iP)*nz(i) ) / densP
 
      alpha = max( sqrt( gamm * ( PRES_hyd(iM) + DPRES_(iM) ) / densM ) + abs(VelM), &
                   sqrt( gamm * ( PRES_hyd(iP) + DPRES_(iP) ) / densP ) + abs(VelP)  )
      tavg_weight = tavg_weight_h * ( abs(nx(i)) + abs(ny(i)) ) + tavg_weight_v * abs(nz(i))

      alph_dens_M(i) = alph_dens_M(i) + tavg_weight * alpha * densM
      alph_dens_P(i) = alph_dens_P(i) + tavg_weight * alpha * densP
    end do
 
    return
  end subroutine atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_dyn

!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_heve_get_delflux_generalhvc( &
    del_flux,                                                      & ! (out)
    QTRC_, MOMX_, MOMY_, MOMZ_, AlphDens_M, AlphDens_P, fct_coef,  & ! (in)
    Gsqrt, GsqrtH, G13, G23, nx, ny, nz, J, Fscale,                & ! (in)
    vmapM, vmapP, iM2Dto3D, FaceIntMat,                            & ! (in)
    lmesh, elem, lmesh2D, elem2D                                   ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  QTRC_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: AlphDens_M(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: AlphDens_P(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: fct_coef(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GsqrtH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)     
    real(RP), intent(in) :: J(elem%Np*lmesh%Ne)    
    real(RP), intent(in) :: Fscale(elem%NfpTot,lmesh%Ne)      
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: iM2Dto3D(elem%NfpTot)
    type(SparseMat), intent(in) :: FaceIntMat
    
    integer :: ke, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: MomFlxP(elem%NfpTot), MomFlxM(elem%NfpTot), alpha(elem%NfpTot)
    real(RP) :: QTRC_P(elem%NfpTot), QTRC_M(elem%NfpTot)
    real(RP) :: GsqrtMOMX_P(elem%NfpTot), GsqrtMOMX_M(elem%NfpTot)
    real(RP) :: GsqrtMOMY_P(elem%NfpTot), GsqrtMOMY_M(elem%NfpTot)
    real(RP) :: GsqrtMOMZ_P(elem%NfpTot), GsqrtMOMZ_M(elem%NfpTot)
    real(RP) :: Gsqrt_P(elem%NfpTot), Gsqrt_M(elem%NfpTot)
    real(RP) :: GsqrtV_P(elem%NfpTot), GsqrtV_M(elem%NfpTot)
    real(RP) :: G13_M(elem%NfpTot), G13_P(elem%NfpTot)
    real(RP) :: G23_M(elem%NfpTot), G23_P(elem%NfpTot)

    real(RP) :: R_M(elem%NfpTot), R_P(elem%NfpTot)
    real(RP) :: numflux(elem%NfpTot)
    real(RP) :: outward_flux_tmp(elem%Nfaces)   
    integer :: f, p, fp, is, ie
    !------------------------------------------------------------------------

    !$omp parallel do private( &
    !$omp ke, iM, iP, ke2D,                                                                     &
    !$omp alpha, MomFlxM, MomFlxP, R_M, R_P, numflux, outward_flux_tmp,                         &
    !$omp f, p, fp,                                                                             &
    !$omp GsqrtMOMX_M, GsqrtMOMX_P, GsqrtMOMY_M, GsqrtMOMY_P, GsqrtMOMZ_M, GsqrtMOMZ_P,         &
    !$omp QTRC_M, QTRC_P, Gsqrt_P, Gsqrt_M, GsqrtV_P, GsqrtV_M, G13_P, G13_M, G23_P, G23_M      )
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

      R_M(:) = fct_coef(iM)
      R_P(:) = fct_coef(iP)  

      QTRC_M(:) = QTRC_(iM)
      QTRC_P(:) = QTRC_(iP)
      GsqrtMOMX_M (:) = Gsqrt_M(:) * MOMX_ (iM)
      GsqrtMOMX_P (:) = Gsqrt_P(:) * MOMX_ (iP)
      GsqrtMOMY_M (:) = Gsqrt_M(:) * MOMY_ (iM)
      GsqrtMOMY_P (:) = Gsqrt_P(:) * MOMY_ (iP)
      GsqrtMOMZ_M (:) = Gsqrt_M(:) * MOMZ_ (iM)
      GsqrtMOMZ_P (:) = Gsqrt_P(:) * MOMZ_ (iP)

      MomFlxM(:) = ( GsqrtMOMX_M(:) * nx(:,ke) + GsqrtMOMY_M(:) * ny(:,ke)                  &
                 + ( ( GsqrtMOMZ_M(:) / GsqrtV_M(:)                                         &
                     + G13_M(:) * GsqrtMOMX_M(:) + G23_M(:) * GsqrtMOMY_M(:) ) * nz(:,ke) ) &
                 ) 
      MomFlxP(:) = ( GsqrtMOMX_P(:) * nx(:,ke) + GsqrtMOMY_P(:) * ny(:,ke)                  &
                 + ( ( GsqrtMOMZ_P(:) / GsqrtV_P(:)                                         &
                     + G13_P(:) * GsqrtMOMX_P(:) + G23_P(:) * GsqrtMOMY_P(:) ) * nz(:,ke) ) &
                 ) 

      numflux(:) =  0.5_RP * ( ( QTRC_P(:) * MomFlxP(:) + QTRC_M(:) * MomFlxM(:) )           &
                              - AlphDENS_P(:,ke) * QTRC_P(:) + AlphDENS_M(:,ke) * QTRC_M(:)  ) 
      ! alpha(:) = max( abs(MomFlxM(:)), abs(MomFlxP(:)) )       
      ! numflux(:) =  0.5_RP * ( ( QTRC_P(:) * MomFlxP(:) + QTRC_M(:) * MomFlxM(:) )    &
      !                         - alpha(:) * (QTRC_P(:) - QTRC_M(:) )                   )  

      call sparsemat_matmul( FaceIntMat, J(iM) * Fscale(:,ke) * numflux(:), outward_flux_tmp )
      do f=1, elem%Nfaces_h
        do p=1, elem%Nfp_h
          fp = p + (f-1)*elem%Nfp_h
          del_flux(fp,ke) = numflux(fp) * 0.5_RP * ( R_P(fp) + R_M(fp) - ( R_P(fp) - R_M(fp) ) * sign( 1.0_RP, outward_flux_tmp(f) ) ) &
                          - QTRC_M(fp) * MomFlxM(fp)    
        end do
      end do 
      do f=1, elem%Nfaces_v
        do p=1, elem%Nfp_v
          fp = p + (f-1)*elem%Nfp_v + elem%Nfaces_h * elem%Nfp_h
          del_flux(fp,ke) = numflux(fp) * 0.5_RP * ( R_P(fp) + R_M(fp) - ( R_P(fp) - R_M(fp) ) * sign( 1.0_RP, outward_flux_tmp(elem%Nfaces_h+f) ) ) &
                          - QTRC_M(fp) * MomFlxM(fp) 
        end do
      end do
    end do

    return
  end subroutine atm_dyn_dgm_trcadvect3d_heve_get_delflux_generalhvc


!OCL SERIAL
  subroutine atm_dyn_dgm_trcadvect3d_heve_get_netOutwardFlux_generalhvc( &
    net_outward_flux,                                                    & ! (out)
    QTRC_, MOMX_, MOMY_, MOMZ_, AlphDens_M, AlphDens_P,                  & ! (in)
    Gsqrt, GsqrtH, G13, G23, nx, ny, nz, J, Fscale,                      & ! (in)
    vmapM, vmapP, iM2Dto3D, FaceIntMat,                                  & ! (in)
    lmesh, elem, lmesh2D, elem2D                                         ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  net_outward_flux(lmesh%Ne)
    real(RP), intent(in) ::  QTRC_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  AlphDens_M(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  AlphDens_P(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GsqrtH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: J(elem%Np*lmesh%Ne)    
    real(RP), intent(in) :: Fscale(elem%NfpTot,lmesh%Ne)    
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: iM2Dto3D(elem%NfpTot)
    type(SparseMat), intent(in) :: FaceIntMat
    
    integer :: ke, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: MomFlxP(elem%NfpTot), MomFlxM(elem%NfpTot), alpha(elem%NfpTot)
    real(RP) :: QTRC_P(elem%NfpTot), QTRC_M(elem%NfpTot)
    real(RP) :: GsqrtMOMX_P(elem%NfpTot), GsqrtMOMX_M(elem%NfpTot)
    real(RP) :: GsqrtMOMY_P(elem%NfpTot), GsqrtMOMY_M(elem%NfpTot)
    real(RP) :: GsqrtMOMZ_P(elem%NfpTot), GsqrtMOMZ_M(elem%NfpTot)
    real(RP) :: Gsqrt_P(elem%NfpTot), Gsqrt_M(elem%NfpTot)
    real(RP) :: GsqrtV_P(elem%NfpTot), GsqrtV_M(elem%NfpTot)
    real(RP) :: G13_M(elem%NfpTot), G13_P(elem%NfpTot)
    real(RP) :: G23_M(elem%NfpTot), G23_P(elem%NfpTot)

    real(RP) :: numflux(elem%NfpTot)
    real(RP) :: outward_flux_tmp(elem%Nfaces)
    !------------------------------------------------------------------------

    !$omp parallel do private( &
    !$omp ke, iM, iP, ke2D,                                                                 &
    !$omp alpha, MomFlxM, MomFlxP, numflux, outward_flux_tmp,                               &
    !$omp GsqrtMOMX_M, GsqrtMOMX_P, GsqrtMOMY_M, GsqrtMOMY_P, GsqrtMOMZ_M, GsqrtMOMZ_P,     &
    !$omp QTRC_M, QTRC_P, Gsqrt_P, Gsqrt_M, GsqrtV_P, GsqrtV_M, G13_P, G13_M, G23_P, G23_M  )
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

      QTRC_M(:) = QTRC_(iM)
      QTRC_P(:) = QTRC_(iP)
      GsqrtMOMX_M (:) = Gsqrt_M(:) * MOMX_ (iM)
      GsqrtMOMX_P (:) = Gsqrt_P(:) * MOMX_ (iP)
      GsqrtMOMY_M (:) = Gsqrt_M(:) * MOMY_ (iM)
      GsqrtMOMY_P (:) = Gsqrt_P(:) * MOMY_ (iP)
      GsqrtMOMZ_M (:) = Gsqrt_M(:) * MOMZ_ (iM)
      GsqrtMOMZ_P (:) = Gsqrt_P(:) * MOMZ_ (iP)

      MomFlxM(:) = ( GsqrtMOMX_M(:) * nx(:,ke) + GsqrtMOMY_M(:) * ny(:,ke)                  &
                 + ( ( GsqrtMOMZ_M(:) / GsqrtV_M(:)                                         &
                     + G13_M(:) * GsqrtMOMX_M(:) + G23_M(:) * GsqrtMOMY_M(:) ) * nz(:,ke) ) &
                 ) 
      MomFlxP(:) = ( GsqrtMOMX_P(:) * nx(:,ke) + GsqrtMOMY_P(:) * ny(:,ke)                  &
                 + ( ( GsqrtMOMZ_P(:) / GsqrtV_P(:)                                         &
                     + G13_P(:) * GsqrtMOMX_P(:) + G23_P(:) * GsqrtMOMY_P(:) ) * nz(:,ke) ) &
                 ) 

      alpha(:) = max( abs(MomFlxM(:)), abs(MomFlxP(:)) )

      numflux(:) = 0.5_RP * ( ( QTRC_P(:) * MomFlxP(:) + QTRC_M(:) * MomFlxM(:) )          &
                            - AlphDENS_P(:,ke) * QTRC_P(:) + AlphDENS_M(:,ke) * QTRC_M(:)  )        
      ! numflux(:) = 0.5_RP * ( ( QTRC_P(:) * MomFlxP(:) + QTRC_M(:) * MomFlxM(:) ) &
      !                       - alpha(:) * ( QTRC_P(:) - QTRC_M(:) )                )

      call sparsemat_matmul( FaceIntMat, J(iM) * Fscale(:,ke) * numflux(:), outward_flux_tmp )
      net_outward_flux(ke) = sum( max( 0.0_RP, outward_flux_tmp(:) ) )
    end do

    return
  end subroutine atm_dyn_dgm_trcadvect3d_heve_get_netOutwardFlux_generalhvc

end module scale_atm_dyn_dgm_trcadvect3d_heve