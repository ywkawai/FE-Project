!> module FElib / Atmosphere / Physics turbulence
!!
!! @par Description
!!      Sub-grid scale turbulnce process for global domain
!!      Smagorinsky-type
!!      The equations are discretized by local DGM. 
!!
!! @author Team SCALE
!!
!! @par Reference
!!  - Brown et al., 1994:
!!    Large-eddy simulaition of stable atmospheric boundary layers with a revised stochastic subgrid model.
!!    Roy. Meteor. Soc., 120, 1485-1512
!!  - Scotti et al., 1993:
!!    Generalized Smagorinsky model for anisotropic grids.
!!    Phys. Fluids A, 5, 2306-2308
!!  - Nishizawa et al., 2015:
!!    Influence of grid aspect ratio on planetary boundary layer turbulence in large-eddy simulations
!!    Geosci. Model Dev., 8, 3393–3419
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_tb_dgm_globalsmg
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    EPS  => CONST_EPS,     &
    GRAV => CONST_GRAV,    &
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  &
    CVdry => CONST_CVdry,  &
    PRES00 => CONST_PRE00, &
    KARMAN  => CONST_KARMAN, &
    RPlanet => CONST_RADIUS  

  use scale_sparsemat  
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
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
  public :: atm_phy_tb_dgm_globalsmg_Init
  public :: atm_phy_tb_dgm_globalsmg_Final
  public :: atm_phy_tb_dgm_globalsmg_cal_grad
  public :: atm_phy_tb_dgm_globalsmg_cal_grad_qtrc
  public :: atm_phy_tb_dgm_globalsmg_cal_tend
  public :: atm_phy_tb_dgm_globalsmg_cal_tend_qtrc

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: cal_del_flux_grad
  private :: cal_del_flux_grad_qtrc
  private :: cal_del_flux

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter   :: OneOverThree  = 1.0_RP / 3.0_RP
  real(RP), private, parameter   :: twoOverThree  = 2.0_RP / 3.0_RP
  real(RP), private, parameter   :: FourOverThree = 4.0_RP / 3.0_RP

  real(RP), private              :: Cs            = 0.13_RP ! Smagorinsky constant (Scotti et al. 1993)
  real(RP), private, parameter   :: PrN           = 0.7_RP  ! Prandtl number in neutral conditions
  real(RP), private, parameter   :: RiC           = 0.25_RP ! critical Richardson number
  real(RP), private, parameter   :: FmC           = 16.0_RP ! fum = sqrt(1 - c*Ri)
  real(RP), private, parameter   :: FhB           = 40.0_RP ! fuh = sqrt(1 - b*Ri)/PrN
  real(RP), private              :: RPrN                    ! 1 / PrN
  real(RP), private              :: RRiC                    ! 1 / RiC
  real(RP), private              :: OnemPrNovRiC            ! PrN / RiC
  
  ! for backscatter
  real(RP), private, parameter   :: CB   = 1.4_RP
  real(RP), private, parameter   :: CBt  = 0.45_RP
  real(RP), private, parameter   :: aN   = 0.47958315233127197_RP ! a_N = sqrt(0.23)
  real(RP), private, parameter   :: atN4 = 0.09_RP                ! a_{\theta N}^4 = 0.3**2
  real(RP), private, parameter   :: C1o  = aN**3
  real(RP), private, parameter   :: D1o  = PrN * atN4 / aN


  real(RP), private              :: filter_fac    = 2.0_RP
  real(RP), private              :: NU_MAX        = 10000.0_RP
  real(RP), private              :: tke_fac       

contains
!OCL SERIAL
  subroutine atm_phy_tb_dgm_globalsmg_Init( mesh )
    implicit none    
    class(MeshBase3D), intent(in) :: mesh

    logical  :: consistent_tke = .true.

    namelist / PARAM_ATMOS_PHY_TB_DGM_GLOBALSMG / &
      Cs,                                   &
      NU_MAX,                               &
      filter_fac,                           &
      consistent_tke
    
    integer :: ierr
    !--------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_dgm_globalsmg_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_TB_dgm_globalsmg_setup",*) 'Smagorinsky-type Eddy Viscocity Model'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_DGM_GLOBALSMG,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_TB_dgm_globalsmg_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_TB_dgm_globalsmg_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_TB_DGM_SMG. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_TB_DGM_GLOBALSMG)

    RPrN         = 1.0_RP / PrN
    RRiC         = 1.0_RP / RiC
    OnemPrNovRiC = ( 1.0_RP - PrN ) * RRiC

    if ( consistent_tke ) then
      tke_fac = 1.0_RP
    else
      tke_fac = 0.0_RP
    end if

    return
  end subroutine atm_phy_tb_dgm_globalsmg_Init

!OCL SERIAL
  subroutine atm_phy_tb_dgm_globalsmg_Final()
    implicit none
    !--------------------------------------------------------------------

    return
  end subroutine atm_phy_tb_dgm_globalsmg_Final

!OCL SERIAL  
  subroutine atm_phy_tb_dgm_globalsmg_cal_grad( &
    T11, T12, T13, T21, T22, T23, T31, T32, T33,                & ! (out)
    dPTdx, dPTdy, dPTdz,                                        & ! (out)
    TKE, Nu, Kh,                                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,    & ! (in)
    PRES, PT,                                                   & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound                                                    ) ! (in)

    use scale_atm_phy_tb_dgm_common, only: &
      atm_phy_tb_dgm_common_calc_lambda
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatVec

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) :: T11(elem%Np,lmesh%NeA), T12(elem%Np,lmesh%NeA), T13(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: T21(elem%Np,lmesh%NeA), T22(elem%Np,lmesh%NeA), T23(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: T31(elem%Np,lmesh%NeA), T32(elem%Np,lmesh%NeA), T33(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: dPTdx(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: dPTdy(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: dPTdz(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: TKE(elem%Np,lmesh%NeA)    
    real(RP), intent(out) :: Nu(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: Kh(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PT(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx, Dy, Dz
    type(SparseMat), intent(in) :: Sx, Sy, Sz
    type(SparseMat), intent(in) :: Lift
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DENS(elem%Np), RDENS(elem%Np), RHOT(elem%Np), Q(elem%Np)
    real(RP) :: DdensDxi(elem%Np,3)
    real(RP) :: DqDxi_(elem%Np,2)
    real(RP) :: DVelDxi(elem%Np,3,3)
    real(RP) :: del_flux_rho (elem%NfpTot,lmesh%Ne,3)
    real(RP) :: del_flux_mom (elem%NfpTot,lmesh%Ne,3,3)
    real(RP) :: del_flux_rhot(elem%NfpTot,lmesh%Ne,3)

    real(RP) :: Ri ! local gradient Richardson number
    real(RP) :: S2 ! (2SijSij)^1/2
    real(RP) :: fm ! factor in eddy viscosity which represents the stability dependence of 
                   ! the Brown et al (1994)'s subgrid model 
    real(RP) :: Pr ! Parandtl number (=Nu/Kh= fm/fh)

    real(RP) :: lambda  (elem%Np,lmesh%Ne) ! basic mixing length
    real(RP) :: lambda_r(elem%Np)          ! characteristic subgrid length scale 
    real(RP) :: E(elem%Np)  ! subgrid kinetic energy 
    real(RP) :: C1(elem%Np) ! factor in the relation with energy disspation rate, lambda_r, and E

    real(RP) :: G11(elem%Np), G12(elem%Np), G22(elem%Np), G33(elem%Np)
    real(RP) :: X2D(elem2D%Np,lmesh2D%Ne), Y2D(elem2D%Np,lmesh2D%Ne)
    real(RP) :: X(elem%Np), Y(elem%Np), Rdel2(elem%Np)

    real(RP) :: S11(elem%Np), S12(elem%Np), S22(elem%Np), S23(elem%Np), S31(elem%Np), S33(elem%Np)
    real(RP) :: DivOvThree(elem%Np), TKEMulTwoOvThree(elem%Np)

    real(RP) :: Sabs_tmp(elem%Np)
    real(RP) :: SIJ(elem%Np,3,3)
    real(RP) :: G_ij(elem%Np,3,3)

    integer :: ke, ke2D
    integer :: p
    integer :: i, j
    !--------------------------------------------------------------------

    call cal_del_flux_grad( del_flux_rho, del_flux_mom, del_flux_rhot,        & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, PT,            & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound )                                                   ! (in)

    call atm_phy_tb_dgm_common_calc_lambda( lambda, & ! (out)
      Cs, filter_fac, lmesh, elem, lmesh2D, elem2D  ) ! (in)
  
    !$omp parallel private( ke, ke2D, &
    !$omp Fx, Fy, Fz, LiftDelFlx,                                &
    !$omp DENS, RHOT, RDENS, Q, DdensDxi, DqDxi_, DVelDxi,            & 
    !$omp S11, S12, S22, S23, S31, S33, DivOvThree, TKEMulTwoOvThree, &
    !$omp p, Ri, S2, fm, Pr, lambda_r, E, C1,                         &
    !$omp Sabs_tmp, i, j, SIJ, G_ij, G11, G12, G22, G33,              &
    !$omp X, Y, Rdel2  )

    !$omp do
    do ke2D = lmesh2D%NeS, lmesh2D%NeE
      X2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,1))
      Y2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,2))
    end do
    !$omp end do

    !$omp do
    do ke=lmesh%NeS, lmesh%NeE
      !---
      ke2d = lmesh%EMap3Dto2D(ke)
      G11(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,1)
      G12(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,2)
      G22(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,2,2)
      G33(:) = 1.0_RP

      G_ij(:,:,:) = 0.0_RP
      G_ij(:,1,1) = lmesh%G_ij(elem%IndexH2Dto3D,ke2d,1,1)
      G_ij(:,2,1) = lmesh%G_ij(elem%IndexH2Dto3D,ke2d,2,1)
      G_ij(:,1,2) = G_ij(:,2,1)
      G_ij(:,2,2) = lmesh%G_ij(elem%IndexH2Dto3D,ke2d,2,2)
      G_ij(:,3,3) = 1.0_RP

      X(:) = X2D(elem%IndexH2Dto3D,ke2d)
      Y(:) = Y2D(elem%IndexH2Dto3D,ke2d)
      Rdel2(:) = 1.0_RP / ( 1.0_RP + X(:)**2 + Y(:)**2 )

      DENS (:) = DENS_hyd(:,ke) + DDENS_(:,ke)
      RDENS(:) = 1.0_RP / DENS(:)
      RHOT(:) = DENS(:) * PT(:,ke)

      ! gradient of density
      call sparsemat_matmul( Dx, DENS, Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rho(:,ke,1), LiftDelFlx )
      DdensDxi(:,1) = lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul( Dy, DENS, Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rho(:,ke,2), LiftDelFlx )
      DdensDxi(:,2) = lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:)

      call sparsemat_matmul( Dz, DENS, Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rho(:,ke,3), LiftDelFlx )
      DdensDxi(:,3) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)
      
      ! gradient of u
      Q(:) = MOMX_(:,ke) * RDENS(:)

      call sparsemat_matmul( Dx, MOMX_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,1,1), LiftDelFlx )
      DqDxi_(:,1) = (   lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1)       &
                         + Rdel2(:) * Y(:) *                                                           &
                           ( 2.0_RP * X(:) * Y(:) * MOMX_(:,ke) - ( 1.0_RP + Y(:)**2 ) * MOMY_(:,ke) ) & ! u^r Gam^1_1r
                        ) * RDENS(:)
      
      call sparsemat_matmul( Dy, MOMX_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2,1), LiftDelFlx )
      DqDxi_(:,2) = (   lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2)      &
                         + Rdel2(:) * Y(:) *                                                        &
                          ( - ( 1.0_RP + Y(:)**2 ) * MOMY_(:,ke) )                                  & ! u^r Gam^1_2r
                        ) * RDENS(:)

      call sparsemat_matmul( Dz, MOMX_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3,1), LiftDelFlx )
      DVelDxi(:,3,1) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)
                                                                                                 ! u^r Gam^1_3r = 0

      DivOvThree(:) = DqDxi_(:,1)
      DVelDxi(:,1,1) = G11(:) * DqDxi_(:,1) + G12(:) * DqDxi_(:,2)
      DVelDxi(:,2,1) = G12(:) * DqDxi_(:,1) + G22(:) * DqDxi_(:,2)

      ! gradient of v
      Q(:) = MOMY_(:,ke) * RDENS(:)

      call sparsemat_matmul( Dx, MOMY_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,1,2), LiftDelFlx )
      DqDxi_(:,1) = (   lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) &
                         + Rdel2(:) * X(:) *                                                   &
                           ( - ( 1.0_RP + X(:)**2 ) * MOMY_(:,ke) )                            & ! u^r Gam^2_1r
                        ) * RDENS(:)

      call sparsemat_matmul( Dy, MOMY_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2,2), LiftDelFlx )
      DqDxi_(:,2) = (   lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) &
                         + Rdel2(:) * X(:) *                                                             &
                           ( - ( 1.0_RP + X(:)**2 ) * MOMX_(:,ke) + 2.0_RP * X(:) * Y(:) * MOMY_(:,ke) ) & ! u^r Gam^2_2r
                        ) * RDENS(:)

      call sparsemat_matmul( Dz, MOMY_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3,2), LiftDelFlx )
      DVelDxi(:,3,2) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)
                                                                                                 ! u^r Gam^2_3r = 0

      DivOvThree(:) = DivOvThree(:) + DqDxi_(:,2)
      DVelDxi(:,1,2) = G11(:) * DqDxi_(:,1) + G12(:) * DqDxi_(:,2)
      DVelDxi(:,2,2) = G12(:) * DqDxi_(:,1) + G22(:) * DqDxi_(:,2)

      ! gradient of w
      Q(:) = MOMZ_(:,ke) * RDENS(:)

      call sparsemat_matmul( Dx, MOMZ_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) *del_flux_mom(:,ke,1,3), LiftDelFlx )
      DqDxi_(:,1) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) ) * RDENS(:)
                                                                                                 ! u^r Gam^3_1r = 0

      call sparsemat_matmul( Dy, MOMZ_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2,3), LiftDelFlx )
      DqDxi_(:,2) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) ) * RDENS(:)
                                                                                                 ! u^r Gam^3_2r = 0

      call sparsemat_matmul( Dz, MOMZ_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3,3), LiftDelFlx )
      DVelDxi(:,3,3) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)
                                                                                                 ! u^r Gam^3_3r = 0

      DivOvThree(:) = DivOvThree(:) + DVelDxi(:,3,3)
      DVelDxi(:,1,3) = G11(:) * DqDxi_(:,1) + G12(:) * DqDxi_(:,2)
      DVelDxi(:,2,3) = G12(:) * DqDxi_(:,1) + G22(:) * DqDxi_(:,2)

      ! gradient of pt
      Q(:) = RHOT(:) * RDENS(:)

      call sparsemat_matmul( Dx, RHOT, Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke,1), LiftDelFlx )
      DqDxi_(:,1) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) ) * RDENS(:)

      call sparsemat_matmul( Dy, RHOT, Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke,2), LiftDelFlx )
      DqDxi_(:,2) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) ) * RDENS(:)

      call sparsemat_matmul( Dz, RHOT, Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke,3), LiftDelFlx )
      dPTdz(:,ke) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)

      dPTdx(:,ke) = G11(:) * DqDxi_(:,1) + G12(:) * DqDxi_(:,2)
      dPTdy(:,ke) = G12(:) * DqDxi_(:,1) + G22(:) * DqDxi_(:,2)

      ! Calculate the component of strain velocity tensor
      S11(:) = DVelDxi(:,1,1)
      S12(:) = 0.5_RP * ( DVelDxi(:,1,2) + DVelDxi(:,2,1) )
      S22(:) = DVelDxi(:,2,2)
      S23(:) = 0.5_RP * ( DVelDxi(:,2,3) + DVelDxi(:,3,2) )
      S31(:) = 0.5_RP * ( DVelDxi(:,1,3) + DVelDxi(:,3,1) )
      S33(:) = DVelDxi(:,3,3)

      SIJ(:,1,1) = S11(:)
      SIJ(:,1,2) = S12(:)
      SIJ(:,1,3) = S31(:)
      SIJ(:,2,1) = S12(:)
      SIJ(:,2,2) = S22(:)
      SIJ(:,2,3) = S23(:)
      SIJ(:,3,1) = S31(:)
      SIJ(:,3,2) = S23(:)
      SIJ(:,3,3) = S33(:)

      Sabs_tmp(:) = 0.0_RP
      do j=1, 3
      do i=1, 3
        Sabs_tmp(:) = Sabs_tmp(:) &
        + SIJ(:,i,j) * (  &
           G_ij(:,i,1) * ( G_ij(:,j,1) * SIJ(:,1,1) + G_ij(:,j,2) * SIJ(:,1,2) + G_ij(:,j,3) * SIJ(:,1,3) ) &
         + G_ij(:,i,2) * ( G_ij(:,j,1) * SIJ(:,2,1) + G_ij(:,j,2) * SIJ(:,2,2) + G_ij(:,j,3) * SIJ(:,2,3) ) &
         + G_ij(:,i,3) * ( G_ij(:,j,1) * SIJ(:,3,1) + G_ij(:,j,2) * SIJ(:,3,2) + G_ij(:,j,3) * SIJ(:,3,3) ) )
      end do
      end do

      ! Caclulate eddy viscosity & eddy diffusivity

      do p=1, elem%Np
        S2 = 2.0_RP * Sabs_tmp(p)
        
        Ri = Grav / PT(p,ke) * dPTdz(p,ke) / max( S2, EPS )

        ! The Stability functions fm and fh are given by the appendix A of Brown et al. (1994). 
        if (Ri < 0.0_RP ) then ! unstable
          fm = sqrt( 1.0_RP - FmC * Ri )
          Nu(p,ke) = lambda(p,ke)**2 * sqrt( S2 ) * fm
          Pr = fm / sqrt( 1.0_RP - FhB * Ri ) * PrN
        else if ( Ri < RiC ) then ! stable
          fm = ( 1.0_RP - Ri * RRiC )**4
          Nu(p,ke) = lambda(p,ke)**2 * sqrt( S2 ) * fm
          Pr = PrN / ( 1.0_RP - OnemPrNovRiC * Ri )
        else ! strongly stable
          fm = 0.0_RP
          Nu(p,ke) = 0.0_RP
          Kh(p,ke) = 0.0_RP
          Pr = 1.0_RP
        end if

        if ( Ri < RiC ) then
          Kh(p,ke) = max( min( Nu(p,ke) / Pr, NU_MAX ), EPS )
          Nu(p,ke) = max( min( Nu(p,ke), NU_MAX ), EPS )
          Pr = Nu(p,ke) / Kh(p,ke)
          lambda_r(p) = lambda(p,ke) * sqrt( fm / sqrt( 1.0_RP - Ri/Pr ) )
        else
          lambda_r(p) = 0.0_RP
        end if
      end do

      ! if ( backscatter ) then
      ! else 
        E (:) = Nu(:,ke)**3 / ( lambda_r(:)**4 + EPS )
        C1(:) = C1o
      ! end if

      ! TKE
      TKE(:,ke) = ( E(:) * lambda_r(:) / C1(:) )**twoOverThree 


      ! Calculate components of parameterized stress tensor
      ! Tij = rho * { 2 Nu * [(G^im u^im + G^jm u^jm)/2 - G^ij/3 D] - G^ij 2/3 * TKE }
      !
      DivOvThree(:) = DivOvThree(:) * OneOverThree
      TKEMulTwoOvThree(:) = twoOverThree * TKE(:,ke) * tke_fac
      do p=1, elem%Np
        T11(p,ke) = DENS(p) * ( 2.0_RP * Nu(p,ke) * ( S11(p) - G11(p) * DivOvThree(p) ) - G11(p) * TKEMulTwoOvThree(p) )
        T12(p,ke) = DENS(p) * ( 2.0_RP * Nu(p,ke) * ( S12(p) - G12(p) * DivOvThree(p) ) - G12(p) * TKEMulTwoOvThree(p) )
        T13(p,ke) = DENS(p) * 2.0_RP * Nu(p,ke) * S31(p)
      end do
      do p=1, elem%Np
        T21(p,ke) = DENS(p) * ( 2.0_RP * Nu(p,ke) * ( S12(p) - G12(p) * DivOvThree(p) ) - G12(p) * TKEMulTwoOvThree(p) )
        T22(p,ke) = DENS(p) * ( 2.0_RP * Nu(p,ke) * ( S22(p) - G22(p) * DivOvThree(p) ) - G22(p) * TKEMulTwoOvThree(p) )
        T23(p,ke) = DENS(p) * 2.0_RP * Nu(p,ke) * S23(p)
      end do
      do p=1, elem%Np      
        T31(p,ke) = DENS(p) * 2.0_RP * Nu(p,ke) * S31(p)
        T32(p,ke) = DENS(p) * 2.0_RP * Nu(p,ke) * S23(p)
        T33(p,ke) = DENS(p) * ( 2.0_RP * Nu(p,ke) * ( S33(p) - G33(p) * DivOvThree(p) ) - G33(p) * TKEMulTwoOvThree(p) )
      end do
    end do
    !$omp end do
    !$omp end parallel
    return
  end subroutine atm_phy_tb_dgm_globalsmg_cal_grad

!OCL SERIAL  
  subroutine atm_phy_tb_dgm_globalsmg_cal_grad_qtrc( &
    dQTdx, dQTdy, dQTdz,                                        & ! (out)
    dRdx, dRdy, dRdz,                                           & ! (inout)
    QTRC, DDENS, DENS_hyd,                                      & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound, cal_grad_dens                                     ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) :: dQTdx(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: dQTdy(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: dQTdz(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: dRdx(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: dRdy(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: dRdz(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: QTRC(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)        
    type(SparseMat), intent(in) :: Dx, Dy, Dz
    type(SparseMat), intent(in) :: Sx, Sy, Sz
    type(SparseMat), intent(in) :: Lift
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)
    logical, intent(in) :: cal_grad_dens

    real(RP) :: G11(elem%Np), G12(elem%Np), G22(elem%Np)
    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DqDxi_(elem%Np,2)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,3)
    real(RP) :: del_flux_rho (elem%NfpTot,lmesh%Ne,3)    
    
    real(RP) :: DENS(elem%Np), RDENS(elem%Np)
    real(RP) :: RHOxQTRC(elem%Np)

    integer :: ke, ke2D
    integer :: p
    !--------------------------------------------------------------------

    call cal_del_flux_grad_qtrc( del_flux, del_flux_rho,                      & ! (out)
      QTRC, DDENS, DENS_hyd,                                                  & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound, cal_grad_dens )                                    ! (in)

    !$omp parallel private( ke, ke2D,   &
    !$omp Fx, Fy, Fz, LiftDelFlx,       &
    !$omp DqDxi_, RHOxQTRC, DENS, RDENS )

    ! Calculate gradient of density
    if ( cal_grad_dens ) then
      !$omp do
      do ke=lmesh%NeS, lmesh%NeE
        DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)

        call sparsemat_matmul( Dx, DENS, Fx )
        call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rho(:,ke,1), LiftDelFlx )
        dRdx(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:)
  
        call sparsemat_matmul( Dy, DENS, Fy )
        call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rho(:,ke,2), LiftDelFlx )
        dRdy(:,ke) = lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:)
  
        call sparsemat_matmul( Dz, DENS, Fz )
        call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rho(:,ke,3), LiftDelFlx )
        dRdz(:,ke) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)        
      end do
    end if

    ! Calculate gradient of tracer
    !$omp do
    do ke=lmesh%NeS, lmesh%NeE
      ke2d = lmesh%EMap3Dto2D(ke)
      G11(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,1)
      G12(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,2)
      G22(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,2,2)      

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      RDENS(:) = 1.0_RP / DENS(:)
      RHOxQTRC(:) = DENS(:) * QTRC(:,ke)

      !---
      call sparsemat_matmul( Dx, RHOxQTRC(:), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,1), LiftDelFlx )
      DqDxi_(:,1) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - QTRC(:,ke) * dRdx(:,ke) ) * RDENS(:)

      call sparsemat_matmul( Dy, RHOxQTRC(:), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,2), LiftDelFlx )
      DqDxi_(:,2) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - QTRC(:,ke) * dRdy(:,ke) ) * RDENS(:)

      call sparsemat_matmul( Dz, RHOxQTRC(:), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,3), LiftDelFlx )
      dQTdz(:,ke) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - QTRC(:,ke) * dRdz(:,ke) ) * RDENS(:)

      dQTdx(:,ke) = G11(:) * DqDxi_(:,1) + G12(:) * DqDxi_(:,2)
      dQTdy(:,ke) = G12(:) * DqDxi_(:,1) + G22(:) * DqDxi_(:,2)
    end do

    !$omp end parallel

    return
  end subroutine atm_phy_tb_dgm_globalsmg_cal_grad_qtrc

!OCL SERIAL  
  subroutine cal_del_flux_grad( del_flux_rho, del_flux_mom, del_flux_rhot,  & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, PT_,           & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound                         ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux_rho(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(out) ::  del_flux_mom(elem%NfpTot*lmesh%Ne,3,3)
    real(RP), intent(out) ::  del_flux_rhot(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: densM, densP
    real(RP) :: del
    real(RP) :: facx, facy, facz

    real(RP) :: MOMZ_P
    !------------------------------------------------------------------------
    
    !$omp parallel do private ( iM, iP,   &
    !$omp densM, densP,                   &
    !$omp del, facx, facy, facz, MOMZ_P   )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iP)

      if ( is_bound(i) ) then
        facx = 1.0_RP
        facy = 1.0_RP 
        facz = 1.0_RP
        MOMZ_P = - MOMZ_(iM)
      else
        ! facx = 1.0_RP - sign(1.0_RP,nx(i))
        ! facy = 1.0_RP - sign(1.0_RP,ny(i))
        ! facz = 1.0_RP - sign(1.0_RP,nz(i))
        facx = 1.0_RP
        facy = 1.0_RP
        facz = 1.0_RP
        MOMZ_P = MOMZ_(iP)
      end if

      del = 0.5_RP * ( densP - densM )
      del_flux_rho(i,1) = facx * del * nx(i)
      del_flux_rho(i,2) = facy * del * ny(i)
      del_flux_rho(i,3) = facz * del * nz(i)

      del = 0.5_RP * ( MOMX_(iP) - MOMX_(iM) )
      del_flux_mom(i,1,1) = facx * del * nx(i)
      del_flux_mom(i,2,1) = facy * del * ny(i)
      del_flux_mom(i,3,1) = facz * del * nz(i)

      del = 0.5_RP * ( MOMY_(iP) - MOMY_(iM) )
      del_flux_mom(i,1,2) = facx * del * nx(i)
      del_flux_mom(i,2,2) = facy * del * ny(i)
      del_flux_mom(i,3,2) = facz * del * nz(i)

      del = 0.5_RP * ( MOMZ_P - MOMZ_(iM) )
      del_flux_mom(i,1,3) = facx * del * nx(i)
      del_flux_mom(i,2,3) = facy * del * ny(i)
      del_flux_mom(i,3,3) = facz * del * nz(i)

      del = 0.5_RP * ( densP * PT_(iP) - densM * PT_(iM) )
      del_flux_rhot(i,1) = facx * del * nx(i)
      del_flux_rhot(i,2) = facy * del * ny(i)
      del_flux_rhot(i,3) = facz * del * nz(i)
    end do

    return
  end subroutine cal_del_flux_grad

!OCL SERIAL  
  subroutine cal_del_flux_grad_qtrc( del_flux, del_flux_rho, & ! (out)
    QTRC_, DDENS_, DENS_hyd_,                                & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound,         & ! (in)
    cal_grad_dens                                            ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(out) ::  del_flux_rho(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(in) ::  QTRC_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: cal_grad_dens
    
    integer :: i, iP, iM
    real(RP) :: del
    real(RP) :: facx, facy, facz
    real(RP) :: densM, densP

    !------------------------------------------------------------------------
    
    !$omp parallel do private ( iM, iP,       &
    !$omp del, facx, facy, facz, densM, densP )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      if ( is_bound(i) ) then
        facx = 1.0_RP
        facy = 1.0_RP 
        facz = 1.0_RP
      else
        ! facx = 1.0_RP - sign(1.0_RP,nx(i))
        ! facy = 1.0_RP - sign(1.0_RP,ny(i))
        ! facz = 1.0_RP - sign(1.0_RP,nz(i))
        facx = 1.0_RP
        facy = 1.0_RP
        facz = 1.0_RP
      end if

      densM = DDENS_(iM) + DENS_hyd_(iM)
      densP = DDENS_(iP) + DENS_hyd_(iP)

      if ( cal_grad_dens ) then
        del = 0.5_RP * ( densP - densM )
        ! del = 0.5_RP * ( sqrt_DENS_Kh_P - sqrt_DENS_Kh_M )
        del_flux_rho(i,1) = facx * del * nx(i)
        del_flux_rho(i,2) = facy * del * ny(i)
        del_flux_rho(i,3) = facz * del * nz(i)
      end if
      
      del = 0.5_RP * ( densP * QTRC_(iP) - densM * QTRC_(iM) )
      del_flux(i,1) = facx * del * nx(i)
      del_flux(i,2) = facy * del * ny(i)
      del_flux(i,3) = facz * del * nz(i)
    end do

    return
  end subroutine cal_del_flux_grad_qtrc

!OCL SERIAL  
  subroutine atm_phy_tb_dgm_globalsmg_cal_tend( &
    MOMX_t, MOMY_t, MOMZ_t, RHOT_t,                             & ! (out)
    T11, T12, T13, T21, T22, T23, T31, T32, T33,                & ! (in)
    dPTdx, dPTdy, dPTdz,                                        & ! (in)
    Nu, Kh,                                                     & ! (in)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,                        & ! (in)
    DENS_hyd, PRES_hyd,  PRES_, PT_,                            & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound                                                    ) ! (in)

    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSVec
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) :: MOMX_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_t(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: T11(elem%Np,lmesh%NeA), T12(elem%Np,lmesh%NeA), T13(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: T21(elem%Np,lmesh%NeA), T22(elem%Np,lmesh%NeA), T23(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: T31(elem%Np,lmesh%NeA), T32(elem%Np,lmesh%NeA), T33(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dPTdx(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dPTdy(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dPTdz(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Nu   (elem%Np,lmesh%NeA) ! Eddy viscosity
    real(RP), intent(in)  :: Kh   (elem%Np,lmesh%NeA) ! Eddy diffusivity
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_ (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_ (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_ (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PT_  (elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx, Dy, Dz
    type(SparseMat), intent(in) :: Sx, Sy, Sz
    type(SparseMat), intent(in) :: Lift
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)

    integer :: ke, ke2d

    real(RP) :: X2D(elem2D%Np,lmesh2D%Ne), Y2D(elem2D%Np,lmesh2D%Ne)
    real(RP) :: X(elem%Np), Y(elem%Np), twoOVdel2(elem%Np)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: GsqrtDENS(elem%Np), RHOT(elem%Np)
    real(RP) :: del_flux_mom(elem%NfpTot,lmesh%Ne,3)
    real(RP) :: del_flux_rhot(elem%NfpTot,lmesh%Ne)
    !--------------------------------------------------------------------

    call cal_del_flux( del_flux_mom, del_flux_rhot,                           & ! (out)
      T11, T12, T13, T21, T22, T23, T31, T32, T33,                            & ! (in)
      dPTdx, dPTdy, dPTdz,                                                    & ! (in)
      Nu, Kh,                                                                 & ! (in)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                & ! (in)
      lmesh%Gsqrt,                                                            & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound )                                                   ! (in)

    !$omp parallel private( ke, ke2D, &
    !$omp Fx, Fy, Fz, LiftDelFlx,          &
    !$omp GsqrtDENS, X, Y, twoOVdel2       )

    !$omp do
    do ke2D = lmesh2D%NeS, lmesh2D%NeE
      X2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,1))
      Y2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,2))
    end do

    !$omp do
    do ke=lmesh%NeS, lmesh%NeE
      ke2d = lmesh%EMap3Dto2D(ke)
      X(:) = X2D(elem%IndexH2Dto3D,ke2d)
      Y(:) = Y2D(elem%IndexH2Dto3D,ke2d)
      twoOVdel2(:) = 2.0_RP / ( 1.0_RP + X(:)**2 + Y(:)**2 )

      GsqrtDENS(:) = lmesh%Gsqrt(:,ke) * ( DENS_hyd(:,ke) + DDENS_(:,ke) )

      ! MOMX
      call sparsemat_matmul( Dx, lmesh%Gsqrt(:,ke) * T11(:,ke), Fx )
      call sparsemat_matmul( Dy, lmesh%Gsqrt(:,ke) * T12(:,ke), Fy )
      call sparsemat_matmul( Dz, lmesh%Gsqrt(:,ke) * T13(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,1), LiftDelFlx )

      MOMX_t(:,ke) = ( lmesh%Escale(:,ke,1,1) * Fx(:)       &
                     + lmesh%Escale(:,ke,2,2) * Fy(:)       &
                     + lmesh%Escale(:,ke,3,3) * Fz(:)       &
                     + LiftDelFlx(:)  ) / lmesh%Gsqrt(:,ke) &
                     + twoOVdel2(:) * Y(:) *                                        &
                     ( X(:) * Y(:) * T11(:,ke) - (1.0_RP + Y(:)**2) * T12(:,ke) )                   
      
      ! MOMY
      call sparsemat_matmul( Dx, lmesh%Gsqrt(:,ke) * T21(:,ke), Fx )
      call sparsemat_matmul( Dy, lmesh%Gsqrt(:,ke) * T22(:,ke), Fy )
      call sparsemat_matmul( Dz, lmesh%Gsqrt(:,ke) * T23(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2), LiftDelFlx )

      MOMY_t(:,ke)  = ( lmesh%Escale(:,ke,1,1) * Fx(:)       &
                      + lmesh%Escale(:,ke,2,2) * Fy(:)       &
                      + lmesh%Escale(:,ke,3,3) * Fz(:)       &
                      + LiftDelFlx(:)  ) / lmesh%Gsqrt(:,ke) &
                      + twoOVdel2(:) * X(:) *                                      &
                        ( - (1.0_RP + X(:)**2) * T21(:,ke) + X(:) * Y(:) * T22(:,ke) )
 
      ! MOMZ
      call sparsemat_matmul( Dx, lmesh%Gsqrt(:,ke) * T31(:,ke), Fx )
      call sparsemat_matmul( Dy, lmesh%Gsqrt(:,ke) * T32(:,ke), Fy )
      call sparsemat_matmul( Dz, lmesh%Gsqrt(:,ke) * T33(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3), LiftDelFlx )

      MOMZ_t(:,ke)  = ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                      + lmesh%Escale(:,ke,2,2) * Fy(:) &
                      + lmesh%Escale(:,ke,3,3) * Fz(:) &
                      + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)

      ! RHOT
      call sparsemat_matmul( Dx, GsqrtDENS(:) * Kh(:,ke) * dPTdx(:,ke), Fx )
      call sparsemat_matmul( Dy, GsqrtDENS(:) * Kh(:,ke) * dPTdy(:,ke), Fy )
      call sparsemat_matmul( Dz, GsqrtDENS(:) * Kh(:,ke) * dPTdz(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke), LiftDelFlx )

      RHOT_t(:,ke)  = ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                      + lmesh%Escale(:,ke,2,2) * Fy(:) &
                      + lmesh%Escale(:,ke,3,3) * Fz(:) &
                      + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)

    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine atm_phy_tb_dgm_globalsmg_cal_tend


!OCL SERIAL  
  subroutine atm_phy_tb_dgm_globalsmg_cal_tend_qtrc( &
    RHOQ_t,                                                     & ! (out)
    dQTdx, dQTdy, dQTdz,                                        & ! (in)
    Kh, DDENS_,DENS_hyd,                                        & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound                                                    ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) :: RHOQ_t(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dQTdx(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dQTdy(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dQTdz(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Kh   (elem%Np,lmesh%NeA) ! Eddy diffusivity
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx, Dy, Dz
    type(SparseMat), intent(in) :: Sx, Sy, Sz
    type(SparseMat), intent(in) :: Lift
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)

    integer :: ke

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: GsqrtDENS(elem%Np), RHOT(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
    !--------------------------------------------------------------------

    call cal_del_flux_qtrc( del_flux,                                                      & ! (out)
      dQTdx, dQTdy, dQTdz, Kh, DDENS_, DENS_hyd,                                           & ! (in)
      lmesh%Gsqrt, lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem, is_bound )                                      ! (in)

    !$omp parallel do private( &
    !$omp Fx, Fy, Fz, LiftDelFlx,            &
    !$omp GsqrtDENS                          )
    do ke=lmesh%NeS, lmesh%NeE
      GsqrtDENS(:) = lmesh%Gsqrt(:,ke) * ( DENS_hyd(:,ke) + DDENS_(:,ke) )

      ! RHOQ
      call sparsemat_matmul( Dx, GsqrtDENS(:) * Kh(:,ke) * dQTdx(:,ke), Fx )
      call sparsemat_matmul( Dy, GsqrtDENS(:) * Kh(:,ke) * dQTdy(:,ke), Fy )
      call sparsemat_matmul( Dz, GsqrtDENS(:) * Kh(:,ke) * dQTdz(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx )

      RHOQ_t(:,ke) = ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                     + lmesh%Escale(:,ke,2,2) * Fy(:) &
                     + lmesh%Escale(:,ke,3,3) * Fz(:) &
                     + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
    end do

    return
  end subroutine atm_phy_tb_dgm_globalsmg_cal_tend_qtrc

!-- private --------------------------------------------------------

!OCL SERIAL  
  subroutine cal_del_flux( del_flux_mom, del_flux_rhot,                & ! (out)
    T11, T12, T13, T21, T22, T23, T31, T32, T33,                       & ! (in)
    dPTdx, dPTdy, dPTdz,                                               & ! (in)
    Nu, Kh,                                                            & ! (in)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,           & ! (in)
    Gsqrt, nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound             ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux_mom (elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(out) ::  del_flux_rhot(elem%NfpTot*lmesh%Ne) 
    real(RP), intent(in)  :: T11(elem%Np*lmesh%NeA), T12(elem%Np*lmesh%NeA), T13(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: T21(elem%Np*lmesh%NeA), T22(elem%Np*lmesh%NeA), T23(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: T31(elem%Np*lmesh%NeA), T32(elem%Np*lmesh%NeA), T33(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: dPTdx(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: dPTdy(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: dPTdz(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: Nu   (elem%Np*lmesh%NeA) ! Eddy viscosity
    real(RP), intent(in)  :: Kh   (elem%Np*lmesh%NeA) ! Eddy diffusivity
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: GsqrtDensM, GsqrtDensP
    real(RP) :: TauM_x, TauP_x
    real(RP) :: TauM_y, TauP_y
    real(RP) :: TauM_z, TauP_z
    real(RP) :: nx_, ny_, nz_
    !------------------------------------------------------------------------

    !$omp parallel do private( iM, iP, &
    !$omp GsqrtDensM, GsqrtDensP,                                        &
    !$omp TauM_x, TauP_x, TauM_y, TauP_y, TauM_z, TauP_z, nx_, ny_, nz_  )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      GsqrtDensM = Gsqrt(iM) * ( DDENS_(iM) + DENS_hyd(iM) )
      GsqrtDensP = Gsqrt(iP) * ( DDENS_(iP) + DENS_hyd(iP) )

      if ( iP > elem%Np * lmesh%Ne .and. abs(nz(i)) > EPS ) then ! Tentative implementation for the treatmnet of lower/upper boundary. 
        nx_ = nx(i)
        ny_ = ny(i)
        nz_ = nz(i)
      else
        ! nx_ = ( 1.0_RP + sign(1.0_RP,nx(i)) ) * nx(i)
        ! ny_ = ( 1.0_RP + sign(1.0_RP,ny(i)) ) * ny(i)
        ! nz_ = ( 1.0_RP + sign(1.0_RP,nz(i)) ) * nz(i)
        nx_ = nx(i)
        ny_ = ny(i)
        nz_ = nz(i)
      end if

      TauM_x = Gsqrt(iM) * ( T11(iM) * nx_ + T12(iM) * ny_ + T13(iM) * nz_ )
      TauP_x = Gsqrt(iP) * ( T11(iP) * nx_ + T12(iP) * ny_ + T13(iP) * nz_ )

      TauM_y = Gsqrt(iM) * ( T21(iM) * nx_ + T22(iM) * ny_ + T23(iM) * nz_ )
      TauP_y = Gsqrt(iP) * ( T21(iP) * nx_ + T22(iP) * ny_ + T23(iP) * nz_ )

      TauM_z = Gsqrt(iM) * ( T31(iM) * nx_ + T32(iM) * ny_ + T33(iM) * nz_ )
      TauP_z = Gsqrt(iP) * ( T31(iP) * nx_ + T32(iP) * ny_ + T33(iP) * nz_ )

      if ( is_bound(i) )  then
        del_flux_mom(i,1) = - TauM_x
        del_flux_mom(i,2) = - TauM_y
        del_flux_mom(i,3) = 0.5_RP * ( TauP_z - TauM_z )
        del_flux_rhot(i)  = - GsqrtDensM * Kh(iM) * ( dPTdx(iM) * nx(i) + dPTdy(iM) * ny(i) + dPTdz(iM) * nz(i) )
      else        
        del_flux_mom(i,1) = 0.5_RP * ( TauP_x - TauM_x )
        del_flux_mom(i,2) = 0.5_RP * ( TauP_y - TauM_y )
        del_flux_mom(i,3) = 0.5_RP * ( TauP_z - TauM_z )
        del_flux_rhot(i)  = 0.5_RP * ( GsqrtDensP * Kh(iP) * ( dPTdx(iP) * nx_ + dPTdy(iP) * ny_ + dPTdz(iP) * nz_ ) &
                                     - GsqrtDensM * Kh(iM) * ( dPTdx(iM) * nx_ + dPTdy(iM) * ny_ + dPTdz(iM) * nz_ ) )
      end if
    end do

    return
  end subroutine cal_del_flux

!OCL SERIAL  
  subroutine cal_del_flux_qtrc( del_flux,                    & ! (out)
    dQTdx, dQTdy, dQTdz,                                     & ! (in)
    Kh, DDENS_, DENS_hyd,                                    & ! (in)
    Gsqrt, nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound   ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne) 
    real(RP), intent(in)  :: dQTdx(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: dQTdy(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: dQTdz(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: Kh   (elem%Np*lmesh%NeA) ! Eddy diffusivity
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: GsqrtDensM, GsqrtDensP
    real(RP) :: nx_, ny_, nz_
    !------------------------------------------------------------------------

    !$omp parallel do private( iM, iP, &
    !$omp GsqrtDensM, GsqrtDensP, nx_, ny_, nz_  )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      GsqrtDensM = Gsqrt(iM) * ( DDENS_(iM) + DENS_hyd(iM) )
      GsqrtDensP = Gsqrt(iP) * ( DDENS_(iP) + DENS_hyd(iP) )

      if ( iP > elem%Np * lmesh%Ne .and. abs(nz(i)) > EPS ) then ! Tentative implementation for the treatmnet of lower/upper boundary. 
        nx_ = nx(i)
        ny_ = ny(i)
        nz_ = nz(i)
      else
        ! nx_ = ( 1.0_RP + sign(1.0_RP,nx(i)) ) * nx(i)
        ! ny_ = ( 1.0_RP + sign(1.0_RP,ny(i)) ) * ny(i)
        ! nz_ = ( 1.0_RP + sign(1.0_RP,nz(i)) ) * nz(i)
        nx_ = nx(i)
        ny_ = ny(i)
        nz_ = nz(i)
      end if

      if ( is_bound(i) )  then
        del_flux(i)  = - GsqrtDensM * Kh(iM) * ( dQTdx(iM) * nx(i) + dQTdy(iM) * ny(i) + dQTdz(iM) * nz(i) )
      else        
        del_flux(i)  = 0.5_RP * ( GsqrtDensP * Kh(iP) * ( dQTdx(iP) * nx_ + dQTdy(iP) * ny_ + dQTdz(iP) * nz_ ) &
                                - GsqrtDensM * Kh(iM) * ( dQTdx(iM) * nx_ + dQTdy(iM) * ny_ + dQTdz(iM) * nz_ ) )
      end if
    end do

    return
  end subroutine cal_del_flux_qtrc

end module scale_atm_phy_tb_dgm_globalsmg
