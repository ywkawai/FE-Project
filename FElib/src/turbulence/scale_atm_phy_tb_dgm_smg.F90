!> module FElib / Atmosphere / Physics turbulence
!!
!! @par Description
!!      Sub-grid scale turbulnce process
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
module scale_atm_phy_tb_dgm_smg
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
  public :: atm_phy_tb_dgm_smg_Init
  public :: atm_phy_tb_dgm_smg_Final
  public :: atm_phy_tb_dgm_smg_cal_grad
  public :: atm_phy_tb_dgm_smg_cal_grad_qtrc
  public :: atm_phy_tb_dgm_smg_cal_tend
  public :: atm_phy_tb_dgm_smg_cal_tend_qtrc

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
  subroutine atm_phy_tb_dgm_smg_Init( mesh )
    implicit none    
    class(MeshBase3D), intent(in) :: mesh

    logical  :: consistent_tke = .true.

    namelist / PARAM_ATMOS_PHY_TB_DGM_SMG / &
      Cs,                                   &
      NU_MAX,                               &
      filter_fac,                           &
      consistent_tke
    
    integer :: ierr
    !--------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_dgm_smg_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_TB_dgm_smg_setup",*) 'Smagorinsky-type Eddy Viscocity Model'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_DGM_SMG,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_TB_dgm_smg_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_TB_dgm_smg_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_TB_DGM_SMG. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_TB_DGM_SMG)

    RPrN         = 1.0_RP / PrN
    RRiC         = 1.0_RP / RiC
    OnemPrNovRiC = ( 1.0_RP - PrN ) * RRiC

    if ( consistent_tke ) then
      tke_fac = 1.0_RP
    else
      tke_fac = 0.0_RP
    end if

    return
  end subroutine atm_phy_tb_dgm_smg_Init

!OCL SERIAL
  subroutine atm_phy_tb_dgm_smg_Final()
    implicit none
    !--------------------------------------------------------------------

    return
  end subroutine atm_phy_tb_dgm_smg_Final

!> Calculate parameterized stress tensor and eddy heat flux with turbulent model
!!
!OCL SERIAL  
  subroutine atm_phy_tb_dgm_smg_cal_grad( &
    T11, T12, T13, T21, T22, T23, T31, T32, T33,                & ! (out)
    DF1, DF2, DF3,                                              & ! (out)
    TKE, Nu, Kh,                                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,    & ! (in)
    PRES, PT,                                                   & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound                                                    ) ! (in)

    use scale_atm_phy_tb_dgm_common, only: &
      atm_phy_tb_dgm_common_calc_lambda
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) :: T11(elem%Np,lmesh%NeA)      !< (1,1) component of stress tensor
    real(RP), intent(out) :: T12(elem%Np,lmesh%NeA)      !< (1,2) component of stress tensor
    real(RP), intent(out) :: T13(elem%Np,lmesh%NeA)      !< (1,3) component of stress tensor
    real(RP), intent(out) :: T21(elem%Np,lmesh%NeA)      !< (2,1) component of stress tensor
    real(RP), intent(out) :: T22(elem%Np,lmesh%NeA)      !< (2,2) component of stress tensor
    real(RP), intent(out) :: T23(elem%Np,lmesh%NeA)      !< (2,3) component of stress tensor
    real(RP), intent(out) :: T31(elem%Np,lmesh%NeA)      !< (3,1) component of stress tensor
    real(RP), intent(out) :: T32(elem%Np,lmesh%NeA)      !< (3,2) component of stress tensor
    real(RP), intent(out) :: T33(elem%Np,lmesh%NeA)      !< (3,3) component of stress tensor
    real(RP), intent(out)  :: DF1(elem%Np,lmesh%NeA)     !< Diffusive heat flux in x1 direction / density
    real(RP), intent(out)  :: DF2(elem%Np,lmesh%NeA)     !< Diffusive heat flux in x2 direction / density
    real(RP), intent(out)  :: DF3(elem%Np,lmesh%NeA)     !< Diffusive heat flux in x3 direction / density
    real(RP), intent(out) :: TKE(elem%Np,lmesh%NeA)      !< Parameterized turbulent kinetic energy
    real(RP), intent(out) :: Nu(elem%Np,lmesh%NeA)       !< Eddy viscosity
    real(RP), intent(out) :: Kh(elem%Np,lmesh%NeA)       !< Eddy diffusivity
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)   !< Density perturbation
    real(RP), intent(in)  :: MOMX_ (elem%Np,lmesh%NeA)   !< Momentum in x1 direction
    real(RP), intent(in)  :: MOMY_ (elem%Np,lmesh%NeA)   !< Momentum in x2 direction
    real(RP), intent(in)  :: MOMZ_ (elem%Np,lmesh%NeA)   !< Momentum in x3 direction
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)   !< Density x potential temperature perturbation
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA) !< Reference pressure in hydrostatic balance
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA) !< Reference density in hydrostatic balance
    real(RP), intent(in)  :: PRES(elem%Np,lmesh%NeA)     !< Pressure
    real(RP), intent(in)  :: PT(elem%Np,lmesh%NeA)       !< Potential temperature
    type(SparseMat), intent(in) :: Dx, Dy, Dz             !< Differential matrix managed by sparse matrix type
    type(SparseMat), intent(in) :: Sx, Sy, Sz             !< Stiffness matrix managed by sparse matrix type
    type(SparseMat), intent(in) :: Lift                   !< Lifting matrix managed by sparse matrix type
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne) !< Flag whether nodes are located at domain boundaries

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DENS(elem%Np), RDENS(elem%Np), RHOT(elem%Np), Q(elem%Np)
    real(RP) :: DdensDxi(elem%Np,3)
    real(RP) :: DVelDxi(elem%Np,3,3)
    real(RP) :: del_flux_rho (elem%NfpTot,lmesh%Ne,3)
    real(RP) :: del_flux_mom (elem%NfpTot,lmesh%Ne,3,3)
    real(RP) :: del_flux_rhot(elem%NfpTot,lmesh%Ne,3)

    real(RP) :: S11(elem%Np), S12(elem%Np), S22(elem%Np), S23(elem%Np), S31(elem%Np), S33(elem%Np)
    real(RP) :: SkkOvThree
    real(RP) :: TKEMulTwoOvThree
    real(RP) :: coef

    real(RP) :: Ri ! local gradient Richardson number
    real(RP) :: S2 ! (2SijSij)^1/2
    real(RP) :: fm ! factor in eddy viscosity which represents the stability dependence of 
                   ! the Brown et al (1994)'s subgrid model 
    real(RP) :: Pr ! Parandtl number (=Nu/Kh= fm/fh)

    real(RP) :: lambda  (elem%Np,lmesh%Ne) ! basic mixing length
    real(RP) :: lambda_r(elem%Np)          ! characteristic subgrid length scale 
    real(RP) :: E(elem%Np)  ! subgrid kinetic energy 
    real(RP) :: C1(elem%Np) ! factor in the relation with energy disspation rate, lambda_r, and E

    integer :: ke
    integer :: p
    !--------------------------------------------------------------------

    call cal_del_flux_grad( del_flux_rho, del_flux_mom, del_flux_rhot,        & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, PT,            & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound )                                                   ! (in)

    call atm_phy_tb_dgm_common_calc_lambda( lambda, & ! (out)
      Cs, filter_fac, lmesh, elem, lmesh2D, elem2D  ) ! (in)
  
    !$omp parallel do private( &
    !$omp Fx, Fy, Fz, LiftDelFlx,                  &
    !$omp DENS, RHOT, RDENS, Q, DdensDxi, DVelDxi, &
    !$omp S11, S12, S22, S23, S31, S33,            &
    !$omp coef, SkkOvThree, TKEMulTwoOvThree,      &
    !$omp p, Ri, S2, fm, Pr, lambda_r, E, C1       )
    do ke=lmesh%NeS, lmesh%NeE
      !---
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
      DVelDxi(:,1,1) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) ) * RDENS(:)

      call sparsemat_matmul( Dy, MOMX_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2,1), LiftDelFlx )
      DVelDxi(:,2,1) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) ) * RDENS(:)

      call sparsemat_matmul( Dz, MOMX_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3,1), LiftDelFlx )
      DVelDxi(:,3,1) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)

      ! gradient of v
      Q(:) = MOMY_(:,ke) * RDENS(:)

      call sparsemat_matmul( Dx, MOMY_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,1,2), LiftDelFlx )
      DVelDxi(:,1,2) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) ) * RDENS(:)

      call sparsemat_matmul( Dy, MOMY_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2,2), LiftDelFlx )
      DVelDxi(:,2,2) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) ) * RDENS(:)

      call sparsemat_matmul( Dz, MOMY_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3,2), LiftDelFlx )
      DVelDxi(:,3,2) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)

      ! gradient of w
      Q(:) = MOMZ_(:,ke) * RDENS(:)

      call sparsemat_matmul( Dx, MOMZ_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,1,3), LiftDelFlx )
      DVelDxi(:,1,3) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) ) * RDENS(:)

      call sparsemat_matmul( Dy, MOMZ_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2,3), LiftDelFlx )
      DVelDxi(:,2,3) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) ) * RDENS(:)

      call sparsemat_matmul( Dz, MOMZ_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3,3), LiftDelFlx )
      DVelDxi(:,3,3) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)

      ! gradient of pt
      Q(:) = RHOT(:) * RDENS(:)

      call sparsemat_matmul( Dx, RHOT, Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke,1), LiftDelFlx )
      DF1(:,ke) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) ) * RDENS(:)

      call sparsemat_matmul( Dy, RHOT, Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke,2), LiftDelFlx )
      DF2(:,ke) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) ) * RDENS(:)

      call sparsemat_matmul( Dz, RHOT, Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke,3), LiftDelFlx )
      DF3(:,ke) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)


      ! Calculate the component of strain velocity tensor
      S11(:) = DVelDxi(:,1,1)
      S12(:) = 0.5_RP * ( DVelDxi(:,1,2) + DVelDxi(:,2,1) )
      S22(:) = DVelDxi(:,2,2)
      S23(:) = 0.5_RP * ( DVelDxi(:,2,3) + DVelDxi(:,3,2) )
      S31(:) = 0.5_RP * ( DVelDxi(:,1,3) + DVelDxi(:,3,1) )
      S33(:) = DVelDxi(:,3,3)

      ! Caclulate eddy viscosity & eddy diffusivity
      
      do p=1, elem%Np
        S2 = 2.0_RP * ( S11(p)**2 + S22(p)**2 + S33(p)**2 ) &
           + 4.0_RP * ( S31(p)**2 + S12(p)**2 + S23(p)**2 )
        
        Ri = Grav / PT(p,ke) * DF3(p,ke) / max( S2, EPS )

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
      
!      Nu(:,ke) = 0.0_RP

      ! if ( backscatter ) then
      ! else 
        E (:) = Nu(:,ke)**3 / ( lambda_r(:)**4 + EPS )
        C1(:) = C1o
      ! end if

      ! TKE
      TKE(:,ke) = ( E(:) * lambda_r(:) / C1(:) )**twoOverThree

      !---

      do p=1, elem%Np
        TKEMulTwoOvThree = twoOverThree * TKE(p,ke) * tke_fac
        SkkOvThree = ( S11(p) + S22(p) + S33(p) ) * OneOverThree
        coef = 2.0_RP * Nu(p,ke)

        T11(p,ke) = DENS(p) * ( coef * ( S11(p) - SkkOvThree ) - TKEMulTwoOvThree )
        T12(p,ke) = DENS(p) * coef * S12(p)
        T13(p,ke) = DENS(p) * coef * S31(p)

        T21(p,ke) = DENS(p) * coef * S12(p)
        T22(p,ke) = DENS(p) * ( coef * ( S22(p) - SkkOvThree ) - TKEMulTwoOvThree )
        T23(p,ke) = DENS(p) * coef * S23(p)

        T31(p,ke) = DENS(p) * coef * S31(p)
        T32(p,ke) = DENS(p) * coef * S23(p)
        T33(p,ke) = DENS(p) * ( coef * ( S33(p) - SkkOvThree ) - TKEMulTwoOvThree )
      end do

      DF1(:,ke) = Kh(:,ke) * DF1(:,ke)
      DF2(:,ke) = Kh(:,ke) * DF2(:,ke)
      DF3(:,ke) = Kh(:,ke) * DF3(:,ke)
    end do
  
    return
  end subroutine atm_phy_tb_dgm_smg_cal_grad

!> Calculate parameterized diffusive mass flux of tracer with turbulent model
!!
!OCL SERIAL  
  subroutine atm_phy_tb_dgm_smg_cal_grad_qtrc( &
    DFQ1, DFQ2, DFQ3,                                           & ! (out)
    dRdx, dRdy, dRdz,                                           & ! (inout)
    Kh, QTRC, DDENS, DENS_hyd,                                  & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound, cal_grad_dens                                     ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) :: DFQ1(elem%Np,lmesh%NeA)      !< Diffusive mass flux in x1 direction / density (Kh dq/dx1)
    real(RP), intent(out) :: DFQ2(elem%Np,lmesh%NeA)      !< Diffusive mass flux in x2 direction / density (Kh dq/dx2)
    real(RP), intent(out) :: DFQ3(elem%Np,lmesh%NeA)      !< Diffusive mass flux in x3 direction / density (Kh dq/dx3)
    real(RP), intent(inout) :: dRdx(elem%Np,lmesh%NeA)    !< Spatial gradient of density in x1 direction
    real(RP), intent(inout) :: dRdy(elem%Np,lmesh%NeA)    !< Spatial gradient of density in x2 direction
    real(RP), intent(inout) :: dRdz(elem%Np,lmesh%NeA)    !< Spatial gradient of density in x3 direction
    real(RP), intent(in)  :: Kh(elem%Np,lmesh%NeA)        !< Eddy diffusivity
    real(RP), intent(in)  :: QTRC(elem%Np,lmesh%NeA)      !< Mass faction of tracer
    real(RP), intent(in)  :: DDENS(elem%Np,lmesh%NeA)     !< Density perturbation
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)  !< Reference desity in hydrostatic state
    type(SparseMat), intent(in) :: Dx, Dy, Dz             !< Differential matrix managed by sparse matrix type
    type(SparseMat), intent(in) :: Sx, Sy, Sz             !< Stiffness matrix managed by sparse matrix type
    type(SparseMat), intent(in) :: Lift                   !< Lifting matrix managed by sparse matrix type
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne) !< Flag whether nodes are located at domain boundaries
    logical, intent(in) :: cal_grad_dens                  !< Flag whether spatial gradients of density are calcuated

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,3)
    real(RP) :: del_flux_rho (elem%NfpTot,lmesh%Ne,3)    
    
    real(RP) :: DENS(elem%Np), RDENS(elem%Np)
    real(RP) :: RHOxQTRC(elem%Np)

    integer :: ke
    integer :: p
    !--------------------------------------------------------------------

    call cal_del_flux_grad_qtrc( del_flux, del_flux_rho,                      & ! (out)
      QTRC, DDENS, DENS_hyd,                                                  & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound, cal_grad_dens )                                    ! (in)

    !$omp parallel private( ke, Fx, Fy, Fz, LiftDelFlx, RHOxQTRC, DENS, RDENS )

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
      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      RDENS(:) = 1.0_RP / DENS(:)
      RHOxQTRC(:) = DENS(:) * QTRC(:,ke)

      !---
      call sparsemat_matmul( Dx, RHOxQTRC(:), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,1), LiftDelFlx )
      DFQ1(:,ke) = Kh(:,ke) * ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - QTRC(:,ke) * dRdx(:,ke) ) * RDENS(:)

      call sparsemat_matmul( Dy, RHOxQTRC(:), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,2), LiftDelFlx )
      DFQ2(:,ke) = Kh(:,ke) * ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - QTRC(:,ke) * dRdy(:,ke) ) * RDENS(:)

      call sparsemat_matmul( Dz, RHOxQTRC(:), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,3), LiftDelFlx )
      DFQ3(:,ke) = Kh(:,ke) * ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - QTRC(:,ke) * dRdz(:,ke) ) * RDENS(:)
    end do

    !$omp end parallel

    return
  end subroutine atm_phy_tb_dgm_smg_cal_grad_qtrc

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
    !------------------------------------------------------------------------
    
    !$omp parallel do private ( iM, iP,   &
    !$omp densM, densP,                   &
    !$omp del, facx, facy, facz           )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iP)

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

      del = 0.5_RP * ( MOMZ_(iP) - MOMZ_(iM) )
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

!> Calculate tendecies with turbulent model
!!
!OCL SERIAL  
  subroutine atm_phy_tb_dgm_smg_cal_tend( &
    MOMX_t, MOMY_t, MOMZ_t, RHOT_t,                             & ! (out)
    T11, T12, T13, T21, T22, T23, T31, T32, T33,                & ! (in)
    DF1, DF2, DF3,                                              & ! (in)
    Nu, Kh,                                                     & ! (in)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,                        & ! (in)
    DENS_hyd, PRES_hyd,  PRES_, PT_,                            & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound                                                    ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) :: MOMX_t(elem%Np,lmesh%NeA)   !< Tendency of momentum in x1 direction with turbulent model
    real(RP), intent(out) :: MOMY_t(elem%Np,lmesh%NeA)   !< Tendency of momentum in x2 direction with turbulent model
    real(RP), intent(out) :: MOMZ_t(elem%Np,lmesh%NeA)   !< Tendency of momentum in x3 direction with turbulent model
    real(RP), intent(out) :: RHOT_t(elem%Np,lmesh%NeA)   !< Tendency of density x potential temperature with turbulent model
    real(RP), intent(in) :: T11(elem%Np,lmesh%NeA)       !< (1,1) component of stress tensor
    real(RP), intent(in) :: T12(elem%Np,lmesh%NeA)       !< (1,2) component of stress tensor
    real(RP), intent(in) :: T13(elem%Np,lmesh%NeA)       !< (1,3) component of stress tensor
    real(RP), intent(in) :: T21(elem%Np,lmesh%NeA)       !< (2,1) component of stress tensor
    real(RP), intent(in) :: T22(elem%Np,lmesh%NeA)       !< (2,2) component of stress tensor
    real(RP), intent(in) :: T23(elem%Np,lmesh%NeA)       !< (2,3) component of stress tensor
    real(RP), intent(in) :: T31(elem%Np,lmesh%NeA)       !< (3,1) component of stress tensor
    real(RP), intent(in) :: T32(elem%Np,lmesh%NeA)       !< (3,2) component of stress tensor
    real(RP), intent(in) :: T33(elem%Np,lmesh%NeA)       !< (3,3) component of stress tensor
    real(RP), intent(in)  :: DF1(elem%Np,lmesh%NeA)      !< Diffusive heat flux in x1 direction / density
    real(RP), intent(in)  :: DF2(elem%Np,lmesh%NeA)      !< Diffusive heat flux in x2 direction / density
    real(RP), intent(in)  :: DF3(elem%Np,lmesh%NeA)      !< Diffusive heat flux in x3 direction / density
    real(RP), intent(in)  :: Nu   (elem%Np,lmesh%NeA)    !< Eddy viscosity
    real(RP), intent(in)  :: Kh   (elem%Np,lmesh%NeA)    !< Eddy diffusivity
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)   !< Density perturbation
    real(RP), intent(in)  :: MOMX_ (elem%Np,lmesh%NeA)   !< Momentum in x1 direction
    real(RP), intent(in)  :: MOMY_ (elem%Np,lmesh%NeA)   !< Momentum in x2 direction
    real(RP), intent(in)  :: MOMZ_ (elem%Np,lmesh%NeA)   !< Momentum in x3 direction
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)   !< Density x potential temperature perturbation
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA) !< Reference pressure in hydrostatic balance
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA) !< Reference density in hydrostatic balance
    real(RP), intent(in)  :: PRES_(elem%Np,lmesh%NeA)    !< Pressure
    real(RP), intent(in)  :: PT_  (elem%Np,lmesh%NeA)    !< Potential temperature
    type(SparseMat), intent(in) :: Dx, Dy, Dz             !< Differential matrix managed by sparse matrix type
    type(SparseMat), intent(in) :: Sx, Sy, Sz             !< Stiffness matrix managed by sparse matrix type
    type(SparseMat), intent(in) :: Lift                   !< Lifting matrix managed by sparse matrix type
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne) !< Flag whether nodes are located at domain boundaries

    integer :: ke

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DENS(elem%Np)
    real(RP) :: del_flux_mom(elem%NfpTot,lmesh%Ne,3)
    real(RP) :: del_flux_rhot(elem%NfpTot,lmesh%Ne)    
    !--------------------------------------------------------------------

    call cal_del_flux( del_flux_mom, del_flux_rhot,                           & ! (out)
      T11, T12, T13, T21, T22, T23, T31, T32, T33,                            & ! (in)
      DF1, DF2, DF3,                                                          & ! (in)
      Nu, Kh,                                                                 & ! (in)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem, is_bound   )                       ! (in)

    !$omp parallel do private( &
    !$omp Fx, Fy, Fz, LiftDelFlx,                    &
    !$omp DENS                                       )
    do ke=lmesh%NeS, lmesh%NeE
      DENS(:) = DENS_hyd(:,ke) + DDENS_(:,ke)

      ! MOMX
      call sparsemat_matmul( Dx, T11(:,ke), Fx )
      call sparsemat_matmul( Dy, T12(:,ke), Fy )
      call sparsemat_matmul( Dz, T13(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,1), LiftDelFlx )

      MOMX_t(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + lmesh%Escale(:,ke,2,2) * Fy(:) &
                   + lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      ! MOMY
      call sparsemat_matmul( Dx, T21(:,ke), Fx )
      call sparsemat_matmul( Dy, T22(:,ke), Fy )
      call sparsemat_matmul( Dz, T23(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2), LiftDelFlx )

      MOMY_t(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + lmesh%Escale(:,ke,2,2) * Fy(:) &
                   + lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      ! MOMZ
      call sparsemat_matmul( Dx, T31(:,ke), Fx )
      call sparsemat_matmul( Dy, T32(:,ke), Fy )
      call sparsemat_matmul( Dz, T33(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3), LiftDelFlx )

      MOMZ_t(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + lmesh%Escale(:,ke,2,2) * Fy(:) &
                   + lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      ! RHOT
      call sparsemat_matmul( Dx, DENS(:) * DF1(:,ke), Fx )
      call sparsemat_matmul( Dy, DENS(:) * DF2(:,ke), Fy )
      call sparsemat_matmul( Dz, DENS(:) * DF3(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke), LiftDelFlx )

      RHOT_t(:,ke)  = ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                      + lmesh%Escale(:,ke,2,2) * Fy(:) &
                      + lmesh%Escale(:,ke,3,3) * Fz(:) &
                      + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)

    end do

    return
  end subroutine atm_phy_tb_dgm_smg_cal_tend

!> Calculate tendecies of tracer density with turbulent model
!!
!OCL SERIAL  
  subroutine atm_phy_tb_dgm_smg_cal_tend_qtrc( &
    RHOQ_t,                                                     & ! (out)
    DFQ1, DFQ2, DFQ3,                                           & ! (in)
    Kh, DDENS_,DENS_hyd,                                        & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound                                                    ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) :: RHOQ_t(elem%Np,lmesh%NeA)    !< Tendency of tracer mass fraction
    real(RP), intent(in)  :: DFQ1(elem%Np,lmesh%NeA)      !< Diffusive mass flux in x1 direction / density (Kh dq/dx1)
    real(RP), intent(in)  :: DFQ2(elem%Np,lmesh%NeA)      !< Diffusive mass flux in x2 direction / density (Kh dq/dx2)
    real(RP), intent(in)  :: DFQ3(elem%Np,lmesh%NeA)      !< Diffusive mass flux in x3 direction / density (Kh dq/dx3)
    real(RP), intent(in)  :: Kh   (elem%Np,lmesh%NeA)     !< Eddy diffusivity
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)    !< Density perturbation
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)  !< Reference pressure in hydrostatic state
    type(SparseMat), intent(in) :: Dx, Dy, Dz             !< Differential matrix managed by sparse matrix type
    type(SparseMat), intent(in) :: Sx, Sy, Sz             !< Stiffness matrix managed by sparse matrix type
    type(SparseMat), intent(in) :: Lift                   !< Lifting matrix managed by sparse matrix type
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne) !< Flag whether nodes are located at domain boundaries

    integer :: ke

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DENS(elem%Np), RHOT(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
    !--------------------------------------------------------------------

    call cal_del_flux_qtrc( del_flux,                                         & ! (out)
      DFQ1, DFQ2, DFQ3,                                                       & ! (in)
      Kh, DDENS_, DENS_hyd,                                                   & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound )                                                   ! (in)

    !$omp parallel do private( &
    !$omp Fx, Fy, Fz, LiftDelFlx,            &
    !$omp DENS                               )
    do ke=lmesh%NeS, lmesh%NeE
      DENS(:) = DENS_hyd(:,ke) + DDENS_(:,ke)

      ! RHOQ
      call sparsemat_matmul( Dx, DENS(:) * DFQ1(:,ke), Fx )
      call sparsemat_matmul( Dy, DENS(:) * DFQ2(:,ke), Fy )
      call sparsemat_matmul( Dz, DENS(:) * DFQ3(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx )

      RHOQ_t(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + lmesh%Escale(:,ke,2,2) * Fy(:) &
                   + lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)
    end do

    return
  end subroutine atm_phy_tb_dgm_smg_cal_tend_qtrc

!-- private --------------------------------------------------------

!OCL SERIAL  
  subroutine cal_del_flux( del_flux_mom, del_flux_rhot,                & ! (out)
    T11, T12, T13, T21, T22, T23, T31, T32, T33,                       & ! (in)
    DF1, DF2, DF3,                                                     & ! (in)
    Nu, Kh,                                                            & ! (in)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,           & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound                    ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux_mom (elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(out) ::  del_flux_rhot(elem%NfpTot*lmesh%Ne) 
    real(RP), intent(in)  :: T11(elem%Np*lmesh%NeA), T12(elem%Np*lmesh%NeA), T13(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: T21(elem%Np*lmesh%NeA), T22(elem%Np*lmesh%NeA), T23(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: T31(elem%Np*lmesh%NeA), T32(elem%Np*lmesh%NeA), T33(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: DF1(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: DF2(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: DF3(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: Nu   (elem%Np*lmesh%NeA) ! Eddy viscosity
    real(RP), intent(in)  :: Kh   (elem%Np*lmesh%NeA) ! Eddy diffusivity
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: densM, densP
    real(RP) :: TauM_x, TauP_x
    real(RP) :: TauM_y, TauP_y
    real(RP) :: TauM_z, TauP_z
    real(RP) :: nx_, ny_, nz_
    !------------------------------------------------------------------------

    !$omp parallel do private( iM, iP, &
    !$omp densM, densP, TauM_x, TauP_x, TauM_y, TauP_y, TauM_z, TauP_z, nx_, ny_, nz_  )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iP)

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

      TauM_x = T11(iM) * nx_ + T12(iM) * ny_ + T13(iM) * nz_
      TauP_x = T11(iP) * nx_ + T12(iP) * ny_ + T13(iP) * nz_

      TauM_y = T21(iM) * nx_ + T22(iM) * ny_ + T23(iM) * nz_
      TauP_y = T21(iP) * nx_ + T22(iP) * ny_ + T23(iP) * nz_

      TauM_z = T31(iM) * nx_ + T32(iM) * ny_ + T33(iM) * nz_
      TauP_z = T31(iP) * nx_ + T32(iP) * ny_ + T33(iP) * nz_

      del_flux_mom(i,1) = 0.5_RP * ( TauP_x - TauM_x )
      del_flux_mom(i,2) = 0.5_RP * ( TauP_y - TauM_y )
      del_flux_mom(i,3) = 0.5_RP * ( TauP_z - TauM_z )
      del_flux_rhot(i)  = 0.5_RP * ( densP * ( DF1(iP) * nx_ + DF2(iP) * ny_ + DF3(iP) * nz_ ) &
                                   - densM * ( DF1(iM) * nx_ + DF2(iM) * ny_ + DF3(iM) * nz_ ) )
    end do

    return
  end subroutine cal_del_flux

!OCL SERIAL  
  subroutine cal_del_flux_qtrc( del_flux,             & ! (out)
    DFQ1, DFQ2, DFQ3,                                 & ! (in)
    Kh, DDENS_, DENS_hyd,                             & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound   ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne) 
    real(RP), intent(in)  :: DFQ1(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: DFQ2(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: DFQ3(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: Kh   (elem%Np*lmesh%NeA) ! Eddy diffusivity
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: densM, densP
    real(RP) :: nx_, ny_, nz_
    !------------------------------------------------------------------------

    !$omp parallel do private( iM, iP, &
    !$omp densM, densP, nx_, ny_, nz_  )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iP)

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

      del_flux(i)  = 0.5_RP * ( densP * ( DFQ1(iP) * nx_ + DFQ2(iP) * ny_ + DFQ3(iP) * nz_ ) &
                              - densM * ( DFQ1(iM) * nx_ + DFQ2(iM) * ny_ + DFQ3(iM) * nz_ ) )
    end do

    return
  end subroutine cal_del_flux_qtrc

end module scale_atm_phy_tb_dgm_smg
