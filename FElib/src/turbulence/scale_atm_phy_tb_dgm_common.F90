!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Physics turbulence / Common
!!
!! @par Description
!!      A coomon modules for atmospheric turbulent parameterization
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_tb_dgm_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    EPS  => CONST_EPS
  
  use scale_sparsemat
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D

  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_phy_tb_dgm_common_get_varinfo
  public :: atm_phy_tb_dgm_common_calc_lambda

  public :: atm_phy_tb_dgm_common_cal_grad_qtrc
  public :: atm_phy_tb_dgm_common_cal_tend
  public :: atm_phy_tb_dgm_common_cal_tend_qtrc
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T13_ID      = 1
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T23_ID      = 2
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T33_ID      = 3
  integer, public, parameter :: ATMOS_PHY_TB_AUX_DIFFFLX3_ID = 4
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T31_ID      = 5
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T32_ID      = 6
  integer, public, parameter :: ATMOS_PHY_TB_AUX_DIFFFLX1_ID = 7
  integer, public, parameter :: ATMOS_PHY_TB_AUX_DIFFFLX2_ID = 8
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T11_ID      = 9
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T12_ID      = 10
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T21_ID      = 11
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T22_ID      = 12
  integer, public, parameter :: ATMOS_PHY_TB_AUX_SCALAR_NUM  = 4
  integer, public, parameter :: ATMOS_PHY_TB_AUX_HVEC_NUM    = 2
  integer, public, parameter :: ATMOS_PHY_TB_AUX_HTENSOR_NUM = 1
  integer, public, parameter :: ATMOS_PHY_TB_AUX_NUM         = 12
  
  integer, public, parameter :: ATMOS_PHY_TB_DIAG_TKE_ID     = 1
  integer, public, parameter :: ATMOS_PHY_TB_DIAG_NU_ID      = 2
  integer, public, parameter :: ATMOS_PHY_TB_DIAG_KH_ID      = 3
  integer, public, parameter :: ATMOS_PHY_TB_DIAG_NUM        = 3

  integer, public, parameter :: ATMOS_PHY_TB_MOMX_t_ID  = 1
  integer, public, parameter :: ATMOS_PHY_TB_MOMY_t_ID  = 2
  integer, public, parameter :: ATMOS_PHY_TB_MOMZ_t_ID  = 3
  integer, public, parameter :: ATMOS_PHY_TB_RHOT_t_ID  = 4
  integer, public, parameter :: ATMOS_PHY_TB_TENDS_NUM1 = 4
  !-

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------
  
  private :: cal_bnd_flux
  private :: cal_bnd_flux_grad_qtrc
  private :: cal_bnd_flux_qtrc

  private :: fact, p1, p2, p3
  
contains
!OCL SERIAL
  subroutine atm_phy_tb_dgm_common_get_varinfo( &
    auxvar_info, diagvar_info, tend_info        )

    implicit none

    type(VariableInfo), intent(out) :: auxvar_info(ATMOS_PHY_TB_AUX_NUM)
    type(VariableInfo), intent(out) :: diagvar_info(ATMOS_PHY_TB_DIAG_NUM)
    type(VariableInfo), intent(out) :: tend_info(ATMOS_PHY_TB_TENDS_NUM1)

    type(VariableInfo) :: ATMOS_PHY_TB_AUX_VINFO(ATMOS_PHY_TB_AUX_NUM)
    DATA ATMOS_PHY_TB_AUX_VINFO / &
      VariableInfo( ATMOS_PHY_TB_AUX_T13_ID, 'TB_T13', 'stress tensor (T13)',                         & ! rho x nu x S : [kg/m3] x [m2/s] x [s-1]
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                                ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T23_ID, 'TB_T23', 'stress tensor (T23)',                         &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                                ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T33_ID, 'TB_T33', 'stress tensor (T33)',                         &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                                ), &
      VariableInfo( ATMOS_PHY_TB_AUX_DIFFFLX3_ID, 'DIFF_FLX3', 'diffusive heat flux (z) / density',   &
                    'm2/s.K/m',  3, 'XYZ',  ''                                                     ), &                                      
      VariableInfo( ATMOS_PHY_TB_AUX_T31_ID, 'TB_T31', 'stress tensor (T31)',                         &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                                ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T32_ID, 'TB_T32', 'stress tensor (T32)',                         &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                                ), &
      VariableInfo( ATMOS_PHY_TB_AUX_DIFFFLX1_ID, 'DIFF_FLX1', 'diffusive heat flux (x) / density',   &
                    'm2/s.K/m',  3, 'XYZ',  ''                                                     ), &
      VariableInfo( ATMOS_PHY_TB_AUX_DIFFFLX2_ID, 'DIFF_FLX2', 'diffusive heat flux (y) / density',   &
                    'm2/s.K/m',  3, 'XYZ',  ''                                                     ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T11_ID, 'TB_T11', 'stress tensor (T11)',                         &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                                ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T12_ID, 'TB_T12', 'stress tensor (T12)',                         &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                                ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T21_ID, 'TB_T21', 'stress tensor (T21)',                         &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                                ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T22_ID, 'TB_T22', 'stress tensor (T22)',                         &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                                )  /

    type(VariableInfo) :: ATMOS_PHY_TB_DIAG_VINFO(ATMOS_PHY_TB_DIAG_NUM)
    DATA ATMOS_PHY_TB_DIAG_VINFO / &
      VariableInfo( ATMOS_PHY_TB_DIAG_TKE_ID, 'TKE', 'SGS turbluence kinetic energy',   &
                    'm2/s2',  3, 'XYZ',  ''                                          ), &
      VariableInfo( ATMOS_PHY_TB_DIAG_NU_ID, 'NU', 'eddy viscosity',                    &
                    'm2/s',  3, 'XYZ',  ''                                           ), &
      VariableInfo( ATMOS_PHY_TB_DIAG_KH_ID, 'KH', 'eddy diffusion',                    &
                    'm2/s',  3, 'XYZ',  ''                                           )  /
                    
    type(VariableInfo) :: ATMOS_PHY_TB_TEND_VINFO(ATMOS_PHY_TB_TENDS_NUM1)
    DATA ATMOS_PHY_TB_TEND_VINFO / &
      VariableInfo( ATMOS_PHY_TB_MOMX_t_ID, 'TB_MOMX_t', 'tendency of x-momentum in TB process',    &
                    'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
      VariableInfo( ATMOS_PHY_TB_MOMY_t_ID, 'TB_MOMY_t', 'tendency of y-momentum in TB process',    &
                    'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
      VariableInfo( ATMOS_PHY_TB_MOMZ_t_ID, 'TB_MOMZ_t', 'tendency of z-momentum in TB process',    &
                    'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
      VariableInfo( ATMOS_PHY_TB_RHOT_t_ID, 'TB_RHOT_t', 'tendency of rho*PT in TB process',        &
                    'kg/m3.K/s', 3, 'XYZ',  ''                                                   )  / 

    !----------------------------------------------------------

    auxvar_info(:) = ATMOS_PHY_TB_AUX_VINFO
    diagvar_info(:) = ATMOS_PHY_TB_DIAG_VINFO
    tend_info(:) = ATMOS_PHY_TB_TEND_VINFO

    return
  end subroutine atm_phy_tb_dgm_common_get_varinfo

!OCL SERIAL  
  subroutine atm_phy_tb_dgm_common_calc_lambda( lambda, & ! (out)
    Cs, filter_fac, lmesh, elem, lmesh2D, elem2D        ) ! (in)

    use scale_const, only: &
      KARMAN  => CONST_KARMAN
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: lambda(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: Cs
    real(RP), intent(in) :: filter_fac

    integer :: ke

    real(RP) :: vol
    real(RP) :: he
    real(RP) :: lambda0
    real(RP) :: Zs(elem2D%Np)
    real(RP) :: Z1(elem2D%Np)
    real(RP) :: dz(elem  %Np)
    real(RP) :: FZ
    real(RP) :: elem_aspect_eff

    real(RP), parameter :: OneOverThree  = 1.0_RP / 3.0_RP
    !--------------------------------------------------------------------

    !$omp parallel do private( &
    !$omp vol, lambda0, Zs, Z1, dz, he, FZ, elem_aspect_eff )
    do ke=lmesh%NeS, lmesh%NeE
      vol = sum( elem%IntWeight_lgl(:) * lmesh%J(:,ke) * lmesh%Gsqrt(:,ke) )
      he = ( vol / ( dble(elem%PolyOrder_h+1)**2 * dble(elem%PolyOrder_v+1) ) )**OneOverThree

      FZ = lmesh%zlev(elem%Colmask(elem%Nnode_v,1),ke) - lmesh%zlev(elem%Colmask(1,1),ke)      
      elem_aspect_eff = fact( FZ/dble(elem%PolyOrder_v+1), sqrt(vol/FZ)/dble(elem%PolyOrder_h+1), sqrt(vol/FZ)/dble(elem%PolyOrder_h+1) )
      lambda0 = elem_aspect_eff * Cs * filter_fac * he

      Zs(:) = lmesh%zlev(elem%Hslice(:,1),lmesh%EMap3Dto2D(ke))
      Z1(:) = lmesh%zlev(elem%Hslice(:,2),lmesh%EMap3Dto2D(ke))
      dz(:) = max( lmesh%zlev(:,ke) - Zs(elem%IndexH2Dto3D(:)), he )

      !lambda(:,ke) = sqrt( 1.0_RP / (1.0_RP / lambda0**2 + 1.0_RP / ( KARMAN * max( dz(:), EPS ) )**2 ) )
      lambda(:,ke) = sqrt( 1.0_RP / (1.0_RP / lambda0**2 + 1.0_RP / ( KARMAN * ( dz(:) + 1.0E-4_RP ) )**2 ) )    
    end do

    return
  end subroutine atm_phy_tb_dgm_common_calc_lambda


!> Calculate parameterized diffusive mass flux of tracer with turbulent model
!!
!OCL SERIAL  
  subroutine atm_phy_tb_dgm_common_cal_grad_qtrc( &
    DFQ1, DFQ2, DFQ3,                                           & ! (out)
    dRdx, dRdy, dRdz,                                           & ! (inout)
    Kh, QTRC, DDENS, DENS_hyd,                                  & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound, cal_grad_dens                                     ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
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
    real(RP) :: bnd_flux(elem%NfpTot,lmesh%Ne,3)
    real(RP) :: bnd_flux_rho (elem%NfpTot,lmesh%Ne,3)    
    
    real(RP) :: DENS(elem%Np), RDENS(elem%Np)
    real(RP) :: RHOxQTRC(elem%Np)

    integer :: ke
    integer :: p
    !--------------------------------------------------------------------

    call cal_bnd_flux_grad_qtrc( bnd_flux, bnd_flux_rho,                      & ! (out)
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
        call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * bnd_flux_rho(:,ke,1), LiftDelFlx )
        dRdx(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:)
  
        call sparsemat_matmul( Dy, DENS, Fy )
        call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * bnd_flux_rho(:,ke,2), LiftDelFlx )
        dRdy(:,ke) = lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:)
  
        call sparsemat_matmul( Dz, DENS, Fz )
        call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * bnd_flux_rho(:,ke,3), LiftDelFlx )
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
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * bnd_flux(:,ke,1), LiftDelFlx )
      DFQ1(:,ke) = Kh(:,ke) * ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - QTRC(:,ke) * dRdx(:,ke) ) * RDENS(:)

      call sparsemat_matmul( Dy, RHOxQTRC(:), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * bnd_flux(:,ke,2), LiftDelFlx )
      DFQ2(:,ke) = Kh(:,ke) * ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - QTRC(:,ke) * dRdy(:,ke) ) * RDENS(:)

      call sparsemat_matmul( Dz, RHOxQTRC(:), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * bnd_flux(:,ke,3), LiftDelFlx )
      DFQ3(:,ke) = Kh(:,ke) * ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - QTRC(:,ke) * dRdz(:,ke) ) * RDENS(:)
    end do

    !$omp end parallel

    return
  end subroutine atm_phy_tb_dgm_common_cal_grad_qtrc

!> Calculate tendecies with turbulent model
!!
!OCL SERIAL  
  subroutine atm_phy_tb_dgm_common_cal_tend( &
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
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
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
    real(RP) :: bnd_flux_mom(elem%NfpTot,lmesh%Ne,3)
    real(RP) :: bnd_flux_rhot(elem%NfpTot,lmesh%Ne)    
    !--------------------------------------------------------------------

    call cal_bnd_flux( bnd_flux_mom, bnd_flux_rhot,                           & ! (out)
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
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * bnd_flux_mom(:,ke,1), LiftDelFlx )

      MOMX_t(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + lmesh%Escale(:,ke,2,2) * Fy(:) &
                   + lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      ! MOMY
      call sparsemat_matmul( Dx, T21(:,ke), Fx )
      call sparsemat_matmul( Dy, T22(:,ke), Fy )
      call sparsemat_matmul( Dz, T23(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * bnd_flux_mom(:,ke,2), LiftDelFlx )

      MOMY_t(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + lmesh%Escale(:,ke,2,2) * Fy(:) &
                   + lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      ! MOMZ
      call sparsemat_matmul( Dx, T31(:,ke), Fx )
      call sparsemat_matmul( Dy, T32(:,ke), Fy )
      call sparsemat_matmul( Dz, T33(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * bnd_flux_mom(:,ke,3), LiftDelFlx )

      MOMZ_t(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + lmesh%Escale(:,ke,2,2) * Fy(:) &
                   + lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      ! RHOT
      call sparsemat_matmul( Dx, DENS(:) * DF1(:,ke), Fx )
      call sparsemat_matmul( Dy, DENS(:) * DF2(:,ke), Fy )
      call sparsemat_matmul( Dz, DENS(:) * DF3(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * bnd_flux_rhot(:,ke), LiftDelFlx )

      RHOT_t(:,ke)  = ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                      + lmesh%Escale(:,ke,2,2) * Fy(:) &
                      + lmesh%Escale(:,ke,3,3) * Fz(:) &
                      + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)

    end do

    return
  end subroutine atm_phy_tb_dgm_common_cal_tend

!> Calculate tendecies of tracer density with turbulent model
!!
!OCL SERIAL  
  subroutine atm_phy_tb_dgm_common_cal_tend_qtrc( &
    RHOQ_t,                                                     & ! (out)
    DFQ1, DFQ2, DFQ3,                                           & ! (in)
    Kh, DDENS_,DENS_hyd,                                        & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound                                                    ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
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

    call cal_bnd_flux_qtrc( del_flux,                                         & ! (out)
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
  end subroutine atm_phy_tb_dgm_common_cal_tend_qtrc

!-- private ---------------------------------------------


!OCL SERIAL  
  subroutine cal_bnd_flux_grad_qtrc( bnd_flux, bnd_flux_rho, & ! (out)
    QTRC_, DDENS_, DENS_hyd_,                                & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound,         & ! (in)
    cal_grad_dens                                            ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    real(RP), intent(out) ::  bnd_flux(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(out) ::  bnd_flux_rho(elem%NfpTot*lmesh%Ne,3)
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
        bnd_flux_rho(i,1) = facx * del * nx(i)
        bnd_flux_rho(i,2) = facy * del * ny(i)
        bnd_flux_rho(i,3) = facz * del * nz(i)
      end if
      
      del = 0.5_RP * ( densP * QTRC_(iP) - densM * QTRC_(iM) )
      bnd_flux(i,1) = facx * del * nx(i)
      bnd_flux(i,2) = facy * del * ny(i)
      bnd_flux(i,3) = facz * del * nz(i)
    end do

    return
  end subroutine cal_bnd_flux_grad_qtrc

!OCL SERIAL  
  subroutine cal_bnd_flux( bnd_flux_mom, bnd_flux_rhot,                & ! (out)
    T11, T12, T13, T21, T22, T23, T31, T32, T33,                       & ! (in)
    DF1, DF2, DF3,                                                     & ! (in)
    Nu, Kh,                                                            & ! (in)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,           & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound                    ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    real(RP), intent(out) ::  bnd_flux_mom (elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(out) ::  bnd_flux_rhot(elem%NfpTot*lmesh%Ne) 
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

      bnd_flux_mom(i,1) = 0.5_RP * ( TauP_x - TauM_x )
      bnd_flux_mom(i,2) = 0.5_RP * ( TauP_y - TauM_y )
      bnd_flux_mom(i,3) = 0.5_RP * ( TauP_z - TauM_z )
      bnd_flux_rhot(i)  = 0.5_RP * ( densP * ( DF1(iP) * nx_ + DF2(iP) * ny_ + DF3(iP) * nz_ ) &
                                   - densM * ( DF1(iM) * nx_ + DF2(iM) * ny_ + DF3(iM) * nz_ ) )
    end do

    return
  end subroutine cal_bnd_flux

!OCL SERIAL  
  subroutine cal_bnd_flux_qtrc( bnd_flux,             & ! (out)
    DFQ1, DFQ2, DFQ3,                                 & ! (in)
    Kh, DDENS_, DENS_hyd,                             & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound   ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    real(RP), intent(out) ::  bnd_flux(elem%NfpTot*lmesh%Ne) 
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

      bnd_flux(i)  = 0.5_RP * ( densP * ( DFQ1(iP) * nx_ + DFQ2(iP) * ny_ + DFQ3(iP) * nz_ ) &
                              - densM * ( DFQ1(iM) * nx_ + DFQ2(iM) * ny_ + DFQ3(iM) * nz_ ) )
    end do

    return
  end subroutine cal_bnd_flux_qtrc

!OCL SERIAL
  elemental function fact(dz, dx, dy)
    implicit none
    real(RP), intent(in) :: dz
    real(RP), intent(in) :: dx
    real(RP), intent(in) :: dy
    real(RP) :: fact ! (out)

    real(RP), parameter :: oot = -1.0_RP/3.0_RP
    real(RP), parameter :: fot =  5.0_RP/3.0_RP
    real(RP), parameter :: eot = 11.0_RP/3.0_RP
    real(RP), parameter :: tof = -3.0_RP/4.0_RP
    real(RP) :: a1, a2, b1, b2, dmax
    !--------------------------------------------------------------------

    dmax = max(dz, dx, dy)
    if ( dz == dmax ) then
       a1 = dx / dmax
       a2 = dy / dmax
    else if ( dx == dmax ) then
       a1 = dz / dmax
       a2 = dy / dmax
    else ! dy == dmax
       a1 = dz / dmax
       a2 = dx / dmax
    end if
    b1 = atan( a1/a2 )
    b2 = atan( a2/a1 )

   fact = 1.736_RP * (a1*a2)**oot &
         * ( 4.0_RP*p1(b1)*a1**oot + 0.222_RP*p2(b1)*a1**fot + 0.077*p3(b1)*a1**eot - 3.0_RP*b1 &
           + 4.0_RP*p1(b2)*a2**oot + 0.222_RP*p2(b2)*a2**fot + 0.077*p3(b2)*a2**eot - 3.0_RP*b2 &
           )**tof
   return
  end function fact
!OCL SERIAL
  elemental function p1(z)
    implicit none
    real(RP), intent(in) :: z
    real(RP) :: p1 ! (out)

    real(RP), parameter :: TwoOverThree  = 2.0_RP / 3.0_RP
    !--------------------------------------------------------------------

    p1 = 2.5_RP * p2(z) - 1.5_RP * sin(z) * cos(z)**TwoOverThree
    return
  end function p1
!OCL SERIAL
  elemental function p2(z)
    implicit none
    real(RP), intent(in) :: z
    real(RP) :: p2 ! (out)
    !--------------------------------------------------------------------

    p2 = 0.986_RP * z + 0.073_RP * z**2 - 0.418_RP * z**3 + 0.120_RP * z**4
    return
  end function p2
!OCL SERIAL
  elemental function p3(z)
    implicit none
    real(RP), intent(in) :: z
    real(RP) :: p3 ! (out)
    !--------------------------------------------------------------------

    p3 = 0.976_RP * z + 0.188_RP * z**2 - 1.169_RP * z**3 + 0.755_RP * z**4 - 0.151_RP * z**5
    return
  end function p3

end module scale_atm_phy_tb_dgm_common
