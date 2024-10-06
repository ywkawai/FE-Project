!> module FElib / Atmosphere / Physics turbulence
!!
!! @par Description
!!      Turbulnce process for DNS
!!
!! @author Team SCALE
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_tb_dgm_dns
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
  public :: atm_phy_tb_dgm_dns_Init
  public :: atm_phy_tb_dgm_dns_Final
  public :: atm_phy_tb_dgm_dns_cal_grad

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: cal_del_flux_grad

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter   :: OneOverThree  = 1.0_RP / 3.0_RP
  real(RP), private, parameter   :: twoOverThree  = 2.0_RP / 3.0_RP
  real(RP), private, parameter   :: FourOverThree = 4.0_RP / 3.0_RP

  real(RP), private              :: DNS_NU        = 1.512E-5_RP !< kinematic viscosity coefficient [m2/s] for air at 20degC
  real(RP), private              :: DNS_MU        = 1.8E-5_RP   !< molecular diffusive coefficient [m2/s] for air at 20degC

contains
!OCL SERIAL
  subroutine atm_phy_tb_dgm_dns_Init( mesh )
    implicit none    
    class(MeshBase3D), intent(in) :: mesh

    logical  :: consistent_tke = .true.

    namelist / PARAM_ATMOS_PHY_TB_DGM_DNS / &
      DNS_NU, DNS_MU
        
    integer :: ierr
    !--------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_dgm_dns_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_TB_dgm_dns_setup",*) 'Eddy Viscocity Model for DNS'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_DGM_DNS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_TB_dgm_dns_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_TB_dgm_dns_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_TB_DGM_DNS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_TB_DGM_DNS)

    return
  end subroutine atm_phy_tb_dgm_dns_Init

!OCL SERIAL
  subroutine atm_phy_tb_dgm_dns_Final()
    implicit none
    !--------------------------------------------------------------------

    return
  end subroutine atm_phy_tb_dgm_dns_Final

!> Calculate parameterized stress tensor and eddy heat flux with turbulent model
!!
!OCL SERIAL  
  subroutine atm_phy_tb_dgm_dns_cal_grad( &
    T11, T12, T13, T21, T22, T23, T31, T32, T33,                & ! (out)
    DF1, DF2, DF3,                                              & ! (out)
    TKE, Nu, Kh,                                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,    & ! (in)
    PRES, PT,                                                   & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
    is_bound                                                    ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
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
    real(RP) :: coef

    integer :: ke
    integer :: p
    !--------------------------------------------------------------------

    call cal_del_flux_grad( del_flux_rho, del_flux_mom, del_flux_rhot,        & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, PT,            & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound )                                                   ! (in)
  
    !$omp parallel do private( &
    !$omp Fx, Fy, Fz, LiftDelFlx,                  &
    !$omp DENS, RHOT, RDENS, Q, DdensDxi, DVelDxi, &
    !$omp S11, S12, S22, S23, S31, S33,            &
    !$omp coef, SkkOvThree,                        &
    !$omp p )
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

      Nu(:,ke) = DNS_NU
      Kh(:,ke) = DNS_MU

      !---

      do p=1, elem%Np
        coef = 2.0_RP * Nu(p,ke)

        T11(p,ke) = DENS(p) * ( coef * ( S11(p) - SkkOvThree )  )
        T12(p,ke) = DENS(p) * coef * S12(p)
        T13(p,ke) = DENS(p) * coef * S31(p)

        T21(p,ke) = DENS(p) * coef * S12(p)
        T22(p,ke) = DENS(p) * ( coef * ( S22(p) - SkkOvThree )  )
        T23(p,ke) = DENS(p) * coef * S23(p)

        T31(p,ke) = DENS(p) * coef * S31(p)
        T32(p,ke) = DENS(p) * coef * S23(p)
        T33(p,ke) = DENS(p) * ( coef * ( S33(p) - SkkOvThree )  )
      end do

      DF1(:,ke) = Kh(:,ke) * DF1(:,ke)
      DF2(:,ke) = Kh(:,ke) * DF2(:,ke)
      DF3(:,ke) = Kh(:,ke) * DF3(:,ke)
    end do
  
    return
  end subroutine atm_phy_tb_dgm_dns_cal_grad

!OCL SERIAL  
  subroutine cal_del_flux_grad( del_flux_rho, del_flux_mom, del_flux_rhot,  & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, PT_,           & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound                         ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
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

!-- private --------------------------------------------------------

end module scale_atm_phy_tb_dgm_dns
