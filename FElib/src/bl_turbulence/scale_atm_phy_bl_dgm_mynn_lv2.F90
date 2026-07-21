!> module FElib / Atmosphere / Physics / boundary layer turbulence
!!
!! @par Description
!! A module to provide routines to calculate turbulent diffusion coefficients based on the Mellor-Yamada-Nakanishi-Niino (MYNN) Level 2 closure.
!!
!! This implementation uses the Level 2 diagnostic closure only. Prognostic turbulent kinetic energy, nonlocal mixing, and the MY Level 2.5 equations are not included.
!!
!! @par Reference
!! Mellor, G. L. and T. Yamada, 1982:
!! Development of a turbulence closure model for geophysical fluid problems.
!! Reviews of Geophysics and Space Physics, 20, 851-875.
!!
!! Nakanishi, M. and H. Niino, 2009:
!! Development of an improved turbulence closure model for the atmospheric boundary layer.
!! Journal of the Meteorological Society of Japan, 87, 895-912.
!!
!! Blackadar, A. K., 1962:
!! The vertical distribution of wind and turbulent exchange in a neutral atmosphere.
!! Journal of Geophysical Research, 67, 3095-3102.
!!
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_bl_dgm_mynn_lv2
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
  public :: atm_phy_bl_dgm_mynn_lv2_Init
  public :: atm_phy_bl_dgm_mynn_lv2_Final
  public :: atm_phy_bl_dgm_mynn_lv2_cal_VViscDiffCoef

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  !- MYNN closure constants -----------

  real(RP), parameter :: A1 = 1.18_RP
  real(RP), parameter :: A2 = 0.665_RP

  real(RP), parameter :: B1 = 24.0_RP
  real(RP), parameter :: B2 = 15.0_RP

  real(RP), parameter :: C1 = 0.137_RP
  real(RP), parameter :: C2 = 0.75_RP
  real(RP), parameter :: C3 = 0.352_RP
  real(RP), parameter :: C5 = 0.2_RP

  real(RP), parameter :: G1 = 0.235_RP

  ! Derived constants
  real(RP) :: G2
  real(RP) :: F2
  real(RP) :: Rf2
  real(RP) :: RFc

  ! Limiter with mixing length
  real(RP) :: L_INF = 100.0_RP
  real(RP) :: L_MIN = 1.0E-6_RP

contains
  !> Initialize a module of MYNN Level 2 PBL turbulence parameterization
!OCL SERIAL
  subroutine atm_phy_bl_dgm_mynn_lv2_Init( mesh )
    implicit none    
    class(MeshBase3D), intent(in) :: mesh

    namelist / PARAM_ATMOS_PHY_BL_DGM_MYNN_LV2 / &
      L_INF, L_MIN
    
    integer :: ierr    
    !--------------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_DGM_MYNN_LV2_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_BL_DGM_MYNN_LV2_setup",*) 'MYNN Level 2 PBL turbulence parameterization'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_BL_DGM_MYNN_LV2,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_BL_DGM_MYNN_LV2_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_BL_DGM_MYNN_LV2_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_BL_DGM_MYNN_LV2. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_BL_DGM_MYNN_LV2)

    !-

    G2 = ( 2.0_RP * A1 * ( 3.0_RP - 2.0_RP * C2 )  &
         + B2 * ( 1.0_RP - C3 )                    &
         ) / B1

    F2 = B1 * ( G1 + G2 ) &
        - 3.0_RP * A1 * ( 1.0_RP - C2 )

    Rf2 = B1 * G1 / F2

    RFc = G1 / ( G1 + G2 )
    return
  end subroutine atm_phy_bl_dgm_mynn_lv2_Init

  !> Finalize a module of MYNN Level 2 PBL turbulence parameterization
!OCL SERIAL
  subroutine atm_phy_bl_dgm_mynn_lv2_Final()
    implicit none    
    !--------------------------------------------------------------------------------
    return
  end subroutine atm_phy_bl_dgm_mynn_lv2_Final

!OCL SERIAL
  subroutine atm_phy_bl_dgm_mynn_lv2_cal_VViscDiffCoef( &
    Nu, Kh, TKE,                                             & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, & ! (in)
    PRES, PT,                                                & ! (in)
    Dz, Lift, lmesh, elem, is_bound  )                         ! (in)
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: Nu(elem%Np,lmesh%NeA)       !< Vertical eddy viscosity
    real(RP), intent(out) :: Kh(elem%Np,lmesh%NeA)       !< Vertical eddy diffusivity
    real(RP), intent(out) :: TKE(elem%Np,lmesh%NeA)      !< Turbulent kinetic energy
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)   !< Density perturbation
    real(RP), intent(in)  :: MOMX_ (elem%Np,lmesh%NeA)   !< Momentum in x1 direction
    real(RP), intent(in)  :: MOMY_ (elem%Np,lmesh%NeA)   !< Momentum in x2 direction
    real(RP), intent(in)  :: MOMZ_ (elem%Np,lmesh%NeA)   !< Momentum in x3 direction
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)   !< Density x potential temperature perturbation
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA) !< Reference density in hydrostatic balance
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA) !< Reference pressure in hydrostatic balance
    real(RP), intent(in)  :: PRES(elem%Np,lmesh%NeA)     !< Pressure
    real(RP), intent(in)  :: PT(elem%Np,lmesh%NeA)       !< Potential temperature
    type(SparseMat), intent(in) :: Dz, Lift
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne) !< Flag whether nodes are located at domain boundaries

    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DENS(elem%Np), RDENS(elem%Np), RHOT(elem%Np), Q(elem%Np)
    real(RP) :: DdensDz(elem%Np), DVelDz(elem%Np,3), DptDz(elem%Np)

    real(RP) :: del_flux_rho (elem%NfpTot,lmesh%Ne)
    real(RP) :: del_flux_mom (elem%NfpTot,lmesh%Ne,3)
    real(RP) :: del_flux_rhot(elem%NfpTot,lmesh%Ne)

    real(RP) :: S2(elem%Np) ! Square of vertical shear of horizontal wind
    real(RP) :: Ri          ! Richardson number
    real(RP) :: Rf(elem%Np) ! Flux Richardson number

    real(RP) :: A2_loc
    real(RP) :: F1
    real(RP) :: Rf1
    real(RP) :: AF12
    real(RP) :: discriminant
    real(RP) :: denom_h
    real(RP) :: denom_m
    real(RP) :: S_M(elem%Np), S_H(elem%Np)     ! Stability functions for momentum and heat

    real(RP) :: mixlen(elem%Np)                    ! Mixing length
    real(RP) :: zsfc(elem%Nnode_h1D**2,lmesh%Ne2D) ! Surface height
    real(RP) :: kz

    integer :: ke, ke2D
    integer :: p, ph, pz
    class(LocalMesh2D), pointer :: lmesh2D
    class(ElementBase2D), pointer :: elem2D

    real(RP), parameter :: EPS_RF   = 1.0e-10_RP
    real(RP), parameter :: EPS_DISC = 1.0e-14_RP    
    !--------------------------------------------------------------------------------

    ! For standard case of NN2009, these parameter are independent of Ri.
    ! If we will implement the K2010 correction in the future, 
    ! these parameter should be calculated in the ke loop because they depend on Ri (i.e., A2_loc = A_2/(1+Ri)).
    A2_loc = A2

    F1 = B1 * ( G1 - C1 )                                    &
       + 2.0_RP * A1 * ( 3.0_RP - 2.0_RP * C2 )              &
       + 3.0_RP * A2_loc * ( 1.0_RP - C3 ) * ( 1.0_RP - C5 )
    
    Rf1 = B1 * ( G1  - C1 ) / F1

    AF12 = A1 * F1 / ( A2_loc * F2 )

    lmesh2D => lmesh%lcmesh2D
    elem2D  => lmesh2D%refElem2D

    !---------------------------------

    call cal_del_flux_grad( del_flux_rho, del_flux_mom, del_flux_rhot,        & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, PT,            & ! (in)
      lmesh%normal_fn(:,:,3), lmesh%vmapM, lmesh%vmapP,                       & ! (in)
      lmesh, elem, is_bound )                                                   ! (in)

    !$omp parallel &
    !$omp private( Fz, LiftDelFlx, DENS, RDENS, RHOT, Q, DdensDz, DVelDz, DptDz, &
    !$omp S2, Ri, Rf, discriminant, denom_m, denom_h, S_M, S_H, mixlen, kz )

    !$omp do
    do ke2D=lmesh2D%NeS, lmesh2D%NeE
      zsfc(:,ke2D) = lmesh%zlev(elem%Hslice(:,1),ke2D)
    end do
    !$omp end do
    !$omp do
    do ke=lmesh%NeS, lmesh%NeE
      DENS (:) = DENS_hyd(:,ke) + DDENS_(:,ke)
      RDENS(:) = 1.0_RP / DENS(:)
      RHOT(:) = DENS(:) * PT(:,ke)

      ! gradient of density
      call sparsemat_matmul( Dz, DENS, Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rho(:,ke), LiftDelFlx )
      DdensDz(:) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      ! gradient of u
      Q(:) = MOMX_(:,ke) * RDENS(:) 
      call sparsemat_matmul( Dz, MOMX_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,1), LiftDelFlx )
      DVelDz(:,1) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDz(:) ) * RDENS(:)

      ! gradient of v
      Q(:) = MOMY_(:,ke) * RDENS(:) 
      call sparsemat_matmul( Dz, MOMY_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2), LiftDelFlx )
      DVelDz(:,2) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDz(:) ) * RDENS(:)
    
      ! gradient of pt
      Q(:) = RHOT(:) * RDENS(:)
      call sparsemat_matmul( Dz, RHOT, Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke), LiftDelFlx )
      DptDz(:) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDz(:) ) * RDENS(:)

     ! Calculate flux Richardson number: Rf
      do p=1, elem%Np
        S2(p) = DVelDz(p,1)**2 + DVelDz(p,2)**2
        Ri = GRAV * DptDz(p) / ( PT(p,ke) * max(S2(p), EPS) )

        discriminant = Ri * Ri &
                     + 2.0_RP * AF12 * ( Rf1 - 2.0_RP * Rf2 ) * Ri &
                     + ( AF12 * Rf1 )**2
        discriminant = max( discriminant, 0.0_RP )

        Rf(p) = 0.5_RP / AF12 * ( Ri + AF12 * Rf1 - sqrt(discriminant) )
        Rf(p) = min( Rf(p), RFc - EPS_RF )
      end do

      ! Calculate stability function: S_M and S_H
      do p=1, elem%Np
        denom_h = max( 1.0_RP - Rf(p), EPS_RF )
        S_H(p) = 3.0_RP * A2_loc * ( G1 + G2 ) * ( RFc - Rf(p) ) / denom_h
        
        denom_m = Rf2 - Rf(p)
        if ( abs(denom_m) < EPS_RF ) then
            denom_m = sign(EPS_RF, denom_m)
        end if
        S_M(p) = S_H(p) * AF12 * ( Rf1 - Rf(p) ) / denom_m

        S_M(p) = max( S_M(p), 0.0_RP )
        S_H(p) = max( S_H(p), 0.0_RP )
      end do

      ! Calculate mixing length

      ke2D = lmesh%EMap3Dto2D(ke)
      do pz=1, elem%Nnode_v
      do ph=1, elem%Nnode_h1D**2
        p = ph + (pz-1)*elem%Nnode_h1D**2
        kz = KARMAN * ( lmesh%zlev(p,ke) - zsfc(ph,ke2D) )

        mixlen(p) = kz * L_INF  / max( kz + L_INF, L_MIN )
      end do
      end do

      ! Calculate eddy viscosity and diffusivity
      do p=1, elem%Np
        ! q^2
        q(p) = B1 * mixlen(p)**2 &
             * S_M(p) * max( 1.0_RP - Rf(p), 0.0_RP ) * S2(p)
        TKE(p,ke) = 0.5_RP * max( q(p), 0.0_RP )

        q(p) = sqrt( 2.0_RP * TKE(p,ke) )
        Kh(p,ke) = mixlen(p) * q(p) * S_H(p)
        Nu(p,ke) = mixlen(p) * q(p) * S_M(p)
      end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine atm_phy_bl_dgm_mynn_lv2_cal_VViscDiffCoef

!-- private --------------------------------------------------------

!OCL SERIAL  
  subroutine cal_del_flux_grad( del_flux_rho, del_flux_mom, del_flux_rhot,  & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, PT_,           & ! (in)
    nz, vmapM, vmapP, lmesh, elem, is_bound                                 ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux_rho(elem%NfpTot*lmesh%Ne)
    real(RP), intent(out) ::  del_flux_mom(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(out) ::  del_flux_rhot(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: densM, densP
    real(RP) :: del
    real(RP) :: facz
    !------------------------------------------------------------------------
    
    !$omp parallel do private ( iM, iP, densM, densP, del, facz )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iP)

      if ( is_bound(i) ) then
        facz = 1.0_RP
      else
        ! facz = 1.0_RP - sign(1.0_RP,nz(i))
        facz = 1.0_RP
      end if

      del = 0.5_RP * ( densP - densM )
      del_flux_rho(i) = facz * del * nz(i)

      del = 0.5_RP * ( MOMX_(iP) - MOMX_(iM) )
      del_flux_mom(i,1) = facz * del * nz(i)

      del = 0.5_RP * ( MOMY_(iP) - MOMY_(iM) )
      del_flux_mom(i,2) = facz * del * nz(i)

      del = 0.5_RP * ( MOMZ_(iP) - MOMZ_(iM) )
      del_flux_mom(i,3) = facz * del * nz(i)

      del = 0.5_RP * ( densP * PT_(iP) - densM * PT_(iM) )
      del_flux_rhot(i) = facz * del * nz(i)
    end do

    return
  end subroutine cal_del_flux_grad
end module scale_atm_phy_bl_dgm_mynn_lv2