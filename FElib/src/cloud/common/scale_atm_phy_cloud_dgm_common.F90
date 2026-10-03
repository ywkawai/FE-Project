!> module FElib / Atmosphere / Physics cloud  / common
!!
!! @par Description
!!      Sedimentation of hydrometeors with cloud process
!!      common subroutines
!!  
!! To preserve nonnegativity in precipitation process, 
!! a limiter proposed by Light and Durran (2016, MWR) is used
!!
!! @par Reference
!!  - Light and Durran 2016: 
!!    Preserving Nonnegativity in Discontinuous Galerkin Approximations to Scalar Transport via Truncation and Mass Aware Rescaling (TMAR).
!!    Monthly Weather Review, 144(12), 4771–4786.
!!
!! @author Yuta Kawai, Team SCALE
!!
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_cloud_dgm_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    UNDEF => CONST_UNDEF8, &
    GRAV => CONST_GRAV,    &
    PRES00 => CONST_PRE00
  use scale_atmos_hydrometeor, only: &
    CV_VAPOR, &
    CP_VAPOR, &
    CV_WATER, &
    CP_WATER, &
    CV_ICE,   &
    CP_ICE,   &
    LHV

  use scale_sparsemat
  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
  use scale_localmesh_3d, only: &
    LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  public :: atm_phy_cloud_dgm_common_sedimentation
  public :: atm_phy_cloud_dgm_common_sedimentation_momentum
  public :: atm_phy_cloud_dgm_common_condensate_removal
  public :: atm_phy_cloud_dgm_common_condensate_removal_momentum
  public :: atm_phy_cloud_dgm_common_negative_fixer

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  private :: atm_phy_cloud_dgm_netOutwardFlux
  private :: atm_phy_cloud_dgm_sedimentation_get_delflux
  private :: atm_phy_cloud_dgm_sedimentation_momentum_get_delflux

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------


contains
  !> Update the variable state through the sedimentation process
!OCL SERIAL
  subroutine atm_phy_cloud_dgm_common_sedimentation( &
    DENS, RHOQ, CPtot, CVtot, RHOE,         & ! (inout)
    FLX_hydro, sflx_rain, sflx_snow, esflx, & ! (inout)
    TEMP, vterm, dt, rnstep,                & ! (in)
    Dz, Lift, nz, vmapM, vmapP, IntWeight,  & ! (in)
    QHA, QLA, QIA, lcmesh, elem             ) ! (in)

    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: QHA                   !< hydrometeor (water + ice)
    real(RP), intent(inout) :: DENS (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: RHOQ (elem%Np,lcmesh%NeZ,lcmesh%Ne2D,QHA)
    real(RP), intent(inout) :: CPtot(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: CVtot(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: RHOE (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: FLX_hydro(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: sflx_rain(elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(inout) :: sflx_snow(elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(inout) :: esflx    (elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(in) :: TEMP (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: vterm(elem%Np,lcmesh%NeZ,lcmesh%Ne2D,QHA)
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: rnstep
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    real(RP), intent(in) :: nz(elem%NfpTot,lcmesh%NeZ,lcmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lcmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lcmesh%NeZ)
    real(RP), intent(in) :: IntWeight(elem%Nfaces,elem%NfpTot)
    integer, intent(in) :: QLA, QIA

    real(RP) :: qflx(elem%Np)
    real(RP) :: eflx(elem%Np)
    real(RP) :: DENS0(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOCP(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOCV(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: NDcoefEuler(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: DzRHOQ(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: DzRHOE(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: dDENS(elem%Np)
    real(RP) :: CP(QHA), CV(QHA)

    real(RP) :: fct_coef(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOQ0, RHOQ1, RHOQ_tmp(elem%Np)
    real(RP) :: netOutwardFlux(lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: del_flux(elem%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,2)

    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: RHOQ_save(elem%Np)

    integer :: ke2D
    integer :: ke_z
    integer :: ke
    integer :: iq

    real(RP) :: Q
    real(RP) :: delz
    !-------------------------------------------------------

    call set_hydrometeor_heat_capacity( QLA, QIA, QHA, & ! (in)
      CP, CV) ! (out)

    !$omp parallel do collapse(2)
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ
      DENS0(:,ke_z,ke2D) = DENS(:,ke_z,ke2D)
      RHOCP(:,ke_z,ke2D) = CPtot(:,ke_z,ke2D) * DENS(:,ke_z,ke2D)
      RHOCV(:,ke_z,ke2D) = CVtot(:,ke_z,ke2D) * DENS(:,ke_z,ke2D)
    end do
    end do

    ! Process each hydrometeor category separately
    do iq = 1, QHA
    
      call atm_phy_cloud_dgm_sedimentation_get_delflux_dq( &
        del_flux(:,:,:,:),                                                      & ! (out)
        DENS0(:,:,:), RHOQ(:,:,:,iq), TEMP(:,:,:), CV(iq), nz(:,:,:), vmapM(:,:), vmapP(:,:), & ! (in)
        lcmesh, elem                                                            ) ! (in)
      
      !$omp parallel do private( &
      !$omp ke2D, ke_z, ke, delz, Fz, LiftDelFlx )
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D
        delz = ( lcmesh%pos_ev(lcmesh%EToV(ke,5),3) - lcmesh%pos_ev(lcmesh%EToV(ke,1),3) ) / dble( elem%Nnode_v )

        call sparsemat_matmul( Dz, RHOQ(:,ke_z,ke2D,iq), Fz )
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,1), LiftDelFlx )
        DzRHOQ(:,ke_z,ke2D) = lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

        call sparsemat_matmul( Dz, RHOQ(:,ke_z,ke2D,iq) * CV(iq) * TEMP(:,ke_z,ke2D), Fz )
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,2), LiftDelFlx )
        DzRHOE(:,ke_z,ke2D) = lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)
        
        NDcoefEuler(:,ke_z,ke2D) = 0.5_RP * delz * abs(vterm(:,ke_z,ke2D,iq))
      end do
      end do

      call atm_phy_cloud_dgm_netOutwardFlux( &
        netOutwardFlux(:,:),                                                  & ! (out)
        RHOQ(:,:,:,iq), vterm(:,:,:,iq), DzRHOQ(:,:,:), NDcoefEuler(:,:,:),   & ! (in) 
        lcmesh%J(:,:), lcmesh%Fscale(:,:),                                    & ! (in)
        nz(:,:,:), vmapM(:,:), vmapP(:,:), lcmesh%VMapM(:,:), IntWeight(:,:), & ! (in)
        lcmesh, elem                                                          ) ! (in)
      
      !$omp parallel do collapse(2) private(ke, Q)
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D

        Q = sum( lcmesh%J(:,ke) * elem%IntWeight_lgl(:) * RHOQ(:,ke_z,ke2D,iq) ) / dt      
        fct_coef(:,ke_z,ke2D) = max( 0.0_RP, min( 1.0_RP, Q / ( netOutwardFlux(ke_z,ke2D) + 1.0E-10_RP ) ) )  
      end do ! end loop for ke_z
      end do ! end loop for ke2D

      call atm_phy_cloud_dgm_sedimentation_get_delflux( &
        del_flux(:,:,:,:),                                           & ! (out)
        DENS0(:,:,:), RHOQ(:,:,:,iq), TEMP(:,:,:), vterm(:,:,:,iq),  & ! (in)
        DzRHOQ(:,:,:), DzRHOE(:,:,:), NDcoefEuler(:,:,:),            & ! (in)
        fct_coef(:,:,:),                                             & ! (in)
        CV(iq), lcmesh%J(:,:), lcmesh%Fscale(:,:), nz(:,:,:),        & ! (in)
        vmapM(:,:), vmapP(:,:), lcmesh%vmapM(:,:), IntWeight(:,:),   & ! (in)
        lcmesh, elem                                                 ) ! (in)

      !$omp parallel do collapse(2) private( &
      !$omp ke2D, ke_z, ke,                            &
      !$omp qflx, eflx, dDENS, RHOQ_tmp, RHOQ0, RHOQ1, &
      !$omp RHOQ_save, &
      !$omp Fz, LiftDelFlx )
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D

        !--- update falling tracer 

        RHOQ_save(:) = RHOQ(:,ke_z,ke2D,iq)
        qflx(:) = vterm(:,ke_z,ke2D,iq) * RHOQ(:,ke_z,ke2D,iq)   &
                - NDcoefEuler(:,ke_z,ke2D) * DzRHOQ(:,ke_z,ke2D)

        call sparsemat_matmul( Dz, qflx(:), Fz )
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,1), LiftDelFlx )
        
        dDENS(:) =  - dt * ( &
          lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )
        RHOQ_tmp(:) = max( 0.0_RP, RHOQ(:,ke_z,ke2D,iq) + dDENS(:) )
!        RHOQ_tmp(:) = RHOQ(:,ke_z,ke2D,iq) + dDENS(:)

        !
        RHOQ0 = sum( lcmesh%Gsqrt(:,ke) * lcmesh%J(:,ke) * elem%IntWeight_lgl(:) * ( RHOQ(:,ke_z,ke2D,iq) + dDENS(:) ) )
        RHOQ1 = sum( lcmesh%Gsqrt(:,ke) * lcmesh%J(:,ke) * elem%IntWeight_lgl(:) * RHOQ_tmp(:)          )

        dDENS(:) = RHOQ0 / ( RHOQ1 + 1.0E-32_RP ) * RHOQ_tmp(:) &
                 - RHOQ(:,ke_z,ke2D,iq)
        RHOQ(:,ke_z,ke2D,iq) = RHOQ(:,ke_z,ke2D,iq) + dDENS(:)


        ! QTRC(iq; iq>QLA+QLI) is not mass tracer, such as number density
        if ( iq > QLA + QIA ) cycle

        FLX_hydro(:,ke_z,ke2D) = FLX_hydro(:,ke_z,ke2D) &
                               + qflx(:) * rnstep
        if ( ke_z == 1 ) then
          if ( iq > QLA ) then ! ice water
              sflx_snow(:,ke2D) = sflx_snow(:,ke2D)               &
                                + qflx(elem%Hslice(:,1)) * rnstep
          else                 ! liquid water
              sflx_rain(:,ke2D) = sflx_rain(:,ke2D)               &
                                + qflx(elem%Hslice(:,1)) * rnstep
          end if
        end if

        !--- update density

        RHOCP(:,ke_z,ke2D) = RHOCP(:,ke_z,ke2D) + CP(iq) * dDENS(:)
        RHOCV(:,ke_z,ke2D) = RHOCV(:,ke_z,ke2D) + CV(iq) * dDENS(:)
        DENS (:,ke_z,ke2D) = DENS(:,ke_z,ke2D) + dDENS(:)
  
        !--- update internal energy   

        eflx(:) = vterm(:,ke_z,ke2D,iq) * RHOQ_save(:) * TEMP(:,ke_z,ke2D) * CV(iq)  &
                - NDcoefEuler(:,ke_z,ke2D) * DzRHOE(:,ke_z,ke2D)
        
        call sparsemat_matmul( Dz, eflx(:), Fz )
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,2), LiftDelFlx )

        RHOE(:,ke_z,ke2D) = RHOE(:,ke_z,ke2D) - dt * ( &
          + lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) &
          + qflx(:) * Grav                                  )

        if ( ke_z == 1 ) then
          esflx(:,ke2D) = esflx(:,ke2D) &
                        + eflx(elem%Hslice(:,1)) * rnstep
        end if

      end do ! end loop for ke_z
      end do ! end loop for ke2D
    end do ! end loop for iq

    !$omp parallel do collapse(2)
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ  
      CPtot(:,ke_z,ke2D) = RHOCP(:,ke_z,ke2D) / DENS(:,ke_z,ke2D)
      CVtot(:,ke_z,ke2D) = RHOCV(:,ke_z,ke2D) / DENS(:,ke_z,ke2D)
    end do
    end do

    return
  end subroutine atm_phy_cloud_dgm_common_sedimentation

!OCL SERIAL
  subroutine atm_phy_cloud_dgm_common_sedimentation_momentum( &
    MOMU_t, MOMV_t, MOMZ_t,                & ! (out)
    DENS, MOMU, MOMV, MOMZ, mflx,          & ! (in)
    Dz, Lift, nz, vmapM, vmapP,            & ! (in)
    lcmesh, elem                           ) ! (in)
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: MOMU_t(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMV_t(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMZ_t(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMU(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMV(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMZ(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: mflx(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    real(RP), intent(in) :: nz(elem%NfpTot,lcmesh%NeZ,lcmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lcmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lcmesh%NeZ)
    
    integer :: ke2D
    integer :: ke_z
    integer :: ke

    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,3)

    real(RP) :: RDENS(elem%Np)
    !-------------------------------------------------------

    call atm_phy_cloud_dgm_sedimentation_momentum_get_delflux( &
      del_flux(:,:,:,:),                                  & ! (out)
      DENS(:,:,:), MOMU(:,:,:), MOMV(:,:,:), MOMZ(:,:,:), & ! (in)
      mflx(:,:,:),                                        & ! (in)
      nz(:,:,:), vmapM(:,:), vmapP(:,:),                  & ! (in)
      lcmesh, elem                                        ) ! (in)

    !$omp parallel do collapse(2) private( &
    !$omp ke2D, ke_z, ke, RDENS, Fz, LiftDelFlx )
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      RDENS(:) = 1.0_RP / DENS(:,ke_z,ke2D)

      call sparsemat_matmul( Dz, mflx(:,ke_z,ke2D) * MOMU(:,ke_z,ke2D) * RDENS(:), Fz )
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,1), LiftDelFlx )
      MOMU_t(:,ke) = - ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )

      call sparsemat_matmul( Dz, mflx(:,ke_z,ke2D) * MOMV(:,ke_z,ke2D) * RDENS(:), Fz )
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,2), LiftDelFlx )
      MOMV_t(:,ke) = - ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )

      call sparsemat_matmul( Dz, mflx(:,ke_z,ke2D) * MOMZ(:,ke_z,ke2D) * RDENS(:), Fz )
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,3), LiftDelFlx )
      MOMZ_t(:,ke) = - ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )
    end do
    end do

    return
  end subroutine atm_phy_cloud_dgm_common_sedimentation_momentum

  !> Calculate a state after the precipitation process
  !!
!OCL SERIAL
  subroutine atm_phy_cloud_dgm_common_condensate_removal( &
    DENS, RHOQ, CPtot, CVtot, RHOE,            & ! (inout)
    sflx_rain, sflx_snow, esflx,               & ! (inout)
    TEMP, dt,                                  & ! (in)
    QHA, QLA, QIA, QS, lcmesh, elem, elem1D,   & ! (in)
    do_evap_precip, tau_evap_precip,           & ! (in, optional)
    RHOQV                                      ) ! (inout, optional)
    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_pres2qsat_liq
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: QHA                   !< hydrometeor (water + ice)
    real(RP), intent(inout) :: DENS (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: RHOQ (elem%Np,lcmesh%NeZ,lcmesh%Ne2D,QHA)
    real(RP), intent(inout) :: CPtot(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: CVtot(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: RHOE (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: sflx_rain(elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(inout) :: sflx_snow(elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(inout) :: esflx    (elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(in) :: TEMP (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: dt
    integer, intent(in) :: QLA, QIA
    integer, intent(in) :: QS        !< Start index for hydrometeor loop
    class(ElementBase1D), intent(in) :: elem1D
    logical, intent(in), optional :: do_evap_precip
    real(RP), intent(in), optional :: tau_evap_precip
    real(RP), intent(inout), optional :: RHOQV(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)

    real(RP) :: RHOCP(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOCV(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)

    real(RP) :: dDENS(elem%Np)
    real(RP) :: dInternalEn(elem%Np)

    real(RP) :: vint_weight(elem%Nnode_h1D**2,elem%Nnode_v)
    real(RP) :: r_dz

    real(RP) :: precip_mass_lc  (elem%Nnode_h1D**2)
    real(RP) :: precip_energy_lc(elem%Nnode_h1D**2)

    real(RP) :: eflx(elem%Np)
    real(RP) :: CP(QHA), CV(QHA)

    integer :: ke2D, ke_z, ke
    integer :: p2D, pz, p
    integer :: iq

    real(RP) :: rdt

    logical :: l_do_evap_precip
    real(RP) :: l_tau_evap_precip
    real(RP) :: relax_factor
    real(RP) :: dens_work, qsat_work
    real(RP) :: dm_evap, drho_evap, precip_energy_evap
    !-------------------------------------------------------

    call set_hydrometeor_heat_capacity( QLA, QIA, QHA, & ! (in)
      CP, CV) ! (out)

    rdt = 1.0_RP / dt

    ! Set local flags for optional evaporation and precipitation parameters

    if ( present(do_evap_precip) ) then
      l_do_evap_precip = do_evap_precip
    else
      l_do_evap_precip = .false.
    end if

    if ( l_do_evap_precip .and. ( .not. present(RHOQV) ) ) then
      LOG_INFO("atm_phy_cloud_dgm_common_condensate_removal",*) "RHOQV must be present when do_evap_precip is true. Check!"
      call PRC_abort
    end if

    if ( l_do_evap_precip ) then
      if ( present(tau_evap_precip) ) then
        l_tau_evap_precip = tau_evap_precip
      else
        LOG_INFO("atm_phy_cloud_dgm_common_condensate_removal",*) "tau_evap_precip must be present when do_evap_precip is true. Check!"
        call PRC_abort
      end if
    end if

    ! Set relaxation factor for evaporation and precipitation
    if ( l_do_evap_precip ) then
      relax_factor = 1.0_RP - exp(- dt / l_tau_evap_precip)
    else
      relax_factor = 0.0_RP
    end if

    ! Calculate initial local values for specific heat and density
    
    !$omp parallel do collapse(2)
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ
      RHOCP(:,ke_z,ke2D) = CPtot(:,ke_z,ke2D) * DENS(:,ke_z,ke2D)
      RHOCV(:,ke_z,ke2D) = CVtot(:,ke_z,ke2D) * DENS(:,ke_z,ke2D)
    end do
    end do

    ! Process each hydrometeor category separately
    do iq = QS, QLA + QIA

      !$omp parallel do private(ke2D,ke_z,ke,p2D,pz,p, &
      !$omp dDENS, vint_weight, r_dz, precip_mass_lc, dInternalEn, precip_energy_lc, &
      !$omp dens_work, qsat_work, dm_evap, drho_evap, precip_energy_evap             )
      do ke2D = 1, lcmesh%Ne2D

        precip_mass_lc(:) = 0.0_RP
        precip_energy_lc(:) = 0.0_RP

        do ke_z = lcmesh%NeZ, 1, -1
          ke = ke2D + (ke_z-1)*lcmesh%Ne2D

          dDENS      (:) = - RHOQ(:,ke_z,ke2D,iq)
          dInternalEn(:) = - RHOQ(:,ke_z,ke2D,iq) * CV(iq) * TEMP(:,ke_z,ke2D)

          RHOQ(:,ke_z,ke2D,iq) = 0.0_RP

          RHOCP(:,ke_z,ke2D) = RHOCP(:,ke_z,ke2D) + CP(iq) * dDENS(:)
          RHOCV(:,ke_z,ke2D) = RHOCV(:,ke_z,ke2D) + CV(iq) * dDENS(:)

          !-- Accumulate condensate and related internal energy vertically

          do pz=1, elem%Nnode_v
          do p2D=1, elem%Nnode_h1D**2
            vint_weight(p2D,pz) = 0.5_RP * elem1D%IntWeight_lgl(pz) &
                                * ( lcmesh%zlev(elem%Colmask(elem%Nnode_v,p2D),ke) - lcmesh%zlev(elem%Colmask(1,p2D),ke) )
          end do
          end do

          do pz=elem%Nnode_v, 1, -1
            do p2D=1, elem%Nnode_h1D**2
              p = elem%Colmask(pz,p2D)
              precip_mass_lc  (p2D) = precip_mass_lc  (p2D) - vint_weight(p2D,pz) * dDENS(p)
              precip_energy_lc(p2D) = precip_energy_lc(p2D) - vint_weight(p2D,pz) * dInternalEn(p)

              !- Evaporation of diagnostic liquid precipitation
              ! Here evaporation is applied only to liquid precipitation.
              ! Ice/snow should preferably use qsat_ice and sublimation latent heat separately.            
              if ( l_do_evap_precip .and. iq <= QLA ) then
                if ( precip_mass_lc(p2D) > 0.0_RP ) then
                  dens_work = DENS(p,ke_z,ke2D) + dDENS(p)
                  call ATMOS_SATURATION_pres2qsat_liq( &
                    TEMP(p,ke_z,ke2D), dens_work, &
                    qsat_work )
                  
                  dm_evap = max(qsat_work - RHOQV(p,ke_z,ke2D) / dens_work, 0.0_RP) * relax_factor &
                          * dens_work * vint_weight(p2D,pz)
                  dm_evap = min(dm_evap, precip_mass_lc(p2D))

                  if ( dm_evap > 0.0_RP ) then
                    r_dz = 1.0_RP / vint_weight(p2D,pz)
                    drho_evap = dm_evap * r_dz

                    precip_energy_evap = dm_evap / precip_mass_lc(p2D) * precip_energy_lc(p2D)

                    precip_mass_lc  (p2D) = precip_mass_lc  (p2D) - dm_evap
                    precip_energy_lc(p2D) = precip_energy_lc(p2D) - precip_energy_evap

                    dDENS      (p) = dDENS      (p) + drho_evap
                    dInternalEn(p) = dInternalEn(p) + precip_energy_evap * r_dz &
                                  - LHV * drho_evap

                    RHOQV(p,ke_z,ke2D) = RHOQV(p,ke_z,ke2D) + drho_evap
                    RHOCP(p,ke_z,ke2D) = RHOCP(p,ke_z,ke2D) + CP_VAPOR * drho_evap
                    RHOCV(p,ke_z,ke2D) = RHOCV(p,ke_z,ke2D) + CV_VAPOR * drho_evap
                  end if
                end if
              end if

            end do
          end do

          !--- Update density and internal energy

          DENS (:,ke_z,ke2D) = DENS(:,ke_z,ke2D) + dDENS(:)
          RHOE(:,ke_z,ke2D) = RHOE(:,ke_z,ke2D) + dInternalEn(:)
        end do ! ke_z loop

        ! Note that the sign is negative for downward flux              
        if ( iq > QLA ) then ! ice water
          sflx_snow(:,ke2D) = sflx_snow(:,ke2D) - precip_mass_lc(:) * rdt
        else                 ! liquid water
          sflx_rain(:,ke2D) = sflx_rain(:,ke2D) - precip_mass_lc(:) * rdt
        end if
        esflx(:,ke2D) = esflx(:,ke2D) - precip_energy_lc(:) * rdt

      end do ! ke2D loop
    end do ! iq loop

    !$omp parallel do collapse(2)
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ  
      CPtot(:,ke_z,ke2D) = RHOCP(:,ke_z,ke2D) / DENS(:,ke_z,ke2D)
      CVtot(:,ke_z,ke2D) = RHOCV(:,ke_z,ke2D) / DENS(:,ke_z,ke2D)
    end do
    end do

    return
  end subroutine atm_phy_cloud_dgm_common_condensate_removal

  !> Calculate a tendency of momentum due to the condensate removal process
  !!
!OCL SERIAL
  subroutine atm_phy_cloud_dgm_common_condensate_removal_momentum( &
    MOMU_t, MOMV_t, MOMZ_t,                & ! (out)
    DENS, MOMU, MOMV, MOMZ, DENS_new,      & ! (in)
    rdt_MP, lcmesh, elem                   ) ! (in)
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: MOMU_t(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMV_t(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMZ_t(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMU(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMV(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMZ(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: DENS_new(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: rdt_MP
    
    integer :: ke2D, ke_z, ke
    real(RP) :: coef(elem%Np)
    !----------------------------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke2D, ke_z, ke, coef )
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      coef(:) = ( DENS_new(:,ke_z,ke2D) / DENS(:,ke_z,ke2D) - 1.0_RP ) * rdt_MP

      MOMU_t(:,ke) = coef(:) * MOMU(:,ke_z,ke2D)
      MOMV_t(:,ke) = coef(:) * MOMV(:,ke_z,ke2D)
      MOMZ_t(:,ke) = coef(:) * MOMZ(:,ke_z,ke2D)
    end do
    end do
    return
  end subroutine atm_phy_cloud_dgm_common_condensate_removal_momentum

!OCL SERIAL
  subroutine atm_phy_cloud_dgm_common_negative_fixer( &
    QTRC, DDENS, PRES,                       &
    CVtot, CPtot, Rtot,                      &
    DENS_hyd, PRES_hyd,                      &
    dt, lmesh, elem, QA, QLA, QIA,           &
    DRHOT                                    )

    use scale_const, only: &
      CVdry => CONST_CVdry,  &
      CPdry => CONST_CPdry,  &
      Rdry => CONST_Rdry

    use scale_tracer, only: &
      TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat
    use scale_localmeshfield_base, only: LocalMeshFieldBaseList
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: QA
    type(LocalMeshFieldBaseList), intent(inout) :: QTRC(QA)
    real(RP), intent(inout) :: DDENS(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: PRES(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: dt
    integer, intent(in) :: QLA, QIA
    real(RP), intent(inout), optional :: DRHOT(elem%Np,lmesh%NeA)

    integer :: ke
    integer :: iq

    real(RP) ::  int_w(elem%Np)

    real(RP) :: DENS(elem%Np)
    real(RP) :: DDENS0(elem%Np)

    real(RP) :: TRCMASS0(elem%Np), TRCMASS1(elem%Np,QA)
    real(RP) :: MASS0_elem, MASS1_elem
    real(RP) :: IntEn0_elem, IntEn_elem

    real(RP) :: QTRC_tmp(elem%Np,QA), Qdry(elem%Np)
    real(RP) :: CVtot_old(elem%Np), CPtot_old(elem%Np), Rtot_old(elem%Np)
    real(RP) :: InternalEn(elem%Np), InternalEn0(elem%Np), TEMP(elem%Np)
    real(RP) :: RHOT_hyd(elem%Np)

#ifdef SINGLE
    real(RP), parameter :: TRC_EPS = 1E-32_RP
#else
    real(RP), parameter :: TRC_EPS = 1E-128_RP
#endif
    !------------------------------------------------

    !$omp parallel do private( &
    !$omp ke, iq, DENS, DDENS0, InternalEn, InternalEn0, TEMP, QTRC_tmp,       &
    !$omp TRCMASS0, TRCMASS1, MASS0_elem, MASS1_elem, IntEn0_elem, IntEn_elem, &
    !$omp Qdry, CVtot_old, CPtot_old, Rtot_old,                                &
    !$omp int_w, RHOT_hyd )
    do ke = lmesh%NeS, lmesh%NeE

      do iq = 1, QA
        QTRC_tmp(:,iq) = QTRC(iq)%ptr%val(:,ke)
      end do
      call ATMOS_THERMODYN_specific_heat( & 
        elem%Np, 1, elem%Np, QA,                                           & ! (in)
        QTRC_tmp, TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! (in)
        Qdry, Rtot_old, CVtot_old, CPtot_old                               ) ! (out)

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      DDENS0(:) = DDENS(:,ke)
      
      ! RHOT_hyd(:) = PRES00 / Rdry * ( PRES_hyd(:,ke) / PRES00 )**( CVdry / CPdry )
      ! ( Internal energy ) = Cvtot * RHO * T = CVtot * RHO * ( PT * EXNER )      
      ! InternalEn0(:) = CVtot_old(:) * ( RHOT_hyd(:) + DRHOT(:,ke) ) &
      !                * ( Rtot_old(:) * ( RHOT_hyd(:) + DRHOT(:,ke) ) / PRES00 )**( Rtot_old(:) / CVtot_old(:) ) 
      !InternalEn0(:) = CVtot_old(:) * PRES(:,ke) / Rtot_old(:)
      InternalEn0(:) = CVtot(:,ke) * PRES(:,ke) / Rtot(:,ke)
      !TEMP(:) = InternalEn0(:) / ( DENS(:) * CVtot_old(:) )
      TEMP(:) = InternalEn0(:) / ( DENS(:) * CVtot(:,ke) )
      InternalEn(:) = InternalEn0(:)

      int_w(:) = lmesh%Gsqrt(:,ke) * lmesh%J(:,ke) * elem%IntWeight_lgl(:)
      do iq = 1, 1 + QLA + QIA  
        TRCMASS0(:) = DENS(:) * QTRC_tmp(:,iq)
        TRCMASS1(:,iq) = max( TRC_EPS, TRCMASS0(:) )

        MASS0_elem = sum( int_w(:) * TRCMASS0(:)    )
        MASS1_elem = sum( int_w(:) * TRCMASS1(:,iq) )
        TRCMASS1(:,iq) = max(MASS0_elem, 0.0E0_RP) / MASS1_elem * TRCMASS1(:,iq)
!        TRCMASS1(:,iq) = MASS0_elem / MASS1_elem * TRCMASS1(:,iq)

        DDENS(:,ke) = DDENS(:,ke) + ( TRCMASS1(:,iq) - TRCMASS0(:) )
        InternalEn(:) = InternalEn(:) + ( TRCMASS1(:,iq) - TRCMASS0(:) ) * TRACER_CV(iq) * TEMP(:)
      end do

      !--

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      do iq = 1, QA
        QTRC_tmp(:,iq) = TRCMASS1(:,iq) / DENS(:)
        QTRC(iq)%ptr%val(:,ke) = QTRC_tmp(:,iq)
      end do

      IntEn0_elem = sum( int_w(:) * InternalEn0(:) )
      IntEn_elem  = sum( int_w(:) * InternalEn (:) )
      InternalEn(:) = IntEn0_elem / IntEn_elem * InternalEn(:)

      call ATMOS_THERMODYN_specific_heat( &
        elem%Np, 1, elem%Np, QA,                                           & ! (in)
        QTRC_tmp, TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! (in)
        Qdry, Rtot(:,ke), CVtot(:,ke), CPtot(:,ke)                         ) ! (out)

      InternalEn(:) = InternalEn(:) - ( DDENS(:,ke) - DDENS0(:) ) * Grav * lmesh%zlev(:,ke)
      PRES(:,ke) = InternalEn(:) * Rtot(:,ke) / CVtot(:,ke)

      if ( present(DRHOT) ) then
        DRHOT(:,ke) = PRES00 / Rtot(:,ke) * ( PRES(:,ke) / PRES00 )**( CVtot(:,ke) / CPtot(:,ke) ) &
                    - PRES00 / Rdry * ( PRES_hyd(:,ke) / PRES00 )**( CVdry / CPdry )
      end if
    end do

    return
  end subroutine atm_phy_cloud_dgm_common_negative_fixer

!- private ---------------------------------------------------------------------

!OCL SERIA
  subroutine set_hydrometeor_heat_capacity( QLA, QIA, QHA, &
    CP, CV )
    implicit none
    integer, intent(in) :: QLA, QIA
    integer, intent(in) :: QHA
    real(RP), intent(out) :: CP(QHA)
    real(RP), intent(out) :: CV(QHA)

    integer :: iq
    !----------------------------------------

    do iq = 1, QHA
      if ( iq > QLA + QIA ) then
        CP(iq) = UNDEF 
        CV(iq) = UNDEF
      else if ( iq > QLA ) then ! ice water
        CP(iq) = CP_ICE
        CV(iq) = CV_ICE
      else                      ! liquid water
        CP(iq) = CP_WATER
        CV(iq) = CV_WATER
      end if
    end do

    return
  end subroutine set_hydrometeor_heat_capacity

!OCL SERIAL
  subroutine atm_phy_cloud_dgm_netOutwardFlux( &
    net_outward_flux,                     &
    RHOQ_, vterm_,                        &
    DzRHOQ_, NDcoefEuler_,                &
    J, Fscale,                            &
    nz, vmapM, vmapP, vmapM3D, IntWeight, &
    lmesh, elem                           )
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: net_outward_flux(lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: RHOQ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: vterm_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: DzRHOQ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: NDcoefEuler_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: J(elem%Np*lmesh%Ne)
    real(RP), intent(in) :: Fscale(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM3D(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: IntWeight(elem%Nfaces,elem%NfpTot)

    real(RP) :: numflux(elem%NfpTot)
    real(RP) :: outward_flux_tmp(elem%Nfaces)    
    real(RP) :: alpha(elem%NfpTot)
    real(RP) :: velM(elem%NfpTot), velP(elem%NfpTot)
    real(RP) :: RHOQ_M(elem%NfpTot), RHOQ_P(elem%NfpTot)

    integer :: ke
    integer :: ke_z, ke2D
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: iM3D(elem%NfpTot)
    !------------------------------------------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke2D, ke_z, ke, iM3D, iM, iP,      &
    !$omp RHOQ_M, RHOQ_P, velM, velP, alpha, &
    !$omp numflux, outward_flux_tmp          )
    do ke2D=1, lmesh%Ne2D
    do ke_z=1, lmesh%NeZ
      ke = ke2D + (ke_z-1)*lmesh%Ne2D

      iM3D(:) = vmapM3D(:,ke)
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      RHOQ_M(:) = RHOQ_(iM(:),ke2D)
      RHOQ_P(:) = RHOQ_(iP(:),ke2D)
      velM(:) = vterm_(iM(:),ke2D) * nz(:,ke_z,ke2D)
      velP(:) = vterm_(iP(:),ke2D) * nz(:,ke_z,ke2D)
      alpha(:) = nz(:,ke_z,ke2D)**2 * max( abs(velM(:)), abs(velP(:)) )

      where (nz(:,ke_z,ke2D) > 1.0E-10 .and. iP(:) == iM(:) )
        velP(:) = - velM(:)
      end where      

      numflux(:) = 0.5_RP * (  RHOQ_P(:) * velP(:) + RHOQ_M(:) * velM(:)                                                        &
        - ( NDcoefEuler_(iP(:),ke2D) * DzRHOQ_(iP(:),ke2D) + NDcoefEuler_(iM(:),ke2D) * DzRHOQ_(iM(:),ke2D) ) * nz(:,ke_z,ke2D) &
        - alpha(:) * ( RHOQ_P(:) - RHOQ_M(:) )                                                                                  )

      outward_flux_tmp(:) = matmul( IntWeight(:,:), J(iM3D(:)) * Fscale(:,ke) * numflux(:) )
      net_outward_flux(ke_z,ke2D) = sum( max( 0.0_RP, outward_flux_tmp(:) ) )
    end do
    end do

    return
  end subroutine atm_phy_cloud_dgm_netOutwardFlux

!OCL SERIAL
  subroutine atm_phy_cloud_dgm_sedimentation_get_delflux_dq( &
    del_flux,                                             &
    DENS_, RHOQ_,TEMP_, CV, nz, vmapM, vmapP,             &
    lmesh, elem                                           )
    
    implicit none
    
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: del_flux(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D,2)
    real(RP), intent(in) :: DENS_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: RHOQ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: TEMP_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: CV
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)

    integer :: ke
    integer :: ke_z, ke2D
    real(RP) :: RHOQ_P(elem%NfpTot), RHOQ_M(elem%NfpTot)
    integer :: iM(elem%NfpTot), iP(elem%NfpTot)
    !-----------------------------------------

    !$omp parallel do collapse(2) private(       &
    !$omp ke2D, ke_z, ke, iM, iP, RHOQ_M, RHOQ_P )
    do ke2D=1, lmesh%Ne2D
    do ke_z=1, lmesh%NeZ
      ke = ke2D + (ke_z-1)*lmesh%Ne2D
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      RHOQ_M(:) = RHOQ_(iM(:),ke2D) 
      RHOQ_P(:) = RHOQ_(iP(:),ke2D) 
      del_flux(:,ke_z,ke2D,1) = 0.5_RP * ( RHOQ_P(:) - RHOQ_M(:) ) * nz(:,ke_z,ke2D)
      del_flux(:,ke_z,ke2D,2) = 0.5_RP * CV * ( RHOQ_P(:) * TEMP_(iP(:),ke2D) - RHOQ_M(:) * TEMP_(iM(:),ke2D) ) * nz(:,ke_z,ke2D)
    end do
    end do
    
    return
  end subroutine atm_phy_cloud_dgm_sedimentation_get_delflux_dq


!OCL SERIAL
  subroutine atm_phy_cloud_dgm_sedimentation_get_delflux( &
    del_flux,                                          &
    DENS_, RHOQ_, TEMP_, vterm_,                       &
    DzRHOQ_, DzRHOE_, NDcoefEuler_,                    &
    fct_coef_, CV,                                     &
    J, Fscale, nz, vmapM, vmapP, vmapM3D, IntWeight,   &
    lmesh, elem                                        )

    implicit none
    
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: del_flux(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D,2)
    real(RP), intent(in) :: DENS_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: RHOQ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: TEMP_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: vterm_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: DzRHOQ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: DzRHOE_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: NDcoefEuler_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: fct_coef_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: CV
    real(RP), intent(in) :: J(elem%Np*lmesh%Ne)
    real(RP), intent(in) :: Fscale(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM3D(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: IntWeight(elem%Nfaces,elem%NfpTot)

    integer :: ke
    integer :: ke_z, ke2D
    integer :: f, p, fp
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: alpha(elem%NfpTot)
    real(RP) :: velM(elem%NfpTot), velP(elem%NfpTot)
    real(RP) :: RHOQ_M(elem%NfpTot), RHOQ_P(elem%NfpTot)
    real(RP) :: TEMP_M(elem%NfpTot), TEMP_P(elem%NfpTot)

    integer :: iM3D(elem%NfpTot)
    real(RP) :: R_M(elem%NfpTot), R_P(elem%NfpTot)
    real(RP) :: NDcoef_M(elem%NfpTot), NDcoef_P(elem%NfpTot)
    real(RP) :: numflux   (elem%NfpTot)
    real(RP) :: numflux_ei(elem%NfpTot)
    real(RP) :: outward_flux_tmp(elem%Nfaces)
    real(RP) :: R
    !-----------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke2D, ke_z, ke, iM3D, iM, iP,                 &
    !$omp R_M, R_P, RHOQ_M, RHOQ_P, TEMP_M, TEMP_P,     &
    !$omp velM, velP, alpha, numflux, numflux_ei,       &
    !$omp NDcoef_M, NDcoef_P,                           &
    !$omp outward_flux_tmp,                             &
    !$omp f, p, fp, R                                   )
    do ke2D=1, lmesh%Ne2D
    do ke_z=1, lmesh%NeZ
      ke = ke2D + (ke_z-1)*lmesh%Ne2D

      iM3D(:) = vmapM3D(:,ke)
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      R_M(:) = fct_coef_(iM(:),ke2D)
      R_P(:) = fct_coef_(iP(:),ke2D) 
      RHOQ_M(:) = RHOQ_(iM(:),ke2D)
      RHOQ_P(:) = RHOQ_(iP(:),ke2D)
      TEMP_M(:) = TEMP_(iM(:),ke2D)
      TEMP_P(:) = TEMP_(iP(:),ke2D)

      velM(:) = vterm_(iM(:),ke2D) * nz(:,ke_z,ke2D)
      velP(:) = vterm_(iP(:),ke2D) * nz(:,ke_z,ke2D)
      alpha(:) = nz(:,ke_z,ke2D)**2 * max( abs(velM(:)), abs(velP(:)) )

      where (nz(:,ke_z,ke2D) > 1.0E-10 .and. iP(:) == iM(:) )
        velP(:) = - velM(:)
      end where      

      NDcoef_M(:) = NDcoefEuler_(iM(:),ke2D)
      NDcoef_P(:) = NDcoefEuler_(iP(:),ke2D)

      numflux(:) = 0.5_RP * (  RHOQ_P(:) * velP(:) + RHOQ_M(:) * velM(:)                              &
        - ( NDcoef_P(:) * DzRHOQ_(iP(:),ke2D) + NDcoef_M(:) * DzRHOQ_(iM(:),ke2D) ) * nz(:,ke_z,ke2D) &
        - alpha(:) * ( RHOQ_P(:) - RHOQ_M(:) )                                                        )

      numflux_ei(:) = 0.5_RP * ( CV * ( RHOQ_P(:) * TEMP_P(:) * velP(:) + RHOQ_M(:) * TEMP_M(:) * velM(:) ) &
        - ( NDcoef_P(:) * DzRHOE_(iP(:),ke2D) + NDcoef_M(:) * DzRHOE_(iM(:),ke2D) ) * nz(:,ke_z,ke2D)       &
        - alpha(:) * CV * ( RHOQ_P(:) * TEMP_P(:) - RHOQ_M(:) * TEMP_M(:) )                                 )


      del_flux(:,ke_z,ke2D,1) = 0.0_RP
      del_flux(:,ke_z,ke2D,2) = 0.0_RP  
      outward_flux_tmp(:) = matmul( IntWeight(:,:), J(iM3D(:)) * Fscale(:,ke) * numflux(:) )
      do f=1, elem%Nfaces_v
      do p=1, elem%Nfp_v
        fp = p + (f-1)*elem%Nfp_v + elem%Nfaces_h * elem%Nfp_h
        R = 0.5_RP * ( R_P(fp) + R_M(fp) - ( R_P(fp) - R_M(fp) ) * sign( 1.0_RP, outward_flux_tmp(elem%Nfaces_h+f) ) )
        del_flux(fp,ke_z,ke2D,1) = numflux   (fp) * R  &
                                 - RHOQ_M(fp) * velM(fp)                                  &
                                 + NDcoef_M(fp) * DzRHOQ_(iM(fp),ke2D) * nz(fp,ke_z,ke2D)
        del_flux(fp,ke_z,ke2D,2) = numflux_ei(fp) * R  &
                                 - RHOQ_M(fp) * velM(fp) * CV * TEMP_M(fp)                &
                                 + NDcoef_M(fp) * DzRHOE_(iM(fp),ke2D) * nz(fp,ke_z,ke2D)
      end do
      end do                     

    end do
    end do

    return
  end subroutine atm_phy_cloud_dgm_sedimentation_get_delflux

!OCL SERIAL
  subroutine atm_phy_cloud_dgm_sedimentation_momentum_get_delflux( &
    del_flux,                          & ! (out)
    DENS_, MOMU_, MOMV_, MOMZ_, mflx_, & ! (in)
    nz, vmapM, vmapP, lmesh, elem      ) ! (in)

    implicit none
    
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: del_flux(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D,3)
    real(RP), intent(in) :: DENS_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: MOMU_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: MOMV_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: MOMZ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: mflx_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)

    integer :: ke_z, ke2D
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: alpha(elem%NfpTot)
    real(RP) :: densM(elem%NfpTot), densP(elem%NfpTot)
    real(RP) :: VelM(elem%NfpTot), VelP(elem%NfpTot)
    !-----------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke2D, ke_z, iM, iP, alpha,       &
    !$omp densM, densP, VelM, VelP         )
    do ke2D=1, lmesh%Ne2D
    do ke_z=1, lmesh%NeZ
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      densM(:) = DENS_(iM(:),ke2D)
      densP(:) = DENS_(iP(:),ke2D)
      VelM(:) = mflx_(iM(:),ke2D) * nz(:,ke_z,ke2D) / densM(:)
      VelP(:) = mflx_(iP(:),ke2D) * nz(:,ke_z,ke2D) / densP(:)
      alpha(:) = nz(:,ke_z,ke2D)**2 * max( abs(VelM(:)), abs(VelP(:)) )

      del_flux(:,ke_z,ke2D,1) = 0.5_RP * ( &
        + ( MOMU_(iP(:),ke2D) * VelP - MOMU_(iM(:),ke2D) * VelM ) &
        - alpha(:) * ( MOMU_(iP(:),ke2D) - MOMU_(iM(:),ke2D) )    )
      
      del_flux(:,ke_z,ke2D,2) = 0.5_RP * ( &
        + ( MOMV_(iP(:),ke2D) * VelP - MOMV_(iM(:),ke2D) * VelM ) &
        - alpha(:) * ( MOMV_(iP(:),ke2D) - MOMV_(iM(:),ke2D) )    )

      del_flux(:,ke_z,ke2D,3) = 0.5_RP * ( &
        + ( MOMZ_(iP(:),ke2D) * VelP - MOMZ_(iM(:),ke2D) * VelM ) &
        - alpha(:) * ( MOMZ_(iP(:),ke2D) - MOMZ_(iM(:),ke2D) )    )
    end do
    end do

    return
  end subroutine atm_phy_cloud_dgm_sedimentation_momentum_get_delflux

end module scale_atm_phy_cloud_dgm_common
