!> module FElib / Atmosphere / Physics cumulus parameterization
!!
!! @par Description
!! A module providing a moist convective adjustment scheme.
!!
!! Saturated convective instability is diagnosed from the vertical gradient of saturated moist static energy in nearly saturated layers.
!! The diagnosed layer is adjusted toward a moist-neutral temperature profile subject to column-integrated moist-static-energy and available-water constraints. 
!! Condensed water is immediately removed as precipitation.
!!
!! @author Yuta Kawai, Team SCALE
!!
!! @par Reference
!!      Manabe, S., J. Smagorinsky, and R. F. Strickler (1965):
!!      Simulated Climatology of a General Circulation Model with a
!!      Hydrologic Cycle. Monthly Weather Review, 93, 769-798.
!!
!! 
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_cp_dgm_mconv_adjustment
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,   &
    LHV0 => CONST_LHV0,   &
    CPdry => CONST_CPDRY, &
    CVdry => CONST_CVdry, &
    Rdry => CONST_Rdry,   &
    Rvap => CONST_Rvap,   &
    PRES0 => CONST_PRE00, &
    EPS => CONST_EPS
  use scale_atmos_hydrometeor, only: &
    CV_VAPOR, &
    CP_VAPOR
  use scale_atmos_saturation, only: &
    ATMOS_SATURATION_psat_liq,     &
    ATMOS_SATURATION_pres2qsat_liq

  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
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
  !++ Public type & procedure
  !
  public :: atm_phy_cp_dgm_mconv_adjustment_setup
  public :: atm_phy_cp_dgm_mconv_adjustment_calc_tendency
  public :: atm_phy_cp_dgm_mconv_adjustment_finalize

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  
  !* MCA paramters

  ! Iteration limits
  integer, parameter :: MCA_MAX_ADJUST_ITER = 200
  integer, parameter :: MCA_MAX_ENERGY_ITER = 50
  integer, parameter :: MCA_MAX_WATER_ITER  = 60

  ! Column-energy root solve
  real(RP), parameter :: MCA_COLUMN_ENERGY_RTOL   = 5.0E-8_RP
  real(RP), parameter :: MCA_COLUMN_ENERGY_ATOL   = 1.0_RP
  real(RP), parameter :: MCA_BOUNDARY_ENERGY_RTOL = 1.0E-7_RP

  ! Local saturated-MSE solve for moist-adiabat integration
  real(RP), parameter :: MCA_LOCAL_HMSE_RTOL = 1.0E-11_RP
  real(RP), parameter :: MCA_LOCAL_HMSE_ATOL = 1.0E-8_RP
  real(RP), parameter :: MCA_TEMP_TOL        = 1.0E-8_RP

  ! Water-feasibility constraint
  real(RP), parameter :: MCA_WATER_RTOL = 1.0E-12_RP
  real(RP), parameter :: MCA_WATER_ATOL = 1.0E-12_RP

  ! Saturation and instability thresholds
  real(RP), parameter :: MCA_RH_FORCED_SATURATION = 0.999_RP
  real(RP), parameter :: MCA_RH_TRIGGER           = 0.999_RP

  ! Thresholds for suppressing numerical-noise detection
  real(RP), parameter :: MCA_HMSE_GRAD_TOL      = 1.0E-3_RP
  real(RP), parameter :: MCA_HMSE_DIFF_TOL      = 1.0_RP
  real(RP), parameter :: MCA_MIN_UNSTABLE_DEPTH = 0.0_RP
  real(RP), parameter :: MCA_Z_TOL              = 100.0_RP * epsilon(1.0_RP)

  ! Moist-neutral-profile solver status
  integer, parameter :: MCA_PROFILE_SUCCESS                 = 0
  integer, parameter :: MCA_PROFILE_NO_ENERGY_ROOT          = 1
  integer, parameter :: MCA_PROFILE_ADIABAT_FAILURE         = 2
  integer, parameter :: MCA_PROFILE_ENERGY_MAXITER          = 3
  integer, parameter :: MCA_PROFILE_NO_ADJUSTED_NODE        = 4
  integer, parameter :: MCA_PROFILE_NO_WATER_FEASIBLE_STATE = 5

contains
  !> Setup a module for moist convective adjustment scheme
  subroutine atm_phy_cp_dgm_mconv_adjustment_setup()
    implicit none
    !----------------------------------------------------
    return
  end subroutine atm_phy_cp_dgm_mconv_adjustment_setup

  !> Calculate tendencies with moist convective adjustment scheme
!OCL SERIAL
  subroutine atm_phy_cp_dgm_mconv_adjustment_calc_tendency( &
    DENS_t, RHOT_t, RHOQV_t, SFLX_RAIN,                & ! (out)
    DDENS, DRHOT, QV, PT, PRES, DENS_hyd, Rtot, CPtot, & ! (in)
    dtsec, lmesh, elem, elem1D                         ) ! (in)
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(ElementBase1D), intent(in) :: elem1D
    real(RP), intent(out) :: DENS_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOQV_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: SFLX_RAIN(elem%Nnode_h1D**2,lmesh%Ne2DA)
    real(RP), intent(in) :: DDENS(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: QV(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PT(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PRES(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: dtsec

    integer :: ke, ke_xy, ke_z
    integer :: ph, pz, p
    real(RP) :: dens_z(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: pres_z(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: zlev_z(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: temp_z(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: pott_z(elem%Nnode_v,lmesh%NeZ) ! work array
    real(RP) :: qvap_z(elem%Nnode_v,lmesh%NeZ) ! work array
    real(RP) :: qsat
    real(RP) :: rh_z
    real(RP) :: rtot_, cptot_, qdry

    real(RP) :: elem_width_z
    real(RP) :: int_weight(elem%Nnode_v,lmesh%NeZ)

    real(RP) :: pott_ini(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: qvap_ini(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: dens_ini(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: rhoqvap_ini(elem%Nnode_v,lmesh%NeZ)

    integer :: iter_expand
    integer :: MCA_MAX_MASK_EXPAND
    integer :: lbase, ltop, lbase_try, ltop_try
    real(RP) :: zbase_try, ztop_try
    logical :: mask_expanded

    logical :: is_unstable
    logical  :: unstable_core_mask(elem%Nnode_v,lmesh%NeZ)    
    logical :: forced_saturation_mask(elem%Nnode_v,lmesh%NeZ)
    logical :: adjustment_mask(elem%Nnode_v,lmesh%NeZ)
    
    integer :: iter_adj
    real(RP) :: temp_adj(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: qvap_adj(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: rhoqvap_adj(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: rhoprecip_adj(elem%Nnode_v,lmesh%NeZ)    ! Water-vapor mass density diagnosed to be removed as precipitation during the current adjustment iteration [kg m-3]
    real(RP) :: rhoprecip_accum(elem%Nnode_v,lmesh%NeZ)  ! Accumulated vapor mass density removed as precipitation [kg m-3]
    real(RP) :: rhoqdry_adj(elem%Nnode_v,lmesh%NeZ)

    real(RP), allocatable :: zlev_diag(:)
    integer :: nlev_diag

    logical :: profile_converged
    integer :: profile_status
    logical :: do_adjustment
    logical :: adjustment_converged

    integer :: Hslice_b(elem%Nnode_h1D**2), Hslice_t(elem%Nnode_h1D**2)

    logical :: debug_flag

    integer  :: lactive_base
    integer  :: lactive_top
    logical  :: active_region_initialized    
    !----------------------------------------------------

    nlev_diag = lmesh%NeZ * (elem%Nnode_v-1) + 1
    allocate( zlev_diag(nlev_diag) )
    
    Hslice_b(:) = elem%Hslice(:,1)
    Hslice_t(:) = elem%Hslice(:,elem%Nnode_v)

    MCA_MAX_MASK_EXPAND = nlev_diag

    !$omp parallel do collapse(2) private(ke,p, &
    !$omp dens_z, pres_z, temp_z, zlev_z, pott_z, qvap_z, qsat, rh_z, rtot_,cptot_,                 &
    !$omp dens_ini, pott_ini, qvap_ini, rhoqvap_ini,                                                &
    !$omp zlev_diag,                                                                                &
    !$omp lbase, ltop, lbase_try, ltop_try, zbase_try, ztop_try, iter_expand, mask_expanded,        &
    !$omp is_unstable, unstable_core_mask, forced_saturation_mask, adjustment_mask,                 &
    !$omp temp_adj, qdry, qvap_adj, rhoqdry_adj, rhoqvap_adj, rhoprecip_adj, rhoprecip_accum,       &
    !$omp profile_converged, profile_status, do_adjustment, adjustment_converged,                   &
    !$omp int_weight,elem_width_z, debug_flag, lactive_base, lactive_top, active_region_initialized )
    do ke_xy=1, lmesh%Ne2D
    do ph=1, elem%Nnode_h1D**2
    
      ! if ( lmesh%PRC_myrank == 1 .and. ph == 8 .and. ke_xy == 8 ) then
      !   debug_flag = .true.
      ! else
        debug_flag = .false.
      ! end if

      !- Extract vertical 1D DG column

      do ke_z=1, lmesh%NeZ
        ke = ke_xy + (ke_z-1)*lmesh%Ne2D

        elem_width_z = lmesh%zlev(Hslice_t(ph),ke) - lmesh%zlev(Hslice_b(ph),ke)
        do pz=1, elem%Nnode_v
          p = ph + (pz-1)*elem%Nnode_h1D**2
          dens_z(pz,ke_z) = DDENS(p,ke) + DENS_hyd(p,ke)
          pres_z(pz,ke_z) = PRES(p,ke)
          zlev_z(pz,ke_z) = lmesh%zlev(p,ke)
          int_weight(pz,ke_z) = elem1D%IntWeight_lgl(pz) * 0.5_RP * elem_width_z

          temp_z(pz,ke_z) = PRES(p,ke)/ ( dens_z(pz,ke_z) * Rtot(p,ke) )
          qvap_z(pz,ke_z) = QV(p,ke)
          pott_z(pz,ke_z) = PT(p,ke)
        end do
      end do

      !- Set initial state

      do ke_z=1, lmesh%NeZ
        dens_ini(:,ke_z) = dens_z(:,ke_z)
        pott_ini(:,ke_z) = pott_z(:,ke_z)
        qvap_ini(:,ke_z) = qvap_z(:,ke_z)
        rhoqvap_ini(:,ke_z) = dens_ini(:,ke_z) * qvap_ini(:,ke_z)

        temp_adj(:,ke_z) = temp_z(:,ke_z)
        rhoqvap_adj(:,ke_z) = rhoqvap_ini(:,ke_z)
        rhoprecip_accum(:,ke_z) = 0.0_RP
      end do

      !* Iteratively diagnose and remove remaining saturated convective instability.  ***************************************
      !  A single profile adjustment may generate a new unstable layer adjacent to or outside the previously adjusted region.

      adjustment_converged = .false.
      do_adjustment        = .false.

      active_region_initialized = .false.
      lactive_base = 0
      lactive_top  = 0
      forced_saturation_mask(:,:) = .false.

      do iter_adj=1, MCA_MAX_ADJUST_ITER

        qvap_adj(:,:) = rhoqvap_adj(:,:) / dens_z(:,:)

        ! Reconstruct the saturation constraint from the current state.
        forced_saturation_mask(:,:) = .false.        

        !- Diagnose unstable layers 

        call diagnose_convective_layer( lbase, ltop, is_unstable, zlev_diag,     & ! (out)
          temp_adj, qvap_adj, pres_z, zlev_z, lmesh, elem, nlev_diag, debug_flag ) ! (in)

        if ( debug_flag ) then
          write(*,*) "ke_xy=", ke_xy, "ph=", ph, "iter_adj=", iter_adj
          write(*,*) "   lbase, ltop, is_unstable = ", lbase, ltop, is_unstable
        end if

        if ( .not. is_unstable ) then
          adjustment_converged = .true.
          exit
        end if
        do_adjustment = .true.

        ! Maintain the smallest diagnostic-level interval containing all unstable or expanded regions encountered during this column adjustment. 
        if ( .not. active_region_initialized ) then
          lactive_base = lbase
          lactive_top  = ltop
          active_region_initialized = .true.
        else if ( lbase <= lactive_top + 1 .and. &
                  ltop  >= lactive_base - 1 ) then
          ! The new unstable layer overlaps or directly touches the previously active connected layer.
          lactive_base = min(lactive_base, lbase)
          lactive_top  = max(lactive_top,  ltop)        
        else
          ! Union of the previously adjusted region and the newly diagnosed unstable region.
          lactive_base = min(lactive_base, lbase)
          lactive_top  = max(lactive_top,  ltop)
        end if

        if ( debug_flag ) then
          write(*,*) "   persistent active region=", lactive_base, lactive_top
        end if        

        ! Mark nearly saturated DG nodes inside the diagnosed unstable core. 
        ! These nodes will be constrained to saturation in each trial moist-neutral profile.
        call make_dg_vertical_range_mask( &
          unstable_core_mask,                                 & ! (out)
          zlev_z, zlev_diag(lbase), zlev_diag(ltop), & ! (in)
          elem%Nnode_v, lmesh%NeZ                    ) ! (in)        

        do ke_z = 1, lmesh%NeZ
        do pz = 1, elem%Nnode_v
          if ( unstable_core_mask(pz,ke_z) ) then
            call ATMOS_SATURATION_pres2qsat_liq( temp_adj(pz,ke_z), pres_z(pz,ke_z),  &
              qsat ) ! (out)

            rh_z = qvap_adj(pz,ke_z) / max(qsat, EPS)
            forced_saturation_mask(pz,ke_z) = rh_z >= MCA_RH_FORCED_SATURATION
          end if
        end do
        end do

        !- Construct a moist-neutral profile.
        !  If no water-feasible and energy-conserving solution exists, expand the adjustment layer and retry.        

        lbase_try = lactive_base
        ltop_try  = lactive_top

        profile_converged     = .false.
        profile_status            = MCA_PROFILE_ENERGY_MAXITER

        do iter_expand=0, MCA_MAX_MASK_EXPAND
          zbase_try = zlev_diag(lbase_try)
          ztop_try  = zlev_diag(ltop_try)
          
          !- Construct a mask for the adjustment layer
          call make_dg_vertical_range_mask( adjustment_mask,     & ! (out)
            zlev_z, zbase_try, ztop_try, elem%Nnode_v, lmesh%NeZ ) ! (in)

          if ( debug_flag ) then
            write(*,*) "* Adjustment-profile attempt: iter_expand=", iter_expand
            write(*,*) "   lactive_base, lactive_top=", lactive_base, lactive_top
            write(*,*) "   lbase_try, ltop_try=", lbase_try, ltop_try
            write(*,*) "   zbase_try, ztop_try=", zbase_try, ztop_try
            write(*,*) "   saturation, mixing nodes=", count(forced_saturation_mask), count(adjustment_mask)
          end if

          !- Solve for a water-feasible and energy-conserving saturated profile in the adjustment layer

          call solve_moist_neutral_adjustment( &
            temp_adj, qvap_adj, rhoqvap_adj, rhoprecip_adj, profile_converged, profile_status,       & ! (out)
            qvap_z, pres_z, zlev_z, dens_z, temp_z,                                                  & ! (in)
            int_weight, adjustment_mask, forced_saturation_mask, elem%Nnode_v, lmesh%NeZ, debug_flag ) ! (in)

          if ( profile_converged ) then
            ! The region added by water-feasibility expansion also becomes part of the persistent active region.
            lactive_base = min(lactive_base, lbase_try)
            lactive_top  = max(lactive_top,  ltop_try)

            if ( debug_flag ) then
              write(*,*) "Moist-neutral profile converged: iter_expand=", iter_expand, "final: lbase, ltop=", lbase_try, ltop_try, ", zbase, ztop=", zbase_try, ztop_try

              ! write(*,*) " check:"
              ! call diagnose_convective_layer( lbase, ltop, is_unstable, zlev_diag,   & ! (out)
              !   temp_aj, qvap_aj, pres_z, zlev_z, lmesh, elem, nlev_diag, debug_flag, .true., lbase_try, ltop_try ) ! (in)
            end if
            exit            
          end if

          select case ( profile_status )
          case ( MCA_PROFILE_NO_ENERGY_ROOT )
            call expand_adjustment_layer( lbase_try, ltop_try, mask_expanded, nlev_diag )
            if ( .not. mask_expanded ) then
              exit
            end if
          case default
            exit
          end select

        end do ! end loop for iter_expand

        if ( .not. profile_converged ) then
          LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "iter_adj=", iter_adj, "ph=", ph, "ke_xy=", ke_xy, "iter_expand=", iter_expand
          select case ( profile_status )
          case ( MCA_PROFILE_NO_ENERGY_ROOT )
            LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "No energy-conserving root exists in the water-feasible temperature interval, even after adjustment-layer expansion."
            LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "original lbase,ltop=", lbase, ltop, "final lbase,ltop=", lbase_try, ltop_try
            do_adjustment = .false.
            exit
            ! call PRC_abort
          case ( MCA_PROFILE_NO_WATER_FEASIBLE_STATE )
            LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "No water-feasible trial moist-neutral profile exists."
            call PRC_abort
          case ( MCA_PROFILE_ADIABAT_FAILURE )
            LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "Moist-adiabat profile integration failed."
            call PRC_abort
          case ( MCA_PROFILE_ENERGY_MAXITER )
            LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "Energy root was bracketed, but bisection did not converge."
            call PRC_abort
          case ( MCA_PROFILE_NO_ADJUSTED_NODE )
            LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "No adjusted DG node was found."
            call PRC_abort
          case default
            LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "Unknown moist-neutral-profile construction error."
            call PRC_abort
          end select
        end if
      
        ! Accumulate the precipitation density from the current iteration
        rhoprecip_accum(:,:) = rhoprecip_accum(:,:) + rhoprecip_adj(:,:)

        ! Update the state for the next iteration
        temp_z(:,:) = temp_adj(:,:)
        qvap_z(:,:) = qvap_adj(:,:)
      
      end do ! End loop for iteration
  
      if ( do_adjustment .and. ( .not. adjustment_converged ) ) then
        LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "Moist convective adjustment did not converge: ph=", ph, "ke_xy=", ke_xy
        call PRC_abort
      end if
      
      
      !- Convert the adjusted state to DG tendencies

      if ( do_adjustment ) then
        do ke_z=1, lmesh%NeZ
          ke = ke_xy + (ke_z-1)*lmesh%Ne2D

          dens_z(:,ke_z) = dens_ini(:,ke_z) - rhoprecip_accum(:,ke_z)
          qvap_z(:,ke_z) = rhoqvap_adj(:,ke_z) / dens_z(:,ke_z)

          do pz=1, elem%Nnode_v
            qdry = 1.0_RP - qvap_z(pz,ke_z)
            rtot_ =  Rdry  * qdry +     Rvap * qvap_z(pz,ke_z)
            cptot_ = CPdry * qdry + CP_VAPOR * qvap_z(pz,ke_z)

            pres_z(pz,ke_z) = dens_z(pz,ke_z) * rtot_ * temp_z(pz,ke_z)
            pott_z(pz,ke_z) = temp_adj(pz,ke_z) * ( PRES0 / pres_z(pz,ke_z) )**( rtot_ / cptot_ )
          end do

          do pz=1, elem%Nnode_v
            p = ph + (pz-1)*elem%Nnode_h1D**2
            DENS_t (p,ke) = ( dens_z(pz,ke_z) - dens_ini(pz,ke_z) ) / dtsec
            RHOT_t (p,ke) = ( dens_z(pz,ke_z) * pott_z(pz,ke_z) - dens_ini(pz,ke_z) * pott_ini(pz,ke_z) ) / dtsec
            RHOQV_t(p,ke) = ( dens_z(pz,ke_z) * qvap_z(pz,ke_z) - rhoqvap_ini(pz,ke_z) ) / dtsec
          end do
        end do

        SFLX_RAIN(ph,ke_xy) = sum( rhoprecip_accum(:,:) * int_weight(:,:) ) / dtsec

      else
        do ke_z=1, lmesh%NeZ
          ke = ke_xy + (ke_z-1)*lmesh%Ne2D
          do pz=1, elem%Nnode_v
            p = ph + (pz-1)*elem%Nnode_h1D**2
            DENS_t (p,ke) = 0.0_RP
            RHOT_t (p,ke) = 0.0_RP
            RHOQV_t(p,ke) = 0.0_RP
          end do
        end do
        SFLX_RAIN(ph,ke_xy) = 0.0_RP
      end if

    end do ! end loop for ph
    end do ! end loop for ke_xy

    return
  end subroutine atm_phy_cp_dgm_mconv_adjustment_calc_tendency

  !> Finalize a module for moist convective adjustment scheme
  subroutine atm_phy_cp_dgm_mconv_adjustment_finalize()
    implicit none
    !----------------------------------------------------
    return
  end subroutine atm_phy_cp_dgm_mconv_adjustment_finalize


!- private subroutines --------------------------

!OCL SERIAL
  subroutine build_diagnostic_profile_from_dg( &
    hmse_diag, rh_diag, zlev_diag,             & ! (out)
    hmse_z, rh_z, zlev_z, lmesh, elem, nlev    ) ! (in)
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: nlev
    real(RP), intent(out) :: hmse_diag(nlev)
    real(RP), intent(out) :: rh_diag(nlev)
    real(RP), intent(out) :: zlev_diag(nlev)
    real(RP), intent(in) :: hmse_z(elem%Nnode_v,lmesh%NeZ)
    real(RP), intent(in) :: rh_z(elem%Nnode_v,lmesh%NeZ)
    real(RP), intent(in) :: zlev_z(elem%Nnode_v,lmesh%NeZ)

    integer :: pz, ke_z
    integer :: l
    integer :: Nnode_v
    !----------------------------------------------------

    Nnode_v = elem%Nnode_v
    l = 1

    hmse_diag(l) = hmse_z(1,1)
    rh_diag(l) = rh_z(1,1)
    zlev_diag(l) = zlev_z(1,1)
    do ke_z=1, lmesh%NeZ
      do pz=2, Nnode_v-1
        l = l + 1
        hmse_diag(l) = hmse_z(pz,ke_z)
        rh_diag(l) = rh_z(pz,ke_z)
        zlev_diag(l) = zlev_z(pz,ke_z)
      end do
      l = l + 1
      if ( ke_z < lmesh%NeZ ) then
        hmse_diag(l) = 0.5_RP * ( hmse_z(Nnode_v,ke_z) + hmse_z(1,ke_z+1) )
        rh_diag(l) = 0.5_RP * ( rh_z(Nnode_v,ke_z) + rh_z(1,ke_z+1) )
        zlev_diag(l) = 0.5_RP * ( zlev_z(Nnode_v,ke_z) + zlev_z(1,ke_z+1) )
      else
        hmse_diag(l) = hmse_z(Nnode_v,ke_z)
        rh_diag(l) = rh_z(Nnode_v,ke_z)
        zlev_diag(l) = zlev_z(Nnode_v,ke_z)
      end if
    end do
    return
  end subroutine build_diagnostic_profile_from_dg

!OCL SERIAL
  subroutine diagnose_convective_layer( lbase, ltop, is_unstable, zlev_diag, & ! (out)
    temp_z, qvap_z, pres_z, zlev_z, lmesh, elem, nlev,  debug_flag )           ! (in)
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: nlev
    integer, intent(out) :: lbase
    integer, intent(out) :: ltop
    logical, intent(out) :: is_unstable
    real(RP), intent(out) :: zlev_diag(nlev)
    real(RP), intent(in) :: temp_z(elem%Nnode_v,lmesh%NeZ)
    real(RP), intent(in) :: qvap_z(elem%Nnode_v,lmesh%NeZ)
    real(RP), intent(in) :: pres_z(elem%Nnode_v,lmesh%NeZ)
    real(RP), intent(in) :: zlev_z(elem%Nnode_v,lmesh%NeZ)
    logical, intent(in) :: debug_flag
    
    integer :: pz, ke_z
    integer :: l

    real(RP) :: hmse_sat_z(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: rh_z(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: qsat, cptot

    real(RP) :: rh(nlev)
    real(RP) :: hmse_sat(nlev)

    real(RP) :: dz
    real(RP) :: dhmse_sat, dhmse_sat_dz
    real(RP) :: rh_lyr
    logical :: unstable_pair
    logical :: is_inside

    integer :: Nnode_v
    !------------------------------------

    Nnode_v = elem%Nnode_v

    do ke_z=1, lmesh%NeZ
    do pz=1, Nnode_v

      call ATMOS_SATURATION_pres2qsat_liq( temp_z(pz,ke_z), pres_z(pz,ke_z), &
        qsat ) ! (out)

      rh_z(pz,ke_z) = qvap_z(pz,ke_z) / max(qsat, EPS)

      cptot = CPdry * (1.0_RP - qsat) + CP_VAPOR * qsat
      hmse_sat_z(pz,ke_z) = cptot * temp_z(pz,ke_z) + GRAV * zlev_z(pz,ke_z) + LHV0 * qsat
    end do
    end do

    call build_diagnostic_profile_from_dg( hmse_sat, rh, zlev_diag, &
      hmse_sat_z, rh_z, zlev_z, &
      lmesh, elem, nlev )

    !- Diagnose unstable layers based on the saturated MSE profile and relative humidity

    is_unstable = .false.
    is_inside   = .false.
    lbase = 0; ltop = 0

    do l=1, nlev-1
      dz = zlev_diag(l+1) - zlev_diag(l)
      dhmse_sat = hmse_sat(l+1) - hmse_sat(l)
      dhmse_sat_dz = dhmse_sat / dz
      rh_lyr = 0.5_RP * ( rh(l) + rh(l+1) )
      ! write(*,*) "l=", l, "zlev(l)=", zlev(l), "zlev(l+1)=", zlev(l+1), &
      !            "rh_lyr=", rh_lyr, "dhmse_sat_dz=", dhmse_sat_dz
      unstable_pair =                                           &
            ( min(rh(l), rh(l+1)) >= MCA_RH_TRIGGER )           &
        .and. ( dhmse_sat < -MCA_HMSE_DIFF_TOL )                &
        .and. ( dhmse_sat_dz < -MCA_HMSE_GRAD_TOL )


      if ( unstable_pair ) then
        if ( .not. is_inside ) then
          lbase = l
          is_inside = .true.
        end if
        ltop = l + 1
      else if ( is_inside ) then
       exit
      end if
    end do

    if ( is_inside ) then
      if ( zlev_diag(ltop) - zlev_diag(lbase) >= MCA_MIN_UNSTABLE_DEPTH ) is_unstable = .true.
    end if

    if ( debug_flag .and. is_unstable ) then
      write(*,*) "Diagnose convective layer:"
      write(*,*) "hmse_sat=", hmse_sat(lbase:ltop)
      write(*,*) "rh=", rh(lbase:ltop)
      write(*,*) "hmse_sat_dg=", hmse_sat_z(:,:)
      write(*,*) "rh_dg=", rh_z(:,:)
    end if

    return
  end subroutine diagnose_convective_layer

!OCL SERIAL
  subroutine make_dg_vertical_range_mask( adjustment_mask, & ! (out)
    zlev, zbase, ztop, npz, nez                            ) ! (in)
    implicit none
    integer, intent(in) :: npz
    integer, intent(in) :: nez
    logical, intent(out) :: adjustment_mask(npz,nez)
    real(RP), intent(in) :: zlev(npz,nez)
    real(RP), intent(in) :: zbase
    real(RP), intent(in) :: ztop

    integer :: ke_z, pz
    !----------------------------------------------------

    adjustment_mask(:,:) = .false.
    do ke_z = 1, nez
    do pz = 1, npz
      if ( zlev(pz,ke_z) >= zbase - MCA_Z_TOL .and. &
           zlev(pz,ke_z) <= ztop  + MCA_Z_TOL ) then
        adjustment_mask(pz,ke_z) = .true.
      end if
    end do
    end do    
    return
  end subroutine make_dg_vertical_range_mask

  !> Expand the adjustment layer by one diagnostic level above and below, if possible
  !! Numerical fallback: enlarge the adjustment region by one diagnostic level on each available side. This is not a physically based entrainment closure.  
!OCL SERIAL
  subroutine expand_adjustment_layer( &
    lbase, ltop, expanded,            & ! (inout,out)
    nlev       ) ! (in)
    implicit none

    integer, intent(inout) :: lbase
    integer, intent(inout) :: ltop
    logical, intent(out) :: expanded
    integer, intent(in) :: nlev

    logical :: can_expand_below
    logical :: can_expand_above
    !------------------------------------------------------------

    can_expand_below = lbase > 1
    can_expand_above = ltop  < nlev

    expanded = .false.

    ! Symmetrically expand by one diagnostic level whenever possible.
    if ( can_expand_below ) then
      lbase = lbase - 1
      expanded = .true.
    end if

    if ( can_expand_above ) then
      ltop = ltop + 1
      expanded = .true.
    end if
    return
  end subroutine expand_adjustment_layer

  !> Solve for a trial moist-neutral profile satisfying: 1. available-water constraint, 2. column-integrated moist-static-energy conservation.
  !!
  !! The pressure and density profiles are held fixed during the solve.
  !! Nodes selected by forced_saturation_mask are set to saturation.
  !! At other adjusted nodes, the original vapor mixing ratio is retained unless it exceeds saturation.  
  !!
!OCL SERIAL
  subroutine solve_moist_neutral_adjustment( &
    temp_adj, qvap_adj, rhoqvap_aj, rhoprecip_aj, converged, profile_status, & ! (out)
    qvap, pres, zlev, dens, temp,                                            & ! (in)
    int_weight, mix_mask, saturation_mask, npz, nez, debug_flag              ) ! (in)

    implicit none
    integer, intent(in) :: npz
    integer, intent(in) :: nez
    real(RP), intent(out) :: temp_adj(npz,nez)
    real(RP), intent(out) :: qvap_adj(npz,nez)
    real(RP), intent(out) :: rhoqvap_aj(npz,nez)
    real(RP), intent(out) :: rhoprecip_aj(npz,nez)
    logical, intent(out) :: converged
    integer, intent(out) :: profile_status
    real(RP), intent(in) :: qvap(npz,nez)
    real(RP), intent(in) :: pres(npz,nez)
    real(RP), intent(in) :: zlev(npz,nez)
    real(RP), intent(in) :: dens(npz,nez)
    real(RP), intent(in) :: temp(npz,nez)
    real(RP), intent(in) :: int_weight(npz,nez)
    logical, intent(in) :: mix_mask(npz,nez)
    logical, intent(in) :: saturation_mask(npz,nez)
    logical, intent(in) :: debug_flag

    ! Initial column-integrated quantities
    real(RP) :: energy_target
    real(RP) :: water_mass_ini

    ! Trial-profile quantities
    real(RP) :: energy_trial
    real(RP) :: water_trial
    real(RP) :: residual

    ! Energy root bracket
    real(RP) :: tbase_lo, tbase_mid, tbase_hi
    real(RP) :: residual_lo, residual_hi
    real(RP) :: energy_lo, energy_hi
    real(RP) :: water_lo, water_hi

    ! Water-feasible upper-bound search
    real(RP) :: twater_lo, twater_mid, twater_hi
    real(RP) :: tbase_water_max
    real(RP) :: water_mid
    real(RP) :: energy_dummy

    ! Work arrays
    real(RP) :: temp_trial(npz,nez)
    real(RP) :: qvap_trial(npz,nez)

    ! Precipitation diagnostics
    real(RP) :: water_mass_adj
    real(RP) :: precip_mass_local
    real(RP) :: precip_mass_total
    real(RP) :: positive_cond_mass
    real(RP) :: precip_scale

    ! Tolerances
    real(RP) :: energy_tol
    real(RP) :: boundary_energy_tol    
    real(RP) :: water_tol

    integer :: ke_z, pz
    integer :: iter
    integer :: iter_water

    logical :: found_base
    logical :: adiabat_profile_ok

    logical :: water_lo_feasible
    logical :: water_hi_feasible

    real(RP), parameter :: MCA_TBASE_RANGE = 40.0_RP
    !------------------------------------------------------------

    converged = .false.
    profile_status = MCA_PROFILE_ENERGY_MAXITER

    ! Initialize the adjusted profile with the original state

    temp_adj(:,:) = temp(:,:)
    qvap_adj(:,:) = qvap(:,:)
    rhoqvap_aj(:,:) = dens(:,:) * qvap(:,:)
    rhoprecip_aj(:,:) = 0.0_RP

    !- Calculate initial column-integrated MSE

    call integ_masked_column_mse( energy_target, & ! (out)
      temp, qvap, zlev, dens,           & ! (in)
      int_weight, mix_mask, npz, nez    ) ! (in)

    energy_tol = MCA_COLUMN_ENERGY_ATOL + MCA_COLUMN_ENERGY_RTOL * abs(energy_target)
    boundary_energy_tol = MCA_COLUMN_ENERGY_ATOL + MCA_BOUNDARY_ENERGY_RTOL * abs(energy_target)

    !- Find temperature at the lowest adjusted node

    tbase_mid = 0.0_RP
    found_base = .false.

    do ke_z = 1, nez
      do pz = 1, npz
        if ( mix_mask(pz,ke_z) ) then
          tbase_mid = temp(pz,ke_z)
          found_base = .true.
          exit
        end if
      end do
      if ( found_base ) exit
    end do

    if ( .not. found_base ) then
      profile_status = MCA_PROFILE_NO_ADJUSTED_NODE
      if ( debug_flag ) then
        write(*,*) "construct_moist_neutral_profile: No adjusted DG node was found."
      end if
      return
    end if

    !-  Calculate the total water mass in the column before adjustment

    water_mass_ini = 0.0_RP

    do ke_z = 1, nez
    do pz = 1, npz
      if ( mix_mask(pz,ke_z) ) then
        water_mass_ini = water_mass_ini + int_weight(pz,ke_z) * dens(pz,ke_z) * qvap(pz,ke_z)
      end if
    end do
    end do

    water_tol = MCA_WATER_ATOL + MCA_WATER_RTOL * max(abs(water_mass_ini),1.0_RP)    

    !- Step 1:
    ! Determine the upper limit of the water-feasible base temperature.
    ! A trial moist-neutral profile is water-feasible when its vapor mass does not exceed the initial vapor mass in the adjustment region.

    twater_lo = tbase_mid - MCA_TBASE_RANGE
    twater_hi = tbase_mid + MCA_TBASE_RANGE

    ! Evaluate the water mass at the lower bound of the temperature range

    call evaluate_trial_moist_neutral_profile( temp_trial, qvap_trial, water_lo, energy_dummy, adiabat_profile_ok, & ! (out)
      twater_lo, temp, qvap, dens, pres, zlev, int_weight, mix_mask, saturation_mask, nez, npz, debug_flag         ) ! (in)

    if ( adiabat_profile_ok ) then
      water_lo_feasible = ( water_lo <= water_mass_ini + water_tol )
    else
      profile_status = MCA_PROFILE_ADIABAT_FAILURE
      if ( debug_flag ) then
        write(*,*) "construct_moist_neutral_profile: Moist-adiabat integration failed. tbase_lo=", twater_lo
      end if
      return
    end if

    if ( .not. water_lo_feasible ) then
      profile_status = MCA_PROFILE_NO_WATER_FEASIBLE_STATE

      if ( debug_flag ) then
        write(*,*) "No water-feasible saturated profile exists."
        write(*,*) "  tbase_lo=", twater_lo, "water_mass_lo=", water_lo, "water_mass_ini =", water_mass_ini
      end if
      return
    end if  
    
    ! Evaluate the water mass at the upper bound of the temperature range

    call evaluate_trial_moist_neutral_profile( temp_trial, qvap_trial, water_hi, energy_dummy, adiabat_profile_ok, & ! (out)
      twater_hi, temp, qvap, dens, pres, zlev, int_weight, mix_mask, saturation_mask, nez, npz, debug_flag         ) ! (in)

    if ( adiabat_profile_ok ) then
      water_hi_feasible = ( water_hi <= water_mass_ini + water_tol )
    else
      water_hi_feasible = .false.
      water_hi = huge(1.0_RP)

      if ( debug_flag ) then
        write(*,*) "High-temperature trial profile is not admissible."
        write(*,*) "  Treating it as the infeasible upper bracket: twater_hi=", twater_hi 
      end if
    end if

    !-

    if ( debug_flag ) then
      write(*,*) "Water feasibility at initial bounds:"
      write(*,*) "  t_lo, water_lo  =", twater_lo, water_lo
      write(*,*) "  t_hi, water_hi  =", twater_hi, water_hi
      write(*,*) "  lower feasible  =", water_lo_feasible
      write(*,*) "  upper feasible  =", water_hi_feasible      
      write(*,*) "  water_mass_ini  =", water_mass_ini
      write(*,*) "  water_tol       =", water_tol
    end if


    if ( water_hi_feasible ) then
      ! If the high side is still feasible, the chosen initial temperature range does not yet contain the water or thermodynamic upper boundary.
      tbase_water_max = twater_hi
    else    

      ! Search for the maximum base temperature for which the moist-adiabat construction succeeds 
      ! and the resulting vapor mass remains feasible.      
      ! The lower endpoint is always:
      !   - thermodynamically admissible
      !   - water feasible
      !
      ! The upper endpoint is either:
      !   - water infeasible, or
      !   - thermodynamically inadmissible

      do iter_water = 1, MCA_MAX_WATER_ITER
        twater_mid = 0.5_RP * (twater_lo + twater_hi)

        call evaluate_trial_moist_neutral_profile( temp_trial, qvap_trial, water_mid, energy_dummy, adiabat_profile_ok, & ! (out)
          twater_mid, temp, qvap, dens, pres, zlev, int_weight, mix_mask, saturation_mask, nez, npz, debug_flag         ) ! (in)

        if ( .not. adiabat_profile_ok ) then
          ! Thermodynamically inadmissible midpoint:
          ! move the upper boundary downward.
          twater_hi = twater_mid
          water_hi  = huge(1.0_RP)

          if ( debug_flag ) then
            write(*,*) "water iter=", iter_water
            write(*,*) "  inadmissible midpoint=", twater_mid, "  new bracket=", twater_lo, twater_hi
          end if

          if ( abs(twater_hi - twater_lo) <= MCA_TEMP_TOL ) exit

          cycle
        end if

        if ( water_mid <= water_mass_ini + water_tol ) then
          ! Admissible and water feasible.
          twater_lo = twater_mid
          water_lo  = water_mid
        else
          ! Admissible but water infeasible.
          twater_hi = twater_mid
          water_hi  = water_mid
        end if

        if ( debug_flag ) then
          write(*,*) "water iter=", iter_water, ", tbase=", twater_lo, twater_mid, twater_hi, ", water=", water_lo, water_mid, water_hi, ", target=", water_mass_ini
        end if

        if ( abs(twater_hi - twater_lo) <= MCA_TEMP_TOL ) exit
        if ( abs(water_mid - water_mass_ini) <= water_tol ) exit
      end do

      ! The lower endpoint is maintained on the valid and water-feasible side.      
      tbase_water_max = twater_lo      
    end if


    !- Step 2:
    ! Evaluate the energy residual at both ends of the water-feasible interval.

    tbase_lo = tbase_mid - MCA_TBASE_RANGE
    tbase_hi = tbase_water_max

    !-
    call evaluate_trial_moist_neutral_profile( temp_trial, qvap_trial, water_lo, energy_lo, adiabat_profile_ok, & ! (out)
      tbase_lo, temp, qvap, dens, pres, zlev, int_weight, mix_mask, saturation_mask, nez, npz, debug_flag       ) ! (in)

    if ( .not. adiabat_profile_ok ) then
      profile_status = MCA_PROFILE_ADIABAT_FAILURE
      return
    end if

    residual_lo = energy_lo - energy_target

    !-
    call evaluate_trial_moist_neutral_profile( temp_trial, qvap_trial, water_hi, energy_hi, adiabat_profile_ok, & ! (out)
      tbase_hi, temp, qvap, dens, pres, zlev, int_weight, mix_mask, saturation_mask, nez, npz, debug_flag       ) ! (in)

    if ( .not. adiabat_profile_ok ) then
      profile_status = MCA_PROFILE_ADIABAT_FAILURE

      if ( debug_flag ) then
        write(*,*) "Internal inconsistency:"
        write(*,*) "  tbase_water_max should be admissible but evaluation failed."
        write(*,*) "  tbase_hi=", tbase_hi
      end if
      return
    end if

    residual_hi = energy_hi - energy_target

    if ( debug_flag ) then
      write(*,*) "Energy-root feasibility check:"
      write(*,*) "  energy_target       =", energy_target
      write(*,*) "  energy_tol          =", energy_tol
      write(*,*) "  tbase_lo            =", tbase_lo
      write(*,*) "  energy_lo           =", energy_lo
      write(*,*) "  residual_lo         =", residual_lo
      write(*,*) "  tbase_water_max     =", tbase_hi
      write(*,*) "  energy_water_max    =", energy_hi
      write(*,*) "  residual_water_max  =", residual_hi
      write(*,*) "  water_at_upper      =", water_hi
      write(*,*) "  water_mass_ini      =", water_mass_ini
    end if    

    ! Check whether either endpoint is already an energy root.
    if ( abs(residual_lo) <= boundary_energy_tol ) then

      call evaluate_trial_moist_neutral_profile( temp_adj, qvap_adj, water_mass_adj, energy_trial, adiabat_profile_ok, & ! (out)
        tbase_lo, temp, qvap, dens, pres, zlev, int_weight, mix_mask, saturation_mask, nez, npz, debug_flag            ) ! (in)

      if ( .not. adiabat_profile_ok ) then
        profile_status = MCA_PROFILE_ADIABAT_FAILURE
        return
      end if
      converged = .true.
      profile_status = MCA_PROFILE_SUCCESS

    else if ( abs(residual_hi) <= boundary_energy_tol ) then

      call evaluate_trial_moist_neutral_profile( temp_adj, qvap_adj, water_mass_adj, energy_trial, adiabat_profile_ok, & ! (out)
        tbase_hi, temp, qvap, dens, pres, zlev, int_weight, mix_mask, saturation_mask, nez, npz, debug_flag            ) ! (in)

      if ( .not. adiabat_profile_ok ) then
        profile_status = MCA_PROFILE_ADIABAT_FAILURE
        return
      end if

      converged = .true.
      profile_status = MCA_PROFILE_SUCCESS

    else if ( residual_lo * residual_hi > 0.0_RP ) then
      profile_status = MCA_PROFILE_NO_ENERGY_ROOT
      if ( debug_flag ) then
        write(*,*) "No energy root exists in the water-feasible interval."
        write(*,*) "  tbase_lo           =", tbase_lo
        write(*,*) "  tbase_water_max    =", tbase_hi
        write(*,*) "  residual_lo        =", residual_lo
        write(*,*) "  residual_upper     =", residual_hi
        write(*,*) "  water_mass_ini     =", water_mass_ini
        write(*,*) "  water_mass_upper   =", water_hi
      end if
      return
    end if

    !- Step 3:
    ! If the endpoint check found a valid bracket, solve the column-energy constraint by bisection.
    
    if ( .not. converged ) then
      do iter=1, MCA_MAX_ENERGY_ITER

        tbase_mid = 0.5_RP * ( tbase_lo + tbase_hi )

        call evaluate_trial_moist_neutral_profile( temp_adj, qvap_adj, water_trial, energy_trial, adiabat_profile_ok, & ! (out)
          tbase_mid, temp, qvap, dens, pres, zlev, int_weight, mix_mask, saturation_mask, nez, npz, debug_flag        ) ! (in)

        if ( .not. adiabat_profile_ok ) then
          profile_status = MCA_PROFILE_ADIABAT_FAILURE

          if ( debug_flag ) then
            write(*,*) "Moist-adiabat integration failed during"
            write(*,*) "energy-root search: iteration=", iter, ", tbase=", tbase_mid
          end if
          return
        end if

        ! The energy search must remain inside the previously determined water-feasible interval.

        if ( water_trial > water_mass_ini + water_tol ) then
          profile_status = MCA_PROFILE_NO_ENERGY_ROOT
          if ( debug_flag ) then
            write(*,*) "Energy bisection entered a water-infeasible state."
            write(*,*) "  tbase          =", tbase_mid
            write(*,*) "  water_trial    =", water_trial
            write(*,*) "  water_mass_ini =", water_mass_ini
          end if
          return
        end if

        residual = energy_trial - energy_target

        if ( debug_flag ) then
          write(*,*) "energy iter=", iter
          write(*,*) "  tbase=", tbase_lo, tbase_mid, tbase_hi
          write(*,*) "  energy_trial=", energy_trial, ", energy_target=", energy_target
          write(*,*) "  residual=", residual
          write(*,*) "  water_trial=", water_trial
        end if

        if ( abs(residual) <= energy_tol ) then
          water_mass_adj = water_trial
          converged = .true.
          profile_status = MCA_PROFILE_SUCCESS
          exit
        end if

        ! General bisection update based on the residual signs.

        if ( residual_lo * residual <= 0.0_RP ) then
          tbase_hi   = tbase_mid
          residual_hi = residual
        else
          tbase_lo   = tbase_mid
          residual_lo = residual
        end if
      
      end do ! end loop for iteration
    end if

    if ( .not. converged ) then
      profile_status = MCA_PROFILE_ENERGY_MAXITER

      if ( debug_flag ) then
        write(*,*) "Energy bisection did not converge."
        write(*,*) "  tbase_lo=", tbase_lo, ", tbase_hi=", tbase_hi
        write(*,*) "  residual_lo=", residual_lo, ", residual_hi=", residual_hi
        write(*,*) "  energy_tol=", energy_tol
      end if
      return
    end if

    ! Step 4:
    ! Convert the converged vapor mixing ratio to density form and diagnose precipitation.

    rhoqvap_aj(:,:) = dens(:,:) * qvap_adj(:,:)

    precip_mass_total = max( 0.0_RP, water_mass_ini - water_mass_adj )
    positive_cond_mass = 0.0_RP
    do ke_z=1, nez
    do pz=1, npz
      if ( mix_mask(pz,ke_z) ) then
        precip_mass_local = max( 0.0_RP, dens(pz,ke_z) * qvap(pz,ke_z) - rhoqvap_aj(pz,ke_z) )

        positive_cond_mass = positive_cond_mass + precip_mass_local * int_weight(pz,ke_z)
      end if
    end do
    end do

    ! Distribute the column-integrated removed water mass in proportion to positive local vapor-density reductions.

    if ( positive_cond_mass > EPS ) then
      precip_scale = precip_mass_total / positive_cond_mass
      do ke_z=1, nez
      do pz=1, npz
        if ( mix_mask(pz,ke_z) ) then
          precip_mass_local = max( 0.0_RP, dens(pz,ke_z) * qvap(pz,ke_z) - rhoqvap_aj(pz,ke_z) )
          rhoprecip_aj(pz,ke_z) = precip_scale * precip_mass_local
        end if
      end do
      end do
    end if

    return
  end subroutine solve_moist_neutral_adjustment

!OCL SERIAL
  subroutine evaluate_trial_moist_neutral_profile( &
    temp_work, qvap_work, water_mass, energy, profile_ok,  & ! (out)
    tbase, temp, qvap, dens, pres, zlev, int_weight,       & ! (in)
    adj_mask, forced_saturation_mask, nez, npz, debug_flag ) ! (in)
    implicit none

    integer, intent(in) :: nez
    integer, intent(in) :: npz
    real(RP), intent(out) :: temp_work(npz,nez)
    real(RP), intent(out) :: qvap_work(npz,nez)
    real(RP), intent(out) :: water_mass
    real(RP), intent(out) :: energy
    logical,  intent(out) :: profile_ok
    real(RP), intent(in)  :: tbase
    real(RP), intent(in)  :: temp(npz,nez)
    real(RP), intent(in)  :: qvap(npz,nez)
    real(RP), intent(in)  :: pres(npz,nez)
    real(RP), intent(in)  :: zlev(npz,nez)
    real(RP), intent(in)  :: dens(npz,nez)
    real(RP), intent(in)  :: int_weight(npz,nez)
    logical,  intent(in)  :: adj_mask(npz,nez)
    logical,  intent(in)  :: forced_saturation_mask(npz,nez)
    logical,  intent(in)  :: debug_flag

    real(RP) :: qsat_work(npz,nez)

    integer :: ke_z, pz
    !------------------------------------------------

    temp_work(:,:) = temp(:,:)
    qvap_work(:,:) = qvap(:,:)
    qsat_work(:,:) = 0.0_RP

    !- Evaluate the reference temperature profile over mix_mask.

    call build_moist_adiabat_temperature( &
      temp_work, qsat_work, profile_ok,                 & ! (inout,out)
      tbase, pres, zlev, adj_mask, npz, nez, debug_flag ) ! (in)

    if ( .not. profile_ok ) then
      water_mass = huge(1.0_RP)
      energy     = huge(1.0_RP)
      return
    end if

    !- Humidity construction:
    ! 1. Nearly saturated nodes in the diagnosed unstable core:
    !      qv = qsat.
    ! 2. Other nodes inside the adjustment region:
    !      retain the input qv unless the trial temperature produces supersaturation.
    ! 3. Outside the adjustment region:
    !      retain the input temperature and qv.

    do ke_z = 1, nez
    do pz = 1, npz
      if ( adj_mask(pz,ke_z) ) then
        if ( forced_saturation_mask(pz,ke_z) ) then
          qvap_work(pz,ke_z) = qsat_work(pz,ke_z)
        else
          qvap_work(pz,ke_z) = min(qvap(pz,ke_z), qsat_work(pz,ke_z))
        end if
      end if
      qvap_work(pz,ke_z) = max(0.0_RP, qvap_work(pz,ke_z))      
    end do
    end do    

    !- Column water over the full mixing region.

    water_mass = 0.0_RP
    do ke_z = 1, nez
    do pz = 1, npz
      if ( adj_mask(pz,ke_z) ) then
        water_mass = water_mass + int_weight(pz,ke_z) * dens(pz,ke_z) * qvap_work(pz,ke_z)
      end if
    end do
    end do

    !- Column MSE over the full mixing region.

    call integ_masked_column_mse( energy, & ! (out)
      temp_work, qvap_work, zlev, dens,   & ! (in)
      int_weight, adj_mask, npz, nez      ) ! (in)

    return
  end subroutine evaluate_trial_moist_neutral_profile

  ! The column-MSE constraint is imposed while holding the density and pressure profiles fixed. 
  ! Thermodynamic composition includes dry air and water vapor only. 
  ! Condensate sensible heat and the energy carried away by precipitation are neglected.  
!OCL SERIAL
  subroutine integ_masked_column_mse( energy, & ! (out)
    temp, qv, zlev, dens,                     & ! (in)
    int_weight, integ_mask, npz, nez          ) ! (in)
    implicit none
    integer, intent(in) :: npz
    integer, intent(in) :: nez
    real(RP), intent(out) :: energy
    real(RP), intent(in) :: temp(npz,nez)
    real(RP), intent(in) :: qv(npz,nez)
    real(RP), intent(in) :: zlev(npz,nez)
    real(RP), intent(in) :: dens(npz,nez)
    real(RP), intent(in) :: int_weight(npz,nez)
    logical, intent(in) :: integ_mask(npz,nez)

    integer :: ke_z, pz
    real(RP) :: qdry
    real(RP) :: cptot
    !------------------------------------------------

    energy = 0.0_RP

    do ke_z = 1, nez
    do pz = 1, npz
      if ( integ_mask(pz,ke_z) ) then
        qdry = 1.0_RP - qv(pz,ke_z) ! - qcon_(pz_,ke_z_)

        cptot = CPdry    * qdry         &
              + CP_VAPOR * qv(pz,ke_z)
                                        ! Liquid-water sensible heat is neglected.
        
        energy = energy + int_weight(pz,ke_z) * dens(pz,ke_z) * &
          ( cptot * temp(pz,ke_z) + GRAV * zlev(pz,ke_z)        &
          + LHV0 * qv(pz,ke_z)                                  )
      end if
    end do
    end do
    return
  end subroutine integ_masked_column_mse

  !----
!OCL SERIAL
  subroutine build_moist_adiabat_temperature( &
    temp_profile, qsat_profile, profile_converged, & ! (out)
    tbase, pres, zlev, mask, npz, nez, debug_flag  ) ! (in)
    implicit none
    integer, intent(in) :: npz
    integer, intent(in) :: nez
    real(RP), intent(inout) :: temp_profile(npz,nez)
    real(RP), intent(inout) :: qsat_profile(npz,nez)
    logical, intent(out) :: profile_converged
    real(RP), intent(in) :: tbase
    real(RP), intent(in) :: pres(npz,nez)
    real(RP), intent(in) :: zlev(npz,nez)
    logical, intent(in) :: mask(npz,nez)
    logical, intent(in) :: debug_flag

    integer :: ke_z, pz

    real(RP) :: t_prev, p_prev, z_prev
    real(RP) :: t_now, p_now, z_now

    logical :: started
    logical :: is_converge_local
    !----------------------------------------------------

    started = .false.
    profile_converged = .true.

    do ke_z = 1, nez
    do pz = 1, npz
      if ( .not. mask(pz,ke_z) ) cycle

      z_now = zlev(pz,ke_z)
      p_now = pres(pz,ke_z)

      if ( .not. started ) then
        t_now = tbase
        started = .true.
      ! else if ( abs(z_now-z_prev) <= MCA_Z_TOL ) then
      !   t_now = t_prev
      else
        call solve_next_moist_adiabat_temp( t_now, is_converge_local, & ! (out)
          t_prev, p_prev, z_prev, p_now, z_now                        ) ! (in)
        
        if ( .not. is_converge_local ) then
          profile_converged = .false.
          return
        end if
      end if

      temp_profile(pz,ke_z) = t_now
      call ATMOS_SATURATION_pres2qsat_liq( t_now, p_now, & ! (in)
        qsat_profile(pz,ke_z)                            ) ! (out)

      t_prev = t_now 
      p_prev = p_now
      z_prev = z_now
    end do
    end do

    return
  end subroutine build_moist_adiabat_temperature

  !> Solve the temperature at (pres1,zlev1) such that the saturated MSE equals that at (temp0,pres0,zlev0).  
!OCL SERIAL
  subroutine solve_next_moist_adiabat_temp( temp1, converged, & ! (out)
    temp0, pres0, zlev0, pres1, zlev1                         ) ! (in)
    implicit none
    real(RP), intent(out) :: temp1
    logical,  intent(out) :: converged
    real(RP), intent(in)  :: temp0
    real(RP), intent(in)  :: pres0
    real(RP), intent(in)  :: zlev0
    real(RP), intent(in)  :: pres1
    real(RP), intent(in)  :: zlev1

    real(RP) :: temp_lo, temp_mid, temp_hi
    real(RP) :: f_lo, f_mid, f_hi
    real(RP) :: hs_target, hs_mid
    real(RP) :: temp_tol,  energy_tol

    integer :: iter

    integer, parameter :: MAX_ITER = 60
    real(RP), parameter :: TEMP_RANGE = 40.0_RP
    !----------------------------------------------------

    converged = .false.

    call saturated_mse_point( temp0, pres0, zlev0, &
      hs_target )

    temp_lo = temp0 - TEMP_RANGE
    temp_hi = temp0 + TEMP_RANGE

    call saturated_mse_point( temp_lo, pres1, zlev1, &
      hs_mid )
    f_lo = hs_mid - hs_target

    call saturated_mse_point( temp_hi, pres1, zlev1, &
      hs_mid )
    f_hi = hs_mid - hs_target

    if ( f_lo * f_hi > 0.0_RP ) then
      return
    end if

    energy_tol = MCA_LOCAL_HMSE_ATOL + MCA_LOCAL_HMSE_RTOL * abs(hs_target)
    temp_tol   = MCA_TEMP_TOL

    do iter = 1, MAX_ITER
      temp_mid = 0.5_RP * ( temp_lo + temp_hi )

      call saturated_mse_point( temp_mid, pres1, zlev1, &
        hs_mid )

      f_mid = hs_mid - hs_target

      if ( f_lo * f_mid <= 0.0_RP ) then
        temp_hi = temp_mid
        f_hi    = f_mid
      else
        temp_lo = temp_mid
        f_lo    = f_mid
      end if

      if (     ( abs(temp_hi - temp_lo) <= temp_tol ) &
          .or. ( abs(f_mid) <= energy_tol )           ) then
        temp1 = temp_mid
        converged = .true.
        return
      end if
      
    end do

    temp1 = temp_mid
    return
  contains
    subroutine saturated_mse_point( &
      temp, pres, zlev, hmse_sat )
      implicit none
      real(RP), intent(in)  :: temp
      real(RP), intent(in)  :: pres
      real(RP), intent(in)  :: zlev
      real(RP), intent(out) :: hmse_sat

      real(RP) :: qsat
      real(RP) :: cptot
      !------------------------------------------------
      call ATMOS_SATURATION_pres2qsat_liq( &
        temp, pres, qsat )

      cptot = CPdry * (1.0_RP-qsat) + CP_VAPOR * qsat

      hmse_sat = cptot * temp + GRAV * zlev  + LHV0 * qsat
      return
    end subroutine saturated_mse_point
  end subroutine solve_next_moist_adiabat_temp

end module scale_atm_phy_cp_dgm_mconv_adjustment

