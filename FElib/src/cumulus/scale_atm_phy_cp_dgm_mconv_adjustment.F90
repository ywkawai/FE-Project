!> module FElib / Atmosphere / Physics cumulus parameterization
!!
!! @par Description
!!      A module to provide a moist convective adjustment scheme for cumulus parameterization in atmospheric model
!!  
!! @author Yuta Kawai, Team SCALE
!!
!! @par Reference
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
  !
  integer, parameter :: MCA_MAX_ADJUST_ITER = 10
  integer, parameter :: MCA_MAX_ENERGY_ITER = 40
  integer, parameter :: MCA_MAX_LOCAL_ITER  = 20  

  real(RP), parameter :: MCA_RH_TRIGGER  = 1.00_RP
  real(RP), parameter :: MCA_ENERGY_RTOL = 1.0E-10_RP
  real(RP), parameter :: MCA_TEMP_TOL    = 1.0E-8_RP

  ! Parameters to avoid detection of numerical noise
  real(RP), parameter :: MCA_HMSE_GRAD_TOL      = 0.0_RP  
  real(RP), parameter :: MCA_MIN_UNSTABLE_DEPTH = 0.0_RP
  real(RP), parameter :: MCA_Z_TOL = 100.0_RP * epsilon(1.0_RP)

contains
  subroutine atm_phy_cp_dgm_mconv_adjustment_setup()
    implicit none
    !----------------------------------------------------
    return
  end subroutine atm_phy_cp_dgm_mconv_adjustment_setup

!OCL SERIAL
  subroutine atm_phy_cp_dgm_mconv_adjustment_calc_tendency( &
    DENS_t, RHOT_t, RHOQV_t, SFLX_RAIN, &
    DDENS, DRHOT, QV, PT, PRES, DENS_hyd, Rtot, CPtot,                          &
    dtsec, lmesh, elem, elem1D                                                  )    
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
    real(RP) :: rtot_(elem%Nnode_v)
    real(RP) :: cptot_(elem%Nnode_v)

    real(RP) :: elem_width_z
    real(RP) :: int_weight(elem%Nnode_v,lmesh%NeZ)

    real(RP) :: pott_ini(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: qvap_ini(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: dens_ini(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: rhoqvap_ini(elem%Nnode_v,lmesh%NeZ)

    integer :: lbase, ltop
    logical :: is_unstable
    logical :: adjust_mask(elem%Nnode_v,lmesh%NeZ)

    integer :: iter_adj
    real(RP) :: temp_aj(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: qdry_aj(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: qvap_aj(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: rhoqvap_aj(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: rhoprecip_aj(elem%Nnode_v,lmesh%NeZ)    ! Condensate water density [kg/m^3] with each iteration, which is removed by precipitation
    real(RP) :: rhoprecip_accum(elem%Nnode_v,lmesh%NeZ)
    real(RP) :: rhoqdry_aj(elem%Nnode_v,lmesh%NeZ)

    integer :: l
    real(RP), allocatable :: pott_diag(:)
    real(RP), allocatable :: qvap_diag(:)
    real(RP), allocatable :: pres_diag(:)
    real(RP), allocatable :: temp_diag(:)
    real(RP), allocatable :: zlev_diag(:)
    real(RP) :: qdry_diag, Rtot_diag, CPtot_diag 
    integer :: nlev_diag

    logical :: is_profile_converged
    logical :: do_adjustment
    logical :: adjustment_converged

    integer :: Hslice_b(elem%Nnode_h1D**2), Hslice_t(elem%Nnode_h1D**2)

    real(RP) :: Rvap

    integer, parameter :: MAX_ITER_ADJ = 10
    !----------------------------------------------------

    nlev_diag = lmesh%NeZ * (elem%Nnode_v-1) + 1
    allocate( pott_diag(nlev_diag), qvap_diag(nlev_diag), pres_diag(nlev_diag), &
              temp_diag(nlev_diag), zlev_diag(nlev_diag) )
    
    Hslice_b(:) = elem%Hslice(:,1)
    Hslice_t(:) = elem%Hslice(:,elem%Nnode_v)

    Rvap = CP_VAPOR - CV_VAPOR

    !$omp parallel do collapse(2) private(ke,p, &
    !$omp dens_z, pres_z, temp_z, zlev_z, pott_z, qvap_z, rtot_,cptot_,                            &
    !$omp dens_ini, pott_ini, qvap_ini, rhoqvap_ini,                                               &
    !$omp pott_diag, qvap_diag, pres_diag, temp_diag, zlev_diag, qdry_diag, Rtot_diag, CPtot_diag, &
    !$omp lbase, ltop, is_unstable, adjust_mask,                                                   &
    !$omp temp_aj, qdry_aj, qvap_aj, rhoqdry_aj, rhoqvap_aj, rhoprecip_aj, rhoprecip_accum,        &
    !$omp is_profile_converged, do_adjustment, adjustment_converged,                               &
    !$omp int_weight,elem_width_z )
    do ke_xy=1, lmesh%Ne2D
    do ph=1, elem%Nnode_h1D**2
    
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

        rhoprecip_accum(:,ke_z) = 0.0_RP
      end do

      !** Loop for iteration of the moist convective adjustment

      adjustment_converged = .false.
      do_adjustment        = .false.
      do iter_adj=1, MAX_ITER_ADJ

        !- Build a vertical profile which is continuous at vertical element interfaces from DG solution

        call build_diag_zprofile_from_dg( &
          pott_diag, qvap_diag, pres_diag, zlev_diag, & ! (out)
          pott_z, qvap_z, pres_z, zlev_z,             & ! (in)
          lmesh, elem, nlev_diag                      ) ! (in)
        
        do l=1, nlev_diag
          qdry_diag = 1.0_RP - qvap_diag(l)
          CPtot_diag = CPdry * qdry_diag + CP_VAPOR * qvap_diag(l)
          Rtot_diag  =  Rdry * qdry_diag +     Rvap * qvap_diag(l)
          temp_diag(l) = pott_diag(l) * ( pres_diag(l) / PRES0 )**( Rtot_diag / CPtot_diag )
        end do

        !- Diagnose unstable layers 

        call diagnose_convective_layer( lbase, ltop, is_unstable, & ! (out)
          temp_diag, pott_diag, qvap_diag, pres_diag, zlev_diag,  & ! (in)
          nlev_diag )                                               ! (in)

        ! write(*,*) "-- ph=", ph, "ke_xy=", ke_xy, "iter_adj=", iter_adj
        ! write(*,*) "   lbase, ltop, is_unstable = ", lbase, ltop, is_unstable

        if ( .not. is_unstable ) then
          adjustment_converged = .true.
          exit
        end if

        do_adjustment = .true.

        !- Generate a mask for DG vertical nodes where adjustment is needed

        call make_dg_adjustment_mask( adjust_mask,                           & ! (out)
          zlev_z, zlev_diag(lbase), zlev_diag(ltop), elem%Nnode_v, lmesh%NeZ ) ! (in)

        !- Construct a moist neutral profile in the unstable layer

        call construct_moist_neutral_profile( &
          temp_aj, qvap_aj,rhoqvap_aj, rhoprecip_aj, is_profile_converged, & ! (out)
          qvap_z, pres_z, zlev_z, dens_z, temp_z,                          & ! (in)
          int_weight, adjust_mask, elem%Nnode_v, lmesh%NeZ                 ) ! (in)

        if ( .not. is_profile_converged ) then
          LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "Moist neutral profile construction did not converge. Aborting adjustment."
          LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "iter_adj=", iter_adj, "ph=", ph, "ke_xy=", ke_xy
          call PRC_abort
        end if
        
        !- Immediate precipitation removal

        rhoprecip_accum(:,:) = rhoprecip_accum(:,:) + rhoprecip_aj(:,:)

        rhoqdry_aj(:,:) = dens_z(:,:) - rhoqvap_aj(:,:) - rhoprecip_aj(:,:)
        dens_z(:,:) = dens_z(:,:) - rhoprecip_aj(:,:)

        qvap_aj(:,:) = rhoqvap_aj(:,:) / dens_z(:,:)
        qdry_aj(:,:) = rhoqdry_aj(:,:) / dens_z(:,:)

        !- Use the post-precipitation state in the next iteration

        do ke_z=1, lmesh%NeZ
          pres_z(:,ke_z) = ( Rdry * rhoqdry_aj(:,ke_z) + Rvap * rhoqvap_aj(:,ke_z) ) * temp_aj(:,ke_z)
          temp_z(:,ke_z) = temp_aj(:,ke_z)
          qvap_z(:,ke_z) = qvap_aj(:,ke_z)

          rtot_ (:) =  Rdry * qdry_aj(:,ke_z) +     Rvap * qvap_aj(:,ke_z)
          cptot_(:) = CPdry * qdry_aj(:,ke_z) + CP_VAPOR * qvap_aj(:,ke_z)
          pott_z(:,ke_z) = temp_aj(:,ke_z) * ( PRES0 / pres_z(:,ke_z) )**( rtot_(:) / cptot_(:) )
        end do
      end do ! end loop for iteration
  
      if ( .not. adjustment_converged ) then
        LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "Moist convective adjustment did not converge."
        LOG_INFO("atm_phy_cp_dgm_mconv_adjustment_calc_tendency",*) "ph=", ph, "ke_xy=", ke_xy
        call PRC_abort
      end if
      
      
      !- Convert the adjusted state to DG tendencies

      if ( do_adjustment ) then
        do ke_z=1, lmesh%NeZ
          ke = ke_xy + (ke_z-1)*lmesh%Ne2D
          do pz=1, elem%Nnode_v
            p = ph + (pz-1)*elem%Nnode_h1D**2
            DENS_t (p,ke) = ( dens_z(pz,ke_z) - dens_ini(pz,ke_z) ) / dtsec
            RHOT_t (p,ke) = ( dens_z(pz,ke_z) * pott_z(pz,ke_z) - dens_ini(pz,ke_z) * pott_ini(pz,ke_z) ) / dtsec
            RHOQV_t(p,ke) = ( dens_z(pz,ke_z) * qvap_z(pz,ke_z) - rhoqvap_ini(pz,ke_z) ) / dtsec
          end do
        end do
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
      end if

      SFLX_RAIN(ph,ke_xy) = sum( rhoprecip_accum(:,:) * int_weight(:,:) ) / dtsec

    end do ! end loop for ph
    end do ! end loop for ke_xy

    return
  end subroutine atm_phy_cp_dgm_mconv_adjustment_calc_tendency

  subroutine atm_phy_cp_dgm_mconv_adjustment_finalize()
    implicit none
    !----------------------------------------------------
    return
  end subroutine atm_phy_cp_dgm_mconv_adjustment_finalize


!- private subroutines
!OCL SERIAL
  subroutine build_diag_zprofile_from_dg( &
    pott_diag, qvap_diag, pres_diag, zlev_diag, &
    pott_z, qvap_z, pres_z, zlev_z,             &
    lmesh, elem, nlev )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: nlev
    real(RP), intent(out) :: pott_diag(nlev)
    real(RP), intent(out) :: qvap_diag(nlev)
    real(RP), intent(out) :: pres_diag(nlev)
    real(RP), intent(out) :: zlev_diag(nlev)
    real(RP), intent(in) :: pott_z(elem%Nnode_v,lmesh%NeZ)
    real(RP), intent(in) :: qvap_z(elem%Nnode_v,lmesh%NeZ)
    real(RP), intent(in) :: pres_z(elem%Nnode_v,lmesh%NeZ)
    real(RP), intent(in) :: zlev_z(elem%Nnode_v,lmesh%NeZ)

    integer :: pz, ke_z
    integer :: l
    integer :: Nnode_v
    !----------------------------------------------------

    Nnode_v = elem%Nnode_v
    l = 1

    pott_diag(l) = pott_z(1,1)
    qvap_diag(l) = qvap_z(1,1)
    pres_diag(l) = pres_z(1,1)
    zlev_diag(l) = zlev_z(1,1)
    do ke_z=1, lmesh%NeZ
      do pz=2, Nnode_v-1
        l = l + 1
        pott_diag(l) = pott_z(pz,ke_z)
        qvap_diag(l) = qvap_z(pz,ke_z)
        pres_diag(l) = pres_z(pz,ke_z)
        zlev_diag(l) = zlev_z(pz,ke_z)
      end do
      l = l + 1
      if ( ke_z < lmesh%NeZ ) then
        pott_diag(l) = 0.5_RP * ( pott_z(Nnode_v,ke_z) + pott_z(1,ke_z+1) )
        qvap_diag(l) = 0.5_RP * ( qvap_z(Nnode_v,ke_z) + qvap_z(1,ke_z+1) )
        pres_diag(l) = 0.5_RP * ( pres_z(Nnode_v,ke_z) + pres_z(1,ke_z+1) )
        zlev_diag(l) = 0.5_RP * ( zlev_z(Nnode_v,ke_z) + zlev_z(1,ke_z+1) )
      else
        pott_diag(l) = pott_z(Nnode_v,ke_z)
        qvap_diag(l) = qvap_z(Nnode_v,ke_z)
        pres_diag(l) = pres_z(Nnode_v,ke_z)
        zlev_diag(l) = zlev_z(Nnode_v,ke_z)
      end if
    end do
    return
  end subroutine build_diag_zprofile_from_dg

!OCL SERIAL
  subroutine diagnose_convective_layer( lbase, ltop, is_unstable, &
    temp, pott, qvap, pres, zlev, nlev )
    implicit none
    integer, intent(out) :: lbase
    integer, intent(out) :: ltop
    logical, intent(out) :: is_unstable
    integer, intent(in) :: nlev
    real(RP), intent(in) :: temp(nlev)
    real(RP), intent(in) :: pott(nlev)
    real(RP), intent(in) :: qvap(nlev)
    real(RP), intent(in) :: pres(nlev)
    real(RP), intent(in) :: zlev(nlev)
    
    integer :: l
    real(RP) :: qsat(nlev)
    real(RP) :: rh(nlev)
    real(RP) :: hmse_sat(nlev)
    real(RP) :: cptot

    real(RP) :: dz
    real(RP) :: dhmse_sat_dz
    real(RP) :: rh_lyr
    logical :: unstable_pair
    logical :: is_inside
    !------------------------------------

    call ATMOS_SATURATION_pres2qsat_liq( &
      nlev, 1, nlev, & ! (in)
      temp, pres,    & ! (in)
      qsat           ) ! (out)

    do l=1, nlev
      rh(l) = qvap(l) / max(qsat(l), EPS)

      cptot = CPdry * (1.0_RP - qsat(l)) + CP_VAPOR * qsat(l)
      hmse_sat(l) = cptot * temp(l) + GRAV * zlev(l) + LHV0 * qsat(l)
    end do
    
    is_unstable = .false.
    is_inside   = .false.
    lbase = 0; ltop = 0

    do l=1, nlev-1
      dz = zlev(l+1) - zlev(l)
      dhmse_sat_dz = ( hmse_sat(l+1) - hmse_sat(l) ) / dz
      rh_lyr = 0.5_RP * ( rh(l) + rh(l+1) )
      ! write(*,*) "l=", l, "zlev(l)=", zlev(l), "zlev(l+1)=", zlev(l+1), &
      !            "rh_lyr=", rh_lyr, "dhmse_sat_dz=", dhmse_sat_dz
      unstable_pair =     ( rh_lyr >= MCA_RH_TRIGGER )           &
                    .and. ( dhmse_sat_dz < - MCA_HMSE_GRAD_TOL )
    
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
      if ( zlev(ltop) - zlev(lbase) >= MCA_MIN_UNSTABLE_DEPTH ) is_unstable = .true.
    end if
    return
  end subroutine diagnose_convective_layer

!OCL SERIAL
  subroutine make_dg_adjustment_mask( adjust_mask, & ! (out)
    zlev, zbase, ztop, npz, nez                    ) ! (in)
    implicit none
    integer, intent(in) :: npz
    integer, intent(in) :: nez
    logical, intent(out) :: adjust_mask(npz,nez)
    real(RP), intent(in) :: zlev(npz,nez)
    real(RP), intent(in) :: zbase
    real(RP), intent(in) :: ztop

    integer :: ke_z, pz
    !----------------------------------------------------

    adjust_mask(:,:) = .false.
    do ke_z = 1, nez
    do pz = 1, npz
      if ( zlev(pz,ke_z) >= zbase - MCA_Z_TOL .and. &
           zlev(pz,ke_z) <= ztop  + MCA_Z_TOL ) then
        adjust_mask(pz,ke_z) = .true.
      end if
    end do
    end do
    return
  end subroutine make_dg_adjustment_mask

!OCL SERIAL
  subroutine construct_moist_neutral_profile( &
    temp_aj, qvap_aj, rhoqvap_aj, rhoprecip_aj,   &
    converged,                                    &
    qvap, pres, zlev, dens, temp,                 &
    int_weight, adjust_mask, npz, nez             )

    implicit none
    integer, intent(in) :: npz
    integer, intent(in) :: nez
    real(RP), intent(out) :: temp_aj(npz,nez)
    real(RP), intent(out) :: qvap_aj(npz,nez)
    real(RP), intent(out) :: rhoqvap_aj(npz,nez)
    real(RP), intent(out) :: rhoprecip_aj(npz,nez)
    logical, intent(out) :: converged
    real(RP), intent(in) :: qvap(npz,nez)
    real(RP), intent(in) :: pres(npz,nez)
    real(RP), intent(in) :: zlev(npz,nez)
    real(RP), intent(in) :: dens(npz,nez)
    real(RP), intent(in) :: temp(npz,nez)
    real(RP), intent(in) :: int_weight(npz,nez)
    logical, intent(in) :: adjust_mask(npz,nez)

    real(RP) :: energy_target
    real(RP) :: energy_aj
    real(RP) :: residual

    real(RP) :: tbase_lo
    real(RP) :: tbase_hi
    real(RP) :: tbase_mid

    integer :: ke_z, pz
    integer :: iter

    real(RP) :: qsat_aj(npz,nez)

    real(RP) :: water_mass_ini
    real(RP) :: water_mass_aj
    real(RP) :: precip_mass_total
    real(RP) :: precip_mass_local
    real(RP) :: positive_cond_mass
    real(RP) :: precip_scale

    logical :: found_base
    real(RP), parameter :: MCA_TBASE_RANGE = 40.0_RP 

    logical :: is_converged_ref_profile
    logical :: water_feasible
    logical :: residual_evaluated    
    !----------------------------------------------------

    converged = .false.

    ! Initialize the adjusted profile with the original state

    temp_aj(:,:) = temp(:,:)
    qvap_aj(:,:) = qvap(:,:)
    rhoqvap_aj(:,:) = dens(:,:) * qvap(:,:)
    rhoprecip_aj(:,:) = 0.0_RP

    !- Calculate initial column-integrated MSE

    call column_mse( energy_target,     & ! (out)
      temp, qvap, zlev, dens,           & ! (in)
      int_weight, adjust_mask, npz, nez ) ! (in)

    !- Find temperature at the lowest adjusted node

    tbase_mid = 0.0_RP
    found_base = .false.
    do ke_z = 1, nez
      do pz = 1, npz
        if ( adjust_mask(pz,ke_z) ) then
          tbase_mid = temp(pz,ke_z)
          found_base = .true.
          exit
        end if
      end do
      if ( found_base ) exit
    end do
    if ( .not. found_base ) then
      write(*,*) "Error: No adjusted nodes found in the column. Aborting adjustment."
      return
    end if

    !-  Calculate the total water mass in the column before adjustment

    water_mass_ini = 0.0_RP
    do ke_z = 1, nez
    do pz = 1, npz
      if ( adjust_mask(pz,ke_z) ) then
        water_mass_ini = water_mass_ini + int_weight(pz,ke_z) * dens(pz,ke_z) * qvap(pz,ke_z)
      end if
    end do
    end do

    !- Set the initial bounds for the base temperature search

    tbase_lo = tbase_mid - MCA_TBASE_RANGE
    tbase_hi = tbase_mid + MCA_TBASE_RANGE

    water_feasible = .false.
    residual_evaluated = .false.

    !-
    do iter=1, MCA_MAX_ENERGY_ITER

      tbase_mid = 0.5_RP * ( tbase_lo + tbase_hi )

      ! Evaluate the moist adiabat temperature profile for the current base temperature

      temp_aj(:,:) = temp(:,:)
      call evaluate_moist_adiabat_temperature( &
        temp_aj, qsat_aj, is_converged_ref_profile,   & ! (out)
        tbase_mid, pres, zlev, adjust_mask, npz, nez  ) ! (in)

      if ( .not. is_converged_ref_profile ) then
        converged = .false.
        exit
      end if

      ! Set the adjusted vapor to stauration

      qvap_aj(:,:) = qvap(:,:)
      where ( adjust_mask(:,:) )
        qvap_aj(:,:) = qsat_aj(:,:)
      end where

      water_mass_aj = 0.0_RP
      do ke_z=1, nez
      do pz=1, npz
        if ( adjust_mask(pz,ke_z) ) then
          water_mass_aj = water_mass_aj + int_weight(pz,ke_z) * dens(pz,ke_z) * qvap_aj(pz,ke_z)
        end if
      end do
      end do

      ! Reject the current profile becuase it requires more water than is available in the column.
      if ( water_mass_aj > water_mass_ini ) then
        ! write(*,*) "water_mass_ini=", water_mass_ini, "water_mass_aj=", water_mass_aj, "tbase_mid=", tbase_mid
        tbase_hi = tbase_mid
        cycle
      end if
      water_feasible = .true.

      call column_mse( energy_aj,         & ! (out)
        temp_aj, qvap_aj, zlev, dens,     & ! (in)
        int_weight, adjust_mask, npz, nez ) ! (in)

      residual = energy_aj - energy_target
      residual_evaluated = .true.
      ! write(*,*) "iter=", iter, "tbase=", tbase_lo, tbase_mid, tbase_hi,"energy_aj=", energy_aj, "energy_target=", energy_target, "residual=", residual

      if ( abs(residual) <= MCA_ENERGY_RTOL * max(abs(energy_target),1.0_RP) ) then
        converged = .true.
        exit
      end if

      if ( residual > 0.0_RP ) then
        tbase_hi = tbase_mid
      else
        tbase_lo = tbase_mid
      end if
    end do ! end loop for iteration

    if ( .not. converged ) then
      write(*,*) "Error: Moist neutral profile construction did not converge."
      if ( residual_evaluated ) then
        write(*,*) "  residual = ", residual
      else if ( .not. water_feasible ) then
        write(*,*) "  No water-feasible saturated profile was evaluated."
      end if
      return
    end if

    ! Convert the converged state to density variables

    rhoqvap_aj(:,:) = dens(:,:) * qvap_aj(:,:)

    precip_mass_total = max( 0.0_RP, water_mass_ini - water_mass_aj )
    positive_cond_mass = 0.0_RP
    do ke_z=1, nez
    do pz=1, npz
      if ( adjust_mask(pz,ke_z) ) then
        precip_mass_local = max( 0.0_RP, dens(pz,ke_z) * qvap(pz,ke_z) - rhoqvap_aj(pz,ke_z) )

        positive_cond_mass = positive_cond_mass + precip_mass_local * int_weight(pz,ke_z)
      end if
    end do
    end do

    ! Distribute precipitation mass over locally condensing nodes

    if ( positive_cond_mass > EPS ) then
      precip_scale = precip_mass_total / positive_cond_mass
      do ke_z=1, nez
      do pz=1, npz
        if ( adjust_mask(pz,ke_z) ) then
          precip_mass_local = max( 0.0_RP, dens(pz,ke_z) * qvap(pz,ke_z) - rhoqvap_aj(pz,ke_z) )
          rhoprecip_aj(pz,ke_z) = precip_scale * precip_mass_local
        end if
      end do
      end do
    end if

    return
  contains
    ! The column MSE constraint is imposed on the fixed-density, vapor-only adjusted state. 
    ! Energy carried away by precipitation is neglected in the present simplified adjustment scheme.  
    subroutine column_mse( energy,          &
      temp_, qv_, zlev_, dens_,             &
      int_weight_, adjust_mask_, npz_, nez_ )
      implicit none
      integer, intent(in) :: npz_
      integer, intent(in) :: nez_
      real(RP), intent(out) :: energy
      real(RP), intent(in) :: temp_(npz_,nez_)
      real(RP), intent(in) :: qv_(npz_,nez_)
      real(RP), intent(in) :: zlev_(npz_,nez_)
      real(RP), intent(in) :: dens_(npz_,nez_)
      real(RP), intent(in) :: int_weight_(npz_,nez_)
      logical, intent(in) :: adjust_mask_(npz_,nez_)

      integer :: ke_z_, pz_
      real(RP) :: qdry_
      real(RP) :: cptot_
      !------------------------------------------------

      energy = 0.0_RP

      do ke_z_ = 1, nez_
      do pz_ = 1, npz_
        if ( adjust_mask_(pz_,ke_z_) ) then
          qdry_ = 1.0_RP - qv_(pz_,ke_z_) !- qcon_(pz_,ke_z_)

          cptot_ = CPdry    * qdry_ &
                 + CP_VAPOR * qv_(pz_,ke_z_)
                                              ! Liquid-water sensible heat is neglected.
          
          energy = energy + int_weight_(pz_,ke_z_) * dens_(pz_,ke_z_) * &
            ( cptot_ * temp_(pz_,ke_z_)            &
            + GRAV * zlev_(pz_,ke_z_)              &
            + LHV0 * qv_(pz_,ke_z_)                )
        end if
      end do
      end do
      return
    end subroutine column_mse
  end subroutine construct_moist_neutral_profile

!OCL SERIAL
  subroutine evaluate_moist_adiabat_temperature( &
    temp_ref, qsat_ref, is_converged_ref_profile, & ! (out)
    tbase, pres, zlev, mask, npz, nez             ) ! (in)
    implicit none
    integer, intent(in) :: npz
    integer, intent(in) :: nez
    real(RP), intent(inout) :: temp_ref(npz,nez)
    real(RP), intent(inout) :: qsat_ref(npz,nez)
    logical, intent(out) :: is_converged_ref_profile
    real(RP), intent(in) :: tbase
    real(RP), intent(in) :: pres(npz,nez)
    real(RP), intent(in) :: zlev(npz,nez)
    logical, intent(in) :: mask(npz,nez)

    integer :: ke_z, pz

    real(RP) :: t_prev, p_prev, z_prev
    real(RP) :: t_now, p_now, z_now

    logical :: started
    logical :: is_converge_local
    !----------------------------------------------------

    started = .false.
    is_converged_ref_profile = .true.

    do ke_z = 1, nez
    do pz = 1, npz
      if ( .not. mask(pz,ke_z) ) cycle

      z_now = zlev(pz,ke_z)
      p_now = pres(pz,ke_z)

      if ( .not. started ) then
        t_now = tbase
        started = .true.
      else if ( abs(z_now-z_prev) <= MCA_Z_TOL ) then
        t_now = t_prev
      else
        call advance_moist_adiabat( t_now, is_converge_local, & ! (out)
          t_prev, p_prev, z_prev, p_now, z_now                ) ! (in)
        
        if ( .not. is_converge_local ) then
          is_converged_ref_profile = .false.
          return
        end if
      end if

      temp_ref(pz,ke_z) = t_now
      call ATMOS_SATURATION_pres2qsat_liq( &
        t_now, p_now,          & ! (in)
        qsat_ref(pz,ke_z)      ) ! (out)

      t_prev = t_now 
      p_prev = p_now
      z_prev = z_now
    end do
    end do

    return
  end subroutine evaluate_moist_adiabat_temperature

!OCL SERIAL
  subroutine advance_moist_adiabat( temp1, converged,  & ! (out)
    temp0, pres0, zlev0, pres1, zlev1                  ) ! (in)
    implicit none
    real(RP), intent(out) :: temp1
    logical,  intent(out) :: converged
    real(RP), intent(in)  :: temp0
    real(RP), intent(in)  :: pres0
    real(RP), intent(in)  :: zlev0
    real(RP), intent(in)  :: pres1
    real(RP), intent(in)  :: zlev1
    integer :: iter

    real(RP) :: dz
    real(RP) :: p_mid, t_mid
    real(RP) :: tguess, tnew
    real(RP) :: qsat_mid
    real(RP) :: Gamma_m
    real(RP) :: error

    ! Local temperature iteration tolerances
    real(RP), parameter :: MCA_TEMP_ATOL = 1.0E-8_RP  ! [K]
    real(RP), parameter :: MCA_TEMP_RTOL = 1.0E-10_RP
    real(RP), parameter :: MCA_PRES_RTOL = 1.0E-12_RP
    !----------------------------------------------------

    converged = .false.
    dz = zlev1 - zlev0

    if ( abs(dz) <= MCA_Z_TOL ) then
      temp1 = temp0
      converged = .true.
      return
    end if

    !- Evaluate the midpoint pressure for the moist adiabatic lapse rate calculation

    if ( abs(pres1 - pres0) <= MCA_PRES_RTOL * max(pres0, pres1) ) then
      p_mid = 0.5_RP * ( pres0 + pres1 )
    else
      p_mid = ( pres1 - pres0 ) / log( pres1 / pres0 )
    end if

    !- Calculate an initial guess

    
    call ATMOS_SATURATION_pres2qsat_liq( &
      temp0, p_mid,      & ! (in)
      qsat_mid           ) ! (out)
    
    call moist_adiabatic_lapse_rate( Gamma_m, & ! (out)
      temp0, qsat_mid )                         ! (in)

    tguess = temp0 - Gamma_m * dz
    tnew   = tguess

    !----------------------------------------------
    ! Iteration to solve the implicit equation for moist adiabatic lapse rate
    ! T1^* = T0 - Gamma_m ( (T0 + T1^*)/2, pmid ) * dz
    !----------------------------------------------

    do iter = 1, MCA_MAX_LOCAL_ITER
      
      t_mid = 0.5_RP * ( temp0 + tguess )

      call ATMOS_SATURATION_pres2qsat_liq( t_mid, p_mid, & ! (in)
                                           qsat_mid      ) ! (out)

      call moist_adiabatic_lapse_rate( Gamma_m, & ! (out)
        t_mid, qsat_mid )                         ! (in)

        tnew = temp0 - Gamma_m * dz

        error = abs(tnew - tguess)
        if (  error <= MCA_TEMP_ATOL                                 &
                       + MCA_TEMP_RTOL * max(abs(tnew), abs(tguess)) ) then
          temp1 = tnew
          converged = .true.
          return
        end if

        tguess = tnew
        ! If you want to stabilize the fixed-point iteration,  
        !   tguess = ( 1.0_RP - MCA_LOCAL_RELAX ) * tguess + MCA_LOCAL_RELAX * tnew
    end do

    temp1 = tnew
    return
  contains
    subroutine moist_adiabatic_lapse_rate( GammaM, &
      temp, qsat )
      implicit none
      real(RP), intent(out) :: GammaM
      real(RP), intent(in) :: temp
      real(RP), intent(in) :: qsat

      real(RP) :: Rvap
      real(RP) :: Rtot, Cptot
      real(RP) :: EpsV
      real(RP) :: fac, A
      !------------------------------------------------
      Rvap = CP_VAPOR - CV_VAPOR
      Rtot = Rdry * ( 1.0_RP - qsat ) + Rvap * qsat
      Cptot = CPdry * ( 1.0_RP - qsat ) + CP_VAPOR * qsat
      EpsV = Rdry / Rvap
      fac = 1.0_RP + ( 1.0_RP - EpsV ) / EpsV * qsat
      A = LHV0 + ( CP_VAPOR - CPdry ) * temp

      GammaM = GRAV &
             * ( 1.0_RP + A * qsat * fac / ( Rtot * temp ) ) &
             / ( Cptot + A * LHV0 * qsat * fac / ( Rvap * temp**2 ) )
      return
    end subroutine moist_adiabatic_lapse_rate
  end subroutine advance_moist_adiabat
end module scale_atm_phy_cp_dgm_mconv_adjustment

