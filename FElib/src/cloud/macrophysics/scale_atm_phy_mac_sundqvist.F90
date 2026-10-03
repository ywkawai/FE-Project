!-------------------------------------------------------------------------------
!> module atmosphere / physics / cloud microphysics / Sundqvist-type cloud scheme
!!
!! @par Description
!!   Simple prognostic cloud-water scheme for coarse-resolution GCM simulations.
!!
!!   Prognostic water tracers: water vapor and cloud liquid water
!!   Cloud fraction is diagnosed using a Sundqvist-type RH relation.
!!   Cloud condensation / evaporation: QV <-> QC
!!   Precipitation: QC -> diagnostic precipitation flux
!!
!!   Precipitating water is NOT retained as a prognostic QR tracer.
!!   Generated precipitation falls through the atmospheric column within one physics call. 
!!   Part of the precipitation may re-evaporate.
!!
!! @author Yuta Kawai, Team SCALE
!!
!! @Reference
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_mac_sundqvist
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use scale_const, only: &
    EPS => CONST_EPS
  use scale_atmos_hydrometeor, only: &
    LHV,      &
    CP_VAPOR, &
    CP_WATER, &
    CV_VAPOR, &
    CV_WATER
  use scale_atmos_saturation, only: &
    SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !++ Public procedures
  !
  public :: ATMOS_PHY_mac_sundqvist_setup
  public :: ATMOS_PHY_mac_sundqvist_finalize
  public :: ATMOS_PHY_mac_sundqvist_adjustment
  public :: ATMOS_PHY_mac_sundqvist_cloud_fraction
  public :: ATMOS_PHY_mac_sundqvist_incloud_water
  public :: ATMOS_PHY_mac_sundqvist_effective_radius

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !
  
  integer, private, parameter :: QA_CLD = 3

  integer, parameter, public :: ATMOS_PHY_mac_sundqvist_ntracers = QA_CLD
  integer, parameter, public :: ATMOS_PHY_mac_sundqvist_nwaters  = 2
  integer, parameter, public :: ATMOS_PHY_mac_sundqvist_nices    = 0

  character(len=H_SHORT), parameter, public :: &
       ATMOS_PHY_mac_sundqvist_tracer_names(QA_CLD) = (/ &
       'QV', &
       'QC', &
       'QR'  /)

  character(len=H_MID), parameter, public :: &
       ATMOS_PHY_mac_sundqvist_tracer_descriptions(QA_CLD) = (/ &
       'Ratio of water vapor mass to total mass        ', &
       'Ratio of cloud liquid mass to total mass       ', &
       'Ratio of precipitating water mass to total mass'  /)

  character(len=H_SHORT), parameter, public :: &
       ATMOS_PHY_mac_sundqvist_tracer_units(QA_CLD) = (/ &
       'kg/kg', &
       'kg/kg', &
       'kg/kg'  /)


  !-----------------------------------------------------------------------------
  !++ Private parameters & variables
  !

  ! Tracer indices
  integer, private, parameter :: I_QV = 1
  integer, private, parameter :: I_QC = 2
  integer, private, parameter :: I_QR = 3

  !- Scheme parameters

  real(RP), private :: RH_CRIT         = 0.80_RP
  real(RP), private :: TAU_COND        = 600.0_RP   !< Relaxation time scale for subgrid condensation [s]
  real(RP), private :: TAU_EVAP_CLOUD  = 600.0_RP   !< Relaxation time scale for cloud-water evaporation [s]
  real(RP), private :: QC_CRIT         = 5.0E-4_RP  !< Characteristic in-cloud QC for precipitation formation [kg/kg]
  real(RP), private :: TAU_PRECIP      = 1000.0_RP  !< Precipitation conversion timescale [s]
  real(RP), private :: TAU_EVAP_PRECIP = 1800.0_RP  !< Evaporation time scale for falling precipitation [s]
  real(RP), private :: RE_CLOUD        = 8.0E-6_RP  !< Fixed effective radius of cloud liquid water [m]


  real(RP), private, parameter :: CLDFRAC_MIN = 1.0E-6_RP

contains

  !> Setup a module to configure the Sundqvist-type cloud scheme
  subroutine ATMOS_PHY_mac_sundqvist_setup()
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_ATMOS_PHY_MAC_sundqvist / &
      RH_CRIT,         &
      TAU_COND,        &
      TAU_EVAP_CLOUD,  &
      QC_CRIT,         &
      TAU_PRECIP,      &
      TAU_EVAP_PRECIP, &
      RE_CLOUD

    integer :: ierr
    !----------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MAC_sundqvist_setup",*) 'Setup Sundqvist-type simple cloud scheme'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MAC_sundqvist,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MAC_sundqvist_setup",*) 'Not found namelist. Default used.'
    else if( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MAC_sundqvist_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MAC_sundqvist. Check!'
       call PRC_abort
    end if
    LOG_NML(PARAM_ATMOS_PHY_MAC_sundqvist)

    return
  end subroutine ATMOS_PHY_MAC_sundqvist_setup

  subroutine ATMOS_PHY_MAC_sundqvist_finalize()
    implicit none
    !--------------------------------------------------
    return
  end subroutine ATMOS_PHY_MAC_sundqvist_finalize

  !> Sundqvist-type cloud adjustment
  !!
  !! QV <-> QC
  !! QC -> diagnostic precipitation
  !! PRECIP_GEN is NOT a prognostic rain-water tracer.
  !!
  subroutine ATMOS_PHY_MAC_sundqvist_adjustment( &
     KA, KS, KE, IA, IS, IE, JA, JS, JE,         & ! (in)
     DENS, PRES, dt,                             & ! (in)
     TEMP, QTRC, CPtot, CVtot,                   & ! (inout)
     EVAPORATE, RHOE_t                           ) ! (out)

    implicit none
    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, IS, IE
    integer,  intent(in) :: JA, JS, JE
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: PRES(KA,IA,JA)
    real(RP), intent(in) :: dt
    real(RP), intent(inout) :: TEMP (KA,IA,JA)
    real(RP), intent(inout) :: QTRC (KA,IA,JA,QA_CLD)
    real(RP), intent(inout) :: CPtot(KA,IA,JA)
    real(RP), intent(inout) :: CVtot(KA,IA,JA)
    real(RP), intent(out) :: EVAPORATE(KA,IA,JA)
    real(RP), intent(out) :: RHOE_t   (KA,IA,JA)

    real(RP) :: dens_
    real(RP) :: temp_work, qv_work, qc_work
    real(RP) :: cvtot_work, cptot_work

    real(RP) :: qsat
    real(RP) :: rh
    real(RP) :: cf

    real(RP) :: qv_target
    real(RP) :: qc_incloud

    real(RP) :: dq_cond
    real(RP) :: dq_evap
    real(RP) :: dq_precip

    real(RP) :: dq_cloud
    real(RP) :: dq_sat    ! Additional condensation by supersaturation correction

    real(RP) :: qv_t_phase
    real(RP) :: cp_t, cv_t

    real(RP) :: cptot_new
    real(RP) :: cvtot_new

    real(RP) :: auto_factor
    real(RP) :: rdt

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('ATMOS_PHY_MAC_sundqvist_adjustment', 3)

    rdt = 1.0_RP / dt

   !$omp parallel do collapse(2)                            &
   !$omp private(k, dens_,temp_work,qv_work,qc_work,        &
   !$omp         cptot_work,cvtot_work,                     &
   !$omp         qsat,rh,cf,qv_target,qc_incloud,           &
   !$omp         dq_cond,dq_evap,dq_precip,dq_cloud,dq_sat, &
   !$omp         qv_t_phase,cp_t,cv_t,cptot_new,cvtot_new,  &
   !$omp         auto_factor )
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dens_ = DENS(k,i,j)

       temp_work = TEMP(k,i,j)
       qv_work = QTRC(k,i,j,I_QV)
       qc_work = QTRC(k,i,j,I_QC)

       cvtot_work = CVtot(k,i,j)
       cptot_work = CPtot(k,i,j)

       !- 1. Saturation specific humidity and RH       

       call SATURATION_dens2qsat_liq( temp_work, dens_, &
          qsat )

       rh = qv_work / qsat

       !- 2. Sundqvist-type cloud fraction
       !       C = 0                             RH <= RHcrit
       !       C = 1 - sqrt[(1-RH)/(1-RHcrit)]   otherwise
       !       C = 1                             RH >= 1
       if ( rh <= RH_CRIT ) then
          cf = 0.0_RP
       else if ( rh >= 1.0_RP ) then
          cf = 1.0_RP
       else
          cf = 1.0_RP - sqrt( (1.0_RP-rh) / max(1.0_RP-RH_CRIT,EPS) )
       endif

       !- 3. Target grid-mean vapor
       ! NOTE: This is a SCALE-DG prototype closure, which is NOT the complete original Sundqvist condensation formulation.
       !
       ! Simple subgrid interpretation:
       !  cloudy fraction: RH = 1
       !  clear fraction : RH = RH_CRIT

       qv_target = ( cf + ( 1.0_RP - cf ) * RH_CRIT ) * qsat

       !- 4. Sundqvist-type condensation QV -> QC

       if ( qv_work > qv_target ) then
         dq_cond = ( qv_work - qv_target ) / TAU_COND
         dq_cond = min(dq_cond, qv_work * rdt)
       else
         dq_cond = 0.0_RP
       end if

       !- 5. Sundqvist-type cloud evaporation QC -> QV
       ! This is the process that allows a fractional-cloud grid box to remain subsaturated in the grid-mean sense.

       if ( qv_work < qv_target .and. qc_work > 0.0_RP ) then
         dq_evap = ( qv_target - qv_work ) / TAU_EVAP_CLOUD
         dq_evap = min(dq_evap, qc_work * rdt)
       else
         dq_evap = 0.0_RP
       end if

       !- 6. Apply Sundqvist condensation / evaporation to working state
       dq_cloud = ( dq_cond - dq_evap ) * dt

       qv_work = qv_work - dq_cloud
       qc_work = qc_work + dq_cloud
       
       cptot_work = cptot_work + ( CP_WATER - CP_VAPOR ) * dq_cloud
       cvtot_work = cvtot_work + ( CV_WATER - CV_VAPOR ) * dq_cloud
       temp_work = ( temp_work * CVtot(k,i,j) + LHV * dq_cloud ) / cvtot_work

       !- 7. One-sided supersaturation correction
       
       call mac_sundqvist_supersat_correct( dens_, & ! (in)
          temp_work, qv_work, qc_work,               & ! (inout)
          cptot_work, cvtot_work,                    & ! (inout)
          dq_sat                                     ) ! (out)

       !- 8. In-cloud condensate

       if ( cf > CLDFRAC_MIN ) then
          qc_incloud = qc_work / cf
       else
          qc_incloud = 0.0_RP
       end if

       !- 9. Diagnostic precipitation formation, which removes water from QC

       if ( qc_work > 0.0_RP ) then
         ! Smooth Sundqvist-type autoconversion
         auto_factor = 1.0_RP - exp( -( qc_incloud / max(QC_CRIT, EPS) )**2 )
         dq_precip = qc_work / TAU_PRECIP * auto_factor
          
         dq_precip = min( dq_precip, qc_work * rdt )
       else
         dq_precip = 0.0_RP
       end if

       qc_work = qc_work - dq_precip * dt
       
       !- 10. Final state update

       QTRC(k,i,j,I_QV) = qv_work
       QTRC(k,i,j,I_QC) = qc_work
       QTRC(k,i,j,I_QR) = dq_precip * dt

       TEMP(k,i,j) = temp_work
       CPtot(k,i,j) = cptot_work
       CVtot(k,i,j) = cvtot_work

       qv_t_phase = - dq_cond + dq_evap &
                    - dq_sat * rdt
      
       RHOE_t(k,i,j) = dens_ * ( - LHV * qv_t_phase )
    end do
    end do
    end do

    call PROF_rapend('ATMOS_PHY_MAC_sundqvist_adjustment', 3)
    return
  end subroutine ATMOS_PHY_MAC_sundqvist_adjustment

  !> Diagnose cloud fraction from final atmospheric state
!OCL SERIAL
  subroutine ATMOS_PHY_mac_sundqvist_cloud_fraction( &
    KA, KS, KE, IA, IS, IE, JA, JS, JE,             &
    DENS, TEMP, QV,                                 &
    CLDFRAC                                         )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: TEMP(KA,IA,JA)
    real(RP), intent(in) :: QV  (KA,IA,JA)
    real(RP), intent(out) :: CLDFRAC(KA,IA,JA)

    real(RP) :: qsat
    real(RP) :: rh
    real(RP) :: arg

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(qsat, rh, arg) collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
      call SATURATION_dens2qsat_liq( TEMP(k,i,j), DENS(k,i,j), &
        qsat )

       rh = QV(k,i,j) / max(qsat,EPS)

       if ( rh <= RH_CRIT ) then
         CLDFRAC(k,i,j) = 0.0_RP
       else if ( rh >= 1.0_RP ) then
         CLDFRAC(k,i,j) = 1.0_RP
       else
         arg = (1.0_RP - rh) / max(1.0_RP - RH_CRIT, EPS)
         CLDFRAC(k,i,j) = 1.0_RP - sqrt(arg)
       endif
    end do
    end do
    end do

    return
  end subroutine ATMOS_PHY_mac_sundqvist_cloud_fraction

  !> Calculate in-cloud liquid water
  !! QC_INCLOUD = QC / cloud_fraction
!OCL SERIAL
  subroutine ATMOS_PHY_mac_sundqvist_incloud_water( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE,            &
       QC, CLDFRAC,                                   &
       QC_INCLOUD                                     )

    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: QC     (KA,IA,JA)
    real(RP), intent(in) :: CLDFRAC(KA,IA,JA)
    real(RP), intent(out) :: QC_INCLOUD(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
      if ( CLDFRAC(k,i,j) > CLDFRAC_MIN ) then
        QC_INCLOUD(k,i,j) = QC(k,i,j) / CLDFRAC(k,i,j)
      else
        QC_INCLOUD(k,i,j) = 0.0_RP
      end if
    end do
    end do
    end do

    return
  end subroutine ATMOS_PHY_mac_sundqvist_incloud_water


  !> Calculate the effective radius
  !!
  !! Fixed liquid-cloud effective radius.
  !! Same simple philosophy as current Kessler implementation.
!OCL SERIAL
  subroutine ATMOS_PHY_mac_sundqvist_effective_radius( &
    KA, KS, KE, IA, IS, IE, JA, JS, JE,               &
    Re                                                )
    use scale_atmos_hydrometeor, only: &
      N_HYD, &
      I_HC
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(out) :: Re(KA,IA,JA,N_HYD)

    real(RP), parameter :: m2cm = 100.0_RP
    !---------------------------------------------------------------------------

    !$omp parallel workshare
    Re(:,:,:,:) = 0.0_RP
    Re(:,:,:,I_HC) = RE_CLOUD * m2cm
    !$omp end parallel workshare
    return
  end subroutine ATMOS_PHY_mac_sundqvist_effective_radius


!--- private subroutines --------------------------------------------------

  !> One-sided supersaturation correction
  !! Remove grid-mean supersaturation by additional condensation.
  !!
  !! Only QV -> QC is allowed.
  !! Cloud-water evaporation QC -> QV is NEVER performed here.
  !!
!OCL SERIAL
  subroutine mac_sundqvist_supersat_correct( &
    dens,                            & ! (in)
    temp, qv, qc,                    & ! (inout)
    cptot, cvtot,                    & ! (inout)
    dq_sat                           ) ! (out)

    use scale_atmos_saturation, only: &
      SATURATION_dqs_dtem_dens_liq    => ATMOS_SATURATION_dqs_dtem_dens_liq     
    implicit none
    real(RP), intent(in) :: dens
    real(RP), intent(inout) :: temp
    real(RP), intent(inout) :: qv
    real(RP), intent(inout) :: qc
    real(RP), intent(inout) :: cptot
    real(RP), intent(inout) :: cvtot
    real(RP), intent(out) :: dq_sat

    real(RP) :: temp0
    real(RP) :: qv0, qc0
    real(RP) :: cptot0, cvtot0

    real(RP) :: dq, dq_new

    real(RP) :: temp_work
    real(RP) :: qv_work
    real(RP) :: cvtot_work

    real(RP) :: qsat0
    real(RP) :: qsat_work
    real(RP) :: dqsat_dtemp

    real(RP) :: residual

    real(RP) :: dres_ddq
    real(RP) :: dtemp_ddq
    real(RP) :: dcv

    integer :: iter

    integer, parameter :: MAX_ITER = 8
    real(RP), parameter :: Q_TOL = 1.0E-12_RP
    !---------------------------------------------------------------------------

    temp0  = temp
    qv0    = qv
    qc0    = qc
    cptot0 = cptot
    cvtot0 = cvtot

    dcv = CV_WATER - CV_VAPOR

    ! Saturation state before correction

    call SATURATION_dens2qsat_liq( &
        temp0, dens,              &
        qsat0                     )

    ! One-sided correction: subsaturated states are left untouched.
    if ( qv0 <= qsat0 ) return

    ! Initial guess
    dq = min(qv0 - qsat0, qv0)

    !- Newton iteration : F = qv,0 - dq - qsat(T(dq),rho) = 0

    do iter = 1, MAX_ITER
      qv_work = qv0 - dq

      cvtot_work = cvtot0 + dcv * dq
      temp_work = ( temp0 * cvtot0 + LHV * dq ) / cvtot_work

      call SATURATION_dens2qsat_liq( temp_work, dens, &
        qsat_work )

      residual = qv_work - qsat_work

      if ( abs(residual) <= Q_TOL ) exit

      call SATURATION_dqs_dtem_dens_liq( temp_work, dens, &
        dqsat_dtemp )

      ! dT / d(dq) when assumuing moist internal energy is conserved
      dtemp_ddq = ( LHV - temp_work * dcv ) / cvtot_work
      ! dF / d(dq)
      dres_ddq = - 1.0_RP - dqsat_dtemp * dtemp_ddq

      ! Newton step
      dq_new = dq - residual / dres_ddq

      ! Condensation cannot be negative and cannot exceed available vapor.
      dq = min(max(dq_new,0.0_RP),qv0)
    end do

    dq_sat = dq

    !- Final corrected state

    qv = qv0 - dq_sat
    qc = qc0 + dq_sat

    cptot = cptot0 + ( CP_WATER - CP_VAPOR ) * dq_sat
    cvtot = cvtot0 + ( CV_WATER - CV_VAPOR ) * dq_sat
    temp = ( temp0 * cvtot0 + LHV * dq_sat ) / cvtot
    return
  end subroutine mac_sundqvist_supersat_correct

end module scale_atm_phy_mac_sundqvist