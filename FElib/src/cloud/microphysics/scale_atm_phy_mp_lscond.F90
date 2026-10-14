!> module FElib / Atmosphere / Physics cloud microphysics / Large-scale condensation
!!
!! @par Description
!!     A module to provide a large-scale condensation scheme
!!
!! @par Reference
!!
!! @author Yuta Kawai, Team SCALE
!!
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_mp_lscond
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort 
  use scale_const, only: &
    Rdry => CONST_Rdry,    &
    Rvap => CONST_Rvap,    &
    CPdry => CONST_CPdry,  & 
    CVdry => CONST_CVdry,  & 
    PRES00 => CONST_PRE00, &
    Grav => CONST_GRAV
  use scale_atmos_hydrometeor, only: &
    CP_VAPOR, &
    CP_WATER, &
    CP_ICE,   &
    CV_VAPOR, &
    CV_WATER, &
    CV_ICE,   &
    LHV
  
  use scale_element_base, only: &
    ElementBase1D, ElementBase3D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  public :: ATMOS_PHY_MP_lscond_setup
  public :: ATMOS_PHY_MP_lscond_adjustment

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !
  integer, private, parameter :: QA_MP  = 2

  integer,                parameter, public :: ATMOS_PHY_MP_lscond_ntracers = QA_MP
  integer,                parameter, public :: ATMOS_PHY_MP_lscond_nwaters = 1
  integer,                parameter, public :: ATMOS_PHY_MP_lscond_nices = 0
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_MP_lscond_tracer_names(QA_MP) = (/ &
       'QV', &
       'QC'  /)
  character(len=H_MID)  , parameter, public :: ATMOS_PHY_MP_lscond_tracer_descriptions(QA_MP) = (/ &
       'Ratio of Water Vapor mass to total mass (Specific humidity)', &
       'Ratio of Cloud Water mass to total mass                    '  /)
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_MP_lscond_tracer_units(QA_MP) = (/ &
       'kg/kg',  &
       'kg/kg'   /)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !


  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter   :: I_QV = 1
  integer,  private, parameter   :: I_QC = 2

  integer,  private, parameter   :: I_hyd_QC =  1

  logical,  private              :: flag_liquid = .true.     ! warm rain
  logical,  private              :: couple_aerosol = .false. ! Consider CCN effect ?

  ! real(RP), private, parameter   :: re_qc =  8.E-6_RP        ! effective radius for cloud water

contains

  !-----------------------------------------------------------------------------
  !> Setup a module for a large-scale condensation scheme
  !!
  subroutine ATMOS_PHY_MP_lscond_setup
    use scale_prc, only: &
       PRC_abort
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_lscond_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_MP_lscond_setup",*) 'large-scale condensation'

    if( couple_aerosol ) then
       LOG_ERROR("ATMOS_PHY_MP_lscond_setup",*) 'MP_aerosol_couple should be .false. for large-scale condensation type MP!'
       call PRC_abort
    endif

    return
  end subroutine ATMOS_PHY_MP_lscond_setup

  !> Calculate a state after the saturation process
  !!
!OCL SERIAL
  subroutine ATMOS_PHY_MP_lscond_adjustment( &
    KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! (in)
    DENS, PRES,                         & ! (in)
    dt,                                 & ! (in)
    TEMP, QTRC, CPtot, CVtot,           & ! (inout)
    RHOE_t, EVAPORATE                   ) ! (out)

    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_pres2qsat_liq
    use scale_const, only: &
      CL => CONST_CL
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: PRES(KA,IA,JA)
    real(RP), intent(in) :: dt
    real(RP), intent(inout) :: TEMP(KA,IA,JA) 
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP)
    real(RP), intent(inout) :: CPtot(KA,IA,JA) 
    real(RP), intent(inout) :: CVtot(KA,IA,JA) 
    real(RP), intent(out) :: RHOE_t(KA,IA,JA)
    real(RP), intent(out) :: EVAPORATE(KA,IA,JA)

    real(RP) :: Qsat(KA,IA,JA)
    real(RP) :: d_qc(KA)
    real(RP) :: d_cv(KA)
    real(RP) :: d_cp(KA)
    real(RP) :: d_e(KA)
    real(RP) :: cptot_(KA)
    real(RP) :: cvtot_(KA)
    real(RP) :: coef1(KA)
    real(RP) :: coef2

    integer :: i, j

    real(RP) :: rdt
    !----------------------------------------------

    call ATMOS_SATURATION_pres2qsat_liq( &
      KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! (in)
      TEMP, PRES,                         & ! (in)
      Qsat                                ) ! (out)


    coef2 = Rvap / LHV
    rdt = 1.0_RP / dt

    !$omp parallel do private(i,j, d_qc, d_cv, d_cp, d_e, cptot_, cvtot_, coef1) collapse(2)
    do j=JS, JE
    do i=IS, IE
  !    coef1 = LHV**2 / ( CVDry * Rvap )
      coef1(:) = LHV * ( LHV + ( CV_VAPOR - CV_WATER ) * TEMP(:,i,j) ) / ( CVtot(:,i,j) * Rvap )

      d_qc(:) = max( ( QTRC(:,i,j,I_QV) - Qsat(:,i,j) ) / ( 1.0_RP + coef1 / TEMP(:,i,j)**2 * ( 1.0_RP - coef2 * TEMP(:,i,j) ) * Qsat(:,i,j) ), &
                     0.0_RP )

      QTRC(:,i,j,I_QV) = QTRC(:,i,j,I_QV) - d_qc(:)
      QTRC(:,i,j,I_QC) = QTRC(:,i,j,I_QC) + d_qc(:)

      d_cp(:) = - CP_VAPOR * d_qc(:) &
                + CP_WATER * d_qc(:)
      cptot_(:) = CPtot(:,i,j) + d_cp(:)

      d_cv(:) = - CV_VAPOR * d_qc(:) &
                + CV_WATER * d_qc(:)
      cvtot_(:) = CVtot(:,i,j) + d_cv(:)

      d_e(:) = LHV * d_qc(:)
      RHOE_t(:,i,j) = DENS(:,i,j) * d_e(:) * rdt


      TEMP(:,i,j) = ( TEMP(:,i,j) * CVtot(:,i,j) + d_e(:) ) / cvtot_(:)
      CPtot(:,i,j) = cptot_(:)
      CVtot(:,i,j) = cvtot_(:)
    end do
    end do

    return
  end subroutine ATMOS_PHY_MP_lscond_adjustment  
end module scale_atm_phy_mp_lscond
