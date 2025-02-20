!-------------------------------------------------------------------------------
!> module atmosphere / physics / surface / simple
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!          Constant bulk coefficient
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_sf_bulk_simple
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
    PRES00 => CONST_PRE00, &
    CpDry => CONST_CPdry,  &
    EPSvap => CONST_EPSvap

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_simple_setup
  public :: ATMOS_PHY_SF_simple_flux

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: ATMOS_PHY_SF_U_maxM      =  100.0_RP ! maximum limit of absolute velocity for momentum [m/s]
  real(RP), private            :: ATMOS_PHY_SF_U_minM      =    0.0_RP ! minimum limit of absolute velocity for momentum [m/s]

  real(RP), private            :: ATMOS_PHY_SF_Const_Cm    = 0.0011_RP ! constant bulk coefficient for momentum    [NIL]
  real(RP), private            :: ATMOS_PHY_SF_Const_Ch    = 0.0044_RP ! constant bulk coefficient for heat        [NIL]
  real(RP), private            :: ATMOS_PHY_SF_Const_Ce    = 0.0044_RP ! constant bulk coefficient for evaporation [NIL]

  real(RP), private :: ATMOS_PHY_SF_BULK_beta = 1.0_RP ! evaporation efficiency (0-1)

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_simple_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_ATMOS_PHY_SF_SIMPLE / &
       ATMOS_PHY_SF_U_minM,         &
       ATMOS_PHY_SF_Const_Cm,       &
       ATMOS_PHY_SF_Const_Ch,       &
       ATMOS_PHY_SF_Const_Ce,       &
       ATMOS_PHY_SF_BULK_beta

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SF_simple_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_SF_simple_setup",*) 'Simple flux'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_SIMPLE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_SF_simple_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_SF_simple_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_SF_SIMPLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_SF_SIMPLE)

    return
  end subroutine ATMOS_PHY_SF_simple_setup

  !-----------------------------------------------------------------------------
  !> Constant flux
  subroutine ATMOS_PHY_SF_simple_flux( &
       IA, IS, IE, JA, JS, JE, &
       ATM_W, ATM_U, ATM_V, ATM_TEMP, ATM_PRES, ATM_QV, &
       SFC_DENS, SFC_TEMP, SFC_PRES,                    &
       ATM_Z1,                                          &
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_LH,     &
       SFLX_QV,                                         &
       U10, V10                                         )
    use scale_const, only: &
       PI    => CONST_PI
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_atmos_saturation, only: &
       SATURATION_psat_all => ATMOS_SATURATION_psat_all
    use scale_time, only: &
       TIME_NOWSEC
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: ATM_W   (IA,JA) ! velocity w  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in) :: ATM_U   (IA,JA) ! velocity u  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in) :: ATM_V   (IA,JA) ! velocity v  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in) :: ATM_TEMP(IA,JA) ! temperature at the lowermost layer (cell center) [K]
    real(RP), intent(in) :: ATM_PRES(IA,JA) ! pressure    at the lowermost layer (cell center) [Pa]
    real(RP), intent(in) :: ATM_QV  (IA,JA) ! qv          at the lowermost layer (cell center) [kg/kg]
    real(RP), intent(in) :: SFC_DENS(IA,JA) ! density     at the surface atmosphere [kg/m3]
    real(RP), intent(in) :: SFC_TEMP(IA,JA) ! tempertire  at the surface atmosphere [K]
    real(RP), intent(in) :: SFC_PRES(IA,JA) ! pressure    at the surface atmosphere [Pa]
    real(RP), intent(in) :: ATM_Z1  (IA,JA) ! height of the lowermost grid from surface (cell center) [m]

    real(RP), intent(out) :: SFLX_MW(IA,JA) ! surface flux for z-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out) :: SFLX_MU(IA,JA) ! surface flux for x-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out) :: SFLX_MV(IA,JA) ! surface flux for y-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out) :: SFLX_SH(IA,JA) ! surface flux for sensible heat (area center)   [J/m2/s]
    real(RP), intent(out) :: SFLX_LH(IA,JA) ! surface flux for latent   heat (area center)   [J/m2/s]
    real(RP), intent(out) :: SFLX_QV(IA,JA) ! surface flux for qv            (area center)   [kg/m2/s]
    real(RP), intent(out) :: U10    (IA,JA) ! velocity u        at 10m height
    real(RP), intent(out) :: V10    (IA,JA) ! velocity v        at 10m height

    real(RP) :: ATM_Uabs(IA,JA) ! absolute velocity at z1 [m/s]
    real(RP) :: R10

    real(RP) :: SFC_PSAT (IA,JA) ! saturatad water vapor pressure [Pa]
    real(RP) :: LHV(IA,JA)

    real(RP) :: SFC_QSAT      ! saturatad water vapor mixing ratio [kg/kg]
    real(RP) :: SFC_QV(IA,JA) ! water vapor mixing ratio [kg/kg]

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / surface flux / simple'

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       ATM_Uabs(i,j) = min( ATMOS_PHY_SF_U_maxM, max( ATMOS_PHY_SF_U_minM, &
            sqrt( ATM_W(i,j)**2 + ATM_U(i,j)**2 + ATM_V(i,j)**2 ) ) ) ! at cell center
    enddo
    enddo

    !-----< momentum >-----

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       SFLX_MW(i,j) = - ATMOS_PHY_SF_Const_Cm * ATM_Uabs(i,j) * SFC_DENS(i,j) * ATM_W(i,j)
       SFLX_MU(i,j) = - ATMOS_PHY_SF_Const_Cm * ATM_Uabs(i,j) * SFC_DENS(i,j) * ATM_U(i,j)
       SFLX_MV(i,j) = - ATMOS_PHY_SF_Const_Cm * ATM_Uabs(i,j) * SFC_DENS(i,j) * ATM_V(i,j)
    enddo
    enddo

    !-----< heat & mass flux >-----
    call SATURATION_psat_all( IA, IS, IE, JA, JS, JE, &
                              SFC_TEMP(:,:), & ! [IN]
                              SFC_PSAT(:,:)  ) ! [OUT]

    call HYDROMETEOR_LHV( &
         IA, IS, IE, JA, JS, JE, &
         SFC_TEMP(:,:), & ! [IN]
         LHV(:,:)       ) ! [OUT]

    !$omp parallel do private( SFC_QSAT )
    do j = JS, JE
    do i = IS, IE
      SFC_QSAT = EPSvap * SFC_PSAT(i,j) / ( SFC_PRES(i,j) - ( 1.0_RP-EPSvap ) * SFC_PSAT(i,j) )
      SFC_QV(i,j) = ( 1.0_RP - ATMOS_PHY_SF_BULK_beta ) * ATM_QV(i,j) + ATMOS_PHY_SF_BULK_beta * SFC_QSAT

      SFLX_SH(i,j) = ATMOS_PHY_SF_Const_Ch * ATM_Uabs(i,j) * SFC_DENS(i,j) * CpDry * ( SFC_TEMP(i,j) - ATM_TEMP(i,j) )
      SFLX_LH(i,j) = ATMOS_PHY_SF_Const_Ce * ATM_Uabs(i,j) * SFC_DENS(i,j) * LHV(i,j) * ( SFC_QV(i,j) - ATM_QV(i,j) )

      SFLX_QV(i,j) = SFLX_LH(i,j) / LHV(i,j)
    enddo
    enddo

    !-----< U10, V10 >-----

    !$omp parallel do &
    !$omp private(R10)
    do j = JS, JE
    do i = IS, IE
       R10 = 10.0_RP / ATM_Z1(i,j)

       U10   (i,j) = R10 * ATM_U(i,j)
       V10   (i,j) = R10 * ATM_V(i,j)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_simple_flux

end module scale_atm_phy_sf_bulk_simple
