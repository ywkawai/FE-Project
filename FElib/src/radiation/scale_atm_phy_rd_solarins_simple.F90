!> module FElib / Atmosphere / Physics radiation / Solar insolation / Simple gray-radiation scheme
!!
!! @par Description
!!      A module for idealized solar insolation scheme
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_rd_solarins_simple
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_const, only: &
    UNDEF => CONST_UNDEF, &
    PI => CONST_PI

  !-----------------------------------------------------------------------------
  implicit none
  private  
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  public :: atm_phy_rd_solarins_simple_setup
  public :: atm_phy_rd_solarins_simple_update
  public :: atm_phy_rd_solarins_simple_get

  !-----------------------------------------------------------------------------
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
  integer :: SOLARINS_TYPE_ID

  integer, parameter :: SOLARINS_SIMPLE_TYPE_ID_CONST    = 1
  integer, parameter :: SOLARINS_SIMPLE_TYPE_ANNUAL_MEAN = 2

  real(RP) :: CONST_FLUX            = 340.0_RP
  integer :: ORBIT_REFERENCE_YEAR   = 2000
  integer :: ANNUAL_MEAN_YEAR       = 2000
  integer :: ANNUAL_SAMPLES_PER_DAY = 1
  logical :: CACHE_ANNUAL_MEAN      = .true.

  integer, parameter :: DATE_YEAR   = 1
  integer, parameter :: DATE_MONTH  = 2
  integer, parameter :: DATE_DAY    = 3
  integer, parameter :: DATE_HOUR   = 4
  integer, parameter :: DATE_MINUTE = 5
  integer, parameter :: DATE_SECOND = 6
  integer, parameter :: DATE_SIZE = 6

contains
  subroutine atm_phy_rd_solarins_simple_setup()
    use scale_atmos_solarins, only: &
        ATMOS_SOLARINS_setup    
    implicit none
    character(len=H_SHORT) :: SOLARINS_TYPE = 'ANNUAL_MEAN'
    namelist / PARAM_ATMOS_SOLARINS_SIMPLE / &
        SOLARINS_TYPE,                    &
        CONST_FLUX,                       &
        ORBIT_REFERENCE_YEAR,             &
        ANNUAL_MEAN_YEAR,                 &
        ANNUAL_SAMPLES_PER_DAY,           &
        CACHE_ANNUAL_MEAN

    integer :: ierr
    !------------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_SOLARINS_SIMPLE_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_SOLARINS_SIMPLE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_RD_SOLARINS_SIMPLE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_RD_SOLARINS_SIMPLE_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_SOLARINS_SIMPLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_SOLARINS_SIMPLE)

    select case (SOLARINS_TYPE)
    case ('CONST')
       SOLARINS_TYPE_ID = SOLARINS_SIMPLE_TYPE_ID_CONST
    case ('ANNUAL_MEAN')
       SOLARINS_TYPE_ID = SOLARINS_SIMPLE_TYPE_ANNUAL_MEAN
    case default
       LOG_ERROR("ATMOS_PHY_RD_SOLARINS_SIMPLE_setup",*) 'Not appropriate SOLARINS_TYPE. Check!'
       call PRC_abort
    end select

    !- Initializes orbital parameters, the vernal-equinox reference date, and related internal state.
    !  Here, longitude and latitude are irrelevant because this module uses only the orbital state returned by ecliptic_longitude.
    call ATMOS_SOLARINS_setup( 0.0_RP, 0.0_RP, ORBIT_REFERENCE_YEAR )

    return
  end subroutine atm_phy_rd_solarins_simple_setup

!OCL SERIAL
  subroutine atm_phy_rd_solarins_simple_update()
  end subroutine atm_phy_rd_solarins_simple_update

!OCL SERIAL
  subroutine atm_phy_rd_solarins_simple_get( solins, cosSZA, &
    lat, Np )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: solins(Np)
    real(RP), intent(out) :: cosSZA(Np)
    real(RP), intent(in) :: lat(Np)
    !------------------------------------------------------------------------------

    select case (SOLARINS_TYPE_ID)
    case (SOLARINS_SIMPLE_TYPE_ID_CONST)
      solins(:) = CONST_FLUX
      cosSZA(:) = UNDEF
    case (SOLARINS_SIMPLE_TYPE_ANNUAL_MEAN)
      call annual_mean_insol( solins, & ! (out)
        lat, Np )                       ! (in)
      cosSZA(:) = UNDEF
    end select
    return
  end subroutine atm_phy_rd_solarins_simple_get

  !- Private subroutines ------------------------------------
  subroutine annual_mean_insol( solins, &
    lat, Np )
    use scale_atmos_solarins, only: &
        ATMOS_SOLARINS_ecliptic_longitude
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: solins(Np)
    real(RP), intent(in) :: lat(Np)

    integer :: date(DATE_SIZE)    
    integer :: month, day
    integer :: ndays_month
    integer :: nsample

    real(RP) :: re_factor
    real(RP) :: sin_decl
    real(RP) :: cos_decl
    real(RP) :: hour_angle
    integer, parameter :: OFFSET_YEAR = 0

    integer :: i

    real(RP) :: sample_second

    real(RP) :: daily_solins(Np)
    real(RP) :: sum_flux(Np)
    !------------------------------------------------------------------------------

    date(DATE_YEAR)   = ANNUAL_MEAN_YEAR
    nsample = 1
    sum_flux(:) = 0.0_RP

    do month=1, 12
      ndays_month = days_in_month( ANNUAL_MEAN_YEAR, month )
      date(DATE_MONTH) = month
      do day=1, ndays_month
        date(DATE_DAY) = day
        do i=1, ANNUAL_SAMPLES_PER_DAY
          ! Midpoint sampling within each day.
          sample_second = 86400.0_RP                                       &
               * ( real(i,RP) - 0.5_RP ) / real(ANNUAL_SAMPLES_PER_DAY,RP)
        
          call seconds_to_hms( sample_second,                        & ! (in)
               date(DATE_HOUR), date(DATE_MINUTE), date(DATE_SECOND) ) ! (out)
        
          call ATMOS_SOLARINS_ecliptic_longitude( &
            re_factor, sin_decl, cos_decl, hour_angle, & ! (out)
            date, OFFSET_YEAR                          ) ! (in)
        
          call daily_mean_insol( daily_solins,     & ! (out)
            lat, re_factor, sin_decl, cos_decl, Np ) ! (in)
    
          sum_flux(:) = sum_flux(:) + daily_solins(:)
          nsample = nsample + 1
        end do
      end do
    end do

    solins(:) = sum_flux(:) / real(nsample, RP)
    return
  end subroutine annual_mean_insol

!OCL SERIAL
  subroutine daily_mean_insol( solins, &
    lat, re_factor, sin_decl, cos_decl, Np )
    use scale_atmos_solarins, only: &
        ATMOS_SOLARINS_CONSTANT
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: solins(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: re_factor
    real(RP), intent(in) :: sin_decl
    real(RP), intent(in) :: cos_decl

    real(RP) :: sin_lat
    real(RP) :: cos_lat
    real(RP) :: cos_h0
    real(RP) :: h0

    real(RP), parameter :: POLE_EPS = 100.0_RP * epsilon(1.0_RP)

    integer :: i
    !------------------------------------------------------------------------------

    !$omp parallel do private(sin_lat, cos_lat, cos_h0, h0)
    do i=1, Np
      sin_lat = sin(lat(i))
      cos_lat = cos(lat(i))
      if ( abs(cos_lat) <= POLE_EPS ) then
        if ( sin_lat * sin_decl > 0.0_RP ) then
          h0 = PI
        else
          h0 = 0.0_RP
        end if
      else if ( abs(cos_decl) <= POLE_EPS ) then
        if ( sin_lat * sin_decl > 0.0_RP ) then
          h0 = PI
        else
          h0 = 0.0_RP
        end if
      else
        cos_h0 = - ( sin_lat * sin_decl ) / ( cos_lat * cos_decl )
        if ( cos_h0 >= 1.0_RP ) then ! Polar night
          h0 = 0.0_RP
        else if ( cos_h0 <= -1.0_RP ) then ! Polar day
          h0 = PI
        else
          h0 = acos(cos_h0)
        end if
      end if

      solins(i) = ATMOS_SOLARINS_CONSTANT * re_factor / PI             &
            * ( h0 * sin_lat * sin_decl + sin(h0) * cos_lat * cos_decl )
      
      solins(i) = max( solins(i), 0.0_RP )
    end do
    return
  end subroutine daily_mean_insol

  !> Number of days in a month
  pure function days_in_month( YEAR, MONTH ) result(NDAYS)
    implicit none
    integer, intent(in) :: YEAR
    integer, intent(in) :: MONTH

    integer :: NDAYS
    integer, parameter :: NDAYS_NORMAL(12) = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
    !---------------------------------------------------------------------------

    NDAYS = NDAYS_NORMAL(MONTH)
    if ( MONTH == 2 .and. is_leap_year(YEAR) ) NDAYS = 29
    return
  end function days_in_month

  !> Gregorian leap-year test
!OCL SERIAL
  pure function is_leap_year( YEAR ) result(IS_LEAP)
    implicit none
    integer, intent(in) :: YEAR
    logical :: IS_LEAP
    !---------------------------------------------------------------------------

    IS_LEAP = mod(YEAR,4) == 0

    if ( mod(YEAR,100) == 0 ) IS_LEAP = .false.
    if ( mod(YEAR,400) == 0 ) IS_LEAP = .true.
    return
  end function is_leap_year  

  !> Convert seconds from start of day to hour/minute/second
!OCL SERIAL
  pure subroutine seconds_to_hms( DAY_SECOND, &
    HOUR, MINUTE, SECOND                      )
    implicit none
    real(RP), intent(in)  :: DAY_SECOND
    integer, intent(out) :: HOUR
    integer, intent(out) :: MINUTE
    integer, intent(out) :: SECOND

    integer :: total_second
    !---------------------------------------------------------------------------

    total_second = int(DAY_SECOND)
    total_second = max(0, min(total_second, 86399))

    HOUR = total_second / 3600
    total_second = total_second - HOUR * 3600

    MINUTE = total_second / 60
    SECOND = total_second - MINUTE * 60
    return
  end subroutine seconds_to_hms  

end module scale_atm_phy_rd_solarins_simple
