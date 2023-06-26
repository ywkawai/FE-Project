!-------------------------------------------------------------------------------
!> module Utility for mktopo
!!
!! @par Description
!!          subroutines useful to prepare topography data
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_mktopo_util
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_const, only: &
    PI => CONST_PI

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: mktopoutil_barocwave_global_JW2006_calc_topo

contains

  !-----------------------------------------------------------------------------
  !> Calculate a topography for a global test case of idealized baroclinic wave in Jablonowski and Williamson (2006)
  !!
  !! @par Reference
  !! - Jablonowski, C., & Williamson, D. L., 2006: A baroclinic instability test case for atmospheric model dynamical cores. Quarterly Journal of the Royal Meteorological Society: A journal of the atmospheric sciences, applied meteorology and physical oceanography, 132(621C), 2943-2975.
  !!
!OCL SERIAL
  subroutine mktopoutil_barocwave_global_JW2006_calc_topo( topo, &
    U0, ETA0, lat, Np )

    use scale_const, only: &
      OHM => CONST_OHM,      &
      Grav => CONST_GRAV,    &
      RPlanet => CONST_RADIUS

    implicit none

    integer, intent(in) :: Np
    real(RP), intent(out) :: topo(Np)
    real(RP), intent(in) :: U0        !< The value of η at a reference level (position of the jet)
    real(RP), intent(in) :: ETA0      !< The parameter associated with zonal jet maximum amplitude  [m/s]
    real(RP), intent(in) :: lat(Np)   !< latitude [radian]

    real(RP) :: sin_lat(Np)
    real(RP) :: sin_lat_pow_6(Np)
    real(RP) :: cos_lat(Np)

    real(RP ) :: tmp
    !-------------------------------------------

    sin_lat(:) = sin(lat(:))
    sin_lat_pow_6(:) = sin_lat(:)**6
    cos_lat(:) = cos(lat(:))

    tmp = cos( (1.0_RP - ETA0) * 0.5_RP * PI )
    tmp = U0 * tmp * sqrt(tmp)

    ! Calc horizontal variation of geopotential height
    topo(:) = tmp * & 
      (   tmp * ( - 2.0_RP * sin_lat_pow_6(:) * ( cos_lat(:)**2 + 1.0_RP / 3.0_RP ) + 10.0_RP / 63.0_RP )          &
        + RPlanet * OHM * ( 8.0_RP / 5.0_RP * cos_lat(:)**3 * ( sin_lat(:)**2 + 2.0_RP / 3.0_RP ) - 0.25_RP * PI ) &
      ) / Grav

    return
  end subroutine mktopoutil_barocwave_global_JW2006_calc_topo

end module mod_mktopo_util