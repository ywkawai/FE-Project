!-------------------------------------------------------------------------------
!> module common / Coriolis parameter 
!!
!! @par Description
!!      Setup coriolis parameter (for regional model)
!!
!!
!! @par Reference
!!
!! @author Yuta Kawai, Team SCALE
!!
#include "scaleFElib.h"
module scale_coriolis_param
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_const, only: &
    PI => CONST_PI,      &
    OHM => CONST_OHM
  
  use scale_io
  use scale_prc

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !  
  public :: get_coriolis_parameter

contains

!> Get Coriolis parameter
!OCL SERIAL
  subroutine get_coriolis_parameter( &
    coriolis,                        &
    COLIORIS_type, Np,               &
    y, f0, beta, y0,                 &
    lat                              )

    implicit none

    integer, intent(in) :: Np                 !< Array size
    real(RP), intent(out) :: coriolis(Np)     !< Array storing Coriolis parameter
    character(*), intent(in) :: COLIORIS_type !< Type of Coriolis parameter [PLANE / SPHERE / NONE]
    real(RP), intent(in), optional :: y(Np)   !< Y-coordinate which is used when COLIORIS_type=PLANE
    real(RP), intent(in), optional :: f0      !< Value of Coriolis parameter at reference latitude which is used when COLIORIS_type=PLANE
    real(RP), intent(in), optional :: beta    !< Rossby parameter at reference latitude which is used when COLIORIS_type=PLANE
    real(RP), intent(in), optional :: y0      !< Distance from reference latitude
    real(RP), intent(in), optional :: lat(Np) !< Latitude

    integer :: i
    !------------------------------------------

    if ( trim(COLIORIS_type) == 'PLANE' ) then
      if ( ( .not. present(f0) ) .or. ( .not. present(beta) )    &
            .or. ( .not. present(y) ) .or. ( .not. present(y0) ) ) then
        LOG_ERROR('get_coriolis_parameter',*) 'If COLIORIS_type is set to PLANE, f0 and beta must be passed. Check!'
        call PRC_abort
      end if
      !$omp parallel do
      do i=1, Np
        coriolis(i) = f0 + beta * ( y(i) - y0 )
      end do
    else if ( trim(COLIORIS_type) == 'SPHERE' ) then
      if ( .not. present(lat) ) then
        LOG_ERROR('get_coriolis_parameter',*) 'If COLIORIS_type is set to SPHERE, lat must be passed. Check!'
        call PRC_abort
      end if
      !$omp parallel do
      do i=1, Np
        coriolis(i) =  2.0_RP * OHM * sin( lat(i) )
      end do
    else if ( trim(COLIORIS_type) == 'NONE' ) then
      !$omp parallel do
      do i=1, Np
        coriolis(i) =  0.0_RP
      end do      
    else
      LOG_ERROR('get_coriolis_parameter',*) 'Unexpected COLIORIS_type is specified. Check! COLIORIS_type=', COLIORIS_type
      call PRC_abort
    end if  

    return
  end subroutine get_coriolis_parameter

end module scale_coriolis_param