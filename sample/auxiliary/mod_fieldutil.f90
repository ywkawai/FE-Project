!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_fieldutil
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
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
  !++ Public procedures
  !

  public :: fieldutil_get_profile1d
  public :: fieldutil_get_upwind_pos1d
  public :: fieldutil_get_profile2d

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !


  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------
contains  

!--- 1d -----------------------------------------------------

function fieldutil_get_upwind_pos1d(pos, ADV_VEL, nowtime, dom_min, dom_max) result(upos)
  real(RP), intent(in) :: pos(:)
  real(RP), intent(in) :: ADV_VEL
  real(RP), intent(in) :: nowtime
  real(RP), intent(in) :: dom_min, dom_max
  real(RP) :: upos(size(pos))

  integer :: period
  !-------

  period = ADV_VEL*nowtime/(dom_max - dom_min)
  upos(:) = pos(:) - (ADV_VEL*nowtime - dble(period)*(dom_max - dom_min))
  where (upos < dom_min)
    upos = dom_max + (upos - dom_min)
  end where
end function fieldutil_get_upwind_pos1d

!> Get 1D data whose value is set based on specified pattern. 
!
! Assume that range of domain is -1 <= x <= 1.
! If the profile_name is 'sin', param1 indicates the wavenumber. 
! If the profile_name is 'cosbell', param1 indicates the half of width. 
! If the profile_name is 'top-hat', param1 indicates the half of width. 
!
function fieldutil_get_profile1d(profile_name, x, param1) result(profile)
  use scale_const, only: PI => CONST_PI 
  implicit none

  character(*), intent(in) :: profile_name
  real(RP), intent(in) :: x(:)
  real(RP), intent(in) :: param1

  real(RP) :: profile(size(x))
  !--------------------------------

  profile(:) = 0.0_RP

  select case(profile_name)
  case ('sin')
    profile(:) = sin( param1*PI*x(:) )
  case ('cosbell')
    where( abs(x) <= param1 )
      profile(:) = (1.0_RP + cos(PI*x/param1))*0.5_RP
    end where
  case ('top-hat')
    where( abs(x) <= param1 )
      profile(:) = 1.0_RP
    end where
  end select
  return
end function fieldutil_get_profile1d

!-- 2d--------------------------------------------------------------

!> Get 2D data whose value is set based on specified pattern. 
!
! Assume that range of domain is -1 <= x <= 1.
! If the profile_name is 'sin', param1 and param2 indicate the wavenumbers in x- and y- directions, respectively. 
! If the profile_name is 'cosbell', param1 indicates the half of width. The shape is isotropic about the center of domain. 
! If the profile_name is 'top-hat', param1 and param2 indicate the half of width in x- and y- directions, respectively. 
!
function fieldutil_get_profile2d(profile_name, x, y, param1, param2) result(profile)  
  implicit none
  !--------------------------------

  character(*), intent(in) :: profile_name
  real(RP), intent(in) :: x(:)
  real(RP), intent(in) :: y(:)
  real(RP), intent(in) :: param1
  real(RP), intent(in) :: param2
  real(RP) :: profile(size(x))

  real(RP) :: dist(size(x))
  !------------------------------------------------------------------------

  profile(:) = 0.0_RP
  dist(:) = sqrt(x(:)**2 + y(:)**2)

  select case(profile_name)
  case ('sin')
    profile(:) = sin( param1*PI*x(:) ) * sin( param2*PI*y(:) )
  case ('cosbell')
    where( dist <= param1 )
      profile(:) = (1.0_RP + cos(PI*dist(:)/param1))*0.5_RP
    end where
  case ('top-hat')
    where( abs(x) <= param1 .and. abs(y) <= param2 )
      profile(:) = 1.0_RP
    end where
  end select

  return
end function fieldutil_get_profile2d

!-------------------------------

end module mod_fieldutil
