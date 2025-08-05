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
module scale_cbinding_util
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use iso_c_binding
  use scale_precision
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  type, abstract, public :: CBindingBase
  contains
  end type CBindingBase
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !  
  public :: cbinding_util_logical_c2f
contains
  function cbinding_util_logical_c2f( flag_c ) result( flag_f )
    implicit none
    logical(c_bool), intent(in) :: flag_c
    logical :: flag_f
    !------------------------------
    flag_f = ( flag_c .eqv. .true._c_bool )
    return
  end function cbinding_util_logical_c2f

end module scale_cbinding_util