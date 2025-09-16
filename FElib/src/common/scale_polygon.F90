!> module common / Polygon
!!
!! @par Description
!!      A module to provide utilities for polygon
!!
!! @par Reference
!!
!! @author Yuta Kawai, Team SCALE
!!
#include "scaleFElib.h"
module scale_polygon
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  public :: polygon_inpoly

contains

  !> Check whether the point is located inside a polygon
!OCL SERIAL
  function polygon_inpoly( pt_x, pt_y, num_node, v_x, v_y ) result(ret)
    implicit none
    real(RP), intent(in) :: pt_x
    real(RP), intent(in) :: pt_y
    integer, intent(in) :: num_node
    real(RP), intent(in) :: v_x(num_node)
    real(RP), intent(in) :: v_y(num_node)
    logical :: ret

    integer :: wn
    integer :: i, ii
    !------------------------------------------

    wn = 0
    do i=1, num_node
      ii = mod(i, num_node) + 1
      if ( v_y(i) <= pt_y .and. pt_y < v_y(ii)) then
        if( pt_x < v_x(i) + (pt_y - v_y(i)) * (v_x(ii) - v_x(i))/(v_y(ii) - v_y(i)) ) then
          wn = wn + 1
        end if
      else if ( v_y(i) > pt_y .and. v_y(ii) <= pt_y ) then
        if( pt_x < v_x(i) + (pt_y - v_y(i)) * (v_x(ii) - v_x(i))/(v_y(ii) - v_y(i)) ) then
          wn = wn - 1
        end if
      end if
    end do

    if (wn == 0) then
      ret = .false.
    else
      ret = .true.
    end if

    return
  end function polygon_inpoly

end module scale_polygon