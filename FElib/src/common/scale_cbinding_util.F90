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
  public :: cbinding_util_logical_f2c
  interface cbinding_util_ArrayPtr_c2f
    module procedure cbinding_util_DArray1DPtr_c2f
    module procedure cbinding_util_DArray2DPtr_c2f
  end interface
  public :: cbinding_util_ArrayPtr_c2f
  public :: cbinding_util_string_c2f

  interface
     function c_strlen(str) bind(C, name="strlen")
       import :: c_ptr, c_size_t
       type(c_ptr), value :: str
       integer(c_size_t) :: c_strlen
     end function c_strlen
  end interface
  
contains
  function cbinding_util_logical_c2f( flag_c ) result( flag_f )
    implicit none
    logical(c_bool), intent(in) :: flag_c
    logical :: flag_f
    !------------------------------
    flag_f = ( flag_c .eqv. .true._c_bool )
    return
  end function cbinding_util_logical_c2f

  function cbinding_util_logical_f2c( flag_f ) result( flag_c )
    implicit none
    logical, intent(in) :: flag_f
    logical(c_bool) :: flag_c
    !------------------------------
    flag_c = flag_f
    return
  end function cbinding_util_logical_f2c

  function cbinding_util_DArray1DPtr_c2f( array_c, Nx ) result( array_f )
    implicit none
    type(c_ptr), value :: array_c
    integer(c_int), value :: Nx
    real(c_double), pointer :: array_f(:)
    !------------------------------
    call c_f_pointer( array_c, array_f, [Nx] )
    return
  end function cbinding_util_Darray1DPtr_c2f

  function cbinding_util_DArray2DPtr_c2f( array_c, Nx, Ny ) result( array_f )
    implicit none
    type(c_ptr), value :: array_c
    integer(c_int), value :: Nx, Ny
    real(c_double), pointer :: array_f(:,:)
    !------------------------------
    call c_f_pointer( array_c, array_f, [Nx, Ny] )
    return
  end function cbinding_util_Darray2DPtr_c2f

  function cbinding_util_string_c2f(cptr) result(fstr)
    implicit none
    type(c_ptr), intent(in) :: cptr
    character(:), allocatable :: fstr

    character(kind=c_char), pointer :: cstr(:)
    integer(c_size_t) :: n
    integer :: i
    !----------------------------

    n = c_strlen(cptr)
    call c_f_pointer(cptr, cstr, [n])

    allocate(character(len=n) :: fstr)
    do i = 1, n
       fstr(i:i) = transfer(cstr(i), ' ')
    end do
    return
  end function cbinding_util_string_c2f
end module scale_cbinding_util