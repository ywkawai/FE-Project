!-------------------------------------------------------------------------------
!> Module common / unit test (assertion)
!!
!! @par Description
!!      A module to provide assertion routines for unit tests.
!!      Each assertion increments the counters of assertions and failures.
!!      When an assertion fails, the detail is written to the standard output.
!!      The program is not stopped immediately (unless abort_on_failure is set to .true.);
!!      the total number of failures is checked in UNITTEST_final (scale_unittest).
!!
!!      This module depends only on scale_precision, scale_io, and scale_prc,
!!      so that it can also be used in sample programs without UNITTEST_init.
!!
!! @author Yuta Kawai, Team SCALE
!!
#include "scaleFElib.h"
module scale_unittest_assert
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !

  !> Check that integer values are equal to the answer
  !! (A 1D array can also be compared with a 2D or 3D answer array in array element order)
  interface UNITTEST_assert_equal
    module procedure assert_equal_int_0d
    module procedure assert_equal_int_1d
    module procedure assert_equal_int_2d
    module procedure assert_equal_int_3d
    module procedure assert_equal_int_1d_2d
    module procedure assert_equal_int_1d_3d
  end interface
  public :: UNITTEST_assert_equal

  !> Check that real (or complex) values are close to the answer
  !! (pass if |val - ans| <= atol + rtol * |ans| for all elements)
  interface UNITTEST_assert_close
    module procedure assert_close_real_0d
    module procedure assert_close_real_1d
    module procedure assert_close_real_2d
    module procedure assert_close_real_3d
    module procedure assert_close_cmplx_1d
  end interface
  public :: UNITTEST_assert_close

  public :: UNITTEST_assert_true
  public :: UNITTEST_assert_less
  public :: UNITTEST_assert_greater

  public :: UNITTEST_assert_set_verbose
  public :: UNITTEST_assert_set_abort_on_failure
  public :: UNITTEST_assert_get_num_assert
  public :: UNITTEST_assert_get_num_failure

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  integer, private :: num_assert  = 0
  integer, private :: num_failure = 0
  logical, private :: verbose          = .false.
  logical, private :: abort_on_failure = .false.

  integer, private, parameter :: MAX_NUM_PRINT_MISMATCH = 5  !< Maximum number of mismatched elements to be printed
  integer, private, parameter :: MAX_NUM_PRINT_VALS     = 32 !< Maximum number of elements to be printed in verbose mode

  !> Default relative tolerance used in UNITTEST_assert_close
  real(RP), public, parameter :: UNITTEST_RTOL_DEFAULT = 100.0_RP * epsilon(1.0_RP)

contains

  !-- Settings and counters -----------------------------------------

  subroutine UNITTEST_assert_set_verbose( flag )
    implicit none
    logical, intent(in) :: flag
    !------------------------------------
    verbose = flag
    return
  end subroutine UNITTEST_assert_set_verbose

  subroutine UNITTEST_assert_set_abort_on_failure( flag )
    implicit none
    logical, intent(in) :: flag
    !------------------------------------
    abort_on_failure = flag
    return
  end subroutine UNITTEST_assert_set_abort_on_failure

  function UNITTEST_assert_get_num_assert() result(n)
    implicit none
    integer :: n
    !------------------------------------
    n = num_assert
    return
  end function UNITTEST_assert_get_num_assert

  function UNITTEST_assert_get_num_failure() result(n)
    implicit none
    integer :: n
    !------------------------------------
    n = num_failure
    return
  end function UNITTEST_assert_get_num_failure

  !-- Logical / scalar comparisons ----------------------------------

  !> Check that the condition is true
  subroutine UNITTEST_assert_true( cond, label, msg )
    implicit none
    logical, intent(in) :: cond
    character(len=*), intent(in) :: label
    character(len=*), intent(in), optional :: msg
    !------------------------------------

    if ( cond ) then
      call report_pass( label )
    else
      if ( present(msg) ) then
        call report_failure( label, msg )
      else
        call report_failure( label, 'The condition is false.' )
      end if
      call abort_if_required()
    end if
    return
  end subroutine UNITTEST_assert_true

  !> Check that val < bound (NaN is treated as a failure)
  subroutine UNITTEST_assert_less( val, bound, label )
    implicit none
    real(RP), intent(in) :: val
    real(RP), intent(in) :: bound
    character(len=*), intent(in) :: label

    character(len=H_MID) :: msg
    !------------------------------------

    if ( val < bound ) then
      call report_pass( label, val )
    else
      write(msg,'(a,es15.7,a,es15.7)') 'val=', val, ' is not less than bound=', bound
      call report_failure( label, msg )
      call abort_if_required()
    end if
    return
  end subroutine UNITTEST_assert_less

  !> Check that val > bound (NaN is treated as a failure)
  subroutine UNITTEST_assert_greater( val, bound, label )
    implicit none
    real(RP), intent(in) :: val
    real(RP), intent(in) :: bound
    character(len=*), intent(in) :: label

    character(len=H_MID) :: msg
    !------------------------------------

    if ( val > bound ) then
      call report_pass( label, val )
    else
      write(msg,'(a,es15.7,a,es15.7)') 'val=', val, ' is not greater than bound=', bound
      call report_failure( label, msg )
      call abort_if_required()
    end if
    return
  end subroutine UNITTEST_assert_greater

  !-- Integer comparisons -------------------------------------------

  subroutine assert_equal_int_0d( val, ans, label, idx )
    implicit none
    integer, intent(in) :: val
    integer, intent(in) :: ans
    character(len=*), intent(in) :: label
    integer, intent(in), optional :: idx
    !------------------------------------
    call assert_equal_int_core( (/ val /), (/ ans /), 1, label, idx )
    return
  end subroutine assert_equal_int_0d

  subroutine assert_equal_int_1d( vals, ans, label, idx )
    implicit none
    integer, intent(in) :: vals(:)
    integer, intent(in) :: ans(:)
    character(len=*), intent(in) :: label
    integer, intent(in), optional :: idx
    !------------------------------------
    if ( .not. check_size( size(vals), size(ans), label ) ) return
    call assert_equal_int_core( vals, ans, size(vals), label, idx )
    return
  end subroutine assert_equal_int_1d

  subroutine assert_equal_int_2d( vals, ans, label, idx )
    implicit none
    integer, intent(in) :: vals(:,:)
    integer, intent(in) :: ans(:,:)
    character(len=*), intent(in) :: label
    integer, intent(in), optional :: idx
    !------------------------------------
    if ( .not. check_size( size(vals), size(ans), label ) ) return
    call assert_equal_int_core( reshape(vals, (/ size(vals) /)), reshape(ans, (/ size(ans) /)), &
                                size(vals), label, idx )
    return
  end subroutine assert_equal_int_2d

  subroutine assert_equal_int_3d( vals, ans, label, idx )
    implicit none
    integer, intent(in) :: vals(:,:,:)
    integer, intent(in) :: ans(:,:,:)
    character(len=*), intent(in) :: label
    integer, intent(in), optional :: idx
    !------------------------------------
    if ( .not. check_size( size(vals), size(ans), label ) ) return
    call assert_equal_int_core( reshape(vals, (/ size(vals) /)), reshape(ans, (/ size(ans) /)), &
                                size(vals), label, idx )
    return
  end subroutine assert_equal_int_3d

  !> Compare a 1D array with a 2D answer array (elements are compared in array element order)
  subroutine assert_equal_int_1d_2d( vals, ans, label, idx )
    implicit none
    integer, intent(in) :: vals(:)
    integer, intent(in) :: ans(:,:)
    character(len=*), intent(in) :: label
    integer, intent(in), optional :: idx
    !------------------------------------
    if ( .not. check_size( size(vals), size(ans), label ) ) return
    call assert_equal_int_core( vals, reshape(ans, (/ size(ans) /)), size(vals), label, idx )
    return
  end subroutine assert_equal_int_1d_2d

  !> Compare a 1D array with a 3D answer array (elements are compared in array element order)
  subroutine assert_equal_int_1d_3d( vals, ans, label, idx )
    implicit none
    integer, intent(in) :: vals(:)
    integer, intent(in) :: ans(:,:,:)
    character(len=*), intent(in) :: label
    integer, intent(in), optional :: idx
    !------------------------------------
    if ( .not. check_size( size(vals), size(ans), label ) ) return
    call assert_equal_int_core( vals, reshape(ans, (/ size(ans) /)), size(vals), label, idx )
    return
  end subroutine assert_equal_int_1d_3d

  subroutine assert_equal_int_core( vals, ans, n, label, idx )
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: vals(n)
    integer, intent(in) :: ans(n)
    character(len=*), intent(in) :: label
    integer, intent(in), optional :: idx

    integer :: i, nmis, nprint
    character(len=H_LONG) :: label_
    !------------------------------------

    label_ = make_label( label, idx )
    nmis = count( vals(:) /= ans(:) )

    if ( nmis == 0 ) then
      call report_pass( label_ )
      if ( verbose ) write(*,*) '   val=', vals(1:min(n,MAX_NUM_PRINT_VALS))
    else
      call report_failure( label_, 'Integer values are inconsistent with the answer.' )
      write(*,'(a,i0,a,i0)') '   number of mismatched elements: ', nmis, ' / ', n
      nprint = 0
      do i=1, n
        if ( vals(i) /= ans(i) ) then
          write(*,'(a,i0,a,i0,a,i0)') '   i=', i, ': val=', vals(i), ', ans=', ans(i)
          nprint = nprint + 1
          if ( nprint >= MAX_NUM_PRINT_MISMATCH ) exit
        end if
      end do
      call abort_if_required()
    end if
    return
  end subroutine assert_equal_int_core

  !-- Real / complex comparisons ------------------------------------

  subroutine assert_close_real_0d( val, ans, label, rtol, atol, idx )
    implicit none
    real(RP), intent(in) :: val
    real(RP), intent(in) :: ans
    character(len=*), intent(in) :: label
    real(RP), intent(in), optional :: rtol
    real(RP), intent(in), optional :: atol
    integer, intent(in), optional :: idx
    !------------------------------------
    call assert_close_real_core( (/ val /), (/ ans /), 1, label, rtol, atol, idx )
    return
  end subroutine assert_close_real_0d

  subroutine assert_close_real_1d( vals, ans, label, rtol, atol, idx )
    implicit none
    real(RP), intent(in) :: vals(:)
    real(RP), intent(in) :: ans(:)
    character(len=*), intent(in) :: label
    real(RP), intent(in), optional :: rtol
    real(RP), intent(in), optional :: atol
    integer, intent(in), optional :: idx
    !------------------------------------
    if ( .not. check_size( size(vals), size(ans), label ) ) return
    call assert_close_real_core( vals, ans, size(vals), label, rtol, atol, idx )
    return
  end subroutine assert_close_real_1d

  subroutine assert_close_real_2d( vals, ans, label, rtol, atol, idx )
    implicit none
    real(RP), intent(in) :: vals(:,:)
    real(RP), intent(in) :: ans(:,:)
    character(len=*), intent(in) :: label
    real(RP), intent(in), optional :: rtol
    real(RP), intent(in), optional :: atol
    integer, intent(in), optional :: idx
    !------------------------------------
    if ( .not. check_size( size(vals), size(ans), label ) ) return
    call assert_close_real_core( reshape(vals, (/ size(vals) /)), reshape(ans, (/ size(ans) /)), &
                                 size(vals), label, rtol, atol, idx )
    return
  end subroutine assert_close_real_2d

  subroutine assert_close_real_3d( vals, ans, label, rtol, atol, idx )
    implicit none
    real(RP), intent(in) :: vals(:,:,:)
    real(RP), intent(in) :: ans(:,:,:)
    character(len=*), intent(in) :: label
    real(RP), intent(in), optional :: rtol
    real(RP), intent(in), optional :: atol
    integer, intent(in), optional :: idx
    !------------------------------------
    if ( .not. check_size( size(vals), size(ans), label ) ) return
    call assert_close_real_core( reshape(vals, (/ size(vals) /)), reshape(ans, (/ size(ans) /)), &
                                 size(vals), label, rtol, atol, idx )
    return
  end subroutine assert_close_real_3d

  subroutine assert_close_cmplx_1d( vals, ans, label, rtol, atol, idx )
    implicit none
    complex(RP), intent(in) :: vals(:)
    complex(RP), intent(in) :: ans(:)
    character(len=*), intent(in) :: label
    real(RP), intent(in), optional :: rtol
    real(RP), intent(in), optional :: atol
    integer, intent(in), optional :: idx

    real(RP) :: diff(size(vals)), ans_abs(size(vals))
    !------------------------------------
    if ( .not. check_size( size(vals), size(ans), label ) ) return
    diff(:)    = abs( vals(:) - ans(:) )
    ans_abs(:) = abs( ans(:) )
    call assert_close_core( diff, ans_abs, size(vals), label, rtol, atol, idx )
    return
  end subroutine assert_close_cmplx_1d

  subroutine assert_close_real_core( vals, ans, n, label, rtol, atol, idx )
    implicit none
    integer, intent(in) :: n
    real(RP), intent(in) :: vals(n)
    real(RP), intent(in) :: ans(n)
    character(len=*), intent(in) :: label
    real(RP), intent(in), optional :: rtol
    real(RP), intent(in), optional :: atol
    integer, intent(in), optional :: idx
    !------------------------------------

    call assert_close_core( abs(vals(:) - ans(:)), abs(ans(:)), n, label, rtol, atol, idx, &
                            vals, ans )
    return
  end subroutine assert_close_real_core

  !> Core routine of UNITTEST_assert_close
  !! Pass if diff <= atol + rtol * ans_abs for all elements (NaN is treated as a failure)
  subroutine assert_close_core( diff, ans_abs, n, label, rtol, atol, idx, &
                                vals, ans )
    implicit none
    integer, intent(in) :: n
    real(RP), intent(in) :: diff(n)
    real(RP), intent(in) :: ans_abs(n)
    character(len=*), intent(in) :: label
    real(RP), intent(in), optional :: rtol
    real(RP), intent(in), optional :: atol
    integer, intent(in), optional :: idx
    real(RP), intent(in), optional :: vals(n)
    real(RP), intent(in), optional :: ans(n)

    real(RP) :: rtol_, atol_
    logical :: is_ok(n)
    integer :: i, nmis, nprint
    character(len=H_LONG) :: label_
    !------------------------------------

    rtol_ = UNITTEST_RTOL_DEFAULT
    if ( present(rtol) ) rtol_ = rtol
    atol_ = 0.0_RP
    if ( present(atol) ) atol_ = atol

    label_ = make_label( label, idx )
    do i=1, n
      is_ok(i) = ( diff(i) <= atol_ + rtol_ * ans_abs(i) ) ! NaN gives .false.
    end do
    nmis = count( .not. is_ok(:) )

    if ( nmis == 0 ) then
      call report_pass( label_, maxval(diff(:)) )
      if ( verbose .and. present(vals) ) write(*,*) '   val=', vals(1:min(n,MAX_NUM_PRINT_VALS))
    else
      call report_failure( label_, 'Values are not close to the answer.' )
      write(*,'(a,i0,a,i0,a,es12.4,a,es12.4,a)') '   number of mismatched elements: ', nmis, ' / ', n, &
        ' (rtol=', rtol_, ', atol=', atol_, ')'
      nprint = 0
      do i=1, n
        if ( .not. is_ok(i) ) then
          if ( present(vals) .and. present(ans) ) then
            write(*,'(a,i0,a,es24.16,a,es24.16,a,es12.4)') '   i=', i, ': val=', vals(i), ', ans=', ans(i), &
              ', |diff|=', diff(i)
          else
            write(*,'(a,i0,a,es12.4,a,es12.4)') '   i=', i, ': |diff|=', diff(i), ', |ans|=', ans_abs(i)
          end if
          nprint = nprint + 1
          if ( nprint >= MAX_NUM_PRINT_MISMATCH ) exit
        end if
      end do
      call abort_if_required()
    end if
    return
  end subroutine assert_close_core

  !-- Private utilities ---------------------------------------------

  function check_size( n_vals, n_ans, label ) result(is_ok)
    implicit none
    integer, intent(in) :: n_vals
    integer, intent(in) :: n_ans
    character(len=*), intent(in) :: label
    logical :: is_ok

    character(len=H_MID) :: msg
    !------------------------------------
    is_ok = ( n_vals == n_ans )
    if ( .not. is_ok ) then
      write(msg,'(a,i0,a,i0)') 'Array sizes are inconsistent: size(val)=', n_vals, ', size(ans)=', n_ans
      call report_failure( label, msg )
      call abort_if_required()
    end if
    return
  end function check_size

  function make_label( label, idx ) result(label_)
    implicit none
    character(len=*), intent(in) :: label
    integer, intent(in), optional :: idx
    character(len=H_LONG) :: label_
    !------------------------------------
    if ( present(idx) ) then
      write(label_,'(a,a,i0,a)') trim(label), ' (idx=', idx, ')'
    else
      label_ = label
    end if
    return
  end function make_label

  subroutine report_pass( label, err )
    implicit none
    character(len=*), intent(in) :: label
    real(RP), intent(in), optional :: err
    !------------------------------------
    num_assert = num_assert + 1
    if ( verbose ) then
      if ( present(err) ) then
        write(*,'(a,a,a,es12.4)') '[PASS] ', trim(label), ' : ', err
      else
        write(*,'(a,a)') '[PASS] ', trim(label)
      end if
    end if
    return
  end subroutine report_pass

  subroutine report_failure( label, msg )
    implicit none
    character(len=*), intent(in) :: label
    character(len=*), intent(in) :: msg
    !------------------------------------
    num_assert  = num_assert + 1
    num_failure = num_failure + 1
    write(*,'(a,a,a,a)') '[FAIL] ', trim(label), ' : ', trim(msg)
    if ( IO_L .and. IO_FID_LOG > 0 .and. IO_FID_LOG /= IO_FID_STDOUT ) then
      write(IO_FID_LOG,'(a,a,a,a)') '[FAIL] ', trim(label), ' : ', trim(msg)
    end if
    return
  end subroutine report_failure

  subroutine abort_if_required()
    implicit none
    !------------------------------------
    if ( abort_on_failure ) then
      write(*,*) 'Abort since abort_on_failure is .true.'
      call PRC_abort
    end if
    return
  end subroutine abort_if_required

end module scale_unittest_assert
