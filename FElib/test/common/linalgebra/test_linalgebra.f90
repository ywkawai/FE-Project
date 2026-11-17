#include "scalelib.h"
program test_linalgebra
  use scale_precision
  use scale_prc
  use scale_io
  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: &
    UNITTEST_assert_equal, UNITTEST_assert_less
  implicit none

  real(RP), allocatable :: A(:,:)
  real(RP), allocatable :: A_batch(:,:,:)
  real(RP), allocatable :: A_bnd(:,:)
  real(RP), allocatable :: B(:)
  real(RP), allocatable :: B_batch(:,:)
  real(RP), allocatable :: X_ans(:)
  real(RP), allocatable :: X_batch_ans(:,:)
  integer, allocatable :: ipiv_ans(:)
  integer :: N
  integer :: KL, KU
  !-------------------------------

  write(*,*) "Test LinAlgebra module ..."
  call init()

  call test_SolveLinEq()
  call test_SolveLinEq_bndmat()

  call final()

contains
  subroutine init()
    implicit none
    !-----------------------------------------

    integer :: i, j
    !-----------------------------------

    call UNITTEST_init( "test_linalgebra", conf_fname="test.conf" )
    !----

    KL = 1; KU = 2
    N = 4
    allocate( A(N,N), B(N) )
    allocate( A_bnd(2*KL+KU+1,N) )

    A(1,:) = (/ -0.23_RP, 2.54_RP, -3.66_RP, 0.0_RP /)
    A(2,:) = (/ -6.98_RP, 2.46_RP, -2.73_RP, -2.13_RP /)
    A(3,:) = (/ 0.0_RP, 2.56_RP, 2.46_RP, 4.07_RP /)
    A(4,:) = (/ 0.0_RP, 0.0_RP, -4.78_RP, -3.82_RP /)
    B(:) = (/ 4.42_RP, 27.13_RP, -6.14_RP, 10.50_RP /)

    do i=1, N
    do j=max(i-KL,1), min(i+KU,N)
      A_bnd(KL+KU+1+i-j,j) = A(i,j)
    end do
    end do

    allocate( X_ans(N) )
    X_ans(:) = (/ -2.0_RP, 3.0_RP, 1.0_RP, -4.0_RP /)

    allocate( ipiv_ans(N) )
    ipiv_ans(:) = (/ 2, 3, 3, 4 /)

    !---------------------------
    allocate( A_batch(3,N,N), B_batch(3,N) )
    allocate( X_batch_ans(3,N) )

    A_batch(1,:,:) = A(:,:)
    B_batch(1,:) = B(:)

    A_batch(2,1,:) = (/ 2.0_RP, 1.0_RP, 3.0_RP, 4.0_RP /)
    A_batch(2,2,:) = (/ 1.0_RP, 2.0_RP, 1.0_RP, 3.0_RP /)
    A_batch(2,3,:) = (/ 3.0_RP, 1.0_RP, 2.0_RP, 1.0_RP /)
    A_batch(2,4,:) = (/ 4.0_RP, 3.0_RP, 2.0_RP, 1.0_RP /)
    B_batch(2,:) = (/ 29.0_RP, 20.0_RP, 15.0_RP, 20.0_RP /)

    A_batch(3,1,:) = (/ 1.0_RP, 2.0_RP, 3.0_RP, 4.0_RP /)
    A_batch(3,2,:) = (/ 0.0_RP, 1.0_RP, 0.0_RP, 2.0_RP /)
    A_batch(3,3,:) = (/ 3.0_RP, 0.0_RP, 2.0_RP, 1.0_RP /)
    A_batch(3,4,:) = (/ 1.0_RP, 1.0_RP, 1.0_RP, 1.0_RP /)
    B_batch(3,:) = (/ 13.0_RP, 1.0_RP, 12.0_RP, 6.0_RP /)

    X_batch_ans(1,:) = X_ans(:)
    X_batch_ans(2,:) = (/ 1.0_RP, 2.0_RP, 3.0_RP, 4.0_RP /)
    X_batch_ans(3,:) = (/ 2.0_RP, 1.0_RP, 3.0_RP, 0.0_RP /)

    return
  end subroutine init

  subroutine test_SolveLinEq
    use scale_linalgebra, only: linalgebra_SolveLinEq
    implicit none

    real(RP) :: x(N)
    !-------------------------------

    call linalgebra_SolveLinEq( A, b, x )
    call check_ans( x, 'test_SolveLinEq' )
    return
  end subroutine test_SolveLinEq

  subroutine test_SolveLinEq_bndmat
    use scale_linalgebra, only: linalgebra_SolveLinEq_BndMat
    implicit none

    real(RP) :: A_bnd_(2*KL+KU+1,N)
    real(RP) :: x(N)
    integer :: ipiv(N)
    !-------------------------------

    !--  The case where lapack is used
    A_bnd_(:,:) = A_bnd(:,:)
    x(:) = b(:)
    call linalgebra_SolveLinEq_BndMat( A_bnd_, x, ipiv, N, KL, KU, 1, use_lapack=.true. )
    call check_ans( x, 'test_SolveLinEq_bndmat (lapack is used)' )
    call check_ipiv( ipiv, 'IPIV' )

    !-- The case where lapack is not used
    A_bnd_(:,:) = A_bnd(:,:)
    x(:) = b(:)
    call linalgebra_SolveLinEq_BndMat( A_bnd_, x, ipiv, N, KL, KU, 1, use_lapack=.false. )
    call check_ans( x, 'test_SolveLinEq_bndmat (lapack is not used)' )
    call check_ipiv( ipiv, 'IPIV' )

    return
  end subroutine test_SolveLinEq_bndmat

  subroutine check_ans( x, lbl, x_ans_ )
    implicit none
    real(RP), intent(in) :: x(N)
    character(len=*), intent(in) :: lbl
    real(RP), intent(in), optional :: x_ans_(N)

    real(RP) :: l2error
    real(RP), parameter :: EPS = 4.0E-15_RP
    !----------------------------------------------

    if ( present(x_ans_) ) then
      l2error = sqrt( sum( (x(:) - X_ans_)**2 ) ) / real(N, kind=RP)
    else
      l2error = sqrt( sum( (x(:) - X_ans)**2 ) ) / real(N, kind=RP)
    end if
    call UNITTEST_assert_less( l2error, EPS, lbl )

    return
  end subroutine check_ans

  subroutine check_ipiv( ipiv, lbl )
    implicit none
    integer, intent(in) :: ipiv(N)
    character(len=*), intent(in) :: lbl

    !----------------------------------------------

    call UNITTEST_assert_equal( ipiv(:), ipiv_ans(:), lbl )
    return
  end subroutine check_ipiv

  subroutine final()
    implicit none
    !------------------------------

    call UNITTEST_final()

    return
  end subroutine final

end program test_linalgebra