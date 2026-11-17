#include "scalelib.h"
program test_sparsemat
  use scale_precision
  use scale_io
  use scale_const, only: &
    EPS => CONST_EPS
  use scale_prc
  use scale_sparsemat
  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: UNITTEST_assert_close
  implicit none

  integer, parameter :: N = 5 
  real(RP), save :: A(N,N)
  type(SparseMat) :: Acsr
  type(SparseMat) :: Aell
  real(RP), save :: x(N), b(N), b_ans(N)
  integer :: i, j
  !----------------------------------------

  !* Initialize
  call init()

  !* Print original matrix
  write(*,*) "Print: A="
  do i=1, N
    write(*,'(a,5I4,a)') '|', int(A(i,:)), '|'
  end do
  write(*,*) "-----------------------------"

  !* Check the storage 
  write(*,*) "- Print buffer information (CSR format)"
  call Acsr%Print()
  call check_storage( Acsr, "GetVal (CSR format)" )

  write(*,*) "- Print buffer information (ELL format)"
  call Aell%Print()
  call check_storage( Aell, "GetVal (ELL format)" )
  
  !* Check SpMV operation

  x(:) = 1.0_RP
  b_ans(:) = matmul(A,x)

  call sparsemat_matmul(Acsr, x, b)
  call UNITTEST_assert_close( b(:), b_ans(:), "SpMV (CSR format)", rtol=0.0_RP, atol=EPS )

  call sparsemat_matmul(Aell, x, b)
  call UNITTEST_assert_close( b(:), b_ans(:), "SpMV (ELL format)", rtol=0.0_RP, atol=EPS )

  !* Finalize

  call final()

contains  
  !> Check consistency of storage data with the original matrix (only for nonzero elements)
  subroutine check_storage( Amat, lbl )
    implicit none
    type(SparseMat), intent(in) :: Amat
    character(len=*), intent(in) :: lbl

    real(RP) :: A_stored(N,N)
    integer :: i_, j_
    !-----------------------------------

    do j_=1, N
    do i_=1, N
      if ( abs(A(i_,j_)) > EPS ) then
        A_stored(i_,j_) = Amat%GetVal(i_,j_)
      else
        A_stored(i_,j_) = A(i_,j_)
      end if
    end do
    end do
    call UNITTEST_assert_close( A_stored(:,:), A(:,:), lbl, rtol=0.0_RP, atol=EPS )

    return
  end subroutine check_storage

  subroutine init()
    implicit none

    !-----------------------------------

    call UNITTEST_init( "test_sparsemat", conf_fname="test.conf", allow_noconf=.false. )

    !-
    A(1,:) = (/ 1.0_RP, 3.0_RP, 0.0_RP, 0.0_RP, 0.0_RP /)
    A(2,:) = (/ 1.0_RP, 2.0_RP, 5.0_RP, 0.0_RP, 0.0_RP /)
    A(3,:) = (/ 4.0_RP, 1.0_RP, 3.0_RP, 0.0_RP, 0.0_RP /)
    A(4,:) = (/ 0.0_RP, 3.0_RP, 7.0_RP, 4.0_RP, 0.0_RP /)
    A(5,:) = (/ 1.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 5.0_RP /)

    call Acsr%Init( A, EPS=EPS, storage_format='CSR' )
    call Aell%Init( A, EPS=EPS, storage_format='ELL' )

    return
  end subroutine init

  subroutine final()
    implicit none
    !-----------------------------------
    call Acsr%Final()
    call Aell%Final()
    call UNITTEST_final()
    return
  end subroutine final

end program test_sparsemat