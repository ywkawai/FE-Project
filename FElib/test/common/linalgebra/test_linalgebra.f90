#include "scalelib.h"
program test_linalgebra
  use scale_precision
  use scale_prc
  use scale_io
  implicit none

  real(RP), allocatable :: A(:,:)
  real(RP), allocatable :: A_bnd(:,:)
  real(RP), allocatable :: B(:)
  real(RP), allocatable :: X_ans(:)
  integer, allocatable :: ipiv_ans(:)
  integer :: N
  integer :: KL, KU
  !-------------------------------

  write(*,*) "Test LinAlgebra module ..."
  call init()

  call test_SolveLinEq()
  call test_SolveLinEq_bndmat()
  call test_SolveLinEq_vi_linkernel()

  write(*,*) "Test of LinAlgebra module has been succeeded!"
  call final()

contains
  subroutine init()
    implicit none
    !-----------------------------------------

    integer :: comm, myrank, nprocs
    logical :: ismaster  

    integer :: i, j
    !-----------------------------------

    call PRC_MPIstart( comm )
    
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]

    ! setup scale_io                                                                                                                                                                                 
    call IO_setup( "test_linalgebra", "test.conf", allow_noconf = .true. )
  
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
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

  subroutine test_SolveLinEq_vi_linkernel
    use scale_atm_dyn_dgm_hevi_common_linalgebra, only: &
      lin_ludecomp => atm_dyn_dgm_hevi_common_linalgebra_ludecomp, &
      lin_solver => atm_dyn_dgm_hevi_common_linalgebra_solve
    implicit none

    real(RP) :: D(2,N,N)
    real(RP) :: x(2,N)
    integer :: ipiv(2,N)
    integer :: v
    !-------------------------------
    do v=1, 2
      D(v,:,:) = A(:,:)
      x(v,:) = b(:)
    end do
    call lin_ludecomp( D, ipiv, 2, N )
    call lin_solver( D, x, ipiv, 2, N )
    do v=1, 2
      call check_ans( x(v,:), 'test_vi_linkernel' )
    end do
    return
  end subroutine test_SolveLinEq_vi_linkernel

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

  subroutine check_ans( x, lbl )
    implicit none
    real(RP), intent(in) :: x(N)
    character(len=*), intent(in) :: lbl

    real(RP) :: l2error
    real(RP), parameter :: EPS = 4.0E-15_RP
    !----------------------------------------------

    l2error = sqrt( sum( (x(:) - X_ans)**2 ) ) / real(N, kind=RP)
    write(*,*) "* ", lbl, ": x=", x(:)
    write(*,*) "l2error=", l2error

    if ( l2error < EPS ) then
        write(*,*) "=> OK!"
    else
        write(*,*) "=> The error is too large. Check!"
        call PRC_abort
    end if

    return
  end subroutine check_ans

  subroutine check_ipiv( ipiv, lbl )
    implicit none
    integer, intent(in) :: ipiv(N)
    character(len=*), intent(in) :: lbl

    real(RP) :: l2error
    !----------------------------------------------

    write(*,*) "* ", lbl, ": IPIV=", ipiv(:)
    if ( sum( ipiv(:) - ipiv_ans(:) ) == 0 ) then
        write(*,*) "=> OK!"
    else
        write(*,*) "=> IPIV is wrong. Check!"
        call PRC_abort
    end if
    return
  end subroutine check_ipiv

  subroutine final()
    implicit none
    !------------------------------

    call PRC_MPIfinish()

    return
  end subroutine final

end program test_linalgebra