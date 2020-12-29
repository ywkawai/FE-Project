#include "scalelib.h"
program test_sparsemat
  use scale_precision
  use scale_io
  use scale_const, only: &
    EPS => CONST_EPS
  use scale_prc
  use scale_sparsemat
  implicit none

  integer, parameter :: N = 5 
  real(RP) :: A(N,N)
  type(SparseMat) :: Acsr
  type(SparseMat) :: Aell
  real(RP) :: x(N), b(N), b_ans(N)
  integer :: i, j
  !----------------------------------------

  !* Initialize
  write(*,*) "- Start test_sparsemat .."
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
  write(*,*) "- Check consistency of storage data (CSR format)..."
  do j=1, N
  do i=1, N
    if ( abs(A(i,j)) > EPS ) then
      if ( abs(Acsr%GetVal(i,j) - A(i,j)) > EPS ) then
        write(*,'(a,2i3,a,2f12.5)') "GetVal: Buffer data is inconsitent to original matrix data! i,j=", i, j, &
          ": A(i,j), storage data=", A(i,j), Acsr%GetVal(i,j)
        call PRC_abort
      end if
    end if
  end do
  end do
  write(*,*) "OK!"

  write(*,*) "- Print buffer information (ELL format)"
  call Aell%Print()
  write(*,*) "- Check consistency of storage data (ELL format)..."
  do j=1, N
  do i=1, N
    if ( abs(A(i,j)) > EPS ) then
      if ( abs(Aell%GetVal(i,j) - A(i,j)) > EPS ) then
        write(*,'(a,2i3,a,2f12.5)') "GetVal: Buffer data is inconsitent to original matrix data! i,j=", i, j, &
          ": A(i,j), storage data'=", A(i,j), Aell%GetVal(i,j)
        call PRC_abort
      end if
    end if
  end do
  end do
  write(*,*) "OK!"
  
  !* Check SpMV operation

  write(*,*) "- Check the result of SpMV operation.."

  x(:) = 1.0_RP

  b_ans(:) = matmul(A,x)
  write(*,'(a,5i5)') "* Answer b=", int(b_ans(:))

  call sparsemat_matmul(Acsr, x, b)
  write(*,'(a,5i5)') "* CSR format: b=", int(b(:))
  do i=1, N
    if ( abs(b(i) - b_ans(i)) > EPS ) then
      write(*,'(a,i3,a,2f12.5)') "SpMV (CSR format): the result is invalid! i=", i, &
        ": b, b_answer=", b(i), b_ans(i)
    call PRC_abort
    end if  
  end do
  write(*,*) "OK!"

  call sparsemat_matmul(Aell, x, b)
  write(*,'(a,5i5)') "* ELL format: b=", int(b(:))
  do i=1, N
    if ( abs(b(i) - b_ans(i)) > EPS ) then
      write(*,'(a,i3,a,2f12.5)') "SpMV (ELL format): the result is invalid! i=", i, &
        ": b, b_answer=", b(i), b_ans(i)
    call PRC_abort
    end if  
  end do
  write(*,*) "OK!"


  !* Finalize

  call final()
  write(*,*) 'test_sparsemat has been succeeded!'

contains  
  subroutine init()
    implicit none

    integer :: comm, myrank, nprocs
    logical :: ismaster  
    !-----------------------------------

    call PRC_MPIstart( comm )
    
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]

    ! setup scale_io
    call IO_setup( "test_sparsemat", "test.conf", allow_noconf = .false. )

    ! setup log
    call IO_LOG_setup( myrank, ismaster )   

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
    call PRC_MPIfinish()
    return
  end subroutine final

end program test_sparsemat