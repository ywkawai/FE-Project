!> module common / Linear algebra
!!
!! @par Description
!!      A module to provide utilities for linear algebra
!!
!! @par Reference
!!
!! @author Team SCALE
!!
#include "scaleFElib.h"
module scale_linalgebra
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_sparsemat, only: sparsemat

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !  

  public :: linalgebra_LU
  public :: linalgebra_inv
  interface linalgebra_SolveLinEq
    module procedure linalgebra_SolveLinEq_b1D
    module procedure linalgebra_SolveLinEq_b2D
  end interface
  public :: linalgebra_SolveLinEq
  public :: linalgebra_SolveLinEq_GMRES


  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  
  !-----------------------------------------------------------------------------

  ! Private procedure
  !

  ! Private variable
  !

contains

!OCL SERIAL
  function linalgebra_inv(A) result(Ainv)

    real(RP), intent(in) :: A(:,:)
    real(RP) :: Ainv(size(A,1),size(A,2))

    real(RP) :: work(size(A,1))
    integer :: ipiv(size(A,1))
    integer :: n, info

    !--------------------------------------------------------------------------- 

    Ainv(:,:) = A
    n = size(A,1)

    call DGETRF(n, n, Ainv, n, ipiv, info)
    if (info /=0 ) then
      LOG_ERROR("linalgebra_inv",*) "Matrix is singular"
      call PRC_abort
    end if

    call DGETRI(n, Ainv, n, ipiv, work, n, info)
    if (info /=0 ) then
      LOG_ERROR("linalgebra_inv",*) "Matrix inversion is failed"
      call PRC_abort
    end if

    return
  end function linalgebra_inv

!OCL SERIAL
  subroutine linalgebra_SolveLinEq_b1D(A, b, x)

    real(RP), intent(in) :: A(:,:)
    real(RP), intent(in) :: b(:)
    real(RP), intent(out) :: x(size(b))

    real(RP) :: A_lu(size(A,1),size(A,2))
    integer :: ipiv(size(A,1))
    integer :: n, info

    !--------------------------------------------------------------------------- 

    A_lu = A
    n = size(A,1)

    call DGETRF(n, n, A_lu, n, ipiv, info)
    if (info /=0 ) then
      LOG_ERROR("linalgebra_SolveLinEq",*) "Matrix is singular"
      call PRC_abort
    end if

    x(:) = b
    call DGETRS('N', n, 1, A_lu, n, ipiv, x, n, info)
    if (info /=0 ) then
      LOG_ERROR("linalgebra_SolveLinEq",*)  "Matrix inversion is failed"
      call PRC_abort
    end if

    return
  end subroutine linalgebra_SolveLinEq_b1D

!OCL SERIAL
  subroutine linalgebra_SolveLinEq_b2D(A, b, x)

    real(RP), intent(in) :: A(:,:)
    real(RP), intent(in) :: b(:,:)
    real(RP), intent(out) :: x(size(b,1),size(b,2))

    real(RP) :: A_lu(size(A,1),size(A,2))
    integer :: ipiv(size(A,1))
    integer :: n, info

    !--------------------------------------------------------------------------- 

    A_lu = A
    n = size(A,1)

    call DGETRF(n, n, A_lu, n, ipiv, info)
    if (info /=0 ) then
      LOG_ERROR("linalgebra_SolveLinEq",*) "Matrix is singular"
      call PRC_abort
    end if

    x(:,:) = b
    call DGETRS('N', n, size(b,2), A_lu, n, ipiv, x, n, info)
    if (info /=0 ) then
      LOG_ERROR("linalgebra_SolveLinEq",*)  "Matrix inversion is failed"
      call PRC_abort
    end if

    return
  end subroutine linalgebra_SolveLinEq_b2D

!OCL SERIAL
  subroutine linalgebra_LU(A_lu, ipiv)

    real(RP), intent(inout) :: A_lu(:,:)
    integer, intent(out) :: ipiv(size(A_lu,1))

    integer :: n
    integer :: info

    !--------------------------------------------------------------------------- 
    n = size(A_lu,1)

    call DGETRF(n, n, A_lu, n, ipiv, info)
    if (info /=0 ) then
      LOG_ERROR("linalgebra_LU",*)  "Matrix is singular"
      call PRC_abort
    end if

    return
  end subroutine linalgebra_LU

  subroutine linalgebra_SolveLinEq_GMRES(A, b, x, m, restart_num, CONV_CRIT)

    use scale_sparsemat, only: sparsemat_matmul
    implicit none
    
    type(sparsemat), intent(in) :: A
    real(RP), intent(in) :: b(:)
    real(RP), intent(inout) :: x(:)
    integer, intent(in) :: m
    integer, intent(in) ::restart_num
    real(RP), intent(in) :: CONV_CRIT

    real(RP) :: x0(size(x))
    real(RP) :: w(size(x))
    real(RP) :: v(size(x),m+1)
    real(RP) :: z(size(x),m+1)
    real(RP) :: Av(size(x))
    real(RP) :: g(size(x)+1)
    real(RP) :: r0_l2, r0_l2_new
    real(RP) :: h(m+1,m)
    real(RP) :: r(m,m)
    real(RP) :: c(m), s(m)
    real(RP) :: y(m)
    real(RP) :: tmp, tmp1, tmp2
    integer :: i, j
    integer :: restart_i

    type(sparsemat) :: M_PC

    !--------------------------------------------------------------------------- 

    call PreCondStep_ILU0_constructmat(A, M_PC)
    x0(:) = x

    do restart_i=1, restart_num

      if (restart_i == 1) then
       call SparseMat_matmul(A, x, Av)
       v(:,1) = b - Av
       r0_l2 = sqrt(sum(v(:,1)**2))
      end if
      v(:,1) = v(:,1)/r0_l2
      g(1) = r0_l2

      do j=1, m

        call PreCondStep_ILU0_solve(M_PC, v(:,j), z(:,j))

        call SparseMat_matmul(A, z(:,j), w)

        do i=1, j
          h(i,j) = sum(w*v(:,i))
          w(:) = w - h(i,j)*v(:,i)
        end do
        h(j+1,j) = sqrt(sum(w**2))
        v(:,j+1) = w(:)/h(j+1,j)

        r(1,j) = h(1,j)
        do i=1, j-1
          tmp1 =  c(i)*r(i,j) + s(i)*h(i+1,j)
          tmp2 = -s(i)*r(i,j) + c(i)*h(i+1,j)
          r(i,j) = tmp1
          r(i+1,j) = tmp2
        end do

        tmp = sqrt(r(j,j)**2 + h(j+1,j)**2)
        c(j) = r(j,j)/tmp
        s(j) = h(j+1,j)/tmp

        g(j+1) = -s(j)*g(j)
        g(j) = c(j)*g(j)
        r(j,j) = c(j)*r(j,j) + s(j)*h(j+1,j)
      end do

      y(m) = g(m)/r(m,m)
      do i=m-1,1,-1
        y(i) = (g(i) - sum(r(i,i+1:m)*y(i+1:m)))/r(i,i)
      end do

      x(:) = x + matmul(Z(:,1:m),y)

      !write(*,*) "x:", x
      !write(*,*) "x0:", x0

      call SparseMat_matmul(A, x, Av)
      v(:,1) = b - Av
      r0_l2_new = sqrt(sum(v(:,1)**2))

      !write(*,*) "r:", r0_l2, "->", r0_l2_new
      if (r0_l2_new < CONV_CRIT) exit
      r0_l2 = r0_l2_new
    end do

    return
  end subroutine linalgebra_SolveLinEq_GMRES

  subroutine PreCondStep_PtJacobi(A, b, x)
    implicit none
    type(sparsemat), intent(in) :: A
    real(RP), intent(in) :: b(:)
    real(RP), intent(out) :: x(size(b))

    integer :: n
    !--------------------------------------------------------------------------- 

    do n=1, size(x)
      x(n) = b(n) / A%GetVal(n, n)
    end do

    return
  end subroutine PreCondStep_PtJacobi

  subroutine PreCondStep_ILU0_constructmat(A, M)
    use scale_const, only: &
      EPS => CONST_EPS
    implicit none

    type(SparseMat), intent(in) :: A
    type(SparseMat), intent(inout) :: M

    integer :: i, j, k, n
    real(RP) :: Mij, M_ik, M_kk
    !--------------------------------------------------------------------------- 

    n = A%rowPtrSize-1

    M = A
    do i=2, n
      do k=1, i-1
        M_ik = M%GetVal(i,k)
        M_kk = M%GetVal(k,k)
        if ( abs(M_ik) < EPS .and. abs(M_kk) < EPS ) then
          M_ik = M_ik / M_kk
          call M%ReplaceVal( i, k,  M_ik )

          do j=k+1, n
            Mij = M%GetVal(i,j)
            if ( abs(Mij) < EPS ) then
              call M%ReplaceVal( i, j, Mij - M%GetVal(i,k) * M%GetVal(k,j) )
            end if
          end do
        end if
      end do
    end do

    return
  end subroutine PreCondStep_ILU0_constructmat

  subroutine PreCondStep_ILU0_solve(M, b, x)
    use scale_sparsemat, only: &
      SPARSEMAT_STORAGE_TYPEID_CSR
    implicit none

    type(SparseMat), intent(in) :: M
    real(RP), intent(in) :: b(:)
    real(RP), intent(out) :: x(size(b))

    integer :: i
    integer :: j
    integer :: n

    integer :: j1, j2
    !--------------------------------------------------------------------------- 

    if ( M%GetStorageFormatId() /= SPARSEMAT_STORAGE_TYPEID_CSR ) then
      LOG_ERROR("linalgebra_PreCondStep_ILU0_solve",*)  "The strorge type of specified sparse matrix is not supported. Check!"
      call PRC_abort
    end if

    n = size(x)

    x(1) = b(1)
    do i=2, n
      j1 = M%rowPtr(i)
      j2 = M%rowPtr(i+1)-1
      x(i) = b(i)
      do j=j1, j2
        if (M%colIdx(j) <= i-1) &
          x(i) = x(i) - M%val(j)*x(M%colIdx(j))
      end do
    end do

    x(n) = x(n) / M%GetVal(n,n)
    do i=n-1, 1, -1
      j1 = M%rowPtr(i)
      j2 = M%rowPtr(i+1)-1
      do j=j1, j2
        if (M%colIdx(j) >= i+1) &
          x(i) = x(i) - M%val(j)*x(M%colIdx(j))
      end do
      x(i) = x(i)/M%GetVal(i,i)
    end do

    return
  end subroutine PreCondStep_ILU0_solve

end module scale_linalgebra
