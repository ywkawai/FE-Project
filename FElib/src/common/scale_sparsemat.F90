#include "scaleFElib.h"
module scale_sparsemat
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_const, only:    &
    UNDEF8 => CONST_UNDEF8, &
    CONST_EPS
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  

  type, public :: sparsemat
    integer :: M, N
    real(RP), allocatable :: val(:)
    integer, allocatable :: colInd(:)

    ! for CSR format
    integer, allocatable :: rowPtr(:)
    integer :: rowPtrSize

    ! for ELL format
    integer :: col_size

    integer :: storage_format_id
  contains
    procedure :: Init => sparsemat_Init
    procedure :: Final => sparsemat_Final
    procedure :: Print => sparsemat_Print
    procedure :: ReplaceVal => sparsemat_ReplceVal
    procedure :: GetVal => sparsemat_GetVal
  end type sparsemat

  public :: sparsemat_matmul

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: SPARSEMAT_STORAGE_TYPEID_CSR = 1
  integer, public, parameter :: SPARSEMAT_STORAGE_TYPEID_ELL = 2

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: sparsemat_matmul_CSR
  private :: sparsemat_matmul_ELL

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  
  !-----------------------------------------------------------------------------

contains
  subroutine sparsemat_Init( this, mat, &
      EPS, storage_format )
    implicit none

    class(SparseMat), intent(inout) :: this
    real(RP), intent(in) :: mat(:,:)
    real(RP), optional, intent(in) :: EPS
    character(len=*), optional, intent(in) :: storage_format

    integer :: i
    integer :: j
    integer :: l

    integer :: val_counter
    integer :: rowptr_counter
    
    real(DP) :: tmp_val(size(mat)+1)
    integer :: tmp_colInd(size(mat)+1)
    integer :: tmp_rowptr(0:size(mat,1)+1)

    real(RP) :: EPS_
    character(len=H_MID) :: storage_format_ = 'CSR' 

    integer :: row_nonzero_counter(size(mat,1))
    integer :: col_size_l

    !--------------------------------------------------------------------------- 
 
    this%M = size(mat,1)
    this%N = size(mat,2)

    val_counter     = 1
    rowptr_counter  = 1
    tmp_rowptr(1)   = 1
    tmp_val(:)      = 0.0_RP

    if (present(EPS)) then
      EPS_ = EPS
    else
      EPS_ = CONST_EPS
    end if

    if (present(storage_format)) then
      storage_format_ = storage_format
    end if

    select case (storage_format_)
    case ('CSR')
      this%storage_format_id = SPARSEMAT_STORAGE_TYPEID_CSR

      do i=1, this%M
        do j=1, this%N
          if ( abs(mat(i,j)) > EPS_ ) then
            tmp_val(val_counter) = mat(i,j)
            tmp_colInd(val_counter) = j
            val_counter = val_counter + 1
          end if
        end do
        rowptr_counter = rowptr_counter + 1
        tmp_rowptr(rowptr_counter) = val_counter
      end do
    case ('ELL')
      this%storage_format_id = SPARSEMAT_STORAGE_TYPEID_ELL

      do i=1, this%M
        row_nonzero_counter(i) = 0
        do j=1, this%N
          if ( abs(mat(i,j)) > EPS_ ) &
            row_nonzero_counter(i) = row_nonzero_counter(i) + 1
        end do
      end do
      this%col_size = maxval(row_nonzero_counter(:))
      val_counter   = this%M * this%col_size

      tmp_val(:)    = 0.0_RP
      tmp_colInd(:) = 1
      do i=1, this%M
        col_size_l = 0
        do j=1,this%N
          if ( abs(mat(i,j)) > EPS_ ) then
              col_size_l = col_size_l + 1
              l = i+(col_size_l-1)*this%M
              tmp_val   (l) = mat(i,j)   
              tmp_colInd(l) = j
          end if
        end do
      end do
    case default
      LOG_ERROR("SparseMat_Init",*) 'Not appropriate names of storage format. Check!', trim(storage_format_)
      call PRC_abort
    end select

    allocate( this%val(val_counter) )
    allocate( this%colInd(val_counter) )
    this%val(:)    = tmp_val   (1:val_counter)
    this%colInd(:) = tmp_colInd(1:val_counter)

    select case (storage_format_)
    case ('CSR')
      allocate( this%rowPtr(rowptr_counter) )
      this%rowPtr(:) = tmp_rowptr(1:rowptr_counter)
      this%rowPtrSize = rowptr_counter
    end select

    !write(*,*) "--- Mat ------"
    !write(*,*) "shape:", shape(mat)
    !write(*,*) "size:", size(mat)
    ! do j=1, size(mat,2)
    !    write(*,*) mat(:,j)
    ! end do
    !write(*,*) "--- Compressed Mat ------"
    !write(*,*) "shape (Mat, colInd, rowPtr):", &
    !  & shape(this%val), shape(this%colInd), shape(this%rowPtr)
    ! write(*,*) "val:", this%val(:)
    ! write(*,*) "colInd:", this%colInd(:)
    ! write(*,*) "rowPtr:", this%rowPtr(:)

    return
  end subroutine sparsemat_Init

  subroutine sparsemat_Final(this)
    implicit none
    class(SparseMat), intent(inout) :: this

    !--------------------------------------------------------------------------- 

    deallocate( this%val )
    deallocate( this%colInd )

    select case( this%storage_format_id )
    case( SPARSEMAT_STORAGE_TYPEID_CSR )
      deallocate( this%rowPtr )
    end select

    return
  end subroutine sparsemat_Final

  function sparsemat_GetVal(A, i, j) result(v)
    class(sparsemat), intent(in) :: A
    integer, intent(in) :: i, j

    real(DP) ::v
    integer :: n
    integer :: jj
    !--------------------------------------------------------------------------- 

    v = UNDEF8
    select case( A%storage_format_id )
    case( SPARSEMAT_STORAGE_TYPEID_CSR )
      do n=A%rowPtr(i),A%rowPtr(i+1)-1
        if (A%colInd(n) == j) then
          v = A%val(n); exit
        end if
      end do
    case ( SPARSEMAT_STORAGE_TYPEID_ELL )
      do jj=1, A%col_size
        n = i + (jj-1)*A%M
        if (A%colInd(n) == j) then
          v = A%val(n); exit
        end if
      end do
    end select

    return
  end function sparsemat_GetVal

  subroutine sparsemat_ReplceVal(A, i, j, v)
    
    class(sparsemat), intent(inout) :: A
    integer, intent(in) :: i, j
    real(RP), intent(in) ::v

    integer :: n
    integer :: jj
    !--------------------------------------------------------------------------- 

    select case ( A%storage_format_id )
    case ( SPARSEMAT_STORAGE_TYPEID_CSR )
      do n=A%rowPtr(i),A%rowPtr(i+1)-1
        if (A%colInd(n) == j) then
          A%val(n) = v; exit
        end if
      end do
    case ( SPARSEMAT_STORAGE_TYPEID_ELL )
      do jj=1, A%col_size
        n = i + (jj-1)*A%M
        if (A%colInd(n) == j) then
          A%val(n) = v; exit
        end if
      end do
    end select

    return
  end subroutine sparsemat_ReplceVal

  subroutine sparsemat_print(A)
    implicit none
    class(sparsemat), intent(in) :: A

    real(RP) :: row_val(A%N)
    integer :: p
    integer :: j1, j2
    !--------------------------------------------------------------------------- 

    write(*,*) "-- print matrix:"
    select case ( A%storage_format_id )
    case ( SPARSEMAT_STORAGE_TYPEID_CSR )
      write(*,*) "rowPtr:", A%rowPtr(:)
    end select
    write(*,*) "colInd:", A%colInd(:)
    write(*,*) "val:"

    select case ( A%storage_format_id )
    case ( SPARSEMAT_STORAGE_TYPEID_CSR )
      j1 = A%rowPtr(1)
      do p=1, A%rowPtrSize-1
        j2 = A%rowPtr(p+1)
        row_val(:) = 0.0_RP
        row_val(A%colInd(j1:j2-1)) = A%val(j1:j2-1)
        write(*,*) row_val(:)
        j1 = j2
      end do
    case ( SPARSEMAT_STORAGE_TYPEID_ELL )
      do p=1, A%M
        row_val(:) = 0.0_RP
        do j1=1, A%col_size
          j2 = p + (j1-1)*A%M
         row_val(A%colInd(j2)) = A%val(j2)
        end do
        write(*,*) row_val(:)
      end do
    end select

    return
  end subroutine sparsemat_print

!OCL SERIAL  
  subroutine sparsemat_matmul(A, b, c)
    implicit none

    type(sparsemat), intent(in) :: A
    real(RP), intent(in ) :: b(:)
    real(RP), intent(out) :: c(size(b))

    !--------------------------------------------------------------------------- 

    select case( A%storage_format_id )
    case( SPARSEMAT_STORAGE_TYPEID_CSR )
      call sparsemat_matmul_CSR( A, b, c )
    case( SPARSEMAT_STORAGE_TYPEID_ELL )
      call sparsemat_matmul_ELL( A, b, c )
    end select

    return
  end subroutine sparsemat_matmul

!--- private ----------------------------------------------

!OCL SERIAL
  subroutine sparsemat_matmul_CSR(A, b, c)
    implicit none

    type(sparsemat), intent(in) :: A
    real(RP), intent(in ) :: b(:)
    real(RP), intent(out) :: c(size(b))

    integer :: p
    integer :: j1, j2, j

    !--------------------------------------------------------------------------- 

    !call mkl_dcsrgemv( 'N', A%rowPtrSize-1, A%val, A%rowPtr, A%colInd, b, c)
    j1 = A%rowPtr(1)
    do p=1, A%rowPtrSize-1
       j2 = A%rowPtr(p+1) 
       c(p) = 0.0_RP
       do j=j1, j2-1
          c(p) = c(p) + A%val(j) * b(A%colInd(j))
       end do
       j1 = j2
    end do

    return
  end subroutine sparsemat_matmul_CSR

!OCL SERIAL
  subroutine sparsemat_matmul_ELL(A, b, c)
    implicit none

    type(sparsemat), intent(in) :: A
    real(RP), intent(in ) :: b(:)
    real(RP), intent(out) :: c(size(b))

    integer :: k, kk, i, ii
    integer :: j_ptr
    integer :: N

    !--------------------------------------------------------------------------- 

    N = size(b)

    c(:) = 0.0_RP
    do k=1, A%col_size
      kk = N * (k-1)
      do i=1, N
        j_ptr = kk + i        
        ii = A%colInd(j_ptr)
        c(i) = c(i) + A%val(j_ptr) * b(ii)
      end do
    end do

    return
  end subroutine sparsemat_matmul_ELL

end module scale_sparsemat

