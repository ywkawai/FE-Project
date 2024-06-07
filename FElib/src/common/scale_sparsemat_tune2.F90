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
    integer :: M                      !< Number of row of original matrix
    integer :: N                      !< Number of column of original matrix
    integer :: nnz                    !< Number of nonzero of the matrix
    real(RP), allocatable :: val(:)   !< Array storing the values of the nonzeros
    integer, allocatable :: colIdx(:) !< Array saving the column indices of the nonzeros

    ! for CSR format
    integer, allocatable :: rowPtr(:) !< Array saving the start and end pointers of the nonzeros of the rows.
    integer :: rowPtrSize             !< Size of rowPtr

    ! for ELL format
    integer :: col_size

    integer, private :: storage_format_id      !< Number of row of original matrix
  contains
    procedure, public :: Init => sparsemat_Init
    procedure, public :: Final => sparsemat_Final
    procedure, public :: Print => sparsemat_Print
    procedure, public :: ReplaceVal => sparsemat_ReplceVal
    procedure, public :: GetVal => sparsemat_GetVal
    procedure, public :: GetStorageFormatId => sparsemat_get_storage_format_id
  end type sparsemat

  interface sparsemat_matmul
    module procedure sparsemat_matmul1
    module procedure sparsemat_matmul2
    module procedure sprasemat_matmul3
    module procedure sprasemat_matmul4
  end interface
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
  private :: sparsemat_matmul_CSR_1
  private :: sparsemat_matmul_CSR_2
  private :: sparsemat_matmul_ELL_1
  private :: sparsemat_matmul_ELL_2

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
      EPS_ = CONST_EPS * 500.0_RP ! ~ 1x10^-13
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
      tmp_colInd(:) = -1
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

    this%nnz = val_counter
    allocate( this%val(val_counter) )
    allocate( this%colIdx(val_counter) )
    this%val(:)    = tmp_val   (1:val_counter)
    this%colIdx(:) = tmp_colInd(1:val_counter)

    select case (storage_format_)
    case ('CSR')
      allocate( this%rowPtr(rowptr_counter) )
      this%rowPtr(:) = tmp_rowptr(1:rowptr_counter)
      this%rowPtrSize = rowptr_counter
    case('ELL')
      do l=2, val_counter
        if (this%colIdx(l) == -1) then
          this%colIdx(l) = this%colIdx(l-1)
        end if
      end do
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
    deallocate( this%colIdx )

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
        if (A%colIdx(n) == j) then
          v = A%val(n); exit
        end if
      end do
    case ( SPARSEMAT_STORAGE_TYPEID_ELL )
      do jj=1, A%col_size
        n = i + (jj-1)*A%M
        if (A%colIdx(n) == j) then
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
        if (A%colIdx(n) == j) then
          A%val(n) = v; exit
        end if
      end do
    case ( SPARSEMAT_STORAGE_TYPEID_ELL )
      do jj=1, A%col_size
        n = i + (jj-1)*A%M
        if (A%colIdx(n) == j) then
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
    write(*,'(a,i5,a,i5)') "orginal matrix shape:", A%M, 'x', A%N
    write(*,'(a,i5)') "size of compressed matrix:", A%nnz
    select case ( A%storage_format_id )
    case ( SPARSEMAT_STORAGE_TYPEID_CSR )
      write(*,*) "rowPtr:", A%rowPtr(:)
    end select
    write(*,*) "colInd:", A%colIdx(:)
    write(*,*) "val:"

    select case ( A%storage_format_id )
    case ( SPARSEMAT_STORAGE_TYPEID_CSR )
      j1 = A%rowPtr(1)
      do p=1, A%rowPtrSize-1
        j2 = A%rowPtr(p+1)
        row_val(:) = 0.0_RP
        row_val(A%colIdx(j1:j2-1)) = A%val(j1:j2-1)
        write(*,*) row_val(:)
        j1 = j2
      end do
    case ( SPARSEMAT_STORAGE_TYPEID_ELL )
      do p=1, A%M
        row_val(:) = 0.0_RP
        do j1=1, A%col_size
          j2 = p + (j1-1)*A%M
         row_val(A%colIdx(j2)) = A%val(j2)
        end do
        write(*,*) row_val(:)
      end do
    end select

    return
  end subroutine sparsemat_print

  function sparsemat_get_storage_format_id( A ) result(id)
    implicit none
    class(sparsemat), intent(in) :: A
    real(RP) :: id
    !------------------------------------------------------

    id = A%storage_format_id
    return 
  end function sparsemat_get_storage_format_id

!OCL SERIAL  
  subroutine sparsemat_matmul1(A, b, c)
    implicit none

    type(sparsemat), intent(in) :: A
    real(RP), intent(in ) :: b(:)
    real(RP), intent(out) :: c(A%M)

    !--------------------------------------------------------------------------- 

    select case( A%storage_format_id )
    case( SPARSEMAT_STORAGE_TYPEID_CSR )
      call sparsemat_matmul_CSR_1( A%val, A%colIdx, A%rowPtr, b, c, &
        A%M, A%N, A%nnz, A%rowPtrSize                          )
    case( SPARSEMAT_STORAGE_TYPEID_ELL )
      call sparsemat_matmul_ELL_1( A%val, A%colIdx, b, c,           &
        A%M, A%N, A%nnz, A%col_size                            )
    end select

    return
  end subroutine sparsemat_matmul1

!OCL SERIAL  
  subroutine sparsemat_matmul2(A, b, c)
    implicit none

    type(sparsemat), intent(in) :: A
    real(RP), intent(in ) :: b(:,:)
    real(RP), intent(out) :: c(size(b, 1),A%M)

    !--------------------------------------------------------------------------- 
    select case( A%storage_format_id )
    case( SPARSEMAT_STORAGE_TYPEID_CSR )
      call sparsemat_matmul_CSR_2( A%val, A%colIdx, A%rowPtr, b, c, &
        A%M, A%N, A%nnz, A%rowPtrSize, size(b,1)               )
    case( SPARSEMAT_STORAGE_TYPEID_ELL )
      call sparsemat_matmul_ELL_2( A%val, A%colIdx, b, c,           &
        A%M, A%N, A%nnz, A%col_size, size(b,1)                 )
    end select

    return
  end subroutine sparsemat_matmul2

!OCL SERIAL
  subroutine sprasemat_matmul3(A, b, c, Nvec)
    implicit none

    type(sparsemat), intent(in) :: A
    integer, intent(in) :: Nvec
    real(RP), intent(in) :: b(:,:)
    real(RP), intent(out) :: c(size(b, 2), Nvec)

    integer :: kvec
    integer :: M
    integer :: col_size
    ! integer :: col_Ind(A%nnz)
    ! real(RP) :: A_(A%nnz)

    integer :: k, kk, i
    integer :: j_ptr
    
    !---Note: Only for ELL format
    M = A%M
    col_size = A%col_size
    
    c(:,:) = 0.0_RP

    do k=1, col_size
      kk = M * (k-1)
      do i=1, M
        j_ptr = kk + i        
        c(i, 1) = c(i, 1) + A%val(j_ptr) * b(1, A%colIdx(j_ptr))
        c(i, 2) = c(i, 2) + A%val(j_ptr) * b(2, A%colIdx(j_ptr))
        c(i, 3) = c(i, 3) + A%val(j_ptr) * b(3, A%colIdx(j_ptr))
        c(i, 4) = c(i, 4) + A%val(j_ptr) * b(4, A%colIdx(j_ptr))
        c(i, 5) = c(i, 5) + A%val(j_ptr) * b(5, A%colIdx(j_ptr))
        c(i, 6) = c(i, 6) + A%val(j_ptr) * b(6, A%colIdx(j_ptr))
      end do
    end do

    return
  end subroutine

!OCL SERIAL
  subroutine sprasemat_matmul4(A, b, c, Nvec, DUMMY)
    implicit none

    type(sparsemat), intent(in) :: A
    integer, intent(in) :: Nvec, DUMMY
    real(RP), intent(in) :: b(:,:)
    real(RP), intent(out) :: c(size(b, 2), Nvec)

    integer :: kvec
    integer :: M
    integer :: col_size
    ! integer :: col_Ind(A%nnz)
    ! real(RP) :: A_(A%nnz)

    integer :: k, kk, i
    integer :: j_ptr
    
    !---Note: Only for ELL format
    M = A%M
    col_size = A%col_size
    
    c(:,:) = 0.0_RP

    do k=1, col_size
      kk = M * (k-1)
      do i=1, M
        j_ptr = kk + i        
        c(i, 1) = c(i, 1) + A%val(j_ptr) * b(1, A%colIdx(j_ptr))
        c(i, 2) = c(i, 2) + A%val(j_ptr) * b(2, A%colIdx(j_ptr))
        c(i, 3) = c(i, 3) + A%val(j_ptr) * b(3, A%colIdx(j_ptr))
        c(i, 4) = c(i, 4) + A%val(j_ptr) * b(4, A%colIdx(j_ptr))
        c(i, 5) = c(i, 5) + A%val(j_ptr) * b(5, A%colIdx(j_ptr))
        c(i, 6) = c(i, 6) + A%val(j_ptr) * b(6, A%colIdx(j_ptr))
        c(i, 7) = c(i, 7) + A%val(j_ptr) * b(7, A%colIdx(j_ptr))
      end do
    end do

    return
  end subroutine
!--- private ----------------------------------------------

!OCL SERIAL
  subroutine sparsemat_matmul_CSR_1(A, col_Ind, rowPtr, b, c, M, N, buf_size, rowPtr_size)
    implicit none

    integer, intent(in) :: M
    integer, intent(in) :: N
    integer, intent(in) :: buf_size
    integer, intent(in) :: rowPtr_size
    real(RP), intent(in) :: A(buf_size)
    integer, intent(in) :: col_Ind(buf_size)
    integer, intent(in) :: rowPtr(rowPtr_size)
    real(RP), intent(in ) :: b(N)
    real(RP), intent(out) :: c(M)

    integer :: p, j
    integer :: j1, j2

    !--------------------------------------------------------------------------- 

    !call mkl_dcsrgemv( 'N', rowPtr_size-1, A, rowPtr, col_Ind, b, c)
    j1 = rowPtr(1)
    do p=1, rowPtr_size-1
       j2 = rowPtr(p+1) 
       c(p) = 0.0_RP
       do j=j1, j2-1
          c(p) = c(p) + A(j) * b(col_Ind(j))
       end do
       j1 = j2
    end do

    return
  end subroutine sparsemat_matmul_CSR_1

!OCL SERIAL
  subroutine sparsemat_matmul_CSR_2(A, col_Ind, rowPtr, b, c, M, N, buf_size, rowPtr_size, NQ)
    implicit none

    integer, intent(in) :: M
    integer, intent(in) :: N
    integer, intent(in) :: NQ
    integer, intent(in) :: buf_size
    integer, intent(in) :: rowPtr_size
    real(RP), intent(in) :: A(buf_size)
    integer, intent(in) :: col_Ind(buf_size)
    integer, intent(in) :: rowPtr(rowPtr_size)
    real(RP), intent(in ) :: b(NQ,N)
    real(RP), intent(out) :: c(NQ,M)

    integer :: p, j
    integer :: j1, j2

    !--------------------------------------------------------------------------- 

    !call mkl_dcsrgemv( 'N', rowPtr_size-1, A, rowPtr, col_Ind, b, c)
    j1 = rowPtr(1)
    c(:,:) = 0.0_RP
    do p=1, rowPtr_size-1
      j2 = rowPtr(p+1) 
      do j=j1, j2-1
        c(:,p) = c(:,p) + A(j) * b(:,col_Ind(j))
      end do
      j1 = j2
    end do

    return
  end subroutine sparsemat_matmul_CSR_2

!OCL SERIAL
  subroutine sparsemat_matmul_ELL_1(A, col_Ind, b, c, M, N, buf_size, col_size)
    implicit none

    integer, intent(in) :: M
    integer, intent(in) :: N
    integer, intent(in) :: buf_size
    integer, intent(in) :: col_size
    real(RP), intent(in) :: A(buf_size)
    integer, intent(in) :: col_Ind(buf_size)
    real(RP), intent(in ) :: b(N)
    real(RP), intent(out) :: c(M)

    integer :: k, kk, i
    integer :: j_ptr
    !--------------------------------------------------------------------------- 

    c(:) = 0.0_RP
    do k=1, col_size
      kk = M * (k-1)
      do i=1, M
        j_ptr = kk + i        
        c(i) = c(i) + A(j_ptr) * b(col_Ind(j_ptr))
      end do
    end do

    return
  end subroutine sparsemat_matmul_ELL_1

!OCL SERIAL
  subroutine sparsemat_matmul_ELL_2(A, col_Ind, b, c, M, N, buf_size, col_size, NQ)
    implicit none

    integer, intent(in) :: M
    integer, intent(in) :: N
    integer, intent(in) :: buf_size
    integer, intent(in) :: col_size
    integer, intent(in) :: NQ
    real(RP), intent(in) :: A(buf_size)
    integer, intent(in) :: col_Ind(buf_size)
    real(RP), intent(in ) :: b(NQ,N)
    real(RP), intent(out) :: c(NQ,M)

    integer :: k, kk, i
    integer :: j_ptr
    !--------------------------------------------------------------------------- 

    c(:,:) = 0.0_RP
    do k=1, col_size
      kk = M * (k-1)
      do i=1, M
        j_ptr = kk + i        
        c(:,i) = c(:,i) + A(j_ptr) * b(:,col_Ind(j_ptr))
      end do
    end do

    return
  end subroutine sparsemat_matmul_ELL_2

end module scale_sparsemat

