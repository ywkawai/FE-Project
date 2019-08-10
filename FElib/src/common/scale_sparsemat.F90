#include "scalelib.h"
module scale_sparsemat
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision  
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  

  type, public :: sparsemat
     real(RP), allocatable :: val(:)
     integer, allocatable :: colInd(:)
     integer, allocatable :: rowPtr(:)
     integer :: rowPtrSize
  contains
    procedure :: Init => sparsemat_Init
    procedure :: Final => sparsemat_Final
    procedure :: Print => sparsemat_Print
  end type sparsemat

  public :: sparsemat_matmul
  public :: sparsemat_GetVal
  public :: sparsemat_SetVal


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

  
contains
  subroutine sparsemat_Init(this, mat, EPS)
 
    class(SparseMat), intent(inout) :: this
    real(RP), intent(in) :: mat(:,:)
    real(RP), optional, intent(in) :: EPS

    integer :: i
    integer :: j

    integer :: val_counter
    integer :: rowptr_counter
    
    real(DP) :: tmp_val(size(mat)+1)
    integer :: tmp_colInd(size(mat)+1)
    integer :: tmp_rowptr(0:size(mat,1)+1)
    logical :: is_zero_row
    real(RP) :: EPS_ = 1.0E-15_RP
    
    !--------------------------------------------------------------------------- 
 
    val_counter    = 1
    rowptr_counter = 1
    tmp_rowptr(1)  = 1
    tmp_val(:)     = 0.0_RP

    if (present(EPS)) then
      EPS_ = EPS
    end if

    do i=1, size(mat,1)
      do j=1, size(mat,2)
        if ( abs(mat(i,j)) > EPS_ ) then
            tmp_val(val_counter) = mat(i,j)
            tmp_colInd(val_counter) = j
            val_counter = val_counter + 1
        end if
      end do

      rowptr_counter = rowptr_counter + 1
      tmp_rowptr(rowptr_counter) = val_counter
    end do
    
    allocate( this%val(val_counter) )
    allocate( this%colInd(val_counter) )
    allocate( this%rowPtr(rowptr_counter) )

    this%val(:) = tmp_val(1:val_counter)
    this%colInd(:) = tmp_colInd(1:val_counter)
    this%rowPtr(:) = tmp_rowptr(1:rowptr_counter)
    this%rowPtrSize = rowptr_counter
    
    ! write(*,*) "--- Mat ------"
    ! write(*,*) "shape:", shape(mat)
    ! write(*,*) "size:", size(mat)
    ! do j=1, size(mat,2)
    !    write(*,*) mat(:,j)
    ! end do
    ! write(*,*) "--- Compressed Mat ------"
    ! write(*,*) "shape (Mat, colInd, rowPtr):", &
    !   & shape(this%val), shape(this%colInd), shape(this%rowPtr)
    ! write(*,*) "val:", this%val(:)
    ! write(*,*) "colInd:", this%colInd(:)
    ! write(*,*) "rowPtr:", this%rowPtr(:)

  end subroutine sparsemat_Init

  subroutine sparsemat_Final(this)

    class(SparseMat), intent(inout) :: this

    !--------------------------------------------------------------------------- 

    deallocate( this%val )
    deallocate( this%colInd )
    deallocate( this%rowPtr )
    
  end subroutine sparsemat_Final

  function sparsemat_GetVal(A, i, j) result(v)

    type(sparsemat), intent(in) :: A
    integer, intent(in) :: i, j

    real(DP) ::v
    integer :: n

    !--------------------------------------------------------------------------- 

    v = 0.0_RP
    do n=A%rowPtr(i),A%rowPtr(i+1)-1
      if (A%colInd(n) == j) v = A%val(n)
    end do

  end function sparsemat_GetVal

  subroutine sparsemat_SetVal(A, i, j, v)
    
    type(sparsemat), intent(inout) :: A
    integer, intent(in) :: i, j
    real(RP), intent(in) ::v

    integer :: n

    !--------------------------------------------------------------------------- 

    do n=A%rowPtr(i),A%rowPtr(i+1)-1
      if (A%colInd(n) == j) then
        A%val(n) = v
      end if
    end do

  end subroutine sparsemat_SetVal

  subroutine sparsemat_print(A)
    class(sparsemat), intent(in) :: A

    real(RP) :: row_val(A%rowPtrSize-1)
    integer :: p
    integer :: j1, j2, j

    !--------------------------------------------------------------------------- 

    write(*,*) "-- print matrix:"
    write(*,*) "rowPtr:", A%rowPtr(:)
    write(*,*) "val:"

    j1 = A%rowPtr(1)
    do p=1, A%rowPtrSize-1
       j2 = A%rowPtr(p+1)
       row_val(:) = 0d0
       row_val(A%colInd(j1:j2-1)) = A%val(j1:j2-1)
       write(*,*) row_val(:)
       j1 = j2
    end do
  
  end subroutine sparsemat_print

  subroutine sparsemat_matmul(A, b, c)

    class(sparsemat), intent(in) :: A
    real(RP), intent(in) :: b(:)
    real(RP), intent(out) :: c(:)

    integer :: p
    integer :: j1, j2, j

    !--------------------------------------------------------------------------- 

!!$    call mkl_dcsrgemv( 'N', A%rowPtrSize-1, A%val, A%rowPtr, A%colInd, b, c)
    j1 = A%rowPtr(1)
    do p=1, A%rowPtrSize-1
       j2 = A%rowPtr(p+1) 
       c(p) = 0d0
       do j=j1, j2-1
          c(p) = c(p) + A%val(j)*b(A%colInd(j))
       end do
       j1 = j2
    end do
        
  end subroutine sparsemat_matmul

end module scale_sparsemat

