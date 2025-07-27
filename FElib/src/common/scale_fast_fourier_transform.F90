!-------------------------------------------------------------------------------
!> module common / sparsemat
!!
!! @par Description
!!          module to treat sparse matrix and the associated operations
!!
!! @author Yuta Kawai, Team SCALE
!!
!! -
!! @par Reference
!!  - Bluestein 1970:
!!    A linear filtering approach to the computation of discrete Fourier transform
!!    IEEE Transactions on Audio and Electroacoustics, 18(4), 451-455. 
!<
#include "scaleFElib.h"
module scale_fast_fourier_transform
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_const, only:    &
    UNDEF8 => CONST_UNDEF8, &
    PI => CONST_PI
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  type, public :: FastFourierTransform1D
    integer :: N
    integer :: M
    complex(RP), allocatable :: W_forward(:)
    complex(RP), allocatable :: W_backward(:)

    ! for N /= 2^n    
    logical :: use_bluestein
    complex(RP), allocatable :: W1_bluestein(:)
    complex(RP), allocatable :: W2_bluestein(:)

  contains
    procedure :: Init => FastFourierTransform1D_Init
    procedure :: Final => FastFourierTransform1D_Final
    procedure :: Forward => FastFourierTransform1D_forward
    procedure :: Backward => FastFourierTransform1D_backward
  end type FastFourierTransform1D

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
!OCL SERIAL
  subroutine FastFourierTransform1D_Init( this, N )
    implicit none
    class(FastFourierTransform1D), intent(inout) :: this
    integer, intent(in) :: N

    complex(RP) :: JJ
    integer :: k

    real(RP) :: angle0
    real(RP) :: angle

    integer :: M
    integer :: MM

    complex(RP) :: w1, w2
    !-----------------------------------

    JJ = cmplx( 0.0_RP, 1.0_RP )

    this%N = N
    if ( is_power_of_two(N) ) then
      this%use_bluestein = .false.
      MM = N / 2
      LOG_INFO("FastFourierTransform1D_Init",*) "use_bluestein=F"
      LOG_INFO("FastFourierTransform1D_Init",*) "MM=", MM
    else
      this%use_bluestein = .true.
      M = 1
      do while ( M < 2*N-1 )
        M = M * 2
      end do
      MM = M / 2
      LOG_INFO("FastFourierTransform1D_Init",*) "use_bluestein=T"
      LOG_INFO("FastFourierTransform1D_Init",*) "M=", M, "MM=", MM
    end if
    this%M = M

    allocate( this%W_forward(MM), this%W_backward(MM) )

    angle0 = 2.0_RP * PI / real(2*MM, kind=RP)
    !$omp parallel do private(angle)
    do k=0, MM-1
      angle = angle0 * real(k, kind=RP)
      this%W_forward (k+1) = exp(- JJ * angle )
      this%W_backward(k+1) = exp(  JJ * angle )
    end do

    if ( this%use_bluestein ) then
      allocate( this%W1_bluestein(N), this%W2_bluestein(N) )
      !$omp parallel do private(angle)
      do k=0, N-1
        angle = PI * real(k*k,kind=RP) / real(N, kind=RP)
        this%W1_bluestein(k+1) = exp(- JJ * angle )
        this%W2_bluestein(k+1) = exp(  JJ * angle )
      end do
    end if

    return
  end subroutine FastFourierTransform1D_Init

!OCL SERIAL
  subroutine FastFourierTransform1D_forward( this, q, s )
    implicit none
    class(FastFourierTransform1D), intent(inout) :: this
    real(RP), intent(in) :: q(this%N)
    complex(RP), intent(out) :: s(this%N)

    complex(RP) :: q_(this%N)
    integer :: i

    real(RP) :: scaling
    !-----------------------------------

    scaling = 1.0_RP / real( this%N )
    do i=1, this%N
        q_(i) = cmplx(q(i), 0.0_RP, kind=RP) * scaling
    end do
    if ( this%use_bluestein ) then
      call fft1d_bluestein( s, &
        q_, this%W_forward, this%W_backward, this%W1_bluestein, this%W2_bluestein, &
        this%N, this%M )
    else
      call fft1d_norecursive_core( s, &
        q_, this%W_forward, this%N )
      ! call fft1d_recursive_core( s, &
      !   q_, this%W_forward, this%N, this%N )

    end if
    return
  end subroutine FastFourierTransform1D_forward

!OCL SERIAL
  subroutine FastFourierTransform1D_backward( this, s, q )
    implicit none
    class(FastFourierTransform1D), intent(inout) :: this
    complex(RP), intent(in) :: s(this%N)
    real(RP), intent(out) :: q(this%N)

    complex(RP) :: q_(this%N)
    integer :: i
    !-----------------------------------

    if ( this%use_bluestein ) then
      call fft1d_bluestein( q_, &
        s, this%W_forward, this%W_backward, this%W2_bluestein, this%W1_bluestein, &
        this%N, this%M )
    else
      call fft1d_norecursive_core( q_, &
        s, this%W_backward, this%N )
    end if  

    do i=1, this%N
        q(i) = real(q_(i), kind=RP)
    end do
    return
  end subroutine FastFourierTransform1D_backward

!OCL SERIAL
  subroutine FastFourierTransform1D_Final( this )
    implicit none
    class(FastFourierTransform1D), intent(inout) :: this
    !-----------------------------------

    deallocate( this%W_forward, this%W_backward)

    if ( this%use_bluestein ) then
      deallocate( this%W1_bluestein, this%W2_bluestein )
    end if
    return
  end subroutine FastFourierTransform1D_Final

!- private -------------
!OCL SERIAL
  subroutine fft1d_bluestein( x_out, x_in, W_forward, W_backward, W1, W2, N, M )
    implicit none
    integer, intent(in) :: N, M
    complex(RP), intent(out) :: x_out(N)
    complex(RP), intent(in) :: x_in(N)
    complex(RP), intent(in) :: W_forward(M/2)  !< twiddle factor for forward transformation
    complex(RP), intent(in) :: W_backward(M/2) !< twiddle factor for backward transformation
    complex(RP), intent(in) :: W1(N)
    complex(RP), intent(in) :: W2(N)

    complex(RP) :: a(M), b(M)
    complex(RP) :: AA(M), BB(M), CC(M)
    integer :: i
    real(RP) :: scaling
    !---------------------------

    !- Chirp pre-processing
    a(:) = 0.0_RP
    b(:) = 0.0_RP
    do i=0, N-1
      a(i+1) = x_in(i+1) * W1(i+1)
      b(i+1) = W2(i+1)
    end do
    ! Fill symmetric part of b
    do i=1, N-1
      b(M-i+1) = W2(i+1)
    end do

    !- Convolution C = A * B 

    call fft1d_norecursive_core( AA, &
      a, W_forward, M )
    call fft1d_norecursive_core( BB, &
      b, W_forward, M )

    call fft1d_norecursive_core( CC, &
      AA(:) * BB(:), W_backward, M )
    
    !- Chirp post-processing
    scaling = 1.0_RP / real( M )    
    do i=1, N
      x_out(i) = scaling * CC(i) * W1(i)
    end do
    return
  end subroutine fft1d_bluestein

!OCL SERIAL
  subroutine fft1d_norecursive_core( x_out, x_in, W, N )
    implicit none
    integer, intent(in) :: N
    complex(RP), intent(out) :: x_out(N)
    complex(RP), intent(in) :: x_in(N)
    complex(RP), intent(in) :: W(N/2) !< twiddle factor

    integer :: rev
    integer :: nbits
    complex(RP) :: tmp(N)

    integer :: i, j
    integer :: k
    integer :: m
    integer :: ind

    integer :: stage
    integer :: step

    integer :: s
    complex(RP) :: t, u
    !---------------------------

    !- bit-reversal permutation
    tmp(:) = x_in(:)

    nbits = int( log(real(N,kind=RP)) / log(2.0_RP) + 0.5_RP )
    do i=0, N-1
      rev = 0
      do j=0, nbits-1
        if ( i/2**j - 2*(i/2**(j+1)) == 1 ) rev = rev + 2**(nbits-j-1)
      end do
      x_out(rev+1) = tmp(i+1)
    end do

    !- Iterative butterfly using W
    m = 1
    do stage=1, nbits
      step = 2*m
      s = N / step
      do k=0, N-1, step
        do j=0, m-1
          ind = j*s
          t = W(ind+1) * x_out(k+j+m+1)
          u = x_out(k+j+1)
          x_out(k+j+1  ) = u + t
          x_out(k+j+1+m) = u - t
        end do
      end do
      m = step
    end do
    return
  end subroutine fft1d_norecursive_core

!OCL SERIAL
  recursive subroutine fft1d_recursive_core( x_out, x_in, W, N, N0 )
    implicit none
    integer, intent(in) :: N
    complex(RP), intent(out) :: x_out(N)
    complex(RP), intent(in) :: x_in(N)
    complex(RP), intent(in) :: W(N0/2) !< twiddle factor
    integer, intent(in) :: N0

    complex(RP) :: even(N/2), odd(N/2)
    complex(RP) :: E(N/2), O(N/2)

    integer :: k
    integer :: s
    complex(RP) :: twiddle
    !---------------------------

    if ( N==1 ) then
      x_out(1) = x_in(1)
      return
    end if

    even(:) = x_in(1:N:2)
    odd (:) = x_in(2:N:2)

    call fft1d_recursive_core( E, &
      even, W, N/2, N0 )
    call fft1d_recursive_core( O, &
       odd, W, N/2, N0 )

    s = N0 / N
    do k=0, N/2-1
      x_out(k+1    ) = E(k+1) + W(s*k+1) * O(k+1)
      x_out(k+1+N/2) = E(k+1) - W(s*k+1) * O(k+1)
    end do
    return
  end subroutine fft1d_recursive_core

!OCL SERIAL
  logical function is_power_of_two(N)
    implicit none
    integer, intent(in) :: N
    !--------------
    is_power_of_two = (N > 0) .and. (iand(N, N-1) == 0)
    return
  end function is_power_of_two  
end module scale_fast_fourier_transform
