#include "scaleFElib.h"
program test_fast_fourier_transform
  use scale_precision
  use scale_prc
  use scale_io
  use scale_const, only: &
    PI => CONST_PI

  use scale_fast_fourier_transform, only: &
    FastFourierTransform1D
  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: UNITTEST_assert_less
  
  implicit none

  integer, parameter :: N1 = 256
  integer, parameter :: N2 = 240
  !-------------------------------

  call init()

  !****
  call test_FFT( N1 )
  call test_FFT( N2 )

  call final()

contains
  !---------------------------------------------------------------------------------
  subroutine test_FFT( N )
    implicit none
    integer, intent(in) :: N

    type(FastFourierTransform1D) :: fft

    real(RP) :: dx
    real(RP) :: x(N)

    real(RP) :: g_q(N)
    complex(RP) :: s_q(N)
    real(RP) :: g_q_ans(N)
    complex(RP) :: s_q_ans(N)

    integer :: i
    character(len=H_MID) :: lbl
    !-------------------------------

    dx = 1.0_RP / real(N, kind=RP)
    call fft%Init( N )

    !-
    do i=1, N
      x(i) = real(i-1, kind=RP) * dx
      g_q_ans(i) = 2.0_RP + cos(4.0_RP*PI*x(i)) + 2.0_RP * sin( 6.0_RP * PI * x(i)) 
    end do
    s_q_ans(:) = 0.0_RP
    s_q_ans(1  ) = cmplx(2.0_RP, 0.0_RP, kind=RP)
    s_q_ans(1+2) = cmplx(0.5_RP, 0.0_RP, kind=RP)
    s_q_ans(1+3) = cmplx(0.0_RP, -1.0_RP, kind=RP)
    s_q_ans(N-2) = cmplx(0.0_RP, 1.0_RP, kind=RP)
    s_q_ans(N-1) = cmplx(0.5_RP, 0.0_RP, kind=RP)

    call fft%Forward( g_q_ans, s_q )
    call fft%Backward( s_q, g_q )
    ! do i=1, N
    !   LOG_INFO('test_FFT',*) i, real(s_q(i)), aimag(s_q(i)), ":", abs(s_q(i))
    ! end do

    write(lbl,'(a,i5)') 'FFT forward N=', N
    call check_ans_cmplx( s_q, trim(lbl), s_q_ans, N )
    write(lbl,'(a,i5)') 'FFT backward N=', N
    call check_ans( g_q, trim(lbl), g_q_ans, N )

    call fft%Final()
    return
  end subroutine test_FFT
  subroutine init()
    implicit none

    !------------------------------

    call UNITTEST_init( "test_fast_fourier_transform", conf_fname="test.conf" )
    return
  end subroutine init

  subroutine final()
    implicit none
    !------------------------------

    call UNITTEST_final()
    return
  end subroutine final

  subroutine check_ans( x, lbl, x_ans_, N )
    implicit none
    integer, intent(in) :: N
    real(RP), intent(in) :: x(N)
    character(len=*), intent(in) :: lbl
    real(RP), intent(in) :: x_ans_(N)

    real(RP) :: l2error
    real(RP), parameter :: EPS = 4.0E-15_RP
    !----------------------------------------------

    l2error = sqrt( sum( (x(:) - X_ans_(:))**2 ) ) / real(N, kind=RP)
    call UNITTEST_assert_less( l2error, EPS, lbl )
    return
  end subroutine check_ans
  subroutine check_ans_cmplx( x, lbl, x_ans_, N )
    implicit none
    integer, intent(in) :: N
    complex(RP), intent(in) :: x(N)
    character(len=*), intent(in) :: lbl
    complex(RP), intent(in) :: x_ans_(N)

    real(RP) :: l2error
    real(RP), parameter :: EPS = 4.0E-15_RP
    !----------------------------------------------

    l2error = sqrt( sum( abs(x(:) - X_ans_(:))**2 ) ) / real(N, kind=RP)
    call UNITTEST_assert_less( l2error, EPS, lbl )
    return
  end subroutine check_ans_cmplx
end program test_fast_fourier_transform
