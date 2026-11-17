#include "scalelib.h"
program test_timeint_rk
  use scale_precision
  use scale_prc
  use scale_io
  use scale_const, only: &
    PI => CONST_PI

  use scale_timeint_rk, only: &
    timeint_rk
  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: UNITTEST_assert_true
  
  implicit none

  real(RP), parameter :: end_time = 2.0_RP
  real(RP), parameter :: dt1 = 0.02_RP
  integer, parameter :: save_error_dstep1 = 1
  real(RP), parameter :: dt2 = 0.0025_RP
  integer, parameter :: save_error_dstep2 = 8

  integer, parameter :: error_array_size = 100
  real(RP), save :: answer(error_array_size)
  real(RP), save :: error1_mem(error_array_size)
  real(RP), save :: error2_mem(error_array_size)

  !-------------------------------

  call init()

  !****

  ! explicit schemes
  call check_tscheme_ex( 'ERK_1s1o', 0.98_RP )
  call check_tscheme_ex( 'ERK_4s4o', 3.98_RP )
  call check_tscheme_ex( 'ERK_SSP_2s2o', 1.98_RP )
  call check_tscheme_ex( 'ERK_SSP_3s3o', 2.98_RP )
  call check_tscheme_ex( 'ERK_SSP_4s3o', 2.98_RP )
  call check_tscheme_ex( 'ERK_SSP_5s3o_2N2*', 2.98_RP )
  call check_tscheme_ex( 'ERK_SSP_10s4o_2N', 3.98_RP )

  ! IMEX schemes
  call check_tscheme_imex( 'IMEX_ARK232', 1.98_RP )
  call check_tscheme_imex( 'IMEX_ARK324', 2.98_RP )

  ! Check conversion betweem Shu-Osher and Butcher forms
  call check_ShuOsher2Butcher()

  call final()

contains
  subroutine check_tscheme_ex( temporal_scheme_name, check_ord )
    character(len=*), intent(in) :: temporal_scheme_name
    real(RP), intent(in) :: check_ord
    !------------------------------------------------

    call calc_sol_ex( temporal_scheme_name, dt1, answer, error1_mem, save_error_dstep1 )
    call calc_sol_ex( temporal_scheme_name, dt2, answer, error2_mem, save_error_dstep2 )
    call check_terror( temporal_scheme_name, error1_mem, error2_mem, dt1 / dt2, check_ord  )
  
    return
  end subroutine check_tscheme_ex

  subroutine check_tscheme_imex( temporal_scheme_name, check_ord )
    character(len=*), intent(in) :: temporal_scheme_name
    real(RP), intent(in) :: check_ord
    !------------------------------------------------

    call calc_sol_imex( temporal_scheme_name, dt1, answer, error1_mem, save_error_dstep1 )
    call calc_sol_imex( temporal_scheme_name, dt2, answer, error2_mem, save_error_dstep2 )
    call check_terror( temporal_scheme_name, error1_mem, error2_mem, dt1 / dt2, check_ord  )
  
    return
  end subroutine check_tscheme_imex

  subroutine check_terror( tname, error_dt1, error_dt2, dt_ratio, check_ord )
    implicit none
    character(*), intent(in) :: tname
    real(RP), intent(in) :: error_dt1(error_array_size)
    real(RP), intent(in) :: error_dt2(error_array_size)
    real(RP), intent(in) :: dt_ratio
    real(RP), intent(in) :: check_ord

    real(RP) :: error1, error2
    real(RP) :: conv_rate
    !------------------------------------------------

    error1 = sqrt(sum(error_dt1(:)**2))
    error2 = sqrt(sum(error_dt2(:)**2))

    conv_rate = log(error1 / error2) / log(dt_ratio)
    write(*,'(a,a,f5.2,a,e12.5,a,e12.5)') trim(tname), ": conv rate=", conv_rate, &
      ', error1=', error1, ", error2=", error2

    call UNITTEST_assert_true( conv_rate >= check_ord, trim(tname)//': convergence rate', &
      'Convergence rate is too low.' )

    return
  end subroutine check_terror

  subroutine calc_sol_ex( tint_type_name, dt, &
      ans_mem, error_mem, save_error_dstep )
    implicit none

    character(*), intent(in) :: tint_type_name
    real(RP), intent(in) :: dt
    real(RP), intent(out) :: ans_mem(error_array_size)
    real(RP), intent(out) :: error_mem(error_array_size)
    integer, intent(in) :: save_error_dstep

    type(timeint_rk) :: tint

    integer :: n
    integer :: nstep
    integer :: rkstage
    integer :: tintbuf_ind

    real(RP) :: u(1,1)
    integer, parameter :: ID_U = 1
    real(RP) :: v(1,1)
    integer, parameter :: ID_V = 2
    real(RP) :: omg

    integer :: err_count
    !-----------------------------------

    !-
    omg = 2.0_RP * PI
    nstep = end_time / dt
    u = 1.0_RP
    v = 0.0_RP
    err_count = 0

    call tint%Init( tint_type_name, &
      dt, 2, 2, (/ 1, 1 /) )
    
    !$acc data copyin( u, v )
    do n = 1, nstep
      do rkstage = 1, tint%nstage
        tintbuf_ind = tint%tend_buf_indmap(rkstage)

        !$acc serial present( tint%tend_buf2D_ex, u, v )
        tint%tend_buf2D_ex(1,1,ID_U,tintbuf_ind) &
          = + omg * v(1,1)
        tint%tend_buf2D_ex(1,1,ID_V,tintbuf_ind) &
          = - omg * u(1,1)
        !$acc end serial

        call tint%Advance( rkstage, u, ID_U, &
          1, 1, 1, 1 )
        call tint%Advance( rkstage, v, ID_V, &
          1, 1, 1, 1 )
      end do
      
      if ( mod(n,save_error_dstep) == 0 ) then
        !$acc update host( u )
        err_count = err_count + 1
        ans_mem(err_count) = cos(omg * dble(n)*dt)
        error_mem(err_count) = u(1,1) - ans_mem(err_count)
        ! write(*,'(4f12.5)') n*dt, u, ans_mem(err_count), error_mem(err_count) 
      end if
    end do
    !$acc end data

    !-
    call tint%Final()

    return
  end subroutine calc_sol_ex

  subroutine calc_sol_imex( tint_type_name, dt, &
    ans_mem, error_mem, save_error_dstep )

    implicit none

    character(*), intent(in) :: tint_type_name
    real(RP), intent(in) :: dt
    real(RP), intent(out) :: ans_mem(error_array_size)
    real(RP), intent(out) :: error_mem(error_array_size)
    integer, intent(in) :: save_error_dstep

    type(timeint_rk) :: tint

    integer :: n
    integer :: nstep
    integer :: rkstage
    integer :: tintbuf_ind
    real(RP) :: implicit_fac

    real(RP) :: u(1,1)
    integer, parameter :: ID_U = 1
    real(RP) :: v(1,1)
    integer, parameter :: ID_V = 2
    real(RP) :: omg
    real(RP) :: omgg
    real(RP) :: r

    real(RP) :: ui
    real(RP) :: vi
    real(RP) :: coef
    
    integer :: err_count
    !-----------------------------------

    !-
    omg  = 2.0_RP * PI
    r    = 0.1_RP
    omgg = sqrt( omg**2 - r**2 / 4.0_RP )

    nstep = end_time / dt
    u = 1.0_RP
    v = 0.0_RP
    err_count = 0

    call tint%Init( tint_type_name, &
      dt, 2, 2, (/ 1, 1 /) )
    
    !$acc data copyin( u, v )
    do n = 1, nstep
      do rkstage = 1, tint%nstage
        tintbuf_ind = tint%tend_buf_indmap(rkstage)
        implicit_fac = tint%Get_implicit_diagfac(rkstage)

        !--
        coef = implicit_fac * omg

        if ( abs(implicit_fac) > 0.0_RP ) then
          !$acc serial present( tint%tend_buf2D_im, u, v )
          ui = ( u(1,1) + coef * v(1,1) ) / ( 1.0_RP + coef**2 )
          vi = ( v(1,1) - coef * u(1,1) ) / ( 1.0_RP + coef**2 ) 
          tint%tend_buf2D_im(1,1,ID_U,tintbuf_ind) = ( ui - u(1,1) ) / implicit_fac
          tint%tend_buf2D_im(1,1,ID_V,tintbuf_ind) = ( vi - v(1,1) ) / implicit_fac
          !$acc end serial
        else
          !$acc serial present( tint%tend_buf2D_im, u, v )
          tint%tend_buf2D_im(1,1,ID_U,tintbuf_ind) &
            = + omg * v(1,1)
          tint%tend_buf2D_im(1,1,ID_V,tintbuf_ind) &
            = - omg * u(1,1)
          !$acc end serial
        end if

        call tint%StoreImplicit( rkstage, u, ID_U, &
          1, 1, 1, 1 )
        call tint%StoreImplicit( rkstage, v, ID_V, &
          1, 1, 1, 1 )

        !--
        
        !$acc serial present( tint%tend_buf2D_ex, u, v )
        tint%tend_buf2D_ex(1,1,ID_U,tintbuf_ind) &
          = - r * u(1,1)
        tint%tend_buf2D_ex(1,1,ID_V,tintbuf_ind) &
          = 0.0_RP
        !$acc end serial

        call tint%Advance( rkstage, u, ID_U, &
          1, 1, 1, 1 )
        call tint%Advance( rkstage, v, ID_V, &
          1, 1, 1, 1 )
      end do

      if ( mod(n,save_error_dstep) == 0 ) then
        !$acc update host( u )

        err_count = err_count + 1
        ans_mem(err_count) = &
            ( cos(omgg * dble(n)*dt) - r/(2.0_RP * omgg) * sin(omgg * dble(n)*dt) ) &
          * exp(- 0.5_RP * r * dble(n)*dt)
        error_mem(err_count) = u(1,1) - ans_mem(err_count)
!        write(*,'(4f12.5)') n*dt, u, ans_mem(err_count), error_mem(err_count) 
      end if
    end do
    !$acc end data
    
    !-
    call tint%Final()

    return
  end subroutine calc_sol_imex

  subroutine check_ShuOsher2Butcher()
    implicit none

    type(timeint_rk) :: tint
    real(RP) :: A_3s3o(3,3), b_3s3o(3)
    real(RP) :: A_4s3o(4,4), b_4s3o(4)
    !------------------------------

    !* ERK_SSP_3s3o

    call tint%Init( 'ERK_SSP_3s3o', &
      1.0_RP, 1, 1, (/ 1, 1 /) )

    A_3s3o(:,:) = 0.0_RP
    A_3s3o(2,1) = 1.0_RP
    A_3s3o(3,1:2) = 1.0_RP / 4.0_RP
    b_3s3o(:) = (/ 1.0_RP, 1.0_RP, 4.0_RP /) / 6.0_RP
    call check_ButcherMat( 'ERK_SSP_3s3o', tint, A_3s3o, b_3s3o )

    call tint%Final()

    !* ERK_SSP_4s3o

    call tint%Init( 'ERK_SSP_4s3o', &
      1.0_RP, 1, 1, (/ 1, 1 /) )

    A_4s3o(:,:) = 0.0_RP
    A_4s3o(2,1) = 1.0_RP / 2.0_RP
    A_4s3o(3,1:2) = 1.0_RP / 2.0_RP
    A_4s3o(4,1:3) = 1.0_RP / 6.0_RP
    b_4s3o(:) = (/ 1.0_RP, 1.0_RP, 1.0_RP, 3.0_RP /) / 6.0_RP
    call check_ButcherMat( 'ERK_SSP_4s3o', tint, A_4s3o, b_4s3o )

    call tint%Final()

    return
  end subroutine check_ShuOsher2Butcher

  subroutine check_ButcherMat( tint_label, &
    tint, a_ans, b_ans )
    implicit none

    character(len=*), intent(in) :: tint_label
    type(timeint_rk), intent(in) :: tint
    real(RP), intent(in) :: a_ans(tint%nstage,tint%nstage)
    real(RP), intent(in) :: b_ans(tint%nstage)
    !---------------------------------

    ! Conversion from Shu-Osher form into Butcher form
    call UNITTEST_assert_true( &
      sqrt( sum( ( a_ans(:,:) - tint%coef_a_ex(:,:) )**2 ) ) <= 1.0E-14_RP * dble(tint%nstage**2), &
      'check_ShuOsher2Butcher: '//trim(tint_label)//': a', 'Fail to convert Shu-Osher form into Butcher form.' )
    call UNITTEST_assert_true( &
      sqrt( sum( ( b_ans(:) - tint%coef_b_ex(:) )**2 ) ) <= 1.0E-14_RP * dble(tint%nstage), &
      'check_ShuOsher2Butcher: '//trim(tint_label)//': b', 'Fail to convert Shu-Osher form into Butcher form.' )

    return
  end subroutine check_ButcherMat
  

  !---------------------------------------------------------------------------------

  subroutine init()
    implicit none

    !------------------------------

    call UNITTEST_init( "test_timeint_rk", conf_fname="test.conf" )

    return
  end subroutine init

  subroutine final()
    implicit none
    !------------------------------

    call UNITTEST_final()

    return
  end subroutine final
  
end program test_timeint_rk
