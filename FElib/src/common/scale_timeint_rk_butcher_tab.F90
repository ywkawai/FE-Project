!-------------------------------------------------------------------------------
!> module common / Runge-Kutta scheme 
!!
!! @par Description
!!      Butcher tableau
!!
!! @author Team SCALE
!<
#include "scaleFElib.h"
module scale_timeint_rk_butcher_tab
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision  
  use scale_io
  use scale_prc
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  public :: timeint_rk_butcher_tab_get_info
  public :: timeint_rk_butcher_tab_get

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
  subroutine timeint_rk_butcher_tab_get_info( rk_scheme_name, &
    nstage, tend_buf_size, low_storage_flag, imex_flag        )
    implicit none
    character(len=*), intent(in) :: rk_scheme_name
    integer, intent(out) :: nstage
    integer, intent(out) :: tend_buf_size
    logical, intent(out) :: low_storage_flag
    logical, intent(out) :: imex_flag
    !-----------------------------------------------------------

    select case(rk_scheme_name)
    case('Euler')
      nstage = 1
      tend_buf_size = 1
      low_storage_flag = .true.
      imex_flag = .false.    
    case('RK_4')
      nstage = 4
      tend_buf_size = 1
      low_storage_flag = .false.
      imex_flag = .false.
    case('RK_TVD_2')
      nstage = 2
      tend_buf_size = 1
      low_storage_flag = .true.
      imex_flag = .false.
    case('RK_TVD_3')
      nstage = 3
      tend_buf_size = 1
      low_storage_flag = .true.
      imex_flag = .false.
    case('ARK_232')
      nstage = 3
      tend_buf_size = 3
      low_storage_flag = .false.
      imex_flag = .true.      
    case('ARK_324')
      nstage = 4
      tend_buf_size = 4
      low_storage_flag = .false.
      imex_flag = .true.        
    case('HEVI_debug')
      nstage = 1
      tend_buf_size = 1
      low_storage_flag = .false.
      imex_flag = .true.         
    case default
      LOG_ERROR("timeint_rk_butcher_tab_get_info",*) trim(rk_scheme_name)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine timeint_rk_butcher_tab_get_info

  subroutine timeint_rk_butcher_tab_get( &
    rk_scheme_name, nstage, imex_flag,         & ! (in)
    coef_a_ex, coef_b_ex, coef_c_ex,           & ! (in)
    coef_a_im, coef_b_im, coef_c_im,           & ! (in)
    tend_buf_indmap                            ) ! (out)

    implicit none
    character(len=*), intent(in) :: rk_scheme_name
    integer, intent(in) :: nstage
    logical, intent(in) :: imex_flag
    real(RP), intent(out) :: coef_a_ex(nstage,nstage)
    real(RP), intent(out) :: coef_b_ex(nstage)
    real(RP), intent(out) :: coef_c_ex(nstage)
    real(RP), intent(out) :: coef_a_im(nstage,nstage)
    real(RP), intent(out) :: coef_b_im(nstage)
    real(RP), intent(out) :: coef_c_im(nstage)
    integer, intent(out) :: tend_buf_indmap(nstage)

    real(RP) :: alp, gam, del
    !-----------------------------------------------------------

    coef_a_ex(:,:) = 0.0_RP
    coef_b_ex(:)   = 0.0_RP
    coef_c_ex(:)   = 0.0_RP

    if (imex_flag) then    
      coef_a_im(:,:) = 0.0_RP
      coef_b_im(:) = 0.0_RP
      coef_c_im(:) = 0.0_RP
    end if

    select case(rk_scheme_name)
    case('Euler')
      coef_a_ex(1,1) = 1.0_RP
      coef_b_ex(:) = (/ 1.0_RP /)
      tend_buf_indmap(:) = 1      
    case('RK_4')      
      coef_a_ex(2,1) = 0.5_RP
      coef_a_ex(3,2) = 0.5_RP
      coef_a_ex(4,3) = 1.0_RP
      coef_b_ex(:) = (/ 1.0_RP, 2.0_RP, 2.0_RP, 1.0_RP /)/6.0_RP
      tend_buf_indmap(:) = 1
    case('RK_TVD_2')
      coef_a_ex(1,1) = 1.0_RP
      coef_a_ex(2,2) = 0.5_RP
      tend_buf_indmap(:) = 1
    case('RK_TVD_3')
      coef_a_ex(1,1) = 1.0_RP
      coef_a_ex(2,2) = 1.0_RP/4.0_RP
      coef_a_ex(3,3) = 2.0_RP/3.0_RP
      tend_buf_indmap(:) = 1
    case('ARK_232')
      alp = (3.0_RP + 2.0_RP*sqrt(2.0_RP))/6.0_RP
      gam = 1.0_RP - 1.0_RP/sqrt(2.0_RP)
      del = 1.0_RP/(2.0_RP*sqrt(2.0_RP))
      
      coef_a_ex(2,1) = 2.0_RP*gam
      coef_a_ex(3,:) = (/ 1.0_RP - alp, alp, 0.0_RP /)
      coef_b_ex(:) = (/ del, del, gam /)

      coef_a_im(2,:) = (/ gam, gam, 0.0_RP /)
      coef_a_im(3,:) = (/ del, del, gam /)
      coef_b_im(:) = coef_b_ex(:)
      
      tend_buf_indmap(:) = (/ 1, 2, 3 /)   
    case('ARK_324')
      
      coef_a_ex(2,1) = 1767732205903.0_RP/2027836641118.0_RP
      coef_a_ex(3,1) = 5535828885825.0_RP/10492691773637.0_RP
      coef_a_ex(3,2) = 788022342437.0_RP/10882634858940.0_RP
      coef_a_ex(4,1) = 6485989280629.0_RP/16251701735622.0_RP
      coef_a_ex(4,2) = -4246266847089.0_RP/9704473918619.0_RP
      coef_a_ex(4,3) = 10755448449292.0_RP/10357097424841.0_RP
      coef_b_ex(1) = 1471266399579.0_RP/7840856788654.0_RP
      coef_b_ex(2) = -4482444167858.0_RP/7529755066697.0_RP
      coef_b_ex(3) = 11266239266428.0_RP/11593286722821.0_RP
      coef_b_ex(4) = 1767732205903.0_RP/4055673282236.0_RP

      coef_a_im(2,1:2) = 1767732205903.0_RP/4055673282236.0_RP
      coef_a_im(3,1) = 2746238789719.0_RP/10658868560708.0_RP
      coef_a_im(3,2) = -640167445237.0_RP/6845629431997.0_RP
      coef_a_im(3,3) = 1767732205903.0_RP/4055673282236.0_RP
      coef_a_im(4,:) = coef_b_ex(:)
      coef_b_im(:)   = coef_b_ex(:)
      
      tend_buf_indmap(:) = (/ 1, 2, 3, 4 /)   

    case ('HEVI_debug')
      coef_a_ex(1,:) = (/ 1.0_RP /)
      coef_a_im(1,:) = (/ 0.0_RP /)
      coef_b_ex(:) = (/ 1.0_RP /)
      coef_b_im(:) = (/ 1.0_RP /)

      tend_buf_indmap(:) = 1.0_RP
    end select

    return
  end subroutine timeint_rk_butcher_tab_get

end module scale_timeint_rk_butcher_tab
