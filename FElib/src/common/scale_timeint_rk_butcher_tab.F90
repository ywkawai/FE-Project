!-------------------------------------------------------------------------------
!> module common / Runge-Kutta scheme 
!!
!! @par Description
!!      Butcher tableau
!!
!! Notation rules:
!!   * ERK schemes
!!    ERK(_abbrev1_)MsPo(_abbrev2) : explicit RK schme with M-stage and P-th order.
!!      ERK_1s1o (, Euler ) : Euler scheme
!!      ERK_4s4o (, RK4   ) : classical 4 stage and 4th-order  RK scheme
!!      ERK_SSP_MsPo : strong stability preserving (SSP) (M,P) RK scheme
!!
!!   * IMEX shcemes (we assume that ERK for non-stiff part and  DIRK for stiff part)
!!     IMEX_abbrev : abbreviation of scheme name follows the table in Vogl et al. (2019). 
!!
!! @author Team SCALE
!!
!! @par Reference
!!  - Vogl et al. 2019:
!!    Evaluation of implicit-explicit additive Runge-Kutta integrators for the HOMME-NH dynamical core. 
!!    Journal of Advances in Modeling Earth Systems, 11, 4228–4244. 
!!  - Higueras and Roldan, 2019:
!!    New third order low-storage SSP explicit Runge–Kutta methods. 
!!    Journal of Scientific Computing, 79(3), 1882-1906.
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
!OCL SERIAL
  subroutine timeint_rk_butcher_tab_get_info( rk_scheme_name, &
    nstage, tend_buf_size, low_storage_flag, imex_flag        )
    implicit none
    character(len=*), intent(in) :: rk_scheme_name
    integer, intent(out) :: nstage
    integer, intent(out) :: tend_buf_size
    logical, intent(out) :: low_storage_flag
    logical, intent(out) :: imex_flag
    !-----------------------------------------------------------

    low_storage_flag = .false.

    select case(rk_scheme_name)
    case( 'ERK_1s1o', 'ERK_Euler' )
      nstage = 1
      tend_buf_size = 1
      imex_flag = .false.    
    case( 'ERK_4s4o', 'ERK_RK4' )
      nstage = 4
      tend_buf_size = 1
      imex_flag = .false.
    case( 'ERK_SSP_2s2o' )
      nstage = 2
      tend_buf_size = 1
      low_storage_flag = .true.
      imex_flag = .false.
    case( 'ERK_SSP_3s3o' ) ! SSP coefficient: 1, effective SSP coefficient Ceff: 1/3 
      nstage = 3
      tend_buf_size = 1
      low_storage_flag = .true.
      imex_flag = .false.
    case( 'ERK_SSP_4s3o' )
      nstage = 4
      tend_buf_size = 1
      low_storage_flag = .true.
      imex_flag = .false.
    case( 'ERK_SSP_5s3o_2N2*' ) ! Higueras and Roldan, 2018: New third order low-storage SSP explicit Runge–Kutta methods 
      nstage = 5
      tend_buf_size = 1
      low_storage_flag = .true.
      imex_flag = .false.
    case ( 'ERK_SSP_10s4o_2N' ) ! Ketcheson, 2008: HIGHLY EFFICIENT STRONG STABILITY-PRESERVING RUNGE–KUTTA METHODS WITH LOW-STORAGE IMPLEMENTATIONS, SIAM J. SCI. COMPUT. 
                                ! SSP coefficient: 6, effective SSP coefficient Ceff: 0.6
      nstage = 10
      tend_buf_size = 1
      low_storage_flag = .true.
      imex_flag = .false.
    case( 'IMEX_ARK232' ) ! Giraldo et al. (2013): Implicit-Explicit Formulations of a Three-Dimensional Nonhydrostatic Unified Model
                          ! of the Atmosphere (NUMA), SIAM J. Sci. Comp.
      nstage = 3
      tend_buf_size = 3
      imex_flag = .true.      
    case( 'IMEX_ARK324' ) ! Kennedy and Carpenter, 2003: Additive Runge-Kutta schemes for convection-diffusion-reaction equations, Appl. Numer. Math.
      nstage = 4
      tend_buf_size = 4
      imex_flag = .true.
    case default
      LOG_ERROR("timeint_rk_butcher_tab_get_info",*) trim(rk_scheme_name)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine timeint_rk_butcher_tab_get_info

!OCL SERIAL
  subroutine timeint_rk_butcher_tab_get( &
    rk_scheme_name, nstage, imex_flag,         & ! (in)
    coef_a_ex, coef_b_ex, coef_c_ex,           & ! (in)
    coef_sig_ex, coef_gam_ex,                  & ! (in)
    coef_a_im, coef_b_im, coef_c_im,           & ! (in)
    tend_buf_indmap                            ) ! (out)

    implicit none
    character(len=*), intent(in) :: rk_scheme_name
    integer, intent(in) :: nstage
    logical, intent(in) :: imex_flag
    real(RP), intent(out) :: coef_a_ex(nstage,nstage)
    real(RP), intent(out) :: coef_b_ex(nstage)
    real(RP), intent(out) :: coef_c_ex(nstage)
    real(RP), intent(out) :: coef_sig_ex(nstage+1,nstage)
    real(RP), intent(out) :: coef_gam_ex(nstage+1,nstage)
    real(RP), intent(out) :: coef_a_im(nstage,nstage)
    real(RP), intent(out) :: coef_b_im(nstage)
    real(RP), intent(out) :: coef_c_im(nstage)
    integer, intent(out) :: tend_buf_indmap(nstage)

    integer :: n
    real(RP) :: alp, gam, del
    logical :: call_ShuOsher2Butcher
    
    !-----------------------------------------------------------

    coef_a_ex(:,:) = 0.0_RP
    coef_b_ex(:)   = 0.0_RP
    coef_c_ex(:)   = 0.0_RP

    coef_sig_ex(:,:) = 0.0_RP
    coef_gam_ex(:,:) = 0.0_RP

    if (imex_flag) then    
      coef_a_im(:,:) = 0.0_RP
      coef_b_im(:) = 0.0_RP
      coef_c_im(:) = 0.0_RP
    end if

    call_ShuOsher2Butcher = .false.

    !---

    select case(rk_scheme_name)
    case( 'ERK_1s1o', 'ERK_Euler' )

      coef_a_ex(1,1) = 1.0_RP
      coef_b_ex(:) = (/ 1.0_RP /)

      tend_buf_indmap(:) = 1
      
    case( 'ERK_4s4o', 'ERK_RK4' )

      coef_a_ex(2,1) = 0.5_RP
      coef_a_ex(3,2) = 0.5_RP
      coef_a_ex(4,3) = 1.0_RP
      coef_b_ex(:) = (/ 1.0_RP, 2.0_RP, 2.0_RP, 1.0_RP /)/6.0_RP

      tend_buf_indmap(:) = 1
    
    case( 'ERK_SSP_2s2o' )

      coef_sig_ex(2,1  ) = 1.0_RP
      coef_sig_ex(3,1:2) = (/ 1.0_RP, 1.0_RP /) / 2.0_RP
      coef_gam_ex(2,1  ) = 1.0_RP
      coef_gam_ex(3,2  ) = 1.0_RP / 2.0_RP

      tend_buf_indmap(:) = 1
      call_ShuOsher2Butcher = .true.
    
    case( 'ERK_SSP_3s3o' ) 

      coef_sig_ex(2,1  ) = 1.0_RP
      coef_sig_ex(3,1:2) = (/ 3.0_RP, 1.0_RP /) / 4.0_RP
      coef_sig_ex(4,1:3) = (/ 1.0_RP, 0.0_RP, 2.0_RP /) / 3.0_RP
      coef_gam_ex(2,1  ) = 1.0_RP
      coef_gam_ex(3,2  ) = 1.0_RP / 4.0_RP
      coef_gam_ex(4,3  ) = 2.0_RP / 3.0_RP

      tend_buf_indmap(:) = 1
      call_ShuOsher2Butcher = .true.

    case( 'ERK_SSP_4s3o' )

      coef_sig_ex(2,1  ) = 1.0_RP
      coef_sig_ex(3,1:2) = (/ 0.0_RP, 1.0_RP /)
      coef_sig_ex(4,1:3) = (/ 2.0_RP, 0.0_RP, 1.0_RP /) / 3.0_RP
      coef_sig_ex(5,1:4) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 1.0_RP /)
      coef_gam_ex(2,1  ) = 0.5_RP
      coef_gam_ex(3,2  ) = 0.5_RP
      coef_gam_ex(4,3  ) = 1.0_RP / 6.0_RP
      coef_gam_ex(5,4  ) = 0.5_RP

      tend_buf_indmap(:) = 1
      call_ShuOsher2Butcher = .true.

    case( 'ERK_SSP_5s3o_2N2*' )

      coef_sig_ex(2,1  ) = 1.0_RP
      coef_sig_ex(3,1:2) = (/ 0.0_RP, 1.0_RP /)
      coef_sig_ex(4,1:3) = (/ 0.682342861037239_RP, 0.0_RP, 0.317657138962761_RP /)
      coef_sig_ex(5,1:4) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 1.0_RP /)
      coef_sig_ex(6,1:5) = (/ 0.045230974482400_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.954769025517600_RP  /)
      coef_gam_ex(2,1  ) = 0.465388589249323_RP
      coef_gam_ex(3,2  ) = 0.465388589249323_RP
      coef_gam_ex(4,3  ) = 0.124745797313998_RP
      coef_gam_ex(5,4  ) = 0.465388589249323_RP
      coef_gam_ex(6,5  ) = 0.154263303748666_RP

      tend_buf_indmap(:) = 1
      call_ShuOsher2Butcher = .true.

    case( 'ERK_SSP_10s4o_2N' )

      do n=1, 4
        coef_sig_ex(n+1,n) = 1.0_RP
        coef_gam_ex(n+1,n) = 1.0_RP / 6.0_RP
      end do
      coef_sig_ex(6,1:5) = (/ 3.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 2.0_RP /) / 5.0_RP
      coef_gam_ex(6,5  ) = 1.0_RP / 15.0_RP
      do n=6, 9
        coef_sig_ex(n+1,n) = 1.0_RP
        coef_gam_ex(n+1,n) = 1.0_RP / 6.0_RP
      end do
      coef_sig_ex(11,1:10) = (/ 0.2_RP, 0.0_RP, 0.0_RP, 0.0_RP, 1.8_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 3.0_RP /) * 0.2_RP
      coef_gam_ex(11,1:10) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 1.8_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 3.0_RP /) / 30.0_RP

      tend_buf_indmap(:) = 1
      call_ShuOsher2Butcher = .true.

    case('IMEX_ARK232')

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
    
    case('IMEX_ARK324')

      coef_a_ex(2,1) = 1767732205903.0_RP/2027836641118.0_RP
      coef_a_ex(3,1) = 5535828885825.0_RP/10492691773637.0_RP
      coef_a_ex(3,2) =  788022342437.0_RP/10882634858940.0_RP
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

    case default
      LOG_ERROR("timeint_rk_butcher_tab_get",*) trim(rk_scheme_name)//' is not supported. Check!'
      call PRC_abort

    end select
    
    !--
    if ( call_ShuOsher2Butcher ) then  
      call ShuOsher2Butcher( &
        coef_sig_ex, coef_gam_ex, nstage, & ! (in)
        coef_a_ex, coef_b_ex              ) ! (out)
    end if

    !--
    coef_c_ex(:) = 0.0_RP
    do n=1, nstage
      coef_c_ex(n) = sum( coef_a_ex(n,1:n) )
    end do
    if ( imex_flag ) then
      coef_c_ex(:) = 0.0_RP
      do n=1, nstage
        coef_c_ex(n) = sum( coef_a_ex(n,1:n) )
      end do
    end if

    return
  end subroutine timeint_rk_butcher_tab_get

  !- private --------------------------------------------

  !> Convert Shu-Osher matrix to Butcher matrix
  !! For detail of the strategy, see Higueras and Roldan (2019).
!OCL SERIAL
  subroutine ShuOsher2Butcher( &
    coef_sig, coef_gam, nstage, &
    coef_a, coef_b )

    use scale_linalgebra, only: &
      linalgebra_SolveLinEq
    implicit none

    integer, intent(in) :: nstage
    real(RP), intent(in) :: coef_sig(nstage+1,nstage)
    real(RP), intent(in) :: coef_gam(nstage+1,nstage)
    real(RP), intent(out) :: coef_a(nstage,nstage)
    real(RP), intent(out) :: coef_b(nstage)

    integer :: n
    real(RP) :: L(nstage,nstage)
    !-----------------------------------------------------------

    ! L=I-SIG
    L(:,:) = - coef_sig(1:nstage,:)
    do n=1, nstage
     L(n,n) = 1.0_RP + L(n,n)
    end do
    ! A=L^-1 GAM 
    call linalgebra_SolveLinEq( L(:,:), coef_gam(1:nstage,:), coef_a(:,:) )

    ! b = GAM_N,: + SIG_N,: * A
    do n=1, nstage
      coef_b(n) = coef_gam(nstage,n) + sum( coef_sig(nstage,:) * coef_a(:,n) )
    end do 
  
    return
  end subroutine ShuOsher2Butcher

end module scale_timeint_rk_butcher_tab
