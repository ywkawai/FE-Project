!-------------------------------------------------------------------------------
!> module FElib / element/ SIAC filter
!!
!! @par Description
!!      A module to provide a Smoothness-Increasing Accuracy-Increasing (SIAC) filter
!!
!! @par Reference
!!  - Cockburn et al. 2003:
!!    Enhanced accuracy by post-processing for finite element methods for hyperbolic equations.
!!    Mathematics of Computation, 72(242), 577–506. 
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_element_SIACfilter
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_element_base, only: &
    ElementBase1D
  use scale_element_line, only: &
    LineElement
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !

  !-----------------------------------------------------------------------------
  !
  !++ Public type
  !
  type, public :: SIAC_filter
    integer :: spline_ord
    integer :: spline_num
    integer :: spline_r

    integer :: Npts_per_elem !< Number of sampling points per element
    
    integer :: KernelHalfW   !< Half width of kernel function 

    integer :: NintGLPt
    real(RP), allocatable :: IntrpMat(:,:,:,:)
    real(RP), allocatable :: int_x(:)
    real(RP), allocatable :: int_w(:)

    real(RP), allocatable :: kernel_func_coef(:)
    real(RP), allocatable :: kernel_func(:,:,:,:)

    real(RP), allocatable :: x_pts_per_elem(:)
  contains
    procedure :: Init => SIAC_filter_Init
    procedure :: Final => SIAC_filter_Final
    procedure :: Apply1D => SIAC_filter_apply1D
    procedure :: Get_kernel_func => SIAC_filter_get_kernel_func
  end type SIAC_filter

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  private :: central_Bspline

contains
!> Initialize a object for SIAC filter
!!
!! @param r the number of B-splines - 1
!! @param l the order of B-spline
!OCL SERIAL
  subroutine SIAC_filter_Init( this, &
    r, l, &
    x_pts_per_elem, elem1D )
    use scale_polynominal, only: &
      Polynominal_GenLagrangePoly, &
      Polynominal_GenLegendrePoly
    implicit none
    class(SIAC_filter), intent(inout) :: this
    integer, intent(in) :: r
    integer, intent(in) :: l
    real(RP), intent(in) :: x_pts_per_elem(:)
    class(LineElement), intent(in) :: elem1D

    integer :: i
    integer :: p
    integer :: m

    real(RP), allocatable :: IntrpMat_dummy(:,:)
    real(RP) :: x0, x1, x2
    real(RP), allocatable :: int_x_tmp(:)

    real(RP), allocatable :: P1D_ori(:,:)
    !--------------------------------------

    this%spline_r = r
    this%spline_num = r + 1
    this%spline_ord = l
    
    this%KernelHalfW = ceiling(0.5_RP * real(r+l,kind=RP))

    this%Npts_per_elem = size(x_pts_per_elem)

    !---
    allocate( this%x_pts_per_elem(this%Npts_per_elem) )
    this%x_pts_per_elem(:) = x_pts_per_elem(:)

    this%NintGLPt = ceiling( ( r + l ) / 2.0_RP )
    allocate( IntrpMat_dummy(this%NintGLPt,elem1D%Np) )
    allocate( this%int_x(this%NintGLPt), this%int_w(this%NintGLPt) )

    IntrpMat_dummy(:,:) = elem1D%GenIntGaussLegendreIntrpMat( this%NintGLPt, this%int_w, this%int_x )

    allocate( this%IntrpMat(this%NintGLPt,elem1D%Np,2,this%Npts_per_elem) )
    allocate( int_x_tmp(this%NintGLPt) )
    allocate( P1D_ori(this%NintGLPt,elem1D%Np) )
    do i=1, this%Npts_per_elem
      x0 = - 1.0_RP; x2 = 1.0_RP
      x1 = x_pts_per_elem(i)

      int_x_tmp(:) = x0 + 0.5_RP * (x1 - x0) * ( 1.0_RP + this%int_x(:) )
      P1D_ori(:,:) = Polynominal_GenLegendrePoly( elem1D%PolyOrder, int_x_tmp(:) )
      do p=1, elem1D%Np
        P1D_ori(:,p) = P1D_ori(:,p) * sqrt(real(p-1,kind=RP) + 0.5_RP)
      end do
      this%IntrpMat(:,:,1,i) = matmul( P1D_ori, elem1D%invV )
      ! this%IntrpMat(:,:,1,i) = Polynominal_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, int_x_tmp(:) )
      
      int_x_tmp(:) = x1 + 0.5_RP * (x2 - x1) * ( 1.0_RP + this%int_x(:) )
      P1D_ori(:,:) = Polynominal_GenLegendrePoly( elem1D%PolyOrder, int_x_tmp(:) )
      do p=1, elem1D%Np
        P1D_ori(:,p) = P1D_ori(:,p) * sqrt(real(p-1,kind=RP) + 0.5_RP)
      end do
      this%IntrpMat(:,:,2,i) = matmul( P1D_ori, elem1D%invV )
      ! this%IntrpMat(:,:,2,i) = Polynominal_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, int_x_tmp(:) )
    end do

    !---
    allocate( this%kernel_func_coef(0:r) )
    call calculate_kernel_func_coef( this%kernel_func_coef, &
        r, l, this%int_x, this%int_w, this%NintGLPt )
    
    allocate( this%kernel_func(this%NintGLPt,2,-this%KernelHalfW:this%KernelHalfW,this%Npts_per_elem) )

    do i=1, this%Npts_per_elem
      call construct_kernel_func( this%kernel_func(:,:,:,i), &
        r, l, this%kernel_func_coef, this%KernelHalfW, &
        this%int_x, this%NintGLPt, x_pts_per_elem(i) )
    end do

    LOG_INFO("SIAC_filter_Init",*) "r, l=", r, l        
    LOG_INFO("SIAC_filter_Init",*) "KernelFunc coef:", this%kernel_func_coef
    do i=1, this%Npts_per_elem
      LOG_INFO("SIAC_filter_Init",*) "--- KernelFunc xi=", x_pts_per_elem(i)
      do m=-this%KernelHalfW,this%KernelHalfW
        LOG_INFO("SIAC_filter_Init",*) this%kernel_func(:,1,m,i), ":", this%kernel_func(:,2,m,i)
      end do
    end do

    do i=1, this%Npts_per_elem
      LOG_INFO("SIAC_filter_Init",*) "--- InterpMat xi=", x_pts_per_elem(i)
      LOG_INFO("SIAC_filter_Init",*) "L", this%IntrpMat(1,:,1,i)
      LOG_INFO("SIAC_filter_Init",*) "L", this%IntrpMat(2,:,1,i)
      LOG_INFO("SIAC_filter_Init",*) "R", this%IntrpMat(1,:,2,i)
      LOG_INFO("SIAC_filter_Init",*) "R", this%IntrpMat(2,:,2,i)
    end do

    return
  end subroutine SIAC_filter_Init

!OCL SERIAL
  subroutine SIAC_filter_Final( this )
    implicit none
    class(SIAC_filter), intent(inout) :: this
    !--------------------------------------

    deallocate( this%IntrpMat, this%int_x, this%int_w )
    deallocate( this%kernel_func )

    deallocate( this%x_pts_per_elem )

    return
  end subroutine SIAC_filter_Final

!OCL SERIAL
  subroutine SIAC_filter_apply1D( this, filtered_q, &
    q, Np1D, Ne, Nmesh, NmeshHalo )
    implicit none
    class(SIAC_filter), intent(in) :: this
    integer, intent(in) :: Np1D
    integer, intent(in) :: Ne
    integer, intent(in) :: Nmesh
    integer, intent(in) :: NmeshHalo
    real(RP), intent(out) :: filtered_q(this%Npts_per_elem,Ne)
    real(RP), intent(in) :: q(Np1D,Ne*Nmesh)
    !--------------------------------------------

    call SIAC_filter_apply_core( filtered_q, &
      q, Np1D, Ne, Nmesh, NmeshHalo, this%Npts_per_elem, &
      this%kernel_func, this%x_pts_per_elem, this%IntrpMat, this%int_w,       &
      this%KernelHalfW, this%NintGLPt )

    return
  end subroutine SIAC_filter_apply1D

!OCL SERIAL
  subroutine SIAC_filter_apply_core( filtered_q,      &
    q, Np1D, Ne, Nmesh, NmeshHalo, Npts_per_elem,     &
    kernel_func, xi, IntrpMat, int_w, HalfW, NintGLPt )
    implicit none
    integer, intent(in) :: Np1D
    integer, intent(in) :: Ne
    integer, intent(in) :: Nmesh
    integer, intent(in) :: NmeshHalo
    integer, intent(in) :: Npts_per_elem 
    real(RP), intent(out) :: filtered_q(Npts_per_elem,Ne)
    real(RP), intent(in) :: q(Np1D,Ne*Nmesh)
    integer, intent(in) :: NintGLpt
    real(RP), intent(in) :: kernel_func(NintGLpt,2,-HalfW:HalfW,Npts_per_elem)
    real(RP), intent(in) :: xi(Npts_per_elem)
    real(RP), intent(in) :: IntrpMat(NintGLpt,Np1D,2,Npts_per_elem)
    real(RP), intent(in) :: int_w(NintGLpt)
    integer, intent(in) :: HalfW

    integer :: p
    integer :: m
    
    integer :: ke_os
    integer :: ke, kee

    real(RP) :: q_intrp(NintGLpt,2)
    real(RP) :: int_w2(NintGLpt,2)
    real(RP) :: tmp
    !--------------------------------------------

    ke_os = Ne * NmeshHalo
    !$omp parallel do collapse(2) private(kee, tmp, m, q_intrp, int_w2)
    do ke=1, Ne
    do p=1, Npts_per_elem
      kee = ke + ke_os
      tmp = 0.0_RP

      int_w2(:,1) = ( xi(p) + 1.0_RP ) * int_w(:)
      int_w2(:,2) = ( 1.0_RP - xi(p) ) * int_w(:)
      do m=-HalfW, HalfW
        q_intrp(:,1) = matmul( IntrpMat(:,:,1,p), q(:,kee+m) )
        q_intrp(:,2) = matmul( IntrpMat(:,:,2,p), q(:,kee+m) )
        tmp = tmp + sum( int_w2(:,:) * q_intrp(:,:) * kernel_func(:,:,m,p) )
      end do
      filtered_q(p,ke) = tmp * 0.25_RP
    end do
    end do

    return
  end subroutine SIAC_filter_apply_core  

!OCL SERIAL
  subroutine SIAC_filter_get_kernel_func( this, kernel_func, Np, &
    x )
    implicit none
    class(SIAC_filter), intent(in) :: this
    integer, intent(in) :: Np
    real(RP), intent(in) :: x(Np)
    real(RP), intent(out) :: kernel_func(Np)

    integer :: i
    integer :: gam
    real(RP) :: kernel_func_tmp
    real(RP) :: x_gam
    real(RP) :: psi
    !--------------------------------------

    do i=1, Np
      kernel_func_tmp = 0.0_RP
      do gam=0, this%spline_r
        x_gam = - 0.5_RP * real(this%spline_r, kind=RP) + gam

        psi = central_Bspline(x(i) - x_gam, this%spline_ord)
        kernel_func_tmp = kernel_func_tmp + psi * this%kernel_func_coef(gam)
      end do
      kernel_func(i) = kernel_func_tmp
    end do
    return
  end subroutine SIAC_filter_get_kernel_func

!-- private --
!OCL SERAIL
  subroutine construct_kernel_func( kernel_func, &
    r, l, kernel_coef, halfW, int_xi, NintPts, xi )
    use scale_linalgebra, only: linalgebra_SolveLinEq
    implicit none
    integer, intent(in) :: halfW
    integer, intent(in) :: NintPts
    real(RP), intent(out) :: kernel_func(NintPts,2,-halfW:halfW)
    integer, intent(in) :: r
    integer, intent(in) :: l !< l=k+1
    real(RP), intent(in) :: kernel_coef(0:r)
    real(RP), intent(in) :: int_xi(NintPts)
    real(RP), intent(in ):: xi

    real(RP) :: x_gam
    integer :: gam

    real(RP) :: psi(NintPts,2)
    real(RP) :: y1, y2

    integer :: i
    integer :: m

    real(RP) :: kernel_func_tmp(NintPts,2)
    real(RP) :: x0, x1, x2
    !--------------------------------------------------------

    do i=-halfW, halfW
      kernel_func_tmp(:,:) = 0.0_RP
      do gam=0, r
        x_gam = - 0.5_RP * real(r, kind=RP) + gam
        do m=1, NintPts
          x0 = i - 0.5_RP; x2 = i + 0.5_RP
          x1 = x0 + 0.5_RP * ( xi + 1.0_RP )

          y1 = x0 + 0.5_RP * ( ( x1 - x0 ) * ( int_xi(m) + 1.0_RP ) - xi ) &
               - x_gam
          y2 = x1 + 0.5_RP * ( ( x2 - x1 ) * ( int_xi(m) + 1.0_RP ) - xi ) &
               - x_gam
          psi(m,1) = central_Bspline(y1,l)
          psi(m,2) = central_Bspline(y2,l)
        end do
        kernel_func_tmp(:,:) = kernel_func_tmp(:,:) + kernel_coef(gam) * psi(:,:)
      end do
      kernel_func(:,:,i) = kernel_func_tmp(:,:)
    end do
    return
  end subroutine construct_kernel_func

!OCL SERAIL
  subroutine calculate_kernel_func_coef( coef, &
    r, l, int_xi, int_w, NintPts )
    use scale_linalgebra, only: linalgebra_SolveLinEq
    use scale_polynominal, only: Polynominal_GenLegendrePoly
    implicit none
    integer, intent(in) :: r
    integer, intent(in) :: l !< l=k+1
    integer, intent(in) :: NintPts
    real(RP), intent(in) :: int_xi(NintPts)
    real(RP), intent(in) :: int_w(NintPts)
    real(RP), intent(out) :: coef(0:r)
    real(RP) :: gam
    real(RP) :: x_gam

    integer :: i, j
    integer :: ii, jj
    integer :: m, mm

    real(RP) :: LinMat(0:r,0:r)
    real(RP) :: b(0:r)
    real(RP) :: coef_(0:r)

    real(RP) :: x_knots(0:l)
    real(RP) :: x0, x1

    real(RP) :: x_intrp(NintPts)
    real(RP) :: psi(NintPts)

    real(RP) :: int_tmp
    real(RP) :: int_w2(NintPts)
    real(RP) :: int_coef

    real(RP) :: P(NintPts,0:r)
    real(RP) :: zero(1)
    real(RP) :: P_x0(1,0:r)

    real(RP) :: scale_s(0:r)
    real(RP) :: scale_r(0:r)
    !---------------------------------------------

    do i=0, l
      x_knots(i) = - 0.5_RP * real(l, kind=RP) + real(i, kind=RP)
    end do
    !$omp parallel do collapse(2) private(i, ii, j, jj, &
    !$omp gam, x_gam, m, int_tmp, x0, x1, int_coef, x_intrp, int_w2, mm, psi, P)
    do jj=0, r
      do ii=0, r
        i = ii; j = jj
        gam = j
        x_gam = - 0.5_RP * real(r, kind=RP) + gam

        int_tmp = 0.0_RP
        do m=0, l-1
          x0 = x_knots(m); x1 = x_knots(m+1)
          int_coef = 0.5_RP * ( x1 - x0 )

          x_intrp(:) = x0 + int_coef * ( 1.0_RP + int_xi(:) )
          int_w2(:) = int_coef * int_w(:)

          do mm=1, NintPts
            psi(mm) = central_Bspline( x_intrp(mm), l )
          end do

          ! int_tmp = int_tmp &
          !   + sum( int_w2(:) * psi(:) * ( x_intrp(:) - x_gam )**i )          
          P(:,0:i) = Polynominal_GenLegendrePoly( i, x_intrp(:) - x_gam )
          int_tmp = int_tmp &
            + sum( int_w2(:) * psi(:) * P(:,i) )

        !   if (i==1) then
        !       LOG_INFO("SIAC_filter_calc_coef",*) m, x0, x1, "Psi: ", psi(:)
        !       LOG_INFO("SIAC_filter_calc_coef",*) m, x0, x1, "x+x_gam: ", x_intrp(:) + x_gam
        !   end if
        end do
        LinMat(ii,jj) = int_tmp
      end do
    end do

    ! do i=0, r
    !   LOG_INFO("SIAC_filter_calc_coef",*) "LinMat: ",  LinMat(i,:)
    ! end do

    ! b(:) = 0.0_RP
    ! b(0) = 1.0_RP
    zero(:) = 0.0_RP
    P_x0(:,:) = Polynominal_GenLegendrePoly( r, zero(:) )
    b(:) = P_x0(1,:)

    do j=0, r
      scale_s(j) = sqrt(sum(LinMat(:,j)**2))
      if ( scale_s(j) /= 0.0_RP ) then
        LinMat(:,j) = LinMat(:,j) / scale_s(j)
      end if
    end do
    do i=0, r
      scale_r(i) = sqrt(sum(LinMat(i,:)**2))
      if ( scale_r(i) /= 0.0_RP ) then
        LinMat(i,:) = LinMat(i,:) / scale_r(i)
        b(i) = b(i) / scale_r(i)
      end if
    end do
    call linalgebra_SolveLinEq( LinMat, b, coef_ )
    do j=0, r
      coef(j) = coef_(j) / scale_s(j)
    end do

    return
  end subroutine calculate_kernel_func_coef

!OCL SERIAL
  recursive function central_Bspline( x, k ) result(b)
    implicit none
    real(RP), intent(in) :: x
    integer, intent(in) :: k
    real(RP) :: b

    real(RP) :: coef1, coef2
    !--------------------------------------

    if ( k==1 ) then
      if ( - 0.5_RP < x .and. x <= 0.5_RP ) then
        b = 1.0_RP
      else
        b = 0.0_RP
      end if
    else
      coef1 = central_Bspline( x + 0.5_RP, k-1 )
      coef2 = central_Bspline( x - 0.5_RP, k-1 )
      b = ( (  0.5_RP * k + x ) * coef1 &
          + (  0.5_RP * k - x ) * coef2 ) / real(k-1, kind=RP)
    end if
    return
  end function central_Bspline
end module scale_element_SIACfilter