#include "scaleFElib.h"
module scale_polynominal
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_const, only: &
    PI => CONST_PI
  use scale_io
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !  

  public :: Polynominal_GenLegendrePoly
  public :: Polynominal_GenDLegendrePoly

  public :: Polynominal_GenGaussLobattoPt
  public :: Polynominal_GenGaussLobattoPtIntWeight

  public :: Polynominal_GenGaussLegendrePt
  public :: Polynominal_GenGaussLegendrePtIntWeight

  public :: Polynominal_GenLagrangePoly
  public :: Polynominal_GenDLagrangePoly_lglpt
  
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
  !> A function to obtain the values of Lagrange basis functions which are evaluated over aribitary points
  !!
  function Polynominal_GenLagrangePoly(Nord, x_lgl, x) result(l)
    implicit none
    
    integer, intent(in) :: Nord
    real(RP), intent(in) :: x_lgl(Nord+1)
    real(RP), intent(in) :: x(:)
    real(RP) :: l(size(x), Nord+1)

    integer :: n, i
    real(RP) :: P_lgl(Nord+1,Nord+1)    
    real(RP) :: P(size(x),Nord+1)
    real(RP) :: Pr(size(x),Nord+1)

    !---------------------------------------------------------------------------

    P_lgl(:,:)  = Polynominal_GenLegendrePoly(Nord, x_lgl)
    P (:,:)  = Polynominal_GenLegendrePoly(Nord, x)
    Pr(:,:) = Polynominal_GenDLegendrePoly(Nord, x, P)

    do n=1, Nord+1
      do i=1, size(x)
        if ( abs(x(i)-x_lgl(n)) < 1E-15_RP ) then
          l(i,n) = 1.0_RP
        else
          l(i,n) = &
            (x(i) - 1.0_RP)*(x(i) + 1.0_RP)*Pr(i,Nord+1)              &
            / (dble(Nord*(Nord+1))*P_lgl(n,Nord+1)*(x(i) - x_lgl(n)))
        end if
      end do
    end do

    return
  end function Polynominal_GenLagrangePoly

  !> A function to obtain the differential values of Lagrange basis functions which are evaluated over aribitary points
  !!
  function Polynominal_GenDLagrangePoly_lglpt(Nord, x_lgl) result(lr)
    implicit none

    integer, intent(in) :: Nord
    real(RP), intent(in) :: x_lgl(Nord+1)
    real(RP) :: lr(Nord+1, Nord+1)

    integer :: n, k
    real(RP) :: P(Nord+1,Nord+1)

    !---------------------------------------------------------------------------

    P(:,:)  = Polynominal_GenLegendrePoly(Nord, x_lgl)

    do n=1, Nord+1
      do k=1, Nord+1
        if (k==1 .and. n==1) then
          lr(k,n) = - 0.25_RP*dble(Nord*(Nord+1))
        else if (k==Nord+1 .and. n==Nord+1) then
          lr(k,n) = + 0.25_RP*dble(Nord*(Nord+1))
        else if (k==n) then
          lr(k,n) = 0.0_RP
        else
          lr(k,n) = P(n,Nord+1)/(P(k,Nord+1)*(x_lgl(n) - x_lgl(k)))
        end if
      end do
    end do
    
    return
  end function Polynominal_GenDLagrangePoly_lglpt

  !> A function to obtain the values of Legendre polynominals which are evaluated at aribitary points. 
  !!
  function Polynominal_GenLegendrePoly(Nord, x) result(P)
    implicit none

    integer, intent(in) :: Nord
    real(RP), intent(in) :: x(:)
    real(RP) :: P(size(x), Nord+1)

    integer :: n

    !---------------------------------------------------------------------------

    if (Nord < 0) then
       LOG_ERROR("Polynominal_GenLegendrePoly",*)   "Nord must be larger than 0."
    end if

    P(:,1) = 1.0_RP
    if (Nord==0) return

    P(:,2) = x(:)
    do n=2, Nord
      P(:,n+1) = ( dble(2*n-1)*x(:)*P(:,n) - dble(n-1)*P(:,n-1))/dble(n)
    end do
    
    return    
  end function Polynominal_GenLegendrePoly

  !> A function to obtain differential values of Legendre polynominals which are evaluated at aribitary points. 
  !! 
  function Polynominal_GenDLegendrePoly(Nord, x, P) result(GradP)
    implicit none
 
    integer, intent(in) :: Nord
    real(DP), intent(in) :: x(:)
    real(DP), intent(in) :: P(:,:)
    real(DP) :: GradP(size(x), Nord+1)

    integer :: n

    !---------------------------------------------------------------------------

    if (Nord < 0) then
      LOG_ERROR("Polynominal_GenDLegendrePoly",*)   "Nord must be larger than 0."
    end if

    GradP(:,1) = 0.0_RP
    if (Nord == 0) return

    GradP(:,2) = 1.0_RP
    do n=2, Nord
      GradP(:,n+1) = 2.0_RP*x(:)*GradP(:,n) - GradP(:,n-1) + P(:,n)
    end do

    return    
  end function Polynominal_GenDLegendrePoly

  !> A function to calcuate the Legendre-Gauss-Lobtatto (LGL) points.
  !!
  function Polynominal_GenGaussLobattoPt(Nord) result(pts)
    implicit none

    integer, intent(in) :: Nord
    real(RP) :: pts(Nord+1)

    integer :: N1
    integer :: i

    real(RP) :: xold
    real(RP) :: x(Nord+1)
    real(RP), parameter :: EPS = 1E-16_RP
    real(RP) :: P(1,Nord+1)

    !---------------------------------------------------------------------------

    N1 = Nord+1

    x(1) = -1.0_RP; x(N1) = 1.0_RP
    do i=2, Nord
      x(i) = - cos( PI * dble(i-1)/dble(Nord) )
    end do

    xold = 10.0_RP
    do i=2, Nord
      do while(abs(x(i) - xold) > EPS)
        xold = x(i)
        
        P(:,:) = Polynominal_GenLegendrePoly(Nord, (/ x(i) /))
        x(i) = xold - (x(i) * P(1,N1) - P(1,Nord))/(dble(N1) * P(1,N1))
      end do
    end do

    pts(:) = x(:)
    return   
  end function Polynominal_GenGaussLobattoPt

  !> A function to calcuate the Gauss-Lobbato weights. 
  !!  
  function Polynominal_GenGaussLobattoPtIntWeight(Nord) result(int_weight_lgl)
    implicit none

    integer, intent(in) :: Nord
    real(RP) :: int_weight_lgl(Nord+1)

    real(RP) :: lglPts1D(Nord+1)
    real(RP) :: P1D_ori(Nord+1, Nord+1)
    !---------------------------------------------------------------------------

    lglPts1D(:)      = polynominal_genGaussLobattoPt( Nord )
    P1D_ori(:,:)     = polynominal_genLegendrePoly( Nord, lglPts1D )

    int_weight_lgl(:) = 2.0_RP/(dble(Nord*(Nord+1))*P1D_ori(:,Nord+1)**2)

    return
  end function Polynominal_GenGaussLobattoPtIntWeight

  !> A function to calcuate the Gauss-Legendre points.
  !!
  function Polynominal_GenGaussLegendrePt(Nord) result(pts)
    implicit none

    integer, intent(in) :: Nord
    real(RP) :: pts(Nord)

    integer :: i
    integer :: N1
    real(RP) :: xold
    real(RP) :: x(Nord)
    real(RP), parameter :: EPS = 1E-15_RP
    real(RP) :: P(1,Nord+1)

    !---------------------------------------------------------------------------

    do i=1, Nord
      x(i) = - cos( 0.5_RP*PI * dble(4*i-1)/dble(2*Nord+1) )
    end do

    N1 = Nord+1
    xold = 10.0_RP
    do i=1, Nord
      do while(abs(x(i) - xold) > EPS)
        xold = x(i)
        
        P(:,:) = Polynominal_GenLegendrePoly(Nord, (/ x(i) /))
        x(i) = xold - P(1,N1) * (xold**2 - 1.0_RP)/(dble(Nord)*(xold * P(1,N1) - P(1,Nord)))
      end do
    end do

    pts(:) = x(:)

    return   
  end function Polynominal_GenGaussLegendrePt

  !> A function to calcuate the Gauss-Legendre weights. 
  !!  
  function Polynominal_GenGaussLegendrePtIntWeight(Nord) result(int_weight_lgl)
    implicit none

    integer, intent(in) :: Nord
    real(RP) :: int_weight_lgl(Nord)

    real(RP) :: glPts1D(Nord)
    real(RP) :: P1D_ori(Nord, Nord+1)
    real(RP) :: dP1D_ori(Nord, Nord+1)
    !---------------------------------------------------------------------------

    glPts1D(:)    = polynominal_genGaussLegendrePt( Nord )
    P1D_ori(:,:)    = Polynominal_GenLegendrePoly( Nord, glPts1D )   
    dP1D_ori(:,:) = polynominal_genDLegendrePoly( Nord, glPts1D, P1D_ori )

    int_weight_lgl(:) = 2d0/( (1.0_RP - glPts1D(:)**2) * dP1D_ori(:,Nord+1)**2 )

    return
  end function Polynominal_GenGaussLegendrePtIntWeight

end module scale_polynominal
