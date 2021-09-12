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
  public :: Polynominal_GenLegendrePoly_sub
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
!OCL SERIAL
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
!OCL SERIAL
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
!OCL SERIAL
  subroutine Polynominal_GenLegendrePoly_sub(Nord, x, P)
    implicit none

    integer, intent(in) :: Nord
    real(RP), intent(in) :: x(:)
    real(RP), intent(out) :: P(size(x), Nord+1)

    integer :: n

    !---------------------------------------------------------------------------

    if (Nord < 0) then
       LOG_ERROR("Polynominal_GenLegendrePoly",*)   "Nord must be larger than 0."
    end if

    P(:,1) = 1.0_RP
    if (Nord==0) return

    P(:,2) = x(:)
    do n=2, Nord
      P(:,n+1) = ( dble(2*n-1)*x(:)*P(:,n) - dble(n-1)*P(:,n-1) )/dble(n)
    end do
    
    return    
  end subroutine Polynominal_GenLegendrePoly_sub

  !> A function to obtain the values of Legendre polynominals which are evaluated at aribitary points. 
  !!
!OCL SERIAL
  function Polynominal_GenLegendrePoly(Nord, x) result(P)
    implicit none

    integer, intent(in) :: Nord
    real(RP), intent(in) :: x(:)
    real(RP) :: P(size(x), Nord+1)
    !---------------------------------------------------------------------------

    call Polynominal_GenLegendrePoly_sub( Nord, x(:), & ! (in)
       P(:,:)                                         ) ! (out)
    return    
  end function Polynominal_GenLegendrePoly

  !> A function to obtain differential values of Legendre polynominals which are evaluated at aribitary points. 
  !! 
!OCL SERIAL
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
!OCL SERIAL
  function Polynominal_GenGaussLobattoPt(Nord) result(pts)
    implicit none

    integer, intent(in) :: Nord
    real(RP) :: pts(Nord+1)

    integer :: N1
    !---------------------------------------------------------------------------

    pts(1) = -1.0_RP; pts(Nord+1) = 1.0_RP
    if (Nord==1) return

    call gen_JacobiGaussQuadraturePts( 1, 1, Nord-2, pts(2:Nord) )
    return   
  end function Polynominal_GenGaussLobattoPt

  !> A function to calcuate the Gauss-Lobbato weights. 
  !!  
!OCL SERIAL
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
!OCL SERIAL
  function Polynominal_GenGaussLegendrePt(Nord) result(pts)
    implicit none

    integer, intent(in) :: Nord
    real(RP) :: pts(Nord)
    !---------------------------------------------------------------------------

    call gen_JacobiGaussQuadraturePts( 0, 0, Nord-1, pts(:) )
    return   
  end function Polynominal_GenGaussLegendrePt

  !> A function to calcuate the Gauss-Legendre weights. 
  !!  
!OCL SERIAL
  function Polynominal_GenGaussLegendrePtIntWeight(Nord) result(int_weight_gl)
    implicit none

    integer, intent(in) :: Nord
    real(RP) :: int_weight_gl(Nord)

    real(RP) :: glPts1D(Nord)
    real(RP) :: P1D_ori(Nord, Nord+1)
    real(RP) :: dP1D_ori(Nord, Nord+1)
    !---------------------------------------------------------------------------

    glPts1D(:)    = polynominal_genGaussLegendrePt( Nord )
    P1D_ori(:,:)    = Polynominal_GenLegendrePoly( Nord, glPts1D )   
    dP1D_ori(:,:) = polynominal_genDLegendrePoly( Nord, glPts1D, P1D_ori )

    int_weight_gl(:) = 2.0_RP / ( (1.0_RP - glPts1D(:)**2) * dP1D_ori(:,Nord+1)**2 )

    return
  end function Polynominal_GenGaussLegendrePtIntWeight

  !- private -------------------------------

  !> Calculate the N'th-order Gauss quadrature points and weights associated the Jacobi polynomial of type (alpja,beta).
!OCL SERIAL
  subroutine gen_JacobiGaussQuadraturePts( alpha, beta, N, &
      x )

    implicit none
    integer, intent(in) :: alpha
    integer, intent(in) :: beta
    integer, intent(in) :: N
    real(RP), intent(out) :: x(N+1)

    integer :: i
    real(DP) :: d(N+1), e(N)
    real(DP) :: work(2*(N+1)-2), z(N+1,N+1)
    real(DP) :: h1(N+1)
    integer :: info
    !--------------------------------------------------------------

    if (N==0) then
      x(1) = - dble(alpha - beta) / dble(alpha + beta + 2)
      return
    end if

    do i=0, N
      h1(i+1) = dble( 2*i + alpha + beta )
    end do

    do i=1, N+1
      d(i) = - dble(alpha**2 - beta**2) / (h1(i) * (h1(i) + 2D0))
    end do
    do i=1, N
      e(i) = 2D0 / (h1(i) + 2D0) &
           * sqrt( dble(i * (i + alpha + beta) * (i + alpha) * (i + beta)) &
                   / ((h1(i) + 1D0) * (h1(i) + 3D0))                       )
    end do
    if ( dble(alpha + beta) < 1D-16) d(1) = 0.0_RP
    
    call dstev( 'Vectors', N+1, d, e, z, N+1, work, info )
    x(:) = d(:) 

    return
  end subroutine gen_JacobiGaussQuadraturePts

end module scale_polynominal
