!-------------------------------------------------------------------------------
!> module common / Coordinate conversion with a cubed-sphere 
!!
!! @par Description
!!      Coordinate conversion with a cubed-sphere 
!!
!! @author Team SCALE
!!
#include "scaleFElib.h"
module scale_cubedsphere_cnv
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_const, only: &
    PI => CONST_PI,      &
    EPS => CONST_EPS
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prc

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  public :: CubedSphereCnv_CS2LonLatCoord
  public :: CubedSphereCnv_CS2CartCoord
  public :: CubedSphereCnv_CS2LonLatVec
  public :: CubedSphereCnv_LonLat2CSPos
  public :: CubedSphereCnv_LonLat2CSVec
  public :: CubedSphereCnv_GetMetric

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !  
  !-----------------------------------------------------------------------------


contains

  subroutine CubedSphereCnv_CS2LonLatCoord( &
    panelID, alpha, beta, Np, radius,        & ! (in)
    lon, lat )                                 ! (out)
    implicit none

    integer, intent(in) :: panelID
    integer, intent(in) :: Np
    real(RP), intent(in) :: alpha(Np)
    real(RP), intent(in) :: beta (Np)
    real(RP), intent(in) :: radius
    real(RP), intent(out) :: lon(Np)
    real(RP), intent(out) :: lat(Np)

    integer :: p
    !-----------------------------------------------------------------------------

    select case( panelID )
    case( 1, 2, 3, 4 )
      !$omp parallel 
      !$omp do
      do p=1, Np
        lon(p) = alpha(p) + 0.5_RP * PI * dble(panelID - 1)
        lat(p) = atan( tan( beta(p) ) * cos( alpha(p) ) )
      end do
      !$omp workshare
      where( lon(:) < 0.0_RP )
        lon(:) = lon(:) + 2.0_RP * PI
      end where
      !$omp end workshare
      !$omp end parallel
    case(5)
      !$omp parallel
      !$omp do
      do p=1, Np
        lon(p) = - atan( tan( alpha(p) ) / tan( sign(max(abs(beta(p)),EPS), beta(p)) ) )
        lat(p) = + atan( 1.0_RP / sqrt( max(tan(alpha(p))**2 + tan(beta(p))**2, EPS) ) )
      end do
      !$omp end do
      !$omp workshare
      where( beta(:) >= 0.0_RP )
        lon(:) = lon(:) + PI
      end where
      where( alpha < 0.0_RP .and. beta(:) < 0.0_RP )
        lon(:) = lon(:) + 2.0_RP * PI
      end where
      !$omp end workshare
      !$omp end parallel
    case(6)
      !$omp parallel
      !$omp do
      do p=1, Np
        lon(p) = + atan( tan( alpha(p) ) / tan( sign(max(abs(beta(p)),EPS), beta(p)) ) )
        lat(p) = - atan( 1.0_RP / sqrt( max(tan(alpha(p))**2 + tan(beta(p))**2, EPS) ) )
      end do
      !$omp end do
      !$omp workshare
      where( beta(:) < 0.0_RP )
        lon(:) = lon(:) + PI
      end where
      where( alpha(:) < 0.0_RP .and. beta(:) >= 0.0_RP )
        lon(:) = lon(:) + 2.0_RP * PI
      end where
      !$omp end workshare
      !$omp end parallel
    case default
      LOG_ERROR("CubedSphereUtil_CS2LonLatCoord",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    return
  end subroutine CubedSphereCnv_CS2LonLatCoord

  !
  !
  subroutine CubedSphereCnv_LonLat2CSVec( &
    panelID, alpha, beta, Np,  radius,     & ! (in)
    VecLon, VecLat,                        & ! (in)
    VecAlpha, VecBeta                      ) ! (out)

    implicit none

    integer, intent(in) :: panelID
    integer, intent(in) :: Np
    real(RP), intent(in) :: alpha(Np)
    real(RP), intent(in) :: beta (Np)
    real(RP), intent(in) :: radius
    real(DP), intent(in) :: VecLon(Np)
    real(DP), intent(in) :: VecLat(Np)
    real(RP), intent(out) :: VecAlpha(Np)
    real(RP), intent(out) :: VecBeta (Np)

    integer :: p
    real(RP) :: X ,Y, del2
    real(RP) :: s
    !-----------------------------------------------------------------------------

    select case( panelID )
    case( 1, 2, 3, 4 )
      !$omp parallel do private( X, Y, del2 )
      do p=1, Np
        X = tan( alpha(p) )
        Y = tan( beta (p) )
        del2 = 1.0_RP + X**2 + Y**2

        VecAlpha(p) = VecLon(p) / radius
        VecBeta (p)  = ( X * Y * VecLon(p) + del2 / sqrt( 1.0_RP + X**2 ) * VecLat(p) ) &
                     / ( radius * (1.0_RP + Y**2) )
      end do
    case ( 5 )
      s = 1.0_RP
    case ( 6 )
      s = -1.0_RP
    case default
      LOG_ERROR("CubedSphereUtil_LonLat2CSVec",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    select case( panelID )
    case( 5, 6 )
      !$omp parallel do private( X, Y, del2 )
      do p=1, Np
        X = tan( alpha(p) )
        Y = tan( beta (p) )
        del2 = 1.0_RP + X**2 + Y**2

        VecAlpha(p) = (- Y * VecLon(p) - del2 * X / sqrt( max(del2 - 1.0_RP,EPS)) * VecLat(p)) * s &
                    / ( radius * ( 1.0_RP + X**2 ) )
        VecBeta (p) = (  X * VecLon(p) - del2 * Y / sqrt( max(del2 - 1.0_RP,EPS)) * VecLat(p)) * s &
                    / ( radius * ( 1.0_RP + Y**2 ) )
      end do
    end select

    return
  end subroutine CubedSphereCnv_LonLat2CSVec

  subroutine CubedSphereCnv_CS2LonLatVec( &
    panelID, alpha, beta, Np, radius,     & ! (in)
    VecAlpha, VecBeta,                    & ! (in)
    VecLon, VecLat                        ) ! (out)

    use scale_const, only: &
      EPS => CONST_EPS
    implicit none

    integer, intent(in) :: panelID
    integer, intent(in) :: Np
    real(RP), intent(in) :: alpha(Np)
    real(RP), intent(in) :: beta (Np)
    real(RP), intent(in) :: radius
    real(DP), intent(in) :: VecAlpha(Np)
    real(DP), intent(in) :: VecBeta (Np)
    real(RP), intent(out) :: VecLon(Np)
    real(RP), intent(out) :: VecLat(Np)

    integer :: p
    real(RP) :: X ,Y, del2
    real(RP) :: s
    !-----------------------------------------------------------------------------

    select case( panelID )
    case( 1, 2, 3, 4 )
      !$omp parallel do private( X, Y, del2 )
      do p=1, Np
        X = tan( alpha(p) )
        Y = tan( beta (p) )
        del2 = 1.0_RP + X**2 + Y**2

        VecLon(p) = VecAlpha(p) * radius
        VecLat(p) = ( - X * Y * VecAlpha(p) + ( 1.0_RP + Y**2 ) * VecBeta(p) ) &
                  * radius * sqrt( 1.0_RP + X**2 ) / del2 
      end do
    case ( 5 )
      s = radius
    case ( 6 )
      s = - radius
    case default
      LOG_ERROR("CubedSphereUtil_LonLat2CSVec",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    select case( panelID )
    case( 5, 6 )
      !$omp parallel do private( X, Y, del2 )
      do p=1, Np
        X = tan( alpha(p) )
        Y = tan( beta (p) )
        del2 = 1.0_RP + X**2 + Y**2

        VecLon(p) = (- Y * ( 1.0 + X**2 ) * VecAlpha(p) + X * ( 1.0_RP + Y**2 ) * VecBeta(p) ) &
                  * s / ( X**2 + Y**2 + EPS )
        VecLat(p) = (- X * ( 1.0 + X**2 ) * VecAlpha(p) - Y * ( 1.0_RP + Y**2 ) * VecBeta(p) ) &
                  * s / ( del2 * sqrt( X**2 + Y**2 ) + EPS )
      end do
    end select

    return
  end subroutine CubedSphereCnv_CS2LonLatVec

  subroutine CubedSphereCnv_CS2CartCoord( &
    panelID, alpha, beta, Np, radius,     & ! (in)
    X, Y, Z                               ) ! (out)

    implicit none
    integer, intent(in) :: panelID
    integer, intent(in) :: Np
    real(RP), intent(in) :: alpha(Np)
    real(RP), intent(in) :: beta (Np)
    real(RP), intent(in) :: radius
    real(RP), intent(out) :: X(Np)
    real(RP), intent(out) :: Y(Np)
    real(RP), intent(out) :: Z(Np)

    integer :: p
    real(DP) :: a
    real(RP) :: x1, x2
    real(RP) :: fac
    !-----------------------------------------------------------------------------

    a = 1.0_RP / sqrt(3.0_RP)

    select case(panelID)
    case(1)
      !$omp parallel do private(x1, x2, fac)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = 1.0_RP / sqrt( a**2 + x1**2 + x2**2 )
        X(p) = fac * a
        Y(p) = fac * x1
        Z(p) = fac * x2
      end do
    case(2)
      !$omp parallel do private(x1, x2, fac)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = 1.0_RP / sqrt( a**2 + x1**2 + x2**2 )
        X(p) = - fac * x1 
        Y(p) =   fac * a
        Z(p) =   fac * x2
      end do
    case(3)
      !$omp parallel do private(x1, x2, fac)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = 1.0_RP / sqrt( a**2 + x1**2 + x2**2 )
        X(p) = - fac * a
        Y(p) = - fac * x1 
        Z(p) =   fac * x2
      end do
    case(4)
      !$omp parallel do private(x1, x2, fac)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = 1.0_RP / sqrt( a**2 + x1**2 + x2**2 )
        X(p) =   fac * x1 
        Y(p) = - fac * a
        Z(p) =   fac * x2
      end do
    case(5)
      !$omp parallel do private(x1, x2, fac)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = 1.0_RP / sqrt( a**2 + x1**2 + x2**2 )
        X(p) = - fac * x2
        Y(p) =   fac * x1
        Z(p) =   fac * a
      end do
    case(6)
      !$omp parallel do private(x1, x2, fac)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = 1.0_RP / sqrt( a**2 + x1**2 + x2**2 )
        X(p) =   fac * x2
        Y(p) =   fac * x1
        Z(p) = - fac * a
      end do
    case default
      LOG_ERROR("CubedSphereUtil_CS2CartCoord",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    return
  end subroutine CubedSphereCnv_CS2CartCoord

  subroutine CubedSphereCnv_LonLat2CSPos(  &
    panelID, lon, lat, Np,                 & ! (in)
    alpha, beta                            ) ! (out)

    implicit none

    integer, intent(in) :: panelID
    integer, intent(in) :: Np
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(out) ::alpha(Np)
    real(RP), intent(out) ::beta (Np)

    integer :: p
    real(RP) :: tan_lat
    real(RP) :: lon_(Np)
    !-----------------------------------------------------------------------------

    select case(panelID)
    case ( 1, 2, 3, 4 )
      !$omp parallel
      if ( panelID == 1 ) then
        !$omp workshare
        where (lon(:) >= 2.0_RP * PI - 0.25_RP * PI )
          lon_(:) = lon(:) - 2.0_RP * PI
        elsewhere
          lon_(:) = lon(:)
        end where
        !$omp end workshare
      else
        !$omp workshare
        where (lon(:) < 0.0_RP )
          lon_(:) = lon(:) + 2.0_RP * PI
        elsewhere
          lon_(:) = lon(:)
        end where    
        !$omp end workshare
      end if
      !$omp do
      do p=1, Np
        alpha(p) = lon_(p) - 0.5_RP * PI * ( dble(panelID) - 1.0_RP )
        beta (p) = atan( tan(lat(p)) / cos(alpha(p)) )
      end do
      !$omp end parallel
    case ( 5 )
      !$omp parallel private(tan_lat)   
      !$omp do
      do p=1, Np
        tan_lat = tan(lat(p)) 
        alpha(p) = + atan( sin(lon(p)) / tan_lat )
        beta (p) = - atan( cos(lon(p)) / tan_lat )
      end do
      !$omp end parallel
    case ( 6 )
      !$omp parallel private(tan_lat)    
      !$omp do
      do p=1, Np
        tan_lat = tan(lat(p)) 
        alpha(p) = - atan( sin(lon(p)) / tan_lat )
        beta (p) = - atan( cos(lon(p)) / tan_lat )
      end do
      !$omp end parallel
    case default
      LOG_ERROR("CubedSphereUtil_LonLat2CSPos",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    return
  end subroutine CubedSphereCnv_LonLat2CSPos

  subroutine CubedSphereCnv_GetMetric( &
    alpha, beta, Np, radius,            & ! (in)
    G_ij, GIJ, Gsqrt                    ) ! (out)

    implicit none

    integer, intent(in) :: Np
    real(RP), intent(in) :: alpha(Np)
    real(RP), intent(in) :: beta(Np)
    real(RP), intent(in) :: radius
    real(RP), intent(out) :: G_ij(Np,2,2)
    real(RP), intent(out) :: GIJ (Np,2,2)
    real(RP), intent(out) :: Gsqrt(Np)

    real(RP) :: X, Y
    real(RP) :: r2
    real(RP) :: fac
    real(RP) :: G_ij_(2,2)
    real(RP) :: GIJ_ (2,2)

    integer :: p
    real(RP) :: OnePlusX2, OnePlusY2
    !-----------------------------------------------------------------------------

    !$omp parallel do private( &
    !$omp X, Y, r2, OnePlusX2, OnePlusY2, fac,    &
    !$omp G_ij_, GIJ_                             )
    do p=1, Np
      X = tan(alpha(p))
      Y = tan(beta (p))
      r2 = 1.0_RP + X**2 + Y**2
      OnePlusX2 = 1.0_RP + X**2
      OnePlusY2 = 1.0_RP + Y**2

      fac = OnePlusX2 * OnePlusY2 * ( radius / r2 )**2
      G_ij_(1,1) =   fac * OnePlusX2
      G_ij_(1,2) = - fac * (X * Y)
      G_ij_(2,1) = - fac * (X * Y)
      G_ij_(2,2) =   fac * OnePlusY2
      G_ij(p,:,:) = G_ij_(:,:)

      Gsqrt(p) = radius**2 * OnePlusX2 * OnePlusY2 / ( r2 * sqrt(r2) )

      fac = 1.0_RP / Gsqrt(p)**2
      GIJ_(1,1) =   fac * G_ij_(2,2)
      GIJ_(1,2) = - fac * G_ij_(1,2)
      GIJ_(2,1) = - fac * G_ij_(2,1)
      GIJ_(2,2) =   fac * G_ij_(1,1)
      GIJ(p,:,:) = GIJ_(:,:)
    end do

    return
  end subroutine CubedSphereCnv_GetMetric

end module scale_cubedsphere_cnv
