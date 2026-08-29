!-------------------------------------------------------------------------------
!> Module common / Coordinate conversion with cubed-sphere projection
!!
!! @par Description
!!      A module to provide coordinate conversions with an equiangular gnomonic cubed-sphere projection
!!
!!
!! @par Reference
!!  - Nair et al. 2015:
!!    A Discontinuous Galerkin Transport Scheme on the Cubed Sphere. 
!!    Monthly Weather Review, 133, 814–828. 
!!    (Appendix A)
!!  - Yin et al. 2017:
!!    Parallel numerical simulation of the thermal convection in the Earth’s outer core on the cubed-sphere. 
!!    Geophysical Journal International, 209, 1934–1954. 
!!    (Appendix A)
!!
!! @author Yuta Kawai, Team SCALE
!!
#include "scaleFElib.h"
module scale_cubedsphere_coord_cnv
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_const, only: &
    PI => CONST_PI,         &
    EPS => CONST_EPS,       &
    RPlanet => CONST_RADIUS
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prc

#ifdef _OPENACC
    use openacc, only: acc_async_sync
#endif

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  public :: CubedSphereCoordCnv_CS2LonLatPos
  public :: CubedSphereCoordCnv_CS2LonLatVec
  public :: CubedSphereCoordCnv_LonLat2CSPos
  public :: CubedSphereCoordCnv_LonLat2CSVec
  public :: CubedSphereCoordCnv_CS2CartPos
  public :: CubedSphereCoordCnv_Cart2CSVec  
  public :: CubedSphereCoordCnv_GetMetric

  interface CubedSphereCoordCnv_CS2LocalOrthVec_alpha
    module procedure :: CubedSphereCoordCnv_CS2LocalOrthVec_alpha_0
    module procedure :: CubedSphereCoordCnv_CS2LocalOrthVec_alpha_1
  end interface
  public :: CubedSphereCoordCnv_CS2LocalOrthVec_alpha
  
  interface CubedSphereCoordCnv_LocalOrth2CSVec_alpha
    module procedure :: CubedSphereCoordCnv_LocalOrth2CSVec_alpha_0
    module procedure :: CubedSphereCoordCnv_LocalOrth2CSVec_alpha_1
  end interface
  public :: CubedSphereCoordCnv_LocalOrth2CSVec_alpha

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !  
  !-----------------------------------------------------------------------------


contains
  !> Calculate longitude and latitude coordinates from local coordinates using the central angles in an equiangular gnomonic cubed-sphere projection
  !!
!OCL SERIAL
  subroutine CubedSphereCoordCnv_CS2LonLatPos( &
    panelID, alpha, beta, gam, Np,             & ! (in)
    lon, lat )                                   ! (out)
    use scale_geographic_coord_cnv, only: &
      GeographicCoordCnv_orth_to_geo_pos
    implicit none

    integer, intent(in) :: panelID     !< Panel ID of cubed-sphere coordinates
    integer, intent(in) :: Np          !< Array size
    real(RP), intent(in) :: alpha(Np)  !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: beta (Np)  !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: gam(Np)    !< A factor of r/a
    real(RP), intent(out) :: lon(Np)   !< Longitude coordinate [rad]
    real(RP), intent(out) :: lat(Np)   !< Latitude coordinate [rad]

    real(RP) :: CartPos(Np,3)
    real(RP) :: GeogPos(Np,3)

    integer :: p
    !-----------------------------------------------------------------------------

    select case( panelID )
    case( 1, 2, 3, 4 )
      !$omp parallel do
      !$acc parallel loop present(alpha, beta, lon, lat)
      do p=1, Np
        lon(p) = alpha(p) + 0.5_RP * PI * dble(panelID - 1)
        lat(p) = atan( tan( beta(p) ) * cos( alpha(p) ) )

        if ( lon(p) < 0.0_RP ) lon(p) = lon(p) + 2.0_RP * PI
      end do
      !$omp end parallel do
    case(5, 6)
      !$acc data create(CartPos, GeogPos) present(alpha, beta, gam, lon, lat)
      call CubedSphereCoordCnv_CS2CartPos( panelID, alpha, beta, gam, Np, & ! (in)
        CartPos(:,1), CartPos(:,2), CartPos(:,3)                          ) ! (out)
      
      call GeographicCoordCnv_orth_to_geo_pos( CartPos(:,:), Np, & ! (in)
        GeogPos(:,:) ) ! (out)
      
      !$omp parallel do
      !$acc parallel loop
      do p=1, Np
        lon(p) = GeogPos(p,1)
        lat(p) = GeogPos(p,2)

        if ( alpha(p) < 0.0_RP ) lon(p) = lon(p) + 2.0_RP * PI
      end do
      !$omp end parallel do
      !$acc end data
    case default
      LOG_ERROR("CubedSphereCoordCnv_CS2LonLatPos",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    return
  end subroutine CubedSphereCoordCnv_CS2LonLatPos

  !> Convert the components of a vector in local coordinates with an equiangular gnomonic cubed-sphere projection to those in longitude and latitude coordinates
  !!
!OCL SERIAL
  subroutine CubedSphereCoordCnv_CS2LonLatVec( &
    panelID, alpha, beta, gam, Np,         & ! (in)
    VecAlpha, VecBeta,                     & ! (in)
    VecLon, VecLat,                        & ! (out)
    lat, gpu_async_id                      ) ! (in, optional)

    use scale_const, only: &
      EPS => CONST_EPS
    implicit none

    integer, intent(in) :: panelID       !< Panel ID of cubed-sphere coordinates
    integer, intent(in) :: Np            !< Array size
    real(RP), intent(in) :: alpha(Np)    !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: beta (Np)    !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: gam(Np)      !< A factor of RPlanet / r
    real(DP), intent(in) :: VecAlpha(Np) !< A component of vector in the alpha-coordinate
    real(DP), intent(in) :: VecBeta (Np) !< A component of vector in the beta-coordinate
    real(RP), intent(out) :: VecLon(Np)  !< A component of vector in the longitude-coordinate
    real(RP), intent(out) :: VecLat(Np)  !< A component of vector in the latitude-coordinate
    real(RP), intent(in), optional :: lat(Np) !< latitude [rad]
    integer, intent(in), optional :: gpu_async_id !< GPU asynchronous ID

    integer :: p
    real(RP) :: X ,Y, del2
    real(RP) :: s

    real(RP) :: radius
    real(RP) :: cos_Lat(Np)

    integer :: gpu_async_id_
    !-----------------------------------------------------------------------------

#ifdef _OPENACC
    if (present(gpu_async_id)) then
      gpu_async_id_ = gpu_async_id
    else
      gpu_async_id_ = acc_async_sync
    end if
#endif

    !$acc data create(cos_Lat) present(alpha, beta, gam, VecLon, VecLat, VecAlpha, VecBeta)

    if (present(lat)) then
      !$omp parallel do
      !$acc parallel loop
      do p=1, Np
        cos_Lat(p) = cos(lat(p))
      end do
    end if

    select case( panelID )
    case( 1, 2, 3, 4 )
      if (.not. present(lat)) then
        !$omp parallel do
        !$acc parallel loop async(gpu_async_id_)
        do p=1, Np
          cos_Lat(p) = cos( atan( tan( beta(p) ) * cos( alpha(p) ) ) )
        end do
      end if

      !$omp parallel do private( X, Y, del2, radius )
      !$acc parallel loop async(gpu_async_id_)
      do p=1, Np
        X = tan( alpha(p) )
        Y = tan( beta (p) )
        del2 = 1.0_RP + X**2 + Y**2
        radius = RPlanet * gam(p)

        VecLon(p) = VecAlpha(p) * cos_Lat(p) * radius
        VecLat(p) = ( - X * Y * VecAlpha(p) + ( 1.0_RP + Y**2 ) * VecBeta(p) ) &
                  * radius * sqrt( 1.0_RP + X**2 ) / del2 
      end do
    case ( 5 )
      s = 1.0_RP
    case ( 6 )
      s = - 1.0_RP
    case default
      LOG_ERROR("CubedSphereCoordCnv_CS2LonLatVec",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    select case( panelID )
    case( 5, 6 )
      !$omp parallel do private( X, Y, del2, radius )
      !$acc parallel loop async(gpu_async_id_)
      do p=1, Np
        X = tan( alpha(p) )
        Y = tan( beta (p) )
        del2 = 1.0_RP + X**2 + Y**2
        radius = s * RPlanet * gam(p)

        if (.not. present(lat)) then
          cos_Lat(p) =  cos( atan( sign(1.0_RP, s) / max( sqrt( X**2 + Y**2 ), EPS ) ) )
        end if

        VecLon(p) = (- Y * ( 1.0_RP + X**2 ) * VecAlpha(p) + X * ( 1.0_RP + Y**2 ) * VecBeta(p) ) &
                  * radius / max( X**2 + Y**2, EPS ) * cos_Lat(p)
        VecLat(p) = (- X * ( 1.0_RP + X**2 ) * VecAlpha(p) - Y * ( 1.0_RP + Y**2 ) * VecBeta(p) ) &
                  * radius / ( del2 * ( max( sqrt( X**2 + Y**2 ), EPS ) ) )
      end do
    end select

    !$acc end data

    return
  end subroutine CubedSphereCoordCnv_CS2LonLatVec

  !> Calculate local coordinates using the central angles in an equiangular gnomonic cubed-sphere projection from longitude and latitude coordinates
  !!
!OCL SERIAL
  subroutine CubedSphereCoordCnv_LonLat2CSPos(  &
    panelID, lon, lat, Np,                 & ! (in)
    alpha, beta                            ) ! (out)

    implicit none
    integer, intent(in) :: panelID      !< Panel ID of cubed-sphere coordinates
    integer, intent(in) :: Np           !< Array size
    real(RP), intent(in) :: lon(Np)     !< Longitude coordinate [rad]
    real(RP), intent(in) :: lat(Np)     !< Latitude coordinate [rad]
    real(RP), intent(out) ::alpha(Np)   !< Local coordinate using the central angles [rad]
    real(RP), intent(out) ::beta (Np)   !< Local coordinate using the central angles [rad]

    integer :: p
    real(RP) :: tan_lat
    real(RP) :: lon_(Np)
    !-----------------------------------------------------------------------------

    !$acc data create(lon_) present(lon, lat, alpha, beta)

    select case(panelID)
    case ( 1, 2, 3, 4 )
      !$omp parallel
      if ( panelID == 1 ) then
        !$omp do
        !$acc parallel loop
        do p=1, Np
          if ( lon(p) >= 2.0_RP * PI - 0.25_RP * PI ) then
            lon_(p) = lon(p) - 2.0_RP * PI
          else
            lon_(p) = lon(p)
          end if
        end do
      else
        !$omp do
        !$acc parallel loop
        do p=1, Np
          if ( lon(p) < 0.0_RP ) then
            lon_(p) = lon(p) + 2.0_RP * PI
          else
            lon_(p) = lon(p)
          end if
        end do
      end if
      !$omp do
      !$acc parallel loop
      do p=1, Np
        alpha(p) = lon_(p) - 0.5_RP * PI * ( dble(panelID) - 1.0_RP )
        beta (p) = atan( tan(lat(p)) / cos(alpha(p)) )
      end do
      !$omp end parallel
    case ( 5 )
      !$omp parallel private(tan_lat)   
      !$omp do
      !$acc parallel loop
      do p=1, Np
        tan_lat = tan(lat(p)) 
        alpha(p) = + atan( sin(lon(p)) / tan_lat )
        beta (p) = - atan( cos(lon(p)) / tan_lat )
      end do
      !$omp end parallel
    case ( 6 )
      !$omp parallel private(tan_lat)    
      !$omp do
      !$acc parallel loop
      do p=1, Np
        tan_lat = tan(lat(p)) 
        alpha(p) = - atan( sin(lon(p)) / tan_lat )
        beta (p) = - atan( cos(lon(p)) / tan_lat )
      end do
      !$omp end parallel
    case default
      LOG_ERROR("CubedSphereCoordCnv_LonLat2CSPos",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    !$acc end data
    return
  end subroutine CubedSphereCoordCnv_LonLat2CSPos

  !> Convert the components of a vector in longitude and latitude coordinates to those in local coordinates with an equiangular gnomonic cubed-sphere projection
  !!
!OCL SERIAL
  subroutine CubedSphereCoordCnv_LonLat2CSVec( &
    panelID, alpha, beta, gam, Np,         & ! (in)
    VecLon, VecLat,                        & ! (in)
    VecAlpha, VecBeta,                     & ! (out)
    lat, gpu_async_id )                      ! (in, optional)

    implicit none

    integer, intent(in) :: panelID        !< Panel ID of cubed-sphere coordinates
    integer, intent(in) :: Np             !< Array size
    real(RP), intent(in) :: alpha(Np)     !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: beta (Np)     !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: gam(Np)       !< A factor of RPlanet / r
    real(DP), intent(in) :: VecLon(Np)    !< A component of vector in the longitude-coordinate
    real(DP), intent(in) :: VecLat(Np)    !< A component of vector in the latitude-coordinate
    real(DP), intent(out) :: VecAlpha(Np) !< A component of vector in the alpha-coordinate
    real(DP), intent(out) :: VecBeta (Np) !< A component of vector in the beta-coordinate
    real(RP), intent(in), optional :: lat(Np)     !< latitude [rad]
    integer, intent(in), optional :: gpu_async_id !< GPU asynchronous ID for OpenACC operations

    integer :: p
    real(RP) :: X ,Y, del2
    real(RP) :: s

    real(RP) :: radius
    real(RP) :: cos_Lat(Np)
    real(RP) :: VecLon_ov_cosLat

    integer :: gpu_async_id_
    !-----------------------------------------------------------------------------
    
#ifdef _OPENACC
    if (present(gpu_async_id)) then
      gpu_async_id_ = gpu_async_id
    else
      gpu_async_id_ = acc_async_sync
    end if
#endif

    !$acc data create(cos_Lat) present(alpha, beta, gam, VecLon, VecLat, VecAlpha, VecBeta)

    if (present(lat)) then
      !$omp parallel do
      !$acc parallel loop async(gpu_async_id_)
      do p=1, Np
        cos_Lat(p) = cos(lat(p))
      end do
    end if

    select case( panelID )
    case( 1, 2, 3, 4 )
      if (.not. present(lat)) then
        !$omp parallel do
        !$acc parallel loop async(gpu_async_id_)
        do p=1, Np
          cos_Lat(p) = cos( atan( tan( beta(p) ) * cos( alpha(p) ) ) )
        end do
      end if
  
      !$omp parallel do private( X, Y, del2, radius, VecLon_ov_cosLat )
      !$acc parallel loop async(gpu_async_id_)
      do p=1, Np
        X = tan( alpha(p) )
        Y = tan( beta (p) )
        del2 = 1.0_RP + X**2 + Y**2
        radius = RPlanet * gam(p)
        VecLon_ov_cosLat = VecLon(p) / cos_Lat(p)

        VecAlpha(p) = VecLon_ov_cosLat / radius
        VecBeta (p)  = ( X * Y * VecLon_ov_cosLat + del2 / sqrt( 1.0_RP + X**2 ) * VecLat(p) ) &
                     / ( radius * (1.0_RP + Y**2) )
      end do
    case ( 5 )
      s = 1.0_RP
    case ( 6 )
      s = -1.0_RP
    case default
      LOG_ERROR("CubedSphereCoordCnv_LonLat2CSVec",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    select case( panelID )
    case( 5, 6 )
      !$omp parallel do private( X, Y, del2, radius, VecLon_ov_cosLat )
      !$acc parallel loop async(gpu_async_id_)
      do p=1, Np
        X = tan( alpha(p) )
        Y = tan( beta (p) )
        del2 = 1.0_RP + X**2 + Y**2
        radius = s * RPlanet * gam(p)

        if (.not. present(lat)) then
          cos_Lat(p) =  cos( atan( s / sqrt(max(X**2 + Y**2, EPS)) ) )
        end if
        VecLon_ov_cosLat = VecLon(p) / cos_Lat(p)

        VecAlpha(p) = (- Y * VecLon_ov_cosLat - del2 * X / sqrt( max(del2 - 1.0_RP,EPS)) * VecLat(p)) &
                    / ( radius * ( 1.0_RP + X**2 ) )
        VecBeta (p) = (  X * VecLon_ov_cosLat - del2 * Y / sqrt( max(del2 - 1.0_RP,EPS)) * VecLat(p)) &
                    / ( radius * ( 1.0_RP + Y**2 ) )
      end do
    end select

    !$acc end data

    return
  end subroutine CubedSphereCoordCnv_LonLat2CSVec

  !> Calculate the Cartesian coordinates from local coordinates using the central angles in an equiangular gnomonic cubed-sphere projection
  !!
!OCL SERIAL
  subroutine CubedSphereCoordCnv_CS2CartPos( &
    panelID, alpha, beta, gam, Np,           & ! (in)
    X, Y, Z                                  ) ! (out)

    implicit none
    integer, intent(in) :: panelID     !< Panel ID of cubed-sphere coordinates
    integer, intent(in) :: Np          !< Array size
    real(RP), intent(in) :: alpha(Np)  !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: beta (Np)  !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: gam(Np)    !< A factor of r/a
    real(RP), intent(out) :: X(Np)     !< x-coordinate in the Cartesian coordinate
    real(RP), intent(out) :: Y(Np)     !< y-coordinate in the Cartesian coordinate
    real(RP), intent(out) :: Z(Np)     !< z-coordinate in the Cartesian coordinate

    integer :: p
    real(RP) :: x1, x2, fac

    !-----------------------------------------------------------------------------

    select case(panelID)
    case(1)
      !$omp parallel do private(x1, x2, fac)
      !$acc parallel loop present(alpha, beta, gam, X, Y, Z)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = RPlanet * gam(p) / sqrt( 1.0_RP + x1**2 + x2**2 )
        X(p) = fac
        Y(p) = fac * x1
        Z(p) = fac * x2
      end do
    case(2)
      !$omp parallel do private(x1, x2, fac)
      !$acc parallel loop present(alpha, beta, gam, X, Y, Z)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = RPlanet * gam(p) / sqrt( 1.0_RP + x1**2 + x2**2 )
        X(p) = - fac * x1 
        Y(p) =   fac
        Z(p) =   fac * x2
      end do
    case(3)
      !$omp parallel do private(x1, x2, fac)
      !$acc parallel loop present(alpha, beta, gam, X, Y, Z)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = RPlanet * gam(p) / sqrt( 1.0_RP + x1**2 + x2**2 )
        X(p) = - fac
        Y(p) = - fac * x1 
        Z(p) =   fac * x2
      end do
    case(4)
      !$omp parallel do private(x1, x2, fac)
      !$acc parallel loop present(alpha, beta, gam, X, Y, Z)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = RPlanet * gam(p) / sqrt( 1.0_RP + x1**2 + x2**2 )
        X(p) =   fac * x1 
        Y(p) = - fac
        Z(p) =   fac * x2
      end do
    case(5)
      !$omp parallel do private(x1, x2, fac)
      !$acc parallel loop present(alpha, beta, gam, X, Y, Z)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = RPlanet * gam(p) / sqrt( 1.0_RP + x1**2 + x2**2 )
        X(p) = - fac * x2
        Y(p) =   fac * x1
        Z(p) =   fac
      end do
    case(6)
      !$omp parallel do private(x1, x2, fac)
      !$acc parallel loop present(alpha, beta, gam, X, Y, Z)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        fac = RPlanet * gam(p) / sqrt( 1.0_RP + x1**2 + x2**2 )
        X(p) =   fac * x2
        Y(p) =   fac * x1
        Z(p) = - fac
      end do
    case default
      LOG_ERROR("CubedSphereCoordCnv_CS2CartPos",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    return
  end subroutine CubedSphereCoordCnv_CS2CartPos

  !> Convert the components of a vector in local coordinates with an equiangular gnomonic cubed-sphere projection to those in the Cartesian coordinates
  !!
!OCL SERIAL
  subroutine CubedSphereCoordCnv_Cart2CSVec( &
    panelID, alpha, beta, gam, Np,         & ! (in)
    Vec_x, Vec_y, Vec_z,                   & ! (in)
    VecAlpha, VecBeta                      ) ! (out)

    implicit none

    integer, intent(in) :: panelID        !< Panel ID of cubed-sphere coordinates
    integer, intent(in) :: Np             !< Array size
    real(RP), intent(in) :: alpha(Np)     !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: beta (Np)     !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: gam(Np)       !< A factor of RPlanet / r
    real(DP), intent(in) :: Vec_x(Np)     !< A component of vector in the x-coordinate with the Cartesian coordinate
    real(DP), intent(in) :: Vec_y(Np)     !< A component of vector in the y-coordinate with the Cartesian coordinate
    real(DP), intent(in) :: Vec_z(Np)     !< A component of vector in the z-coordinate with the Cartesian coordinate
    real(RP), intent(out) :: VecAlpha(Np) !< A component of vector in the alpha-coordinate
    real(RP), intent(out) :: VecBeta (Np) !< A component of vector in the beta-coordinate

    integer :: p
    real(RP) :: x1, x2, fac
    real(RP) :: r_sec2_alpha, r_sec2_beta
    !-----------------------------------------------------------------------------

    select case( panelID )
    case(1)
      !$omp parallel do private(x1, x2, r_sec2_alpha, r_sec2_beta, fac)
      !$acc parallel loop present(alpha, beta, gam, Vec_x, Vec_y, Vec_z, VecAlpha, VecBeta)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        r_sec2_alpha = cos(alpha(p))**2
        r_sec2_beta  = cos(beta (p))**2
        fac = sqrt( 1.0_RP + x1**2 + x2**2 ) / ( RPlanet * gam(p) )

        VecAlpha(p) = fac * r_sec2_alpha * ( - x1  * Vec_x(p) + Vec_y(p) )
        VecBeta (p) = fac * r_sec2_beta  * ( - x2  * Vec_x(p) + Vec_z(p) )
      end do      
    case(2)
      !$omp parallel do private(x1, x2, r_sec2_alpha, r_sec2_beta, fac)
      !$acc parallel loop present(alpha, beta, gam, Vec_x, Vec_y, Vec_z, VecAlpha, VecBeta)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        r_sec2_alpha = cos(alpha(p))**2
        r_sec2_beta  = cos(beta (p))**2
        fac = sqrt( 1.0_RP + x1**2 + x2**2 ) / ( RPlanet * gam(p) )

        VecAlpha(p) = - fac * r_sec2_alpha * ( Vec_x(p) + x1 * Vec_y(p) )
        VecBeta (p) =   fac * r_sec2_beta  * ( - x2  * Vec_y(p) + Vec_z(p) )
      end do
    case(3)
      !$omp parallel do private(x1, x2, r_sec2_alpha, r_sec2_beta, fac)
      !$acc parallel loop present(alpha, beta, gam, Vec_x, Vec_y, Vec_z, VecAlpha, VecBeta)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        r_sec2_alpha = cos(alpha(p))**2
        r_sec2_beta  = cos(beta (p))**2
        fac = sqrt( 1.0_RP + x1**2 + x2**2 ) / ( RPlanet * gam(p) )

        VecAlpha(p) = fac * r_sec2_alpha * ( x1  * Vec_x(p) - Vec_y(p) )
        VecBeta (p) = fac * r_sec2_beta  * ( x2  * Vec_x(p) + Vec_z(p) )
      end do
    case(4)
      !$omp parallel do private(x1, x2, r_sec2_alpha, r_sec2_beta, fac)
      !$acc parallel loop present(alpha, beta, gam, Vec_x, Vec_y, Vec_z, VecAlpha, VecBeta)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        r_sec2_alpha = cos(alpha(p))**2
        r_sec2_beta  = cos(beta (p))**2
        fac = sqrt( 1.0_RP + x1**2 + x2**2 ) / ( RPlanet * gam(p) )

        VecAlpha(p) = fac * r_sec2_alpha * ( Vec_x(p) + x1 * Vec_y(p) )
        VecBeta (p) = fac * r_sec2_beta  * ( x2 * Vec_y(p) + Vec_z(p) )
      end do
    case ( 5 )
      !$omp parallel do private(x1, x2, r_sec2_alpha, r_sec2_beta, fac)
      !$acc parallel loop present(alpha, beta, gam, Vec_x, Vec_y, Vec_z, VecAlpha, VecBeta)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        r_sec2_alpha = cos(alpha(p))**2
        r_sec2_beta  = cos(beta (p))**2
        fac = sqrt( 1.0_RP + x1**2 + x2**2 ) / ( RPlanet * gam(p) )

        VecAlpha(p) =   fac * r_sec2_alpha * ( Vec_y(p) - x1 * Vec_z(p) )
        VecBeta (p) = - fac * r_sec2_beta  * ( Vec_x(p) + x2 * Vec_z(p) )
      end do
    case ( 6 )
      !$omp parallel do private(x1, x2, r_sec2_alpha, r_sec2_beta, fac)
      !$acc parallel loop present(alpha, beta, gam, Vec_x, Vec_y, Vec_z, VecAlpha, VecBeta)
      do p=1, Np
        x1 = tan( alpha(p) )
        x2 = tan( beta (p) )
        r_sec2_alpha = cos(alpha(p))**2
        r_sec2_beta  = cos(beta (p))**2
        fac = sqrt( 1.0_RP + x1**2 + x2**2 ) / ( RPlanet * gam(p) )

        VecAlpha(p) = fac * r_sec2_alpha * ( Vec_y(p) + x1 * Vec_z(p) )
        VecBeta (p) = fac * r_sec2_beta  * ( Vec_x(p) + x2 * Vec_z(p) )
      end do
    case default
      LOG_ERROR("CubedSphereCoordCnv_Cart2CSVec",'(a,i2,a)') "panelID ", panelID, " is invalid. Check!"
      call PRC_abort
    end select

    return
  end subroutine CubedSphereCoordCnv_Cart2CSVec

!OCL SERIAL  
  subroutine CubedSphereCoordCnv_CS2LocalOrthVec_alpha_0( &
    alpha, beta, radius, Np,                      & ! (in)
    VecOrth1, VecOrth2                            ) ! (inout)

    implicit none
    
    integer, intent(in) :: Np
    real(RP), intent(in) :: alpha(Np)
    real(RP), intent(in) :: beta (Np)
    real(RP), intent(in) :: radius(Np)
    real(RP), intent(inout) :: VecOrth1(Np)
    real(RP), intent(inout) :: VecOrth2(Np)

    integer :: p
    real(RP) :: x1, x2, del
    real(RP) :: fac
    real(RP) :: tmp
    !-----------------------------------------------------------------------------

    !$omp parallel do private(x1, x2, del, fac, tmp)
    !$acc parallel loop present(alpha, beta, radius, VecOrth1, VecOrth2)
    do p=1, Np
      x1 = tan( alpha(p) )
      x2 = tan( beta (p) )
      del = sqrt(1.0_RP + x1**2 + x2**2)
      fac = radius(p) * sqrt(1.0_RP + x1**2) / del**2

      tmp = VecOrth1(p)
      VecOrth1(p) = fac * del * tmp
      VecOrth2(p) = fac * ( - x1 * x2 * tmp + (1.0_RP + x2**2) * VecOrth2(p) )
    end do

    return
  end subroutine CubedSphereCoordCnv_CS2LocalOrthVec_alpha_0

!OCL SERIAL  
  subroutine CubedSphereCoordCnv_CS2LocalOrthVec_alpha_1( &
    alpha, beta, radius, Np,                      & ! (in)
    VecAlpha, VecBeta,                            & ! (in)
    VecOrth1, VecOrth2                            ) ! (out)

    implicit none
    
    integer, intent(in) :: Np
    real(RP), intent(in) :: alpha(Np)
    real(RP), intent(in) :: beta (Np)
    real(RP), intent(in) :: radius(Np)
    real(RP), intent(in) :: VecAlpha(Np)
    real(RP), intent(in) :: VecBeta (Np)
    real(RP), intent(out) :: VecOrth1(Np)
    real(RP), intent(out) :: VecOrth2(Np)

    integer :: p
    !-----------------------------------------------------------------------------

    !$omp parallel do
    !$acc parallel loop present(VecAlpha, VecBeta, VecOrth1, VecOrth2)
    do p=1, Np
      VecOrth1(p) = VecAlpha(p)
      VecOrth2(p) = VecBeta(p)
    end do
    !$omp end parallel do

    call CubedSphereCoordCnv_CS2LocalOrthVec_alpha_0( &
      alpha, beta, radius, Np,                      & ! (in)
      VecOrth1, VecOrth2                            ) ! (inout)    

    return
  end subroutine CubedSphereCoordCnv_CS2LocalOrthVec_alpha_1

!OCL SERIAL  
  subroutine CubedSphereCoordCnv_LocalOrth2CSVec_alpha_0( &
    alpha, beta, radius, Np,                      & ! (in)
    VecAlpha, VecBeta                             ) ! (inout)

    implicit none
    
    integer, intent(in) :: Np
    real(RP), intent(in) :: alpha(Np)
    real(RP), intent(in) :: beta (Np)
    real(RP), intent(in) :: radius(Np)
    real(RP), intent(inout) :: VecAlpha(Np)
    real(RP), intent(inout) :: VecBeta (Np)    

    integer :: p
    real(RP) :: x1, x2, del
    real(RP) :: fac
    real(RP) :: tmp
    !-----------------------------------------------------------------------------

    !$omp parallel do private(x1, x2, del, fac, tmp)
    !$acc parallel loop present(alpha, beta, radius, VecAlpha, VecBeta)
    do p=1, Np
      x1 = tan( alpha(p) )
      x2 = tan( beta (p) )
      del = sqrt(1.0_RP + x1**2 + x2**2)
      fac = del / ( radius(p) * ( 1.0_RP + x2**2 ) * sqrt( 1.0_RP + x1**2 ) )

      tmp = VecAlpha(p)
      VecAlpha(p) = fac * ( 1.0_RP + x2**2 ) * tmp
      VecBeta (p) = fac * ( x1 * x2 * tmp + del * VecBeta(p) )
    end do

    return
  end subroutine CubedSphereCoordCnv_LocalOrth2CSVec_alpha_0

!OCL SERIAL  
  subroutine CubedSphereCoordCnv_LocalOrth2CSVec_alpha_1( &
    alpha, beta, radius, Np,                      & ! (in)
    VecOrth1, VecOrth2,                           & ! (in)
    VecAlpha, VecBeta                             ) ! (out)
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(in) :: alpha(Np)
    real(RP), intent(in) :: beta (Np)
    real(RP), intent(in) :: radius(Np)
    real(RP), intent(in) :: VecOrth1(Np)
    real(RP), intent(in) :: VecOrth2(Np)
    real(RP), intent(out) :: VecAlpha(Np)
    real(RP), intent(out) :: VecBeta (Np)    

    integer :: p
    !-----------------------------------------------------------------------------

    !$omp parallel do
    !$acc parallel loop present(VecOrth1, VecOrth2, VecAlpha, VecBeta)
    do p=1, Np
      VecAlpha(p) = VecOrth1(p)
      VecBeta (p) = VecOrth2(p)
    end do
    !$omp end parallel do

    call CubedSphereCoordCnv_LocalOrth2CSVec_alpha_0( &
      alpha, beta, radius, Np,                        & ! (in)
      VecAlpha, VecBeta                               ) ! (inout)    

    return
  end subroutine CubedSphereCoordCnv_LocalOrth2CSVec_alpha_1

  !> Calculate the metrics associated with an equiangular gnomonic cubed-sphere projection to those in longitude and latitude coordinates
  !!
!OCL SERIAL  
  subroutine CubedSphereCoordCnv_GetMetric( &
    alpha, beta, Np, radius,            & ! (in)
    G_ij, GIJ, Gsqrt                    ) ! (out)

    implicit none

    integer, intent(in) :: Np             !< Array size
    real(RP), intent(in) :: alpha(Np)     !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: beta (Np)     !< Local coordinate using the central angles [rad]
    real(RP), intent(in) :: radius        !< Planetary radius
    real(RP), intent(out) :: G_ij(Np,2,2) !< Horizontal covariant metric tensor
    real(RP), intent(out) :: GIJ (Np,2,2) !< Horizontal contravariant metric tensor
    real(RP), intent(out) :: Gsqrt(Np)    !< Horizontal Jacobian

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
    !$acc parallel loop private(G_ij_, GIJ_) present(alpha, beta, G_ij, GIJ, Gsqrt)
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
  end subroutine CubedSphereCoordCnv_GetMetric

end module scale_cubedsphere_coord_cnv
