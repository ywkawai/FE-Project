!> module common / Coordinate conversion with a geographic coordinate
!!
!! @par Description
!!      Coordinate conversion with a geographic coordinate
!!
!! @author Team SCALE
!!
#include "scaleFElib.h"
module scale_geographic_coord_cnv
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
  public :: GeographicCoordCnv_orth_to_geo_pos
  public :: GeographicCoordCnv_geo_to_orth_pos

  public :: GeographicCoordCnv_orth_to_geo_vec
  public :: GeographicCoordCnv_geo_to_orth_vec

  public :: GeographicCoordCnv_rotateX
  public :: GeographicCoordCnv_rotateY
  public :: GeographicCoordCnv_rotateZ

contains
!OCL SERIAL
  subroutine GeographicCoordCnv_orth_to_geo_pos( orth_p, Np, &
    geo_p )
    
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(in) :: orth_p(Np,3)
    real(RP), intent(out) :: geo_p(Np,3)

    integer :: p
    !---------------------------------------------------------  

    !$omp parallel private(p)
    !$omp do
    do p=1, Np
      geo_p(p,3) = sqrt( orth_p(p,1)**2 + orth_p(p,2)**2 + orth_p(p,3)**2 )
      geo_p(p,2) = asin( orth_p(p,3) / geo_p(p,3) )
      geo_p(p,1) = atan( orth_p(p,2) / orth_p(p,1) )
    end do
    !$omp do
    do p=1, Np
      if ( geo_p(p,1) <= 0.0_RP .and. orth_p(p,1) < 0.0_RP ) then
        geo_p(p,1) = geo_p(p,1) + PI
      else if ( geo_p(p,1) >= 0.0_RP .and. orth_p(p,1) < 0.0_RP ) then
        geo_p(p,1) = geo_p(p,1) - PI
      end if
      if ( orth_p(p,1) == 0.0_RP .and. orth_p(p,2) == 0.0_RP ) then
        geo_p(p,1) = 0.0_RP
      else if ( orth_p(p,1) == 0.0_RP ) then
        geo_p(p,1) = sign(1.0_RP, orth_p(p,2)) * 0.5_RP * PI
      end if   
    end do
    !$omp end parallel

    return
  end subroutine GeographicCoordCnv_orth_to_geo_pos

!OCL SERIAL
  subroutine GeographicCoordCnv_geo_to_orth_pos( geo_p, Np, &
    orth_p )

    implicit none
    integer, intent(in) :: Np
    real(RP), intent(in) :: geo_p(Np,3)
    real(RP), intent(out) :: orth_p(Np,3)

    integer :: p
    !---------------------------------------------------------  

    !$omp parallel do private(p)
    do p=1, Np
      orth_p(p,1) = geo_p(p,3) * cos(geo_p(p,2)) * cos(geo_p(p,1))
      orth_p(p,2) = geo_p(p,3) * cos(geo_p(p,2)) * sin(geo_p(p,1))
      orth_p(p,3) = geo_p(p,3) * sin(geo_p(p,2))
    end do

    return
  end subroutine GeographicCoordCnv_geo_to_orth_pos

!OCL SERIAL
  subroutine GeographicCoordCnv_orth_to_geo_vec( orth_v, geo_p, Np, &
    geo_v ) 
    implicit none

    integer, intent(in) :: Np
    real(RP), intent(in) :: orth_v(Np,3)
    real(RP), intent(in) :: geo_p(Np,3)
    real(RP), intent(out) :: geo_v(Np,3)

    integer :: p
    real(RP) :: sin_geo_p1
    real(RP) :: cos_geo_p1
    real(RP) :: sin_geo_p2
    real(RP) :: cos_geo_p2    
    !---------------------------------------------------------  

    !$omp parallel do private(p, cos_geo_p2)
    do p=1, Np
      sin_geo_p1 = sin(geo_p(p,1))
      cos_geo_p1 = cos(geo_p(p,1))
      sin_geo_p2 = sin(geo_p(p,2))
      cos_geo_p2 = cos(geo_p(p,2))

      geo_v(p,3) =   orth_v(p,1) * cos_geo_p2 * cos_geo_p1  &
                   + orth_v(p,2) * cos_geo_p2 * sin_geo_p1 &
                   + orth_v(p,3) * sin_geo_p2
  
      geo_v(p,2) = - orth_v(p,1) * sin_geo_p2 * cos_geo_p1 &
                   - orth_v(p,2) * sin_geo_p2 * sin_geo_p1 &
                   + orth_v(p,3) * cos_geo_p2
    
      geo_v(p,1) = - orth_v(p,1) * sin_geo_p1 &
                   + orth_v(p,2) * cos_geo_p1
    end do

    return
  end subroutine GeographicCoordCnv_orth_to_geo_vec

!OCL SERIAL
  subroutine GeographicCoordCnv_geo_to_orth_vec( geo_v, geo_p, Np, &
    orth_v ) 
    implicit none

    integer, intent(in) :: Np
    real(RP), intent(in) :: geo_v(Np,3)
    real(RP), intent(in) :: geo_p(Np,3)
    real(RP), intent(out) :: orth_v(Np,3)

    integer :: p
    real(RP) :: sin_geo_p1
    real(RP) :: cos_geo_p1
    real(RP) :: sin_geo_p2
    real(RP) :: cos_geo_p2    
    !---------------------------------------------------------  

    !$omp parallel do private(p, sin_geo_p1, cos_geo_p1, sin_geo_p2, cos_geo_p2)
    do p=1, Np
      sin_geo_p1 = sin(geo_p(p,1))
      cos_geo_p1 = cos(geo_p(p,1))
      sin_geo_p2 = sin(geo_p(p,2))
      cos_geo_p2 = cos(geo_p(p,2))

      orth_v(p,1) =   geo_v(p,3) * cos_geo_p1 * cos_geo_p2 &
                    - geo_v(p,2) * cos_geo_p1 * sin_geo_p2 &
                    - geo_v(p,1) * sin_geo_p1

      orth_v(p,2) =   geo_v(p,3) * sin_geo_p1 * cos_geo_p2 &
                    - geo_v(p,2) * sin_geo_p1 * sin_geo_p2 &
                    + geo_v(p,1) * cos_geo_p1

      orth_v(p,3) = geo_v(p,3) * sin_geo_p2 + geo_v(p,2) * cos_geo_p2                  
    end do

    return
  end subroutine GeographicCoordCnv_geo_to_orth_vec

!OCL SERIAL
  subroutine GeographicCoordCnv_rotateX( pos_vec, angle, &
    rotated_pos_vec )

    implicit none
    real(RP), intent(in) :: pos_vec(3)
    real(RP), intent(in) :: angle
    real(RP), intent(out) :: rotated_pos_vec(3)

    real(RP) :: sin_angle, cos_angle
    !---------------------------------------------------------  

    sin_angle = sin(angle)
    cos_angle = cos(angle)

    rotated_pos_vec(1) = pos_vec(1)
    rotated_pos_vec(2) = cos_angle * pos_vec(2) - sin_angle * pos_vec(3)
    rotated_pos_vec(3) = sin_angle * pos_vec(2) + cos_angle * pos_vec(3)

    return
  end subroutine GeographicCoordCnv_rotateX

!OCL SERIAL
  subroutine GeographicCoordCnv_rotateY( pos_vec, angle, &
    rotated_pos_vec )
    
    implicit none
    real(RP), intent(in) :: pos_vec(3)
    real(RP), intent(in) :: angle
    real(RP), intent(out) :: rotated_pos_vec(3)

    real(RP) :: sin_angle, cos_angle
    !---------------------------------------------------------  

    sin_angle = sin(angle)
    cos_angle = cos(angle)

    rotated_pos_vec(1) = cos_angle * pos_vec(1) + sin_angle * pos_vec(3)
    rotated_pos_vec(2) = pos_vec(2)
    rotated_pos_vec(3) = - sin_angle * pos_vec(1) + cos_angle * pos_vec(3)

    return
  end subroutine GeographicCoordCnv_rotateY

!OCL SERIAL
  subroutine GeographicCoordCnv_rotateZ( pos_vec, angle, &
    rotated_pos_vec )

    implicit none
    real(RP), intent(in) :: pos_vec(3)
    real(RP), intent(in) :: angle
    real(RP), intent(out) :: rotated_pos_vec(3)

    real(RP) :: sin_angle, cos_angle
    !---------------------------------------------------------  

    sin_angle = sin(angle)
    cos_angle = cos(angle)

    rotated_pos_vec(1) = cos_angle * pos_vec(1) - sin_angle * pos_vec(2)
    rotated_pos_vec(2) = sin_angle * pos_vec(1) + cos_angle * pos_vec(2)
    rotated_pos_vec(3) = pos_vec(3)

    return
  end subroutine GeographicCoordCnv_rotateZ
  
end  module scale_geographic_coord_cnv
