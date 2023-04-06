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
  subroutine GeographicCoordCnv_orth_to_geo_pos( geo_p, &
    orth_p, Np )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: geo_p(3,Np)
    real(RP), intent(in) :: orth_p(3,Np)

    integer :: p
    !---------------------------------------------------------  

    !$omp parallel private(p)
    !$omp do
    do p=1, Np
      geo_p(3,p) = sqrt( orth_p(1,p)**2 + orth_p(2,p)**2 + orth_p(3,p)**2 )
      geo_p(2,p) = asin( orth_p(3,p) / geo_p(3,p) )
      geo_p(1,p) = atan( orth_p(2,p) / orth_p(1,p) )
    end do
    !$omp do
    do p=1, Np
      if ( geo_p(1,p) < 0 .and. orth_p(1,p) < 0 ) then
        geo_p(1,p) = geo_p(1,p) + PI
      else if ( geo_p(1,p) > 0 .and. orth_p(1,p) < 0 ) then
        geo_p(1,p) = geo_p(1,p) - PI
      end if      
    end do
    !$omp end parallel

    return
  end subroutine GeographicCoordCnv_orth_to_geo_pos

!OCL SERIAL
  subroutine GeographicCoordCnv_geo_to_orth_pos( orth_p, &
    geo_p, Np )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: orth_p(3,Np)
    real(RP), intent(in) :: geo_p(3,Np)

    integer :: p
    !---------------------------------------------------------  

    !$omp parallel do private(p)
    do p=1, Np
      orth_p(1,p) = geo_p(3,p) * cos(geo_p(2,p)) * cos(geo_p(1,p))
      orth_p(2,p) = geo_p(3,p) * cos(geo_p(2,p)) * sin(geo_p(1,p))
      orth_p(3,p) = geo_p(3,p) * sin(geo_p(2,p))
    end do

    return
  end subroutine GeographicCoordCnv_geo_to_orth_pos

!OCL SERIAL
  subroutine GeographicCoordCnv_orth_to_geo_vec( geo_v, &
    orth_v, geo_p, Np ) 
    implicit none

    integer, intent(in) :: Np
    real(RP), intent(out) :: geo_v(3,Np)
    real(RP), intent(in) :: orth_v(3,Np)
    real(RP), intent(in) :: geo_p(3,Np)
    integer :: p

    real(RP) :: sin_geo_p1
    real(RP) :: cos_geo_p1
    real(RP) :: sin_geo_p2
    real(RP) :: cos_geo_p2    
    !---------------------------------------------------------  

    !$omp parallel do private(p, cos_geo_p2)
    do p=1, Np
      sin_geo_p1 = sin(geo_p(1,p))
      cos_geo_p1 = cos(geo_p(1,p))
      sin_geo_p2 = sin(geo_p(2,p))
      cos_geo_p2 = cos(geo_p(2,p))

      geo_v(3,p) =   orth_v(1,p) * cos_geo_p2 * cos_geo_p1  &
                   + orth_v(2,p) * cos_geo_p2 * sin_geo_p1 &
                   + orth_v(3,p) * sin_geo_p2
  
      geo_v(2,p) = - orth_v(1,p) * sin_geo_p2 * cos_geo_p1 &
                   - orth_v(2,p) * sin_geo_p2 * sin_geo_p1 &
                   + orth_v(3,p) * cos_geo_p2
    
      geo_v(1,p) = - orth_v(1,p) * sin_geo_p1 &
                   + orth_v(2,p) * cos_geo_p1
    end do

    return
  end subroutine GeographicCoordCnv_orth_to_geo_vec

!OCL SERIAL
  subroutine GeographicCoordCnv_geo_to_orth_vec( orth_v, &
    geo_v, geo_p, Np ) 
    implicit none

    integer, intent(in) :: Np
    real(RP), intent(out) :: orth_v(3,Np)
    real(RP), intent(in) :: geo_v(3,Np)
    real(RP), intent(in) :: geo_p(3,Np)

    integer :: p
    real(RP) :: sin_geo_p1
    real(RP) :: cos_geo_p1
    real(RP) :: sin_geo_p2
    real(RP) :: cos_geo_p2    
    !---------------------------------------------------------  

    !$omp parallel do private(p, sin_geo_p1, cos_geo_p1, sin_geo_p2, cos_geo_p2)
    do p=1, Np
      sin_geo_p1 = sin(geo_p(1,p))
      cos_geo_p1 = cos(geo_p(1,p))
      sin_geo_p2 = sin(geo_p(2,p))
      cos_geo_p2 = cos(geo_p(2,p))

      orth_v(1,p) =   geo_v(3,p) * cos_geo_p1 * cos_geo_p2 &
                    - geo_v(2,p) * cos_geo_p1 * sin_geo_p2 &
                    - geo_v(1,p) * sin_geo_p1

      orth_v(2,p) =   geo_v(3,p) * sin_geo_p1 * cos_geo_p2 &
                    - geo_v(2,p) * sin_geo_p1 * sin_geo_p2 &
                    + geo_v(1,p) * cos_geo_p1

      orth_v(3,p) = geo_v(3,p) * sin_geo_p2 + geo_v(2,p) * cos_geo_p2                  
    end do

    return
  end subroutine GeographicCoordCnv_geo_to_orth_vec

!OCL SERIAL
  subroutine GeographicCoordCnv_rotateX( rotated_pos_vec, &
    pos_vec, angle )
    implicit none
    real(RP), intent(out) :: rotated_pos_vec(3)
    real(RP), intent(in) :: pos_vec(3)
    real(RP), intent(in) :: angle

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
  subroutine GeographicCoordCnv_rotateY( rotated_pos_vec, &
    pos_vec, angle )
    implicit none
    real(RP), intent(out) :: rotated_pos_vec(3)
    real(RP), intent(in) :: pos_vec(3)
    real(RP), intent(in) :: angle

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
  subroutine GeographicCoordCnv_rotateZ( rotated_pos_vec, &
    pos_vec, angle )
    implicit none
    real(RP), intent(out) :: rotated_pos_vec(3)
    real(RP), intent(in) :: pos_vec(3)
    real(RP), intent(in) :: angle

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
