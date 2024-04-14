!-------------------------------------------------------------------------------
!> module Utility for mktopo
!!
!! @par Description
!!          subroutines useful to prepare topography data
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_mktopo_util
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_const, only: &
    PI => CONST_PI
    
  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement    
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_meshfield_base, only: MeshField2D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: mktopoutil_barocwave_global_JW2006_calc_topo
  public :: mktopoutil_gen_GPMat
  public :: mktopoutil_GalerkinProjection

contains

  !-----------------------------------------------------------------------------
  !> Calculate a topography for a global test case of idealized baroclinic wave in Jablonowski and Williamson (2006)
  !!
  !! @par Reference
  !! - Jablonowski, C., & Williamson, D. L., 2006: A baroclinic instability test case for atmospheric model dynamical cores. Quarterly Journal of the Royal Meteorological Society: A journal of the atmospheric sciences, applied meteorology and physical oceanography, 132(621C), 2943-2975.
  !!
!OCL SERIAL
  subroutine mktopoutil_barocwave_global_JW2006_calc_topo( topo, &
    U0, ETA0, lat, Np )

    use scale_const, only: &
      OHM => CONST_OHM,      &
      Grav => CONST_GRAV,    &
      RPlanet => CONST_RADIUS

    implicit none

    integer, intent(in) :: Np
    real(RP), intent(out) :: topo(Np)
    real(RP), intent(in) :: U0        !< The value of η at a reference level (position of the jet)
    real(RP), intent(in) :: ETA0      !< The parameter associated with zonal jet maximum amplitude  [m/s]
    real(RP), intent(in) :: lat(Np)   !< latitude [radian]

    real(RP) :: sin_lat(Np)
    real(RP) :: sin_lat_pow_6(Np)
    real(RP) :: cos_lat(Np)

    real(RP ) :: tmp
    !-------------------------------------------

    sin_lat(:) = sin(lat(:))
    sin_lat_pow_6(:) = sin_lat(:)**6
    cos_lat(:) = cos(lat(:))

    tmp = cos( (1.0_RP - ETA0) * 0.5_RP * PI )
    tmp = U0 * tmp * sqrt(tmp)

    ! Calc horizontal variation of geopotential height
    topo(:) = tmp * & 
      (   tmp * ( - 2.0_RP * sin_lat_pow_6(:) * ( cos_lat(:)**2 + 1.0_RP / 3.0_RP ) + 10.0_RP / 63.0_RP )          &
        + RPlanet * OHM * ( 8.0_RP / 5.0_RP * cos_lat(:)**3 * ( sin_lat(:)**2 + 2.0_RP / 3.0_RP ) - 0.25_RP * PI ) &
      ) / Grav

    return
  end subroutine mktopoutil_barocwave_global_JW2006_calc_topo

!OCL SERIAL
  subroutine mktopoutil_gen_GPMat( GPMat, &
    elem_intrp, elem )
    implicit none

    class(ElementBase2D), intent(in) :: elem_intrp
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: GPMat(elem%Np,elem_intrp%Np)

    integer :: p1, p2, p_
    integer :: p_intrp

    real(RP) :: InvV_intrp(elem%Np,elem_intrp%Np)
    !---------------------------------------------

    InvV_intrp(:,:) = 0.0_RP
    do p2=1, elem%PolyOrder+1
    do p1=1, elem%PolyOrder+1
      p_ = p1 + (p2-1)*(elem%PolyOrder + 1) 
      p_intrp = p1 + (p2-1)*(elem_intrp%PolyOrder + 1)
      InvV_intrp(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do
    GPMat(:,:) = matmul(elem%V, InvV_intrp)

    return
  end subroutine mktopoutil_gen_GPMat

!OCL SERIAL  
  subroutine mktopoutil_GalerkinProjection( q, &
    func, IntrpPolyOrder_h,                    &
    lcmesh2D, elem                             )
  
  implicit none
  class(LocalMesh2D), intent(in) :: lcmesh2D
  class(ElementBase2D), intent(in) :: elem
  real(RP), intent(out) :: q(elem%Np,lcmesh2D%NeA)
  integer, intent(in) :: IntrpPolyOrder_h

  interface
    subroutine func( q_intrp, &
        x, y, elem_intrp   )
      import ElementBase2D
      import RP
      class(ElementBase2D), intent(in) :: elem_intrp
      real(RP), intent(out) :: q_intrp(elem_intrp%Np)
      real(RP), intent(in) :: x(elem_intrp%Np)
      real(RP), intent(in) :: y(elem_intrp%Np)
    end subroutine func
  end interface

  type(QuadrilateralElement) :: elem_intrp
  real(RP), allocatable :: x_intrp(:,:), y_intrp(:,:)
  real(RP) :: vx(elem%Nv), vy(elem%Nv)

  real(RP), allocatable :: IntrpMat(:,:)
  real(RP), allocatable :: q_intrp(:)

  integer :: ke
  !-----------------------------------------------

  call elem_intrp%Init( IntrpPolyOrder_h, .false. )

  allocate( IntrpMat(elem%Np,elem_intrp%Np) )
  call mktopoutil_gen_GPMat( IntrpMat, elem_intrp, elem )

  allocate( x_intrp(elem_intrp%Np,lcmesh2D%Ne), y_intrp(elem_intrp%Np,lcmesh2D%Ne) )
  allocate( q_intrp(elem_intrp%Np) )

  !$omp parallel do private(vx, vy)
  do ke=lcmesh2D%NeS, lcmesh2D%NeE
    vx(:) = lcmesh2D%pos_ev(lcmesh2D%EToV(ke,:),1)
    vy(:) = lcmesh2D%pos_ev(lcmesh2D%EToV(ke,:),2)
    x_intrp(:,ke) = vx(1) + 0.5_RP * ( elem_intrp%x1(:) + 1.0_RP ) * ( vx(2) - vx(1) ) 
    y_intrp(:,ke) = vy(1) + 0.5_RP * ( elem_intrp%x2(:) + 1.0_RP ) * ( vy(4) - vy(1) )
  end do

  !$omp parallel do private( q_intrp )
  do ke=lcmesh2D%NeS, lcmesh2D%NeE

    call func( q_intrp,                 & ! (out)
      x_intrp(:,ke), y_intrp(:,ke),     & ! (in)
      elem_intrp                        ) ! (in)
    
    ! Perform Galerkin projection
    q(:,ke) = matmul( IntrpMat, q_intrp )
  end do

  call elem_intrp%Final()

  return
end subroutine mktopoutil_GalerkinProjection

end module mod_mktopo_util