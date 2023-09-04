!-------------------------------------------------------------------------------
!> module Utility for mkinit
!!
!! @par Description
!!          subroutines useful to prepare initial data
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_mkinit_util
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc

  use scale_const, only: &
    PI => CONST_PI

  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D  
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  
    
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: mkinitutil_gen_GPMat 
  public :: mkinitutil_gen_Vm1Mat 
  public :: mkinitutil_calc_cosinebell
  public :: mkinitutil_calc_cosinebell_global
  public :: mkinitutil_GalerkinProjection
  public :: mkinitutil_GalerkinProjection_global

contains

!OCL SERIAL
  subroutine mkinitutil_gen_GPMat( GPMat, &
    elem_intrp, elem )
    implicit none

    class(ElementBase3D), intent(in) :: elem_intrp
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: GPMat(elem%Np,elem_intrp%Np)

    integer :: p1, p2, p3, p_
    integer :: p_intrp

    real(RP) :: InvV_intrp(elem%Np,elem_intrp%Np)
    !---------------------------------------------

    InvV_intrp(:,:) = 0.0_RP
    do p3=1, elem%PolyOrder_v+1
    do p2=1, elem%PolyOrder_h+1
    do p1=1, elem%PolyOrder_h+1
      p_ = p1 + (p2-1)*(elem%PolyOrder_h + 1) + (p3-1)*(elem%PolyOrder_h + 1)**2
      p_intrp = p1 + (p2-1)*(elem_intrp%PolyOrder_h + 1) + (p3-1)*(elem_intrp%PolyOrder_h + 1)**2
      InvV_intrp(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do
    end do
    GPMat(:,:) = matmul(elem%V, InvV_intrp)

    return
  end subroutine mkinitutil_gen_GPMat

!OCL SERIAL  
  subroutine mkinitutil_gen_Vm1Mat( Vm1Mat, &
    elem_intrp, elem )
    implicit none

    class(ElementBase3D), intent(in) :: elem_intrp
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: Vm1Mat(elem%Np,elem_intrp%Np)

    integer :: p1, p2, p3, p_
    integer :: p_intrp

    real(RP) :: InvV_intrpVm1(elem%Np,elem_intrp%Np)
    !---------------------------------------------

    InvV_intrpVm1(:,:) = 0.0_RP
    do p3=1, elem%PolyOrder_v
    do p2=1, elem%PolyOrder_h+1
    do p1=1, elem%PolyOrder_h+1
      p_ = p1 + (p2-1)*(elem%PolyOrder_h + 1) + (p3-1)*(elem%PolyOrder_h + 1)**2
      p_intrp = p1 + (p2-1)*(elem_intrp%PolyOrder_h + 1) + (p3-1)*(elem_intrp%PolyOrder_h + 1)**2
      InvV_intrpVm1(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do
    end do
    Vm1Mat(:,:) = matmul(elem%V, InvV_intrpVm1)

    return
  end subroutine mkinitutil_gen_Vm1Mat

  !>
  !! Calculate the distribution function of a cosine bell in regional domain
  !! 
  !! If the vertical dependence is considered, specify z_func_type and z_func_params.
  !!  For z_func_type = 'sin', the values of z_func_params is
  !!   1:  the vertical model, 2: the half of wavelength
  !! 
!OCL SERIAL  
  subroutine mkinitutil_calc_cosinebell( &
    q,                                           &
    qmax, rx, ry, rz, xc, yc, zc,                &
    x, y, z, lcmesh3D, elem,                     &
    IntrpPolyOrder_h, IntrpPolyOrder_v,          &
    z_func_type, z_func_params, cosbell_exponent )

    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: q(elem%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: qmax
    real(RP), intent(in) :: rx, ry, rz
    real(RP), intent(in) :: xc, yc, zc
    real(RP), intent(in) :: x(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh3D%Ne)
    integer, intent(in) :: IntrpPolyOrder_h
    integer, intent(in) :: IntrpPolyOrder_v
    character(len=*), optional, intent(in) :: z_func_type
    real(RP), optional, intent(in) :: z_func_params(:)
    integer, intent(in), optional :: cosbell_exponent  !< parameter to ensure 2*cosbell_exponent-1 continuous derivatives

    integer :: ke

    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: x_intrp(:,:), y_intrp(:,:), z_intrp(:,:)
    real(RP), allocatable :: r_intrp(:)
    real(RP), allocatable :: z_func(:,:)
    real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)

    real(RP), allocatable :: IntrpMat(:,:)
    real(RP), allocatable :: q_intrp(:)

    integer :: exponent    
    !-----------------------------------------------

    if ( present(cosbell_exponent) ) then
      exponent = cosbell_exponent
    else
      exponent = 1
    end if

    call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, .false. )

    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    call mkinitutil_gen_GPMat( IntrpMat, elem_intrp, elem )

    allocate( x_intrp(elem_intrp%Np,lcmesh3D%Ne), y_intrp(elem_intrp%Np,lcmesh3D%Ne), z_intrp(elem_intrp%Np,lcmesh3D%Ne) )
    allocate( z_func(elem_intrp%Np,lcmesh3D%Ne))
    allocate( r_intrp(elem_intrp%Np) )
    allocate( q_intrp(elem_intrp%Np) )

    !$omp parallel private( &
    !$omp q_intrp, vx, vy, vz,               &
    !$omp x_intrp, y_intrp, z_intrp, r_intrp )

    !$omp do
    do ke=lcmesh3D%NeS, lcmesh3D%NeE

      vx(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),1)
      vy(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),2)
      vz(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),3)
      x_intrp(:,ke) = vx(1) + 0.5_RP * ( elem_intrp%x1(:) + 1.0_RP ) * ( vx(2) - vx(1) ) 
      y_intrp(:,ke) = vy(1) + 0.5_RP * ( elem_intrp%x2(:) + 1.0_RP ) * ( vy(4) - vy(1) )
      z_intrp(:,ke) = vz(1) + 0.5_RP * ( elem_intrp%x3(:) + 1.0_RP ) * ( vz(5) - vz(1) )

      z_func(:,ke) = 1.0_RP
    end do
    !$omp end do

    ! Calculate the vertical function
    if ( present(z_func_type) ) then
      select case(z_func_type)
      case ('sin')
        !$omp do
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          z_func(:,ke) = sin( z_func_params(1) * PI * z_intrp(:,ke) / z_func_params(2) )
        end do
        !$omp end do
      end select
    end if

    !$omp do
    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      r_intrp(:) = sqrt( &
          ( (x_intrp(:,ke) - xc) / rx )**2 &
        + ( (y_intrp(:,ke) - yc) / ry )**2 &
        + ( (z_intrp(:,ke) - zc) / rz )**2 )
      
      where( r_intrp(:) <= 1.0_RP ) 
        q_intrp(:) = qmax * ( 0.5_RP * (1.0_RP + cos( PI * r_intrp(:) ) ) )**exponent
      elsewhere
        q_intrp(:) = 0.0_RP
      end where
      q(:,ke) = matmul(IntrpMat, q_intrp(:) * z_func(:,ke))
    end do
    !$omp end do
    !$omp end parallel

    call elem_intrp%Final()

    return
  end subroutine mkinitutil_calc_cosinebell

  !>
  !! Calculate the distribution function of a cosine bell in global domain
  !! 
  !! If the vertical dependence is considered, specify z_func_type and z_func_params.
  !!  For z_func_type = 'sin', the values of z_func_params is
  !!   1:  the vertical model, 2: the half of wavelength
  !! 
!OCL SERIAL
  subroutine mkinitutil_calc_cosinebell_global( &
    q,                                           & 
    qmax, rh, lonc, latc, rplanet,               &
    x, y, z, lcmesh3D, elem,                     &
    IntrpPolyOrder_h, IntrpPolyOrder_v,          &
    z_func_type, z_func_params, cosbell_exponent )

    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatPos
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: q(elem%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: qmax
    real(RP), intent(in) :: rh
    real(RP), intent(in) :: lonc, latc
    real(RP), intent(in) :: rplanet
    real(RP), intent(in) :: x(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh3D%Ne)
    integer, intent(in) :: IntrpPolyOrder_h
    integer, intent(in) :: IntrpPolyOrder_v
    character(len=*), optional, intent(in) :: z_func_type
    real(RP), optional, intent(in) :: z_func_params(:)
    integer, intent(in), optional :: cosbell_exponent      !< parameter to ensure 2*cosbell_exponent-1 continuous derivatives

    integer :: ke

    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: x_intrp(:,:), y_intrp(:,:), z_intrp(:,:)
    real(RP), allocatable :: gam_intrp(:,:)
    real(RP), allocatable :: lon_intrp(:,:), lat_intrp(:,:)
    real(RP), allocatable :: z_func(:,:)
    real(RP), allocatable :: r_intrp(:)
    real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)

    real(RP), allocatable :: IntrpMat(:,:)
    real(RP), allocatable :: q_intrp(:)

    integer :: exponent

    !-----------------------------------------------

    if ( present(cosbell_exponent) ) then
      exponent = cosbell_exponent
    else
      exponent = 1
    end if

    call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, .false. )

    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    call mkinitutil_gen_GPMat( IntrpMat, elem_intrp, elem )

    allocate( x_intrp(elem_intrp%Np,lcmesh3D%Ne), y_intrp(elem_intrp%Np,lcmesh3D%Ne), z_intrp(elem_intrp%Np,lcmesh3D%Ne) )
    allocate( gam_intrp(elem_intrp%Np,lcmesh3D%Ne) )
    allocate( lon_intrp(elem_intrp%Np,lcmesh3D%Ne), lat_intrp(elem_intrp%Np,lcmesh3D%Ne) )
    allocate( z_func(elem_intrp%Np,lcmesh3D%Ne) )
    allocate( r_intrp(elem_intrp%Np) )
    allocate( q_intrp(elem_intrp%Np) )


    !$omp parallel do private(vx, vy, vz)
    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      vx(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),1)
      vy(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),2)
      vz(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),3)
      x_intrp(:,ke) = vx(1) + 0.5_RP * ( elem_intrp%x1(:) + 1.0_RP ) * ( vx(2) - vx(1) ) 
      y_intrp(:,ke) = vy(1) + 0.5_RP * ( elem_intrp%x2(:) + 1.0_RP ) * ( vy(4) - vy(1) )
      z_intrp(:,ke) = vz(1) + 0.5_RP * ( elem_intrp%x3(:) + 1.0_RP ) * ( vz(5) - vz(1) )
      gam_intrp(:,ke) = 1.0_RP

      z_func(:,ke) = 1.0_RP
    end do

    call CubedSphereCoordCnv_CS2LonLatPos( lcmesh3D%panelID, x_intrp, y_intrp, gam_intrp, & ! (in)
      elem_intrp%Np * lcmesh3D%Ne,                                                        & ! (in)
      lon_intrp(:,:), lat_intrp(:,:) )                                                      ! (out)

    ! Calculate the vertical function
    if ( present(z_func_type) ) then
      select case(z_func_type)
      case ('sin')
        !$omp parallel do
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          z_func(:,ke) = sin( z_func_params(1) * PI * z_intrp(:,ke) / z_func_params(2) )
        end do
      end select
    end if

    !$omp parallel do private( r_intrp, q_intrp )
    do ke=lcmesh3D%NeS, lcmesh3D%NeE

      ! Calculate the horizontal function
      r_intrp(:) = rplanet / rh * acos( sin(latc) * sin(lat_intrp(:,ke)) + cos(latc) * cos(lat_intrp(:,ke)) * cos(lon_intrp(:,ke) - lonc) )
      where( r_intrp(:) <= 1.0_RP ) 
        q_intrp(:) = qmax * ( 0.5_RP * (1.0_RP + cos( PI * r_intrp(:) ) ) )**exponent
      elsewhere
        q_intrp(:) = 0.0_RP
      end where

      ! Perform Galerkin projection
      q(:,ke) = matmul(IntrpMat, q_intrp(:) * z_func(:,ke))
    end do

    call elem_intrp%Final()

    return
  end subroutine mkinitutil_calc_cosinebell_global


!>
!! Apply the Galerkin projection to the user-defined function 
!! 
!OCL SERIAL  
  subroutine mkinitutil_GalerkinProjection( q, &
    func, IntrpPolyOrder_h, IntrpPolyOrder_v,  &
    lcmesh3D, elem                             )
  
  implicit none
  class(LocalMesh3D), intent(in) :: lcmesh3D
  class(ElementBase3D), intent(in) :: elem
  real(RP), intent(out) :: q(elem%Np,lcmesh3D%NeA)
  integer, intent(in) :: IntrpPolyOrder_h
  integer, intent(in) :: IntrpPolyOrder_v

  interface
    subroutine func( q_intrp, &
        x, y, z, elem_intrp   )
      import ElementBase3D
      import RP
      class(ElementBase3D), intent(in) :: elem_intrp
      real(RP), intent(out) :: q_intrp(elem_intrp%Np)
      real(RP), intent(in) :: x(elem_intrp%Np)
      real(RP), intent(in) :: y(elem_intrp%Np)
      real(RP), intent(in) :: z(elem_intrp%Np)
    end subroutine func
  end interface

  type(HexahedralElement) :: elem_intrp
  real(RP), allocatable :: x_intrp(:,:), y_intrp(:,:), z_intrp(:,:)
  real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)

  real(RP), allocatable :: IntrpMat(:,:)
  real(RP), allocatable :: q_intrp(:)

  integer :: ke
  !-----------------------------------------------

  call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, .false. )

  allocate( IntrpMat(elem%Np,elem_intrp%Np) )
  call mkinitutil_gen_GPMat( IntrpMat, elem_intrp, elem )

  allocate( x_intrp(elem_intrp%Np,lcmesh3D%Ne), y_intrp(elem_intrp%Np,lcmesh3D%Ne), z_intrp(elem_intrp%Np,lcmesh3D%Ne) )
  allocate( q_intrp(elem_intrp%Np) )

  !$omp parallel do private(vx, vy, vz)
  do ke=lcmesh3D%NeS, lcmesh3D%NeE
    vx(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),1)
    vy(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),2)
    vz(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),3)
    x_intrp(:,ke) = vx(1) + 0.5_RP * ( elem_intrp%x1(:) + 1.0_RP ) * ( vx(2) - vx(1) ) 
    y_intrp(:,ke) = vy(1) + 0.5_RP * ( elem_intrp%x2(:) + 1.0_RP ) * ( vy(4) - vy(1) )
    z_intrp(:,ke) = vz(1) + 0.5_RP * ( elem_intrp%x3(:) + 1.0_RP ) * ( vz(5) - vz(1) )
  end do

  !$omp parallel do private( q_intrp )
  do ke=lcmesh3D%NeS, lcmesh3D%NeE

    call func( q_intrp,                                & ! (out)
      x_intrp(:,ke), y_intrp(:,ke), z_intrp(:,ke),     & ! (in)
      elem_intrp                                       ) ! (in)
    
    ! Perform Galerkin projection
    q(:,ke) = matmul( IntrpMat, q_intrp )
  end do

  call elem_intrp%Final()

  return
end subroutine mkinitutil_GalerkinProjection

!>
!! Apply the Galerkin projection to the user-defined function 
!! (for global model)
!OCL SERIAL  
  subroutine mkinitutil_GalerkinProjection_global( q, &
      func, IntrpPolyOrder_h, IntrpPolyOrder_v,       &
      lcmesh3D, elem, rplanet )

    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatPos
    
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: q(elem%Np,lcmesh3D%NeA)
    integer, intent(in) :: IntrpPolyOrder_h
    integer, intent(in) :: IntrpPolyOrder_v
    real(RP), intent(in) :: rplanet

    interface
      subroutine func( q_intrp,             &
          lon, lat, z, elem_intrp, rplanet_ )
        import ElementBase3D
        import RP
        class(ElementBase3D), intent(in) :: elem_intrp
        real(RP), intent(out) :: q_intrp(elem_intrp%Np)
        real(RP), intent(in) :: lon(elem_intrp%Np)
        real(RP), intent(in) :: lat(elem_intrp%Np)
        real(RP), intent(in) :: z(elem_intrp%Np)
        real(RP), intent(in) :: rplanet_
      end subroutine func
    end interface

    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: x_intrp(:,:), y_intrp(:,:), z_intrp(:,:)
    real(RP), allocatable :: gam_intrp(:,:)
    real(RP), allocatable :: lon_intrp(:,:), lat_intrp(:,:)
    real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)

    real(RP), allocatable :: IntrpMat(:,:)
    real(RP), allocatable :: q_intrp(:)

    integer :: ke
    !-----------------------------------------------

    call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, .false. )

    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    call mkinitutil_gen_GPMat( IntrpMat, elem_intrp, elem )

    allocate( x_intrp(elem_intrp%Np,lcmesh3D%Ne), y_intrp(elem_intrp%Np,lcmesh3D%Ne), z_intrp(elem_intrp%Np,lcmesh3D%Ne) )
    allocate( gam_intrp(elem_intrp%Np,lcmesh3D%Ne) )
    allocate( lon_intrp(elem_intrp%Np,lcmesh3D%Ne), lat_intrp(elem_intrp%Np,lcmesh3D%Ne) )
    allocate( q_intrp(elem_intrp%Np) )

    !$omp parallel do private(vx, vy, vz)
    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      vx(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),1)
      vy(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),2)
      vz(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),3)
      x_intrp(:,ke) = vx(1) + 0.5_RP * ( elem_intrp%x1(:) + 1.0_RP ) * ( vx(2) - vx(1) ) 
      y_intrp(:,ke) = vy(1) + 0.5_RP * ( elem_intrp%x2(:) + 1.0_RP ) * ( vy(4) - vy(1) )
      z_intrp(:,ke) = vz(1) + 0.5_RP * ( elem_intrp%x3(:) + 1.0_RP ) * ( vz(5) - vz(1) )

      gam_intrp(:,ke) = 1.0_RP
    end do

    call CubedSphereCoordCnv_CS2LonLatPos( lcmesh3D%panelID, x_intrp, y_intrp, gam_intrp, & ! (in)
      elem_intrp%Np * lcmesh3D%Ne,                                                        & ! (in)
      lon_intrp(:,:), lat_intrp(:,:) )                                                      ! (out)

    !$omp parallel do private( q_intrp )
    do ke=lcmesh3D%NeS, lcmesh3D%NeE

      call func( q_intrp,                                & ! (out)
        lon_intrp(:,ke), lat_intrp(:,ke), z_intrp(:,ke), & ! (in)
        elem_intrp, rplanet                              ) ! (in)
      
      ! Perform Galerkin projection
      q(:,ke) = matmul( IntrpMat, q_intrp )
    end do

    call elem_intrp%Final()
  
    return
  end subroutine mkinitutil_GalerkinProjection_global

end module mod_mkinit_util
