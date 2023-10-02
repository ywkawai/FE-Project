#include "scalelib.h"
module mod_spectral_transform
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
    PRC_abort
  use scale_const, only: &
    PI => CONST_PI, &
    RPlanet => CONST_RADIUS
  
  use scale_mesh_cubedspheredom3d, only: &
    MeshCubedSphereDom3D
  use scale_element_base, only: &
    ElementBase2D
  implicit none
  private

  public :: spectral_tranform

contains
!OCL SERIAL
  subroutine spectral_tranform( g_var, lon, lat, Gsqrt, J,     &
    mesh3D_list, varNum, levelNum, elem2D, Ne2D, Mt, mesh_num, &
    s_var )
    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD  
    implicit none
    integer, intent(in) :: Ne2D
    integer, intent(in) :: mesh_num
    integer, intent(in) :: Mt
    integer, intent(in) :: varNum
    integer, intent(in) :: levelNum
    class(MeshCubedSphereDom3D), intent(in) :: mesh3D_list(mesh_num)
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(in) :: g_var(varNum,levelNum,elem2D%Np,Ne2D,mesh_num)
    real(RP), intent(in) :: lon(elem2D%Np,Ne2D,mesh_num)
    real(RP), intent(in) :: lat(elem2D%Np,Ne2D,mesh_num)
    real(RP), intent(in) :: Gsqrt(elem2D%Np,Ne2D,mesh_num)
    real(RP), intent(in) :: J(elem2D%Np,Ne2D,mesh_num)
    real(RP), intent(out) :: s_var(varNum,levelNum,0:Mt,0:Mt,2)

    integer :: mesh_id
    integer :: m, l

    real(RP) :: s_local(varNum,levelNum,0:Mt,0:Mt,2)

    real(RP) :: fact0
    integer :: ierr
    !----------------------------------------------------------

    !$omp parallel do collapse(2)
    do m=0, Mt
      do l=0, Mt
        s_local(:,:,l,m,1) = 0.0_RP
        s_local(:,:,l,m,2) = 0.0_RP
      end do
    end do

    do mesh_id=1, mesh_num
      call spectral_inv_tranform_local( &
        g_var(:,:,:,:,mesh_id), lon(:,:,mesh_id), lat(:,:,mesh_id),   &
        Gsqrt(:,:,mesh_id), J(:,:,mesh_id), elem2D%IntWeight_lgl,     &
        varNum, levelNum, elem2D%Np, Ne2D, Mt, &
        s_local(:,:,:,:,1), s_local(:,:,:,:,2) )
    end do

    ! global sum
    call MPI_AllReduce( s_local, s_var, varNum * levelNum * (Mt+1)**2 * 2, &
      MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr  )
        
    !-- Nomarization

    !$omp parallel do private(m,l,fact0)
    do m=0, Mt
    do l=0, Mt
      if ( l < m ) cycle

      fact0 = sqrt( 1.0_RP / (2.0_RP * PI) ) / RPlanet**2
      s_var(:,:,l,m,1) = fact0 * s_var(:,:,l,m,1)
      s_var(:,:,l,m,2) = fact0 * s_var(:,:,l,m,2)
    end do
    end do

    return
  end subroutine spectral_tranform

!OCL SERIAL
  subroutine spectral_inv_tranform_local( g_var, lon, lat, Gsqrt, J, intw, &
    varNum, levelNum, Np, Ne2D, Mt, &
    s_r_local, s_i_local )

    integer, intent(in) :: varNum
    integer, intent(in) :: levelNum
    integer, intent(in) :: Np
    integer, intent(in) :: Ne2D
    integer, intent(in) :: Mt
    real(RP), intent(in) :: g_var(varNum*levelNum,Np,Ne2D)
    real(RP), intent(in) :: lon(Np,Ne2D)
    real(RP), intent(in) :: lat(Np,Ne2D)
    real(RP), intent(in) :: Gsqrt(Np,Ne2D)
    real(RP), intent(in) :: J(Np,Ne2D)
    real(RP), intent(in) :: intw(Np)
    real(RP), intent(inout) :: s_r_local(varNum*levelNum,0:Mt,0:Mt)
    real(RP), intent(inout) :: s_i_local(varNum*levelNum,0:Mt,0:Mt)

    integer :: m, l

    integer :: ke2D
    integer :: p

    real(RP) :: mu(Np), cos_lat(Np)
    real(RP) :: Pm_l(-1:Mt,0:Mt)
    real(RP) :: sign_
    real(RP) :: cos_m_lon
    real(RP) :: sin_m_lon
    real(RP) :: coef0

    real(RP) :: fact0(0:Mt)
    integer :: i
    !----------------------------------------------------------

    fact0(:) = 1.0_RP
    fact0(0) = 1.0_RP / sqrt(2.0_RP)

    !$omp parallel do private(i,sign_)
    do l=1, Mt
      ! (2l-1)!!
      ! do i=2*l-1, 2, -2
      !   fact0(l) = fact0(l) * dble(i)
      ! end do
      do i=1, l
        fact0(l) = fact0(l) * dble(2*i+1)/dble(2*i)
      end do
      sign_ = (-1)**(mod(l,2)) 
      fact0(l) = sign_ * sqrt(0.5_RP * fact0(l))
    end do

    do ke2D=1, Ne2D
      mu(:) = sin(lat(:,ke2D))
      cos_lat(:) = cos(lat(:,ke2D))

      do p=1, Np
        !-----
        coef0 = Gsqrt(p,ke2D) * J(p,ke2D) * intw(p)

        !$omp parallel private(l, m, cos_m_lon, sin_m_lon, sign_)
        !$omp do
        do l=0, Mt
          Pm_l(:,l) = 0.0_RP
          ! P^m_m = [ 1/2 * (2m+1)!!/(2m)!! ]^1/2 * (1-mu^2)^m/2
          !
          Pm_l(l,l) = fact0(l) * cos_lat(p)**l
        end do
        !$omp end do

        !$omp do
        do m=0, Mt
          cos_m_lon = cos( mod(dble(m) * lon(p,ke2D), 2.0_RP * PI) )
          sin_m_lon = sin( mod(dble(m) * lon(p,ke2D), 2.0_RP * PI) )

          do l=0, Mt
            if ( m > l ) cycle
            if ( l /= m ) then
              ! 
              Pm_l(l,m) =  ( &
                mu(p) * Pm_l(l-1,m) &
              - sqrt( dble((l-1)**2-m**2)/dble(4*(l-1)**2-1) ) * Pm_l(l-2,m)  &
              ) / sqrt( dble(l**2-m**2)/dble(4*l**2-1) )
            end if

            s_r_local(:,l,m) = s_r_local(:,l,m) &
              + coef0 * Pm_l(l,m) * cos_m_lon * g_var(:,p,ke2D)
            s_i_local(:,l,m) = s_i_local(:,l,m) &
              + coef0 * Pm_l(l,m) * sin_m_lon * g_var(:,p,ke2D)
          end do
        end do
        !$omp end do
        !$omp end parallel
      end do
    end do

    return    
  end subroutine spectral_inv_tranform_local

end module mod_spectral_transform
