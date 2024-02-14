#include "scaleFElib.h"
module scale_matrix_vec_kernel
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_const, only:    &
    UNDEF8 => CONST_UNDEF8, &
    CONST_EPS
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  

  public :: matrix_vec_kernel_dx_8
  public :: matrix_vec_kernel_dy_8
  public :: matrix_vec_kernel_dz_8
  public :: matrix_vec_kernel_DxDyDz_8  
  public :: matrix_vec_kernel_Lift_8
  public :: matrix_vec_kernel_VPM1_8

contains
!OCL SERIAL
  subroutine matrix_vec_kernel_dx_8( A, B, C )
    implicit none
    real(RP), intent(in) :: A(8,8)
    real(RP), intent(in) :: B(8,8**2)
    real(RP), intent(out) :: C(8,8**2)

    integer :: i, j
    !-------------------------------

    do j=1, 8**2
    do i=1, 8
      C(i,j) = a(i,1) * b(1,j) + a(i,2) * b(2,j) + a(i,3) * b(3,j) + a(i,4) * b(4,j) &
             + a(i,5) * b(5,j) + a(i,6) * b(6,j) + a(i,7) * b(7,j) + a(i,8) * b(8,j)
    end do
    end do

    return
  end subroutine matrix_vec_kernel_dx_8

  !OCL SERIAL
  subroutine matrix_vec_kernel_dy_8( B, A, C )
    implicit none
    real(RP), intent(in) :: A(8,8,8)
    real(RP), intent(in) :: B(8,8)
    real(RP), intent(out) :: C(8,8,8)

    integer :: i, j, k
    integer :: kj
    !-------------------------------

    do k=1, 8
    do j=1, 8
    do i=1, 8
      C(i,j,k) = a(i,1,k) * b(1,j) + a(i,2,k) * b(2,j) + a(i,3,k) * b(3,j) + a(i,4,k) * b(4,j) &
               + a(i,5,k) * b(5,j) + a(i,6,k) * b(6,j) + a(i,7,k) * b(7,j) + a(i,8,k) * b(8,j)
    end do
    end do
    end do

    return
  end subroutine matrix_vec_kernel_dy_8
  
  !OCL SERIAL
  subroutine matrix_vec_kernel_dz_8( B, A, C )
    implicit none
    real(RP), intent(in) :: A(8**2,8)
    real(RP), intent(in) :: B(8,8)
    real(RP), intent(out) :: C(8**2,8)

    integer :: i, j
    !-------------------------------

    do j=1, 8
    do i=1, 8**2
      C(i,j) = a(i,1) * b(1,j) + a(i,2) * b(2,j) + a(i,3) * b(3,j) + a(i,4) * b(4,j) &
             + a(i,5) * b(5,j) + a(i,6) * b(6,j) + a(i,7) * b(7,j) + a(i,8) * b(8,j)
    end do
    end do

    return
  end subroutine matrix_vec_kernel_dz_8

!OCL SERIAL
  subroutine matrix_vec_kernel_DxDyDz_8( D, Dtr, q1, q2, q3, dx,  dy,  dz )
    implicit none
    real(RP), intent(in) :: D(8,8)
    real(RP), intent(in) :: Dtr(8,8)
!    real(RP), intent(in) :: q1(8,8,8)
    real(RP), intent(in) :: q1(8,8**2)
    real(RP), intent(in) :: q2(8,8,8)
    real(RP), intent(in) :: q3(8,8,8)
!    real(RP), intent(out) :: grad3D_x(8,8,8)
    real(RP), intent(out) :: dx(8,8**2)
    real(RP), intent(out) :: dy(8,8,8)
    real(RP), intent(out) :: dz(8,8,8)

    integer :: i, j, k
    integer :: ii, jj, ij
    !-------------------------------

    ! do k=1, 8
    !   do j=1, 8
    !   do i=1, 8
    !     jj = j + (k-1)*8
    !     grad3D_x(i,j,k) = D(i,1) * q1(1,jj) + D(i,2) * q1(2,jj) + D(i,3) * q1(3,jj) + D(i,4) * q1(4,jj) &
    !                     + D(i,5) * q1(5,jj) + D(i,6) * q1(6,jj) + D(i,7) * q1(7,jj) + D(i,8) * q1(8,jj)
    !   end do
    !   end do
    ! end do
    do jj=1, 8**2
      do i=1, 8
        dx(i,jj) = D(i,1) * q1(1,jj) + D(i,2) * q1(2,jj) + D(i,3) * q1(3,jj) + D(i,4) * q1(4,jj) &
                 + D(i,5) * q1(5,jj) + D(i,6) * q1(6,jj) + D(i,7) * q1(7,jj) + D(i,8) * q1(8,jj)
      end do
    end do
    ! do k=1, 8
    ! do j=1, 8
    !   do i=1, 8
    !     grad3D_x(i,j,k) = D(i,1) * q1(1,j,k) + D(i,2) * q1(2,j,k) + D(i,3) * q1(3,j,k) + D(i,4) * q1(4,j,k) &
    !                     + D(i,5) * q1(5,j,k) + D(i,6) * q1(6,j,k) + D(i,7) * q1(7,j,k) + D(i,8) * q1(8,j,k)
    !   end do
    ! end do
    ! end do

    do k=1, 8
      do j=1, 8
      do i=1, 8
        dy(i,j,k) = q2(i,1,k) * Dtr(1,j) + q2(i,2,k) * Dtr(2,j) + q2(i,3,k) * Dtr(3,j) + q2(i,4,k) * Dtr(4,j) &
                  + q2(i,5,k) * Dtr(5,j) + q2(i,6,k) * Dtr(6,j) + q2(i,7,k) * Dtr(7,j) + q2(i,8,k) * Dtr(8,j)
      end do
      end do
    end do

    do k=1, 8
    do j=1, 8
    do i=1, 8      
      dz(i,j,k) = q3(i,j,1) * Dtr(1,k) + q3(i,j,2) * Dtr(2,k) + q3(i,j,3) * Dtr(3,k) + q3(i,j,4) * Dtr(4,k) &
                + q3(i,j,5) * Dtr(5,k) + q3(i,j,6) * Dtr(6,k) + q3(i,j,7) * Dtr(7,k) + q3(i,j,8) * Dtr(8,k)
    end do
    end do
    end do

    return
  end subroutine matrix_vec_kernel_DxDyDz_8

!OCL SERIAL
  subroutine matrix_vec_kernel_div3D_8( D, Dtr, q1, q2, &
    E11, E22, E33, Gsqrt, grad3D )
    implicit none
    real(RP), intent(in) :: D(8,8)
    real(RP), intent(in) :: Dtr(8,8)
    real(RP), intent(in) :: q1(8,8**2)
    real(RP), intent(in) :: q2(8,8,8)
    real(RP), intent(in) :: E11(8,8**2)
    real(RP), intent(in) :: E22(8,8**2)
    real(RP), intent(in) :: E33(8,8**2)
    real(RP), intent(in) :: Gsqrt(8,8**2)
    real(RP), intent(out) :: grad3D(8,8**2)

    ! real(RP), intent(in) :: q1(8,8,8)
    ! real(RP), intent(in) :: q2(8,8,8)
    ! real(RP), intent(in) :: E11(8,8,8)
    ! real(RP), intent(in) :: E22(8,8,8)
    ! real(RP), intent(in) :: E33(8,8,8)
    ! real(RP), intent(in) :: Gsqrt(8,8,8)
    ! real(RP), intent(out) :: grad3D(8,8,8)

    integer :: i, j, k
    integer :: ii, jj
    !-------------------------------

    do jj=1, 8**2
    do i=1, 8
      grad3D(i,jj) = &
            E11(i,jj) * ( &
              D(i,1) * q1(1,jj) + D(i,2) * q1(2,jj) + D(i,3) * q1(3,jj) + D(i,4) * q1(4,jj) &
            + D(i,5) * q1(5,jj) + D(i,6) * q1(6,jj) + D(i,7) * q1(7,jj) + D(i,8) * q1(8,jj) )
    end do
    end do
    do k=1, 8
    do j=1, 8
      jj = j + (k-1)*8
      do i=1, 8
        grad3D(i,jj) = ( grad3D(i,jj) &        
            + E22(i,jj) * ( &
                q2(i,1,k) * Dtr(1,j) + q2(i,2,k) * Dtr(2,j) + q2(i,3,k) * Dtr(3,j) + q2(i,4,k) * Dtr(4,j)   &
              + q2(i,5,k) * Dtr(5,j) + q2(i,6,k) * Dtr(6,j) + q2(i,7,k) * Dtr(7,j) + q2(i,8,k) * Dtr(8,j) ) &
            + E33(i,jj) * ( & 
                q2(i,j,1) * Dtr(1,k) + q2(i,j,2) * Dtr(2,k) + q2(i,j,3) * Dtr(3,k) + q2(i,j,4) * Dtr(4,k)   &
              + q2(i,j,5) * Dtr(5,k) + q2(i,j,6) * Dtr(6,k) + q2(i,j,7) * Dtr(7,k) + q2(i,j,8) * Dtr(8,k) ) &
          ) / Gsqrt(i,jj)
      end do
    end do
    end do

    return
  end subroutine matrix_vec_kernel_div3D_8

!OCL SERIAL
  subroutine matrix_vec_kernel_div3D_8_2( D, Dtr, q1, q2, &
    Escale, Gsqrt, Fx, Fy, Fz, grad3D )
    implicit none
    real(RP), intent(in) :: D(8,8)
    real(RP), intent(in) :: Dtr(8,8)
    real(RP), intent(in) :: q1(8,8**2)
    real(RP), intent(in) :: q2(8,8,8)
    real(RP), intent(in) :: Escale(3,8,8**2)
    real(RP), intent(in) :: Gsqrt(8,8**2)
    real(RP), intent(out) :: Fx(8,8**2)
    real(RP), intent(out) :: Fy(8,8**2)
    real(RP), intent(out) :: Fz(8,8**2)
    real(RP), intent(out) :: grad3D(8,8**2)

    ! real(RP), intent(in) :: q1(8,8,8)
    ! real(RP), intent(in) :: q2(8,8,8)
    ! real(RP), intent(in) :: E11(8,8,8)
    ! real(RP), intent(in) :: E22(8,8,8)
    ! real(RP), intent(in) :: E33(8,8,8)
    ! real(RP), intent(in) :: Gsqrt(8,8,8)
    ! real(RP), intent(out) :: grad3D(8,8,8)

    integer :: i, j, k
    integer :: ii, jj
    !-------------------------------

    do jj=1, 8**2
    do i=1, 8
      Fx(i,jj) =  &
              D(i,1) * q1(1,jj) + D(i,2) * q1(2,jj) + D(i,3) * q1(3,jj) + D(i,4) * q1(4,jj) &
            + D(i,5) * q1(5,jj) + D(i,6) * q1(6,jj) + D(i,7) * q1(7,jj) + D(i,8) * q1(8,jj)
    end do
    end do
    do k=1, 8
    do j=1, 8
      jj = j + (k-1)*8
      do i=1, 8
        Fy(i,jj) = &
                q2(i,1,k) * Dtr(1,j) + q2(i,2,k) * Dtr(2,j) + q2(i,3,k) * Dtr(3,j) + q2(i,4,k) * Dtr(4,j)   &
              + q2(i,5,k) * Dtr(5,j) + q2(i,6,k) * Dtr(6,j) + q2(i,7,k) * Dtr(7,j) + q2(i,8,k) * Dtr(8,j)
      end do 
      do i=1, 8
        Fz(i,jj) = &
                q2(i,j,1) * Dtr(1,k) + q2(i,j,2) * Dtr(2,k) + q2(i,j,3) * Dtr(3,k) + q2(i,j,4) * Dtr(4,k)   &
              + q2(i,j,5) * Dtr(5,k) + q2(i,j,6) * Dtr(6,k) + q2(i,j,7) * Dtr(7,k) + q2(i,j,8) * Dtr(8,k) 
      end do
    end do
    end do

    do jj=1, 8**2
    do i=1, 8
      grad3D(i,jj) = - ( Escale(1,i,jj) * Fx(i,jj) + Escale(2,i,jj) * Fy(i,jj) + Escale(3,i,jj) * Fz(i,jj) ) / Gsqrt(i,jj)
    end do
    end do

    return
  end subroutine matrix_vec_kernel_div3D_8_2

!OCL SERIAL
  subroutine matrix_vec_kernel_Lift_8( &
    Lift, bnd_flux, LiftDFlx )
    implicit none
    real(RP), intent(in) :: Lift(8,8,8,6)
    real(RP), intent(in) :: bnd_flux(8,8,6)
    real(RP), intent(out) :: LiftDFlx(8,8,8)

    integer :: i, j, k
    !-------------------------------------------

    do k=1, 8
    do j=1, 8
    do i=1, 8
      LiftDFlx(i,j,k) = &
          Lift(i,j,k,1) * bnd_flux(i,k,1) &
        + Lift(i,j,k,2) * bnd_flux(j,k,2) &
        + Lift(i,j,k,3) * bnd_flux(i,k,3) &
        + Lift(i,j,k,4) * bnd_flux(j,k,4) &
        + Lift(i,j,k,5) * bnd_flux(i,j,5) &
        + Lift(i,j,k,6) * bnd_flux(i,j,6)
    end do
    end do
    end do  

    return
  end subroutine matrix_vec_kernel_Lift_8

!OCL SERIAL
  subroutine matrix_vec_kernel_VPM1_8( &
      VPM1_tr, q0, q )
    implicit none
    real(RP), intent(in) :: VPM1_tr(8,8)
    real(RP), intent(in) :: q0(8**2,8)
    real(RP), intent(out) :: q(8**2,8)
    integer :: ij
    integer :: k
    !--------------------------------

    do k=1, 8
    do ij=1, 8**2
      q(ij,k) = &
              q0(ij,1) * VPM1_tr(1,k) + q0(ij,2) * VPM1_tr(2,k) + q0(ij,3) * VPM1_tr(3,k) + q0(ij,4) * VPM1_tr(4,k)   &
            + q0(ij,5) * VPM1_tr(5,k) + q0(ij,6) * VPM1_tr(6,k) + q0(ij,7) * VPM1_tr(7,k) + q0(ij,8) * VPM1_tr(8,k) 
    end do
  end do

    return
  end subroutine matrix_vec_kernel_VPM1_8

end module scale_matrix_vec_kernel
