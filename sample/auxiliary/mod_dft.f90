!-------------------------------------------------------------------------------
!> module sample / auxiliary
!!
!! @par Description
!!      A module to provide discrete Fourier transform
!!
!! @par Reference
!!
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_dft
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_const, only: &
    PI => CONST_PI

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  public :: DFT_dft1d
  public :: DFT_idft1d
  public :: DFT_interp1d

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------
contains  
!OCL SERIAL
  subroutine DFT_dft1d(q_out, q_in, nmax, mmax)
    implicit none
    integer, intent(in) :: nmax, mmax
    complex(RP), intent(out) :: q_out(-mmax/2:mmax/2)
    real(RP), intent(in) :: q_in(0:nmax)

    complex(RP) :: phi
    complex(RP), parameter :: i = (0.0, 1.0)
    integer :: j, k
    !--------------------------------------------

    phi = -2.0_RP*i*PI/(nmax+1)
    !$omp parallel do private(k)
    do j=-mmax/2, mmax/2
      q_out(j) = q_in(0)
      do k=1, nmax
        q_out(j) = q_out(j) + q_in(k) * exp( phi*real(j*k,kind=RP) )
      end do
      q_out(j) = q_out(j) / real(nmax+1, kind=RP)
    end do
    return
  end subroutine DFT_dft1d

!OCL SERIAL
  subroutine DFT_idft1d(q_out, q_in, nmax, mmax)
    implicit none
    integer, intent(in) :: nmax, mmax
    real(RP), intent(out) :: q_out(0:nmax)
    complex(RP), intent(in) :: q_in(-mmax/2:mmax/2)

    complex(RP) :: phi
    complex(RP), parameter :: i = (0.0, 1.0)
    integer :: j, k
    !--------------------------------------------

    phi = 2.0_RP*i*PI/(nmax+1)
    !$omp parallel do private(k)
    do j=0, nmax
      q_out(j) = 0.0_RP
      do k=-mmax/2, mmax/2
        q_out(j) = q_out(j) + real(q_in(k) * exp( phi*real(j*k,kind=RP) ), kind=RP)
      end do
    end do
    return
  end subroutine DFT_idft1d

!OCL SERIAL
  subroutine DFT_interp1d(q_out, q_in, x, mmax)
    implicit none
    integer, intent(in) :: mmax
    real(RP), intent(in) :: x(:)
    real(RP), intent(out) :: q_out(size(x))
    complex(RP), intent(in) :: q_in(-mmax/2:mmax/2)

    complex(RP) :: phi
    complex(RP), parameter :: i = (0.0, 1.0)
    integer :: j, k
    !--------------------------------------------

    phi = 2.0_RP*i*PI
    !$omp parallel do private(k)
    do j=1, size(x)
      q_out(j) = 0.0_RP
      do k=-mmax/2, mmax/2
        q_out(j) = q_out(j) + real(q_in(k) * exp( phi*x(j)*real(k,kind=RP) ), kind=RP)
      end do
    end do
    return
  end subroutine DFT_interp1d

end module mod_dft
