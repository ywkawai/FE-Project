!-------------------------------------------------------------------------------
!> module common / GMRES
!!
!! @par Description
!!      A module to provide a iterative solver for system of linear equations using generalized minimal residual method (GMRES)
!!
!! @par Reference
!!  - Y. Saad and M.H. Schultz, 1986:
!!     GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems.
!!     SIAM J. Sci. Stat. Comput., 7:856–869
!!
!! @author Yuta Kawai, Team SCALE
!!
#include "scaleFElib.h"
module scale_gmres
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_prc
 
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type
  !  

  type, public :: GMRES
    real(RP) :: EPS
    real(RP) :: EPS0
    integer :: N
    integer :: m
    integer :: m_out
    real(RP), allocatable :: v(:,:)
    real(RP), allocatable :: hj(:)
    real(RP), allocatable :: g(:)
    real(RP), allocatable :: r(:,:)
    real(RP), allocatable :: co(:)
    real(RP), allocatable :: si(:)
    real(RP), allocatable :: y(:)
  contains 
    procedure, public :: Init => GMRES_Init
    procedure, public :: Iterate_pre => GMRES_Iterate_pre
    procedure, public :: Iterate_step_j => GMRES_Iterate_step_j
    procedure, public :: Iterate_post => GMRES_Iterate_post
    procedure, public :: Final => GMRES_Final
  end type

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

  ! Private procedure
  !

  ! Private variable
  !

contains

subroutine GMRES_Init( this, N, m, eps, eps0 )

  implicit none
  class(GMRES), intent(inout) :: this
  integer, intent(in) :: N
  integer, intent(in) :: m
  real(RP), intent(in) :: eps
  real(RP), intent(in) :: eps0

  !--------------------------------------

  this%N = N
  this%m = m
  this%EPS = eps
  this%EPS0 = eps0

  allocate( this%v(N,m+1) )
  allocate( this%hj(m+1) )
  allocate( this%g(m+1) )
  allocate( this%r(m+1,m) )
  allocate( this%co(m), this%si(m) )
  allocate( this%y(m) )

  return
end subroutine GMRES_Init

!OCL SERIAL
subroutine GMRES_Iterate_pre( this, b, w0, is_converged )
  implicit none
  class(GMRES), intent(inout) :: this
  real(RP), intent(in) :: b(this%N)
  real(RP), intent(in) :: w0(this%N)
  logical, intent(out) :: is_converged
  
  integer :: i
  !--------------------------------------

  this%g(1) = 0.0_RP
  do i=1, this%N
    this%v(i,1) = b(i) - w0(i)
    this%g(1) = this%g(1) + this%v(i,1)**2
  end do
  this%g(1) = sqrt(this%g(1))

  if ( this%g(1) < this%EPS0 * this%N ) then
    is_converged = .true.
    return
  else
    is_converged = .false.
  end if

  do i=1, this%N
    this%v(i,1) = this%v(i,1) / this%g(1)
  end do

  this%m_out = min(this%m, this%N)
  return
end subroutine GMRES_Iterate_pre

!OCL SERIAL
subroutine GMRES_Iterate_step_j( this, j, wj, is_converged )
  implicit none
  class(GMRES), intent(inout) :: this
  integer, intent(in) :: j
  real(RP), intent(inout) :: wj(this%N)
  logical, intent(out) :: is_converged
  
  integer :: i
  real(RP) :: tmp1, tmp2
  !--------------------------------------

  do i=1, j
    this%hj(i) = sum( wj(:) * this%v(:,i) )
    wj(:) = wj(:) - this%hj(i) * this%v(:,i)
  end do
  this%hj(j+1) = sqrt( sum(wj(:)**2) )

  ! if ( abs(this%hj(j+1)) < this%EPS0 ) then
  !   this%m_out = j
  !   is_converged = .true.
  !   return
  ! else
    this%v(:,j+1) = wj(:) / this%hj(j+1)
  ! end if

  this%r(1,j) = this%hj(1)
  do i=1, j-1
    tmp1 =   this%co(i) * this%r(i,j) + this%si(i) * this%hj(i+1)
    tmp2 = - this%si(i) * this%r(i,j) + this%co(i) * this%hj(i+1)
    this%r(i  ,j) = tmp1
    this%r(i+1,j) = tmp2
  end do

  tmp1 = 1.0_RP / sqrt(this%r(j,j)**2 + this%hj(j+1)**2)
  this%co(j) = tmp1 * this%r(j,j)
  this%si(j) = tmp1 * this%hj(j+1)

  this%g(j+1) = - this%si(j) * this%g(j)
  this%g(j)   =   this%co(j) * this%g(j)

  this%r(j,j) = this%co(j) * this%r(j,j) + this%si(j) * this%hj(j+1)
  this%r(j+1,j) = 0.0_RP

  if ( abs(this%g(j+1)) < this%EPS ) then
    this%m_out = j
    is_converged = .true.
  else
    is_converged = .false.
  end if

  return
end subroutine GMRES_Iterate_step_j

!OCL SERIAL
subroutine GMRES_Iterate_post( this, x )
  implicit none
  class(GMRES), intent(inout) :: this
  real(RP), intent(inout) :: x(this%N)
  
  integer :: i, j
  !--------------------------------------

  do j=this%m_out, 1, -1
    this%y(j) = this%g(j)
    do i=j+1, this%m_out
      this%y(j) = this%y(j) - this%r(j,i) * this%y(i)
    end do
    this%y(j) = this%y(j) / this%r(j,j)
  end do

  x(:) = x(:) + matmul(this%v(:,1:this%m_out), this%y(1:this%m_out))

  return
end subroutine GMRES_Iterate_post

subroutine GMRES_Final( this )
  implicit none
  class(GMRES), intent(inout) :: this
  !--------------------------------------

  if ( allocated(this%v) ) then
    deallocate( this%v )
    deallocate( this%hj )
    deallocate( this%g )
    deallocate( this%r )
    deallocate( this%co, this%si )
    deallocate( this%y )
  end if

  return
end subroutine GMRES_Final

end module scale_gmres
