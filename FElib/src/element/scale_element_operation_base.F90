
!-------------------------------------------------------------------------------
!> module FElib / Element / Operation / Base
!!
!! @par Description
!!           A base module for providing mathematical operations
!!
!! @author Yuta Kawai, Xuanzhengbo Ren, and Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_element_operation_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_sparsemat, only: &
    SparseMat
  use scale_element_base, only: &
    ElementBase3D, &
    ElementBase3D_Init, ElementBase3D_Final
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  type, public, abstract :: ElementOperationBase3D
    type(SparseMat), pointer :: Dx_sm
    type(SparseMat), pointer :: Dy_sm
    type(SparseMat), pointer :: Dz_sm
    type(SparseMat), pointer :: Lift_sm
    integer :: Np1D
  contains
    procedure(ElementOperationBase_Dx), deferred, public :: Dx
    procedure(ElementOperationBase_Dy), deferred, public :: Dy
    procedure(ElementOperationBase_Dz), deferred, public :: Dz
  end type ElementOperationBase3D
  interface
    subroutine ElementOperationBase_Dx( this, vec_in, vec_out )
      import ElementOperationBase3D
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%Np1D,this%Np1D**2)
      real(RP), intent(out) :: vec_out(this%Np1D,this%Np1D**2)
    end subroutine ElementOperationBase_Dx

    subroutine ElementOperationBase_Dy( this, vec_in, vec_out )
      import ElementOperationBase3D
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%Np1D,this%Np1D,this%Np1D)
      real(RP), intent(out) :: vec_out(this%Np1D,this%Np1D,this%Np1D)
    end subroutine ElementOperationBase_Dy

    subroutine ElementOperationBase_Dz( this, vec_in, vec_out )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%Np1D,this%Np1D,this%Np1D)
      real(RP), intent(out) :: vec_out(this%Np1D,this%Np1D,this%Np1D)
    end subroutine ElementOperationBase_Dz    
  end interface
contains
end module scale_element_operation_base