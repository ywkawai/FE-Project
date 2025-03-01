
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
    class(ElementBase3D), pointer :: elem3D 
  contains
    procedure(ElementOperationBase_Final), deferred, public :: Final
    procedure(ElementOperationBase_Dx), deferred, public :: Dx
    procedure(ElementOperationBase_Dy), deferred, public :: Dy
    procedure(ElementOperationBase_Dz), deferred, public :: Dz
    procedure(ElementOperationBase_Lift), deferred, public :: Lift
    procedure(ElementOperationBase_DxDyDzLift), deferred, public :: DxDyDzLift
    procedure(ElementOperationBase_Div), deferred, public :: Div
    procedure(ElementOperationBase_VFilterPM1), deferred, public :: VFilterPM1
  end type ElementOperationBase3D
  interface
    subroutine ElementOperationBase_Final( this )
      import ElementOperationBase3D
      class(ElementOperationBase3D), intent(inout) :: this
    end subroutine ElementOperationBase_Final

    subroutine ElementOperationBase_Dx( this, vec_in, vec_out )
      import ElementOperationBase3D
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%elem3D%Np)
      real(RP), intent(out) :: vec_out(this%elem3D%Np)
    end subroutine ElementOperationBase_Dx

    subroutine ElementOperationBase_Dy( this, vec_in, vec_out )
      import ElementOperationBase3D
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%elem3D%Np)
      real(RP), intent(out) :: vec_out(this%elem3D%Np)
    end subroutine ElementOperationBase_Dy

    subroutine ElementOperationBase_Dz( this, vec_in, vec_out )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%elem3D%Np)
      real(RP), intent(out) :: vec_out(this%elem3D%Np)
    end subroutine ElementOperationBase_Dz

    subroutine ElementOperationBase_Lift( this, vec_in, vec_out )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%elem3D%NfpTot)
      real(RP), intent(out) :: vec_out(this%elem3D%Np)
    end subroutine ElementOperationBase_Lift

    subroutine ElementOperationBase_DxDyDzLift( this, vec_in, vec_in_lift, vec_out_dx, vec_out_dy, vec_out_dz, vec_out_lift )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%elem3D%Np)
      real(RP), intent(in) :: vec_in_lift(this%elem3D%NfpTot)
      real(RP), intent(out) :: vec_out_dx(this%elem3D%Np)
      real(RP), intent(out) :: vec_out_dy(this%elem3D%Np)
      real(RP), intent(out) :: vec_out_dz(this%elem3D%Np)
      real(RP), intent(out) :: vec_out_lift(this%elem3D%Np)
    end subroutine ElementOperationBase_DxDyDzLift

    subroutine ElementOperationBase_Div( this, vec_in_x, vec_in_y, vec_in_z, vec_in_lift, Escale, Gsqrt, sign_, &
      vec_out_dx, vec_out_dy, vec_out_dz, vec_out_lift, vec_out )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in_x(this%elem3D%Np)
      real(RP), intent(in) :: vec_in_y(this%elem3D%Np)
      real(RP), intent(in) :: vec_in_z(this%elem3D%Np)
      real(RP), intent(in) :: vec_in_lift(this%elem3D%NfpTot)
      real(RP), intent(in) :: Escale(3,this%elem3D%Np)
      real(RP), intent(in) :: Gsqrt(this%elem3D%Np)
      real(RP), intent(in) :: sign_
      real(RP), intent(out) :: vec_out_dx(this%elem3D%Np)
      real(RP), intent(out) :: vec_out_dy(this%elem3D%Np)
      real(RP), intent(out) :: vec_out_dz(this%elem3D%Np)
      real(RP), intent(out) :: vec_out_lift(this%elem3D%Np)
      real(RP), intent(out) :: vec_out(this%elem3D%Np)
    end subroutine ElementOperationBase_Div
    
    subroutine ElementOperationBase_VFilterPM1( this, vec_in, vec_out )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%elem3D%Np)
      real(RP), intent(out) :: vec_out(this%elem3D%Np)
    end subroutine ElementOperationBase_VFilterPM1    
  end interface
contains
end module scale_element_operation_base