
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
    procedure(ElementOperationBase_Div_var5), deferred, public :: Div_var5
    procedure(ElementOperationBase_VFilterPM1), deferred, public :: VFilterPM1
    procedure(ElementOperationBase_Setup_ModalFilter), deferred, public :: Setup_ModalFilter
    procedure(ElementOperationBase_Setup_ModalFilter_tracer), deferred, public :: Setup_ModalFilter_tracer
    procedure(ElementOperationBase_ModalFilter_tracer), deferred, public :: ModalFilter_tracer
    procedure(ElementOperationBase_ModalFilter_var5), deferred, public :: ModalFilter_var5
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

    subroutine ElementOperationBase_Div( this, vec_in_x, vec_in_y, vec_in_z, vec_in_lift, &
      vec_out_dx, vec_out_dy, vec_out_dz, vec_out_lift )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in_x(this%elem3D%Np)
      real(RP), intent(in) :: vec_in_y(this%elem3D%Np)
      real(RP), intent(in) :: vec_in_z(this%elem3D%Np)
      real(RP), intent(in) :: vec_in_lift(this%elem3D%NfpTot)
      real(RP), intent(out) :: vec_out_dx(this%elem3D%Np)
      real(RP), intent(out) :: vec_out_dy(this%elem3D%Np)
      real(RP), intent(out) :: vec_out_dz(this%elem3D%Np)
      real(RP), intent(out) :: vec_out_lift(this%elem3D%Np)
    end subroutine ElementOperationBase_Div
    
    subroutine ElementOperationBase_Div_var5( this, vec_in, vec_in_lift, &
      vec_out_d )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%elem3D%Np,3,5)
      real(RP), intent(in) :: vec_in_lift(this%elem3D%NfpTot,5)
      real(RP), intent(out) :: vec_out_d(this%elem3D%Np,4,5)
    end subroutine ElementOperationBase_Div_var5

    subroutine ElementOperationBase_VFilterPM1( this, vec_in, vec_out )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%elem3D%Np)
      real(RP), intent(out) :: vec_out(this%elem3D%Np)
    end subroutine ElementOperationBase_VFilterPM1    

    subroutine ElementOperationBase_Setup_ModalFilter( this, &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h, &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(inout) :: this
      real(RP), intent(in) :: MF_ETAC_h
      real(RP), intent(in) :: MF_ALPHA_h
      integer, intent(in) :: MF_ORDER_h
      real(RP), intent(in) :: MF_ETAC_v
      real(RP), intent(in) :: MF_ALPHA_v
      integer, intent(in) :: MF_ORDER_v
    end subroutine ElementOperationBase_Setup_ModalFilter

    subroutine ElementOperationBase_Setup_ModalFilter_tracer( this, &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h, &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(inout) :: this
      real(RP), intent(in) :: MF_ETAC_h
      real(RP), intent(in) :: MF_ALPHA_h
      integer, intent(in) :: MF_ORDER_h
      real(RP), intent(in) :: MF_ETAC_v
      real(RP), intent(in) :: MF_ALPHA_v
      integer, intent(in) :: MF_ORDER_v
    end subroutine ElementOperationBase_Setup_ModalFilter_tracer

    subroutine ElementOperationBase_ModalFilter_tracer( this, vec_in, vec_work, vec_out )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%elem3D%Np)
      real(RP), intent(out) :: vec_work(this%elem3D%Np)
      real(RP), intent(out) :: vec_out(this%elem3D%Np)
    end subroutine ElementOperationBase_ModalFilter_tracer

    subroutine ElementOperationBase_ModalFilter_var5( this, vec_in, vec_work, vec_out )
      import ElementOperationBase3D      
      import RP
      class(ElementOperationBase3D), intent(in) :: this
      real(RP), intent(in) :: vec_in(this%elem3D%Np,5)
      real(RP), intent(out) :: vec_work(this%elem3D%Np)
      real(RP), intent(out) :: vec_out(this%elem3D%Np,5)
    end subroutine ElementOperationBase_ModalFilter_var5

  end interface
contains
end module scale_element_operation_base