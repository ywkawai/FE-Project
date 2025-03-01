
!-------------------------------------------------------------------------------
!> module FElib / Element / Operation with arbitary elements
!!
!! @par Description
!!           A module for providing mathematical operations with arbitary elements using a module for SpMV
!!
!! @author Yuta Kawai, Xuanzhengbo Ren, and Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_element_operation_general

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_sparsemat, only: &
    SparseMat, sparsemat_matmul

  use scale_element_base, only: &
    ElementBase3D, &
    ElementBase3D_Init, ElementBase3D_Final

  use scale_element_operation_base, only: ElementOperationBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  type, public, extends(ElementOperationBase3D) :: ElementOperationGenral
    type(sparsemat), pointer :: Dx_sm
    type(sparsemat), pointer :: Dy_sm
    type(sparsemat), pointer :: Dz_sm
    type(sparsemat), pointer :: Lift_sm
    real(RP), allocatable :: IntrpMat_VPOrdM1(:,:)
  contains
    procedure, public :: Init => element_operation_general_Init
    procedure, public :: Final => element_operation_general_Final
    procedure, public :: Dx => element_operation_general_Dx
    procedure, public :: Dy => element_operation_general_Dy
    procedure, public :: Dz => element_operation_general_Dz
    procedure, public :: Lift => element_operation_general_Lift
    procedure, public :: DxDyDzLift => element_operation_general_DxDyDzLift
    procedure, public :: Div => element_operation_general_Div
    procedure, public :: VFilterPM1 => element_operation_general_VFilterPM1
  end type ElementOperationGenral
  
contains

  !> Initialization
  !!
  !OCL SERIAL
  subroutine element_operation_general_Init( this, elem3D, &
      Dx, Dy, Dz, Lift )
    implicit none
    class(ElementOperationGenral), intent(inout) :: this
    class(ElementBase3D), intent(in), target :: elem3D
    type(SparseMat), intent(in), target :: Dx
    type(SparseMat), intent(in), target :: Dy
    type(SparseMat), intent(in), target :: Dz
    type(SparseMat), intent(in), target :: Lift

    integer :: p1, p2, p_
    real(RP) :: invV_VPOrdM1(elem3D%Np,elem3D%Np)
    !----------------------------------------------------------

    this%elem3D => elem3D
    this%Dx_sm => Dx
    this%Dy_sm => Dy
    this%Dz_sm => Dz
    this%Lift_sm => Lift

    !--
    allocate( this%IntrpMat_VPOrdM1(elem3D%Np,elem3D%Np) )

    InvV_VPOrdM1(:,:) = elem3D%invV(:,:)
    do p2=1, elem3D%Nnode_h1D
    do p1=1, elem3D%Nnode_h1D
      p_ = p1 + (p2-1)*elem3D%Nnode_h1D + (elem3D%Nnode_v-1)*elem3D%Nnode_h1D**2
      InvV_VPOrdM1(p_,:) = 0.0_RP
    end do
    end do
    this%IntrpMat_VPOrdM1(:,:) = matmul(elem3D%V, invV_VPOrdM1)

    return
  end subroutine element_operation_general_Init

  !> Finalization
  !!
  !OCL SERIAL
  subroutine element_operation_general_Final( this )
    implicit none
    class(ElementOperationGenral), intent(inout) :: this
    !----------------------------------------------------------

    nullify( this%elem3D )
    nullify( this%Dx_sm, this%Dy_sm, this%Dz_sm, this%Lift_sm )

    return
  end subroutine element_operation_general_Final

!> Calculate the differential in x-direction
!!
!OCL SERIAL
  subroutine element_operation_general_Dx( this, vec_in, vec_out )
    implicit none
    class(ElementOperationGenral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np)
    real(RP), intent(out) :: vec_out(this%elem3D%Np)
    !----------------------------------------------------------
    call sparsemat_matmul( this%Dx_sm, vec_in, vec_out )
    return
  end subroutine element_operation_general_Dx

!> Calculate the differential in y-direction
!!
!OCL SERIAL
  subroutine element_operation_general_Dy( this, vec_in, vec_out )
    implicit none
    class(ElementOperationGenral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np)
    real(RP), intent(out) :: vec_out(this%elem3D%Np)
    !----------------------------------------------------------
    call sparsemat_matmul( this%Dy_sm, vec_in, vec_out )
    return
  end subroutine element_operation_general_Dy

!> Calculate the differential in z-direction
!!
!OCL SERIAL
  subroutine element_operation_general_Dz( this, vec_in, vec_out )
    implicit none
    class(ElementOperationGenral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np)
    real(RP), intent(out) :: vec_out(this%elem3D%Np)
    !----------------------------------------------------------
    call sparsemat_matmul( this%Dz_sm, vec_in, vec_out )
    return
  end subroutine element_operation_general_Dz  

!> Calculate the differential in z-direction
!!
!OCL SERIAL
  subroutine element_operation_general_Lift( this, vec_in, vec_out )
    implicit none
    class(ElementOperationGenral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np)
    real(RP), intent(out) :: vec_out(this%elem3D%Np)
    !----------------------------------------------------------
    call sparsemat_matmul( this%Lift_sm, vec_in, vec_out )
    return
  end subroutine element_operation_general_Lift 

!> Calculate the 3D gradient
!!
!OCL SERIAL
  subroutine element_operation_general_DxDyDzLift( this, vec_in, vec_in_lift, vec_out_dx, vec_out_dy, vec_out_dz, vec_out_lift )
    implicit none
    class(ElementOperationGenral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np)
    real(RP), intent(in) :: vec_in_lift(this%elem3D%NfpTot)
    real(RP), intent(out) :: vec_out_dx(this%elem3D%Np)
    real(RP), intent(out) :: vec_out_dy(this%elem3D%Np)
    real(RP), intent(out) :: vec_out_dz(this%elem3D%Np)
    real(RP), intent(out) :: vec_out_lift(this%elem3D%Np)
    !----------------------------------------------------------

    call sparsemat_matmul( this%Dx_sm, vec_in, vec_out_dx )
    call sparsemat_matmul( this%Dy_sm, vec_in, vec_out_dy )
    call sparsemat_matmul( this%Dz_sm, vec_in, vec_out_dz )
    call sparsemat_matmul( this%Lift_sm, vec_in_lift, vec_out_lift )
    return
  end subroutine element_operation_general_DxDyDzLift

!> Calculate the 3D gradient
!!
!OCL SERIAL
  subroutine element_operation_general_Div( this, vec_in_x, vec_in_y, vec_in_z, vec_in_lift, Escale, Gsqrt, sign_, &
    vec_out_dx, vec_out_dy, vec_out_dz, vec_out_lift, vec_out )
    implicit none
    class(ElementOperationGenral), intent(in) :: this
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

    integer :: p
    !---------------------------------------------------------------

    call sparsemat_matmul( this%Dx_sm, vec_in_x, vec_out_dx )
    call sparsemat_matmul( this%Dy_sm, vec_in_y, vec_out_dy )
    call sparsemat_matmul( this%Dz_sm, vec_in_z, vec_out_dz )
    call sparsemat_matmul( this%Lift_sm, vec_in_lift, vec_out_lift )

    do p=1, this%elem3D%Np
      vec_out(p) = sign_ * ( &
                   Escale(1,p) * vec_out_dx(p) + Escale(2,p) * vec_out_dy(p) + Escale(3,p) * vec_out_dz(p) &
                 + vec_out_lift(p) ) / Gsqrt(p)
    end do

    return
  end subroutine element_operation_general_Div    

  subroutine element_operation_general_VFilterPM1( this, vec_in, vec_out )
    implicit none
    class(ElementOperationGenral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np)
    real(RP), intent(out) :: vec_out(this%elem3D%Np)
    !---------------------------------------------------------------

    vec_out(:) = matmul( this%IntrpMat_VPOrdM1, vec_in(:) )
    return
  end subroutine element_operation_general_VFilterPM1 

end module scale_element_operation_general