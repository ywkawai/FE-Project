
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

  use scale_element_base, only: &
    ElementBase3D, &
    ElementBase3D_Init, ElementBase3D_Final

  use scale_sparsemat, only: &
    SparseMat, sparsemat_matmul
  use scale_element_modalfilter, only: ModalFilter  

  use scale_element_operation_base, only: ElementOperationBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  type, public, extends(ElementOperationBase3D) :: ElementOperationGeneral
    type(SparseMat), pointer :: Dx_sm
    type(SparseMat), pointer :: Dy_sm
    type(SparseMat), pointer :: Dz_sm
    type(SparseMat), pointer :: Lift_sm
    real(RP), allocatable :: IntrpMat_VPOrdM1(:,:)

    type(ModalFilter) :: MFilter
    type(ModalFilter) :: MFilter_tracer
  contains
    procedure, public :: Init => element_operation_general_Init
    procedure, public :: Final => element_operation_general_Final
    procedure, public :: Dx => element_operation_general_Dx
    procedure, public :: Dy => element_operation_general_Dy
    procedure, public :: Dz => element_operation_general_Dz
    procedure, public :: Lift => element_operation_general_Lift
    procedure, public :: DxDyDzLift => element_operation_general_DxDyDzLift
    procedure, public :: Div => element_operation_general_Div
    procedure, public :: Div_var5 => element_operation_general_Div_var5
    procedure, public :: VFilterPM1 => element_operation_general_VFilterPM1
    !-
    procedure, public :: Setup_ModalFilter => element_operation_general_Setup_ModalFilter
    procedure, public :: Setup_ModalFilter_tracer => element_operation_general_Setup_ModalFilter_tracer
    procedure, public :: ModalFilter_tracer => element_operation_general_ModalFilter_tracer
    procedure, public :: ModalFilter_var5 => element_operation_general_ModalFilter_var5
  end type ElementOperationGeneral
  
  interface ElementOperationGeneral_Generate_VPOrdM1
    module procedure element_operation_general_generate_VPOrdM1
  end interface
  public :: ElementOperationGeneral_Generate_VPOrdM1
  
contains

  !> Initialization
  !!
!OCL SERIAL
  subroutine element_operation_general_Init( this, elem3D, &
      Dx, Dy, Dz, Lift )
    implicit none
    class(ElementOperationGeneral), intent(inout) :: this
    class(ElementBase3D), intent(in), target :: elem3D
    type(SparseMat), intent(in), target :: Dx
    type(SparseMat), intent(in), target :: Dy
    type(SparseMat), intent(in), target :: Dz
    type(SparseMat), intent(in), target :: Lift
    !----------------------------------------------------------

    this%elem3D => elem3D
    this%Dx_sm => Dx
    this%Dy_sm => Dy
    this%Dz_sm => Dz
    this%Lift_sm => Lift

    !--
    allocate( this%IntrpMat_VPOrdM1(elem3D%Np,elem3D%Np) )
    call element_operation_general_generate_VPOrdM1( this%IntrpMat_VPOrdM1, &
      elem3D )

    return
  end subroutine element_operation_general_Init

  !> Generate a vertical filter matrix to remove the highest mode
  !!
!OCL SERIAL
  subroutine element_operation_general_generate_VPOrdM1( IntrpMat_VPOrdM1, &
    elem3D )
    implicit none
    class(ElementBase3D), intent(in), target :: elem3D
    real(RP), intent(out) :: IntrpMat_VPOrdM1(elem3D%Np,elem3D%Np)

    integer :: p1, p2, p_
    real(RP) :: invV_VPOrdM1(elem3D%Np,elem3D%Np)
    !----------------------------------------------------------

    InvV_VPOrdM1(:,:) = elem3D%invV(:,:)
    do p2=1, elem3D%Nnode_h1D
    do p1=1, elem3D%Nnode_h1D
      p_ = p1 + (p2-1)*elem3D%Nnode_h1D + (elem3D%Nnode_v-1)*elem3D%Nnode_h1D**2
      InvV_VPOrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem3D%V, invV_VPOrdM1)

    return
  end subroutine element_operation_general_generate_VPOrdM1

  !> Setup modal filter
  !!
!OCL SERIAL
  subroutine element_operation_general_Setup_ModalFilter( this, &
    MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h, &
    MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v )

    implicit none
    class(ElementOperationGeneral), intent(inout) :: this
    real(RP), intent(in) :: MF_ETAC_h
    real(RP), intent(in) :: MF_ALPHA_h
    integer, intent(in) :: MF_ORDER_h
    real(RP), intent(in) :: MF_ETAC_v
    real(RP), intent(in) :: MF_ALPHA_v
    integer, intent(in) :: MF_ORDER_v
    !--------------------------------------------------------

    call setup_ModalFilter( this%MFilter, &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,         &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v,         &
      this%elem3D%PolyOrder_h, this%elem3D%PolyOrder_v  )
  
    return
  end subroutine element_operation_general_Setup_ModalFilter  

  !> Setup modal filter for tracer
  !!
!OCL SERIAL
  subroutine element_operation_general_Setup_ModalFilter_tracer( this, &
    MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h, &
    MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v )
    implicit none
    class(ElementOperationGeneral), intent(inout) :: this
    real(RP), intent(in) :: MF_ETAC_h
    real(RP), intent(in) :: MF_ALPHA_h
    integer, intent(in) :: MF_ORDER_h
    real(RP), intent(in) :: MF_ETAC_v
    real(RP), intent(in) :: MF_ALPHA_v
    integer, intent(in) :: MF_ORDER_v
    !--------------------------------------------------------

    call setup_ModalFilter( this%MFilter_tracer, &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,         &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v,         &
      this%elem3D%PolyOrder_h, this%elem3D%PolyOrder_v  )
    
    return
  end subroutine element_operation_general_Setup_ModalFilter_tracer
  

  !> Finalization
  !!
  !OCL SERIAL
  subroutine element_operation_general_Final( this )
    implicit none
    class(ElementOperationGeneral), intent(inout) :: this
    !----------------------------------------------------------

    nullify( this%elem3D )
    nullify( this%Dx_sm, this%Dy_sm, this%Dz_sm, this%Lift_sm )

    deallocate( this%IntrpMat_VPOrdM1 )

    return
  end subroutine element_operation_general_Final

!> Calculate the differential in x-direction
!!
!OCL SERIAL
  subroutine element_operation_general_Dx( this, vec_in, vec_out )
    implicit none
    class(ElementOperationGeneral), intent(in) :: this
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
    class(ElementOperationGeneral), intent(in) :: this
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
    class(ElementOperationGeneral), intent(in) :: this
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
    class(ElementOperationGeneral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%NfpTot)
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
    class(ElementOperationGeneral), intent(in) :: this
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
  subroutine element_operation_general_Div( this, vec_in, vec_in_lift, &
    vec_out )
    implicit none
    class(ElementOperationGeneral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np,3)
    real(RP), intent(in) :: vec_in_lift(this%elem3D%NfpTot)
    real(RP), intent(out) :: vec_out(this%elem3D%Np,4)
    !---------------------------------------------------------------

    call sparsemat_matmul( this%Dx_sm, vec_in(:,1), vec_out(:,1) )
    call sparsemat_matmul( this%Dy_sm, vec_in(:,2), vec_out(:,2) )
    call sparsemat_matmul( this%Dz_sm, vec_in(:,3), vec_out(:,3) )
    call sparsemat_matmul( this%Lift_sm, vec_in_lift, vec_out(:,4) )
    return
  end subroutine element_operation_general_Div    


!> Calculate the 3D gradient
!!
!OCL SERIAL
  subroutine element_operation_general_Div_var5( this, vec_in, vec_in_lift, &
    vec_out_d )
    implicit none
    class(ElementOperationGeneral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np,3,5)
    real(RP), intent(in) :: vec_in_lift(this%elem3D%NfpTot,5)
    real(RP), intent(out) :: vec_out_d(this%elem3D%Np,4,5)

    integer :: iv
    !---------------------------------------------------------------

    do iv=1, 5
      call sparsemat_matmul( this%Dx_sm, vec_in(:,1,iv), vec_out_d(:,1,iv) )
      call sparsemat_matmul( this%Dy_sm, vec_in(:,2,iv), vec_out_d(:,2,iv) )
      call sparsemat_matmul( this%Dz_sm, vec_in(:,3,iv), vec_out_d(:,3,iv) )
    end do
    do iv=1, 5
      call sparsemat_matmul( this%Lift_sm, vec_in_lift(:,iv), vec_out_d(:,4,iv) )
    end do
    return
  end subroutine element_operation_general_Div_var5

!OCL SERIAL
  subroutine element_operation_general_VFilterPM1( this, vec_in, vec_out )
    implicit none
    class(ElementOperationGeneral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np)
    real(RP), intent(out) :: vec_out(this%elem3D%Np)
    !---------------------------------------------------------------

    call matmul_( this%IntrpMat_VPOrdM1, vec_in, this%elem3D%Np, &
      vec_out )    
    return
  end subroutine element_operation_general_VFilterPM1 
!--
!OCL SERIAL
  subroutine matmul_( IntrpMat_VPOrdM1, vec_in_, Np, vec_out_ )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(in) :: IntrpMat_VPOrdM1(Np,Np)
    real(RP), intent(in) :: vec_in_(Np)
    real(RP), intent(out) :: vec_out_(Np)
    !-------------------------------------------
    vec_out_(:) = matmul( IntrpMat_VPOrdM1(:,:), vec_in_(:) )
    return
  end subroutine matmul_

!OCL SERIAL
  subroutine element_operation_general_ModalFilter_tracer( this, vec_in, vec_work, vec_out )
    implicit none
    class(ElementOperationGeneral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np)
    real(RP), intent(out) :: vec_work(this%elem3D%Np)
    real(RP), intent(out) :: vec_out(this%elem3D%Np)

    integer :: ii, kk, Np
    real(RP) :: Mik
    !---------------------------------------------

    Np = this%elem3D%Np 
    vec_out(:) = 0.0_RP

    do ii=1, Np
    do kk=1, Np
      Mik = this%MFilter_tracer%FilterMat(ii,kk)
      vec_out(ii) = vec_out(ii) + Mik * vec_in(kk)
    end do
    end do

    return
  end subroutine element_operation_general_ModalFilter_tracer

!OCL SERIAL
  subroutine element_operation_general_ModalFilter_var5( this, vec_in, vec_work, vec_out )
    implicit none
    class(ElementOperationGeneral), intent(in) :: this
    real(RP), intent(in) :: vec_in(this%elem3D%Np,5)
    real(RP), intent(out) :: vec_work(this%elem3D%Np)
    real(RP), intent(out) :: vec_out(this%elem3D%Np,5)

    integer :: ii, kk, Np
    real(RP) :: Mik
    !---------------------------------------------

    Np = this%elem3D%Np 
    vec_out(:,:) = 0.0_RP

    do ii=1, Np
    do kk=1, Np
      Mik = this%MFilter%FilterMat(ii,kk)

      vec_out(ii,1) = vec_out(ii,1) + Mik * vec_in(kk,1)
      vec_out(ii,2) = vec_out(ii,2) + Mik * vec_in(kk,2)
      vec_out(ii,3) = vec_out(ii,3) + Mik * vec_in(kk,3)
      vec_out(ii,4) = vec_out(ii,4) + Mik * vec_in(kk,4)
      vec_out(ii,5) = vec_out(ii,5) + Mik * vec_in(kk,5)
    end do
    end do

    return
  end subroutine element_operation_general_ModalFilter_var5

!- private -

!OCL SERIAL
  subroutine setup_ModalFilter( MFilter, &
    MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h, &
    MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v, &
    PolyOrder_h, PolyOrder_v )

    use scale_element_hexahedral, only: HexahedralElement
    implicit none

    class(ModalFilter), intent(inout) :: MFilter
    real(RP), intent(in) :: MF_ETAC_h
    real(RP), intent(in) :: MF_ALPHA_h
    integer, intent(in) :: MF_ORDER_h
    real(RP), intent(in) :: MF_ETAC_v
    real(RP), intent(in) :: MF_ALPHA_v
    integer, intent(in) :: MF_ORDER_v
    integer, intent(in) :: PolyOrder_h
    integer, intent(in) :: PolyOrder_v

    type(HexahedralElement) :: elem3D
    !--------------------------------------------------------

    call elem3D%Init( PolyOrder_h, PolyOrder_v, .false. )

    call MFilter%Init( &
      elem3D,                              & ! (in)
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   & ! (in)
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    ) ! (in)
    
    call elem3D%Final()
    return
  end subroutine setup_ModalFilter

end module scale_element_operation_general