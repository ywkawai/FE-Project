!-------------------------------------------------------------------------------
!> module FElib / element/ ModalFilter
!!
!! @par Description
!!      Modal filter
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_element_modalfilter
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
  
  use scale_element_line, only: LineElement
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !

  !-----------------------------------------------------------------------------
  !
  !++ Public type
  !
  type, public :: ModalFilter
    real(RP), allocatable :: FilterMat(:,:)
  contains
    procedure :: Init_line => ModalFilter_Init_line
    procedure :: Init_quadrilateral => ModalFilter_Init_quadrilateral
    procedure :: Init_hexahedral => ModalFilter_Init_hexahedral
    generic :: Init => Init_line, Init_quadrilateral, Init_hexahedral
    procedure :: Final => ModalFilter_Final
  end type ModalFilter

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  private :: get_exp_filter

contains
  subroutine ModalFilter_Init_line( this,  & ! (inout)
    elem,                                  & ! (in)
    etac, alpha, ord,                      & ! (in)
    tend_flag                              ) ! (in)

    implicit none
    class(ModalFilter), intent(inout) :: this
    class(LineElement), intent(in) :: elem
    real(RP), intent(in) :: etac
    real(RP), intent(in) :: alpha
    integer, intent(in) :: ord
    logical, intent(in), optional :: tend_flag

    real(RP) :: filter1D(elem%Np)
    integer :: p
    logical :: tend_flag_
    !----------------------------------------------------
    
    tend_flag_ = .false.
    if ( present(tend_flag) ) tend_flag_ =  tend_flag

    call get_exp_filter( filter1D,               & ! (out)
      etac, alpha, ord, elem%Np, elem%PolyOrder, & ! (in)
      tend_flag_ )                                 ! (in)

    allocate( this%FilterMat(elem%Np,elem%Np) )
    this%FilterMat(:,:) = 0.0_RP
    do p=1, elem%Np
      this%FilterMat(p,p) = filter1D(p)
    end do
    this%FilterMat(:,:) = matmul(this%FilterMat, elem%invV)
    this%FilterMat(:,:) = matmul(elem%V, this%FilterMat)
    
    return
  end subroutine ModalFilter_Init_line

  subroutine ModalFilter_Init_quadrilateral( this,   & ! (inout)
    elem,                                            & ! (in)
    etac, alpha, ord,                                & ! (in)
    tend_flag                                        ) ! (in)

    implicit none
    class(ModalFilter), intent(inout) :: this
    class(QuadrilateralElement), intent(in) :: elem
    real(RP), intent(in) :: etac
    real(RP), intent(in) :: alpha
    integer, intent(in) :: ord
    logical, intent(in), optional :: tend_flag

    real(RP) :: filter1D(elem%Nfp)
    integer :: p1, p2
    integer :: l
    logical :: tend_flag_
    !----------------------------------------------------
    
    tend_flag_ = .false.
    if ( present(tend_flag) ) tend_flag_ =  tend_flag

    call get_exp_filter( filter1D,                & ! (out)
      etac, alpha, ord, elem%Nfp, elem%PolyOrder, & ! (in)
      tend_flag_ )                                  ! (in)
    
    allocate( this%FilterMat(elem%Np,elem%Np) )
    this%FilterMat(:,:) = 0.0_RP
    do p2=1, elem%Nfp
    do p1=1, elem%Nfp
      l = p1 + (p2-1)*elem%Nfp
      this%FilterMat(l,l) = filter1D(p1) * filter1D(p2)
    end do  
    end do
    this%FilterMat(:,:) = matmul(this%FilterMat, elem%invV)
    this%FilterMat(:,:) = matmul(elem%V, this%FilterMat)
    
    return
  end subroutine ModalFilter_Init_quadrilateral

  subroutine ModalFilter_Init_hexahedral( this,  & ! (inout)
    elem,                                        & ! (in)
    etac_h, alpha_h, ord_h,                      & ! (in)
    etac_v, alpha_v, ord_v,                      & ! (in)
    tend_flag                                    ) ! (in)

    implicit none
    class(ModalFilter), intent(inout) :: this
    class(HexahedralElement), intent(in) :: elem
    real(RP), intent(in) :: etac_h
    real(RP), intent(in) :: alpha_h
    integer, intent(in) :: ord_h
    real(RP), intent(in) :: etac_v
    real(RP), intent(in) :: alpha_v
    integer, intent(in) :: ord_v
    logical, intent(in), optional :: tend_flag

    real(RP) :: filter1D_h(elem%Nnode_h1D)
    real(RP) :: filter1D_v(elem%Nnode_v)
    !--- for 3 direction
    !--- should create 1D_h and 1D_v, just for simplicity
    type(LineElement) :: elem1D
    !--- end
    integer :: p1, p2, p3
    integer :: l
    logical :: tend_flag_
    !----------------------------------------------------
    
    tend_flag_ = .false.
    if ( present(tend_flag) ) tend_flag_ =  tend_flag

    call get_exp_filter( filter1D_h,                            & ! (out)
      etac_h, alpha_h, ord_h, elem%Nnode_h1D, elem%PolyOrder_h, & ! (in)
      tend_flag_ )                                                ! (in)
    
    call get_exp_filter( filter1D_v,                           & ! (out)
      etac_v, alpha_v, ord_V, elem%Nnode_v, elem%PolyOrder_v,  & ! (in)
      tend_flag_ )                                               ! (in)

    ! -- original method
    !allocate( this%FilterMat(elem%Np,elem%Np) )
    !this%FilterMat(:,:) = 0.0_RP

    !do p3=1, elem%Nnode_v
    !do p2=1, elem%Nnode_h1D
    !do p1=1, elem%Nnode_h1D
    !  l = p1 + (p2-1)*elem%Nnode_h1D + (p3-1)*elem%Nnode_h1D**2
    !  this%FilterMat(l,l) = filter1D_h(p1) * filter1D_h(p2) * filter1D_v(p3)
    !end do  
    !end do
    !end do
    !this%FilterMat(:,:) = matmul(this%FilterMat, elem%invV)
    !this%FilterMat(:,:) = matmul(elem%V, this%FilterMat)

    !--- init a line element
    call elem1D%Init(elem%PolyOrder_h, .false.)

    allocate( this%FilterMat(elem1D%Np,elem1D%Np) )
    this%FilterMat(:,:) = 0.0_RP

    !--- assgin 1d filter to the matrix
    do p1=1, elem1D%Np
      this%FilterMat(p1,p1) = filter1D_h(p1)
    end do

    this%FilterMat(:,:) = matmul(this%FilterMat, elem1D%invV)
    this%FilterMat(:,:) = matmul(elem1D%V, this%FilterMat)

    !-- clean up

    call elem1D%Final()

    return
  end subroutine ModalFilter_Init_hexahedral

  subroutine ModalFilter_Final( this )
    implicit none

    class(ModalFilter), intent(inout) :: this
    !--------------------------------------------

    if( allocated(this%FilterMat) ) deallocate( this%FilterMat )
    
    return
  end subroutine ModalFilter_Final

!-- private --------------------------------------------------

  subroutine get_exp_filter( filter, &  
    etac, alpha, ord, Np, polyOrder, &
    tend_flag )

    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: filter(Np)
    real(RP), intent(in) :: etac  
    real(RP), intent(in) :: alpha
    integer , intent(in) :: ord
    integer , intent(in) :: polyOrder
    logical, intent(in) :: tend_flag

    integer :: p
    real(RP) :: eta
    !-----------------------------------------

    if ( tend_flag ) then
      filter(:) = 0.0_RP
    else
      filter(:) = 1.0_RP
    end if

    do p=1, Np
      eta = dble(p-1)/dble(polyOrder)
      if ( eta >  etac .and. p /= 1) then
        filter(p) = - alpha * ( ((eta - etac)/(1.0_RP - etac))**ord )
        if ( .not. tend_flag ) filter(p) = exp( filter(p) )
      end if
    end do

    return
  end subroutine get_exp_filter

end module scale_element_modalfilter
