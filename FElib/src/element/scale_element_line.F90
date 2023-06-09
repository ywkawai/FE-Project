!> module FElib / Element / line
!!
!! @par Description
!!           A module for a line finite element
!!
!! @author Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_element_line

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_element_base, only: &
    ElementBase1D, &
    ElementBase1D_Init, ElementBase1D_Final
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  type, public, extends(ElementBase1D) :: LineElement
  contains
    procedure :: Init => LineElement_Init
    procedure :: Final => LineElement_Final
    procedure :: GenIntGaussLegendreIntrpMat => LineElement_gen_IntGaussLegendreIntrpMat    
  end type LineElement

contains
!OCL SERIAL
  subroutine LineElement_Init( &
    elem, elemOrder,           &
    LumpedMassMatFlag )
    
    implicit none

    class(LineElement), intent(inout) :: elem
    integer, intent(in) :: elemOrder
    logical, intent(in) :: LumpedMassMatFlag

    !-----------------------------------------------------------------------------
    
    elem%PolyOrder = elemOrder
    elem%Nv = 2
    elem%Np = elemOrder + 1
    elem%Nfp = 1
    elem%Nfaces = 2
    elem%NfpTot = elem%Nfp*elem%Nfaces
       
    call ElementBase1D_Init(elem, LumpedMassMatFlag)
    call construct_Element(elem)

    return
  end subroutine LineElement_Init

!OCL SERIAL  
  subroutine LineElement_Final(elem)
    implicit none

    class(LineElement), intent(inout) :: elem
    !-----------------------------------------------------------------------------

    call ElementBase1D_Final(elem)

    return
  end subroutine LineElement_Final

!OCL SERIAL  
  subroutine construct_Element(elem)

    use scale_linalgebra, only: linalgebra_inv
    use scale_polynominal, only: &
      polynominal_genGaussLobattoPt, Polynominal_GenGaussLobattoPtIntWeight,   &
      polynominal_genLegendrePoly, Polynominal_genDLegendrePoly,               &
      polynominal_genLagrangePoly, polynominal_genDLagrangePoly_lglpt

    implicit none

    type(LineElement), intent(inout) :: elem

    integer :: nodes(elem%Np)

    real(RP) :: lglPts1D(elem%Np)
    real(DP) :: intWeight_lgl1DPts(elem%Np)

    real(RP) :: P1D_ori(elem%Np, elem%Np)
    real(RP) :: DP1D_ori(elem%Np, elem%Np)
    real(RP) :: DLagr1D(elem%Np, elem%Np)
    real(RP) :: Emat(elem%Np, elem%Nfp*elem%Nfaces)
    real(RP) :: MassEdge(elem%Nfp, elem%Nfp)

    integer :: i
    integer :: p1
    integer :: n, l, f
    integer :: Nord

    !-----------------------------------------------------------------------------

    lglPts1D(:)      = polynominal_genGaussLobattoPt( elem%PolyOrder )

    P1D_ori(:,:)     = polynominal_genLegendrePoly( elem%PolyOrder, lglPts1D )
    DP1D_ori(:,:) = polynominal_genDLegendrePoly( elem%PolyOrder, lglPts1D, P1D_ori )
    DLagr1D(:,:) = polynominal_GenDLagrangePoly_lglpt(elem%PolyOrder, lglPts1D)
    
    !* Preparation 
    
    do i=1, elem%Np
        nodes(i) = i
    end do

    ! Set the mask to extract the values at faces
    
    elem%Fmask(:,1) = 1
    elem%Fmask(:,2) = elem%Np

    !* Set the coordinates of LGL points, and the Vandermonde and differential matricies

    elem%Dx1(:,:) = 0.0_RP
  
    do n=1, elem%Np      
      !* Set the coordinates of LGL points
      elem%x1(n) = lglPts1D(n)

      !* Set the Vandermonde and differential matricies
      do l=1, elem%Np
        elem%V(n,l) = P1D_ori(n,l) * sqrt(dble(l-1) + 0.5_RP)
        elem%Dx1(n,l) = DLagr1D(l,n)
      end do
    end do
    elem%invV(:,:) = linAlgebra_inv(elem%V)
    
    !* Set the weights at LGL points to integrate over element
    elem%IntWeight_lgl(:) = Polynominal_GenGaussLobattoPtIntWeight(elem%PolyOrder)

    !* Set the mass matrix

    if (elem%IsLumpedMatrix()) then
      elem%invM(:,:) = 0.0_RP
      elem%M(:,:)    = 0.0_RP
      do i=1, elem%Np
        elem%M(i,i) = elem%IntWeight_lgl(i)
        elem%invM(i,i) = 1.0_RP/elem%IntWeight_lgl(i)
      end do      
    else
      elem%invM(:,:) = matmul(elem%V, transpose(elem%V))
      elem%M(:,:) = linAlgebra_inv( elem%invM )
    end if

    !* Set the stiffness matrix
    elem%Sx1(:,:) = transpose(matmul( elem%M, elem%Dx1))
    elem%Sx1(:,:) = matmul( elem%invM, elem%Sx1 )

    !* Set the lift matrix

    Emat(:,:) = 0.0_RP  
    do f=1, elem%Nfaces
      MassEdge(:,:) = 0.0_RP
      do l=1, elem%Nfp
        MassEdge(l,l) = 1.0_RP
      end do  
      Emat(elem%Fmask(:,f), (f-1)*elem%Nfp+1:f*elem%Nfp) = MassEdge
    end do
    elem%Lift(:,:) = matmul( elem%invM, Emat )

    return
  end subroutine construct_Element

!OCL SERIAL  
  function LineElement_gen_IntGaussLegendreIntrpMat( this, IntrpPolyOrder, &
    intw_intrp, x_intrp ) result(IntrpMat)

    use scale_polynominal, only: &
      Polynominal_genLegendrePoly,    &
      Polynominal_GenGaussLegendrePt, &
      Polynominal_GenGaussLegendrePtIntWeight
    
    implicit none

    class(LineElement), intent(in) :: this
    integer, intent(in) :: IntrpPolyOrder
    real(RP), intent(out), optional :: intw_intrp(IntrpPolyOrder)
    real(RP), intent(out), optional :: x_intrp(IntrpPolyOrder)
    real(RP) :: IntrpMat(IntrpPolyOrder,this%Np)

    real(RP) :: r_int1D_i(IntrpPolyOrder)
    real(RP) :: r_int1Dw_i(IntrpPolyOrder)
    real(RP) :: P_int1D_ori(IntrpPolyOrder,this%PolyOrder+1)
    real(RP) :: Vint(IntrpPolyOrder,this%PolyOrder+1)

    integer :: p1, p1_
    !-----------------------------------------------------

    r_int1D_i(:) = Polynominal_GenGaussLegendrePt( IntrpPolyOrder )
    r_int1Dw_i(:) = Polynominal_GenGaussLegendrePtIntWeight( IntrpPolyOrder )
    P_int1D_ori(:,:) = Polynominal_GenLegendrePoly( this%PolyOrder, r_int1D_i)

    do p1_=1, IntrpPolyOrder
      if (present(intw_intrp)) intw_intrp(p1_) = r_int1Dw_i(p1_)
      if (present(x_intrp)) x_intrp(p1_) = r_int1D_i(p1_)
      do p1=1, this%Np
        Vint(p1_,p1) =  P_int1D_ori(p1_,p1) * sqrt(real(p1-1,kind=RP) + 0.5_RP)
      end do
    end do
    IntrpMat(:,:) = matmul(Vint, this%invV)

    return
  end function LineElement_gen_IntGaussLegendreIntrpMat

end module scale_element_Line