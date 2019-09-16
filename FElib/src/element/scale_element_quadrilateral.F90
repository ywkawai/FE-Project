#include "scaleFElib.h"
module scale_element_quadrilateral

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_element_base, only: &
    elementbase2D, &
    ElementBase2D_Init, ElementBase2D_Final
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  type, public, extends(elementbase2D) :: QuadrilateralElement
  contains
    procedure :: Init => QuadrilateralElement_Init
    procedure :: Final => QuadrilateralElement_Final
    procedure :: GenIntGaussLegendreIntrpMat => QuadrilateralElement_gen_IntGaussLegendreIntrpMat
  end type QuadrilateralElement

contains

  subroutine QuadrilateralElement_Init( &
      elem, elemOrder,             &
      DumpedMassMatFlag )
    
    class(QuadrilateralElement), intent(inout) :: elem
    integer, intent(in) :: elemOrder
    logical, intent(in) :: DumpedMassMatFlag

    !-----------------------------------------------------------------------------
    
    elem%PolyOrder = elemOrder
    elem%Nv = 4
    elem%Np = (elemOrder + 1)**2
    elem%Nfp = elemOrder + 1
    elem%Nfaces = 4
    elem%NfpTot = elem%Nfp*elem%Nfaces
       
    call ElementBase2D_Init(elem)
    call construct_Element(elem, DumpedMassMatFlag)

  end subroutine QuadrilateralElement_Init

  subroutine QuadrilateralElement_Final(elem)
    class(QuadrilateralElement), intent(inout) :: elem
    !-----------------------------------------------------------------------------

    call ElementBase2D_Final(elem)
  end subroutine QuadrilateralElement_Final

  subroutine construct_Element(elem, DumpedMassMatFlag)

    use scale_linalgebra, only: linalgebra_inv
    use scale_polynominal, only: &
    polynominal_genGaussLobattoPt, Polynominal_GenGaussLobattoPtIntWeight,   &
    polynominal_genLegendrePoly, Polynominal_genDLegendrePoly,               &
    polynominal_genLagrangePoly, polynominal_genDLagrangePoly_lglpt

    type(QuadrilateralElement), intent(inout) :: elem
    logical, intent(in) :: DumpedMassMatFlag

    integer :: nodes_ij(elem%Nfp, elem%Nfp)

    real(RP) :: lglPts1D(elem%Nfp)
    real(DP) :: intWeight_lgl1DPts(elem%Nfp)

    real(RP) :: P1D_ori(elem%Nfp, elem%Nfp)
    real(RP) :: DP1D_ori(elem%Nfp, elem%Nfp)
    real(RP) :: DLagr1D(elem%Nfp, elem%Nfp)
    real(RP) :: V1D(elem%Nfp, elem%Nfp)
    real(RP) :: Emat(elem%Np, elem%Nfp*elem%Nfaces)
    real(RP) :: MassEdge(elem%Nfp, elem%Nfp)

    integer :: i, j
    integer :: p1, p2
    integer :: n, l, f
    integer :: Nord

    !-----------------------------------------------------------------------------

    lglPts1D(:)      = polynominal_genGaussLobattoPt( elem%PolyOrder )

    P1D_ori(:,:)     = polynominal_genLegendrePoly( elem%PolyOrder, lglPts1D )
    DP1D_ori(:,:) = polynominal_genDLegendrePoly( elem%PolyOrder, lglPts1D, P1D_ori )
    DLagr1D(:,:) = polynominal_GenDLagrangePoly_lglpt(elem%PolyOrder, lglPts1D)
    
    !* Preparation 
    
    do j=1, elem%Nfp
    do i=1, elem%Nfp
        nodes_ij(i,j) = i + (j-1)*elem%Nfp
    end do
    end do

    ! Set the mask to extract the values at faces
    
    elem%Fmask(:,1) = nodes_ij(:,1)
    elem%Fmask(:,2) = nodes_ij(elem%Nfp,:)
    elem%Fmask(:,3) = nodes_ij(:,elem%Nfp)
    elem%Fmask(:,4) = nodes_ij(1,:)

    !* Set the coordinates of LGL points, and the Vandermonde and differential matricies

    elem%Dx1(:,:) = 0.0_RP
    elem%Dx2(:,:) = 0.0_RP

    do j=1, elem%Nfp
    do i=1, elem%Nfp
      n = i + (j-1)*elem%Nfp
      
      !* Set the coordinates of LGL points
      elem%x1(n) = lglPts1D(i)
      elem%x2(n) = lglPts1D(j)

      !* Set the Vandermonde and differential matricies
      do p2=1, elem%Nfp
      do p1=1, elem%Nfp
        l = p1 + (p2 - 1)*elem%Nfp
        elem%V(n,l) = (P1D_ori(i,p1)*P1D_ori(j,p2))                       &
                      * sqrt((dble(p1-1) + 0.5_RP)*(dble(p2-1) + 0.5_RP))
  
        if(p2==j) elem%Dx1(n,l) = DLagr1D(p1,i)
        if(p1==i) elem%Dx2(n,l) = Dlagr1D(p2,j)
      end do
      end do
    end do
    end do
    elem%invV(:,:) = linAlgebra_inv(elem%V)
    
    !* Set the weights at LGL points to integrate over element

    intWeight_lgl1DPts(:) = Polynominal_GenGaussLobattoPtIntWeight(elem%PolyOrder)
        
    do j=1, elem%Nfp
    do i=1, elem%Nfp
        l = i + (j - 1)*elem%Nfp
        elem%IntWeight_lgl(l) = &
          intWeight_lgl1DPts(i) * intWeight_lgl1DPts(j)
    end do
    end do

    !* Set the mass matrix

    if (DumpedMassMatFlag) then
      elem%invM(:,:) = 0.0_RP
      elem%M(:,:)    = 0.0_RP
      do j=1, elem%Nfp
      do i=1, elem%Nfp
        l = i + (j - 1)*elem%Nfp
        elem%M(l,l) = elem%IntWeight_lgl(l)
        elem%invM(l,l) = 1.0_RP/elem%IntWeight_lgl(l)
      end do
      end do      
    else
      elem%invM(:,:) = matmul(elem%V, transpose(elem%V))
      elem%M(:,:) = linAlgebra_inv( elem%invM )
    end if

    !* Set the stiffness matrix
    elem%Sx1(:,:) = transpose(matmul( elem%M, elem%Dx1))
    elem%Sx1(:,:) = matmul( elem%invM, elem%Sx1 )
    
    elem%Sx2(:,:) = transpose(matmul( elem%M, elem%Dx2))
    elem%Sx2(:,:) = matmul( elem%invM, elem%Sx2 )

    !* Set the lift matrix

    do n=1, elem%Nfp
      V1D(:,n) = P1D_ori(:,n)*sqrt(dble(n-1) + 0.5_RP)
    end do

    Emat(:,:) = 0.0_RP
    do f=1, elem%Nfaces

      if (DumpedMassMatFlag) then
        MassEdge = 0.0_RP
        do l=1, elem%Nfp
          MassEdge(l,l) = intWeight_lgl1DPts(l)
        end do  
      else
        MassEdge(:,:) = linalgebra_inv(matmul(V1D, transpose(V1D)))
      end if

      Emat(elem%Fmask(:,f), (f-1)*elem%Nfp+1:f*elem%Nfp) = MassEdge
    end do
    elem%Lift(:,:) = matmul( elem%invM, Emat )
  
  end subroutine construct_Element

  function QuadrilateralElement_gen_IntGaussLegendreIntrpMat( this, IntrpPolyOrder, &
    intw_intrp, x_intrp, y_intrp ) result(IntrpMat)

    use scale_polynominal, only: &
      Polynominal_genLegendrePoly,    &
      Polynominal_GenGaussLegendrePt, &
      Polynominal_GenGaussLegendrePtIntWeight
    
    implicit none

    class(QuadrilateralElement), intent(in) :: this
    integer, intent(in) :: IntrpPolyOrder
    real(RP), intent(out), optional :: intw_intrp(IntrpPolyOrder**2)
    real(RP), intent(out), optional :: x_intrp(IntrpPolyOrder**2)
    real(RP), intent(out), optional :: y_intrp(IntrpPolyOrder**2)
    real(RP) :: IntrpMat(IntrpPolyOrder**2,this%Np)

    real(RP) :: r_int1D_i(IntrpPolyOrder)
    real(RP) :: r_int1Dw_i(IntrpPolyOrder)
    real(RP) :: P_int1D_ori(IntrpPolyOrder,this%PolyOrder+1)
    real(RP) :: Vint(IntrpPolyOrder**2,(this%PolyOrder+1)**2)

    integer :: p1, p2, p1_, p2_
    integer :: n_, l_
    !-----------------------------------------------------

    r_int1D_i(:) = Polynominal_GenGaussLegendrePt( IntrpPolyOrder )
    r_int1Dw_i(:) = Polynominal_GenGaussLegendrePtIntWeight( IntrpPolyOrder )
    P_int1D_ori(:,:) = Polynominal_GenLegendrePoly( this%PolyOrder, r_int1D_i)

    do p2_=1, IntrpPolyOrder
    do p1_=1, IntrpPolyOrder
      n_= p1_ + (p2_-1)*IntrpPolyOrder
      if (present(intw_intrp)) intw_intrp(n_) = r_int1Dw_i(p1_) * r_int1Dw_i(p2_)
      if (present(x_intrp)) x_intrp(n_) = r_int1D_i(p1_)
      if (present(y_intrp)) y_intrp(n_) = r_int1D_i(p2_)
      
      do p2=1, this%Nfp
      do p1=1, this%Nfp
        l_ = p1 + (p2-1)*this%Nfp
        Vint(n_,l_) =  P_int1D_ori(p1_,p1) * sqrt(real(p1-1,kind=RP) + 0.5_RP) &
                     * P_int1D_ori(p2_,p2) * sqrt(real(p2-1,kind=RP) + 0.5_RP)
      end do
      end do
    end do
    end do
    IntrpMat(:,:) = matmul(Vint, this%invV)
  end function QuadrilateralElement_gen_IntGaussLegendreIntrpMat

end module scale_element_quadrilateral