!> module FElib / Element / hexahedron
!!
!! @par Description
!!           A module for a hexahedral finite element
!!
!! @author Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_element_hexahedral

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
  type, public, extends(ElementBase3D) :: HexahedralElement
  contains
    procedure :: Init => HexhedralElement_Init
    procedure :: Final => HexhedralElement_Final
    procedure :: GenIntGaussLegendreIntrpMat => HexhedralElement_gen_IntGaussLegendreIntrpMat
  end type HexahedralElement

contains
!OCL SERIAL
  subroutine HexhedralElement_Init( &
      elem, elemOrder_h, elemOrder_v,      &
      LumpedMassMatFlag )
    implicit none
    
    class(HexahedralElement), intent(inout) :: elem
    integer, intent(in) :: elemOrder_h
    integer, intent(in) :: elemOrder_v
    logical, intent(in) :: LumpedMassMatFlag

    !-----------------------------------------------------------------------------
    
    elem%PolyOrder_h = elemOrder_h
    elem%PolyOrder_v = elemOrder_v
    elem%Nnode_h1D = elemOrder_h + 1
    elem%Nnode_v = elemOrder_v + 1

    elem%Nv = 8
    elem%Nfaces_h = 4
    elem%Nfaces_v = 2  
    elem%Nfaces = elem%Nfaces_h + elem%Nfaces_v

    elem%Nfp_h = elem%Nnode_h1D*elem%Nnode_v
    elem%Nfp_v = elem%Nnode_h1D**2
    elem%NfpTot = elem%Nfp_h*elem%Nfaces_h + elem%Nfp_v*elem%Nfaces_v

    elem%Np = elem%Nfp_v * elem%Nnode_v
    
    call ElementBase3D_Init(elem, LumpedMassMatFlag)
    call construct_Element(elem)

    return
  end subroutine HexhedralElement_Init

!OCL SERIAL  
  subroutine HexhedralElement_Final(elem)
    implicit none

    class(HexahedralElement), intent(inout) :: elem
    !-----------------------------------------------------------------------------

    call ElementBase3D_Final(elem)

    return
  end subroutine HexhedralElement_Final

!OCL SERIAL  
  subroutine construct_Element(elem)

    use scale_linalgebra, only: linalgebra_inv
    use scale_polynominal, only: &
      polynominal_genGaussLobattoPt, Polynominal_GenGaussLobattoPtIntWeight,   &
      polynominal_genLegendrePoly, Polynominal_genDLegendrePoly,               &
      polynominal_genLagrangePoly, polynominal_genDLagrangePoly_lglpt
    use scale_element_quadrilateral, only: QuadrilateralElement

    implicit none
    
    type(HexahedralElement), intent(inout) :: elem

    integer :: nodes_ijk(elem%Nnode_h1D, elem%Nnode_h1D, elem%Nnode_v)

    real(RP) :: lglPts1D_h(elem%Nnode_h1D)
    real(RP) :: lglPts1D_v(elem%Nnode_v)

    real(DP) :: intWeight_lgl1DPts_h(elem%Nnode_h1D)
    real(DP) :: intWeight_lgl1DPts_v(elem%Nnode_v)

    real(RP) :: P1D_ori_h(elem%Nnode_h1D, elem%Nnode_h1D)
    real(RP) :: P1D_ori_v(elem%Nnode_v, elem%Nnode_v)
    real(RP) :: DP1D_ori_h(elem%Nnode_h1D, elem%Nnode_h1D)
    real(RP) :: DP1D_ori_v(elem%Nnode_v, elem%Nnode_v)
    real(RP) :: DLagr1D_h(elem%Nnode_h1D, elem%Nnode_h1D)
    real(RP) :: DLagr1D_v(elem%Nnode_v, elem%Nnode_v)
    real(RP) :: V2D_h(elem%Nfp_h, elem%Nfp_h)
    real(RP) :: V2D_v(elem%Nfp_v, elem%Nfp_v)
    real(RP) :: Emat(elem%Np, elem%NfpTot)
    real(RP) :: MassEdge_h(elem%Nfp_h, elem%Nfp_h)
    real(RP) :: MassEdge_v(elem%Nfp_v, elem%Nfp_v)

    integer :: i, j, k
    integer :: p1, p2, p3
    integer :: n, l, f
    integer :: Nord
    integer :: is, ie

    integer :: f_h, f_v
    integer :: fp, fp_h1, fp_h2, fp_v

    type(QuadrilateralElement) :: elem2D
    !-----------------------------------------------------------------------------

    lglPts1D_h(:)      = polynominal_genGaussLobattoPt( elem%PolyOrder_h )
    P1D_ori_h(:,:)     = polynominal_genLegendrePoly( elem%PolyOrder_h, lglPts1D_h )
    DP1D_ori_h(:,:) = polynominal_genDLegendrePoly( elem%PolyOrder_h, lglPts1D_h, P1D_ori_h )
    DLagr1D_h(:,:) = polynominal_GenDLagrangePoly_lglpt( elem%PolyOrder_h, lglPts1D_h )

    lglPts1D_v(:)      = polynominal_genGaussLobattoPt( elem%PolyOrder_v )
    P1D_ori_v(:,:)     = polynominal_genLegendrePoly( elem%PolyOrder_v, lglPts1D_v )
    DP1D_ori_v(:,:) = polynominal_genDLegendrePoly( elem%PolyOrder_v, lglPts1D_v, P1D_ori_v )
    DLagr1D_v(:,:) = polynominal_GenDLagrangePoly_lglpt( elem%PolyOrder_v, lglPts1D_v )
    
    !* Preparation 

    do k=1, elem%Nnode_v   
    do j=1, elem%Nnode_h1D
    do i=1, elem%Nnode_h1D
        nodes_ijk(i,j,k) = i + (j-1)*elem%Nnode_h1D + (k-1)*elem%Nnode_h1D**2
    end do
    end do
    end do
    
    ! Set the mask to extract the values at faces
    
    elem%Fmask_h(:,1) = reshape(nodes_ijk(:,1,:), (/ elem%Nfp_h /))
    elem%Fmask_h(:,2) = reshape(nodes_ijk(elem%Nnode_h1D,:,:), (/ elem%Nfp_h /))
    elem%Fmask_h(:,3) = reshape(nodes_ijk(:,elem%Nnode_h1D,:), (/ elem%Nfp_h /))
    elem%Fmask_h(:,4) = reshape(nodes_ijk(1,:,:), (/ elem%Nfp_h /))

    elem%Fmask_v(:,1) = reshape(nodes_ijk(:,:,1), (/ elem%Nfp_v /))
    elem%Fmask_v(:,2) = reshape(nodes_ijk(:,:,elem%Nnode_v), (/ elem%Nfp_v /))
    
    !- ColMask

    do j=1, elem%Nnode_h1D
    do i=1, elem%Nnode_h1D
      n = i + (j-1)*elem%Nnode_h1D
      elem%Colmask(:,n) = nodes_ijk(i,j,:)
    end do
    end do
    
    != Hslice

    do k=1, elem%Nnode_v
      elem%Hslice(:,k) = reshape(nodes_ijk(:,:,k), (/ elem%Nfp_v /))
    end do

    !- IndexH2Dto3D

    do k=1, elem%Nnode_v    
    do j=1, elem%Nnode_h1D
    do i=1, elem%Nnode_h1D
      n = i + (j-1)*elem%Nnode_h1D + (k-1)*elem%Nnode_h1D**2
      elem%IndexH2Dto3D(n) = nodes_ijk(i,j,1)
    end do
    end do    
    end do

    !- IndexH2Dto3D_bnd

    call elem2D%Init( elem%PolyOrder_h, .false. )

    do f_h=1, 4
      do fp_v=1, elem%Nnode_v
      do fp_h1=1, elem%Nnode_h1D
        fp = fp_h1 + (fp_v-1)*elem%Nnode_h1D + (f_h-1)*elem%Nfp_h
        elem%IndexH2Dto3D_bnd(fp) = elem2D%Fmask(fp_h1,f_h)
      end do  
      end do
    end do
    do f_v=1, 2
      do fp_h2=1, elem%Nnode_h1D
      do fp_h1=1, elem%Nnode_h1D
        fp = fp_h1 + (fp_h2-1)*elem%Nnode_h1D    &
           + (f_v-1) * elem%Nfp_v                &
           + 4 * elem%Nnode_h1D * elem%Nnode_v
          elem%IndexH2Dto3D_bnd(fp) = fp_h1 + (fp_h2-1)*elem%Nnode_h1D
      end do  
      end do
    end do

    call elem2D%Final()

    !- IndexZ1Dto3D

    do k=1, elem%Nnode_v    
    do j=1, elem%Nnode_h1D
    do i=1, elem%Nnode_h1D
      n = i + (j-1)*elem%Nnode_h1D + (k-1)*elem%Nnode_h1D**2
      elem%IndexZ1Dto3D(n) = k
    end do
    end do    
    end do
    
    !* Set the coordinates of LGL points, and the Vandermonde and differential matricies

    elem%Dx1(:,:) = 0.0_RP
    elem%Dx2(:,:) = 0.0_RP
    elem%Dx3(:,:) = 0.0_RP

    do k=1, elem%Nnode_v
    do j=1, elem%Nnode_h1D
    do i=1, elem%Nnode_h1D
      n = i + (j-1)*elem%Nnode_h1D + (k-1)*elem%Nnode_h1D**2
      
      !* Set the coordinates of LGL points
      elem%x1(n) = lglPts1D_h(i)
      elem%x2(n) = lglPts1D_h(j)
      elem%x3(n) = lglPts1D_v(k)

      !* Set the Vandermonde and differential matricies
      do p3=1, elem%Nnode_v
      do p2=1, elem%Nnode_h1D
      do p1=1, elem%Nnode_h1D
        l = p1 + (p2-1)*elem%Nnode_h1D + (p3-1)*elem%Nnode_h1D**2
        elem%V(n,l) = (P1D_ori_h(i,p1)*P1D_ori_h(j,p2)*P1D_ori_v(k,p3))                         &
                      * sqrt((dble(p1-1) + 0.5_DP)*(dble(p2-1) + 0.5_DP)*(dble(p3-1) + 0.5_DP))
  
        if(p2==j .and. p3==k) elem%Dx1(n,l) = DLagr1D_h(p1,i)
        if(p1==i .and. p3==k) elem%Dx2(n,l) = DLagr1D_h(p2,j)
        if(p1==i .and. p2==j) elem%Dx3(n,l) = DLagr1D_v(p3,k)
      end do
      end do
      end do
    end do
    end do
    end do
    elem%invV(:,:) = linAlgebra_inv(elem%V)
    
    !* Set the weights at LGL points to integrate over element

    intWeight_lgl1DPts_h(:) = Polynominal_GenGaussLobattoPtIntWeight(elem%PolyOrder_h)
    intWeight_lgl1DPts_v(:) = Polynominal_GenGaussLobattoPtIntWeight(elem%PolyOrder_v)

    do k=1, elem%Nnode_v
    do j=1, elem%Nnode_h1D
    do i=1, elem%Nnode_h1D
        l = i + (j-1)*elem%Nnode_h1D + (k-1)*elem%Nnode_h1D**2
        elem%IntWeight_lgl(l) = &
          intWeight_lgl1DPts_h(i) * intWeight_lgl1DPts_h(j) * intWeight_lgl1DPts_v(k)
    end do
    end do
    end do

    !* Set the mass matrix

    if (elem%IsLumpedMatrix()) then
      elem%invM(:,:) = 0.0_RP
      elem%M(:,:)    = 0.0_RP
      do k=1, elem%Nnode_v
      do j=1, elem%Nnode_h1D
      do i=1, elem%Nnode_h1D
        l = i + (j-1)*elem%Nnode_h1D + (k-1)*elem%Nnode_h1D**2
        elem%M(l,l) = elem%IntWeight_lgl(l)
        elem%invM(l,l) = 1.0_DP/elem%IntWeight_lgl(l)
      end do
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

    elem%Sx3(:,:) = transpose(matmul( elem%M, elem%Dx3))
    elem%Sx3(:,:) = matmul( elem%invM, elem%Sx3 )

    !* Set the lift matrix

    do k=1, elem%Nnode_v
    do i=1, elem%Nnode_h1D
      n = i + (k-1)*elem%Nnode_h1D
      do p3=1, elem%Nnode_v
      do p1=1, elem%Nnode_h1D
        l = p1 + (p3-1)*elem%Nnode_h1D
        V2D_h(n,l) =   P1D_ori_h(i,p1)*P1D_ori_v(k,p3) &
                     * sqrt( (dble(p1-1) + 0.5_DP)*(dble(p3-1) + 0.5_DP) )
      end do
      end do
    end do
    end do
    do j=1, elem%Nnode_h1D
    do i=1, elem%Nnode_h1D
      n = i + (j-1)*elem%Nnode_h1D
      do p2=1, elem%Nnode_h1D
      do p1=1, elem%Nnode_h1D
        l = p1 + (p2-1)*elem%Nnode_h1D
        V2D_v(n,l) =   P1D_ori_h(i,p1)*P1D_ori_h(j,p2) &
                     * sqrt( (dble(p1-1) + 0.5_DP)*(dble(p2-1) + 0.5_DP) )
      end do
      end do
    end do
    end do
  
    !--

    Emat(:,:) = 0.0_RP
    do f=1, elem%Nfaces_h
      if (elem%IsLumpedMatrix()) then
        MassEdge_h(:,:) = 0.0_RP
        do k=1, elem%Nnode_v
        do i=1, elem%Nnode_h1D
          l = i + (k-1)*elem%Nnode_h1D
          MassEdge_h(l,l) = intWeight_lgl1DPts_h(i) * intWeight_lgl1DPts_v(k)
        end do
        end do
      else
        MassEdge_h(:,:) = linalgebra_inv(matmul(V2D_h, transpose(V2D_h)))
      end if

      is = (f-1)*elem%Nfp_h + 1
      ie = is + elem%Nfp_h - 1
      Emat(elem%Fmask_h(:,f), is:ie) = MassEdge_h
    end do

    do f=1, elem%Nfaces_v
      if (elem%IsLumpedMatrix()) then
        MassEdge_v(:,:) = 0.0_RP
        do j=1, elem%Nnode_h1D
        do i=1, elem%Nnode_h1D
          l = i + (j-1)*elem%Nnode_h1D
          MassEdge_v(l,l) = intWeight_lgl1DPts_h(i) * intWeight_lgl1DPts_h(j)
        end do
        end do
      else
        MassEdge_v(:,:) = linalgebra_inv(matmul(V2D_v, transpose(V2D_v)))
      end if

      is = elem%Nfaces_h*elem%Nfp_h + (f-1)*elem%Nfp_v + 1
      ie = is + elem%Nfp_v - 1
      Emat(elem%Fmask_v(:,f), is:ie) = MassEdge_v
    end do

    elem%Lift(:,:) = matmul( elem%invM, Emat )

    return
  end subroutine construct_Element

!OCL SERIAL  
  function HexhedralElement_gen_IntGaussLegendreIntrpMat( this, IntrpPolyOrder, &
    intw_intrp, x_intrp, y_intrp, z_intrp ) result(IntrpMat)

    use scale_polynominal, only: &
      Polynominal_genLegendrePoly,    &
      Polynominal_GenGaussLegendrePt, &
      Polynominal_GenGaussLegendrePtIntWeight
    
    implicit none

    class(HexahedralElement), intent(in) :: this
    integer, intent(in) :: IntrpPolyOrder
    real(RP), intent(out), optional :: intw_intrp(IntrpPolyOrder**3)
    real(RP), intent(out), optional :: x_intrp(IntrpPolyOrder**3)
    real(RP), intent(out), optional :: y_intrp(IntrpPolyOrder**3)
    real(RP), intent(out), optional :: z_intrp(IntrpPolyOrder**3)
    real(RP) :: IntrpMat(IntrpPolyOrder**3,this%Np)

    real(RP) :: r_int1D_i(IntrpPolyOrder)
    real(RP) :: r_int1Dw_i(IntrpPolyOrder)
    real(RP) :: P_int1D_ori_h(IntrpPolyOrder,this%Nnode_h1D)
    real(RP) :: P_int1D_ori_v(IntrpPolyOrder,this%Nfp_v)
    real(RP) :: Vint(IntrpPolyOrder**3,this%Np)

    integer :: p1, p2, p3, p1_, p2_, p3_
    integer :: n_, l_, m_
    !-----------------------------------------------------

    r_int1D_i(:) = Polynominal_GenGaussLegendrePt( IntrpPolyOrder )
    r_int1Dw_i(:) = Polynominal_GenGaussLegendrePtIntWeight( IntrpPolyOrder )
    P_int1D_ori_h(:,:) = Polynominal_GenLegendrePoly( this%PolyOrder_h, r_int1D_i)
    P_int1D_ori_v(:,:) = Polynominal_GenLegendrePoly( this%PolyOrder_v, r_int1D_i)

    do p3_=1, IntrpPolyOrder
    do p2_=1, IntrpPolyOrder
    do p1_=1, IntrpPolyOrder
      n_= p1_ + (p2_-1)*IntrpPolyOrder + (p3_-1)*IntrpPolyOrder**2
      if (present(intw_intrp)) intw_intrp(n_) = r_int1Dw_i(p1_) * r_int1Dw_i(p2_) * r_int1Dw_i(p3_)
      if (present(x_intrp)) x_intrp(n_) = r_int1D_i(p1_)
      if (present(y_intrp)) y_intrp(n_) = r_int1D_i(p2_)
      if (present(z_intrp)) z_intrp(n_) = r_int1D_i(p3_)

      do p3=1, this%Nnode_v
      do p2=1, this%Nnode_h1D
      do p1=1, this%Nnode_h1D
        l_ = p1 + (p2-1)*this%Nnode_h1D + (p3-1)*this%Nnode_h1D**2
        Vint(n_,l_) =  P_int1D_ori_h(p1_,p1) * sqrt(dble(p1-1) + 0.5_DP) &
                     * P_int1D_ori_h(p2_,p2) * sqrt(dble(p2-1) + 0.5_DP) &
                     * P_int1D_ori_v(p3_,p3) * sqrt(dble(p3-1) + 0.5_DP)
      end do
      end do
      end do
    end do
    end do
    end do
    IntrpMat(:,:) = matmul(Vint, this%invV)

    return
  end function HexhedralElement_gen_IntGaussLegendreIntrpMat

end module scale_element_hexahedral