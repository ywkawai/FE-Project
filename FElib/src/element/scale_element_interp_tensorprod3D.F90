!-------------------------------------------------------------------------------
!> module FElib / Element / Interpolation with 3D tensor product elements
!!
!! @par Description
!!           A module for providing interpolation operations assuming a 3D tensor product element with (p+1)^3 DOF
!!
!! @author Yuta Kawai, Xuanzhengbo Ren, and Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_element_interp_tensorprod3d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort 

  use scale_element_base, only: ElementBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public :: ElementInterpTensorProd3D
    type(ElementBase3D), pointer :: elem  
    integer :: PolyOrder_h_intrp
    integer :: NnodeH1D_GL
    integer :: Np3D_hGLvGLL
    integer :: Np3D
    real(RP), allocatable :: IntrpMat1D_gll2gl(:,:)
    real(RP), allocatable :: IntrpMat1D_gll2gl_tr(:,:)
    real(RP), allocatable :: IntW1D_gl(:)

    real(RP),  allocatable :: V1D_GL(:,:)
    real(RP),  allocatable :: inv_V1D_GL(:,:)
    real(RP), allocatable :: IntrpMat1D_gl2gll(:,:)
    real(RP), allocatable :: IntrpMat1D_gl2gll_tr(:,:)

  contains
    procedure :: Init => ElementInterpTensorProd3D_Init
    procedure :: Final => ElementInterpTensorProd3D_Final
    procedure :: do_GLL2GL => ElementInterpTensorProd3D_gll2gl
    procedure :: do_GL2GLL => ElementInterpTensorProd3D_gl2gll
  end type ElementInterpTensorProd3D

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
  private :: intrp_core
  
contains
!OCL SERIAL
  subroutine ElementInterpTensorProd3D_Init( this, elem, PolyOrder_h_ip )
    use scale_linalgebra, only: linalgebra_inv
    use scale_polynominal, only: &
      Polynominal_genLegendrePoly,    &
      Polynominal_GenGaussLegendrePt, &
      Polynominal_GenGaussLegendrePtIntWeight
    use scale_element_line, only: LineElement
    use scale_prc
    implicit none
    class(ElementInterpTensorProd3D), intent(inout) :: this
    type(ElementBase3D), intent(in), target :: elem
    integer, intent(in) :: Polyorder_h_ip

    real(RP) :: r_int1D_i(Polyorder_h_ip+1)
    real(RP) :: r_int1Dw_i(Polyorder_h_ip+1)
    real(RP) :: P_int1D_ori_h(Polyorder_h_ip+1,elem%PolyOrder_h+1)
    real(RP) :: Vint1D(Polyorder_h_ip+1,elem%PolyOrder_h+1)

    real(RP) :: P1D_gl_h(Polyorder_h_ip+1,Polyorder_h_ip+1)
    real(RP) :: P1D_gll_h(elem%Nnode_h1D,Polyorder_h_ip+1)
    real(RP) :: Vint1D_gl(elem%Nnode_h1D,Polyorder_h_ip+1)

    integer :: p1, p1_
    integer :: n_, l_

    type(LineElement) :: elem1D
    !----------------------------------------

    this%elem => elem
    this%PolyOrder_h_intrp = Polyorder_h_ip
    this%NnodeH1D_GL = Polyorder_h_ip + 1
    this%Np3D_hGLvGLL = this%NnodeH1D_GL**2 * elem%Nnode_v
    this%Np3D = elem%Np

    allocate( this%IntrpMat1D_gll2gl(this%NnodeH1D_GL,elem%Nnode_h1D) )
    allocate( this%IntrpMat1D_gll2gl_tr(elem%Nnode_h1D,this%NnodeH1D_GL) )
    allocate( this%IntW1D_gl(this%NnodeH1D_GL) )

    !-
    call elem1D%Init( elem%PolyOrder_h, .false. )
    r_int1D_i(:) = Polynominal_GenGaussLegendrePt( this%NnodeH1D_GL )
    r_int1Dw_i(:) = Polynominal_GenGaussLegendrePtIntWeight(this%NnodeH1D_GL )
    P_int1D_ori_h(:,:) = Polynominal_GenLegendrePoly( elem%PolyOrder_h, r_int1D_i )

    do p1_=1,this%NnodeH1D_GL
      n_ = p1_
      this%IntW1D_gl(n_) = r_int1D_i(p1_)
      do p1=1, elem%Nnode_h1D
        l_ = p1
        Vint1D(n_,l_) = P_int1D_ori_h(p1_,p1) * sqrt( real(p1-1,kind=RP) + 0.5_RP )
      end do
    end do
    this%IntrpMat1D_gll2gl(:,:) = matmul( Vint1D, elem1D%invV )
    this%IntrpMat1D_gll2gl_tr(:,:) = transpose(this%IntrpMat1D_gll2gl)

    !-
    allocate( this%V1D_GL(this%NnodeH1D_GL,Polyorder_h_ip+1) )
    allocate( this%inv_V1D_GL(Polyorder_h_ip+1,this%NnodeH1D_GL) )
    allocate( this%IntrpMat1D_gl2gll(elem%Nnode_h1D,this%NnodeH1D_GL) )
    allocate( this%IntrpMat1D_gl2gll_tr(this%NnodeH1D_GL,elem%Nnode_h1D) )

    P1D_gl_h(:,:) = Polynominal_GenLegendrePoly( Polyorder_h_ip, r_int1D_i )
    P1D_gll_h(:,:) = Polynominal_GenLegendrePoly( Polyorder_h_ip, elem1D%x1 )

    do p1_=1, this%NnodeH1D_GL
      n_ = p1_
      do p1=1, Polyorder_h_ip+1
        l_ = p1
        this%V1D_GL(n_,l_) = P1D_gl_h(p1_,p1) * sqrt( real(p1-1,kind=RP) + 0.5_RP )
      end do
    end do
    this%inv_V1D_GL(:,:) = linalgebra_inv( this%V1D_GL(:,:) )

    Vint1D_gl(:,:) = 0.0_RP
    do p1_=1,this%elem%Nnode_h1D
      n_ = p1_
      this%IntW1D_gl(n_) = r_int1D_i(p1)
      do p1=1, elem%PolyOrder_h+1
        l_ = p1
        Vint1D_gl(n_,l_) = P1D_gll_h(p1_,p1) * sqrt( real(p1-1,kind=RP) + 0.5_RP )
      end do
    end do
    this%IntrpMat1D_gl2gll(:,:) = matmul( Vint1D_gl, this%inv_V1D_GL )
    this%IntrpMat1D_gl2gll_tr(:,:) = transpose(this%IntrpMat1D_gl2gll(:,:))

    !-
    call elem1D%Final()
    return
  end subroutine ElementInterpTensorProd3D_Init

!OCL SERIAL
  subroutine ElementInterpTensorProd3D_Final( this )
    implicit none
    class(ElementInterpTensorProd3D), intent(inout) :: this
    !----------------------------------------
    deallocate( this%IntrpMat1D_gll2gl, this%IntrpMat1D_gll2gl_tr )
    deallocate( this%IntW1D_gl )

    deallocate( this%V1D_GL, this%inv_V1D_GL )
    deallocate( this%IntrpMat1D_gl2gll, this%IntrpMat1D_gl2gll_tr )
    return
  end subroutine ElementInterpTensorProd3D_Final

!OCL SERIAL
  subroutine ElementInterpTensorProd3D_gll2gl( this, q, q_intrp )
    implicit none
    class(ElementInterpTensorProd3D), intent(in) :: this
    real(RP), intent(in) :: q(this%Np3D)
    real(RP), intent(out) :: q_intrp(this%Np3D_hGLvGLL)
    !----------------------------------------

    call intrp_core( q_intrp, &
      q, this%elem%Nnode_h1D, this%NnodeH1D_GL, this%elem%Nnode_v, &
      this%IntrpMat1D_gll2gl_tr )
    return
  end subroutine ElementInterpTensorProd3D_gll2gl

!OCL SERIAL
  subroutine ElementInterpTensorProd3D_gl2gll( this, q, q_intrp )
    implicit none
    class(ElementInterpTensorProd3D), intent(in) :: this
    real(RP), intent(in) :: q(this%Np3D_hGLvGLL)
    real(RP), intent(out) :: q_intrp(this%Np3D)
    !----------------------------------------

    call intrp_core( q_intrp, &
      q, this%NnodeH1D_GL, this%elem%Nnode_h1D, this%elem%Nnode_v, &
      this%IntrpMat1D_gl2gll_tr )
    return
  end subroutine ElementInterpTensorProd3D_gl2gll

!---------
!OCL SERIAL
  subroutine intrp_core( q_intrp, &
    q, Nnode_h1D, Nnode_h1D_ip, Nnode_v, IntrpMat1D_tr )
    integer, intent(in) :: Nnode_h1D
    integer, intent(in) :: Nnode_h1D_ip
    integer, intent(in) :: Nnode_v
    real(RP), intent(out) :: q_intrp(Nnode_h1D_ip,Nnode_h1D_ip,Nnode_v)
    real(RP), intent(in) :: q(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: IntrpMat1D_tr(Nnode_h1D,Nnode_h1D_ip)

    integer :: p1, p2, p3
    integer :: pp
    real(RP) :: tmp
    real(RP) :: q_tmp(Nnode_h1D_ip,Nnode_h1D)
    !--------------------------------

    do p3=1, Nnode_v
      do p2=1, Nnode_h1D
      do p1=1, Nnode_h1D_ip
        tmp = 0.0_RP
        do pp=1, Nnode_h1D
          tmp = tmp + IntrpMat1D_tr(pp,p1) * q(pp,p2,p3)
        end do
        q_tmp(p1,p2)= tmp
      end do
      end do

      q_intrp(:,:,p3) = 0.0_RP
      do p2=1, Nnode_h1D_ip
      do pp=1, Nnode_h1D
        do p1=1, Nnode_h1D_ip
          q_intrp(p1,p2,p3) = q_intrp(p1,p2,p3) + IntrpMat1D_tr(pp,p2) * q_tmp(p1,pp)
        end do
      end do
      end do
    end do
    return
  end subroutine intrp_core

end module scale_element_interp_tensorprod3d
