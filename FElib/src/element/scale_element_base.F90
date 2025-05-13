!> module FElib / Element / Base
!!
!! @par Description
!!           A base module for finite element
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_element_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision  
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !- Base

  !> Derived type representing an arbitrary finite element
  type, public :: ElementBase
    integer :: Np                       !< Number of nodes within an element
    integer :: Nfaces                   !< Number of faces
    integer :: NfpTot                   !< Total number of nodes on faces
    integer :: Nv                       !< Number of vertices with an element
    logical, private :: LumpedMatFlag   !< Flag whether the lumped mass matrix is used

    real(RP), allocatable :: V(:,:)             !< The Vandermonde matrix (V) whose size is Np x Np 
    real(RP), allocatable :: invV(:,:)          !< Inversion of the Vandermonde matrix (V^-1) whose size is Np x Np
    real(RP), allocatable :: M(:,:)             !< Mass matrix (M) whose size is Np x Np
    real(RP), allocatable :: invM(:,:)          !< Inversion of the mass matrix (M^-1) whose size Np x NP
    real(RP), allocatable :: Lift(:,:)          !< Lifting matrix with element boundary integrals whose size Np x NfpTot
    real(RP), allocatable :: IntWeight_lgl(:)   !< Weights of gaussian quadrature with the LGL nodes
  contains
    procedure :: IsLumpedMatrix => ElementBase_isLumpedMatrix
  end type ElementBase
  public :: ElementBase_construct_MassMat
  public :: ElementBase_construct_StiffMat
  public :: ElementBase_construct_LiftMat

  !- 1D

  !> Derived type representing a 1D reference element
  type, public, extends(ElementBase) :: ElementBase1D
    integer :: PolyOrder                 !< Polynomial order
    integer :: Nfp                       !< Number of nodes on an element face
    integer, allocatable :: Fmask(:,:)   !< Array saving indices to extract nodal values on the faces

    real(RP), allocatable :: x1(:)       !< Array saving x1-coordinate of nodes within the reference element

    real(RP), allocatable :: Dx1(:,:)    !< Elementwise differential matrix for the x1-coordinate direction (Dx1 = M^-1 Sx1)

    real(RP), allocatable :: Sx1(:,:)    !< Elementwise stiffness matrix for the x1-coordinate direction
  end type ElementBase1D
  
  public :: ElementBase1D_Init
  public :: ElementBase1D_Final

  !- 2D

  !> Derived type representing a 2D reference element
  type, public, extends(ElementBase) :: ElementBase2D
    integer :: PolyOrder               !< Polynomial order
    integer :: Nfp                     !< Number of nodes on an element face
    integer, allocatable :: Fmask(:,:) !< Array saving indices to extract nodal values on the faces

    real(RP), allocatable :: x1(:)     !< Array saving x1-coordinate of nodes within the reference element
    real(RP), allocatable :: x2(:)     !< Array saving x2-coordinate of nodes within the reference element

    real(RP), allocatable :: Dx1(:,:)  !< Elementwise differential matrix for the x1-coordinate direction (Dx1 = M^-1 Sx1)
    real(RP), allocatable :: Dx2(:,:)  !< Elementwise differential matrix for the x2-coordinate direction (Dx2 = M^-1 Sx2)

    real(RP), allocatable :: Sx1(:,:)  !< Elementwise stiffness matrix for the x1-coordinate direction
    real(RP), allocatable :: Sx2(:,:)  !< Elementwise stiffness matrix for the x2-coordinate direction
  end type ElementBase2D

  public :: ElementBase2D_Init
  public :: ElementBase2D_Final

  interface
    function ElementBase2D_GenIntGaussLegendreIntrpMat(this, IntrpPolyOrder, &
      intw_intrp, x_intrp, y_intrp ) result(IntrpMat)
  
      import ElementBase2D
      import RP
      class(ElementBase2D), intent(in) :: this
      integer, intent(in) :: IntrpPolyOrder
      real(RP), intent(out), optional :: intw_intrp(IntrpPolyOrder**2)
      real(RP), intent(out), optional :: x_intrp(IntrpPolyOrder**2)
      real(RP), intent(out), optional :: y_intrp(IntrpPolyOrder**2)
      real(RP) :: IntrpMat(IntrpPolyOrder**2,this%Np)
    end function ElementBase2D_GenIntGaussLegendreIntrpMat
  end interface

  !- 3D

  !> Derived type representing a 3D reference element
  type, public, extends(ElementBase) :: ElementBase3D    
    integer :: PolyOrder_h               !< Polynomial order with the horizontal direction
    integer :: Nnode_h1D                 !< Number of nodes along the horizontal coordinate
    integer :: Nfaces_h                  !< Number of nodes on an horizontal face of the reference element
    integer :: Nfp_h                     !< Number of horizontal faces of the reference element
    integer, allocatable :: Fmask_h(:,:) !< Array saving indices to extract nodal values on the horizontal faces

    integer :: PolyOrder_v               !< Polynomial order with the vertical direction
    integer :: Nnode_v                   !< Number of nodes along the vertical coordinate
    integer :: Nfaces_v                  !< Number of nodes on an vertical face of the reference element
    integer :: Nfp_v                     !< Number of vertical faces of the reference element
    integer, allocatable :: Fmask_v(:,:) !< Number of vertical faces of the reference element

    integer, allocatable :: Colmask(:,:)        !< Array saving indices to extract nodal values on the vertical columns
    integer, allocatable :: Hslice(:,:)         !< Array saving indices to extract nodal values on the horizontal plane
    integer, allocatable :: IndexH2Dto3D(:)     !< Array saving indices to expand 2D horizontal nodal values into the 3D nodal values
    integer, allocatable :: IndexH2Dto3D_bnd(:) !< Array saving indices to expand 2D horizontal nodal values into the 3D nodal values on element faces
    integer, allocatable :: IndexZ1Dto3D(:)     !< Array saving indices to expand 1D vertical nodal values into the 3D nodal values

    real(RP), allocatable :: x1(:) !< Array saving x1-coordinate of nodes within the reference element
    real(RP), allocatable :: x2(:) !< Array saving x2-coordinate of nodes within the reference element
    real(RP), allocatable :: x3(:) !< Array saving x3-coordinate of nodes within the reference element
    
    real(RP), allocatable :: Dx1(:,:) !< Elementwise differential matrix for the x1-coordinate direction (Dx1 = M^-1 Sx1)
    real(RP), allocatable :: Dx2(:,:) !< Elementwise differential matrix for the x2-coordinate direction (Dx2 = M^-1 Sx2)
    real(RP), allocatable :: Dx3(:,:) !< Elementwise differential matrix for the x3-coordinate direction (Dx3 = M^-1 Sx3)
  
    real(RP), allocatable :: Sx1(:,:) !< Elementwise stiffness matrix for the x1-coordinate direction
    real(RP), allocatable :: Sx2(:,:) !< Elementwise stiffness matrix for the x2-coordinate direction    
    real(RP), allocatable :: Sx3(:,:) !< Elementwise stiffness matrix for the x3-coordinate direction
  end type ElementBase3D
  
  public :: ElementBase3D_Init
  public :: ElementBase3D_Final

  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  ! 

  private :: ElementBase_Init
  private :: ElementBase_Final

contains
  !-- Base Element ------------------------------------------------------------------------------

!> Initialize a base object to manage a reference element
!OCL SERIAL
  subroutine ElementBase_Init( elem, lumpedmat_flag )
    implicit none

    class(ElementBase), intent(inout) :: elem
    logical, intent(in) :: lumpedmat_flag
    !-----------------------------------------------------------------------------

    allocate( elem%M(elem%Np, elem%Np) )
    allocate( elem%invM(elem%Np, elem%Np) )
    allocate( elem%V(elem%Np, elem%Np) )
    allocate( elem%invV(elem%Np, elem%Np) )
    allocate( elem%Lift(elem%Np, elem%NfpTot) )    
    
    allocate( elem%IntWeight_lgl(elem%Np) )  

    elem%LumpedMatFlag = lumpedmat_flag

    return
  end subroutine ElementBase_Init

!> Finalize a base object to manage a reference element
!OCL SERIAL
  subroutine ElementBase_Final( elem )
    implicit none

    class(ElementBase), intent(inout) :: elem
    !-----------------------------------------------------------------------------
    
    if ( allocated(elem%M) ) then
      deallocate( elem%M )
      deallocate( elem%invM )
      deallocate( elem%V )
      deallocate( elem%invV )
      deallocate( elem%Lift )
      deallocate( elem%IntWeight_lgl )
    end if

    return
  end subroutine ElementBase_Final

!> Get a flag whether the lumped mass matrix is used
!OCL SERIAL
  function ElementBase_isLumpedMatrix( elem ) result(lumpedmat_flag)
    implicit none
    class(ElementBase), intent(in) :: elem
    logical :: lumpedmat_flag
    !---------------------------------------------

    lumpedmat_flag = elem%LumpedMatFlag
    return
  end function ElementBase_isLumpedMatrix


!> Construct mass matrix
!!  M^-1 = V V^T
!!  M = ( M^-1 )^-1
!OCL SERIAL
  subroutine ElementBase_construct_MassMat( V, Np, &
    MassMat, invMassMat )
    use scale_linalgebra, only: linalgebra_inv
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(in) :: V(Np,Np)
    real(RP), intent(out) :: MassMat(Np,Np)
    real(RP), intent(out), optional :: invMassMat(Np,Np)

    real(RP) :: tmpMat(Np,Np)
    real(RP) :: invM(Np,Np)
    !------------------------------------

    tmpMat(:,:) = transpose(V)
    invM(:,:) = matmul( V, tmpMat )
    MassMat(:,:) = linAlgebra_inv( invM )

    if ( present(invMassMat) ) invMassMat(:,:) = invM(:,:)
    return
  end subroutine ElementBase_construct_MassMat

!> Construct stiffness matrix
!!  StiffMat_i = M^-1 ( M D_xi )^T
!OCL SERIAL
  subroutine ElementBase_construct_StiffMat( MassMat, invMassMat, DMat, Np, &
      StiffMat )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(in) :: MassMat(Np,Np)
    real(RP), intent(in) :: invMassMat(Np,Np)
    real(RP), intent(in) :: DMat(Np,Np)
    real(RP), intent(out) :: StiffMat(Np,Np)

    real(RP) :: tmpMat1(Np,Np)
    real(RP) :: tmpMat2(Np,Np)
    !------------------------------------

    tmpMat1(:,:) = matmul( MassMat, DMat )
    tmpMat2(:,:) = transpose( tmpMat1 )
    StiffMat(:,:) = matmul( invMassMat, tmpMat2 )

    return
  end subroutine ElementBase_construct_StiffMat

!> Construct stiffness matrix
!!  StiffMat_i = M^-1 ( M D_xi )^T
!OCL SERIAL
  subroutine ElementBase_construct_LiftMat( invM, EMat, Np, NfpTot, &
      LiftMat )
    implicit none
    integer, intent(in) :: Np
    integer, intent(in) :: NfpTot
    real(RP), intent(in) :: invM(Np,Np)
    real(RP), intent(in) :: EMat(Np,NfpTot)
    real(RP), intent(out) :: LiftMat(Np,NfpTot)
    !------------------------------------

    LiftMat(:,:) = matmul( invM, EMat )
    return
  end subroutine ElementBase_construct_LiftMat 

  !-- 1D Element ------------------------------------------------------------------------------

!> Initialize an object to manage a 1D reference element
!OCL SERIAL
  subroutine ElementBase1D_Init( elem, lumpedmat_flag )
    implicit none

    class(ElementBase1D), intent(inout) :: elem
    logical, intent(in) :: lumpedmat_flag    
    !-----------------------------------------------------------------------------

    call ElementBase_Init( elem, lumpedmat_flag )

    allocate( elem%x1(elem%Np) )  
    allocate( elem%Fmask(elem%Nfp, elem%Nfaces) )
    
    allocate( elem%Dx1(elem%Np, elem%Np) )
    allocate( elem%Sx1(elem%Np, elem%Np) )

    return
  end subroutine ElementBase1D_Init

!> Finalize an object to manage a 1D reference element
!OCL SERIAL
  subroutine ElementBase1D_Final( elem )
    implicit none

    class(ElementBase1D), intent(inout) :: elem
    !-----------------------------------------------------------------------------

    if ( allocated( elem%x1 ) )  then
      deallocate( elem%x1 )
      deallocate( elem%Dx1 )
      deallocate( elem%Sx1 )
      deallocate( elem%Fmask )
    end if

    call ElementBase_Final( elem )

    return
  end subroutine ElementBase1D_Final

  !-- 2D Element ------------------------------------------------------------------------------

!> Initialize an object to manage a 2D reference element
!OCL SERIAL
  subroutine ElementBase2D_Init( elem, lumpedmat_flag )
    implicit none

    class(ElementBase2D), intent(inout) :: elem
    logical, intent(in) :: lumpedmat_flag    
    !-----------------------------------------------------------------------------

    call ElementBase_Init( elem, lumpedmat_flag )

    allocate( elem%x1(elem%Np), elem%x2(elem%Np) )
    allocate( elem%Fmask(elem%Nfp, elem%Nfaces) )

    allocate( elem%Dx1(elem%Np, elem%Np), elem%Dx2(elem%Np, elem%Np) )
    allocate( elem%Sx1(elem%Np, elem%Np), elem%Sx2(elem%Np, elem%Np) )

    return
  end subroutine ElementBase2D_Init

!> Finalize an object to manage a 2D reference element
!OCL SERIAL
  subroutine ElementBase2D_Final( elem )
    implicit none

    class(ElementBase2D), intent(inout) :: elem
    !-----------------------------------------------------------------------------

    if ( allocated( elem%x1 ) )  then
      deallocate( elem%x1, elem%x2 )
      deallocate( elem%Fmask )

      deallocate( elem%Dx1, elem%Dx2 )
      deallocate( elem%Sx1, elem%Sx2 )
    end if

    call ElementBase_Final( elem )

    return
  end subroutine ElementBase2D_Final

  !-- 3D Element ------------------------------------------------------------------------------
  
!> Initialize an object to manage a 3D reference element
!OCL SERIAL
  subroutine ElementBase3D_Init( elem, lumpedmat_flag )
    implicit none

    class(ElementBase3D), intent(inout) :: elem
    logical, intent(in) :: lumpedmat_flag
    !-----------------------------------------------------------------------------

    call ElementBase_Init( elem, lumpedmat_flag )

    allocate( elem%x1(elem%Np), elem%x2(elem%Np), elem%x3(elem%Np) )
    allocate( elem%Dx1(elem%Np, elem%Np), elem%Dx2(elem%Np, elem%Np), elem%Dx3(elem%Np, elem%Np) )
    allocate( elem%Sx1(elem%Np, elem%Np), elem%Sx2(elem%Np, elem%Np), elem%Sx3(elem%Np, elem%Np) )
    allocate( elem%Fmask_h(elem%Nfp_h, elem%Nfaces_h), elem%Fmask_v(elem%Nfp_v, elem%Nfaces_v) )
    allocate( elem%Colmask(elem%Nnode_v,elem%Nfp_v))
    allocate( elem%Hslice(elem%Nfp_v,elem%Nnode_v) )
    allocate( elem%IndexH2Dto3D(elem%Np) )
    allocate( elem%IndexH2Dto3D_bnd(elem%NfpTot) )
    allocate( elem%IndexZ1Dto3D(elem%Np) )
  
    return
  end subroutine ElementBase3D_Init

!> Finalize an object to manage a 3D reference element
!OCL SERIAL
  subroutine ElementBase3D_Final( elem )
    implicit none

    class(ElementBase3D), intent(inout) :: elem
    !-----------------------------------------------------------------------------

    if ( allocated( elem%x1 ) )  then
      deallocate( elem%x1, elem%x2, elem%x3 )
      deallocate( elem%Dx1, elem%Dx2, elem%Dx3 )
      deallocate( elem%Sx1, elem%Sx2, elem%Sx3 )
      deallocate( elem%Fmask_h, elem%Fmask_v )
      deallocate( elem%Colmask, elem%Hslice )
      deallocate( elem%IndexH2Dto3D, elem%IndexH2Dto3D_bnd )
      deallocate( elem%IndexZ1Dto3D )
    end if

    call ElementBase_Final( elem )    

    return
  end subroutine ElementBase3D_Final

end module scale_element_base
