!> module FElib/ Element / Base
!!
!! @par Description
!!           A base module for finite element
!!
!! @author Team SCALE
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

  type, public :: ElementBase
    integer :: Np
    integer :: Nfaces
    integer :: NfpTot
    integer :: Nv
    logical, private :: LumpedMatFlag

    real(RP), allocatable :: V(:,:)
    real(RP), allocatable :: invV(:,:)
    real(RP), allocatable :: M(:,:)
    real(RP), allocatable :: invM(:,:)
    real(RP), allocatable :: Lift(:,:)
    real(RP), allocatable :: Filter(:,:)
    real(RP), allocatable :: IntWeight_lgl(:)
  contains
    procedure :: IsLumpedMatrix => ElementBase_isLumpedMatrix
  end type ElementBase
  
  !- 1D

  type, public, extends(ElementBase) :: ElementBase1D
    integer :: PolyOrder
    integer :: Nfp
    integer, allocatable :: Fmask(:,:)

    real(RP), allocatable :: x1(:) 

    real(RP), allocatable :: Dx1(:,:)

    real(RP), allocatable :: Sx1(:,:) 
  end type ElementBase1D
  
  public :: ElementBase1D_Init
  public :: ElementBase1D_Final

  !- 2D

  type, public, extends(ElementBase) :: ElementBase2D
    integer :: PolyOrder
    integer :: Nfp
    integer, allocatable :: Fmask(:,:)

    real(RP), allocatable :: x1(:)
    real(RP), allocatable :: x2(:)    

    real(RP), allocatable :: Dx1(:,:)
    real(RP), allocatable :: Dx2(:,:)

    real(RP), allocatable :: Sx1(:,:)
    real(RP), allocatable :: Sx2(:,:) 
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

  type, public, extends(ElementBase) :: ElementBase3D    
    integer :: PolyOrder_h
    integer :: Nnode_h1D
    integer :: Nfaces_h
    integer :: Nfp_h
    integer, allocatable :: Fmask_h(:,:)

    integer :: PolyOrder_v
    integer :: Nnode_v
    integer :: Nfaces_v
    integer :: Nfp_v
    integer, allocatable :: Fmask_v(:,:)

    integer, allocatable :: Colmask(:,:)
    integer, allocatable :: Hslice(:,:)
    integer, allocatable :: IndexH2Dto3D(:)
    integer, allocatable :: IndexH2Dto3D_bnd(:)    
    integer, allocatable :: IndexZ1Dto3D(:)

    real(RP), allocatable :: x1(:)
    real(RP), allocatable :: x2(:)
    real(RP), allocatable :: x3(:)
    
    real(RP), allocatable :: Dx1(:,:)
    real(RP), allocatable :: Dx2(:,:)    
    real(RP), allocatable :: Dx3(:,:)
  
    real(RP), allocatable :: Sx1(:,:)
    real(RP), allocatable :: Sx2(:,:)    
    real(RP), allocatable :: Sx3(:,:)
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
    allocate( elem%Filter(elem%Np, elem%Np))

    elem%LumpedMatFlag = lumpedmat_flag

    return
  end subroutine ElementBase_Init

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
      deallocate( elem%Filter )
    end if

    return
  end subroutine ElementBase_Final

!OCL SERIAL
  function ElementBase_isLumpedMatrix( elem ) result(lumpedmat_flag)
    implicit none
    class(ElementBase), intent(in) :: elem
    logical :: lumpedmat_flag
    !---------------------------------------------

    lumpedmat_flag = elem%LumpedMatFlag
    return
  end function ElementBase_isLumpedMatrix

  !--------------------------------------------------------------------------------
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
