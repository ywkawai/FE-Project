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
  type, public :: elementbase
    integer :: Np
    integer :: Nfaces
    integer :: NfpTot
    integer :: Nv

    real(RP), allocatable :: V(:,:)
    real(RP), allocatable :: invV(:,:)
    real(RP), allocatable :: M(:,:)
    real(RP), allocatable :: invM(:,:)
    real(RP), allocatable :: Lift(:,:)
    real(RP), allocatable :: IntWeight_lgl(:)    
  end type elementbase
  
  !-

  type, public, extends(elementbase) :: elementbase1D
    integer :: PolyOrder
    integer :: Nfp
    integer, allocatable :: Fmask(:,:)

    real(RP), allocatable :: x1(:) 

    real(RP), allocatable :: Dx1(:,:)

    real(RP), allocatable :: Sx1(:,:) 
  end type elementbase1D
  
  public :: ElementBase1D_Init
  public :: ElementBase1D_Final

  !-

  type, public, extends(elementbase) :: elementbase2D
    integer :: PolyOrder
    integer :: Nfp
    integer, allocatable :: Fmask(:,:)

    real(RP), allocatable :: x1(:)
    real(RP), allocatable :: x2(:)    

    real(RP), allocatable :: Dx1(:,:)
    real(RP), allocatable :: Dx2(:,:)

    real(RP), allocatable :: Sx1(:,:)
    real(RP), allocatable :: Sx2(:,:) 
  end type elementbase2D

  public :: ElementBase2D_Init
  public :: ElementBase2D_Final

  !-

  type, public, extends(elementbase) :: elementbase3D    
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
  end type elementbase3D

  public :: ElementBase3D_Init
  public :: ElementBase3D_Final

  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  ! 

contains
  subroutine ElementBase1D_Init( elem )
    class(ElementBase1D), intent(inout) :: elem

    !-----------------------------------------------------------------------------

    allocate( elem%x1(elem%Np) )
    allocate( elem%M(elem%Np, elem%Np) )
    allocate( elem%invM(elem%Np, elem%Np) )
    allocate( elem%V(elem%Np, elem%Np) )
    allocate( elem%invV(elem%Np, elem%Np) )
    allocate( elem%Dx1(elem%Np, elem%Np) )
    allocate( elem%Sx1(elem%Np, elem%Np) )
    allocate( elem%IntWeight_lgl(elem%Np) )    
    allocate( elem%Lift(elem%Np, elem%NfpTot) )    
    allocate( elem%Fmask(elem%Nfp, elem%Nfaces) )
    
  end subroutine ElementBase1D_Init

  subroutine ElementBase1D_Final( elem )
    class(ElementBase1D), intent(inout) :: elem

    !-----------------------------------------------------------------------------

    deallocate( elem%x1 )
    deallocate( elem%M )
    deallocate( elem%invM )
    deallocate( elem%V )
    deallocate( elem%invV )
    deallocate( elem%Dx1 )
    deallocate( elem%Sx1 )
    deallocate( elem%Lift )
    deallocate( elem%IntWeight_lgl )
    deallocate( elem%Fmask )

  end subroutine ElementBase1D_Final

  subroutine ElementBase2D_Init( elem )
    class(ElementBase2D), intent(inout) :: elem

    !-----------------------------------------------------------------------------

    allocate( elem%x1(elem%Np), elem%x2(elem%Np) )
    allocate( elem%M(elem%Np, elem%Np) )
    allocate( elem%invM(elem%Np, elem%Np) )
    allocate( elem%V(elem%Np, elem%Np) )
    allocate( elem%invV(elem%Np, elem%Np) )
    allocate( elem%Dx1(elem%Np, elem%Np), elem%Dx2(elem%Np, elem%Np) )
    allocate( elem%Sx1(elem%Np, elem%Np), elem%Sx2(elem%Np, elem%Np) )
    allocate( elem%IntWeight_lgl(elem%Np) )    
    allocate( elem%Lift(elem%Np, elem%NfpTot) )    
    allocate( elem%Fmask(elem%Nfp, elem%Nfaces) )

  end subroutine ElementBase2D_Init

  subroutine ElementBase2D_Final( elem )
    class(ElementBase2D), intent(inout) :: elem

    !-----------------------------------------------------------------------------

    deallocate( elem%x1, elem%x2 )
    deallocate( elem%M )
    deallocate( elem%invM )
    deallocate( elem%V )
    deallocate( elem%invV )
    deallocate( elem%Dx1, elem%Dx2 )
    deallocate( elem%Sx1, elem%Sx2 )
    deallocate( elem%Lift )
    deallocate( elem%IntWeight_lgl )
    deallocate( elem%Fmask )

  end subroutine ElementBase2D_Final

  subroutine ElementBase3D_Init( elem )
    class(ElementBase3D), intent(inout) :: elem

    !-----------------------------------------------------------------------------

    allocate( elem%x1(elem%Np), elem%x2(elem%Np), elem%x3(elem%Np) )
    allocate( elem%M(elem%Np, elem%Np) )
    allocate( elem%invM(elem%Np, elem%Np) )
    allocate( elem%V(elem%Np, elem%Np) )
    allocate( elem%invV(elem%Np, elem%Np) )
    allocate( elem%Dx1(elem%Np, elem%Np), elem%Dx2(elem%Np, elem%Np), elem%Dx3(elem%Np, elem%Np) )
    allocate( elem%Sx1(elem%Np, elem%Np), elem%Sx2(elem%Np, elem%Np), elem%Sx3(elem%Np, elem%Np) )
    allocate( elem%IntWeight_lgl(elem%Np) )    
    allocate( elem%Lift(elem%Np, elem%NfpTot) )    
    allocate( elem%Fmask_h(elem%Nfp_h, elem%Nfaces_h), elem%Fmask_v(elem%Nfp_v, elem%Nfaces_v) )
    allocate( elem%Colmask(elem%Nnode_v,elem%Nfp_h))
    allocate( elem%Hslice(elem%Nfp_v,elem%Nnode_v) )
    allocate( elem%IndexH2Dto3D(elem%Np) )
    allocate( elem%IndexZ1Dto3D(elem%Np) )
        
  end subroutine ElementBase3D_Init

  subroutine ElementBase3D_Final( elem )
    class(ElementBase3D), intent(inout) :: elem

    !-----------------------------------------------------------------------------

    deallocate( elem%x1, elem%x2, elem%x3 )
    deallocate( elem%M )
    deallocate( elem%invM )
    deallocate( elem%V )
    deallocate( elem%invV )
    deallocate( elem%Dx1, elem%Dx2, elem%Dx3 )
    deallocate( elem%Sx1, elem%Sx2, elem%Sx3 )
    deallocate( elem%Lift )
    deallocate( elem%IntWeight_lgl )
    deallocate( elem%Fmask_h, elem%Fmask_v )
    deallocate( elem%Colmask, elem%Hslice )
    deallocate( elem%IndexH2Dto3D, elem%IndexZ1Dto3D )

  end subroutine ElementBase3D_Final

end module scale_element_base
