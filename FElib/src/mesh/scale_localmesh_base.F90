#include "scaleFElib.h"
module scale_localmesh_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: elementbase

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  
  type, public :: LocalMeshBase
    integer :: Ne
    integer :: NeS
    integer :: NeE
    integer :: NeA
    integer :: Nv

    class(ElementBase), pointer :: refElem
    
    real(RP), allocatable :: pos_ev(:,:)
    real(RP), allocatable :: pos_en(:,:,:)
    real(RP), allocatable :: normal_fn(:,:,:)

    real(RP), allocatable :: sJ(:,:)
    real(RP), allocatable :: J(:,:)
    
    real(RP), allocatable :: Escale(:,:,:,:)
    real(RP), allocatable :: Fscale(:,:)

    integer, allocatable :: EToV(:,:)
    integer, allocatable :: EToE(:,:)
    integer, allocatable :: EToF(:,:)
    integer, allocatable :: VMapM(:,:)
    integer, allocatable :: VMapP(:,:)
    integer, allocatable :: MapM(:,:)
    integer, allocatable :: MapP(:,:)

    integer, allocatable :: BCType(:,:)
    integer, allocatable :: MapB(:)
    integer, allocatable :: VMapB(:)
    ! integer, allocatable :: MapD(:)
    ! integer, allocatable :: VMapD(:)
    ! integer, allocatable :: MapN(:)
    ! integer, allocatable :: VMapN(:)

    integer :: tileID
    integer :: panelID
    integer :: PRC_myrank
    integer :: lcdomID
  
    real(RP), allocatable :: G_ij(:,:,:,:) !< The covariant component of metric tensor with horizontal general curvilinear coordinate 
    real(RP), allocatable :: GIJ(:,:,:,:)  !< The contravariant component of metric tensor with horizontal general curvilinear coordinate 
    real(RP), allocatable :: Gsqrt(:,:)    !< The Jacobian of 3D transformation in the computational coordinate (=GsqrtH * GsqrtV)
  end type LocalMeshBase

  public :: LocalMeshBase_Init
  public :: LocalMeshBase_Final

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: BCTYPE_INTERIOR  = 0
  integer, public, parameter :: BCTYPE_UNDEFBC   = 1
  integer, public, parameter :: BCTYPE_PERIODIC  = 2
  integer, public, parameter :: BCTYPE_DIRCHLET  = 3
  integer, public, parameter :: BCTYPE_NEUMAN    = 4
  integer, public, parameter :: BCTYPE_SHORELINE = 5

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
  subroutine LocalMeshBase_Init( this, lcdomID, refElem, ndim, myrank )
    
    use scale_prc, only: PRC_myrank
    implicit none

    class(LocalMeshBase), intent(inout) :: this
    integer, intent(in) :: lcdomID
    class(ElementBase), intent(in), target :: refElem
    integer, intent(in) :: ndim
    integer, intent(in), optional :: myrank
    !-----------------------------------------------------------------------------

    this%lcdomID = lcdomID
    this%refElem => refElem

    if (present(myrank)) then
      this%PRC_myrank = myrank
    else
      this%PRC_myrank = PRC_myrank
    end if

    return
  end subroutine LocalMeshBase_Init

  subroutine LocalMeshBase_Final( this, is_generated )
    implicit none
    
    class(LocalMeshBase), intent(inout) :: this
    logical, intent(in) :: is_generated
    !-----------------------------------------------------------------------------

    if ( is_generated ) then
      deallocate( this%pos_ev, this%pos_en, this%normal_fn )
      deallocate( this%sJ, this%J )

      deallocate( this%Escale, this%Fscale )

      deallocate( this%EToV, this%EToE, this%EToF )

      if ( allocated(this%VMapM) ) then
        deallocate( this%VMapM, this%VMapP, this%MapM, this%MapP )
      end if
      if ( allocated(this%VMapB) ) then
        deallocate( this%BCType )
        deallocate( this%VMapB, this%MapB )
      end if
      if ( allocated(this%G_ij) ) then
        deallocate( this%G_ij, this%GIJ )
        deallocate( this%Gsqrt )
      end if
    end if
    
    return
  end subroutine LocalMeshBase_Final
end module scale_localmesh_base