!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_mesh_bndinfo
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  type, public :: MeshBndInfo
    integer, allocatable :: list(:)
    real(RP), allocatable :: val(:)
    character(len=H_SHORT) :: tag
  contains
    procedure :: Init => MeshBndInfo_Init
    procedure :: Final => MeshBndInfo_Final
    procedure, private :: Set_by_ID => MeshBndInfo_set_by_ID
    procedure, private :: Set_by_name => MeshBndInfo_set_by_name
    generic :: Set => Set_by_ID, Set_by_name
  end type MeshBndInfo
  
  public :: BndType_NameToID

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=*), public , parameter :: BND_TYPE_NOSPEC_NAME  = 'NONSPEC'
  integer, public, parameter :: BND_TYPE_NOSPEC_ID              = 0
  character(len=*), public, parameter :: BND_TYPE_PERIODIC_NAME = 'PERIODIC'
  integer, public, parameter  :: BND_TYPE_PERIODIC_ID           = 1
  character(len=*), public, parameter  :: BND_TYPE_SLIP_NAME    = 'SLIP'
  integer, public, parameter :: BND_TYPE_SLIP_ID                = 2
  character(len=*), public, parameter  :: BND_TYPE_NOSLIP_NAME  = 'NOSLIP'
  integer, public, parameter:: BND_TYPE_NOSLIP_ID               = 3
  character(len=*), public , parameter :: BND_TYPE_ADIABAT_NAME = 'ADIABATIC'
  integer, public, parameter :: BND_TYPE_ADIABAT_ID             = 4
  character(len=*), public , parameter :: BND_TYPE_FIXVAL_NAME  = 'FIXVAL'
  integer, public, parameter :: BND_TYPE_FIXVAL_ID              = 5

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------
 
contains
!OCL SERIAL
  subroutine MeshBndInfo_Init(this, list_size, tag)
    use scale_const, only: &
      UNDEF8 => CONST_UNDEF8
    implicit none
    class(MeshBndInfo), intent(inout) :: this
    integer, intent(in) :: list_size
    character(*), optional, intent(in) :: tag

    integer :: i
    !------------------------------------------------------

    allocate( this%list(list_size) )
    allocate( this%val(list_size) )

    !$omp parallel do
    do i = 1, list_size
      this%list(i) = BND_TYPE_NOSPEC_ID
      this%val(i)  = UNDEF8
    end do

    if (present(tag)) then
      this%tag = tag
    else
      this%tag = ''
    end if
    
    return
  end subroutine MeshBndInfo_Init

!OCL SERIAL
  subroutine MeshBndInfo_Final(this)
    implicit none
    class(MeshBndInfo), intent(inout) :: this
    !------------------------------------------------------

    if (allocated(this%list)) deallocate(this%list)
    if (allocated(this%val)) deallocate(this%val)
    
    return
  end subroutine MeshBndInfo_Final  

!OCL SERIAL
  subroutine MeshBndInfo_set_by_ID(this, is, ie, bnd_type_id, val)
    implicit none
    class(MeshBndInfo), intent(inout) :: this
    integer, intent(in) :: is
    integer, intent(in) :: ie
    integer, intent(in) :: bnd_type_id
    real(RP), intent(in), optional :: val
    !------------------------------------------------------

    this%list(is:ie) = bnd_type_id
    if ( present(val) ) then
      this%val(is:ie) = val
    end if

    return
  end subroutine MeshBndInfo_set_by_Id

!OCL SERIAL
  subroutine MeshBndInfo_set_by_name(this, is, ie, bnd_type_name, val)
    implicit none
    class(MeshBndInfo), intent(inout) :: this
    integer, intent(in) :: is
    integer, intent(in) :: ie
    character(*), intent(in) :: bnd_type_name
    real(RP), intent(in), optional :: val

    integer :: bnd_type_id
    !------------------------------------------------------
    bnd_type_id = BndType_NameToID(bnd_type_name)
    call MeshBndInfo_set_by_ID(this, is, ie, bnd_type_id, val)

    return
  end subroutine MeshBndInfo_set_by_name

!OCL SERIAL
  function BndType_NameToID(bnd_type_name) result(bnd_type_id)
    implicit none
    character(*), intent(in) :: bnd_type_name

    integer :: bnd_type_id
    !------------------------------------------------------

    select case(trim(bnd_type_name))
    case (BND_TYPE_NOSPEC_NAME)
      bnd_type_id = BND_TYPE_NOSPEC_ID
    case (BND_TYPE_PERIODIC_NAME)
      bnd_type_id = BND_TYPE_PERIODIC_ID
    case (BND_TYPE_SLIP_NAME)
      bnd_type_id = BND_TYPE_SLIP_ID
    case (BND_TYPE_NOSLIP_NAME)
      bnd_type_id = BND_TYPE_NOSLIP_ID
    case (BND_TYPE_ADIABAT_NAME)
      bnd_type_id = BND_TYPE_ADIABAT_ID
    case (BND_TYPE_FIXVAL_NAME)
      bnd_type_id = BND_TYPE_FIXVAL_ID      
    case default
      LOG_ERROR('BndType_NameToID ',*) trim(bnd_type_name) // ' is not supported. Check!'
      call PRC_abort
    end select

    return
  end function BndType_NameToID 

end module scale_mesh_bndinfo