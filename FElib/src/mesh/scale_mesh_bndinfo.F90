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
    character(len=H_SHORT) :: tag
  contains
    procedure :: Init => MeshBndInfo_Init
    procedure :: Final => MeshBndInfo_Final
    procedure, private :: Set_by_ID => MeshBndInfo_set_by_ID
    procedure, private :: Set_by_name => MeshBndInfo_set_by_ID
    generic :: Set => Set_by_ID, Set_by_name
  end type MeshBndInfo
  
  public :: BndType_NameToID

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=*), public , parameter :: BND_TYPE_NOSPEC_NAME  = 'NONSPEC'
  integer, public :: BND_TYPE_NOSPEC_ID                         = 0
  character(len=*), public, parameter :: BND_TYPE_PERIODIC_NAME = 'PERIODIC'
  integer, public, parameter  :: BND_TYPE_PERIODIC_ID           = 1
  character(len=*), public, parameter  :: BND_TYPE_SLIP_NAME    = 'SLIP'
  integer, public :: BND_TYPE_SLIP_ID                           = 2
  character(len=*), public, parameter  :: BND_TYPE_NOSLIP_NAME  = 'NOSLIP'
  integer, public :: BND_TYPE_NOSLIP_ID                         = 3
  character(len=*), public , parameter :: BND_TYPE_ADIABAT_NAME = 'ADIABATIC'
  integer, public :: BND_TYPE_ADIABAT_ID                        = 5

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------
 
contains
  subroutine MeshBndInfo_Init(this, list_size, tag)
    implicit none
    class(MeshBndInfo), intent(inout) :: this
    integer, intent(in) :: list_size
    character(*), optional, intent(in) :: tag
    !------------------------------------------------------

    allocate( this%list(list_size) )
    if (present(tag)) then
      this%tag = tag
    else
      this%tag = ''
    end if
    
    return
  end subroutine MeshBndInfo_Init

  subroutine MeshBndInfo_Final(this)
    implicit none
    class(MeshBndInfo), intent(inout) :: this
    !------------------------------------------------------

    if (allocated(this%list)) deallocate(this%list)
    
    return
  end subroutine MeshBndInfo_Final  

  subroutine MeshBndInfo_set_by_ID(this, is, ie, bnd_type_id)
    implicit none
    class(MeshBndInfo), intent(inout) :: this
    integer, intent(in) :: is
    integer, intent(in) :: ie
    integer, intent(in) :: bnd_type_id
    !------------------------------------------------------

    this%list(is:ie) = bnd_type_id
    return
  end subroutine MeshBndInfo_set_by_Id

  subroutine MeshBndInfo_set_by_name(this, is, ie, bnd_type_name)
    implicit none
    class(MeshBndInfo), intent(inout) :: this
    integer, intent(in) :: is
    integer, intent(in) :: ie
    character(*), intent(in) :: bnd_type_name

    integer :: bnd_type_id
    !------------------------------------------------------
    bnd_type_id = BndType_NameToID(bnd_type_name)
    call MeshBndInfo_set_by_ID(this, is, ie, bnd_type_id)

    return
  end subroutine MeshBndInfo_set_by_name

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
    case default
      LOG_ERROR('BndType_NameToID ',*) trim(bnd_type_name) // ' is not supported. Check!'
      call PRC_abort
    end select

    return
  end function BndType_NameToID 

end module scale_mesh_bndinfo