#include "scaleFElib.h"
module scale_model_var_manager
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io

  use scale_file_history, only: FILE_HISTORY_reg

  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField1D, MeshField2D, MeshField3D
  use scale_linkedlist, only: &
    LinkedList
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase, MeshFieldContainer
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: VariableInfo
    integer                :: keyID
    character(len=H_SHORT) :: NAME
    character(len=H_MID)   :: DESC
    character(len=H_SHORT) :: UNIT
    integer                :: ndims
    character(len=H_SHORT) :: dim_type
    character(len=H_MID)   :: STDNAME
  end type VariableInfo

  type, public :: ModelVarManager
    type(LinkedList) :: list
    class(MeshFieldCommBase), pointer :: ptr_comm
    type(MeshFieldContainer), allocatable :: comm_list(:)  
  contains
    procedure, public :: Init => ModelVarManager_Init
    procedure, public :: Final => ModelVarManager_Final
    procedure, private :: Regist1D => ModelVarManager_regist1D
    procedure, private :: Regist2D => ModelVarManager_regist2D
    procedure, private :: Regist3D => ModelVarManager_regist3D
    generic, public :: Regist => Regist1D, Regist2D, Regist3D
    procedure, public :: Get => ModelVarManager_Get

    !--
    procedure, public :: MeshFieldComm_Prepair => ModelVarManager_meshfiled_comm_prepare
    procedure, public :: MeshFieldComm_Exchange => ModelVarManager_meshfiled_comm_exchange
  end type ModelVarManager


  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !------------------

contains
  subroutine ModelVarManager_Init( this )
    implicit none
    class(ModelVarManager), intent(inout) :: this

    !------------------------------------------------
    call this%list%Init()
    nullify( this%ptr_comm )

    return
  end subroutine ModelVarManager_Init

  subroutine ModelVarManager_Final( this )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    !------------------------------------------------

    call this%list%Traverse( field_release )

    call this%list%Final()

    nullify( this%ptr_comm )
    if ( allocated(this%comm_list) ) deallocate(this%comm_list)

    return
  end subroutine ModelVarManager_Final
  subroutine field_release( key, pField, done)
    implicit none
    class(*), intent(in) :: key
    class(*), pointer    :: pField
    logical, intent(out) :: done

    select type( pField )
    type is ( MeshField1D )
      call pField%Final()
    type is ( MeshField2D )
      call pField%Final() 
    type is ( MeshField3D )
      call pField%Final()            
    end select

    return
  end subroutine field_release

  subroutine ModelVarManager_Regist1D( this,      &
      varinfo, mesh, field, reg_file_history_flag )
    
    use scale_mesh_base1d, only: Meshbase1d
    implicit none

    class(ModelVarManager), intent(inout) :: this
    type(VariableInfo), intent(in) :: varinfo
    class(MeshBase1D), intent(in) :: mesh
    class(MeshField1D), intent(inout), target :: field
    logical, intent(in) :: reg_file_history_flag

    class(*), pointer :: ptr_field    
    !------------------------------------------------

    call field%Init( varinfo%NAME, varinfo%UNIT, mesh )
    if (reg_file_history_flag) then
      call FILE_HISTORY_reg( varinfo%NAME, varinfo%DESC, varinfo%UNIT, field%hist_id, dim_type=varinfo%dim_type)
    end if

    ptr_field => field
    call this%list%AddByPointer( varinfo%keyID, ptr_field )

    return
  end subroutine ModelVarManager_Regist1D

  subroutine ModelVarManager_Regist2D( this,    &
    varinfo, mesh, field, reg_file_history_flag )
    use scale_mesh_base2d, only: Meshbase2d
    implicit none

    class(ModelVarManager), intent(inout) :: this
    type(VariableInfo), intent(in) :: varinfo
    class(MeshBase2D), intent(in) :: mesh
    class(MeshField2D), intent(inout), target :: field
    logical, intent(in) :: reg_file_history_flag

    class(*), pointer :: ptr_field    
    !------------------------------------------------

    call field%Init( varinfo%NAME, varinfo%UNIT, mesh )
    if (reg_file_history_flag) then
      call FILE_HISTORY_reg( varinfo%NAME, varinfo%DESC, varinfo%UNIT, field%hist_id, dim_type=varinfo%dim_type)
    end if

    ptr_field => field
    call this%list%AddByPointer( varinfo%keyID, ptr_field )

    return
  end subroutine ModelVarManager_Regist2D

  subroutine ModelVarManager_Regist3D( this,    &
    varinfo, mesh, field, reg_file_history_flag )
    use scale_mesh_base3d, only: MeshBase3D
    implicit none

    class(ModelVarManager), intent(inout) :: this
    type(VariableInfo), intent(in) :: varinfo
    class(MeshBase3D), intent(in) :: mesh
    class(MeshField3D), intent(inout), target :: field
    logical, intent(in) :: reg_file_history_flag

    class(*), pointer :: ptr_field    
    !------------------------------------------------

    call field%Init( varinfo%NAME, varinfo%UNIT, mesh )
    if (reg_file_history_flag) then
      call FILE_HISTORY_reg( field%varname, varinfo%DESC, field%unit, field%hist_id, &
                             dim_type=varinfo%dim_type )
    end if

    ptr_field => field
    call this%list%AddByPointer( varinfo%keyID, ptr_field )

    return
  end subroutine ModelVarManager_Regist3D

  subroutine ModelVarManager_Get( this, keyID, pField )
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    class(ModelVarManager), intent(inout) :: this
    integer, intent(in) :: keyID
    class(MeshFieldBase), pointer, intent(out) :: pField

    class(*), pointer :: ptr_field    
    !------------------------------------------------

    call this%list%Get(keyID, ptr_field)

    nullify( pField )
    select type( ptr_field )
    type is (MeshField1D)
      pField => ptr_field
    type is (MeshField2D)
      pField => ptr_field
    type is (MeshField3D)
      pField => ptr_field          
    end select
    
    return
  end subroutine ModelVarManager_Get

  subroutine ModelVarManager_meshfiled_comm_prepare( this, &
      comm,  fields )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    class(MeshFieldCommBase), target, intent(in) :: comm
    class(MeshFieldBase), target, intent(in) :: fields(:)

    integer :: v
    integer :: nFields
    class(MeshFieldBase), pointer :: ptr_field
    !------------------------------------------------

    this%ptr_comm => comm

    nFields = size(fields)
    allocate( this%comm_list(nFields) )

    do v = 1, size(fields)
      ptr_field => fields(v)
      select type( ptr_field )
      type is (MeshField1D)
        this%comm_list(v)%field1d => ptr_field
      type is (MeshField2D)
        this%comm_list(v)%field2d => ptr_field
      type is (MeshField3D)
        this%comm_list(v)%field3d => ptr_field               
      end select
    end do

    return
  end subroutine ModelVarManager_meshfiled_comm_prepare

  subroutine ModelVarManager_meshfiled_comm_exchange( this )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    !------------------------------------------------
    
    call this%ptr_comm%Put( this%comm_list, 1 )
    call this%ptr_comm%Exchange()
    call this%ptr_comm%Get(this%comm_list, 1 )

    return
  end subroutine ModelVarManager_meshfiled_comm_exchange

end module scale_model_var_manager