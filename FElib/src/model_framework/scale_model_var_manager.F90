#include "scaleFElib.h"
module scale_model_var_manager
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_meshfield_base, only: &
    MeshField1D, MeshField2D, MeshField3D
  use scale_linkedlist, only: &
    LinkedList
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
   
  type, public :: ModelVarManager
    type(LinkedList) :: list
  contains
    procedure, public :: Init => ModelVarManager_Init
    procedure, public :: Final => ModelVarManager_Final
    procedure, private :: Regist1D => ModelVarManager_regist1D
    procedure, private :: Regist2D => ModelVarManager_regist1D
    procedure, private :: Regist3D => ModelVarManager_regist1D
    generic, public :: Regist => Regist1D, Regist2D, Regist3D
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
    return
  end subroutine ModelVarManager_Init

  subroutine ModelVarManager_Final( this )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    !------------------------------------------------

    call this%list%Traverse( field_release )
    call this%list%Final()
    
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

  subroutine ModelVarManager_regist1D( this, key_id, field )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    integer, intent(in) :: key_id
    class(MeshField1D), target, intent(in) :: field

    class(*), pointer :: ptr_field
    !------------------------------------------------

    ptr_field => field
    call this%list%AddByPointer( key_id, ptr_field )
    return
  end subroutine ModelVarManager_regist1D

  subroutine ModelVarManager_regist2D( this, key_id, field )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    integer, intent(in) :: key_id
    class(MeshField2D), target, intent(in) :: field

    class(*), pointer :: ptr_field
    !------------------------------------------------

    ptr_field => field
    call this%list%AddByPointer( key_id, ptr_field ) 
    return
  end subroutine ModelVarManager_regist2D

  subroutine ModelVarManager_regist3D( this, key_id, field )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    integer, intent(in) :: key_id
    class(MeshField3D), target, intent(in) :: field

    class(*), pointer :: ptr_field
    !------------------------------------------------

    ptr_field => field
    call this%list%AddByPointer( key_id, ptr_field )
    return
  end subroutine ModelVarManager_regist3D

end module scale_model_var_manager