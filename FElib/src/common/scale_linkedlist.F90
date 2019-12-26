#include "scaleFElib.h"
module scale_linkedlist
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  
  !---
  type, abstract, public :: LinkedListKey
  contains
    procedure(key_equal_func), private, deferred :: key_equal
    generic, public :: operator(==) => key_equal
  end type LinkedListKey

  abstract interface
    pure elemental logical function key_equal_func( item1, item2 )
      import LinkedListKey
      class(LinkedListKey), intent(in) :: item1
      class(LinkedListKey), intent(in) :: item2
    end function key_equal_func
  end interface

  !---
  type :: Node
    private
    class(*), allocatable :: key
    class(*), pointer     :: value => null()
    logical :: destroy_on_delete
    type(node), pointer :: next => null()
    type(node), pointer :: previous => null()
  contains
    procedure, public :: Init => Node_Init
    procedure, public :: Final => Node_Final
    procedure, public :: GetData => Node_get_data
  end type


  type, public :: LinkedList
    type(node), pointer, private :: head
    type(node), pointer, private :: tail
    integer, private :: counter
  contains
    procedure, public :: Init  => LinkedList_Init
    procedure, public :: Final => LinkedList_Final
    procedure, public :: AddByPointer => LinkedList_Add_by_pointer
    procedure, public :: RemoveByPointer => LinkedList_remove_by_pointer
    procedure, public :: Remove => LinkedList_remove_by_key
    procedure, public :: AddByClone => LinkedList_add_by_clone
    procedure, public :: GetNode => LinkedList_get_node
    procedure, public :: Get => LinkedList_get_data
    procedure, public :: TraverseList =>  LinkedList_traverse_list
    procedure, public :: Traverse =>  LinkedList_traverse

    procedure, private :: keysEqual => LinkedList_keys_eqaul
  end type 

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !------------------

  abstract interface
    subroutine iterator_func(this, done)
      import :: node
      implicit none
      type(node), pointer  :: this
      logical, intent(out) :: done
    end subroutine iterator_func

    subroutine key_iterator(key, value, done)
      implicit none
      class(*), intent(in) :: key
      class(*), pointer    :: value
      logical, intent(out) :: done
    end subroutine key_iterator   
  end interface

contains
  subroutine Node_Init( this, key, value, previous_node, destroy_on_delete ) 
    implicit none
    class(Node) :: this
    class(*), intent(in) :: key
    class(*), pointer, intent(in) :: value
    type(node), pointer :: previous_node
    logical, intent(in) :: destroy_on_delete
    !---------------------------------------------

    allocate( this%key, source=key )
    this%value => value
    this%previous => previous_node
    this%destroy_on_delete = destroy_on_delete

    return
  end subroutine Node_Init

  subroutine Node_Final( this ) 
    implicit none
    class(Node), intent(inout) :: this
    !---------------------------------------------

    if( allocated( this%key ) ) deallocate( this%key )
    if ( this%destroy_on_delete ) then
      if (associated( this%value) ) deallocate( this%value )
    end if

    return
  end subroutine Node_Final

  subroutine Node_get_data( this, value ) 
    implicit none
    class(Node) :: this
    class(*), pointer, intent(out) :: value
    !---------------------------------------------
    
    if ( associated(this%value) )then
      value => this%value
    else
      LOG_ERROR(" Node_get_data",*)   "The pointer into data of node is not associated. Check!"      
    end if

    return
  end subroutine Node_get_data

  !-------------------------------------------------------

  subroutine LinkedList_Init( this )
    implicit none
    class(LinkedList), intent(inout) :: this
    !---------------------------------------------

    this%counter = 0
    nullify( this%head, this%tail )

    return
  end subroutine LinkedList_Init

  subroutine LinkedList_Final( this )
    implicit none
    class(LinkedList), intent(inout) :: this
    !---------------------------------------------

    this%counter = 0
    if ( associated(this%head) ) call destroy_subsequent_nodes( this%head )

    return
  end subroutine LinkedList_Final

  recursive subroutine destroy_subsequent_nodes( this_node )
    implicit none
    type(Node), pointer, intent(inout) :: this_node

    !---------------------------------------------
    if ( associated(this_node) ) then
      call this_node%Final()
      call destroy_subsequent_nodes( this_node%next )
      nullify( this_node%previous )
      deallocate( this_node )
      nullify( this_node )
    end if

    return
  end subroutine destroy_subsequent_nodes


  function LinkedList_has_key( this, key ) result(has_key)
    implicit none
    class(LinkedList), intent(inout) :: this
    class(*), intent(in)             :: key
    logical                          :: has_key
    !---------------------------------------------

    has_key = .false.
    call this%TraverseList( key_search )
    return
  contains
    subroutine key_search(ptr, done)
      implicit none
      type(node), pointer :: ptr
      logical, intent(out) :: done
      !-------------------------------------------------
      has_key = this%keysEqual(ptr%key, key)
      done = has_key
      return
    end subroutine key_search
  end function LinkedList_has_key

  subroutine LinkedList_traverse_list( &
    this, iterator )

    implicit none
    class(LinkedList), intent(inout) :: this
    procedure(iterator_func) :: iterator

    type(node), pointer :: ptr
    logical :: done
    !----------------------------------------------------
    
    done = .false.
    ptr => this%head

    do 
      if (associated(ptr)) then
        call iterator( ptr, done )
        if (done) exit
        ptr => ptr%next
      else
        exit
      end if
    end do

    return
  end subroutine LinkedList_traverse_list

  subroutine LinkedList_traverse( &
    this, iterator )

    implicit none
    class(LinkedList), intent(inout) :: this
    procedure(key_iterator) :: iterator

    !----------------------------------------------------
    
    call this%TraverseList( key_iterator_wrapper )
    return

  contains
    subroutine key_iterator_wrapper( this_node, done )
      implicit none
      type(node), pointer :: this_node
      logical, intent(out) :: done
      !----------------------------------------------------

      call iterator( this_node%key, this_node%value, done )
      return
    end subroutine key_iterator_wrapper

  end subroutine LinkedList_traverse

  subroutine LinkedList_add_by_pointer( &
      this, key, value, destroy_on_delete )
    
    implicit none
    class(LinkedList), intent(inout) :: this
    class(*), intent(in) :: key
    class(*), pointer, intent(in) :: value
    logical, intent(in), optional :: destroy_on_delete

    type(node), pointer :: pNode
    logical :: destroy_on_delete_ = .false.
    !---------------------------------------------

    select type (key)
    type is (integer)
    type is (character(len=*))
    class is (LinkedListKey)
    class default
      LOG_ERROR("LinkedList_add_by_pointer",*)   "The type of key is invalid. Check!"
    end select

    call this%GetNode( key, pNode )
    if ( associated(pNode) ) call this%RemoveByPointer( pNode )

    allocate( pNode )
    if (present(destroy_on_delete)) then
      destroy_on_delete_ = destroy_on_delete
    end if
    call pNode%Init( key, value, this%tail, destroy_on_delete_ )

    if ( associated(this%tail) ) then
      this%tail%next => pNode
    else
      this%head => pNode
    end if
    this%tail => pNode
    this%counter = this%counter + 1

    return
  end subroutine LinkedList_add_by_pointer

  subroutine LinkedList_add_by_clone( this, key, value )
    implicit none
    class(LinkedList), intent(inout) :: this
    class(*), intent(in) :: key
    class(*), intent(in) :: value

    class(*), pointer :: ptr_value
    !---------------------------------------------

    allocate(ptr_value, source=value)
    call this%AddByPointer( key, ptr_value, destroy_on_delete=.true. )

    return
  end subroutine LinkedList_add_by_clone

  subroutine LinkedList_remove_by_key( this, key )
    implicit none
    class(LinkedList), intent(inout) :: this
    class(*), intent(in) :: key

    type(Node), pointer :: pNode
    !---------------------------------------------

    call this%GetNode( key, pNode )
    call this%RemoveByPointer( pNode )

    return
  end subroutine LinkedList_remove_by_key

  subroutine LinkedList_remove_by_pointer( this, pNode )
    implicit none
    class(LinkedList), intent(inout) :: this
    type(Node), pointer :: pNode

    logical :: has_next
    logical :: has_previous
    !---------------------------------------------

    if ( associated(pNode) ) then
      has_next = associated( pNode%next )
      has_previous = associated( pNode%previous )

      if ( has_next .and. has_previous ) then
        pNode%previous%next => pNode%next
        pNode%next%previous => pNode%previous
      else if (       has_next .and. .not. has_previous ) then
        this%head => pNode%next
        nullify( this%head%previous )        
      else if ( .not. has_next .and.       has_previous ) then
        this%tail => pNode%previous
        nullify( this%tail%next )
      else if ( .not. (has_next .or. has_previous) ) then
        nullify( this%head, this%tail )
      end if

      call pNode%Final()
      deallocate( pNode )
      nullify( pNode )

      this%counter = this%counter - 1
    end if

    return
  end subroutine LinkedList_remove_by_pointer

  subroutine LinkedList_get_node( this, key, ptr_node )
    implicit none
    class(LinkedList), intent(in) :: this
    class(*), intent(in) :: key
    type(node), pointer, intent(out) :: ptr_node

    type(node), pointer :: ptr
    !---------------------------------------------

    nullify( ptr_node )
    
    ptr => this%head
    do 
      if ( associated(ptr) ) then
        if (this%keysEqual(ptr%key, key)) then
          ptr_node => ptr
          return
        end if
        ptr => ptr%next
      else
        return
      end if
    end do

  end subroutine LinkedList_get_node

  subroutine LinkedList_get_data( this, key, ptr_value )
    implicit none
    class(LinkedList), intent(in) :: this
    class(*), intent(in) :: key
    class(*), pointer, intent(out) :: ptr_value

    type(node), pointer :: ptr
    !---------------------------------------------

    call this%GetNode( key, ptr )
    if ( associated(ptr) ) then
      ptr_value => ptr%value
    else
      nullify( ptr_value )
    end if

    return
  end subroutine LinkedList_get_data

  function LinkedList_keys_eqaul( this, key1, key2 ) result(is_keys_eqaul)
    implicit none
    class(LinkedList), intent(in) :: this
    class(*), intent(in) :: key1 
    class(*), intent(in) :: key2
    logical :: is_keys_eqaul

    !---------------------------------------------

    is_keys_eqaul = .false.

    if ( same_type_as(key1, key2) ) then
      select type (key1)
      class is (LinkedListKey)
        select type(key2)
        class is (LinkedListKey)
          is_keys_eqaul = (key1 == key2)
        end select      
      type is (integer)
        select type(key2)
        type is (integer)
          is_keys_eqaul = (key1 == key2)
        end select         
      type is (character(len=*))
        select type(key2)
        type is (character(len=*))
          is_keys_eqaul = (key1 == key2)
        end select    
      end select
    end if

    return
  end function LinkedList_keys_eqaul

end module scale_linkedlist