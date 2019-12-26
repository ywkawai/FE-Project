#include "scalelib.h"
program test_linkedlist
  use scale_precision
  use scale_io
  use scale_linkedlist
  implicit none

  type(LinkedList) :: list
  character(len=H_MID), target :: item_s1
  character(len=H_MID), target :: item_s2
  integer :: i
  class(*), pointer :: ptr_item 
  !-------------------------------

  write(*,*) "test LinkedList..."
  call list%Init()

  write(*,*) "Add item by pointer ... (1, .., 10) "
  do i=1, 5
    item_s1 = gen_item(i)
    ptr_item => item_s1
    call list%AddByPointer( i, ptr_item )
  end do
  write(*,*) "Add item by clone ... (1, .., 10) "
  do i=6, 10
    item_s2 = gen_item(i)
    ptr_item => item_s2
    call list%AddByClone( i, ptr_item )
  end do
  call list%Traverse( print_items )

  write(*,*) "Remove item by key (key=8)"
  call list%Remove(8)
  call list%Traverse( print_items )

  write(*,*) "Add item by clone (key=11)"
  item_s2 = gen_item(11)
  ptr_item => item_s2  
  call list%AddByClone(11, item_s2)
  call list%Traverse( print_items )

  call list%Final()
  write(*,*) "test_linkedlist has been succeeded!"


contains
  function gen_item( id ) result(item_s)
    integer, intent(in) :: id
    character(len=H_MID) :: item_s
    !------------------------------------------

    write(item_s,'(a,i3.3)') "item_", id   
    return
  end function gen_item

  subroutine print_items( key, val, done )
    implicit none
    class(*), intent(in) :: key
    class(*), pointer    :: val
    logical, intent(out) :: done

    integer :: key_i
    character(len=H_MID) :: val_s
    !------------------------------------------
    select type(key)
    type is (integer)
      key_i = key
    end select
    select type(val)
    type is (character(len=*))
      val_s = val
    end select

    write(*,'(a,i2,a,a)') "key=", key_i, ": val=", trim(val_s)
    done = .false.

    return
  end subroutine print_items

end program test_linkedlist