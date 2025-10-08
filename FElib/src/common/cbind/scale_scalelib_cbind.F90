#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_scalelib_cbind
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  subroutine scalelib_init( &
       app_name ) bind(C, name="CSCALElib_Init")
    use scale, only: SCALE_init    
    use scale_cbinding_util, only: &
      string_c2f => cbinding_util_string_c2f

    implicit none
    type(c_ptr), value :: app_name !< application name

    character(:), allocatable :: name
    !------------------------------------------------------------

    name = string_c2f(app_name)
    call SCALE_init( name )
    return
  end subroutine scalelib_init

  subroutine scalelib_final() bind(C, name="CSCALElib_Final")
    use scale, only: SCALE_finalize
    implicit none
    !------------------------------------------------------------
    call SCALE_finalize()
    return
  end subroutine scalelib_final

end module scale_scalelib_cbind