!-------------------------------------------------------------------------------
!> module ATMOSPHERE driver
!!
!! @par Description
!!          Atmosphere module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use mod_atmos_component
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, extends(scale_model_component) :: AtmosComponent
  contains
    procedure, public :: setup => AtmosComponent_setup
    procedure, public :: calc_tendency => AtmosComponent_calc_tendency
    procedure, public :: update => AtmosComponent_update
    procedure, public :: finalize => AtmosComponent_finalize
  end type

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  subroutine AtmosComponent_setup( this )
    implicit none

    class(AtmosComponent), intent(inout) :: this

    !----------------------------------------
  end subroutine AtmosComponent_setup
  

  subroutine AtmosComponent_calc_tendency( this )
    implicit none
    
    class(AtmosComponent), intent(inout) :: this

    !----------------------------------------
  end subroutine AtmosComponent_calc_tendency

  subroutine AtmosComponent_update( this )
    implicit none
    
    class(AtmosComponent), intent(inout) :: this

    !----------------------------------------
  end subroutine AtmosComponent_update

  subroutine AtmosComponent_finalize( this )
    implicit none
    
    class(AtmosComponent), intent(inout) :: this

    !----------------------------------------
  end subroutine AtmosComponent_finalize

end module mod_atmos_driver