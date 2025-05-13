!-------------------------------------------------------------------------------
!> module USER_base
!!
!! @par Description
!!          User defined module (base class)
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_user_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use mod_atmos_component, only: AtmosComponent
  use mod_experiment, only: Experiment

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  type, public :: UserBase
    logical :: USER_do !< do user step?
  contains
    procedure :: mkinit_base => USER_base_mkinit
    generic :: mkinit => mkinit_base
    procedure :: mkfinal => USER_base_mkfinal
    procedure :: setup_base => USER_base_setup
    generic :: setup => setup_base
    procedure :: final => USER_base_final
    procedure :: calc_tendency => USER_base_calc_tendency
    procedure :: update_pre => USER_base_update_pre
    procedure :: update => USER_base_update
  end type UserBase

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
  subroutine USER_base_mkinit( this, atm, exp )
    implicit none
    class(UserBase), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm
    class(Experiment), intent(inout) :: exp
    !------------------------------------------

    call exp%SetInitCond( atm%mesh,                        &
      atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager, &
      atm%vars%QTRCVARS_manager                            )
    
    return
  end subroutine USER_base_mkinit

  subroutine USER_base_mkfinal( this )
    implicit none
    class(UserBase), intent(inout) :: this
    !------------------------------------------
    return
  end subroutine USER_base_mkfinal

  subroutine USER_base_setup( this, atm, user_do )
    implicit none
    class(UserBase), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm
    logical, intent(in) :: user_do
    !------------------------------------------
    
    this%USER_do = user_do
    
    return
  end subroutine USER_base_setup

  subroutine USER_base_final( this )
    implicit none
    class(UserBase), intent(inout) :: this
    !------------------------------------------
    return
  end subroutine USER_base_final
  
  subroutine USER_base_calc_tendency( this, atm )
    implicit none
    class(UserBase), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_base_calc_tendency

  subroutine USER_base_update_pre( this, atm )
    implicit none
    class(UserBase), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------
    return
  end subroutine USER_base_update_pre

  subroutine USER_base_update( this, atm )
    implicit none
    class(UserBase), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------
    return
  end subroutine USER_base_update

end module mod_user_base
