#include "scaleFElib.h"
module scale_model_component
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use scale_time_manager, only: &
    TIME_manager_component

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  
  type, abstract, public :: ModelComponent
    character(len=H_SHORT), private :: name
    logical, private :: is_activated = .false.
    type(TIME_manager_component) :: time_manager
  contains
    procedure(ModelComponent_setup), deferred, public :: setup
    procedure(ModelComponent_calc_tendency), deferred, public :: calc_tendency
    procedure(ModelComponent_update), deferred, public :: update
    procedure(ModelComponent_finalize), deferred, public :: finalize 

    procedure, public :: ModelComponent_Init    
    procedure, public :: IsActivated => ModelComponent_isActivated
    procedure, public :: GetComponentName => ModelComponent_getCompName
  end type ModelComponent

  interface
    subroutine ModelComponent_setup( this )
      import ModelComponent
      class(ModelComponent), intent(inout), target :: this
    end subroutine ModelComponent_setup

    subroutine ModelComponent_calc_tendency( this )
      import ModelComponent
      class(ModelComponent), intent(inout) :: this
    end subroutine ModelComponent_calc_tendency

    subroutine ModelComponent_update( this )
      import ModelComponent
      class(ModelComponent), intent(inout) :: this
    end subroutine ModelComponent_update

    subroutine ModelComponent_finalize( this )
      import ModelComponent
      class(ModelComponent), intent(inout) :: this
    end subroutine ModelComponent_finalize    
  end interface

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
  subroutine ModelComponent_Init( this, name, is_activated )
    implicit none

    class(ModelComponent), intent(inout) :: this
    character(len=*), intent(in) :: name
    logical, intent(in) :: is_activated 
    !---------------------------------------------------
    
    this%name = name
    this%is_activated = is_activated

    return
  end subroutine ModelComponent_Init

  function ModelComponent_isActivated( this ) result( is_activated )
    implicit none

    class(ModelComponent), intent(in) :: this
    logical :: is_activated 
    !---------------------------------------------------

    is_activated = this%is_activated
    return
  end function ModelComponent_isActivated

  function ModelComponent_getCompName( this ) result( model_name )
    implicit none

    class(ModelComponent), intent(in) :: this
    character(len=H_SHORT) :: model_name
    !---------------------------------------------------

    model_name = this%name
    return
  end function ModelComponent_getCompName

end module scale_model_component