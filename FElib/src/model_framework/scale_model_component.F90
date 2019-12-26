#include "scaleFElib.h"
module scale_model_component
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use scale_model_component_proc
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  
  type, abstract, public :: ModelComponent
    character(len=H_SHORT) :: name
  contains
    procedure(ModelComponent_setup), deferred, public :: setup
    procedure(ModelComponent_calc_tendency), deferred, public :: calc_tendency
    procedure(ModelComponent_update), deferred, public :: update
    procedure(ModelComponent_finalize), deferred, public :: finalize    
  end type ModelComponent

  interface
    subroutine ModelComponent_setup( this )
      import ModelComponent
      class(ModelComponent), intent(inout) :: this
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


end module scale_model_component