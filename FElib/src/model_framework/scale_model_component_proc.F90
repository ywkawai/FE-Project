#include "scaleFElib.h"
module scale_model_component_proc
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  
  type, abstract, public :: ModelComponentProc
    character(len=H_SHORT) :: name

    real(DP) :: dtsec
  contains
    procedure(ModelComponentProc_setup), deferred, public :: setup
    procedure(ModelComponentProc_calc_tendency), deferred, public :: calc_tendency
    procedure(ModelComponentProc_update), deferred, public :: update
    procedure(ModelComponentProc_finalize), deferred, public :: finalize    
  end type ModelComponentProc

  interface
    subroutine ModelComponentProc_setup( this )
      import ModelComponentProc
      class(ModelComponentProc), intent(inout) :: this
    end subroutine ModelComponentProc_setup

    subroutine ModelComponentProc_calc_tendency( this )
      import ModelComponentProc
      class(ModelComponentProc), intent(inout) :: this
    end subroutine ModelComponentProc_calc_tendency

    subroutine ModelComponentProc_update( this )
      import ModelComponentProc
      class(ModelComponentProc), intent(inout) :: this
    end subroutine ModelComponentProc_update

    subroutine ModelComponentProc_finalize( this )
      import ModelComponentProc
      class(ModelComponentProc), intent(inout) :: this
    end subroutine ModelComponentProc_finalize    
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


end module scale_model_component_proc