#include "scaleFElib.h"
module scale_model_component_proc
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use scale_time_manager, only: &
    TIME_manager_component
  use scale_model_mesh_manager, only: &
    ModelMeshBase
  use scale_model_var_manager, only: &
    ModelVarManager
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  
  type, abstract, public :: ModelComponentProc
    character(len=H_SHORT) :: name
    logical, private :: is_activated = .false.    
    integer :: tm_process_id
  contains
    procedure(ModelComponentProc_setup), deferred, public :: setup
    procedure(ModelComponentProc_calc_tendency), deferred, public :: calc_tendency
    procedure(ModelComponentProc_update), deferred, public :: update
    procedure(ModelComponentProc_finalize), deferred, public :: finalize  
    
    procedure, public :: ModelComponentProc_Init    
    procedure, public :: IsActivated => ModelComponentProc_isActivated    
  end type ModelComponentProc

  interface
    subroutine ModelComponentProc_setup( this, model_mesh, tm_parent_comp )
      import ModelComponentProc
      import ModelMeshBase
      import TIME_manager_component
      class(ModelComponentProc), intent(inout) :: this
      class(ModelMeshBase), target, intent(in) :: model_mesh
      class(TIME_manager_component), intent(inout) :: tm_parent_comp
    end subroutine ModelComponentProc_setup

    subroutine ModelComponentProc_calc_tendency( this, model_mesh, prgvars_list, trcvars_list, auxvars_list, forcing_list, is_update )
      import ModelComponentProc
      import ModelMeshBase
      import ModelVarManager      
      class(ModelComponentProc), intent(inout) :: this
      class(ModelMeshBase), intent(in) :: model_mesh
      class(ModelVarManager), intent(inout) :: prgvars_list
      class(ModelVarManager), intent(inout) :: trcvars_list  
      class(ModelVarManager), intent(inout) :: auxvars_list         
      class(ModelVarManager), intent(inout) :: forcing_list 
      logical, intent(in) :: is_update
    end subroutine ModelComponentProc_calc_tendency

    subroutine ModelComponentProc_update( this, model_mesh, prgvars_list, trcvars_list, auxvars_list, forcing_list, is_update )
      import ModelComponentProc
      import ModelMeshBase
      import ModelVarManager
      class(ModelComponentProc), intent(inout) :: this
      class(ModelMeshBase), intent(in) :: model_mesh
      class(ModelVarManager), intent(inout) :: prgvars_list
      class(ModelVarManager), intent(inout) :: trcvars_list      
      class(ModelVarManager), intent(inout) :: auxvars_list
      class(ModelVarManager), intent(inout) :: forcing_list
      logical, intent(in) :: is_update
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

contains
  subroutine ModelComponentProc_Init( this, name, is_activated )
    implicit none

    class(ModelComponentProc), intent(inout) :: this
    character(len=*), intent(in) :: name
    logical, intent(in) :: is_activated 
    !---------------------------------------------------
    
    this%name = name
    this%is_activated = is_activated

    return
  end subroutine ModelComponentProc_Init

  function ModelComponentProc_isActivated( this ) result( is_activated )
    implicit none

    class(ModelComponentProc), intent(in) :: this
    logical :: is_activated 
    !---------------------------------------------------

    is_activated = this%is_activated
    return
  end function ModelComponentProc_isActivated

end module scale_model_component_proc