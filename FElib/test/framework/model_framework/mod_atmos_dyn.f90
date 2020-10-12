#include "scalelib.h"
module mod_atmos_dyn
  use scale_precision
  use scale_prc
  use scale_io

  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d

  use scale_model_component_proc, only: &
    ModelComponentProc

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_time_manager, only: TIME_manager_component

  implicit none
  private

  type, extends(ModelComponentProc), public :: AtmosDynProc
  contains
    procedure, public :: setup => AtmosDyn_setup 
    procedure, public :: calc_tendency => AtmosDyn_calc_tendency
    procedure, public :: update => AtmosDyn_update
    procedure, public :: finalize => AtmosDyn_finalize
  end type AtmosDynProc


contains

subroutine AtmosDyn_setup( this, model_mesh, tm_parent_comp )
  implicit none
  class(AtmosDynProc), intent(inout) :: this
  class(ModelMeshBase), target, intent(in) :: model_mesh
  class(TIME_manager_component), intent(inout) :: tm_parent_comp

  !--------------------------------------------------

  LOG_INFO('AtmosDyn_setup',*)
  this%name = 'AtmosDyn'

  return  
end subroutine AtmosDyn_setup

subroutine AtmosDyn_calc_tendency( this, model_mesh, prgvars_list, auxvars_list  )
  implicit none
  class(AtmosDynProc), intent(inout) :: this
  class(ModelMeshBase), intent(in) :: model_mesh
  class(ModelVarManager), intent(inout) :: prgvars_list
  class(ModelVarManager), intent(inout) :: auxvars_list  
  !--------------------------------------------------

  LOG_INFO('AtmosDyn_tendency',*)

  return  
end subroutine AtmosDyn_calc_tendency

subroutine AtmosDyn_update( this, model_mesh, prgvars_list, auxvars_list  )
  implicit none
  class(AtmosDynProc), intent(inout) :: this
  class(ModelMeshBase), intent(in) :: model_mesh
  class(ModelVarManager), intent(inout) :: prgvars_list
  class(ModelVarManager), intent(inout) :: auxvars_list  
  !--------------------------------------------------

  LOG_INFO('AtmosDyn_update',*)

  return  
end subroutine AtmosDyn_update

subroutine AtmosDyn_finalize( this )
  implicit none
  class(AtmosDynProc), intent(inout) :: this
  !--------------------------------------------------

  LOG_INFO('AtmosDyn_finalize',*)

  return  
end subroutine AtmosDyn_finalize  

end module mod_atmos_dyn