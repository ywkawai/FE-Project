#include "scalelib.h"
module mod_atmos_component
  use scale_precision
  use scale_prc
  use scale_io

  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d

  use scale_localmeshfield_base, only: LocalMeshField1D
  use scale_meshfield_base, only: MeshField1D

  use scale_model_component, only: &
    ModelComponent
  use scale_model_component_proc, only: &
    ModelComponentProc  
  
  use mod_atmos_vars, only: AtmosVars
  use mod_atmos_mesh, only: AtmosMesh
  use mod_atmos_dyn, only: AtmosDynProc

  implicit none
  private

  type, extends(ModelComponent), public :: AtmosComponent
    type(AtmosVars) :: vars
    type(AtmosMesh) :: mesh

    type(AtmosDynProc) :: dyn_proc

  contains
    procedure, public :: setup => Atmos_setup 
    procedure, public :: calc_tendency => Atmos_calc_tendency
    procedure, public :: update => Atmos_update
    procedure, public :: finalize => Atmos_finalize
  end type AtmosComponent


contains

subroutine Atmos_setup( this )
  implicit none
  class(AtmosComponent), intent(inout) :: this

  integer :: n
  !--------------------------------------------------

  LOG_INFO('AtmosComponent_setup',*)
  this%name = 'Atmos'
  
  call this%mesh%Init()
  call this%vars%Init( this%mesh )

  !-------------------------------------
  call this%dyn_proc%setup()

  return
end subroutine Atmos_setup

subroutine Atmos_calc_tendency( this )
  implicit none
  class(AtmosComponent), intent(inout) :: this

  !--------------------------------------------------

  LOG_INFO('AtmosComponent_calc_tendency',*)

  call this%dyn_proc%calc_tendency()

  return  
end subroutine Atmos_calc_tendency

subroutine Atmos_update( this )
  implicit none
  class(AtmosComponent), intent(inout) :: this
  !--------------------------------------------------

  LOG_INFO('AtmosComponent_update',*)

  call this%dyn_proc%update()

  return  
end subroutine Atmos_update

subroutine Atmos_finalize( this )
  implicit none
  class(AtmosComponent), intent(inout) :: this
  !--------------------------------------------------

  LOG_INFO('AtmosComponent_finalize',*)

  call this%dyn_proc%finalize()
  call this%vars%Final()
  call this%mesh%Final()

  return  
end subroutine Atmos_finalize

end module mod_atmos_component