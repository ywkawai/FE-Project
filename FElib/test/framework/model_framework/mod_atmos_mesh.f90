#include "scalelib.h"
module mod_atmos_mesh
  use scale_precision
  use scale_prc
  use scale_io

  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d
  use scale_model_mesh_manager, only: ModelMesh1D
  
  implicit none
  private

  
  type, extends(ModelMesh1D), public :: AtmosMesh
    type(MeshLineDom1D) :: mesh
    type(LineElement) :: element
  contains
    procedure :: Init => AtmosMesh_Init
    procedure :: Final => AtmosMesh_Final
  end type AtmosMesh

  integer, parameter :: DEFAULT_NLocalMeshPerPrc = 1
  real(RP), parameter :: DEFAULT_dom_xmin = -1.0_RP
  real(RP), parameter :: DEFAULT_dom_xmax = +1.0_RP
  integer, parameter :: DEFAULT_NeGX = 4
  integer, parameter :: DEFAULT_POLYORDER = 2

contains

subroutine AtmosMesh_Init( this )
  implicit none
  class(AtmosMesh), intent(inout) :: this

  !--------------------------------------------------

  LOG_INFO('AtmosMesh_Init',*)
  
  call this%element%Init( DEFAULT_PolyOrder, .true.)

  call this%mesh%Init( &
    DEFAULT_NeGX,                              &
    DEFAULT_dom_xmin, DEFAULT_dom_xmax,        &
    this%element, DEFAULT_NLocalMeshPerPrc )
  
  call this%mesh%Generate()

  return
end subroutine AtmosMesh_Init

subroutine AtmosMesh_Final( this )
  implicit none
  class(AtmosMesh), intent(inout) :: this
  !--------------------------------------------------

  LOG_INFO('AtmosMesh_Final',*)

  call this%mesh%Final()
  call this%element%Final()

  return
end subroutine AtmosMesh_Final

end module mod_atmos_mesh