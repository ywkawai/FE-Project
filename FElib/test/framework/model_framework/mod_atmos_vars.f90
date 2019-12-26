#include "scalelib.h"
module mod_atmos_vars
  use scale_precision
  use scale_prc
  use scale_io

  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d

  use scale_localmeshfield_base, only: LocalMeshField1D
  use scale_meshfield_base, only: MeshField1D

  use scale_model_var_manager, only: &
    ModelVarManager  
  
  use mod_atmos_mesh, only: AtmosMesh

  implicit none
  private

  type, public :: AtmosVars
    type(MeshField1D) :: V
    type(MeshField1D) :: T
    type(ModelVarManager) :: vars_list
  contains
    procedure :: Init => AtmosVars_Init
    procedure :: Final => AtmosVars_Final
  end type AtmosVars

  integer, public, parameter :: ATMOSVARS_V_ID = 1
  integer, public, parameter :: ATMOSVARS_T_ID = 1

contains

subroutine AtmosVars_Init( this, mesh )
  implicit none
  class(AtmosVars), intent(inout) :: this
  class(AtmosMesh), intent(in) :: mesh

  integer :: n
  !--------------------------------------------------

  LOG_INFO('AtmosVars_Init',*)

  call this%vars_list%Init()

  call this%V%Init("V", "m/s", mesh%mesh)
  call this%vars_list%Regist( ATMOSVARS_V_ID, this%V )

  call this%T%Init("T", "K", mesh%mesh)
  call this%vars_list%Regist( ATMOSVARS_T_ID, this%T )
  
  return
end subroutine AtmosVars_Init

subroutine AtmosVars_Final( this )
  implicit none
  class(AtmosVars), intent(inout) :: this

  class(ModelVarManager), pointer :: ptr_base
  !--------------------------------------------------

  LOG_INFO('AtmosVars_Final',*)

  call this%vars_list%Final()

  return
end subroutine AtmosVars_Final

!-----------------------------
elemental function get_field_val(tileID, k, p) result(val)
  implicit none
  integer, intent(in) :: tileID, k, p
  real(RP) :: val
  !----------------------------
  val = tileID*1000000 + k*1000 + p

  return
end function get_field_val

end module mod_atmos_vars