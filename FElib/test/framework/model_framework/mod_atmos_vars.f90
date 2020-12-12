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
    ModelVarManager, VariableInfo
  
  use mod_atmos_mesh, only: AtmosMesh

  implicit none
  private

  type, public :: AtmosVars
    type(MeshField1D) :: V
    type(MeshField1D) :: T
    type(MeshField1D) :: P
    type(MeshField1D) :: F_V
    type(MeshField1D) :: F_T
    type(ModelVarManager) :: prgvars_list
    type(ModelVarManager) :: auxvars_list
    type(ModelVarManager) :: forcing_list
  contains
    procedure :: Init => AtmosVars_Init
    procedure :: Final => AtmosVars_Final
  end type AtmosVars

  integer, public, parameter :: ATMOSVARS_V_ID = 1
  integer, public, parameter :: ATMOSVARS_T_ID = 2
  
  integer, public, parameter :: ATMOSVARS_P_ID = 1

  integer, public, parameter :: ATMOSVARS_F_V_ID = 1
  integer, public, parameter :: ATMOSVARS_F_T_ID = 2

contains

subroutine AtmosVars_Init( this, mesh )
  implicit none
  class(AtmosVars), intent(inout) :: this
  class(AtmosMesh), intent(in) :: mesh

  integer :: n
  type(VariableInfo) :: VINFO_PRG(2)
  DATA VINFO_PRG / &
  VariableInfo( ATMOSVARS_V_ID, "V", "velocity",    "m/s", 1, "X", ""), &
  VariableInfo( ATMOSVARS_T_ID, "T", "temperature", "K"  , 1, "X", "")  /
  
  type(VariableInfo) :: VINFO_AUX(1)
  DATA VINFO_AUX / &
    VariableInfo( ATMOSVARS_P_ID, "P", "pressure",    "Pa", 1, "X", "") /

  type(VariableInfo) :: VINFO_FORCING(2)
  DATA VINFO_FORCING / &
    VariableInfo( ATMOSVARS_F_V_ID, "F_V", "forcing_V", "m/s2", 1, "X", ""), &
    VariableInfo( ATMOSVARS_F_T_ID, "F_T", "forcing_T", "K/s",  1, "X", "") /

  !--------------------------------------------------

  LOG_INFO('AtmosVars_Init',*)

  call this%prgvars_list%Init()
  call this%prgvars_list%Regist( VINFO_PRG(ATMOSVARS_V_ID), mesh%mesh, this%V, .false. )
  call this%prgvars_list%Regist( VINFO_PRG(ATMOSVARS_T_ID), mesh%mesh, this%T, .false. )

  call this%auxvars_list%Init()
  call this%auxvars_list%Regist( VINFO_AUX(ATMOSVARS_P_ID), mesh%mesh, this%P, .false. )

  call this%forcing_list%Init()
  call this%forcing_list%Regist( VINFO_FORCING(ATMOSVARS_F_V_ID), mesh%mesh, this%F_T, .false. )
  call this%forcing_list%Regist( VINFO_FORCING(ATMOSVARS_F_T_ID), mesh%mesh, this%F_V, .false. )
  
  return
end subroutine AtmosVars_Init

subroutine AtmosVars_Final( this )
  implicit none
  class(AtmosVars), intent(inout) :: this

  class(ModelVarManager), pointer :: ptr_base
  !--------------------------------------------------

  LOG_INFO('AtmosVars_Final',*)

  call this%prgvars_list%Final()
  call this%auxvars_list%Final()
  call this%forcing_list%Final()

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