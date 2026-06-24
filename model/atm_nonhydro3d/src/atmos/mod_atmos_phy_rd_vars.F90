!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Radiation
!!
!! @par Description
!!          Container for variables with radiation component in atmospheric model
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_rd_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_element_base, only: ElementBase3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: &
    MeshBase3D,                              &
    DIMTYPE_XYZ  => MeshBase3D_DIMTYPEID_XYZ
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField2D, MeshField3D

  use scale_file_restart_meshfield, only: &
    FILE_restart_meshfield_component
  
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  
  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo
  use scale_model_mesh_manager, only: ModelMeshBase
    
  use mod_atmos_mesh, only: AtmosMesh

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  !> Derived type to manage variables with radiation component in atmospheric model
  type, public :: AtmosPhyRdVars
    type(MeshField3D), allocatable :: tends(:)     !< Array of tendency variables
    type(ModelVarManager) :: tends_manager         !< Object to manage tendencies


    integer :: TENDS_NUM_TOT                        !< Number of tendency variables with cloud microphysics
  contains
    procedure :: Init => AtmosPhyRdVars_Init
    procedure :: Final => AtmosPhyRdVars_Final
    procedure :: History => AtmosPhyRdVars_history
  end type AtmosPhyRdVars
  
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  integer, public, parameter :: ATMOS_PHY_RD_RHOH_ID      = 1
  integer, public, parameter :: ATMOS_PHY_RD_TENDS_NUM1   = 1

  type(VariableInfo), public :: ATMOS_PHY_RD_TEND_VINFO(ATMOS_PHY_RD_TENDS_NUM1)
  DATA ATMOS_PHY_RD_TEND_VINFO / &
    VariableInfo( ATMOS_PHY_RD_RHOH_ID, 'RD_RHOH', 'diabatic heating rate in RD process',              &
                  'J/kg/s',   3, 'XYZ',  ''                                                           )  /                   

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
contains
  !> Setup an object to manage variables with a radiation component  
!OCL SERIAL
  subroutine AtmosPhyRdVars_Init( this, model_mesh )

    use scale_tracer, only: &
      TRACER_NAME, TRACER_DESC, TRACER_UNIT
    use scale_file_history, only: &
      FILE_HISTORY_reg

    implicit none
    class(AtmosPhyRdVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D

    !----------------------------------------------------

    LOG_INFO('AtmosPhyRdVars_Init',*)

    this%TENDS_NUM_TOT = ATMOS_PHY_RD_TENDS_NUM1

    !- Initialize auxiliary and diagnostic variables

    nullify( atm_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select
    mesh3D => atm_mesh%ptr_mesh

    call mesh3D%GetMesh2D( mesh2D )

    !----

    call this%tends_manager%Init()
    allocate( this%tends(this%TENDS_NUM_TOT) )

    return
  end subroutine AtmosPhyRdVars_Init

  !> Finalize an object to manage variables with radiation component  
!OCL SERIAL
  subroutine AtmosPhyRdVars_Final( this )
    implicit none
    class(AtmosPhyRdVars), intent(inout) :: this
    !----------------------------------------------------

    LOG_INFO('AtmosPhyMpVars_Final',*)

    call this%tends_manager%Final()
    deallocate( this%tends )

    return
  end subroutine AtmosPhyRdVars_Final

!OCL SERIAL
  subroutine AtmosPhyRdVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosPhyRdVars), intent(inout) :: this
    !----------------------------------------------------
    return
  end subroutine AtmosPhyRdVars_history

end module mod_atmos_phy_rd_vars