!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_sfc_vars
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
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField2D, MeshField3D

  use scale_file_restart_meshfield, only: &
    FILE_restart_meshfield_component
  use scale_file_common_meshfield, only: &
    DIMTYPE_XYZ  => FILE_COMMON_MESHFILED3D_DIMTYPEID_XYZ
  
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
  type, public :: AtmosPhySfcVars
    type(MeshField2D), allocatable :: SFC_FLX(:)
    type(ModelVarManager) :: SFCFLX_manager
  contains
    procedure :: Init => AtmosPhySfcVars_Init
    procedure :: Final => AtmosPhySfcVars_Final
  end type AtmosPhySfcVars

  integer, public, parameter :: ATMOS_PHY_SF_SFLX_MU_ID  = 1
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_MV_ID  = 2
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_MW_ID  = 3
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_SH_ID  = 4
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_LH_ID  = 5
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_NUM    = 5

  type(VariableInfo), public :: ATMOS_PHY_SF_SFLX_VINFO(ATMOS_PHY_SF_SFLX_NUM)
  DATA ATMOS_PHY_SF_SFLX_VINFO / &
    VariableInfo( ATMOS_PHY_SF_SFLX_MU_ID, 'SFLX_MU', 'x-momentum flux',    &
                  'm/s*kg/m2/s',  2, 'XY',  ''                           ), &
    VariableInfo( ATMOS_PHY_SF_SFLX_MV_ID, 'SFLX_MV', 'y-momentum flux',    &
                  'm/s*kg/m2/s',  2, 'XY',  ''                           ), &
    VariableInfo( ATMOS_PHY_SF_SFLX_MW_ID, 'SFLX_MW', 'z-momentum flux',    &
                  'm/s*kg/m2/s',  2, 'XY',  ''                           ), &
    VariableInfo( ATMOS_PHY_SF_SFLX_SH_ID, 'SFLX_SH', 'sensible heat flux', &
                  'J/m2/s',  2, 'XY',  ''                                ), &
    VariableInfo( ATMOS_PHY_SF_SFLX_LH_ID, 'SFLX_LH', 'latent heat flux',   &
                  'J/m2/s',  2, 'XY',  ''                                )  / 

  public :: AtmosPhySfcVars_GetLocalMeshFields

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains
  subroutine AtmosPhySfcVars_Init( this, model_mesh )
    implicit none
    class(AtmosPhySfcVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    integer :: v
    integer :: n
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D
    !--------------------------------------------------

    LOG_INFO('AtmosPhySfcVars_Init',*)

    !- Initialize auxiliary and diagnostic variables

    nullify( atm_mesh )
    select type(model_mesh)
    type is (AtmosMesh)
      atm_mesh => model_mesh
    end select
    
    mesh3D => atm_mesh%mesh
    call atm_mesh%mesh%GetMesh2D( mesh2D )

    !----
    call this%SFCFLX_manager%Init()
    allocate( this%SFC_FLX(ATMOS_PHY_SF_SFLX_NUM) )

    reg_file_hist = .false.    
    do v = 1, ATMOS_PHY_SF_SFLX_NUM
      call this%SFCFLX_manager%Regist(            &
        ATMOS_PHY_SF_SFLX_VINFO(v), mesh2D,       & ! (in) 
        this%SFC_FLX(v), reg_file_hist            ) ! (out)
      
      do n = 1, atm_mesh%mesh%LOCAL_MESH_NUM
        this%SFC_FLX(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    return
  end subroutine AtmosPhySfcVars_Init

  subroutine AtmosPhySfcVars_Final( this )
    implicit none
    class(AtmosPhySfcVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosPhySfcVars_Final',*)

    call this%SFCFLX_manager%Final()

    return
  end subroutine AtmosPhySfcVars_Final


  subroutine AtmosPhySfcVars_GetLocalMeshFields( domID, mesh, sflx_list, &
    SFLX_MU, SFLX_MV, SFLX_MW, SFLX_SH, SFLX_LH,                         &
    lcmesh3D                                                             &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: sflx_list
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_MU
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_MV
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_MW
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_SH
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_LH
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call sflx_list%Get(ATMOS_PHY_SF_SFLX_MU_ID, field)
    call field%GetLocalMeshField(domID, SFLX_MU)

    call sflx_list%Get(ATMOS_PHY_SF_SFLX_MV_ID, field)
    call field%GetLocalMeshField(domID, SFLX_MV)

    call sflx_list%Get(ATMOS_PHY_SF_SFLX_MW_ID, field)
    call field%GetLocalMeshField(domID, SFLX_MW)

    call sflx_list%Get(ATMOS_PHY_SF_SFLX_SH_ID, field)
    call field%GetLocalMeshField(domID, SFLX_SH)

    call sflx_list%Get(ATMOS_PHY_SF_SFLX_LH_ID, field)
    call field%GetLocalMeshField(domID, SFLX_LH)
    !---
    
    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosPhySfcVars_GetLocalMeshFields

end  module mod_atmos_phy_sfc_vars
