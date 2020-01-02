!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_dyn_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_element_base, only: ElementBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base2d, only: MeshBase2D  
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfield_base, only: MeshField3D, MeshField2D
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
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
  type, public :: AtmosDynVars
    type(MeshField2D), allocatable :: AUX_VARS2D(:)
    type(ModelVarManager) :: AUXVARS_manager
  contains
    procedure :: Init => AtmosDynVars_Init
    procedure :: Final => AtmosDynVars_Final
    procedure :: History => AtmosDynVars_History
  end type AtmosDynVars

  public AtmosDynVars_GetLocalMeshFields

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, parameter :: ATMOS_DYN_AUXVARS2D_NUM          = 1
  integer, public, parameter :: ATMOS_DYN_AUXVARS2D_CORIOLIS_ID  = 1

  type(VariableInfo), public :: ATMOS_DYN_AUXVARS2D_VINFO(ATMOS_DYN_AUXVARS2D_NUM)
  DATA ATMOS_DYN_AUXVARS2D_VINFO / &
    VariableInfo( ATMOS_DYN_AUXVARS2D_CORIOLIS_ID, 'CORIOLIS', 'coriolis parameter',  &
                  's-1',  2, 'XY',  ''                                             )  / 
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains
  subroutine AtmosDynVars_Init( this, model_mesh )
    implicit none
    class(AtmosDynVars), target, intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh

    integer :: v
    integer :: n
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    !--------------------------------------------------

    LOG_INFO('AtmosDynVars_Init',*)

    !- Initialize auxiliary and diagnostic variables

    nullify( atm_mesh )
    select type(model_mesh)
    type is (AtmosMesh)
      atm_mesh => model_mesh
    end select
    

    call this%AUXVARS_manager%Init()
    allocate( this%AUX_VARS2D(ATMOS_DYN_AUXVARS2D_NUM) )

    call atm_mesh%mesh%GetMesh2D( mesh2D )

    reg_file_hist = .false.    
    do v = 1, ATMOS_DYN_AUXVARS2D_NUM
      call this%AUXVARS_manager%Regist(                     &
        ATMOS_DYN_AUXVARS2D_VINFO(v), mesh2D,               & ! (in) 
        this%AUX_VARS2D(v), reg_file_hist                   ) ! (out)
      
      do n = 1, atm_mesh%mesh%LOCAL_MESH_NUM
        this%AUX_VARS2D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    return
  end subroutine AtmosDynVars_Init

  subroutine AtmosDynVars_Final( this )
    implicit none
    class(AtmosDynVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosDynVars_Final',*)

    call this%AUXVARS_manager%Final()

    return
  end subroutine AtmosDynVars_Final

  subroutine AtmosDynVars_History( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosDynVars), intent(in) :: this
    !-------------------------------------------------------------------------

    return
  end subroutine AtmosDynVars_History

  subroutine AtmosDynVars_GetLocalMeshFields( domID, mesh, auxvars_list, &
    Coriolis,                                                            &
    lcmesh3D                                                             &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: Coriolis
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call auxvars_list%Get(ATMOS_DYN_AUXVARS2D_CORIOLIS_ID, field)
    call field%GetLocalMeshField(domID, Coriolis)
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
  end subroutine AtmosDynVars_GetLocalMeshFields

end module mod_atmos_dyn_vars