!-------------------------------------------------------------------------------
!> module Ocean / Dynamics
!!
!! @par Description
!!          Container for variables with dynamics component
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_ocean_dyn_vars
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
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D  
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfield_base, only: MeshField3D, MeshField2D
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use scale_meshfieldcomm_base, only: MeshFieldContainer
  
  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo
  use scale_model_meshbase_manager, only: ModelMeshBase

  use mod_ocean_mesh, only: OceanMesh
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  !> Derived type to manage variables with oceanic dynamics component
  type, public :: OceanDynVars
    type(MeshField2D), allocatable :: AUX_VARS2D(:) !< Array of 2D auxiliary variables
    type(ModelVarManager) :: AUXVARS2D_manager      !< Object to manage 2D auxiliary variables
  contains
    procedure :: Init => OceanDynVars_Init
    procedure :: Final => OceanDynVars_Final
    procedure :: History => OceanDynVars_History
  end type OceanDynVars

  public :: OceanDynAuxVars_GetLocalMeshFields
  !public :: OceanDynVars_GetLocalMeshFields_analysis

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-
  integer, public, parameter :: OCEAN_DYN_AUXVARS2D_NUM          = 1
  integer, public, parameter :: OCEAN_DYN_AUXVARS2D_CORIOLIS_ID  = 1

  type(VariableInfo), public :: OCEAN_DYN_AUXVARS2D_VINFO(OCEAN_DYN_AUXVARS2D_NUM)
  DATA OCEAN_DYN_AUXVARS2D_VINFO / &
    VariableInfo( OCEAN_DYN_AUXVARS2D_CORIOLIS_ID, 'CORIOLIS', 'coriolis parameter',  &
                  's-1',  2, 'XY',  ''                                             )  /   
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains

!> Setup an object to manage variables with oceanic dynamics component
!!
!! @param model_mesh Object to manage computational mesh of oceanic model 
!!
!OCL SERIAL
  subroutine OceanDynVars_Init( this, model_mesh )
    use scale_localmeshfield_base, only: LOCAL_MESHFIELD_TYPE_NODES_FACEVAL
    implicit none
    class(OceanDynVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    integer :: v
    integer :: n
    logical :: reg_file_hist

    class(OceanMesh), pointer :: ocean_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D

    !--------------------------------------------------

    LOG_INFO('OceanDynVars_Init',*)
    
    nullify( ocean_mesh )
    select type(model_mesh)
    class is (OceanMesh)
      ocean_mesh => model_mesh
    end select
    mesh3D => ocean_mesh%ptr_mesh
    
    call mesh3D%GetMesh2D( mesh2D )

    !- Initialize 2D auxiliary variables

    call this%AUXVARS2D_manager%Init()
    allocate( this%AUX_VARS2D(OCEAN_DYN_AUXVARS2D_NUM) )

    reg_file_hist = .false.    
    do v = 1, OCEAN_DYN_AUXVARS2D_NUM
      call this%AUXVARS2D_manager%Regist(                   &
        OCEAN_DYN_AUXVARS2D_VINFO(v), mesh2D,               & ! (in) 
        this%AUX_VARS2D(v), reg_file_hist                   ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%AUX_VARS2D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    return
  end subroutine OceanDynVars_Init

!> Finalize an object to manage variables with oceanic dynamics component
  subroutine OceanDynVars_Final( this )
    implicit none
    class(OceanDynVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('OceanDynVars_Final',*)

    call this%AUXVARS2D_manager%Final()
    deallocate( this%AUX_VARS2D )
    
    !call this%ANALYSISVARS_manager%Final()

    return
  end subroutine OceanDynVars_Final

  !> Output history files for variables with oceanic dynamics component
  subroutine OceanDynVars_History( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(OceanDynVars), intent(in) :: this
    !-------------------------------------------------------------------------
    return
  end subroutine OceanDynVars_History

  subroutine OceanDynAuxVars_GetLocalMeshFields( domID, mesh, auxvars_list, &
    Coriolis,                                                               &
    lcmesh3D                                                                &
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
    call auxvars_list%Get(OCEAN_DYN_AUXVARS2D_CORIOLIS_ID, field)
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
  end subroutine OceanDynAuxVars_GetLocalMeshFields
end module mod_ocean_dyn_vars