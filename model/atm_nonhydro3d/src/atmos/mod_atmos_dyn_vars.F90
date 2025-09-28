!-------------------------------------------------------------------------------
!> module ATMOSPHERE dynamics
!!
!! @par Description
!!          Container for variables with dynamics component
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
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
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D  
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfield_base, only: MeshField3D, MeshField2D
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use scale_meshfieldcomm_base, only: MeshFieldContainer
  
  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo
  use scale_model_meshbase_manager, only: ModelMeshBase

  use mod_atmos_mesh, only: AtmosMesh
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  !> Derived type to manage variables with atmospheric dynamics component
  type, public :: AtmosDynVars
    type(MeshField2D), allocatable :: AUX_VARS2D(:) !< Array of 2D auxiliary variables
    type(ModelVarManager) :: AUXVARS2D_manager      !< Object to manage 2D auxiliary variables
  contains
    procedure :: Init => AtmosDynVars_Init
    procedure :: Final => AtmosDynVars_Final
    procedure :: History => AtmosDynVars_History
  end type AtmosDynVars

  public :: AtmosDynAuxVars_GetLocalMeshFields
  !public :: AtmosDynVars_GetLocalMeshFields_analysis

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-
  integer, public, parameter :: ATMOS_DYN_AUXVARS2D_NUM          = 1
  integer, public, parameter :: ATMOS_DYN_AUXVARS2D_CORIOLIS_ID  = 1

  type(VariableInfo), public :: ATMOS_DYN_AUXVARS2D_VINFO(ATMOS_DYN_AUXVARS2D_NUM)
  DATA ATMOS_DYN_AUXVARS2D_VINFO / &
    VariableInfo( ATMOS_DYN_AUXVARS2D_CORIOLIS_ID, 'CORIOLIS', 'coriolis parameter',  &
                  's-1',  2, 'XY',  ''                                             )  / 


  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVARS_NUM          = 6
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t        = 1
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t_advX   = 2
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t_advY   = 3
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t_advZ   = 4
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t_lift   = 5
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t_buoy   = 6

  ! type(VariableInfo), public :: ATMOS_DYN_ANALYSISVAR_VINFO(ATMOS_DYN_ANALYSISVARS_NUM)
  ! DATA ATMOS_DYN_ANALYSISVAR_VINFO / &
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t, 'MOMZ_t', 'MOMZ_t',    &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ), & 
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t_advX, 'MOMZ_t_advx', 'MOMZ_t_advX',  &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ), & 
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t_advY, 'MOMZ_t_advy', 'MOMZ_t_advY',  &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ), & 
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t_advZ, 'MOMZ_t_advz', 'MOMZ_t_advZ',  &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ), & 
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t_lift, 'MOMZ_t_lift', 'MOMZ_t_lift',  &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ), & 
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t_buoy, 'MOMZ_t_buoy', 'MOMZ_t_buo',    &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ) & 
  ! / 
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains

!> Setup an object to manage variables with atmospheric dynamics component
!OCL SERIAL
  subroutine AtmosDynVars_Init( this, model_mesh )
    use scale_localmeshfield_base, only: LOCAL_MESHFIELD_TYPE_NODES_FACEVAL
    implicit none
    class(AtmosDynVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    integer :: v
    integer :: n
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D

    !--------------------------------------------------

    LOG_INFO('AtmosDynVars_Init',*)

    nullify( atm_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select
    mesh3D => atm_mesh%ptr_mesh
    
    call mesh3D%GetMesh2D( mesh2D )

    !- Initialize 2D auxiliary variables

    call this%AUXVARS2D_manager%Init()
    allocate( this%AUX_VARS2D(ATMOS_DYN_AUXVARS2D_NUM) )

    reg_file_hist = .false.    
    do v = 1, ATMOS_DYN_AUXVARS2D_NUM
      call this%AUXVARS2D_manager%Regist(                   &
        ATMOS_DYN_AUXVARS2D_VINFO(v), mesh2D,               & ! (in) 
        this%AUX_VARS2D(v), reg_file_hist                   ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%AUX_VARS2D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    !-
    ! call this%ANALYSISVARS_manager%Init()
    ! allocate( this%ANALYSIS_VARS3D(ATMOS_DYN_ANALYSISVARS_NUM) )

    ! reg_file_hist = .true.
    ! do v = 1, ATMOS_DYN_ANALYSISVARS_NUM
    !   call this%ANALYSISVARS_manager%Regist( &
    !     ATMOS_DYN_ANALYSISVAR_VINFO(v), atm_mesh%mesh, &
    !     this%ANALYSIS_VARS3D(v), reg_file_hist )
    ! end do

    return
  end subroutine AtmosDynVars_Init

!> Finalize an object to manage variables with atmospheric dynamics component
  subroutine AtmosDynVars_Final( this )
    implicit none
    class(AtmosDynVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosDynVars_Final',*)

    call this%AUXVARS2D_manager%Final()
    deallocate( this%AUX_VARS2D )
    
    !call this%ANALYSISVARS_manager%Final()

    return
  end subroutine AtmosDynVars_Final

  subroutine AtmosDynVars_History( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosDynVars), intent(in) :: this

    ! integer :: v
    ! integer :: hst_id
    !-------------------------------------------------------------------------

    ! do v = 1, ATMOS_DYN_ANALYSISVARS_NUM
    !   hst_id = this%ANALYSIS_VARS3D(v)%hist_id
    !   if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%ANALYSIS_VARS3D(v) )
    ! end do
    
    return
  end subroutine AtmosDynVars_History

  subroutine AtmosDynAuxVars_GetLocalMeshFields( domID, mesh, auxvars_list, &
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
  end subroutine AtmosDynAuxVars_GetLocalMeshFields

  ! subroutine AtmosDynVars_GetLocalMeshFields_analysis( domID, mesh, analysis_list, &
  !   MOMZ_t, MOMZ_t_advx, MOMZ_t_advY, MOMZ_t_advZ, MOMZ_t_lift, MOMZ_t_buoy,       &
  !   lcmesh3D                                                             &
  !   )

  !   use scale_mesh_base, only: MeshBase
  !   use scale_meshfield_base, only: MeshFieldBase
  !   implicit none

  !   integer, intent(in) :: domID
  !   class(MeshBase), intent(in) :: mesh
  !   class(ModelVarManager), intent(inout) :: analysis_list
  !   class(LocalMeshFieldBase), pointer, intent(out) :: &
  !     MOMZ_t, MOMZ_t_advx, MOMZ_t_advY, MOMZ_t_advZ, MOMZ_t_lift, MOMZ_t_buoy
  !   class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

  !   class(MeshFieldBase), pointer :: field   
  !   class(LocalMeshBase), pointer :: lcmesh
  !   !-------------------------------------------------------

  !   !--
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t)
  !   !---
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t_advx, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t_advx)
  !   !---
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t_advy, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t_advy)
  !   !---
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t_advz, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t_advz)
  !   !---
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t_lift, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t_lift)
  !   !---
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t_buoy, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t_buoy)

    
  !   if (present(lcmesh3D)) then
  !     call mesh%GetLocalMesh( domID, lcmesh )
  !     nullify( lcmesh3D )

  !     select type(lcmesh)
  !     type is (LocalMesh3D)
  !       if (present(lcmesh3D)) lcmesh3D => lcmesh
  !     end select
  !   end if

  !   return
  ! end subroutine AtmosDynVars_GetLocalMeshFields_analysis


end module mod_atmos_dyn_vars