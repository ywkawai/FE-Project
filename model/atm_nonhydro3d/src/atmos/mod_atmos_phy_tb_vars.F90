!-------------------------------------------------------------------------------
!> module Atmosphere / Physics turbulence
!!
!! @par Description
!!          Container for variables with turbulence model
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_tb_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_tracer, only: QA

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
    
  use scale_atm_phy_tb_dgm_common, only: &
    ATMOS_PHY_TB_AUX_NUM, ATMOS_PHY_TB_AUX_SCALAR_NUM, ATMOS_PHY_TB_AUX_HVEC_NUM, &
    ATMOS_PHY_TB_AUX_HTENSOR_NUM,                                                    &
    ATMOS_PHY_TB_DIAG_NUM,                                                           &
    ATMOS_PHY_TB_TENDS_NUM1, &
    TB_MOMX_t_VID => ATMOS_PHY_TB_MOMX_t_ID,TB_MOMY_t_VID => ATMOS_PHY_TB_MOMY_t_ID, &
    TB_MOMZ_t_VID => ATMOS_PHY_TB_MOMZ_t_ID, TB_RHOT_t_VID => ATMOS_PHY_TB_RHOT_t_ID

  use mod_atmos_mesh, only: AtmosMesh

  !-----------------------------------------------------------------------------
  implicit none
  private


  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: AtmosPhyTbVars
    type(MeshField3D), allocatable :: tends(:)
    type(ModelVarManager) :: tends_manager

    type(MeshField3D), allocatable :: auxvars(:)
    type(ModelVarManager) :: auxvars_manager
    integer :: auxvars_commid

    type(MeshField3D), allocatable :: diagvars(:)
    type(ModelVarManager) :: diagvars_manager

    integer :: TENDS_NUM_TOT 
  contains
    procedure :: Init => AtmosPhyTbVars_Init
    procedure :: Final => AtmosPhyTbVars_Final
    procedure :: History => AtmosPhyTbVars_history
  end type AtmosPhyTbVars

  public :: AtmosPhyTbVars_GetLocalMeshFields_tend

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains
!OCL SERIAL
  subroutine AtmosPhyTbVars_Init( this, model_mesh )

    use scale_tracer, only: &
      TRACER_NAME, TRACER_DESC, TRACER_UNIT    
    use scale_atm_phy_tb_dgm_common, only: &
      atm_phy_tb_dgm_common_setup_variables
    
    implicit none
    class(AtmosPhyTbVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    integer :: iv
    integer :: iq
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D
    !--------------------------------------------------

    LOG_INFO('AtmosPhyTbVars_Init',*)

    this%TENDS_NUM_TOT = ATMOS_PHY_TB_TENDS_NUM1 + QA

    !- Initialize auxiliary and diagnostic variables

    nullify( atm_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select
    mesh3D => atm_mesh%ptr_mesh
    
    call mesh3D%GetMesh2D( mesh2D )

    !- Get variable information
    ! call atm_phy_tb_dgm_common_get_varinfo( auxvar_info, diagvar_info, tbtend_info ) ! (out)

    !- Initialize objects to variables with turbulent model

    call this%tends_manager%Init()
    call this%auxvars_manager%Init()
    call this%diagvars_manager%Init()

    allocate( this%diagvars(ATMOS_PHY_TB_DIAG_NUM) )
    allocate( this%tends(this%TENDS_NUM_TOT) )
    allocate( this%auxvars(ATMOS_PHY_TB_AUX_NUM) )

    call atm_phy_tb_dgm_common_setup_variables( &
      this%tends, this%auxvars, this%diagvars,                         & ! (inout)
      this%tends_manager, this%auxvars_manager, this%diagvars_manager, & ! (inout)
      this%TENDS_NUM_TOT, mesh3D )                                       ! (in)

    !- Setup communication

    call atm_mesh%Create_communicator( &
      ATMOS_PHY_TB_AUX_SCALAR_NUM, ATMOS_PHY_TB_AUX_HVEC_NUM, & ! (in)
      ATMOS_PHY_TB_AUX_HTENSOR_NUM,                           & ! (in)
      this%auxvars_manager,            & ! (inout)
      this%auxvars(:),                 & ! (in)
      this%auxvars_commid              ) ! (out)

    return
  end subroutine AtmosPhyTbVars_Init

!OCL SERIAL
  subroutine AtmosPhyTbVars_Final( this )
    implicit none
    class(AtmosPhyTbVars), intent(inout) :: this
    !--------------------------------------------------

    LOG_INFO('AtmosPhyTbVars_Final',*)

    call this%tends_manager%Final()
    deallocate( this%tends )

    call this%auxvars_manager%Final()
    deallocate( this%auxvars )

    call this%diagvars_manager%Final()
    deallocate( this%diagvars )

    return
  end subroutine AtmosPhyTbVars_Final

!OCL SERIAL
  subroutine AtmosPhyTbVars_GetLocalMeshFields_tend( domID, mesh, tb_tends_list, &
    tb_MOMX_t, tb_MOMY_t, tb_MOMZ_t, tb_RHOT_t, tb_RHOQ_t,                       &
    lcmesh3D                                                                     &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: tb_tends_list
    class(LocalMeshFieldBase), pointer, intent(out) :: tb_MOMX_t
    class(LocalMeshFieldBase), pointer, intent(out) :: tb_MOMY_t
    class(LocalMeshFieldBase), pointer, intent(out) :: tb_MOMZ_t
    class(LocalMeshFieldBase), pointer, intent(out) :: tb_RHOT_t
    type(LocalMeshFieldBaseList), intent(out), optional :: tb_RHOQ_t(:)
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh

    integer :: iq
    !-------------------------------------------------------

    !--
    call tb_tends_list%Get(TB_MOMX_t_VID, field)
    call field%GetLocalMeshField(domID, tb_MOMX_t)

    call tb_tends_list%Get(TB_MOMY_t_VID, field)
    call field%GetLocalMeshField(domID, tb_MOMY_t)

    call tb_tends_list%Get(TB_MOMZ_t_VID, field)
    call field%GetLocalMeshField(domID, tb_MOMZ_t)

    call tb_tends_list%Get(TB_RHOT_t_VID, field)
    call field%GetLocalMeshField(domID, tb_RHOT_t)

    !---
    if ( present(tb_RHOQ_t) ) then
      do iq = 1, size(tb_RHOQ_t)
        call tb_tends_list%Get(ATMOS_PHY_TB_TENDS_NUM1 + iq, field)
        call field%GetLocalMeshField(domID, tb_RHOQ_t(iq)%ptr)
      end do    
    end if

    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosPhyTbVars_GetLocalMeshFields_tend

!OCL SERIAL
  subroutine AtmosPhyTbVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosPhyTbVars), intent(inout) :: this
  
    integer :: v
    integer :: hst_id
    type(MeshField3D) :: tmp_field
    class(MeshBase3D), pointer :: mesh3D
    !-------------------------------------------------------------------------

    mesh3D => this%auxvars(1)%mesh

    do v = 1, this%TENDS_NUM_TOT
      hst_id = this%tends(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%tends(v) )
    end do

    do v = 1, ATMOS_PHY_TB_AUX_NUM
      hst_id = this%auxvars(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%auxvars(v) )
    end do

    do v = 1, ATMOS_PHY_TB_DIAG_NUM
      hst_id = this%diagvars(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%diagvars(v) )
    end do

    return
  end subroutine AtmosPhyTbVars_history

end  module mod_atmos_phy_tb_vars
