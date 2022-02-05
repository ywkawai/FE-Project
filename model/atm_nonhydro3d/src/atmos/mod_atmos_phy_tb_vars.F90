!-------------------------------------------------------------------------------
!> module Atmosphere / Physics turbulence
!!
!! @par Description
!!          Container for mod_atmos_phy_mp
!!
!! @author Team SCALE
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

  use scale_element_base, only: ElementBase3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: &
    MeshBase3D,                              &
    DIMTYPE_XYZ  => MeshBase3D_DIMTYPEID_XYZ
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
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
  type, public :: AtmosPhyTbVars
    type(MeshField3D), allocatable :: tends(:)
    type(ModelVarManager) :: tends_manager

    type(MeshField3D), allocatable :: auxvars(:)
    type(ModelVarManager) :: auxvars_manager
    integer :: auxvars_commid

  contains
    procedure :: Init => AtmosPhyTbVars_Init
    procedure :: Final => AtmosPhyTbVars_Final
    procedure :: History => AtmosPhyTbVars_history
  end type AtmosPhyTbVars

  integer, public, parameter :: ATMOS_PHY_TB_MOMX_t_ID  = 1
  integer, public, parameter :: ATMOS_PHY_TB_MOMY_t_ID  = 2
  integer, public, parameter :: ATMOS_PHY_TB_MOMZ_t_ID  = 3
  integer, public, parameter :: ATMOS_PHY_TB_RHOT_t_ID  = 4
  integer, public, parameter :: ATMOS_PHY_TB_RHOQ_t_ID  = 5
  integer, public, parameter :: ATMOS_PHY_TB_TENDS_NUM  = 5

  type(VariableInfo), public :: ATMOS_PHY_TB_TEND_VINFO(ATMOS_PHY_TB_TENDS_NUM)
  DATA ATMOS_PHY_TB_TEND_VINFO / &
    VariableInfo( ATMOS_PHY_TB_MOMX_t_ID, 'TB_MOMX_t', 'tendency of x-momentum in TB process',    &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
    VariableInfo( ATMOS_PHY_TB_MOMY_t_ID, 'TB_MOMY_t', 'tendency of y-momentum in TB process',    &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
    VariableInfo( ATMOS_PHY_TB_MOMZ_t_ID, 'TB_MOMZ_t', 'tendency of z-momentum in TB process',    &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
    VariableInfo( ATMOS_PHY_TB_RHOT_t_ID, 'TB_RHOT_t', 'tendency of rho*PT in TB process',        &
                  'kg/m3.K/s', 3, 'XYZ',  ''                                                   ), &
    VariableInfo( ATMOS_PHY_TB_RHOQ_t_ID, 'TB_RHOQ_t', 'tendency of rho*QTRC in TB process',      &
                  'kg/m3/s',   3, 'XYZ',  ''                                                   )  / 


  integer, public, parameter :: ATMOS_PHY_TB_AUX_S11_ID   = 1
  integer, public, parameter :: ATMOS_PHY_TB_AUX_S12_ID   = 2
  integer, public, parameter :: ATMOS_PHY_TB_AUX_S22_ID   = 3
  integer, public, parameter :: ATMOS_PHY_TB_AUX_S23_ID   = 4
  integer, public, parameter :: ATMOS_PHY_TB_AUX_S31_ID   = 5
  integer, public, parameter :: ATMOS_PHY_TB_AUX_S33_ID   = 6
  integer, public, parameter :: ATMOS_PHY_TB_AUX_TKE_ID   = 7
  integer, public, parameter :: ATMOS_PHY_TB_AUX_DPTDX_ID = 8
  integer, public, parameter :: ATMOS_PHY_TB_AUX_DPTDY_ID = 9
  integer, public, parameter :: ATMOS_PHY_TB_AUX_DPTDZ_ID = 10
  integer, public, parameter :: ATMOS_PHY_TB_AUX_NU_ID    = 11
  integer, public, parameter :: ATMOS_PHY_TB_AUX_KH_ID    = 12
  integer, public, parameter :: ATMOS_PHY_TB_AUX_NUM      = 12
  type(VariableInfo), public :: ATMOS_PHY_TB_AUX_VINFO(ATMOS_PHY_TB_AUX_NUM)
  DATA ATMOS_PHY_TB_AUX_VINFO / &
    VariableInfo( ATMOS_PHY_TB_AUX_S11_ID, 'S11', 'strain velocity tensor (S11)',    &
                  's-1',  3, 'XYZ',  ''                                           ), &
    VariableInfo( ATMOS_PHY_TB_AUX_S12_ID, 'S12', 'strain velocity tensor (S12)',    &
                  's-1',  3, 'XYZ',  ''                                           ), &
    VariableInfo( ATMOS_PHY_TB_AUX_S22_ID, 'S22', 'strain velocity tensor (S22)',    &
                  's-1',  3, 'XYZ',  ''                                           ), &
    VariableInfo( ATMOS_PHY_TB_AUX_S23_ID, 'S23', 'strain velocity tensor (S23)',    &
                  's-1',  3, 'XYZ',  ''                                           ), &
    VariableInfo( ATMOS_PHY_TB_AUX_S31_ID, 'S31', 'strain velocity tensor (S31)',    &
                  's-1',  3, 'XYZ',  ''                                           ), &
    VariableInfo( ATMOS_PHY_TB_AUX_S33_ID, 'S33', 'strain velocity tensor (S33)',    &
                  's-1',  3, 'XYZ',  ''                                           ), &
    VariableInfo( ATMOS_PHY_TB_AUX_TKE_ID, 'TKE', 'SGS turbluence kinetic energy',   &
                  'm2/s2',  3, 'XYZ',  ''                                         ), &
    VariableInfo( ATMOS_PHY_TB_AUX_DPTDX_ID, 'DPTDX', 'gradient of PT (x)',          &
                  'K/m',  3, 'XYZ',  ''                                           ), &
    VariableInfo( ATMOS_PHY_TB_AUX_DPTDY_ID, 'DPTDY', 'gradient of PT (y)',          &
                  'K/m',  3, 'XYZ',  ''                                           ), &
    VariableInfo( ATMOS_PHY_TB_AUX_DPTDZ_ID, 'DPTDZ', 'gradient of PT (z)',          &
                  'K/m',  3, 'XYZ',  ''                                           ), &
    VariableInfo( ATMOS_PHY_TB_AUX_NU_ID, 'NU', 'eddy viscosity',                    &
                  'm2/s',  3, 'XYZ',  ''                                          ), &
    VariableInfo( ATMOS_PHY_TB_AUX_KH_ID, 'KH', 'eddy diffusion',                    &
                  'm2/s',  3, 'XYZ',  ''                                          )  /


  public :: AtmosPhyTbVars_GetLocalMeshFields_tend
  public :: AtmosPhyTbVars_GetLocalMeshFields_aux

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains
  subroutine AtmosPhyTbVars_Init( this, model_mesh )
    implicit none
    class(AtmosPhyTbVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    integer :: v
    integer :: n
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D
    !--------------------------------------------------

    LOG_INFO('AtmosPhyTbVars_Init',*)

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

    !-
    allocate( this%tends(ATMOS_PHY_TB_TENDS_NUM) )

    reg_file_hist = .true.    
    do v = 1, ATMOS_PHY_TB_TENDS_NUM
      call this%tends_manager%Regist(           &
        ATMOS_PHY_TB_TEND_VINFO(v), mesh3D,     & ! (in) 
        this%tends(v), reg_file_hist            ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%tends(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    !-
    allocate( this%auxvars(ATMOS_PHY_TB_AUX_NUM) )

    reg_file_hist = .true.    
    do v = 1, ATMOS_PHY_TB_AUX_NUM
      call this%auxvars_manager%Regist(         &
        ATMOS_PHY_TB_AUX_VINFO(v), mesh3D,      & ! (in) 
        this%auxvars(v), reg_file_hist          ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%auxvars(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    call atm_mesh%Create_communicator( &
      ATMOS_PHY_TB_AUX_NUM, 0,         & ! (in)
      this%auxvars_manager,            & ! (inout)
      this%auxvars(:),                 & ! (in)
      this%auxvars_commid              ) ! (out)

    return
  end subroutine AtmosPhyTbVars_Init

  subroutine AtmosPhyTbVars_Final( this )
    implicit none
    class(AtmosPhyTbVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosPhyTbVars_Final',*)

    call this%tends_manager%Final()
    call this%auxvars_manager%Final()

    return
  end subroutine AtmosPhyTbVars_Final


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
    class(LocalMeshFieldBase), pointer, intent(out) :: tb_RHOQ_t
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call tb_tends_list%Get(ATMOS_PHY_TB_MOMX_t_ID, field)
    call field%GetLocalMeshField(domID, tb_MOMX_t)

    call tb_tends_list%Get(ATMOS_PHY_TB_MOMY_t_ID, field)
    call field%GetLocalMeshField(domID, tb_MOMY_t)

    call tb_tends_list%Get(ATMOS_PHY_TB_MOMZ_t_ID, field)
    call field%GetLocalMeshField(domID, tb_MOMZ_t)

    call tb_tends_list%Get(ATMOS_PHY_TB_RHOT_t_ID, field)
    call field%GetLocalMeshField(domID, tb_RHOT_t)

    call tb_tends_list%Get(ATMOS_PHY_TB_RHOQ_t_ID, field)
    call field%GetLocalMeshField(domID, tb_RHOQ_t)
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
  end subroutine AtmosPhyTbVars_GetLocalMeshFields_tend

  subroutine AtmosPhyTbVars_GetLocalMeshFields_aux( domID, mesh, tb_aux_list, &
    S11, S12, S22, S23, S31, S33, TKE,                                        &
    dPTdx, dPTdy, dPTdz,                                                      &
    Nu, Kh,                                                                   &
    lcmesh3D                                                                  )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: tb_aux_list
    class(LocalMeshFieldBase), pointer, intent(out) :: S11, S12, S22, S23, S31, S33, TKE
    class(LocalMeshFieldBase), pointer, intent(out) :: dPTdx, DPTdy, DPTdz
    class(LocalMeshFieldBase), pointer, intent(out) :: Nu, Kh
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_S11_ID, field)
    call field%GetLocalMeshField(domID, S11)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_S12_ID, field)
    call field%GetLocalMeshField(domID, S12)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_S22_ID, field)
    call field%GetLocalMeshField(domID, S22)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_S23_ID, field)
    call field%GetLocalMeshField(domID, S23)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_S31_ID, field)
    call field%GetLocalMeshField(domID, S31)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_S33_ID, field)
    call field%GetLocalMeshField(domID, S33)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_TKE_ID, field)
    call field%GetLocalMeshField(domID, TKE)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_DPTDX_ID, field)
    call field%GetLocalMeshField(domID, dPTdx)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_DPTDY_ID, field)
    call field%GetLocalMeshField(domID, dPTdy)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_DPTDZ_ID, field)
    call field%GetLocalMeshField(domID, dPTdz)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_NU_ID, field)
    call field%GetLocalMeshField(domID, Nu)

    call tb_aux_list%Get(ATMOS_PHY_TB_AUX_KH_ID, field)
    call field%GetLocalMeshField(domID, Kh)

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
  end subroutine AtmosPhyTbVars_GetLocalMeshFields_aux

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

    do v = 1, ATMOS_PHY_TB_TENDS_NUM
      hst_id = this%tends(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%tends(v) )
    end do

    do v = 1, ATMOS_PHY_TB_AUX_NUM
      hst_id = this%auxvars(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%auxvars(v) )
    end do

    return
  end subroutine AtmosPhyTbVars_history

end  module mod_atmos_phy_tb_vars
