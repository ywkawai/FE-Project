!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cloud Microphysics
!!
!! @par Description
!!          Container for mod_atmos_phy_mp
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_mp_vars
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
  type, public :: AtmosPhyMpVars
    type(MeshField3D), allocatable :: tends(:)
    type(ModelVarManager) :: tends_manager

    integer :: QS
    integer :: QE
    integer :: QA

    integer :: TENDS_NUM_TOT 
  contains
    procedure :: Init => AtmosPhyMpVars_Init
    procedure :: Final => AtmosPhyMpVars_Final
  end type AtmosPhyMpVars

  public :: AtmosPhyMpVars_GetLocalMeshFields_tend
  public :: AtmosPhyMpVars_GetLocalMeshFields_tend_qtrc

  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !

  integer, public, parameter :: ATMOS_PHY_MP_DENS_t_ID  = 1  
  integer, public, parameter :: ATMOS_PHY_MP_MOMX_t_ID  = 2
  integer, public, parameter :: ATMOS_PHY_MP_MOMY_t_ID  = 3
  integer, public, parameter :: ATMOS_PHY_MP_MOMZ_t_ID  = 4
  integer, public, parameter :: ATMOS_PHY_MP_RHOT_t_ID  = 5
  integer, public, parameter :: ATMOS_PHY_MP_RHOH_ID    = 6  
  integer, public, parameter :: ATMOS_PHY_MP_TENDS_NUM1 = 6 

  type(VariableInfo), public :: ATMOS_PHY_MP_TEND_VINFO(ATMOS_PHY_MP_TENDS_NUM1)
  DATA ATMOS_PHY_MP_TEND_VINFO / &
    VariableInfo( ATMOS_PHY_MP_DENS_t_ID, 'MP_DENS_t', 'tendency of x-momentum in MP process',    &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
    VariableInfo( ATMOS_PHY_MP_MOMX_t_ID, 'MP_MOMX_t', 'tendency of x-momentum in MP process',    &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
    VariableInfo( ATMOS_PHY_MP_MOMY_t_ID, 'MP_MOMY_t', 'tendency of y-momentum in MP process',    &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
    VariableInfo( ATMOS_PHY_MP_MOMZ_t_ID, 'MP_MOMZ_t', 'tendency of z-momentum in MP process',    &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
    VariableInfo( ATMOS_PHY_MP_RHOT_t_ID, 'MP_RHOT_t', 'tendency of rho*PT in MP process',        &
                  'kg/m3.K/s', 3, 'XYZ',  ''                                                   ), &
    VariableInfo( ATMOS_PHY_MP_RHOH_ID  ,   'mp_RHOH', 'diabatic heating rate in MP process',     &
                  'J/kg/s',   3, 'XYZ',  ''                                                    )  / 

  type(VariableInfo), public, allocatable :: ATMOS_PHY_MP_TEND_VINFO_Q(:)                  


  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains
  subroutine AtmosPhyMpVars_Init( this, model_mesh, &
    QS_MP, QE_MP, QA_MP )
    implicit none
    class(AtmosPhyMpVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    integer, intent(in) :: QS_MP
    integer, intent(in) :: QE_MP
    integer, intent(in) :: QA_MP

    integer :: v
    integer :: n
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D
    !--------------------------------------------------

    LOG_INFO('AtmosPhyMpVars_Init',*)

    this%QS = QS_MP
    this%QE = QE_MP
    this%QA = QA_MP
    this%TENDS_NUM_TOT = ATMOS_PHY_MP_TENDS_NUM1 + QE_MP - QS_MP + 1

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
    allocate( this%tends(this%TENDS_NUM_TOT) )

    reg_file_hist = .false.    
    do v = 1, this%TENDS_NUM_TOT
      call this%tends_manager%Regist(           &
        ATMOS_PHY_MP_TEND_VINFO(v), mesh3D,     & ! (in) 
        this%tends(v), reg_file_hist            ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%tends(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    return
  end subroutine AtmosPhyMpVars_Init

  subroutine AtmosPhyMpVars_Final( this )
    implicit none
    class(AtmosPhyMpVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosPhyMpVars_Final',*)

    call this%tends_manager%Final()
    
    return
  end subroutine AtmosPhyMpVars_Final


  subroutine AtmosPhyMpVars_GetLocalMeshFields_tend( domID, mesh, mp_tends_list, &
    mp_DENS_t, mp_MOMX_t, mp_MOMY_t, mp_MOMZ_t, mp_RHOT_t, mp_RHOQ_t,            &
    lcmesh3D                                                                     &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: mp_tends_list
    class(LocalMeshFieldBase), pointer, intent(out) :: mp_DENS_t    
    class(LocalMeshFieldBase), pointer, intent(out) :: mp_MOMX_t
    class(LocalMeshFieldBase), pointer, intent(out) :: mp_MOMY_t
    class(LocalMeshFieldBase), pointer, intent(out) :: mp_MOMZ_t
    class(LocalMeshFieldBase), pointer, intent(out) :: mp_RHOT_t
    class(LocalMeshFieldBase), pointer, intent(out) :: mp_RHOQ_t
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call mp_tends_list%Get(ATMOS_PHY_MP_DENS_t_ID, field)
    call field%GetLocalMeshField(domID, mp_DENS_t)

    call mp_tends_list%Get(ATMOS_PHY_MP_MOMX_t_ID, field)
    call field%GetLocalMeshField(domID, mp_MOMX_t)

    call mp_tends_list%Get(ATMOS_PHY_MP_MOMY_t_ID, field)
    call field%GetLocalMeshField(domID, mp_MOMY_t)

    call mp_tends_list%Get(ATMOS_PHY_MP_MOMZ_t_ID, field)
    call field%GetLocalMeshField(domID, mp_MOMZ_t)

    call mp_tends_list%Get(ATMOS_PHY_MP_RHOT_t_ID, field)
    call field%GetLocalMeshField(domID, mp_RHOT_t)

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
  end subroutine AtmosPhyMpVars_GetLocalMeshFields_tend

  subroutine AtmosPhyMpVars_GetLocalMeshFields_tend_qtrc( domID, mesh, mp_tends_list, &
    QTRC_ID, QS_MP, mp_QTRC_t                                                         )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: mp_tends_list
    integer, intent(in) :: QTRC_ID
    integer, intent(in) :: QS_MP
    class(LocalMeshFieldBase), pointer, intent(out) :: mp_QTRC_t

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    call mp_tends_list%Get(ATMOS_PHY_MP_RHOT_t_ID + QTRC_ID - QS_MP + 1, field)
    call field%GetLocalMeshField(domID, mp_QTRC_t)

    return
  end subroutine AtmosPhyMpVars_GetLocalMeshFields_tend_qtrc

end  module mod_atmos_phy_mp_vars
