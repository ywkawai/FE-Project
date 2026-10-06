!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / cloud macrophysics component
!!
!! @par Description
!!          Container for variables with cloud macrophysics component in atmospheric model
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_mac_vars
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

  !> Derived type to manage variables with cloud component in atmospheric model
  type, public :: AtmosPhyMacVars
    type(MeshField3D), allocatable :: tends(:)     !< Array of tendency variables with cloud component
    type(ModelVarManager) :: tends_manager         !< Object to manage tendencies with cloud component

    type(MeshField2D), allocatable :: auxvars2D(:) !< Array of 2D auxiliary variables with cloud component
    type(ModelVarManager) :: auxvars2D_manager     !< Object to manage 2D auxiliary variables with cloud component

    integer :: QS      !< Start index of tracer variables with cloud component
    integer :: QE      !< End index of tracer variables with cloud component
    integer :: QA      !< Number of tracer variables with cloud component

    integer :: TENDS_NUM_TOT                        !< Number of tendency variables with cloud component
  contains
    procedure :: Init => AtmosPhyMacVars_Init
    procedure :: Final => AtmosPhyMacVars_Final
    procedure :: Setup => AtmosPhyMacVars_Setup
    procedure :: History => AtmosPhyMacVars_history
  end type AtmosPhyMacVars

  public :: AtmosPhyMacVars_GetLocalMeshFields_tend
  public :: AtmosPhyMacVars_GetLocalMeshFields_sfcflx

  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  integer, public, parameter :: ATMOS_PHY_MAC_DENS_t_ID    = 1  
  integer, public, parameter :: ATMOS_PHY_MAC_MOMX_t_ID    = 2
  integer, public, parameter :: ATMOS_PHY_MAC_MOMY_t_ID    = 3
  integer, public, parameter :: ATMOS_PHY_MAC_MOMZ_t_ID    = 4
  integer, public, parameter :: ATMOS_PHY_MAC_RHOT_t_ID    = 5
  integer, public, parameter :: ATMOS_PHY_MAC_RHOH_ID      = 6 
  integer, public, parameter :: ATMOS_PHY_MAC_EVAPORATE_ID = 7
  integer, public, parameter :: ATMOS_PHY_MAC_TENDS_NUM1   = 7 

  type(VariableInfo), public :: ATMOS_PHY_MAC_TEND_VINFO(ATMOS_PHY_MAC_TENDS_NUM1)
  DATA ATMOS_PHY_MAC_TEND_VINFO / &
    VariableInfo( ATMOS_PHY_MAC_DENS_t_ID, 'MAC_DENS_t', 'tendency of x-momentum in CLD process',           &
                   'kg/m3/s',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_MAC_MOMX_t_ID, 'MAC_MOMX_t', 'tendency of x-momentum in CLD process',           &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_MAC_MOMY_t_ID, 'MAC_MOMY_t', 'tendency of y-momentum in CLD process',           &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_MAC_MOMZ_t_ID, 'MAC_MOMZ_t', 'tendency of z-momentum in CLD process',           &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_MAC_RHOT_t_ID, 'MAC_RHOT_t', 'tendency of rho*PT in CLD process',               &
                  'kg/m3.K/s', 3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_MAC_RHOH_ID  , 'MAC_RHOH', 'diabatic heating rate in CLD process',              &
                  'J/kg/s',   3, 'XYZ',  ''                                                           ), & 
    VariableInfo( ATMOS_PHY_MAC_EVAPORATE_ID, 'MAC_EVAPORATE', 'number concentration of evaporated cloud', &
                  'm-3'   ,   3, 'XYZ',  ''                                                           )  /                   

  integer, public, parameter :: ATMOS_PHY_MAC_AUX2D_SFLX_RAIN_ID   = 1
  integer, public, parameter :: ATMOS_PHY_MAC_AUX2D_SFLX_SNOW_ID   = 2
  integer, public, parameter :: ATMOS_PHY_MAC_AUX2D_SFLX_ENGI_ID   = 3
  integer, public, parameter :: ATMOS_PHY_MAC_AUX2D_NUM            = 3

  type(VariableInfo), public :: ATMOS_PHY_MAC_AUX2D_VINFO(ATMOS_PHY_MAC_AUX2D_NUM)
  DATA ATMOS_PHY_MAC_AUX2D_VINFO / &
    VariableInfo( ATMOS_PHY_MAC_AUX2D_SFLX_RAIN_ID, 'MAC_SFLX_RAIN', 'precipitation flux (liquid) in cld process',    &
                  'kg/m2/s',  2, 'XY',  ''                                                                      ), &
    VariableInfo( ATMOS_PHY_MAC_AUX2D_SFLX_SNOW_ID, 'MAC_SFLX_SNOW', 'precipitation flux (solid) in cld process',     &
                  'kg/m2/s',  2, 'XY',  ''                                                                      ), &
    VariableInfo( ATMOS_PHY_MAC_AUX2D_SFLX_ENGI_ID, 'MAC_SFLX_ENGI', 'internal energy flux flux in cld process',      &
                  'J/m2/s',   2, 'XY',  ''                                                                      )  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
contains
  !> Setup an object to manage variables with a cloud component  
!OCL SERIAL
  subroutine AtmosPhyMacVars_Init( this, model_mesh, &
    QS_mac, QE_mac, QA_mac )
    implicit none
    class(AtmosPhyMacVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    integer, intent(in) :: QS_mac
    integer, intent(in) :: QE_mac
    integer, intent(in) :: QA_mac
    !---------------------------------------------------

    LOG_INFO('AtmosPhyMacVars_Init',*)

    this%QS = QS_mac
    this%QE = QE_mac
    this%QA = QA_mac
    this%TENDS_NUM_TOT = ATMOS_PHY_MAC_TENDS_NUM1 + QE_mac - QS_mac + 1
    return
  end subroutine AtmosPhyMacVars_Init

  !> Setup variable objects with cloud macrophysics component in atmospheric model
  subroutine AtmosPhyMacVars_Setup( this, model_mesh )
    use scale_tracer, only: &
      TRACER_NAME, TRACER_DESC, TRACER_UNIT
    use scale_file_history, only: &
      FILE_HISTORY_reg
    implicit none
    class(AtmosPhyMacVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    integer :: iv
    integer :: iq
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D

    type(VariableInfo) :: qtrc_tp_vinfo_tmp
    type(VariableInfo) :: qtrc_vterm_vinfo_tmp
    !----------------------------------------------------


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

    reg_file_hist = .true.    
    do iv = 1, ATMOS_PHY_MAC_TENDS_NUM1
      call this%tends_manager%Regist( &
        ATMOS_PHY_MAC_TEND_VINFO(iv), mesh3D, &
        this%tends(iv), reg_file_hist,        &
        fill_zero=.true.                      )
    end do

    qtrc_tp_vinfo_tmp%ndims    = 3
    qtrc_tp_vinfo_tmp%dim_type = 'XYZ'
    qtrc_tp_vinfo_tmp%STDNAME  = ''
    
    do iq = 1, this%QA
      iv = ATMOS_PHY_MAC_TENDS_NUM1 + iq 
      qtrc_tp_vinfo_tmp%keyID = iv
      qtrc_tp_vinfo_tmp%NAME  = 'MAC_'//trim(TRACER_NAME(this%QS+iq-1))//'_t'
      qtrc_tp_vinfo_tmp%DESC  = 'tendency of rho*'//trim(TRACER_NAME(this%QS+iq-1))//' in cloud macrophysics process'
      qtrc_tp_vinfo_tmp%UNIT  = 'kg/m3/s'

      reg_file_hist = .true.
      call this%tends_manager%Regist( &
        qtrc_tp_vinfo_tmp, mesh3D,     & 
        this%tends(iv), reg_file_hist, &
        fill_zero=.true.               ) 
    end do    

    !--
    
    call this%auxvars2D_manager%Init()
    allocate( this%auxvars2D(ATMOS_PHY_MAC_AUX2D_NUM) )

    reg_file_hist = .true.    
    do iv = 1, ATMOS_PHY_MAC_AUX2D_NUM
      call this%auxvars2D_manager%Regist( &
        ATMOS_PHY_MAC_AUX2D_VINFO(iv), mesh2D, &
        this%auxvars2D(iv), reg_file_hist,     &
        fill_zero=.true.                       ) 
    end do

    return
  end subroutine AtmosPhyMacVars_Setup

  !> Finalize an object to manage variables with cloud macrophysics component in atmospheric model
!OCL SERIAL
  subroutine AtmosPhyMacVars_Final( this )
    implicit none
    class(AtmosPhyMacVars), intent(inout) :: this
    !----------------------------------------------------

     LOG_INFO('AtmosPhyMacVars_Final',*)

     call this%tends_manager%Final()
     deallocate( this%tends )

     call this%auxvars2D_manager%Final()
     deallocate( this%auxvars2D )

    return
  end subroutine AtmosPhyMacVars_Final

!OCL SERIAL
  subroutine AtmosPhyMacVars_GetLocalMeshFields_tend( domID, mesh, MAC_tends_list,   &
    MAC_DENS_t, MAC_MOMX_t, MAC_MOMY_t, MAC_MOMZ_t, MAC_RHOT_t, MAC_RHOH, MAC_EVAP,  &
    MAC_RHOQ_t,                                                                      &
    lcmesh3D                                                                         &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: MAC_tends_list
    class(LocalMeshFieldBase), pointer, intent(out) :: MAC_DENS_t    
    class(LocalMeshFieldBase), pointer, intent(out) :: MAC_MOMX_t
    class(LocalMeshFieldBase), pointer, intent(out) :: MAC_MOMY_t
    class(LocalMeshFieldBase), pointer, intent(out) :: MAC_MOMZ_t
    class(LocalMeshFieldBase), pointer, intent(out) :: MAC_RHOT_t
    class(LocalMeshFieldBase), pointer, intent(out) :: MAC_RHOH
    class(LocalMeshFieldBase), pointer, intent(out) :: MAC_EVAP
    type(LocalMeshFieldBaseList), intent(out) :: MAC_RHOQ_t(:)
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh

    integer :: iq
    !-------------------------------------------------------

    !--
    call MAC_tends_list%Get(ATMOS_PHY_MAC_DENS_t_ID, field)
    call field%GetLocalMeshField(domID, MAC_DENS_t)

    call MAC_tends_list%Get(ATMOS_PHY_MAC_MOMX_t_ID, field)
    call field%GetLocalMeshField(domID, MAC_MOMX_t)

    call MAC_tends_list%Get(ATMOS_PHY_MAC_MOMY_t_ID, field)
    call field%GetLocalMeshField(domID, MAC_MOMY_t)

    call MAC_tends_list%Get(ATMOS_PHY_MAC_MOMZ_t_ID, field)
    call field%GetLocalMeshField(domID, MAC_MOMZ_t)

    call MAC_tends_list%Get(ATMOS_PHY_MAC_RHOT_t_ID, field)
    call field%GetLocalMeshField(domID, MAC_RHOT_t)

    call MAC_tends_list%Get(ATMOS_PHY_MAC_RHOH_ID, field)
    call field%GetLocalMeshField(domID, MAC_RHOH)

    call MAC_tends_list%Get(ATMOS_PHY_MAC_EVAPORATE_ID, field)
    call field%GetLocalMeshField(domID, MAC_EVAP)

    !---
    do iq = 1, size(MAC_RHOQ_t)
      call MAC_tends_list%Get(ATMOS_PHY_MAC_TENDS_NUM1 + iq, field)
      call field%GetLocalMeshField(domID, MAC_RHOQ_t(iq)%ptr)
    end do

    
    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosPhyMacVars_GetLocalMeshFields_tend

!OCL SERIAL
  subroutine AtmosPhyMacVars_GetLocalMeshFields_sfcflx( domID, mesh, sfcflx_list, &
    SFLX_rain, SFLX_snow, SFLX_engi                                              )
    
    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: sfcflx_list
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_rain
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_snow
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_engi

    class(MeshFieldBase), pointer :: field
    !-------------------------------------------------------

    call sfcflx_list%Get(ATMOS_PHY_MAC_AUX2D_SFLX_RAIN_ID, field)
    call field%GetLocalMeshField(domID, SFLX_rain)

    call sfcflx_list%Get(ATMOS_PHY_MAC_AUX2D_SFLX_SNOW_ID, field)
    call field%GetLocalMeshField(domID, SFLX_snow)

    call sfcflx_list%Get(ATMOS_PHY_MAC_AUX2D_SFLX_engi_ID, field)
    call field%GetLocalMeshField(domID, SFLX_engi)
    
    return
  end subroutine AtmosPhyMacVars_GetLocalMeshFields_sfcflx

!OCL SERIAL
  subroutine AtmosPhyMacVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosPhyMacVars), intent(inout) :: this

    integer :: v
    integer :: hst_id
    !----------------------------------------------------

    do v=1, this%TENDS_NUM_TOT
      hst_id = this%tends(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%tends(v) )
    end do

    do v=1, ATMOS_PHY_MAC_AUX2D_NUM
      hst_id = this%auxvars2D(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%auxvars2D(v) )
    end do
    return
  end subroutine AtmosPhyMacVars_history

end module mod_atmos_phy_mac_vars