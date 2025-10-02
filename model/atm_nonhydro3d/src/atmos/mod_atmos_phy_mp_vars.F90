!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Cloud Microphysics
!!
!! @par Description
!!          Container for variables with cloud microphysics component
!!
!! @author Yuta Kawai, Team SCALE
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

  !> Derived type to manage variables with cloud microphysics component  
  type, public :: AtmosPhyMpVars
    type(MeshField3D), allocatable :: tends(:)     !< Array of tendency variables
    type(ModelVarManager) :: tends_manager         !< Object to manage tendencies

    type(MeshField2D), allocatable :: auxvars2D(:) !< Array of 2D auxiliary variables
    type(ModelVarManager) :: auxvars2D_manager     !< Object to manage 2D auxiliary variables

    integer :: QS      !< Start index of tracer variables with cloud microphysics
    integer :: QE      !< End index of tracer variables with cloud microphysics
    integer :: QA      !< Number of tracer variables with cloud microphysics

    integer :: TENDS_NUM_TOT                        !< Number of tendency variables with cloud microphysics
    integer, allocatable :: vterm_hist_id(:)
    type(MeshField3D), allocatable :: vterm_hist(:)
  contains
    procedure :: Init => AtmosPhyMpVars_Init
    procedure :: Final => AtmosPhyMpVars_Final
    procedure :: History => AtmosPhyMpVars_history
  end type AtmosPhyMpVars

  public :: AtmosPhyMpVars_GetLocalMeshFields_tend
  public :: AtmosPhyMpVars_GetLocalMeshFields_sfcflx

  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !

  integer, public, parameter :: ATMOS_PHY_MP_DENS_t_ID    = 1  
  integer, public, parameter :: ATMOS_PHY_MP_MOMX_t_ID    = 2
  integer, public, parameter :: ATMOS_PHY_MP_MOMY_t_ID    = 3
  integer, public, parameter :: ATMOS_PHY_MP_MOMZ_t_ID    = 4
  integer, public, parameter :: ATMOS_PHY_MP_RHOT_t_ID    = 5
  integer, public, parameter :: ATMOS_PHY_MP_RHOH_ID      = 6 
  integer, public, parameter :: ATMOS_PHY_MP_EVAPORATE_ID = 7
  integer, public, parameter :: ATMOS_PHY_MP_TENDS_NUM1   = 7 

  type(VariableInfo), public :: ATMOS_PHY_MP_TEND_VINFO(ATMOS_PHY_MP_TENDS_NUM1)
  DATA ATMOS_PHY_MP_TEND_VINFO / &
    VariableInfo( ATMOS_PHY_MP_DENS_t_ID, 'MP_DENS_t', 'tendency of x-momentum in MP process',           &
                   'kg/m3/s',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_MP_MOMX_t_ID, 'MP_MOMX_t', 'tendency of x-momentum in MP process',           &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_MP_MOMY_t_ID, 'MP_MOMY_t', 'tendency of y-momentum in MP process',           &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_MP_MOMZ_t_ID, 'MP_MOMZ_t', 'tendency of z-momentum in MP process',           &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_MP_RHOT_t_ID, 'MP_RHOT_t', 'tendency of rho*PT in MP process',               &
                  'kg/m3.K/s', 3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_MP_RHOH_ID  , 'mp_RHOH', 'diabatic heating rate in MP process',              &
                  'J/kg/s',   3, 'XYZ',  ''                                                           ), & 
    VariableInfo( ATMOS_PHY_MP_EVAPORATE_ID, 'mp_EVAPORATE', 'number concentration of evaporated cloud', &
                  'm-3'   ,   3, 'XYZ',  ''                                                           )  /                   


  integer, public, parameter :: ATMOS_PHY_MP_AUX2D_SFLX_RAIN_ID   = 1
  integer, public, parameter :: ATMOS_PHY_MP_AUX2D_SFLX_SNOW_ID   = 2
  integer, public, parameter :: ATMOS_PHY_MP_AUX2D_SFLX_ENGI_ID   = 3
  integer, public, parameter :: ATMOS_PHY_MP_AUX2D_NUM            = 3

  type(VariableInfo), public :: ATMOS_PHY_MP_AUX2D_VINFO(ATMOS_PHY_MP_AUX2D_NUM)
  DATA ATMOS_PHY_MP_AUX2D_VINFO / &
    VariableInfo( ATMOS_PHY_MP_AUX2D_SFLX_RAIN_ID, 'MP_SFLX_RAIN', 'precipitation flux (liquid) in MP process',    &
                  'kg/m2/s',  2, 'XY',  ''                                                                      ), &
    VariableInfo( ATMOS_PHY_MP_AUX2D_SFLX_SNOW_ID, 'MP_SFLX_SNOW', 'precipitation flux (solid) in MP process',     &
                  'kg/m2/s',  2, 'XY',  ''                                                                      ), &
    VariableInfo( ATMOS_PHY_MP_AUX2D_SFLX_ENGI_ID, 'MP_SFLX_ENGI', 'internal energy flux flux in MP process',      &
                  'J/m2/s',   2, 'XY',  ''                                                                      )  /

  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains

!> Setup an object to manage variables with a cloud microphysics component  
!OCL SERIAL
  subroutine AtmosPhyMpVars_Init( this, model_mesh, &
    QS_MP, QE_MP, QA_MP )

    use scale_tracer, only: &
      TRACER_NAME, TRACER_DESC, TRACER_UNIT
    use scale_file_history, only: &
      FILE_HISTORY_reg

    implicit none
    class(AtmosPhyMpVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    integer, intent(in) :: QS_MP
    integer, intent(in) :: QE_MP
    integer, intent(in) :: QA_MP

    integer :: iv
    integer :: iq
    integer :: n
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D

    type(VariableInfo) :: qtrc_tp_vinfo_tmp
    type(VariableInfo) :: qtrc_vterm_vinfo_tmp
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
    allocate( this%tends(this%TENDS_NUM_TOT) )

    reg_file_hist = .true.    
    do iv = 1, ATMOS_PHY_MP_TENDS_NUM1
      call this%tends_manager%Regist(           &
        ATMOS_PHY_MP_TEND_VINFO(iv), mesh3D,    &
        this%tends(iv), reg_file_hist           )
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%tends(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    qtrc_tp_vinfo_tmp%ndims    = 3
    qtrc_tp_vinfo_tmp%dim_type = 'XYZ'
    qtrc_tp_vinfo_tmp%STDNAME  = ''
    
    do iq = 1, this%QA
      iv = ATMOS_PHY_MP_TENDS_NUM1 + iq 
      qtrc_tp_vinfo_tmp%keyID = iv
      qtrc_tp_vinfo_tmp%NAME  = 'MP_'//trim(TRACER_NAME(this%QS+iq-1))//'_t'
      qtrc_tp_vinfo_tmp%DESC  = 'tendency of rho*'//trim(TRACER_NAME(this%QS+iq-1))//' in MP process'
      qtrc_tp_vinfo_tmp%UNIT  = 'kg/m3/s'

      reg_file_hist = .true.
      call this%tends_manager%Regist( &
        qtrc_tp_vinfo_tmp, mesh3D,              & 
        this%tends(iv), reg_file_hist           ) 
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%tends(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do    

    !--
    
    allocate( this%vterm_hist_id(this%QS+1:this%QE) )
    allocate( this%vterm_hist   (this%QS+1:this%QE) )

    qtrc_vterm_vinfo_tmp%ndims    = 3
    qtrc_vterm_vinfo_tmp%dim_type = 'XYZ'
    qtrc_vterm_vinfo_tmp%STDNAME  = ''

    do iq = this%QS+1, this%QE
      qtrc_vterm_vinfo_tmp%NAME = 'Vterm_'//trim(TRACER_NAME(this%QS+iq-1))
      qtrc_vterm_vinfo_tmp%DESC = 'terminal velocity of '//trim(TRACER_NAME(this%QS+iq-1))
      qtrc_vterm_vinfo_tmp%UNIT = 'm/s'
      call FILE_HISTORY_reg( qtrc_vterm_vinfo_tmp%NAME, qtrc_vterm_vinfo_tmp%DESC, qtrc_vterm_vinfo_tmp%UNIT, &
        this%vterm_hist_id(iq), dim_type='XYZ'                                              )
      if ( this%vterm_hist_id(iq) > 0 ) call this%vterm_hist(iq)%Init( qtrc_vterm_vinfo_tmp%NAME, qtrc_vterm_vinfo_tmp%UNIT, mesh3D )
    end do

    !--
    
    call this%auxvars2D_manager%Init()
    allocate( this%auxvars2D(ATMOS_PHY_MP_AUX2D_NUM) )

    reg_file_hist = .true.    
    do iv = 1, ATMOS_PHY_MP_AUX2D_NUM
      call this%auxvars2D_manager%Regist( &
        ATMOS_PHY_MP_AUX2D_VINFO(iv), mesh2D,    & ! (in) 
        this%auxvars2D(iv), reg_file_hist        ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%auxvars2D(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do
    
    return
  end subroutine AtmosPhyMpVars_Init

!> Finalize an object to manage variables with cloud microphysics component  
!OCL SERIAL
  subroutine AtmosPhyMpVars_Final( this )
    implicit none
    class(AtmosPhyMpVars), intent(inout) :: this

    integer :: iq
    !--------------------------------------------------

    LOG_INFO('AtmosPhyMpVars_Final',*)

    call this%tends_manager%Final()
    deallocate( this%tends )

    call this%auxvars2D_manager%Final()
    deallocate( this%auxvars2D )

    do iq = this%QS+1, this%QE
      if ( this%vterm_hist_id(iq) > 0 ) call this%vterm_hist(iq)%Final()
    end do
    deallocate( this%vterm_hist_id )

    return
  end subroutine AtmosPhyMpVars_Final

!OCL SERIAL
  subroutine AtmosPhyMpVars_GetLocalMeshFields_tend( domID, mesh, mp_tends_list, &
    mp_DENS_t, mp_MOMX_t, mp_MOMY_t, mp_MOMZ_t, mp_RHOT_t, mp_RHOH, mp_EVAP,     &
    mp_RHOQ_t,                                                                   &
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
    class(LocalMeshFieldBase), pointer, intent(out) :: mp_RHOH
    class(LocalMeshFieldBase), pointer, intent(out) :: mp_EVAP
    type(LocalMeshFieldBaseList), intent(out) :: mp_RHOQ_t(:)
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh

    integer :: iq
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

    call mp_tends_list%Get(ATMOS_PHY_MP_RHOH_ID, field)
    call field%GetLocalMeshField(domID, mp_RHOH)

    call mp_tends_list%Get(ATMOS_PHY_MP_EVAPORATE_ID, field)
    call field%GetLocalMeshField(domID, mp_EVAP)

    !---
    do iq = 1, size(mp_RHOQ_t)
      call mp_tends_list%Get(ATMOS_PHY_MP_TENDS_NUM1 + iq, field)
      call field%GetLocalMeshField(domID, mp_RHOQ_t(iq)%ptr)
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
  end subroutine AtmosPhyMpVars_GetLocalMeshFields_tend

!OCL SERIAL
  subroutine AtmosPhyMpVars_GetLocalMeshFields_sfcflx( domID, mesh, sfcflx_list, &
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

    call sfcflx_list%Get(ATMOS_PHY_MP_AUX2D_SFLX_RAIN_ID, field)
    call field%GetLocalMeshField(domID, SFLX_rain)

    call sfcflx_list%Get(ATMOS_PHY_MP_AUX2D_SFLX_SNOW_ID, field)
    call field%GetLocalMeshField(domID, SFLX_snow)

    call sfcflx_list%Get(ATMOS_PHY_MP_AUX2D_SFLX_engi_ID, field)
    call field%GetLocalMeshField(domID, SFLX_engi)
    
    return
  end subroutine AtmosPhyMpVars_GetLocalMeshFields_sfcflx

!OCL SERIAL
  subroutine AtmosPhyMpVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosPhyMpVars), intent(inout) :: this
  
    integer :: v
    integer :: iq
    integer :: hst_id
    type(MeshField3D) :: tmp_field
    class(MeshBase3D), pointer :: mesh3D
    !-------------------------------------------------------------------------

    mesh3D => this%tends(1)%mesh

    do v = 1, this%TENDS_NUM_TOT
      hst_id = this%tends(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%tends(v) )
    end do

    do v = 1, ATMOS_PHY_MP_AUX2D_NUM
      hst_id = this%auxvars2D(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%auxvars2D(v) )
    end do

    do iq = this%QS+1, this%QE
      hst_id = this%vterm_hist_id(iq)
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%vterm_hist(iq) )
    end do

    return
  end subroutine AtmosPhyMpVars_history

end  module mod_atmos_phy_mp_vars
