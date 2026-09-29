!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Cumulus Parameterization
!!
!! @par Description
!!          Container for variables with cumulus parameterization component in atmospheric model
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_cp_vars
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

  !> Derived type to manage variables with cumulus parameterization component in atmospheric model
  type, public :: AtmosPhyCpVars
    type(MeshField3D), allocatable :: tends(:)     !< Array of tendency variables
    type(ModelVarManager) :: tends_manager         !< Object to manage tendencies

    type(MeshField2D), allocatable :: auxvars2D(:) !< Array of 2D auxiliary variables
    type(ModelVarManager) :: auxvars2D_manager     !< Object to manage 2D auxiliary variables

    integer :: TENDS_NUM_TOT                        !< Number of tendency variables with cumulus parameterization
  contains
    procedure :: Init => AtmosPhyCpVars_Init
    procedure :: Final => AtmosPhyCpVars_Final
    procedure :: Setup => AtmosPhyCpVars_Setup
    procedure :: History => AtmosPhyCpVars_history
  end type AtmosPhyCpVars

  public :: AtmosPhyCpVars_GetLocalMeshFields_tend
  public :: AtmosPhyCPVars_GetLocalMeshFields_sfcflx

  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  integer, public, parameter :: ATMOS_PHY_CP_DENS_t_ID    = 1
  integer, public, parameter :: ATMOS_PHY_CP_RHOW_t_ID    = 2
  integer, public, parameter :: ATMOS_PHY_CP_RHOT_t_ID    = 3
  integer, public, parameter :: ATMOS_PHY_CP_RHOQV_t_ID   = 4
  integer, public, parameter :: ATMOS_PHY_CP_TENDS_NUM1   = 4 

  type(VariableInfo), public :: ATMOS_PHY_CP_TEND_VINFO(ATMOS_PHY_CP_TENDS_NUM1)
  DATA ATMOS_PHY_CP_TEND_VINFO / &
    VariableInfo( ATMOS_PHY_CP_DENS_t_ID, 'CP_DENS_t', 'tendency of density in CP process',          &
                  'kg/m3/s',  3, 'XYZ',  ''                                                       ), &
    VariableInfo( ATMOS_PHY_CP_RHOW_t_ID, 'CP_RHOW_t', 'tendency of water vapor in CP process',      &
                  'kg/m3/s',  3, 'XYZ',  ''                                                       ), &
    VariableInfo( ATMOS_PHY_CP_RHOT_t_ID, 'CP_RHOT_t', 'tendency of rho*PT in CP process',           &
                  'kg/m3.K/s', 3, 'XYZ',  ''                                                      ), &
    VariableInfo( ATMOS_PHY_CP_RHOQV_t_ID, 'CP_RHOQV_t', 'tendency of rho*QV in CP process',         &
                  'kg/m3/s', 3, 'XYZ',  ''                                                        )  /

  integer, public, parameter :: ATMOS_PHY_CP_AUX2D_SFLX_RAIN_ID   = 1
  integer, public, parameter :: ATMOS_PHY_CP_AUX2D_SFLX_SNOW_ID   = 2
  integer, public, parameter :: ATMOS_PHY_CP_AUX2D_SFLX_ENGI_ID   = 3
  integer, public, parameter :: ATMOS_PHY_CP_AUX2D_NUM            = 3

  type(VariableInfo), public :: ATMOS_PHY_CP_AUX2D_VINFO(ATMOS_PHY_CP_AUX2D_NUM)
  DATA ATMOS_PHY_CP_AUX2D_VINFO / &
    VariableInfo( ATMOS_PHY_CP_AUX2D_SFLX_RAIN_ID, 'CP_SFLX_RAIN', 'precipitation flux (liquid) in CP process',    &
                  'kg/m2/s',  2, 'XY',  ''                                                                      ), &
    VariableInfo( ATMOS_PHY_CP_AUX2D_SFLX_SNOW_ID, 'CP_SFLX_SNOW', 'precipitation flux (solid) in CP process',     &
                  'kg/m2/s',  2, 'XY',  ''                                                                      ), &
    VariableInfo( ATMOS_PHY_CP_AUX2D_SFLX_ENGI_ID, 'CP_SFLX_ENGI', 'internal energy flux flux in CP process',      &
                  'J/m2/s',   2, 'XY',  ''                                                                      )  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
contains
  !> Setup an object to manage variables with a cumulus parameterization component  
!OCL SERIAL
  subroutine AtmosPhyCpVars_Init( this, model_mesh )
    implicit none
    class(AtmosPhyCpVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    !----------------------------------------------------

    LOG_INFO('AtmosPhyCpVars_Init',*)
    return
  end subroutine AtmosPhyCpVars_Init

  !> Setup variable objects with cumulus parameterization component
  subroutine AtmosPhyCpVars_Setup( this, model_mesh )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       HYD_NAME
    use scale_tracer, only: &
      TRACER_NAME, TRACER_DESC, TRACER_UNIT
    use scale_file_history, only: &
      FILE_HISTORY_reg
    implicit none
    class(AtmosPhyCpVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    integer :: iv
    integer :: iq
    integer :: n
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D

    type(VariableInfo) :: qtrc_tp_vinfo_tmp
    type(VariableInfo) :: qtrc_vterm_vinfo_tmp
    !----------------------------------------------------


    this%TENDS_NUM_TOT = ATMOS_PHY_CP_TENDS_NUM1 + N_HYD

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
    do iv = 1, ATMOS_PHY_CP_TENDS_NUM1
      call this%tends_manager%Regist(           &
        ATMOS_PHY_CP_TEND_VINFO(iv), mesh3D,    &
        this%tends(iv), reg_file_hist           )
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%tends(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    qtrc_tp_vinfo_tmp%ndims    = 3
    qtrc_tp_vinfo_tmp%dim_type = 'XYZ'
    qtrc_tp_vinfo_tmp%STDNAME  = ''
    
    do iq = 1, N_HYD
      iv = ATMOS_PHY_CP_TENDS_NUM1 + iq 
      qtrc_tp_vinfo_tmp%keyID = iv
      qtrc_tp_vinfo_tmp%NAME  = 'CP_'//trim(HYD_NAME(iq))//'_t'
      qtrc_tp_vinfo_tmp%DESC  = 'tendency of rho*'//trim(HYD_NAME(iq))//' in CP process'
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
    
    call this%auxvars2D_manager%Init()
    allocate( this%auxvars2D(ATMOS_PHY_CP_AUX2D_NUM) )

    reg_file_hist = .true.    
    do iv = 1, ATMOS_PHY_CP_AUX2D_NUM
      call this%auxvars2D_manager%Regist( &
        ATMOS_PHY_CP_AUX2D_VINFO(iv), mesh2D,    & ! (in) 
        this%auxvars2D(iv), reg_file_hist        ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%auxvars2D(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    return
  end subroutine AtmosPhyCpVars_Setup

  !> Finalize an object to manage variables with cumulus parameterization component  
!OCL SERIAL
  subroutine AtmosPhyCpVars_Final( this )
    implicit none
    class(AtmosPhyCpVars), intent(inout) :: this
    !----------------------------------------------------

     LOG_INFO('AtmosPhyCpVars_Final',*)

     call this%tends_manager%Final()
     deallocate( this%tends )

     call this%auxvars2D_manager%Final()
     deallocate( this%auxvars2D )

    return
  end subroutine AtmosPhyCpVars_Final

!OCL SERIAL
  subroutine AtmosPhyCpVars_GetLocalMeshFields_tend( domID, mesh, bl_tends_list, &
    cp_DENS_t, cp_RHOT_t, cp_RHOQv_t,                                             &
    lcmesh3D                                                                     &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: bl_tends_list
    class(LocalMeshFieldBase), pointer, intent(out) :: cp_DENS_t
    class(LocalMeshFieldBase), pointer, intent(out) :: cp_RHOT_t
    class(LocalMeshFieldBase), pointer, intent(out) :: cp_RHOQv_t
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh

    integer :: iq
    !-------------------------------------------------------

    !--
    call bl_tends_list%Get(ATMOS_PHY_CP_DENS_t_ID, field)
    call field%GetLocalMeshField(domID, cp_DENS_t)

    call bl_tends_list%Get(ATMOS_PHY_CP_RHOT_t_ID, field)
    call field%GetLocalMeshField(domID, cp_RHOT_t)

    call bl_tends_list%Get(ATMOS_PHY_CP_RHOQv_t_ID, field)
    call field%GetLocalMeshField(domID, cp_RHOQv_t)

    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosPhyCpVars_GetLocalMeshFields_tend

!OCL SERIAL
  subroutine AtmosPhyCPVars_GetLocalMeshFields_sfcflx( domID, mesh, sfcflx_list, &
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

    call sfcflx_list%Get(ATMOS_PHY_CP_AUX2D_SFLX_RAIN_ID, field)
    call field%GetLocalMeshField(domID, SFLX_rain)

    call sfcflx_list%Get(ATMOS_PHY_CP_AUX2D_SFLX_SNOW_ID, field)
    call field%GetLocalMeshField(domID, SFLX_snow)

    call sfcflx_list%Get(ATMOS_PHY_CP_AUX2D_SFLX_engi_ID, field)
    call field%GetLocalMeshField(domID, SFLX_engi)
    
    return
  end subroutine AtmosPhyCPVars_GetLocalMeshFields_sfcflx

!OCL SERIAL
  subroutine AtmosPhyCpVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosPhyCpVars), intent(inout) :: this

    integer :: v
    integer :: hst_id
    !----------------------------------------------------

    do v=1, this%TENDS_NUM_TOT
      hst_id = this%tends(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%tends(v) )
    end do

    do v=1, ATMOS_PHY_CP_AUX2D_NUM
      hst_id = this%auxvars2D(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%auxvars2D(v) )
    end do
    return
  end subroutine AtmosPhyCpVars_history

end module mod_atmos_phy_cp_vars