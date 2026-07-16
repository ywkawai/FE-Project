!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Radiation
!!
!! @par Description
!!          Container for variables with radiation component in atmospheric model
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_rd_vars
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
    MeshBase3D,                               &
    DIMTYPE_XYZ  => MeshBase3D_DIMTYPEID_XYZ, &
    DIMTYPE_XY   => MeshBase3D_DIMTYPEID_XY
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField2D, MeshField3D
  use scale_meshfieldcomm_base, only: MeshFieldContainer

  use scale_file_restart_meshfield, only: &
    FILE_restart_meshfield_component
  
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

  !> Derived type to manage variables with radiation component in atmospheric model
  type, public :: AtmosPhyRdVars
    type(MeshField3D), allocatable :: tends(:)     !< Array of tendency variables
    type(ModelVarManager) :: tends_manager         !< Object to manage tendencies

    type(MeshField2D), allocatable :: auxvars2D(:) !< Array of 2D auxiliary variables
    type(ModelVarManager) :: auxvars2D_manager     !< Object to manage 2D auxiliary variables

    type(MeshField2D), allocatable :: sflx_rad_dn(:) !< Array of downward radiative fluxes at the surface

    integer :: TENDS_NUM_TOT                       !< Number of tendency variables with cloud microphysics

    type(FILE_restart_meshfield_component) restart_file !< Restart file for variables with radiation component
  contains
    procedure :: Init => AtmosPhyRdVars_Init
    procedure :: Final => AtmosPhyRdVars_Final
    procedure :: History => AtmosPhyRdVars_history
    procedure :: Read_restart_file => AtmosPhyRdVars_Read_restart_file
    procedure :: Write_restart_file_prep => AtmosPhyRdVars_Write_restart_file_prep
    procedure :: Write_restart_file => AtmosPhyRdVars_Write_restart_file
    procedure :: Write_restart_file_post => AtmosPhyRdVars_Write_restart_file_post
  end type AtmosPhyRdVars
  
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  integer, public, parameter :: ATMOS_PHY_RD_RHOH_ID      = 1 !< ID of diabatic heating rate in RD process
  integer, public, parameter :: ATMOS_PHY_RD_TENDS_NUM1   = 1

  type(VariableInfo), public :: ATMOS_PHY_RD_TEND_VINFO(ATMOS_PHY_RD_TENDS_NUM1)
  DATA ATMOS_PHY_RD_TEND_VINFO / &
    VariableInfo( ATMOS_PHY_RD_RHOH_ID, 'RD_RHOH', 'diabatic heating rate in RD process',              &
                  'J/kg/s',   3, 'XYZ',  ''                                                         )  /

  !-- 2D Auxiliary variables for radiation

  ! Radiative fluxes at the surface
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_SFLX_LW_up_ID   = 1  !< ID of upward longwave surface flux
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_SFLX_LW_dn_ID   = 2  !< ID of downward longwave surface flux
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_SFLX_SW_up_ID   = 3  !< ID of upward shortwave surface flux
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_SFLX_SW_dn_ID   = 4  !< ID of downward shortwave surface flux
  integer, public, parameter :: ATMOS_PHY_RD_RESTART_VAR_NUM       = 4  !< Number of restart variables for radiation component
  ! Radiative fluxes at the top of the model
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_TOMFLX_LW_up_ID = 5  !< ID of upward longwave flux at the top of the model
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_TOMFLX_LW_dn_ID = 6  !< ID of downward longwave flux at the top of the model
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_TOMFLX_SW_up_ID = 7  !< ID of upward shortwave flux at the top of the model
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_TOMFLX_SW_dn_ID = 8  !< ID of downward shortwave flux at the top of the model
  !
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_SOLINS_ID       = 9  !< ID of solar insolation flux at the top of the model
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_COSSZA_ID       = 10 !< ID of cosine of solar zenith angle
  integer, public, parameter :: ATMOS_PHY_RD_AUX2D_NUM             = 10

  type(VariableInfo), public :: ATMOS_PHY_RD_AUX2D_VINFO(ATMOS_PHY_RD_AUX2D_NUM)
  DATA ATMOS_PHY_RD_AUX2D_VINFO / &
    VariableInfo( ATMOS_PHY_RD_AUX2D_SFLX_LW_up_ID  , 'RD_SFLX_LW_up'  , 'upward longwave surface flux'                   , 'W/m2', 2, 'XY', '' ), &
    VariableInfo( ATMOS_PHY_RD_AUX2D_SFLX_LW_dn_ID  , 'RD_SFLX_LW_dn'  , 'downward longwave surface flux'                 , 'W/m2', 2, 'XY', '' ), &
    VariableInfo( ATMOS_PHY_RD_AUX2D_SFLX_SW_up_ID  , 'RD_SFLX_SW_up'  , 'upward shortwave surface flux'                  , 'W/m2', 2, 'XY', '' ), &
    VariableInfo( ATMOS_PHY_RD_AUX2D_SFLX_SW_dn_ID  , 'RD_SFLX_SW_dn'  , 'downward shortwave surface flux'                , 'W/m2', 2, 'XY', '' ), &
    VariableInfo( ATMOS_PHY_RD_AUX2D_TOMFLX_LW_up_ID, 'RD_TOMFLX_LW_up', 'upward longwave flux at the top of the model'   , 'W/m2', 2, 'XY', '' ), &
    VariableInfo( ATMOS_PHY_RD_AUX2D_TOMFLX_LW_dn_ID, 'RD_TOMFLX_LW_dn', 'downward longwave flux at the top of the model' , 'W/m2', 2, 'XY', '' ), &
    VariableInfo( ATMOS_PHY_RD_AUX2D_TOMFLX_SW_up_ID, 'RD_TOMFLX_SW_up', 'upward shortwave flux at the top of the model'  , 'W/m2', 2, 'XY', '' ), &
    VariableInfo( ATMOS_PHY_RD_AUX2D_TOMFLX_SW_dn_ID, 'RD_TOMFLX_SW_dn', 'downward shortwave flux at the top of the model', 'W/m2', 2, 'XY', '' ), &
    VariableInfo( ATMOS_PHY_RD_AUX2D_SOLINS_ID      , 'RD_SolINS'      , 'solar insolation flux at the top of the model'  , 'W/m2', 2, 'XY', '' ), &
    VariableInfo( ATMOS_PHY_RD_AUX2D_COSSZA_ID      , 'RD_cosSZA'      , 'cosine of solar zenith angle'                   ,    '1', 2, 'XY', '' )  /
  
  !-
  integer, public, parameter :: ATMOS_PHY_RD_SFLX_dn_SW_dir = 1  !< ID of downward direct shortwave flux at the surface
  integer, public, parameter :: ATMOS_PHY_RD_SFLX_dn_SW_dif = 2  !< ID of downward diffuse shortwave flux at the surface
  integer, public, parameter :: ATMOS_PHY_RD_SFLX_dn_LW     = 3  !< ID of downward longwave flux at the surface
  integer, public, parameter :: ATMOS_PHY_RD_SFLX_dn_NUM    = 3  !< Number of downward radiative fluxes
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
contains
  !> Setup an object to manage variables with a radiation component  
!OCL SERIAL
  subroutine AtmosPhyRdVars_Init( this, model_mesh )

    use scale_tracer, only: &
      TRACER_NAME, TRACER_DESC, TRACER_UNIT
    use scale_file_history, only: &
      FILE_HISTORY_reg

    implicit none
    class(AtmosPhyRdVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D

    integer :: iv
    integer :: iq
    integer :: n
    logical :: reg_file_hist
    !----------------------------------------------------

    LOG_INFO('AtmosPhyRdVars_Init',*)

    this%TENDS_NUM_TOT = ATMOS_PHY_RD_TENDS_NUM1

    !- Initialize auxiliary and diagnostic variables

    nullify( atm_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select
    mesh3D => atm_mesh%ptr_mesh

    call mesh3D%GetMesh2D( mesh2D )

    !----

    !-
    call this%tends_manager%Init()
    allocate( this%tends(this%TENDS_NUM_TOT) )

    reg_file_hist = .true.
    do iv=1, ATMOS_PHY_RD_TENDS_NUM1
      call this%tends_manager%Regist( &
        ATMOS_PHY_RD_TEND_VINFO(iv), mesh3D, & ! (in) 
        this%tends(iv), reg_file_hist        ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%tends(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    !- 
    call this%auxvars2D_manager%Init()
    allocate( this%auxvars2D(ATMOS_PHY_RD_AUX2D_NUM) )

    reg_file_hist = .true.
    do iv=1, ATMOS_PHY_RD_AUX2D_NUM
      call this%auxvars2D_manager%Regist( &
        ATMOS_PHY_RD_AUX2D_VINFO(iv), mesh2D, & ! (in) 
        this%auxvars2D(iv), reg_file_hist     ) ! (out)
      
      do n = 1, mesh2D%LOCAL_MESH_NUM
        this%auxvars2D(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    !- Setup restart file for variables with radiation component
    call atm_mesh%Setup_restartfile( this%restart_file, &
      ATMOS_PHY_RD_RESTART_VAR_NUM                      )
    
    return
  end subroutine AtmosPhyRdVars_Init

  !> Finalize an object to manage variables with radiation component  
!OCL SERIAL
  subroutine AtmosPhyRdVars_Final( this )
    implicit none
    class(AtmosPhyRdVars), intent(inout) :: this
    !----------------------------------------------------

    LOG_INFO('AtmosPhyRdVars_Final',*)

    call this%tends_manager%Final()
    deallocate( this%tends )

    call this%auxvars2D_manager%Final()
    deallocate( this%auxvars2D )

    return
  end subroutine AtmosPhyRdVars_Final

  !> Write history data of variables with radiation component
!OCL SERIAL
  subroutine AtmosPhyRdVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosPhyRdVars), intent(inout) :: this

    integer :: v
    integer :: hst_id
    !----------------------------------------------------

    do v = 1, this%TENDS_NUM_TOT
      hst_id = this%tends(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%tends(v) )
    end do

    do v = 1, ATMOS_PHY_RD_AUX2D_NUM
      hst_id = this%auxvars2D(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%auxvars2D(v) )
    end do

    return
  end subroutine AtmosPhyRdVars_history

!> Read data with atmospheric variables with radiation component from restart file
!!
!OCL SERIAL
  subroutine AtmosPhyRdVars_Read_restart_file( this )
    use scale_meshfieldcomm_base, only: MeshFieldContainer
    implicit none
    
    class(AtmosPhyRdVars), intent(inout), target :: this

    integer :: iv
    !--------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("AtmosPhyRdVars_Read_restart_file",*) 'Open restart file (ATMOS/Physics radiation) '

    !- Open restart file
    call this%restart_file%Open()

    !- Read restart file

    do iv=ATMOS_PHY_RD_AUX2D_SFLX_LW_up_ID, ATMOS_PHY_RD_AUX2D_SFLX_SW_dn_ID
      call this%restart_file%Read_var( DIMTYPE_XY, this%auxvars2D(iv)%varname, &
        this%auxvars2D(iv) )
    end do

    !- Close restart file
    LOG_INFO("AtmosPhyRdVars_Read_restart_file",*) 'Close restart file (ATMOS/Physics radiation) '
    call this%restart_file%Close()
    return
  end subroutine AtmosPhyRdVars_Read_restart_file

!> Write data with atmospheric variables with radiation component to restart file
!!
!OCL SERIAL
  subroutine AtmosPhyRdVars_Write_restart_file_prep( this )
    implicit none
    class(AtmosPhyRdVars), intent(inout) :: this

    integer :: iv, rf_vid
    !---------------------------------------

    LOG_NEWLINE
    LOG_INFO("AtmosPhyRdVars_Write_restart_file_prep",*) 'Define variables with ATMOS/Physics radiation'

    !- Create restart file
    call this%restart_file%Create()
    call PRC_mpibarrier()

    !- Define variables

    do iv=ATMOS_PHY_RD_AUX2D_SFLX_LW_up_ID, ATMOS_PHY_RD_AUX2D_SFLX_SW_dn_ID
      rf_vid = iv
      call this%restart_file%Def_var( this%auxvars2D(iv),     &
        ATMOS_PHY_RD_AUX2D_VINFO(iv)%DESC, rf_vid, DIMTYPE_XY )
    end do
    return
  end subroutine AtmosPhyRdVars_Write_restart_file_prep

!> Write data with atmospheric variables with radiation component to restart file
!!
!OCL SERIAL
  subroutine AtmosPhyRdVars_Write_restart_file( this )
    implicit none
    class(AtmosPhyRdVars), intent(inout) :: this

    integer :: iv, rf_vid 
    !---------------------------------------
    !- Write restart file
    do iv=ATMOS_PHY_RD_AUX2D_SFLX_LW_up_ID, ATMOS_PHY_RD_AUX2D_SFLX_SW_dn_ID
      rf_vid = iv
      call this%restart_file%Write_var(rf_vid, this%auxvars2D(iv) )
    end do
    return
  end subroutine AtmosPhyRdVars_Write_restart_file

  !> Close restart file with atmospheric radiation component
  !!
!OCL SERIAL
  subroutine AtmosPhyRdVars_Write_restart_file_post( this )
    implicit none
    class(AtmosPhyRdVars), intent(inout) :: this
    !---------------------------------------
    LOG_INFO("AtmosPhyRdVars_Write_restart_file_post",*) 'Close restart file (ATMOS/Physics radiation) '

    call this%restart_file%Close()
    return
  end subroutine AtmosPhyRdVars_Write_restart_file_post

end module mod_atmos_phy_rd_vars