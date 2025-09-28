!-------------------------------------------------------------------------------
!> module ATMOSPHERE physics / surface process
!!
!! @par Description
!!          Container for variables with surface model
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_sfc_vars
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

  !> Derived type to manage variables with a surface component
  !!
  type, public :: AtmosPhySfcVars
    type(MeshField2D), allocatable :: SFC_VARS(:)
    type(ModelVarManager) :: SFCVARS_manager

    type(MeshField2D), allocatable :: SFC_FLX(:)
    type(ModelVarManager) :: SFCFLX_manager

    integer :: SFCFLX_NUM_TOT
  contains
    procedure :: Init => AtmosPhySfcVars_Init
    procedure :: Final => AtmosPhySfcVars_Final
    procedure :: History => AtmosPhySfcVars_history
  end type AtmosPhySfcVars

  integer, public, parameter :: ATMOS_PHY_SF_SVAR_TEMP_ID  = 1
  integer, public, parameter :: ATMOS_PHY_SF_SVAR_NUM      = 1
  type(VariableInfo), public :: ATMOS_PHY_SF_SVAR_VINFO(ATMOS_PHY_SF_SVAR_NUM)
  DATA ATMOS_PHY_SF_SVAR_VINFO / &
    VariableInfo( ATMOS_PHY_SF_SVAR_NUM, 'SFC_TEMP', 'surface skin temperature',    &
                  'K',  2, 'XY',  ''                                             )  /

  integer, public, parameter :: ATMOS_PHY_SF_SFLX_MU_ID  = 1
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_MV_ID  = 2
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_MW_ID  = 3
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_SH_ID  = 4
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_LH_ID  = 5
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_QV_ID  = 6
  integer, public, parameter :: ATMOS_PHY_SF_SFLX_NUM1   = 6

  type(VariableInfo), public :: ATMOS_PHY_SF_SFLX_VINFO(ATMOS_PHY_SF_SFLX_NUM1)
  DATA ATMOS_PHY_SF_SFLX_VINFO / &
    VariableInfo( ATMOS_PHY_SF_SFLX_MU_ID, 'SFLX_MU', 'x-momentum flux',    &
                  'm/s*kg/m2/s',  2, 'XY',  ''                           ), &
    VariableInfo( ATMOS_PHY_SF_SFLX_MV_ID, 'SFLX_MV', 'y-momentum flux',    &
                  'm/s*kg/m2/s',  2, 'XY',  ''                           ), &
    VariableInfo( ATMOS_PHY_SF_SFLX_MW_ID, 'SFLX_MW', 'z-momentum flux',    &
                  'm/s*kg/m2/s',  2, 'XY',  ''                           ), &
    VariableInfo( ATMOS_PHY_SF_SFLX_SH_ID, 'SFLX_SH', 'sensible heat flux', &
                  'J/m2/s',  2, 'XY',  ''                                ), &
    VariableInfo( ATMOS_PHY_SF_SFLX_LH_ID, 'SFLX_LH', 'latent heat flux',   &
                  'J/m2/s',  2, 'XY',  ''                                ), &
    VariableInfo( ATMOS_PHY_SF_SFLX_QV_ID, 'SFLX_QV', 'water vapor flux',   &
                  'kg/m2/s',  2, 'XY',  ''                               )  / 

  public :: AtmosPhySfcVars_GetLocalMeshFields

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: ATMOS_PHY_SFC_DEFAULT_SFC_TEMP = 300.0_RP

contains
  subroutine AtmosPhySfcVars_Init( this, model_mesh )
    implicit none
    class(AtmosPhySfcVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    integer :: v
    integer :: n
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D

    namelist / PARAM_ATMOS_PHY_SFC_VARS / &
      ATMOS_PHY_SFC_DEFAULT_SFC_TEMP
      
    integer :: ierr
    !--------------------------------------------------

    LOG_INFO('AtmosPhySfcVars_Init',*)

    this%SFCFLX_NUM_TOT = ATMOS_PHY_SF_SFLX_NUM1

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SFC_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_SFC_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_SFC_vars_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_SFC_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_SFC_VARS)

    !- Initialize auxiliary and diagnostic variables

    nullify( atm_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select
    mesh3D => atm_mesh%ptr_mesh

    call mesh3D%GetMesh2D( mesh2D )

    !----
    call this%SFCVARS_manager%Init()
    allocate( this%SFC_VARS(ATMOS_PHY_SF_SVAR_NUM) )

    reg_file_hist = .true.    
    do v = 1, ATMOS_PHY_SF_SVAR_NUM
      call this%SFCVARS_manager%Regist( &
        ATMOS_PHY_SF_SVAR_VINFO(v), mesh2D,               & 
        this%SFC_VARS(v), reg_file_hist, fill_zero=.true. )
    end do

    do n=1, mesh2D%LOCAL_MESH_NUM
      this%SFC_VARS(ATMOS_PHY_SF_SVAR_TEMP_ID)%local(n)%val(:,:) = ATMOS_PHY_SFC_DEFAULT_SFC_TEMP
    end do

    !----
    call this%SFCFLX_manager%Init()
    allocate( this%SFC_FLX(this%SFCFLX_NUM_TOT) )

    reg_file_hist = .true.    
    do v = 1, this%SFCFLX_NUM_TOT
      call this%SFCFLX_manager%Regist( &
        ATMOS_PHY_SF_SFLX_VINFO(v), mesh2D,              & 
        this%SFC_FLX(v), reg_file_hist, fill_zero=.true. )
    end do

    return
  end subroutine AtmosPhySfcVars_Init

  subroutine AtmosPhySfcVars_Final( this )
    implicit none
    class(AtmosPhySfcVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosPhySfcVars_Final',*)

    call this%SFCVARS_manager%Final()
    deallocate( this%SFC_VARS )

    call this%SFCFLX_manager%Final()
    deallocate( this%SFC_FLX )

    return
  end subroutine AtmosPhySfcVars_Final

!OCL SERIAL
  subroutine AtmosPhySfcVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosPhySfcVars), intent(inout) :: this
  
    integer :: v
    integer :: iq
    integer :: hst_id
    !-------------------------------------------------------------------------

    do v = 1, ATMOS_PHY_SF_SVAR_NUM
      hst_id = this%SFC_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%SFC_VARS(v) )
    end do

    do v = 1, this%SFCFLX_NUM_TOT
      hst_id = this%SFC_FLX(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%SFC_FLX(v) )
    end do

    return
  end subroutine AtmosPhySfcVars_history

  subroutine AtmosPhySfcVars_GetLocalMeshFields( domID, mesh, svars_list, sflx_list, &
    SFC_TEMP, SFLX_MU, SFLX_MV, SFLX_MW, SFLX_SH, SFLX_LH, SFLX_QV,              &
    lcmesh3D                                                                     &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: svars_list
    class(ModelVarManager), intent(inout) :: sflx_list
    class(LocalMeshFieldBase), pointer, intent(out) :: SFC_TEMP
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_MU
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_MV
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_MW
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_SH
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_LH
    class(LocalMeshFieldBase), pointer, intent(out) :: SFLX_QV
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call svars_list%Get(ATMOS_PHY_SF_SVAR_TEMP_ID, field)
    call field%GetLocalMeshField(domID, SFC_TEMP)

    call sflx_list%Get(ATMOS_PHY_SF_SFLX_MU_ID, field)
    call field%GetLocalMeshField(domID, SFLX_MU)

    call sflx_list%Get(ATMOS_PHY_SF_SFLX_MV_ID, field)
    call field%GetLocalMeshField(domID, SFLX_MV)

    call sflx_list%Get(ATMOS_PHY_SF_SFLX_MW_ID, field)
    call field%GetLocalMeshField(domID, SFLX_MW)

    call sflx_list%Get(ATMOS_PHY_SF_SFLX_SH_ID, field)
    call field%GetLocalMeshField(domID, SFLX_SH)

    call sflx_list%Get(ATMOS_PHY_SF_SFLX_LH_ID, field)
    call field%GetLocalMeshField(domID, SFLX_LH)

    call sflx_list%Get(ATMOS_PHY_SF_SFLX_QV_ID, field)
    call field%GetLocalMeshField(domID, SFLX_QV)
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
  end subroutine AtmosPhySfcVars_GetLocalMeshFields

end  module mod_atmos_phy_sfc_vars
