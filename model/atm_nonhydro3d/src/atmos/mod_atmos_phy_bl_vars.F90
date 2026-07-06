!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Planetary Boundary Layer Turbulence
!!
!! @par Description
!!          Container for variables with planetary boundary layer (PBL) turbulence parameterization in atmospheric model
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_bl_vars
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

  !> Derived type to manage variables with planetary boundary layer (PBL) turbulence parameterization component in atmospheric model
  type, public :: AtmosPhyBlVars
    type(MeshField3D), allocatable :: tends(:)     !< Array of tendency variables
    type(ModelVarManager) :: tends_manager         !< Object to manage tendencies

    integer :: QS      !< Start index of tracer variables with PBL turbulence parameterization
    integer :: QE      !< End index of tracer variables with PBL turbulence parameterization
    integer :: QA      !< Number of tracer variables with PBL turbulence parameterization

    integer :: TENDS_NUM_TOT                        !< Number of tendency variables with PBL turbulence parameterization
  contains
    procedure :: Init => AtmosPhyBlVars_Init
    procedure :: Final => AtmosPhyBlVars_Final
    procedure :: History => AtmosPhyBlVars_history
  end type AtmosPhyBlVars

  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  integer, public, parameter :: ATMOS_PHY_BL_RHOU_t_ID    = 1
  integer, public, parameter :: ATMOS_PHY_BL_RHOV_t_ID    = 2
  integer, public, parameter :: ATMOS_PHY_BL_RHOT_t_ID    = 3
  integer, public, parameter :: ATMOS_PHY_BL_TENDS_NUM1   = 3 

  type(VariableInfo), public :: ATMOS_PHY_BL_TEND_VINFO(ATMOS_PHY_BL_TENDS_NUM1)
  DATA ATMOS_PHY_BL_TEND_VINFO / &
    VariableInfo( ATMOS_PHY_BL_RHOU_t_ID, 'BL_RHOU_t', 'tendency of x-momentum in BL process',           &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_BL_RHOV_t_ID, 'BL_RHOV_t', 'tendency of y-momentum in BL process',           &
                  'kg/m2/s2',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_BL_RHOT_t_ID, 'BL_RHOT_t', 'tendency of rho*PT in BL process',               &
                  'kg/m3.K/s', 3, 'XYZ',  ''                                                          )  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
contains
  !> Setup an object to manage variables with a planetary boundary layer (PBL) turbulence parameterization component  
!OCL SERIAL
  subroutine AtmosPhyBlVars_Init( this, model_mesh, &
    QS_BL, QE_BL, QA_BL )

    use scale_tracer, only: &
      TRACER_NAME, TRACER_DESC, TRACER_UNIT
    use scale_file_history, only: &
      FILE_HISTORY_reg

    implicit none
    class(AtmosPhyBlVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    integer, intent(in) :: QS_BL
    integer, intent(in) :: QE_BL
    integer, intent(in) :: QA_BL

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

    LOG_INFO('AtmosPhyBlVars_Init',*)

    this%QS = QS_BL
    this%QE = QE_BL
    this%QA = QA_BL
    this%TENDS_NUM_TOT = ATMOS_PHY_BL_TENDS_NUM1 + QE_BL - QS_BL + 1

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
    do iv = 1, ATMOS_PHY_BL_TENDS_NUM1
      call this%tends_manager%Regist(           &
        ATMOS_PHY_BL_TEND_VINFO(iv), mesh3D,    &
        this%tends(iv), reg_file_hist           )
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%tends(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    qtrc_tp_vinfo_tmp%ndims    = 3
    qtrc_tp_vinfo_tmp%dim_type = 'XYZ'
    qtrc_tp_vinfo_tmp%STDNAME  = ''
    
    do iq = 1, this%QA
      iv = ATMOS_PHY_BL_TENDS_NUM1 + iq 
      qtrc_tp_vinfo_tmp%keyID = iv
      qtrc_tp_vinfo_tmp%NAME  = 'BL_'//trim(TRACER_NAME(this%QS+iq-1))//'_t'
      qtrc_tp_vinfo_tmp%DESC  = 'tendency of rho*'//trim(TRACER_NAME(this%QS+iq-1))//' in BL process'
      qtrc_tp_vinfo_tmp%UNIT  = 'kg/m3/s'

      reg_file_hist = .true.
      call this%tends_manager%Regist( &
        qtrc_tp_vinfo_tmp, mesh3D,              & 
        this%tends(iv), reg_file_hist           ) 
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%tends(iv)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do    

    return
  end subroutine AtmosPhyBlVars_Init

  !> Finalize an object to manage variables with planetary boundary layer (PBL) turbulence parameterization component  
!OCL SERIAL
  subroutine AtmosPhyBlVars_Final( this )
    implicit none
    class(AtmosPhyBlVars), intent(inout) :: this
    !----------------------------------------------------
    return
  end subroutine AtmosPhyBlVars_Final

!OCL SERIAL
  subroutine AtmosPhyBlVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosPhyBlVars), intent(inout) :: this
    !----------------------------------------------------
    return
  end subroutine AtmosPhyBlVars_history

end module mod_atmos_phy_bl_vars