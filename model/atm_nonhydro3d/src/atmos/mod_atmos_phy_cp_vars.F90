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

    integer :: TENDS_NUM_TOT                        !< Number of tendency variables with cumulus parameterization
  contains
    procedure :: Init => AtmosPhyCpVars_Init
    procedure :: Final => AtmosPhyCpVars_Final
    procedure :: History => AtmosPhyCpVars_history
  end type AtmosPhyCpVars

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
    VariableInfo( ATMOS_PHY_CP_DENS_t_ID, 'CP_DENS_t', 'tendency of density in CP process',           &
                  'kg/m3',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_CP_RHOW_t_ID, 'CP_RHOW_t', 'tendency of water vapor in CP process',       &
                  'kg/m3',  3, 'XYZ',  ''                                                          ), &
    VariableInfo( ATMOS_PHY_CP_RHOT_t_ID, 'CP_RHOT_t', 'tendency of rho*PT in CP process',            &
                  'kg/m3.K/s', 3, 'XYZ',  ''                                                       ), &
    VariableInfo( ATMOS_PHY_CP_RHOQV_t_ID, 'CP_RHOQV_t', 'tendency of rho*QV in CP process',          &
                  'kg/m3', 3, 'XYZ',  ''                                                           )  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
contains
  !> Setup an object to manage variables with a cumulus parameterization component  
!OCL SERIAL
  subroutine AtmosPhyCpVars_Init( this, model_mesh )
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

    LOG_INFO('AtmosPhyCpVars_Init',*)

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

    return
  end subroutine AtmosPhyCpVars_Init

  !> Finalize an object to manage variables with cumulus parameterization component  
!OCL SERIAL
  subroutine AtmosPhyCpVars_Final( this )
    implicit none
    class(AtmosPhyCpVars), intent(inout) :: this
    !----------------------------------------------------
    return
  end subroutine AtmosPhyCpVars_Final

!OCL SERIAL
  subroutine AtmosPhyCpVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosPhyCpVars), intent(inout) :: this
    !----------------------------------------------------
    return
  end subroutine AtmosPhyCpVars_history

end module mod_atmos_phy_cp_vars