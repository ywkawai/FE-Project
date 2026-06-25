!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Radiation
!!
!! @par Description
!!          Module for radiation component in atmospheric model
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_rd
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof
  use scale_const, only: &
    UNDEF8 => CONST_UNDEF8

  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase, &
    ElementBase1D, ElementBase2D, ElementBase3D

  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField3D
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_component_proc, only:  ModelComponentProc

  use mod_atmos_phy_rd_vars, only: AtmosPhyRdVars

  use mod_atmos_vars_container, only: &
    AtmosVarsContainer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage a component of radiation in atmospheric model
  !!
  type, extends(ModelComponentProc), public :: AtmosPhyRd
    integer :: RD_TYPEID         !< Type id of radiation scheme
    type(AtmosPhyRdVars) :: vars !< Object to manage variables with radiation

    integer :: atm_var_container_typeid     !< Type ID of variable container for radiation

  contains
    procedure :: setup => AtmosPhyRd_setup
    procedure :: calc_tendency => AtmosPhyRd_calc_tendency
    procedure :: update => AtmosPhyRd_update
    procedure :: finalize => AtmosPhyRd_finalize
  end type AtmosPhyRd

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
contains

!> Setup a component of radiation in atmospheric model
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param tm_parent_comp Object to mange a temporal scheme in a parent component
!!
  subroutine AtmosPhyRd_setup( this, model_mesh, tm_parent_comp )
    use mod_atmos_mesh, only: AtmosMesh
    use scale_time_manager, only: TIME_manager_component
    use mod_atmos_vars, only: ATM_VARS_CONTAINER_PRIMARY_ID
    implicit none
    class(AtmosPhyRd), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8 !< Timestep for radiation
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  !< Unit of timestep

    character(len=H_MID) :: RD_TYPE = 'NONE'                 !< Type of a radiation scheme
    integer :: atm_var_container_typeid

    namelist /PARAM_ATMOS_PHY_RD/ &
      TIME_DT,             &
      TIME_DT_UNIT,        &
      RD_TYPE,             &
      atm_var_container_typeid

    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D

    integer :: ierr
    !-----------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_setup",*) 'Setup'

    atm_var_container_typeid = ATM_VARS_CONTAINER_PRIMARY_ID

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_RD_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_PHY_RD_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_RD. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_RD)

    this%atm_var_container_typeid = atm_var_container_typeid

    !--- Register this component in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_RD', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                       ! (out) 

    !--- Set the type of radiation scheme
    
    select case ( RD_TYPE )
    case default
      LOG_ERROR("ATMOS_PHY_RD_setup",*) 'Not appropriate RD_TYPE. Check!'
      call PRC_abort
    end select


    !- Initialize the variables 
    call this%vars%Init( model_mesh )

    return
  end subroutine AtmosPhyRd_setup

!> Calculate tendencies associated with radiation in atmospheric model
!!
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param prgvars_list Object to manage prognostic variables with atmospheric dynamical core
!! @param trcvars_list Object to manage auxiliary variables 
!! @param forcing_list Object to manage forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL
  subroutine AtmosPhyRd_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )
    implicit none
    class(AtmosPhyRd), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !------------------------------------------------------------------------
    return
  end subroutine AtmosPhyRd_calc_tendency

!> Update variables in a component of radiation in atmospheric model
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param prgvars_list Object to manage prognostic variables with atmospheric dynamical core
!! @param trcvars_list Object to manage auxiliary variables 
!! @param forcing_list Object to manage forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL  
  subroutine AtmosPhyRd_update( this, model_mesh, &
    prgvars_list, trcvars_list,                   &
    auxvars_list, forcing_list, is_update         )  
    
    implicit none
    class(AtmosPhyRd), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------
    return
  end subroutine AtmosPhyRd_update

!> Finalize a component of radiation in atmospheric model
!!
!OCL SERIAL  
  subroutine AtmosPhyRd_finalize( this )
    implicit none
    class(AtmosPhyRd), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    ! select case ( this%RD_TYPEID )
    ! case( RD_TYPEID_LSCOND )
    ! end select

    call this%vars%Final()
    return
  end subroutine AtmosPhyRd_finalize

!- private ------------------------------------------------

end module mod_atmos_phy_rd
