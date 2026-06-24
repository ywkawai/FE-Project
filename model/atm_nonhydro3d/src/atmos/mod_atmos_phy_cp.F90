!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Cumulus Parameterization
!!
!! @par Description
!!          Module for cumulus parameterization
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_cp
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

  use mod_atmos_phy_cp_vars, only: AtmosPhyCpVars

  use mod_atmos_vars_container, only: &
    AtmosVarsContainer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage a component of cumulus parameterization in atmospheric model
  !!
  type, extends(ModelComponentProc), public :: AtmosPhyCp
    integer :: CP_TYPEID         !< Type id of cumulus parameterization scheme
    type(AtmosPhyCpVars) :: vars !< Object to manage variables with cumulus parameterization

    integer :: atm_var_container_typeid     !< Type ID of variable container for cumulus parameterization

  contains
    procedure :: setup => AtmosPhyCp_setup
    procedure :: calc_tendency => AtmosPhyCp_calc_tendency
    procedure :: update => AtmosPhyCp_update
    procedure :: finalize => AtmosPhyCp_finalize
  end type AtmosPhyCp

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

!> Setup a component of cumulus parameterization in atmospheric model
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param tm_parent_comp Object to mange a temporal scheme in a parent component
!!
  subroutine AtmosPhyCp_setup( this, model_mesh, tm_parent_comp )
    use scale_atmos_hydrometeor, only: &
       N_HYD    
    use mod_atmos_mesh, only: AtmosMesh
    use scale_time_manager, only: TIME_manager_component
    use mod_atmos_vars, only: ATM_VARS_CONTAINER_PRIMARY_ID
    implicit none
    class(AtmosPhyCp), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8 !< Timestep for cumulus parameterization
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  !< Unit of timestep

    character(len=H_MID) :: CP_TYPE = 'NONE'              !< Type of a cumulus parameterization scheme
    integer :: atm_var_container_typeid

    namelist /PARAM_ATMOS_PHY_CP/ &
      TIME_DT,             &
      TIME_DT_UNIT,        &
      CP_TYPE,             &
      atm_var_container_typeid

    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D

    integer :: ierr
    !-----------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_CP_setup",*) 'Setup'

    atm_var_container_typeid = ATM_VARS_CONTAINER_PRIMARY_ID

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_CP_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_PHY_CP_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_CP. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_CP)
 
    !- Get atmospheric mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select

    !--- Register this component in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_CP', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                       ! (out) 

    !--- Set the type of cumulus parameterization

    select case( CP_TYPE )
    case default
      LOG_ERROR("ATMOS_PHY_CP_setup",*) 'Not appropriate cumulus parameterization type. Check!'
      call PRC_abort
    end select

    !- Initialize the variables 
    call this%vars%Init( model_mesh )

    return
  end subroutine AtmosPhyCp_setup

!> Calculate tendencies associated with cumulus parameterization in atmospheric model
!!
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param prgvars_list Object to manage prognostic variables with atmospheric dynamical core
!! @param trcvars_list Object to manage auxiliary variables 
!! @param forcing_list Object to manage forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL
  subroutine AtmosPhyCp_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )
    implicit none
    class(AtmosPhyCp), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !------------------------------------------------------------------------
    return
  end subroutine AtmosPhyCp_calc_tendency

!> Update variables in a component of cumulus parameterization in atmospheric model
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param prgvars_list Object to manage prognostic variables with atmospheric dynamical core
!! @param trcvars_list Object to manage auxiliary variables 
!! @param forcing_list Object to manage forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL  
  subroutine AtmosPhyCp_update( this, model_mesh, &
    prgvars_list, trcvars_list,                   &
    auxvars_list, forcing_list, is_update         )  
    
    implicit none
    class(AtmosPhyCp), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------
    return
  end subroutine AtmosPhyCp_update

!> Finalize a component of cumulus parameterization in atmospheric model
!!
!OCL SERIAL  
  subroutine AtmosPhyCp_finalize( this )
    implicit none
    class(AtmosPhyCp), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    ! select case ( this%CP_TYPEID )
    ! case( CP_TYPEID_LSCOND )
    ! end select

    call this%vars%Final()
    return
  end subroutine AtmosPhyCp_finalize

!- private ------------------------------------------------

end module mod_atmos_phy_cp
