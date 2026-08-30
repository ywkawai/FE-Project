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

  use scale_element_line, only: LineElement

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

    real(RP) :: dtsec !< Timestep for cumulus parameterization

    type(LineElement) :: v_elem1D
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
  integer, parameter :: CP_TYPEID_MCONV_ADJUSTMENT = 1 !< Type ID of a moist convective adjustment scheme
  
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

    use scale_atm_phy_cp_dgm_mconv_adjustment, only: &
      atm_phy_cp_dgm_mconv_adjustment_setup
    implicit none
    class(AtmosPhyCp), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8 !< Timestep for cumulus parameterization
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  !< Unit of timestep

    character(len=H_MID) :: CP_TYPE = 'NONE'                 !< Type of a cumulus parameterization scheme
    integer :: atm_var_container_typeid                      !< Type ID of variable container for cumulus parameterization

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
 
    this%atm_var_container_typeid = atm_var_container_typeid
    
    !- Get atmospheric mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select

    !--- Register this component in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_CP', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                       ! (out) 

    this%dtsec = tm_parent_comp%process_list(this%tm_process_id)%dtsec

    !--- Set the type of cumulus parameterization

    select case( CP_TYPE )
    case( 'MOIST_CONV_ADJUSTMENT' )
      this%CP_TYPEID = CP_TYPEID_MCONV_ADJUSTMENT
      call atm_phy_cp_dgm_mconv_adjustment_setup()
    case default
      LOG_ERROR("ATMOS_PHY_CP_setup",*) 'Not appropriate cumulus parameterization type. Check!'
      call PRC_abort
    end select

    !- Initialize the variables 
    call this%vars%Init( model_mesh )

    !-
    call this%v_elem1D%Init( atm_mesh%ptr_mesh%refElem3D%PolyOrder_v, .false. ) 

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
    use scale_tracer, only: &
      QA
    use scale_atm_phy_cp_dgm_mconv_adjustment, only: &
      atm_phy_cp_dgm_mconv_adjustment_calc_tendency
    
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,     &
      AtmosVars_GetLocalMeshPhyAuxVars,  &
      AtmosVars_GetLocalMeshQTRCVarList, & 
      AtmosVars_GetLocalMeshPhyTends,    &
      AtmosVars_GetLocalMeshQTRC_Qv
    use mod_atmos_phy_cp_vars, only: &
      AtmosPhyCpVars_GetLocalMeshFields_tend, &
      SFLX_RAIN_ID => ATMOS_PHY_CP_AUX2D_SFLX_RAIN_ID, &
      SFLX_ENGI_ID => ATMOS_PHY_CP_AUX2D_SFLX_ENGI_ID
    implicit none
    class(AtmosPhyCp), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    class(MeshBase), pointer :: mesh
    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lcmesh

    integer :: n
    integer :: ke
    integer :: iq

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT, QV
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: PRES, PT

    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_P
    class(LocalMeshFieldBase), pointer :: RHOQ_tp
    class(LocalMeshFieldBase), pointer :: cp_DENS_t, cp_RHOT_t, cp_RHOQv_t
    !------------------------------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_PROGRESS(*) 'atmosphere / physics / cumulus parameterization' 

    call model_mesh%GetModelMesh( mesh )
    select type(mesh)
    class is (MeshBase3D)
      mesh3D => mesh
    end select

    !-
    if ( is_update ) then
      call PROF_rapstart( 'ATM_CP_tendency', 2)
      do n=1, mesh3D%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshPrgVars( n,    &
          mesh, prgvars_list, auxvars_list,       &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,         &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
          lcmesh                                  )
        
        call AtmosVars_GetLocalMeshPhyAuxVars( n,  &
          mesh, auxvars_list,                      &
          PRES, PT )

        call AtmosVars_GetLocalMeshQTRC_Qv( n, &
          mesh, trcvars_list, forcing_list,    &
          QV )

        call AtmosPhyCpVars_GetLocalMeshFields_tend( n, &
          mesh, this%vars%tends_manager,            &
          cp_DENS_t, cp_RHOT_t, cp_RHOQv_t          )

        select case( this%CP_TYPEID )
        case( CP_TYPEID_MCONV_ADJUSTMENT )
          call atm_phy_cp_dgm_mconv_adjustment_calc_tendency( &
            cp_DENS_t%val, cp_RHOT_t%val, cp_RHOQv_t%val,     & ! (out)
            this%vars%auxvars2D(SFLX_RAIN_ID)%local(n)%val,   & ! (out)
            this%vars%auxvars2D(SFLX_ENGI_ID)%local(n)%val,   & ! (out)
            DDENS%val, DRHOT%val, QV%val, PT%val, PRES%val,   & ! (in)
            DENS_hyd%val, Rtot%val, CPtot%val, this%dtsec,    & ! (in)
            lcmesh, lcmesh%refElem3D, this%v_elem1D           ) ! (in)
        end select
      end do

      call PROF_rapend( 'ATM_CP_tendency', 2)
    end if

    call PROF_rapstart('ATM_PHY_CP_add_tend', 2)
    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPhyTends( n,        &
        mesh, forcing_list,                          &
        DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, &
        RHOH_p )

      call AtmosVars_GetLocalMeshQTRC_Qv( n, &
        mesh, trcvars_list, forcing_list,    &
        QV, RHOQ_tp )
      
      call AtmosPhyCpVars_GetLocalMeshFields_tend( n,  &
        mesh, this%vars%tends_manager,                 &
        cp_DENS_t, cp_RHOT_t, cp_RHOQv_t,              &
        lcmesh                                         )
      
      !$omp parallel private(ke)
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
        DENS_tp%val(:,ke) = DENS_tp%val(:,ke) + cp_DENS_t%val(:,ke)
        RHOT_tp%val(:,ke) = RHOT_tp%val(:,ke) + cp_RHOT_t%val(:,ke)
        RHOQ_tp%val(:,ke) = RHOQ_tp%val(:,ke) + cp_RHOQv_t%val(:,ke)
      end do
      !$omp end parallel
    end do
    call PROF_rapend('ATM_PHY_CP_add_tend', 2)

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
    use scale_atm_phy_cp_dgm_mconv_adjustment, only: &
      atm_phy_cp_dgm_mconv_adjustment_finalize
    implicit none
    class(AtmosPhyCp), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    select case ( this%CP_TYPEID )
    case( CP_TYPEID_MCONV_ADJUSTMENT )
      call atm_phy_cp_dgm_mconv_adjustment_finalize()
    end select

    call this%vars%Final()
    call this%v_elem1D%Final()
    return
  end subroutine AtmosPhyCp_finalize

!- private ------------------------------------------------

end module mod_atmos_phy_cp
