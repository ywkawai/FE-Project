!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / cloud macrophysics component
!!
!! @par Description
!!          Module for cloud macrophysics component
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_mac
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

  use mod_atmos_phy_mac_vars, only: AtmosPhyMacVars

  use mod_atmos_vars_container, only: &
    AtmosVarsContainer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage a component of cloud macrophysics model in atmospheric model
  !!
  type, extends(ModelComponentProc), public :: AtmosPhyMac
    integer :: mac_TYPEID         !< Type id of cloud macrophysics model
    type(AtmosPhyMacVars) :: vars !< Object to manage variables with cloud macrophysics model

    integer :: atm_var_container_typeid     !< Type ID of variable container for cloud macrophysics model

    logical :: do_precipitation    !< Apply sedimentation (precipitation)?  
    logical :: evap_precip         !< Apply evaporation of precipitation?
    real(RP) :: tau_evap_precip    !< Timescale for evaporation of precipitation
    real(RP) :: dtsec !< Timestep for cloud macrophysics model

    type(LineElement) :: v_elem1D

    logical :: CP_flag
    class(MeshFieldBase), pointer :: ptr_CP_RHOT_t
  contains
    procedure :: setup => AtmosPhyMac_setup
    procedure :: calc_tendency => AtmosPhyMac_calc_tendency
    procedure :: update => AtmosPhyMac_update
    procedure :: finalize => AtmosPhyMac_finalize
    procedure, public :: Set_CP_tends_manager => AtmosPhyMac_set_CP_tends_manager
    procedure, private :: calc_tendency_core => AtmosPhyMac_calc_tendency_core
  end type AtmosPhyMac

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
  integer, parameter :: mac_TYPEID_SUNDQVIST = 1 !< Type ID of the Sundqvist-type cloud model
  
contains

  !> Setup a component of cloud model in atmospheric model
  !!
  !! @param model_mesh Object to manage computational mesh of atmospheric model 
  !! @param tm_parent_comp Object to manage a temporal scheme in a parent component
  !!
  subroutine AtmosPhyMac_setup( this, model_mesh, tm_parent_comp )
    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_regist
    use mod_atmos_mesh, only: AtmosMesh
    use scale_time_manager, only: TIME_manager_component
    use mod_atmos_vars, only: ATM_VARS_CONTAINER_PRIMARY_ID

    use scale_atm_phy_mac_sundqvist, only: &
      ATMOS_PHY_mac_sundqvist_setup, &
      ATMOS_PHY_mac_sundqvist_ntracers,            &
      ATMOS_PHY_mac_sundqvist_nwaters,             &
      ATMOS_PHY_mac_sundqvist_nices,               &
      ATMOS_PHY_mac_sundqvist_tracer_names,        &
      ATMOS_PHY_mac_sundqvist_tracer_descriptions, &
      ATMOS_PHY_mac_sundqvist_tracer_units


    implicit none
    class(AtmosPhyMac), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8 !< Timestep for cloud model
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  !< Unit of timestep

    character(len=H_MID) :: mac_TYPE = 'NONE'                !< Type of a cloud model scheme

    logical :: do_precipitation     !< Flag whether sedimentation (precipitation) is applied
    logical :: evap_precip          !< Apply evaporation of precipitation?
    real(RP) :: tau_evap_precip     !< Timescale for evaporation of precipitation [s]

    integer :: atm_var_container_typeid !< Type ID of variable container for mac model

    namelist /PARAM_ATMOS_PHY_MAC/ &
      TIME_DT,                &
      TIME_DT_UNIT,           &
      mac_TYPE,               &
      do_precipitation,       &
      evap_precip,            &
      tau_evap_precip,        &
      atm_var_container_typeid
    
    integer :: ierr

    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D

    integer :: QS_mac, QE_mac, QA_mac
    !-----------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_phy_mac_setup",*) 'Setup'

    do_precipitation = .true.
    evap_precip      = .true.
    tau_evap_precip  = 3600.0_RP

    atm_var_container_typeid = ATM_VARS_CONTAINER_PRIMARY_ID

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MAC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_phy_mac_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_phy_mac_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MAC. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MAC)
 
    this%do_precipitation = do_precipitation
    this%evap_precip      = evap_precip
    this%tau_evap_precip  = tau_evap_precip

    this%atm_var_container_typeid = atm_var_container_typeid
    
    !- Get atmospheric mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select

    !--- Register this component in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_phy_mac', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                       ! (out) 

    this%dtsec = tm_parent_comp%process_list(this%tm_process_id)%dtsec

    !--- Set the type of mac model

    select case( mac_TYPE )
    case( 'SUNDQVIST' )
      this%mac_TYPEID = mac_TYPEID_SUNDQVIST

      call ATMOS_HYDROMETEOR_regist( &
           ATMOS_PHY_mac_SUNDQVIST_nwaters,                & ! (in)
           ATMOS_PHY_mac_SUNDQVIST_nices,                  & ! (in)
           ATMOS_PHY_mac_SUNDQVIST_tracer_names(:),        & ! (in)
           ATMOS_PHY_mac_SUNDQVIST_tracer_descriptions(:), & ! (in)
           ATMOS_PHY_mac_SUNDQVIST_tracer_units(:),        & ! (in)
           QS_mac                                          ) ! (out)

      QA_mac = ATMOS_PHY_mac_SUNDQVIST_ntracers
    case default
      LOG_ERROR("ATMOS_phy_mac_setup",*) 'Not appropriate mac model type. Check!'
      call PRC_abort
    end select

    this%do_precipitation  = do_precipitation

    this%atm_var_container_typeid = atm_var_container_typeid

    !- Initialize the variables
    QE_mac = QS_mac + QA_mac - 1
    call this%vars%Init( model_mesh, QS_mac, QE_mac, QA_mac )

    !-
    call this%v_elem1D%Init( atm_mesh%ptr_mesh%refElem3D%PolyOrder_v, .false. ) 

    !- Setup a module for cloud modules

    select case( this%mac_TYPEID )
    case( mac_TYPEID_SUNDQVIST )
      call ATMOS_PHY_mac_SUNDQVIST_setup()
    end select

    return
  end subroutine AtmosPhyMac_setup

  subroutine AtmosPhyMac_set_CP_tends_manager( this, CP_tends_manager )
    use mod_atmos_phy_cp_vars, only: &
      CP_RHOT_t_ID => ATMOS_PHY_CP_RHOT_t_ID    
    implicit none
    class(AtmosPhyMac), intent(inout) :: this
    class(ModelVarManager), intent(inout) :: CP_tends_manager
    !--------------------------------------------------

    this%CP_flag = .true.
    call CP_tends_manager%Get( CP_RHOT_t_ID, this%ptr_CP_RHOT_t )
    return
  end subroutine AtmosPhyMac_set_CP_tends_manager

  !> Calculate tendencies associated with cloud macrophysics model in atmospheric model
  !!
  !!
  !! @param model_mesh Object to manage computational mesh of atmospheric model 
  !! @param prgvars_list Object to manage prognostic variables with atmospheric dynamical core
  !! @param trcvars_list Object to manage auxiliary variables 
  !! @param forcing_list Object to manage forcing terms
  !! @param is_update Flag to specify whether the tendencies are updated in this call
  !!
!OCL SERIAL
  subroutine AtmosPhyMac_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )
    use scale_tracer, only: &
      QA
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,     &
      AtmosVars_GetLocalMeshPhyAuxVars,  &
      AtmosVars_GetLocalMeshQTRCVarList, & 
      AtmosVars_GetLocalMeshPhyTends
    use mod_atmos_phy_mac_vars, only: &
      AtmosPhyMacVars_GetLocalMeshFields_tend,   &
      AtmosPhyMacVars_GetLocalMeshFields_sfcflx, &
      SFLX_RAIN_ID => ATMOS_phy_mac_AUX2D_SFLX_RAIN_ID, &
      SFLX_ENGI_ID => ATMOS_phy_mac_AUX2D_SFLX_ENGI_ID
    implicit none
    class(AtmosPhyMac), intent(inout) :: this
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

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    type(LocalMeshFieldBaseList) :: QTRC(this%vars%QS:this%vars%QE)
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT


    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_P
    type(LocalMeshFieldBaseList) :: RHOQ_tp(QA)
    class(LocalMeshFieldBase), pointer :: mac_DENS_t, mac_MOMX_t, mac_MOMY_t, mac_MOMZ_t, mac_RHOT_t, mac_RHOH, mac_EVAP
    type(LocalMeshFieldBaseList) :: mac_RHOQ_t(this%vars%QS:this%vars%QE)
    class(LocalMeshFieldBase), pointer :: SFLX_rain, SFLX_snow, SFLX_engi

    class(LocalMeshFieldBase), pointer :: CP_RHOT_t
    logical, allocatable :: CP_mask(:,:)
    !------------------------------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_PROGRESS(*) 'atmosphere / physics / cloud model' 

    call model_mesh%GetModelMesh( mesh )
    select type(mesh)
    class is (MeshBase3D)
      mesh3D => mesh
    end select

    !-
    if ( is_update ) then
      call PROF_rapstart( 'ATM_mac_tendency', 2)

      do n=1, mesh3D%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshPrgVars( n,    &
          mesh, prgvars_list, auxvars_list,       &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,         &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
          lcmesh                                  )
        
        call AtmosVars_GetLocalMeshQTRCVarList( n,  &
          mesh, trcvars_list, this%vars%QS, QTRC(:) )

        call AtmosVars_GetLocalMeshPhyAuxVars( n,  &
          mesh, auxvars_list,                      &
          PRES, PT )

        call AtmosPhyMacVars_GetLocalMeshFields_tend( n, &
          mesh, this%vars%tends_manager,                                                  &
          mac_DENS_t, mac_MOMX_t, mac_MOMY_t, mac_MOMZ_t, mac_RHOT_t, mac_RHOH, mac_EVAP, &
          mac_RHOQ_t, lcmesh                                                              )

        call AtmosPhyMacVars_GetLocalMeshFields_sfcflx( n, &
          mesh, this%vars%auxvars2D_manager,               &
          SFLX_rain, SFLX_snow, SFLX_engi                  )     

        !-
        allocate( CP_mask(lcmesh%refElem3D%Np,lcmesh%Ne) )
        if ( this%CP_flag ) then
          call this%ptr_CP_RHOT_t%GetLocalMeshField( n, CP_RHOT_t )
          where ( abs(CP_RHOT_t%val(:,lcmesh%NeS:lcmesh%NeE)) > 1.0E-12_RP )
            CP_mask(:,:) = .true.
          elsewhere
            CP_mask(:,:) = .false.
          end where
        end if

        call this%calc_tendency_core( &
          mac_DENS_t%val, mac_MOMX_t%val, mac_MOMY_t%val, mac_MOMZ_t%val, mac_RHOQ_t,  & ! (out)
          mac_RHOH%val, mac_EVAP%val, SFLX_rain%val, SFLX_snow%val, SFLX_engi%val,     & ! (out)
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, PT%val, QTRC, PRES%val,             & ! (in)
          PRES_hyd%val, DENS_hyd%val, Rtot%val, CVtot%val, CPtot%val,                  & ! (in)
          CP_mask, model_mesh%DOptrMat(3), model_mesh%LiftOptrMat,                     & ! (in)
          lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D,        & ! (in)
          this%v_elem1D ) ! (in)

        deallocate( CP_mask )
      end do

      call PROF_rapend( 'ATM_mac_tendency', 2)
    end if

    call PROF_rapstart('ATM_phy_mac_add_tend', 2)
    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPhyTends( n,        &
        mesh, forcing_list,                          &
        DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, &
        RHOH_p, RHOQ_tp  )

      call AtmosPhyMacVars_GetLocalMeshFields_tend( n, &
        mesh, this%vars%tends_manager,                                                  &
        mac_DENS_t, mac_MOMX_t, mac_MOMY_t, mac_MOMZ_t, mac_RHOT_t, mac_RHOH, mac_EVAP, &
        mac_RHOQ_t, lcmesh                                                              )
      
      !$omp parallel private(ke, iq)
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
        DENS_tp%val(:,ke) = DENS_tp%val(:,ke) + mac_DENS_t%val(:,ke)
        MOMX_tp%val(:,ke) = MOMX_tp%val(:,ke) + mac_MOMX_t%val(:,ke)
        MOMY_tp%val(:,ke) = MOMY_tp%val(:,ke) + mac_MOMY_t%val(:,ke)
        MOMZ_tp%val(:,ke) = MOMZ_tp%val(:,ke) + mac_MOMZ_t%val(:,ke)
        RHOH_p %val(:,ke) = RHOH_p %val(:,ke) + mac_RHOH  %val(:,ke)
      end do
      !$omp end do
      !$omp do collapse(2)
      do iq = this%vars%QS, this%vars%QE
      do ke = lcmesh%NeS, lcmesh%NeE
        RHOQ_tp(iq)%ptr%val(:,ke) = RHOQ_tp(iq)%ptr%val(:,ke)  &
                                  + mac_RHOQ_t(iq)%ptr%val(:,ke)
      end do
      end do 
      !$omp end do
      !$omp end parallel
    end do
    call PROF_rapend('ATM_phy_mac_add_tend', 2)

    return
  end subroutine AtmosPhyMac_calc_tendency

  !> Update variables in a component of mac model in atmospheric model
  !!
  !! @param model_mesh Object to manage computational mesh of atmospheric model 
  !! @param prgvars_list Object to manage prognostic variables with atmospheric dynamical core
  !! @param trcvars_list Object to manage auxiliary variables 
  !! @param forcing_list Object to manage forcing terms
  !! @param is_update Flag to speicfy whether the tendencies are updated in this call
  !!
!OCL SERIAL  
  subroutine AtmosPhyMac_update( this, model_mesh, &
    prgvars_list, trcvars_list,                   &
    auxvars_list, forcing_list, is_update         )  
    
    implicit none
    class(AtmosPhyMac), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------
    return
  end subroutine AtmosPhyMac_update

!> Finalize a component of mac model in atmospheric model
!!
!OCL SERIAL  
  subroutine AtmosPhyMac_finalize( this )
    use scale_atm_phy_mac_sundqvist, only: &
      ATMOS_PHY_mac_sundqvist_finalize
    implicit none
    class(AtmosPhyMac), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    select case ( this%mac_TYPEID )
    case( mac_TYPEID_SUNDQVIST )
      call ATMOS_PHY_mac_sundqvist_finalize()
    end select

    call this%vars%Final()
    call this%v_elem1D%Final()
    return
  end subroutine AtmosPhyMac_finalize

!- private ------------------------------------------------

  !> Calculate tendencies associated with cloud microphysics and precipitation for each local mesh
!OCL SERIAL
  subroutine AtmosPhyMac_calc_tendency_core( this, &
    DENS_t_mac, RHOU_t_mac, RHOV_t_mac, MOMZ_t_mac, RHOQ_t_mac,  & ! (out)
    RHOH_mac, EVAPORATE, SFLX_rain, SFLX_snow, SFLX_ENGI,   & ! (out)
    DDENS, RHOU, RHOV, MOMZ, PT, QTRC,                      & ! (in)
    PRES, PRES_hyd, DENS_hyd,                               & ! (in)
    Rtot, CVtot, CPtot,                                     & ! (in)
    CP_mask,                                                & ! (in)
    Dz, Lift,                                               & ! (in)
    lcmesh, elem3D, lcmesh2D, elem2D, elem_v1D              ) ! (in)

    use scale_const, only: &
      PRE00 => CONST_PRE00
    use scale_const, only: &
      CVdry => CONST_CVdry, &
      CPdry => CONST_CPdry
    use scale_atmos_hydrometeor, only: & 
      LHF,                             &
      QHA, QLA, QIA
    use scale_sparsemat, only: SparseMat

    use scale_atm_phy_cloud_dgm_common, only: &
      atm_phy_cloud_dgm_common_condensate_removal,         &
      atm_phy_cloud_dgm_common_condensate_removal_momentum
    implicit none

    class(AtmosPhyMac), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    class(ElementBase1D), intent(in) :: elem_v1D
    real(RP), intent(out) :: DENS_t_mac(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: RHOU_t_mac(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: RHOV_t_mac(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMZ_t_mac(elem3D%Np,lcmesh%NeA)
    type(LocalMeshFieldBaseList), intent(inout) :: RHOQ_t_mac(this%vars%QS:this%vars%QE)
    real(RP), intent(out) :: RHOH_mac(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: EVAPORATE(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: SFLX_rain(elem2D%Np,lcmesh2D%NeA) 
    real(RP), intent(out) :: SFLX_snow(elem2D%Np,lcmesh2D%NeA) 
    real(RP), intent(out) :: SFLX_ENGI(elem2D%Np,lcmesh2D%NeA) 
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: RHOU(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: RHOV(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PT  (elem3D%Np,lcmesh%NeA)
    type(LocalMeshFieldBaseList), intent(in) :: QTRC(this%vars%QS:this%vars%QE)
    real(RP), intent(in) :: PRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot (elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA)
    logical, intent(in) :: CP_mask(elem3D%Np,lcmesh%Ne)
    class(SparseMat), intent(in) :: Dz
    class(SparseMat), intent(in) :: Lift

    real(RP) :: DENS (elem3D%Np,lcmesh%NeA)
    real(RP) :: RHOE_t(elem3D%Np,lcmesh%NeA)

    real(RP) :: CPtot_t(elem3D%Np,lcmesh%NeA)
    real(RP) :: CVtot_t(elem3D%Np,lcmesh%NeA)

    real(RP) :: RHOQ(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D,this%vars%QS+1:this%vars%QE)
    real(RP) :: RHOQ2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D,this%vars%QS+1:this%vars%QE)
    real(RP) :: RHOQV2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: DENS0(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: DENS2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOU2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOV2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: MOMZ2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: TEMP2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: CPtot2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: CVtot2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOE (elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOE2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)

    real(RP) :: CP_t(elem3D%Np), CV_t(elem3D%Np)

    integer :: iq
    integer :: iq_QV
    integer :: domid
    integer :: ke
    integer :: ke2D, ke_z

    real(RP) :: rdt_mac

    integer :: QS_remove
    !--------------------------------------------------

    rdt_mac = 1.0_RP / this%dtsec

    iq_QV = this%vars%QS

    !$omp parallel do private(ke)
    do ke = lcmesh%NeS, lcmesh%NeE
      DENS(:,ke) = DENS_hyd(:,ke) + DDENS(:,ke)
    end do    

    !- Calculate tendencies of cloud processes ----------------------

    select case( this%mac_TYPEID )
    case( mac_TYPEID_SUNDQVIST )
      call calc_tendency_Sundqvist( this, &
        RHOQ_t_mac, CPtot_t, CVtot_t, RHOE_t, EVAPORATE, & ! (out)
        DENS, QTRC, PRES, DENS_hyd, Rtot, CVtot, CPtot,  & ! (in)
        CP_mask, rdt_mac, lcmesh, elem3D )                 ! (in)
    end select

    !$omp parallel do
    do ke = lcmesh%NeS, lcmesh%NeE
      RHOH_mac(:,ke) = RHOE_t(:,ke) &
        - ( CPtot_t(:,ke) + log( PRES(:,ke) / PRE00 ) * ( CVtot(:,ke) / CPtot(:,ke) * CPtot(:,ke) - CVtot(:,ke) ) ) &
        * PRES(:,ke) / Rtot(:,ke)
    end do

    !- Calculate precipitation processes if enabled ----------------------

    if ( this%do_precipitation ) then
      !$omp parallel private(ke,ke2D,ke_z,iq)
      !$omp do collapse(2)
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D
        DENS0(:,ke_z,ke2D) = DENS_hyd(:,ke) + DDENS(:,ke)
        DENS2(:,ke_z,ke2D) = DENS_hyd(:,ke) + DDENS(:,ke)

        RHOU2(:,ke_z,ke2D) = RHOU(:,ke)
        RHOV2(:,ke_z,ke2D) = RHOV(:,ke)
        MOMZ2(:,ke_z,ke2D) = MOMZ(:,ke)

        TEMP2(:,ke_z,ke2D) = PRES(:,ke) / ( DENS0(:,ke_z,ke2D) * Rtot(:,ke) )
        CPtot2(:,ke_z,ke2D) = CPtot(:,ke)
        CVtot2(:,ke_z,ke2D) = CVtot(:,ke)
        RHOE(:,ke_z,ke2D) = TEMP2(:,ke_z,ke2D) * DENS0(:,ke_z,ke2D) * CVtot(:,ke)
        RHOE2(:,ke_z,ke2D) = RHOE(:,ke_z,ke2D)

        RHOQV2(:,ke_z,ke2D) = DENS(:,ke) * QTRC(iq_QV)%ptr%val(:,ke) &
                            + RHOQ_t_mac(iq_QV)%ptr%val(:,ke) * this%dtsec
      end do
      end do
      !$omp end do
      !$omp do collapse(3)
      do iq = this%vars%QS+1, this%vars%QE
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D
        RHOQ(:,ke_z,ke2D,iq) = DENS(:,ke) * QTRC(iq)%ptr%val(:,ke) &
                             + RHOQ_t_mac(iq)%ptr%val(:,ke) * this%dtsec
        RHOQ2(:,ke_z,ke2D,iq) = RHOQ(:,ke_z,ke2D,iq)
      end do
      end do
      end do
      !$omp end do
      !$omp workshare
      SFLX_rain(:,:) = 0.0_RP
      SFLX_snow(:,:) = 0.0_RP
      SFLX_ENGI(:,:) = 0.0_RP
      !$omp end workshare      
      !$omp end parallel

      !- Remove condensation water

      select case( this%mac_TYPEID )
      case( mac_TYPEID_SUNDQVIST )
        QS_remove = 2
      end select

      call atm_phy_cloud_dgm_common_condensate_removal( &
        DENS2, RHOQ2, CPtot2, CVtot2, RHOE2,              & ! (inout)
        SFLX_rain, SFLX_snow, SFLX_ENGI,                  & ! (inout)
        TEMP2, this%dtsec,                                & ! (in)
        this%vars%QE - this%vars%QS, QLA, QIA, QS_remove, & ! (in)
        lcmesh, elem3D, elem_v1D,                         & ! (in)
        this%evap_precip, this%tau_evap_precip,           & ! (in)
        RHOQV2                                            ) ! (inout)

      !$omp parallel private(ke2D, ke_z, iq, ke, CP_t, CV_t)
      !$omp workshare
      SFLX_ENGI(:,:) = SFLX_ENGI(:,:) - SFLX_snow(:,:) * LHF ! moist internal energy
      !$omp end workshare
      !$omp do collapse(2)
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D
        DENS_t_mac(:,ke) = ( DENS2(:,ke_z,ke2D) - DENS0(:,ke_z,ke2D) ) * rdt_mac

        CP_t(:) = ( CPtot2(:,ke_z,ke2D) - CPtot(:,ke) ) * rdt_mac
        CV_t(:) = ( CVtot2(:,ke_z,ke2D) - CVtot(:,ke) ) * rdt_mac
        RHOH_mac(:,ke) = RHOH_mac(:,ke) &
          + ( RHOE2(:,ke_z,ke2D) - RHOE(:,ke_z,ke2D) ) * rdt_mac &
          - ( CP_t(:) &
              + log( PRES(:,ke) / PRE00 ) * ( CVtot(:,ke) / CPtot(:,ke) * CP_t(:) - CV_t(:) ) &
            ) * PRES(:,ke) / Rtot(:,ke)

        if ( this%evap_precip ) then
          RHOQ_t_mac(iq_QV)%ptr%val(:,ke) = RHOQ_t_mac(iq_QV)%ptr%val(:,ke) &
            + ( RHOQV2(:,ke_z,ke2D) - DENS(:,ke) * QTRC(iq_QV)%ptr%val(:,ke) ) * rdt_mac
        end if
      end do
      end do
      !$omp end do
      !$omp do collapse(2)
      do iq = this%vars%QS+1, this%vars%QE
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D
        RHOQ_t_mac(iq)%ptr%val(:,ke) = RHOQ_t_mac(iq)%ptr%val(:,ke)       &
          + ( RHOQ2(:,ke_z,ke2D,iq) - RHOQ(:,ke_z,ke2D,iq) ) * rdt_mac
      end do
      end do
      end do
      !$omp end do
      !$omp end parallel

      !- Sedimentation of momentum

      call atm_phy_cloud_dgm_common_condensate_removal_momentum( &
        RHOU_t_mac, RHOV_t_mac, MOMZ_t_mac,                    & ! (out)
        DENS0, RHOU2, RHOV2, MOMZ2, DENS2,                     & ! (in)
        rdt_mac, lcmesh, elem3D                                ) ! (in)
    end if

    return
  end subroutine AtmosPhyMac_calc_tendency_core


  !> Calculate tendencies of cloud microphysics processes with Kessler scheme
  !!
!OCL SERIAL
  subroutine calc_tendency_Sundqvist( this, &
    RHOQ_t_mac, CPtot_t, CVtot_t, RHOE_t, EVAPORATE, & ! (out)
    DENS, QTRC, PRES, DENS_hyd, Rtot, CVtot, CPtot,  & ! (in)
    CP_mask,                                         & ! (in)
    rdt_mac, lcmesh, elem3D )                          ! (in)

    use scale_atmos_hydrometeor, only: &
       LHV    
    use scale_atm_phy_mac_sundqvist, only: &
      ATMOS_PHY_mac_sundqvist_adjustment
    implicit none

    class(AtmosPhyMac), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    type(LocalMeshFieldBaseList), intent(inout) :: RHOQ_t_mac(this%vars%QS:this%vars%QE)
    real(RP), intent(out) :: CPtot_t(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: CVtot_t(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: RHOE_t(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: EVAPORATE(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS(elem3D%Np,lcmesh%NeA)
    type(LocalMeshFieldBaseList), intent(in) :: QTRC(this%vars%QS:this%vars%QE)
    real(RP), intent(in) :: PRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot (elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA)
    logical, intent(in) :: CP_mask(elem3D%Np,lcmesh%Ne)
    real(RP), intent(in) :: rdt_mac

    real(RP) :: TEMP1(elem3D%Np,lcmesh%NeA)
    real(RP) :: CPtot1(elem3D%Np,lcmesh%NeA)
    real(RP) :: CVtot1(elem3D%Np,lcmesh%NeA)
    real(RP) :: QTRC1(elem3D%Np,lcmesh%NeA,this%vars%QS:this%vars%QE)

    integer :: ke
    integer :: iq

    real(RP) :: RHOQ_t(elem3D%Np)
    real(RP) :: RHOQ_pri(elem3D%Np)
    real(RP) :: RHOQ_t_cor(elem3D%Np)
    real(RP) :: RHOQV_t(elem3D%Np,lcmesh%Ne)

    real(RP) :: QR_tmp(elem3D%Np,lcmesh%NeA)

    real(RP) :: sw(elem3D%Np,lcmesh%Ne)
    !------------------------------------------------------------

    !$omp parallel private(ke, iq)
    !$omp do
    do ke = lcmesh%NeS, lcmesh%NeE
      TEMP1(:,ke) = PRES(:,ke) / ( DENS(:,ke) * Rtot(:,ke) )
      CPtot1(:,ke) = CPtot(:,ke)
      CVtot1(:,ke) = CVtot(:,ke)
      RHOQV_t(:,ke) = 0.0_RP

      sw(:,ke) = 1.0_RP
    end do
    if ( this%CP_flag ) then
      !$omp do
      do ke = lcmesh%NeS, lcmesh%NeE
        where (CP_mask(:,ke))
          sw(:,ke) = 0.0_RP
        end where
      end do
    end if
    !$omp do collapse(2)
    do iq = this%vars%QS, this%vars%QE
    do ke = lcmesh%NeS, lcmesh%NeE
      QTRC1(:,ke,iq) = QTRC(iq)%ptr%val(:,ke)
    end do
    end do
    !$omp end parallel

    call ATMOS_PHY_mac_sundqvist_adjustment( &
      elem3D%Np, 1, elem3D%Np, lcmesh%NeA, lcmesh%NeS, lcmesh%NeE, 1, 1, 1,  & ! (in)
      DENS, PRES, this%dtsec,                                                & ! (in)
      TEMP1, QTRC1, CPtot1, CVtot1,                                          & ! (inout)
      EVAPORATE, RHOE_t                                                      ) ! (out)
  
    !$omp parallel private(ke, iq, RHOQ_t, RHOQ_t_cor, RHOQ_pri)
    !$omp do collapse(2)
    do ke = lcmesh%NeS, lcmesh%NeE
    do iq = this%vars%QS+1, this%vars%QE      
      RHOQ_t(:) = ( QTRC1(:,ke,iq) - QTRC(iq)%ptr%val(:,ke) ) * DENS(:,ke) * rdt_mac

      RHOQ_pri(:) = DENS(:,ke) * QTRC(iq)%ptr%val(:,ke)
      RHOQ_t_cor(:) = max( RHOQ_t(:), - RHOQ_pri(:) * rdt_mac ) * sw(:,ke)

      RHOQV_t(:,ke) = RHOQV_t(:,ke) - RHOQ_t_cor(:)
      RHOE_t(:,ke) = ( RHOE_t(:,ke) + LHV * ( RHOQ_t_cor(:) - RHOQ_t(:) ) ) * sw(:,ke)
      RHOQ_t_mac(iq)%ptr%val(:,ke) = RHOQ_t_cor(:)
    end do  
    end do
    !$omp end do
    !$omp do
    do ke = lcmesh%NeS, lcmesh%NeE
      RHOQ_t_mac(this%vars%QS)%ptr%val(:,ke) = RHOQV_t(:,ke)

      CPtot_t(:,ke) = ( CPtot1(:,ke) - CPtot(:,ke) ) * rdt_mac * sw(:,ke)
      CVtot_t(:,ke) = ( CVtot1(:,ke) - CVtot(:,ke) ) * rdt_mac * sw(:,ke)
    end do
    !$omp end do
    !$omp end parallel
    return
  end subroutine calc_tendency_Sundqvist

end module mod_atmos_phy_mac
