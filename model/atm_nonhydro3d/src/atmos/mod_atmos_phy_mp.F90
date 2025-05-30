!-------------------------------------------------------------------------------
!> module ATMOSPHERE physics / cloud microphysics
!!
!! @par Description
!!          Module for cloud microphysics
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_mp
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
  use scale_tracer, only: QA

  use scale_sparsemat, only: sparsemat
  use scale_element_hexahedral, only: HexahedralElement

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
  
  use mod_atmos_phy_mp_vars, only: AtmosPhyMpVars


  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage a component of cloud microphysics
  !!
  type, extends(ModelComponentProc), public :: AtmosPhyMp
    integer :: MP_TYPEID         !< Type id of cloud microphysics scheme
    type(AtmosPhyMpVars) :: vars !< A object to manage variables with cloud microphysics

    logical, private :: do_precipitation    !< Apply sedimentation (precipitation)?
    logical, private :: do_negative_fixer   !< Apply negative fixer?
    real(RP), private :: limit_negative     !< Abort if abs(fixed negative value) > abs(MP_limit_negative)
    integer, private :: ntmax_sedimentation !< Number of time step for sedimentation
    real(RP), private :: max_term_vel       !< Terminal velocity for calculate dt of sedimentation
    real(RP), private :: cldfrac_thleshold  !< Threshold for cloud fraction

    real(DP), private :: dtsec
    integer, private :: nstep_sedmientation
    real(RP), private :: rnstep_sedmientation
    real(DP), private :: dtsec_sedmientation

    type(sparsemat) :: Dz, Lift
    type(HexahedralElement) :: elem
  contains
    procedure, public :: setup => AtmosPhyMp_setup 
    procedure, public :: calc_tendency => AtmosPhyMp_calc_tendency
    procedure, public :: update => AtmosPhyMp_update
    procedure, public :: finalize => AtmosPhyMp_finalize
    procedure, private :: calc_tendency_core => AtmosPhyMp_calc_tendency_core
  end type AtmosPhyMp

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------

  integer, parameter :: MP_TYPEID_KESSLER  = 1
  integer, parameter :: MP_TYPEID_TOMITA08 = 2
  integer, parameter :: MP_TYPEID_SN14     = 3

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
    
  logical :: Flag_LT = .false. ! Tentative

contains

!> Setup a component of cloud microphysics
!!
!! @param model_mesh a object to manage computational mesh of atmospheric model 
!! @param tm_parent_comp a object to mange a temporal scheme in a parent component
!!
  subroutine AtmosPhyMp_setup( this, model_mesh, tm_parent_comp )
    use scale_const, only: &
      EPS => CONST_EPS
    use scale_tracer, only: &
      TRACER_regist
    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_regist
    use scale_atmos_phy_mp_kessler, only: &
      ATMOS_PHY_MP_kessler_setup,               &
      ATMOS_PHY_MP_KESSLER_ntracers,            &
      ATMOS_PHY_MP_KESSLER_nwaters,             &
      ATMOS_PHY_MP_KESSLER_nices,               &
      ATMOS_PHY_MP_KESSLER_tracer_names,        &
      ATMOS_PHY_MP_KESSLER_tracer_descriptions, &
      ATMOS_PHY_MP_KESSLER_tracer_units
    use scale_atmos_phy_mp_tomita08, only: &
      ATMOS_PHY_MP_tomita08_setup,               &
      ATMOS_PHY_MP_tomita08_ntracers,            &
      ATMOS_PHY_MP_tomita08_nwaters,             &
      ATMOS_PHY_MP_tomita08_nices,               &
      ATMOS_PHY_MP_tomita08_tracer_names,        &
      ATMOS_PHY_MP_tomita08_tracer_descriptions, &
      ATMOS_PHY_MP_tomita08_tracer_units
    use scale_atmos_phy_mp_sn14, only: &
      ATMOS_PHY_MP_sn14_setup,               &
      ATMOS_PHY_MP_sn14_ntracers,            &
      ATMOS_PHY_MP_sn14_nwaters,             &
      ATMOS_PHY_MP_sn14_nices,               &
      ATMOS_PHY_MP_sn14_tracer_names,        &
      ATMOS_PHY_MP_sn14_tracer_descriptions, &
      ATMOS_PHY_MP_sn14_tracer_units      
    !------------------------------
    use scale_time_manager, only: TIME_manager_component    
    use mod_atmos_mesh, only: AtmosMesh

    implicit none
    class(AtmosPhyMp), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8 !< Timestep for cloud microcloud physics
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  !< Unit of timestep

    character(len=H_MID) :: MP_TYPE = 'KESSLER'              !< Type of cloud microcloud physics scheme

    logical :: do_precipitation     !< Flag whether sedimentation (precipitation) is applied
    logical :: do_negative_fixer    !< Flag whether negative fixer is applied
    real(RP) :: limit_negative
    integer :: ntmax_sedimentation  
    real(RP) :: max_term_vel        
    real(RP) :: cldfrac_thleshold

    namelist /PARAM_ATMOS_PHY_MP/ &
      TIME_DT,             &
      TIME_DT_UNIT,        &
      MP_TYPE,             &
      do_precipitation,    &
      do_negative_fixer,   &
      limit_negative,      &
      ntmax_sedimentation, &
      max_term_vel,        &
      cldfrac_thleshold
    
    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D

    integer :: ierr

    integer :: QS_MP, QE_MP, QA_MP
    integer :: QS2

    integer :: nstep_max
    !--------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_setup",*) 'Setup'

    cldfrac_thleshold = EPS

    do_precipitation    = .true.
    do_negative_fixer   = .true.
    limit_negative      = 0.1_RP
    ntmax_sedimentation = 1
    max_term_vel        = 10.0_RP
    
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_MP_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_PHY_MP_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP)

    this%cldfrac_thleshold = cldfrac_thleshold
    this%do_precipitation  = do_precipitation
    this%do_negative_fixer = do_negative_fixer
    this%ntmax_sedimentation = ntmax_sedimentation
    this%max_term_vel        = max_term_vel

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_setup",*) 'Enable negative fixer?                    : ', this%do_negative_fixer
    LOG_INFO("ATMOS_PHY_MP_setup",*) 'Value limit of negative fixer (abs)       : ', abs(this%limit_negative)
    LOG_INFO("ATMOS_PHY_MP_setup",*) 'Enable sedimentation (precipitation)?     : ', this%do_precipitation

    !- get mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select

    !--- Register this component in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_MP', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                       ! (out) 

    this%dtsec = tm_parent_comp%process_list(this%tm_process_id)%dtsec
    nstep_max = 0 ! ceiling( this%dtsec * this%max_term_vel / maxval( CDZ) )
    this%ntmax_sedimentation = max( ntmax_sedimentation, nstep_max )

    this%nstep_sedmientation = ntmax_sedimentation
    this%rnstep_sedmientation = 1.0_RP / real(ntmax_sedimentation,kind=RP)
    this%dtsec_sedmientation  = this%dtsec * this%rnstep_sedmientation
  
    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_setup",*) 'Enable negative fixer?                    : ', this%do_negative_fixer
    !LOG_INFO("ATMOS_PHY_MP_setup",*) 'Value limit of negative fixer (abs)       : ', abs(MP_limit_negative)
    LOG_INFO("ATMOS_PHY_MP_setup",*) 'Enable sedimentation (precipitation)?     : ', this%do_precipitation
    LOG_INFO("ATMOS_PHY_MP_setup",*) 'Timestep of sedimentation is divided into : ', this%ntmax_sedimentation, 'step'
    LOG_INFO("ATMOS_PHY_MP_setup",*) 'DT of sedimentation                       : ', this%dtsec_sedmientation, '[s]'

    !--- Set the type of microphysics scheme

    select case( MP_TYPE )
    case ('KESSLER')
      this%MP_TYPEID = MP_TYPEID_KESSLER
      call ATMOS_HYDROMETEOR_regist( &
           ATMOS_PHY_MP_KESSLER_nwaters,                & ! (in)
           ATMOS_PHY_MP_KESSLER_nices,                  & ! (in)
           ATMOS_PHY_MP_KESSLER_tracer_names(:),        & ! (in)
           ATMOS_PHY_MP_KESSLER_tracer_descriptions(:), & ! (in)
           ATMOS_PHY_MP_KESSLER_tracer_units(:),        & ! (in)
           QS_MP                                        ) ! (out)
      QA_MP = ATMOS_PHY_MP_KESSLER_ntracers
    case ( 'TOMITA08' )
      this%MP_TYPEID = MP_TYPEID_TOMITA08
      call ATMOS_HYDROMETEOR_regist( &
            ATMOS_PHY_MP_TOMITA08_nwaters,                & ! (in)
            ATMOS_PHY_MP_TOMITA08_nices,                  & ! (in)
            ATMOS_PHY_MP_TOMITA08_tracer_names(:),        & ! (in)
            ATMOS_PHY_MP_TOMITA08_tracer_descriptions(:), & ! (in)
            ATMOS_PHY_MP_TOMITA08_tracer_units(:),        & ! (in)
            QS_MP                                         ) ! (out)
      QA_MP = ATMOS_PHY_MP_TOMITA08_ntracers
    case( 'SN14' )
      this%MP_TYPEID = MP_TYPEID_SN14
      call ATMOS_HYDROMETEOR_regist( &
            ATMOS_PHY_MP_SN14_nwaters,                   & ! (in)
            ATMOS_PHY_MP_SN14_nices,                     & ! (in)
            ATMOS_PHY_MP_SN14_tracer_names(1:6),         & ! (in)
            ATMOS_PHY_MP_SN14_tracer_descriptions(1:6),  & ! (in)
            ATMOS_PHY_MP_SN14_tracer_units(1:6),         & ! (in)
            QS_MP                                        ) ! (out)
      
      call TRACER_regist( QS2,                           & ! (out)
            5,                                           & ! (in)
            ATMOS_PHY_MP_SN14_tracer_names(7:11),        & ! (in)
            ATMOS_PHY_MP_SN14_tracer_descriptions(7:11), & ! (in)
            ATMOS_PHY_MP_SN14_tracer_units(7:11)         ) ! (in)
      QA_MP = ATMOS_PHY_MP_SN14_ntracers
    case default
      LOG_ERROR("ATMOS_PHY_MP_setup",*) 'Not appropriate names of MP_TYPE in namelist PARAM_ATMOS_PHY_MP. Check!'
      call PRC_abort
    end select

    QE_MP = QS_MP + QA_MP - 1

    !- initialize the variables 
    call this%vars%Init( model_mesh, QS_MP, QE_MP, QA_MP )     

    !- Setup a module for microcloud physics

    select case( this%MP_TYPEID )
    case( MP_TYPEID_KESSLER )
      call ATMOS_PHY_MP_kessler_setup()
    case( MP_TYPEID_TOMITA08 )
      lcmesh3D => atm_mesh%ptr_mesh%lcmesh_list(1)
      elem3D => lcmesh3D%refElem3D
      call ATMOS_PHY_MP_tomita08_setup( &
        elem3D%Nnode_v*lcmesh3D%NeZ, 1, elem3D%Nnode_v*lcmesh3D%NeZ, elem3D%Nnode_h1D**2, 1, elem3D%Nnode_h1D**2, lcmesh3D%Ne2D, 1, lcmesh3D%Ne2D, &
        Flag_LT )    
    case( MP_TYPEID_SN14 )
    end select

    !
    call this%elem%Init( atm_mesh%ptr_mesh%refElem3D%PolyOrder_h, atm_mesh%ptr_mesh%refElem3D%PolyOrder_v, .true. )
    call this%Dz%Init( this%elem%Dx3, storage_format='ELL' )
    call this%Lift%Init( this%elem%Lift, storage_format='ELL' )

    return
  end subroutine AtmosPhyMp_setup

!> Calculate tendencies associated with cloud microphysics
!!
!!
!! @param model_mesh a object to manage computational mesh of atmospheric model 
!! @param prgvars_list a object to mange prognostic variables with atmospheric dynamical core
!! @param trcvars_list a object to mange auxiliary variables 
!! @param forcing_list a object to mange forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL
  subroutine AtmosPhyMp_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,     &
      AtmosVars_GetLocalMeshPhyAuxVars,  &
      AtmosVars_GetLocalMeshQTRCVarList, & 
      AtmosVars_GetLocalMeshPhyTends
    use mod_atmos_phy_mp_vars, only:           &
      AtmosPhyMPVars_GetLocalMeshFields_tend,  &
      AtmosPhyMpVars_GetLocalMeshFields_sfcflx
    
    implicit none
    
    class(AtmosPhyMp), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    class(MeshBase), pointer :: mesh
    class(LocalMesh3D), pointer :: lcmesh

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    type(LocalMeshFieldBaseList) :: QTRC(this%vars%QS:this%vars%QE)
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p
    type(LocalMeshFieldBaseList) :: RHOQ_tp(QA)
    class(LocalMeshFieldBase), pointer :: mp_DENS_t, mp_MOMX_t, mp_MOMY_t, mp_MOMZ_t, mp_RHOT_t, mp_RHOH, mp_EVAP
    type(LocalMeshFieldBaseList) :: mp_RHOQ_t(this%vars%QS:this%vars%QE)
    class(LocalMeshFieldBase), pointer :: SFLX_rain, SFLX_snow, SFLX_engi
    
    integer :: n
    integer :: ke
    integer :: iq
    !--------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_PROGRESS(*) 'atmosphere / physics / cloud microphysics'

    call model_mesh%GetModelMesh( mesh )

    do n=1, mesh%LOCAL_MESH_NUM
      call PROF_rapstart('ATM_PHY_MP_get_localmesh_ptr', 2)         
      call AtmosVars_GetLocalMeshPrgVars( n,  &
        mesh, prgvars_list, auxvars_list,       &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,         &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
        lcmesh                                  )
      
      call AtmosVars_GetLocalMeshQTRCVarList( n,  &
        mesh, trcvars_list, this%vars%QS, QTRC(:) )
      
      call AtmosVars_GetLocalMeshPhyAuxVars( n,  &
        mesh, auxvars_list,                      &
        PRES, PT )
      
      call AtmosVars_GetLocalMeshPhyTends( n,        &
        mesh, forcing_list,                          &
        DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, &
        RHOH_p, RHOQ_tp  )

      call AtmosPhyMpVars_GetLocalMeshFields_tend( n, &
        mesh, this%vars%tends_manager,                                           &
        mp_DENS_t, mp_MOMX_t, mp_MOMY_t, mp_MOMZ_t, mp_RHOT_t, mp_RHOH, mp_EVAP, &
        mp_RHOQ_t                                                                )

      call AtmosPhyMpVars_GetLocalMeshFields_sfcflx( n, &
        mesh, this%vars%auxvars2D_manager,              &
        SFLX_rain, SFLX_snow, SFLX_engi                 )     

      call PROF_rapend('ATM_PHY_MP_get_localmesh_ptr', 2)   

      call PROF_rapstart('ATM_PHY_MP_cal_tend', 2)
      if (is_update) then
        call this%calc_tendency_core( &
          mp_DENS_t%val, mp_MOMX_t%val, mp_MOMY_t%val, mp_MOMZ_t%val, mp_RHOQ_t,  & ! (out)
          mp_RHOH%val, mp_EVAP%val, SFLX_rain%val, SFLX_snow%val, SFLX_ENGI%val,  & ! (out)
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, PT%val, QTRC,                  & ! (in)
          PRES%val, PRES_hyd%val, DENS_hyd%val, Rtot%val, CVtot%val, CPtot%val,   & ! (in)
          model_mesh%DOptrMat(3), model_mesh%LiftOptrMat,                         & ! (in)
          lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D    ) ! (in)
      end if
      
      !$omp parallel private(ke, iq)
      !$omp do
      do ke = lcmesh%NeS, lcmesh%NeE
        DENS_tp%val(:,ke) = DENS_tp%val(:,ke) + mp_DENS_t%val(:,ke)
        MOMX_tp%val(:,ke) = MOMX_tp%val(:,ke) + mp_MOMX_t%val(:,ke)
        MOMY_tp%val(:,ke) = MOMY_tp%val(:,ke) + mp_MOMY_t%val(:,ke)
        MOMZ_tp%val(:,ke) = MOMZ_tp%val(:,ke) + mp_MOMZ_t%val(:,ke)
        RHOH_p %val(:,ke) = RHOH_p %val(:,ke) + mp_RHOH  %val(:,ke)
      end do
      !$omp end do
      !$omp do collapse(2)
      do iq = this%vars%QS, this%vars%QE
      do ke = lcmesh%NeS, lcmesh%NeE
        RHOQ_tp(iq)%ptr%val(:,ke) = RHOQ_tp(iq)%ptr%val(:,ke)  &
                                  + mp_RHOQ_t(iq)%ptr%val(:,ke)
      end do
      end do 
      !$omp end do
      !$omp end parallel

      call PROF_rapend('ATM_PHY_MP_cal_tend', 2)
    end do

    return  
  end subroutine AtmosPhyMp_calc_tendency

!> Update variables in a component of cloud microphysics
!!
!! @param model_mesh a object to manage computational mesh of atmospheric model 
!! @param prgvars_list a object to mange prognostic variables with atmospheric dynamical core
!! @param trcvars_list a object to mange auxiliary variables 
!! @param forcing_list a object to mange forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL  
  subroutine AtmosPhyMp_update( this, model_mesh, &
    prgvars_list, trcvars_list,                   &
    auxvars_list, forcing_list, is_update         )  
    
    implicit none
    class(AtmosPhyMp), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------

    return
  end subroutine AtmosPhyMp_update

!> Finalize a component of cloud microphysics
!!
!OCL SERIAL  
  subroutine AtmosPhyMp_finalize( this )
    implicit none
    class(AtmosPhyMp), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    select case ( this%MP_TYPEID )
    case( MP_TYPEID_KESSLER )
    end select

    call this%vars%Final()

    return
  end subroutine AtmosPhyMp_finalize
  
!- private ------------------------------------------------

!OCL SERIAL
  subroutine AtmosPhyMp_calc_tendency_core( this, &
    DENS_t_MP, RHOU_t_MP, RHOV_t_MP, MOMZ_t_MP, RHOQ_t_MP, & ! (out)
    RHOH_MP, EVAPORATE, SFLX_rain, SFLX_snow, SFLX_ENGI,   & ! (out)
    DDENS, RHOU, RHOV, MOMZ, PT, QTRC,                     & ! (in)
    PRES, PRES_hyd, DENS_hyd, Rtot, CVtot, CPtot,          & ! (in)
    Dz, Lift,                                              & ! (in)
    lcmesh, elem3D, lcmesh2D, elem2D                       ) ! (in)

    use scale_const, only: &
      PRE00 => CONST_PRE00
    use scale_const, only: &
      CVdry => CONST_CVdry, &
      CPdry => CONST_CPdry
    use scale_atmos_hydrometeor, only: & 
      LHF,                             &
      QHA, QLA, QIA
    use scale_atmos_phy_mp_kessler, only:    &
      ATMOS_PHY_MP_KESSLER_terminal_velocity
    use scale_atmos_phy_mp_tomita08, only: &
      ATMOS_PHY_MP_TOMITA08_terminal_velocity
    use scale_atmos_phy_mp_sn14, only: &
      ATMOS_PHY_MP_SN14_tendency, &
      ATMOS_PHY_MP_SN14_terminal_velocity
    use scale_sparsemat, only: SparseMat
    use scale_atm_phy_mp_dgm_common, only:  &
      atm_phy_mp_dgm_common_gen_intweight,         &
      atm_phy_mp_dgm_common_precipitation,         &
      atm_phy_mp_dgm_common_precipitation_momentum

    implicit none

    class(AtmosPhyMp), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: DENS_t_MP(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: RHOU_t_MP(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: RHOV_t_MP(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMZ_t_MP(elem3D%Np,lcmesh%NeA)
    type(LocalMeshFieldBaseList), intent(inout) :: RHOQ_t_MP(this%vars%QS:this%vars%QE)
    real(RP), intent(out) :: RHOH_MP(elem3D%Np,lcmesh%NeA)
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
    class(SparseMat), intent(in) :: Dz
    class(SparseMat), intent(in) :: Lift

    real(RP) :: DENS (elem3D%Np,lcmesh%NeA)

    real(RP) :: RHOE_t(elem3D%Np,lcmesh%NeA)

    real(RP) :: CPtot_t(elem3D%Np,lcmesh%NeA)
    real(RP) :: CVtot_t(elem3D%Np,lcmesh%NeA)

    real(RP) :: vterm(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D,this%vars%QS+1:this%vars%QE)
    real(RP) :: DENS0(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)    
    real(RP) :: DENS2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: REF_DENS(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOU2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOV2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: MOMZ2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: TEMP2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: PRES2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: CPtot2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: CVtot2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOE (elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOE2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOQ (elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D,this%vars%QS+1:this%vars%QE)
    real(RP) :: RHOQ2(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D,this%vars%QS+1:this%vars%QE)
    real(RP) :: RHOQ2_tmp(elem3D%Nnode_v,lcmesh%NeZ,this%vars%QS+1:this%vars%QE)
    real(RP) :: DENS2_tmp(elem3D%Nnode_v,lcmesh%NeZ)
    real(RP) :: TEMP2_tmp(elem3D%Nnode_v,lcmesh%NeZ)
    real(RP) :: PRES2_tmp(elem3D%Nnode_v,lcmesh%NeZ)
    real(RP) :: REF_DENS_tmp(elem3D%Nnode_v,lcmesh%NeZ)
    real(RP) :: vterm_tmp(elem3D%Nnode_v,lcmesh%NeZ,this%vars%QS+1:this%vars%QE)
    real(RP) :: FLX_hydro(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: CP_t(elem3D%Np), CV_t(elem3D%Np)

    integer :: iq
    integer :: domid
    integer :: ke
    integer :: ke2D, ke_z
    integer :: p2D, pv1D, p
    integer :: ColMask(elem3D%Nnode_v)

    integer :: step
    real(RP) :: rdt_MP

    integer :: vmapM(elem3D%NfpTot,lcmesh%NeZ)
    integer :: vmapP(elem3D%NfpTot,lcmesh%NeZ)
    real(RP) :: IntWeight(elem3D%Nfaces,elem3D%NfpTot)
    real(RP) :: nz(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D)
    !----------------------------------------------

    rdt_MP = 1.0_RP / this%dtsec
    domid  = lcmesh%lcdomID

    call lcmesh%GetVmapZ1D( vmapM, vmapP ) ! (out)
    
    call atm_phy_mp_dgm_common_gen_intweight( IntWeight, & ! (out)
      lcmesh                                             ) ! (in)
    
    !$omp parallel do private(ke)
    do ke = lcmesh%NeS, lcmesh%NeE
      DENS(:,ke) = DENS_hyd(:,ke) + DDENS(:,ke)
    end do

    select case( this%MP_TYPEID )
    case( MP_TYPEID_KESSLER )
      call calc_tendency_Kessler( this, &
        RHOQ_t_MP, CPtot_t, CVtot_t, RHOE_t, EVAPORATE,  & ! (out)
        DENS, QTRC, PRES, DENS_hyd, Rtot, CVtot, CPtot,  & ! (in)
        rdt_MP, lcmesh, elem3D )                           ! (in)

    case( MP_TYPEID_TOMITA08 )
      call calc_tendency_Tomita08( this, &
        RHOQ_t_MP, CPtot_t, CVtot_t, RHOE_t, EVAPORATE,  & ! (out)
        DENS, QTRC, PRES, DENS_hyd, Rtot, CVtot, CPtot,  & ! (in)
        rdt_MP, lcmesh, elem3D )                           ! (in)

    case( MP_TYPEID_SN14 )
      call calc_tendency_SN14( this, &
        RHOQ_t_MP, CPtot_t, CVtot_t, RHOE_t, EVAPORATE,  & ! (out)
        DENS, QTRC, PRES, DENS_hyd, Rtot, CVtot, CPtot,  & ! (in)
        rdt_MP, lcmesh, elem3D )                           ! (in)

    end select

    !$omp parallel do
    do ke = lcmesh%NeS, lcmesh%NeE
      RHOH_MP(:,ke) = RHOE_t(:,ke) &
        - ( CPtot_t(:,ke) + log( PRES(:,ke) / PRE00 ) * ( CVtot(:,ke) / CPtot(:,ke) * CPtot_t(:,ke) - CVtot_t(:,ke) ) ) &
        * PRES(:,ke) / Rtot(:,ke)
    end do

    if ( this%do_precipitation ) then

      !$omp parallel private(ke,ke2D,ke_z,iq)
      !$omp do collapse(2)
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D

        DENS0(:,ke_z,ke2D) = DENS(:,ke)        
        DENS2(:,ke_z,ke2D) = DENS(:,ke)
        REF_DENS(:,ke_z,ke2D) = DENS_hyd(:,ke)
        RHOU2(:,ke_z,ke2D) = RHOU(:,ke)
        RHOV2(:,ke_z,ke2D) = RHOV(:,ke)
        MOMZ2(:,ke_z,ke2D) = MOMZ(:,ke)
        TEMP2(:,ke_z,ke2D) = PRES(:,ke) / ( DENS(:,ke) * Rtot(:,ke) )
        PRES2(:,ke_z,ke2D) = PRES(:,ke)
        CPtot2(:,ke_z,ke2D) = CPtot(:,ke)
        CVtot2(:,ke_z,ke2D) = CVtot(:,ke)
        RHOE (:,ke_z,ke2D) = TEMP2(:,ke_z,ke2D) * CVtot(:,ke) * DENS(:,ke)
        RHOE2(:,ke_z,ke2D) = RHOE(:,ke_z,ke2D)

        nz(:,ke_z,ke2D) = lcmesh%normal_fn(:,ke,3)

        FLX_hydro(:,ke_z,ke2D) = 0.0_RP
      end do
      end do
      !$omp end do
      !$omp do collapse(3)
      do iq = this%vars%QS+1, this%vars%QE
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D
        RHOQ (:,ke_z,ke2D,iq) = DENS2(:,ke_z,ke2D) * QTRC(iq)%ptr%val(:,ke) &
                              + RHOQ_t_MP(iq)%ptr%val(:,ke) * this%dtsec
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

      do iq = this%vars%QS+1, this%vars%QE
        if ( this%vars%vterm_hist_id(iq) > 0 ) then
          !$omp parallel do
          do ke = lcmesh%NeS, lcmesh%NeE
            this%vars%vterm_hist(iq)%local(domid)%val(:,ke) = 0.0_RP
          end do
        end if
      end do

      do step = 1, this%nstep_sedmientation

        !- Calculate terminal velocity 

          do ke2D = 1, lcmesh%Ne2D
            !$omp parallel do private( &
            !$omp p2D, ke_z, iq, ColMask,                                   &
            !$omp DENS2_tmp, REF_DENS_tmp, RHOQ2_tmp, TEMP2_tmp, PRES2_tmp, &
            !$omp vterm_tmp ) 
            do p2D = 1, elem3D%Nnode_h1D**2
              ColMask(:) = elem3D%Colmask(:,p2D)
              do ke_z = 1, lcmesh%NeZ
                DENS2_tmp   (:,ke_z) = DENS2   (ColMask(:),ke_z,ke2D)
                REF_DENS_tmp(:,ke_z) = REF_DENS(ColMask(:),ke_z,ke2D)
                TEMP2_tmp   (:,ke_z) = TEMP2   (ColMask(:),ke_z,ke2D)
                PRES2_tmp   (:,ke_z) = PRES2   (ColMask(:),ke_z,ke2D)
              end do
              do iq = this%vars%QS+1, this%vars%QE
              do ke_z = 1, lcmesh%NeZ            
                RHOQ2_tmp(:,ke_z,iq) = RHOQ2(ColMask(:),ke_z,ke2D,iq)
              end do
              end do

              select case( this%MP_TYPEID )
              case( MP_TYPEID_KESSLER )
                call ATMOS_PHY_MP_KESSLER_terminal_velocity( &
                  elem3D%Nnode_v*lcmesh%NeZ, 1, elem3D%Nnode_v*lcmesh%NeZ, & ! (in)
                  DENS2_tmp(:,:), RHOQ2_tmp(:,:,:), REF_DENS_tmp(:,:),     & ! (in)
                  vterm_tmp(:,:,:)                                         ) ! (out)
              case( MP_TYPEID_TOMITA08 )
                call ATMOS_PHY_MP_TOMITA08_terminal_velocity( &
                  elem3D%Nnode_v*lcmesh%NeZ, 1, elem3D%Nnode_v*lcmesh%NeZ, & ! (in)
                  DENS2_tmp(:,:), TEMP2_tmp(:,:), RHOQ2_tmp(:,:,:),        & ! (in)
                  vterm_tmp(:,:,:)                                         ) ! (out)
              case( MP_TYPEID_SN14 )
                call ATMOS_PHY_MP_SN14_terminal_velocity( &
                  elem3D%Nnode_v*lcmesh%NeZ, 1, elem3D%Nnode_v*lcmesh%NeZ,          & ! (in)
                  DENS2_tmp(:,:), TEMP2_tmp(:,:), RHOQ2_tmp(:,:,:), PRES2_tmp(:,:), & ! (in)
                  vterm_tmp(:,:,:)                                                  ) ! (out)
              end select   

              do iq = this%vars%QS+1, this%vars%QE
              do ke_z = 1, lcmesh%NeZ            
                vterm(ColMask(:),ke_z,ke2D,iq) = vterm_tmp(:,ke_z,iq)
              end do
              end do
            end do
          end do

        
        do iq = this%vars%QS+1, this%vars%QE
          if ( this%vars%vterm_hist_id(iq) > 0 ) then
            !$omp parallel do collapse(2) private(ke2D, ke_z, ke)
            do ke2D = 1, lcmesh%Ne2D
            do ke_z = 1, lcmesh%NeZ
              ke = ke2D + (ke_z-1)*lcmesh%Ne2D
              this%vars%vterm_hist(iq)%local(domid)%val(:,ke) = &
                this%vars%vterm_hist(iq)%local(domid)%val(:,ke) + vterm(:,ke_z,ke2D,iq) * this%rnstep_sedmientation
            end do
            end do
          end if
        end do 

        !- precipiation of hydrometers

        call atm_phy_mp_dgm_common_precipitation( &
          DENS2, RHOQ2, CPtot2, CVtot2, RHOE2,                 & ! (inout)
          FLX_hydro, SFLX_rain, SFLX_snow, SFLX_ENGI,          & ! (inout)
          TEMP2, vterm,                                        & ! (in)
          this%dtsec_sedmientation, this%rnstep_sedmientation, & ! (in)
!          Dz, Lift, nz, vmapM, vmapP, IntWeight,               & ! (in)
          this%Dz, this%Lift, nz, vmapM, vmapP, IntWeight,               & ! (in)
          this%vars%QE - this%vars%QS, QLA, QIA,               & ! (in)
          lcmesh, elem3D )

        !$omp parallel do private(ke2D, ke_z, ke) collapse(2)
        do ke2D = 1, lcmesh%Ne2D
        do ke_z = 1, lcmesh%NeZ
          ke = ke2D + (ke_z-1)*lcmesh%Ne2D
          TEMP2(:,ke_z,ke2D) = RHOE2(:,ke_z,ke2D) / ( DENS2(:,ke_z,ke2D) * CVtot2(:,ke_z,ke2D) )
        end do
        end do 
        
      end do ! end loop for step

      !$omp parallel private(ke2D, ke_z, iq, ke, CP_t, CV_t)
      !$omp workshare
      SFLX_ENGI(:,:) = SFLX_ENGI(:,:) - SFLX_snow(:,:) * LHF ! moist internal energy
      !$omp end workshare
      !$omp do collapse(2)
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D
        DENS_t_MP(:,ke) = ( DENS2(:,ke_z,ke2D) - DENS(:,ke) ) * rdt_MP

        CP_t(:) = ( CPtot2(:,ke_z,ke2D) - CPtot(:,ke) ) * rdt_MP
        CV_t(:) = ( CVtot2(:,ke_z,ke2D) - CVtot(:,ke) ) * rdt_MP
        RHOH_MP(:,ke) = RHOH_MP(:,ke) &
          + ( RHOE2(:,ke_z,ke2D) - RHOE(:,ke_z,ke2D) ) * rdt_MP &
          - ( CP_t(:) &
              + log( PRES(:,ke) / PRE00 ) * ( CVtot(:,ke) / CPtot(:,ke) * CP_t(:) - CV_t(:) ) &
            ) * PRES(:,ke) / Rtot(:,ke) 
      end do
      end do
      !$omp end do
      !$omp do collapse(2)
      do iq = this%vars%QS+1, this%vars%QE
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D
        RHOQ_t_MP(iq)%ptr%val(:,ke) = RHOQ_t_MP(iq)%ptr%val(:,ke)       &
          + ( RHOQ2(:,ke_z,ke2D,iq) - RHOQ(:,ke_z,ke2D,iq) ) * rdt_MP
      end do
      end do
      end do
      !$omp end do
      !$omp end parallel

      !- precipiation of momentum

      call atm_phy_mp_dgm_common_precipitation_momentum( &
        RHOU_t_MP, RHOV_t_MP, MOMZ_t_MP,          & ! (out)
        DENS0, RHOU2, RHOV2, MOMZ2, FLX_hydro,    & ! (in)
!        Dz, Lift, nz, vmapM, vmapP,              & ! (in)
        this%Dz, this%Lift, nz, vmapM, vmapP,     & ! (in)
        lcmesh, elem3D                            ) ! (in)

    end if ! end if do_precipitation

    return
  end subroutine AtmosPhyMp_calc_tendency_core

!OCL SERIAL
  subroutine calc_tendency_Kessler( this, &
    RHOQ_t_MP, CPtot_t, CVtot_t, RHOE_t, EVAPORATE,  & ! (out)
    DENS, QTRC, PRES, DENS_hyd, Rtot, CVtot, CPtot,  & ! (in)
    rdt_MP, lcmesh, elem3D )                           ! (in)

    use scale_atmos_phy_mp_kessler, only:    &
      ATMOS_PHY_MP_KESSLER_adjustment
    implicit none

    class(AtmosPhyMp), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    type(LocalMeshFieldBaseList), intent(inout) :: RHOQ_t_MP(this%vars%QS:this%vars%QE)
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
    real(RP), intent(in) :: rdt_MP

    real(RP) :: TEMP1(elem3D%Np,lcmesh%NeA)
    real(RP) :: CPtot1(elem3D%Np,lcmesh%NeA)
    real(RP) :: CVtot1(elem3D%Np,lcmesh%NeA)
    real(RP) :: QTRC1(elem3D%Np,lcmesh%NeA,this%vars%QS:this%vars%QE)

    integer :: ke
    integer :: iq
    !------------------------------------------------------------

    !$omp parallel private(ke, iq)
    !$omp do
    do ke = lcmesh%NeS, lcmesh%NeE
      TEMP1(:,ke) = PRES(:,ke) / ( DENS(:,ke) * Rtot(:,ke) )
      CPtot1(:,ke) = CPtot(:,ke)
      CVtot1(:,ke) = CVtot(:,ke)
    end do
    !$omp end do
    !$omp do collapse(2)
    do iq = this%vars%QS, this%vars%QE
    do ke = lcmesh%NeS, lcmesh%NeE
      QTRC1(:,ke,iq) = QTRC(iq)%ptr%val(:,ke)
    end do
    end do
    !$omp end do
    !$omp end parallel

    call ATMOS_PHY_MP_kessler_adjustment( &
      elem3D%Np, 1, elem3D%Np, lcmesh%NeA, lcmesh%NeS, lcmesh%NeE, 1, 1, 1,  & ! (in)
      DENS, PRES, this%dtsec,                                                & ! (in)
      TEMP1, QTRC1, CPtot1, CVtot1,                                          & ! (inout)
      RHOE_t, EVAPORATE                                                      ) ! (out)
  
    !$omp parallel private(ke, iq)
    !$omp do collapse(2)
    do iq = this%vars%QS, this%vars%QE
    do ke = lcmesh%NeS, lcmesh%NeE
      RHOQ_t_MP(iq)%ptr%val(:,ke) = ( QTRC1(:,ke,iq) - QTRC(iq)%ptr%val(:,ke) ) * DENS(:,ke) * rdt_MP
    end do  
    end do
    !$omp end do
    !$omp do
    do ke = lcmesh%NeS, lcmesh%NeE
      CPtot_t(:,ke) = ( CPtot1(:,ke) - CPtot(:,ke) ) * rdt_MP
      CVtot_t(:,ke) = ( CVtot1(:,ke) - CVtot(:,ke) ) * rdt_MP
    end do
    !$omp end do
    !$omp end parallel
    return
  end subroutine calc_tendency_Kessler

!OCL SERIAL
  subroutine calc_tendency_Tomita08( this, &
    RHOQ_t_MP, CPtot_t, CVtot_t, RHOE_t, EVAPORATE,  & ! (out)
    DENS, QTRC, PRES, DENS_hyd, Rtot, CVtot, CPtot,  & ! (in)
    rdt_MP, lcmesh, elem3D )                           ! (in)

    use scale_const, only: &
      UNDEF => CONST_UNDEF
    use scale_atmos_phy_mp_tomita08, only: &
      ATMOS_PHY_MP_TOMITA08_adjustment
    implicit none

    class(AtmosPhyMp), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    type(LocalMeshFieldBaseList), intent(inout) :: RHOQ_t_MP(this%vars%QS:this%vars%QE)
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
    real(RP), intent(in) :: rdt_MP

    real(RP) :: DENS_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: PRES_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: CCN_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: TEMP1_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: CPtot1_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: CVtot1_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: QTRC1_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D,this%vars%QS:this%vars%QE)
    real(RP) :: RHOE_t_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: EVAPORATE_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)

    integer :: ke, ke_xy, ke_z
    integer :: p,p_xy, p_z
    integer :: iq
    !------------------------------------------------------------

    !$omp parallel do collapse(2) private(p_z,p_xy,iq, ke,p)
    do ke_z=1, lcmesh%NeZ
    do ke_xy=1, lcmesh%Ne2D
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D
      do p_z=1, elem3D%Nnode_v
      do p_xy=1, elem3D%Nnode_h1D**2
        p = p_xy + (p_z-1)*elem3D%Nnode_h1D**2
        DENS_z(p_z,ke_z,p_xy,ke_xy) = DENS(p,ke)
        PRES_z(p_z,ke_z,p_xy,ke_xy) = PRES(p,ke)
        TEMP1_z(p_z,ke_z,p_xy,ke_xy) = PRES(p,ke) / ( DENS(p,ke) * Rtot(p,ke) )
        CPtot1_z(p_z,ke_z,p_xy,ke_xy) = CPtot(p,ke)
        CVtot1_z(p_z,ke_z,p_xy,ke_xy) = CVtot(p,ke)

        CCN_z(p_z,ke_z,p_xy,ke_xy) = UNDEF
      end do
      end do
      do iq = this%vars%QS, this%vars%QE
        do p_z=1, elem3D%Nnode_v
        do p_xy=1, elem3D%Nnode_h1D**2
          p = p_xy + (p_z-1)*elem3D%Nnode_h1D**2
          QTRC1_z(p_z,ke_z,p_xy,ke_xy,iq) = QTRC(iq)%ptr%val(p,ke)
        end do
        end do
      end do
    end do
    end do

    call ATMOS_PHY_MP_Tomita08_adjustment( &
      elem3D%Nnode_v*lcmesh%NeZ, 1, elem3D%Nnode_v*lcmesh%NeZ, elem3D%Nnode_h1D**2, 1, elem3D%Nnode_h1D**2, lcmesh%Ne2D, 1, lcmesh%Ne2D, & ! (in)
      DENS_z, PRES_z, CCN_z, this%dtsec,       & ! (in)
      TEMP1_z, QTRC1_z, CPtot1_z, CVtot1_z,    & ! (inout)
      RHOE_t_z, EVAPORATE_z                    ) ! (out)
  
    !$omp parallel do collapse(2) private(p_z,p_xy,iq, ke,p)
    do ke_z=1, lcmesh%NeZ
    do ke_xy=1, lcmesh%Ne2D
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D
      do p_z=1, elem3D%Nnode_v
      do p_xy=1, elem3D%Nnode_h1D**2
        p = p_xy + (p_z-1)*elem3D%Nnode_h1D**2
        CPtot_t(p,ke) = ( CPtot1_z(p_z,ke_z,p_xy,ke_xy) - CPtot(p,ke) ) * rdt_MP
        CVtot_t(p,ke) = ( CVtot1_z(p_z,ke_z,p_xy,ke_xy) - CVtot(p,ke) ) * rdt_MP
        RHOE_t(p,ke) = RHOE_t_z(p_z,ke_z,p_xy,ke_xy)
        EVAPORATE(p,ke) = EVAPORATE_z(p_z,ke_z,p_xy,ke_xy)
      end do
      end do
      do iq = this%vars%QS, this%vars%QE
        do p_z=1, elem3D%Nnode_v
        do p_xy=1, elem3D%Nnode_h1D**2
          p = p_xy + (p_z-1)*elem3D%Nnode_h1D**2
          RHOQ_t_MP(iq)%ptr%val(p,ke) = ( QTRC1_z(p_z,ke_z,p_xy,ke_xy,iq) - QTRC(iq)%ptr%val(p,ke) ) * DENS(p,ke) * rdt_MP
        end do
        end do
      end do
    end do
    end do

    return
  end subroutine calc_tendency_Tomita08

!OCL SERIAL
  subroutine calc_tendency_SN14( this, &
    RHOQ_t_MP, CPtot_t, CVtot_t, RHOE_t, EVAPORATE,  & ! (out)
    DENS, QTRC, PRES, DENS_hyd, Rtot, CVtot, CPtot,  & ! (in)
    rdt_MP, lcmesh, elem3D )                           ! (in)

    use scale_const, only: &
      UNDEF => CONST_UNDEF
    use scale_atmos_phy_mp_sn14, only: &
       ATMOS_PHY_MP_sn14_tendency
    implicit none

    class(AtmosPhyMp), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    type(LocalMeshFieldBaseList), intent(inout) :: RHOQ_t_MP(this%vars%QS:this%vars%QE)
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
    real(RP), intent(in) :: rdt_MP

    real(RP) :: DENS_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: PRES_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: CCN_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: TEMP1_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: CPtot1_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: CVtot1_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: QTRC1_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D,this%vars%QS:this%vars%QE)
    real(RP) :: RHOQ_t_MP_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D,this%vars%QS:this%vars%QE)
    real(RP) :: RHOE_t_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: EVAPORATE_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: CPtot_t_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: CVtot_t_z(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)

    integer :: ke, ke_xy, ke_z
    integer :: p,p_xy, p_z
    integer :: iq
    !------------------------------------------------------------

    !$omp parallel do collapse(2) private(p_z,p_xy,iq, ke,p)
    do ke_z=1, lcmesh%NeZ
    do ke_xy=1, lcmesh%Ne2D
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D
      do p_z=1, elem3D%Nnode_v
      do p_xy=1, elem3D%Nnode_h1D**2
        p = p_xy + (p_z-1)*elem3D%Nnode_h1D**2
        DENS_z(p_z,ke_z,p_xy,ke_xy) = DENS(p,ke)
        PRES_z(p_z,ke_z,p_xy,ke_xy) = PRES(p,ke)
        TEMP1_z(p_z,ke_z,p_xy,ke_xy) = PRES(p,ke) / ( DENS(p,ke) * Rtot(p,ke) )
        CPtot1_z(p_z,ke_z,p_xy,ke_xy) = CPtot(p,ke)
        CVtot1_z(p_z,ke_z,p_xy,ke_xy) = CVtot(p,ke)

        CCN_z(p_z,ke_z,p_xy,ke_xy) = UNDEF
      end do
      end do
      do iq = this%vars%QS, this%vars%QE
        do p_z=1, elem3D%Nnode_v
        do p_xy=1, elem3D%Nnode_h1D**2
          p = p_xy + (p_z-1)*elem3D%Nnode_h1D**2
          QTRC1_z(p_z,ke_z,p_xy,ke_xy,iq) = QTRC(iq)%ptr%val(p,ke)
        end do
        end do
      end do
    end do
    end do

    ! call  ATMOS_PHY_MP_sn14_tendency( &
    !   elem3D%Nnode_v*lcmesh%NeZ, 1, elem3D%Nnode_v*lcmesh%NeZ, elem3D%Nnode_h1D**2, 1, elem3D%Nnode_h1D**2, lcmesh%Ne2D, 1, lcmesh%Ne2D, & ! (in)
    !   DENS_z, W_z, QTRC1_z, PRES_z, TEMP_z, Qdry_z, CPtot1_z, CVtot1_z, CCN_z, this%dtsec, REAL_CZ, REAL_FZ,  & ! (in)
    !   RHOQ_t_MP_z, RHOE_t_z, CPtot_t_z, CVtot_t_z, EVAPORATE_z )                                                ! (out)
  

    !$omp parallel do collapse(2) private(p_z,p_xy,iq, ke,p)
    do ke_z=1, lcmesh%NeZ
    do ke_xy=1, lcmesh%Ne2D
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D
      do p_z=1, elem3D%Nnode_v
      do p_xy=1, elem3D%Nnode_h1D**2
        p = p_xy + (p_z-1)*elem3D%Nnode_h1D**2
        CPtot_t(p,ke) = CPtot_t_z(p_z,ke_z,p_xy,ke_xy)
        CVtot_t(p,ke) = CVtot_t_z(p_z,ke_z,p_xy,ke_xy)
      end do
      end do
      do iq = this%vars%QS, this%vars%QE
        do p_z=1, elem3D%Nnode_v
        do p_xy=1, elem3D%Nnode_h1D**2
          p = p_xy + (p_z-1)*elem3D%Nnode_h1D**2
          RHOQ_t_MP(iq)%ptr%val(p,ke) = RHOQ_t_MP_z(p_z,ke_z,p_xy,ke_xy,iq)
        end do
        end do
      end do
    end do
    end do

    return
  end subroutine calc_tendency_SN14  
end module mod_atmos_phy_mp