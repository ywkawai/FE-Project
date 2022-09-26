!-------------------------------------------------------------------------------
!> module ATMOSPHERE dynamical process
!!
!! @par Description
!!          Module for atmosphere dynamical process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_dyn
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

  use scale_sparsemat, only: SparseMat
  use scale_timeint_rk, only: TimeInt_RK
  
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
    
  use scale_meshfield_base, only: MeshFieldBase
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_component_proc, only:  ModelComponentProc

  use scale_atm_dyn_dgm_nonhydro3d_heve, only: &
    atm_dyn_dgm_nonhydro3d_heve_Init,          &
    atm_dyn_dgm_nonhydro3d_heve_Final,         &
    atm_dyn_dgm_nonhydro3d_heve_cal_tend
    
  use scale_atm_dyn_dgm_nonhydro3d_hevi, only: &
    atm_dyn_dgm_nonhydro3d_hevi_Init,          &
    atm_dyn_dgm_nonhydro3d_hevi_Final,         &
    atm_dyn_dgm_nonhydro3d_hevi_cal_tend,      &
    atm_dyn_dgm_nonhydro3d_hevi_cal_vi

  use scale_atm_dyn_dgm_nonhydro3d_splitform_heve, only: &
    atm_dyn_dgm_nonhydro3d_heve_splitform_Init,          &
    atm_dyn_dgm_nonhydro3d_heve_splitform_Final,         &
    atm_dyn_dgm_nonhydro3d_heve_splitform_cal_tend

  use scale_atm_dyn_dgm_nonhydro3d_splitform_hevi, only: &
    atm_dyn_dgm_nonhydro3d_hevi_splitform_Init,          &
    atm_dyn_dgm_nonhydro3d_hevi_splitform_Final,         &
    atm_dyn_dgm_nonhydro3d_hevi_splitform_cal_tend,      &
    atm_dyn_dgm_nonhydro3d_hevi_splitform_cal_vi    
  
  use scale_atm_dyn_dgm_globalnonhydro3d_heve, only: &
    atm_dyn_dgm_globalnonhydro3d_heve_Init,          &
    atm_dyn_dgm_globalnonhydro3d_heve_Final,         &
    atm_dyn_dgm_globalnonhydro3d_heve_cal_tend

  use scale_atm_dyn_dgm_globalnonhydro3d_etot_heve, only: &
    atm_dyn_dgm_globalnonhydro3d_etot_heve_Init,          &
    atm_dyn_dgm_globalnonhydro3d_etot_heve_Final,         &
    atm_dyn_dgm_globalnonhydro3d_etot_heve_cal_tend

  use scale_atm_dyn_dgm_globalnonhydro3d_hevi, only:  &
    atm_dyn_dgm_globalnonhydro3d_hevi_Init,           &
    atm_dyn_dgm_globalnonhydro3d_hevi_Final,          &
    atm_dyn_dgm_globalnonhydro3d_hevi_cal_tend,       &
    atm_dyn_dgm_globalnonhydro3d_hevi_cal_vi

  use scale_atm_dyn_dgm_globalnonhydro3d_etot_hevi, only:  &
    atm_dyn_dgm_globalnonhydro3d_etot_hevi_Init,           &
    atm_dyn_dgm_globalnonhydro3d_etot_hevi_Final,          &
    atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_tend,       &
    atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_vi

  use scale_atm_dyn_dgm_nonhydro3d_numdiff, only: &
    atm_dyn_dgm_nonhydro3d_numdiff_Init,          &
    atm_dyn_dgm_nonhydro3d_numdiff_Final
  
  use scale_element_modalfilter, only: ModalFilter

  use scale_atm_dyn_dgm_trcadvect3d_heve, only: &
    atm_dyn_dgm_trcadvect3d_heve_Init,                 &
    atm_dyn_dgm_trcadvect3d_heve_Final,                &    
    atm_dyn_dgm_trcadvect3d_heve_calc_fct_coef,        &
    atm_dyn_dgm_trcadvect3d_heve_cal_tend,             &
    atm_dyn_dgm_trcadvect3d_TMAR,                      &
    atm_dyn_dgm_trcadvect3d_save_massflux,             &
    atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_advtest

  use mod_atmos_mesh, only: AtmosMesh
  use mod_atmos_vars, only: &
    AtmosVars_GetLocalMeshPrgVar,        &
    AtmosVars_GetLocalMeshPrgVars,       &
    AtmosVars_GetLocalMeshQTRCVar,       &
    AtmosVars_GetLocalMeshQTRCPhyTend,   &
    ATMOS_PROGVARS_NUM,                  &
    DDENS_ID => ATMOS_PROGVARS_DDENS_ID, &
    DRHOT_ID => ATMOS_PROGVARS_DRHOT_ID, &
    EnTot_ID => ATMOS_PROGVARS_EnTot_ID, &
    THERM_ID => ATMOS_PROGVARS_THERM_ID, &
    MOMX_ID  => ATMOS_PROGVARS_MOMX_ID,  &
    MOMY_ID  => ATMOS_PROGVARS_MOMY_ID,  &
    MOMZ_ID  => ATMOS_PROGVARS_MOMZ_ID
  use mod_atmos_dyn_bnd, only: AtmosDynBnd
  use mod_atmos_dyn_vars, only: &
    AtmosDynVars,                                &
    AtmosDynAuxVars_GetLocalMeshFields,          &
    AtmosDynMassFlux_GetLocalMeshFields,         &
    TRCQ_ID    => ATMOS_DYN_TRCVARS3D_TRCADV_ID, &
    TRCDDENS_ID => ATMOS_DYN_TRCVARS3D_DENS_ID

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  abstract interface    
    subroutine atm_dyn_nonhydro3d_cal_tend_ex( &
      DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, CORIOLIS,          & ! (in)
      Rtot, CVtot, CPtot,                                                         & ! (in)
      SL_flag, wdamp_tau, wdamp_height, hveldamp_flag,                            & ! (in)
      Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )

      import RP
      import LocalMesh3D
      import elementbase3D
      import LocalMesh2D
      import elementbase2D
      import SparseMat
      implicit none

      class(LocalMesh3D), intent(in) :: lmesh
      class(elementbase3D), intent(in) :: elem
      class(LocalMesh2D), intent(in) :: lmesh2D
      class(elementbase2D), intent(in) :: elem2D
      type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift
      real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMY_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: RHOT_dt(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
      real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
      logical, intent(in)   :: SL_flag
      real(RP), intent(in)  :: wdamp_tau
      real(RP), intent(in)  :: wdamp_height
      logical, intent(in) :: hveldamp_flag
    end subroutine atm_dyn_nonhydro3d_cal_tend_ex
  end interface

  abstract interface    
    subroutine atm_dyn_nonhydro3d_cal_vi( &
      DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,             & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, & ! (in)
      DDENS0_, MOMX0_, MOMY0_, MOMZ0_, DRHOT0_,                & ! (in) 
      Rtot, CVtot, CPtot,                                      & ! (in)
      Dz, Lift,                                                & ! (in)
      modalFilterFlag, VModalFilter,                           & ! (in)
      impl_fac, dt,                                            & ! (in)
      lmesh, elem, lmesh2D, elem2D )
  
      import RP
      import LocalMesh3D
      import ElementBase3D
      import LocalMesh2D
      import ElementBase2D
      import ModalFilter
      import SparseMat
      implicit none
  
      class(LocalMesh3D), intent(in) :: lmesh
      class(elementbase3D), intent(in) :: elem
      class(LocalMesh2D), intent(in) :: lmesh2D
      class(elementbase2D), intent(in) :: elem2D
      real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMY_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: RHOT_dt(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DDENS0_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMX0_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMY0_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMZ0_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DRHOT0_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
      class(SparseMat), intent(in) :: Dz, Lift
      logical, intent(in) :: modalFilterFlag
      class(ModalFilter), intent(in) :: VModalFilter
      real(RP), intent(in) :: impl_fac
      real(RP), intent(in) :: dt
    end subroutine atm_dyn_nonhydro3d_cal_vi
  end interface

  type, extends(ModelComponentProc), public :: AtmosDyn
    integer :: EQS_TYPEID
    type(TimeInt_RK), allocatable :: tint(:)
    type(TimeInt_RK), allocatable :: tint_qtrc(:)
    type(AtmosDynBnd) :: boundary_cond
    type(AtmosDynVars) :: dyn_vars

    procedure (atm_dyn_nonhydro3d_cal_vi), pointer, nopass :: cal_vi => null()
    procedure (atm_dyn_nonhydro3d_cal_tend_ex), pointer, nopass :: cal_tend_ex => null()

    ! explicit numerical diffusion
    logical :: CALC_NUMDIFF_FLAG
    integer  :: ND_LAPLACIAN_NUM
    real(RP) :: ND_COEF_H
    real(RP) :: ND_COEF_V

    ! element-wise modal filter
    logical :: MODALFILTER_FLAG
    type(ModalFilter) :: modal_filter_3d
    type(ModalFilter) :: modal_filter_tracer_3d
    type(ModalFilter) :: modal_filter_v1D

    ! sponge layer
    logical :: SPONGELAYER_FLAG
    real(RP) :: wdamp_tau
    real(RP) :: wdamp_height
    logical  :: hvel_damp_flag

    ! tracer advection
    logical :: ONLY_TRACERADV_FLAG
    logical :: TRACERADV_disable_limiter
    logical :: TRACERADV_MODALFILTER_FLAG
    type(SparseMat) :: FaceIntMat    

    ! 
    logical :: ENTOT_CONSERVE_SCHEME_FLAG

  contains
    procedure, public :: setup => AtmosDyn_setup 
    procedure, public :: calc_tendency => AtmosDyn_calc_tendency
    procedure, public :: update => AtmosDyn_update
    procedure, public :: finalize => AtmosDyn_finalize
  end type AtmosDyn

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------
  
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_HEVE             = 1
  integer, public, parameter :: EQS_TYPEID_GLOBALNONHYD3D_HEVE       = 2
  integer, public, parameter :: EQS_TYPEID_GLOBALNONHYD3D_HEVE_ENTOT = 3
  integer, public, parameter :: EQS_TYPEID_GLOBALNONHYD3D_HEVI       = 4
  integer, public, parameter :: EQS_TYPEID_GLOBALNONHYD3D_HEVI_ENTOT = 5
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_SPLITFORM_HEVE   = 6
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_HEVI             = 7
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_SPLITFORM_HEVI   = 8


  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: cal_numfilter_tend
  private :: add_phy_tend

  private :: setup_modalfilter
  private :: setup_numdiff
  private :: setup_spongelayer
  private :: setup_coriolis_parameter

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------

contains

!OCL SERIAL
  subroutine AtmosDyn_setup( this, model_mesh, tm_parent_comp )
    use mod_atmos_mesh, only: AtmosMesh
    use mod_atmos_vars, only: ATMOS_PROGVARS_NUM
    use scale_time_manager, only: TIME_manager_component

    implicit none

    class(AtmosDyn), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    character(len=H_MID) :: EQS_TYPE             = "NONHYDRO3D_HEVE"
    character(len=H_SHORT) :: TINTEG_TYPE        = 'ERK_SSP_3s3o'
    character(len=H_SHORT) :: TINTEG_TYPE_TRACER = 'ERK_SSP_3s3o'    
    real(DP) :: TIME_DT                          = UNDEF8
    character(len=H_SHORT) :: TIME_DT_UNIT       = 'SEC'  
    
    logical :: MODALFILTER_FLAG           = .false.
    logical :: NUMDIFF_FLAG               = .false.
    logical :: SPONGELAYER_FLAG           = .false.
    logical :: ONLY_TRACERADV_FLAG        = .false.
    logical :: TRACERADV_DISABLE_LIMITER  = .false.
    logical :: TRACERADV_MODALFILTER_FLAG = .false.

    namelist / PARAM_ATMOS_DYN /       &
      EQS_TYPE,                        &
      TINTEG_TYPE,                     &
      TINTEG_TYPE_TRACER,              &      
      TIME_DT,                         &
      TIME_DT_UNIT,                    &
      MODALFILTER_FLAG,                &
      NUMDIFF_FLAG,                    &
      SPONGELAYER_FLAG,                &
      ONLY_TRACERADV_FLAG,             &
      TRACERADV_DISABLE_LIMITER,       &
      TRACERADV_MODALFILTER_FLAG
    
    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh
    class(LocalMeshBase), pointer :: ptr_lcmesh
    class(ElementBase3D), pointer :: elem3D
    integer :: n
    real(DP) :: dtsec

    class(MeshBase3D), pointer :: mesh3D

    integer :: ierr
    !--------------------------------------------------

    if (.not. this%IsActivated()) return
    LOG_INFO('AtmosDyn_setup',*)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN)
    
    !- get mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select
    mesh3D => atm_mesh%ptr_mesh

    !- Setup the temporal integrator

    call tm_parent_comp%Regist_process( 'ATMOS_DYN', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                    ! (out)

    dtsec = tm_parent_comp%process_list(this%tm_process_id)%dtsec
    
    allocate( this%tint(ptr_mesh%LOCAL_MESH_NUM) )
    allocate( this%tint_qtrc(ptr_mesh%LOCAL_MESH_NUM) )

    do n = 1, ptr_mesh%LOCAL_MESH_NUM
      call ptr_mesh%GetLocalMesh( n, ptr_lcmesh )

      call this%tint(n)%Init( TINTEG_TYPE, dtsec, ATMOS_PROGVARS_NUM, 2, &
        (/ ptr_mesh%refElem%Np, ptr_lcmesh%NeA /) )
      
      call this%tint_qtrc(n)%Init( TINTEG_TYPE_TRACER, dtsec, 1, 2, &
        (/ ptr_mesh%refElem%Np, ptr_lcmesh%NeA /) )        
    end do

    !- initialize an object to manage boundary conditions
    call this%boundary_cond%Init()
    call this%boundary_cond%SetBCInfo( ptr_mesh )

    !- initialize the variables 
    call this%dyn_vars%Init( model_mesh )

    call setup_coriolis_parameter( this%dyn_vars, atm_mesh )

    !- Initialize a module for 3D dynamical core 

    this%ENTOT_CONSERVE_SCHEME_FLAG = .false.
    select case(EQS_TYPE)
    case("NONHYDRO3D_HEVE")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_HEVE
      call atm_dyn_dgm_nonhydro3d_heve_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_nonhydro3d_heve_cal_tend
      this%cal_vi => null()
    case("GLOBALNONHYDRO3D_HEVE")
      this%EQS_TYPEID = EQS_TYPEID_GLOBALNONHYD3D_HEVE
      call atm_dyn_dgm_globalnonhydro3d_heve_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_globalnonhydro3d_heve_cal_tend
      this%cal_vi => null()
    case("GLOBALNONHYDRO3D_HEVE_ENTOT")
      this%EQS_TYPEID = EQS_TYPEID_GLOBALNONHYD3D_HEVE_ENTOT
      call atm_dyn_dgm_globalnonhydro3d_etot_heve_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_globalnonhydro3d_etot_heve_cal_tend
      this%cal_vi => null()
      this%ENTOT_CONSERVE_SCHEME_FLAG = .true.
    case("GLOBALNONHYDRO3D_HEVI")
      this%EQS_TYPEID = EQS_TYPEID_GLOBALNONHYD3D_HEVI
      call atm_dyn_dgm_globalnonhydro3d_hevi_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_globalnonhydro3d_hevi_cal_tend
      this%cal_vi => atm_dyn_dgm_globalnonhydro3d_hevi_cal_vi
    case("GLOBALNONHYDRO3D_HEVI_ENTOT")
      this%EQS_TYPEID = EQS_TYPEID_GLOBALNONHYD3D_HEVI_ENTOT
      call atm_dyn_dgm_globalnonhydro3d_etot_hevi_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_tend
      this%cal_vi => atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_vi
      this%ENTOT_CONSERVE_SCHEME_FLAG = .true.
    case("NONHYDRO3D_SPLITFORM_HEVE")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_SPLITFORM_HEVE
      call atm_dyn_dgm_nonhydro3d_heve_splitform_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_nonhydro3d_heve_splitform_cal_tend
      this%cal_vi => null()
    case("NONHYDRO3D_HEVI")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_HEVI
      call atm_dyn_dgm_nonhydro3d_hevi_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_nonhydro3d_hevi_cal_tend
      this%cal_vi => atm_dyn_dgm_nonhydro3d_hevi_cal_vi      
    case("NONHYDRO3D_SPLITFORM_HEVI")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_SPLITFORM_HEVI
      call atm_dyn_dgm_nonhydro3d_hevi_splitform_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_nonhydro3d_hevi_splitform_cal_tend
      this%cal_vi => atm_dyn_dgm_nonhydro3d_hevi_splitform_cal_vi
    case default
      LOG_ERROR("ATMOS_DYN_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
      call PRC_abort
    end select    

    !- Initialize a module for tracer equations
    call atm_dyn_dgm_trcadvect3d_heve_Init( mesh3D, this%FaceIntMat ) 

    !- Setup the numerical diffusion
    this%CALC_NUMDIFF_FLAG = NUMDIFF_FLAG
    if( NUMDIFF_FLAG ) call setup_numdiff( this, atm_mesh )

    !- Setup the modal filter
    this%MODALFILTER_FLAG = MODALFILTER_FLAG
    if ( MODALFILTER_FLAG ) call setup_modalfilter( this, atm_mesh, 'dyn' )

    !- Setup the sponge layer
    this%SPONGELAYER_FLAG = SPONGELAYER_FLAG
    if ( SPONGELAYER_FLAG ) call setup_spongelayer( this, atm_mesh, dtsec )

    !- Setup flags associated with tracer advection
    this%ONLY_TRACERADV_FLAG = ONLY_TRACERADV_FLAG
    this%TRACERADV_disable_limiter = TRACERADV_DISABLE_LIMITER
    this%TRACERADV_MODALFILTER_FLAG = TRACERADV_MODALFILTER_FLAG
    if ( TRACERADV_MODALFILTER_FLAG ) call setup_modalfilter( this, atm_mesh, 'tracer' )

    return
  end subroutine AtmosDyn_setup

!OCL SERIAL  
  subroutine AtmosDyn_calc_tendency( this, model_mesh, prgvars_list, trcvars_list, auxvars_list, forcing_list, is_update )
    implicit none
    
    class(AtmosDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------
    if (.not. this%IsActivated()) return
    !LOG_INFO('AtmosDyn_tendency',*)

    return  
  end subroutine AtmosDyn_calc_tendency

!OCL SERIAL
  subroutine AtmosDyn_update( this, model_mesh, prgvars_list, trcvars_list, auxvars_list, forcing_list, is_update )
    use scale_tracer, only: &
      QA, TRACER_ADVC, TRACER_NAME
    use scale_atm_dyn_dgm_modalfilter, only: &
      atm_dyn_dgm_modalfilter_apply,       &
      atm_dyn_dgm_tracer_modalfilter_apply
    use scale_const, only: &
      GRAV => CONST_GRAV,  &
      Rdry => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    implicit none

    class(AtmosDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    integer :: rkstage
    integer :: tintbuf_ind

    class(MeshBase), pointer :: mesh
    class(LocalMesh3D), pointer :: lcmesh
    integer :: n
    integer :: ke, ke2D, p

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: Coriolis
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: ALPH_DENS_M_tavg, ALPH_DENS_P_tavg, MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg
    class(LocalMeshFieldBase), pointer :: QTRC
    class(LocalMeshFieldBase), pointer :: RHOQ_tp
    class(LocalMeshFieldBase), pointer :: ThermodynVar


!    class(LocalMeshFieldBase), pointer :: MOMZ_t, MOMZ_t_advx, MOMZ_t_advY, MOMZ_t_advZ, MOMZ_t_lift, MOMZ_t_buoy
    integer :: v
    integer :: iq
    integer :: nRKstage
    real(RP) :: implicit_fac
    real(RP) :: dt
    real(RP) :: dttmp_trc

    real(RP) :: tavg_coef_MFLXZ(this%tint(1)%nstage)
    !--------------------------------------------------
    
    call PROF_rapstart( 'ATM_DYN_update', 1)   

    call model_mesh%GetModelMesh( mesh )
    
    if ( this%ONLY_TRACERADV_FLAG ) then
      nRKstage = -1
    else
      nRKstage = this%tint(1)%nstage
    end if

    !-
    do rkstage=1, nRKstage

      if ( this%ENTOT_CONSERVE_SCHEME_FLAG ) then
        do n=1, mesh%LOCAL_MESH_NUM
          call AtmosVars_GetLocalMeshPrgVars( n, &
            mesh, prgvars_list, auxvars_list,                               &
            DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
            DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh                  )
  
          call cal_DRHOT2Entot( this%dyn_vars%EnTot%local(n)%val, &
            DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,             &
            PRES_hyd%val, DENS_hyd%val, CPtot%val, CVtot%val, Rtot%val,     &
            lcmesh, lcmesh%refElem3D                                        ) 
        end do
      end if

      if (this%tint(1)%imex_flag) then        
        do n=1, mesh%LOCAL_MESH_NUM
          call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
          call AtmosVars_GetLocalMeshPrgVars( n, &
            mesh, prgvars_list, auxvars_list,                               &
            DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
            DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh                  )
          call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)   

          if ( this%ENTOT_CONSERVE_SCHEME_FLAG ) then
            call this%dyn_vars%EnTot%GetLocalMeshField( n, ThermodynVar )
          else
            ThermodynVar => DRHOT
          end if
          if (rkstage==1) then
            call this%tint(n)%StoreVar0( DDENS%val, DDENS_ID,    &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
            call this%tint(n)%StoreVar0( MOMX%val, MOMX_ID,      &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
            call this%tint(n)%StoreVar0( MOMY%val, MOMY_ID,      &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
            call this%tint(n)%StoreVar0( MOMZ%val, MOMZ_ID,      &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
            call this%tint(n)%StoreVar0( ThermodynVar%val, THERM_ID, &
                      1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE   )
          end if

          call PROF_rapstart( 'ATM_DYN_cal_vi', 2)
          implicit_fac = this%tint(n)%Get_implicit_diagfac(rkstage)
          tintbuf_ind = this%tint(n)%tend_buf_indmap(rkstage)
          dt = this%tint(n)%Get_deltime()
          call this%cal_vi( &
            this%tint(n)%tend_buf2D_im(:,:,DDENS_ID,tintbuf_ind),                   & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,MOMX_ID ,tintbuf_ind),                   & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,MOMY_ID ,tintbuf_ind),                   & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,MOMZ_ID ,tintbuf_ind),                   & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,THERM_ID,tintbuf_ind),                   & ! (out)
            DDENS%val, MOMX%val, MOMY%val, MOMZ%val,                                & ! (in)
            ThermodynVar%val,                                                       & ! (in)
            DENS_hyd%val, PRES_hyd%val,                                             & ! (in)
            this%tint(n)%var0_2D(:,:,DDENS_ID), this%tint(n)%var0_2D(:,:,MOMX_ID),  & ! (in)
            this%tint(n)%var0_2D(:,:,MOMY_ID ), this%tint(n)%var0_2D(:,:,MOMZ_ID),  & ! (in)
            this%tint(n)%var0_2D(:,:,THERM_ID ),                                    & ! (in)
            Rtot%val, CVtot%val, CPtot%val,                                         & ! (in)
            model_mesh%DOptrMat(3), model_mesh%LiftOptrMat,                         & ! (in)
            this%MODALFILTER_FLAG, this%modal_filter_v1D,                           & ! (in)
            implicit_fac, dt,                                                       & ! (in)
            lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D    ) ! (in)
          
          call PROF_rapend( 'ATM_DYN_cal_vi', 2)  
          
          call PROF_rapstart( 'ATM_DYN_store_impl', 2)      
          call this%tint(n)%StoreImplicit( rkstage, DDENS%val, DDENS_ID,  &
                             1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
          
          call this%tint(n)%StoreImplicit( rkstage, MOMX%val, MOMX_ID,    &
                             1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
          
          call this%tint(n)%StoreImplicit( rkstage, MOMY%val, MOMY_ID,    &
                             1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )

          call this%tint(n)%StoreImplicit( rkstage, MOMZ%val, MOMZ_ID,    &
                             1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )

          call this%tint(n)%StoreImplicit( rkstage, ThermodynVar%val, THERM_ID, &
                             1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )

          if (this%ENTOT_CONSERVE_SCHEME_FLAG) then
            call cal_EnTot2DRHOT( DRHOT%val, &
              DDENS%val, MOMX%val, MOMY%val, MOMZ%val, this%dyn_vars%EnTot%local(n)%val,  &
              PRES_hyd%val, DENS_hyd%val, CPtot%val, CVtot%val, Rtot%val,                 &
              lcmesh, lcmesh%refElem3D                                                    ) 
          end if
          call PROF_rapend( 'ATM_DYN_store_impl', 2) 
        end do
      end if

      !* Exchange halo data
      call PROF_rapstart( 'ATM_DYN_exchange_prgv', 2)
      call prgvars_list%MeshFieldComm_Exchange()
      call PROF_rapend( 'ATM_DYN_exchange_prgv', 2)
  
      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
        call AtmosVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list,           &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,             &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot,     &
          lcmesh                                      )
        call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)         
        
        !* Apply boundary conditions
        call PROF_rapstart( 'ATM_DYN_applyBC_prgv', 2)
        call this%boundary_cond%ApplyBC_PROGVARS_lc( n,                                & ! (in)
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                          & ! (inout)
          DENS_hyd%val, PRES_hyd%val,                                                  & ! (in)
          lcmesh%Gsqrt(:,:), lcmesh%GsqrtH(:,:), lcmesh%GI3(:,:,1), lcmesh%GI3(:,:,2), & ! (in)
          lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3),   & ! (in)
          lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB,                                    & ! (in)
          lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D         ) ! (in)
        call PROF_rapend( 'ATM_DYN_applyBC_prgv', 2)
      end do


      do n=1, mesh%LOCAL_MESH_NUM
        tintbuf_ind = this%tint(n)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
        call AtmosVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list,        &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,          &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot,  &
          lcmesh                                   )
        
        call AtmosDynAuxVars_GetLocalMeshFields( n,      &
          mesh, this%dyn_vars%AUXVARS2D_manager,         &
          Coriolis )

        if ( this%ENTOT_CONSERVE_SCHEME_FLAG ) then
          call this%dyn_vars%EnTot%GetLocalMeshField( n, ThermodynVar )
        else
          ThermodynVar => DRHOT
        end if          
        call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)

        call PROF_rapstart( 'ATM_DYN_update_caltend_ex', 2)
        call this%cal_tend_ex( &
          this%tint(n)%tend_buf2D_ex(:,:,DDENS_ID,tintbuf_ind),                   &
          this%tint(n)%tend_buf2D_ex(:,:,MOMX_ID ,tintbuf_ind),                   &
          this%tint(n)%tend_buf2D_ex(:,:,MOMY_ID ,tintbuf_ind),                   &
          this%tint(n)%tend_buf2D_ex(:,:,MOMZ_ID ,tintbuf_ind),                   &
          this%tint(n)%tend_buf2D_ex(:,:,DRHOT_ID,tintbuf_ind),                   &
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                     &
          DENS_hyd%val, PRES_hyd%val,                                             &
          Coriolis%val,                                                           &
          Rtot%val, CVtot%val, CPtot%val,                                         & 
          this%SPONGELAYER_FLAG, this%wdamp_tau, this%wdamp_height,               &
          this%hvel_damp_flag,                                                    &
          model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3), &
          model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3), &
          model_mesh%LiftOptrMat,                                                 &
          lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D ) 
        call PROF_rapend( 'ATM_DYN_update_caltend_ex', 2)

        call PROF_rapstart( 'ATM_DYN_update_add_tp', 2)
        call add_phy_tend( &
          this, this%tint(n)%tend_buf2D_ex(:,:,:,tintbuf_ind), & ! (inout)
          DRHOT%val, PRES_hyd%val,                             & ! (in)
          Rtot%val, CVtot%val, CPtot%val,                      & ! (in)
          forcing_list,                                        & ! (in)
          mesh, n, lcmesh, lcmesh%refElem3D                    ) ! (in)
        call PROF_rapend( 'ATM_DYN_update_add_tp', 2)

        if ( QA > 0 ) then
          call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
          call AtmosDynMassFlux_GetLocalMeshFields( n, mesh, this%dyn_vars,            &
            ALPH_DENS_M_tavg, ALPH_DENS_P_tavg, MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg  )
          call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)  

          call PROF_rapstart( 'ATM_DYN_tavg_mflx', 2)
          if (this%tint(1)%imex_flag) then
            tavg_coef_MFLXZ(:) = this%tint(n)%coef_b_im(:)
          else
            tavg_coef_MFLXZ(:) = this%tint(n)%coef_b_ex(:)
          end if
          call atm_dyn_dgm_trcadvect3d_save_massflux( &
            MFLX_x_tavg%val, MFLX_y_tavg%val, MFLX_z_tavg%val,                 & ! (inout)
            ALPH_DENS_M_tavg%face_val, ALPH_DENS_P_tavg%face_val,              & ! (inout)
            DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                & ! (in)
            DENS_hyd%val, PRES_hyd%val,                                        & ! (in)
            Rtot%val, CVtot%val, CPtot%val,                                    & ! (in)
            lcmesh, lcmesh%refElem3D,                                          & ! (in)
            rkstage, this%tint(n)%coef_b_ex(rkstage), tavg_coef_MFLXZ(rkstage) ) ! (in)
          call PROF_rapend( 'ATM_DYN_tavg_mflx', 2)
        end if
        
        call PROF_rapstart( 'ATM_DYN_update_advance', 2)      
        call this%tint(n)%Advance( rkstage, DDENS%val, DDENS_ID, &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        call this%tint(n)%Advance( rkstage, MOMX%val, MOMX_ID,   &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        call this%tint(n)%Advance( rkstage, MOMY%val, MOMY_ID,   &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )

        call this%tint(n)%Advance( rkstage, MOMZ%val, MOMZ_ID,   &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )

        call this%tint(n)%Advance( rkstage, ThermodynVar%val, THERM_ID, &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE        )

        if (this%ENTOT_CONSERVE_SCHEME_FLAG) then
          call cal_EnTot2DRHOT( DRHOT%val, &
            DDENS%val, MOMX%val, MOMY%val, MOMZ%val, this%dyn_vars%EnTot%local(n)%val, &
            PRES_hyd%val, DENS_hyd%val, CPtot%val, CVtot%val, Rtot%val,                &
            lcmesh, lcmesh%refElem3D                                                   ) 
        end if
        call PROF_rapend( 'ATM_DYN_update_advance', 2)
      end do
    end do


    !-- Tracer advection (prepair) ------------------------------------------------
    

    if ( QA > 0 ) then
      call PROF_rapstart( 'ATM_DYN_qtracer', 2)

      if ( this%ONLY_TRACERADV_FLAG ) then
        do n=1, mesh%LOCAL_MESH_NUM
          call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 3) 
          call AtmosVars_GetLocalMeshPrgVars( n, &
            mesh, prgvars_list, auxvars_list,       &
            DDENS, MOMX, MOMY, MOMZ, DRHOT,         &
            DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
            lcmesh                                  )
                       
          call AtmosDynMassFlux_GetLocalMeshFields( n, mesh, this%dyn_vars,            &
            ALPH_DENS_M_tavg, ALPH_DENS_P_tavg, MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg  )
          call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 3)
          !$omp parallel do
          do ke=lcmesh%NeS, lcmesh%NeE
            MFLX_x_tavg%val(:,ke) = MOMX%val(:,ke)
            MFLX_y_tavg%val(:,ke) = MOMY%val(:,ke)
            MFLX_z_tavg%val(:,ke) = MOMZ%val(:,ke)
            this%tint(n)%var0_2D(:,ke,DDENS_ID) = DDENS%val(:,ke)           
          end do
          call atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_advtest( &
            ALPH_DENS_M_tavg%face_val, ALPH_DENS_P_tavg%face_val,                      & ! (inout)
            DDENS%val, MOMX%val, MOMY%val, MOMz%val, DENS_hyd%val,                     & ! (in)
            lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), & ! (in)
            lcmesh%VMapM, lcmesh%VMapP, lcmesh, lcmesh%refElem3D                       ) ! (in)
        end do
      end if

      !* Exchange halo data of mass flux

      call PROF_rapstart( 'ATM_DYN_exchange_mflx', 3)
      call this%dyn_vars%AUXTRC_FLUX_VAR3D_manager%MeshFieldComm_Exchange()
      call PROF_rapend( 'ATM_DYN_exchange_mflx', 3)

      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 3) 
        call AtmosVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list,       &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,         &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
          lcmesh                                  )
                   
        call AtmosDynMassFlux_GetLocalMeshFields( n, mesh, this%dyn_vars,            &
          ALPH_DENS_M_tavg, ALPH_DENS_P_tavg, MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg  )
        call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 3)

        call PROF_rapstart( 'ATM_DYN_applyBC_mflux', 3)
        !$omp parallel do
        do ke=lcmesh%NeS, lcmesh%NeE
          this%dyn_vars%TRCVARS3D(TRCDDENS_ID)%local(n)%val(:,ke) = DDENS%val(:,ke)
        end do

        call this%boundary_cond%ApplyBC_PROGVARS_lc( n,                                & ! (in)
          this%dyn_vars%TRCVARS3D(TRCDDENS_ID)%local(n)%val(:,:),                      & ! (inout)
          MFLX_x_tavg%val, MFLX_y_tavg%val, MFLX_z_tavg%val, DRHOT%val,                & ! (inout)
          DENS_hyd%val, PRES_hyd%val,                                                  & ! (in)
          lcmesh%Gsqrt(:,:), lcmesh%GsqrtH(:,:), lcmesh%GI3(:,:,1), lcmesh%GI3(:,:,2), & ! (in)
          lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3),   & ! (in)
          lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB,                                    & ! (in)
          lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D         ) ! (in)
        call PROF_rapend( 'ATM_DYN_applyBC_mflux', 3)
      end do

      call PROF_rapend( 'ATM_DYN_qtracer', 2)
    end if

    !-- modal filter  -----------------------------------------------------------

    if ( this%MODALFILTER_FLAG ) then
      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
        call AtmosVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list,                               &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh                  )
        call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)

        call PROF_rapstart( 'ATM_DYN_update_modalfilter', 2)
        call atm_dyn_dgm_modalfilter_apply(  & 
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val, & ! (inout)
          lcmesh, lcmesh%refElem3D, this%modal_filter_3d,     & ! (in)
!          do_weight_Gsqrt = .false.                          ) ! (in)          
          do_weight_Gsqrt = .true.                            ) ! (in)
        call PROF_rapend( 'ATM_DYN_update_modalfilter', 2)
      end do
    end if  

    !-- Tracer advection ------------------------------------------------

    if ( QA > 0 ) then
      call PROF_rapstart( 'ATM_DYN_qtracer', 2)

      do iq=1, QA
        do n=1, mesh%LOCAL_MESH_NUM
          call PROF_rapstart( 'ATM_DYN_get_localmesh_qtrc', 3)
          call AtmosVars_GetLocalMeshQTRCVar( n,       &
            mesh, trcvars_list, iq,                    &
            QTRC, lcmesh                               )
          !$omp parallel do
          do ke=lcmesh%NeS, lcmesh%NeE
            this%dyn_vars%TRCVARS3D(TRCQ_ID)%local(n)%val(:,ke) = QTRC%val(:,ke)
          end do            
          call PROF_rapend( 'ATM_DYN_get_localmesh_qtrc', 3)            
        end do

        do rkstage=1, this%tint_qtrc(1)%nstage

          if ( TRACER_ADVC(iq) ) then
            call PROF_rapstart( 'ATM_DYN_exchange_qtrc', 3)
            call this%dyn_vars%TRCVAR3D_manager%MeshFieldComm_Exchange()
            call PROF_rapend( 'ATM_DYN_exchange_qtrc', 3)

            do n=1, mesh%LOCAL_MESH_NUM
              call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 3)

              call AtmosVars_GetLocalMeshPrgVar( n, &
                mesh, prgvars_list, auxvars_list,   &
                DDENS_ID, DDENS, DENS_hyd           )  

              call AtmosDynMassFlux_GetLocalMeshFields( n, mesh, this%dyn_vars,            &
                ALPH_DENS_M_tavg, ALPH_DENS_P_tavg, MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg  )
        
              call AtmosVars_GetLocalMeshQTRCPhyTend( n, &
                mesh, forcing_list, iq,                  &
                RHOQ_tp                                  )
  
              call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 3)              

              !---
              call PROF_rapstart( 'ATM_DYN_calc_fct_coef', 3)  

              dt = this%tint_qtrc(n)%Get_deltime()
              dttmp_trc = dt * this%tint_qtrc(n)%coef_gam_ex(rkstage+1,rkstage) &
                            / this%tint_qtrc(n)%coef_sig_ex(rkstage+1,rkstage)
              call atm_dyn_dgm_trcadvect3d_heve_calc_fct_coef( &
                this%dyn_vars%AUX_TRCVARS3D(1)  %local(n)%val,                             & ! (out)
                this%dyn_vars%TRCVARS3D(TRCQ_ID)%local(n)%val,                             & ! (in)
                MFLX_x_tavg%val, MFLX_y_tavg%val, MFLX_z_tavg%val, RHOQ_tp%val,            & ! (in)
                ALPH_DENS_M_tavg%face_val, ALPH_DENS_P_tavg%face_val,                      & ! (in)
                DENS_hyd%val(:,:), this%dyn_vars%TRCVARS3D(TRCDDENS_ID)%local(n)%val(:,:), & ! (in)
                this%tint(n)%var0_2D(:,:,DDENS_ID),                                        & ! (in)
                this%tint_qtrc(n)%coef_c_ex(rkstage), dttmp_trc,                           & ! (in) 
                model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),    & ! (in)
                model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3),    & ! (in)
                model_mesh%LiftOptrMat, this%FaceIntMat,                                   & ! (in)
                lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D,      & ! (in)
                this%TRACERADV_disable_limiter                                             ) ! (in)

              call PROF_rapend( 'ATM_DYN_calc_fct_coef', 3)  
            end do

            call PROF_rapstart( 'ATM_DYN_exchange_qtrc', 3)
            call this%dyn_vars%AUXTRCVAR3D_manager%MeshFieldComm_Exchange()
            call PROF_rapend( 'ATM_DYN_exchange_qtrc', 3)
          end if
          do n=1, mesh%LOCAL_MESH_NUM
            tintbuf_ind = this%tint_qtrc(n)%tend_buf_indmap(rkstage)

            call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 3)              
            call AtmosDynMassFlux_GetLocalMeshFields( n, mesh, this%dyn_vars,            &
              ALPH_DENS_M_tavg, ALPH_DENS_P_tavg, MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg  )
            
            call AtmosVars_GetLocalMeshPrgVar( n, &
              mesh, prgvars_list, auxvars_list,   &
              DDENS_ID, DDENS, DENS_hyd           )

            call AtmosVars_GetLocalMeshQTRCPhyTend( n, &
              mesh, forcing_list, iq,                  &
              RHOQ_tp                                  )
            call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 3)

            call PROF_rapstart( 'ATM_DYN_update_caltend_ex_qtrc', 3)  
            if ( TRACER_ADVC(iq) ) then 
              call atm_dyn_dgm_trcadvect3d_heve_cal_tend( &        
                this%tint_qtrc(n)%tend_buf2D_ex(:,:,1,tintbuf_ind),                     & ! (out)
                this%dyn_vars%TRCVARS3D(TRCQ_ID)%local(n)%val,                          & ! (in)
                MFLX_x_tavg%val, MFLX_y_tavg%val, MFLX_z_tavg%val,                      & ! (in)
                this%dyn_vars%alphaDensM%local(n)%face_val,                             & ! (in)
                this%dyn_vars%alphaDensP%local(n)%face_val,                             & ! (in)
                this%dyn_vars%AUX_TRCVARS3D(1)%local(n)%val,                            & ! (in)
                RHOQ_tp%val,                                                            & ! (in) 
                model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3), & ! (in)
                model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3), & ! (in)
                model_mesh%LiftOptrMat, this%FaceIntMat,                                & ! (in)
                lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D    ) ! (in)
            else
              !$omp parallel do
              do ke=lcmesh%NeS, lcmesh%NeE
                this%tint_qtrc(n)%tend_buf2D_ex(:,ke,1,tintbuf_ind) = RHOQ_tp%val(:,ke)
              end do
            end if
            call PROF_rapend( 'ATM_DYN_update_caltend_ex_qtrc', 3)

            call PROF_rapstart( 'ATM_DYN_update_advance_qtrc', 3)                
            call this%tint_qtrc(n)%Advance_trcvar( &
              rkstage, this%dyn_vars%TRCVARS3D(1)%local(n)%val, 1,    &
              1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE,           &
              this%dyn_vars%TRCVARS3D(TRCDDENS_ID)%local(n)%val(:,:), &
              this%tint(n)%var0_2D(:,:,DDENS_ID), DENS_hyd%val(:,:)   ) 
            call PROF_rapend( 'ATM_DYN_update_advance_qtrc', 3)

            if ( rkstage == this%tint_qtrc(1)%nstage         &
                .and. this%TRACERADV_MODALFILTER_FLAG        ) then
              call PROF_rapstart( 'ATM_DYN_update_qtrc_modalfilter', 3)
              call atm_dyn_dgm_tracer_modalfilter_apply( &
                this%dyn_vars%TRCVARS3D(TRCQ_ID)%local(n)%val(:,:),                        & ! (inout)
                DENS_hyd%val(:,:), this%dyn_vars%TRCVARS3D(TRCDDENS_ID)%local(n)%val(:,:), & ! (in)
                lcmesh, lcmesh%refElem3D, this%modal_filter_tracer_3d                      ) ! (in)
              call PROF_rapend( 'ATM_DYN_update_qtrc_modalfilter', 3)
            end if
            if ( TRACER_ADVC(iq)                             &
              .and. rkstage == this%tint_qtrc(1)%nstage      &
              .and. ( .not. this%TRACERADV_disable_limiter ) ) then
              call PROF_rapstart( 'ATM_DYN_update_qtrc_TMAR', 3)             
              call atm_dyn_dgm_trcadvect3d_TMAR( &
                this%dyn_vars%TRCVARS3D(TRCQ_ID)%local(n)%val(:,:),                        & ! (inout)
                DENS_hyd%val(:,:), this%dyn_vars%TRCVARS3D(TRCDDENS_ID)%local(n)%val(:,:), & ! (in)
                lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D       ) ! (in)
              call PROF_rapend( 'ATM_DYN_update_qtrc_TMAR', 3)  
            end if

          end do
        end do

        do n=1, mesh%LOCAL_MESH_NUM
          call PROF_rapstart( 'ATM_DYN_get_localmesh_qtrc', 3)
          call AtmosVars_GetLocalMeshQTRCVar( n,       &
            mesh, trcvars_list, iq,                    &
            QTRC, lcmesh                               )
          call AtmosVars_GetLocalMeshPrgVar( n, &
            mesh, prgvars_list, auxvars_list,   &
            DDENS_ID, DDENS, DENS_hyd           )
          call PROF_rapend( 'ATM_DYN_get_localmesh_qtrc', 3)            

         !$omp parallel do
          do ke=lcmesh%NeS, lcmesh%NeE
            QTRC%val(:,ke) = ( DENS_hyd%val(:,ke) + this%dyn_vars%TRCVARS3D(TRCDDENS_ID)%local(n)%val(:,ke) ) &
                           / ( DENS_hyd%val(:,ke) + DDENS%val(:,ke) )                                         &
                           * this%dyn_vars%TRCVARS3D(TRCQ_ID)%local(n)%val(:,ke)
          end do            
        end do 

      end do ! end do for iq

      call PROF_rapend( 'ATM_DYN_qtracer', 2)      
    end if

    !-- numerical diffusion for dynamical variables -----------------------------

    if ( this%CALC_NUMDIFF_FLAG ) then
      
      call PROF_rapstart( 'ATM_DYN_numfilter', 2)
      call prgvars_list%MeshFieldComm_Exchange()

      do v = 1, ATMOS_PROGVARS_NUM
        call cal_numfilter_tend( this, model_mesh, prgvars_list, auxvars_list, v )
      end do

      do n=1, mesh%LOCAL_MESH_NUM
        dt = this%tint(n)%Get_deltime()

        call AtmosVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list,       &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,         &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
          lcmesh                                  )
        !$omp parallel do
        do ke=1, lcmesh%Ne
         DDENS%val(:,ke) = DDENS%val(:,ke) + dt * this%tint(n)%tend_buf2D_ex(:,ke,DDENS_ID,1)
         MOMX %val(:,ke) = MOMX %val(:,ke) + dt * this%tint(n)%tend_buf2D_ex(:,ke,MOMX_ID ,1)
         MOMY %val(:,ke) = MOMY %val(:,ke) + dt * this%tint(n)%tend_buf2D_ex(:,ke,MOMY_ID ,1)
         MOMZ %val(:,ke) = MOMZ %val(:,ke) + dt * this%tint(n)%tend_buf2D_ex(:,ke,MOMZ_ID ,1)
         DRHOT%val(:,ke) = DRHOT%val(:,ke) + dt * this%tint(n)%tend_buf2D_ex(:,ke,DRHOT_ID,1)
        end do
      end do

      call PROF_rapend( 'ATM_DYN_numfilter', 2)
    end if

    !---------------------------
    call PROF_rapend( 'ATM_DYN_update', 1)

    return  
  end subroutine AtmosDyn_update

!OCL SERIAL
  subroutine AtmosDyn_finalize( this )
    implicit none
    class(AtmosDyn), intent(inout) :: this

    integer :: n
    !--------------------------------------------------
    if (.not. this%IsActivated()) return
    LOG_INFO('AtmosDyn_finalize',*)

    select case(this%EQS_TYPEID)
    case(EQS_TYPEID_NONHYD3D_HEVE)
      call atm_dyn_dgm_nonhydro3d_heve_Final()
    case(EQS_TYPEID_GLOBALNONHYD3D_HEVE)
      call atm_dyn_dgm_globalnonhydro3d_heve_Final()
    case(EQS_TYPEID_GLOBALNONHYD3D_HEVI)
      call atm_dyn_dgm_globalnonhydro3d_hevi_Final()
    case(EQS_TYPEID_NONHYD3D_HEVI)  
      call atm_dyn_dgm_nonhydro3d_hevi_Final()
    case(EQS_TYPEID_NONHYD3D_SPLITFORM_HEVI)  
      call atm_dyn_dgm_nonhydro3d_hevi_splitform_Final()     
    end select 

    call atm_dyn_dgm_trcadvect3d_heve_Final()

    if (this%CALC_NUMDIFF_FLAG) then
      call atm_dyn_dgm_nonhydro3d_numdiff_Final()
    end if

    if (this%MODALFILTER_FLAG) then
      call this%modal_filter_3d%Final()
      if ( associated(this%cal_vi) ) call this%modal_filter_v1D%Final()
    end if
    
    do n = 1, size(this%tint)
      call this%tint(n)%Final()
    end do
    deallocate( this%tint )

    do n = 1, size(this%tint_qtrc)
      call this%tint_qtrc(n)%Final()
    end do
    deallocate( this%tint_qtrc )    

    call this%boundary_cond%Final()
    call this%dyn_vars%Final()

    return  
  end subroutine AtmosDyn_finalize  

  !--- private ---------------

!OCL SERIAL
  subroutine add_phy_tend( this,      & ! (in)
    dyn_tends,                        & ! (inout)
    DRHOT, PRES_hyd,                  & ! (in)
    Rtot, CVtot, CPtot,               & ! (in)
    phytends_list,                    & ! (in)
    mesh, domID, lcmesh, elem3D       ) ! (in)

    use scale_const, only: &
      Rdry => CONST_Rdry,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPhyTends

    implicit none

    class(AtmosDyn), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(inout) :: dyn_tends(elem3D%Np,lcmesh%NeA,ATMOS_PROGVARS_NUM)
    real(RP), intent(in) :: DRHOT(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA)
    class(ModelVarManager), intent(inout) :: phytends_list
    class(MeshBase), intent(in) :: mesh
    integer, intent(in) :: domID

    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p

    integer :: ke
    integer :: iq

    real(RP) :: RHOT(elem3D%Np)
    real(RP) :: EXNER(elem3D%Np)

    real(RP) :: rgamm    
    real(RP) :: rP0
    real(RP) :: P0ovR
    !---------------------------------------------------------------------------------

    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    P0ovR = PRES00 / Rdry 

    call AtmosVars_GetLocalMeshPhyTends( domID, mesh, phytends_list, & ! (in)
      DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p            ) ! (out)

    !$omp parallel do          &
    !$omp private( RHOT, EXNER )
    do ke=lcmesh%NeS, lcmesh%NeE
      RHOT(:) = P0ovR * (PRES_hyd(:,ke) * rP0)**rgamm + DRHOT(:,ke)
      EXNER(:) = ( Rtot(:,ke) * rP0 * RHOT(:) )**( Rtot(:,ke) / CVtot(:,ke) )

      dyn_tends(:,ke,DDENS_ID) = dyn_tends(:,ke,DDENS_ID) + DENS_tp%val(:,ke)
      dyn_tends(:,ke,MOMX_ID ) = dyn_tends(:,ke,MOMX_ID ) + MOMX_tp%val(:,ke)
      dyn_tends(:,ke,MOMY_ID ) = dyn_tends(:,ke,MOMY_ID ) + MOMY_tp%val(:,ke)
      dyn_tends(:,ke,MOMZ_ID ) = dyn_tends(:,ke,MOMZ_ID ) + MOMZ_tp%val(:,ke)
      dyn_tends(:,ke,DRHOT_ID) = dyn_tends(:,ke,DRHOT_ID) + RHOT_tp%val(:,ke) &
                               + RHOH_p %val(:,ke) / ( CPtot(:,ke) * EXNER(:) )
    end do

    return
  end subroutine add_phy_tend

!OCL SERIAL
  subroutine cal_numfilter_tend( this, model_mesh, prgvars_list, auxvars_list, varid )

    use mod_atmos_dyn_vars, only: &
      AtmosDynAuxVars_GetLocalMeshFields,     &
      AtmosDynNumDiffFlux_GetLocalMeshFields, &
      AtmosDynNumDiffTend_GetLocalMeshFields
    
    use scale_atm_dyn_dgm_nonhydro3d_numdiff, only: &
      atm_dyn_dgm_nonhydro3d_numdiff_tend,          &
      atm_dyn_dgm_nonhydro3d_numdiff_cal_laplacian, &
      atm_dyn_dgm_nonhydro3d_numdiff_cal_flx

    implicit none
          
    class(AtmosDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    integer, intent(in) :: varid

    class(LocalMeshFieldBase), pointer :: var
    class(LocalMeshFieldBase), pointer :: ND_flx_x, ND_flx_y, ND_flx_z
    class(LocalMeshFieldBase), pointer :: ND_lapla_h, ND_lapla_v
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot

    class(MeshBase), pointer :: mesh
    class(LocalMesh3D), pointer :: lcmesh
    integer :: n
    integer :: nd_itr
    real(RP) :: nd_sign
    logical :: dens_weight_flag
    logical, allocatable :: is_bound(:,:)

    !-----------------------------------------

    nd_sign = (-1)**(mod(this%ND_LAPLACIAN_NUM+1,2))
    dens_weight_flag = (varid /= DDENS_ID)

    call model_mesh%GetModelMesh( mesh )

    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVar( n, mesh, prgvars_list, auxvars_list, &
        varid, var,                                                           &
        DENS_hyd, PRES_hyd, lcmesh                                            )
      call AtmosVars_GetLocalMeshPrgVars( n, mesh, prgvars_list, auxvars_list, &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                                        &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot                                 )
      call AtmosDynNumDiffFlux_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_FLUX_manager, &
        ND_flx_x, ND_flx_y, ND_flx_z )
      
      allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
      call this%boundary_cond%ApplyBC_numdiff_even_lc( var%val, is_bound, varid, n, &
        MOMX%val, MOMY%val, MOMZ%val, DENS_hyd%val, PRES_hyd%val,                    &
        lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3),   &
        lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )
      
      call atm_dyn_dgm_nonhydro3d_numdiff_cal_flx( ND_flx_x%val, ND_flx_y%val, ND_flx_z%val, &
        var%val, var%val, DDENS%val, DENS_hyd%val,                                       &
        model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),          &
        model_mesh%LiftOptrMat,                                                          &
        lcmesh, lcmesh%refElem3D, is_bound, dens_weight_flag ) 

      deallocate( is_bound )
    end do

    !* Exchange halo data
    call this%dyn_vars%NUMDIFF_FLUX_manager%MeshFieldComm_Exchange()

    do nd_itr=1, this%ND_LAPLACIAN_NUM-1
      do n = 1, mesh%LOCAL_MESH_NUM
        call AtmosDynNumDiffFlux_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_FLUX_manager, &
          ND_flx_x, ND_flx_y, ND_flx_z, lcmesh)
        call AtmosDynNumDiffTend_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_TEND_manager, &
          ND_lapla_h, ND_lapla_v )
          
        allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
        call this%boundary_cond%ApplyBC_numdiff_odd_lc( &
          ND_flx_x%val, ND_flx_y%val, ND_flx_z%val, is_bound, varid, n,              &
          lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), &
          lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )

        call atm_dyn_dgm_nonhydro3d_numdiff_cal_laplacian( ND_lapla_h%val, ND_lapla_v%val, &
          ND_flx_x%val, ND_flx_y%val, ND_flx_z%val,                                        &
          model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),          &
          model_mesh%LiftOptrMat,                                                          &
          lcmesh, lcmesh%refElem3D, is_bound )
        
        deallocate( is_bound )
      end do
      !* Exchange halo data
      call this%dyn_vars%NUMDIFF_TEND_manager%MeshFieldComm_Exchange()

      do n = 1, mesh%LOCAL_MESH_NUM
        call AtmosDynNumDiffFlux_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_FLUX_manager, &
          ND_flx_x, ND_flx_y, ND_flx_z, lcmesh)
        call AtmosDynNumDiffTend_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_TEND_manager, &
          ND_lapla_h, ND_lapla_v )    
        call AtmosVars_GetLocalMeshPrgVars( n, mesh, prgvars_list, auxvars_list, &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                                        &
          DENS_hyd, PRES_hyd, Rtot, CPtot, CVtot                                 )
          
        allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
        call this%boundary_cond%ApplyBC_numdiff_even_lc( &
          ND_lapla_h%val, is_bound, varid, n,                                        &
          MOMX%val, MOMY%val, MOMZ%val, DENS_hyd%val, PRES_hyd%val,                  &
          lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), &
          lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )
          
        call atm_dyn_dgm_nonhydro3d_numdiff_cal_flx( ND_flx_x%val, ND_flx_y%val, ND_flx_z%val, &
          ND_lapla_h%val, ND_lapla_v%val, DDENS%val, DENS_hyd%val,                             &
          model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),              &
          model_mesh%LiftOptrMat,                                                              &
          lcmesh, lcmesh%refElem3D, is_bound, .false. ) 

        deallocate( is_bound )
      end do
      !* Exchange halo data
      call this%dyn_vars%NUMDIFF_FLUX_manager%MeshFieldComm_Exchange()
    end do
    
    do n = 1, mesh%LOCAL_MESH_NUM

      call AtmosDynNumDiffFlux_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_FLUX_manager, &
        ND_flx_x, ND_flx_y, ND_flx_z, lcmesh)  
      call AtmosVars_GetLocalMeshPrgVar( n, mesh, prgvars_list, auxvars_list,  &
        varid, var,                                                            &
        DENS_hyd, PRES_hyd, lcmesh                                             )
      call AtmosVars_GetLocalMeshPrgVar( n, mesh, prgvars_list, auxvars_list,  &
        DDENS_ID, DDENS                                                        )

      allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
      call this%boundary_cond%ApplyBC_numdiff_odd_lc(                              &
        ND_flx_x%val, ND_flx_y%val, ND_flx_z%val, is_bound, varid, n,              &
        lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), &
        lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )

      call atm_dyn_dgm_nonhydro3d_numdiff_tend( this%tint(n)%tend_buf2D_ex(:,:,varid,1),  &
        ND_flx_x%val, ND_flx_y%val, ND_flx_z%val,                                         &
        DDENS%val, DENS_hyd%val, nd_sign * this%ND_COEF_H, nd_sign * this%ND_COEF_V,      &
        model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),           &
        model_mesh%LiftOptrMat,                                                           &
        lcmesh, lcmesh%refElem3D, is_bound, dens_weight_flag ) 

      deallocate( is_bound )
    end do

    return
  end subroutine cal_numfilter_tend

!OCL SERIAL
  subroutine cal_DRHOT2Entot( EnTot,  &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,         &
    PRES_hyd, DENS_hyd, CPtot, CVtot, Rtot, &
    lcmesh, elem3D )

    use scale_const, only: &
      GRAV => CONST_GRAV,  &
      Rdry => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00

    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: EnTot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA)

    integer :: ke, ke2D

    real(RP) :: DENS(elem3D%Np)
    real(RP) :: mom_u1(elem3D%Np), mom_u2(elem3D%Np)
    real(RP) :: RHOT(elem3D%Np), PRES(elem3D%Np)
    !---------------------------------------------------------------

    !$omp parallel do private( ke2D, DENS, mom_u1, mom_u2, RHOT, PRES )
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      mom_u1(:) = lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,1,1) * MOMX(:,ke) + lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,1) * MOMY(:,ke)
      mom_u2(:) = lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,1) * MOMX(:,ke) + lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,2) * MOMY(:,ke)

      RHOT(:) = PRES00 / Rdry * ( PRES_hyd(:,ke) / PRES00 )**(CvDry/CpDry) + DRHOT(:,ke)
      PRES(:) = PRES00 * ( Rtot(:,ke) / PRES00 * RHOT(:) )**( CPtot(:,ke) / CVtot(:,ke) ) 
  
      EnTot(:,ke) = &
          Grav * DENS(:) * lcmesh%zlev(:,ke)                                                     &
        + CVtot(:,ke) * PRES(:) / Rtot(:,ke)                                                     &
        + 0.5_RP * ( MOMX(:,ke) * mom_u1(:) + MOMY(:,ke) * mom_u2(:) + MOMZ(:,ke)**2 ) / DENS(:)
    end do
    
    return
  end subroutine cal_DRHOT2Entot

!OCL SERIAL
  subroutine cal_EnTot2DRHOT( DRHOT, &
    DDENS, MOMX, MOMY, MOMZ, EnTot,         &
    PRES_hyd, DENS_hyd, CPtot, CVtot, Rtot, &
    lcmesh, elem3D )

    use scale_const, only: &
      GRAV => CONST_GRAV,  &
      Rdry => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: DRHOT(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: EnTot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA)

    integer :: ke, ke2D

    real(RP) :: DENS(elem3D%Np)
    real(RP) :: mom_u1(elem3D%Np), mom_u2(elem3D%Np)
    real(RP) :: PRES(elem3D%Np)
    !---------------------------------------------------------------

    !$omp parallel do private( ke2D, DENS, mom_u1, mom_u2, PRES )
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      mom_u1(:) = lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,1,1) * MOMX(:,ke) + lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,1) * MOMY(:,ke)
      mom_u2(:) = lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,1) * MOMX(:,ke) + lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,2) * MOMY(:,ke)

      PRES(:) = (  EnTot(:,ke) - Grav * DENS(:) * lcmesh%zlev(:,ke)                                        &
                  - 0.5_RP * ( MOMX(:,ke) * mom_u1(:) + MOMY(:,ke) * mom_u2(:) + MOMZ(:,ke)**2 ) / DENS(:) &
                ) * Rtot(:,ke) / CVtot(:,ke)
      DRHOT(:,ke) = PRES00 / Rtot(:,ke) * ( PRES(:) / PRES00 )**( CVtot(:,ke) / CPtot(:,ke) ) &
                  - PRES00 / Rdry * ( PRES_hyd(:,ke) / PRES00 )**(CvDry/CpDry)
    end do
    
    return
  end subroutine cal_EnTot2DRHOT

  !-- Setup modal filter
!OCL SERIAL
  subroutine setup_modalfilter( this, atm_mesh, read_type )
    implicit none

    class(AtmosDyn), target, intent(inout) :: this
    class(AtmosMesh), target, intent(in) :: atm_mesh
    character(len=*), intent(in) :: read_type

    real(RP) :: MF_ETAC_h  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_h = 36.0_RP
    integer  :: MF_ORDER_h = 16
    real(RP) :: MF_ETAC_v  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_v = 36.0_RP
    integer  :: MF_ORDER_v = 16

    namelist /PARAM_ATMOS_DYN_MODALFILTER/ &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    

    namelist /PARAM_ATMOS_DYN_TRACER_MODALFILTER/ &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    

    integer :: ierr
    character(len=H_SHORT) :: lbl_readtype
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    select case(read_type)
    case ('dyn')
      read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_MODALFILTER,iostat=ierr)
      lbl_readtype = ''
    case ('tracer')
      read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_TRACER_MODALFILTER,iostat=ierr)
      lbl_readtype = 'TRACER_'
    end select
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_setup_modalfilter",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_setup_modalfilter",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_'//trim(lbl_readtype)//'MODALFILTER. Check!'
      call PRC_abort
    endif

    select case(read_type)
    case ('dyn')
      LOG_NML(PARAM_ATMOS_DYN_MODALFILTER)

      if ( .not. associated( this%cal_vi ) ) then
        call atm_mesh%Construct_ModalFilter3D( &
          this%modal_filter_3d,                & ! (inout)
          MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   & ! (in)
          MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    ) ! (in)
      else
        call atm_mesh%Construct_ModalFilterHV( &
          this%modal_filter_3d, this%modal_filter_v1D, & ! (inout)
          MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,           & ! (in)
          MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v            ) ! (in)
      end if  
    case ('tracer')
      LOG_NML(PARAM_ATMOS_DYN_TRACER_MODALFILTER)

      call atm_mesh%Construct_ModalFilter3D( &
        this%modal_filter_tracer_3d,         & ! (inout)
        MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   & ! (in)
        MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    ) ! (in)
    end select

    return
  end subroutine setup_modalfilter

  !-- Setup explicit numerical diffusion
!OCL SERIAL
  subroutine setup_numdiff( this, atm_mesh )
    implicit none

    class(AtmosDyn), target, intent(inout) :: this
    class(AtmosMesh), target, intent(in) :: atm_mesh

    integer ::  ND_LAPLACIAN_NUM = 1
    real(RP) :: ND_COEF_h        = 0.0_RP
    real(RP) :: ND_COEF_v        = 0.0_RP

    namelist /PARAM_ATMOS_DYN_NUMDIFF/ &
      ND_LAPLACIAN_NUM,                &
      ND_COEF_h, ND_COEF_v

    integer :: ierr
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_NUMDIFF,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_setup_numdiff",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_setup_numdiff",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_NUMDIFF. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN_NUMDIFF)      

    this%ND_LAPLACIAN_NUM = ND_LAPLACIAN_NUM
    this%ND_COEF_H = ND_COEF_h
    this%ND_COEF_v = ND_COEF_v
    call atm_dyn_dgm_nonhydro3d_numdiff_Init( atm_mesh%ptr_mesh )

    return
  end subroutine setup_numdiff

  !-- Setup sponge layer
!OCL SERIAL
  subroutine setup_spongelayer( this, atm_mesh, dtsec )
    use scale_mesh_cubedom3d, only: MeshCubeDom3D
    use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
    implicit none

    class(AtmosDyn), target, intent(inout) :: this
    class(AtmosMesh), target, intent(in) :: atm_mesh
    real(RP), intent(in) :: dtsec

    real(RP) :: SL_WDAMP_TAU        = -1.0_RP ! the maximum tau for Rayleigh damping of w [s]
    real(RP) :: SL_WDAMP_HEIGHT     = -1.0_RP ! the height to start apply Rayleigh damping [m]
    integer  :: SL_WDAMP_LAYER      = -1      ! the vertical number of finite element to start apply Rayleigh damping [num]
    logical  :: SL_HORIVELDAMP_FLAG = .false. ! Is the horizontal velocity damped? 
    
    namelist /PARAM_ATMOS_DYN_SPONGELAYER/ &
      SL_WDAMP_TAU,                        &                
      SL_WDAMP_HEIGHT,                     &
      SL_WDAMP_LAYER,                      &
      SL_HORIVELDAMP_FLAG
    
    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D
  
    integer :: NeGZ
    integer :: ierr
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_SPONGELAYER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_setup_spongelayer",*) 'Not found namelist. Default used.'
    else if( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_setup_spongelayer",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_SPONGELAYER. Check!'
      call PRC_abort
    end if
    LOG_NML(PARAM_ATMOS_DYN_SPONGELAYER)

    this%wdamp_tau    = SL_WDAMP_TAU 
    this%wdamp_height = SL_WDAMP_HEIGHT

    mesh3D => atm_mesh%ptr_mesh
    lcmesh3D => mesh3D%lcmesh_list(1)
    elem3D => lcmesh3D%refElem3D

    select type(mesh3D)
    type is (MeshCubeDom3D)
      NeGZ = mesh3D%NeGZ
    type is (MeshCubedSphereDom3D)
      NeGZ = mesh3D%NeGZ
    end select

    if ( SL_WDAMP_LAYER > NeGZ ) then
      LOG_ERROR("ATMOS_DYN_setup_spongelayer",*) 'SL_wdamp_layer should be less than total of vertical elements (NeGZ). Check!'
      call PRC_abort
    else if( SL_WDAMP_LAYER > 0 ) then
      this%wdamp_height = lcmesh3D%pos_en(1,1+(SL_WDAMP_LAYER-1)*lcmesh3D%NeX*lcmesh3D%NeY,3)
    end if
    if ( this%wdamp_tau < 0.0_RP ) then
      this%wdamp_tau = dtsec * 10.0_RP
    else if ( this%wdamp_tau < dtsec ) then
      LOG_ERROR("ATMOS_DYN_setup_spongelayer",*) 'SL_wdamp_tau should be larger than TIME_DT (ATMOS_DYN). Check!'
      call PRC_abort
    end if
    
    this%hvel_damp_flag = SL_HORIVELDAMP_FLAG

    return
  end subroutine setup_spongelayer

  !-- Setup Coriolis parameter

!OCL SERIAL
  subroutine setup_coriolis_parameter( this, atm_mesh )

    use scale_coriolis_param, only: get_coriolis_parameter
    use scale_mesh_cubedom3d, only: MeshCubeDom3D
    implicit none

    class(AtmosDynVars), target, intent(inout) :: this
    class(AtmosMesh), target, intent(in) :: atm_mesh

    class(LocalMeshFieldBase), pointer :: coriolis
    class(LocalMesh3D), pointer :: lcmesh3D
    class(LocalMesh2D), pointer :: lcmesh2D
    integer :: n

    character(len=H_SHORT) :: CORIOLIS_type ! type of coriolis force: 'PLANE', 'SPHERE'
    real(RP) :: CORIOLIS_f0   = 0.0_RP
    real(RP) :: CORIOLIS_beta = 0.0_RP
    real(RP) :: CORIOLIS_y0

    namelist /PARAM_ATMOS_DYN_CORIOLIS/ &
      CORIOLIS_type,                         &                
      CORIOLIS_f0, CORIOLIS_beta, CORIOLIS_y0
        
    class(MeshBase3D), pointer :: mesh3D
    class(MeshCubeDom3D), pointer :: meshCube
    integer :: ierr
    !---------------------------------------------------------------

    mesh3D => atm_mesh%ptr_mesh

    CORIOLIS_type = 'NONE'

    select type(mesh3D)
    type is (MeshCubeDom3D)
      CORIOLIS_y0 = 0.5_RP*(mesh3D%ymax_gl +  mesh3D%ymin_gl)
    end select

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_CORIOLIS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_setup_coriolis",*) 'Not found namelist. Default used.'
    else if( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_setup_coriolis",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_CORIOLIS. Check!'
      call PRC_abort
    end if
    LOG_NML(PARAM_ATMOS_DYN_CORIOLIS)

    do n = 1, mesh3D%LOCAL_MESH_NUM
      call AtmosDynAuxVars_GetLocalMeshFields( n, mesh3D, this%AUXVARS2D_manager, &
        coriolis, lcmesh3D )
      lcmesh2D => lcmesh3D%lcmesh2D

      call get_coriolis_parameter( &
        coriolis%val(:,lcmesh2D%NeS:lcmesh2D%NeE),                       & ! (out)
        CORIOLIS_type, lcmesh2D%refElem2D%Np * lcmesh2D%Ne,              & ! (in)
        lcmesh2D%pos_en(:,:,2), CORIOLIS_f0, CORIOLIS_beta, CORIOLIS_y0, & ! (in)
        lcmesh3D%lat2D(:,:)                                              ) ! (in)
    end do

    return
  end subroutine setup_coriolis_parameter

!--------

!   subroutine cal_MOMZ_tend( &
!     MOMZ_t, MOMZ_t_advx, MOMZ_t_advY, MOMZ_t_advZ, MOMZ_t_lift, MOMZ_t_buoy,     & ! (out)
!      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                    & ! (in)
!      Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )
 
!      use scale_element_base
!      use scale_sparsemat
!      use scale_const, only: &
!       GRAV => CONST_GRAV,  &
!       Rdry => CONST_Rdry,  &
!       CPdry => CONST_CPdry, &
!       CVdry => CONST_CVdry, &
!       PRES00 => CONST_PRE00
!      use scale_atm_dyn_nonhydro3d, only: IntrpMat_VPOrdM1
!      implicit none
 
!      class(LocalMesh3D), intent(in) :: lmesh
!      class(elementbase3D), intent(in) :: elem
!      class(LocalMesh2D), intent(in) :: lmesh2D
!      class(elementbase2D), intent(in) :: elem2D
!      type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift
!      real(RP), intent(out) :: MOMZ_t(elem%Np,lmesh%NeA)
!      real(RP), intent(out) :: MOMZ_t_advx(elem%Np,lmesh%NeA)
!      real(RP), intent(out) :: MOMZ_t_advy(elem%Np,lmesh%NeA)
!      real(RP), intent(out) :: MOMZ_t_advz(elem%Np,lmesh%NeA)
!      real(RP), intent(out) :: MOMZ_t_lift(elem%Np,lmesh%NeA)
!      real(RP), intent(out) :: MOMZ_t_buoy(elem%Np,lmesh%NeA)

!      real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
 
!      real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
!      real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
!      real(RP) :: dens_(elem%Np), RHOT_(elem%Np), dpres_(elem%Np)
!      real(RP) :: pres_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np)
 
!      integer :: ke
!      !------------------------------------------------------------------------
 
!      call cal_del_flux_dyn( del_flux,                                          & ! (out)
!        DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                & ! (in)
!        lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
!        lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
!        lmesh, elem )                                                             ! (in)
  
!      !-----
!      !$omp parallel do private(RHOT_,pres_,dpres_,dens_,u_,v_,w_,Fx,Fy,Fz,LiftDelFlx)
!      do ke = lmesh%NeS, lmesh%NeE
!        !--
 
!        RHOT_(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
!        pres_(:) = PRES00 * (Rdry*RHOT_(:)/PRES00)**(CPdry/Cvdry)
!        dpres_(:) = pres_(:) - PRES_hyd(:,ke)
!        dens_(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
 
!        u_(:) = MOMX_(:,ke)/dens_(:)
!        v_(:) = MOMY_(:,ke)/dens_(:)
!        w_(:) = MOMZ_(:,ke)/dens_(:)
 
!        !-- MOMZ
!        call sparsemat_matmul(Dx, u_(:)*MOMZ_(:,ke), Fx)
!        call sparsemat_matmul(Dy, v_(:)*MOMZ_(:,ke), Fy)
!        call sparsemat_matmul(Dz, w_(:)*MOMZ_(:,ke), Fz)
!        MOMZ_t_advx(:,ke) = - lmesh%Escale(:,ke,1,1) * Fx(:)
!        MOMZ_t_advy(:,ke) = - lmesh%Escale(:,ke,2,2) * Fy(:)
!        MOMZ_t_advz(:,ke) = - lmesh%Escale(:,ke,3,3) * Fz(:)

!        call sparsemat_matmul(Dz, dpres_(:), Fz)
!        MOMZ_t_buoy(:,ke) = - lmesh%Escale(:,ke,3,3) * Fz(:) &
!                            - matmul(IntrpMat_VPOrdM1, DDENS_(:,ke)) * Grav

!        call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke), LiftDelFlx)
!        MOMZ_t_lift(:,ke) = - LiftDelFlx(:)

!        MOMZ_t(:,ke) = MOMZ_t_advx(:,ke) +  MOMZ_t_advy(:,ke) +  MOMZ_t_advz(:,ke) &
!                     + MOMZ_t_lift(:,ke) + MOMZ_t_buoy(:,ke)
!      end do
 
!      return
!  end subroutine cal_MOMZ_tend


!  subroutine cal_del_flux_dyn( del_flux, &
!    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,   &
!    nx, ny, nz, vmapM, vmapP, lmesh, elem )

!    use scale_const, only: &
!     GRAV => CONST_GRAV,  &
!     Rdry => CONST_Rdry,  &
!     CPdry => CONST_CPdry, &
!     CVdry => CONST_CVdry, &
!     PRES00 => CONST_PRE00
  
!    implicit none

!    class(LocalMesh3D), intent(in) :: lmesh
!    class(elementbase3D), intent(in) :: elem  
!    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
!    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
!    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
!    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
!    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
!    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
!    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
!    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
!    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
!    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
!    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
!    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
!    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
   
!    integer :: i, iP, iM
!    real(RP) :: VelP, VelM, alpha
!    real(RP) :: uM, uP, vM, vP, wM, wP, presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P
!    real(RP) :: gamm, rgamm
!    !------------------------------------------------------------------------

!    gamm = CpDry/CvDry
!    rgamm = CvDry/CpDry

!    !$omp parallel do private( &
!    !$omp iM, iP, uM, VelP, VelM, alpha, &
!    !$omp uP, vM, vP, wM, wP, presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P)
!    do i=1, elem%NfpTot*lmesh%Ne
!      iM = vmapM(i); iP = vmapP(i)

!      rhot_hyd_M = PRES00/Rdry * (PRES_hyd(iM)/PRES00)**rgamm
!      rhot_hyd_P = PRES00/Rdry * (PRES_hyd(iP)/PRES00)**rgamm
     
!      rhotM = rhot_hyd_M + DRHOT_(iM)
!      presM = PRES00 * (Rdry*rhotM/PRES00)**gamm
!      dpresM = presM - PRES_hyd(iM)*abs(nz(i))

!      rhotP = rhot_hyd_P + DRHOT_(iP) 
!      presP = PRES00 * (Rdry*rhotP/PRES00)**gamm
!      dpresP = presP - PRES_hyd(iP)*abs(nz(i))

!      densM = DDENS_(iM) + DENS_hyd(iM)
!      densP = DDENS_(iP) + DENS_hyd(iP)

!      VelM = (MOMX_(iM)*nx(i) + MOMY_(iM)*ny(i) + MOMZ_(iM)*nz(i))/densM
!      VelP = (MOMX_(iP)*nx(i) + MOMY_(iP)*ny(i) + MOMZ_(iP)*nz(i))/densP

!      alpha = max( sqrt(gamm*presM/densM) + abs(VelM), sqrt(gamm*presP/densP) + abs(VelP)  )


!      del_flux(i) = 0.5_RP*(                &
!                    ( MOMZ_(iP)*VelP - MOMZ_(iM)*VelM)   &
!                    + ( dpresP - dpresM )*nz(i)          &                    
!                    - alpha*(MOMZ_(iP) - MOMZ_(iM))      )
!    end do

!    return
!  end subroutine cal_del_flux_dyn

end module mod_atmos_dyn