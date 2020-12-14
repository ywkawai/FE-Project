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
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_component_proc, only:  ModelComponentProc

  use scale_atm_dyn_nonhydro3d_heve, only: &
    atm_dyn_nonhydro3d_heve_Init,              &
    atm_dyn_nonhydro3d_heve_Final,             &
    atm_dyn_nonhydro3d_heve_cal_tend

  use scale_atm_dyn_nonhydro3d_hevi, only: &
    atm_dyn_nonhydro3d_hevi_Init,              &
    atm_dyn_nonhydro3d_hevi_Final,             &
    atm_dyn_nonhydro3d_hevi_cal_tend,          &
    atm_dyn_nonhydro3d_hevi_cal_vi

  use scale_atm_dyn_nonhydro3d_splitform_hevi, only: &
    atm_dyn_nonhydro3d_hevi_splitform_Init,          &
    atm_dyn_nonhydro3d_hevi_splitform_Final,         &
    atm_dyn_nonhydro3d_hevi_splitform_cal_tend,      &
    atm_dyn_nonhydro3d_hevi_splitform_cal_vi    
  
  use scale_atm_dyn_nonhydro3d_numdiff, only: &
    atm_dyn_nonhydro3d_numdiff_Init,          &
    atm_dyn_nonhydro3d_numdiff_Final
  
  use scale_element_modalfilter, only: ModalFilter

  use mod_atmos_mesh, only: AtmosMesh
  use mod_atmos_vars, only: &
    AtmosVars_GetLocalMeshPrgVar,        &
    AtmosVars_GetLocalMeshPrgVars,       &
    ATMOS_PROGVARS_NUM,                  &
    DDENS_ID => ATMOS_PROGVARS_DDENS_ID, &
    DRHOT_ID => ATMOS_PROGVARS_DRHOT_ID, &
    MOMX_ID  => ATMOS_PROGVARS_MOMX_ID,  &
    MOMY_ID  => ATMOS_PROGVARS_MOMY_ID,  &
    MOMZ_ID  => ATMOS_PROGVARS_MOMZ_ID
  use mod_atmos_dyn_bnd, only: AtmosDynBnd
  use mod_atmos_dyn_vars, only: &
    AtmosDynVars,                      &
    AtmosDynAuxVars_GetLocalMeshFields

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
      SL_flag, wdamp_tau, wdamp_height,                                           & ! (in)
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
      logical, intent(in)   :: SL_flag
      real(RP), intent(in)  :: wdamp_tau
      real(RP), intent(in)  :: wdamp_height
  
    end subroutine atm_dyn_nonhydro3d_cal_tend_ex
  end interface

  abstract interface    
    subroutine atm_dyn_nonhydro3d_cal_vi( &
      DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,             & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, & ! (in)
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
    type(ModalFilter) :: modal_filter_v1D

    ! sponge layer
    logical :: SPONGELAYER_FLAG
    real(RP) :: wdamp_tau
    real(RP) :: wdamp_height

  contains
    procedure, public :: setup => AtmosDyn_setup 
    procedure, public :: calc_tendency => AtmosDyn_calc_tendency
    procedure, public :: update => AtmosDyn_update
    procedure, public :: finalize => AtmosDyn_finalize
  end type AtmosDyn

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------
  
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_HEVE           = 1
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_HEVI           = 2
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_SPLITFORM_HEVI = 3


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

  subroutine AtmosDyn_setup( this, model_mesh, tm_parent_comp )
    use mod_atmos_mesh, only: AtmosMesh
    use mod_atmos_vars, only: ATMOS_PROGVARS_NUM
    use scale_time_manager, only: TIME_manager_component

    implicit none

    class(AtmosDyn), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    character(len=H_MID) :: EQS_TYPE    = "NONHYDRO3D_HEVE"
    character(len=H_SHORT) :: TINTEG_TYPE = 'RK_TVD_3'
    real(DP) :: TIME_DT                             = UNDEF8
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  
    
    logical  :: MODALFILTER_FLAG = .false.
    logical :: NUMDIFF_FLAG      = .false.
    logical :: SPONGELAYER_FLAG  = .false.

    character(len=H_SHORT) :: coriolis_type = 'PLANE'   ! type of coriolis force: 'PLANE', 'SPHERE'
    real(RP) :: coriolis_f0         = 0.0_RP
    real(RP) :: coriolis_beta       = 0.0_RP
    real(RP) :: coriolis_y0         = UNDEF8            ! default is domain center    

    namelist / PARAM_ATMOS_DYN /       &
      EQS_TYPE,                               &
      TINTEG_TYPE,                            &
      TIME_DT,                                &
      TIME_DT_UNIT,                           &
      MODALFILTER_FLAG,                       &
      NUMDIFF_FLAG,                           &
      SPONGELAYER_FLAG,                       &
      CORIOLIS_TYPE,                          &
      CORIOLIS_f0, CORIOLIS_beta, CORIOLIS_y0
    
    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh
    class(LocalMeshBase), pointer :: ptr_lcmesh
    class(ElementBase3D), pointer :: elem3D
    integer :: n
    real(DP) :: dtsec

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
    type is (AtmosMesh)
      atm_mesh => model_mesh
    end select

    !- Setup the temporal integrator

    call tm_parent_comp%Regist_process( 'ATMOS_DYN', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                    ! (out)

    dtsec = tm_parent_comp%process_list(this%tm_process_id)%dtsec
    
    allocate( this%tint(ptr_mesh%LOCAL_MESH_NUM) )
    do n = 1, ptr_mesh%LOCAL_MESH_NUM
      call ptr_mesh%GetLocalMesh( n, ptr_lcmesh )
      call this%tint(n)%Init( TINTEG_TYPE, dtsec, ATMOS_PROGVARS_NUM, 2, &
        (/ ptr_mesh%refElem%Np, ptr_lcmesh%NeA /) )
    end do

    !- initialize an object to manage boundary conditions
    call this%boundary_cond%Init()
    call this%boundary_cond%SetBCInfo( ptr_mesh )

    !- initialize the variables 
    call this%dyn_vars%Init( model_mesh )
    call setup_coriolis_parameter( this%dyn_vars, atm_mesh, CORIOLIS_type, CORIOLIS_f0, CORIOLIS_beta, CORIOLIS_y0 )

    !- Initialize a module for 3D dynamical core 

    select case(EQS_TYPE)
    case("NONHYDRO3D_HEVE")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_HEVE
      call atm_dyn_nonhydro3d_heve_Init( atm_mesh%mesh )
      this%cal_tend_ex => atm_dyn_nonhydro3d_heve_cal_tend
      this%cal_vi => null()
    case("NONHYDRO3D_HEVI")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_HEVI
      call atm_dyn_nonhydro3d_hevi_Init( atm_mesh%mesh )
      this%cal_tend_ex => atm_dyn_nonhydro3d_hevi_cal_tend
      this%cal_vi => atm_dyn_nonhydro3d_hevi_cal_vi
    case("NONHYDRO3D_SPLITFORM_HEVI")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_SPLITFORM_HEVI
      call atm_dyn_nonhydro3d_hevi_splitform_Init( atm_mesh%mesh )
      this%cal_tend_ex => atm_dyn_nonhydro3d_hevi_splitform_cal_tend
      this%cal_vi => atm_dyn_nonhydro3d_hevi_splitform_cal_vi
    case default
      LOG_ERROR("ATMOS_DYN_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
      call PRC_abort
    end select    

    !- Setup the numerical diffusion
    this%CALC_NUMDIFF_FLAG = NUMDIFF_FLAG
    if( NUMDIFF_FLAG ) call setup_numdiff( this, atm_mesh )

    !- Setup the modal filter
    this%MODALFILTER_FLAG = MODALFILTER_FLAG
    if ( MODALFILTER_FLAG ) call setup_modalfilter( this, atm_mesh )

    !- Setup the sponge layer
    this%SPONGELAYER_FLAG = SPONGELAYER_FLAG
    if ( SPONGELAYER_FLAG ) call setup_spongelayer( this, atm_mesh, dtsec )

    return
  end subroutine AtmosDyn_setup

  subroutine AtmosDyn_calc_tendency( this, model_mesh, prgvars_list, auxvars_list, forcing_list, is_update )
    implicit none
    
    class(AtmosDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------
    if (.not. this%IsActivated()) return
    !LOG_INFO('AtmosDyn_tendency',*)

    return  
  end subroutine AtmosDyn_calc_tendency

  subroutine AtmosDyn_update( this, model_mesh, prgvars_list, auxvars_list, forcing_list, is_update )

    use scale_atm_dyn_modalfilter, only: &
      atm_dyn_modalfilter_apply

    implicit none

    class(AtmosDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    integer :: rkstage
    integer :: tintbuf_ind

    class(MeshBase), pointer :: mesh
    class(MeshBase2D), pointer :: mesh2D    
    class(LocalMesh3D), pointer :: lcmesh
    integer :: n
    integer :: ke

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: Coriolis

!    class(LocalMeshFieldBase), pointer :: MOMZ_t, MOMZ_t_advx, MOMZ_t_advY, MOMZ_t_advZ, MOMZ_t_lift, MOMZ_t_buoy
    integer :: v
    real(RP) :: implicit_fac
    real(RP) :: dt
    character(len=H_SHORT) :: labl
    !--------------------------------------------------
    
    call PROF_rapstart( 'ATM_DYN_update', 1)   

    call model_mesh%GetModelMesh( mesh )

    !-
    do rkstage=1, this%tint(1)%nstage

      if (this%tint(1)%imex_flag) then        
        do n=1, mesh%LOCAL_MESH_NUM
          call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
          call AtmosVars_GetLocalMeshPrgVars( n, &
            mesh, prgvars_list, auxvars_list,                               &
            DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
            DENS_hyd, PRES_hyd, lcmesh                                      )
          call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)   

          call PROF_rapstart( 'ATM_DYN_cal_vi', 2)
          implicit_fac = this%tint(n)%Get_implicit_diagfac(rkstage)
          tintbuf_ind = this%tint(n)%tend_buf_indmap(rkstage)
          dt = this%tint(n)%Get_deltime()
          call this%cal_vi( &
            this%tint(n)%tend_buf2D_im(:,:,DDENS_ID,tintbuf_ind),    & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,MOMX_ID ,tintbuf_ind),    & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,MOMY_ID ,tintbuf_ind),    & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,MOMZ_ID ,tintbuf_ind),    & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,DRHOT_ID,tintbuf_ind),    & ! (out)
            DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                     & ! (in)
            DENS_hyd%val, PRES_hyd%val,                                             & ! (in)
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

          call this%tint(n)%StoreImplicit( rkstage, DRHOT%val, DRHOT_ID,  &
                             1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
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
          mesh, prgvars_list, auxvars_list,                               &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
          DENS_hyd, PRES_hyd, lcmesh                                      )
        call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)         
        
        !* Apply boundary conditions
        call PROF_rapstart( 'ATM_DYN_applyBC_prgv', 2)
        call this%boundary_cond%ApplyBC_PROGVARS_lc( n,                              & ! (in)
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                        & ! (inout)
          DENS_hyd%val, PRES_hyd%val,                                                & ! (in)
          lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), & ! (in)
          lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D         ) ! (in)
        call PROF_rapend( 'ATM_DYN_applyBC_prgv', 2)
      end do


      do n=1, mesh%LOCAL_MESH_NUM
        tintbuf_ind = this%tint(n)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
        call AtmosVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list,    &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,      &
          DENS_hyd, PRES_hyd, lcmesh           )
        
        call AtmosDynAuxVars_GetLocalMeshFields( n,      &
          mesh, this%dyn_vars%AUXVARS2D_manager,         &
          Coriolis )
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
          this%SPONGELAYER_FLAG, this%wdamp_tau, this%wdamp_height,               &
          model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3), &
          model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3), &
          model_mesh%LiftOptrMat,                                                 &
          lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D ) 
        call PROF_rapend( 'ATM_DYN_update_caltend_ex', 2)

        call PROF_rapstart( 'ATM_DYN_update_add_tp', 2)
        call add_phy_tend( &
          this, this%tint(n)%tend_buf2D_ex(:,:,:,tintbuf_ind), & ! (inout)
          DRHOT%val, PRES_hyd%val, forcing_list,               & ! (in)
          mesh, n, lcmesh, lcmesh%refElem3D                    ) ! (in)
        call PROF_rapend( 'ATM_DYN_update_add_tp', 2)

        call PROF_rapstart( 'ATM_DYN_update_advance', 2)      
        call this%tint(n)%Advance( rkstage, DDENS%val, DDENS_ID, &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        call this%tint(n)%Advance( rkstage, MOMX%val, MOMX_ID,   &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        call this%tint(n)%Advance( rkstage, MOMY%val, MOMY_ID,   &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )

        call this%tint(n)%Advance( rkstage, MOMZ%val, MOMZ_ID,   &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )

        call this%tint(n)%Advance( rkstage, DRHOT%val, DRHOT_ID, &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call PROF_rapend( 'ATM_DYN_update_advance', 2)

        !------------------------------------------------------------------------------
      end do
    end do

    !-- modal filter
    if ( this%MODALFILTER_FLAG ) then
      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
        call AtmosVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list,                               &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
          DENS_hyd, PRES_hyd, lcmesh                                      )
        call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)

        call PROF_rapstart( 'ATM_DYN_update_expfilter', 2)
        call atm_dyn_modalfilter_apply(                       & ! (inout)
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val, & ! (in)
          lcmesh, lcmesh%refElem3D, this%modal_filter_3d      ) ! (in)
        call PROF_rapend( 'ATM_DYN_update_expfilter', 2)
      end do
    end if

    !-- numerical diffusion
    if ( this%CALC_NUMDIFF_FLAG ) then
      
      call PROF_rapstart( 'ATM_DYN_numfilter', 1)
      call prgvars_list%MeshFieldComm_Exchange()

      do v = 1, ATMOS_PROGVARS_NUM
        call cal_numfilter_tend( this, model_mesh, prgvars_list, auxvars_list, v )
      end do

      do n=1, mesh%LOCAL_MESH_NUM
        dt = this%tint(n)%Get_deltime()

        call AtmosVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list,    &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,      &
          DENS_hyd, PRES_hyd, lcmesh           )
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

    call PROF_rapend( 'ATM_DYN_update', 1)   


    return  
  end subroutine AtmosDyn_update

  subroutine AtmosDyn_finalize( this )
    implicit none
    class(AtmosDyn), intent(inout) :: this

    integer :: n
    !--------------------------------------------------
    if (.not. this%IsActivated()) return
    LOG_INFO('AtmosDyn_finalize',*)

    select case(this%EQS_TYPEID)
    case(EQS_TYPEID_NONHYD3D_HEVE)
      call atm_dyn_nonhydro3d_heve_Final()
    case(EQS_TYPEID_NONHYD3D_HEVI)  
      call atm_dyn_nonhydro3d_hevi_Final()
    case(EQS_TYPEID_NONHYD3D_SPLITFORM_HEVI)  
      call atm_dyn_nonhydro3d_hevi_splitform_Final()     
    end select 

    if (this%CALC_NUMDIFF_FLAG) then
      call atm_dyn_nonhydro3d_numdiff_Final()
    end if

    if (this%MODALFILTER_FLAG) then
      call this%modal_filter_3d%Final()
      if ( associated(this%cal_vi) ) call this%modal_filter_v1D%Final()
    end if
    
    do n = 1, size(this%tint)
      call this%tint(n)%Final()
    end do
    deallocate( this%tint )

    call this%boundary_cond%Final()
    call this%dyn_vars%Final()

    return  
  end subroutine AtmosDyn_finalize  

  !--- private ---------------

  subroutine add_phy_tend( this,      & ! (in)
    dyn_tends,                        & ! (inout)
    DRHOT, PRES_hyd,                  & ! (in)
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
    class(ModelVarManager), intent(inout) :: phytends_list
    class(MeshBase), intent(in) :: mesh
    integer, intent(in) :: domID

    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOH_p
    integer :: ke

    real(RP) :: RHOT(elem3D%Np)
    real(RP) :: EXNER(elem3D%Np)
    !---------------------------------------------------------------------------------

    call AtmosVars_GetLocalMeshPhyTends( domID, mesh, phytends_list, & ! (in)
      DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOH_p                     ) ! (out)

    !$omp parallel do          &
    !$Omp private( RHOT, EXNER )
    do ke=lcmesh%NeS, lcmesh%NeE
      RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT(:,ke)
      EXNER(:) = (Rdry*RHOT(:)/PRES00)**(Rdry/Cvdry)

      dyn_tends(:,ke,DDENS_ID) = dyn_tends(:,ke,DDENS_ID) + DENS_tp%val(:,ke)
      dyn_tends(:,ke,MOMX_ID ) = dyn_tends(:,ke,MOMX_ID ) + MOMX_tp%val(:,ke)
      dyn_tends(:,ke,MOMY_ID ) = dyn_tends(:,ke,MOMY_ID ) + MOMY_tp%val(:,ke)
      dyn_tends(:,ke,MOMZ_ID ) = dyn_tends(:,ke,MOMZ_ID ) + MOMZ_tp%val(:,ke)
      dyn_tends(:,ke,DRHOT_ID) = dyn_tends(:,ke,DRHOT_ID) + RHOH_p %val(:,ke) / ( CpDry * EXNER(:) )
    end do

    return
  end subroutine add_phy_tend

  subroutine cal_numfilter_tend( this, model_mesh, prgvars_list, auxvars_list, varid )

    use mod_atmos_dyn_vars, only: &
      AtmosDynAuxVars_GetLocalMeshFields,     &
      AtmosDynNumDiffFlux_GetLocalMeshFields, &
      AtmosDynNumDiffTend_GetLocalMeshFields
    
    use scale_atm_dyn_nonhydro3d_numdiff, only: &
      atm_dyn_nonhydro3d_numdiff_tend,          &
      atm_dyn_nonhydro3d_numdiff_cal_laplacian, &
      atm_dyn_nonhydro3d_numdiff_cal_flx

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

    class(MeshBase), pointer :: mesh
    class(LocalMesh3D), pointer :: lcmesh
    integer :: n
    integer :: ke
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
        DENS_hyd, PRES_hyd                                                     )
      call AtmosDynNumDiffFlux_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_FLUX_manager, &
        ND_flx_x, ND_flx_y, ND_flx_z )
      
      allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
      call this%boundary_cond%ApplyBC_numdiff_even_lc( var%val, is_bound, varid, n, &
        MOMX%val, MOMY%val, MOMZ%val, DENS_hyd%val, PRES_hyd%val,                    &
        lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3),   &
        lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )
      
      call atm_dyn_nonhydro3d_numdiff_cal_flx( ND_flx_x%val, ND_flx_y%val, ND_flx_z%val, &
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

        call atm_dyn_nonhydro3d_numdiff_cal_laplacian( ND_lapla_h%val, ND_lapla_v%val, &
          ND_flx_x%val, ND_flx_y%val, ND_flx_z%val,                                    &
          model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),      &
          model_mesh%LiftOptrMat,                                                      &
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
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                                       &
          DENS_hyd, PRES_hyd                                                    )
          
        allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
        call this%boundary_cond%ApplyBC_numdiff_even_lc( &
          ND_lapla_h%val, is_bound, varid, n,                                        &
          MOMX%val, MOMY%val, MOMZ%val, DENS_hyd%val, PRES_hyd%val,                  &
          lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), &
          lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )
          
        call atm_dyn_nonhydro3d_numdiff_cal_flx( ND_flx_x%val, ND_flx_y%val, ND_flx_z%val, &
          ND_lapla_h%val, ND_lapla_v%val, DDENS%val, DENS_hyd%val,                         &
          model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),          &
          model_mesh%LiftOptrMat,                                                          &
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

      call atm_dyn_nonhydro3d_numdiff_tend( this%tint(n)%tend_buf2D_ex(:,:,varid,1),  &
        ND_flx_x%val, ND_flx_y%val, ND_flx_z%val,                                     &
        DDENS%val, DENS_hyd%val, nd_sign * this%ND_COEF_H, nd_sign * this%ND_COEF_V,  &
        model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),       &
        model_mesh%LiftOptrMat,                                                       &
        lcmesh, lcmesh%refElem3D, is_bound, dens_weight_flag ) 

      deallocate( is_bound )
    end do

    return
  end subroutine cal_numfilter_tend

  !-- Setup modal filter
  subroutine setup_modalfilter( this, atm_mesh )
    implicit none

    class(AtmosDyn), target, intent(inout) :: this
    class(AtmosMesh), target, intent(in) :: atm_mesh

    real(RP) :: MF_ETAC_h  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_h = 36.0_RP
    integer  :: MF_ORDER_h = 16
    real(RP) :: MF_ETAC_v  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_v = 36.0_RP
    integer  :: MF_ORDER_v = 16

    namelist /PARAM_ATMOS_DYN_MODALFILTER/ &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    

    integer :: ierr
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_MODALFILTER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_setup_modalfilter",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_setup_modalfilter",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_MODALFILTER. Check!'
      call PRC_abort
    endif
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

    return
  end subroutine setup_modalfilter

  !-- Setup explicit numerical diffusion

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
    call atm_dyn_nonhydro3d_numdiff_Init( atm_mesh%mesh )

    return
  end subroutine setup_numdiff

  !-- Setup sponge layer
  subroutine setup_spongelayer( this, atm_mesh, dtsec )
    implicit none

    class(AtmosDyn), target, intent(inout) :: this
    class(AtmosMesh), target, intent(in) :: atm_mesh
    real(RP), intent(in) :: dtsec

    real(RP) :: SL_WDAMP_TAU    = -1.0_RP ! the maximum tau for Rayleigh damping of w [s]
    real(RP) :: SL_WDAMP_HEIGHT = -1.0_RP ! the height to start apply Rayleigh damping [m]
    integer  :: SL_WDAMP_LAYER  = -1      ! the vertical number of finite element to start apply Rayleigh damping [num]
    
    namelist /PARAM_ATMOS_DYN_SPONGELAYER/ &
      SL_WDAMP_TAU,                        &                
      SL_WDAMP_HEIGHT,                     &
      SL_WDAMP_LAYER
    
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D
  
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

    lcmesh3D => atm_mesh%mesh%lcmesh_list(1)
    elem3D => lcmesh3D%refElem3D

    if ( SL_WDAMP_LAYER > atm_mesh%mesh%NeGZ ) then
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
    
    return
  end subroutine setup_spongelayer

  !-- Setup Coriolis parameter

  subroutine setup_coriolis_parameter( this, atm_mesh, &
    COLIORIS_type, f0, beta, y0_ )

    use scale_const, only: &
      OHM => CONST_OHM
    implicit none

    class(AtmosDynVars), target, intent(inout) :: this
    class(AtmosMesh), target, intent(in) :: atm_mesh
    character(*), intent(in) :: COLIORIS_type
    real(RP), intent(in) :: f0, beta, y0_

    class(LocalMeshFieldBase), pointer :: coriolis
    class(LocalMesh3D), pointer :: lcmesh3D
    class(LocalMesh2D), pointer :: lcmesh2D
    integer :: n, ke
    real(RP) :: y0
    !---------------------------------------------------------------

    if (y0_ == UNDEF8) then
      y0 = 0.5_RP*(atm_mesh%mesh%ymax_gl +  atm_mesh%mesh%ymin_gl)
    else
      y0 = y0_
    end if

    do n = 1, atm_mesh%mesh%LOCAL_MESH_NUM
      call AtmosDynAuxVars_GetLocalMeshFields( n, atm_mesh%mesh, this%AUXVARS2D_manager, &
        coriolis, lcmesh3D )
      lcmesh2D => lcmesh3D%lcmesh2D

      if ( trim(COLIORIS_type) == 'PLANE' ) then
        do ke=1, lcmesh2D%Ne
          coriolis%val(:,ke) = f0 + beta * (lcmesh2D%pos_en(:,ke,2) - y0)
        end do
      else if ( trim(COLIORIS_type) == 'SPHERE' ) then
        do ke=1, lcmesh2D%Ne
          coriolis%val(:,ke) = 2.0_RP * OHM * sin(lcmesh3D%lat2D(:,ke))
        end do
      else
        LOG_ERROR('AtmosDyn_set_colioris_parameter',*) 'Unexpected COLIORIS_type is specified. Check! COLIORIS_type=', COLIORIS_type
        call PRC_abort
      end if  
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