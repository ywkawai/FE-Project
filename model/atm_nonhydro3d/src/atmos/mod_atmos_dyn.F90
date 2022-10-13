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

  use scale_atm_dyn_dgm_driver_nonhydro3d, only: AtmDynDGMDriver_nonhydro3d
  use scale_atm_dyn_dgm_nonhydro3d_numdiff, only: AtmDyn_Nonhydro3D_Numdiff
  
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
    AtmosVars_GetLocalMeshQTRCPhyTend
  use mod_atmos_dyn_vars, only: &
    AtmosDynVars,                                &
    AtmosDynAuxVars_GetLocalMeshFields,          &
    AtmosDynMassFlux_GetLocalMeshFields,         &
    TRCQ_ID    => ATMOS_DYN_TRCVARS3D_TRCADV_ID    

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  type, extends(ModelComponentProc), public :: AtmosDyn
    integer :: EQS_TYPEID
    type(TimeInt_RK), allocatable :: tint_qtrc(:)

    type(AtmDynDGMDriver_nonhydro3d) :: dyncore_driver
    type(AtmosDynVars) :: dyn_vars

    ! explicit numerical diffusion
    logical :: CALC_NUMDIFF_FLAG
    type(AtmDyn_Nonhydro3D_Numdiff) :: numdiff

    ! element-wise modal filter
    logical :: MODALFILTER_FLAG
    type(ModalFilter) :: modal_filter_tracer_3d

    ! tracer advection
    logical :: ONLY_TRACERADV_FLAG
    logical :: TRACERADV_disable_limiter
    logical :: TRACERADV_MODALFILTER_FLAG
    type(SparseMat) :: FaceIntMat    

  contains
    procedure, public :: setup => AtmosDyn_setup 
    procedure, public :: calc_tendency => AtmosDyn_calc_tendency
    procedure, public :: update => AtmosDyn_update
    procedure, public :: finalize => AtmosDyn_finalize
  end type AtmosDyn

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: setup_modalfilter
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
    
    allocate( this%tint_qtrc(ptr_mesh%LOCAL_MESH_NUM) )

    do n = 1, ptr_mesh%LOCAL_MESH_NUM
      call ptr_mesh%GetLocalMesh( n, ptr_lcmesh )

      call this%tint_qtrc(n)%Init( TINTEG_TYPE_TRACER, dtsec, 1, 2, &
        (/ ptr_mesh%refElem%Np, ptr_lcmesh%NeA /) )        
    end do


    !- initialize the variables 
    call this%dyn_vars%Init( model_mesh )

    call setup_coriolis_parameter( this%dyn_vars, atm_mesh )

    !- Initialize a module for 3D dynamical core 
    call this%dyncore_driver%Init( EQS_TYPE, &
      TINTEG_TYPE, dtsec,                        &
      SPONGELAYER_FLAG, MODALFILTER_FLAG, mesh3D )

    !- Initialize a module for tracer equations
    call atm_dyn_dgm_trcadvect3d_heve_Init( mesh3D, this%FaceIntMat ) 

    !- Setup the numerical diffusion
    this%CALC_NUMDIFF_FLAG = NUMDIFF_FLAG
    if (this%CALC_NUMDIFF_FLAG) call this%numdiff%Init( atm_mesh, dtsec )

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

    use scale_meshfield_base, only: &
      MeshField2D, MeshField3D
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRGVAR_NUM, &
      PRGVAR_DDENS_ID

    use mod_atmos_dyn_vars, only: &
      ATMOS_DYN_AUXVARS2D_CORIOLIS_ID,             &
      TRCDDENS_ID => ATMOS_DYN_TRCVARS3D_DENS_ID,  &
      TRCDDENS0_ID => ATMOS_DYN_TRCVARS3D_DENS_ID, &
      MASSFLX_X_TAVG => ATMOS_DYN_MASSFLX_X_ID,    &
      MASSFLX_Y_TAVG => ATMOS_DYN_MASSFLX_Y_ID,    &
      MASSFLX_Z_TAVG => ATMOS_DYN_MASSFLX_Z_ID
      
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

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, THERM
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot

    class(LocalMeshFieldBase), pointer :: ALPH_DENS_M_tavg, ALPH_DENS_P_tavg, MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg    
    class(LocalMeshFieldBase), pointer :: QTRC
    class(LocalMeshFieldBase), pointer :: RHOQ_tp

    integer :: v
    integer :: iq
    real(RP) :: implicit_fac
    real(RP) :: dt
    real(RP) :: dttmp_trc

    class(MeshBase3D), pointer :: mesh3D
    !--------------------------------------------------
    
    call PROF_rapstart( 'ATM_DYN_update', 1)   

    call model_mesh%GetModelMesh( mesh )
    select type(mesh)
    class is (MeshBase3D)
      mesh3D => mesh
    end select

    if ( .not. this%ONLY_TRACERADV_FLAG ) then

      call this%dyncore_driver%Update( &
        prgvars_list, auxvars_list, forcing_list,                                  & ! (inout)
        this%dyn_vars%AUX_TRCVARS3D(TRCDDENS_ID),                                  & ! (inout)
        this%dyn_vars%AUX_TRCVARS3D(TRCDDENS0_ID),                                 & ! (inout)
        this%dyn_vars%AUXTRC_FLUX_VARS3D(MASSFLX_X_TAVG),                          & ! (inout)
        this%dyn_vars%AUXTRC_FLUX_VARS3D(MASSFLX_Y_TAVG),                          & ! (inout)
        this%dyn_vars%AUXTRC_FLUX_VARS3D(MASSFLX_Z_TAVG),                          & ! (inout)
        this%dyn_vars%alphaDensM, this%dyn_vars%alphaDensP,                        & ! (inout)
        this%dyn_vars%AUX_VARS2D(ATMOS_DYN_AUXVARS2D_CORIOLIS_ID),                 & ! (in)
        model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),    & ! (in)
        model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3),    & ! (in)
        model_mesh%LiftOptrMat, mesh3D                                             ) ! (in)
            
    end if

    !-- Tracer advection (prepair) ------------------------------------------------

    if ( QA > 0 ) then
      call PROF_rapstart( 'ATM_DYN_qtracer', 2)

      if ( this%ONLY_TRACERADV_FLAG ) then
        do n=1, mesh%LOCAL_MESH_NUM
          call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 3) 
          call AtmosVars_GetLocalMeshPrgVars( n, &
            mesh, prgvars_list, auxvars_list,       &
            DDENS, MOMX, MOMY, MOMZ, THERM,         &
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
            this%dyn_vars%AUX_TRCVARS3D(TRCDDENS0_ID)%local(n)%val(:,ke) = DDENS%val(:,ke) 
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
          DDENS, MOMX, MOMY, MOMZ, THERM,         &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
          lcmesh                                  )
                   
        call AtmosDynMassFlux_GetLocalMeshFields( n, mesh, this%dyn_vars,            &
          ALPH_DENS_M_tavg, ALPH_DENS_P_tavg, MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg  )
        call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 3)

        call PROF_rapstart( 'ATM_DYN_applyBC_mflux', 3)
        call this%dyncore_driver%boundary_cond%ApplyBC_PROGVARS_lc( n,                 & ! (in)
          this%dyn_vars%TRCVARS3D(TRCDDENS_ID)%local(n)%val(:,:),                      & ! (inout)
          MFLX_x_tavg%val, MFLX_y_tavg%val, MFLX_z_tavg%val, THERM%val,                & ! (inout)
          DENS_hyd%val, PRES_hyd%val,                                                  & ! (in)
          lcmesh%Gsqrt(:,:), lcmesh%GsqrtH(:,:), lcmesh%GI3(:,:,1), lcmesh%GI3(:,:,2), & ! (in)
          lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3),   & ! (in)
          lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB,                                    & ! (in)
          lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D         ) ! (in)
        call PROF_rapend( 'ATM_DYN_applyBC_mflux', 3)
      end do

      call PROF_rapend( 'ATM_DYN_qtracer', 2)
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
                PRGVAR_DDENS_ID, DDENS, DENS_hyd    )  

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
                DENS_hyd%val(:,:), this%dyn_vars%TRCVARS3D(TRCDDENS_ID)%local(n)%val,      & ! (in)
                this%dyn_vars%TRCVARS3D(TRCDDENS0_ID)%local(n)%val,                        & ! (in)
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
              PRGVAR_DDENS_ID, DDENS, DENS_hyd    )

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
              rkstage, this%dyn_vars%TRCVARS3D(1)%local(n)%val, 1,     &
              1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE,            &
              this%dyn_vars%TRCVARS3D(TRCDDENS_ID)%local(n)%val(:,:),  &
              this%dyn_vars%TRCVARS3D(TRCDDENS0_ID)%local(n)%val(:,:), &
              DENS_hyd%val(:,:)   ) 
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
            PRGVAR_DDENS_ID, DDENS, DENS_hyd    )
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
      call this%numdiff%Apply( prgvars_list, &
        auxvars_list, this%dyncore_driver%boundary_cond, &
        model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),    & ! (in)
        model_mesh%LiftOptrMat, mesh3D                                             ) ! (in)
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

    call this%dyncore_driver%Final()
    call atm_dyn_dgm_trcadvect3d_heve_Final()

    if (this%CALC_NUMDIFF_FLAG) call this%numdiff%Final()

    do n = 1, size(this%tint_qtrc)
      call this%tint_qtrc(n)%Final()
    end do
    deallocate( this%tint_qtrc )    

    call this%dyn_vars%Final()

    return  
  end subroutine AtmosDyn_finalize  

  !--- private ---------------

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

    namelist /PARAM_ATMOS_DYN_TRACER_MODALFILTER/ &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    

    integer :: ierr
    character(len=H_SHORT) :: lbl_readtype
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    select case(read_type)
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
    case ('tracer')
      LOG_NML(PARAM_ATMOS_DYN_TRACER_MODALFILTER)

      call atm_mesh%Construct_ModalFilter3D( &
        this%modal_filter_tracer_3d,         & ! (inout)
        MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   & ! (in)
        MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    ) ! (in)
    end select

    return
  end subroutine setup_modalfilter

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