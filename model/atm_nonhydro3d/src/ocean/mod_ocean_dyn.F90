!-------------------------------------------------------------------------------
!> module Ocean Dynamics
!!
!! @par Description
!!          Module for oceanic dynamical process
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_ocean_dyn
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

!  use mod_ocean_dyn_dgm_driver, only: OceanDynDGMDriver_hydro3d
  
  use mod_ocean_mesh, only: OceanMesh
  use mod_ocean_vars, only: &
    OceanVars_GetLocalMeshPrgVars
  use mod_ocean_dyn_vars, only: &
    OceanDynVars,                      &
    OceanDynAuxVars_GetLocalMeshFields

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage a component of oceanic dynamics
  !!
  type, extends(ModelComponentProc), public :: OceanDyn
!    type(OceanDynDGMDriver_hydro3d) :: dyncore_driver !< A driver object to manage a oceanic dynamical core

    type(OceanDynVars) :: dyn_vars                     !< An object to manage variables in a component of oceanic dynamics
    integer :: eqs_type                                !< Type of governing equations: 0=NONE, 1=SLAB

    class(MeshBase3D), pointer :: mesh3D
    real(RP) :: dtsec                                  !< Timestep for a oceanic dynamical core
  contains
    procedure, public :: setup => OceanDyn_setup 
    procedure, public :: calc_tendency => OceanDyn_calc_tendency
    procedure, public :: update => OceanDyn_update
    procedure, public :: finalize => OceanDyn_finalize
  end type OceanDyn

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: setup_coriolis_parameter

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------

  integer :: OCEAN_DYN_EQS_FIXED_INIT_TEMP = 0
  integer :: OCEAN_DYN_EQS_SLAB            = 1

contains

!> Setup an object to manage a component of oceanic dynamics
!!
!! @param model_mesh Object to manage computational mesh of oceanic model
!! @param tm_parent_comp Object to mange a temporal scheme in a parent component
!!
!OCL SERIAL
  subroutine OceanDyn_setup( this, model_mesh, tm_parent_comp )
    use mod_ocean_mesh, only: OceanMesh
    use scale_time_manager, only: TIME_manager_component

    implicit none

    class(OceanDyn), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    character(len=H_MID) :: EQS_TYPE             = "SLAB"            !< Type of governing equations: FIXED_TEMP or SLAB
    character(len=H_SHORT) :: TINTEG_TYPE        = 'ERK_SSP_3s3o'    !< Type of temporal scheme for a dry dynamical core
    character(len=H_SHORT) :: TINTEG_TYPE_TRACER = 'ERK_SSP_3s3o'    !< Type of temporal scheme for tracer advection equations
    real(DP) :: TIME_DT                          = UNDEF8            !< Timestep for a atmospheric dynamical core
    character(len=H_SHORT) :: TIME_DT_UNIT       = 'SEC'             !< Unit of timestep
    
!    logical :: MODALFILTER_FLAG           = .false. !< Flag to set whether a modal filtering is used
    
    namelist / PARAM_OCEAN_DYN /       &
      EQS_TYPE,                        &
      TINTEG_TYPE,                     &
      TINTEG_TYPE_TRACER,              &      
      TIME_DT,                         &
      TIME_DT_UNIT                    
!      MODALFILTER_FLAG,                
      
    class(OceanMesh), pointer     :: ocean_mesh
    class(MeshBase), pointer      :: ptr_mesh
    class(LocalMeshBase), pointer :: ptr_lcmesh
    class(ElementBase3D), pointer :: elem3D
    integer :: n
    real(DP) :: dtsec

    integer :: ierr
    !--------------------------------------------------

    if (.not. this%IsActivated()) return
    LOG_INFO('OceanDyn_setup',*)

    EQS_TYPE = 'SLAB'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_DYN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("OceanDyn_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("OceanDyn_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_DYN. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_DYN)
    
    !-
    select case (trim(EQS_TYPE))
    case ('FIXED_INIT_TEMP')
      this%eqs_type = OCEAN_DYN_EQS_FIXED_INIT_TEMP
    case ('SLAB')
      this%eqs_type = OCEAN_DYN_EQS_SLAB
    case default
      LOG_INFO('OceanDyn_setup',*) 'Not appropriate EQS_TYPE. Check!'
      call PRC_abort
    end select

    !- get mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    class is (OceanMesh)
      ocean_mesh => model_mesh
    end select
    this%mesh3D => ocean_mesh%ptr_mesh

    !- Setup the temporal integrator

    call tm_parent_comp%Regist_process( 'OCEAN_DYN', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                    ! (out)

    this%dtsec = tm_parent_comp%process_list(this%tm_process_id)%dtsec

    !- initialize the variables 
    call this%dyn_vars%Init( model_mesh )

    call setup_coriolis_parameter( this%dyn_vars, ocean_mesh )

    !- Initialize a module for 3D dynamical core 
    ! call this%dyncore_driver%Init( EQS_TYPE, &
    !   TINTEG_TYPE, dtsec,                    &
    !   MODALFILTER_FLAG,                      &
    !   atm_mesh )

    return
  end subroutine OceanDyn_setup


!> Calculate tendencies associated with oceanic dynamics
!!
!! Because the tendencies with oceanic dynamical cores are treated in OceanDyn_update,
!! no calculation is performed in this subroutine.
!!
!! @param model_mesh Object to manage computational mesh of oceanic model 
!! @param prgvars_list Object to mange prognostic variables with oceanic dynamical core
!! @param trcvars_list Object to mange auxiliary variables 
!! @param forcing_list Object to mange forcing terms
!! @param is_update Flag to specify whether the tendencies are updated in this call
!!
!OCL SERIAL  
  subroutine OceanDyn_calc_tendency( this, model_mesh, prgvars_list, trcvars_list, auxvars_list, forcing_list, is_update )
    implicit none
    
    class(OceanDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    !--------------------------------------------------

    if (.not. this%IsActivated()) return
    !LOG_INFO('OceanDyn_tendency',*)
    return  
  end subroutine OceanDyn_calc_tendency


!> Update variables with a component of oceanic dynamics
!! The tendencies with oceanic dynamical cores are evaluated and the prognostic variables is updated.
!!
!! @param model_mesh Object to manage computational mesh of oceanic model 
!! @param prgvars_list Object to mange prognostic variables with oceanic dynamical core
!! @param trcvars_list Object to mange auxiliary variables 
!! @param forcing_list Object to mange forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL
  subroutine OceanDyn_update( this, model_mesh, prgvars_list, trcvars_list, auxvars_list, forcing_list, is_update )
    use scale_meshfield_base, only: &
      MeshField2D, MeshField3D
    use mod_ocean_vars, only: &
      PRGVAR_THERM_ID,     &
      PHYTEND_RHOH_ID
    implicit none

    class(OceanDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    class(MeshBase), pointer :: mesh
    class(MeshBase3D), pointer :: mesh3D

    integer :: idom
    class(LocalMesh3D), pointer :: lmesh
    class(LocalMesh2D), pointer :: lmesh2D

    class(MeshField3D), pointer :: THERM, RHOH
    !--------------------------------------------------
  
    if ( this%eqs_type == OCEAN_DYN_EQS_FIXED_INIT_TEMP ) then
      LOG_INFO('OceanDyn_update',*) 'No update for fixed initial temperature'
      return
    end if

    call PROF_rapstart( 'OCEAN_DYN_update', 1)

    call model_mesh%GetModelMesh( mesh )
    select type(mesh)
    class is (MeshBase3D)
      mesh3D => mesh
    end select


    call PROF_rapstart( 'OCEAN_DYN_core', 2)
    call prgvars_list%Get3D( PRGVAR_THERM_ID, THERM )
    call forcing_list%Get3D( PHYTEND_RHOH_ID, RHOH )
    do idom=1, mesh3D%LOCAL_MESH_NUM
      lmesh => mesh3D%lcmesh_list(idom)
      lmesh2D => lmesh%lcmesh2D
      call ocn_dyn_update_lc( THERM%local(idom)%val, &
        RHOH%local(idom)%val,                                          &
        this%dtsec, lmesh, lmesh%refElem3D, lmesh2D, lmesh2D%refElem2D )

    end do

    call PROF_rapend( 'OCEAN_DYN_core', 2)

    !---------------------------
    call PROF_rapend( 'OCEAN_DYN_update', 1)

    return  
  end subroutine OceanDyn_update

!> Finalize an object to manage a component of oceanic dynamics
!!
!OCL SERIAL
  subroutine OceanDyn_finalize( this )
    implicit none
    class(OceanDyn), intent(inout) :: this

    integer :: n
    !--------------------------------------------------
    if (.not. this%IsActivated()) return
    LOG_INFO('OceanDyn_finalize',*)

    ! call this%dyncore_driver%Final()
    call this%dyn_vars%Final()

    return  
  end subroutine OceanDyn_finalize

  !--- private ---------------

  !> Setup Coriolis parameter
!OCL SERIAL
  subroutine setup_coriolis_parameter( this, ocn_mesh )

    use scale_coriolis_param, only: get_coriolis_parameter
    use scale_mesh_cubedom3d, only: MeshCubeDom3D
    implicit none

    class(OceanDynVars), target, intent(inout) :: this
    class(OceanMesh), target, intent(in) :: ocn_mesh

    class(LocalMeshFieldBase), pointer :: coriolis
    class(LocalMesh3D), pointer :: lcmesh3D
    class(LocalMesh2D), pointer :: lcmesh2D
    integer :: n

    character(len=H_SHORT) :: CORIOLIS_type !< Type of coriolis force: 'PLANE', 'SPHERE'
    real(RP) :: CORIOLIS_f0   = 0.0_RP
    real(RP) :: CORIOLIS_beta = 0.0_RP
    real(RP) :: CORIOLIS_y0

    namelist /PARAM_OCEAN_DYN_CORIOLIS/ &
      CORIOLIS_type,                         &                
      CORIOLIS_f0, CORIOLIS_beta, CORIOLIS_y0
        
    class(MeshBase3D), pointer :: mesh3D
    class(MeshCubeDom3D), pointer :: meshCube
    integer :: ierr
    !---------------------------------------------------------------

    mesh3D => ocn_mesh%ptr_mesh

    CORIOLIS_type = 'NONE'

    select type(mesh3D)
    type is (MeshCubeDom3D)
      CORIOLIS_y0 = 0.5_RP*(mesh3D%ymax_gl +  mesh3D%ymin_gl)
    end select

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_DYN_CORIOLIS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("OCEAN_DYN_setup_coriolis",*) 'Not found namelist. Default used.'
    else if( ierr > 0 ) then !--- fatal error
      LOG_ERROR("OCEAN_DYN_setup_coriolis",*) 'Not appropriate names in namelist PARAM_OCEAN_DYN_CORIOLIS. Check!'
      call PRC_abort
    end if
    LOG_NML(PARAM_OCEAN_DYN_CORIOLIS)

    do n = 1, mesh3D%LOCAL_MESH_NUM
      call OceanDynAuxVars_GetLocalMeshFields( n, mesh3D, this%AUXVARS2D_manager, &
        coriolis, lcmesh3D )
      lcmesh2D => lcmesh3D%lcmesh2D

      call get_coriolis_parameter( &
        coriolis%val(:,lcmesh2D%NeS:lcmesh2D%NeE),                       & ! (out)
        CORIOLIS_type, lcmesh2D%refElem2D%Np * lcmesh2D%Ne,              & ! (in)
        lcmesh2D%pos_en(:,:,2), CORIOLIS_f0, CORIOLIS_beta, CORIOLIS_y0, & ! (in)
        lcmesh3D%lat2D                                                   ) ! (in)
      
      !$acc update device( coriolis%val )
    end do

    return
  end subroutine setup_coriolis_parameter

!-- private
!OCL SERIAL
  subroutine ocn_dyn_update_lc( THERM, &
    RHOH, dt, lmesh, elem, lmesh2D, elem2D )
    use scale_const, only: &
      DWATR => CONST_DWATR
    use scale_atmos_hydrometeor, only: &
       CV_WATER      
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(inout) :: THERM(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: RHOH(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: dt

    integer :: ke
    real(RP) :: coef
    !----------------------------------------------------------------------

    coef = dt / ( CV_WATER * DWATR )
    !$omp parallel do
    do ke=lmesh%NeS, lmesh%NeE
      THERM(:,ke) = THERM(:,ke) + RHOH(:,ke) * coef
    end do
    return
  end subroutine ocn_dyn_update_lc
end module mod_ocean_dyn