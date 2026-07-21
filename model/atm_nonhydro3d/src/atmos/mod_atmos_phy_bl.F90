!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Planetary Boundary Layer
!!
!! @par Description
!!          Module for planetary boundary layer (PBL) turbulence parameterization
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_bl
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

  use scale_atm_dyn_dgm_bnd, only: AtmDynBnd

  use mod_atmos_phy_bl_vars, only: AtmosPhyBlVars
  use mod_atmos_vars_container, only: &
    AtmosVarsContainer
    

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage a component of planetary boundary layer (PBL) turbulence parameterization in atmospheric model
  !!
  type, extends(ModelComponentProc), public :: AtmosPhyBl
    integer :: BL_TYPEID         !< Type id of PBL turbulence parameterization scheme
    type(AtmosPhyBlVars) :: vars !< Object to manage variables with PBL turbulence parameterization

    integer :: atm_var_container_typeid     !< Type ID of variable container for PBL turbulence parameterization

    real(RP) :: dtsec !< Timestep for PBL turbulence parameterization

    type(LineElement) :: v_elem1D
    type(AtmDynBnd), pointer :: dyn_bnd  !< Pointer to object for treating boundary conditions with atmospheric dynamics  
    real(RP) :: C_IP                     !< Parameter for symmetric interior penalty method in DGM
    logical :: use_delta_form            !< Flag to use delta form in the vertical implicit time integration of PBL scheme
  contains
    procedure :: setup => AtmosPhyBl_setup
    procedure :: calc_tendency => AtmosPhyBl_calc_tendency
    procedure :: update => AtmosPhyBl_update
    procedure :: finalize => AtmosPhyBl_finalize
    procedure, public :: SetDynBC => AtmosPhyBl_SetDynBC
  end type AtmosPhyBl

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

  integer, parameter :: BL_TYPEID_MYNN_LEVEL2 = 1 !< Type ID of MYNN Level 2 PBL scheme
contains

!> Setup a component of planetary boundary layer (PBL) turbulence parameterization in atmospheric model
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param tm_parent_comp Object to mange a temporal scheme in a parent component
!!
  subroutine AtmosPhyBl_setup( this, model_mesh, tm_parent_comp )
    use scale_tracer, only: QA
    use mod_atmos_mesh, only: AtmosMesh
    use scale_time_manager, only: TIME_manager_component

    use scale_atm_phy_bl_dgm_mynn_lv2, only: &
      atm_phy_bl_dgm_mynn_lv2_Init
    use mod_atmos_vars, only: ATM_VARS_CONTAINER_PRIMARY_ID
    implicit none
    class(AtmosPhyBl), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8 !< Timestep for PBL turbulence parameterization
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  !< Unit of timestep

    character(len=H_MID) :: BL_TYPE = 'NONE'      !< Type of a PBL turbulence parameterization scheme
    integer :: atm_var_container_typeid           
    real(RP) :: C_IP = 1.0_RP                     !< Parameter for symmetric interior penalty method in DGM
    logical :: use_delta_form = .false.           !< Flag to use delta form in the vertical implicit time integration of PBL scheme

    namelist /PARAM_ATMOS_PHY_BL/ &
      TIME_DT,                  &
      TIME_DT_UNIT,             &
      BL_TYPE,                  &
      atm_var_container_typeid, &
      C_IP,                     &
      use_delta_form

    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D

    integer :: ierr

    integer :: QS_BL, QE_BL, QA_BL
    !-----------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_BL_setup",*) 'Setup'

    atm_var_container_typeid = ATM_VARS_CONTAINER_PRIMARY_ID

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_BL,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_BL_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_PHY_BL_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_BL. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_BL)
 
    this%atm_var_container_typeid = atm_var_container_typeid
    this%C_IP = C_IP
    this%use_delta_form = use_delta_form
        
    !- Get atmospheric mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select

    !--- Register this component in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_BL', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                       ! (out) 
    
    this%dtsec = tm_parent_comp%process_list(this%tm_process_id)%dtsec

    !--- Set the type of PBL turbulence parameterization

    select case( BL_TYPE )
    case( 'MYNN_LEVEL2' )
      this%BL_TYPEID = BL_TYPEID_MYNN_LEVEL2
      QS_BL = QA
      QA_BL = 0

      call atm_phy_bl_dgm_mynn_lv2_Init( atm_mesh%ptr_mesh )
    case default
      LOG_ERROR("ATMOS_PHY_BL_setup",*) 'Not appropriate PBL turbulence parameterization type. Check!'
      call PRC_abort
    end select

    QE_BL = QS_BL + QA_BL - 1

    !- Initialize the variables 
    call this%vars%Init( model_mesh, QS_BL, QE_BL, QA_BL )

    !-
    call this%v_elem1D%Init( atm_mesh%ptr_mesh%refElem3D%PolyOrder_v, .false. ) 

    return
  end subroutine AtmosPhyBl_setup

!> Calculate tendencies associated with PBL turbulence parameterization in atmospheric model
!!
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param prgvars_list Object to manage prognostic variables with atmospheric dynamical core
!! @param trcvars_list Object to manage auxiliary variables 
!! @param forcing_list Object to manage forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL
  subroutine AtmosPhyBl_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )
    use scale_tracer, only: QA
    use scale_atm_phy_bl_dgm_mynn_lv2, only: &
      atm_phy_bl_dgm_mynn_lv2_cal_VViscDiffCoef
    use scale_atm_phy_bl_dgm_common, only: &
      atm_phy_bl_dgm_common_calc_tendency

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,     &
      AtmosVars_GetLocalMeshPhyAuxVars,  &
      AtmosVars_GetLocalMeshQTRCVarList, & 
      AtmosVars_GetLocalMeshPhyTends      
    use mod_atmos_phy_bl_vars, only: &
      AtmosPhyBLVars_GetLocalMeshFields_tend, &
      RHOU_tp_ID => ATMOS_PHY_BL_RHOU_t_ID, &
      RHOV_tp_ID => ATMOS_PHY_BL_RHOV_t_ID, &
      RHOT_tp_ID => ATMOS_PHY_BL_RHOT_t_ID, &
      TKE_ID => ATMOS_PHY_BL_DIAG_TKE_ID,   &
      NU_ID => ATMOS_PHY_BL_DIAG_NU_ID,     &
      KH_ID => ATMOS_PHY_BL_DIAG_KH_ID

    implicit none
    class(AtmosPhyBl), intent(inout) :: this
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
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: PRES, PT

    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_P
    type(LocalMeshFieldBaseList) :: RHOQ_tp(QA)
    class(LocalMeshFieldBase), pointer :: bl_RHOU_t, bl_RHOV_t, bl_RHOT_t

    type DYN_BNDInfo
      logical, allocatable :: is_bound(:,:)
    end type
    type(DYN_BNDInfo), allocatable :: bnd_info(:)
    !------------------------------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_PROGRESS(*) 'atmosphere / physics / planetary boundary layer' 

    call model_mesh%GetModelMesh( mesh )
    select type(mesh)
    class is (MeshBase3D)
      mesh3D => mesh
    end select

    !-
    if ( is_update ) then
      call PROF_rapstart( 'ATM_BL_tendency', 2)

      allocate( bnd_info(mesh3D%LOCAL_MESH_NUM) )

      do n=1, mesh3D%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshPrgVars( n,    &
          mesh, prgvars_list, auxvars_list,       &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,         &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
          lcmesh                                  )
        
        call AtmosVars_GetLocalMeshPhyAuxVars( n,  &
          mesh, auxvars_list,                      &
          PRES, PT )

        call AtmosPhyBLVars_GetLocalMeshFields_tend( n,  &
          mesh, this%vars%tends_manager,                 &
          bl_RHOU_t, bl_RHOV_t, bl_RHOT_t                )          
        
        !-
        allocate( bnd_info(n)%is_bound(lcmesh%refElem3D%NfpTot,lcmesh%Ne) )
        call this%dyn_bnd%Inquire_bound_flag(  bnd_info(n)%is_bound, & ! (out)
          n, lcmesh%VMapM, lcmesh%VMapP, lcmesh%VMapB,               & ! (in)
          lcmesh, lcmesh%refElem3D                                   ) ! (in)

        select case( this%BL_TYPEID )
        case( BL_TYPEID_MYNN_LEVEL2 )
          call atm_phy_bl_dgm_mynn_lv2_cal_VViscDiffCoef( &
            this%vars%diagvars(NU_ID)%local(n)%val,             & ! (out)
            this%vars%diagvars(KH_ID)%local(n)%val,             & ! (out)
            this%vars%diagvars(TKE_ID)%local(n)%val,            & ! (out)
            DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val, & ! (in)
            DENS_hyd%val, PRES_hyd%val, PRES%val, PT%val,       & ! (in)
            model_mesh%DOptrMat(3), model_mesh%LiftOptrMat,     & ! (in)
            lcmesh, lcmesh%refElem3D, bnd_info(n)%is_bound      ) ! (in)
        end select

        call atm_phy_bl_dgm_common_calc_tendency( &
          bl_RHOU_t%val, bl_RHOV_t%val, bl_RHOT_t%val,    & ! (out)
          DDENS%val, MOMX%val, MOMY%val, DRHOT%val,       & ! (in)
          PT%val, DENS_hyd%val, PRES_hyd%val,             & ! (in)
          this%vars%diagvars(NU_ID)%local(n)%val,         & ! (in)
          this%vars%diagvars(KH_ID)%local(n)%val,         & ! (in)
          model_mesh%element3D_operation,                 & ! (in)
          this%C_IP, this%dtsec,                          & ! (in)
          lcmesh, lcmesh%refElem3D, this%v_elem1D,        & ! (in)
          bnd_info(n)%is_bound, this%use_delta_form       ) ! (in)
      
      end do

      do n=1, mesh3D%LOCAL_MESH_NUM
        deallocate( bnd_info(n)%is_bound )
      end do
      call PROF_rapend( 'ATM_BL_tendency', 2)
    end if

    call PROF_rapstart('ATM_PHY_BL_add_tend', 2)
    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPhyTends( n,        &
        mesh, forcing_list,                          &
        DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, &
        RHOH_p, RHOQ_tp  )

      call AtmosPhyBLVars_GetLocalMeshFields_tend( n,  &
        mesh, this%vars%tends_manager,                 &
        bl_RHOU_t, bl_RHOV_t, bl_RHOT_t,               &
        lcmesh                                         )

      !$omp parallel private(ke, iq)
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
        MOMX_tp%val(:,ke) = MOMX_tp%val(:,ke) + bl_RHOU_t%val(:,ke)
        MOMY_tp%val(:,ke) = MOMY_tp%val(:,ke) + bl_RHOV_t%val(:,ke)
        RHOT_tp%val(:,ke) = RHOT_tp%val(:,ke) + bl_RHOT_t%val(:,ke)
      end do
      !$omp end parallel
    end do
    call PROF_rapend('ATM_PHY_BL_add_tend', 2)


    return
  end subroutine AtmosPhyBl_calc_tendency

!> Update variables in a component of PBL turbulence parameterization in atmospheric model
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param prgvars_list Object to manage prognostic variables with atmospheric dynamical core
!! @param trcvars_list Object to manage auxiliary variables 
!! @param forcing_list Object to manage forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL  
  subroutine AtmosPhyBl_update( this, model_mesh, &
    prgvars_list, trcvars_list,                   &
    auxvars_list, forcing_list, is_update         )  
    
    implicit none
    class(AtmosPhyBl), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------
    return
  end subroutine AtmosPhyBl_update

!> Finalize a component of PBL turbulence parameterization in atmospheric model
!!
!OCL SERIAL  
  subroutine AtmosPhyBl_finalize( this )
    use scale_atm_phy_bl_dgm_mynn_lv2, only: &
      atm_phy_bl_dgm_mynn_lv2_Final
    implicit none
    class(AtmosPhyBl), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    select case ( this%BL_TYPEID )
    case( BL_TYPEID_MYNN_LEVEL2 )
      call atm_phy_bl_dgm_mynn_lv2_Final()
    end select

    call this%vars%Final()
    call this%v_elem1D%Final()
    return
  end subroutine AtmosPhyBl_finalize

!> Set boundary conditions to PBL component in atmospheric model
!!
!! @param dyn_bnd Object to manage boundary conditions of dynamical core
!!
!OCL SERIAL
  subroutine AtmosPhyBl_setDynBC( this, dyn_bnd )
    implicit none
    class(AtmosPhyBl), intent(inout) :: this
    type(AtmDynBnd), intent(in), target :: dyn_bnd
    !--------------------------------------------------

    this%dyn_bnd => dyn_bnd

    return
  end subroutine AtmosPhyBl_setDynBC
!- private ------------------------------------------------

end module mod_atmos_phy_bl
