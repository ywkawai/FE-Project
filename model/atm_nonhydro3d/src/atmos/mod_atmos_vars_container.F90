!-------------------------------------------------------------------------------
!> module Atmosphere / Variables
!!
!! @par Description
!!          Container for atmospheric variables
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_vars_container
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_tracer, only: QA

  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField2D, MeshField3D
  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    PRGVAR_NUM, AUXVAR_NUM, PHYTEND_NUM1 => PHYTEND_NUM,                                        &
    PRGVAR_DDENS_ID, PRGVAR_THERM_ID, PRGVAR_MOMZ_ID, PRGVAR_MOMX_ID, PRGVAR_MOMY_ID,           &
    AUXVAR_DENSHYDRO_ID, AUXVAR_PRESHYDRO_ID, AUXVAR_THERMHYDRO_ID, AUXVAR_PRESHYDRO_REF_ID,    &
    AUXVAR_Rtot_ID, AUXVAR_CPtot_ID, AUXVAR_CVtot_ID,                                           &
    AUXVAR_PRES_ID, AUXVAR_PT_ID, AUXVAR_Qdry_ID,                                               &
    PHYTEND_DENS_ID, PHYTEND_MOMX_ID, PHYTEND_MOMY_ID, PHYTEND_MOMZ_ID, PHYTEND_RHOT_ID,        &
    PHYTEND_RHOH_ID

  use mod_atmos_mesh, only: AtmosMesh  
  use mod_atmos_phy_preproc, only: &
    PhysPreProcBase,        &
    PhysPreProcNone,        &
    PhysPreProcModalFilter, &
    PhysPreProcGlobalFilter
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  !> Derived type to manage a set of variables (prognostic variables, tracer variables, auxiliary variables, and tendencies of physical processes) with atmospheric component
  type, public :: AtmosVarsContainer
    class(AtmosMesh), pointer :: mesh              !< Pointer to an object to manage mesh for atmospheric component

    !- prognostic variables
    type(MeshField3D), allocatable :: PROG_VARS(:) !< Array of 3D prognostic variables
    type(ModelVarManager) :: PROGVARS_manager      !< Object to manage 3D prognostic variables
    integer :: PROG_VARS_commID

    !- tracer variables
    type(MeshField3D), allocatable :: QTRC_VARS(:) !< Array of 3D tracer variables
    type(ModelVarManager) :: QTRCVARS_manager      !< Object to manage 3D tracer variables
    integer :: QTRC_VARS_commID 

    !- auxiliary variables    
    type(MeshField3D), allocatable :: AUX_VARS(:)  !< Array of 3D auxiliary variables
    type(ModelVarManager) :: AUXVARS_manager       !< Object to manage 3D auxiliary variables
    integer :: AUX_VARS_commID

    !- auxiliary variables (2D)
    type(MeshField2D), allocatable :: AUX_VARS2D(:) !< Array of 2D auxiliary variables
    type(ModelVarManager) :: AUXVARS2D_manager      !< Object to manage 3D auxiliary variables

    !- Tendency with physics
    type(MeshField3D), allocatable :: PHY_TEND(:)  !< Array of tendency variables with physics
    type(ModelVarManager) :: PHYTENDS_manager      !< Object to manage tendency variables with physics
    integer :: PHYTENDS_commID
    integer :: PHYTEND_NUM_TOT                     !< Total number of tendency variables with physics

    !-
    integer :: container_type
    integer :: phy_preproc_operation_type
    class(PhysPreProcBase), allocatable :: phy_preproc

  contains
    procedure :: Init => AtmosVarsContainer_Init
    procedure :: Final => AtmosVarsContainer_Final
    procedure :: Calc_diagnostics => AtmosVarsContainer_calculateDiagnostics
    procedure :: Calc_diagVar => AtmosVarsContainer_CalcDiagvar
    procedure :: Calc_SpecificHeat => AtmosVarsContainer_calc_specific_heat
    procedure :: Preproc_operation_for_phys => AtmosVarsContainer_physics_preoperation

    procedure, private :: Setup_phys_preoperation
  end type AtmosVarsContainer

  integer, parameter, public :: ATM_VARS_CONTAINER_PRIMARY_ID = 1

  public :: AtmosVars_GetLocalMeshPrgVar
  public :: AtmosVars_GetLocalMeshPrgVars
  public :: AtmosVars_GetLocalMeshSfcVar
  public :: AtmosVars_GetLocalMeshQTRCVar  
  public :: AtmosVars_GetLocalMeshQTRCVarList
  public :: AtmosVars_GetLocalMeshQTRC_Qv
  public :: AtmosVars_GetLocalMeshPhyAuxVars
  public :: AtmosVars_GetLocalMeshPhyTends
  public :: AtmosVars_GetLocalMeshQTRCPhyTend

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  private :: vars_calc_diagnoseVar_lc

  ! Surface variables
  
  integer, public, parameter :: ATMOS_AUXVARS2D_PREC_ID      = 1
  integer, public, parameter :: ATMOS_AUXVARS2D_PREC_ENGI_ID = 2
  integer, public, parameter :: ATMOS_AUXVARS2D_NUM          = 2

  type(VariableInfo), public :: ATMOS_AUXVARS2D_VINFO(ATMOS_AUXVARS2D_NUM)
  DATA ATMOS_AUXVARS2D_VINFO / &
    VariableInfo( ATMOS_AUXVARS2D_PREC_ID     ,      'PREC', 'surface precipitaion flux'        , 'kg/m2/s', 2, 'XY', 'precipitation_flux'  ), &
    VariableInfo( ATMOS_AUXVARS2D_PREC_ENGI_ID, 'PREC_ENGI', 'internal energy of precipitation' ,    'J/m2', 2, 'XY', ''  )                    /

  !
  integer, parameter :: PHY_PREOPERATION_TYPEID_NONE         = 0
  integer, parameter :: PHY_PREOPERATION_TYPEID_MODALFILTER  = 1
  integer, parameter :: PHY_PREOPERATION_TYPEID_GLOBALFILTER = 2

contains

  !> Initalize an object to manage a set of variables with atmospheric component
!OCL SERIAL
  subroutine AtmosVarsContainer_Init( this, &
    container_type, phy_preproc_file_basename,     &
    atm_mesh )
    use scale_const, only: &
       UNDEF => CONST_UNDEF    
    use scale_tracer, only: &
      TRACER_NAME, TRACER_DESC, TRACER_UNIT
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRGVAR_SCALAR_NUM, PRGVAR_HVEC_NUM,          &
      atm_dyn_dgm_nonhydro3d_common_setup_variables
    
    implicit none
    class(AtmosVarsContainer), intent(inout) :: this
    integer, intent(in) :: container_type
    character(len=*), intent(in) :: phy_preproc_file_basename
    class(AtmosMesh), target, intent(inout) :: atm_mesh

    integer :: iv
    integer :: idom

    type(VariableInfo) :: prgvar_info(PRGVAR_NUM)
    logical :: do_setup_phytend
    logical :: do_setup_auxvar2D
    logical :: reg_file_hist

    class(MeshBase3D), pointer :: mesh3D
    class(MeshBase2D), pointer :: mesh2D
    !---------------------------------------------------------
    
    if ( container_type == ATM_VARS_CONTAINER_PRIMARY_ID ) then
      do_setup_auxvar2D = .true.
      do_setup_phytend = .true.
      reg_file_hist = .true.
    else
      do_setup_auxvar2D = .false.
      do_setup_phytend = .false.
      reg_file_hist = .false.
    end if

    this%mesh => atm_mesh
    mesh3D => atm_mesh%ptr_mesh
    call mesh3D%GetMesh2D( mesh2D )

    !-
    call this%PROGVARS_manager%Init()
    call this%QTRCVARS_manager%Init()
    call this%AUXVARS_manager%Init()

    allocate( this%PROG_VARS(PRGVAR_NUM) )
    allocate( this%QTRC_VARS(0:QA) )
    allocate( this%AUX_VARS(AUXVAR_NUM) )
    !$acc enter data create( this%PROG_VARS, this%QTRC_VARS, this%AUX_VARS )    

    if ( do_setup_phytend ) then
      call this%PHYTENDS_manager%Init()

      this%PHYTEND_NUM_TOT = PHYTEND_NUM1 + max(1,QA)
      allocate( this%PHY_TEND(this%PHYTEND_NUM_TOT) )
      !$acc enter data create( this%PHY_TEND )
    end if

    call atm_dyn_dgm_nonhydro3d_common_setup_variables( &
      this%PROG_VARS, this%QTRC_VARS, this%AUX_VARS, this%PHY_TEND,                              & ! (inout)
      this%PROGVARS_manager, this%QTRCVARS_manager, this%AUXVARS_manager, this%PHYTENDS_manager, & ! (inout)
      reg_file_hist, do_setup_phytend, this%PHYTEND_NUM_TOT, mesh3D,                             & ! (in)
      prgvar_info ) ! (out)
 
    ! Setup communicator
    
    call atm_mesh%Create_communicator( &
      PRGVAR_SCALAR_NUM, PRGVAR_HVEC_NUM, 0,              & ! (in)
      this%PROGVARS_manager,                              & ! (inout)
      this%PROG_VARS(:),                                  & ! (in)
      this%PROG_VARS_commID                               ) ! (out)
        

    if ( QA > 0 ) then
      call atm_mesh%Create_communicator( &
        QA, 0, 0,                        & ! (in)
        this%QTRCVARS_manager,           & ! (inout)
        this%QTRC_VARS(1:QA),            & ! (in)
        this%QTRC_VARS_commID            ) ! (out)
    end if

    call atm_mesh%Create_communicator( &
      AUXVAR_NUM, 0, 0,                & ! (in)
      this%AUXVARS_manager,            & ! (inout)
      this%AUX_VARS(:),                & ! (in)
      this%AUX_VARS_commID             ) ! (out)

    ! Output list of prognostic variables
      
    if ( container_type == ATM_VARS_CONTAINER_PRIMARY_ID ) then
      LOG_NEWLINE
      LOG_INFO("ATMOS_vars_setup",*) 'List of prognostic variables (ATMOS) '
      LOG_INFO_CONT('(1x,A,A24,A,A48,A,A12,A)') &
                '      |', 'VARNAME                 ','|', &
                'DESCRIPTION                                     ', '[', 'UNIT        ', ']'
      do iv = 1, PRGVAR_NUM
        LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
        'NO.',iv,'|',prgvar_info(iv)%NAME,'|', prgvar_info(iv)%DESC,'[', prgvar_info(iv)%UNIT,']'
      end do
      do iv = 1, QA
        LOG_INFO_CONT('(1x,A,I3,A,A24,A,A48,A,A12,A)') &
        'NO.',PRGVAR_NUM+iv,'|',TRACER_NAME(iv),'|', TRACER_DESC(iv),'[', TRACER_UNIT(iv),']'
      end do
      LOG_NEWLINE
    end if
    
    !- Initialize 2D auxiliary variables
    if ( do_setup_auxvar2D ) then
      call this%AUXVARS2D_manager%Init()
      allocate( this%AUX_VARS2D(ATMOS_AUXVARS2D_NUM) )
      !$acc enter data create( this%AUX_VARS2D )
      
      do iv = 1, ATMOS_AUXVARS2D_NUM
        call this%AUXVARS2D_manager%Regist(    &
          ATMOS_AUXVARS2D_VINFO(iv), mesh2D,   & ! (in) 
          this%AUX_VARS2D(iv),                 & ! (inout)
          reg_file_hist                        ) ! (in)
        do idom=1, mesh2D%LOCAL_MESH_NUM
          this%AUX_VARS2D(iv)%local(idom)%val(:,:) = UNDEF
        end do
      end do
    end if

    !- Setup preoperation before physics
    this%container_type = container_type
    call this%Setup_phys_preoperation( phy_preproc_file_basename, mesh3D, mesh3D%refElem3D )

    return
  end subroutine AtmosVarsContainer_Init

  !> Finalize an object to manage a set of variables with atmospheric component
!OCL SERIAL
  subroutine AtmosVarsContainer_Final( this )
    implicit none
    class(AtmosVarsContainer), intent(inout) :: this
    !----------------------------------------

    !$acc exit data delete( this%PROG_VARS, this%QTRC_VARS, this%AUX_VARS )

    call this%PROGVARS_manager%Final()
    deallocate( this%PROG_VARS )
    
    call this%QTRCVARS_manager%Final()
    deallocate( this%QTRC_VARS )

    call this%AUXVARS_manager%Final()
    deallocate( this%AUX_VARS )

    if ( this%container_type == ATM_VARS_CONTAINER_PRIMARY_ID ) then
      !$acc exit data delete( this%PHY_TEND, this%AUX_VARS2D )
      
      call this%PHYTENDS_manager%Final()
      deallocate( this%PHY_TEND )

      call this%AUXVARS2D_manager%Final()
      deallocate( this%AUX_VARS2D )
    else
      ! Finalize preoperation before physics
      call this%phy_preproc%Final()
    end if

    return
  end subroutine AtmosVarsContainer_Final

  !> Calculate diagnostic variables from prognostic variables and auxiliary variables
!OCL SERIAL  
  subroutine AtmosVarsContainer_CalculateDiagnostics( this )    
    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CPdry => CONST_CPdry,    &
      CVdry => CONST_CVdry,    &
      PRES00 => CONST_PRE00
    use scale_tracer, only: &
      TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP    
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat
    implicit none
    class(AtmosVarsContainer), intent(inout), target :: this

    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    integer :: varid

    class(MeshField3D), pointer :: field
    class(ElementBase3D), pointer :: elem3D

    type(LocalMeshFieldBaseList) :: QTRC(QA)
    !-------------------------------------------------------

    ! Calculate specific heat
    call this%Calc_SpecificHeat()
    
    ! Calculate diagnostic variables
    do varid=AUXVAR_THERMHYDRO_ID+1, AUXVAR_PT_ID
      field => this%AUX_VARS(varid)
      do n=1, field%mesh%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshQTRCVarList( n, &
          field%mesh, this%QTRCVARS_manager,       &
          1, QTRC, lcmesh3D )

        elem3D => lcmesh3D%refElem3D

        call vars_calc_diagnoseVar_lc( &
          field%varname, field%local(n)%val,                       &
          this%PROG_VARS(PRGVAR_DDENS_ID)%local(n)%val,            &
          this%PROG_VARS(PRGVAR_MOMX_ID)%local(n)%val,             &
          this%PROG_VARS(PRGVAR_MOMY_ID)%local(n)%val,             &
          this%PROG_VARS(PRGVAR_MOMZ_ID)%local(n)%val,             &
          this%AUX_VARS(AUXVAR_PRES_ID)%local(n)%val,              &
          this%AUX_VARS(AUXVAR_QDRY_ID)%local(n)%val,              &
          QTRC,                                                    &
          this%AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val,         & 
          this%AUX_VARS(AUXVAR_PRESHYDRO_ID)%local(n)%val,         &
          this%AUX_VARS(AUXVAR_Rtot_ID )%local(n)%val,             & 
          this%AUX_VARS(AUXVAR_CVtot_ID)%local(n)%val,             & 
          this%AUX_VARS(AUXVAR_CPtot_ID)%local(n)%val,             & 
          lcmesh3D, lcmesh3D%refElem3D )
      end do
    end do
    !$acc wait(1)

    return
  end subroutine AtmosVarsContainer_CalculateDiagnostics

  !> Calcualte a diagnostic variable from prognostic variables and auxiliary variables
!OCL SERIAL
  subroutine AtmosVarsContainer_CalcDiagvar( this, field_name, field_work ) 
    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CPdry => CONST_CPdry,    &
      CVdry => CONST_CVdry,    &
      PRES00 => CONST_PRE00
    use scale_tracer, only: &
      TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP    
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat

    implicit none
    class(AtmosVarsContainer), intent(inout) :: this
    character(*), intent(in) :: field_name         !< Name of the diagnostic variable to be calculated
    type(MeshField3D), intent(inout) :: field_work !< A work space for the diagnostic variable to be calculated

    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    integer :: ke

    type(MeshField3D) :: field_work_UVmet(2)
    logical :: is_UVmet
    integer :: UVmet_i

    type(LocalMeshFieldBaseList) :: QTRC(QA)
    !--------------------------------------------------

    is_UVmet = .false.
    if ( field_name == 'Umet' ) then
      is_UVmet = .true.; UVmet_i = 1
    else if ( field_name == 'Vmet' ) then      
      is_UVmet = .true.; UVmet_i = 2
    end if

    field_work%varname = field_name

    do n=1, field_work%mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshQTRCVarList( n, &
        field_work%mesh, this%QTRCVARS_manager,  &
        1, QTRC, lcmesh3D )
      
      if ( .not. is_UVmet ) then
        call vars_calc_diagnoseVar_lc( field_name, field_work%local(n)%val,  &
          this%PROG_VARS(PRGVAR_DDENS_ID)%local(n)%val,               &
          this%PROG_VARS(PRGVAR_MOMX_ID)%local(n)%val,                &
          this%PROG_VARS(PRGVAR_MOMY_ID)%local(n)%val,                &
          this%PROG_VARS(PRGVAR_MOMZ_ID)%local(n)%val,                &
          this%AUX_VARS(AUXVAR_PRES_ID)%local(n)%val,                 &
          this%AUX_VARS(AUXVAR_QDRY_ID)%local(n)%val,                 &
          QTRC,                                                       &
          this%AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val,            & 
          this%AUX_VARS(AUXVAR_PRESHYDRO_ID)%local(n)%val,            &
          this%AUX_VARS(AUXVAR_Rtot_ID )%local(n)%val,                & 
          this%AUX_VARS(AUXVAR_CVtot_ID)%local(n)%val,                & 
          this%AUX_VARS(AUXVAR_CPtot_ID)%local(n)%val,                & 
          lcmesh3D, lcmesh3D%refElem3D )
      else
        call field_work_UVmet(1)%Init( 'Umet', '', field_work%mesh )
        call field_work_UVmet(2)%Init( 'Vmet', '', field_work%mesh )
        call this%mesh%Calc_UVmet( &
          this%PROG_VARS(PRGVAR_MOMX_ID), this%PROG_VARS(PRGVAR_MOMY_ID), & ! (in)
          field_work_UVmet(1), field_work_UVmet(2)                        ) ! (inout)
        !$omp parallel do
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          field_work%local(n)%val(:,ke) = field_work_UVmet(UVmet_i)%local(n)%val(:,ke) &
            / ( this%AUX_VARS (AUXVAR_DENSHYDRO_ID)%local(n)%val(:,ke)                 &
              + this%PROG_VARS(PRGVAR_DDENS_ID   )%local(n)%val(:,ke)                  )
        end do
        call field_work_UVmet(1)%Final()
        call field_work_UVmet(2)%Final()
      end if
    end do
    !$acc wait(1)
    return
  end subroutine AtmosVarsContainer_CalcDiagvar

  !> Calculate specific heat coefficients which are weighted by the mass of dry air and tracers
!OCL SERIAL  
  subroutine AtmosVarsContainer_calc_specific_heat( this )
    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CPdry => CONST_CPdry,    &
      CVdry => CONST_CVdry,    &
      PRES00 => CONST_PRE00
    use scale_tracer, only: &
      TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP    
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat
    implicit none
    class(AtmosVarsContainer), intent(inout), target :: this

    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    integer :: ke, p
    integer :: iq

    class(ElementBase3D), pointer :: elem3D

    real(RP), allocatable :: q_tmp(:,:)
    class(LocalMeshFieldBase), pointer :: Qdry, Rtot, CVtot, CPtot
    type(LocalMeshFieldBaseList) :: QTRC(QA)    
    !-------------------------------------------------------

    mesh3D => this%AUX_VARS(1)%mesh

    ! Calculate specific heat
    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(n)
      elem3D => lcmesh3D%refElem3D
      !$acc enter data attach(lcmesh3D, elem3D)

      allocate( q_tmp(elem3D%Np,QA) )

      call AtmosVars_GetLocalMeshQTRCVarList( n, &
        mesh3D, this%QTRCVARS_manager,  &
        1, QTRC, lcmesh3D )
#ifdef _OPENACC
      do iq=1, QA
        !$acc enter data attach(QTRC(iq)%ptr)
      end do
#endif        

      Qdry => this%AUX_VARS(AUXVAR_QDRY_ID )%local(n)
      Rtot => this%AUX_VARS(AUXVAR_Rtot_ID )%local(n)
      CVtot => this%AUX_VARS(AUXVAR_CVtot_ID)%local(n)
      CPtot => this%AUX_VARS(AUXVAR_CPtot_ID)%local(n)
      !$acc enter data attach(Qdry, Rtot, CVtot, CPtot)

      !$omp parallel do private(ke, iq, q_tmp)
      !$acc parallel loop gang private(q_tmp) present(Qdry%val, Rtot%val, CVtot%val, CPtot%val, TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP, lcmesh3D,elem3D)
      do ke = lcmesh3D%NeS, lcmesh3D%NeE
        do iq = 1, QA
        !$acc loop vector
        do p=1, elem3D%Np
          q_tmp(p,iq) = QTRC(iq)%ptr%val(p,ke)
        end do
        end do

        call ATMOS_THERMODYN_specific_heat( &
          elem3D%Np, 1, elem3D%Np, QA,                                       & ! (in)
          q_tmp, TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP,                & ! (in)
          Qdry%val(:,ke), Rtot%val(:,ke), CVtot%val(:,ke), CPtot%val(:,ke)   ) ! (out)
      end do

      deallocate(q_tmp)
    end do
    return
  end subroutine AtmosVarsContainer_calc_specific_heat

  !> Perform preprocessing operation for variables used in physics
!OCL SERIAL  
  subroutine AtmosVarsContainer_physics_preoperation( this, &
    container_ori, dyncore )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_driver_nonhydro3d, only: AtmDynDGMDriver_nonhydro3d

    use scale_atmos_hydrometeor, only: &
      QLA, QIA
    use scale_atm_phy_mp_dgm_common, only: &
      atm_phy_mp_dgm_common_negative_fixer    
    implicit none
    class(AtmosVarsContainer), intent(inout), target :: this
    class(AtmosVarsContainer), intent(in), target :: container_ori
    class(AtmDynDGMDriver_nonhydro3d), intent(in) :: dyncore

    integer :: n
    integer :: iv, iq

    character(len=H_SHORT) :: varname_ori
    character(len=H_SHORT) :: typeid_s

    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: trcid_list(QA)
    type(LocalMeshFieldBaseList) :: lc_qtrc(QA)
    !----------------------------------------------

    if ( this%container_type == ATM_VARS_CONTAINER_PRIMARY_ID ) return

    LOG_INFO("AtmosVarsContainer_physics_preoperation",*) "..."
    
    call this%phy_preproc%Operate( this%PROG_VARS, this%QTRC_VARS, this%AUX_VARS,  & ! (inout)
      container_ori%PROG_VARS, container_ori%QTRC_VARS, container_ori%AUX_VARS ) ! (in)

    call this%Calc_SpecificHeat()

    call dyncore%update_therm_hyd( this%AUXVARS_manager )
    call dyncore%calc_pressure( this%AUX_VARS(AUXVAR_PRES_ID), &
      this%PROGVARS_manager, this%AUXVARS_manager              )

    !- Apply negative fixer for prognostic variables and auxiliary variables which are used in physics
    do iq=1, QA 
      trcid_list(iq) = iq
    end do

    mesh3D => this%AUX_VARS(1)%mesh
    do n=1, mesh3D%LOCAL_MESH_NUM
      call this%QTRCVARS_manager%GetLocalMeshFieldList( trcid_list, n, lc_qtrc )

      lcmesh3D => mesh3D%lcmesh_list(n)
      call atm_phy_mp_dgm_common_negative_fixer( &
        lc_qtrc, this%PROG_VARS(PRGVAR_DDENS_ID)%local(n)%val, this%AUX_VARS(AUXVAR_PRES_ID)%local(n)%val, &
        this%AUX_VARS(AUXVAR_CVtot_ID)%local(n)%val, this%AUX_VARS(AUXVAR_CPtot_ID)%local(n)%val, &
        this%AUX_VARS(AUXVAR_Rtot_ID)%local(n)%val, &        
        this%AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val, this%AUX_VARS(AUXVAR_PRESHYDRO_ID)%local(n)%val, &
        1.0_RP, lcmesh3D, lcmesh3D%refElem3D, QA, QLA, QIA )
    end do

    !-
    call this%Calc_diagnostics()

    !- Tentative : Register the variables for file history

    write(typeid_s,'(I2.2)') this%container_type
    varname_ori = this%PROG_VARS(PRGVAR_DDENS_ID)%varname
    this%PROG_VARS(PRGVAR_DDENS_ID)%varname = trim(varname_ori) // "_"//trim(typeid_s)
    call FILE_HISTORY_meshfield_in( this%PROG_VARS(PRGVAR_DDENS_ID), "" )
    this%PROG_VARS(PRGVAR_DDENS_ID)%varname = varname_ori

    varname_ori = this%AUX_VARS(AUXVAR_PRES_ID)%varname
    this%AUX_VARS(AUXVAR_PRES_ID)%varname = trim(varname_ori) // "_"//trim(typeid_s)
    call FILE_HISTORY_meshfield_in( this%AUX_VARS(AUXVAR_PRES_ID), "" )
    this%AUX_VARS(AUXVAR_PRES_ID)%varname = varname_ori

    do iv=1, 3
      varname_ori = this%QTRC_VARS(iv)%varname
      this%QTRC_VARS(iv)%varname = trim(varname_ori) // "_"//trim(typeid_s)
      call FILE_HISTORY_meshfield_in( this%QTRC_VARS(iv), "" )
      this%QTRC_VARS(iv)%varname = varname_ori
    end do

    return
  end subroutine AtmosVarsContainer_physics_preoperation

  !----  Getter ---------------------------------------------------------------------------

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshPrgVar( domID, mesh, prgvars_list, auxvars_list, &
     varid,                                                                         &
     var, DENS_hyd, PRES_hyd, lcmesh3D                                              )
   
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    integer, intent(in) :: varid
    class(LocalMeshFieldBase), pointer, intent(out) :: var
    class(LocalMeshFieldBase), pointer, intent(out), optional :: DENS_hyd, PRES_hyd
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call prgvars_list%Get(varid, field)
    call field%GetLocalMeshField(domID, var)

    if (present(DENS_hyd)) then
      call auxvars_list%Get(AUXVAR_DENSHYDRO_ID, field)
      call field%GetLocalMeshField(domID, DENS_hyd)
    end if
    if (present(PRES_hyd)) then
      call auxvars_list%Get(AUXVAR_PRESHYDRO_ID, field)
      call field%GetLocalMeshField(domID, PRES_hyd)
    end if

    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshPrgVar

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshPrgVars( domID, mesh, prgvars_list, auxvars_list, &
    DDENS, MOMX, MOMY, MOMZ, THERM,                                                  &
    DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot,                                          &
    lcmesh3D                                                                         )
    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: DDENS, MOMX, MOMY, MOMZ, THERM
    class(LocalMeshFieldBase), pointer, intent(out) :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer, intent(out) :: Rtot, CVtot, CPtot
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call prgvars_list%Get(PRGVAR_DDENS_ID, field)
    call field%GetLocalMeshField(domID, DDENS)

    call prgvars_list%Get(PRGVAR_MOMX_ID, field)
    call field%GetLocalMeshField(domID, MOMX)
    
    call prgvars_list%Get(PRGVAR_MOMY_ID, field)
    call field%GetLocalMeshField(domID, MOMY)

    call prgvars_list%Get(PRGVAR_MOMZ_ID, field)
    call field%GetLocalMeshField(domID, MOMZ)

    call prgvars_list%Get(PRGVAR_THERM_ID, field)
    call field%GetLocalMeshField(domID, THERM)
  
    !--
    call auxvars_list%Get(AUXVAR_DENSHYDRO_ID, field)
    call field%GetLocalMeshField(domID, DENS_hyd)

    call auxvars_list%Get(AUXVAR_PRESHYDRO_ID, field)
    call field%GetLocalMeshField(domID, PRES_hyd)

    call auxvars_list%Get(AUXVAR_Rtot_ID, field)
    call field%GetLocalMeshField(domID, Rtot)

    call auxvars_list%Get(AUXVAR_CVtot_ID, field)
    call field%GetLocalMeshField(domID, CVtot)

    call auxvars_list%Get(AUXVAR_CPtot_ID, field)
    call field%GetLocalMeshField(domID, CPtot)

    !---
    
    if ( present(lcmesh3D) ) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshPrgVars

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshSfcVar( domID, mesh, auxvars2D_list, &
    PREC, PREC_ENGI, lcmesh2D                                           )

   implicit none
   integer, intent(in) :: domID
   class(MeshBase), intent(in) :: mesh
   class(ModelVarManager), intent(inout) :: auxvars2D_list
   class(LocalMeshFieldBase), pointer, intent(out) :: PREC, PREC_ENGI
   class(LocalMesh2D), pointer, intent(out), optional :: lcmesh2D

   class(MeshFieldBase), pointer :: field
   class(LocalMeshBase), pointer :: lcmesh
   !-------------------------------------------------------

   !--
   call auxvars2D_list%Get(ATMOS_AUXVARS2D_PREC_ID, field)
   call field%GetLocalMeshField(domID, PREC)

   call auxvars2D_list%Get(ATMOS_AUXVARS2D_PREC_ENGI_ID, field)
   call field%GetLocalMeshField(domID, PREC_ENGI)

   if (present(lcmesh2D)) then
     call mesh%GetLocalMesh( domID, lcmesh )
     nullify( lcmesh2D )

     select type(lcmesh)
     type is (LocalMesh2D)
       if (present(lcmesh2D)) lcmesh2D => lcmesh
     end select
   end if

   return
 end subroutine AtmosVars_GetLocalMeshSfcVar

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshQTRCVar( domID, mesh, trcvars_list,  &
    varid,                                                              &
    var, lcmesh3D                                                       )

   implicit none
   integer, intent(in) :: domID
   class(MeshBase), intent(in) :: mesh
   class(ModelVarManager), intent(inout) :: trcvars_list
   integer, intent(in) :: varid
   class(LocalMeshFieldBase), pointer, intent(out) :: var
   class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

   class(MeshFieldBase), pointer :: field
   class(LocalMeshBase), pointer :: lcmesh
   !-------------------------------------------------------

   !--
   call trcvars_list%Get(varid, field)
   call field%GetLocalMeshField(domID, var)

   if (present(lcmesh3D)) then
     call mesh%GetLocalMesh( domID, lcmesh )
     nullify( lcmesh3D )

     select type(lcmesh)
     type is (LocalMesh3D)
       if (present(lcmesh3D)) lcmesh3D => lcmesh
     end select
   end if

   return
  end subroutine AtmosVars_GetLocalMeshQTRCVar

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshQTRC_Qv( domID, mesh, trcvars_list, forcing_list, &
    var, var_tp, lcmesh3D                                               )

    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_dry, &
      I_QV    
   implicit none
   integer, intent(in) :: domID
   class(MeshBase), intent(in) :: mesh
   class(ModelVarManager), intent(inout) :: trcvars_list
   class(ModelVarManager), intent(inout) :: forcing_list
   class(LocalMeshFieldBase), pointer, intent(out) :: var
   class(LocalMeshFieldBase), pointer, intent(out), optional :: var_tp
   class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

   class(MeshFieldBase), pointer :: field
   class(LocalMeshBase), pointer :: lcmesh

   integer :: iq, tend_iq
   !-------------------------------------------------------

   !--

   if ( ATMOS_HYDROMETEOR_dry ) then
     iq = 0; tend_iq = PHYTEND_NUM1+1
   else
     iq = I_QV; tend_iq = PHYTEND_NUM1 + I_QV
   end if

   call trcvars_list%Get(iq, field)
   call field%GetLocalMeshField(domID, var)

   call forcing_list%Get(tend_iq, field)
   if (present(var_tp)) call field%GetLocalMeshField(domID, var_tp)

   if (present(lcmesh3D)) then
     call mesh%GetLocalMesh( domID, lcmesh )
     nullify( lcmesh3D )

     select type(lcmesh)
     type is (LocalMesh3D)
       if (present(lcmesh3D)) lcmesh3D => lcmesh
     end select
   end if

   return
  end subroutine AtmosVars_GetLocalMeshQTRC_Qv

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshQTRCVarList( domID, mesh, trcvars_list,  &
    varid_s,                                                                &
    var_list, lcmesh3D                                                      )

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: trcvars_list
    integer, intent(in) :: varid_s
    type(LocalMeshFieldBaseList), intent(out) :: var_list(:)
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh

    integer :: iq
    !-------------------------------------------------------

    !--
    do iq = varid_s, varid_s + size(var_list) - 1
      call trcvars_list%Get(iq, field)
      call field%GetLocalMeshField(domID, var_list(iq-varid_s+1)%ptr)
    end do
    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshQTRCVarList

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshPhyAuxVars( domID, mesh, phyauxvars_list, &
    PRES, PT,                                                                &
    lcmesh3D                                                                 )

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: phyauxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: PRES, PT
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--    
    call phyauxvars_list%Get(AUXVAR_PRES_ID, field)
    call field%GetLocalMeshField(domID, PRES)

    call phyauxvars_list%Get(AUXVAR_PT_ID, field)
    call field%GetLocalMeshField(domID, PT)

    !---
    
    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshPhyAuxVars

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshPhyTends( domID, mesh, phytends_list,  &
    DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p,                  &
    RHOQ_tp,                                                              &
    lcmesh3D                                                              )

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: phytends_list
    class(LocalMeshFieldBase), pointer, intent(out) :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp
    class(LocalMeshFieldBase), pointer, intent(out) :: RHOH_p
    type(LocalMeshFieldBaseList), intent(inout), optional :: RHOQ_tp(QA)
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh

    integer :: iq
    !-------------------------------------------------------

    !--
    call phytends_list%Get(PHYTEND_DENS_ID, field)
    call field%GetLocalMeshField(domID, DENS_tp)

    call phytends_list%Get(PHYTEND_MOMX_ID, field)
    call field%GetLocalMeshField(domID, MOMX_tp)
    
    call phytends_list%Get(PHYTEND_MOMY_ID, field)
    call field%GetLocalMeshField(domID, MOMY_tp)

    call phytends_list%Get(PHYTEND_MOMZ_ID, field)
    call field%GetLocalMeshField(domID, MOMZ_tp)

    call phytends_list%Get(PHYTEND_RHOT_ID, field)
    call field%GetLocalMeshField(domID, RHOT_tp)

    call phytends_list%Get(PHYTEND_RHOH_ID, field)
    call field%GetLocalMeshField(domID, RHOH_p)

    if ( present(RHOQ_tp) ) then
      do iq = 1, QA
        call phytends_list%Get(PHYTEND_NUM1+iq, field)
        call field%GetLocalMeshField(domID, RHOQ_tp(iq)%ptr)  
      end do
    end if

    !---
    if ( present(lcmesh3D) ) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if ( present(lcmesh3D) ) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshPhyTends

!OCL SERIAL
  subroutine AtmosVars_GetLocalMeshQTRCPhyTend( domID, mesh, phytends_list,  &
    qtrcid,                                                                  &
    RHOQ_tp                                                                  )

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: phytends_list
    integer, intent(in) :: qtrcid
    class(LocalMeshFieldBase), pointer, intent(out) :: RHOQ_tp

    class(MeshFieldBase), pointer :: field
    !-------------------------------------------------------

    call phytends_list%Get(PHYTEND_NUM1 + qtrcid, field)
    call field%GetLocalMeshField(domID, RHOQ_tp)

    return
  end subroutine AtmosVars_GetLocalMeshQTRCPhyTend  
  
!-- private -----------------------------------------------------------------------
!OCL SERIAL
  subroutine setup_phys_preoperation( this, phy_preproc_file_basename, mesh3D, elem3D )
    use scale_prc, only: PRC_isMaster
    implicit none
    class(AtmosVarsContainer), intent(inout) :: this
    character(len=*), intent(in) :: phy_preproc_file_basename
    class(MeshBase3D), intent(in) :: mesh3D
    class(ElementBase3D), intent(in) :: elem3D

    integer :: fid
    character(len=H_SHORT) :: PHY_PREPROC_OPERATION_TYPE = 'None' ! ModalFilter / GlobalFilter
    real(RP) :: MF_ALPHA_h
    integer :: MF_ORDER_h
    real(RP) :: MF_ALPHA_v
    integer :: MF_ORDER_v

    character(len=H_SHORT) :: GLFilterOptrType
    character(len=H_SHORT) :: GLFilterShape
    real(RP) :: GLFilterWidthFac
    integer :: Nnode_h1D_reconst
    namelist / PARAM_ATMOS_VARS_CONTAINER / &
      PHY_PREPROC_OPERATION_TYPE, &
      MF_ALPHA_h, MF_ORDER_h, &
      MF_ALPHA_v, MF_ORDER_v, &
      GLFilterOptrType, GLFilterShape, GLFilterWidthFac, Nnode_h1D_reconst

    character(len=H_LONG) :: cnfname
    integer :: ierr

    class(PhysPreProcBase), pointer :: phys_pp_ptr
    !----------------------------------------

    LOG_INFO("ATMOS_vars_container/setup_phys_preoperation",*) "Setting up physical pre-operation for container type ", this%container_type
    
    if ( this%container_type == ATM_VARS_CONTAINER_PRIMARY_ID ) then
      this%phy_preproc_operation_type = PHY_PREOPERATION_TYPEID_NONE
      return
    end if

    !--
    MF_ALPHA_h = 0.0_RP; MF_ORDER_h = 16
    MF_ALPHA_v = 0.0_RP; MF_ORDER_v = 16

    GLFilterOptrType = ''
    GLFilterShape    = 'GAUSSIAN'
    GLFilterWidthFac = 1.0_RP

    Nnode_h1D_reconst = -1

    write(cnfname,'(A,I2.2,A)') trim(phy_preproc_file_basename), this%container_type, ".conf"
    fid = IO_CNF_open(trim(cnfname), PRC_isMaster)

    rewind(fid)
    read(fid,nml=PARAM_ATMOS_VARS_CONTAINER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_vars_container/setup_phys_preoperation",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_vars_container/setup_phys_preoperation",*) 'Invalid names in namelist PARAM_ATMOS_VARS_CONTAINER. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_VARS_CONTAINER)
    close(fid)

    !--
    select case( PHY_PREPROC_OPERATION_TYPE ) 
    case( 'None' )
      this%phy_preproc_operation_type = PHY_PREOPERATION_TYPEID_NONE
      allocate( PhysPreProcNone :: this%phy_preproc )
    case( 'ModalFilter' )
      this%phy_preproc_operation_type = PHY_PREOPERATION_TYPEID_MODALFILTER 
      allocate( PhysPreProcModalFilter :: this%phy_preproc )
    case( 'GlobalFilter' )
      this%phy_preproc_operation_type = PHY_PREOPERATION_TYPEID_GLOBALFILTER
      allocate( PhysPreProcGlobalFilter :: this%phy_preproc )
    case default
      LOG_ERROR("ATMOS_vars_container/setup_phys_preoperation",*) 'Unsupported PHY_PREPROC_OPERATION_TYPE is specified. Check!', PHY_PREPROC_OPERATION_TYPE
      call PRC_abort
    end select

    select type( phys_pp_ptr => this%phy_preproc )
    class is (PhysPreProcNone)
      call phys_pp_ptr%Init( mesh3D )
    class is (PhysPreProcModalFilter)
      call phys_pp_ptr%Init( MF_ALPHA_h, MF_ORDER_h, MF_ALPHA_v, MF_ORDER_v, mesh3D )
    class is (PhysPreProcGlobalFilter)
      if ( Nnode_h1D_reconst < 0 ) Nnode_h1D_reconst = elem3D%Nnode_h1D
      call phys_pp_ptr%Init( GLFilterOptrType, &
        GLFilterShape, GLFilterWidthFac, &
        Nnode_h1D_reconst,               &
        mesh3D )
    end select

    return
  end subroutine setup_phys_preoperation

!OCL SERIAL
  subroutine vars_calc_diagnoseVar_lc( field_name, var_out,  &
    DDENS_, MOMX_, MOMY_, MOMZ_, PRES_, QDRY_, QTRC,         &
    DENS_hyd, PRES_hyd, Rtot, CVtot, CPTot,                  &
    lcmesh, elem )

    use scale_const, only: &
      GRAV => CONST_GRAV,   &
      Rdry => CONST_Rdry,   &
      Rvap => CONST_Rvap,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    use scale_tracer, only: &
      TRACER_inq_id,          &
      TRACER_CV, TRACER_ENGI0
    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_dry      
    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_psat_liq
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    character(*), intent(in) :: field_name
    real(RP), intent(out) :: var_out(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: QDRY_(elem%Np,lcmesh%NeA)
    type(LocalMeshFieldBaseList), intent(in) :: QTRC(QA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot (elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem%Np,lcmesh%NeA)

    integer :: ke, ke2D
    integer :: p
    integer :: iq
    real(RP) :: DENS
    real(RP) :: mom_u1, mom_u2, G_11, G_12, G_22
    real(RP) :: TEMP(elem%Np), PSAT(elem%Np)

    integer :: iq_QV
    integer :: IndexH2Dto3D(elem%Np)

    integer :: Ne, Np
    !-------------------------------------------------------------------------

    Ne = lcmesh%Ne; Np = elem%Np

    select case(trim(field_name))
    case('DENS')
      !$omp parallel do
      !$acc parallel loop collapse(2) present(DDENS_, DENS_hyd, var_out) async(1)
      do ke=1, Ne
      do p=1, Np
        var_out(p,ke) = DDENS_(p,ke) + DENS_hyd(p,ke)
      end do
      end do
    
    case('U')
      !$omp parallel do private (DENS)
      !$acc parallel loop collapse(2) present(DDENS_, DENS_hyd, MOMX_, var_out) async(1)
      do ke=1, Ne
      do p=1, Np
        DENS = DDENS_(p,ke) + DENS_hyd(p,ke)
        var_out(p,ke) = MOMX_(p,ke) / DENS
      end do
      end do
    
    case('V')
      !$omp parallel do private (DENS)
      !$acc parallel loop collapse(2) present(DDENS_, DENS_hyd, MOMY_, var_out) async(1)
      do ke=1, Ne
      do p=1, Np
        DENS = DDENS_(p,ke) + DENS_hyd(p,ke)
        var_out(p,ke) = MOMY_(p,ke) / DENS
      end do
      end do  
          
    case('W')
      !$omp parallel do private (DENS)
      !$acc parallel loop collapse(2) present(DDENS_, DENS_hyd, MOMZ_, var_out) async(1)
      do ke=1, Ne
      do p=1, Np
        DENS = DDENS_(p,ke) + DENS_hyd(p,ke)
        var_out(p,ke) = MOMZ_(p,ke) / DENS
      end do
      end do
    
    case ( 'PRES' )
    case('PRES_diff')  
      !$omp parallel do
      !$acc parallel loop collapse(2) present(PRES_, PRES_hyd, var_out) async(1)
      do ke=1, Ne
      do p=1, Np
        var_out(p,ke) = PRES_(p,ke) - PRES_hyd(p,ke)
      end do
      end do
    
    case('T')
      !$omp parallel do
      !$acc parallel loop collapse(2) present(PRES_, Rtot, DDENS_, DENS_hyd, var_out) async(1)
      do ke=1, Ne
      do p=1, Np
        var_out(p,ke) = PRES_(p,ke) / (Rtot(p,ke) * (DDENS_(p,ke) + DENS_hyd(p,ke)) )
      end do
      end do
    
    case('T_diff')
      !$omp parallel do
      !$acc parallel loop collapse(2) present(PRES_, Rtot, DDENS_, DENS_hyd, PRES_hyd, var_out) async(1)
      do ke=1, Ne
      do p=1, Np
        var_out(p,ke) = PRES_(p,ke) / ( Rtot(p,ke) * (DDENS_(p,ke) + DENS_hyd(p,ke)) ) &
                      - PRES_hyd(p,ke) / ( Rdry * DENS_hyd(p,ke) )
      end do
      end do
    
    case('PT')
      !$omp parallel do private(DENS)
      !$acc parallel loop collapse(2) present(PRES_, Rtot, CVtot, CPtot, DDENS_, DENS_hyd, var_out) async(1)
      do ke=1, Ne
      do p=1, Np
        DENS = DDENS_(p,ke) + DENS_hyd(p,ke)
        var_out(p,ke) = PRES_(p,ke) / (Rtot(p,ke) * DENS ) * ( PRES00 / PRES_(p,ke) )**( Rtot(p,ke) / CPtot(p,ke) )
      end do 
      end do
    
    case('PT_diff')
      !$omp parallel do private( DENS )
      !$acc parallel loop collapse(2) present(PRES_, Rtot, CVtot, CPtot, DDENS_, DENS_hyd, PRES_hyd, var_out) async(1)
      do ke=1, Ne
      do p=1, Np
        DENS = DDENS_(p,ke) + DENS_hyd(p,ke)
        var_out(p,ke) = PRES_(p,ke) / (Rtot(p,ke) * DENS ) * ( PRES00 / PRES_(p,ke) )**( Rtot(p,ke) / CPtot(p,ke) ) &
                      - PRES00/Rdry * (PRES_hyd(p,ke)/PRES00)**(CVdry/CPdry) / DENS_hyd(p,ke)
      end do 
      end do
    
    case( 'RH', 'RHL' )
      if ( ATMOS_HYDROMETEOR_dry ) then
        !$omp parallel do
        !$acc parallel loop collapse(2) present(var_out) async(1)
        do ke=1, Ne
        do p=1, Np
          var_out(p,ke) = 0.0_RP
        end do
        end do
      else
        call TRACER_inq_id( "QV", iq_QV )

        !$omp parallel do private (TEMP, PSAT)
        !$acc parallel loop gang private(TEMP, PSAT) present(PRES_, Rtot, DDENS_, DENS_hyd, var_out) async(1)
        do ke=1, Ne
          !$acc loop vector
          do p=1, Np
            TEMP(p) = PRES_(p,ke) / (Rtot(p,ke) * (DDENS_(p,ke) + DENS_hyd(p,ke)) )
#ifdef _OPENACC
            call ATMOS_SATURATION_psat_liq( TEMP(p), & ! (in)
              PSAT(p)                                ) ! (out)
#endif
          end do

#ifndef _OPENACC
          call ATMOS_SATURATION_psat_liq( &
            Np, 1, Np, TEMP(:),     & ! (in)
            PSAT(:)                           ) ! (out)
#endif

          !$acc loop vector
          do p=1, Np
            var_out(p,ke) = ( DDENS_(p,ke) + DENS_hyd(p,ke) ) * QTRC(iq_QV)%ptr%val(p,ke) &
                          / PSAT(p) * Rvap * TEMP(p) * 100.0_RP
          end do
        end do 
      end if
    
    case('ENGK')
      IndexH2Dto3D(:) = elem%IndexH2Dto3D(:)

      !$omp parallel do private (ke2D, DENS, mom_u1, mom_u2, G_11, G_12, G_22)
      !$acc parallel loop collapse(2) present(DDENS_, DENS_hyd, MOMX_, MOMY_, MOMZ_, var_out, lcmesh,elem) copyin(IndexH2Dto3D) async(1)
      do ke=1, Ne
      do p=1, Np
        ke2D = lcmesh%EMap3Dto2D(ke)
        
        DENS = DDENS_(p,ke) + DENS_hyd(p,ke)
        G_11 = lcmesh%G_ij(IndexH2Dto3D(p),ke2D,1,1)
        G_12 = lcmesh%G_ij(IndexH2Dto3D(p),ke2D,1,2)
        G_22 = lcmesh%G_ij(IndexH2Dto3D(p),ke2D,2,2)

        mom_u1 = G_11 * MOMX_(p,ke) + G_12 * MOMY_(p,ke)
        mom_u2 = G_12 * MOMX_(p,ke) + G_22 * MOMY_(p,ke)

        var_out(p,ke) = 0.5_RP * ( MOMX_(p,ke) * mom_u1 + MOMY_(p,ke) * mom_u2 + MOMZ_(p,ke)**2 ) / DENS
      end do
      end do
    
    case('ENGP')
      !$omp parallel do private (DENS)
      !$acc parallel loop collapse(2) present(DDENS_, DENS_hyd, var_out, lcmesh) async(1)
      do ke=1, Ne
      do p=1, Np
        DENS = DDENS_(p,ke) + DENS_hyd(p,ke)
        var_out(p,ke) = DENS * Grav * lcmesh%zlev(p,ke)
      end do
      end do
    
    case('ENGI')
      !$omp parallel do private (DENS, iq)
      !$acc parallel loop collapse(2) present(DDENS_, DENS_hyd, PRES_, Rtot, QDRY_, var_out) async(1)
      do ke=1, Ne
      do p=1, Np
        DENS = DDENS_(p,ke) + DENS_hyd(p,ke)
        var_out(p,ke) = QDRY_(p,ke) * PRES_(p,ke) / Rtot(p,ke) * CVdry
        !$acc loop seq
        do iq = 1, QA
          var_out(p,ke) = var_out(p,ke) &
            + QTRC(iq)%ptr%val(p,ke) * ( PRES_(p,ke) / Rtot(p,ke) * TRACER_CV(iq) + DENS * TRACER_ENGI0(iq) )
        end do
      end do
      end do
    
    case('ENGT')
      IndexH2Dto3D(:) = elem%IndexH2Dto3D(:)

      !$omp parallel do private (ke2D, DENS, mom_u1, mom_u2, iq, G_11, G_12, G_22)
      !$acc parallel loop collapse(2) present(DDENS_, DENS_hyd, MOMX_, MOMY_, MOMZ_, PRES_, Rtot, QDRY_, var_out, lcmesh, elem) copyin(IndexH2Dto3D) async(1)
      do ke=1, Ne
      do p=1, Np
        ke2D = lcmesh%EMap3Dto2D(ke)

        DENS = DDENS_(p,ke) + DENS_hyd(p,ke)

        G_11 = lcmesh%G_ij(IndexH2Dto3D(p),ke2D,1,1)
        G_12 = lcmesh%G_ij(IndexH2Dto3D(p),ke2D,1,2)
        G_22 = lcmesh%G_ij(IndexH2Dto3D(p),ke2D,2,2)
        mom_u1 = G_11 * MOMX_(p,ke) + G_12 * MOMY_(p,ke)
        mom_u2 = G_12 * MOMX_(p,ke) + G_22 * MOMY_(p,ke)

        ! ENGI
        var_out(p,ke) = QDRY_(p,ke) * PRES_(p,ke) / Rtot(p,ke) * CVdry
        !$acc loop seq
        do iq = 1, QA
          var_out(p,ke) = var_out(p,ke) &
            + QTRC(iq)%ptr%val(p,ke) * ( PRES_(p,ke) / Rtot(p,ke) * TRACER_CV(iq) + DENS * TRACER_ENGI0(iq) )
        end do
        ! ENGT
        var_out(p,ke) = &
            0.5_RP * ( MOMX_(p,ke) * mom_u1 + MOMY_(p,ke) * mom_u2 + MOMZ_(p,ke)**2 ) / DENS & ! ENGK       
          + var_out(p,ke)                                                                    & ! ENGI
          + DENS * Grav * lcmesh%pos_en(p,ke,3)                                                ! ENGP
      end do
      end do
    
    case default
      LOG_ERROR("AtmosVars_calc_diagnoseVar_lc",*) 'The name of diagnostic variable is not suported. Check!', field_name
      call PRC_abort
    
    end select

    return
  end subroutine vars_calc_diagnoseVar_lc
end module mod_atmos_vars_container
  