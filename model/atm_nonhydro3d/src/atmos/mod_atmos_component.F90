!-------------------------------------------------------------------------------
!> module Atmosphere component
!!
!! @par Description
!!          Atmosphere component module
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_component
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc

  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList
  use scale_model_component, only: ModelComponent
  use scale_model_mesh_manager, only: ModelMesh3D

  use mod_atmos_mesh, only: AtmosMesh
  use mod_atmos_mesh_rm, only: AtmosMeshRM
  use mod_atmos_mesh_gm, only: AtmosMeshGM

  use mod_atmos_vars, only: &
    AtmosVars, AtmosVarsContainer

  use mod_atmos_dyn, only: AtmosDyn
  use mod_atmos_phy_sfc, only: AtmosPhySfc
  use mod_atmos_phy_tb , only: AtmosPhyTb
  use mod_atmos_phy_mp , only: AtmosPhyMp

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage atmospheric component
  type, extends(ModelComponent), public :: AtmosComponent
    type(AtmosVars) :: vars             !< Object to mange variables with atmospheric component

    character(len=H_SHORT) :: mesh_type !< Type name of mesh atmospheric component
    class(AtmosMesh), pointer :: mesh   !< Pointer of mesh atmospheric component
    type(AtmosMeshRM) :: mesh_rm        !< Object to manage mesh for the case of regional mode
    type(AtmosMeshGM) :: mesh_gm        !< Object to manage mesh for the case of global mode

    type(AtmosDyn) :: dyn_proc          !< Object to manage dynamical process
    type(AtmosPhySfc) :: phy_sfc_proc   !< Object to manage surface process
    type(AtmosPhyTb ) :: phy_tb_proc    !< Object to manage sub-grid scale turbulence process
    type(AtmosPhyMp ) :: phy_mp_proc    !< Object to manage cloud microphysics process

  contains
    procedure, public :: setup => Atmos_setup 
    procedure, public :: setup_vars => Atmos_setup_vars    
    procedure, public :: calc_tendency => Atmos_calc_tendency
    procedure, public :: update => Atmos_update
    procedure, public :: set_surface => Atmos_set_surface
    procedure, public :: finalize => Atmos_finalize
  end type AtmosComponent

  !-----------------------------------------------------------------------------
  !
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
  !-----------------------------------------------------------------------------
contains

!> Setup an object to mange atmospheric component
!OCL SERIAL
  subroutine Atmos_setup( this )
    use scale_const, only: &
      UNDEF8 => CONST_UNDEF8
    
    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_dry,   &
      ATMOS_HYDROMETEOR_regist    
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_setup
    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_setup

    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_setup  
    use scale_time_manager, only: &
      TIME_manager_Regist_component

    implicit none
    
    class(AtmosComponent), intent(inout), target :: this

    logical :: ACTIVATE_FLAG = .true. !< Flag whether atmospheric component is activated

    real(DP) :: TIME_DT                             = UNDEF8  !< Timestep value of atmospheric component
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'   !< Timestep unit of atmospheric component
    real(DP) :: TIME_DT_RESTART                     = UNDEF8  !< Timestep value when outputting restart file for atmospheric component
    character(len=H_SHORT) :: TIME_DT_RESTART_UNIT  = 'SEC'   !< Timestep unit when outputting restart file for atmospheric component

    logical :: ATMOS_DYN_DO    = .true.  !< Flag whether dynamics process is considered
    logical :: ATMOS_PHY_SF_DO = .false. !< Flag whether surface process is considered
    logical :: ATMOS_PHY_TB_DO = .false. !< Flag whether SGS turbulent process is considered
    logical :: ATMOS_PHY_MP_DO = .false. !< Flag whether cloud microphysics process is considered
    character(len=H_SHORT) :: ATMOS_MESH_TYPE = 'REGIONAL'  !< Name of mesh type for atmospheric component ('REGIONAL' or 'GLOBAL')

    logical :: ATMOS_USE_QV    = .false. !< Flag whether QV is used although cloud microphysics is not considered


    namelist / PARAM_ATMOS / &
      ACTIVATE_FLAG,         &
      TIME_DT,               &
      TIME_DT_UNIT,          &
      TIME_DT_RESTART,       &
      TIME_DT_RESTART_UNIT,  &
      ATMOS_MESH_TYPE,       &
      ATMOS_DYN_DO,          &
      ATMOS_PHY_SF_DO,       &
      ATMOS_PHY_TB_DO,       &
      ATMOS_PHY_MP_DO,       &    
      ATMOS_USE_QV
    
    integer :: ierr
    !--------------------------------------------------
    call PROF_rapstart( 'ATM_setup', 1)
    LOG_INFO('AtmosComponent_setup',*) 'Atmosphere model components '
    
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATM_setup",*) 'Not appropriate names in namelist PARAM_ATMOS. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS)

    !************************************************
    call this%ModelComponent_Init('ATMOS', ACTIVATE_FLAG )
    if ( .not. ACTIVATE_FLAG ) return

    !- Setup time manager

    call this%time_manager%Init( this%GetComponentName(), &
      TIME_DT, TIME_DT_UNIT,                              &
      TIME_DT_RESTART, TIME_DT_RESTART_UNIT               ) 
    
    call TIME_manager_Regist_component( this%time_manager )
    
    !- Setup mesh & file I/O for atmospheric component

    this%mesh_type = ATMOS_MESH_TYPE
    select case( this%mesh_type )
    case('REGIONAL')
      call this%mesh_rm%Init()
      call FILE_HISTORY_meshfield_setup( mesh3d_=this%mesh_rm%mesh )
      this%mesh => this%mesh_rm
    case('GLOBAL')
      call this%mesh_gm%Init()
      call FILE_HISTORY_meshfield_setup( meshCubedSphere3D_=this%mesh_gm%mesh )
      this%mesh => this%mesh_gm
    case default
      LOG_ERROR("ATM_setup",*) 'Unsupported type of mesh is specified. Check!', this%mesh_type
      call PRC_abort    
    end select
    
    !- setup common tools for atmospheric model

    call ATMOS_THERMODYN_setup
    call ATMOS_SATURATION_setup

    !- Setup each processes in atmospheric model ------------------------------------

    !- Setup the module for atmosphere / physics / surface
    call this%phy_sfc_proc%ModelComponentProc_Init( 'AtmosPhysSfc', ATMOS_PHY_SF_DO )
    call this%phy_sfc_proc%setup( this%mesh, this%time_manager )
    
    !- Setup the module for atmosphere / physics / cloud microphysics
    call this%phy_mp_proc%ModelComponentProc_Init( 'AtmosPhysMp', ATMOS_PHY_MP_DO )
    call this%phy_mp_proc%setup( this%mesh, this%time_manager )

    !-- Regist qv if needed
    if ( ATMOS_HYDROMETEOR_dry .and. ATMOS_USE_QV ) then
      LOG_INFO("ATMOS_setup",*) "Regist QV"
      call ATMOS_HYDROMETEOR_regist( 0, 0,                                 & ! (in)
        (/'QV'/),                                                          & ! (in)
        (/'Ratio of Water Vapor mass to total mass (Specific humidity)'/), & ! (in)
        (/'kg/kg'/),                                                       & ! (in)
        this%phy_mp_proc%vars%QS                                           ) ! (out)
      
      this%phy_mp_proc%vars%QA = 1
      this%phy_mp_proc%vars%QE = this%phy_mp_proc%vars%QS    
    end if

    !- Setup the module for atmosphere / dynamics 
    call this%dyn_proc%ModelComponentProc_Init( 'AtmosDyn', ATMOS_DYN_DO )
    call this%dyn_proc%setup( this%mesh, this%time_manager )

    !- Setup the module for atmosphere / physics / turbulence
    call this%phy_tb_proc%ModelComponentProc_Init( 'AtmosPhysTb', ATMOS_PHY_TB_DO )
    call this%phy_tb_proc%setup( this%mesh, this%time_manager )
    call this%phy_tb_proc%SetDynBC( this%dyn_proc%dyncore_driver%boundary_cond )

    LOG_NEWLINE
    LOG_INFO('AtmosComponent_setup',*) 'Finish setup of each atmospheric components.'

    call PROF_rapend( 'ATM_setup', 1)

    return
  end subroutine Atmos_setup

!> Setup variables with the atmospheric component
!OCL SERIAL
  subroutine Atmos_setup_vars( this )
    implicit none

    class(AtmosComponent), intent(inout) :: this
    !----------------------------------------------------------

    call PROF_rapstart( 'ATM_setup_vars', 1)

    call this%vars%Init( this%mesh )
    call this%vars%Regist_physvar_manager( &
      this%phy_mp_proc%vars%auxvars2D_manager )

    call this%vars%Setup_container( this%phy_mp_proc%coarsend_dynvars_typeID, this%mesh )
    call PROF_rapend( 'ATM_setup_vars', 1)   

    return
  end subroutine Atmos_setup_vars

!> Calculate tendencies with the atmospheric component
!OCL SERIAL
  subroutine Atmos_calc_tendency( this, force )
    use scale_tracer, only: QA
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PHYTEND_NUM1 => PHYTEND_NUM, &
      DENS_tp => PHYTEND_DENS_ID,  &
      MOMX_tp => PHYTEND_MOMX_ID,  &
      MOMY_tp => PHYTEND_MOMY_ID,  &
      MOMZ_tp => PHYTEND_MOMZ_ID,  &
      RHOT_tp =>  PHYTEND_RHOT_ID, &
      RHOH_p => PHYTEND_RHOH_ID    
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPhyTends

    implicit none
    class(AtmosComponent), intent(inout), target :: this
    logical, intent(in) :: force


    class(MeshBase), pointer :: mesh
    class(LocalMesh3D), pointer :: lcmesh
    type(LocalMeshFieldBaseList) :: tp_list(PHYTEND_NUM1)
    type(LocalMeshFieldBaseList) :: tp_qtrc(QA)

    integer :: tm_process_id
    logical :: is_update
    integer :: n
    integer :: v
    integer :: iq
    integer :: ke

    type(AtmosVarsContainer), pointer :: vars_container
    !------------------------------------------------------------------
    
    call PROF_rapstart( 'ATM_tendency', 1)
    !LOG_INFO('AtmosComponent_calc_tendency',*)

    call this%mesh%GetModelMesh( mesh )  

    !########## Get Surface Boundary from coupler ##########
    
    
    !########## calculate tendency ##########

    !* Exchange halo data ( for physics )
    call PROF_rapstart( 'ATM_exchange_prgv', 2)
    call this%vars%container%PROGVARS_manager%MeshFieldComm_Exchange()
    if ( QA > 0 ) call this%vars%container%QTRCVARS_manager%MeshFieldComm_Exchange()
    call PROF_rapend( 'ATM_exchange_prgv', 2)

    ! reset tendencies of physics

    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPhyTends( n, mesh, this%vars%container%PHYTENDS_manager,  & ! (in)
        tp_list(DENS_tp)%ptr, tp_list(MOMX_tp)%ptr, tp_list(MOMY_tp)%ptr,                  & ! (out)
        tp_list(MOMZ_tp)%ptr, tp_list(RHOT_tp)%ptr, tp_list(RHOH_p )%ptr, tp_qtrc,         & ! (out)
        lcmesh                                                                             ) ! (out)
      
      !$omp parallel private(v,iq,ke)
      !$omp do collapse(2)
      do v=1, PHYTEND_NUM1
      do ke=lcmesh%NeS, lcmesh%NeE
        tp_list(v)%ptr%val(:,ke) = 0.0_RP
      end do
      end do
      !$omp do collapse(2)
      do iq=1, QA
      do ke=lcmesh%NeS, lcmesh%NeE
        tp_qtrc(iq)%ptr%val(:,ke) = 0.0_RP
      end do
      end do
      !$omp end do
      !$omp end parallel
    end do

    ! Cloud Microphysics

    if ( this%phy_mp_proc%IsActivated() ) then
      call PROF_rapstart('ATM_Microphysics', 1)
      tm_process_id = this%phy_mp_proc%tm_process_id
      is_update = this%time_manager%Do_process(tm_process_id) .or. force

      vars_container => this%vars%container_list(this%phy_mp_proc%coarsend_dynvars_typeID)
      call this%phy_mp_proc%calc_tendency( &
          this%mesh, vars_container%PROGVARS_manager, vars_container%QTRCVARS_manager, &
          vars_container%AUXVARS_manager, vars_container%PHYTENDS_manager, is_update   )
      call PROF_rapend('ATM_Microphysics', 1)
    end if
    
    ! Radiation


    ! Turbulence

    if ( this%phy_tb_proc%IsActivated() ) then
      call PROF_rapstart('ATM_Turbulence', 1)
      tm_process_id = this%phy_tb_proc%tm_process_id
      is_update = this%time_manager%Do_process(tm_process_id) .or. force

      vars_container => this%vars%container
      call this%phy_tb_proc%calc_tendency( &
        this%mesh, vars_container%PROGVARS_manager, vars_container%QTRCVARS_manager, &
        vars_container%AUXVARS_manager, vars_container%PHYTENDS_manager, is_update   )
      call PROF_rapend('ATM_Turbulence', 1)
    end if

    ! Cumulus


!    if ( .not. CPL_sw ) then    
    
    ! Surface flux
    
    if ( this%phy_sfc_proc%IsActivated() ) then
      call PROF_rapstart('ATM_SurfaceFlux', 1)
      tm_process_id = this%phy_sfc_proc%tm_process_id
      is_update = this%time_manager%Do_process(tm_process_id) .or. force

      vars_container => this%vars%container
      call this%phy_sfc_proc%calc_tendency( &
        this%mesh, vars_container%PROGVARS_manager, vars_container%QTRCVARS_manager, &
        vars_container%AUXVARS_manager, vars_container%PHYTENDS_manager, is_update   )
      call PROF_rapend('ATM_SurfaceFlux', 1)
    end if
    
    ! Planetary Boundary layer

!   end if

    call PROF_rapend( 'ATM_tendency', 1)
    return  
  end subroutine Atmos_calc_tendency

!> Update variables with the atmospheric component
!OCL SERIAL
  subroutine Atmos_update( this )
    implicit none
    class(AtmosComponent), intent(inout), target :: this
    
    integer :: tm_process_id
    logical :: is_update
    integer :: inner_itr

    type(AtmosVarsContainer), pointer :: vars_container
    !--------------------------------------------------
    call PROF_rapstart( 'ATM_update', 1)

    !########## Dynamics ########## 

    vars_container => this%vars%container ! container_type = 1

    if ( this%dyn_proc%IsActivated() ) then
      call PROF_rapstart('ATM_Dynamics', 1)
      tm_process_id = this%dyn_proc%tm_process_id
      is_update = this%time_manager%Do_process( tm_process_id )

      LOG_PROGRESS(*) 'atmosphere / dynamics'   
      do inner_itr=1, this%time_manager%Get_process_inner_itr_num( tm_process_id )
        call this%dyn_proc%update( &
          this%mesh, vars_container%PROGVARS_manager, vars_container%QTRCVARS_manager, &
          vars_container%AUXVARS_manager, vars_container%PHYTENDS_manager, is_update   )
      end do
      call PROF_rapend('ATM_Dynamics', 1)
    end if
    
    !########## Calculate diagnostic variables ##########  

    call vars_container%Calc_diagnostics()
    call vars_container%AUXVARS_manager%MeshFieldComm_Exchange()

    !########## Adjustment ##########
    ! Microphysics
    ! Aerosol
    ! Lightning

    !########## Reference State ###########

    !#### Check values #################################
    call this%vars%Check()

    call PROF_rapend('ATM_update', 1)
    return  
  end subroutine Atmos_update

!OCL SERIAL
  subroutine Atmos_set_surface( this )
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshSfcVar
    use mod_atmos_phy_mp_vars, only: &
      AtmosPhyMpVars_GetLocalMeshFields_sfcflx
    implicit none
    class(AtmosComponent), intent(inout) :: this

    class(MeshBase), pointer :: mesh
    class(MeshBase2D), pointer :: mesh2D
    class(LocalMesh2D), pointer :: lcmesh
    integer :: n
    integer :: ke

    class(LocalMeshFieldBase), pointer :: PREC, PREC_ENGI
    class(LocalMeshFieldBase), pointer :: SFLX_rain_MP, SFLX_snow_MP, SFLX_ENGI_MP

    !--------------------------------------------------

    call PROF_rapstart( 'ATM_sfc_exch', 1)

    call this%mesh%GetModelMesh( mesh )  
    select type(mesh)
    class is (MeshBase3D)
      call mesh%GetMesh2D( mesh2D )
    end select

    !- sum of rainfall from mp and cp

    do n=1, mesh2D%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshSfcVar( n, &
        mesh2D, this%vars%container%AUXVARS2D_manager, & ! (in)
        PREC, PREC_ENGI, lcmesh                        ) ! (out)

      !$omp parallel do private(ke)
      do ke=lcmesh%NeS, lcmesh%NeE
        PREC     %val(:,ke) = 0.0_RP
        PREC_ENGI%val(:,ke) = 0.0_RP
      end do

      if ( this%phy_mp_proc%IsActivated() ) then
        call AtmosPhyMpVars_GetLocalMeshFields_sfcflx( n, &
          mesh2D, this%phy_mp_proc%vars%auxvars2D_manager, & ! (in)
          SFLX_rain_MP, SFLX_snow_MP, SFLX_ENGI_MP         ) ! (out)

        !$omp parallel do private(ke)
        do ke=lcmesh%NeS, lcmesh%NeE
          PREC     %val(:,ke) = PREC     %val(:,ke) + SFLX_rain_MP%val(:,ke) + SFLX_snow_MP%val(:,ke)
          PREC_ENGI%val(:,ke) = PREC_ENGI%val(:,ke) + SFLX_ENGI_MP%val(:,ke)
        end do
      end if
   end do

    call PROF_rapend( 'ATM_sfc_exch', 1)

    return
  end subroutine Atmos_set_surface

!> Finalize an object to manage the atmospheric component
!OCL SERIAL
  subroutine Atmos_finalize( this )
    implicit none
    class(AtmosComponent), intent(inout) :: this
    !--------------------------------------------------

    LOG_INFO('AtmosComponent_finalize',*)

    if ( .not. this%IsActivated() ) return

    call this%dyn_proc%finalize()
    call this%phy_sfc_proc%finalize()
    call this%phy_tb_proc%finalize()
    call this%phy_mp_proc%finalize()
    
    call this%vars%Final()

    select case( this%mesh_type )
    case('REGIONAL')
      call this%mesh_rm%Final()
    case('GLOBAL')
      call this%mesh_gm%Final()
    end select  
    this%mesh => null()

    call this%time_manager%Final()

    return  
  end subroutine Atmos_finalize

end module mod_atmos_component