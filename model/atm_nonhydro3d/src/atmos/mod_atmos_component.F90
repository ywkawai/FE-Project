!-------------------------------------------------------------------------------
!> module ATMOSPHERE component
!!
!! @par Description
!!          Atmosphere component module
!!
!! @author Team SCALE
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
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_model_component, only: ModelComponent

  use mod_atmos_vars, only: AtmosVars
  use mod_atmos_mesh, only: AtmosMesh

  use mod_atmos_dyn, only: AtmosDyn
  use mod_atmos_phy_sfc, only: AtmosPhySfc
  use mod_atmos_phy_tb , only: AtmosPhyTb

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, extends(ModelComponent), public :: AtmosComponent
    type(AtmosVars) :: vars
    type(AtmosMesh) :: mesh

    type(AtmosDyn) :: dyn_proc
    type(AtmosPhySfc) :: phy_sfc_proc
    type(AtmosPhyTb ) :: phy_tb_proc
  contains
    procedure, public :: setup => Atmos_setup 
    procedure, public :: calc_tendency => Atmos_calc_tendency
    procedure, public :: update => Atmos_update
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
subroutine Atmos_setup( this )
  use scale_const, only: &
    UNDEF8 => CONST_UNDEF8
  use scale_atmos_hydrometeor, only: &
    ATMOS_HYDROMETEOR_setup

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_setup  
  use scale_time_manager, only: &
    TIME_manager_Regist_component

  implicit none

  class(AtmosComponent), intent(inout) :: this

  logical :: ACTIVATE_FLAG = .true.

  real(DP) :: TIME_DT                             = UNDEF8
  character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'
  real(DP) :: TIME_DT_RESTART                     = UNDEF8
  character(len=H_SHORT) :: TIME_DT_RESTART_UNIT  = 'SEC'

  logical :: ATMOS_DYN_DO    = .true.
  logical :: ATMOS_PHY_SF_DO = .false.
  logical :: ATMOS_PHY_TB_DO = .false.

  namelist / PARAM_ATMOS / &
    ACTIVATE_FLAG,         &
    TIME_DT,               &
    TIME_DT_UNIT,          &
    TIME_DT_RESTART,       &
    TIME_DT_RESTART_UNIT,  &
    ATMOS_DYN_DO,          &
    ATMOS_PHY_SF_DO,       &
    ATMOS_PHY_TB_DO

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
  
  !- Setup mesh

  call this%mesh%Init()

  !- Setup file I/O for atmospheric component

  call FILE_HISTORY_meshfield_setup( mesh3d_=this%mesh%mesh )  
  
  !- Setup variables

  call this%vars%Init( this%mesh )
  
  !- Setup each processes in atmospheric model ------------------------------------

  !- Setup index for atmospheric tracers
  call ATMOS_HYDROMETEOR_setup()
    
  !- Setup the module for atmosphere / dynamics 
  call this%dyn_proc%ModelComponentProc_Init( 'AtmosDyn', ATMOS_DYN_DO )
  call this%dyn_proc%setup( this%mesh, this%time_manager )

  !- Setup the module for atmosphere / physics / surface
  call this%phy_sfc_proc%ModelComponentProc_Init( 'AtmosPhysSfc', ATMOS_PHY_SF_DO )
  call this%phy_sfc_proc%setup( this%mesh, this%time_manager )

  !- Setup the module for atmosphere / physics / turbulence
  call this%phy_tb_proc%ModelComponentProc_Init( 'AtmosPhysTb', ATMOS_PHY_TB_DO )
  call this%phy_tb_proc%setup( this%mesh, this%time_manager )
  call this%phy_tb_proc%SetDynBC( this%dyn_proc%boundary_cond )
  
  call PROF_rapend( 'ATM_setup', 1)
  return
end subroutine Atmos_setup

subroutine Atmos_calc_tendency( this )

  use mod_atmos_vars, only: &
    AtmosVars_GetLocalMeshPhyTends, &
    ATMOS_PHYTEND_NUM,                &
    DENS_tp => ATMOS_PHYTEND_DENS_ID, &
    MOMX_tp => ATMOS_PHYTEND_MOMX_ID, &
    MOMY_tp => ATMOS_PHYTEND_MOMY_ID, &
    MOMZ_tp => ATMOS_PHYTEND_MOMZ_ID, &
    RHOT_tp => ATMOS_PHYTEND_RHOT_ID, &
    RHOH_p  => ATMOS_PHYTEND_RHOH_ID

  implicit none
  class(AtmosComponent), intent(inout) :: this


  class(MeshBase), pointer :: mesh
  class(LocalMesh3D), pointer :: lcmesh
  type :: PhysTendPtrList
    class(LocalMeshFieldBase), pointer :: ptr => null()
  end type PhysTendPtrList
  type(PhysTendPtrList) :: tp_list(ATMOS_PHYTEND_NUM)

  integer :: tm_process_id
  logical :: is_update
  integer :: n
  integer :: v
  integer :: ke
  !------------------------------------------------------------------
  
  call PROF_rapstart( 'ATM_tendency', 1)
  !LOG_INFO('AtmosComponent_calc_tendency',*)

  !########## calculate tendency ##########

  !* Exchange halo data ( for physics )
  call PROF_rapstart( 'ATM_exchange_prgv', 2)
  call this%vars%PROGVARS_manager%MeshFieldComm_Exchange()
  call PROF_rapend( 'ATM_exchange_prgv', 2)

  ! reset tendencies of physics

  call this%mesh%GetModelMesh( mesh )  
  do n=1, mesh%LOCAL_MESH_NUM
    call AtmosVars_GetLocalMeshPhyTends( n, mesh, this%vars%PHYTENDS_manager ,  & ! (in)
      tp_list(DENS_tp)%ptr, tp_list(MOMX_tp)%ptr, tp_list(MOMY_tp)%ptr,         & ! (out)
      tp_list(MOMZ_tp)%ptr, tp_list(RHOT_tp)%ptr, tp_list(RHOH_p )%ptr, lcmesh  ) ! (out)
    
    !$omp parallel do collapse(2)
    do v=1, ATMOS_PHYTEND_NUM
    do ke=lcmesh%NeS, lcmesh%NeE
      tp_list(v)%ptr%val(:,ke) = 0.0_RP
    end do
    end do
  end do

  ! Surface flux
  if ( this%phy_sfc_proc%IsActivated() ) then
    call PROF_rapstart('ATM_SurfaceFlux', 1)
    tm_process_id = this%phy_sfc_proc%tm_process_id
    is_update = this%time_manager%Do_process(tm_process_id)
    call this%phy_sfc_proc%calc_tendency( &
        this%mesh, this%vars%PROGVARS_manager, this%vars%AUXVARS_manager, this%vars%PHYTENDS_manager, is_update )
    call PROF_rapend('ATM_SurfaceFlux', 1)
  end if
  
  ! Turbulence
  if ( this%phy_tb_proc%IsActivated() ) then
    call PROF_rapstart('ATM_Turbulence', 1)
    tm_process_id = this%phy_tb_proc%tm_process_id
    is_update = this%time_manager%Do_process(tm_process_id)
    call this%phy_tb_proc%calc_tendency( &
        this%mesh, this%vars%PROGVARS_manager, this%vars%AUXVARS_manager, this%vars%PHYTENDS_manager, is_update )
    call PROF_rapend('ATM_Turbulence', 1)
  end if

  call PROF_rapend( 'ATM_tendency', 1)
  return  
end subroutine Atmos_calc_tendency

subroutine Atmos_update( this )
  implicit none
  class(AtmosComponent), intent(inout) :: this
  
  integer :: tm_process_id
  logical :: is_update
  integer :: inner_itr
  !--------------------------------------------------
  call PROF_rapstart( 'ATM_update', 1)

  !########## Dynamics ########## 

  if ( this%dyn_proc%IsActivated() ) then
    call PROF_rapstart('ATM_Dynamics', 1)
    tm_process_id = this%dyn_proc%tm_process_id
    is_update = this%time_manager%Do_process( tm_process_id )

    LOG_PROGRESS(*) 'atmosphere / dynamics'   
    do inner_itr=1, this%time_manager%Get_process_inner_itr_num( tm_process_id )
      call this%dyn_proc%update( &
        this%mesh, this%vars%PROGVARS_manager, this%vars%AUXVARS_manager, this%vars%PHYTENDS_manager, is_update )
    end do
    call PROF_rapend('ATM_Dynamics', 1)
  end if
  
  !########## Calculate diagnostic variables ##########  

  call this%vars%Clac_diagnostics()
  call this%vars%AUXVARS_manager%MeshFieldComm_Exchange()

  !########## Adjustment ##########
  ! Microphysics
  ! Aerosol
  ! Lightning

  !########## Reference State ###########

  !#### Check values #################################
  call this%vars%Check()

  call PROF_rapend( 'ATM_update', 1)
  return  
end subroutine Atmos_update

subroutine Atmos_finalize( this )
  implicit none
  class(AtmosComponent), intent(inout) :: this
  !--------------------------------------------------

  LOG_INFO('AtmosComponent_finalize',*)

  if ( .not. this%IsActivated() ) return

  call this%dyn_proc%finalize()
  call this%phy_sfc_proc%finalize()
  call this%phy_tb_proc%finalize()
  
  call this%vars%Final()
  call this%mesh%Final()
  call this%time_manager%Final()

  return  
end subroutine Atmos_finalize

end module mod_atmos_component