!-------------------------------------------------------------------------------
!> module global shallow water component
!!
!! @par Description
!!          global shallow water component module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_globalsw_component
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc

  use scale_mesh_base, only: MeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_model_component, only: ModelComponent

  use mod_sw_vars, only: SWVars
  use mod_sw_mesh, only: SWMesh
  use mod_sw_dyn, only: SWDyn

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, extends(ModelComponent), public :: GlobalSWComponent
   type(SWMesh) :: mesh
   type(SWVars) :: vars
   type(SWDyn) :: dyn_proc
  contains
    procedure, public :: setup => SW_setup 
    procedure, public :: calc_tendency => SW_calc_tendency
    procedure, public :: update => SW_update
    procedure, public :: finalize => SW_finalize
  end type GlobalSWComponent

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
subroutine SW_setup( this )
  use scale_const, only: &
    UNDEF8 => CONST_UNDEF8
  use scale_atmos_hydrometeor, only: &
    ATMOS_HYDROMETEOR_setup

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_setup  
  use scale_time_manager, only: &
    TIME_manager_Regist_component

  implicit none

  class(GlobalSWComponent), intent(inout), target :: this

  logical :: ACTIVATE_FLAG = .true.

  real(DP) :: TIME_DT                             = UNDEF8
  character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'
  real(DP) :: TIME_DT_RESTART                     = UNDEF8
  character(len=H_SHORT) :: TIME_DT_RESTART_UNIT  = 'SEC'

  logical :: GlobalSW_DYN_DO    = .true.

  namelist / PARAM_GLOBALSW / &
    ACTIVATE_FLAG,         &
    TIME_DT,               &
    TIME_DT_UNIT,          &
    TIME_DT_RESTART,       &
    TIME_DT_RESTART_UNIT,  &
    GlobalSW_DYN_DO

  integer :: ierr
  !--------------------------------------------------
  call PROF_rapstart( 'GlobalSW_setup', 1)
  LOG_INFO('GlobalShallowWater_setup',*) 'Global shallow water model components '
  
  !--- read namelist
  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_GLOBALSW,iostat=ierr)
  if( ierr < 0 ) then !--- missing
     LOG_INFO("GlobalShallowWater_setup",*) 'Not found namelist. Default used.'
  elseif( ierr > 0 ) then !--- fatal error
     LOG_ERROR("GlobalShallowWater_setup",*) 'Not appropriate names in namelist PARAM_GLOBALSW. Check!'
     call PRC_abort
  endif
  LOG_NML(PARAM_GLOBALSW)

  !************************************************
  call this%ModelComponent_Init('GlobalSW', ACTIVATE_FLAG )
  if ( .not. ACTIVATE_FLAG ) return

  !- Setup time manager

  call this%time_manager%Init( this%GetComponentName(), &
    TIME_DT, TIME_DT_UNIT,                              &
    TIME_DT_RESTART, TIME_DT_RESTART_UNIT               ) 
  
  call TIME_manager_Regist_component( this%time_manager )
  
  !- Setup mesh

  call this%mesh%Init()

  !- Setup file I/O for atmospheric component

  call FILE_HISTORY_meshfield_setup( meshcubedsphere2D_=this%mesh%mesh )  
  
  !- Setup variables

  call this%vars%Init( this%mesh )
  
  !- Setup each processes in atmospheric model ------------------------------------

  !- Setup index for atmospheric tracers
  call ATMOS_HYDROMETEOR_setup()
    
  !- Setup the module for atmosphere / dynamics 
  call this%dyn_proc%ModelComponentProc_Init( 'SWDyn', GlobalSW_DYN_DO )
  call this%dyn_proc%setup( this%mesh, this%time_manager )

  call PROF_rapend( 'GlobalSW_setup', 1)
  return
end subroutine SW_setup

subroutine SW_update( this )
  implicit none
  class(GlobalSWComponent), intent(inout) :: this

  integer :: tm_process_id
  logical :: is_update
  integer :: inner_itr  
  !------------------------------------------------------------------------
  
  !- ATMOS
  if ( this%dyn_proc%IsActivated() ) then
    call PROF_rapstart('GlobalSW_dynamics', 1)  
    tm_process_id = this%dyn_proc%tm_process_id
    is_update = this%time_manager%Do_process( tm_process_id )

    LOG_PROGRESS(*) 'shallow water / dynamics'   
    do inner_itr=1, this%time_manager%Get_process_inner_itr_num( tm_process_id )
      call this%dyn_proc%update( &
        this%mesh, this%vars%PROGVARS_manager, this%vars%AUXVARS_manager, &
        this%vars%PHYTENDS_manager, is_update                             )
    end do
    call PROF_rapend('GlobalSW_dynamics', 1)  
  end if

  !########## Calculate diagnostic variables ##########  
  call this%vars%Clac_diagnostics( this%mesh )
  call this%vars%AUXVARS_manager%MeshFieldComm_Exchange()
    
  !#### Check values #################################
  call this%vars%Check()

  return
end subroutine SW_update

subroutine SW_calc_tendency( this )
  implicit none
  class(GlobalSWComponent), intent(inout) :: this
  !------------------------------------------------------------------------
  
  call PROF_rapstart( 'GlobalSW_tendency', 1)

  !########## calculate tendency ##########

  !* Exchange halo data ( for physics )
  call PROF_rapstart( 'ATM_exchange_prgv', 2)
  call this%vars%PROGVARS_manager%MeshFieldComm_Exchange()
  call PROF_rapend( 'ATM_exchange_prgv', 2)

  call PROF_rapend( 'GlobalSW_tendency', 1)

  return
end subroutine SW_calc_tendency

subroutine SW_finalize( this )
  implicit none
  class(GlobalSWComponent), intent(inout) :: this
  !--------------------------------------------------

  LOG_INFO('GlobalSWComponent_finalize',*)

  if ( .not. this%IsActivated() ) return

  call this%dyn_proc%finalize()
  
  call this%vars%Final()
  call this%mesh%Final()
  call this%time_manager%Final()

  return  
end subroutine SW_finalize

end module mod_globalsw_component