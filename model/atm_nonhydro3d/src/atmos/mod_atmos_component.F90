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

  use scale_model_component, only: ModelComponent

  use mod_atmos_vars, only: AtmosVars
  use mod_atmos_mesh, only: AtmosMesh
  use mod_atmos_dyn, only: AtmosDyn
  
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

  logical :: ATMOS_DYN_DO  = .true.

  namelist / PARAM_ATMOS / &
    ACTIVATE_FLAG,         &
    TIME_DT,               &
    TIME_DT_UNIT,          &
    TIME_DT_RESTART,       &
    TIME_DT_RESTART_UNIT,  &
    ATMOS_DYN_DO

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

  !- Setup the module for atmosphere / dynamics 
  call this%dyn_proc%ModelComponentProc_Init('AtmosDyn', ATMOS_DYN_DO )
  call this%dyn_proc%setup( this%mesh, this%time_manager )
  
  call PROF_rapend( 'ATM_setup', 1)
  return
end subroutine Atmos_setup

subroutine Atmos_calc_tendency( this )
  implicit none
  class(AtmosComponent), intent(inout) :: this

  !--------------------------------------------------
  call PROF_rapstart( 'ATM_tendency', 1)
  !LOG_INFO('AtmosComponent_calc_tendency',*)

  call this%dyn_proc%calc_tendency( this%mesh, this%vars%PROGVARS_manager, this%vars%AUXVARS_manager )

  call PROF_rapend( 'ATM_tendency', 1)
  return  
end subroutine Atmos_calc_tendency

subroutine Atmos_update( this )

  implicit none
  class(AtmosComponent), intent(inout) :: this
  
  integer :: tm_process_id
  integer :: inner_itr
  !--------------------------------------------------
  call PROF_rapstart( 'ATM_update', 1)

  !########## Dynamics ########## 

  if ( this%dyn_proc%IsActivated() ) then
    tm_process_id = this%dyn_proc%tm_process_id
    if ( this%time_manager%Do_process( tm_process_id ) ) then
      do inner_itr=1, this%time_manager%Get_process_inner_itr_num( tm_process_id )
        call this%dyn_proc%update( this%mesh, this%vars%PROGVARS_manager, this%vars%AUXVARS_manager )
      end do
    end if
  end if
  
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
  call this%vars%Final()
  call this%mesh%Final()
  call this%time_manager%Final()

  return  
end subroutine Atmos_finalize

end module mod_atmos_component