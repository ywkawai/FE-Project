!-------------------------------------------------------------------------------
!> Module common / time
!!
!! @par Description
!!          Module to manage time with temporal integration
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_time_manager
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_const, only: &
    UNDEF8 => CONST_UNDEF8  
  use scale_precision
  use scale_prc
  use scale_io
  
  use scale_time, only: &
    TIME_STARTDAYSEC,           &
    TIME_NOWDATE,               &
    TIME_NOWDAY,                &
    TIME_NOWDAYSEC,             &               
    TIME_NOWSEC,                &
    TIME_NOWSUBSEC,             &
    TIME_NOWSTEP,               &   
    TIME_NSTEP,                 &
    TIME_DTSEC,                 &
    TIME_OFFSET_YEAR,           &
    TIME_DTSEC_WALLCLOCK_CHECK, &
    TIME_DSTEP_WALLCLOCK_CHECK

  
  use scale_calendar, only: &
    CALENDAR_date2daysec,    &
    CALENDAR_daysec2date,    &
    CALENDAR_combine_daysec, &
    CALENDAR_adjust_daysec,  &
    CALENDAR_unit2sec,       &
    CALENDAR_CFunits2sec,    &
    CALENDAR_date2char

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !  
  public :: TIME_manager_Init
  public :: TIME_manager_Final
  public :: TIME_manager_Regist_component
  public :: TIME_manager_checkstate
  public :: TIME_manager_advance
  public :: TIME_manager_report_timeintervals

  
  type :: TIME_manager_process
    real(DP) :: dtsec
    integer :: dstep
    integer :: res_step
    logical :: do_step
    integer :: inner_itr_num
    character(len=H_SHORT) :: process_name
  contains
    procedure, public :: Init => TIME_manager_process_Init
    procedure, public :: Check_state => TIME_manager_process_checkstate
    procedure, public :: Final => TIME_manager_process_Final
  end type TIME_manager_process

  integer, private, parameter :: TIME_MANAGER_PROCESS_MAX_NUM = 10
  type, public :: TIME_manager_component
    real(DP) :: dtsec
    integer :: dstep
    integer :: res_step
    logical :: do_step
    !-
    real(DP) :: dtsec_restart
    integer :: dstep_restart
    integer :: res_step_restart
    logical :: do_restart
    character(len=H_SHORT) :: comp_name
    !--
    type(TIME_manager_process) :: process_list(TIME_MANAGER_PROCESS_MAX_NUM)
    integer :: process_num
  contains
    procedure, public :: Init => TIME_manager_component_Init
    procedure, public :: Regist_process => TIME_manager_component_Regist_process
    procedure, public :: Check_state => TIME_manager_component_checkstate
    procedure, public :: Do_process => TIME_manager_component_do_process
    procedure, public :: Get_process_inner_itr_num => TIME_manager_component_get_process_inner_itr_num
    procedure, public :: Final => TIME_manager_component_Final
  end type TIME_manager_component

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: TIME_STARTDATE(6) = (/ -999, 1, 1, 0, 0, 0 /)
  real(DP), public :: TIME_STARTMS     = 0.0_DP !< [millisec]
  real(DP), public :: TIME_STARTSUBSEC
  integer, public :: TIME_STARTDAY
  real(DP), public :: TIME_STARTSEC

  integer, public  :: TIME_ENDDATE(6)
  real(DP), public :: TIME_ENDSUBSEC
  integer, public  :: TIME_ENDDAY
  real(DP), public :: TIME_ENDSEC

  real(DP), public :: TIME_DURATIONSEC

  logical, public :: TIME_DOresume
  logical, public :: TIME_DOend

  real(DP), public :: TIME_DTSEC_RESUME
  integer, public :: TIME_DSTEP_RESUME

  public :: TIME_STARTDAYSEC
  public :: TIME_NOWDATE,  TIME_NOWSUBSEC, TIME_NOWDAY,  TIME_NOWDAYSEC, TIME_NOWSEC
  public :: TIME_NOWSTEP, TIME_NSTEP
  public :: TIME_DTSEC
  public :: TIME_OFFSET_YEAR

  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  
  !-----------------------------------------------------------------------------

  integer, private :: TIME_RES_RESUME
  
  real(DP), private :: TIME_WALLCLOCK_START             ! Start time of wall clock             [sec]
  real(DP), private :: TIME_WALLCLOCK_LIMIT   = -1.0_DP ! Elapse time limit of wall clock time [sec]
  real(DP), private :: TIME_WALLCLOCK_SAFE    =  0.9_DP ! Safety coefficient for elapse time limit
  real(DP), private :: TIME_WALLCLOCK_safelim           ! TIME_WALLCLOCK_LIMIT * TIME_WALLCLOCK_SAFE

  real(DP), private, parameter :: eps = 1.E-6_DP !> epsilon for timesec


  integer, private, parameter :: TIME_MANAGER_COMPONENT_MAX_NUM = 20
  type TIME_manager_component_ptr
    type(TIME_manager_component), pointer :: ptr
  end type TIME_manager_component_ptr
  type(TIME_manager_component_ptr), private :: time_manager_comp_ptr_list(TIME_MANAGER_COMPONENT_MAX_NUM)
  integer, private :: TIME_MANAGER_COMPONENT_num

  logical :: setup_tinteg_flag

contains

  subroutine TIME_manager_Init( &
      setup_TimeIntegration,    &
      restart_in_basename       )

    use scale_file, only: &
      FILE_Get_Attribute

    implicit none

    logical, intent(in), optional :: setup_TimeIntegration
    character(len=*), intent(in), optional :: restart_in_basename
    
    real(DP)               :: TIME_DURATION                = UNDEF8
    character(len=H_SHORT) :: TIME_DURATION_UNIT           = "SEC"
    real(DP)               :: TIME_DT                      = UNDEF8
    character(len=H_SHORT) :: TIME_DT_UNIT                 = "SEC"

    real(DP)               :: TIME_DT_RESUME               = UNDEF8
    character(len=H_SHORT) :: TIME_DT_RESUME_UNIT          = ""

    real(DP)               :: TIME_DT_WALLCLOCK_CHECK      = UNDEF8
    character(len=H_SHORT) :: TIME_DT_WALLCLOCK_CHECK_UNIT = ""

    namelist /PARAM_TIME/ &
       TIME_STARTDATE,               &
       TIME_STARTMS,                 &
       TIME_DURATION,                &
       TIME_DURATION_UNIT,           &
       TIME_DT,                      &
       TIME_DT_UNIT,                 &
       TIME_DT_RESUME,               &
       TIME_DT_RESUME_UNIT,          &
       TIME_DT_WALLCLOCK_CHECK,      &
       TIME_DT_WALLCLOCK_CHECK_UNIT, &
       TIME_WALLCLOCK_LIMIT,         &
       TIME_WALLCLOCK_SAFE

    integer              :: dateday
    real(DP)             :: datesec
    real(DP)             :: cftime(1)
    character(len=H_MID) :: cfunits
   
    character(len=27) :: startchardate
    character(len=27) :: endchardate
   
    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("TIME_manager_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TIME,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("TIME_manager_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("TIME_manager_setup",*) 'Not appropriate names in namelist PARAM_TIME. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_TIME)

    if ( present(setup_TimeIntegration) ) then
      setup_tinteg_flag = setup_TimeIntegration
    else
      setup_tinteg_flag = .true.
    end if

    if ( setup_tinteg_flag ) then
      if ( TIME_DT == UNDEF8 ) then
        LOG_ERROR("TIME_manager_setup",*) 'Not found TIME_DT. STOP.'
        call PRC_abort
      endif
      if ( TIME_DURATION == UNDEF8 ) then
        LOG_ERROR("TIME_manager_setup",*) 'Not found TIME_DURATION. STOP.'
        call PRC_abort
      endif
            
      if ( TIME_DT_RESUME == UNDEF8 ) then
        TIME_DT_RESUME = TIME_DURATION
      endif
      if ( TIME_DT_RESUME_UNIT == '' ) then
          LOG_INFO_CONT(*) 'Not found TIME_DT_RESUME_UNIT.        TIME_DURATION_UNIT is used.'
          TIME_DT_RESUME_UNIT = TIME_DURATION_UNIT
      endif
    end if

    !--- calculate time
    if ( TIME_STARTDATE(1) == -999 ) then
      if ( present(restart_in_basename) ) then
        if ( restart_in_basename /= '' ) then ! read start time from the restart data
          call FILE_Get_Attribute( RESTART_IN_BASENAME, & ! [IN]
                                   "global",            & ! [IN]
                                   'time_start',        & ! [IN]
                                   cftime(:),           & ! [OUT]
                                   rankid = PRC_myrank, & ! [IN]
                                   single = .false.     ) ! [IN]

          call FILE_Get_Attribute( RESTART_IN_BASENAME, & ! [IN]
                                   "global",            & ! [IN]
                                   'time_units',        & ! [IN]
                                   cfunits,             & ! [OUT]
                                   rankid = PRC_myrank, & ! [IN]
                                   single = .false.     ) ! [IN]

         dateday = 0
         datesec = CALENDAR_CFunits2sec( cftime(1), cfunits, 0 )

         call CALENDAR_adjust_daysec( dateday, datesec )

         call CALENDAR_daysec2date( TIME_STARTDATE,   & ! [OUT]
                                    TIME_STARTSUBSEC, & ! [OUT]
                                    dateday,          & ! [IN]
                                    datesec,          & ! [IN]
                                    0                 ) ! [IN]
        end if
      else
         TIME_STARTDATE   = (/ 0, 1, 1, 0, 0, 0 /)
         TIME_STARTSUBSEC = 0.0_DP
      endif
    else
      TIME_STARTSUBSEC = TIME_STARTMS * 1.E-3_DP
    endif

    TIME_OFFSET_YEAR = TIME_STARTDATE(1)
    
    call CALENDAR_date2daysec( TIME_STARTDAY, TIME_STARTSEC,                        & ! [OUT]
                               TIME_STARTDATE(:), TIME_STARTMS,  TIME_OFFSET_YEAR   ) ! [IN]

    call CALENDAR_date2char( startchardate,                        & ! [OUT]
                             TIME_STARTDATE(:), TIME_STARTSUBSEC   ) ! [IN]
  
    TIME_STARTDAYSEC  = CALENDAR_combine_daysec( TIME_STARTDAY, TIME_STARTSEC )

    TIME_NOWDATE(:)   = TIME_STARTDATE(:)
    TIME_NOWSUBSEC    = TIME_STARTMS
    TIME_NOWDAY       = TIME_STARTDAY
    TIME_NOWSEC       = TIME_STARTSEC
    TIME_NOWDAYSEC    = CALENDAR_combine_daysec( TIME_NOWDAY, TIME_NOWSEC )

    TIME_ENDDAY       = TIME_STARTDAY
    if (setup_tinteg_flag) then
      call CALENDAR_unit2sec( TIME_DURATIONSEC, TIME_DURATION, TIME_DURATION_UNIT )
      TIME_ENDSEC = TIME_STARTSEC + TIME_DURATIONSEC
    else
      TIME_ENDSEC = TIME_STARTSEC
    end if

    call CALENDAR_adjust_daysec( TIME_ENDDAY, TIME_ENDSEC ) ! [INOUT]

    call CALENDAR_daysec2date( TIME_ENDDATE(:), & ! [OUT]
                               TIME_ENDSUBSEC,  & ! [OUT]
                               TIME_ENDDAY,     & ! [IN]
                               TIME_ENDSEC,     & ! [IN]
                               TIME_OFFSET_YEAR ) ! [IN]

    call CALENDAR_date2char( endchardate,     & ! [OUT]
                             TIME_ENDDATE(:), & ! [IN]
                             TIME_ENDSUBSEC   ) ! [IN]

    LOG_NEWLINE
    LOG_INFO("TIME_manager_setup",*) 'Global date / time setting '
    LOG_INFO_CONT('(1x,A,A)') 'START Date     : ', startchardate
    LOG_INFO_CONT('(1x,A,A)') 'END   Date     : ', endchardate

    if (setup_tinteg_flag) then
      call CALENDAR_unit2sec( TIME_DTSEC, TIME_DT, TIME_DT_UNIT )
      TIME_NSTEP   = int( TIME_DURATIONSEC / TIME_DTSEC )
      TIME_NOWSTEP = 1

      call CALENDAR_unit2sec( TIME_DTSEC_RESUME, TIME_DT_RESUME, TIME_DT_RESUME_UNIT )
      TIME_DSTEP_RESUME = nint( TIME_DTSEC_RESUME / TIME_DTSEC )
      TIME_RES_RESUME = TIME_DSTEP_RESUME - 1
    else
      TIME_DTSEC = 1.0_RP
    end if

    !--
    TIME_MANAGER_COMPONENT_num = 0
    do n=1, TIME_MANAGER_COMPONENT_MAX_NUM
      nullify( time_manager_comp_ptr_list(n)%ptr )
    end do
    
    !--
    TIME_WALLCLOCK_START = PRC_MPItime()

    if ( TIME_WALLCLOCK_LIMIT > 0.0_DP ) then
      LOG_NEWLINE
      LOG_INFO("TIME_manager_setup",*) 'Wall clock time limit of execution is specified.'

      if ( TIME_DT_WALLCLOCK_CHECK == UNDEF8 ) then
        LOG_INFO_CONT(*) 'Not found TIME_DT_WALLCLOCK_CHECK. TIME_DT is used.'
        TIME_DTSEC_WALLCLOCK_CHECK = TIME_DTSEC
     else
        if ( TIME_DT_WALLCLOCK_CHECK_UNIT == '' ) then
           LOG_INFO_CONT(*) 'Not found TIME_DT_WALLCLOCK_CHECK_UNIT. TIME_DURATION_UNIT is used.'
           TIME_DT_WALLCLOCK_CHECK_UNIT = TIME_DURATION_UNIT
        endif
        call CALENDAR_unit2sec( TIME_DTSEC_WALLCLOCK_CHECK, TIME_DT_WALLCLOCK_CHECK, TIME_DT_WALLCLOCK_CHECK_UNIT )
        TIME_DTSEC_WALLCLOCK_CHECK = max( TIME_DTSEC_WALLCLOCK_CHECK, TIME_DTSEC )
     endif

     TIME_DSTEP_WALLCLOCK_CHECK = int( TIME_DTSEC_WALLCLOCK_CHECK / TIME_DTSEC )

     TIME_WALLCLOCK_SAFE    = max( min( TIME_WALLCLOCK_SAFE, 1.0_DP ), 0.0_DP )
     TIME_WALLCLOCK_safelim = TIME_WALLCLOCK_LIMIT * TIME_WALLCLOCK_SAFE

     LOG_INFO_CONT('(1x,A,F10.1,A)')      'This job stops after ', TIME_WALLCLOCK_safelim, ' seconds.'
     LOG_INFO_CONT('(1x,A,F10.3,A,I8,A)') 'Time interval for check     : ', TIME_DTSEC_WALLCLOCK_CHECK, &
                                          ' (step interval=', TIME_DSTEP_WALLCLOCK_CHECK, ')'
    end if

    return
  end subroutine TIME_manager_Init

  subroutine TIME_manager_Final()
    implicit none

    integer :: n
    !---------------------------------------------------------------------------

    do n=1, TIME_MANAGER_COMPONENT_num
      nullify( time_manager_comp_ptr_list(n)%ptr )
    end do
    TIME_MANAGER_COMPONENT_num = 0

    return
  end subroutine TIME_manager_Final  

  subroutine TIME_manager_advance()
    implicit none
    !---------------------------------------------------------------------------

    TIME_NOWSTEP = TIME_NOWSTEP + 1
    TIME_NOWDAY  = TIME_STARTDAY
    TIME_NOWSEC  = TIME_STARTSEC + real(TIME_NOWSTEP-1,kind=DP) * TIME_DTSEC
    
    ! reallocate day & sub-day
    call CALENDAR_adjust_daysec( TIME_NOWDAY, TIME_NOWSEC ) ! [INOUT]
  
    call CALENDAR_daysec2date( TIME_NOWDATE(:), TIME_NOWSUBSEC,           & ! [OUT]
                               TIME_NOWDAY, TIME_NOWSEC, TIME_OFFSET_YEAR ) ! [IN]
  
    TIME_NOWDAYSEC = CALENDAR_combine_daysec( TIME_NOWDAY, TIME_NOWSEC )
  
    if (TIME_NOWSTEP > TIME_NSTEP) then
      TIME_DOend = .true.
    end if

    return
  end subroutine TIME_manager_advance
  
  subroutine TIME_manager_checkstate()
    use scale_calendar, only: CALENDAR_date2char    
    implicit none

    integer :: n
    type(TIME_manager_component), pointer :: tm_comp

    real(DP)          :: WALLCLOCK_elapse
    character(len=27) :: nowchardate
    logical           :: TO_STDOUT
    !---------------------------------------------------------------------------

    do n=1, TIME_MANAGER_COMPONENT_num
      tm_comp => time_manager_comp_ptr_list(n)%ptr
      call tm_comp%Check_state()
    end do
 
    TIME_DOresume   = .false.
    TIME_RES_RESUME = TIME_RES_RESUME + 1

    if (TIME_RES_RESUME == TIME_DSTEP_RESUME) then
      TIME_DOresume   = .true.
      TIME_RES_RESUME =  0
    end if

    !----

    TO_STDOUT = .false.
    if ( IO_STEP_TO_STDOUT > 0 ) then
       if( mod(TIME_NOWSTEP-1,IO_STEP_TO_STDOUT) == 0 ) TO_STDOUT = .true.
    endif

    call CALENDAR_date2char( nowchardate,     & ! [OUT]
                             TIME_NOWDATE(:), & ! [IN]
                             TIME_NOWSUBSEC   ) ! [IN]

    WALLCLOCK_elapse = PRC_MPItime() - TIME_WALLCLOCK_START

    LOG_NEWLINE
    if ( TIME_WALLCLOCK_LIMIT > 0.0_DP ) then
       LOG_PROGRESS('(1x,2A,2(A,I7),2(A,F10.1))') 'TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP, &
                                                                 ' WCLOCK:', WALLCLOCK_elapse, '/', TIME_WALLCLOCK_safelim
       if ( PRC_UNIVERSAL_IsMaster .AND. TO_STDOUT ) then ! universal master node
          write(*,'(1x,2A,2(A,I7),2(A,F10.1))') 'TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP, &
                                                ' WCLOCK:', WALLCLOCK_elapse, '/', TIME_WALLCLOCK_safelim
       endif
    else
       LOG_PROGRESS('(1x,2A,2(A,I7),A,F10.1)') 'TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP, &
                                                              ' WCLOCK:', WALLCLOCK_elapse
       if ( PRC_UNIVERSAL_IsMaster .AND. TO_STDOUT ) then ! universal master node
          write(*,'(1x,2A,2(A,I7),A,F10.1)') 'TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP, &
                                             ' WCLOCK:', WALLCLOCK_elapse
       endif
    endif

    return
  end subroutine TIME_manager_checkstate

  subroutine TIME_manager_report_timeintervals()    
    implicit none

    integer :: n, p
    type(TIME_manager_component), pointer :: tm_comp
    type(TIME_manager_process), pointer :: tm_process
    !-------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("TIME_manager_report_timeintervals",*) 'Time interval for each component (sec.)'

    do n=1, TIME_MANAGER_COMPONENT_num
      tm_comp => time_manager_comp_ptr_list(n)%ptr
      LOG_INFO_CONT(*) trim(tm_comp%comp_name)
      do p=1, tm_comp%process_num
        tm_process => tm_comp%process_list(p)
        LOG_INFO_CONT('(1x,A,F10.3)')        trim(tm_process%process_name)//' (time)             : ', tm_process%dtsec
        if (tm_process%inner_itr_num > 1) then
          LOG_INFO_CONT('(1x,A,I10)') ' (step)             : ', tm_process%inner_itr_num 
        end if
        LOG_INFO_CONT('(1x,A,I8,A)')  ' (step interval=', tm_process%dstep,    ')' 
      end do
    end do

    LOG_NEWLINE
    LOG_INFO_CONT(*)                     'Time interval for restart (sec.)'

    do n=1, TIME_MANAGER_COMPONENT_num
      tm_comp => time_manager_comp_ptr_list(n)%ptr
      LOG_INFO_CONT('(1x,A,F10.3,A,I8,A)') trim(tm_comp%comp_name)//' variables       : ', tm_comp%dtsec_restart, &
                                         ' (step interval=', tm_comp%dstep_restart, ')'
    end do
    LOG_INFO_CONT('(1x,A,F10.3,A,I8,A)') 'Resume                      : ', TIME_DTSEC_RESUME, &
    ' (step interval=', TIME_DSTEP_RESUME,        ')'

    return
  end subroutine TIME_manager_report_timeintervals

  subroutine TIME_manager_Regist_component( tmanager_comp )
    implicit none

    type(TIME_manager_component), intent(in), target :: tmanager_comp
    !------------------------------------------------------------------------

    TIME_MANAGER_COMPONENT_num = TIME_MANAGER_COMPONENT_num + 1
    if (TIME_MANAGER_COMPONENT_num > TIME_MANAGER_COMPONENT_MAX_NUM) then
      LOG_ERROR("TIME_manager_regist_component",*) 'The number of TIME_manager_component registered exceeds ', &
        TIME_MANAGER_COMPONENT_num, TIME_MANAGER_COMPONENT_MAX_NUM
      call PRC_abort
    end if

    time_manager_comp_ptr_list(TIME_MANAGER_COMPONENT_num)%ptr => tmanager_comp
    return
  end subroutine TIME_manager_Regist_component

!---
  subroutine TIME_manager_component_Init( this, comp_name, &
    dt, dt_unit, dt_restart, dt_restart_unit )
    implicit none

    class(TIME_manager_component), intent(inout) :: this
    character(*), intent(in) :: comp_name
    real(DP), intent(in) :: dt
    character(*), intent(in) :: dt_unit
    real(DP), intent(in) :: dt_restart
    character(*), intent(in) :: dt_restart_unit

    real(RP) :: start_sec, end_sec
    real(RP) :: absday, absdaysec, abssec
    !------------------------------------------------------------------------

    this%process_num = 0
    this%res_step = 0
    this%res_step_restart = 0
    this%comp_name = comp_name

    !--
    if (.not. setup_tinteg_flag) return

    if (dt == UNDEF8) then
      LOG_INFO_CONT(*) 'Not found TIME_DT_'//trim(comp_name)//'. TIME_DTSEC is used.'
      this%dtsec = TIME_DTSEC
    else 
      call CALENDAR_unit2sec( this%dtsec, dt, dt_unit )
    end if
    
    this%dstep = nint( this%dtsec / TIME_DTSEC )

    if ( abs(this%dtsec - real(this%dstep,kind=DP)*TIME_DTSEC) > eps ) then
      LOG_ERROR("TIME_manager_component_Init",*) 'delta t('//trim(comp_name)//') must be a multiple of delta t ', &
        this%dtsec, real(this%dstep,kind=DP)*TIME_DTSEC
      call PRC_abort
    end if

    !--

    if (dt_restart == UNDEF8) then
      LOG_INFO_CONT(*) 'Not found TIME_DT_'//trim(comp_name)//'_RESTART.       TIME_DURATION is used.'
      this%dtsec_restart = TIME_DURATIONSEC
    else 
      call CALENDAR_unit2sec( this%dtsec_restart, dt_restart, dt_restart_unit )
    end if

    this%dstep_restart = nint( this%dtsec_restart / TIME_DTSEC )

    if ( abs(this%dtsec_restart - real(this%dstep_restart,kind=DP)*TIME_DTSEC) > eps ) then
      LOG_ERROR("TIME_manager_component_Init",*) 'delta t('//trim(comp_name)//'_RESTART) must be a multiple of delta t ', &
        this%dtsec_restart, real(this%dstep_restart,kind=DP)*TIME_DTSEC
      call PRC_abort
    end if

    !--

    return
  end subroutine TIME_manager_component_Init

  subroutine TIME_manager_component_checkstate( this )    
    implicit none
    class(TIME_manager_component), intent(inout), target :: this

    integer :: n
    type(TIME_manager_process), pointer :: tm_process
    !--------------------------------------------------

    this%do_step = .false.      

    if (this%process_num > 0) then
      do n=1, this%process_num
        tm_process => this%process_list(n)
        call tm_process%Check_state()
        if (tm_process%do_step) this%do_step = .true.
      end do  
    else
      this%res_step = this%res_step + 1
      if ( this%res_step == this%dstep ) then
        this%do_step = .true.
        this%res_step = 0
      end if  
    end if

    !-
    this%do_restart = .false.

    this%res_step_restart = this%res_step_restart + 1
    if ( this%res_step_restart == this%dstep_restart ) then
      this%do_restart = .true.
      this%res_step_restart = 0
    end if  

    return
  end subroutine TIME_manager_component_checkstate

  subroutine TIME_manager_component_Regist_process( this, &
      process_name, dt, dt_unit,                          &
      tm_process_id )    
    implicit none

    class(TIME_manager_component), intent(inout), target :: this
    character(*), intent(in) :: process_name
    real(DP), intent(in) :: dt
    character(*), intent(in) :: dt_unit
    integer, intent(out) :: tm_process_id

    type(TIME_manager_process), pointer :: tm_process
    !--------------------------------------------------

    this%process_num = this%process_num + 1
    tm_process_id = this%process_num 

    if( this%process_num > TIME_MANAGER_PROCESS_MAX_NUM) then
      LOG_ERROR("TIME_manager_component_Regist_process",*) 'The number of TIME_manager_process registered exceeds ', &
        this%process_num, TIME_MANAGER_PROCESS_MAX_NUM
      call PRC_abort
    end if

    tm_process => this%process_list(this%process_num)
    call tm_process%Init( process_name, dt, dt_unit )

    return
  end subroutine TIME_manager_component_Regist_process

  subroutine TIME_manager_component_Final( this )    
    implicit none
    class(TIME_manager_component), intent(inout) :: this

    integer :: n
    !--------------------------------------------------

    do n=1, this%process_num
      call this%process_list(n)%Final()
    end do
    this%process_num = 0

    return
  end subroutine TIME_manager_component_Final

  function TIME_manager_component_do_process( this, tm_process_id ) result(do_step)
    implicit none
    class(TIME_manager_component), intent(inout) :: this
    integer, intent(in) :: tm_process_id
    logical :: do_step
    !--------------------------------------------------
    do_step = this%process_list(tm_process_id)%do_step
    return
  end function TIME_manager_component_do_process

  function TIME_manager_component_get_process_inner_itr_num( this, tm_process_id ) result(itr_num)
    implicit none
    class(TIME_manager_component), intent(inout) :: this
    integer, intent(in) :: tm_process_id
    integer :: itr_num
    !--------------------------------------------------
    itr_num = this%process_list(tm_process_id)%inner_itr_num
    return
  end function TIME_manager_component_get_process_inner_itr_num

  !---

  subroutine TIME_manager_process_Init( this, process_name, &
      dt, dt_unit )    
    implicit none

    class(TIME_manager_process), intent(inout) :: this
    character(*), intent(in) :: process_name
    real(DP), intent(in) :: dt
    character(*), intent(in) :: dt_unit
    !--------------------------------------------------

    this%res_step = 0
    this%process_name = process_name

    if (.not. setup_tinteg_flag) return

    if (dt == UNDEF8) then
      LOG_INFO_CONT(*) 'Not found TIME_DT_'//trim(process_name)//'. TIME_DTSEC is used.'
      this%dtsec = TIME_DTSEC
    else 
      call CALENDAR_unit2sec( this%dtsec, dt, dt_unit )
    end if

    this%inner_itr_num = max( nint(TIME_DTSEC/this%dtsec), 1 )
    this%dstep = nint( this%dtsec / TIME_DTSEC * dble(this%inner_itr_num) )

    if ( abs(  real(this%inner_itr_num, kind=DP)*this%dtsec         &
             - real(this%dstep         ,kind=DP)*TIME_DTSEC          ) > eps ) then
      LOG_ERROR("TIME_manager_process_Init",*) 'delta t('//trim(process_name)//') must be a multiple of delta t ', &
        this%dtsec, real(this%dstep,kind=DP)*TIME_DTSEC
      call PRC_abort
    end if

    return
  end subroutine TIME_manager_process_Init  


  subroutine TIME_manager_process_Final( this )    
    implicit none
    class(TIME_manager_process), intent(inout) :: this
    !--------------------------------------------------

    return
  end subroutine TIME_manager_process_Final

  subroutine TIME_manager_process_checkstate( this )    
    implicit none
    class(TIME_manager_process), intent(inout) :: this
    !--------------------------------------------------

    this%do_step = .false.
    this%res_step = this%res_step + 1
    if ( this%res_step == this%dstep ) then
      this%do_step = .true.
      this%res_step = 0
    end if

    return
  end subroutine TIME_manager_process_checkstate

!--

end module scale_time_manager