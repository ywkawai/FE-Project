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
    TIME_STARTDAYSEC, &
    TIME_NOWDATE,     &
    TIME_NOWDAY,      &
    TIME_NOWDAYSEC,   &           
    TIME_NOWSEC,      &
    TIME_NOWSTEP,     &
    TIME_NOWSUBSEC,   &
    TIME_NSTEP,       &
    TIME_DTSEC,       &
    TIME_OFFSET_YEAR
  
  use scale_calendar, only: &
    CALENDAR_date2daysec,    &
    CALENDAR_daysec2date,    &
    CALENDAR_combine_daysec, &
    CALENDAR_adjust_daysec,  &
    CALENDAR_unit2sec

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

  
  type :: TIME_manager_process
    real(DP) :: dtsec
    integer :: dstep
    integer :: res_step
    logical :: do_step
    integer :: inner_itr_num
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
  integer, public :: TIME_STARTDATE(6) = (/ 0, 1, 1, 0, 0, 0 /)
  real(DP), public :: TIME_STARTMS     = 0.0_DP !< [millisec]
  integer, public :: TIME_STARTDAY
  real(DP), public :: TIME_STARTSEC

  integer, public  :: TIME_ENDDATE(6)
  real(DP), public :: TIME_ENDMS
  integer, public  :: TIME_ENDDAY
  real(DP), public :: TIME_ENDSEC

  logical, public :: TIME_DOresume
  logical, public :: TIME_DOend

  real(DP), public :: TIME_DTSEC_RESUME

  ! integer, public :: TIME_DTSEC_RESUME
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

contains

  subroutine TIME_manager_Init()
    implicit none
    
    real(DP)               :: TIME_DURATION                = UNDEF8
    character(len=H_SHORT) :: TIME_DURATION_UNIT           = "SEC"
    real(DP)               :: TIME_DT                      = UNDEF8
    character(len=H_SHORT) :: TIME_DT_UNIT                 = "SEC"
    real(DP) :: TIME_DURATIONSEC

    real(DP)               :: TIME_DT_RESUME               = UNDEF8
    character(len=H_SHORT) :: TIME_DT_RESUME_UNIT          = ""

    namelist /PARAM_TIME/ &
       TIME_STARTDATE,               &
       TIME_STARTMS,                 &
       TIME_DURATION,                &
       TIME_DURATION_UNIT,           &
       TIME_DT,                      &
       TIME_DT_UNIT,                 &
       TIME_DT_RESUME,               &
       TIME_DT_RESUME_UNIT       

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
        
    !--
    TIME_OFFSET_YEAR = TIME_STARTDATE(1)
    
    call CALENDAR_date2daysec( TIME_STARTDAY, TIME_STARTSEC,                       & ! [OUT]
                              TIME_STARTDATE(:), TIME_STARTMS,  TIME_OFFSET_YEAR   ) ! [IN]

    TIME_STARTDAYSEC  = CALENDAR_combine_daysec( TIME_STARTDAY, TIME_STARTSEC )
    TIME_NOWDATE(:)   = TIME_STARTDATE(:)
    TIME_NOWSUBSEC    = TIME_STARTMS

    call CALENDAR_unit2sec( TIME_DURATIONSEC, TIME_DURATION, TIME_DURATION_UNIT )
    TIME_ENDSEC = TIME_STARTSEC + TIME_DURATIONSEC

    call CALENDAR_unit2sec( TIME_DTSEC, TIME_DT, TIME_DT_UNIT )
    TIME_NSTEP   = int( TIME_DURATIONSEC / TIME_DTSEC )
    TIME_NOWSTEP = 1

    call CALENDAR_unit2sec( TIME_DTSEC_RESUME, TIME_DT_RESUME, TIME_DT_RESUME_UNIT )
    TIME_DSTEP_RESUME = nint( TIME_DTSEC_RESUME / TIME_DTSEC )
    TIME_RES_RESUME = TIME_DSTEP_RESUME - 1

    !--
    TIME_MANAGER_COMPONENT_num = 0
    do n=1, TIME_MANAGER_COMPONENT_MAX_NUM
      nullify( time_manager_comp_ptr_list(n)%ptr )
    end do
    
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

    character(len=27) :: nowchardate
    integer :: n
    type(TIME_manager_component), pointer :: tm_comp

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

    if (mod(TIME_NOWSTEP-1, 1000) == 0) then
      call CALENDAR_date2char( nowchardate,                       & ! [OUT]
                              TIME_NOWDATE(:), TIME_NOWSUBSEC     ) ! [IN]    
      LOG_PROGRESS('(1x,2A,2(A,I7),A,F10.1)') 'TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP
    end if

    return
  end subroutine TIME_manager_checkstate

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
    !------------------------------------------------------------------------

    !--

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
      this%dtsec_restart = TIME_ENDSEC - TIME_STARTSEC
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
    this%res_step = 0
    this%res_step_restart = 0

    this%process_num = 0

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

    if (dt == UNDEF8) then
      LOG_INFO_CONT(*) 'Not found TIME_DT_'//trim(process_name)//'. TIME_DTSEC is used.'
      this%dtsec = TIME_DTSEC
    else 
      call CALENDAR_unit2sec( this%dtsec, dt, dt_unit )
    end if

    this%inner_itr_num = max( nint(TIME_DTSEC/this%dtsec), 1 )
    this%dstep = nint( this%dtsec / TIME_DTSEC * this%inner_itr_num )

    if ( abs(  real(this%inner_itr_num, kind=DP)*this%dtsec         &
             - real(this%dstep         ,kind=DP)*TIME_DTSEC          ) > eps ) then
      LOG_ERROR("TIME_manager_process_Init",*) 'delta t('//trim(process_name)//') must be a multiple of delta t ', &
        this%dtsec, real(this%dstep,kind=DP)*TIME_DTSEC
      call PRC_abort
    end if

    this%res_step = 0

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

end module scale_time_manager