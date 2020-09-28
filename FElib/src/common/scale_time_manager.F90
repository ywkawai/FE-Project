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
  public :: TIME_manager_checkstate
  public :: TIME_manager_advance

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
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  
  !-----------------------------------------------------------------------------

  integer, private :: TIME_RES_RESUME

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
    
    return
  end subroutine TIME_manager_Init
  
  subroutine TIME_manager_Final()
    implicit none
    !---------------------------------------------------------------------------

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
    !---------------------------------------------------------------------------

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

end module scale_time_manager