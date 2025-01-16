#include "scalelib.h"
program test_time_manager
  use scale_precision
  use scale_prc
  use scale_io
  use scale_time_manager
  use scale_calendar, only: CALENDAR_date2char    

  implicit none

  type(TIME_manager_component) :: tm_atm
  integer :: TM_ATMOS_DYN_ID
  integer :: TM_ATMOS_PHYS_CP_ID
  integer :: TM_ATMOS_PHYS_RD_ID
  type(TIME_manager_component) :: tm_ocn
  character(len=27) :: nowchardate

  !-------------------------------

  write(*,*) 'Start test_time_manager..'
  call init()

  !****
  do
    call TIME_manager_checkstate()
    if (TIME_DOresume) then
      call restart_read()
    end if

    call TIME_manager_advance()
    call CALENDAR_date2char( nowchardate, & ! [OUT]
      TIME_NOWDATE(:), TIME_NOWSUBSEC     ) ! [IN]    
    LOG_PROGRESS('(1x,2A,2(A,I7),A,F10.1)') 'TIME: ', nowchardate,' STEP:',TIME_NOWSTEP, '/', TIME_NSTEP


    if ( tm_atm%do_step ) call atmos_update()
    if ( tm_ocn%do_step ) call ocn_update()

    call restart_write()

    if (TIME_DOend) exit
    if( IO_L ) call flush(IO_FID_LOG)
  end do
  !****

  call final()
  write(*,*) 'test_time_manager has been succeeded!'

contains  
  subroutine init()
    use scale_const, only: &
      UNDEF8 => CONST_UNDEF8
    use scale_calendar, only: CALENDAR_setup
    implicit none

    integer :: comm, myrank, nprocs
    logical :: ismaster

    real(DP) :: TIME_DT
    real(DP) :: TIME_DT_RESTART
    character(len=H_SHORT) :: TIME_DT_UNIT
    character(len=H_SHORT) :: TIME_DT_RESTART_UNIT

    namelist /PARAM_ATMOS/ &
      TIME_DT, TIME_DT_UNIT,                 &
      TIME_DT_RESTART, TIME_DT_RESTART_UNIT 
    namelist /PARAM_ATMOS_DYN/     &
      TIME_DT, TIME_DT_UNIT 
    namelist /PARAM_ATMOS_PHYS_CP/ &
      TIME_DT, TIME_DT_UNIT      
    namelist /PARAM_ATMOS_PHYS_RD/ &
      TIME_DT, TIME_DT_UNIT      
    
    namelist /PARAM_OCN/ &
      TIME_DT, TIME_DT_UNIT,                 &
      TIME_DT_RESTART, TIME_DT_RESTART_UNIT 
    namelist /PARAM_OCN_DYN/ &
      TIME_DT, TIME_DT_UNIT
    integer :: ierr

    character(len=H_LONG) :: cnf_fname              ! config file for launcher
    !-------------------------------------------------

    call PRC_MPIstart( comm )
    
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]

    ! setup scale_io
    cnf_fname = IO_ARG_getfname( ismaster )
    call IO_setup( "test_time_manager", cnf_fname, allow_noconf = .false. )

    ! setup log
    call IO_LOG_setup( myrank, ismaster )   

    ! setup calendar & initial time
    call CALENDAR_setup

    call TIME_manager_Init

    ! setup ATMOS ------------------

    TIME_DT = UNDEF8
    TIME_DT_UNIT = 'SEC'
    TIME_DT_RESTART = UNDEF8
    TIME_DT_RESTART_UNIT = 'SEC'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_ATMOS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS)
    call tm_atm%Init( 'ATMOS', TIME_DT, TIME_DT_UNIT, TIME_DT_RESTART, TIME_DT_RESTART_UNIT )
    call TIME_manager_Regist_component( tm_atm )

    ! setup ATMOS_DYN ------------------

    TIME_DT = UNDEF8
    TIME_DT_UNIT = 'SEC'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN)
    call tm_atm%Regist_process( 'ATMOS_DYN', TIME_DT, TIME_DT_UNIT, & ! (in)
                                TM_ATMOS_DYN_ID                     ) ! (out)


    ! setup ATMOS_PHYS_CP ------------------

    TIME_DT = UNDEF8
    TIME_DT_UNIT = 'SEC'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHYS_CP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_ATMOS_PHYS_CP. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHYS_CP)
    call tm_atm%Regist_process( 'ATMOS_PHYS_CP', TIME_DT, TIME_DT_UNIT, & ! (in)
                                TM_ATMOS_PHYS_CP_ID                     ) ! (out)

    ! setup ATMOS_PHYS_CP ------------------

    TIME_DT = UNDEF8
    TIME_DT_UNIT = 'SEC'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHYS_RD,iostat=ierr)
    if( ierr < 0 ) then !--- missing
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_ATMOS_PHYS_RD. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHYS_RD)
    call tm_atm%Regist_process( 'ATMOS_PHYS_RD', TIME_DT, TIME_DT_UNIT, & ! (in)
                                TM_ATMOS_PHYS_RD_ID                     ) ! (out)
                                                            
    ! setup OCN ------------------

    TIME_DT = UNDEF8
    TIME_DT_UNIT = 'SEC'
    TIME_DT_RESTART = UNDEF8
    TIME_DT_RESTART_UNIT = 'SEC'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_ATMOS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCN)
    call tm_ocn%Init( 'OCN', TIME_DT, TIME_DT_UNIT, TIME_DT_RESTART, TIME_DT_RESTART_UNIT )
    call TIME_manager_Regist_component( tm_ocn )

    return
  end subroutine init

  subroutine atmos_update()
    implicit none

    integer :: itr
    !---------------------

    LOG_INFO("amots_update",*) ''

    if ( tm_atm%Do_process(TM_ATMOS_DYN_ID) ) then
      do itr=1, tm_atm%process_list(TM_ATMOS_DYN_ID)%inner_itr_num
        LOG_INFO("amots_update_dyn",*) "itr=", itr
      end do
    end if

    if ( tm_atm%Do_process(TM_ATMOS_PHYS_CP_ID) ) then
      do itr=1, tm_atm%process_list(TM_ATMOS_PHYS_CP_ID)%inner_itr_num
        LOG_INFO("amots_update_phys_cp",*) "itr=", itr
      end do
    end if

    if ( tm_atm%Do_process(TM_ATMOS_PHYS_RD_ID) ) then
      do itr=1, tm_atm%process_list(TM_ATMOS_PHYS_RD_ID)%inner_itr_num
        LOG_INFO("amots_update_phys_rd",*) "itr=", itr
      end do
    end if

    return
  end subroutine atmos_update

  subroutine ocn_update()
    implicit none

    LOG_INFO("ocn_update",*) ''
    return
  end subroutine ocn_update

  subroutine restart_read()
    implicit none

    LOG_INFO("resume",*) 'Read restart file..'

    return
  end subroutine restart_read

  subroutine restart_write()
    implicit none

    if( tm_atm%do_restart ) then
      LOG_INFO("atm",*) 'Write restart file..'      
    end if
    if( tm_ocn%do_restart ) then
      LOG_INFO("ocn",*) 'Write restart file..'      
    end if

    return
  end subroutine restart_write
  
  subroutine final()
    implicit none
    !-------------------------------------------------

    call TIME_manager_Final()
    call tm_atm%Final()
    call tm_ocn%Final()

    call PRC_MPIfinish()

    return
  end subroutine final
end program test_time_manager