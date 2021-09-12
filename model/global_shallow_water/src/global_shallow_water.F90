!-------------------------------------------------------------------------------
!> Program global shallow water (a launcher of main routine)
!!
!! @par Description
!!         global shallow water model
!!         which is discretized by a discontinuous galerkin method. 
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
program global_shallow_water
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate
  use scale_file_monitor_meshfield, only: &
    FILE_monitor_meshfield_write  
  
  use scale_time_manager, only: &
    TIME_manager_checkstate, TIME_manager_advance,          &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP, TIME_NSTEP, &
    TIME_DOresume, TIME_DOend, TIME_DTSEC

  use scale_file_history_meshfield, only:   &
    FILE_HISTORY_meshfield_write

  use mod_globalsw_component, only: &
    GlobalSWComponent

  use mod_user, only: &
    USER_update, USER_calc_tendency

  implicit none

  integer :: nowstep
  character(len=H_MID) :: timelabel

  type(GlobalSWComponent) :: swmodel
  !-------------------------------------------------------

  call init()
  
  LOG_NEWLINE
  LOG_PROGRESS(*) 'START TIMESTEP'  
  call PROF_setprefx('MAIN')
  call PROF_rapstart('Loop', 0)

  do
    !*******************************************

    ! report current time
    call TIME_manager_checkstate()

    if (TIME_DOresume) then
      ! set state from restart files
      call restart_read()
      call FILE_HISTORY_meshfield_write()
    end if

    !* Advance time *********************************
    call TIME_manager_advance()
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )

    if ( swmodel%IsActivated() .and. swmodel%time_manager%do_step ) then
      call swmodel%update()
    end if

    !- USER
    call USER_update

    !* restart and monitor output *******************
    if ( swmodel%IsActivated() ) call swmodel%vars%Monitor()
    call restart_write
    call FILE_MONITOR_meshfield_write('MAIN', TIME_NOWSTEP)

    !* calc tendencies and diagnostices *************
    if ( swmodel%IsActivated() .and. swmodel%time_manager%do_step ) then
      call swmodel%calc_tendency()
    end if

    !- USER 
    call USER_calc_tendency
  
    !* output history files *************************

    if ( swmodel%IsActivated() ) call swmodel%vars%History()

    call FILE_HISTORY_meshfield_write()

    !*******************************************    
    if (TIME_DOend) exit
    
    if( IO_L ) call flush(IO_FID_LOG)
  end do

  call PROF_rapend('Loop', 0)

  LOG_PROGRESS(*) 'END TIMESTEP'
  LOG_NEWLINE

  !########## Finalize ##########
  call final()

contains
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine restart_read()
    implicit none
    !--------------------------------------------

    call swmodel%vars%Read_restart_file( swmodel%mesh )

    !- History & Monitor 
    call swmodel%vars%History()
    call swmodel%vars%Monitor()

    return
  end subroutine restart_read

  subroutine restart_write
    implicit none    
    !----------------------------------------

    if ( swmodel%isActivated() .and. swmodel%time_manager%do_restart) then
      call swmodel%vars%Write_restart_file()
    end if

    return
  end subroutine restart_write

  subroutine update()
    use mod_user, only: USER_update
    implicit none

    integer :: tm_process_id
    logical :: is_update
    integer :: inner_itr  
    !------------------------------------------------------------------------
    
    !- ATMOS
    if ( swmodel%dyn_proc%IsActivated() ) then
      call PROF_start('GlobalSW_dynamics', 1)  
      tm_process_id = swmodel%dyn_proc%tm_process_id
      is_update = swmodel%time_manager%Do_process( tm_process_id )
  
      LOG_PROGRESS(*) 'shallow water / dynamics'   
      do inner_itr=1, swmodel%time_manager%Get_process_inner_itr_num( tm_process_id )
        call swmodel%dyn_proc%update( &
          swmodel%mesh, swmodel%vars%PROGVARS_manager, swmodel%vars%AUXVARS_manager, &
          swmodel%vars%PHYTENDS_manager, is_update                                   )
      end do
      call PROF_rapend('GlobalSW_dynamics', 1)  
    end if

    !########## Calculate diagnostic variables ##########  
    call swmodel%vars%Clac_diagnostics( swmodel%mesh )
    call swmodel%vars%AUXVARS_manager%MeshFieldComm_Exchange()
      
    !#### Check values #################################
    call swmodel%vars%Check()

    call USER_update()

    return
  end subroutine update

  subroutine calc_tendency()
    use mod_user, only: USER_calc_tendency
    implicit none
    !------------------------------------------------------------------------
    
    call PROF_rapstart( 'GlobalSW_tendency', 1)
  
    !########## calculate tendency ##########
  
    !* Exchange halo data ( for physics )
    call PROF_rapstart( 'ATM_exchange_prgv', 2)
    call swmodel%vars%PROGVARS_manager%MeshFieldComm_Exchange()
    call PROF_rapend( 'ATM_exchange_prgv', 2)
  
    !- ATMOS

    !- USER
    call USER_calc_tendency()

    call PROF_rapend( 'GlobalSW_tendency', 1)

    return
  end subroutine calc_tendency

  subroutine init()
    use scale_const, only: CONST_setup
    use scale_calendar, only: CALENDAR_setup 
    use scale_random, only: RANDOM_setup
    use scale_time_manager, only: TIME_manager_Init
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_setup

    use scale_time_manager, only:           &
      TIME_manager_Init,                    &
      TIME_manager_report_timeintervals
    use scale_meshfield_statistics, only:   &
      MeshField_statistics_setup
    use scale_file_restart_meshfield, only: &
      restart_file,                         &
      FILE_restart_meshfield_setup
    use scale_file_monitor_meshfield, only: &
      FILE_monitor_meshfield_setup  

    
    use mod_user, only: USER_setup
    implicit none

    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr

    character(len=H_MID) :: conf_name
    !------------------------------------------------------------------------

    call PRC_MPIstart( comm )

    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    call get_command_argument(1, conf_name)
    call IO_setup( "global_shallow_water", trim(conf_name) )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
        
    ! setup profiler
    call PROF_setup
    call PROF_setprefx('INIT')
    call PROF_rapstart( "ALL", 0 )

    ! setup constants
    call CONST_setup

    ! setup calendar & initial time
    call CALENDAR_setup

    ! setup random number
    call RANDOM_setup

    ! setup a module for restart file
    call FILE_restart_meshfield_setup
    call TIME_manager_Init( &
      setup_TimeIntegration = .true.,                   &
      restart_in_basename   =  restart_file%in_basename )

    ! setup statistics
    call MeshField_statistics_setup

    ! setup monitor
    call FILE_monitor_meshfield_setup( TIME_DTSEC )

    ! setup submodels
    call  swmodel%setup()

    ! Setup user module
    call USER_setup( swmodel )

    !---
    call PROF_rapend( "ALL", 0 )
    return
  end subroutine init

  subroutine final()
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    use scale_file_monitor_meshfield, only: &
      FILE_monitor_meshfield_final
    use scale_time_manager, only: &
      TIME_manager_Final   
    implicit none
    !------------------------------------------------------------------------

    call PROF_setprefx('FIN')
    call PROF_rapstart( "ALL", 0 )

    call PROF_rapstart('Monit', 2)
    call FILE_monitor_meshfield_final
    call PROF_rapend  ('Monit', 2)

    call PROF_rapstart('File', 2)
    call FILE_HISTORY_meshfield_finalize()
    call PROF_rapend  ('File', 2)

    ! finalization submodels
    call  swmodel%finalize()

    !-
    call TIME_manager_Final

    call PROF_rapend( "ALL", 0 )
    call PROF_rapreport

    call PRC_MPIfinish()

    return
  end subroutine final

end program global_shallow_water
