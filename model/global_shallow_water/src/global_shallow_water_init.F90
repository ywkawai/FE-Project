!-------------------------------------------------------------------------------
!> Program global shallow water init (a launcher of main routine)
!!
!! @par Description
!!         global shallow water model
!!         which is discretized by a discontinuous Galerkin method. 
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
program global_shallow_water_init
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
  
  use scale_time_manager, only: &
    TIME_manager_checkstate, TIME_manager_advance,          &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP, TIME_NSTEP, &
    TIME_DOresume, TIME_DOend

  use scale_file_history_meshfield, only:   &
    FILE_HISTORY_meshfield_write

  use mod_globalsw_component, only: &
    GlobalSWComponent

  use mod_user, only: &
    USER_mkinit

  implicit none

  integer :: nowstep
  character(len=H_MID) :: timelabel

  type(GlobalSWComponent) :: swmodel
  !-------------------------------------------------------

  call init()
  
  LOG_NEWLINE
  LOG_PROGRESS(*) 'START TIMESTEP'  
  call PROF_setprefx('MAIN')

  !- Execute mkinit
  call PROF_rapstart('mkInit',1)
  call USER_mkinit( swmodel )
  call PROF_rapend  ('mkInit',1)

  !- Output
  call PROF_rapstart('restart',1)
  call restart_write()
  call PROF_rapend  ('restart',1)

  LOG_NEWLINE

  !########## Finalize ##########
  call final()

contains
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine restart_write
    implicit none    
    !----------------------------------------

    if ( swmodel%isActivated() ) call swmodel%vars%Write_restart_file()

    return
  end subroutine restart_write

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
    call IO_setup( "global_shallow_water_init", trim(conf_name) )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
        
    ! setup profiler
    call PROF_setup
    call PROF_setprefx('INIT')
    call PROF_rapstart( "All", 0 )

    ! setup constants
    call CONST_setup

    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init( .false. )

    ! setup random number
    call RANDOM_setup

    ! setup a module for restart file
    call FILE_restart_meshfield_setup

    ! setup submodels
    call  swmodel%setup()

    ! Setup user module
    call USER_setup( swmodel )

    !---
    call PROF_rapend( "All", 0 )
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
    call PROF_rapstart('All', 1)

    call PROF_rapstart('File', 2)
    call FILE_HISTORY_meshfield_finalize()
    call PROF_rapend  ('File', 2)

    ! finalization submodels
    call  swmodel%finalize()

    !-
    call TIME_manager_Final

    call PROF_rapend( "All", 1 )
    call PROF_rapreport

    call PRC_MPIfinish()

    return
  end subroutine final

end program global_shallow_water_init
