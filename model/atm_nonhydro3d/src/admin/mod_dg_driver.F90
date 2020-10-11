!-------------------------------------------------------------------------------
!> module SCALE-DG driver
!!
!! @par Description
!!         
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_dg_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc

  use mod_atmos_component, only: &
    AtmosComponent
  
  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_write
  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate

  use mod_user, only: &
    USER_update, USER_calc_tendency

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
#include "scale-dg.h" 

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: dg_driver
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
  character(len=H_MID), private, parameter :: MODELNAME = "SCALE-DG ver. "//VERSION

  type(AtmosComponent) :: atmos

contains
  subroutine dg_driver( &
    comm_world, intercomm_parent, intercomm_child, cnf_fname )

    use scale_time_manager, only: &
      TIME_manager_checkstate, TIME_manager_advance,      &
      TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP, TIME_NSTEP, &
      TIME_DOresume, TIME_DOend

    implicit none

    integer,          intent(in) :: comm_world
    integer,          intent(in) :: intercomm_parent
    integer,          intent(in) :: intercomm_child
    character(len=*), intent(in) :: cnf_fname

    integer :: myrank
    integer :: fpm_counter
    logical :: ismaster
    logical :: sign_exit
    !---------------------------------------------------------------------------

    !########## Initial setup ##########

    ! setup standard I/O
    call IO_setup( MODELNAME, cnf_fname )

    ! setup MPI
    call PRC_LOCAL_setup( comm_world, & ! [IN]
                          myrank,     & ! [OUT]
                          ismaster    ) ! [OUT]

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
 
    call initialize()

    !###########################################################################
 
    !########## main ##########
    LOG_NEWLINE
    LOG_PROGRESS(*) 'START TIMESTEP'
    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_Loop', 0)

    do

      !*******************************************

      ! report current time
      call TIME_manager_checkstate()
  
      if (TIME_DOresume) then
        ! set state from restart file
        call restart_read()
        ! history & monitor file output        
        call FILE_HISTORY_meshfield_write()
      end if

      !* Advance time *********************************

      call TIME_manager_advance()
      call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )

      !* change to next state *************************

      !- ATMOS
      if ( atmos%IsActivated() .and. atmos%time_manager%do_step) then
        call atmos%update()
      end if

      !- USER
      call USER_update()

      !* restart and monitor output *******************
      call restart_write()

      !* calc tendencies and diagnostices *************

      !- ATMOS 
      if ( atmos%IsActivated() .and. atmos%time_manager%do_step) then
        call atmos%calc_tendency()
      end if

      !- USER 
      call USER_calc_tendency()
  
      !* output history files *************************

      if ( atmos%IsActivated() ) call atmos%vars%History()

      call FILE_HISTORY_meshfield_write()
      
      !*******************************************
      if (TIME_DOend) exit
      
      if( IO_L ) call flush(IO_FID_LOG)
    end do

    call PROF_rapend('Main_Loop', 0)

    LOG_PROGRESS(*) 'END TIMESTEP'
    LOG_NEWLINE

    !########## Finalize ##########
    call finalize()    

    return
  end subroutine dg_driver

  !----------------------------
  
  subroutine initialize()

    use scale_const, only: CONST_setup
    use scale_calendar, only: CALENDAR_setup
    use scale_file_restart_meshfield, only: &
      FILE_restart_meshfield_setup
    use scale_time_manager, only: TIME_manager_Init
    use mod_user, only: USER_setup    
    implicit none

    !----------------------------------------------

    ! namelist compatibility check
    !call ADMIN_versioncheck

    ! setup PROF
    call PROF_setup

    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize', 0)
    
    ! setup constants
    call CONST_setup

    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init

    ! setup a module for restart file
    call FILE_restart_meshfield_setup

    ! setup submodels
    call  atmos%setup()

    call USER_setup( atmos )

    call PROF_rapend('Initialize', 0)

    return
  end subroutine initialize

  subroutine finalize()
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_finalize
    use scale_time_manager, only: TIME_manager_Final   
    implicit none
    
    !----------------------------------------------
    call PROF_setprefx('FIN')
    call PROF_rapstart('All', 1)

    !-
    call FILE_HISTORY_meshfield_finalize()

    ! finalization submodels
    call  atmos%finalize()

    !-
    call TIME_manager_Final()

    call PROF_rapend  ('All', 1)
    call PROF_rapreport()

    return
  end subroutine finalize

  subroutine restart_read()
    implicit none    
    !----------------------------------------

    !- read restart data
    if ( atmos%isActivated() ) then
      call atmos%vars%Read_restart_file()
    end if
      
    !- Calculate the tendencies

    if ( atmos%IsActivated() ) call atmos%calc_tendency()


    !- History & Monitor 

    if ( atmos%isActivated() ) then
      call atmos%vars%History()
      ! call atmos%vars%Monitor()
    end if

    return
  end subroutine restart_read

  subroutine restart_write
    implicit none    
    !----------------------------------------


    if ( atmos%isActivated() .and. atmos%time_manager%do_restart) then
      call atmos%vars%Write_restart_file()
    end if

    return
  end subroutine 

end module mod_dg_driver