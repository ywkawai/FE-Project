!-------------------------------------------------------------------------------
!> module SCALE-DG (a main routine of regional/global model)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          SCALE-DG: Regional / global model with atmospheric dynamical core based on DGM
!!
!! @author Yuta Kawai, Team SCALE
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

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_write
  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate
  use scale_file_monitor_meshfield, only: &
    FILE_monitor_meshfield_write
  
  use mod_atmos_component, only: &
    AtmosComponent
  use mod_ocean_component, only: &
    OceanComponent
  use mod_cpl_component, only: &
    CouplerComponent
  use mod_user, only: &
    User

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
  type(OceanComponent) :: ocean
  type(CouplerComponent) :: coupler
  type(User) :: user_

contains
!OCL SERIAL
  subroutine dg_driver(     &
    comm_world, cnf_fname,  &
    path, add_path          )

    use scale_time_manager, only: &
      TIME_manager_checkstate, TIME_manager_advance,          &
      TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP, TIME_NSTEP, &
      TIME_DOresume, TIME_DOend

    implicit none

    integer,          intent(in) :: comm_world
    character(len=*), intent(in) :: cnf_fname
    character(len=*), intent(in) :: path
    logical,          intent(in) :: add_path

    integer :: myrank
    integer :: fpm_counter
    logical :: ismaster
    logical :: sign_exit
    !---------------------------------------------------------------------------

    !########## Initial setup ##########

    ! setup standard I/O
    if ( add_path .and. path /= "" ) then
      call IO_setup( MODELNAME, trim(path)//cnf_fname, prefix=path )
    else
      call IO_setup( MODELNAME, trim(path)//cnf_fname )
    end if

    ! setup MPI
    call PRC_LOCAL_setup( comm_world, & ! [IN]
                          myrank,     & ! [OUT]
                          ismaster    ) ! [OUT]

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
 
    call initialize

    !###########################################################################
 
    !########## main ##########

#ifdef FIPP
    call fipp_start
#endif

    LOG_NEWLINE
    LOG_PROGRESS(*) 'START TIMESTEP'
    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_Loop', 0)
 
    do

      !*******************************************

      ! report current time
      call TIME_manager_checkstate

      if (TIME_DOresume) then
        ! set state from restart file
        call restart_read
        ! history & monitor file output 
        call FILE_MONITOR_meshfield_write('MAIN', TIME_NOWSTEP)
        call FILE_HISTORY_meshfield_write
      end if   
      
      !* Advance time *********************************

      call TIME_manager_advance
      call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )
      
      !* change to next state *************************

      !- USER
      call user_%update_pre( atmos )

      if ( ocean%IsActivated() .and. ocean%time_manager%do_step ) call ocean%update()
      if ( atmos%IsActivated() .and. atmos%time_manager%do_step ) call atmos%update()
      call user_%update( atmos )

      !* restart and monitor output *******************
      if ( atmos%IsActivated() ) call atmos%vars%Monitor()
      call restart_write
      call FILE_MONITOR_meshfield_write('MAIN', TIME_NOWSTEP)


      !* calc tendencies and diagnostics *************

      if ( atmos  %IsActivated() .and. atmos%time_manager%do_step ) call atmos%calc_tendency( force=.false. )
      if ( ocean  %IsActivated() .and. ocean%time_manager%do_step ) call ocean%calc_tendency( force=.false. )
      if ( coupler%IsActivated() .and. atmos%time_manager%do_step ) call atmos%calc_tendency_from_sflux(force=.false.)
      call user_%calc_tendency( atmos )
  
      !* output history files *************************

      if ( atmos%IsActivated() ) call atmos%vars%History()
      if ( atmos%dyn_proc%IsActivated() ) call atmos%dyn_proc%dyn_vars%History()
      if ( atmos%phy_tb_proc%IsActivated() ) call atmos%phy_tb_proc%vars%History()
      if ( atmos%phy_mp_proc%IsActivated() ) call atmos%phy_mp_proc%vars%History()
      if ( atmos%phy_sfc_proc%IsActivated() ) call atmos%phy_sfc_proc%vars%History()
      if ( atmos%phy_rd_proc%IsActivated() ) call atmos%phy_rd_proc%vars%History()

      if ( ocean%IsActivated() ) call ocean%vars%History()

      call FILE_HISTORY_meshfield_write
      
      !*******************************************
      if (TIME_DOend) exit
      
      if( IO_L ) call flush(IO_FID_LOG)
    end do

    call PROF_rapend('Main_Loop', 0)

    LOG_PROGRESS(*) 'END TIMESTEP'
    LOG_NEWLINE

#ifdef FIPP
    call fipp_stop
#endif

    !########## Finalize ##########
    call finalize

    return
  end subroutine dg_driver

  !----------------------------
  
!OCL SERIAL
  subroutine initialize()

    use scale_const, only: CONST_setup
    use scale_calendar, only: CALENDAR_setup
    use scale_random, only: RANDOM_setup
    use scale_atmos_hydrometeor, only: ATMOS_HYDROMETEOR_setup
    use scale_time_manager, only: TIME_DTSEC

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

    ! setup random number
    call RANDOM_setup

    ! setup tracer index
    call ATMOS_HYDROMETEOR_setup

    ! setup a module for restart file
    call FILE_restart_meshfield_setup
    call TIME_manager_Init( &
      setup_TimeIntegration = .true.,                   &
      restart_in_basename   =  restart_file%in_basename )

    ! setup statistics
    call MeshField_statistics_setup

    ! setup monitor
    call FILE_monitor_meshfield_setup( TIME_DTSEC )

    ! setup sub-models
    call atmos%setup()
    call ocean%setup()
    call coupler%setup()
    call user_%setup( atmos )

    call coupler%evaluate_activation( ocean )
    call atmos%set_coupler( coupler )
    call ocean%set_coupler( coupler )

    call atmos%setup_vars()

    if ( ocean%IsActivated() ) call ocean%setup_vars()
    if ( coupler%IsActivated() ) call coupler%setup_vars( atmos%mesh%ptr_mesh, ocean%mesh%ptr_mesh )
    
    ! report information of time intervals
    call TIME_manager_report_timeintervals

    !----------------------------------------

    call PROF_rapend('Initialize', 0)

    return
  end subroutine initialize

!OCL SERIAL  
  subroutine finalize()
    use scale_file, only: &
      FILE_Close_All
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    use scale_file_restart_meshfield, only: &
      FILE_restart_meshfield_finalize
    use scale_file_monitor_meshfield, only: &
      FILE_monitor_meshfield_final
    use scale_time_manager, only: &
      TIME_manager_Final   
    implicit none
    
    !----------------------------------------------
    call PROF_setprefx('FIN')
    call PROF_rapstart('All', 1)

    call PROF_rapstart('Monit', 2)
    call FILE_monitor_meshfield_final
    call PROF_rapend  ('Monit', 2)

    !-
    call PROF_rapstart('File', 2)
    call FILE_HISTORY_meshfield_finalize
    call FILE_restart_meshfield_finalize
    call PROF_rapend  ('File', 2)

    ! finalization sub-models
    call atmos%finalize()
    call ocean%finalize()
    call coupler%finalize()
    call user_%final()

    !-
    call TIME_manager_Final

    call PROF_rapend  ('All', 1)
    call PROF_rapreport

    return
  end subroutine finalize

!OCL SERIAL
  subroutine restart_read()
    implicit none
    !----------------------------------------

    !- read restart data

    if ( atmos%isActivated() ) then
      call atmos%vars%Read_restart_file( atmos%mesh, atmos%dyn_proc%dyncore_driver )
      if ( atmos%phy_rd_proc%IsActivated() ) call atmos%phy_rd_proc%vars%Read_restart_file()
    end if
    if ( ocean%isActivated() ) then
      call ocean%vars%Read_restart_file( ocean%mesh )
    end if

    !- Setup surface condition

    if ( atmos%isActivated() ) call atmos%set_surface( countup=.false. )
    if ( ocean%isActivated() ) call ocean%set_surface( countup=.false. )

    !- Calculate the tendencies

    if ( atmos%IsActivated() ) call atmos%calc_tendency( force= .true. )
    if ( ocean%IsActivated() ) call ocean%calc_tendency( force= .true. )
    if ( coupler%IsActivated() ) call atmos%calc_tendency_from_sflux(force=.true.)

    call user_%calc_tendency( atmos )

    !- History & Monitor 

    if ( atmos%isActivated() ) then
      call atmos%vars%History()
      if ( atmos%phy_sfc_proc%IsActivated() ) &
        call atmos%phy_sfc_proc%vars%History()
      if ( atmos%phy_tb_proc%IsActivated() )  &
        call atmos%phy_tb_proc%vars%History()
      if ( atmos%phy_mp_proc%IsActivated() )  &
        call atmos%phy_mp_proc%vars%History()
      if ( atmos%phy_rd_proc%IsActivated() )  &
        call atmos%phy_rd_proc%vars%History()
      call atmos%vars%Monitor()
    end if

    if ( ocean%isActivated() ) then
      call ocean%vars%History()
    end if

    return
  end subroutine restart_read

!OCL SERIAL
  subroutine restart_write()
    use scale_file_restart_meshfield, only: &
      restart_file    
    implicit none

    logical :: is_restart_write_atmos
    logical :: is_restart_write_ocean
    !----------------------------------------

    if ( .not. restart_file%flag_output ) return
    
    is_restart_write_atmos = atmos%isActivated() .and. atmos%time_manager%do_restart
    is_restart_write_ocean  = ocean%isActivated()  .and. ocean%time_manager%do_restart

    !- Preprocess
    if ( is_restart_write_atmos ) then
      call atmos%vars%Write_restart_file_prep()
      if ( atmos%phy_rd_proc%IsActivated() ) call atmos%phy_rd_proc%vars%Write_restart_file_prep()
    end if
    if ( is_restart_write_ocean ) then
      call ocean%vars%Write_restart_file_prep()
    end if

    !- Write
    if ( is_restart_write_atmos ) then
      call atmos%vars%Write_restart_file()
      if ( atmos%phy_rd_proc%IsActivated() ) call atmos%phy_rd_proc%vars%Write_restart_file()
    end if
    if ( is_restart_write_ocean ) then
      call ocean%vars%Write_restart_file()
    end if

    !- Postprocess
    if ( is_restart_write_atmos ) then
      call atmos%vars%Write_restart_file_post()
      if ( atmos%phy_rd_proc%IsActivated() ) call atmos%phy_rd_proc%vars%Write_restart_file_post()
    end if
    if ( is_restart_write_ocean ) call ocean%vars%Write_restart_file_post()

    return
  end subroutine restart_write


end module mod_dg_driver