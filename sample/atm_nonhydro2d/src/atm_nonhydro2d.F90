#include "scaleFElib.h"
program atm_nonhydro2d
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
    TIME_manager_advance,                              &
    TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP,            &
    TIME_NOWDAYSEC, TIME_DTSEC, TIME_NSTEP

  implicit none

  integer :: nowstep
  character(len=H_MID) :: timelabel
  !-------------------------------------------------------

  call init()
  
  do nowstep=1, TIME_NSTEP
    !* Advance time
    call TIME_manager_advance()
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP )

    call update()
    call calc_tendency()

    if (mod(nowstep,1000) == 0) then 
      write(timelabel,'(I4.4,I2.2,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') &
        TIME_NOWDATE(1:3), '-', TIME_NOWDATE(4), ":", TIME_NOWDATE(5), ":", TIME_NOWDATE(6)
      LOG_PROGRESS('(A,A)') "time = ", trim(timelabel)
    end if

    if( IO_L ) call flush(IO_FID_LOG)
  end do

  call final()

contains
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine update()
    use mod_atmos_vars, only: ATMOS_VARS_output
    use mod_atmos_dyn_driver, only: ATMOS_DYN_driver
    use mod_user, only: USER_update
    implicit none
    !------------------------------------------------------------------------
    
    !- ATMOS
    call ATMOS_DYN_driver()
    call ATMOS_VARS_output( TIME_NOWDAYSEC )

    !- USER
    call USER_update()

    return
  end subroutine update

  subroutine calc_tendency()
    use mod_user, only: USER_calc_tendency
    implicit none
    !------------------------------------------------------------------------
    
    !- ATMOS

    !- USER
    call USER_calc_tendency()

    return
  end subroutine calc_tendency

  subroutine init()
    use scale_const, only: CONST_setup
    use scale_calendar, only: CALENDAR_setup 
    use scale_time_manager, only: TIME_manager_Init
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_setup

    use mod_atmos_vars, only: &
      ATMOS_VARS_setup, ATMOS_VARS_output
    use mod_atmos_mesh, only: &
      ATMOS_MESH_setup, mesh
    use mod_atmos_bnd, only: &
      ATMOS_bnd_setup, ATMOS_bnd_setBCInfo
    use mod_atmos_dyn_driver, only: &
      ATMOS_DYN_driver_setup
    use mod_user, only: USER_setup
    implicit none

    character(len=H_SHORT) :: exp_name

    namelist /PARAM_ATM_NONHYDRO2D/ &
      exp_name
    
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
    call IO_setup( "atm_nonhydro2d", trim(conf_name) )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
 

    !--- read namelist

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATM_NONHYDRO2D,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_ATM_NONHYDRO2D. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATM_NONHYDRO2D)

    ! setup profiler
    call PROF_setup
    call PROF_rapstart( "total", 0 )
    call PROF_rapstart( "init", 1 )

    ! setup constants
    call CONST_setup

    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init

    !-
    ! Setup a module to manage boundary conditions
    call ATMOS_bnd_setup()

    ! Setup a module to manage mesh
    call ATMOS_MESH_setup()
    
    ! Setup a module to output data
    call FILE_HISTORY_meshfield_setup( mesh2d_=mesh )

    ! Setup a module to manage variables
    call ATMOS_VARS_setup()

    ! Setup a module for dynamics
    call ATMOS_DYN_driver_setup()

    ! Set some informations for boundary conditions
    call ATMOS_bnd_setBCInfo()
    
    call USER_setup()

    !--
    call ATMOS_VARS_output( TIME_NOWDAYSEC )
    LOG_PROGRESS('(A,F13.5,A)') "time=", real(0.0_RP), "[s]"

    !---
    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()
    use scale_time_manager, only: TIME_manager_Final
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_finalize    
    use mod_atmos_vars, only: ATMOS_VARS_finalize
    use mod_atmos_mesh, only: ATMOS_MESH_finalize
    use mod_atmos_bnd, only: ATMOS_bnd_finalize
    use mod_atmos_dyn_driver, only: ATMOS_DYN_driver_finalize  
    implicit none
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )

    call ATMOS_DYN_driver_finalize()
    call ATMOS_bnd_finalize()
    call ATMOS_VARS_finalize()
    call ATMOS_MESH_finalize()   
    
    call FILE_HISTORY_meshfield_finalize()
    call TIME_manager_Final()

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()

    return
  end subroutine final

end program atm_nonhydro2d
