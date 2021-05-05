!-------------------------------------------------------------------------------
!> module SCALE-DG prep
!!
!! @par Description
!!          This program is driver of preprocess tools
!!          1) boundary data (e.g. topography, land use index)
!!          2) initial data for ideal/real test cases
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_dg_prep
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
    USER_setup, &
    USER_mkinit

  use mod_mkinit, only: &
    MKINIT
  
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
  public :: dg_prep
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
  subroutine dg_prep(                     &
    comm_world, cnf_fname, path, add_path )

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
    logical :: ismaster

    logical :: output
    !---------------------------------------------------------------------------

    !########## Initial setup ##########

#ifdef SCALE_DEVELOP
    ! setup standard I/O
    if ( add_path .and. path /= "" ) then
      call IO_setup( MODELNAME, trim(path)//cnf_fname, prefix=path )
    else
#endif
      call IO_setup( MODELNAME, trim(path)//cnf_fname )
#ifdef SCALE_DEVELOP
    end if
#endif

    ! setup MPI
    call PRC_LOCAL_setup( comm_world, & ! [IN]
                          myrank,     & ! [OUT]
                          ismaster    ) ! [OUT]

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
 
    call initialize()

    !###########################################################################
 
    !########## main ##########
    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_prep', 0)

    !- Execute preprocess

    !- Execute mktopo

    !- Re-setup

    !- Execute mkinit
    call PROF_rapstart('MkInit',1)
    call MKINIT( output, &
      atmos%mesh, atmos%vars%PROGVARS_manager, atmos%vars%AUXVARS_manager )
    call USER_mkinit( atmos )
    call PROF_rapend  ('MkInit',1)
    call PROF_rapend('Main_prep', 0)

    !- Output

    ! call TOPOGRAPHY_write

    if (output) then
      call PROF_rapstart('MkInit_restart',1)
      call restart_write()
      call PROF_rapend  ('MkInit_restart',1)
    end if

    !########## Finalize ##########
    call finalize()    

    return
  end subroutine dg_prep

  !----------------------------
  
  subroutine initialize()

    use scale_const, only: CONST_setup
    use scale_calendar, only: CALENDAR_setup
    use scale_random, only: RANDOM_setup

    use scale_file_restart_meshfield, only: &
      FILE_restart_meshfield_setup
    use scale_time_manager, only: TIME_manager_Init

    use mod_user, only: USER_setup    
    use mod_mkinit, only: MKINIT_setup

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
    call TIME_manager_Init( .false. )

    ! setup random number
    call RANDOM_setup

    ! setup a module for restart file
    call FILE_restart_meshfield_setup

    ! setup submodels
    call  atmos%setup()

    ! setup mkinit
    call MKINIT_setup()

    ! setup mod_user
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

  subroutine restart_write
    implicit none    
    !----------------------------------------

    if ( atmos%isActivated() ) call atmos%vars%Write_restart_file()

    return
  end subroutine 

end module mod_dg_prep