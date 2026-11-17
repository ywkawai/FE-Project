!-------------------------------------------------------------------------------
!> Module common / unit test (utility)
!!
!! @par Description
!!      A module to provide common routines for the setup and finalization of unit tests.
!!
!!      UNITTEST_init performs
!!        MPI start -> process setup -> IO setup -> LOG setup
!!        -> (optional) CONST setup -> (optional) CALENDAR & TIME_manager setup.
!!      UNITTEST_final gathers the numbers of assertions and failures over all processes,
!!      reports the summary, and stops the program with an error (PRC_abort) if any assertion fails.
!!
!!      The verbose mode of assertions can be switched with namelist PARAM_UNITTEST:
!!        &PARAM_UNITTEST UNITTEST_VERBOSE = .true. /
!!
!! @author Yuta Kawai, Team SCALE
!!
#include "scaleFElib.h"
module scale_unittest_util
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_unittest_assert, only: &
    UNITTEST_assert_set_verbose,     &
    UNITTEST_assert_get_num_assert,  &
    UNITTEST_assert_get_num_failure

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: UNITTEST_init
  public :: UNITTEST_final

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  character(len=H_MID), private :: unittest_name = ''
  logical, private :: time_manager_initialized = .false.

contains

  !> Setup a unit test
  subroutine UNITTEST_init( test_name, &
      comm, myrank, nprocs, ismaster,  &
      conf_fname, allow_noconf,        &
      use_const, use_time, verbose     )

    use scale_const, only: CONST_setup
    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only: TIME_manager_Init
    implicit none

    character(len=*), intent(in) :: test_name             !< Name of the unit test
    integer, intent(out), optional :: comm                !< Communicator
    integer, intent(out), optional :: myrank              !< My rank
    integer, intent(out), optional :: nprocs              !< Number of processes
    logical, intent(out), optional :: ismaster            !< Master process or not
    character(len=*), intent(in), optional :: conf_fname  !< Name of configuration file (if absent, the first command line argument is used)
    logical, intent(in), optional :: allow_noconf         !< Allow the test to run without configuration file (default: .true.)
    logical, intent(in), optional :: use_const            !< Call CONST_setup (default: .false.)
    logical, intent(in), optional :: use_time             !< Call CALENDAR_setup and TIME_manager_Init (default: .false.)
    logical, intent(in), optional :: verbose              !< Verbose mode of assertions (default: .false. or the value specified in PARAM_UNITTEST)

    integer :: comm_, myrank_, nprocs_
    logical :: ismaster_
    logical :: allow_noconf_

    logical :: UNITTEST_VERBOSE
    namelist / PARAM_UNITTEST / UNITTEST_VERBOSE
    integer :: ierr
    !------------------------------------------------

    unittest_name = test_name
    allow_noconf_ = .true.
    if ( present(allow_noconf) ) allow_noconf_ = allow_noconf

    call PRC_MPIstart( comm_ )

    call PRC_SINGLECOM_setup( comm_, & ! [IN]
      nprocs_, myrank_, ismaster_    ) ! [OUT]

    call PRC_ERRHANDLER_setup( .false., ismaster_ ) ! [IN]

    ! setup scale_io
    if ( present(conf_fname) ) then
      call IO_setup( test_name, conf_fname, allow_noconf = allow_noconf_ )
    else
      call IO_setup( test_name, allow_noconf = allow_noconf_ )
    end if

    ! setup log
    call IO_LOG_setup( myrank_, ismaster_ )

    ! setup the verbose mode of assertions
    UNITTEST_VERBOSE = .false.
    if ( present(verbose) ) UNITTEST_VERBOSE = verbose
    if ( IO_FID_CONF > 0 ) then
      ! Note: iostat is specified for rewind because IO_FID_CONF may be connected to a special file (/dev/null)
      !       when no configuration file is given.
      rewind(IO_FID_CONF, iostat=ierr)
      read(IO_FID_CONF,nml=PARAM_UNITTEST,iostat=ierr)
      if ( ierr > 0 ) then !--- fatal error
        LOG_ERROR("UNITTEST_init",*) 'Not appropriate names in namelist PARAM_UNITTEST. Check!'
        call PRC_abort
      end if
    end if
    call UNITTEST_assert_set_verbose( UNITTEST_VERBOSE )

    ! setup constants
    if ( present(use_const) ) then
      if ( use_const ) call CONST_setup
    end if

    ! setup calendar & initial time
    time_manager_initialized = .false.
    if ( present(use_time) ) then
      if ( use_time ) then
        call CALENDAR_setup
        call TIME_manager_Init
        time_manager_initialized = .true.
      end if
    end if

    LOG_INFO("UNITTEST_init",*) 'Start unit test: ', trim(test_name)

    if ( present(comm)     ) comm     = comm_
    if ( present(myrank)   ) myrank   = myrank_
    if ( present(nprocs)   ) nprocs   = nprocs_
    if ( present(ismaster) ) ismaster = ismaster_

    return
  end subroutine UNITTEST_init

  !> Finalize a unit test
  !! Note: Objects used in the test (mesh, fields, etc.) should be finalized before calling this routine.
  subroutine UNITTEST_final()
    use mpi
    use scale_time_manager, only: TIME_manager_Final
    implicit none

    integer :: num_lc(2), num(2)
    integer :: ierr
    !------------------------------------------------

    if ( time_manager_initialized ) then
      call TIME_manager_Final()
      time_manager_initialized = .false.
    end if

    num_lc(1) = UNITTEST_assert_get_num_assert()
    num_lc(2) = UNITTEST_assert_get_num_failure()
    call MPI_Allreduce( num_lc, num, 2, MPI_INTEGER, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr )

    if ( PRC_IsMaster ) then
      write(*,'(a)') '***********************************************************************'
      write(*,'(a,a,a,i0,a,i0,a)') ' Unit test [', trim(unittest_name), ']: ', &
        num(1), ' assertions, ', num(2), ' failures (summed over all processes)'
      if ( num(2) == 0 ) then
        write(*,'(a)') ' => PASSED'
      else
        write(*,'(a)') ' => FAILED'
      end if
      write(*,'(a)') '***********************************************************************'
    end if
    LOG_INFO("UNITTEST_final",'(a,i0,a,i0)') 'Number of assertions = ', num(1), ', failures = ', num(2)

    if ( num(2) > 0 ) then
      call PRC_abort
    end if

    call PRC_MPIfinish()

    return
  end subroutine UNITTEST_final

end module scale_unittest_util
