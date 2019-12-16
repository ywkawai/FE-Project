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

contains
  subroutine dg_driver( &
    comm_world, intercomm_parent, intercomm_child, cnf_fname )
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
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_setup
    implicit none

    !----------------------------------------------

    ! namelist compatibility check
    !call ADMIN_versioncheck

    ! setup process
    !call PRC_CARTESC_setup

    ! setup PROF
    call PROF_setup

    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize', 0)
    
    ! setup constants
    call CONST_setup

    ! setup calendar & initial time
    call CALENDAR_setup
!    call TIME_manager_Init

    ! setup submodel administrator
    call ATMOS_admin_setup
    
    ! Setup a module to manage mesh
    !call ATMOS_MESH_setup()
    
    ! Setup file I/O
    !call FILE_HISTORY_meshfield_setup( mesh2d_=mesh )

    call PROF_rapend('Initialize', 0)

    return
  end subroutine initialize

  subroutine finalize()
    implicit none
    
    !----------------------------------------------
    call PROF_setprefx('FIN')

    call PROF_rapstart('All', 1)

    call PROF_rapstart('File', 2)
    !call FILE_HISTORY_meshfield_finalize()
    call PROF_rapend('File', 2)

    call PROF_rapend  ('All', 1)

    call PROF_rapreport()

    return
  end subroutine finalize
end module mod_dg_driver