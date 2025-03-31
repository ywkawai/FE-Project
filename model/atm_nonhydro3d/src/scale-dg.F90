!-------------------------------------------------------------------------------
!> Program SCALE-DG (a launcher of main routine)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          numerical model with DGM for regional weather, regional climate, and large-Eddy Simulation (LES)
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
program scaledg
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_dg_launcher
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

  logical               :: EXECUTE_PREPROCESS           = .false. ! execute preprocess tools?
  logical               :: EXECUTE_MODEL                = .true.  ! execute main model?

  !-----------------------------------------------------------

  call launcher( EXECUTE_PREPROCESS, & ! (in)
                 EXECUTE_MODEL       ) ! (in)

  stop
end program scaledg
