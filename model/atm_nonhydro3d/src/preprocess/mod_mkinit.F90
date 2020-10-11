!-------------------------------------------------------------------------------
!> module INITIAL
!!
!! @par Description
!!          subroutines for preparing initial data
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_mkinit
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
  !++ Public procedure
  !
  public :: MKINIT_setup
  public :: MKINIT


  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKINIT_setup
    implicit none

    character(len=H_SHORT) :: initname = 'NONE'

    namelist / PARAM_MKINIT / &
      initname

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("MKINIT_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("MKINIT_setup",*) 'Not appropriate names in namelist PARAM_MKINIT. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT)

  end subroutine MKINIT_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine MKINIT( output )
    implicit none
    logical, intent(out) :: output
    !---------------------------------------------------------------------------

    return
  end subroutine MKINIT
end module mod_mkinit