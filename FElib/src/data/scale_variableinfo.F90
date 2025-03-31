!> module FElib / Data / base
!!
!! @par Description
!!           A module for managing information of variable
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_variableinfo
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: VariableInfo
    integer                :: keyID
    character(len=H_SHORT) :: NAME
    character(len=H_MID)   :: DESC
    character(len=H_SHORT) :: UNIT
    integer                :: ndims
    character(len=H_SHORT) :: dim_type
    character(len=H_MID)   :: STDNAME
  end type VariableInfo

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !------------------

end module scale_variableinfo