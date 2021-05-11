!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_sw_bnd
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,  &
    Rdry => CONST_Rdry,  &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00

  use scale_sparsemat  
  use scale_element_base
  use scale_element_quadrilateral
  use scale_localmesh_2d
  use scale_mesh_rectdom2d

  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D

  use scale_mesh_bndinfo, only: &
    MeshBndInfo

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  type, public :: SWBnd
  contains
    procedure :: Init => SW_bnd_setup
    procedure :: Final => SW_bnd_finalize
  end type
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  logical, save :: is_initialized = .false.

contains
  subroutine SW_bnd_setup( this )
    implicit none
    class(SWBnd), intent(inout) :: this
    !-----------------------------------------------

    is_initialized = .false.
    return
  end subroutine SW_bnd_setup

  subroutine SW_bnd_finalize( this )
    implicit none
    class(SWBnd), intent(inout) :: this
    !--------------------------------------

    if (is_initialized) then
    end if

    is_initialized = .false.
    return
  end subroutine SW_bnd_finalize  
end module mod_sw_bnd