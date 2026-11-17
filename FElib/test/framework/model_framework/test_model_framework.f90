#include "scalelib.h"
program test_model_framework
  use scale_precision
  use scale_prc
  use scale_io
  use scale

  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d

  use scale_localmeshfield_base, only: LocalMeshField1D
  use scale_meshfield_base, only: MeshField1D
  use mod_atmos_component, only: &
    AtmosComponent
  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final

  implicit none

  type(AtmosComponent) :: atmos
  !------------------------------------------------------------------

  call init()
  call final()

contains
  subroutine init()
    implicit none
    !---------

    call UNITTEST_init( "test_model_framework", conf_fname="test.conf", allow_noconf=.false., &
      use_time=.true. )

    !------
    call atmos%setup()

    return
  end subroutine init

  subroutine final()
    implicit none
    !------------------

    call atmos%finalize()
    call UNITTEST_final()

    return
  end subroutine final


end program test_model_framework
