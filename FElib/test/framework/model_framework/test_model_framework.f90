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

  implicit none

  type(AtmosComponent) :: atmos
  !------------------------------------------------------------------

  call init()
  call final()

contains
  subroutine init()
    use scale_calendar, only: CALENDAR_setup

    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: n, k, p

    !---------
    call PRC_MPIstart( comm )
    
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    call IO_setup( "test_model_framework", "test.conf", allow_noconf = .false. )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
    
    !------
    call atmos%setup()

    return
  end subroutine init

  subroutine final()

    call atmos%finalize()
    call PRC_MPIfinish()

    return
  end subroutine final


end program test_model_framework
