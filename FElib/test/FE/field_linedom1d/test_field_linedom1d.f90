#include "scalelib.h"
program test_field1d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale_file_history
  use scale

  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d

  use scale_localmeshfield_base, only: LocalMeshField1D
  use scale_meshfield_base, only: MeshField1D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_1d, only: MeshFieldComm1D

  use scale_time_manager, only: &
    TIME_manager_Init , TIME_manager_Final,            &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,        &
    TIME_DTSEC, TIME_NSTEP 
    
  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_setup,   &
    FILE_HISTORY_meshfield_put,     &
    FILE_HISTORY_meshfield_write,   &
    FILE_HISTORY_meshfield_finalize

  use mod_test_meshfield_util, only: &
    TEST_field_val, TEST_perform_comm, TEST_check_face_node
  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: UNITTEST_assert_equal
  implicit none

  integer, parameter :: NeGX = 4
  integer, parameter :: NLocalMeshPerPrc = 1

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP

  type(LineElement) :: refElem
  integer, parameter :: PolyOrder = 2
  
  type(MeshLineDom1D), target :: mesh

  type(MeshField1D), target :: q
  character(len=H_short), parameter :: q_VARNAME = "q"
  character(len=H_short), parameter :: q_DESC    = "q"
  character(len=H_short), parameter :: q_UNITS   = "K"
  integer, save :: HST_ID(1)

  type(MeshFieldComm1D) :: fields_comm
  type(MeshFieldContainer), save :: field_list(1)  

  type(LocalMesh1D), pointer :: lcmesh
  integer :: n
  
  !------------------------------------------------------------------

  call init()

  !----
  write(*,*) "* Check data communication.."
  call perform_comm()

  do n=1, mesh%LOCAL_MESH_NUM
    lcmesh => mesh%lcmesh_list(n)

    write(*,*) "Check interior & halo data.."
    write(*,*) "tileID=", lcmesh%tileID
    !$acc update host(q%local(n)%val)
    call check_interior_data(n, mesh%lcmesh_list(n), q%local(n)%val, q%local(n)%val)
    call check_halo_data(n, mesh%lcmesh_list(n), q%local(n)%val, q%local(n)%val)
  end do
  !----

  call final()

contains

  subroutine get_field_val(val, lmesh, elem)
    class(LocalMesh1D), intent(in) :: lmesh
    class(LineElement), intent(in) :: elem
    real(RP), intent(out) :: val(elem%Np,lmesh%NeA)

    integer :: tileID, ke, p
    !----------------------------

    tileID = lmesh%tileID
    !$acc parallel loop collapse(2) present(lmesh,val)
    do ke=lmesh%NeS, lmesh%NeE
    do p=1, elem%Np
      val(p,ke) = TEST_field_val(tileID, ke, p)
    end do
    end do
  end subroutine get_field_val

  subroutine init()

    integer :: comm, myrank, nprocs
    logical :: ismaster
    !---------
    call UNITTEST_init( "test_field_linedom1d",          &
      comm=comm, myrank=myrank, nprocs=nprocs, ismaster=ismaster, use_time=.true. )

    !------
    call refElem%Init(PolyOrder, .true.)

    call mesh%Init( &
      NeGX,                              &
      dom_xmin, dom_xmax,                &
      refElem, NLocalMeshPerPrc )
    
    call mesh%Generate()

    !---
    call q%Init( q_VARNAME, q_UNITS, mesh )
    call fields_comm%Init(1, 0, mesh)
    
    call FILE_HISTORY_meshfield_setup( mesh )
    call FILE_HISTORY_reg( q_VARNAME, q_DESC, q_UNITS, HST_ID(1), dim_type='X')

    !---
    do n=1, mesh%LOCAL_MESH_NUM
      call get_field_val( q%local(n)%val, &
        mesh%lcmesh_list(n), refElem)
    end do
    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_write()    
    !---

  end subroutine init

  subroutine final()

    call FILE_HISTORY_meshfield_finalize()

    call q%Final()
    call fields_comm%Final()
    call mesh%Final()
    call refElem%Final()
    
    call UNITTEST_final()

  end subroutine final

  subroutine perform_comm()
    implicit none
    type(MeshFieldContainer) :: field_list(1)
    !---------------------------

    field_list(1)%field1d => q
    call TEST_perform_comm( fields_comm, field_list )

    return
  end subroutine perform_comm

  subroutine check_interior_data(n, lcmesh_, lcfield, lcfield_fl)
    integer, intent(in) :: n
    type(LocalMesh1D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ke, f, p
    integer :: fnid
    real(RP) :: ans(refElem%Np,lcmesh_%NeA)
    !------------------------------

    !$acc data create(ans)
    call get_field_val( ans, &
      lcmesh_, refElem )
    !$acc update host(ans)

    do ke=lcmesh_%NeS, lcmesh_%NeE
      call UNITTEST_assert_equal( int(lcfield(:,ke)), int(ans(:,ke)), 'check_interior_data: q', idx=ke )
      do f=1, refElem%Nfaces
      do p=1, refElem%Nfp
        fnid = p+(f-1)*refElem%Nfp
        call TEST_check_face_node( lcmesh_, lcfield_fl, ke, f, fnid )
      end do
      end do
    end do
    !$acc end data
  end subroutine check_interior_data

  subroutine check_halo_data(n, lcmesh_, lcfield, lcfield_fl)
    implicit none

    integer, intent(in) :: n
    type(LocalMesh1D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ans(2)
    integer :: haloInd_s, haloInd_e

    !------------------------------

    write(*,*) "Check values in halo.."
    write(*,*) "tileID_globalMap=", mesh%tileID_globalMap(:,n)
    write(*,*) "tileFaceID_globalMap=", mesh%tileFaceID_globalMap(:,n)

    haloInd_s = refElem%Np*lcmesh_%Ne + 1
    haloInd_e = haloInd_s + 2 - 1

    ans(1) = TEST_field_val(mesh%tileID_globalMap(1,n), lcmesh_%Ne, refElem%Np) 
    ans(2) = TEST_field_val(mesh%tileID_globalMap(2,n), 1, 1) 
    
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans(:), 'check_halo_data: q_halo' )

  end subroutine check_halo_data


end program test_field1d
