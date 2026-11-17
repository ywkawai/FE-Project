#include "scalelib.h"
program test_field2d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale_file_history
  use scale

  use scale_element_quadrilateral
  use scale_localmesh_2d
  use scale_mesh_rectdom2d

  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  use scale_time_manager, only: &
    TIME_manager_Init , TIME_manager_Final,            &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,            &
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

  integer, parameter :: NeGX = 2
  integer, parameter :: NeGY = 2
  integer, parameter :: NLocalMeshPerPrc = 1

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: dom_ymin = -1.0_RP
  real(RP), parameter :: dom_ymax = +1.0_RP

  type(QuadrilateralElement) :: refElem
  integer, parameter :: PolyOrder = 2
  
  type(MeshRectDom2D), target :: mesh

  type(MeshField2D), target :: q
  character(len=H_short), parameter :: q_VARNAME = "q"
  character(len=H_short), parameter :: q_DESC    = "q"
  character(len=H_short), parameter :: q_UNITS   = "K"
  integer, save :: HST_ID(1)

  type(MeshFieldCommRectDom2D) :: fields_comm
  type(MeshFieldContainer), save :: field_list(1)  

  integer :: n
  type(LocalMesh2D), pointer :: lcmesh

  !---------------
  
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
    class(LocalMesh2D), intent(in) :: lmesh
    class(QuadrilateralElement), intent(in) :: elem
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

    call UNITTEST_init( "test_field_rectdom2d_quadrilateral",          &
      comm=comm, myrank=myrank, nprocs=nprocs, ismaster=ismaster, use_time=.true. )

    !------
    call refElem%Init(PolyOrder, .true.)

    call mesh%Init( &
      NeGX, NeGY,                             &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
      .true., .true.,                         &
      refElem, NLocalMeshPerPrc, 1, 1         )
    
    call mesh%Generate()

    !---
    call q%Init( q_VARNAME, q_UNITS, mesh )
    call fields_comm%Init(1, 0, 0, mesh)
    
    call FILE_HISTORY_meshfield_setup( mesh2d_=mesh )
    call FILE_HISTORY_reg( q_VARNAME, q_DESC, q_UNITS, HST_ID(1), dim_type='XY')

    !---

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

    field_list(1)%field2d => q
    call TEST_perform_comm( fields_comm, field_list )

    return
  end subroutine perform_comm
    
  subroutine check_interior_data(n, lcmesh_, lcfield, lcfield_fl)
    integer, intent(in) :: n
    type(LocalMesh2D), intent(in) :: lcmesh_
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
    type(LocalMesh2D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ans_bx(refElem%Nfp,lcmesh_%NeX)
    integer :: ans_by(refElem%Nfp,lcmesh_%NeY)    
    integer :: haloInd_s, haloInd_e
    integer :: i, j, k, f, fp
    integer :: tileID

    !------------------------------

    tileID = lcmesh_%tileID

    write(*,*) "Check values in halo.."
    write(*,*) "tileID_globalMap=", mesh%tileID_globalMap(:,tileID)
    write(*,*) "tileFaceID_globalMap=", mesh%tileFaceID_globalMap(:,tileID)


    haloInd_s = refElem%Np*lcmesh_%Ne + 1    
    haloInd_e = haloInd_s + refElem%Nfp*lcmesh_%NeX - 1
    do i=1, lcmesh_%NeX
      k = i + (lcmesh_%NeY-1)*lcmesh_%NeX
      do fp=1, refElem%Nfp
        ans_bx(fp,i) = TEST_field_val(mesh%tileID_globalMap(1,tileID), k, refElem%Fmask(fp,3))
      end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_bx(:,:), &
      'check_halo_data: q_halo_south' )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + refElem%Nfp*lcmesh_%NeY - 1
    do j=1, lcmesh_%NeY
      k = 1 + (j-1)*lcmesh_%NeX
      do fp=1, refElem%Nfp
        ans_by(fp,j) = TEST_field_val(mesh%tileID_globalMap(2,tileID), k, refElem%Fmask(fp,4))
      end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_by(:,:), &
      'check_halo_data: q_halo_east' )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + refElem%Nfp*lcmesh_%NeX - 1    
    do i=1, lcmesh_%NeX
      k = i
      do fp=1, refElem%Nfp
        ans_bx(fp,i) = TEST_field_val(mesh%tileID_globalMap(3,tileID), k, refElem%Fmask(fp,1))
      end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_bx(:,:), &
      'check_halo_data: q_halo_north' )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + refElem%Nfp*lcmesh_%NeY - 1    
    do j=1, lcmesh_%NeY
      k = lcmesh_%NeX + (j-1)*lcmesh_%NeX
      do fp=1, refElem%Nfp
        ans_by(fp,j) = TEST_field_val(mesh%tileID_globalMap(4,tileID), k, refElem%Fmask(fp,2))
      end do
    end do    
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_by(:,:), &
      'check_halo_data: q_halo_west' )    

  end subroutine check_halo_data


end program test_field2d
