#include "scalelib.h"
program test_field_cubedspheredom2d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale_file_history
  use scale
  use scale_const, only: &
    RPlanet => CONST_Radius

  use scale_element_quadrilateral
  use scale_localmesh_2d
  use scale_mesh_cubedspheredom2d

  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_cubedspheredom2d, only: MeshFieldCommCubedSphereDom2D

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
    TEST_field_val, TEST_perform_comm, TEST_check_face_node, TEST_get_face_index
  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: UNITTEST_assert_equal
  implicit none

  integer, parameter :: PolyOrder = 2
  integer, parameter :: NeGX = 3
  integer, parameter :: NeGY = 3
  integer, parameter :: NLocalMeshPerPrc = 6

  type(QuadrilateralElement) :: refElem  
  type(MeshCubedSphereDom2D), target :: mesh

  type(MeshField2D), target :: lon
  type(MeshField2D), target :: lat
  type(MeshField2D), target :: q
  character(len=H_short), parameter :: q_VARNAME = "q"
  character(len=H_short), parameter :: q_DESC    = "q"
  character(len=H_short), parameter :: q_UNITS   = "K"
  integer, save :: HST_ID(3)

  type(MeshFieldCommCubedSphereDom2D) :: fields_comm
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
    call check_interior_data(n, mesh%lcmesh_list(n), q%local(n)%val, q%local(n)%val)
    call check_halo_data(n, mesh%lcmesh_list(n), q%local(n)%val, q%local(n)%val)
  end do
  !----

  call final()

contains
!OCL SERIAL
  subroutine init()
    implicit none

    integer :: comm, myrank, nprocs
    logical :: ismaster
    
    integer :: domid
    integer :: ke, p
    !---------

    call UNITTEST_init( "test_field_cubedspheredom2d",          &
      comm=comm, myrank=myrank, nprocs=nprocs, ismaster=ismaster, use_const=.true., use_time=.true. )

    !------
    call refElem%Init(PolyOrder, .true.)

    call mesh%Init( &
      NeGX, NeGY,                                        &
      RPlanet, refElem, NLocalMeshPerPrc, nprocs, myrank )
    
    call mesh%Generate()

    !---
    call q%Init( q_VARNAME, q_UNITS, mesh )
    call lon%Init( "lon", "rad", mesh )
    call lat%Init( "lat", "rad", mesh )
    call fields_comm%Init(1, 0, 0, mesh)
    
    call FILE_HISTORY_meshfield_setup( meshCubedSphere2D_ = mesh )
    call FILE_HISTORY_reg( q_VARNAME, q_DESC, q_UNITS, HST_ID(1), dim_type='XY')
    call FILE_HISTORY_reg( lon%varname, "longitude", lon%unit, HST_ID(2), dim_type='XY')
    call FILE_HISTORY_reg( lat%varname,  "latitude", lat%unit, HST_ID(3), dim_type='XY')

    !---
    do domid=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(domid)

      !$omp parallel do private(p)
      do ke=lcmesh%NeS, lcmesh%NeE
        do p=1, refElem%Np
          q%local(domid)%val(p,ke) = TEST_field_val(lcmesh%tileID, ke, p)
        end do
        lon%local(domid)%val(:,ke) = lcmesh%lon(:,ke)
        lat%local(domid)%val(:,ke) = lcmesh%lat(:,ke)
      end do
    end do

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), lon)
    call FILE_HISTORY_meshfield_put(HST_ID(3), lat)
    call FILE_HISTORY_meshfield_write()    

    !---

    return
  end subroutine init

!OCL SERIAL
  subroutine final()
    implicit none

    call FILE_HISTORY_meshfield_finalize()

    call q%Final()
    call lon%Final()
    call lat%Final()
    
    call fields_comm%Final()
    call mesh%Final()
    call refElem%Final()
    
    call UNITTEST_final()

    return
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
    implicit none

    integer, intent(in) :: n
    type(LocalMesh2D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ke, f, p
    integer :: fnid
    integer :: ans(refElem%Np)
    !------------------------------

    do ke=lcmesh_%NeS, lcmesh_%NeE
      do p=1, refElem%Np
        ans(p) = TEST_field_val(lcmesh_%tileID,ke,p)
      end do

      call UNITTEST_assert_equal( int(lcfield(:,ke)), ans(:), 'check_interior_data: q', idx=ke )
      do f=1, refElem%Nfaces
      do p=1, refElem%Nfp
        fnid = p + (f-1)*refElem%Nfp
        call TEST_check_face_node( lcmesh_, lcfield_fl, ke, f, fnid )
      end do
      end do
    end do

    return
  end subroutine check_interior_data

  !> Check the values in halo 
  !! The halo data are stored in the order of south, east, north, and west boundaries. 
  !! The answer is obtained from the neighbor tile, considering the face ID of the neighbor tile 
  !! and the orientation of the face (tileFaceID_globalMap). 
  subroutine check_halo_data(n, lcmesh_, lcfield, lcfield_fl)
    implicit none

    integer, intent(in) :: n
    type(LocalMesh2D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ans(refElem%Nfp,max(lcmesh_%NeX,lcmesh_%NeY))
    integer :: haloInd_s, haloInd_e
    integer :: i, fph, f
    integer :: ke, fp
    integer :: tileID
    integer :: s_face(4)
    integer :: Ne1D(4)
    character(len=5), parameter :: face_name(4) = (/ 'south', 'east ', 'north', 'west ' /)
    !------------------------------

    tileID = lcmesh_%tileID
    s_face(:) = mesh%tileFaceID_globalMap(:,tileID)
    Ne1D(:) = (/ lcmesh_%NeX, lcmesh_%NeY, lcmesh_%NeX, lcmesh_%NeY /)

    write(*,*) "Check values in halo.."
    write(*,*) "tileID_globalMap=", mesh%tileID_globalMap(:,tileID)
    write(*,*) "tileFaceID_globalMap=", s_face(:)

    haloInd_e = refElem%Np*lcmesh_%Ne
    do f=1, 4
      haloInd_s = haloInd_e + 1
      haloInd_e = haloInd_s + refElem%Nfp*Ne1D(f) - 1
      do i=1, Ne1D(f)
      do fph=1, refElem%Nfp
        call TEST_get_face_index( s_face(f), Ne1D(f), refElem%Nfp, i, fph, &
          ke, fp )
        ans(fph,i) = TEST_field_val(mesh%tileID_globalMap(f,tileID), ke, refElem%Fmask(fp,abs(s_face(f))))
      end do
      end do
      call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans(:,1:Ne1D(f)), &
        'check_halo_data: q_halo_'//trim(face_name(f)), idx=tileID )
    end do

    return
  end subroutine check_halo_data

end program test_field_cubedspheredom2d
