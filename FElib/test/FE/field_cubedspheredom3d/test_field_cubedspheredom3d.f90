#include "scalelib.h"
program test_field_cubedspheredom3d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale_file_history
  use scale
  use scale_const, only: &
    RPlanet => CONST_Radius

  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D

  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_cubedspheredom3d, only: &
    MeshFieldCommCubedSphereDom3D

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
    TEST_field_val, TEST_perform_comm, TEST_check_face_node, TEST_get_face_index
  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: UNITTEST_assert_equal
  implicit none

  integer, parameter :: PolyOrder_h = 2
  integer, parameter :: PolyOrder_v = 1
  integer, parameter :: NeGX = 3
  integer, parameter :: NeGY = 3
  integer, parameter :: NeGZ = 2
  real(RP), parameter :: dom_zmin = 0.0_RP
  real(RP), parameter :: dom_zmax = 100.0_RP
  integer, parameter :: NLocalMeshPerPrc = 6

  type(HexahedralElement) :: refElem  
  type(MeshCubedSphereDom3D), target :: mesh

  type(MeshField3D), target :: q
  character(len=H_short), parameter :: q_VARNAME = "q"
  character(len=H_short), parameter :: q_DESC    = "q"
  character(len=H_short), parameter :: q_UNITS   = "K"
  integer, save :: HST_ID(1)

  type(MeshFieldCommCubedSphereDom3D) :: fields_comm
  type(MeshFieldContainer), save :: field_list(1)  

  integer :: n
  type(LocalMesh3D), pointer :: lcmesh

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

  subroutine init()
    implicit none

    integer :: comm, myrank, nprocs
    logical :: ismaster
    
    integer :: k, p
    !---------

    call UNITTEST_init( "test_field_cubedspheredom3d",          &
      comm=comm, myrank=myrank, nprocs=nprocs, ismaster=ismaster, use_const=.true., use_time=.true. )

    !------
    call refElem%Init(PolyOrder_h, PolyOrder_v, .true.)

    call mesh%Init( &
      NeGX, NeGY, NeGZ, RPlanet, &
      dom_zmin, dom_zmax,        &
      refElem, NLocalMeshPerPrc  )
    
    call mesh%Generate()

    !---
    call q%Init( q_VARNAME, q_UNITS, mesh )
    call fields_comm%Init(1, 0, 0, mesh)
    
    call FILE_HISTORY_meshfield_setup( meshCubedSphere3D_ = mesh )
    call FILE_HISTORY_reg( q_VARNAME, q_DESC, q_UNITS, HST_ID(1), dim_type='XYZ')

    !---

    !---
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE
      do p=1, refElem%Np
        q%local(n)%val(p,k) = TEST_field_val(lcmesh%tileID, k, p)
      end do
      end do
    end do

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_write()    

    !---

    return
  end subroutine init

  subroutine final()
    implicit none

    call FILE_HISTORY_meshfield_finalize()

    call q%Final()
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

    field_list(1)%field3d => q
    call TEST_perform_comm( fields_comm, field_list )

    return
  end subroutine perform_comm
    
  subroutine check_interior_data(n, lcmesh_, lcfield, lcfield_fl)
    implicit none

    integer, intent(in) :: n
    type(LocalMesh3D), intent(in) :: lcmesh_
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
      do f=1, refElem%Nfaces_h
      do p=1, refElem%Nfp_h
        fnid = p + (f-1)*refElem%Nfp_h
        call TEST_check_face_node( lcmesh_, lcfield_fl, ke, f, fnid )
      end do
      end do
      do f=1, refElem%Nfaces_v
      do p=1, refElem%Nfp_v
        fnid = p + refElem%Nfaces_h*refElem%Nfp_h + (f-1)*refElem%Nfp_v
        call TEST_check_face_node( lcmesh_, lcfield_fl, ke, refElem%Nfaces_h+f, fnid )
      end do
      end do
    end do

    return
  end subroutine check_interior_data

  subroutine check_halo_data(n, lcmesh_, lcfield, lcfield_fl)
    implicit none

    integer, intent(in) :: n
    type(LocalMesh3D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ans_bx(refElem%Nfp_h,lcmesh_%NeX,lcmesh_%NeZ)
    integer :: ans_by(refElem%Nfp_h,lcmesh_%NeY,lcmesh_%NeZ)    
    integer :: ans_bz(refElem%Nfp_v,lcmesh_%NeX,lcmesh_%NeY) 
    integer :: haloInd_s, haloInd_e
    integer :: ke, f
    integer :: i, j, k, fp, fph, fpv
    integer :: tileID

    integer :: Ne, NeX, NeY, NeZ
    integer :: Np, Nfp_h, Nfp_v, Nnode_h1D

    integer :: s_face(6)

    !------------------------------

    tileID = lcmesh_%tileID
    Ne  = lcmesh_%Ne
    NeX = lcmesh_%NeX
    NeY = lcmesh_%NeY
    NeZ = lcmesh_%NeZ
    Np = refElem%Np
    Nfp_h = refElem%Nfp_h
    Nfp_v = refElem%Nfp_v
    Nnode_h1D = refElem%Nnode_h1D

    s_face(:) = mesh%tileFaceID_globalMap(:,tileID)
    write(*,*) "Check values in halo.."
    write(*,*) "tileID_globalMap=", mesh%tileID_globalMap(:,tileID)
    write(*,*) "tileFaceID_globalMap=", s_face(:)

    haloInd_s = Np*Ne + 1    
    haloInd_e = haloInd_s + Nfp_h*NeX*NeZ - 1
    do k=1, NeZ
    do i=1, NeX
      do fpv=1, refElem%Nnode_v
      do fph=1, refElem%Nnode_h1D 
        call TEST_get_face_index( s_face(1), NeX, Nnode_h1D, i, fph, &
          ke, fp, k=k, fpv=fpv )
        ans_bx(fph+(fpv-1)*Nnode_h1D,i,k) = &
          TEST_field_val(mesh%tileID_globalMap(1,tileID), ke, refElem%Fmask_h(fp,abs(s_face(1))))
      end do
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_bx(:,:,:), &
      'check_halo_data: q_halo_south' )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + Nfp_h*NeY*NeZ - 1
    do k=1, NeZ
    do j=1, NeY
      do fpv=1, refElem%Nnode_v
      do fph=1, refElem%Nnode_h1D 
        call TEST_get_face_index( s_face(2), NeY, Nnode_h1D, j, fph, &
          ke, fp, k=k, fpv=fpv )
        ans_by(fph+(fpv-1)*Nnode_h1D,j,k) = &
          TEST_field_val(mesh%tileID_globalMap(2,tileID), ke, refElem%Fmask_h(fp,abs(s_face(2))))
      end do
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_by(:,:,:), &
      'check_halo_data: q_halo_east' )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + Nfp_h*NeX*NeZ - 1    
    do k=1, NeZ
    do i=1, NeX
      do fpv=1, refElem%Nnode_v
      do fph=1, refElem%Nnode_h1D 
        call TEST_get_face_index( s_face(3), NeX, Nnode_h1D, i, fph, &
          ke, fp, k=k, fpv=fpv )
        ans_bx(fph+(fpv-1)*Nnode_h1D,i,k) = &
          TEST_field_val(mesh%tileID_globalMap(3,tileID), ke, refElem%Fmask_h(fp,abs(s_face(3))))
      end do
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_bx(:,:,:), &
      'check_halo_data: q_halo_north' )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + Nfp_h*NeY*NeZ - 1
    do k=1, NeZ
    do j=1, NeY
      do fpv=1, refElem%Nnode_v
      do fph=1, refElem%Nnode_h1D 
        call TEST_get_face_index( s_face(4), NeY, Nnode_h1D, j, fph, &
          ke, fp, k=k, fpv=fpv )
        ans_by(fph+(fpv-1)*Nnode_h1D,j,k) = &
          TEST_field_val(mesh%tileID_globalMap(4,tileID), ke, refElem%Fmask_h(fp,abs(s_face(4))))
      end do
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_by(:,:,:), &
      'check_halo_data: q_halo_west' )    

    haloInd_s = haloInd_e + 1   
    haloInd_e = haloInd_s + Nfp_v*NeX*NeY - 1
    do j=1, NeY
    do i=1, NeX
     ke = i + (j-1)*NeX
      do fp=1, Nfp_v
        ans_bz(fp,i,j) = TEST_field_val(mesh%tileID_globalMap(5,tileID), ke, refElem%Fmask_v(fp,1))
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_bz(:,:,:), &
      'check_halo_data: q_halo_bottom' )

    haloInd_s = haloInd_e + 1   
    haloInd_e = haloInd_s + Nfp_v*NeX*NeY - 1
    do j=1, NeY
    do i=1, NeX
      ke = i + (j-1)*NeX + (NeZ-1)*NeX*NeY
      do fp=1, Nfp_v
        ans_bz(fp,i,j) = TEST_field_val(mesh%tileID_globalMap(6,tileID), ke, refElem%Fmask_v(fp,2))
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_bz(:,:,:), &
      'check_halo_data: q_halo_top' ) 

    return
  end subroutine check_halo_data


end program test_field_cubedspheredom3d
