!-------------------------------------------------------------------------------
!> Module for unit tests of FElib / utilities for mesh and field tests
!!
!! @par Description
!!      Common routines used in the unit tests of mesh and field modules (FElib/test/FE/mesh_*, field_*).
!!      This module is compiled in the build directory of each test
!!      (add "mod_test_meshfield_util.o" to OBJS in Makefile of the test).
!!
!! @author Yuta Kawai, Team SCALE
!!
#include "scaleFElib.h"
module mod_test_meshfield_util
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase, MeshFieldContainer
  use scale_unittest_assert, only: &
    UNITTEST_assert_equal, UNITTEST_assert_close

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: TEST_field_val
  public :: TEST_perform_comm
  public :: TEST_check_face_node
  public :: TEST_get_face_index
  public :: TEST_elemID3D
  public :: TEST_nodeID3D

contains

  !> Value of a test field to identify the tile, element, and node
  !! (= tileID*1000000 + ke*1000 + p)
  elemental function TEST_field_val( tileID, ke, p ) result(val)
    implicit none
    integer, intent(in) :: tileID
    integer, intent(in) :: ke
    integer, intent(in) :: p
    real(RP) :: val
    !----------------------------
    !$acc routine seq
    val = tileID*1000000 + ke*1000 + p
    return
  end function TEST_field_val

  !> Perform data communication of fields (Put -> Exchange -> Get)
  subroutine TEST_perform_comm( fields_comm, field_list )
    implicit none
    class(MeshFieldCommBase), intent(inout) :: fields_comm
    type(MeshFieldContainer), intent(inout) :: field_list(:)
    !---------------------------

    write(*,*) " - Put.."
    call fields_comm%Put(field_list, 1)
    write(*,*) " - Exchange.."
    call fields_comm%Exchange()
    write(*,*) " - Get.."
    call fields_comm%Get(field_list, 1)

    return
  end subroutine TEST_perform_comm

  !> Check the data at a face node of element ke
  !! (the field is assumed to be set by TEST_field_val)
  !! - q- (iM) is the value of element ke itself
  !! - If q+ (iP) refers to an interior element, the element must be EToE(ke,f),
  !!   q+ must be unchanged by communication, and the positions of iM and iP must coincide.
  !!   (Values at halo nodes should be checked separately.)
  subroutine TEST_check_face_node( lcmesh, lcfield_fl, ke, f, fnid )
    implicit none
    class(LocalMeshBase), intent(in) :: lcmesh
    real(RP), intent(in) :: lcfield_fl(:)      !< Field data (including halo) as a 1D array
    integer, intent(in) :: ke                  !< Element ID
    integer, intent(in) :: f                   !< Face ID
    integer, intent(in) :: fnid                !< Face node ID

    integer :: Np
    integer :: iM, iP, keM, keP, pM, pP
    real(RP), parameter :: POS_EPS = 1.0E-12_RP
    !------------------------------

    Np = lcmesh%refElem%Np
    iM = lcmesh%VMapM(fnid,ke)
    iP = lcmesh%VMapP(fnid,ke)

    keM = (iM - 1) / Np + 1
    pM  = iM - (keM - 1) * Np
    call UNITTEST_assert_equal( keM, ke, 'check_face_node: element of VMapM', idx=ke )
    call UNITTEST_assert_equal( int(lcfield_fl(iM)), int(TEST_field_val(lcmesh%tileID, keM, pM)), &
      'check_face_node: q-', idx=ke )

    if ( iP > Np * lcmesh%Ne ) return ! halo

    keP = (iP - 1) / Np + 1
    pP  = iP - (keP - 1) * Np
    call UNITTEST_assert_equal( keP, lcmesh%EToE(ke,f), 'check_face_node: element of VMapP', idx=ke )
    call UNITTEST_assert_equal( int(lcfield_fl(iP)), int(TEST_field_val(lcmesh%tileID, keP, pP)), &
      'check_face_node: q+', idx=ke )
    call UNITTEST_assert_close( lcmesh%pos_en(pP,keP,:), lcmesh%pos_en(pM,keM,:), &
      'check_face_node: positions of q- and q+', rtol=POS_EPS, atol=POS_EPS, idx=ke )

    return
  end subroutine TEST_check_face_node

  !> Get the element ID and the node ID (within the face) in the neighbor tile,
  !! which corresponds to the i1D-th element and the fph-th node along a lateral boundary of a tile
  !! s_face: face ID of the neighbor tile connected to the boundary (a negative value means the opposite orientation)
  !! For 3D meshes, k (element index in the vertical direction) and fpv (node index in the vertical direction) are specified.
  !! Note: the neighbor tile is assumed to have Ne1D x Ne1D elements in the horizontal direction.
  subroutine TEST_get_face_index( s_face, Ne1D, Nnode1D, i1D, fph, &
    ke, fp,                                                        &
    k, fpv                                                         )
    implicit none
    integer, intent(in) :: s_face
    integer, intent(in) :: Ne1D
    integer, intent(in) :: Nnode1D
    integer, intent(in) :: i1D
    integer, intent(in) :: fph
    integer, intent(out) :: ke
    integer, intent(out) :: fp
    integer, intent(in), optional :: k
    integer, intent(in), optional :: fpv
    !-----------------------------------

    select case( s_face )
    case(1)
      ke = i1D
      fp = fph
    case(-1)
      ke = Ne1D - i1D + 1
      fp = Nnode1D - fph + 1
    case(2)
      ke = Ne1D + (i1D-1)*Ne1D
      fp = fph
    case(-2)
      ke = Ne1D + (Ne1D-i1D)*Ne1D
      fp = Nnode1D - fph + 1
    case(3)
      ke = i1D + (Ne1D-1)*Ne1D
      fp = fph
    case(-3)
      ke = Ne1D - i1D + 1 + (Ne1D-1)*Ne1D
      fp = Nnode1D - fph + 1
    case(4)
      ke = 1 + (i1D-1)*Ne1D
      fp = fph
    case(-4)
      ke = 1 + (Ne1D-i1D)*Ne1D
      fp = Nnode1D - fph + 1
    case default
      ke = -1; fp = -1
      LOG_ERROR('TEST_get_face_index',*) 'Unexpected s_face=', s_face
      call PRC_abort
    end select

    if ( present(k)   ) ke = ke + (k - 1) * Ne1D**2
    if ( present(fpv) ) fp = fp + (fpv - 1) * Nnode1D

    return
  end subroutine TEST_get_face_index

  !> Element ID in a 3D local mesh
  function TEST_elemID3D( lmesh, i, j, k ) result(eid)
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    integer, intent(in) :: i, j, k
    integer :: eid
    !--------------------------------------
    eid = i + (j-1)*lmesh%NeX + (k-1)*lmesh%NeX*lmesh%NeY
    return
  end function TEST_elemID3D

  !> Node ID in a hexahedral element
  function TEST_nodeID3D( elem, i, j, k ) result(nid)
    implicit none
    class(HexahedralElement), intent(in) :: elem
    integer, intent(in) :: i, j, k
    integer :: nid
    !--------------------------------------
    nid = i + (j-1)*elem%Nnode_h1D + (k-1)*elem%Nnode_h1D**2
    return
  end function TEST_nodeID3D

end module mod_test_meshfield_util
