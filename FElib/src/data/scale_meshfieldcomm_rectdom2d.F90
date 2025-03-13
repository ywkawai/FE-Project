!-------------------------------------------------------------------------------
!> module FElib / Data / Communication 2D rectangle domain
!!
!! @par Description
!!      A module to manage data communication with 2D rectangle domain for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
#include "scaleFElib.h"
module scale_meshfieldcomm_rectdom2d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: &
    ElementBase, ElementBase2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase,                               &
    MeshFieldCommBase_Init, MeshFieldCommBase_Final, &
    MeshFieldCommBase_extract_bounddata,             &
    MeshFieldCommBase_set_bounddata,                 &
    MeshFieldContainer
  use scale_localmesh_2d, only: Localmesh2d
   
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  type, public, extends(MeshFieldCommBase) :: MeshFieldCommRectDom2D
    class(MeshRectDom2D), pointer :: mesh2d
  contains
    procedure, public :: Init => MeshFieldCommRectDom2D_Init
    procedure, public :: Put => MeshFieldCommRectDom2D_put
    procedure, public :: Get => MeshFieldCommRectDom2D_get
    procedure, public :: Exchange => MeshFieldCommRectDom2D_exchange  
    procedure, public :: Final => MeshFieldCommRectDom2D_Final
  end type MeshFieldCommRectDom2D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: push_localsendbuf
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer :: bufsize_per_field
  integer, parameter :: COMM_FACE_NUM = 4

contains
  subroutine MeshFieldCommRectDom2D_Init( this, &
    sfield_num, hvfield_num, htensorfield_num, mesh2d )
    implicit none
    
    class(MeshFieldCommRectDom2D), intent(inout) :: this
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    integer, intent(in) :: htensorfield_num
    class(MeshRectDom2D), intent(in), target :: mesh2d
    
    type(LocalMesh2D), pointer :: lcmesh
    integer :: n
    integer :: Nnode_LCMeshFace(COMM_FACE_NUM,mesh2d%LOCAL_MESH_NUM)
    !-----------------------------------------------------------------------------
    
    this%mesh2d => mesh2d
    lcmesh => mesh2d%lcmesh_list(1)
    bufsize_per_field = 2*(lcmesh%NeX + lcmesh%NeY)*lcmesh%refElem2D%Nfp

    do n=1, this%mesh2d%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      Nnode_LCMeshFace(:,n) = (/ lcmesh%NeX, lcmesh%NeY, lcmesh%NeX, lcmesh%NeY /) * lcmesh%refElem2D%Nfp
    end do

    call MeshFieldCommBase_Init( this, sfield_num, hvfield_num, htensorfield_num, bufsize_per_field, COMM_FACE_NUM, Nnode_LCMeshFace, mesh2d)  
  
    return
  end subroutine MeshFieldCommRectDom2D_Init

  subroutine MeshFieldCommRectDom2D_Final( this )
    implicit none
    
    class(MeshFieldCommRectDom2D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call MeshFieldCommBase_Final( this )

    return
  end subroutine MeshFieldCommRectDom2D_Final

  subroutine MeshFieldCommRectDom2D_put(this, field_list, varid_s)
    implicit none
    class(MeshFieldCommRectDom2D), intent(inout) :: this
    type(MeshFieldContainer), intent(in) :: field_list(:)
    integer, intent(in) :: varid_s
  
    integer :: i
    integer :: n
    type(Localmesh2d), pointer :: lcmesh
    !-----------------------------------------------------------------------------

    do i=1, size(field_list)
    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      call MeshFieldCommBase_extract_bounddata( field_list(i)%field2d%local(n)%val, lcmesh%refElem, lcmesh, & ! (in)
        this%send_buf(:,varid_s+i-1,n) )                                                                  ! (out)
    end do
    end do

    return
  end subroutine MeshFieldCommRectDom2D_put

  subroutine MeshFieldCommRectDom2D_get(this, field_list, varid_s)
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_wait_core
    implicit none
    
    class(MeshFieldCommRectDom2D), intent(inout) :: this
    type(MeshFieldContainer), intent(inout) :: field_list(:)
    integer, intent(in) :: varid_s

    integer :: i
    integer :: n
    type(Localmesh2d), pointer :: lcmesh
    !-----------------------------------------------------------------------------

    !--
    if ( this%call_wait_flag_sub_get ) &
      call MeshFieldCommBase_wait_core( this, this%commdata_list )

    do i=1, size(field_list) 
    do n=1, this%mesh2d%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      call MeshFieldCommBase_set_bounddata( this%recv_buf(:,varid_s+i-1,n), lcmesh%refElem, lcmesh, & !(in)
         field_list(i)%field2d%local(n)%val )                                                         !(out)
    end do
    end do

    return
  end subroutine MeshFieldCommRectDom2D_get

!OCL SERIAL
  subroutine MeshFieldCommRectDom2D_exchange( this, do_wait )
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_exchange_core,  &
      LocalMeshCommData

    implicit none
  
    class(MeshFieldCommRectDom2D), intent(inout), target :: this
    logical, intent(in), optional :: do_wait

    integer :: n, f
    type(LocalMeshCommData), pointer :: commdata
    !-----------------------------------------------------------------------------

    do n=1, this%mesh%LOCAL_MESH_NUM
    do f=1, this%nfaces_comm
      commdata => this%commdata_list(f,n)
      call push_localsendbuf( commdata%send_buf(:,:),             &  ! (inout)
        this%send_buf(:,:,n), commdata%s_faceID, this%is_f(f,n),  &  ! (in)
        commdata%Nnode_LCMeshFace, this%field_num_tot)               ! (in)
    end do
    end do

    !-----------------------

    call MeshFieldCommBase_exchange_core( this, this%commdata_list(:,:), do_wait )

    !---------------------

    return
  end subroutine MeshFieldCommRectDom2D_exchange

!----------------------------

!OCL SERIAL
  subroutine push_localsendbuf( lc_send_buf, send_buf, s_faceID, is, Nnode_LCMeshFace, var_num )
    implicit none

    integer, intent(in) :: var_num
    integer, intent(in) ::  Nnode_LCMeshFace
    real(RP), intent(inout) :: lc_send_buf(Nnode_LCMeshFace,var_num)
    real(RP), intent(in) :: send_buf(bufsize_per_field,var_num)  
    integer, intent(in) :: s_faceID, is
  
    integer :: ie
    !-----------------------------------------------------------------------------

    ie = is + Nnode_LCMeshFace - 1
    if ( s_faceID > 0 ) then
      lc_send_buf(:,:) = send_buf(is:ie,:)   
    else
      lc_send_buf(:,:) = send_buf(ie:is:-1,:)   
    end if 
    
    return
  end subroutine push_localsendbuf

end module scale_meshfieldcomm_rectdom2d