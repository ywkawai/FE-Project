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
    MeshFieldCommBase_extract_bounddata_2,           &
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

  !> Base derived type to manage data communication with 2D rectangle domain
  type, public, extends(MeshFieldCommBase) :: MeshFieldCommRectDom2D
    class(MeshRectDom2D), pointer :: mesh2d  !< Pointer to an object representing 2D rectangular computational mesh
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
  integer, parameter :: COMM_FACE_NUM = 4  !< Number of faces with data communication


contains
!> Initialize an object to manage data communication with 2D rectangle domain
  subroutine MeshFieldCommRectDom2D_Init( this, &
    sfield_num, hvfield_num, htensorfield_num, mesh2d )
    implicit none
    
    class(MeshFieldCommRectDom2D), intent(inout) :: this
    integer, intent(in) :: sfield_num                    !< Number of scalar fields
    integer, intent(in) :: hvfield_num                   !< Number of horizontal vector fields
    integer, intent(in) :: htensorfield_num              !< Number of horizontal vector fields
    class(MeshRectDom2D), intent(in), target :: mesh2d   !< Object to manage a 2D rectangular computational mesh
    
    type(LocalMesh2D), pointer :: lcmesh
    integer :: n
    integer :: Nnode_LCMeshFace(COMM_FACE_NUM,mesh2d%LOCAL_MESH_NUM)
    !-----------------------------------------------------------------------------
    
    this%mesh2d => mesh2d
    lcmesh => mesh2d%lcmesh_list(1)
    this%bufsize_per_field = 2*(lcmesh%NeX + lcmesh%NeY)*lcmesh%refElem2D%Nfp

    do n=1, this%mesh2d%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      Nnode_LCMeshFace(:,n) = (/ lcmesh%NeX, lcmesh%NeY, lcmesh%NeX, lcmesh%NeY /) * lcmesh%refElem2D%Nfp
    end do

    call MeshFieldCommBase_Init( this, sfield_num, hvfield_num, htensorfield_num, this%bufsize_per_field, COMM_FACE_NUM, Nnode_LCMeshFace, mesh2d)  
  
    return
  end subroutine MeshFieldCommRectDom2D_Init

!> Finalize an object to manage data communication with 2D rectangle domain
  subroutine MeshFieldCommRectDom2D_Final( this )
    implicit none
    
    class(MeshFieldCommRectDom2D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call MeshFieldCommBase_Final( this )

    return
  end subroutine MeshFieldCommRectDom2D_Final

!> Put field data into temporary buffers  
  subroutine MeshFieldCommRectDom2D_put(this, field_list, varid_s)
    implicit none
    class(MeshFieldCommRectDom2D), intent(inout) :: this
    type(MeshFieldContainer), intent(in) :: field_list(:)
    integer, intent(in) :: varid_s
  
    integer :: i
    integer :: n
    type(LocalMesh2D), pointer :: lcmesh
    !-----------------------------------------------------------------------------

    ! do i=1, size(field_list)
    ! do n=1, this%mesh%LOCAL_MESH_NUM
    !   lcmesh => this%mesh2d%lcmesh_list(n)
    !   call MeshFieldCommBase_extract_bounddata( field_list(i)%field2d%local(n)%val, lcmesh%refElem, lcmesh, & ! (in)
    !     this%send_buf(:,varid_s+i-1,n) )                                                                      ! (out)
    ! end do
    ! end do
    call MeshFieldCommBase_extract_bounddata_2( &
      field_list, 2, varid_s, this%mesh2d%lcmesh_list, size(this%mesh2d%lcmesh_list(1)%VMapB), & !(in)
      this%send_buf ) ! (out)

    return
  end subroutine MeshFieldCommRectDom2D_put

 !> Extract field data from temporary buffers 
  subroutine MeshFieldCommRectDom2D_get(this, field_list, varid_s)
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_wait_core
    implicit none
    
    class(MeshFieldCommRectDom2D), intent(inout) :: this
    type(MeshFieldContainer), intent(inout) :: field_list(:)
    integer, intent(in) :: varid_s

    integer :: i
    integer :: n
    type(LocalMesh2D), pointer :: lcmesh
    !-----------------------------------------------------------------------------

    !--
    if ( this%call_wait_flag_sub_get ) then
      call MeshFieldCommBase_wait_core( this, this%commdata_list, &
        field_list, 2, varid_s, this%mesh2d%lcmesh_list )
    else
      do i=1, size(field_list) 
      do n=1, this%mesh2d%LOCAL_MESH_NUM
        lcmesh => this%mesh2d%lcmesh_list(n)
        call MeshFieldCommBase_set_bounddata( this%recv_buf(:,varid_s+i-1,n), lcmesh%refElem, lcmesh, & !(in)
          field_list(i)%field2d%local(n)%val )                                                         !(out)
      end do
      end do
    end if
    return
  end subroutine MeshFieldCommRectDom2D_get

!> Exchange field data between neighboring MPI processes
!!
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
      call push_localsendbuf( commdata%send_buf,                  &  ! (inout)
        this%send_buf(:,:,n), commdata%s_faceID, this%is_f(f,n),  &  ! (in)
        commdata%Nnode_LCMeshFace, this%bufsize_per_field,        &  ! (in)
        this%field_num_tot )                                         ! (in)
    end do
    end do
    !$acc wait(1)
    !-----------------------

    call MeshFieldCommBase_exchange_core( this, this%commdata_list, do_wait )

    !---------------------

    return
  end subroutine MeshFieldCommRectDom2D_exchange

!----------------------------

!OCL SERIAL
  subroutine push_localsendbuf( lc_send_buf, send_buf, s_faceID, is, Nnode_LCMeshFace, bufsize_per_field, var_num )
    implicit none

    integer, intent(in) :: var_num
    integer, intent(in) ::  Nnode_LCMeshFace
    integer, intent(in) :: bufsize_per_field
    real(RP), intent(inout) :: lc_send_buf(Nnode_LCMeshFace,var_num)
    real(RP), intent(in) :: send_buf(bufsize_per_field,var_num)  
    integer, intent(in) :: s_faceID, is
  
    integer :: is_, ie, lincrement
    integer :: i, v
    !-----------------------------------------------------------------------------

    if ( s_faceID > 0 ) then
      is_ = is
      ie  = is + Nnode_LCMeshFace - 1
      lincrement = +1          
    else
      is_ = is + Nnode_LCMeshFace - 1
      ie  = is          
      lincrement = -1          
    end if 

#ifdef _OPENACC
    !$acc parallel loop present(lc_send_buf, send_buf) async(1)
    do v=1, var_num
    do i=1, Nnode_LCMeshFace
      lc_send_buf(i,v) = send_buf(is_+(i-1)*lincrement,v)
    end do
    end do
#else
    lc_send_buf(:,:) = send_buf(is_:ie:lincrement,:) 
#endif   
    
    return
  end subroutine push_localsendbuf

end module scale_meshfieldcomm_rectdom2d