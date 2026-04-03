!-------------------------------------------------------------------------------
!> module FElib / Data / Communication 1D
!!
!! @par Description
!!      A module to manage 1D data communication for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
#include "scaleFElib.h"
module scale_meshfieldcomm_1d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: ElementBase, ElementBase1d
  use scale_mesh_base1d, only: MeshBase1D
  use scale_meshfield_base, only: MeshField1D
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase,                               &
    MeshFieldCommBase_Init, MeshFieldCommBase_Final, &
    MeshFieldCommBase_extract_bounddata,             &
    MeshFieldCommBase_set_bounddata,                 &
    MeshFieldContainer
  use scale_localmesh_1d, only: LocalMesh1D
   
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !> Base derived type to manage data communication with 1D domain
  type, public, extends(MeshFieldCommBase) :: MeshFieldComm1D
    class(MeshBase1D), pointer :: mesh1d  !< Pointer to an object representing 1D computational mesh
  contains
    procedure, public :: Init => MeshFieldComm1D_Init
    procedure, public :: Put => MeshFieldComm1D_put
    procedure, public :: Get => MeshFieldComm1D_get
    procedure, public :: Exchange => MeshFieldComm1D_exchange  
    procedure, public :: Final => MeshFieldComm1D_Final
  end type MeshFieldComm1D

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
  integer, parameter :: COMM_FACE_NUM = 2  !< Number of faces with data communication

contains
!> Initialize an object to manage data communication with 1D domain
  subroutine MeshFieldComm1D_Init( this, &
    sfield_num, hvfield_num, mesh1d )

    implicit none
    
    class(MeshFieldComm1D), intent(inout) :: this
    integer, intent(in) :: sfield_num                    !< Number of scalar fields
    integer, intent(in) :: hvfield_num                   !< Number of vector fields
    class(MeshBase1D), intent(in), target :: mesh1d      !< Object to manage a 1D computational mesh

    integer :: n
    integer :: Nnode_LCMeshFace(COMM_FACE_NUM,mesh1d%LOCAL_MESH_NUM)
    !-----------------------------------------------------------------------------
    
    this%mesh1d => mesh1d 
    this%bufsize_per_field =  mesh1d%refElem1D%Nfp * 2

    ! Dummy
    do n=1, this%mesh1d%LOCAL_MESH_NUM
      Nnode_LCMeshFace(:,n) = (/ 1, 1 /)
    end do
    
    call MeshFieldCommBase_Init( this, sfield_num, hvfield_num, 0, this%bufsize_per_field, 2, Nnode_LCMeshFace, mesh1d)  
  
    return
  end subroutine MeshFieldComm1D_Init

!> Finalize an object to manage data communication with 1D domain  
  subroutine MeshFieldComm1D_Final( this )
    implicit none
    
    class(MeshFieldComm1D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call MeshFieldCommBase_Final( this )

    return
  end subroutine MeshFieldComm1D_Final

!> Put field data into temporary buffers
  subroutine MeshFieldComm1D_put(this, field_list, varid_s)
    implicit none
    class(MeshFieldComm1D), intent(inout) :: this
    type(MeshFieldContainer), intent(in) :: field_list(:)
    integer, intent(in) :: varid_s
  
    integer :: i
    integer :: n
    type(LocalMesh1D), pointer :: lcmesh
    !-----------------------------------------------------------------------------
    
    do i=1, size(field_list)
    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh1d%lcmesh_list(n)
      call MeshFieldCommBase_extract_bounddata( field_list(i)%field1d%local(n)%val, lcmesh%refElem, lcmesh, & ! (in)
        this%send_buf(:,varid_s+i-1,n) )                                                                      ! (out)
    end do
    end do

    return
  end subroutine MeshFieldComm1D_put

!> Extract field data from temporary buffers
  subroutine MeshFieldComm1D_get(this, field_list, varid_s)
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_wait_core
    implicit none
    
    class(MeshFieldComm1D), intent(inout) :: this
    type(MeshFieldContainer), intent(inout) :: field_list(:)
    integer, intent(in) :: varid_s

    integer :: i
    integer :: n
    type(LocalMesh1D), pointer :: lcmesh
    !-----------------------------------------------------------------------------
    
    if ( this%call_wait_flag_sub_get ) then
      call MeshFieldCommBase_wait_core( this, this%commdata_list )
    else
      do i=1, size(field_list) 
      do n=1, this%mesh1D%LOCAL_MESH_NUM
        lcmesh => this%mesh1D%lcmesh_list(n)
        call MeshFieldCommBase_set_bounddata( this%recv_buf(:,varid_s+i-1,n), lcmesh%refElem, lcmesh, & !(in)
          field_list(i)%field1d%local(n)%val )                                                          !(out)
      end do
      end do
      !$acc wait(1)
    end if

    return
  end subroutine MeshFieldComm1D_get

!> Exchange field data between neighboring MPI processes
!!
!! @param do_wait Flag whether MPI_waitall is called and move tmp data of LocalMeshCommData object to a recv buffer
!OCL SERIAL
  subroutine MeshFieldComm1D_exchange( this, do_wait )
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_exchange_core,  &
      MeshFieldCommBase_wait_core,      &
      LocalMeshCommData

    implicit none
  
    class(MeshFieldComm1D), intent(inout), target :: this
    logical, intent(in), optional :: do_wait
  
    integer :: n, f
    type(LocalMesh1D), pointer :: mesh
    type(LocalMeshCommData), pointer :: commdata
    !-----------------------------------------------------------------------------

    do n=1, this%mesh%LOCAL_MESH_NUM
    do f=1, this%nfaces_comm      
      commdata => this%commdata_list(f,n)
      call push_localsendbuf( commdata%send_buf,           & ! (inout)
        this%send_buf(:,:,n), commdata%s_faceID, f,        & ! (in)
        commdata%Nnode_LCMeshFace, this%bufsize_per_field, & ! (in)
        this%field_num_tot )                                 ! (in)
    end do
    end do
    !$acc wait(1)
    !-----------------------
    call MeshFieldCommBase_exchange_core( this, this%commdata_list, do_wait )

    return
  end subroutine MeshFieldComm1D_exchange

!----------------------------

  subroutine push_localsendbuf( lc_send_buf, send_buf, s_faceID, f, Nnode_LCMeshFace, bufsize_per_field, var_num )
    implicit none

    integer, intent(in) :: var_num
    integer, intent(in) ::  Nnode_LCMeshFace
    integer, intent(in) :: bufsize_per_field
    real(RP), intent(inout) :: lc_send_buf(Nnode_LCMeshFace,var_num)
    real(RP), intent(in) :: send_buf(bufsize_per_field,var_num)  
    integer, intent(in) :: s_faceID, f
  
    integer :: is, ie, lincrement
    integer :: i, v
    !-----------------------------------------------------------------------------
  
    if ( s_faceID > 0 ) then
      is = 1 + (f-1)*Nnode_LCMeshFace
      ie = is + Nnode_LCMeshFace - 1
      lincrement = +1          
    else
      is   = f*Nnode_LCMeshFace
      ie   = 1 + (f-1)*Nnode_LCMeshFace          
      lincrement = -1          
    end if 
#ifdef _OPENACC
    !$acc parallel loop present(lc_send_buf, send_buf) async(1)
    do v=1, var_num
    do i=1, Nnode_LCMeshFace
      lc_send_buf(i,v) = send_buf(is+(i-1)*lincrement,v)
    end do
    end do
#else
    lc_send_buf(:,:) = send_buf(is:ie:lincrement,:) 
#endif   
    return
  end subroutine push_localsendbuf

end module scale_meshfieldcomm_1d