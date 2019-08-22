#include "scaleFElib.h"
module scale_meshfieldcomm_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: elementbase
  use scale_mesh_base, only: MeshBase
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField1D, MeshField2D
  use scale_localmesh_base, only: LocalMeshBase

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 


  type, abstract, public :: MeshFieldCommBase
    integer :: sfield_num
    integer :: hvfield_num
    integer :: field_num_tot

    class(MeshBase), pointer :: mesh
    real(RP), allocatable :: send_buf(:,:,:)
    real(RP), allocatable :: recv_buf(:,:,:)
    integer :: nfaces_comm

  contains
    procedure(MeshFieldCommBase_put), public, deferred :: Put   
    procedure(MeshFieldCommBase_get), public, deferred :: Get
    procedure(MeshFieldCommBase_exchange), public, deferred :: Exchange
  end type MeshFieldCommBase

  public :: MeshFieldCommBase_Init, MeshFieldCommBase_Final
  public :: MeshFieldCommBase_exchange_core
  public :: MeshFieldCommBase_extract_bounddata
  public :: MeshFieldCommBase_set_bounddata

  type, public :: MeshFieldContainer
    class(MeshField1D), pointer :: field1d
    class(MeshField2D), pointer :: field2d
    !class(MeshFieldBase), pointer :: field3d        
  end type

  type, public :: LocalMeshCommData
    class(LocalMeshBase), pointer :: lcmesh
    real(RP), allocatable :: send_buf(:,:)
    real(RP), allocatable :: recv_buf(:,:)
    integer :: Nnode_LCMeshFace
    integer :: s_tileID
    integer :: s_faceID
    integer :: s_panelID
    integer :: s_rank
    integer :: s_tilelocalID
    integer :: faceID
  contains
    procedure, public :: Init => LocalMeshCommData_Init
    procedure, public :: Final => LocalMeshCommData_Final
    procedure, public :: ISendRecv => LocalMeshCommData_iSendRecv
  end type  
  
  interface
    subroutine MeshFieldCommBase_put(this, field_list, varid_s)
      import MeshFieldCommBase
      import MeshFieldContainer
      class(MeshFieldCommBase), intent(inout) :: this
      type(MeshFieldContainer), intent(in) :: field_list(:)
      integer, intent(in) :: varid_s
    end subroutine MeshFieldCommBase_put

    subroutine MeshFieldCommBase_get(this, field_list, varid_s)
      import MeshFieldCommBase
      import MeshFieldContainer
      class(MeshFieldCommBase), intent(in) :: this
      type(MeshFieldContainer), intent(inout) :: field_list(:)
      integer, intent(in) :: varid_s
    end subroutine MeshFieldCommBase_get

    subroutine MeshFieldComm_get_vec( this, &
      & varid_s, U, V, u1, u2 )

      import MeshFieldCommBase
      import MeshFieldBase      
      class(MeshFieldCommBase), intent(in) :: this
      integer, intent(in) :: varid_s
      class(MeshFieldBase), intent(inout) :: U
      class(MeshFieldBase), intent(inout) :: V
      class(MeshFieldBase), intent(inout), optional :: u1
      class(MeshFieldBase), intent(inout), optional :: u2
    end subroutine MeshFieldComm_get_vec

    subroutine MeshFieldCommBase_Exchange(this)
        import MeshFieldCommBase
        class(MeshFieldCommBase), intent(inout) :: this
    end subroutine MeshFieldCommBase_Exchange

    subroutine MeshFieldComm_Final(this)
      import MeshFieldCommBase  
      class(MeshFieldCommBase), intent(inout) :: this
    end subroutine MeshFieldComm_Final    
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
  subroutine MeshFieldCommBase_Init( this, &
    & sfield_num, hvfield_num, bufsize_per_field, comm_face_num, mesh )

    implicit none
    
    class(MeshFieldCommBase), intent(inout) :: this
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    integer, intent(in) :: bufsize_per_field
    integer, intent(in) :: comm_face_num
    class(Meshbase), intent(in), target :: mesh

    this%mesh => mesh
    this%sfield_num = sfield_num
    this%hvfield_num = hvfield_num
    this%field_num_tot = sfield_num + hvfield_num*2
    this%nfaces_comm = comm_face_num

    if (this%field_num_tot > 0) then
      allocate( this%send_buf(bufsize_per_field, this%field_num_tot, mesh%LOCAL_MESH_NUM) )
      allocate( this%recv_buf(bufsize_per_field, this%field_num_tot, mesh%LOCAL_MESH_NUM) )
    end if 

    return
  end subroutine MeshFieldCommBase_Init

  subroutine MeshFieldCommBase_Final( this )

    implicit none
    
    class(MeshFieldCommBase), intent(inout) :: this

    if (this%field_num_tot > 0) then
      deallocate( this%send_buf, this%recv_buf )
    end if

    return
  end subroutine MeshFieldCommBase_Final

  subroutine MeshFieldCommBase_exchange_core( this, commdata_list )

    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD, PRC_abort, PRC_MPIbarrier
    use mpi, only: &
      !MPI_waitall,    &
      MPI_STATUS_SIZE
    implicit none
  
    class(MeshFieldCommBase), intent(inout) :: this
    type(LocalMeshCommData), intent(inout), target :: commdata_list(this%nfaces_comm, this%mesh%LOCAL_MESH_NUM)
  
    integer :: n, f, k, l, p

    integer :: var_id
    integer :: irs, ire
    integer :: ierr
    integer :: req_counter 
    integer :: request_send(this%nfaces_comm * this%mesh%LOCAL_MESH_NUM)
    integer :: request_recv(this%nfaces_comm * this%mesh%LOCAL_MESH_NUM)
    integer :: stat_send(MPI_STATUS_SIZE, this%nfaces_comm * this%mesh%LOCAL_MESH_NUM)
    integer :: stat_recv(MPI_STATUS_SIZE, this%nfaces_comm * this%mesh%LOCAL_MESH_NUM)
  
    !-----------------------------------------------------------------------------
    
    req_counter = 0
    do n=1, this%mesh%LOCAL_MESH_NUM      
    do f=1, this%nfaces_comm
      call commdata_list(f,n)%ISendRecv( &
        req_counter, request_send(:), request_recv(:), & ! (inout)
        commdata_list(:,:) )                             ! (inout)      
    end do
    end do

    if (req_counter > 0) then
      call MPI_waitall(req_counter, request_recv(1:req_counter), stat_recv(:,1:req_counter), ierr)
      call MPI_waitall(req_counter, request_send(1:req_counter), stat_send(:,1:req_counter), ierr)
    end if
  
    !---------------------

    do n=1, this%mesh%LOCAL_MESH_NUM
      irs = 1
      do f=1, this%nfaces_comm
        ire = irs + commdata_list(f,n)%Nnode_LCMeshFace - 1
        do var_id=1, this%field_num_tot
          this%recv_buf(irs:ire,var_id,n) = commdata_list(f,n)%recv_buf(:,var_id)
        end do
        irs = ire + 1
      end do ! end loop for face
    end do

    return
  end subroutine MeshFieldCommBase_exchange_core

  subroutine MeshFieldCommBase_extract_bounddata(var, refElem, mesh, buf)
    implicit none
  
    class(ElementBase), intent(in) :: refElem
    class(LocalMeshBase), intent(in) :: mesh
    real(RP), intent(in) :: var(refElem%Np * mesh%NeA)
    real(RP), intent(inout) :: buf(size(mesh%VmapB))
  
    buf(:) = var(mesh%VmapB(:))
    return
  end subroutine MeshFieldCommBase_extract_bounddata

  subroutine MeshFieldCommBase_set_bounddata(buf, refElem, mesh, var)
    implicit none
    
    class(ElementBase), intent(in) :: refElem
    class(LocalMeshBase), intent(in) :: mesh
    real(RP), intent(in) :: buf(size(mesh%VmapB))
    real(RP), intent(inout) :: var(refElem%Np * mesh%NeA)

    var(refElem%Np*mesh%NeE+1:refElem%Np*mesh%NeE+size(buf)) = buf(:)
    return
  end subroutine MeshFieldCommBase_set_bounddata  

  !-------------------------------------------
  
  subroutine LocalMeshCommData_Init( this, comm, lcmesh, faceID, Nnode_LCMeshFace )
    class(LocalMeshCommData), intent(inout) :: this
    class(MeshFieldCommBase), intent(in) :: comm
    class(LocalMeshBase), intent(in), target :: lcmesh
    integer, intent(in) :: faceID
    integer, intent(in) :: Nnode_LCMeshFace

    this%lcmesh => lcmesh
    this%Nnode_LCMeshFace = Nnode_LCMeshFace

    allocate( this%send_buf(Nnode_LCMeshFace, comm%field_num_tot))
    allocate( this%recv_buf(Nnode_LCMeshFace, comm%field_num_tot))

    this%s_tileID  = comm%mesh%tileID_globalMap(faceID, lcmesh%tileID)
    this%s_faceID  = comm%mesh%tileFaceID_globalMap(faceID, lcmesh%tileID)
    this%s_panelID = comm%mesh%tilePanelID_globalMap(faceID, lcmesh%tileID)      
    this%s_rank    = comm%mesh%PRCrank_globalMap(this%s_tileID)
    this%s_tilelocalID = comm%mesh%tileID_global2localMap(this%s_tileID)
    this%faceID    = faceID

    return
  end subroutine LocalMeshCommData_Init

  subroutine LocalMeshCommData_iSendRecv( this, &
      req_counter, req_send, req_recv,          &
      lccommdat_list )
    
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD
    use mpi, only: &
      !MPI_isend, MPI_irecv, &
      MPI_DOUBLE_PRECISION
    implicit none

    class(LocalMeshCommData), intent(inout) :: this
    integer, intent(inout) :: req_counter
    integer, intent(inout) :: req_send(:)
    integer, intent(inout) :: req_recv(:)
    type(LocalMeshCommData), intent(inout) :: lccommdat_list(:,:)
  
    integer :: tag
    integer :: bufsize
    integer :: ierr

    !-------------------------------------------

    if ( this%s_rank /= this%lcmesh%PRC_myrank ) then

      req_counter = req_counter + 1

      bufsize = size(this%send_buf)
      tag = 10 * this%s_tileID + abs(this%s_faceID)
      call MPI_isend( this%send_buf(1,1), bufsize, MPI_DOUBLE_PRECISION, &
       this%s_rank, tag, PRC_LOCAL_COMM_WORLD,                           &
       req_send(req_counter), ierr )

      tag = 10 * this%lcmesh%tileID + this%faceID
      bufsize = size(this%recv_buf)
      call MPI_irecv( this%recv_buf(1,1), bufsize, MPI_DOUBLE_PRECISION, &
        this%s_rank, tag, PRC_LOCAL_COMM_WORLD,                          &
        req_recv(req_counter), ierr)
      
    else  !<----- end if s_rank /= mesh%PRC_myrank
      lccommdat_list(abs(this%s_faceID), this%s_tilelocalID)%recv_buf(:,:) &
        = this%send_buf(:,:)
    end if         
    
    return
  end subroutine LocalMeshCommData_isendrecv

  subroutine LocalMeshCommData_Final( this )
    class(LocalMeshCommData), intent(inout) :: this
    
    deallocate( this%send_buf )
    deallocate( this%recv_buf )

    return
  end subroutine LocalMeshCommData_Final

end module scale_meshfieldcomm_base