!-------------------------------------------------------------------------------
!> module FElib / Data / Communication base
!!
!! @par Description
!!      Base module to manage data communication for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
#include "scaleFElib.h"
module scale_meshfieldcomm_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: ElementBase
  use scale_mesh_base, only: MeshBase
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField1D, MeshField2D, MeshField3D
  use scale_localmesh_base, only: &
    LocalMeshBase

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

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
    procedure, public :: SendRecv => LocalMeshCommData_SendRecv
    procedure, public :: PC_Init => LocalMeshCommData_pc_init
  end type

  type, abstract, public :: MeshFieldCommBase
    integer :: sfield_num
    integer :: hvfield_num
    integer :: htensorfield_num
    integer :: field_num_tot
    integer :: nfaces_comm

    class(MeshBase), pointer :: mesh
    real(RP), allocatable :: send_buf(:,:,:)
    real(RP), allocatable :: recv_buf(:,:,:)
    integer, allocatable :: request_send(:)
    integer, allocatable :: request_recv(:)

    type(LocalMeshCommData), allocatable :: commdata_list(:,:)
    integer, allocatable :: is_f(:,:)

    logical :: MPI_pc_flag
    integer, allocatable :: request_pc(:)

    integer :: req_counter
    logical :: call_wait_flag_sub_get
  contains
    procedure(MeshFieldCommBase_put), public, deferred :: Put   
    procedure(MeshFieldCommBase_get), public, deferred :: Get
    procedure(MeshFieldCommBase_exchange), public, deferred :: Exchange
    procedure, public :: Prepare_PC => MeshFieldCommBase_prepare_PC
  end type MeshFieldCommBase

  public :: MeshFieldCommBase_Init, MeshFieldCommBase_Final
  public :: MeshFieldCommBase_exchange_core
  public :: MeshFieldCommBase_wait_core
  public :: MeshFieldCommBase_extract_bounddata
  public :: MeshFieldCommBase_set_bounddata

  type, public :: MeshFieldContainer
    class(MeshField1D), pointer :: field1d
    class(MeshField2D), pointer :: field2d
    class(MeshField3D), pointer :: field3d        
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
      class(MeshFieldCommBase), intent(inout) :: this
      type(MeshFieldContainer), intent(inout) :: field_list(:)
      integer, intent(in) :: varid_s
    end subroutine MeshFieldCommBase_get

    subroutine MeshFieldCommBase_Exchange(this, do_wait)
        import MeshFieldCommBase
        class(MeshFieldCommBase), intent(inout), target :: this
        logical, intent(in), optional :: do_wait
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
!OCL SERIAL
  subroutine MeshFieldCommBase_Init( this, &
    sfield_num, hvfield_num, htensorfield_num, bufsize_per_field, comm_face_num, &
    Nnode_LCMeshFace, mesh )

    implicit none
    
    class(MeshFieldCommBase), intent(inout) :: this
    class(Meshbase), intent(in), target :: mesh
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    integer, intent(in) :: htensorfield_num
    integer, intent(in) :: bufsize_per_field
    integer, intent(in) :: comm_face_num
    integer, intent(in) :: Nnode_LCMeshFace(comm_face_num,mesh%LOCAL_MESH_NUM)

    integer :: n, f
    class(LocalMeshBase), pointer :: lcmesh
    !-----------------------------------------------------------------------------

    this%mesh => mesh
    this%sfield_num = sfield_num
    this%hvfield_num = hvfield_num
    this%htensorfield_num = htensorfield_num
    this%field_num_tot = sfield_num + hvfield_num*2 + htensorfield_num*4
    this%nfaces_comm = comm_face_num

    if (this%field_num_tot > 0) then
      allocate( this%send_buf(bufsize_per_field, this%field_num_tot, mesh%LOCAL_MESH_NUM) )
      allocate( this%recv_buf(bufsize_per_field, this%field_num_tot, mesh%LOCAL_MESH_NUM) )
      allocate( this%request_send(comm_face_num*mesh%LOCAL_MESH_NUM) )
      allocate( this%request_recv(comm_face_num*mesh%LOCAL_MESH_NUM) )

      allocate( this%commdata_list(comm_face_num,mesh%LOCAL_MESH_NUM) )
      allocate( this%is_f(comm_face_num,mesh%LOCAL_MESH_NUM) ) 

      do n=1, mesh%LOCAL_MESH_NUM
        this%is_f(1,n) = 1
        do f=2, this%nfaces_comm
          this%is_f(f,n) = this%is_f(f-1,n) + Nnode_LCMeshFace(f-1,n)
        end do

        call mesh%GetLocalMesh(n, lcmesh)
        do f=1, this%nfaces_comm
          call this%commdata_list(f,n)%Init( this, lcmesh, f, Nnode_LCMeshFace(f,n) )
        end do
      end do  
    end if 

    this%MPI_pc_flag = .false.

    return
  end subroutine MeshFieldCommBase_Init

!OCL SERIAL
  subroutine MeshFieldCommBase_Final( this )
    implicit none
    
    class(MeshFieldCommBase), intent(inout) :: this

    integer :: n, f
    integer :: ireq
    integer :: ierr
    !-----------------------------------------------------------------------------

    if (this%field_num_tot > 0) then
      deallocate( this%send_buf, this%recv_buf )
      deallocate( this%request_send, this%request_recv )

      do n=1, this%mesh%LOCAL_MESH_NUM
      do f=1, this%nfaces_comm
        call this%commdata_list(f,n)%Final()
      end do
      end do     
      deallocate( this%commdata_list, this%is_f )

      if (this%MPI_pc_flag) then
        do ireq=1, this%req_counter
          call MPI_request_free( this%request_pc(ireq), ierr )
        end do           
        deallocate( this%request_pc ) 
      end if
    end if

    return
  end subroutine MeshFieldCommBase_Final

!> Prepare persistent communication
!OCL SERIAL
  subroutine MeshFieldCommBase_prepare_PC( this )
    implicit none    
    class(MeshFieldCommBase), intent(inout) :: this

    integer :: n, f
    !-----------------------------------------------------------------------------

    this%MPI_pc_flag = .true.
    if (this%field_num_tot > 0) then
      allocate( this%request_pc(2*this%nfaces_comm*this%mesh%LOCAL_MESH_NUM) )
    end if

    this%req_counter = 0
    do n=1, this%mesh%LOCAL_MESH_NUM
    do f=1, this%nfaces_comm
      call this%commdata_list(f,n)%PC_Init(     &
        this%req_counter, this%request_pc(:) ) ! (inout)
    end do
    end do     

    return
  end subroutine MeshFieldCommBase_prepare_PC

!> Exchange halo data
!!
!OCL SERIAL
  subroutine MeshFieldCommBase_exchange_core( this, commdata_list, do_wait )
!    use scale_prof
    implicit none
  
    class(MeshFieldCommBase), intent(inout) :: this
    type(LocalMeshCommData), intent(inout), target :: commdata_list(this%nfaces_comm, this%mesh%LOCAL_MESH_NUM)
    logical, intent(in), optional :: do_wait
  
    integer :: n, f    
    integer :: ierr
    logical :: do_wait_
    !-----------------------------------------------------------------------------

    if ( present(do_wait) ) then
      do_wait_ = do_wait
    else
      do_wait_ = .true.
    end if    
    
!    call PROF_rapstart( 'meshfiled_comm_ex_core', 3)
    if ( this%MPI_pc_flag ) then
      call MPI_startall( this%req_counter, this%request_pc(1:this%req_counter), ierr )
    else
      this%req_counter = 0
    end if
    do n=1, this%mesh%LOCAL_MESH_NUM      
    do f=1, this%nfaces_comm
      call commdata_list(f,n)%SendRecv( &
        this%req_counter, this%request_send(:), this%request_recv(:), & ! (inout)
        commdata_list(:,:),                                           & ! (inout) 
        this%MPI_pc_flag )                                              ! (in)  
    end do
    end do
!    call PROF_rapend( 'meshfiled_comm_ex_core', 3)

    if ( do_wait_ ) then
      call MeshFieldCommBase_wait_core( this, commdata_list )
      this%call_wait_flag_sub_get = .false.
    else
      this%call_wait_flag_sub_get = .true.
    end if

    return
  end subroutine MeshFieldCommBase_exchange_core

!> Wait data communication
!!
!OCL SERIAL
  subroutine MeshFieldCommBase_wait_core( this, commdata_list )
    use mpi, only: &
      !MPI_waitall,    &
      MPI_STATUS_SIZE
    use scale_prof
    implicit none
  
    class(MeshFieldCommBase), intent(inout) :: this
    type(LocalMeshCommData), intent(inout), target :: commdata_list(this%nfaces_comm, this%mesh%LOCAL_MESH_NUM)

    integer :: ierr
    integer :: stat_send(MPI_STATUS_SIZE, this%nfaces_comm * this%mesh%LOCAL_MESH_NUM)
    integer :: stat_recv(MPI_STATUS_SIZE, this%nfaces_comm * this%mesh%LOCAL_MESH_NUM)
    integer :: stat_pc(MPI_STATUS_SIZE, 2*this%nfaces_comm * this%mesh%LOCAL_MESH_NUM)

    integer :: n, f
    integer :: var_id
    integer :: irs, ire
    !----------------------------

!    call PROF_rapstart( 'meshfiled_comm_wait_core', 3)
    if ( this%MPI_pc_flag ) then
      if (this%req_counter > 0) then
        call MPI_waitall( this%req_counter, this%request_pc(1:this%req_counter), stat_pc(:,1:this%req_counter), ierr )
      end if
    else
      if ( this%req_counter > 0 ) then
        call MPI_waitall( this%req_counter, this%request_recv(1:this%req_counter), stat_recv(:,1:this%req_counter), ierr )
        call MPI_waitall( this%req_counter, this%request_send(1:this%req_counter), stat_send(:,1:this%req_counter), ierr )
      end if
    end if
!    call PROF_rapend( 'meshfiled_comm_wait_core', 3)
  
    !---------------------

!    call PROF_rapstart( 'meshfiled_comm_wait_post', 3)
    do n=1, this%mesh%LOCAL_MESH_NUM
      !$omp parallel do private(var_id,f,irs,ire)
      do var_id=1, this%field_num_tot
        irs = 1
        do f=1, this%nfaces_comm
          ire = irs + commdata_list(f,n)%Nnode_LCMeshFace - 1
          this%recv_buf(irs:ire,var_id,n) = commdata_list(f,n)%recv_buf(:,var_id)
          irs = ire + 1
        end do ! end loop for face
      end do
    end do
!    call PROF_rapend( 'meshfiled_comm_wait_post', 3)

    return
  end subroutine MeshFieldCommBase_wait_core

!OCL SERIAL
  subroutine MeshFieldCommBase_extract_bounddata(var, refElem, mesh, buf)
    implicit none
  
    class(ElementBase), intent(in) :: refElem
    class(LocalMeshBase), intent(in) :: mesh
    real(RP), intent(in) :: var(refElem%Np * mesh%NeA)
    real(RP), intent(out) :: buf(size(mesh%VmapB))

    integer :: i
    !-----------------------------------------------------------------------------
    !$omp parallel do
!OCL PREFETCH
    do i=1, size(buf)
      buf(i) = var(mesh%VmapB(i))
    end do
    return
  end subroutine MeshFieldCommBase_extract_bounddata

  subroutine MeshFieldCommBase_set_bounddata(buf, refElem, mesh, var)
    implicit none
    
    class(ElementBase), intent(in) :: refElem
    class(LocalMeshBase), intent(in) :: mesh
    real(RP), intent(in) :: buf(size(mesh%VmapB))
    real(RP), intent(inout) :: var(refElem%Np * mesh%NeA)
    !-----------------------------------------------------------------------------

    var(refElem%Np*mesh%NeE+1:refElem%Np*mesh%NeE+size(buf)) = buf(:)
    return
  end subroutine MeshFieldCommBase_set_bounddata  

  !-------------------------------------------
  
  subroutine LocalMeshCommData_Init( this, comm, lcmesh, faceID, Nnode_LCMeshFace )
    implicit none

    class(LocalMeshCommData), intent(inout) :: this
    class(MeshFieldCommBase), intent(in) :: comm
    class(LocalMeshBase), intent(in), target :: lcmesh
    integer, intent(in) :: faceID
    integer, intent(in) :: Nnode_LCMeshFace
    !-----------------------------------------------------------------------------

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

  subroutine LocalMeshCommData_SendRecv( this, &
      req_counter, req_send, req_recv,          &
      lccommdat_list, MPI_pc_flag )
    
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
    logical, intent(in) :: MPI_pc_flag
  
    integer :: tag
    integer :: bufsize
    integer :: ierr
    !-------------------------------------------

    if ( this%s_rank /= this%lcmesh%PRC_myrank &
      .and. (.not. MPI_pc_flag ) ) then

      req_counter = req_counter + 1

      tag = 10 * this%lcmesh%tileID + this%faceID
      bufsize = size(this%recv_buf)
      call MPI_irecv( this%recv_buf(1,1), bufsize, MPI_DOUBLE_PRECISION, &
        this%s_rank, tag, PRC_LOCAL_COMM_WORLD,                          &
        req_recv(req_counter), ierr)

      bufsize = size(this%send_buf)
      tag = 10 * this%s_tileID + abs(this%s_faceID)
      call MPI_isend( this%send_buf(1,1), bufsize, MPI_DOUBLE_PRECISION, &
       this%s_rank, tag, PRC_LOCAL_COMM_WORLD,                           &
       req_send(req_counter), ierr )
      
    else if ( this%s_rank == this%lcmesh%PRC_myrank ) then
      lccommdat_list(abs(this%s_faceID), this%s_tilelocalID)%recv_buf(:,:) &
        = this%send_buf(:,:)
    end if         
    
    return
  end subroutine LocalMeshCommData_sendrecv

  subroutine LocalMeshCommData_pc_init( this, &
      req_counter, req                        )
    
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD
    use mpi, only: &
      !MPI_isend, MPI_irecv, &
      MPI_DOUBLE_PRECISION
    implicit none

    class(LocalMeshCommData), intent(inout) :: this
    integer, intent(inout) :: req_counter
    integer, intent(inout) :: req(:)
  
    integer :: tag
    integer :: bufsize
    integer :: ierr
    !-------------------------------------------

    if ( this%s_rank /= this%lcmesh%PRC_myrank ) then

      req_counter = req_counter + 1

      tag = 100 * this%lcmesh%tileID + this%faceID
      bufsize = size(this%recv_buf)
      call MPI_recv_init( this%recv_buf(1,1), bufsize, MPI_DOUBLE_PRECISION, &
        this%s_rank, tag, PRC_LOCAL_COMM_WORLD,                              &
        req(req_counter), ierr)

      !--
      req_counter = req_counter + 1

      bufsize = size(this%send_buf)
      tag = 100 * this%s_tileID + abs(this%s_faceID)
      call MPI_send_init( this%send_buf(1,1), bufsize, MPI_DOUBLE_PRECISION, &
       this%s_rank, tag, PRC_LOCAL_COMM_WORLD,                               &
       req(req_counter), ierr )
    end if         
    
    return
  end subroutine LocalMeshCommData_pc_init

  subroutine LocalMeshCommData_Final( this )
    implicit none
    
    class(LocalMeshCommData), intent(inout) :: this
    !-----------------------------------------------------------------------------

    deallocate( this%send_buf )
    deallocate( this%recv_buf )

    return
  end subroutine LocalMeshCommData_Final

end module scale_meshfieldcomm_base