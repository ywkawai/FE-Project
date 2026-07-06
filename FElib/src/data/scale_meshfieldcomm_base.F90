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
  use scale_prc, only: &
    PRC_abort

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

  !> Derived type to manage data communication at a face between adjacent local meshes
  type, public :: LocalMeshCommData
    class(LocalMeshBase), pointer :: lcmesh
    real(RP), allocatable :: send_buf(:,:) !< Buffer for sending data
    real(RP), allocatable :: recv_buf(:,:) !< Buffer for receiving data
    integer :: Nnode_LCMeshFace !< Number of nodes at a face
    integer :: s_tileID       !< Destination tile ID
    integer :: s_faceID       !< Destination face ID
    integer :: s_panelID      !< Destination panel ID
    integer :: s_rank         !< Destination MPI rank
    integer :: s_tilelocalID  !< Destination local tile ID 
    integer :: faceID         !< Own face ID
  contains
    procedure, public :: Init => LocalMeshCommData_Init
    procedure, public :: Final => LocalMeshCommData_Final
    procedure, public :: SendRecv => LocalMeshCommData_SendRecv
    procedure, public :: PC_Init_recv => LocalMeshCommData_pc_init_recv
    procedure, public :: PC_Init_send => LocalMeshCommData_pc_init_send
  end type

  !> Base derived type to manage data communication
  type, abstract, public :: MeshFieldCommBase
    integer :: sfield_num        !< Number of scalar fields
    integer :: hvfield_num       !< Number of horizontal vector fields
    integer :: htensorfield_num  !< Number of horizontal tensor fields
    integer :: field_num_tot     !< Total number of fields
    integer :: nfaces_comm       !< Number of faces where halo data is communicated

    integer :: bufsize_per_field !< Buffer size per a field data

    class(MeshBase), pointer :: mesh
    real(RP), allocatable :: send_buf(:,:,:) !< Buffer for sending data
    real(RP), allocatable :: recv_buf(:,:,:) !< Buffer for receiving data
    integer, allocatable :: request_send(:)
    integer, allocatable :: request_recv(:)

    type(LocalMeshCommData), allocatable :: commdata_list(:,:)
    integer, allocatable :: is_f(:,:)

    logical :: MPI_pc_flag                !< Flag whether persistent communication is used
    logical :: use_mpi_pc_fujitsu_ext     !< Flag whether Fujitsu extension routines are used for persistent communication
    integer, allocatable :: request_pc(:)

    integer :: req_counter
    logical :: call_wait_flag_sub_get     !< Flag whether MPI_wait need to be called before getting halo data from recv_buf

    integer :: obj_ind

    !-
    logical :: use_vmap_wide_flag
    integer, allocatable :: VMapB_size(:)
    integer, allocatable :: VMapB2(:)
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
  public :: MeshFieldCommBase_extract_bounddata_2
  public :: MeshFieldCommBase_extract_bounddata_3
  public :: MeshFieldCommBase_set_bounddata
  public :: MeshFieldCommBase_set_bounddata_3

  !> Container to save a pointer of MeshField(1D, 2D, 3D) object
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

  integer            :: obj_ind       = 0
  integer, parameter :: OBJ_INDEX_MAX = 8192

contains

!> Initialize a base object to manage data communication of fields
!!
!! @param sfield_num Number of scalar fields
!! @param hvfield_num Number of horizontal vector fields
!! @param htensorfield_num Number of horizontal tensor fields
!! @param bufsize_per_field Buffer size per a field
!! @param comm_face_num Number of faces on a local mesh which perform data communication
!! @param Nnode_LCMeshFace Array to store the number of nodes  
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

    !$acc enter data create(this)
    this%mesh => mesh
    this%sfield_num = sfield_num
    this%hvfield_num = hvfield_num
    this%htensorfield_num = htensorfield_num
    this%field_num_tot = sfield_num + hvfield_num*2 + htensorfield_num*4
    this%nfaces_comm = comm_face_num
    !$acc update device(this%sfield_num, this%hvfield_num, this%htensorfield_num, this%field_num_tot, this%nfaces_comm)
    !$acc enter data attach(this%mesh)

    if (this%field_num_tot > 0) then
      allocate( this%send_buf(bufsize_per_field, this%field_num_tot, mesh%LOCAL_MESH_NUM) )
      allocate( this%recv_buf(bufsize_per_field, this%field_num_tot, mesh%LOCAL_MESH_NUM) )
      allocate( this%request_send(comm_face_num*mesh%LOCAL_MESH_NUM) )
      allocate( this%request_recv(comm_face_num*mesh%LOCAL_MESH_NUM) )
      !$acc enter data create(this%send_buf, this%recv_buf, this%request_send, this%request_recv)

      allocate( this%commdata_list(comm_face_num,mesh%LOCAL_MESH_NUM) )
      allocate( this%is_f(comm_face_num,mesh%LOCAL_MESH_NUM) ) 
      !$acc enter data create(this%is_f)

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
      !$acc update device(this%is_f)
      !$acc enter data copyin(this%commdata_list)
#ifdef _OPENACC
      do n=1, mesh%LOCAL_MESH_NUM
      do f=1, this%nfaces_comm
        call LocalMeshCommData_SetupForGPU( this%commdata_list(f,n) )
      end do
      end do
#endif
    end if 

    this%MPI_pc_flag = .false.

    obj_ind = obj_ind + 1
    if ( obj_ind > OBJ_INDEX_MAX ) then
      LOG_ERROR("MeshFieldCommBase_Init",*) 'obj_ind > OBJ_INDEX_MAX. Check!'
      call PRC_abort
    end if
    this%obj_ind = obj_ind
    
    return
  end subroutine MeshFieldCommBase_Init

!> Finalize a base object to manage data communication of fields
!!
!OCL SERIAL
  subroutine MeshFieldCommBase_Final( this )
    implicit none
    
    class(MeshFieldCommBase), intent(inout) :: this

    integer :: n, f
    integer :: ireq
    integer :: ierr
    !-----------------------------------------------------------------------------

    if (this%field_num_tot > 0) then
      !$acc exit data delete(this%send_buf, this%recv_buf)
      deallocate( this%send_buf, this%recv_buf )
      !$acc exit data delete(this%request_send, this%request_recv)
      deallocate( this%request_send, this%request_recv )

      do n=1, this%mesh%LOCAL_MESH_NUM
      do f=1, this%nfaces_comm
        call this%commdata_list(f,n)%Final()
      end do
      end do     
      !$acc exit data delete(this%commdata_list, this%is_f)
      deallocate( this%commdata_list, this%is_f )

      if ( this%MPI_pc_flag ) then
        do ireq=1, this%req_counter
          call MPI_request_free( this%request_pc(ireq), ierr )
        end do           
        deallocate( this%request_pc ) 
      end if
    end if
    return
  end subroutine MeshFieldCommBase_Final

!> Prepare persistent communication
!!
!! @param use_mpi_pc_fujitsu_ext Flag whether the extension routines of Fujitsu MPI are used
!!
!OCL SERIAL
  subroutine MeshFieldCommBase_prepare_PC( this, &
    use_mpi_pc_fujitsu_ext )
    implicit none    
    class(MeshFieldCommBase), intent(inout) :: this
    logical, intent(in), optional :: use_mpi_pc_fujitsu_ext

    integer :: n, f
    !-----------------------------------------------------------------------------

    this%MPI_pc_flag = .true.

#ifdef __FUJITSU
    this%use_mpi_pc_fujitsu_ext = .true.
#else
    this%use_mpi_pc_fujitsu_ext = .false.
#endif
    if ( present(use_mpi_pc_fujitsu_ext) ) then
      this%use_mpi_pc_fujitsu_ext = use_mpi_pc_fujitsu_ext
#ifndef __FUJITSU
      if ( use_mpi_pc_fujitsu_ext ) then
        LOG_ERROR("MeshFieldCommBase_prepare_PC",*) 'use_mpi_pc_fujitsu_ext=.true., but Fujitsu MPI is unavailable. Check!'
        call PRC_abort
      end if
#endif
    end if

    if (this%field_num_tot > 0) then
      allocate( this%request_pc(2*this%nfaces_comm*this%mesh%LOCAL_MESH_NUM) )
    end if

    this%req_counter = 0
    do n=1, this%mesh%LOCAL_MESH_NUM
    do f=1, this%nfaces_comm
      call this%commdata_list(f,n)%PC_Init_recv( &
        this%req_counter, this%request_pc(:),     & ! (inout)
        this%obj_ind, this%use_mpi_pc_fujitsu_ext ) ! (in)
    end do
    end do    
    do n=1, this%mesh%LOCAL_MESH_NUM
    do f=1, this%nfaces_comm
      call this%commdata_list(f,n)%PC_Init_send( &
        this%req_counter, this%request_pc(:),     & ! (inout)
        this%obj_ind, this%use_mpi_pc_fujitsu_ext ) ! (in)
    end do
    end do

    return
  end subroutine MeshFieldCommBase_prepare_PC

!> Exchange halo data
!!
!! @param commdata_list Array of LocalMeshCommData objects which manage information and halo data 
!! @param do_wait Flag whether MPI_waitall is called and move tmp data of LocalMeshCommData object to a recv buffer
!OCL SERIAL
  subroutine MeshFieldCommBase_exchange_core( this, commdata_list, do_wait )
#ifdef __FUJITSU
     use mpi_ext, only: &
       FJMPI_prequest_startall
#endif
!    use mpi, only: &
!      MPI_startall
    use scale_prof
    implicit none
  
    class(MeshFieldCommBase), intent(inout) :: this
    type(LocalMeshCommData), intent(inout), target :: commdata_list(this%nfaces_comm, this%mesh%LOCAL_MESH_NUM)
    logical, intent(in), optional :: do_wait
  
    integer :: n, f    
    integer :: ierr
    logical :: do_wait_

    integer :: i
    type(LocalMeshCommData), pointer :: lcommdata
    !-----------------------------------------------------------------------------

    if ( present(do_wait) ) then
      do_wait_ = do_wait
    else
      do_wait_ = .true.
    end if    
    
!    call PROF_rapstart( 'meshfiled_comm_ex_core', 3)
    !
    if ( this%MPI_pc_flag ) then
      !$omp parallel
      !$omp master
      if ( this%use_mpi_pc_fujitsu_ext ) then
#ifdef __FUJITSU
        ! Use Fujitsu MPI extension for persistent communication
        call FJMPI_prequest_startall( this%req_counter, this%request_pc(1:this%req_counter), ierr )
#endif
      else
        call MPI_startall( this%req_counter, this%request_pc(1:this%req_counter), ierr )
      end if
      !$omp end master
      !$omp do collapse(2) private(lcommdata)
      do n=1, this%mesh%LOCAL_MESH_NUM      
      do f=1, this%nfaces_comm
        lcommdata => commdata_list(f,n)
        if ( lcommdata%s_rank == lcommdata%lcmesh%PRC_myrank ) then
          commdata_list(abs(lcommdata%s_faceID), lcommdata%s_tilelocalID)%recv_buf(:,:) &
            = lcommdata%send_buf(:,:)
        end if
      end do 
      end do
      !$omp end parallel
    else
      this%req_counter = 0
      do n=1, this%mesh%LOCAL_MESH_NUM      
      do f=1, this%nfaces_comm
        call commdata_list(f,n)%SendRecv( &
          this%req_counter, this%request_send, this%request_recv, & ! (inout)
          commdata_list                                           ) ! (inout) 
      end do
      end do
      !$acc wait(1)
    end if

!    call PROF_rapend( 'meshfiled_comm_ex_core', 3)

    if ( do_wait_ ) then
      call MeshFieldCommBase_wait_core( this, commdata_list )
      this%call_wait_flag_sub_get = .false.
    else
      this%call_wait_flag_sub_get = .true.
    end if

    return
  end subroutine MeshFieldCommBase_exchange_core

!> Wait data communication and move tmp data of LocalMeshCommData object to a recv buffer
!! 
!! @param commdata_list Array of LocalMeshCommData objects which manage information and halo data 
!OCL SERIAL
  subroutine MeshFieldCommBase_wait_core( this, commdata_list, &
    field_list, dim, varid_s, lcmesh_list )
    use mpi, only: &
!      MPI_waitall,    &
      MPI_STATUS_SIZE
    use scale_prof
    implicit none
  
    class(MeshFieldCommBase), intent(inout) :: this
    type(LocalMeshCommData), intent(inout), target :: commdata_list(this%nfaces_comm, this%mesh%LOCAL_MESH_NUM)
    type(MeshFieldContainer), intent(inout), optional :: field_list(:)
    integer, intent(in), optional :: dim
    integer, intent(in), optional :: varid_s
    class(LocalMeshBase), intent(in), optional, target :: lcmesh_list(:)

    integer :: ierr
    integer :: stat_send(MPI_STATUS_SIZE, this%req_counter)
    integer :: stat_recv(MPI_STATUS_SIZE, this%req_counter)
    integer :: stat_pc(MPI_STATUS_SIZE, this%req_counter)

    integer :: n, f
    integer :: i, var_id
    integer :: irs(this%nfaces_comm,this%mesh%LOCAL_MESH_NUM), ire(this%nfaces_comm,this%mesh%LOCAL_MESH_NUM)

    class(LocalMeshBase), pointer :: lcmesh
    integer :: val_size(this%mesh%LOCAL_MESH_NUM)
    !----------------------------

!   call PROF_rapstart( 'meshfiled_comm_wait_core', 2)
    if ( this%MPI_pc_flag ) then
      if (this%req_counter > 0) then
        call MPI_waitall( this%req_counter, this%request_pc(1:this%req_counter), stat_pc, ierr )
      end if
    else
      if ( this%req_counter > 0 ) then
        call MPI_waitall( this%req_counter, this%request_recv(1:this%req_counter), stat_recv, ierr )
        call MPI_waitall( this%req_counter, this%request_send(1:this%req_counter), stat_send, ierr )
      end if
    end if

#ifdef _OPENACC
    do n=1, this%mesh%LOCAL_MESH_NUM
    do f=1, this%nfaces_comm
        if ( commdata_list(f,n)%s_rank /= commdata_list(f,n)%lcmesh%PRC_myrank ) then
          !$acc update device( commdata_list(f,n)%recv_buf )
        end if
    end do
    end do
#endif

!   call PROF_rapend( 'meshfiled_comm_wait_core', 2)
  
    !---------------------

!   call PROF_rapstart( 'meshfiled_comm_wait_post', 2)

    if ( present(field_list) ) then
      do n=1, this%mesh%LOCAL_MESH_NUM
        lcmesh => lcmesh_list(n)
        val_size(n) = lcmesh%refElem%Np * lcmesh%NeA
        irs(1,n) = lcmesh%refElem%Np * lcmesh%Ne + 1
        do f=1, this%nfaces_comm
          ire(f,n) = irs(f,n) + commdata_list(f,n)%Nnode_LCMeshFace - 1
          if (f<this%nfaces_comm) irs(f+1,n) = ire(f,n) + 1
        end do
      end do
      
      !$omp parallel do private(var_id,n,i,f) collapse(3)
      !$acc parallel loop gang collapse(3) present(field_list, commdata_list) copyin(irs, ire, val_size)
      do n=1, this%mesh%LOCAL_MESH_NUM
      do i=1, size(field_list)
      do f=1, this%nfaces_comm
        var_id = varid_s + i - 1

        if (dim==1) then
          call set_bounddata( field_list(var_id)%field1d%local(n)%val, val_size(n), irs(f,n), ire(f,n), commdata_list(f,n)%recv_buf(:,var_id) )
        else if (dim==2) then
          call set_bounddata( field_list(var_id)%field2d%local(n)%val, val_size(n), irs(f,n), ire(f,n), commdata_list(f,n)%recv_buf(:,var_id) )
        else if (dim==3) then
          call set_bounddata( field_list(var_id)%field3d%local(n)%val, val_size(n), irs(f,n), ire(f,n), commdata_list(f,n)%recv_buf(:,var_id) )
        end if
      end do ! end loop for face
      end do
      end do
#ifdef _OPENACC
      do n=1, this%mesh%LOCAL_MESH_NUM
      do i=1, size(field_list)
      do f=1, this%nfaces_comm
        var_id = varid_s + i - 1
        if (dim==1) then
          !$acc update device( field_list(var_id)%field1d%local(n)%val(irs(f,n):ire(f,n)) )
        else if (dim==2) then
          !$acc update device( field_list(var_id)%field2d%local(n)%val(irs(f,n):ire(f,n)) )
        else if (dim==3) then
          !$acc update device( field_list(var_id)%field3d%local(n)%val(irs(f,n):ire(f,n)) )
        end if
      end do ! end loop for face
      end do
      end do
#endif

    else
      
      do n=1, this%mesh%LOCAL_MESH_NUM
        irs(1,n) = 1
        do f=1, this%nfaces_comm
          ire(f,n) = irs(f,n) + commdata_list(f,n)%Nnode_LCMeshFace - 1
          if (f<this%nfaces_comm) irs(f+1,n) = ire(f,n) + 1
        end do
      end do

      !$omp parallel do private(n,var_id,f) collapse(3)
      !$acc parallel loop gang collapse(3) present(this%recv_buf, commdata_list) copyin(irs, ire)
      do n=1, this%mesh%LOCAL_MESH_NUM
      do var_id=1, this%field_num_tot
      do f=1, this%nfaces_comm
#ifdef _OPENACC
        call set_bounddata( this%recv_buf(:,var_id,n), size(this%recv_buf(:,var_id,n)), irs(f,n), ire(f,n), commdata_list(f,n)%recv_buf(:,var_id) )
#else        
        this%recv_buf(irs(f,n):ire(f,n),var_id,n) = commdata_list(f,n)%recv_buf(:,var_id)
#endif
      end do ! end loop for face
      end do
      end do
    end if
!   call PROF_rapend( 'meshfiled_comm_wait_post', 2)
    return
  contains
!OCL SERIAL
    subroutine set_bounddata( var, IA, irs_, ire_, recv_buf )
      implicit none
      integer, intent(in) :: IA
      real(RP), intent(inout) :: var(IA)
      integer, intent(in) :: irs_, ire_
      real(RP), intent(in) :: recv_buf(ire_-irs_+1)

      integer :: ii
      !-----------------------------
      !$acc routine vector

      !$acc loop vector
      do ii=1, size(recv_buf)
        var(irs_+ii-1) = recv_buf(ii)
      end do
      return
    end subroutine set_bounddata
  end subroutine MeshFieldCommBase_wait_core

!> Extract halo data from data array with MeshField object and set it to the receiving buffer
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
    !$acc parallel loop present(var, mesh%VmapB, buf)
!OCL PREFETCH
    do i=1, size(buf)
      buf(i) = var(mesh%VmapB(i))
    end do
    return
  end subroutine MeshFieldCommBase_extract_bounddata

!> Extract halo data from data array with MeshField object and set it to the receiving buffer
!OCL SERIAL
  subroutine MeshFieldCommBase_extract_bounddata2(var, VMapB, VMapB_size, NpxNeA, buf)
    implicit none
    integer, intent(in) :: VMapB_size
    integer, intent(in) :: NpxNeA
    real(RP), intent(in) :: var(NpxNeA)
    integer, intent(in) :: VMapB(VMapB_size)
    real(RP), intent(out) :: buf(VMapB_size)

    integer :: i
    !-----------------------------------------------------------------------------
    !$omp parallel do
    !$acc parallel loop present(var, VmapB, buf) async(1)
!OCL PREFETCH
    do i=1, size(buf)
      buf(i) = var(VmapB(i))
    end do
    return
  end subroutine MeshFieldCommBase_extract_bounddata2  

!> Extract halo data from data array with MeshField object and set it to the receiving buffer
!!
!! Note: We assume that the size of VmapB is the same for all local meshes. 
!! For the future, this subroutine should be modified to handle the case with different sizes of VmapB among local meshes.
!!
!OCL SERIAL
  subroutine MeshFieldCommBase_extract_bounddata_2(field_list, dim, varid_s, lcmesh_list, vmapB_size, buf)
    implicit none
    type(MeshFieldContainer), intent(in), target :: field_list(:)
    integer, intent(in) :: varid_s
    integer, intent(in) :: dim
    class(LocalMeshBase), intent(in), target :: lcmesh_list(:)
    integer, intent(in) :: vmapB_size                                            !< The size of VmapB (we assume that the size of VmapB is the same for all local meshes)
    real(RP), intent(out) :: buf(vmapB_size,size(field_list),size(lcmesh_list))

    class(LocalMeshBase), pointer :: lcmesh
    integer :: varid
    integer :: i, n

    integer :: field_num
    !-----------------------------------------------------------------------------

    field_num = size(field_list)

    do n=1, size(lcmesh_list)
      lcmesh => lcmesh_list(n)
      i = 1
      do while( i <= field_num )
        varid = varid_s + i - 1
        if ( i+1 <= field_num ) then
          if (dim==1) then
            call extract_bounddata_var2( buf(:,varid,n), buf(:,varid+1,n), field_list(varid)%field1d%local(n)%val, field_list(varid+1)%field1d%local(n)%val, lcmesh, lcmesh%refElem )
          else if(dim==2) then
            call extract_bounddata_var2( buf(:,varid,n), buf(:,varid+1,n), field_list(varid)%field2d%local(n)%val, field_list(varid+1)%field2d%local(n)%val, lcmesh, lcmesh%refElem )
          else if(dim==3) then
            call extract_bounddata_var2( buf(:,varid,n), buf(:,varid+1,n), field_list(varid)%field3d%local(n)%val, field_list(varid+1)%field3d%local(n)%val, lcmesh, lcmesh%refElem )
          end if
          i = i + 2
        else
          if (dim==1) then
            call MeshFieldCommBase_extract_bounddata2( field_list(varid)%field1d%local(n)%val, lcmesh%VMapB, size(lcmesh%VMapB), lcmesh%refElem%Np*lcmesh%NeA,  buf(:,varid,n) )
          else if(dim==2) then
            call MeshFieldCommBase_extract_bounddata2( field_list(varid)%field2d%local(n)%val, lcmesh%VMapB, size(lcmesh%VMapB), lcmesh%refElem%Np*lcmesh%NeA,  buf(:,varid,n) )
          else if(dim==3) then
            call MeshFieldCommBase_extract_bounddata2( field_list(varid)%field3d%local(n)%val, lcmesh%VMapB, size(lcmesh%VMapB), lcmesh%refElem%Np*lcmesh%NeA,  buf(:,varid,n) )
          end if
          i = i + 1
        end if
      end do
      !!$acc update host( buf(:,varid_s:field_num,n) )
    end do
    !$acc wait(1)
    return
  contains
!OCL SERIAL
    subroutine extract_bounddata_var2( buf1_, buf2_, var1, var2, lmesh, elem )
      implicit none
      class(LocalMeshBase), intent(in) :: lmesh
      class(ElementBase), intent(in) :: elem
      real(RP), intent(out) :: buf1_(size(lmesh%VMapB))
      real(RP), intent(out) :: buf2_(size(lmesh%VMapB))
      real(RP), intent(inout) :: var1(elem%Np*lmesh%NeA)
      real(RP), intent(inout) :: var2(elem%Np*lmesh%NeA)

      integer :: ii, b
      !-----------------------------
      !$omp parallel do
      !$acc parallel loop present(buf1_, buf2_, var1, var2, lmesh%vmapB) async(1)
!OCL PREFETCH
      do ii=1, size(buf1_)
        b = lmesh%vmapB(ii)
        buf1_(ii) = var1(b)
        buf2_(ii) = var2(b)
      end do
      return
    end subroutine extract_bounddata_var2
  end subroutine MeshFieldCommBase_extract_bounddata_2

!> Extract halo data from data array with MeshField object and set it to the recieving buffer
!OCL SERIAL
  subroutine MeshFieldCommBase_extract_bounddata_3(field_list, dim, varid_s, lcmesh_list, VMapB2, VMapB2_size, buf)
    implicit none
    class(MeshFieldContainer), intent(in), target :: field_list(:)
    integer, intent(in) :: varid_s
    integer, intent(in) :: dim
    class(LocalMeshBase), intent(in), target :: lcmesh_list(:)
    integer, intent(in) :: VMapB2_size
    integer, intent(in) :: VMapB2(VMapB2_size)
    real(RP), intent(out) :: buf(VMapB2_size,size(field_list),size(lcmesh_list))

    class(LocalMeshBase), pointer :: lcmesh
    integer :: varid
    integer :: i, n

    integer :: field_num
    !-----------------------------------------------------------------------------

    field_num = size(field_list)

    do n=1, size(lcmesh_list)
      lcmesh => lcmesh_list(n)
      i = 1
      do while( i <= field_num )
        varid = varid_s + i - 1
        if ( i+1 <= field_num ) then
          if (dim==1) then
            call extract_bounddata_var2( buf(:,varid,n), buf(:,varid+1,n), &
              VMapB2, VMapB2_size, field_list(varid)%field1d%local(n)%val, field_list(varid+1)%field1d%local(n)%val, lcmesh%refElem%Np * lcmesh%NeA )
          else if(dim==2) then
            call extract_bounddata_var2( buf(:,varid,n), buf(:,varid+1,n), &
              VMapB2, VMapB2_size, field_list(varid)%field2d%local(n)%val, field_list(varid+1)%field2d%local(n)%val, lcmesh%refElem%Np * lcmesh%NeA )
          else if(dim==3) then
            call extract_bounddata_var2( buf(:,varid,n), buf(:,varid+1,n), &
              VMapB2, VMapB2_size, field_list(varid)%field3d%local(n)%val, field_list(varid+1)%field3d%local(n)%val, lcmesh%refElem%Np * lcmesh%NeA )
          end if
          i = i + 2
        else
          if (dim==1) then
            call MeshFieldCommBase_extract_bounddata2( field_list(varid)%field1d%local(n)%val, VMapB2, VMapB2_size, lcmesh%refElem%Np*lcmesh%NeA,  &
              buf(:,varid,n) )
          else if(dim==2) then
            call MeshFieldCommBase_extract_bounddata2( field_list(varid)%field2d%local(n)%val, VMapB2, VMapB2_size, lcmesh%refElem%Np*lcmesh%NeA,  &
              buf(:,varid,n) )
          else if(dim==3) then
            call MeshFieldCommBase_extract_bounddata2( field_list(varid)%field3d%local(n)%val, VMapB2, VMapB2_size, lcmesh%refElem%Np*lcmesh%NeA,  &
              buf(:,varid,n) )
          end if
          i = i + 1
        end if
      end do
    end do
    return
  contains
!OCL SERIAL
    subroutine extract_bounddata_var2( buf1_, buf2_, VMapB, VMapB_size, var1, var2, varSize )
      implicit none
      integer, intent(in) :: vmapB_size
      real(RP), intent(out) :: buf1_(vmapB_size)
      real(RP), intent(out) :: buf2_(vmapB_size)
      integer, intent(in) :: vmapB(vmapB_size)
      integer, intent(in) :: varSize
      real(RP), intent(inout) :: var1(varSize)
      real(RP), intent(inout) :: var2(varSize)

      integer :: ii
      !-----------------------------
      !$omp parallel do
!OCL PREFETCH
      do ii=1, vmapB_size
        buf1_(ii) = var1(vmapB(ii))
        buf2_(ii) = var2(vmapB(ii))
      end do
      return
    end subroutine extract_bounddata_var2
  end subroutine MeshFieldCommBase_extract_bounddata_3

!> Extract halo data from the receiving buffer and set it to data array with MeshField object
  subroutine MeshFieldCommBase_set_bounddata(buf, refElem, mesh, var)
    implicit none
    
    class(ElementBase), intent(in) :: refElem
    class(LocalMeshBase), intent(in) :: mesh
    real(RP), intent(in) :: buf(size(mesh%VmapB))
    real(RP), intent(inout) :: var(refElem%Np * mesh%NeA)

    integer :: ii, iis
    !-----------------------------------------------------------------------------

#ifdef _OPENACC
     iis = refElem%Np*mesh%NeE
     !$acc parallel loop present(buf, var) async(1)
     do ii=1, size(buf)
       var(iis+ii) = buf(ii)
     end do
#else
    var(refElem%Np*mesh%NeE+1:refElem%Np*mesh%NeE+size(buf)) = buf(:)
#endif
    return
  end subroutine MeshFieldCommBase_set_bounddata  

!> Extract halo data from the recieving buffer and set it to data array with MeshField object
  subroutine MeshFieldCommBase_set_bounddata_3(buf, refElem, mesh, vmapB, vmapB_size, var)
    implicit none
    
    class(ElementBase), intent(in) :: refElem
    class(LocalMeshBase), intent(in) :: mesh
    integer, intent(in) :: vmapB_size
    integer, intent(in) :: vmapB(vmapB_size)
    real(RP), intent(in) :: buf(vmapB_size)
    real(RP), intent(inout) :: var(refElem%Np * mesh%NeA)
    !-----------------------------------------------------------------------------

    var(refElem%Np*mesh%NeE+1:refElem%Np*mesh%NeE+size(buf)) = buf(:)
    return
  end subroutine MeshFieldCommBase_set_bounddata_3
  
  !-------------------------------------------
  
  !> Initialize an object to manage data communication of fields for a face on a local mesh
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

    !-
    this%s_tileID  = comm%mesh%tileID_globalMap(faceID, lcmesh%tileID)
    this%s_faceID  = comm%mesh%tileFaceID_globalMap(faceID, lcmesh%tileID)
    this%s_panelID = comm%mesh%tilePanelID_globalMap(faceID, lcmesh%tileID)      
    this%s_rank    = comm%mesh%PRCrank_globalMap(this%s_tileID)
    this%s_tilelocalID = comm%mesh%tileID_global2localMap(this%s_tileID)
    this%faceID    = faceID

    allocate( this%send_buf(Nnode_LCMeshFace, comm%field_num_tot))
    allocate( this%recv_buf(Nnode_LCMeshFace, comm%field_num_tot))
    return
  end subroutine LocalMeshCommData_Init

  subroutine LocalMeshCommData_SetupForGPU( this )
    implicit none
    type(LocalMeshCommData), intent(inout) :: this
    !-----------------------------------------------------------------------------
    !$acc enter data create(this%send_buf, this%recv_buf)
    !$acc enter data attach(this%lcmesh)
    return
  end subroutine LocalMeshCommData_SetupForGPU

  !> Send and receive halo data for a face on a local mesh
  subroutine LocalMeshCommData_SendRecv( this, &
    req_counter, req_send, req_recv,           &
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

      !$acc update host(this%send_buf)

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
#ifdef _OPENACC
      call set_recvbuf_from_sendbuf( lccommdat_list(abs(this%s_faceID), this%s_tilelocalID)%recv_buf, &
        this%send_buf, size(this%send_buf) )
#else      
      lccommdat_list(abs(this%s_faceID), this%s_tilelocalID)%recv_buf(:,:) &
        = this%send_buf(:,:)
#endif
    end if         
    
    return
  contains
    subroutine set_recvbuf_from_sendbuf( recv_buf, send_buf, buf_size )
      implicit none
      integer, intent(in) :: buf_size
      real(RP), intent(out) :: recv_buf(buf_size)
      real(RP), intent(in) :: send_buf(buf_size)

      integer :: i
      !-----------------------------
      !$acc parallel loop present(recv_buf, send_buf) async(1)
      do i=1, buf_size
        recv_buf(i) = send_buf(i)
      end do
      return
    end subroutine set_recvbuf_from_sendbuf
  end subroutine LocalMeshCommData_SendRecv

  !> Initialize persistent communication for sending halo data
  subroutine LocalMeshCommData_pc_init_send( this, &
    req_counter, req, obj_ind_ ,      &
    use_mpi_pc_fujisu_ext             )
    
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD
    use mpi, only: &
      MPI_DOUBLE_PRECISION!, &
!      MPI_send_init
#ifdef __FUJITSU
    use mpi_ext, only: &
      FJMPI_Prequest_send_init
#endif
    implicit none

    class(LocalMeshCommData), intent(inout) :: this
    integer, intent(inout) :: req_counter
    integer, intent(inout) :: req(:)
    integer, intent(in) :: obj_ind_
    logical, intent(in) :: use_mpi_pc_fujisu_ext

    integer :: tag
    integer :: bufsize
    integer :: ierr
    !-------------------------------------------

    if ( this%s_rank /= this%lcmesh%PRC_myrank ) then
      req_counter = req_counter + 1

      bufsize = size(this%send_buf)
!      tag = 1000 * this%s_tileID + 10 * obj_ind_ + abs(this%s_faceID)
      tag = obj_ind_

      if ( use_mpi_pc_fujisu_ext ) then
#ifdef __FUJITSU
        call FJMPI_Prequest_send_init( this%send_buf(1,1), bufsize, MPI_DOUBLE_PRECISION, &
          this%s_rank, tag, PRC_LOCAL_COMM_WORLD,                                         &
          req(req_counter), ierr )
#endif
      else
        call MPI_send_init( this%send_buf(1,1), bufsize, MPI_DOUBLE_PRECISION, &
          this%s_rank, tag, PRC_LOCAL_COMM_WORLD,                              &
          req(req_counter), ierr )
      end if
    end if         

    return
  end subroutine LocalMeshCommData_pc_init_send
  !> Initialize persistent communication for receiving halo data
  subroutine LocalMeshCommData_pc_init_recv( this, &
    req_counter, req, obj_ind_,       &
    use_mpi_pc_fujisu_ext             )
    
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD
    use mpi, only: &
      MPI_DOUBLE_PRECISION!, &
!      MPI_recv_init
#ifdef __FUJITSU
    use mpi_ext, only: &
      FJMPI_Prequest_recv_init
#endif
    implicit none

    class(LocalMeshCommData), intent(inout) :: this
    integer, intent(inout) :: req_counter
    integer, intent(inout) :: req(:)
    integer, intent(in) :: obj_ind_
    logical, intent(in) :: use_mpi_pc_fujisu_ext
  
    integer :: tag
    integer :: bufsize
    integer :: ierr
    !-------------------------------------------

    if ( this%s_rank /= this%lcmesh%PRC_myrank ) then
      req_counter = req_counter + 1

      bufsize = size(this%recv_buf)
!      tag = 1000 * this%lcmesh%tileID + 10 * obj_ind_ + this%faceID
      tag = obj_ind_

      if ( use_mpi_pc_fujisu_ext ) then
#ifdef __FUJITSU
        call FJMPI_Prequest_recv_init( this%recv_buf(1,1), bufsize, MPI_DOUBLE_PRECISION, &
          this%s_rank, tag, PRC_LOCAL_COMM_WORLD,                                         &
          req(req_counter), ierr )
#endif
      else 
        call MPI_recv_init( this%recv_buf(1,1), bufsize, MPI_DOUBLE_PRECISION, &
          this%s_rank, tag, PRC_LOCAL_COMM_WORLD,                              &
          req(req_counter), ierr )
      end if
    end if         

    return
  end subroutine LocalMeshCommData_pc_init_recv

  !> Finalize an object to manage data communication of fields for a face on a local mesh
  subroutine LocalMeshCommData_Final( this )
    implicit none
    
    class(LocalMeshCommData), intent(inout) :: this
    !-----------------------------------------------------------------------------

    !$acc exit data delete(this%send_buf, this%recv_buf)
    deallocate( this%send_buf )
    deallocate( this%recv_buf )
    return
  end subroutine LocalMeshCommData_Final

end module scale_meshfieldcomm_base