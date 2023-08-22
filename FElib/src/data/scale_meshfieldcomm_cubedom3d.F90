#include "scaleFElib.h"
module scale_meshfieldcomm_cubedom3d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: elementbase, elementBase3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_meshfield_base, only: MeshField3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase,                               &
    MeshFieldCommBase_Init, MeshFieldCommBase_Final, &
    MeshFieldCommBase_extract_bounddata,             &
    MeshFieldCommBase_set_bounddata,                 &
    MeshFieldContainer
  use scale_localmesh_3d, only: Localmesh3d
   
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 


  type, public, extends(MeshFieldCommBase) :: MeshFieldCommCubeDom3D
    class(MeshCubeDom3D), pointer :: mesh3d
  contains
    procedure, public :: Init => MeshFieldCommCubeDom3D_Init
    procedure, public :: Put => MeshFieldCommCubeDom3D_put
    procedure, public :: Get => MeshFieldCommCubeDom3D_get
    procedure, public :: Exchange => MeshFieldCommCubeDom3D_exchange  
    procedure, public :: Final => MeshFieldCommCubeDom3D_Final
  end type MeshFieldCommCubeDom3D

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
  integer :: bufsize_per_field

contains
  subroutine MeshFieldCommCubeDom3D_Init( this, &
    sfield_num, hvfield_num, htensorfield_num, mesh3d )

    implicit none
    
    class(MeshFieldCommCubeDom3D), intent(inout) :: this
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    integer, intent(in) :: htensorfield_num
    class(MeshCubeDom3D), intent(in), target :: mesh3d
    
    type(LocalMesh3D), pointer :: lcmesh
    type(ElementBase3D), pointer :: elem
    !-----------------------------------------------------------------------------
    
    this%mesh3d => mesh3d
    lcmesh => mesh3d%lcmesh_list(1)
    elem => lcmesh%refElem3D

    bufsize_per_field =  2*(lcmesh%NeX + lcmesh%NeY)*lcmesh%NeZ*elem%Nfp_h &
                       + 2*lcmesh%NeX*lcmesh%NeY*elem%Nfp_v
    call MeshFieldCommBase_Init( this, sfield_num, hvfield_num, htensorfield_num, bufsize_per_field, 6, mesh3d )  
  
    return
  end subroutine MeshFieldCommCubeDom3D_Init

  subroutine MeshFieldCommCubeDom3D_Final( this )

    implicit none    
    class(MeshFieldCommCubeDom3D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call MeshFieldCommBase_Final( this )

    return
  end subroutine MeshFieldCommCubeDom3D_Final

  subroutine MeshFieldCommCubeDom3D_put(this, field_list, varid_s)
    implicit none
    class(MeshFieldCommCubeDom3D), intent(inout) :: this
    type(MeshFieldContainer), intent(in) :: field_list(:)
    integer, intent(in) :: varid_s
  
    integer :: i
    integer :: n
    type(Localmesh3d), pointer :: lcmesh
    !-----------------------------------------------------------------------------

    do i=1, size(field_list)
    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      call MeshFieldCommBase_extract_bounddata( field_list(i)%field3d%local(n)%val, lcmesh%refElem, lcmesh, & ! (in)
        this%send_buf(:,varid_s+i-1,n) )                                                                      ! (out)
    end do
    end do

    return
  end subroutine MeshFieldCommCubeDom3D_put

  subroutine MeshFieldCommCubeDom3D_get(this, field_list, varid_s)
    implicit none
    
    class(MeshFieldCommCubeDom3D), intent(in) :: this
    type(MeshFieldContainer), intent(inout) :: field_list(:)
    integer, intent(in) :: varid_s

    integer :: i
    integer :: n
    type(Localmesh3d), pointer :: lcmesh
    !-----------------------------------------------------------------------------

    do i=1, size(field_list) 
    do n=1, this%mesh3d%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      call MeshFieldCommBase_set_bounddata( this%recv_buf(:,varid_s+i-1,n), lcmesh%refElem, lcmesh, & !(in)
         field_list(i)%field3d%local(n)%val )                                                         !(out)
    end do
    end do

    return
  end subroutine MeshFieldCommCubeDom3D_get


  subroutine MeshFieldCommCubeDom3D_exchange( this )

    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD, PRC_abort, PRC_MPIbarrier
    
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_exchange_core,  &
      LocalMeshCommData
    use scale_prof
    implicit none
  
    class(MeshFieldCommCubeDom3D), intent(inout) :: this
  
    integer :: n, f
    type(LocalMesh3D), pointer :: lcmesh
    integer :: Nnode_LCMeshFace(this%nfaces_comm)
    integer :: is_f(this%nfaces_comm)
    type(LocalMeshCommData), target :: commdata_list(this%nfaces_comm, this%mesh%LOCAL_MESH_NUM)
    type(LocalMeshCommData), pointer :: commdata
    !-----------------------------------------------------------------------------
    
    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      
      Nnode_LCMeshFace(:) = &
          (/ lcmesh%NeX, lcmesh%NeY, lcmesh%NeX, lcmesh%NeY, 0, 0 /) * lcmesh%NeZ * lcmesh%refElem3D%Nfp_h &
        + (/ 0, 0, 0, 0, 1, 1 /) * lcmesh%NeX*lcmesh%NeY * lcmesh%refElem3D%Nfp_v
      is_f(1) = 1
      do f=2, this%nfaces_comm
        is_f(f) = is_f(f-1) + Nnode_LCMeshFace(f-1)
      end do

      do f=1, this%nfaces_comm
        commdata => commdata_list(f,n)
        call commdata%Init(this, lcmesh, f, Nnode_LCMeshFace(f))

        call push_localsendbuf( commdata%send_buf(:,:),                                                    & ! (inout)
          this%send_buf(:,:,n), commdata%s_faceID, is_f(f), commdata%Nnode_LCMeshFace, this%field_num_tot, & ! (in)
          lcmesh )                                                                                           ! (in)
      end do
    end do

    !-----------------------

    call MeshFieldCommBase_exchange_core(this, commdata_list(:,:))

    !---------------------

    do n=1, this%mesh%LOCAL_MESH_NUM
    do f=1, this%nfaces_comm
      call commdata_list(f,n)%Final()
    end do
    end do 

    return
  end subroutine MeshFieldCommCubeDom3D_exchange

!----------------------------

  subroutine push_localsendbuf( lc_send_buf, send_buf, s_faceID, is, Nnode_LCMeshFace, var_num, lcmesh )
    implicit none

    integer, intent(in) ::  Nnode_LCMeshFace
    integer, intent(in) :: var_num
    type(LocalMesh3D), intent(in) :: lcmesh
    real(RP), intent(inout) :: lc_send_buf(Nnode_LCMeshFace,var_num)
    real(RP), intent(in) :: send_buf(bufsize_per_field,var_num)  
    integer, intent(in) :: s_faceID, is
  
    integer :: ie
    class(ElementBase3D), pointer :: e3D
    !-----------------------------------------------------------------------------

    e3D => lcmesh%refElem3D
    ie = is + Nnode_LCMeshFace - 1
    if ( s_faceID > 0 ) then
      lc_send_buf(:,:) = send_buf(is:ie,:)        
    else if ( -5 < s_faceID .and. s_faceID < 0) then
      call revert_hori( lc_send_buf(:,:), send_buf(is:ie,:), Nnode_LCMeshFace/lcmesh%NeZ, lcmesh%NeZ )      
    end if 

    return
  contains
    subroutine revert_hori(revert, ori, Ne_h1D, NeZ)
      integer, intent(in) :: Ne_h1D, NeZ
      real(RP), intent(out) :: revert(e3D%Nnode_h1D, e3D%Nnode_v, Ne_h1D, NeZ, var_num)
      real(RP), intent(in)  :: ori(e3D%Nnode_h1D, e3D%Nnode_v, Ne_h1D, NeZ, var_num)
      
      integer :: p1, p3, i, k, n
      integer :: i_, p1_
      !-----------------------------------------------------------------------------
      
      do n=1, var_num
      do k=1, NeZ
      do i=1, Ne_h1D
        i_ = Ne_h1D - i + 1
        do p3=1, e3D%Nnode_v
        do p1=1, e3D%Nnode_h1D
          p1_ = e3D%Nnode_h1D - p1 + 1
          revert(p1,p3,i,k,n) = ori(p1_,p3,i_,k,n)
        end do
        end do
      end do
      end do
      end do
    end subroutine revert_hori    
  end subroutine push_localsendbuf

end module scale_meshfieldcomm_cubedom3d