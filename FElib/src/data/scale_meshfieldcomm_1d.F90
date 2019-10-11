#include "scaleFElib.h"
module scale_meshfieldcomm_1d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: elementbase, elementbase1d
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


  type, public, extends(MeshFieldCommBase) :: MeshFieldComm1D
    class(MeshBase1D), pointer :: mesh1d
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
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
  subroutine MeshFieldComm1D_Init( this, &
    sfield_num, hvfield_num, mesh1d )

    implicit none
    
    class(MeshFieldComm1D), intent(inout) :: this
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    class(Meshbase1d), intent(in), target :: mesh1d    
    !-----------------------------------------------------------------------------
    
    this%mesh1d => mesh1d 
    call MeshFieldCommBase_Init( this, sfield_num, hvfield_num, mesh1d%refElem1D%Nfp * 2, 2, mesh1d)  
  
    return
  end subroutine MeshFieldComm1D_Init

  subroutine MeshFieldComm1D_Final( this )
    implicit none
    
    class(MeshFieldComm1D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call MeshFieldCommBase_Final( this )

    return
  end subroutine MeshFieldComm1D_Final

  subroutine MeshFieldComm1D_put(this, field_list, varid_s)
    implicit none
    class(MeshFieldComm1D), intent(inout) :: this
    type(MeshFieldContainer), intent(in) :: field_list(:)
    integer, intent(in) :: varid_s
  
    integer :: i
    integer :: n
    type(LocalMesh1D), pointer :: mesh
    !-----------------------------------------------------------------------------
    
    do i=1, size(field_list)
    do n=1, this%mesh%LOCAL_MESH_NUM
      mesh => this%mesh1d%lcmesh_list(n)
      call MeshFieldCommBase_extract_bounddata( field_list(i)%field1d%local(n)%val, mesh%refElem, mesh, & ! (in)
        this%send_buf(:,varid_s+i-1,n) )                                                                  ! (out)
    end do
    end do

    return
  end subroutine MeshFieldComm1D_put

  subroutine MeshFieldComm1D_get(this, field_list, varid_s)
    implicit none
    
    class(MeshFieldComm1D), intent(in) :: this
    type(MeshFieldContainer), intent(inout) :: field_list(:)
    integer, intent(in) :: varid_s

    integer :: i
    integer :: n
    type(LocalMesh1D), pointer :: mesh
    !-----------------------------------------------------------------------------
    
    do i=1, size(field_list) 
    do n=1, this%mesh1D%LOCAL_MESH_NUM
      mesh => this%mesh1d%lcmesh_list(n)
      call MeshFieldCommBase_set_bounddata( this%recv_buf(:,varid_s+i-1,n), mesh%refElem, mesh, & !(in)
        field_list(i)%field1d%local(n)%val )                                                      !(out)
    end do
    end do

    return
  end subroutine MeshFieldComm1D_get

  subroutine MeshFieldComm1D_exchange( this )

    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD, PRC_abort, PRC_MPIbarrier
    
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_exchange_core,  &
      LocalMeshCommData

    implicit none
  
    class(MeshFieldComm1D), intent(inout) :: this
  
    integer :: n, f
    type(LocalMesh1D), pointer :: mesh
    integer, parameter :: Nnode_LCMeshFace = 1
    type(LocalMeshCommData), target :: commdata_list(this%nfaces_comm, this%mesh%LOCAL_MESH_NUM)
    type(LocalMeshCommData), pointer :: commdata
    !-----------------------------------------------------------------------------
    
    do n=1, this%mesh%LOCAL_MESH_NUM
      mesh => this%mesh1d%lcmesh_list(n)
      do f=1, this%nfaces_comm
        
        commdata => commdata_list(f,n)
        call commdata%Init(this, mesh, f, Nnode_LCMeshFace)

        call push_localsendbuf( commdata%send_buf(:,:),      &  ! (inout)
          this%send_buf(:,:,n), commdata%s_faceID, f,        &  ! (in)
          commdata%Nnode_LCMeshFace, this%field_num_tot)        ! (in)
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
  end subroutine MeshFieldComm1D_exchange

!----------------------------

  subroutine push_localsendbuf( lc_send_buf, send_buf, s_faceID, f, Nnode_LCMeshFace, var_num )
    implicit none

    integer, intent(in) :: var_num
    integer, intent(in) ::  Nnode_LCMeshFace
    real(RP), intent(inout) :: lc_send_buf(Nnode_LCMeshFace,var_num)
    real(RP), intent(in) :: send_buf(Nnode_LCMeshFace*2,var_num)  
    integer, intent(in) :: s_faceID, f
  
    integer :: is, ie, lincrement
    !-----------------------------------------------------------------------------
  
    if ( s_faceID > 0 ) then
      is = 1 + (f-1)*Nnode_LCMeshFace
      ie = is + Nnode_LCMeshFace - 1
      lincrement = +1          
    else
      is   = f*Nnode_LCMeshFace
      ie   = 1 + (f - 1)*Nnode_LCMeshFace          
      lincrement = -1          
    end if 
    lc_send_buf(:,:) = send_buf(is:ie:lincrement,:) 
    
    return
  end subroutine push_localsendbuf

end module scale_meshfieldcomm_1d