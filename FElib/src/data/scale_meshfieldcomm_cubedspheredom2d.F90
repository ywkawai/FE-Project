!-------------------------------------------------------------------------------
!> module FElib / Data / Communication in 2D cubed-sphere domain
!!
!! @par Description
!!      A module to mangage data communication with 2D cubed-sphere domain for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
#include "scaleFElib.h"
module scale_meshfieldcomm_cubedspheredom2d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: elementbase, elementBase2D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
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

  type :: VecCovariantComp
    type(MeshField2D), pointer :: u1 => null()
    type(MeshField2D), pointer :: u2 => null()
  end type    

  type, public, extends(MeshFieldCommBase) :: MeshFieldCommCubedSphereDom2D
    class(MeshCubedSphereDom2D), pointer :: mesh2d
    type(VecCovariantComp), allocatable :: vec_covariant_comp_ptrlist(:)
    integer, allocatable :: Nnode_LCMeshAllFace(:)
  contains
    procedure, public :: Init   => MeshFieldCommCubedSphereDom2D_Init
    procedure, public :: Put    => MeshFieldCommCubedSphereDom2D_put
    procedure, public :: Get    => MeshFieldCommCubedSphereDom2D_get
    procedure, public :: Exchange => MeshFieldCommCubedSphereDom2D_exchange  
    procedure, public :: SetCovariantVec => MeshFieldCommCubedSphereDom2D_set_covariantvec
    procedure, public :: Final => MeshFieldCommCubedSphereDom2D_Final
  end type MeshFieldCommCubedSphereDom2D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer :: bufsize_per_field
  integer, parameter :: COMM_FACE_NUM = 4

contains
  subroutine MeshFieldCommCubedSphereDom2D_Init( this, &
    sfield_num, hvfield_num, htensorfield_num, mesh2d )

    implicit none
    
    class(MeshFieldCommCubedSphereDom2D), intent(inout) :: this
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    integer, intent(in) :: htensorfield_num
    class(MeshCubedSphereDom2D), intent(in), target :: mesh2d
    
    type(LocalMesh2D), pointer :: lcmesh
    integer :: n
    integer :: Nnode_LCMeshFace(COMM_FACE_NUM,mesh2d%LOCAL_MESH_NUM)
    !-----------------------------------------------------------------------------
    
    this%mesh2d => mesh2d
    lcmesh => mesh2d%lcmesh_list(1)
    bufsize_per_field = 2*(lcmesh%NeX + lcmesh%NeY)*lcmesh%refElem2D%Nfp

    allocate( this%Nnode_LCMeshAllFace(mesh2d%LOCAL_MESH_NUM) )
    do n=1, this%mesh2d%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      Nnode_LCMeshFace(:,n) = (/ lcmesh%NeX, lcmesh%NeY, lcmesh%NeX, lcmesh%NeY /) * lcmesh%refElem2D%Nfp
      this%Nnode_LCMeshAllFace(n) = sum(Nnode_LCMeshFace(:,n))
    end do

    call MeshFieldCommBase_Init( this, sfield_num, hvfield_num, htensorfield_num, bufsize_per_field, COMM_FACE_NUM, Nnode_LCMeshFace, mesh2d)  

    if (hvfield_num > 0) then
      allocate( this%vec_covariant_comp_ptrlist(hvfield_num) )
    end if

    return
  end subroutine MeshFieldCommCubedSphereDom2D_Init

  subroutine MeshFieldCommCubedSphereDom2D_Final( this )
    implicit none
    
    class(MeshFieldCommCubedSphereDom2D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    deallocate( this%Nnode_LCMeshAllFace )

    if ( this%hvfield_num > 0 ) then
      deallocate( this%vec_covariant_comp_ptrlist )
    end if

    call MeshFieldCommBase_Final( this )

    return
  end subroutine MeshFieldCommCubedSphereDom2D_Final

  subroutine MeshFieldCommCubedSphereDom2D_set_covariantvec( &
    this, hvfield_ID, u1, u2  )
    implicit none
    class(MeshFieldCommCubedSphereDom2D), intent(inout) :: this
    integer, intent(in) :: hvfield_ID
    type(MeshField2D), intent(in), target :: u1
    type(MeshField2D), intent(in), target :: u2
    !--------------------------------------------------------------

    this%vec_covariant_comp_ptrlist(hvfield_ID)%u1 => u1
    this%vec_covariant_comp_ptrlist(hvfield_ID)%u2 => u2

    return
  end subroutine MeshFieldCommCubedSphereDom2D_set_covariantvec

  subroutine MeshFieldCommCubedSphereDom2D_put(this, field_list, varid_s)
    implicit none
    class(MeshFieldCommCubedSphereDom2D), intent(inout) :: this
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
        this%send_buf(:,varid_s+i-1,n) )                                                                      ! (out)
    end do
    end do

    return
  end subroutine MeshFieldCommCubedSphereDom2D_put

  subroutine MeshFieldCommCubedSphereDom2D_get(this, field_list, varid_s)
    implicit none
    
    class(MeshFieldCommCubedSphereDom2D), intent(in) :: this
    type(MeshFieldContainer), intent(inout) :: field_list(:)
    integer, intent(in) :: varid_s

    integer :: i
    integer :: n
    type(Localmesh2d), pointer :: lcmesh

    integer :: varnum
    integer :: varid_e
    integer :: varid_vec_s
    type(MeshField2D), pointer :: u1, u2
    !-----------------------------------------------------------------------------

    varnum = size(field_list) 

    do i=1, varnum
    do n=1, this%mesh2d%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      call MeshFieldCommBase_set_bounddata( this%recv_buf(:,varid_s+i-1,n), lcmesh%refElem, lcmesh, & !(in)
         field_list(i)%field2d%local(n)%val )  !(out)
    end do
    end do

    varid_e = varid_s + varnum - 1
    if ( varid_e > this%sfield_num ) then
      do i=1, this%hvfield_num

        varid_vec_s = this%sfield_num + 2*i - 1 
        if ( varid_vec_s > varid_e ) exit

        if (       associated(this%vec_covariant_comp_ptrlist(i)%u1 ) &
            .and.  associated(this%vec_covariant_comp_ptrlist(i)%u2 ) ) then
        
          do n=1, this%mesh2d%LOCAL_MESH_NUM
            call set_boundary_data2D_u1u2( &
              this%recv_buf(:,varid_vec_s,n), this%recv_buf(:,varid_vec_s+1,n), & ! (in)
              lcmesh%refElem2D, lcmesh, lcmesh%G_ij,                            & ! (in)
              this%vec_covariant_comp_ptrlist(i)%u1%local(n)%val,               & ! (out)
              this%vec_covariant_comp_ptrlist(i)%u2%local(n)%val                ) ! (out)
          end do
       end if
      end do
    end if

    return
  end subroutine MeshFieldCommCubedSphereDom2D_get

!OCL SERIAL
  subroutine MeshFieldCommCubedSphereDom2D_exchange( this )    
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_exchange_core,  &
      LocalMeshCommData

    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatVec, &
      CubedSphereCoordCnv_LonLat2CSVec
    
    implicit none
  
    class(MeshFieldCommCubedSphereDom2D), intent(inout), target :: this
  
    integer :: n, f
    integer :: varid

    real(RP), allocatable :: fpos2D(:,:)
    real(RP), allocatable :: lcfpos2D(:,:)
    real(RP), allocatable :: unity_fac(:)
    real(RP), allocatable :: tmp_svec2D(:,:)
    
    class(ElementBase2D), pointer :: elem
    type(LocalMesh2D), pointer :: lcmesh
    type(LocalMeshCommData), pointer :: commdata

    integer :: irs, ire    
    !-----------------------------------------------------------------------------

    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      elem => lcmesh%refElem2D

      allocate( fpos2D(this%Nnode_LCMeshAllFace(n),2) )
      call extract_boundary_data2D( lcmesh%pos_en(:,:,1), elem, lcmesh, fpos2D(:,1) )  
      call extract_boundary_data2D( lcmesh%pos_en(:,:,2), elem, lcmesh, fpos2D(:,2) ) 

      irs = 1
      do f=1, this%nfaces_comm
        commdata => this%commdata_list(f,n)
        call push_localsendbuf( commdata%send_buf(:,:),             &  ! (inout)
          this%send_buf(:,:,n), commdata%s_faceID, this%is_f(f,n),  &  ! (in)
          commdata%Nnode_LCMeshFace, this%field_num_tot)               ! (in)
        
        if ( commdata%s_panelID /= lcmesh%panelID ) then
          if ( this%hvfield_num > 0 ) then
            
            allocate( lcfpos2D(commdata%Nnode_LCMeshFace,2), unity_fac(commdata%Nnode_LCMeshFace) )
            allocate( tmp_svec2D(commdata%Nnode_LCMeshFace,2) )

            call push_localsendbuf( lcfpos2D,                                   &
              fpos2D, commdata%s_faceID, this%is_f(f,n), commdata%Nnode_LCMeshFace, 2 )
            unity_fac(:) = 1.0_RP

            ire = irs + commdata%Nnode_LCMeshFace - 1

            do varid=this%sfield_num+1, this%field_num_tot-1,2
              tmp_svec2D(:,1) = commdata%send_buf(:,varid  )
              tmp_svec2D(:,2) = commdata%send_buf(:,varid+1)
              call CubedSphereCoordCnv_CS2LonLatVec( &
                lcmesh%panelID, lcfpos2D(:,1), lcfpos2D(:,2), unity_fac,           &
                commdata%Nnode_LCMeshFace,                                         &
                tmp_svec2D(:,1), tmp_svec2D(:,2),                                  &
                commdata%send_buf(:,varid), commdata%send_buf(:,varid+1)           )
            end do

            deallocate( lcfpos2D, unity_fac, tmp_svec2D )
          end if
        end if
        
        irs = ire + 1
      end do
      deallocate( fpos2D )
    end do

    !-----------------------

    call MeshFieldCommBase_exchange_core(this, this%commdata_list(:,:))

    !-----------------------
  
    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      elem => lcmesh%refElem2D

      allocate( fpos2D(this%Nnode_LCMeshAllFace(n),2) )
      call extract_boundary_data2D( lcmesh%pos_en(:,:,1), elem, lcmesh, fpos2D(:,1) )  
      call extract_boundary_data2D( lcmesh%pos_en(:,:,2), elem, lcmesh, fpos2D(:,2) ) 

      irs = 1
      do f=1, this%nfaces_comm
        commdata => this%commdata_list(f,n)
        ire = irs + commdata%Nnode_LCMeshFace - 1

        if ( commdata%s_panelID /= lcmesh%panelID ) then
          if ( this%hvfield_num > 0 ) then

            allocate( lcfpos2D(commdata%Nnode_LCMeshFace,2), unity_fac(commdata%Nnode_LCMeshFace) )

            call push_localsendbuf( lcfpos2D,                          &
              fpos2D, f, this%is_f(f,n), commdata%Nnode_LCMeshFace, 2 )
            unity_fac(:) = 1.0_RP

            do varid=this%sfield_num+1, this%field_num_tot-1, 2
              call CubedSphereCoordCnv_LonLat2CSVec( &
                lcmesh%panelID, lcfpos2D(:,1), lcfpos2D(:,2), unity_fac(:),          &
                commdata%Nnode_LCMeshFace,                                           &
                commdata%recv_buf(:,varid), commdata%recv_buf(:,varid+1),            &
                this%recv_buf(irs:ire,varid,n), this%recv_buf(irs:ire,varid+1,n)     )
            end do
            deallocate( lcfpos2D, unity_fac )
          end if
        end if

        irs = ire + 1
      end do

      deallocate( fpos2D )
    end do 

    return
  end subroutine MeshFieldCommCubedSphereDom2D_exchange

!----------------------------

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

  subroutine extract_boundary_data2D( var, refElem, mesh, buf )
    implicit none
  
    type(elementbase2D), intent(in) :: refElem
    type(LocalMesh2D), intent(in) :: mesh
    real(DP), intent(in) :: var(refElem%Np * mesh%Ne)
    real(DP), intent(inout) :: buf(refElem%Nfp * mesh%NeX * 4)
    !------------------------------------------------------------
  
    buf(:) = var(mesh%VmapB(:))

    return
  end subroutine extract_boundary_data2D

  subroutine set_boundary_data2D_u1u2( buf_U, buf_V, &
    elem, mesh, G_ij,                             &
    u1, u2)
  
    implicit none
  
    type(ElementBase2D), intent(in) :: elem
    type(LocalMesh2D), intent(in) :: mesh
    real(DP), intent(in) :: buf_U(elem%Nfp * mesh%NeX * 4)
    real(DP), intent(in) :: buf_V(elem%Nfp * mesh%NeX * 4)
    real(DP), intent(in) :: G_ij(elem%Np * mesh%Ne,2,2)
    real(DP), intent(inout) :: u1(elem%Np * mesh%NeA)
    real(DP), intent(inout) :: u2(elem%Np * mesh%NeA)
    !------------------------------------------------------------

    u1(elem%Np*mesh%NeE+1:elem%Np*mesh%NeE+size(buf_U)) &
      = G_ij(mesh%VmapB,1,1) * buf_U(:) + G_ij(mesh%VmapB,1,2) * buf_V(:)
    u2(elem%Np*mesh%NeE+1:elem%Np*mesh%NeE+size(buf_U)) &
      = G_ij(mesh%VmapB,2,1) * buf_U(:) + G_ij(mesh%VmapB,2,2) * buf_V(:)
  
    return
  end subroutine set_boundary_data2D_u1u2

end module scale_meshfieldcomm_cubedspheredom2d