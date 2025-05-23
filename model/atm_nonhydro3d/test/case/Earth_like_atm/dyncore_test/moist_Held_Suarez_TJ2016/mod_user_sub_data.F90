#include "scaleFElib.h"
module mod_user_sub_data

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: ElementBase, ElementBase3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_meshfield_base, only: MeshField3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase,                               &
    MeshFieldCommBase_Init, MeshFieldCommBase_Final, &
    MeshFieldCommBase_extract_bounddata,             &
    MeshFieldCommBase_set_bounddata,                 &
    MeshFieldContainer

  use scale_localmesh_2d, only: Localmesh2D
  use scale_localmesh_3d, only: Localmesh3D

   
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  type :: VecCovariantComp
    type(MeshField3D), pointer :: u1 => null()
    type(MeshField3D), pointer :: u2 => null()
  end type    

  type, public, extends(MeshFieldCommBase) :: MeshFieldCommCubedSphereDom3D_2
    class(MeshCubedSphereDom3D), pointer :: mesh3d
    type(VecCovariantComp), allocatable :: vec_covariant_comp_ptrlist(:)
    integer, allocatable :: Nnode_LCMeshAllFace(:)
  contains
    procedure, public :: Init   => MeshFieldCommCubedSphereDom3D_Init
    procedure, public :: Put    => MeshFieldCommCubedSphereDom3D_put
    procedure, public :: Get    => MeshFieldCommCubedSphereDom3D_get
    procedure, public :: Exchange => MeshFieldCommCubedSphereDom3D_exchange  
    procedure, public :: SetCovariantVec => MeshFieldCommCubedSphereDom3D_set_covariantvec
    procedure, public :: Final => MeshFieldCommCubedSphereDom3D_Final
  end type MeshFieldCommCubedSphereDom3D_2

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
  integer, parameter :: COMM_FACE_NUM = 6

contains
  subroutine MeshFieldCommCubedSphereDom3D_Init( this, &
    sfield_num, hvfield_num, htensorfield_num, mesh3d )

    implicit none
    
    class(MeshFieldCommCubedSphereDom3D_2), intent(inout), target :: this
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    integer, intent(in) :: htensorfield_num
    class(MeshCubedSphereDom3D), intent(in), target :: mesh3d
    
    type(LocalMesh3D), pointer :: lcmesh
    type(ElementBase3D), pointer :: elem
    integer :: n
    integer :: Nnode_LCMeshFace(COMM_FACE_NUM,mesh3d%LOCAL_MESH_NUM)
    !-----------------------------------------------------------------------------
    
    this%mesh3d => mesh3d
    lcmesh => mesh3d%lcmesh_list(1)
    elem => lcmesh%refElem3D

    bufsize_per_field =  2*(lcmesh%NeX + lcmesh%NeY)*lcmesh%NeZ*elem%Nfp_h*elem%Nnode_h1D &
                       + 2*lcmesh%NeX*lcmesh%NeY*elem%Nfp_v

    allocate( this%Nnode_LCMeshAllFace(mesh3d%LOCAL_MESH_NUM) )
    do n=1, this%mesh3d%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      Nnode_LCMeshFace(:,n) = &
          (/ lcmesh%NeX, lcmesh%NeY, lcmesh%NeX, lcmesh%NeY, 0, 0 /) * lcmesh%NeZ * lcmesh%refElem3D%Nfp_h*elem%Nnode_h1D &
        + (/ 0, 0, 0, 0, 1, 1 /) * lcmesh%NeX*lcmesh%NeY * lcmesh%refElem3D%Nfp_v
      this%Nnode_LCMeshAllFace(n) = sum(Nnode_LCMeshFace(:,n))
    end do

    call MeshFieldCommBase_Init( this, sfield_num, hvfield_num, htensorfield_num, bufsize_per_field, COMM_FACE_NUM, Nnode_LCMeshFace, mesh3d )  
  
    if (hvfield_num > 0) then
      allocate( this%vec_covariant_comp_ptrlist(hvfield_num) )
    end if

    return
  end subroutine MeshFieldCommCubedSphereDom3D_Init

  subroutine MeshFieldCommCubedSphereDom3D_Final( this )

    implicit none
    
    class(MeshFieldCommCubedSphereDom3D_2), intent(inout) :: this
    !-----------------------------------------------------------------------------

    if ( this%hvfield_num > 0 ) then
      deallocate( this%vec_covariant_comp_ptrlist )
    end if

    call MeshFieldCommBase_Final( this )

    return
  end subroutine MeshFieldCommCubedSphereDom3D_Final

  subroutine MeshFieldCommCubedSphereDom3D_set_covariantvec( &
    this, hvfield_ID, u1, u2  )
    implicit none
    class(MeshFieldCommCubedSphereDom3D_2), intent(inout) :: this
    integer, intent(in) :: hvfield_ID
    type(MeshField3D), intent(in), target :: u1
    type(MeshField3D), intent(in), target :: u2
    !--------------------------------------------------------------

    this%vec_covariant_comp_ptrlist(hvfield_ID)%u1 => u1
    this%vec_covariant_comp_ptrlist(hvfield_ID)%u2 => u2

    return
  end subroutine MeshFieldCommCubedSphereDom3D_set_covariantvec

  subroutine MeshFieldCommCubedSphereDom3D_put(this, field_list, varid_s)
    implicit none
    class(MeshFieldCommCubedSphereDom3D_2), intent(inout) :: this
    type(MeshFieldContainer), intent(in) :: field_list(:)
    integer, intent(in) :: varid_s
  
    integer :: i
    integer :: n
    type(Localmesh3D), pointer :: lcmesh
    !-----------------------------------------------------------------------------

    do i=1, size(field_list)
    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      ! call extract_bounddata( field_list(i)%field3d%local(n)%val, lcmesh%refElem, lcmesh, & ! (in)
      !   this%send_buf(:,varid_s+i-1,n) )                                                    ! (out)
      call extract_bounddata3D( field_list(i)%field3d%local(n)%val, lcmesh%refElem3D, lcmesh, & ! (in)
        this%send_buf(:,varid_s+i-1,n)   )                                                      ! (out)
    end do
    end do

    return
  end subroutine MeshFieldCommCubedSphereDom3D_put

  subroutine MeshFieldCommCubedSphereDom3D_get(this, field_list, varid_s )
    implicit none
    
    class(MeshFieldCommCubedSphereDom3D_2), intent(inout) :: this
    type(MeshFieldContainer), intent(inout) :: field_list(:)
    integer, intent(in) :: varid_s

    integer :: i
    integer :: n
    integer :: ke
    integer :: ke2D
    type(Localmesh3D), pointer :: lcmesh
    type(ElementBase3D), pointer :: elem

    integer :: varnum
    integer :: varid_e
    integer :: varid_vec_s

    real(RP), allocatable :: G_ij(:,:,:,:)
    !-----------------------------------------------------------------------------

    varnum = size(field_list) 

    do i=1, varnum
    do n=1, this%mesh3d%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      call set_bounddata( this%recv_buf(:,varid_s+i-1,n), lcmesh%refElem, lcmesh, & !(in)
         field_list(i)%field3d%local(n)%val )  !(out)
    end do
    end do

    varid_e = varid_s + varnum - 1
    if ( varid_e > this%sfield_num ) then
      do i=1, this%hvfield_num

        varid_vec_s = this%sfield_num + 2*i - 1 
        if ( varid_vec_s > varid_e ) exit

        if (       associated(this%vec_covariant_comp_ptrlist(i)%u1 ) &
            .and.  associated(this%vec_covariant_comp_ptrlist(i)%u2 ) ) then
        
          do n=1, this%mesh3d%LOCAL_MESH_NUM
            lcmesh => this%mesh3d%lcmesh_list(n)
            elem => lcmesh%refElem3D

            allocate( G_ij(elem%Np,lcmesh%Ne,2,2) )
            !$omp parallel do private(ke2D)
            do ke=lcmesh%NeS, lcmesh%NeE
              ke2D = lcmesh%EMap3Dto2D(ke)
              G_ij(:,ke,1,1) = lcmesh%G_ij(elem%IndexH2Dto3D(:),ke2D,1,1)
              G_ij(:,ke,2,1) = lcmesh%G_ij(elem%IndexH2Dto3D(:),ke2D,2,1)
              G_ij(:,ke,1,2) = lcmesh%G_ij(elem%IndexH2Dto3D(:),ke2D,1,2)
              G_ij(:,ke,2,2) = lcmesh%G_ij(elem%IndexH2Dto3D(:),ke2D,2,2)
            end do

            call set_boundary_data3D_u1u2( &
              this%recv_buf(:,varid_vec_s,n), this%recv_buf(:,varid_vec_s+1,n), & ! (in)
              lcmesh%refElem3D, lcmesh, G_ij(:,:,:,:),                          & ! (in)
              this%vec_covariant_comp_ptrlist(i)%u1%local(n)%val,               & ! (out)
              this%vec_covariant_comp_ptrlist(i)%u2%local(n)%val                ) ! (out)
          end do
       end if
      end do
    end if

    return
  end subroutine MeshFieldCommCubedSphereDom3D_get

!OCL SERIAL
  subroutine MeshFieldCommCubedSphereDom3D_exchange( this, do_wait )
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_exchange_core,  &
      LocalMeshCommData

    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatVec, &
      CubedSphereCoordCnv_LonLat2CSVec
    
    implicit none
  
    class(MeshFieldCommCubedSphereDom3D_2), intent(inout), target :: this
    logical, intent(in), optional :: do_wait

    integer :: n, f
    integer :: varid

    real(RP), allocatable :: fpos3D(:,:)
    real(RP), allocatable :: lcfpos3D(:,:)
    real(RP), allocatable :: unity_fac(:)
    real(RP), allocatable :: tmp_svec3D(:,:)
    real(RP), allocatable :: tmp1_htensor3D(:,:,:)
    real(RP), allocatable :: tmp2_htensor3D(:,:,:)
    
    class(ElementBase3D), pointer :: elem
    type(LocalMesh3D), pointer :: lcmesh
    type(LocalMeshCommData), pointer :: commdata

    integer :: irs, ire    
    !-----------------------------------------------------------------------------

    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      elem => lcmesh%refElem3D

      ! allocate( fpos3D(this%Nnode_LCMeshAllFace(n),2) )
      ! call extract_boundary_data3D( lcmesh%pos_en(:,:,1), elem, lcmesh, fpos3D(:,1) )  
      ! call extract_boundary_data3D( lcmesh%pos_en(:,:,2), elem, lcmesh, fpos3D(:,2) ) 

      do f=1, this%nfaces_comm
        commdata => this%commdata_list(f,n)
        call push_localsendbuf( commdata%send_buf(:,:),            &  ! (inout)
          this%send_buf(:,:,n), commdata%s_faceID, this%is_f(f,n), &  ! (in)
          commdata%Nnode_LCMeshFace, this%field_num_tot,           &  ! (in)
          lcmesh, elem                                             )  ! (in)
        
        ! if ( commdata%s_panelID /= lcmesh%panelID ) then
        !   if ( this%hvfield_num > 0 ) then
            
        !     allocate( lcfpos3D(commdata%Nnode_LCMeshFace,2), unity_fac(commdata%Nnode_LCMeshFace) )
        !     allocate( tmp_svec3D(commdata%Nnode_LCMeshFace,2) )

        !     call push_localsendbuf( lcfpos3D,                                          &
        !       fpos3D, commdata%s_faceID, this%is_f(f,n), commdata%Nnode_LCMeshFace, 2, &
        !       lcmesh, elem )
        !     unity_fac(:) = 1.0_RP

        !     do varid=this%sfield_num+1, this%sfield_num+2*this%hvfield_num-1, 2
        !       tmp_svec3D(:,1) = commdata%send_buf(:,varid  )
        !       tmp_svec3D(:,2) = commdata%send_buf(:,varid+1)
  
        !       call CubedSphereCoordCnv_CS2LonLatVec( &
        !         lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),        &
        !         commdata%Nnode_LCMeshFace,                                               &
        !         tmp_svec3D(:,1), tmp_svec3D(:,2),                                  &
        !         commdata%send_buf(:,varid), commdata%send_buf(:,varid+1)           )

        !     end do
        !     deallocate( lcfpos3D, unity_fac, tmp_svec3D )
        !   end if

        !   if ( this%htensorfield_num > 0 ) then
        !     allocate( lcfpos3D(commdata%Nnode_LCMeshFace,2), unity_fac(commdata%Nnode_LCMeshFace) )
        !     allocate( tmp1_htensor3D(commdata%Nnode_LCMeshFace,2,2) )
        !     allocate( tmp2_htensor3D(commdata%Nnode_LCMeshFace,2,2) )

        !     call push_localsendbuf( lcfpos3D,                                          &
        !       fpos3D, commdata%s_faceID, this%is_f(f,n), commdata%Nnode_LCMeshFace, 2, &
        !       lcmesh, elem )
        !     unity_fac(:) = 1.0_RP
            
        !     do varid=this%sfield_num+2*this%hvfield_num+1, this%field_num_tot-3, 4
        !       tmp1_htensor3D(:,1,1) = commdata%send_buf(:,varid  )
        !       tmp1_htensor3D(:,2,1) = commdata%send_buf(:,varid+1)
        !       tmp1_htensor3D(:,1,2) = commdata%send_buf(:,varid+2)
        !       tmp1_htensor3D(:,2,2) = commdata%send_buf(:,varid+3)
  
        !       call CubedSphereCoordCnv_CS2LonLatVec( &
        !         lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),        &
        !         commdata%Nnode_LCMeshFace,                                         &
        !         tmp1_htensor3D(:,1,1), tmp1_htensor3D(:,2,1),                      &
        !         tmp2_htensor3D(:,1,1), tmp2_htensor3D(:,2,1)                       )
        !       call CubedSphereCoordCnv_CS2LonLatVec( &
        !         lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),        &
        !         commdata%Nnode_LCMeshFace,                                         &
        !         tmp1_htensor3D(:,1,2), tmp1_htensor3D(:,2,2),                      &
        !         tmp2_htensor3D(:,1,2), tmp2_htensor3D(:,2,2)                       )
        !       call CubedSphereCoordCnv_CS2LonLatVec( &
        !         lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),        &
        !         commdata%Nnode_LCMeshFace,                                         &
        !         tmp2_htensor3D(:,1,1), tmp2_htensor3D(:,1,2),                      &
        !         commdata%send_buf(:,varid), commdata%send_buf(:,varid+2)           )
        !       call CubedSphereCoordCnv_CS2LonLatVec( &
        !         lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),        &
        !         commdata%Nnode_LCMeshFace,                                         &
        !         tmp2_htensor3D(:,2,1), tmp2_htensor3D(:,2,2),                      &
        !         commdata%send_buf(:,varid+1), commdata%send_buf(:,varid+3)         )
        !     end do
        !     deallocate( lcfpos3D, unity_fac, tmp1_htensor3D, tmp2_htensor3D )
        !   end if

        ! end if
        
      end do
      ! deallocate( fpos3D )
    end do

    !-----------------------

    call MeshFieldCommBase_exchange_core( this, this%commdata_list(:,:), do_wait )

    !-----------------------
  
    ! do n=1, this%mesh%LOCAL_MESH_NUM
    !   lcmesh => this%mesh3d%lcmesh_list(n)
    !   elem => lcmesh%refElem3D

    !   allocate( fpos3D(this%Nnode_LCMeshAllFace(n),2) )
    !   call extract_boundary_data3D( lcmesh%pos_en(:,:,1), elem, lcmesh, fpos3D(:,1) )  
    !   call extract_boundary_data3D( lcmesh%pos_en(:,:,2), elem, lcmesh, fpos3D(:,2) ) 

    !   irs = 1
    !   do f=1, this%nfaces_comm
    !     commdata => this%commdata_list(f,n)
    !     ire = irs + commdata%Nnode_LCMeshFace - 1

    !     if ( commdata%s_panelID /= lcmesh%panelID ) then
    !       if ( this%hvfield_num > 0 ) then

    !         allocate( lcfpos3D(commdata%Nnode_LCMeshFace,2), unity_fac(commdata%Nnode_LCMeshFace) )
    !         call push_localsendbuf( lcfpos3D,                          &
    !           fpos3D, f, this%is_f(f,n), commdata%Nnode_LCMeshFace, 2, &
    !           lcmesh, elem )
    !         unity_fac(:) = 1.0_RP
            
    !         do varid=this%sfield_num+1, this%sfield_num+2*this%hvfield_num-1, 2
    !           call CubedSphereCoordCnv_LonLat2CSVec( &
    !             lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),          &
    !             commdata%Nnode_LCMeshFace,                                           &
    !             commdata%recv_buf(:,varid), commdata%recv_buf(:,varid+1),            &
    !             this%recv_buf(irs:ire,varid,n), this%recv_buf(irs:ire,varid+1,n)     )
    !         end do
    !         deallocate( lcfpos3D, unity_fac )
    !       end if

    !       if ( this%htensorfield_num > 0 ) then
    !         allocate( lcfpos3D(commdata%Nnode_LCMeshFace,2), unity_fac(commdata%Nnode_LCMeshFace) )
    !         allocate( tmp1_htensor3D(commdata%Nnode_LCMeshFace,2,2) )
    !         unity_fac(:) = 1.0_RP

    !         call push_localsendbuf( lcfpos3D,                          &
    !           fpos3D, f, this%is_f(f,n), commdata%Nnode_LCMeshFace, 2, &
    !           lcmesh, elem )
            
    !         do varid=this%sfield_num+2*this%hvfield_num+1, this%field_num_tot-3, 4
    !           call CubedSphereCoordCnv_LonLat2CSVec( &
    !             lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),          &
    !             commdata%Nnode_LCMeshFace,                                           &
    !             commdata%recv_buf(:,varid), commdata%recv_buf(:,varid+1),            &
    !             tmp1_htensor3D(:,1,1), tmp1_htensor3D(:,2,1)    )
    !           call CubedSphereCoordCnv_LonLat2CSVec( &
    !             lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),          &
    !             commdata%Nnode_LCMeshFace,                                           &
    !             commdata%recv_buf(:,varid+2), commdata%recv_buf(:,varid+3),          &
    !             tmp1_htensor3D(:,1,2), tmp1_htensor3D(:,2,2)    )
    !           call CubedSphereCoordCnv_LonLat2CSVec( &
    !             lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),          &
    !             commdata%Nnode_LCMeshFace,                                           &
    !             tmp1_htensor3D(:,1,1), tmp1_htensor3D(:,1,2),                        &
    !             this%recv_buf(irs:ire,varid,n), this%recv_buf(irs:ire,varid+2,n)     )
    !           call CubedSphereCoordCnv_LonLat2CSVec( &
    !             lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),          &
    !             commdata%Nnode_LCMeshFace,                                           &
    !             tmp1_htensor3D(:,2,1), tmp1_htensor3D(:,2,2),                        &
    !             this%recv_buf(irs:ire,varid+1,n), this%recv_buf(irs:ire,varid+3,n)   )
    !         end do
    !         deallocate( lcfpos3D, unity_fac, tmp1_htensor3D )
    !       end if          
    !     end if

    !     irs = ire + 1
    !   end do

    !   deallocate( fpos3D )
    ! end do 

    return
  end subroutine MeshFieldCommCubedSphereDom3D_exchange

!----------------------------

!OCL SERIAL
  subroutine push_localsendbuf( lc_send_buf, send_buf, s_faceID, is, Nnode_LCMeshFace, var_num, &
      mesh3D, elem3D )
    implicit none

    integer, intent(in) :: var_num
    integer, intent(in) ::  Nnode_LCMeshFace
    real(RP), intent(inout) :: lc_send_buf(Nnode_LCMeshFace,var_num)
    real(RP), intent(in) :: send_buf(bufsize_per_field,var_num)  
    integer, intent(in) :: s_faceID, is
    type(LocalMesh3D), pointer :: mesh3D
    type(ElementBase3D), intent(in) :: elem3D

    integer :: ie
    !-----------------------------------------------------------------------------

    ie = is + Nnode_LCMeshFace - 1
    if ( s_faceID > 0 ) then
      lc_send_buf(:,:) = send_buf(is:ie,:)   
    else
      call revert_hori( lc_send_buf, send_buf(is:ie,:), mesh3D, elem3D, mesh3D%lcmesh2D )
    end if 
    
    return
  contains
    subroutine revert_hori( revert, ori, mesh, e3D, mesh2D )
      implicit none
      type(LocalMesh3D), intent(in) :: mesh
      type(ElementBase3D), intent(in) :: e3D
      type(LocalMesh2D), intent(in) :: mesh2D
      real(RP), intent(out) :: revert(e3D%Nnode_h1D,e3D%Nnode_v,e3D%Nnode_h1D,mesh2D%NeX,mesh%NeZ, var_num)
      real(RP), intent(in)  :: ori(e3D%Nnode_h1D,e3D%Nnode_v,e3D%Nnode_h1D,mesh2D%NeX,mesh%NeZ, var_num)
      
      integer :: p1, p3, ph, i, k, n
      integer :: i_, p1_
      !------------------------------------------------------------------------
      
      do n=1, var_num
      do k=1, mesh%NeZ
      do i=1, mesh2D%NeX
        i_ = mesh2D%NeX - i + 1
        do ph=1, e3D%Nnode_h1D
        do p3=1, e3D%Nnode_v
        do p1=1, e3D%Nnode_h1D
          p1_ = e3D%Nnode_h1D - p1 + 1
          revert(p1,p3,ph,i,k,n) = ori(p1_,p3,ph,i_,k,n)
        end do
        end do
        end do
      end do
      end do
      end do

      return
    end subroutine revert_hori
  end subroutine push_localsendbuf

!OCL SERIAL
  subroutine extract_boundary_data3D(var, refElem, mesh, buf)
    implicit none
  
    class(ElementBase3D), intent(in) :: refElem
    class(LocalMesh3D), intent(in) :: mesh
    real(RP), intent(in) :: var(refElem%Np * mesh%Ne)
    real(RP), intent(inout) :: buf( 2*(mesh%NeX + mesh%NeY)*mesh%NeZ*refElem%Nfp_h*refElem%Nnode_h1D &
                                    + 2*mesh%NeX*mesh%NeY*refElem%Nfp_v )
    !-----------------------------------------------------------------------------

    buf(:) = var(mesh%VmapB(:))
    return
  end subroutine extract_boundary_data3D

!OCL SERIAL
  subroutine extract_bounddata3D( var, elem, mesh, buf )
    implicit none
  
    type(ElementBase3D), intent(in) :: elem
    type(LocalMesh3D), intent(in) :: mesh
    real(DP), intent(in) :: var(elem%Np * mesh%NeA)
    real(DP), intent(inout) :: buf( 2*(mesh%NeX + mesh%NeY)*mesh%NeZ*elem%Nfp_h*elem%Nnode_h1D &
                                   + 2*mesh%NeX*mesh%NeY*elem%Nfp_v )

    integer :: f
    integer :: iso, isso, i, ii
    integer :: ke, p, ph
    !------------------------------------------------------------
  
    !--
    f = 1
    iso = 0; isso = 0
    do ke=1, mesh%NeX * mesh%NeZ
      do ph=1, elem%Nnode_h1D
        do p=1, elem%Nfp_h
          i  = iso  + p                     + (ke-1)*elem%Nfp_h
          ii = isso + p + (ph-1)*elem%Nfp_h + (ke-1)*elem%Nfp_h*elem%Nnode_h1D
          buf(ii) = var(mesh%VMapB(i) + (ph-1)*elem%Nnode_h1D)
        end do
      end do
    end do

    !--
    f = 2
    iso = iso + mesh%NeX*mesh%NeZ * elem%Nfp_h
    isso = isso + mesh%NeX*mesh%NeZ * elem%Nfp_h*elem%Nnode_h1D
    do ke=1, mesh%NeY * mesh%NeZ
      do ph=1, elem%Nnode_h1D
        do p=1, elem%Nfp_h
          i  = iso  + p                     + (ke-1)*elem%Nfp_h
          ii = isso + p + (ph-1)*elem%Nfp_h + (ke-1)*elem%Nfp_h*elem%Nnode_h1D
          buf(ii) = var(mesh%VMapB(i) - (ph-1))
        end do
      end do
    end do

    !--
    f = 3
    iso = iso + mesh%NeY*mesh%NeZ * elem%Nfp_h
    isso = isso + mesh%NeY*mesh%NeZ * elem%Nfp_h*elem%Nnode_h1D
    do ke=1, mesh%NeX * mesh%NeZ
      do ph=1, elem%Nnode_h1D
        do p=1, elem%Nfp_h
          i  = iso  + p                     + (ke-1)*elem%Nfp_h
          ii = isso + p + (ph-1)*elem%Nfp_h + (ke-1)*elem%Nfp_h*elem%Nnode_h1D
          buf(ii) = var(mesh%VMapB(i) - (ph-1)*elem%Nnode_h1D)
        end do
      end do
    end do

    !--
    f = 4
    iso = iso + mesh%NeX*mesh%NeZ * elem%Nfp_h
    isso = isso + mesh%NeX*mesh%NeZ * elem%Nfp_h*elem%Nnode_h1D
    do ke=1, mesh%NeY * mesh%NeZ
      do ph=1, elem%Nnode_h1D
        do p=1, elem%Nfp_h
          i  = iso  + p                     + (ke-1)*elem%Nfp_h
          ii = isso + p + (ph-1)*elem%Nfp_h + (ke-1)*elem%Nfp_h*elem%Nnode_h1D
          buf(ii) = var(mesh%VMapB(i) + (ph-1))
        end do
      end do
    end do

    !-- Upper & lower faces
    do f=1, 2
      iso = 2*(mesh%NeX + mesh%NeY)*mesh%NeZ*elem%Nfp_h + (f-1) * mesh%NeX*mesh%NeY*elem%Nfp_v
      isso = 2*(mesh%NeX + mesh%NeY)*mesh%NeZ*elem%Nfp_h*elem%Nnode_h1D + (f-1) * mesh%NeX*mesh%NeY*elem%Nfp_v
      do ke=1, mesh%NeX * mesh%NeY
        do p=1, elem%Nfp_v
          i  = iso  + p + (ke-1)*elem%Nfp_v
          ii = isso + p + (ke-1)*elem%Nfp_v
          buf(ii) = var(mesh%VMapB(i))
        end do
      end do
    end do

    return
  end subroutine extract_bounddata3D

  subroutine set_bounddata(buf, refElem, mesh, var)
    implicit none
    
    class(ElementBase), intent(in) :: refElem
    class(LocalMesh3D), intent(in) :: mesh
    real(RP), intent(in) :: buf(:)
    real(RP), intent(inout) :: var(refElem%Np * mesh%NeA)
    !-----------------------------------------------------------------------------

    var(refElem%Np*mesh%NeE+1:refElem%Np*mesh%NeE+size(buf)) = buf(:)
    return
  end subroutine set_bounddata  

!OCL SERIAL
  subroutine set_boundary_data3D_u1u2( buf_U, buf_V, &
    elem, mesh, G_ij,                                &
    u1, u2)
  
    implicit none
  
    type(ElementBase3D), intent(in) :: elem
    type(LocalMesh3D), intent(in) :: mesh
    real(DP), intent(in) :: buf_U(2*(mesh%NeX + mesh%NeY)*mesh%NeZ*elem%Nfp_h &
                                  + 2*mesh%NeX*mesh%NeY*elem%Nfp_v)
    real(DP), intent(in) :: buf_V(2*(mesh%NeX + mesh%NeY)*mesh%NeZ*elem%Nfp_h &
                                  + 2*mesh%NeX*mesh%NeY*elem%Nfp_v)
    real(DP), intent(in) :: G_ij(elem%Np * mesh%Ne,2,2)
    real(DP), intent(inout) :: u1(elem%Np * mesh%NeA)
    real(DP), intent(inout) :: u2(elem%Np * mesh%NeA)
    !------------------------------------------------------------

    u1(elem%Np*mesh%NeE+1:elem%Np*mesh%NeE+size(buf_U)) &
      = G_ij(mesh%VmapB,1,1) * buf_U(:) + G_ij(mesh%VmapB,1,2) * buf_V(:)
    u2(elem%Np*mesh%NeE+1:elem%Np*mesh%NeE+size(buf_U)) &
      = G_ij(mesh%VmapB,2,1) * buf_U(:) + G_ij(mesh%VmapB,2,2) * buf_V(:)
  
    return
  end subroutine set_boundary_data3D_u1u2

end module mod_user_sub_data