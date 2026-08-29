!-------------------------------------------------------------------------------
!> module FElib / Data / Communication in 3D cubed-sphere domain
!!
!! @par Description
!!      A module to manage data communication with 3D cubed-sphere domain for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
#include "scaleFElib.h"
module scale_meshfieldcomm_cubedspheredom3d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use scale_element_base, only: ElementBase, ElementBase3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_meshfield_base, only: MeshField3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase,                               &
    MeshFieldCommBase_Init, MeshFieldCommBase_Final, &
    MeshFieldCommBase_extract_bounddata,             &
    MeshFieldCommBase_extract_bounddata2,            &
    MeshFieldCommBase_extract_bounddata_2,           &
    MeshFieldCommBase_extract_bounddata_3,           &
    MeshFieldCommBase_set_bounddata,                 &
    MeshFieldCommBase_set_bounddata_3,               &
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

  !> Derived type to represent covariant vector components
  type :: VecCovariantComp
    type(MeshField3D), pointer :: u1 => null()
    type(MeshField3D), pointer :: u2 => null()
  end type    

  !> Base derived type to manage data communication with 3D cubed-sphere domain
  type, public, extends(MeshFieldCommBase) :: MeshFieldCommCubedSphereDom3D
    class(MeshCubedSphereDom3D), pointer :: mesh3d                          !< Pointer to an object representing 3D cubed-sphere computational mesh
    type(VecCovariantComp), allocatable :: vec_covariant_comp_ptrlist(:)
    integer, allocatable :: Nnode_LCMeshAllFace(:)

    integer :: haloSize_h1D !< Halo size for 1D horizontal direction
    integer :: haloSize_v   !< Halo size for vertical direction  
  contains
    procedure, public :: Init   => MeshFieldCommCubedSphereDom3D_Init
    procedure, public :: Put    => MeshFieldCommCubedSphereDom3D_put
    procedure, public :: Get    => MeshFieldCommCubedSphereDom3D_get
    procedure, public :: Exchange => MeshFieldCommCubedSphereDom3D_exchange  
    procedure, public :: SetCovariantVec => MeshFieldCommCubedSphereDom3D_set_covariantvec
    procedure, public :: Final => MeshFieldCommCubedSphereDom3D_Final
  end type MeshFieldCommCubedSphereDom3D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !
  private :: post_exchange_core
  private :: push_localsendbuf

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, parameter :: COMM_FACE_NUM = 6  !< Number of faces with data communication

contains
!> Initialize an object to manage data communication with 3D cubed-sphere computational mesh
  subroutine MeshFieldCommCubedSphereDom3D_Init( this, &
    sfield_num, hvfield_num, htensorfield_num, mesh3d, &
    haloSize_h1D, haloSize_v )

    use scale_meshutil_cubedsphere3d, only: MeshUtilCubedSphere3D_genPatchBoundaryMap_wide
    implicit none
    
    class(MeshFieldCommCubedSphereDom3D), intent(inout), target :: this
    integer, intent(in) :: sfield_num                         !< Number of scalar fields
    integer, intent(in) :: hvfield_num                        !< Number of horizontal vector fields
    integer, intent(in) :: htensorfield_num                   !< Number of horizontal vector fields
    class(MeshCubedSphereDom3D), intent(in), target :: mesh3d !< Object to manage a 3D cubed-sphere computational mesh
    integer, intent(in), optional :: haloSize_h1D             !< Halo size for 1D horizontal direction
    integer, intent(in), optional :: haloSize_v               !< Halo size for vertical direction

    type(LocalMesh3D), pointer :: lcmesh
    type(ElementBase3D), pointer :: elem
    integer :: n
    integer :: Nnode_LCMeshFace(COMM_FACE_NUM,mesh3d%LOCAL_MESH_NUM)
    !-----------------------------------------------------------------------------
    
    this%mesh3d => mesh3d
    lcmesh => mesh3d%lcmesh_list(1)
    elem => lcmesh%refElem3D

    !-
    if ( present(haloSize_h1D) ) then
      this%haloSize_h1D = haloSize_h1D
    else
      this%haloSize_h1D = 1
    end if
    if ( present(haloSize_v) ) then
      this%haloSize_v = haloSize_v
    else
      this%haloSize_v = 1
    end if

    !-
    allocate( this%Nnode_LCMeshAllFace(mesh3d%LOCAL_MESH_NUM) )    
    allocate( this%VMapB_size(this%mesh3d%LOCAL_MESH_NUM) )

    this%bufsize_per_field =  2*(lcmesh%NeX + lcmesh%NeY)*lcmesh%NeZ*elem%Nfp_h*this%haloSize_h1D &
                            + 2*lcmesh%NeX*lcmesh%NeY*elem%Nfp_v*this%haloSize_v

    do n=1, this%mesh3d%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      Nnode_LCMeshFace(:,n) = &
          (/ lcmesh%NeX, lcmesh%NeY, lcmesh%NeX, lcmesh%NeY, 0, 0 /) * lcmesh%NeZ * lcmesh%refElem3D%Nfp_h*this%haloSize_h1D &
        + (/ 0, 0, 0, 0, 1, 1 /) * lcmesh%NeX*lcmesh%NeY * lcmesh%refElem3D%Nfp_v*this%haloSize_v

      this%Nnode_LCMeshAllFace(n) = sum(Nnode_LCMeshFace(:,n))        
    end do
    !$acc enter data copyin(this%Nnode_LCMeshAllFace)

    call MeshFieldCommBase_Init( this, sfield_num, hvfield_num, htensorfield_num, this%bufsize_per_field, COMM_FACE_NUM, Nnode_LCMeshFace, mesh3d )  
  
    !-
    if ( this%haloSize_h1D > 1 .or. this%haloSize_v > 1) then
      this%use_vmap_wide_flag = .true.
      allocate( this%VMapB2(this%bufsize_per_field) )

      lcmesh => this%mesh3d%lcmesh_list(1)      
      call MeshUtilCubedSphere3D_genPatchBoundaryMap_wide( this%VMapB2, &
        lcmesh%VMapB, this%haloSize_h1D, this%haloSize_v,    &
        lcmesh%NeX, lcmesh%NeY, lcmesh%NeZ,                  &
        elem%Nfp_h, elem%Nfp_v, elem%Nnode_h1D, elem%Nnode_v )
      !$acc enter data copyin(this%VMapB2)
    else
      this%use_vmap_wide_flag = .false.
    end if
    
    do n=1, this%mesh3d%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      if ( this%use_vmap_wide_flag ) then
        this%VMapB_size(n) = size(this%VMapB2)        
      else
        this%VMapB_size(n) = size(lcmesh%VMapB)
      end if
    end do
    !$acc enter data copyin(this%VMapB_size)

    !-
    if (hvfield_num > 0) then
      allocate( this%vec_covariant_comp_ptrlist(hvfield_num) )
    end if

    return
  end subroutine MeshFieldCommCubedSphereDom3D_Init

!> Finalize an object to manage data communication with 3D cubed-sphere computational mesh
  subroutine MeshFieldCommCubedSphereDom3D_Final( this )

    implicit none
    
    class(MeshFieldCommCubedSphereDom3D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    !$acc exit data delete(this%Nnode_LCMeshAllFace, this%VMapB_size)
    deallocate( this%Nnode_LCMeshAllFace, this%VMapB_size )

    if ( this%hvfield_num > 0 ) then
      deallocate( this%vec_covariant_comp_ptrlist )
    end if

    call MeshFieldCommBase_Final( this )
    return
  end subroutine MeshFieldCommCubedSphereDom3D_Final

  !> Register objects to manage covariant vector components
  subroutine MeshFieldCommCubedSphereDom3D_set_covariantvec( &
    this, hvfield_ID, u1, u2  )
    implicit none
    class(MeshFieldCommCubedSphereDom3D), intent(inout) :: this
    integer, intent(in) :: hvfield_ID
    type(MeshField3D), intent(in), target :: u1
    type(MeshField3D), intent(in), target :: u2
    !--------------------------------------------------------------

    this%vec_covariant_comp_ptrlist(hvfield_ID)%u1 => u1
    this%vec_covariant_comp_ptrlist(hvfield_ID)%u2 => u2
    return
  end subroutine MeshFieldCommCubedSphereDom3D_set_covariantvec

!> Put field data into temporary buffers
  subroutine MeshFieldCommCubedSphereDom3D_put(this, field_list, varid_s)
    implicit none
    class(MeshFieldCommCubedSphereDom3D), intent(inout) :: this
    type(MeshFieldContainer), intent(in) :: field_list(:)  !< Array of objects with 3D mesh field
    integer, intent(in) :: varid_s                         !< Start index with variables when field_list(1) is written to buffers for data communication
  
    integer :: i
    integer :: n
    type(Localmesh3D), pointer :: lcmesh
    !-----------------------------------------------------------------------------
    ! call PROF_rapstart( 'comm_put', 2)

    if ( this%use_vmap_wide_flag ) then
      call MeshFieldCommBase_extract_bounddata_3( field_list, 3, varid_s, this%mesh3d%lcmesh_list, this%VMapB2, this%VMapB_size(1), & ! (in)
        this%send_buf ) ! (out)
    else
      call MeshFieldCommBase_extract_bounddata_2( &
        field_list, 3, varid_s, this%mesh3d%lcmesh_list, size(this%mesh3d%lcmesh_list(1)%VMapB), & !(in)
        this%send_buf ) ! (out)
    end if
    ! call PROF_rapend( 'comm_put', 2)
    return
  end subroutine MeshFieldCommCubedSphereDom3D_put

!> Extract field data from temporary buffers
  subroutine MeshFieldCommCubedSphereDom3D_get(this, field_list, varid_s)
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_wait_core
    implicit none
    
    class(MeshFieldCommCubedSphereDom3D), intent(inout) :: this
    type(MeshFieldContainer), intent(inout) :: field_list(:)  !< Array of objects with 3D mesh field
    integer, intent(in) :: varid_s                            !< Start index with variables when field_list(1) is written to buffers for data communication

    integer :: i
    integer :: n
    integer :: ke, ke2D
    integer :: p
    type(Localmesh3D), pointer :: lcmesh
    class(ElementBase3D), pointer :: elem

    integer :: varnum
    integer :: varid_e
    integer :: varid_vec_s

    integer, allocatable :: IndexH2Dto3D(:)
    real(RP), allocatable :: G_ij(:,:,:,:)

    integer :: ke_z
    integer, allocatable :: VMapB(:)
    !-----------------------------------------------------------------------------

    ! call PROF_rapstart( 'comm_get', 2)

    varnum = size(field_list) 

    !--
    if ( this%call_wait_flag_sub_get ) then
      ! This workflow should be reconsidered in near future.
      ! The coordinate convresions for horizontal vector/tensor fields in post_exchange_core are applied for the received data (this%recv_buf).
      ! However, MeshFieldCommBase_wait_core does not store the communication data into this%recv_buf, but directly into the field_list.
      call MeshFieldCommBase_wait_core( this, this%commdata_list, &
        field_list, 3, varid_s, this%mesh3d%lcmesh_list )
      call post_exchange_core( this )
    else
      do i=1, varnum
      do n=1, this%mesh3d%LOCAL_MESH_NUM
        lcmesh => this%mesh3d%lcmesh_list(n)
        if ( this%use_vmap_wide_flag ) then
          call MeshFieldCommBase_set_bounddata_3( this%recv_buf(:,varid_s+i-1,n), lcmesh%refElem, lcmesh, & ! (in)
            this%VMapB2, this%VMapB_size(1),                                                              & ! (in)          
            field_list(i)%field3d%local(n)%val )                                                            !(out)
        else
          call MeshFieldCommBase_set_bounddata( this%recv_buf(:,varid_s+i-1,n), lcmesh%refElem, lcmesh, & !(in)
            field_list(i)%field3d%local(n)%val )                                                          !(out)
        end if
      end do
      end do
    end if

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

            allocate( G_ij(elem%Np,lcmesh%Ne,2,2), IndexH2Dto3D(elem%Np) )
            IndexH2Dto3D(:) = elem%IndexH2Dto3D(:)

            if ( this%use_vmap_wide_flag ) then
              allocate(VMapB(this%VMapB_size(n)))
              VMapB(:) = this%VMapB2(:)
            else
              allocate(VMapB(size(lcmesh%VMapB)))
              VMapB(:) = lcmesh%VMapB(:)
            end if
            !$acc enter data create(G_ij) copyin(VMapB, IndexH2Dto3D) async(1)

            !$omp parallel do private(ke2D)
            !$acc parallel loop collapse(2) async(1)
            do ke=lcmesh%NeS, lcmesh%NeE
            do p=1, elem%Np
              ke2D = lcmesh%EMap3Dto2D(ke)
              G_ij(p,ke,1,1) = lcmesh%G_ij(IndexH2Dto3D(p),ke2D,1,1)
              G_ij(p,ke,2,1) = lcmesh%G_ij(IndexH2Dto3D(p),ke2D,2,1)
              G_ij(p,ke,1,2) = lcmesh%G_ij(IndexH2Dto3D(p),ke2D,1,2)
              G_ij(p,ke,2,2) = lcmesh%G_ij(IndexH2Dto3D(p),ke2D,2,2)
            end do
            end do

            call set_boundary_data3D_u1u2( &
              this%recv_buf(:,varid_vec_s,n), this%recv_buf(:,varid_vec_s+1,n), & ! (in)
              this%bufsize_per_field, lcmesh%refElem3D, lcmesh, G_ij, VMapB,    & ! (in)
              this%vec_covariant_comp_ptrlist(i)%u1%local(n)%val,               & ! (out)
              this%vec_covariant_comp_ptrlist(i)%u2%local(n)%val                ) ! (out)
            !$acc exit data delete(G_ij, VMapB, IndexH2Dto3D) async(1)
            deallocate( G_ij, VMapB, IndexH2Dto3D )
          end do
       end if
      end do
    end if

    !$acc wait(1)
    ! call PROF_rapend( 'comm_get', 2)
    return
  end subroutine MeshFieldCommCubedSphereDom3D_get

  !> Exchange field data between neighboring MPI processes
  !!
  !! @param do_wait Flag whether MPI_waitall is called and move tmp data of LocalMeshCommData object to a recv buffer
!OCL SERIAL
  subroutine MeshFieldCommCubedSphereDom3D_exchange( this, do_wait )
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_exchange_core,  &
      LocalMeshCommData
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatVec, &
      CubedSphereCoordCnv_LonLat2CSVec
    use scale_prc, only: PRC_abort
    implicit none
  
    class(MeshFieldCommCubedSphereDom3D), intent(inout), target :: this
    logical, intent(in), optional :: do_wait

    integer :: n, f
    integer :: varid
    integer :: i

    real(RP), allocatable :: fpos3D(:,:)
    real(RP), allocatable :: lcfpos3D(:,:)
    real(RP), allocatable :: unity_fac(:)
    real(RP), allocatable :: tmp_svec3D(:,:)
    real(RP), allocatable :: tmp1_htensor3D(:,:,:)
    real(RP), allocatable :: tmp2_htensor3D(:,:,:)
    
    class(ElementBase3D), pointer :: elem
    type(LocalMesh3D), pointer :: lcmesh
    type(LocalMeshCommData), pointer :: commdata

    integer, pointer :: VmapB(:)
    !-----------------------------------------------------------------------------

    ! call PROF_rapstart( 'comm_exchange_1', 2)

    if ( present(do_wait) ) then
      if ( do_wait == .false. ) then
        LOG_INFO("MeshFieldCommCubedSphereDom3D_exchange",*) "do_wait=False is not currently supported. Check!"
        call PRC_abort
      end if
    end if

    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      elem => lcmesh%refElem3D

      allocate( fpos3D(this%Nnode_LCMeshAllFace(n),2) )

      if ( this%use_vmap_wide_flag ) then
        VmapB => this%VMapB2
      else
        VmapB => lcmesh%VMapB
      end if

      !$acc data create(fpos3D) copyin(VmapB)

      call MeshFieldCommBase_extract_bounddata2( lcmesh%pos_en(:,:,1), VMapB, size(VMapB), elem%Np*lcmesh%Ne, fpos3D(:,1) )
      call MeshFieldCommBase_extract_bounddata2( lcmesh%pos_en(:,:,2), VMapB, size(VMapB), elem%Np*lcmesh%Ne, fpos3D(:,2) )

      do f=1, this%nfaces_comm
        commdata => this%commdata_list(f,n)
        call push_localsendbuf( commdata%send_buf,                 &  ! (inout)
          this%send_buf(:,:,n), commdata%s_faceID, this%is_f(f,n), &  ! (in)
          commdata%Nnode_LCMeshFace, this%bufsize_per_field,       &  ! (in)
          this%field_num_tot, lcmesh, this%HaloSize_h1D            )  ! (in)

        if ( commdata%s_panelID /= lcmesh%panelID &
            .and. ( this%hvfield_num > 0 .or. this%htensorfield_num > 0 ) ) then
          
          allocate( lcfpos3D(commdata%Nnode_LCMeshFace,2), unity_fac(commdata%Nnode_LCMeshFace) )
          unity_fac(:) = 1.0_RP
          !$acc enter data create(lcfpos3D) copyin(unity_fac) async(1)

          call push_localsendbuf( lcfpos3D,                            &
            fpos3D, commdata%s_faceID, this%is_f(f,n),                 &
            commdata%Nnode_LCMeshFace, this%Nnode_LCMeshAllFace(n), 2, &
            lcmesh, this%HaloSize_h1D )
        end if
          
        if ( commdata%s_panelID /= lcmesh%panelID ) then
          if ( this%hvfield_num > 0 ) then
            allocate( tmp_svec3D(commdata%Nnode_LCMeshFace,2) )
            !$acc enter data create(tmp_svec3D) async(1)

            do varid=this%sfield_num+1, this%sfield_num+2*this%hvfield_num-1, 2
              !$acc parallel loop present(tmp_svec3D, commdata%send_buf) async(1)
              do i=1, size(tmp_svec3D,1)
                tmp_svec3D(i,1) = commdata%send_buf(i,varid  )
                tmp_svec3D(i,2) = commdata%send_buf(i,varid+1)
              end do
  
              call CubedSphereCoordCnv_CS2LonLatVec( &
                lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),  & ! (in)
                commdata%Nnode_LCMeshFace,                                   & ! (in)
                tmp_svec3D(:,1), tmp_svec3D(:,2),                            & ! (in)
                commdata%send_buf(:,varid), commdata%send_buf(:,varid+1),    & ! (out) 
                gpu_async_id=1 ) ! (in)

            end do
            !$acc exit data delete(tmp_svec3D) async(1)
            deallocate( tmp_svec3D )
          end if

          if ( this%htensorfield_num > 0 ) then
            allocate( tmp1_htensor3D(commdata%Nnode_LCMeshFace,2,2) )
            allocate( tmp2_htensor3D(commdata%Nnode_LCMeshFace,2,2) )
            !$acc enter data create(tmp1_htensor3D, tmp2_htensor3D) async(1)
            
            do varid=this%sfield_num+2*this%hvfield_num+1, this%field_num_tot-3, 4
              !$acc parallel loop present(tmp1_htensor3D, commdata%send_buf) async(1)
              do i=1, size(tmp1_htensor3D,1)
                tmp1_htensor3D(i,1,1) = commdata%send_buf(i,varid  )
                tmp1_htensor3D(i,2,1) = commdata%send_buf(i,varid+1)
                tmp1_htensor3D(i,1,2) = commdata%send_buf(i,varid+2)
                tmp1_htensor3D(i,2,2) = commdata%send_buf(i,varid+3)
              end do
  
              call CubedSphereCoordCnv_CS2LonLatVec( &
                lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),        & ! (in)
                commdata%Nnode_LCMeshFace,                                         & ! (in)
                tmp1_htensor3D(:,1,1), tmp1_htensor3D(:,2,1),                      & ! (in)
                tmp2_htensor3D(:,1,1), tmp2_htensor3D(:,2,1),                      & ! (out)
                gpu_async_id=1 ) ! (in)
              call CubedSphereCoordCnv_CS2LonLatVec( &
                lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),        & ! (in)
                commdata%Nnode_LCMeshFace,                                         & ! (in)
                tmp1_htensor3D(:,1,2), tmp1_htensor3D(:,2,2),                      & ! (in)
                tmp2_htensor3D(:,1,2), tmp2_htensor3D(:,2,2),                      & ! (out)
                gpu_async_id=1 ) ! (in)
              call CubedSphereCoordCnv_CS2LonLatVec( &
                lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),        & ! (in)
                commdata%Nnode_LCMeshFace,                                         & ! (in)
                tmp2_htensor3D(:,1,1), tmp2_htensor3D(:,1,2),                      & ! (in)
                commdata%send_buf(:,varid), commdata%send_buf(:,varid+2),          & ! (out)
                gpu_async_id=1 ) ! (in)
              call CubedSphereCoordCnv_CS2LonLatVec( &
                lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac(:),        & ! (in)
                commdata%Nnode_LCMeshFace,                                         & ! (in)
                tmp2_htensor3D(:,2,1), tmp2_htensor3D(:,2,2),                      & ! (in)
                commdata%send_buf(:,varid+1), commdata%send_buf(:,varid+3),        & ! (out)
                gpu_async_id=1 ) ! (in)
            end do
            !$acc exit data delete(tmp1_htensor3D, tmp2_htensor3D) async(1)
            deallocate( tmp1_htensor3D, tmp2_htensor3D )
          end if
        
        end if
        
        if ( allocated(lcfpos3D) ) then
          !$acc exit data delete(lcfpos3D, unity_fac) async(1)
          deallocate( lcfpos3D, unity_fac )
        end if

      end do
      !$acc end data
      deallocate( fpos3D )
    end do
    !$acc wait(1)
    ! call PROF_rapend( 'comm_exchange_1', 2)
    
    !-----------------------

    ! call PROF_rapstart( 'comm_exchange_2', 2)
    call MeshFieldCommBase_exchange_core(this, this%commdata_list, do_wait )
    ! call PROF_rapend( 'comm_exchange_2', 2)

    ! call PROF_rapstart( 'comm_exchange_3', 2)
    if ( .not. this%call_wait_flag_sub_get ) &
      call post_exchange_core( this )

    ! call PROF_rapend( 'comm_exchange_3', 2)

    return
  end subroutine MeshFieldCommCubedSphereDom3D_exchange

!----------------------------

!OCL SERIAL
  subroutine post_exchange_core( this )
    use scale_meshfieldcomm_base, only: &
      LocalMeshCommData
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSVec
    implicit none

    class(MeshFieldCommCubedSphereDom3D), intent(inout), target :: this

    integer :: n, f
    integer :: varid

    real(RP), allocatable :: fpos3D(:,:)
    real(RP), allocatable :: lcfpos3D(:,:)
    real(RP), allocatable :: unity_fac(:)
    real(RP), allocatable :: tmp1_htensor3D(:,:,:)
    
    class(ElementBase3D), pointer :: elem
    type(LocalMesh3D), pointer :: lcmesh
    type(LocalMeshCommData), pointer :: commdata

    integer :: irs, ire
    integer, pointer :: VmapB(:)
    !-----------------------------------------------------------------------------

    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh3d%lcmesh_list(n)
      elem => lcmesh%refElem3D

      allocate( fpos3D(this%Nnode_LCMeshAllFace(n),2) )
      if ( this%use_vmap_wide_flag ) then
        VmapB => this%VMapB2
      else
        VmapB => lcmesh%VMapB
      end if
      !$acc data create(fpos3D) copyin(VMapB)

      call MeshFieldCommBase_extract_bounddata2( lcmesh%pos_en(:,:,1), VMapB, size(VMapB), elem%Np*lcmesh%Ne, fpos3D(:,1) )
      call MeshFieldCommBase_extract_bounddata2( lcmesh%pos_en(:,:,2), VMapB, size(VMapB), elem%Np*lcmesh%Ne, fpos3D(:,2) )

      irs = 1
      do f=1, this%nfaces_comm
        commdata => this%commdata_list(f,n)
        ire = irs + commdata%Nnode_LCMeshFace - 1

        if ( commdata%s_panelID /= lcmesh%panelID &
            .and. ( this%hvfield_num > 0 .or. this%htensorfield_num > 0 ) ) then
          
          allocate( lcfpos3D(commdata%Nnode_LCMeshFace,2), unity_fac(commdata%Nnode_LCMeshFace) )
          unity_fac(:) = 1.0_RP
          !$acc enter data create(lcfpos3D) copyin(unity_fac) async(1)

          call push_localsendbuf( lcfpos3D,                            &
            fpos3D, commdata%s_faceID, this%is_f(f,n),                 &
            commdata%Nnode_LCMeshFace, this%Nnode_LCMeshAllFace(n), 2, &
            lcmesh, this%HaloSize_h1D )
        end if

        if ( commdata%s_panelID /= lcmesh%panelID ) then
          if ( this%hvfield_num > 0 ) then            
            do varid=this%sfield_num+1, this%sfield_num+2*this%hvfield_num-1, 2
              call CubedSphereCoordCnv_LonLat2CSVec( &
                lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac,           & ! (in)
                commdata%Nnode_LCMeshFace,                                         & ! (in)
                commdata%recv_buf(:,varid), commdata%recv_buf(:,varid+1),          & ! (in)
                this%recv_buf(irs:ire,varid,n), this%recv_buf(irs:ire,varid+1,n),  & ! (out) 
                gpu_async_id=1 ) ! (in)
            end do
          end if

          if ( this%htensorfield_num > 0 ) then
            allocate( tmp1_htensor3D(commdata%Nnode_LCMeshFace,2,2) )
            !$acc enter data create(tmp1_htensor3D) async(1)
            
            do varid=this%sfield_num+2*this%hvfield_num+1, this%field_num_tot-3, 4
              call CubedSphereCoordCnv_LonLat2CSVec( &
                lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac,            & ! (in)
                commdata%Nnode_LCMeshFace,                                          & ! (in)
                commdata%recv_buf(:,varid), commdata%recv_buf(:,varid+1),           & ! (in)
                tmp1_htensor3D(:,1,1), tmp1_htensor3D(:,2,1),                       & ! (out)
                gpu_async_id=1 ) ! (in)
              call CubedSphereCoordCnv_LonLat2CSVec( &
                lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac,            & ! (in)
                commdata%Nnode_LCMeshFace,                                          & ! (in)
                commdata%recv_buf(:,varid+2), commdata%recv_buf(:,varid+3),         & ! (in)
                tmp1_htensor3D(:,1,2), tmp1_htensor3D(:,2,2),                       & ! (out)
                gpu_async_id=1 ) ! (in)
              call CubedSphereCoordCnv_LonLat2CSVec( &
                lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac,            & ! (in)
                commdata%Nnode_LCMeshFace,                                          & ! (in)
                tmp1_htensor3D(:,1,1), tmp1_htensor3D(:,1,2),                       & ! (in)
                this%recv_buf(irs:ire,varid,n), this%recv_buf(irs:ire,varid+2,n),   & ! (out)
                gpu_async_id=1 ) ! (in)
              call CubedSphereCoordCnv_LonLat2CSVec( &
                lcmesh%panelID, lcfpos3D(:,1), lcfpos3D(:,2), unity_fac,            & ! (in)
                commdata%Nnode_LCMeshFace,                                          & ! (in)
                tmp1_htensor3D(:,2,1), tmp1_htensor3D(:,2,2),                       & ! (in)
                this%recv_buf(irs:ire,varid+1,n), this%recv_buf(irs:ire,varid+3,n), & ! (out)
                gpu_async_id=1 ) ! (in)
            end do
            !$acc exit data delete(tmp1_htensor3D) async(1)
            deallocate( tmp1_htensor3D )
          end if          
        end if

        irs = ire + 1
        if ( allocated(lcfpos3D) ) then
          !$acc exit data delete(lcfpos3D, unity_fac) async(1)
          deallocate( lcfpos3D, unity_fac )
        end if
      end do

      !$acc end data
      deallocate( fpos3D )
    end do 
    !$acc wait(1)

    return
  end subroutine post_exchange_core

  !> Push temporary buffer of the local send data to the communication buffer
!OCL SERIAL
  subroutine push_localsendbuf( lc_send_buf, &
    send_buf, s_faceID, is, Nnode_LCMeshFace, bufsize_per_field, var_num, &
    lcmesh, haloSize_h1D )
    implicit none

    integer, intent(in) :: var_num
    integer, intent(in) ::  Nnode_LCMeshFace
    integer, intent(in) :: bufsize_per_field
    real(RP), intent(inout) :: lc_send_buf(Nnode_LCMeshFace,var_num)
    real(RP), intent(in) :: send_buf(bufsize_per_field,var_num)  
    integer, intent(in) :: s_faceID, is
    type(LocalMesh3D), pointer :: lcmesh
    integer, intent(in) :: haloSize_h1D

    integer :: vid, i
    integer :: Nnode_v, Nnode_h1D, Nfp_h
    !-----------------------------------------------------------------------------

    if ( s_faceID > 0 ) then
      !$omp parallel do
      !$acc parallel loop collapse(2) present(lc_send_buf, send_buf) async(1)
      do vid=1, var_num
      do i=1, Nnode_LCMeshFace
        lc_send_buf(i,vid) = send_buf(is+i-1,vid)
      end do
      end do
    else if ( -5 < s_faceID .and. s_faceID < 0) then

      Nfp_h = lcmesh%refElem3D%Nfp_h
      Nnode_h1D = lcmesh%refElem3D%Nnode_h1D
      Nnode_v = lcmesh%refElem3D%Nnode_v

      call revert_hori( lc_send_buf, &
        send_buf, is, Nnode_LCMeshFace/(Nfp_h * haloSize_h1D * lcmesh%NeZ), &
        lcmesh%NeZ, haloSize_h1D )      
    end if 
    
    return
  contains
    !> Revert the horizontal ordering of the local send buffer for negative s_faceID
    subroutine revert_hori(revert, ori, is_, Ne_h1D, NeZ, haloSize_h1D_ )
      integer, intent(in) :: is_
      integer, intent(in) :: Ne_h1D, NeZ
      integer, intent(in) :: haloSize_h1D_
      real(RP), intent(out) :: revert(Nnode_h1D,Nnode_v,haloSize_h1D_, Ne_h1D,NeZ, var_num)
      real(RP), intent(in)  :: ori(bufsize_per_field, var_num)
      
      integer :: p1, p3, ph, i, k, n
      integer :: i_, p1_
      integer :: id_src
      !-----------------------------------------------------------------------------
      
      !$acc parallel loop collapse(4) present(revert, ori) async(1)
      do n=1, var_num
      do k=1, NeZ
      do i=1, Ne_h1D
      do ph=1, haloSize_h1D_ 
        i_ = Ne_h1D - i + 1
        do p3=1, Nnode_v
        do p1=1, Nnode_h1D
          p1_ = Nnode_h1D - p1 + 1
          ! source ordering:　(p1_,p3,ph,i_,k)
          id_src = is_ - 1 &
            + p1_ &
            + (p3-1) * Nnode_h1D             &
            + (ph-1) * Nfp_h                 &
            + (i_-1) * Nfp_h * haloSize_h1D_ &
            + (k -1) * Nfp_h * haloSize_h1D_ * Ne_h1D
          
          revert(p1,p3,ph,i,k,n) = ori(id_src,n)
        end do
        end do
      end do
      end do
      end do
      end do
      return
    end subroutine revert_hori    
  end subroutine push_localsendbuf

  subroutine set_boundary_data3D_u1u2( buf_U, buf_V, &
    bufsize_per_field, elem, mesh, G_ij, VMapB,      &
    u1, u2)
  
    implicit none
    integer, intent(in) :: bufsize_per_field
    type(ElementBase3D), intent(in) :: elem
    type(LocalMesh3D), intent(in) :: mesh
    real(DP), intent(in) :: buf_U(bufsize_per_field)
    real(DP), intent(in) :: buf_V(bufsize_per_field)
    real(DP), intent(in) :: G_ij(elem%Np * mesh%Ne,2,2)
    integer, intent(in) :: VMapB(bufsize_per_field)
    real(DP), intent(inout) :: u1(elem%Np * mesh%NeA)
    real(DP), intent(inout) :: u2(elem%Np * mesh%NeA)

    integer :: ii, iis
    integer :: vmap_b
    !------------------------------------------------------------

#ifdef _OPENACC
    iis = elem%Np*mesh%NeE
    !$acc parallel loop present(buf_U, buf_V, G_ij, u1, u2, VMapB)
    do ii=1, bufsize_per_field
      vmap_b = VMapB(ii)
      u1(iis+ii) = G_ij(vmap_b,1,1) * buf_U(ii) + G_ij(vmap_b,1,2) * buf_V(ii)
      u2(iis+ii) = G_ij(vmap_b,2,1) * buf_U(ii) + G_ij(vmap_b,2,2) * buf_V(ii)
    end do
#else
    u1(elem%Np*mesh%NeE+1:elem%Np*mesh%NeE+bufsize_per_field) &
      = G_ij(VMapB,1,1) * buf_U(:) + G_ij(VMapB,1,2) * buf_V(:)
    u2(elem%Np*mesh%NeE+1:elem%Np*mesh%NeE+bufsize_per_field) &
      = G_ij(VMapB,2,1) * buf_U(:) + G_ij(VMapB,2,2) * buf_V(:)
#endif
    return
  end subroutine set_boundary_data3D_u1u2

end module scale_meshfieldcomm_cubedspheredom3d