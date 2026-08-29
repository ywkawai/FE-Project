!-------------------------------------------------------------------------------
!> module FElib / Data / Communication in 2D cubed-sphere domain
!!
!! @par Description
!!      A module to manage data communication with 2D cubed-sphere domain for element-based methods
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

  use scale_element_base, only: ElementBase, ElementBase2D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  use scale_meshfield_base, only: MeshField2D
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
  use scale_localmesh_2d, only: LocalMesh2D
   
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !> Derived type to represent covariant vector components
  type :: VecCovariantComp
    type(MeshField2D), pointer :: u1 => null()
    type(MeshField2D), pointer :: u2 => null()
  end type    

  !> Base derived type to manage data communication with 2D cubed-sphere domain
  type, public, extends(MeshFieldCommBase) :: MeshFieldCommCubedSphereDom2D
    class(MeshCubedSphereDom2D), pointer :: mesh2d                         !< Pointer to an object representing 2D cubed-sphere computational mesh
    type(VecCovariantComp), allocatable :: vec_covariant_comp_ptrlist(:)
    integer, allocatable :: Nnode_LCMeshAllFace(:)

    integer :: haloSize_1D !< Halo size for 1D direction
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
  private :: post_exchange_core
  private :: push_localsendbuf
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, parameter :: COMM_FACE_NUM = 4

contains
!> Initialize an object to manage data communication with 2D cubed-sphere computational mesh
  subroutine MeshFieldCommCubedSphereDom2D_Init( this, &
    sfield_num, hvfield_num, htensorfield_num, mesh2d, &
    haloSize_1D )

    use scale_meshutil_cubedsphere2d, only: MeshUtilCubedSphere2D_genPatchBoundaryMap_wide
    implicit none
    
    class(MeshFieldCommCubedSphereDom2D), intent(inout) :: this
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    integer, intent(in) :: htensorfield_num
    class(MeshCubedSphereDom2D), intent(in), target :: mesh2d
    integer, intent(in), optional :: haloSize_1D
    
    type(LocalMesh2D), pointer :: lcmesh
    type(ElementBase2D), pointer :: elem
    integer :: n
    integer :: Nnode_LCMeshFace(COMM_FACE_NUM,mesh2d%LOCAL_MESH_NUM)
    !-----------------------------------------------------------------------------
    
    this%mesh2d => mesh2d
    lcmesh => mesh2d%lcmesh_list(1)
    elem => lcmesh%refElem2D

    !-
    if ( present(haloSize_1D) ) then
      this%haloSize_1D = haloSize_1D
    else
      this%haloSize_1D = 1
    end if

    !-
    allocate( this%VMapB_size(this%mesh2d%LOCAL_MESH_NUM) )

    this%bufsize_per_field = 2*(lcmesh%NeX + lcmesh%NeY) * elem%Nfp*this%haloSize_1D

    allocate( this%Nnode_LCMeshAllFace(mesh2d%LOCAL_MESH_NUM) )
    do n=1, this%mesh2d%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      Nnode_LCMeshFace(:,n) = (/ lcmesh%NeX, lcmesh%NeY, lcmesh%NeX, lcmesh%NeY /) * lcmesh%refElem2D%Nfp*this%haloSize_1D
      this%Nnode_LCMeshAllFace(n) = sum(Nnode_LCMeshFace(:,n))
    end do
    !$acc enter data copyin(this%Nnode_LCMeshAllFace)

    call MeshFieldCommBase_Init( this, sfield_num, hvfield_num, htensorfield_num, this%bufsize_per_field, COMM_FACE_NUM, Nnode_LCMeshFace, mesh2d)  

    !-
    if ( this%haloSize_1D > 1 ) then
      this%use_vmap_wide_flag = .true.
      allocate( this%VMapB2(this%bufsize_per_field) )

      lcmesh => this%mesh2d%lcmesh_list(1)      
      call MeshUtilCubedSphere2D_genPatchBoundaryMap_wide( this%VMapB2, &
        lcmesh%VMapB, this%haloSize_1D,                      &
        lcmesh%NeX, lcmesh%NeY,                              &
        elem%Nfp )
      !$acc enter data copyin(this%VMapB2)
    else
      this%use_vmap_wide_flag = .false.
    end if

    do n=1, this%mesh2d%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
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
  end subroutine MeshFieldCommCubedSphereDom2D_Init

!> Finalize an object to manage data communication with 2D cubed-sphere computational mesh
  subroutine MeshFieldCommCubedSphereDom2D_Final( this )
    implicit none
    
    class(MeshFieldCommCubedSphereDom2D), intent(inout) :: this
    !-----------------------------------------------------------------------------

    !$acc exit data delete(this%Nnode_LCMeshAllFace, this%VMapB_size)
    deallocate( this%Nnode_LCMeshAllFace, this%VMapB_size )

    if ( this%hvfield_num > 0 ) then
      deallocate( this%vec_covariant_comp_ptrlist )
    end if

    call MeshFieldCommBase_Final( this )

    return
  end subroutine MeshFieldCommCubedSphereDom2D_Final

!> Register objects to manage covariant vector components
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

!> Put field data into temporary buffers
  subroutine MeshFieldCommCubedSphereDom2D_put(this, field_list, varid_s)
    implicit none
    class(MeshFieldCommCubedSphereDom2D), intent(inout) :: this
    type(MeshFieldContainer), intent(in) :: field_list(:)
    integer, intent(in) :: varid_s
    !-----------------------------------------------------------------------------

!    call PROF_rapstart( 'comm_put', 2)

    if ( this%use_vmap_wide_flag ) then
       call MeshFieldCommBase_extract_bounddata_3( field_list, 2, varid_s, this%mesh2d%lcmesh_list, this%VMapB2, this%VMapB_size(1), this%send_buf )
    else
      call MeshFieldCommBase_extract_bounddata_2( &
        field_list, 2, varid_s, this%mesh2d%lcmesh_list, size(this%mesh2d%lcmesh_list(1)%VMapB), & !(in)
        this%send_buf ) ! (out)
    end if

!    call PROF_rapend( 'comm_put', 2)
    return
  end subroutine MeshFieldCommCubedSphereDom2D_put

!> Extract field data from temporary buffers
  subroutine MeshFieldCommCubedSphereDom2D_get(this, field_list, varid_s)
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_wait_core
    implicit none
    
    class(MeshFieldCommCubedSphereDom2D), intent(inout) :: this
    type(MeshFieldContainer), intent(inout) :: field_list(:)
    integer, intent(in) :: varid_s

    integer :: i
    integer :: n
    type(LocalMesh2D), pointer :: lcmesh

    integer :: varnum
    integer :: varid_e
    integer :: varid_vec_s

    integer, allocatable :: VMapB(:)
    !-----------------------------------------------------------------------------
!    call PROF_rapstart( 'comm_get', 2)
    varnum = size(field_list) 

    !--
    if ( this%call_wait_flag_sub_get ) then
      ! This workflow should be reconsidered in near future.
      ! The coordinate convresions for horizontal vector/tensor fields in post_exchange_core are applied for the received data (this%recv_buf).
      ! However, MeshFieldCommBase_wait_core does not store the communication data into this%recv_buf, but directly into the field_list.
      call MeshFieldCommBase_wait_core( this, this%commdata_list )
      call post_exchange_core( this )
    else
      do i=1, varnum
      do n=1, this%mesh2d%LOCAL_MESH_NUM
        lcmesh => this%mesh2d%lcmesh_list(n)
        if ( this%use_vmap_wide_flag ) then
          call MeshFieldCommBase_set_bounddata_3( this%recv_buf(:,varid_s+i-1,n), lcmesh%refElem, lcmesh, & ! (in)
            this%VMapB2, this%VMapB_size(1),                                                              & ! (in)
            field_list(i)%field2d%local(n)%val )                                                            ! (inout)
        else        
          call MeshFieldCommBase_set_bounddata( this%recv_buf(:,varid_s+i-1,n), lcmesh%refElem, lcmesh, & !(in)
            field_list(i)%field2d%local(n)%val )                                                         !(out)
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
        
          do n=1, this%mesh2d%LOCAL_MESH_NUM
            lcmesh => this%mesh2d%lcmesh_list(n)

            if ( this%use_vmap_wide_flag ) then
              allocate(VMapB(this%VMapB_size(n)))
              VMapB(:) = this%VMapB2(:)
            else
              allocate(VMapB(size(lcmesh%VMapB)))
              VMapB(:) = lcmesh%VMapB(:)
            end if
            !$acc enter data copyin(VMapB) async(1)
            call set_boundary_data2D_u1u2( &
              this%recv_buf(:,varid_vec_s,n), this%recv_buf(:,varid_vec_s+1,n),     & ! (in)
              this%bufsize_per_field, lcmesh%refElem2D, lcmesh, lcmesh%G_ij, VMapB, & ! (in)
              this%vec_covariant_comp_ptrlist(i)%u1%local(n)%val,                   & ! (out)
              this%vec_covariant_comp_ptrlist(i)%u2%local(n)%val                    ) ! (out)
            !$acc exit data delete(VMapB) async(1)
            deallocate( VMapB )
          end do
       end if
      end do
    end if
    !$acc wait(1)

!    call PROF_rapend( 'comm_get', 2)
    return
  end subroutine MeshFieldCommCubedSphereDom2D_get

  !> Exchange field data between neighboring MPI processes
  !!
  !! @param do_wait Flag whether MPI_waitall is called and move tmp data of LocalMeshCommData object to a recv buffer
!OCL SERIAL
  subroutine MeshFieldCommCubedSphereDom2D_exchange( this, do_wait )    
    use scale_meshfieldcomm_base, only: &
      MeshFieldCommBase_exchange_core,  &
      LocalMeshCommData
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatVec, &
      CubedSphereCoordCnv_LonLat2CSVec
    use scale_prc, only: PRC_abort
    implicit none
  
    class(MeshFieldCommCubedSphereDom2D), intent(inout), target :: this
    logical, intent(in), optional :: do_wait

    integer :: n, f
    integer :: varid
    integer :: i

    real(RP), allocatable :: fpos2D(:,:)
    real(RP), allocatable :: lcfpos2D(:,:)
    real(RP), allocatable :: unity_fac(:)
    real(RP), allocatable :: tmp_svec2D(:,:)
    
    class(ElementBase2D), pointer :: elem
    type(LocalMesh2D), pointer :: lcmesh
    type(LocalMeshCommData), pointer :: commdata

    integer, pointer :: VmapB(:)
    !-----------------------------------------------------------------------------

    ! call PROF_rapstart( 'comm_exchange', 2)

    if ( present(do_wait) ) then
      if ( .not. do_wait ) then
        LOG_INFO("MeshFieldCommCubedSphereDom2D_exchange",*) "do_wait=False is not currently supported. Check!"
        call PRC_abort
      end if
    end if

    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      elem => lcmesh%refElem2D

      allocate( fpos2D(this%Nnode_LCMeshAllFace(n),2) )

      if ( this%use_vmap_wide_flag ) then
        VmapB => this%VMapB2
      else
        VmapB => lcmesh%VMapB
      end if
      !$acc data create(fpos2D) copyin(VMapB)

      call MeshFieldCommBase_extract_bounddata2( lcmesh%pos_en(:,:,1), VmapB, size(VmapB), elem%Np*lcmesh%Ne, fpos2D(:,1) )
      call MeshFieldCommBase_extract_bounddata2( lcmesh%pos_en(:,:,2), VmapB, size(VmapB), elem%Np*lcmesh%Ne, fpos2D(:,2) )

      do f=1, this%nfaces_comm
        commdata => this%commdata_list(f,n)

        call push_localsendbuf( commdata%send_buf,                  & ! (inout)
          this%send_buf(:,:,n), commdata%s_faceID, this%is_f(f,n),  & ! (in)
          commdata%Nnode_LCMeshFace, this%bufsize_per_field,        & ! (in)
          this%field_num_tot, lcmesh, this%HaloSize_1D              ) ! (in)
        
        if ( commdata%s_panelID /= lcmesh%panelID &
            .and. ( this%hvfield_num > 0 .or. this%htensorfield_num > 0 ) ) then
          
          allocate( lcfpos2D(commdata%Nnode_LCMeshFace,2), unity_fac(commdata%Nnode_LCMeshFace) )
          unity_fac(:) = 1.0_RP
          !$acc enter data create(lcfpos2D) copyin(unity_fac) async(1)

          call push_localsendbuf( lcfpos2D,                            &
            fpos2D, commdata%s_faceID, this%is_f(f,n),                 &
            commdata%Nnode_LCMeshFace, this%Nnode_LCMeshAllFace(n), 2, &
            lcmesh, this%HaloSize_1D )
        end if

        if ( commdata%s_panelID /= lcmesh%panelID ) then
          if ( this%hvfield_num > 0 ) then
            allocate( tmp_svec2D(commdata%Nnode_LCMeshFace,2) )
            !$acc enter data create(tmp_svec2D) async(1)

            do varid=this%sfield_num+1, this%field_num_tot-1,2
              !$acc parallel loop present(tmp_svec2D, commdata%send_buf) async(1)
              do i=1, size(tmp_svec2D,1)
                tmp_svec2D(i,1) = commdata%send_buf(i,varid  )
                tmp_svec2D(i,2) = commdata%send_buf(i,varid+1)
              end do

              call CubedSphereCoordCnv_CS2LonLatVec( &
                lcmesh%panelID, lcfpos2D(:,1), lcfpos2D(:,2), unity_fac,  & ! (in)
                commdata%Nnode_LCMeshFace,                                & ! (in)
                tmp_svec2D(:,1), tmp_svec2D(:,2),                         & ! (in)
                commdata%send_buf(:,varid), commdata%send_buf(:,varid+1), & ! (out)
                gpu_async_id=1 ) ! (in)
            end do
            !$acc exit data delete(tmp_svec2D) async(1)
            deallocate( tmp_svec2D )
          end if
        end if
        

        if ( allocated(lcfpos2D) ) then
          !$acc exit data delete(lcfpos2D, unity_fac) async(1)
          deallocate( lcfpos2D, unity_fac )
        end if        
      end do
      !$acc end data
      deallocate( fpos2D )
    end do
    !$acc wait(1)

    !-----------------------

    call MeshFieldCommBase_exchange_core(this, this%commdata_list, do_wait )
    if ( .not. this%call_wait_flag_sub_get ) &
      call post_exchange_core( this )

    ! call PROF_rapend( 'comm_exchange', 2)

    return
  end subroutine MeshFieldCommCubedSphereDom2D_exchange

!- Private ---------------------------

  subroutine post_exchange_core( this )
    use scale_meshfieldcomm_base, only: LocalMeshCommData
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSVec
    implicit none

    class(MeshFieldCommCubedSphereDom2D), intent(inout), target :: this

    integer :: n, f
    integer :: varid

    real(RP), allocatable :: fpos2D(:,:)
    real(RP), allocatable :: lcfpos2D(:,:)
    real(RP), allocatable :: unity_fac(:)
    
    class(ElementBase2D), pointer :: elem
    type(LocalMesh2D), pointer :: lcmesh
    type(LocalMeshCommData), pointer :: commdata

    integer :: irs, ire    
    integer, pointer :: VmapB(:)
    !-----------------------------------------------------------------------------

    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh2d%lcmesh_list(n)
      elem => lcmesh%refElem2D

      allocate( fpos2D(this%Nnode_LCMeshAllFace(n),2) )
      if ( this%use_vmap_wide_flag ) then
        VmapB => this%VMapB2
      else
        VmapB => lcmesh%VMapB
      end if
      !$acc data create(fpos2D) copyin(VMapB)

      call MeshFieldCommBase_extract_bounddata2( lcmesh%pos_en(:,:,1), VMapB, size(VMapB), elem%Np*lcmesh%Ne, fpos2D(:,1) )
      call MeshFieldCommBase_extract_bounddata2( lcmesh%pos_en(:,:,2), VMapB, size(VMapB), elem%Np*lcmesh%Ne, fpos2D(:,2) )

      irs = 1
      do f=1, this%nfaces_comm
        commdata => this%commdata_list(f,n)
        ire = irs + commdata%Nnode_LCMeshFace - 1

        if ( commdata%s_panelID /= lcmesh%panelID &
            .and. ( this%hvfield_num > 0 .or. this%htensorfield_num > 0 ) ) then
          
          allocate( lcfpos2D(commdata%Nnode_LCMeshFace,2), unity_fac(commdata%Nnode_LCMeshFace) )
          unity_fac(:) = 1.0_RP
          !$acc enter data create(lcfpos2D) copyin(unity_fac) async(1)

          call push_localsendbuf( lcfpos2D,                            &
            fpos2D, commdata%s_faceID, this%is_f(f,n),                 &
            commdata%Nnode_LCMeshFace, this%Nnode_LCMeshAllFace(n), 2, &
            lcmesh, this%HaloSize_1D )
        end if


        if ( commdata%s_panelID /= lcmesh%panelID ) then
          if ( this%hvfield_num > 0 ) then
            do varid=this%sfield_num+1, this%field_num_tot-1, 2
              call CubedSphereCoordCnv_LonLat2CSVec( &
                lcmesh%panelID, lcfpos2D(:,1), lcfpos2D(:,2), unity_fac(:),          &
                commdata%Nnode_LCMeshFace,                                           &
                commdata%recv_buf(:,varid), commdata%recv_buf(:,varid+1),            &
                this%recv_buf(irs:ire,varid,n), this%recv_buf(irs:ire,varid+1,n)     )
            end do
          end if
        end if

        irs = ire + 1
        if ( allocated(lcfpos2D) ) then
          !$acc exit data delete(lcfpos2D, unity_fac) async(1)
          deallocate( lcfpos2D, unity_fac )
        end if
      end do

      !$acc end data
      deallocate( fpos2D )
    end do
    !$acc wait(1)

    return    
  end subroutine post_exchange_core

  !> Push temporary buffer of the local send data to the communication buffer
  subroutine push_localsendbuf( lc_send_buf, &
    send_buf, s_faceID, is, Nnode_LCMeshFace, bufsize_per_field, var_num, &
    lcmesh, haloSize_1D )
    implicit none

    integer, intent(in) :: var_num
    integer, intent(in) :: Nnode_LCMeshFace
    integer, intent(in) :: bufsize_per_field
    real(RP), intent(inout) :: lc_send_buf(Nnode_LCMeshFace,var_num)
    real(RP), intent(in) :: send_buf(bufsize_per_field,var_num)
    integer, intent(in) :: s_faceID, is
    type(LocalMesh2D), intent(in) :: lcmesh
    integer, intent(in) :: haloSize_1D

    integer :: i, v
    integer :: Ne_h1D
    integer :: Nfp
    !-----------------------------------------------------------------------------

    if ( s_faceID > 0 ) then
      !$omp parallel do
      !$acc parallel loop collapse(2) present(lc_send_buf, send_buf) async(1)
      do v=1, var_num
      do i=1, Nnode_LCMeshFace
        lc_send_buf(i,v) = send_buf(is+i-1,v)
      end do
      end do
    else if ( s_faceID < 0 ) then
      Nfp = lcmesh%refElem2D%Nfp
      Ne_h1D = Nnode_LCMeshFace / ( Nfp * haloSize_1D )
      call revert_hori( lc_send_buf, &
        send_buf, is, Ne_h1D, haloSize_1D, Nfp )
    end if

    return
  contains
    !> Revert the horizontal ordering of the local send buffer for negative s_faceID
    subroutine revert_hori( revert, ori, is_, Ne_h1D_, haloSize_1D_, Nfp_ )
      implicit none
      integer, intent(in) :: is_
      integer, intent(in) :: Ne_h1D_
      integer, intent(in) :: haloSize_1D_
      integer, intent(in) :: Nfp_

      real(RP), intent(out) :: revert(Nfp_, haloSize_1D_, Ne_h1D_, var_num)
      real(RP), intent(in) :: ori(bufsize_per_field,var_num)

      integer :: p, ph, i, v
      integer :: p_, i_
      integer :: id_src
      !---------------------------------------------------------------------------

      !$acc parallel loop collapse(4) present(revert, ori) async(1)
      do v=1, var_num
      do i=1, Ne_h1D_
      do ph=1, haloSize_1D_
      do p=1, Nfp_
        ! reverse orientation along the face
        i_ = Ne_h1D_ - i + 1
        p_ = Nfp_ - p + 1
        ! source ordering: (p, ph, i)
        id_src = is_ - 1 &
          + p_ &
          + (ph-1) * Nfp_ &
          + (i_-1) * Nfp_ * haloSize_1D_

        revert(p,ph,i,v) = ori(id_src,v)
      end do
      end do
      end do
      end do
      return
    end subroutine revert_hori
  end subroutine push_localsendbuf

  subroutine set_boundary_data2D_u1u2( buf_U, buf_V, &
    bufsize_per_field, elem, mesh, G_ij, VMapB,      &
    u1, u2)
  
    implicit none
    integer, intent(in) :: bufsize_per_field  
    type(ElementBase2D), intent(in) :: elem
    type(LocalMesh2D), intent(in) :: mesh
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
      = G_ij(VmapB,1,1) * buf_U(:) + G_ij(VmapB,1,2) * buf_V(:)
    u2(elem%Np*mesh%NeE+1:elem%Np*mesh%NeE+bufsize_per_field) &
      = G_ij(VmapB,2,1) * buf_U(:) + G_ij(VmapB,2,2) * buf_V(:)
#endif
    return
  end subroutine set_boundary_data2D_u1u2

end module scale_meshfieldcomm_cubedspheredom2d