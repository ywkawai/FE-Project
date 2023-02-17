!-------------------------------------------------------------------------------
!> module Mesh / Cubic 3D domain
!!
!! @par Description
!!      Mangage mesh data of cubic 3D domain for element-based methods
!!
!! @author Team SCALE
!<
#include "scaleFElib.h"
module scale_mesh_cubedom3d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_io
  use scale_precision
  use scale_prc

  use scale_mesh_base3d, only: &
    MeshBase3D, MeshBase3D_Init, MeshBase3D_Final, &
    MeshBase3D_setGeometricInfo

  use scale_mesh_base2d, only: &
    MeshBase2D, MeshBase2D_Init, MeshBase2D_Final, &
    MeshBase2D_setGeometricInfo

  use scale_localmesh_3d, only: &
    LocalMesh3D
  use scale_element_base, only: elementbase3D
  use scale_element_hexahedral, only: HexahedralElement

  use scale_mesh_rectdom2d, only: &
    MeshRectDom2D, MeshRectDom2D_setupLocalDom
  
  use scale_element_quadrilateral, only: QuadrilateralElement

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(MeshBase3D), public :: MeshCubeDom3D
    integer :: NeGX
    integer :: NeGY
    integer :: NeGZ

    integer :: NprcX
    integer :: NprcY
    integer :: NprcZ
    
    real(RP), public :: xmin_gl, xmax_gl
    real(RP), public :: ymin_gl, ymax_gl    
    real(RP), public :: zmin_gl, zmax_gl

    real(RP), allocatable :: FZ(:)

    integer, allocatable :: rcdomIJK2LCMeshID(:,:,:)

    logical :: isPeriodicX
    logical :: isPeriodicY
    logical :: isPeriodicZ

    type(MeshRectDom2D) :: mesh2D
    type(QuadrilateralElement) :: refElem2D
  contains
    procedure :: Init => MeshCubeDom3D_Init
    procedure :: Final => MeshCubeDom3D_Final
    procedure :: Generate => MeshCubeDom3D_generate
    procedure :: GetMesh2D => MeshCubeDom3D_getMesh2D
  end type MeshCubeDom3D

  public :: MeshCubeDom3D_coord_conv

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
  subroutine MeshCubeDom3D_Init( this,                          &
    NeGX, NeGY, NeGZ,                                           &
    dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
    isPeriodicX, isPeriodicY, isPeriodicZ,                      &
    refElem, NLocalMeshPerPrc, NprcX, NprcY,                    &
    nproc, myrank,                                              &
    FZ                                                          )
    
    implicit none

    class(MeshCubeDom3D), intent(inout) :: this
    integer, intent(in) :: NeGX
    integer, intent(in) :: NeGY
    integer, intent(in) :: NeGZ
    real(RP), intent(in) :: dom_xmin
    real(RP), intent(in) :: dom_xmax
    real(RP), intent(in) :: dom_ymin
    real(RP), intent(in) :: dom_ymax
    real(RP), intent(in) :: dom_Zmin
    real(RP), intent(in) :: dom_zmax
    logical, intent(in) :: isPeriodicX
    logical, intent(in) :: isPeriodicY
    logical, intent(in) :: isPeriodicZ
    type(HexahedralElement), intent(in), target :: refElem
    integer, intent(in) :: NLocalMeshPerPrc
    integer, intent(in) :: NprcX
    integer, intent(in) :: NprcY
    integer, intent(in), optional :: nproc
    integer, intent(in), optional :: myrank
    real(RP), intent(in), optional :: FZ(NeGZ+1)

    integer :: k
    real(RP) :: dz
    !-----------------------------------------------------------------------------
    
    this%NeGX = NeGX
    this%NeGY = NeGY
    this%NeGZ = NeGZ

    this%xmin_gl       = dom_xmin
    this%xmax_gl       = dom_xmax
    this%ymin_gl       = dom_ymin
    this%ymax_gl       = dom_ymax
    this%zmin_gl       = dom_zmin
    this%zmax_gl       = dom_zmax
    this%dom_vol = (this%xmax_gl - this%xmin_gl) * (this%ymax_gl - this%ymin_gl) * (this%zmax_gl - this%zmin_gl)  

    this%isPeriodicX = isPeriodicX
    this%isPeriodicY = isPeriodicY
    this%isPeriodicZ = isPeriodicZ

    this%NprcX = NprcX
    this%NprcY = NprcY
    this%NprcZ = 1

    !- Fz
    allocate( this%FZ(this%NeGZ+1) )
    if ( present(FZ) ) then
      this%FZ(:) = FZ(:)
    else
      this%FZ(1          ) = dom_Zmin
      this%FZ(this%NeGZ+1) = dom_Zmax
      dz = (dom_zmax - dom_zmin) / dble(this%NeGZ)
      do k=2, this%NeGZ
        this%FZ(k) = this%FZ(k-1) + dz
      end do
    end if

    !--
    call MeshBase3D_Init( this, refElem, NLocalMeshPerPrc, 6, &
                          nproc, myrank                       )

    !--- 2D mesh

    call this%refElem2D%Init( this%refElem3D%PolyOrder_h, refElem%IsLumpedMatrix() )
    
    this%mesh2D%isPeriodicX = isPeriodicX
    this%mesh2D%isPeriodicY = isPeriodicY
    this%mesh2D%NprcX = NprcX
    this%mesh2D%NprcY = NprcY

    call MeshBase2D_Init( this%mesh2D, this%refElem2D, NLocalMeshPerPrc,           &
                          nproc, myrank                                            )

    return
  end subroutine MeshCubeDom3D_Init

  subroutine MeshCubeDom3D_Final( this )
    use scale_prc
    implicit none

    class(MeshCubeDom3D), intent(inout) :: this

    integer :: n
    !-----------------------------------------------------------------------------
  
    if (this%isGenerated) then
      if ( allocated(this%rcdomIJK2LCMeshID) ) then
        deallocate( this%rcdomIJK2LCMeshID )
      end if
    else
      if ( allocated( this%FZ ) ) deallocate( this%FZ )
    end if

    call this%mesh2D%Final()
    call this%refElem2D%Final()

    call MeshBase3D_Final(this)

    return
  end subroutine MeshCubeDom3D_Final
  
  subroutine MeshCubeDom3D_getMesh2D( this, ptr_mesh2D )
    implicit none
    class(MeshCubeDom3D), intent(in), target :: this
    class(MeshBase2D), pointer, intent(out) :: ptr_mesh2D
    !-------------------------------------------------------

    ptr_mesh2D => this%mesh2D
    return
  end subroutine MeshCubeDom3D_getMesh2D

  subroutine MeshCubeDom3D_generate( this )
    implicit none

    class(MeshCubeDom3D), intent(inout), target :: this
            
    integer :: n
    integer :: p
    type(LocalMesh3D), pointer :: mesh

    integer :: tileID_table(this%LOCAL_MESH_NUM, this%PRC_NUM)
    integer :: panelID_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer :: pi_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer :: pj_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer :: pk_table(this%LOCAL_MESH_NUM*this%PRC_NUM)

    integer :: TILE_NUM_PER_PANEL
    real(RP) :: delx, dely, delz
    integer :: tileID
    !-----------------------------------------------------------------------------

    TILE_NUM_PER_PANEL = this%LOCAL_MESH_NUM_global / 1
    
    
    !--- Construct the connectivity of patches  (only master node)

    call MesshCubeDom3D_assignDomID( this,    & ! (in)
      tileID_table, panelID_table,            & ! (out)
      pi_table, pj_table, pk_table )            ! (out)
    
    !--- Setup local meshes managed by my process

    do n=1, this%LOCAL_MESH_NUM
      mesh => this%lcmesh_list(n)
      tileID = tileID_table(n, mesh%PRC_myrank+1)

      call MeshCubeDom3D_setupLocalDom( mesh, &
        tileID,  panelID_table(tileID),                                                           &
        pi_table(tileID), pj_table(tileID), pk_table(tileID), this%NprcX, this%NprcY, this%NprcZ, &
        this%xmin_gl, this%xmax_gl, this%ymin_gl, this%ymax_gl, this%zmin_gl, this%zmax_gl,       &
        this%NeGX/this%NprcX, this%NeGY/this%NprcY, this%NeGZ/this%NprcZ, this%FZ(:)              )

      call MeshRectDom2D_setupLocalDom( this%mesh2D%lcmesh_list(n), &
        tileID,  panelID_table(tileID),                             &
        pi_table(tileID), pj_table(tileID), this%NprcX, this%NprcY, &
        this%xmin_gl, this%xmax_gl, this%ymin_gl, this%ymax_gl,     &
        this%NeGX/this%NprcX, this%NeGY/this%NprcY )

      call mesh%SetLocalMesh2D( this%mesh2D%lcmesh_list(n) )
      !---
      ! write(*,*) "** my_rank=", mesh%PRC_myrank
      ! write(*,*) " tileID:", mesh%tileID
      ! write(*,*) " pnlID:", mesh%panelID, "-- i,j (within a panel)=", pi_table(tileID), pj_table(tileID)
      ! write(*,*) " local mesh:", n, "( total", this%LOCAL_MESH_NUM, ")"
      ! write(*,*) " panel_connect:", this%tilePanelID_globalMap(:,mesh%tileID)
      ! write(*,*) " tile_connect:", this%tileID_globalMap(:,mesh%tileID)
      ! write(*,*) " face_connect:", this%tileFaceID_globalMap(:,mesh%tileID)
      ! write(*,*) " domain size"
      ! write(*,*) "   NeX, NeY:", mesh%NeX, mesh%NeY
      ! write(*,*) "   [X], [Y]:",  mesh%xmin, mesh%xmax, ":", mesh%ymin, mesh%ymax
    end do

    ! To set rcdomIJP2LCMeshID, call AssignDomID for 2D mesh
    call this%mesh2D%AssignDomID( & 
      tileID_table, panelID_table,   & ! (out)
      pi_table, pj_table             ) ! (out)

    this%isGenerated = .true.
    this%mesh2D%isGenerated = .true.
    
    deallocate( this%FZ )

    return
  end subroutine MeshCubeDom3D_generate

  !- private ------------------------------------------------------

!OCL SERIAL
  subroutine MeshCubeDom3D_setupLocalDom( lcmesh, &
    tileID, panelID,                                            &
    i, j, k, NprcX, NprcY, NprcZ,                               &
    dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
    NeX, NeY, NeZ,                                              &
    FZ                                                          )

    use scale_prof
    use scale_meshutil_3d, only: &
      MeshUtil3D_genConnectivity,   &
      MeshUtil3D_genCubeDomain,     &
      MeshUtil3D_BuildInteriorMap,  &
      MeshUtil3D_genPatchBoundaryMap

    use scale_localmesh_base, only: BCTYPE_INTERIOR
    implicit none
      
    type(LocalMesh3D), intent(inout) :: lcmesh
    integer, intent(in) :: tileID
    integer, intent(in) :: panelID
    integer, intent(in) :: i, j, k
    integer, intent(in) :: NprcX, NprcY, NprcZ
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_ymin, dom_ymax
    real(RP), intent(in) :: dom_Zmin, dom_zmax
    integer, intent(in) :: NeX, NeY, NeZ
    real(RP), intent(in) :: FZ(NeZ*NprcZ+1)
    
    class(ElementBase3D), pointer :: elem
    real(RP) :: delx, dely
    real(RP) :: FZ_lc(NeZ+1)

    integer :: ii, jj, kk
    integer :: ke
    !-----------------------------------------------------------------------------

    elem => lcmesh%refElem3D
    lcmesh%tileID = tileID
    lcmesh%panelID = panelID
    
    !--
    lcmesh%Ne   = NeX * NeY * NeZ
    lcmesh%Nv  = (NeX + 1)*(NeY + 1)*(NeZ + 1)
    lcmesh%NeS = 1
    lcmesh%NeE = lcmesh%Ne
    lcmesh%NeA = lcmesh%Ne + 2*(NeX + NeY)*NeZ + 2*NeX*NeY

    lcmesh%NeX = NeX
    lcmesh%NeY = NeY
    lcmesh%NeZ = NeZ

    lcmesh%Ne2D  = NeX * NeY
    lcmesh%Ne2DA = NeX * NeY + 2*(NeX + NeY)

    !--
    delx = (dom_xmax - dom_xmin)/dble(NprcX)
    dely = (dom_ymax - dom_ymin)/dble(NprcY)
    FZ_lc(:) = Fz((k-1)*NeZ+1:k*NeZ+1)
    lcmesh%xmin = dom_xmin + (i-1)*delx
    lcmesh%xmax = dom_xmin +  i   *delx
    lcmesh%ymin = dom_ymin + (j-1)*dely
    lcmesh%ymax = dom_ymin +  j   *dely
    lcmesh%zmin = FZ_lc(1)
    lcmesh%zmax = FZ_lc(NeZ+1)
    
    !-
    allocate( lcmesh%pos_ev(lcmesh%Nv,3) )
    allocate( lcmesh%EToV(lcmesh%Ne,elem%Nv) )
    allocate( lcmesh%EToE(lcmesh%Ne,elem%Nfaces) )
    allocate( lcmesh%EToF(lcmesh%Ne,elem%Nfaces) )
    allocate( lcmesh%BCType(lcmesh%refElem%Nfaces,lcmesh%Ne) )
    allocate( lcmesh%VMapM(elem%NfpTot, lcmesh%Ne) )
    allocate( lcmesh%VMapP(elem%NfpTot, lcmesh%Ne) )
    allocate( lcmesh%MapM(elem%NfpTot, lcmesh%Ne) )
    allocate( lcmesh%MapP(elem%NfpTot, lcmesh%Ne) )
    
    allocate( lcmesh%EMap3Dto2D(lcmesh%Ne) )

    lcmesh%BCType(:,:) = BCTYPE_INTERIOR
    
    !----

    call MeshUtil3D_genCubeDomain( lcmesh%pos_ev, lcmesh%EToV,   & ! (out)
        lcmesh%NeX, lcmesh%xmin, lcmesh%xmax,                    & ! (in)
        lcmesh%NeY, lcmesh%ymin, lcmesh%ymax,                    & ! (in) 
        lcmesh%NeZ, lcmesh%zmin, lcmesh%zmax, FZ=FZ_lc           ) ! (in) 
    
    !---
    call MeshBase3D_setGeometricInfo( lcmesh, MeshCubeDom3D_coord_conv, MeshCubeDom3D_calc_normal )

    !---
    call MeshUtil3D_genConnectivity( lcmesh%EToE, lcmesh%EToF, & ! (out)
        lcmesh%EToV, lcmesh%Ne, elem%Nfaces )                    ! (in)

    !---
    call MeshUtil3D_BuildInteriorMap( lcmesh%VmapM, lcmesh%VMapP, lcmesh%MapM, lcmesh%MapP,           & ! (out)
      lcmesh%pos_en, lcmesh%pos_ev, lcmesh%EToE, lcmesh%EtoF, lcmesh%EtoV,                            & ! (in)
      elem%Fmask_h, elem%Fmask_v, lcmesh%Ne, lcmesh%Nv, elem%Np, elem%Nfp_h, elem%Nfp_v, elem%NfpTot, & ! (in)
      elem%Nfaces_h, elem%Nfaces_v, elem%Nfaces )                                                       ! (in)

    call MeshUtil3D_genPatchBoundaryMap( lcmesh%VMapB, lcmesh%MapB, lcmesh%VMapP,                       & !(out)
      lcmesh%pos_en, lcmesh%xmin, lcmesh%xmax, lcmesh%ymin, lcmesh%ymax, lcmesh%zmin, lcmesh%zmax,      & ! (in)
      elem%Fmask_h, elem%Fmask_v, lcmesh%Ne, lcmesh%Nv, elem%Np, elem%Nfp_h, elem%Nfp_v, elem%NfpTot,   & ! (in)
      elem%Nfaces_h, elem%Nfaces_v, elem%Nfaces )                                                         ! (in)
    
    !---
    !$omp parallel do collapse(2) private(ii,ke)
    do kk=1, lcmesh%NeZ
    do jj=1, lcmesh%NeY
    do ii=1, lcmesh%NeX
      ke = ii + (jj-1) * lcmesh%NeX + (kk-1) * lcmesh%NeX * lcmesh%NeY
      lcmesh%EMap3Dto2D(ke) = ii + (jj-1) * lcmesh%NeX
    end do
    end do
    end do

    return
  end subroutine MeshCubeDom3D_setupLocalDom

!OCL SERIAL
  subroutine MesshCubeDom3D_assignDomID( this, &
    tileID_table, panelID_table,               &
    pi_table, pj_table, pk_table )
  
    use scale_meshutil_3d, only: &       
      MeshUtil3D_buildGlobalMap    
    implicit none

    type(MeshCubeDom3D), target, intent(inout) :: this    
    integer, intent(out) :: tileID_table(this%LOCAL_MESH_NUM, this%PRC_NUM)
    integer, intent(out) :: panelID_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer, intent(out) :: pi_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer, intent(out) :: pj_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer, intent(out) :: pk_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    
    integer :: n
    integer :: p
    integer :: tileID
    integer :: is_lc, js_lc, ks_lc
    integer :: ilc_count, jlc_count, klc_count
    integer :: ilc, jlc, klc

    type(LocalMesh3D), pointer :: lcmesh
    !-----------------------------------------------------------------------------
    
    call MeshUtil3D_buildGlobalMap( &
      panelID_table, pi_table, pj_table, pk_table,                                          & ! (out)
      this%tileID_globalMap, this%tileFaceID_globalMap, this%tilePanelID_globalMap,         & ! (out)
      this%LOCAL_MESH_NUM_global, 6, 8,                                                     & ! (in)
      this%isPeriodicX, this%isPeriodicY, this%isPeriodicZ,                                 & ! (in)
      this%NprcX, this%NprcY, this%NprcZ )                                                    ! (in)

    !----

    do p=1, this%PRC_NUM
    do n=1, this%LOCAL_MESH_NUM
      tileID = n + (p-1)*this%LOCAL_MESH_NUM
      lcmesh => this%lcmesh_list(n)
      !-
      tileID_table(n,p)                   = tileID
      this%tileID_global2localMap(tileID) = n
      this%PRCRank_globalMap(tileID)      = p - 1

      !-
      if ( this%PRCRank_globalMap(tileID) == lcmesh%PRC_myrank ) then
        if (n==1) then
          is_lc = pi_table(tileID); ilc_count = 1
          js_lc = pj_table(tileID); jlc_count = 1
          ks_lc = pk_table(tileID); klc_count = 1
        end if
        if(is_lc < pi_table(tileID)) ilc_count = ilc_count + 1
        if(js_lc < pj_table(tileID)) jlc_count = jlc_count + 1
        if(ks_lc < pk_table(tileID)) klc_count = klc_count + 1
      end if 
    end do
    end do

    allocate( this%rcdomIJK2LCMeshID(ilc_count,jlc_count,klc_count) )
    do klc=1, klc_count
    do jlc=1, jlc_count
    do ilc=1, ilc_count
      this%rcdomIJK2LCMeshID(ilc,jlc,klc) = ilc + (jlc - 1)*ilc_count + (klc - 1)*ilc_count*jlc_count
    end do
    end do
    end do

    return
  end subroutine MesshCubeDom3D_assignDomID
  
!OCL SERIAL
  subroutine MeshCubeDom3D_coord_conv( x, y, z, xX, xY, xZ, yX, yY, yZ, zX, zY, zZ, &
    vx, vy, vz, elem )

    implicit none

    type(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: x(elem%Np), y(elem%Np), z(elem%Np)
    real(RP), intent(out) :: xX(elem%Np), xY(elem%Np), xZ(elem%Np)
    real(RP), intent(out) :: yX(elem%Np), yY(elem%Np), yZ(elem%Np)
    real(RP), intent(out) :: zX(elem%Np), zY(elem%Np), zZ(elem%Np)
    real(RP), intent(in) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)

    !-------------------------------------------------
    
    x(:) = vx(1) + 0.5_RP*(elem%x1(:) + 1.0_RP)*(vx(2) - vx(1))
    y(:) = vy(1) + 0.5_RP*(elem%x2(:) + 1.0_RP)*(vy(3) - vy(1))
    z(:) = vz(1) + 0.5_RP*(elem%x3(:) + 1.0_RP)*(vz(5) - vz(1))

    xX(:) = 0.5_RP*(vx(2) - vx(1)) !matmul(refElem%Dx1,mesh%x1(:,n))
    xY(:) = 0.0_RP                 !matmul(refElem%Dx2,mesh%x1(:,n))
    xZ(:) = 0.0_RP                 !matmul(refElem%Dx3,mesh%x1(:,n))
    yX(:) = 0.0_RP                 !matmul(refElem%Dx1,mesh%x2(:,n))
    yY(:) = 0.5_RP*(vy(3) - vy(1)) !matmul(refElem%Dx2,mesh%x2(:,n))
    yZ(:) = 0.0_RP                 !matmul(refElem%Dx3,mesh%x2(:,n))
    zX(:) = 0.0_RP                 !matmul(refElem%Dx1,mesh%x3(:,n))
    zY(:) = 0.0_RP                 !matmul(refElem%Dx2,mesh%x3(:,n))
    zZ(:) = 0.5_RP*(vz(5) - vz(1)) !matmul(refElem%Dx3,mesh%x3(:,n))

    return
  end subroutine MeshCubeDom3D_coord_conv  

!OCL SERIAL
  subroutine MeshCubeDom3D_calc_normal( normal_fn, &
    Escale_f, fid_h, fid_v, elem )

    implicit none

    type(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: normal_fn(elem%NfpTot,3)
    integer, intent(in) :: fid_h(elem%Nfp_h,elem%Nfaces_h)
    integer, intent(in) :: fid_v(elem%Nfp_v,elem%Nfaces_v)        
    real(RP), intent(in) :: Escale_f(elem%NfpTot,3,3)

    integer :: d
    !-------------------------------------------------

    do d=1, 3
      normal_fn(fid_h(:,1),d) = - Escale_f(fid_h(:,1),2,d)
      normal_fn(fid_h(:,2),d) = + Escale_f(fid_h(:,2),1,d)
      normal_fn(fid_h(:,3),d) = + Escale_f(fid_h(:,3),2,d)
      normal_fn(fid_h(:,4),d) = - Escale_f(fid_h(:,4),1,d)

      normal_fn(fid_v(:,1),d) = - Escale_f(fid_v(:,1),3,d)
      normal_fn(fid_v(:,2),d) = + Escale_f(fid_v(:,2),3,d)    
    end do

    return
  end subroutine MeshCubeDom3D_calc_normal 

end module scale_mesh_cubedom3d
