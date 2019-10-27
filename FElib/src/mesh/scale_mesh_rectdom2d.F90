#include "scaleFElib.h"
module scale_mesh_rectdom2d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_mesh_base2d, only: &
    MeshBase2D, MeshBase2D_Init, MeshBase2D_Final, &
    MeshBase2D_setGeometricInfo

  use scale_localmesh_2d, only: &
    LocalMesh2D
  use scale_element_base, only: elementbase2D
  use scale_element_quadrilateral, only: QuadrilateralElement

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(MeshBase2D), public :: MeshRectDom2D
    integer :: NeGX
    integer :: NeGY
    
    real(RP), public :: xmin_gl, xmax_gl
    real(RP), public :: ymin_gl, ymax_gl    
    integer, allocatable :: rcdomIJ2LCMeshID(:,:)

    logical :: isPeriodicX
    logical :: isPeriodicY
  contains
    procedure :: Init => MeshRectDom2D_Init
    procedure :: Final => MeshRectDom2D_Final
    procedure :: Generate => MeshRectDom2D_generate
  end type MeshRectDom2D

  public :: MeshRectDom2D_coord_conv

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
  subroutine MeshRectDom2D_Init(this,       &
    NeGX, NeGY,                             &
    dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
    isPeriodicX, isPeriodicY,               &
    refElem, NLocalMeshPerPrc )
    
    implicit none

    class(MeshRectDom2D), intent(inout) :: this
    integer, intent(in) :: NeGX
    integer, intent(in) :: NeGY
    real(RP), intent(in) :: dom_xmin
    real(RP), intent(in) :: dom_xmax
    real(RP), intent(in) :: dom_ymin
    real(RP), intent(in) :: dom_ymax
    logical, intent(in) :: isPeriodicX
    logical, intent(in) :: isPeriodicY
    type(QuadrilateralElement), intent(in), target :: refElem
    integer, intent(in) :: NLocalMeshPerPrc

    !-----------------------------------------------------------------------------
    
    this%NeGX = NeGX
    this%NeGY = NeGY

    this%xmin_gl       = dom_xmin
    this%xmax_gl       = dom_xmax
    this%ymin_gl       = dom_ymin
    this%ymax_gl       = dom_ymax

    this%isPeriodicX = isPeriodicX
    this%isPeriodicY = isPeriodicY

    call MeshBase2D_Init(this, refElem, NLocalMeshPerPrc)
  end subroutine MeshRectDom2D_Init

  subroutine MeshRectDom2D_Final( this )
    
    class(MeshRectDom2D), intent(inout) :: this

    integer :: n

    !-----------------------------------------------------------------------------
  
    if (this%isGenerated) then
      deallocate( this%rcdomIJ2LCMeshID )
    end if

    call MeshBase2D_Final(this)

  end subroutine MeshRectDom2D_Final
  
  subroutine MeshRectDom2D_generate( this )
    
    implicit none

    class(MeshRectDom2D), intent(inout), target :: this
            
    integer :: n
    integer :: p
    type(LocalMesh2D), pointer :: mesh

    integer :: tileID_table(this%LOCAL_MESH_NUM, this%PRC_NUM)
    integer :: panelID_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer :: pi_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer :: pj_table(this%LOCAL_MESH_NUM*this%PRC_NUM)

    integer :: TILE_NUM_PER_PANEL
    integer :: NprcX, NprcY
    real(RP) :: delx, dely
    integer :: tileID
    
    !-----------------------------------------------------------------------------

    TILE_NUM_PER_PANEL = this%LOCAL_MESH_NUM_global / 1
    
    
    !--- Construct the connectivity of patches  (only master node)

    call MesshRectDom2D_assignDomID( this,    & ! (in)
        & NprcX, NprcY,                & ! (out)
        & tileID_table, panelID_table, & ! (out)
        & pi_table, pj_table )           ! (out)
  
    !--- Setup local meshes managed by my process

    do n=1, this%LOCAL_MESH_NUM
      mesh => this%lcmesh_list(n)
      tileID = tileID_table(n, mesh%PRC_myrank+1)

      call MeshRectDom2D_setupLocalDom( mesh, &
         & tileID,  panelID_table(tileID),                     &
         & pi_table(tileID), pj_table(tileID), NprcX, NprcY,       &
         & this%xmin_gl, this%xmax_gl, this%ymin_gl, this%ymax_gl, &
         & this%NeGX/NprcX, this%NeGY/NprcY )

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

    this%isGenerated = .true.
  end subroutine MeshRectDom2D_generate

  !- private ------------------------------------------------------

  subroutine MeshRectDom2D_setupLocalDom( lcmesh, &
    tileID, panelID,                         &
    i, j, NprcX, NprcY,                      &
    dom_xmin, dom_xmax, dom_ymin, dom_ymax,  &
    NeX, NeY )

    use scale_meshutil_2d, only: &
      MeshUtil2D_genConnectivity,   &
      MeshUtil2D_genRectDomain,     &
      MeshUtil2D_BuildInteriorMap,  &
      MeshUtil2D_genPatchBoundaryMap

    use scale_localmesh_base, only: BCTYPE_INTERIOR
    implicit none
      
    type(LocalMesh2D), intent(inout) :: lcmesh
    integer, intent(in) :: tileID
    integer, intent(in) :: panelID
    integer, intent(in) :: i, j
    integer, intent(in) :: NprcX, NprcY
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_ymin, dom_ymax
    integer, intent(in) ::NeX, NeY
    
    class(ElementBase2D), pointer :: elem
    real(RP) :: delx, dely
    
    !-----------------------------------------------------------------------------

    elem => lcmesh%refElem2D
    lcmesh%tileID = tileID
    lcmesh%panelID = panelID
    
    !--

    lcmesh%Ne  = NeX * NeY
    lcmesh%Nv  = (NeX + 1)*(NeY + 1)
    lcmesh%NeS = 1
    lcmesh%NeE = lcmesh%Ne
    lcmesh%NeA = lcmesh%Ne + 2*(NeX + NeY)

    lcmesh%NeX = NeX
    lcmesh%NeY = NeY

    delx = (dom_xmax - dom_xmin)/dble(NprcX)
    dely = (dom_ymax - dom_ymin)/dble(NprcY)

    lcmesh%xmin = dom_xmin + (i-1)*delx
    lcmesh%xmax = dom_xmin +  i   *delx
    lcmesh%ymin = dom_ymin + (j-1)*dely
    lcmesh%ymax = dom_ymin +  j   *dely

    allocate( lcmesh%pos_ev(lcmesh%Nv,2) )
    allocate( lcmesh%EToV(lcmesh%Ne,elem%Nv) )
    allocate( lcmesh%EToE(lcmesh%Ne,elem%Nfaces) )
    allocate( lcmesh%EToF(lcmesh%Ne,elem%Nfaces) )
    allocate( lcmesh%BCType(lcmesh%refElem%Nfaces,lcmesh%Ne) )
    allocate( lcmesh%VMapM(elem%NfpTot, lcmesh%Ne) )
    allocate( lcmesh%VMapP(elem%NfpTot, lcmesh%Ne) )
    allocate( lcmesh%MapM(elem%NfpTot, lcmesh%Ne) )
    allocate( lcmesh%MapP(elem%NfpTot, lcmesh%Ne) )

    lcmesh%BCType(:,:) = BCTYPE_INTERIOR

    !----

    call MeshUtil2D_genRectDomain( lcmesh%pos_ev, lcmesh%EToV,   & ! (out)
      lcmesh%NeX, lcmesh%xmin, lcmesh%xmax,                      & ! (in)
      lcmesh%NeY, lcmesh%ymin, lcmesh%ymax )                       ! (in)
    
    !---
    call MeshBase2D_setGeometricInfo(lcmesh, MeshRectDom2D_coord_conv, MeshRectDom2D_calc_normal )

    !---

    call MeshUtil2D_genConnectivity( lcmesh%EToE, lcmesh%EToF, & ! (out)
      lcmesh%EToV, lcmesh%Ne, elem%Nfaces )                      ! (in)

    !---
    call MeshUtil2D_BuildInteriorMap( lcmesh%VmapM, lcmesh%VMapP, lcmesh%MapM, lcmesh%MapP, &
      lcmesh%pos_en, lcmesh%pos_ev, lcmesh%EToE, lcmesh%EtoF, lcmesh%EtoV,                  &
      elem%Fmask, lcmesh%Ne, elem%Np, elem%Nfp, elem%Nfaces, lcmesh%Nv )

    call MeshUtil2D_genPatchBoundaryMap( lcmesh%VMapB, lcmesh%MapB, lcmesh%VMapP, &
      lcmesh%pos_en, lcmesh%xmin, lcmesh%xmax, lcmesh%ymin, lcmesh%ymax,          &
      elem%Fmask, lcmesh%Ne, elem%Np, elem%Nfp, elem%Nfaces, lcmesh%Nv)
    
    return
  end subroutine MeshRectDom2D_setupLocalDom

  subroutine MesshRectDom2D_assignDomID( this, &
    NprcX, NprcY,                       &
    tileID_table, panelID_table,        &
    pi_table, pj_table )
  
    use scale_prc, only: PRC_myrank
    use scale_meshutil_2d, only: &       
      MeshUtil2D_buildGlobalMap
    
    implicit none

    type(MeshRectDom2D), intent(inout) :: this    
    integer, intent(out) :: NprcX
    integer, intent(out) :: NprcY
    integer, intent(out) :: tileID_table(this%LOCAL_MESH_NUM, this%PRC_NUM)
    integer, intent(out) :: panelID_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer, intent(out) :: pi_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer, intent(out) :: pj_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    
    integer :: n
    integer :: p
    integer :: tileID
    integer :: is_lc, js_lc
    integer :: ilc_count, jlc_count
    integer :: ilc, jlc
    
    !-----------------------------------------------------------------------------

    call MeshUtil2D_buildGlobalMap( &
      panelID_table, pi_table, pj_table,                                            & ! (out)
      this%tileID_globalMap, this%tileFaceID_globalMap, this%tilePanelID_globalMap, & ! (out)
      this%LOCAL_MESH_NUM_global, this%isPeriodicX, this%isPeriodicY )                ! (in)                                        ! (in)

    NprcX = maxval(pi_table)
    NprcY = maxval(pj_table)
    
    !----
    
    do p=1, this%PRC_NUM
    do n=1, this%LOCAL_MESH_NUM
      tileID = n + (p-1)*this%LOCAL_MESH_NUM
      
      !-
      tileID_table(n,p)                   = tileID
      this%tileID_global2localMap(tileID) = n
      this%PRCRank_globalMap(tileID)      = p - 1

      !-
      if ( this%PRCRank_globalMap(tileID) == PRC_myrank ) then
        if (n==1) then
          is_lc = pi_table(tileID); ilc_count = 1
          js_lc = pj_table(tileID); jlc_count = 1
        end if
        if(ilc_count > pi_table(tileID)) ilc_count = ilc_count + 1
        if(jlc_count > pj_table(tileID)) jlc_count = jlc_count + 1
      end if 
    end do
    end do

    allocate( this%rcdomIJ2LCMeshID(ilc_count,jlc_count) )
    do jlc=1, jlc_count
    do ilc=1, ilc_count
      this%rcdomIJ2LCMeshID(ilc,jlc) = ilc + (jlc - 1)*ilc_count
    end do
    end do

    return
  end subroutine MesshRectDom2D_assignDomID
  
  subroutine MeshRectDom2D_coord_conv( x, y, xr, xs, yr, ys, &
    vx, vy, elem )

    implicit none

    type(elementbase2D), intent(in) :: elem
    real(RP), intent(out) :: x(elem%Np), y(elem%Np)
    real(RP), intent(out) :: xr(elem%Np), xs(elem%Np), yr(elem%Np), ys(elem%Np)
    real(RP), intent(in) :: vx(elem%Nv), vy(elem%Nv)

    !-------------------------------------------------
    
    x(:) = vx(1) + 0.5_RP*(elem%x1(:) + 1.0_RP)*(vx(2) - vx(1))
    y(:) = vy(1) + 0.5_RP*(elem%x2(:) + 1.0_RP)*(vy(3) - vy(1))

    xr(:) = 0.5_RP*(vx(2) - vx(1)) !matmul(refElem%Dx1,mesh%x1(:,n))
    xs(:) = 0.0_RP                 !matmul(refElem%Dx2,mesh%x1(:,n))
    yr(:) = 0.0_RP                 !matmul(refElem%Dx1,mesh%x2(:,n))
    ys(:) = 0.5_RP*(vy(3) - vy(1)) !matmul(refElem%Dx2,mesh%x2(:,n))

    return
  end subroutine MeshRectDom2D_coord_conv  

  subroutine MeshRectDom2D_calc_normal( normal_fn, &
    Escale_f, fid, elem )

    implicit none

    type(elementbase2D), intent(in) :: elem
    real(RP), intent(out) :: normal_fn(elem%NfpTot,2)
    integer, intent(in) :: fid(elem%Nfp,elem%Nfaces)
    real(RP), intent(in) :: Escale_f(elem%NfpTot,2,2)

    integer :: d
    !-------------------------------------------------

    do d=1, 2
      normal_fn(fid(:,1),d) = - Escale_f(fid(:,1),2,d)
      normal_fn(fid(:,2),d) = + Escale_f(fid(:,2),1,d)
      normal_fn(fid(:,3),d) = + Escale_f(fid(:,3),2,d)
      normal_fn(fid(:,4),d) = - Escale_f(fid(:,4),1,d)  
    end do

    return
  end subroutine MeshRectDom2D_calc_normal 

end module scale_mesh_rectdom2d
