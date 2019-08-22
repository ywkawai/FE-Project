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
  use scale_element_quadrial, only: QuadrialElement

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

  public :: MeshRectDom2D_Init, MeshRectDom2D_Final
  public :: MeshRectDom2D_generate
  public :: MeshRectDom2D_setupLocalDom

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
    type(QuadrialElement), intent(in), target :: refElem
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

    call Messh2D_assignDomID( this,    & ! (in)
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

  subroutine MeshRectDom2D_setupLocalDom( mesh, &
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
      
    type(LocalMesh2D), intent(inout) :: mesh
    integer, intent(in) :: tileID
    integer, intent(in) :: panelID
    integer, intent(in) :: i, j
    integer, intent(in) :: NprcX, NprcY
    real(RP) :: dom_xmin, dom_xmax
    real(RP) :: dom_ymin, dom_ymax
    integer, intent(in) ::NeX, NeY
    
    class(ElementBase2D), pointer :: elem
    real(RP) :: delx, dely
    
    !-----------------------------------------------------------------------------

    elem => mesh%refElem2D
    mesh%tileID = tileID
    mesh%panelID = panelID
    
    !--

    mesh%Ne  = NeX * NeY
    mesh%Nv  = (NeX + 1)*(NeY + 1)
    mesh%NeS = 1
    mesh%NeE = mesh%Ne
    mesh%NeA = mesh%Ne + 2*(NeX + NeY)

    mesh%NeX = NeX
    mesh%NeY = NeY

    delx = (dom_xmax - dom_xmin)/dble(NprcX)
    dely = (dom_ymax - dom_ymin)/dble(NprcY)

    mesh%xmin = dom_xmin + (i-1)*delx
    mesh%xmax = dom_xmin +  i   *delx
    mesh%ymin = dom_ymin + (j-1)*dely
    mesh%ymax = dom_ymin +  j   *dely

    allocate(mesh%pos_ev(mesh%Nv,2))
    allocate( mesh%EToV(mesh%Ne,4) )
    allocate( mesh%EToE(mesh%Ne,elem%Nfaces) )
    allocate( mesh%EToF(mesh%Ne,elem%Nfaces) )
    allocate( mesh%BCType(mesh%refElem%Nfaces,mesh%Ne) )
    allocate( mesh%VMapM(elem%NfpTot, mesh%Ne) )
    allocate( mesh%VMapP(elem%NfpTot, mesh%Ne) )
    allocate( mesh%MapM(elem%NfpTot, mesh%Ne) )
    allocate( mesh%MapP(elem%NfpTot, mesh%Ne) )

    mesh%BCType(:,:) = BCTYPE_INTERIOR

    !----

    call MeshUtil2D_genRectDomain( mesh%pos_ev, mesh%EToV,   & ! (out)
        & mesh%NeX, mesh%xmin, mesh%xmax,                    & ! (in)
        & mesh%NeY, mesh%ymin, mesh%ymax )                     ! (in)

    !---
    call MeshBase2D_setGeometricInfo(mesh)
 
    !---

    call MeshUtil2D_genConnectivity( mesh%EToE, mesh%EToF, & ! (out)
        & mesh%EToV, mesh%Ne, elem%Nfaces )                  ! (in)

    !---
    call MeshUtil2D_BuildInteriorMap( mesh%VmapM, mesh%VMapP, mesh%MapM, mesh%MapP,  &
      & mesh%pos_en, mesh%pos_ev, mesh%EToE, mesh%EtoF, mesh%EtoV,                   &
      & elem%Fmask, mesh%Ne, elem%Np, elem%Nfp, elem%Nfaces, mesh%Nv )

    call MeshUtil2D_genPatchBoundaryMap( mesh%VMapB, mesh%MapB, mesh%VMapP, &
      & mesh%pos_en, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax,            &
      & elem%Fmask, mesh%Ne, elem%Np, elem%Nfp, elem%Nfaces, mesh%Nv)
    
    return
  end subroutine MeshRectDom2D_setupLocalDom

  subroutine Messh2D_assignDomID( this, &
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
    
    NprcX = int(sqrt(dble(this%PRC_NUM)))
    NprcY = this%PRC_NUM/NprcX

    call MeshUtil2D_buildGlobalMap( &
      panelID_table, pi_table, pj_table,                                            & ! (out)
      this%tileID_globalMap, this%tileFaceID_globalMap, this%tilePanelID_globalMap, & ! (out)
      this%LOCAL_MESH_NUM_global, this%isPeriodicX, this%isPeriodicY )                ! (in)                                        ! (in)

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
  end subroutine Messh2D_assignDomID
  
end module scale_mesh_rectdom2d
