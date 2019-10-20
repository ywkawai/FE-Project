#include "scaleFElib.h"
module scale_mesh_base1d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_localmesh_1d, only: &
    LocalMesh1D, LocalMesh1D_Init, LocalMesh1D_Final
  
  use scale_mesh_base, only: MeshBase, &
    MeshBase_Init, MeshBase_Final, &
    MeshBase_setGeometricInfo
  
  use scale_element_base, only: elementbase1d
  use scale_element_line, only: LineElement

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, abstract, public, extends(MeshBase) :: Meshbase1d
    type(LocalMesh1D), allocatable :: lcmesh_list(:)
    class(elementbase1D), pointer :: refElem1D
    
    integer :: NeG
    real(RP), public :: xmin_gl, xmax_gl    
  contains
    procedure(Meshbase1d_generate), deferred :: Generate 
  end type Meshbase1d

  interface 
    subroutine Meshbase1d_generate(this)
      import Meshbase1d
      class(Meshbase1d), intent(inout), target :: this
    end subroutine Meshbase1d_generate
  end interface

  public :: Meshbase1d_Init, Meshbase1d_Final
  public :: Meshbase1d_setGeometricInfo, Meshbase1D_assignDomID, MeshBase1D_setupLocalDom

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
  subroutine Meshbase1d_Init( this,         &
    NeG,                                    &
    dom_xmin, dom_xmax,                     &
    refElem, NLocalMeshPerPrc )
    
    use scale_prc, only: PRC_myrank
    implicit none

    class(Meshbase1d), intent(inout) :: this
    integer, intent(in) :: NeG
    real(RP), intent(in) :: dom_xmin
    real(RP), intent(in) :: dom_xmax   
    class(elementbase1d), intent(in), target :: refElem
    integer, intent(in) :: NLocalMeshPerPrc

    integer :: n

    !-----------------------------------------------------------------------------
    
    this%NeG = NeG
    this%xmin_gl       = dom_xmin
    this%xmax_gl       = dom_xmax

    this%refElem1D => refElem
    call MeshBase_Init(this, refElem, NLocalMeshPerPrc, 2)
    
    allocate( this%lcmesh_list(this%LOCAL_MESH_NUM) )
    do n=1, this%LOCAL_MESH_NUM
      call LocalMesh1D_Init( this%lcmesh_list(n), refElem, PRC_myrank )
    end do

  end subroutine Meshbase1d_Init

  subroutine Meshbase1d_Final( this )
    
    class(Meshbase1d), intent(inout) :: this

    integer :: n

    !-----------------------------------------------------------------------------
  
    do n=1, this%LOCAL_MESH_NUM
      call LocalMesh1D_Final( this%lcmesh_list(n) )
    end do
    deallocate( this%lcmesh_list )
    
    call MeshBase_Final(this)
  end subroutine Meshbase1d_Final
  
  subroutine Meshbase1d_setGeometricInfo( mesh )

    implicit none
    
    type(LocalMesh1D), intent(inout) :: mesh

    class(Elementbase1D), pointer :: refElem
    integer :: n
    integer :: f
    integer :: i
    real(RP) :: vx(2)
    real(RP) :: xr(mesh%refElem%Np)
    real(DP) :: Escale(1,1,mesh%refElem%Np)
    integer :: fmask(mesh%refElem%NfpTot)
    integer :: fid(mesh%refElem1D%Nfp,mesh%refElem1D%Nfaces)

  !-----------------------------------------------------------------------------

    call MeshBase_setGeometricInfo(mesh, 1)
    refElem => mesh%refElem1D
    
    fmask(:) = reshape(refElem%Fmask, shape(fmask))
    do f=1, refElem%Nfaces
    do i=1, refElem%Nfp
       fid(i,f) = i + (f-1)*refElem%Nfp
    end do
    end do
  
    do n=1, mesh%Ne

       vx(:) = mesh%pos_ev(mesh%EToV(n,:),1)
       mesh%pos_en(:,n,1) = vx(1) + 0.5_RP*(refElem%x1(:) + 1.0_RP)*(vx(2) - vx(1))
       
       xr(:) = 0.5_RP*(vx(2) - vx(1)) !matmul(refElem%Dr,mesh%x(:,n))

       mesh%J(:,n) = xr
       mesh%Escale(:,n,1,1) =   1.0_RP/mesh%J(:,n)

       !* Face

       !
       !mesh%fx(:,n) = mesh%x(fmask(:),n)
       !mesh%fy(:,n) = mesh%y(fmask(:),n)

       ! Calculate normal vectors
       mesh%normal_fn(fid(:,1),n,1) = - 1.0_RP
       mesh%normal_fn(fid(:,2),n,1) = + 1.0_RP
       mesh%sJ(:,n) = 1.0_RP
       mesh%Fscale(:,n) = mesh%sJ(:,n)/mesh%J(fmask(:),n)
       mesh%Gsqrt(:,n) = 1.0_RP
    end do

  end subroutine Meshbase1d_setGeometricInfo

  subroutine Meshbase1D_assignDomID( this, &
    Nprc,                               &
    tileID_table, panelID_table,        &
    pi_table )
  
    use scale_meshutil_1d, only: &       
      MeshUtil1D_buildGlobalMap

    class(Meshbase1d), intent(inout) :: this    
    integer, intent(out) :: Nprc
    integer, intent(out) :: tileID_table(this%LOCAL_MESH_NUM, this%PRC_NUM)
    integer, intent(out) :: panelID_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer, intent(out) :: pi_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    
    integer :: n
    integer :: p
    integer :: tileID
    
    !-----------------------------------------------------------------------------
    
    Nprc = this%PRC_NUM

    call MeshUtil1D_buildGlobalMap( &
      panelID_table, pi_table,                                                      & ! (out)
      this%tileID_globalMap, this%tileFaceID_globalMap, this%tilePanelID_globalMap, & ! (out)
      this%LOCAL_MESH_NUM_global )                                                    ! (in)

    do p=1, this%PRC_NUM
    do n=1, this%LOCAL_MESH_NUM
      tileID = n + (p-1)*this%LOCAL_MESH_NUM

      tileID_table(n,p)                   = tileID
      this%tileID_global2localMap(tileID) = n
      this%PRCRank_globalMap(tileID)      = p - 1
    end do
    end do

  end subroutine Meshbase1D_assignDomID

  subroutine MeshBase1D_setupLocalDom( mesh,    &
    tileID, panelID,                            &
    i, Nprc,                                    &
    dom_xmin, dom_xmax,                         &
    Ne )

    use scale_meshutil_1d, only: &
      MeshUtil1D_genConnectivity,   &
      MeshUtil1D_genLineDomain,     &
      MeshUtil1D_BuildInteriorMap,  &
      MeshUtil1D_genPatchBoundaryMap

    use scale_localmesh_base, only: BCTYPE_INTERIOR
    implicit none
      
    type(LocalMesh1D), intent(inout) :: mesh
    integer, intent(in) :: tileID
    integer, intent(in) :: panelID
    integer, intent(in) :: i
    integer, intent(in) :: Nprc
    real(RP) :: dom_xmin, dom_xmax
    integer, intent(in) ::Ne
    
    class(ElementBase1D), pointer :: elem
    real(RP) :: delx
    
    !-----------------------------------------------------------------------------

    elem => mesh%refElem1D
    mesh%tileID = tileID
    mesh%panelID = panelID
    
    !--

    mesh%Ne  = Ne
    mesh%Nv  = Ne + 1
    mesh%NeS = 1
    mesh%NeE = mesh%Ne
    mesh%NeA = mesh%Ne + 2

    delx = (dom_xmax - dom_xmin)/dble(Nprc)

    mesh%xmin = dom_xmin + (i-1)*delx
    mesh%xmax = dom_xmin +  i   *delx

    allocate(mesh%pos_ev(mesh%Nv,1))
    allocate( mesh%EToV(mesh%Ne,2) )
    allocate( mesh%EToE(mesh%Ne,elem%Nfaces) )
    allocate( mesh%EToF(mesh%Ne,elem%Nfaces) )
    allocate( mesh%BCType(mesh%refElem%Nfaces,mesh%Ne) )
    allocate( mesh%VMapM(elem%NfpTot, mesh%Ne) )
    allocate( mesh%VMapP(elem%NfpTot, mesh%Ne) )
    allocate( mesh%MapM(elem%NfpTot, mesh%Ne) )
    allocate( mesh%MapP(elem%NfpTot, mesh%Ne) )

    mesh%BCType(:,:) = BCTYPE_INTERIOR

    !----

    call MeshUtil1D_genLineDomain( mesh%pos_ev, mesh%EToV,   & ! (out)
        & mesh%Ne, mesh%xmin, mesh%xmax )                      ! (in)


    !---
    call MeshBase1D_setGeometricInfo( mesh )
    
    !---

    call MeshUtil1D_genConnectivity( mesh%EToE, mesh%EToF, & ! (out)
        & mesh%EToV, mesh%Ne, elem%Nfaces )                  ! (in)

    !---
    call MeshUtil1D_BuildInteriorMap( mesh%VmapM, mesh%VMapP, mesh%MapM, mesh%MapP,  &
      & mesh%pos_en, mesh%pos_ev, mesh%EToE, mesh%EtoF, mesh%EtoV,                   &
      & elem%Fmask, mesh%Ne, elem%Np, elem%Nfp, elem%Nfaces, mesh%Nv )

    call MeshUtil1D_genPatchBoundaryMap( mesh%VMapB, mesh%MapB, mesh%VMapP, &
      & mesh%pos_en, mesh%xmin, mesh%xmax,                                  &
      & elem%Fmask, mesh%Ne, elem%Np, elem%Nfp, elem%Nfaces, mesh%Nv)
    
  end subroutine MeshBase1D_setupLocalDom

end module scale_mesh_base1d