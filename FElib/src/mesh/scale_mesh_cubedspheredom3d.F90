#include "scaleFElib.h"
module scale_mesh_cubedspheredom3d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_io
  use scale_precision

  use scale_mesh_base3d, only: &
    MeshBase3D, MeshBase3D_Init, MeshBase3D_Final,                               &
    MeshBase3D_setGeometricInfo,                                                 &
    MeshBase3D_DIMTYPE_NUM,                                                      &
    MeshBase3D_DIMTYPEID_X, MeshBase3D_DIMTYPEID_Y, MeshBase3D_DIMTYPEID_Z,      &
    MeshBase3D_DIMTYPEID_XYZ, MeshBase3D_DIMTYPEID_ZT, MeshBase3D_DIMTYPEID_XYZT
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_hexahedral, only: HexahedralElement

  use scale_mesh_base2d, only: &
    MeshBase2D, MeshBase2D_Init, MeshBase2D_Final, &
    MeshBase2D_setGeometricInfo
  use scale_mesh_cubedspheredom2d, only: &
    MeshCubedSphereDom2D, MeshCubedSphereDom2D_setupLocalDom
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_element_quadrilateral, only: QuadrilateralElement

  use scale_element_base, only: ElementBase2D, ElementBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(MeshBase3D), public :: MeshCubedSphereDom3D
    integer :: NeGX
    integer :: NeGY
    integer :: NeGZ
    
    real(RP), public :: xmin_gl, xmax_gl
    real(RP), public :: ymin_gl, ymax_gl
    real(RP), public :: zmin_gl, zmax_gl

    real(RP), allocatable :: FZ(:)

    integer, allocatable :: rcdomIJKP2LCMeshID(:,:,:,:)

    real(RP) :: RPlanet

    type(MeshCubedSphereDom2D) :: mesh2D
    type(QuadrilateralElement) :: refElem2D

    logical :: shallow_approx
  contains
    procedure :: Init => MeshCubedSphereDom3D_Init
    procedure :: Final => MeshCubedSphereDom3D_Final
    procedure :: Generate => MeshCubedSphereDom3D_generate
    procedure :: AssignDomID => MesshCubedSphereDom3D_assignDomID
    procedure :: GetMesh2D => MeshCubedSphereDom3D_getMesh2D
    procedure :: Set_geometric_with_vcoord => MeshCubeDom3D_set_geometric_with_vcoord
  end type MeshCubedSphereDom3D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  private :: MeshCubedSphereDom3D_calc_normal
  private :: MeshCubedSphereDom3D_coord_conv
  private :: MeshCubedSphereDom3D_set_metric
  private :: fill_halo_metric

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
!OCL SERIAL
  subroutine MeshCubedSphereDom3D_Init( this, &
    NeGX, NeGY, NeGZ, RPlanet,                &
    dom_zmin, dom_zmax,                       &
    refElem, NLocalMeshPerPrc,                &
    nproc, myrank,                            &
    FZ, shallow_approx                        )
    
    use scale_const, only: &
      PI => CONST_PI
    implicit none

    class(MeshCubedSphereDom3D), intent(inout) :: this
    integer, intent(in) :: NeGX
    integer, intent(in) :: NeGY
    integer, intent(in) :: NeGZ
    real(RP), intent(in) :: RPlanet
    real(RP), intent(in) :: dom_Zmin
    real(RP), intent(in) :: dom_zmax
    type(HexahedralElement), intent(in), target :: refElem
    integer, intent(in) :: NLocalMeshPerPrc
    integer, intent(in), optional :: nproc
    integer, intent(in), optional :: myrank
    real(RP), intent(in), optional :: FZ(NeGZ+1)
    logical, intent(in), optional :: shallow_approx

    integer :: k
    real(RP) :: dz
    !-----------------------------------------------------------------------------

    this%NeGX = NeGX
    this%NeGY = NeGY
    this%NeGZ = NeGZ

    this%xmin_gl = - 0.25_RP * PI 
    this%xmax_gl = + 0.25_RP * PI 
    this%ymin_gl = - 0.25_RP * PI 
    this%ymax_gl = + 0.25_RP * PI 
    this%zmin_gl = dom_zmin
    this%zmax_gl = dom_zmax
    this%RPlanet = RPlanet
    this%dom_vol = 4.0_RP / 3.0_RP * PI * ( ( dom_zmax + RPlanet )**3 - ( dom_zmin + RPlanet )**3 )


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

    !-
    if ( present(shallow_approx) ) then
      this%shallow_approx = shallow_approx
    else
      this%shallow_approx = .true.
    end if

    !--
    call MeshBase3D_Init( this, refElem, NLocalMeshPerPrc, 6, &
      nproc, myrank )

    !---
    call this%refElem2D%Init( this%refElem3D%PolyOrder_h, refElem%IsLumpedMatrix() )
    call this%mesh2D%Init( NeGX, NeGY, RPlanet, this%refElem2D, NLocalMeshPerPrc, &
      nproc, myrank )
      
    !-- Modify the information of dimension for the cubed sphere mesh
    call this%SetDimInfo( MeshBase3D_DIMTYPEID_X, "x", "1", "X-coordinate" )
    call this%SetDimInfo( MeshBase3D_DIMTYPEID_Y, "y", "1", "Y-coordinate" )
    call this%SetDimInfo( MeshBase3D_DIMTYPEID_Z, "z", "m", "Z-coordinate" )
    call this%SetDimInfo( MeshBase3D_DIMTYPEID_XYZ, "xyz", "1", "XYZ-coordinate" )
    call this%SetDimInfo( MeshBase3D_DIMTYPEID_ZT, "zt", "1", "XYZ-coordinate" )
    call this%SetDimInfo( MeshBase3D_DIMTYPEID_XYZT, "xyzt", "1", "XYZ-coordinate" )
  
    return
  end subroutine MeshCubedSphereDom3D_Init

!OCL SERIAL
  subroutine MeshCubedSphereDom3D_Final( this )
    implicit none
    class(MeshCubedSphereDom3D), intent(inout) :: this
    !-----------------------------------------------------------------------------
    
    if (this%isGenerated) then
      if ( allocated(this%rcdomIJKP2LCMeshID) ) then
        deallocate( this%rcdomIJKP2LCMeshID )
      end if
    else
      if ( allocated( this%FZ ) ) deallocate( this%FZ )
    end if

    call this%mesh2D%Final()
    call this%refElem2D%Final()

    call MeshBase3D_Final( this )

    return
  end subroutine MeshCubedSphereDom3D_Final

!OCL SERIAL
  subroutine MeshCubedSphereDom3D_getMesh2D( this, ptr_mesh2D )
    implicit none
    class(MeshCubedSphereDom3D), intent(in), target :: this
    class(MeshBase2D), pointer, intent(out) :: ptr_mesh2D
    !-------------------------------------------------------

    ptr_mesh2D => this%mesh2D
    return
  end subroutine MeshCubedSphereDom3D_getMesh2D

!OCL SERIAL
  subroutine MeshCubedSphereDom3D_generate( this )
    use scale_mesh_cubedspheredom2d, only: &
      MeshCubedSphereDom2D_check_division_params
    implicit none

    class(MeshCubedSphereDom3D), intent(inout), target :: this

    integer :: n
    type(LocalMesh3D), pointer :: mesh

    integer :: tileID_table(this%LOCAL_MESH_NUM, this%PRC_NUM)
    integer :: panelID_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer :: pi_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer :: pj_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer :: pk_table(this%LOCAL_MESH_NUM*this%PRC_NUM)

    integer :: NprcX_lc, NprcY_lc, NprcZ_lc
    integer :: tileID   
    !-----------------------------------------------------------------------------

    call MeshCubedSphereDom2D_check_division_params( &
      NprcX_lc, NprcY_lc,                            &
      this%PRC_NUM, this%LOCAL_MESH_NUM_global,      &
      .true. )

    NprcZ_lc = 1

    !--- Construct the connectivity of patches  (only master node)

    call this%AssignDomID( & 
      NprcX_lc, NprcY_lc, NprcZ_lc,  & ! (in)
      tileID_table, panelID_table,   & ! (out)
      pi_table, pj_table, pk_table   ) ! (out)

    do n=1, this%LOCAL_MESH_NUM
      mesh => this%lcmesh_list(n)
      tileID = tileID_table(n, mesh%PRC_myrank+1)

      call MeshCubedSphereDom3D_setupLocalDom( mesh, &
        tileID,  panelID_table(tileID),                             &
        pi_table(tileID), pj_table(tileID), pk_table(tileID),       &
        NprcX_lc, NprcY_lc, NprcZ_lc,                               &
        this%xmin_gl, this%xmax_gl, this%ymin_gl, this%ymax_gl,     &
        this%zmin_gl, this%zmax_gl, this%RPlanet,                   &
        this%NeGX/NprcX_lc, this%NeGY/NprcY_lc, this%NeGZ/NprcZ_lc, &
        this%FZ(:) )

      call MeshCubedSphereDom2D_setupLocalDom( this%mesh2D%lcmesh_list(n),   &
        tileID,  panelID_table(tileID),                                      &
        pi_table(tileID), pj_table(tileID), NprcX_lc, NprcY_lc,              &
        this%xmin_gl, this%xmax_gl, this%ymin_gl, this%ymax_gl,              &
        this%RPlanet, this%NeGX/NprcX_lc, this%NeGY/NprcY_lc                 )

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

    ! Set lon&lat position and metrics with the cubed sphere mesh
    call MeshCubedSphereDom3D_set_metric( this )

    ! To set rcdomIJP2LCMeshID, call AssignDomID for 2D mesh
    call this%mesh2D%AssignDomID( & 
      NprcX_lc, NprcY_lc,            & ! (in)
      tileID_table, panelID_table,   & ! (out)
      pi_table, pj_table             ) ! (out)

    this%isGenerated = .true.
    this%mesh2D%isGenerated = .true.
    
    deallocate( this%FZ )
    
    return
  end subroutine MeshCubedSphereDom3D_generate

!OCL SERIAL
  subroutine MeshCubeDom3D_set_geometric_with_vcoord(this, lcdomID, GsqrtV_lc, zlev_lc, G13_lc, G23_lc)
    implicit none
    class(MeshCubedSphereDom3D), intent(inout), target :: this
    integer, intent(in) :: lcdomID
    real(RP), intent(in) :: GsqrtV_lc(this%refElem3D%Np,this%lcmesh_list(lcdomID)%NeA)
    real(RP), intent(in) :: zlev_lc(this%refElem3D%Np,this%lcmesh_list(lcdomID)%NeA)
    real(RP), intent(in) :: G13_lc(this%refElem3D%Np,this%lcmesh_list(lcdomID)%NeA)
    real(RP), intent(in) :: G23_lc(this%refElem3D%Np,this%lcmesh_list(lcdomID)%NeA)

    integer :: ke
    class(LocalMesh3D), pointer :: lcmesh
    !-------------------------------------------------------

    lcmesh => this%lcmesh_list(lcdomID)

    !$omp parallel private(ke)
    !$omp do
    do ke=lcmesh%NeS, lcmesh%NeE
      lcmesh%zlev(:,ke) = zlev_lc(:,ke)
    end do
    !$omp do
    do ke=lcmesh%NeS, lcmesh%NeA   
      if ( this%shallow_approx ) then
        lcmesh%gam(:,ke) = 1.0_RP
      else
        lcmesh%gam(:,ke) = 1.0_RP + zlev_lc(:,ke) / this%RPlanet
      end if

      lcmesh%Gsqrt(:,ke) = GsqrtV_lc(:,ke) * lcmesh%gam(:,ke)**2 * lcmesh%Gsqrt(:,ke)
      lcmesh%GI3(:,ke,1) = G13_lc(:,ke)
      lcmesh%GI3(:,ke,2) = G23_lc(:,ke)
    end do
    !$omp end parallel

    return
  end subroutine MeshCubeDom3D_set_geometric_with_vcoord

  !- private ------------------------------

!OCL SERIAL
  subroutine MeshCubedSphereDom3D_setupLocalDom( lcmesh,        &
    tileID, panelID,                                            &
    i, j, k, NprcX, NprcY, NprcZ,                               &
    dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
    planet_radius, NeX, NeY, NeZ,                               &
    FZ )

    use scale_meshutil_cubedsphere3d, only: &
      MeshUtilCubedSphere3D_genConnectivity,   &
      MeshUtilCubedSphere3D_genCubeDomain,     &
      MeshUtilCubedSphere3D_BuildInteriorMap,  &
      MeshUtilCubedSphere3D_genPatchBoundaryMap
    
    use scale_localmesh_base, only: BCTYPE_INTERIOR

    implicit none
      
    type(LocalMesh3D), intent(inout) :: lcmesh
    integer, intent(in) :: tileID
    integer, intent(in) :: panelID
    integer, intent(in) :: i, j, k
    integer, intent(in) :: NprcX, NprcY, NprcZ
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_ymin, dom_ymax
    real(RP), intent(in) :: dom_zmin, dom_zmax
    real(RP), intent(in) :: planet_radius
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
    delx = ( dom_xmax - dom_xmin ) / dble(NprcX)
    dely = ( dom_ymax - dom_ymin ) / dble(NprcY)
    FZ_lc(:) = Fz((k-1)*NeZ+1:k*NeZ+1)
    lcmesh%xmin = dom_xmin + (i-1)*delx
    lcmesh%xmax = dom_xmin +  i   *delx
    lcmesh%ymin = dom_ymin + (j-1)*dely
    lcmesh%ymax = dom_ymin +  j   *dely
    lcmesh%zmin = FZ_lc(1)
    lcmesh%zmax = FZ_lc(NeZ+1)

    !--
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
    call MeshUtilCubedSphere3D_genCubeDomain( lcmesh%pos_ev, lcmesh%EToV, & ! (out)
      lcmesh%NeX, lcmesh%xmin, lcmesh%xmax,                               & ! (in)
      lcmesh%NeY, lcmesh%ymin, lcmesh%ymax,                               & ! (in) 
      lcmesh%NeZ, lcmesh%zmin, lcmesh%zmax, FZ=FZ_lc                      ) ! (in) 
    
    !---
    call MeshBase3D_setGeometricInfo(lcmesh, MeshCubedSphereDom3D_coord_conv, MeshCubedSphereDom3D_calc_normal )
    
    !---
    call MeshUtilCubedSphere3D_genConnectivity( lcmesh%EToE, lcmesh%EToF, & ! (out)
      lcmesh%EToV, lcmesh%Ne, elem%Nfaces )                                 ! (in)
    
    !---
    call MeshUtilCubedSphere3D_BuildInteriorMap( lcmesh%VmapM, lcmesh%VMapP, lcmesh%MapM, lcmesh%MapP,  & ! (out)
      lcmesh%pos_en, lcmesh%pos_ev, lcmesh%EToE, lcmesh%EtoF, lcmesh%EtoV,                              & ! (in)
      elem%Fmask_h, elem%Fmask_v, lcmesh%Ne, lcmesh%Nv, elem%Np, elem%Nfp_h, elem%Nfp_v, elem%NfpTot,   & ! (in)
      elem%Nfaces_h, elem%Nfaces_v, elem%Nfaces ) 

    
    call MeshUtilCubedSphere3D_genPatchBoundaryMap( lcmesh%VMapB, lcmesh%MapB, lcmesh%VMapP,            & !(out)
      lcmesh%pos_en, lcmesh%xmin, lcmesh%xmax, lcmesh%ymin, lcmesh%ymax, lcmesh%zmin, lcmesh%zmax,      & ! (in)
      elem%Fmask_h, elem%Fmask_v, lcmesh%Ne, lcmesh%Nv, elem%Np, elem%Nfp_h, elem%Nfp_v, elem%NfpTot,   & ! (in)
      elem%Nfaces_h, elem%Nfaces_v, elem%Nfaces )  
    
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
  end subroutine MeshCubedSphereDom3D_setupLocalDom

!OCL SERIAL
  subroutine MesshCubedSphereDom3D_assignDomID( this, &
    NprcX_lc, NprcY_lc, NprcZ_lc,                     &
    tileID_table, panelID_table,                      &
    pi_table, pj_table, pk_table )
  
    use scale_meshutil_cubedsphere3d, only: &       
      MeshUtilCubedSphere3D_buildGlobalMap
    
    implicit none

    class(MeshCubedSphereDom3D), target, intent(inout) :: this    
    integer, intent(in) :: NprcX_lc
    integer, intent(in) :: NprcY_lc
    integer, intent(in) :: NprcZ_lc
    integer, intent(out) :: tileID_table(this%LOCAL_MESH_NUM, this%PRC_NUM)
    integer, intent(out) :: panelID_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer, intent(out) :: pi_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer, intent(out) :: pj_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer, intent(out) :: pk_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    
    integer :: n
    integer :: prc
    integer :: tileID
    integer :: is_lc, js_lc, ks_lc, ps_lc
    integer :: ilc_count, jlc_count, klc_count, plc_count
    integer :: ilc, jlc, klc, plc
    
    type(LocalMesh3D), pointer :: lcmesh
    !-----------------------------------------------------------------------------

    call MeshUtilCubedSphere3D_buildGlobalMap( &
      panelID_table, pi_table, pj_table, pk_table,                                  & ! (out)
      this%tileID_globalMap, this%tileFaceID_globalMap, this%tilePanelID_globalMap, & ! (out)
      this%LOCAL_MESH_NUM_global, NprcZ_lc )                                          ! (in)

    !----
    
    do prc=1, this%PRC_NUM
    do n=1, this%LOCAL_MESH_NUM
      tileID = n + (prc-1)*this%LOCAL_MESH_NUM
      lcmesh => this%lcmesh_list(n)
      
      !-
      tileID_table(n,prc)                 = tileID
      this%tileID_global2localMap(tileID) = n
      this%PRCRank_globalMap(tileID)      = prc - 1

      !-
      if ( this%PRCRank_globalMap(tileID) == lcmesh%PRC_myrank ) then
        if (n==1) then
          is_lc = pi_table(tileID); ilc_count = 1
          js_lc = pj_table(tileID); jlc_count = 1
          ks_lc = pk_table(tileID); klc_count = 1
          ps_lc = panelID_table(tileID); plc_count = 1
        end if
        if(is_lc < pi_table(tileID)) ilc_count = ilc_count + 1
        if(js_lc < pj_table(tileID)) jlc_count = jlc_count + 1
        if(ks_lc < pk_table(tileID)) klc_count = klc_count + 1
        if(ps_lc < panelID_table(tileID)) plc_count = plc_count + 1
      end if 
    end do
    end do

    allocate( this%rcdomIJKP2LCMeshID(ilc_count,jlc_count,klc_count,plc_count) )
    do plc=1, plc_count
    do klc=1, klc_count
    do jlc=1, jlc_count
    do ilc=1, ilc_count
      this%rcdomIJKP2LCMeshID(ilc,jlc,klc,plc) = ilc + (jlc - 1)*ilc_count + (klc - 1)*ilc_count*jlc_count &
                                               + (plc-1)*ilc_count*jlc_count*klc_count
    end do
    end do
    end do
    end do

    return
  end subroutine MesshCubedSphereDom3D_assignDomID

!OCL SERIAL
  subroutine MeshCubedSphereDom3D_coord_conv( x, y, z, xX, xY, xZ, yX, yY, yZ, zX, zY, zZ, &
    vx, vy, vz, elem )

    implicit none

    type(ElementBase3D), intent(in) :: elem
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
  end subroutine MeshCubedSphereDom3D_coord_conv  

!OCL SERIAL
  subroutine MeshCubedSphereDom3D_calc_normal( normal_fn, &
    Escale_f, fid_h, fid_v, elem )

    implicit none

    type(ElementBase3D), intent(in) :: elem
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
  end subroutine MeshCubedSphereDom3D_calc_normal 

  !--

!OCL SERIAL
  subroutine MeshCubedSphereDom3D_set_metric( this )
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatPos, &
      CubedSphereCoordCnv_GetMetric
        
    implicit none
    class(MeshCubedSphereDom3D), intent(inout), target :: this

    integer :: n
    integer :: ke, ke2D

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase2D), pointer :: elem2D

    real(RP), allocatable :: gam2D(:,:)
    !----------------------------------------------------

    do n=1, this%mesh2D%LOCAL_MESH_NUM
      lcmesh => this%lcmesh_list(n)
      lcmesh2D => this%mesh2D%lcmesh_list(n)
      elem2D => lcmesh2D%refElem2D

      allocate( gam2D(elem2D%Np,lcmesh2D%Ne) )
      gam2D(:,:) = 1.0_RP

      call CubedSphereCoordCnv_CS2LonLatPos( &
        lcmesh2D%panelID, lcmesh2D%pos_en(:,:,1), lcmesh2D%pos_en(:,:,2), gam2D(:,:), & ! (in)
        lcmesh2D%Ne * elem2D%Np,                                                      & ! (in)
        lcmesh%lon2D(:,:), lcmesh%lat2D(:,:)                                          ) ! (out)

      call CubedSphereCoordCnv_GetMetric( &
        lcmesh2D%pos_en(:,:,1), lcmesh2D%pos_en(:,:,2), elem2D%Np * lcmesh2D%Ne, this%RPlanet, & ! (in)
        lcmesh%G_ij, lcmesh%GIJ, lcmesh%GsqrtH                                                 ) ! (out)

      !$omp parallel do private(ke2D)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)

        if ( this%shallow_approx ) then
          lcmesh%gam(:,ke) = 1.0_RP
        else
          lcmesh%gam(:,ke) = 1.0_RP + lcmesh%pos_en(:,ke,3) / this%RPlanet
        end if
        lcmesh%Gsqrt(:,ke) = lcmesh%GsqrtH(lcmesh%refElem3D%IndexH2Dto3D(:),ke2D)        
      end do

      call fill_halo_metric( lcmesh%Gsqrt, lcmesh%gam, lcmesh%VMapM, lcmesh%VMapP, lcmesh, lcmesh%refElem3D )

      deallocate( gam2D )
      !--
    end do


    return
  end subroutine MeshCubedSphereDom3D_set_metric

!OCL SERIAL
  subroutine fill_halo_metric( Gsqrt, gam, vmapM, vmapP, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    real(RP), intent(inout) :: Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: gam(elem%Np*lmesh%NeA)

    integer :: i, iM, iP
    !------------------------------------------------

    !$omp parallel do private(i, iM, iP)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)
      if ( iP > elem%Np * lmesh%Ne ) then
        Gsqrt(iP) = Gsqrt(iM)
        gam(iP) = gam(iM)
      end if
    end do  
    return
  end subroutine fill_halo_metric

end module scale_mesh_cubedspheredom3d