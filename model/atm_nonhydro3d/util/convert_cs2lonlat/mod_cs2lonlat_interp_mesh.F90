!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_cs2lonlat_interp_mesh
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort
  use scale_const, only: &
    PI => CONST_PI,         &
    RPlanet => CONST_radius
  
  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom2d, only: &
    MeshCubedSphereDom2D,                       &
    MeshCubedSphereDom2D_check_division_params
  use scale_mesh_cubedspheredom3d, only: &
    MeshCubedSphereDom3D
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  public :: interp_mesh_Init
  public :: interp_mesh_Final

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, parameter :: ATMOS_MESH_NLocalMeshPerPrc = 1
  type(MeshRectDom2D), public, target :: out_mesh2D
  type(MeshCubeDom3D), public, target :: out_mesh3D
  type(MeshCubedSphereDom2D), public :: in_csmesh2D
  type(MeshCubedSphereDom3D), public :: in_csmesh3D

  type, public :: NodeMappingInfo
    type(LocalMesh2D), pointer :: lcmesh2D => null()
    type(LocalMesh3D), pointer :: lcmesh3D => null()
    integer, allocatable :: local_domID(:,:)
    integer, allocatable :: lcprc(:,:)
    integer, allocatable :: inCSPanelID(:,:)
    integer, allocatable :: elem_i(:,:)
    integer, allocatable :: elem_j(:,:)
    integer, allocatable :: elem_k(:,:)
    real(RP), allocatable :: elem_x(:,:)
    real(RP), allocatable :: elem_y(:,:)
    real(RP), allocatable :: elem_z(:,:)

    type(MeshCubedSphereDom2D), allocatable :: in_mesh2D_list(:)
    type(MeshCubedSphereDom3D), allocatable :: in_mesh3D_list(:)
    integer, allocatable :: in_tileID_list(:)
  contains
    procedure :: Init_2D => NodeMappingInfo_Init_2D
    procedure :: Init_3D => NodeMappingInfo_Init_3D
    generic :: Init => Init_2D, Init_3D
    procedure :: Final => NodeMappingInfo_Final
  end type
  type(NodeMappingInfo), public :: nodeMap_list(ATMOS_MESH_NLocalMeshPerPrc)

  integer, public :: in_Nprc             = 1
  integer, public :: in_NeGX             = 1
  integer, public :: in_NeGY             = 1
  integer, public :: in_NeGZ             = 1
  integer, public :: in_NLocalMeshPerPrc = 6  
  type(QuadrilateralElement), public :: in_elem2D
  type(HexahedralElement), public :: in_elem3D
  
  logical, public :: is_mesh3D

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private :: out_NprcX        = 1        ! x length of 2D processor topology (output)
  integer, private :: out_NprcY        = 1        ! y length of 2D processor topology (output)
  type(QuadrilateralElement), private :: out_elem2D
  type(HexahedralElement), private :: out_elem3D

  real(RP), private :: dom_zmin        = 0.0_RP
  real(RP), private :: dom_zmax        = 0.0_RP

contains
  subroutine interp_mesh_Init()
    use scale_mesh_base2d, only: &
      MeshBase2D_DIMTYPEID_X, MeshBase2D_DIMTYPEID_Y,   &
      MeshBase2D_DIMTYPEID_XY, MeshBase2D_DIMTYPEID_XYT
    use scale_mesh_base3d, only: &
      MeshBase3D_DIMTYPEID_X, MeshBase3D_DIMTYPEID_Y, MeshBase3D_DIMTYPEID_Z,  &
      MeshBase3D_DIMTYPEID_XYZ, MeshBase3D_DIMTYPEID_XYZT
    use scale_const, only: &
      UNDEF => CONST_UNDEF

    implicit none

  
    integer :: out_NeX          =  1
    integer :: out_NeY          =  1
    integer :: out_NeZ          = -1
    integer :: in_PolyOrder_h   = 1
    integer :: in_PolyOrder_v   = 1
    integer :: out_PolyOrder_h  = 1
    integer :: out_PolyOrder_v  = 1
    real(RP) :: out_dom_zmin
    real(RP) :: out_dom_zmax

    integer, parameter :: FZ_nmax = 1000
    logical :: is_spec_in_FZ
    logical :: is_spec_out_FZ
    real(RP) :: in_FZ(FZ_nmax)
    real(RP) :: out_FZ(FZ_nmax)

    namelist / PARAM_INTERP_MESH / &
      in_Nprc,             &
      in_NeGX,             &
      in_NeGY,             &
      in_NeGZ,             &
      in_NLocalMeshPerPrc, &
      in_PolyOrder_h,      &
      in_PolyOrder_v,      &
      out_NprcX,           &
      out_NeX,             &
      out_NprcY,           &
      out_NeY,             &
      out_NeZ,             &
      out_PolyOrder_h,     &
      out_PolyOrder_v,     &
      dom_zmin, dom_zmax,  &
      out_dom_zmin,        &
      out_dom_zmax,        &
      in_Fz, out_Fz

    
    integer :: ierr
    integer :: k
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("interp_mesh",*) 'Setup'
  
    out_dom_zmin = UNDEF
    out_dom_zmax = UNDEF

    in_FZ(:) = UNDEF
    out_FZ(:) = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_INTERP_MESH,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
        LOG_INFO("interp_mesh",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("interp_mesh",*) 'Not appropriate names in namelist PARAM_INTERP_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_INTERP_MESH)

    if ( out_dom_zmin == UNDEF ) out_dom_zmin = dom_zmin
    if ( out_dom_zmax == UNDEF ) out_dom_zmax = dom_zmax

    is_spec_in_FZ = .true.    
    do k=1, in_NeGZ+1
      if ( in_FZ(k) == UNDEF ) then
        is_spec_in_FZ = .false.; exit
      end if
    end do
    is_spec_out_FZ = .true.    
    do k=1, out_NeZ+1
      if ( out_FZ(k) == UNDEF ) then
        is_spec_out_FZ = .false.; exit
      end if
    end do

    !-    

    call in_elem2D%Init( in_PolyOrder_h, .false. )

    if (out_NeZ < 0 ) then
      is_mesh3D = .false.

      call out_elem2D%Init( out_PolyOrder_h, .false. )
  
      call out_mesh2D%Init( out_NprcX * out_NeX, out_NprcY * out_NeY, &
        0.0_RP, 360.0_RP, -90.0_RP, 90.0_RP,                          &
        .true., .false., out_elem2D, 1,                               &
        out_NprcX, out_NprcY          )

      call out_mesh2D%SetDimInfo( MeshBase2D_DIMTYPEID_X, 'lon', 'degree_east', 'longitude' )
      call out_mesh2D%SetDimInfo( MeshBase2D_DIMTYPEID_Y, 'lat', 'degree_north', 'latitude' )
      call out_mesh2D%SetDimInfo( MeshBase2D_DIMTYPEID_XY, 'lonlat', 'degree', 'longitude,latitude' )
      call out_mesh2D%SetDimInfo( MeshBase2D_DIMTYPEID_XYT, 'lonlatt', 'degree', 'longitude,latitude' )
  
      call out_mesh2D%Generate()  
      call construct_map()
    else
      is_mesh3D = .true.
      call in_elem3D%Init( in_PolyOrder_h, in_PolyOrder_v, .false. )
      call out_elem3D%Init( out_PolyOrder_h, out_PolyOrder_v, .false. )
      
      if ( is_spec_out_FZ ) then
        call out_mesh3D%Init( out_NprcX * out_NeX, out_NprcY * out_NeY, out_NeZ,          &
          0.0_RP, 360.0_RP, -90.0_RP, 90.0_RP, out_dom_zmin, out_dom_zmax,                &
          .true., .false., .false., out_elem3D, ATMOS_MESH_NLocalMeshPerPrc,              &
          out_NprcX, out_NprcY,                                                           &
          FZ=out_FZ(1:out_NeZ+1) )
      else
        call out_mesh3D%Init( out_NprcX * out_NeX, out_NprcY * out_NeY, out_NeZ,          &
          0.0_RP, 360.0_RP, -90.0_RP, 90.0_RP, out_dom_zmin, out_dom_zmax,                &
          .true., .false., .false., out_elem3D, ATMOS_MESH_NLocalMeshPerPrc,              &
          out_NprcX, out_NprcY                                                            )
      end if
      
      call out_mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_X, 'lon', 'degree_east', 'longitude' )
      call out_mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_Y, 'lat', 'degree_north', 'latitude' )
      call out_mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_Z, 'z', 'm', 'altitude' )
      call out_mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_XYZ, 'lonlatz', 'degree', 'longitude,latitude,altitude' )
      call out_mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_XYZT, 'lonlatzt', 'degree', 'longitude,latitude,altitude' )

      call out_mesh3D%Generate()
      call construct_map( in_FZ(1:in_NeGZ+1), is_spec_in_FZ )
    end if

    return
  end subroutine interp_mesh_Init

  subroutine interp_mesh_Final()
    implicit none

    integer :: n
    !-------------------------------------------

    call in_elem2D%Final()

    if ( is_mesh3D ) then
      call in_elem3D%Final()

      call out_mesh3D%Final()
      call out_elem3D%Final()  
    else
      call out_mesh2D%Final()
      call out_elem2D%Final()  
    end if

    do n=1, ATMOS_MESH_NLocalMeshPerPrc
      call nodeMap_list(n)%Final()
    end do
    
    return
  end subroutine interp_mesh_Final

  !- private -------------------------------------

!OCL SERIAL  
  subroutine construct_map( in_Fz, is_spec_in_Fz )
    implicit none

    real(RP), intent(in), optional :: in_Fz(1:in_NeGZ+1)
    logical, intent(in), optional :: is_spec_in_Fz

    integer :: n
    integer :: i, j
    real(RP) :: in_tiles_x(4,in_NLocalMeshPerPrc,in_Nprc)
    real(RP) :: in_tiles_y(4,in_NLocalMeshPerPrc,in_Nprc)
    integer :: in_panel(4,in_NLocalMeshPerPrc,in_Nprc)
    real(RP) :: delx, dely
    type(LocalMesh2D), pointer :: lcmesh2D
    type(LocalMesh3D), pointer :: lcmesh3D

    integer :: NprcX_lc, NprcY_lc
    integer :: tileID_table(in_NLocalMeshPerPrc,in_Nprc)
    integer :: panelID_table(in_NLocalMeshPerPrc*in_Nprc)
    integer :: pi_table(in_NLocalMeshPerPrc*in_Nprc)
    integer :: pj_table(in_NLocalMeshPerPrc*in_Nprc)

    integer :: in_p, in_n
    integer :: tileID
    type(MeshCubedSphereDom2D) :: in_csmesh2D_dummy
    !-------------------------------------------

    !---
    call in_csmesh2D_dummy%Init( in_NeGX, in_NeGY, RPlanet,   &
      in_elem2D, in_NLocalMeshPerPrc, nproc=in_Nprc, myrank=0 )
    
    call MeshCubedSphereDom2D_check_division_params(  &
      NprcX_lc, NprcY_lc,                             & ! (out)
      in_Nprc, in_NLocalMeshPerPrc * in_Nprc          ) ! (in)
    
    call in_csmesh2D_dummy%AssignDomID( NprcX_lc, NprcY_lc, & ! (in)
      tileID_table, panelID_table, pi_table, pj_table       ) ! (out)
    
    call in_csmesh2D_dummy%Final()
    !---

    delx = 0.5_RP * PI / dble(NprcX_lc)
    dely = 0.5_RP * PI / dble(NprcY_lc)

    do in_p=1, in_Nprc
    do in_n=1, in_NLocalMeshPerPrc
      tileID = tileID_table(in_n,in_p)
      i = pi_table(tileID); j = pj_table(tileID)

      in_tiles_x(:,in_n,in_p) = - 0.25_RP * PI + delx * dble( (/ i-1, i, i, i-1 /) )
      in_tiles_y(:,in_n,in_p) = - 0.25_RP * PI + dely * dble( (/ j-1, j-1, j, j /) )

      if (i==NprcX_lc) then
        in_tiles_x(2:3,in_n,in_p) = in_tiles_x(2:3,in_n,in_p) + 1.0E-12_RP * delx       
      end if
      if (j==NprcY_lc) then
        in_tiles_y(3:4,in_n,in_p) = in_tiles_y(3:4,in_n,in_p) + 1.0E-12_RP * dely
      end if
    end do
    end do

    if ( is_mesh3D ) then
      do n=1, out_mesh3D%LOCAL_MESH_NUM
        lcmesh3D => out_mesh3D%lcmesh_list(n)
        call nodeMap_list(n)%Init_3D( lcmesh3D, lcmesh3D%refElem3D, in_tiles_x, in_tiles_y, in_Fz, &
          tileID_table, panelID_table, NprcX_lc, NprcY_lc, is_spec_in_Fz )
      end do
    else
      do n=1, out_mesh2D%LOCAL_MESH_NUM
        lcmesh2D => out_mesh2D%lcmesh_list(n)
        call nodeMap_list(n)%Init_2D( lcmesh2D, lcmesh2D%refElem2D, in_tiles_x, in_tiles_y, &
          tileID_table, panelID_table, NprcX_lc, NprcY_lc )
      end do
    end if

    return
  end subroutine construct_map

!OCL SERIAL
  subroutine NodeMappingInfo_Init_2D( this, &
      lcmesh, elem2D, tile_x, tile_y, tileID_table, panelID_table, &
      NprcX_lc, NprcY_lc )
    use scale_prc
    use scale_cubedsphere_cnv, only: &
      CubedSphereCnv_LonLat2CSPos
    implicit none

    class(NodeMappingInfo), intent(inout), target :: this
    type(LocalMesh2D), intent(in) :: lcmesh
    type(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(in) :: tile_x(4,in_NLocalMeshPerPrc,in_Nprc)
    real(RP), intent(in) :: tile_y(4,in_NLocalMeshPerPrc,in_Nprc)
    integer, intent(in) :: tileID_table(in_NLocalMeshPerPrc,in_Nprc)
    integer, intent(in) :: panelID_table(in_NLocalMeshPerPrc,in_Nprc)
    integer, intent(in) :: NprcX_lc, NprcY_lc

    integer :: lc_domID, prcID, in_prc, in_n
    integer :: i, j
    integer :: p, ke
    integer :: ke_h
    integer :: p_h, p_h_x, p_h_y
    integer :: ke_z, ke_z2, p_z, p_z2
    logical :: is_inside_tile(elem2D%Nfp**2)
    logical :: is_inside_elem
    real(RP) :: in_elem_x(4)
    real(RP) :: in_elem_y(4)
    real(RP) :: delx, dely
    logical :: target_tile_flag(in_NLocalMeshPerPrc,in_Nprc)
    logical :: target_prc_flag(in_Nprc)
    integer :: in_tileID_tmp(in_NLocalMeshPerPrc*in_Nprc)
    integer :: in_prc2lcprc_tmp(in_Nprc)
    integer :: in_lcprc2prc_tmp(in_Nprc)
    integer :: in_tile_num
    integer :: in_prc_num
    type(LocalMesh2D), pointer :: in_lcmesh
    real(RP) :: in_Z0, in_Z1

    real(RP) :: out_lon(1), out_lat(1)
    real(RP) :: out_x(1), out_y(1)
    integer :: out_cspanel
    real(RP) :: out_x0, del_x
    real(RP) :: out_y0, del_y

    integer :: in_NeX, in_NeY
    !-------------------------------------------

    allocate( this%local_domID(elem2D%Nfp**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%lcprc(elem2D%Nfp**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%inCSPanelID(elem2D%Nfp**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_i(elem2D%Nfp**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_j(elem2D%Nfp**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_x(elem2D%Nfp**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_y(elem2D%Nfp**2, lcmesh%NeX*lcmesh%NeY) )

    target_tile_flag(:,:) = .false.
    target_prc_flag(:) = .false.
    in_lcprc2prc_tmp(:) = -1
    in_prc2lcprc_tmp(:) = -1

    call get_panelID( this%inCSPanelID(:,:),      &
      lcmesh%pos_en(:,:,1) * PI / 180.0_RP,       &
      lcmesh%pos_en(:,:,2) * PI / 180.0_RP,       &
      elem2D%Nfp**2 * lcmesh%NeX*lcmesh%NeY       )
    
    in_tile_num = 0
    in_prc_num = 0
    do ke_h=1, lcmesh%NeX * lcmesh%NeY

      do p_h_y=1, elem2D%Nfp
      do p_h_x=1, elem2D%Nfp
        p_h = p_h_x + (p_h_y-1)*elem2D%Nfp
        this%local_domID(p_h,ke_h) = -1
        this%lcprc(p_h,ke_h) = -1

        out_cspanel = this%inCSPanelID(p_h,ke_h)
        out_lon(1) = lcmesh%pos_en(p_h,ke_h,1) * PI / 180.0_RP
        out_lat(1) = lcmesh%pos_en(p_h,ke_h,2) * PI / 180.0_RP
        call CubedSphereCnv_LonLat2CSPos( &
          out_cspanel, out_lon, out_lat, 1, &
          out_x, out_y )

        loop_prc: do in_prc=1, in_Nprc
        do in_n=1, in_NLocalMeshPerPrc
          
          if ( out_cspanel /= panelID_table(in_n,in_prc) ) cycle

          is_inside_tile(p_h) = inpoly( out_x(1), out_y(1),                             &
                                        4, tile_x(:,in_n,in_prc), tile_y(:,in_n,in_prc) )

          if ( is_inside_tile(p_h) ) then
            this%local_domID(p_h,ke_h) = in_n

            if ( .not. target_tile_flag(in_n,in_prc) ) then
              target_tile_flag(in_n,in_prc) = .true.
              in_tile_num = in_tile_num + 1
              in_tileID_tmp(in_tile_num) = tileID_table(in_n,in_prc)
            end if
            if ( .not. target_prc_flag(in_prc) ) then
              target_prc_flag(in_prc) = .true.
              in_prc_num = in_prc_num + 1
              in_prc2lcprc_tmp(in_prc) = in_prc_num
              in_lcprc2prc_tmp(in_prc_num) = in_prc
            end if
            this%lcprc(p_h,ke_h) = in_prc2lcprc_tmp(in_prc)

            exit loop_prc
          end if

        end do
        end do loop_prc

      end do
      end do

    end do


    allocate( this%in_tileID_list(in_tile_num) )
    this%in_tileID_list(:) = in_tileID_tmp(1:in_tile_num)

    in_NeX = in_NeGX / NprcX_lc
    in_NeY = in_NeGY / NprcY_lc

    do ke_h=1, lcmesh%NeX * lcmesh%NeY
      do p_h_y=1, elem2D%Nfp
      do p_h_x=1, elem2D%Nfp

        p_h = p_h_x + (p_h_y-1)*elem2D%Nfp

        lc_domID = this%local_domID(p_h,ke_h)
        prcID = in_lcprc2prc_tmp( this%lcprc(p_h,ke_h) )

        out_cspanel = this%inCSPanelID(p_h,ke_h)
        out_lon(1) = lcmesh%pos_en(p_h,ke_h,1) * PI / 180.0_RP
        out_lat(1) = lcmesh%pos_en(p_h,ke_h,2) * PI / 180.0_RP
        call CubedSphereCnv_LonLat2CSPos( &
          out_cspanel, out_lon, out_lat, 1, &
          out_x, out_y )

        this%elem_i(p_h,ke_h) = -1
        this%elem_j(p_h,ke_h) = -1
        if ( lc_domID > 0 .and. prcID > 0) then
          delx = ( tile_x(2,lc_domID,prcID) - tile_x(1,lc_domID,prcID) ) / dble(in_NeX)
          dely = ( tile_y(4,lc_domID,prcID) - tile_y(1,lc_domID,prcID) ) / dble(in_NeY) 
  
          loop_ne: do j=1, in_NeY
          do i=1, in_NeX        
            in_elem_x(:) = tile_x(1,lc_domID,prcID) + delx * dble( (/ i-1, i, i, i-1 /) )
            in_elem_y(:) = tile_y(1,lc_domID,prcID) + dely * dble( (/ j-1, j-1, j, j /) )
            is_inside_elem  = inpoly( out_x(1), out_y(1),           &
                                      4, in_elem_x(:), in_elem_y(:) )
            if (i==in_NeX) then
              in_elem_x(2:3) = in_elem_x(2:3) + 1.0E-12_RP * delx       
            end if
            if (j==in_NeY) then
              in_elem_y(3:4) = in_elem_y(3:4) + 1.0E-12_RP * dely
            end if
                                                            
            if (is_inside_elem) then
              this%elem_i(p_h,ke_h) = i
              this%elem_j(p_h,ke_h) = j
              this%elem_x(p_h,ke_h) = out_x(1)
              this%elem_y(p_h,ke_h) = out_y(1)

              exit loop_ne
            end if
          end do  
          end do loop_ne
        end if
      end do
      end do
    end do

    !-- prepair mesh for input data --------------------------------

    allocate( this%in_mesh2D_list(in_prc_num) )
    do i=1, in_prc_num
      call this%in_mesh2D_list(i)%Init( in_NeGX, in_NeGY,  &
        RPlanet, in_elem2D, in_NLocalMeshPerPrc,           &
        nproc=in_Nprc, myrank=in_lcprc2prc_tmp(i)-1        )
      
      call this%in_mesh2D_list(i)%Generate()
    end do

    return
  end subroutine NodeMappingInfo_Init_2D

!OCL SERIAL
  subroutine NodeMappingInfo_Init_3D( this, &
    lcmesh, elem3D, tile_x, tile_y, in_Fz,  &
    tileID_table, panelID_table,            &
    NprcX_lc, NprcY_lc, is_spec_in_FZ )
    use scale_prc
    use scale_cubedsphere_cnv, only: &
      CubedSphereCnv_LonLat2CSPos
    implicit none

    class(NodeMappingInfo), intent(inout), target :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    type(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(in) :: tile_x(4,in_NLocalMeshPerPrc,in_Nprc)
    real(RP), intent(in) :: tile_y(4,in_NLocalMeshPerPrc,in_Nprc)
    real(RP), intent(in) :: in_Fz(1:in_NeGZ+1)
    integer, intent(in) :: tileID_table(in_NLocalMeshPerPrc,in_Nprc)
    integer, intent(in) :: panelID_table(in_NLocalMeshPerPrc,in_Nprc)
    integer, intent(in) :: NprcX_lc, NprcY_lc
    logical, intent(in) :: is_spec_in_FZ

    integer :: lc_domID, prcID, in_prc, in_n
    integer :: i, j
    integer :: p, ke
    integer :: ke_h
    integer :: p_h, p_h_x, p_h_y
    integer :: ke_z, ke_z2, p_z, p_z2
    logical :: is_inside_tile(elem3D%Nnode_h1D**2)
    logical :: is_inside_elem
    real(RP) :: in_elem_x(4)
    real(RP) :: in_elem_y(4)
    real(RP) :: delx, dely
    logical :: target_tile_flag(in_NLocalMeshPerPrc,in_Nprc)
    logical :: target_prc_flag(in_Nprc)
    integer :: in_tileID_tmp(in_NLocalMeshPerPrc*in_Nprc)
    integer :: in_prc2lcprc_tmp(in_Nprc)
    integer :: in_lcprc2prc_tmp(in_Nprc)
    integer :: in_tile_num
    integer :: in_prc_num
    type(LocalMesh3D), pointer :: in_lcmesh
    real(RP) :: in_Z0, in_Z1
    integer :: in_ke3D

    real(RP) :: out_lon(1), out_lat(1)
    real(RP) :: out_x(1), out_y(1)
    integer :: out_cspanel
    real(RP) :: out_x0, del_x
    real(RP) :: out_y0, del_y

    integer :: in_NeX, in_NeY

    real(RP) :: out_lon2D(elem3D%Nnode_h1D**2,lcmesh%NeX*lcmesh%NeY)
    real(RP) :: out_lat2D(elem3D%Nnode_h1D**2,lcmesh%NeX*lcmesh%NeY)
    !-------------------------------------------

    allocate( this%local_domID(elem3D%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%lcprc(elem3D%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%inCSPanelID(elem3D%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_i(elem3D%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_j(elem3D%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_k(elem3D%Np, lcmesh%NeX*lcmesh%NeY*lcmesh%NeZ) )
    allocate( this%elem_x(elem3D%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_y(elem3D%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_z(elem3D%Np, lcmesh%NeX*lcmesh%NeY*lcmesh%NeZ) )

    target_tile_flag(:,:) = .false.
    target_prc_flag(:)    = .false.
    in_lcprc2prc_tmp(:)   = -1
    in_prc2lcprc_tmp(:)   = -1

    !$omp parallel do
    do ke_h=1, lcmesh%NeX * lcmesh%NeY
      out_lon2D(:,ke_h) = lcmesh%pos_en(elem3D%Hslice(:,1),ke_h,1) * PI / 180.0_RP
      out_lat2D(:,ke_h) = lcmesh%pos_en(elem3D%Hslice(:,1),ke_h,2) * PI / 180.0_RP
    end do
    call get_panelID( this%inCSPanelID(:,:),      &
      out_lon2D, out_lat2D,                       &
      elem3D%Nnode_h1D**2 * lcmesh%NeX*lcmesh%NeY )

    in_tile_num = 0
    in_prc_num = 0
    do ke_h=1, lcmesh%NeX * lcmesh%NeY

      do p_h_y=1, elem3D%Nnode_h1D
      do p_h_x=1, elem3D%Nnode_h1D
        p_h = p_h_x + (p_h_y-1)*elem3D%Nnode_h1D
        this%local_domID(p_h,ke_h) = -1
        this%lcprc(p_h,ke_h)       = -1

        out_cspanel = this%inCSPanelID(p_h,ke_h)
        out_lon(1) = lcmesh%pos_en(p_h,ke_h,1) * PI / 180.0_RP
        out_lat(1) = lcmesh%pos_en(p_h,ke_h,2) * PI / 180.0_RP
        call CubedSphereCnv_LonLat2CSPos(   &
          out_cspanel, out_lon, out_lat, 1, &
          out_x, out_y )

        loop_prc: do in_prc=1, in_Nprc
        do in_n=1, in_NLocalMeshPerPrc
          
          if ( out_cspanel /= panelID_table(in_n,in_prc) ) cycle

          is_inside_tile(p_h) = inpoly( out_x(1), out_y(1),                             &
                                        4, tile_x(:,in_n,in_prc), tile_y(:,in_n,in_prc) )

          if ( is_inside_tile(p_h) ) then
            this%local_domID(p_h,ke_h) = in_n

            if ( .not. target_tile_flag(in_n,in_prc) ) then
              target_tile_flag(in_n,in_prc) = .true.
              in_tile_num = in_tile_num + 1
              in_tileID_tmp(in_tile_num) = tileID_table(in_n,in_prc)
            end if
            if ( .not. target_prc_flag(in_prc) ) then
              target_prc_flag(in_prc) = .true.
              in_prc_num = in_prc_num + 1
              in_prc2lcprc_tmp(in_prc) = in_prc_num
              in_lcprc2prc_tmp(in_prc_num) = in_prc
            end if
            this%lcprc(p_h,ke_h) = in_prc2lcprc_tmp(in_prc)

            exit loop_prc
          end if
        end do
        end do loop_prc
      end do
      end do

    end do

    allocate( this%in_tileID_list(in_tile_num) )
    this%in_tileID_list(:) = in_tileID_tmp(1:in_tile_num)

    in_NeX = in_NeGX / NprcX_lc
    in_NeY = in_NeGY / NprcY_lc

    do ke_h=1, lcmesh%NeX * lcmesh%NeY
      do p_h_y=1, elem3D%Nnode_h1D
      do p_h_x=1, elem3D%Nnode_h1D

        p_h = p_h_x + (p_h_y-1)*elem3D%Nnode_h1D

        lc_domID = this%local_domID(p_h,ke_h)
        prcID = in_lcprc2prc_tmp( this%lcprc(p_h,ke_h) )

        out_cspanel = this%inCSPanelID(p_h,ke_h)
        out_lon(1) = lcmesh%pos_en(p_h,ke_h,1) * PI / 180.0_RP
        out_lat(1) = lcmesh%pos_en(p_h,ke_h,2) * PI / 180.0_RP
        call CubedSphereCnv_LonLat2CSPos( &
          out_cspanel, out_lon, out_lat, 1, &
          out_x, out_y )
  
        this%elem_i(p_h,ke_h) = -1
        this%elem_j(p_h,ke_h) = -1
        if ( lc_domID > 0 .and. prcID > 0 ) then
          delx = ( tile_x(2,lc_domID,prcID) - tile_x(1,lc_domID,prcID) ) / dble(in_NeX)
          dely = ( tile_y(4,lc_domID,prcID) - tile_y(1,lc_domID,prcID) ) / dble(in_NeY) 

          loop_ne: do j=1, in_NeY
          do i=1, in_NeX        
            in_elem_x(:) = tile_x(1,lc_domID,prcID) + delx * dble( (/ i-1, i, i, i-1 /) )
            in_elem_y(:) = tile_y(1,lc_domID,prcID) + dely * dble( (/ j-1, j-1, j, j /) )

            if (i==in_NeX) then
              in_elem_x(2:3) = in_elem_x(2:3) + 1.0E-12_RP * delx       
            end if
            if (j==in_NeY) then
              in_elem_y(3:4) = in_elem_y(3:4) + 1.0E-12_RP * dely
            end if

            is_inside_elem  = inpoly( out_x(1), out_y(1),           &
                                      4, in_elem_x(:), in_elem_y(:) )
                            
            if (is_inside_elem) then
              this%elem_i(p_h,ke_h) = i
              this%elem_j(p_h,ke_h) = j
              this%elem_x(p_h,ke_h) = out_x(1)
              this%elem_y(p_h,ke_h) = out_y(1)
              exit loop_ne
            end if
          end do  
          end do loop_ne
        end if

      end do
      end do
    end do

    !-- prepair mesh for input data --------------------------------

    allocate( this%in_mesh3D_list(in_prc_num) )
    do i=1, in_prc_num
      if (is_spec_in_FZ) then
        call this%in_mesh3D_list(i)%Init( in_NeGX, in_NeGY, in_NeGZ,    &
          RPlanet, dom_zmin, dom_zmax, in_elem3D, in_NLocalMeshPerPrc,  &
          nproc=in_Nprc, myrank=in_lcprc2prc_tmp(i)-1,                  &
          FZ=in_Fz                                                      )
      else
        call this%in_mesh3D_list(i)%Init( in_NeGX, in_NeGY, in_NeGZ,    &
          RPlanet, dom_zmin, dom_zmax, in_elem3D, in_NLocalMeshPerPrc,  &
          nproc=in_Nprc, myrank=in_lcprc2prc_tmp(i)-1                   )
      end if
      call this%in_mesh3D_list(i)%Generate()
    end do

    !--
    !$omp parallel private( &
    !$omp ke_h, p_h, in_lcmesh, in_prc, in_n,            &
    !$omp ke_z, ke_z2, p_z, ke, p, in_ke3D, in_Z0, in_Z1 )

    !$omp workshare
    this%elem_k(:,:) = -1
    !$omp end workshare

    !$omp do collapse(2)
    do ke_h=1, lcmesh%NeX * lcmesh%NeY
    do p_h=1, elem3D%Nnode_h1D**2

      in_n   = this%local_domID(p_h,ke_h)
      in_prc = this%lcprc(p_h,ke_h)

      if ( in_n > 0 .and. in_prc > 0 .and.                           &
           this%elem_i(p_h,ke_h) > 0 .and. this%elem_j(p_h,ke_h) > 0 ) then
        
        in_lcmesh => this%in_mesh3D_list(in_prc)%lcmesh_list(in_n)
        do ke_z=1, lcmesh%NeZ
        do p_z=1, out_elem3D%Nnode_v
          ke = ke_h + (ke_z-1)*lcmesh%NeX*lcmesh%NeY
          p = p_h + (p_z-1)*elem3D%Nnode_h1D**2
          do ke_z2=1, in_NeGZ
            in_ke3D = this%elem_i(p_h,ke_h) + (this%elem_j(p_h,ke_h)-1)*in_NeX &
                    + (ke_z2-1)*in_NeX*in_NeY
            in_Z0 = in_lcmesh%pos_ev(in_lcmesh%EToV(in_ke3D,1),3)
            in_Z1 = in_lcmesh%pos_ev(in_lcmesh%EToV(in_ke3D,5),3)
            if ( in_Z0 <= lcmesh%pos_en(p,ke,3) .and. lcmesh%pos_en(p,ke,3) <= in_Z1 ) then
              this%elem_k(p,ke) = ke_z2
              this%elem_z(p,ke) = lcmesh%pos_en(p,ke,3)
              exit
            end if
          end do
        end do
        end do
      end if

    end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine NodeMappingInfo_Init_3D

!OCL SERIAL  
  subroutine NodeMappingInfo_Final( this )
    implicit none
    class(NodeMappingInfo), intent(inout) :: this

    integer :: i
    !-------------------------------------------

    if ( allocated(this%local_domID) ) then
      deallocate( this%local_domID, this%lcprc )
      deallocate( this%inCSPanelID )
      deallocate( this%elem_i, this%elem_j )
      if ( allocated(this%elem_k) ) deallocate( this%elem_k )
      deallocate( this%in_tileID_list )

      if (is_mesh3D) then
        do i=1, size(this%in_mesh3D_list)
          call this%in_mesh3D_list(i)%Final()
        end do
        deallocate( this%in_mesh3D_list )
      else
        do i=1, size(this%in_mesh2D_list)
          call this%in_mesh2D_list(i)%Final()
        end do
        deallocate( this%in_mesh2D_list )
      end if
    end if

    return
  end subroutine NodeMappingInfo_Final

  !> Check whether the point is located inside a polyngon
!OCL SERIAL
  function inpoly( pt_x, pt_y, num_node, v_x, v_y ) result(ret)
    implicit none
    real(RP), intent(in) :: pt_x
    real(RP), intent(in) :: pt_y
    integer, intent(in) :: num_node
    real(RP), intent(in) :: v_x(num_node)
    real(RP), intent(in) :: v_y(num_node)
    logical :: ret

    integer :: wn
    integer :: i, ii
    !------------------------------------------

    wn = 0
    do i=1, num_node
      ii = mod(i, num_node) + 1
      if ( v_y(i) <= pt_y .and. pt_y < v_y(ii)) then
        if( pt_x < v_x(i) + (pt_y - v_y(i)) * (v_x(ii) - v_x(i))/(v_y(ii) - v_y(i)) ) then
          wn = wn + 1
        end if
      else if ( v_y(i) > pt_y .and. v_y(ii) <= pt_y ) then
        if( pt_x < v_x(i) + (pt_y - v_y(i)) * (v_x(ii) - v_x(i))/(v_y(ii) - v_y(i)) ) then
          wn = wn - 1
        end if
      end if
    end do

    if (wn == 0) then
      ret = .false.
    else
      ret = .true.
    end if

    return
  end function inpoly

!OCL SERIAL
  subroutine get_panelID( panelID, lon, lat, Np )
    use scale_cubedsphere_cnv, only: &
      CubedSphereCnv_LonLat2CSPos
    implicit none

    integer, intent(in) :: Np
    integer, intent(out) :: panelID(Np)
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)

    integer :: p
    integer :: pnl
    real(RP) :: alph(Np), beta(Np)
    real(RP) :: lon_(Np), lat_(Np)
    real(RP), parameter :: EPS = 1.0E-64_RP
    !------------------------------------------

    !$omp parallel do
    do p=1, Np
      if ( abs(cos(lon(p))) < EPS ) then
        lon_(p) = lon(p) + EPS
      else
        lon_(p) = lon(p)
      end if
      if ( abs( lat(p) - 0.5_RP * PI ) < EPS ) then
        lat_(p) = lat(p) - sign(EPS, lat(p))
      else
        lat_(p) = lat(p)
      end if
    end do

    panelID(:) = -1
    do pnl=1, 6
      call CubedSphereCnv_LonLat2CSPos( pnl, lon_, lat_, Np, &
        alph, beta )
    
      select case(pnl)
      case (5)
        where ( lat(:) > 0.0_RP .and. abs(alph(:)) <= 0.25_RP * PI .and. abs(beta(:)) <= 0.25_RP * PI  )
          panelID(:) = pnl
        end where
      case (6)
        where ( lat(:) < 0.0_RP .and. abs(alph(:)) <= 0.25_RP * PI .and. abs(beta(:)) <= 0.25_RP * PI  )
          panelID(:) = pnl
        end where
      case default 
        where ( abs(alph(:)) <= 0.25_RP * PI .and. abs(beta(:)) <= 0.25_RP * PI  )
          panelID(:) = pnl
        end where
      end select     
    end do

    do p=1, Np
      if (panelID(p) <  0) then
        LOG_ERROR("get_panelID",*) 'Fail to search a panel ID of cubed sphere grid!'
        write(*,*) "p=", p, ": (lon,lat)=", lon(p), lat(p), "(alpha,beta)=", alph(p), beta(p) 
        call PRC_abort
      end if 
    end do

    return
  end subroutine get_panelID

end module mod_cs2lonlat_interp_mesh