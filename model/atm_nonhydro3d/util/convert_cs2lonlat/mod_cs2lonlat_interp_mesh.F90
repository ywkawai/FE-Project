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
  
  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedspheredom2d, only: &
    MeshCubedSphereDom2D,                       &
    MeshCubedSphereDom2D_check_division_params

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
  type(MeshCubedSphereDom2D), public :: in_csmesh2D

  type, public :: NodeMappingInfo
    type(LocalMesh2D), pointer :: lcmesh2D
    integer, allocatable :: local_domID(:,:)
    integer, allocatable :: lcprc(:,:)
    integer, allocatable :: inCSPanelID(:,:)
    integer, allocatable :: elem_i(:,:)
    integer, allocatable :: elem_j(:,:)
    integer, allocatable :: elem_k(:,:)
    real(RP), allocatable :: elem_x(:,:)
    real(RP), allocatable :: elem_y(:,:)

    type(MeshCubedSphereDom2D), allocatable :: in_mesh2D_list(:)
    integer, allocatable :: in_tileID_list(:)
  contains
    procedure :: Init => NodeMappingInfo_Init
    procedure :: Final => NodeMappingInfo_Final
  end type
  type(NodeMappingInfo), public :: nodeMap_list(ATMOS_MESH_NLocalMeshPerPrc)

  integer, public :: in_Nprc             = 1
  integer, public :: in_NeGX             = 1
  integer, public :: in_NeGY             = 1
  integer, public :: in_NLocalMeshPerPrc = 6  
  type(QuadrilateralElement), public :: in_elem2D
  
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

contains
  subroutine interp_mesh_Init()
    implicit none

  
    integer :: out_NeX          = 1
    integer :: out_NeY          = 1
    integer :: in_PolyOrder_h   = 1
    integer :: out_PolyOrder_h  = 1

    namelist / PARAM_INTERP_MESH / &
      in_Nprc,             &
      in_NeGX,             &
      in_NeGY,             &
      in_NLocalMeshPerPrc, &
      in_PolyOrder_h,      &
      out_NprcX,           &
      out_NeX,             &
      out_NprcY,           &
      out_NeY,             &
      out_PolyOrder_h
    
    integer :: ierr
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("interp_mesh",*) 'Setup'
  
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

    !-    
    call in_elem2D%Init( in_PolyOrder_h, .false. )
    call out_elem2D%Init( out_PolyOrder_h, .false. )

    call out_mesh2D%Init( out_NprcX * out_NeX, out_NprcY * out_NeY, &
       0.0_RP, 2.0_RP*PI, -0.5_RP*PI, 0.5_RP*PI,                    &
       .true., .false., out_elem2D, 1,                              &
       out_NprcX, out_NprcY          )
    
    call out_mesh2D%Generate()

    !-

    call construct_map()


    return
  end subroutine interp_mesh_Init

  subroutine interp_mesh_Final()
    implicit none

    integer :: n
    !-------------------------------------------

!    call in_mesh%Final()
    call in_elem2D%Final()

    call out_mesh2D%Final()
    call out_elem2D%Final()

    do n=1, ATMOS_MESH_NLocalMeshPerPrc
      call nodeMap_list(n)%Final()
    end do
    
    return
  end subroutine interp_mesh_Final

  subroutine interp_mesh_Interpolate_field()
    implicit none
    !-------------------------------------------

    call in_elem2D%Final()

    call out_mesh2D%Final()
    call out_elem2D%Final()

    return
  end subroutine interp_mesh_Interpolate_field

  !- private -------------------------------------

!OCL SERIAL  
  subroutine construct_map()
    implicit none

    integer :: n
    integer :: i, j
    real(RP) :: in_tiles_x(4,in_NLocalMeshPerPrc,in_Nprc)
    real(RP) :: in_tiles_y(4,in_NLocalMeshPerPrc,in_Nprc)
    integer :: in_panel(4,in_NLocalMeshPerPrc,in_Nprc)
    real(RP) :: delx, dely
    type(LocalMesh2D), pointer :: lcmesh

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
    call in_csmesh2D_dummy%Init( in_NeGX, in_NeGY, RPlanet, &
      in_elem2D, in_NLocalMeshPerPrc, nproc=in_Nprc         )
    
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

    do n=1, out_mesh2D%LOCAL_MESH_NUM
      lcmesh => out_mesh2D%lcmesh_list(n)
      call nodeMap_list(n)%Init( lcmesh, lcmesh%refElem2D, in_tiles_x, in_tiles_y, &
        tileID_table, panelID_table, NprcX_lc, NprcY_lc )
    end do
    
    return
  end subroutine construct_map

!OCL SERIAL
  subroutine NodeMappingInfo_Init( this, &
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
      lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), &
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
        out_lon(1) = lcmesh%pos_en(p_h,ke_h,1)
        out_lat(1) = lcmesh%pos_en(p_h,ke_h,2)
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
        out_lon(1) = lcmesh%pos_en(p_h,ke_h,1)
        out_lat(1) = lcmesh%pos_en(p_h,ke_h,2)
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
  end subroutine NodeMappingInfo_Init

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
      deallocate( this%in_tileID_list )

      do i=1, size(this%in_mesh2D_list)
        call this%in_mesh2D_list(i)%Final()
      end do   
      deallocate( this%in_mesh2D_list )
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

    do p=1, Np
      if ( abs(lon(p)) < EPS ) then
        lon_(p) = lon(p) + EPS
      else
        lon_(p) = lon(p)
      end if
      if ( abs(lat(p)) < EPS ) then
        lat_(p) = lat(p) + EPS
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