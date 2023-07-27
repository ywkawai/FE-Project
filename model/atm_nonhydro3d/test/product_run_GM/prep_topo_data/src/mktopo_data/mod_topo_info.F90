#include "scalelib.h"
module mod_topo_info
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
    PRC_abort
  use scale_const, only: &
    PI => CONST_PI

  use netcdf

  use scale_element_base, only: ElementBase2D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D

  use mod_topo_nodemap, only: &
    TOPO_nodemap_info
  
  implicit none
  private

  type, public :: TOPO_info
    real(RP), allocatable :: lon(:)
    real(RP), allocatable :: lat(:)
    real(RP), allocatable :: topo(:,:)
    
    integer :: nlon, nlat
  contains
    procedure :: Init => TOPO_info_Init
    procedure :: Final => TOPO_info_Final
    procedure :: Read_data => TOPO_info_read_data
  end type

  public :: TOPO_info_list_setup
  public :: TOPO_info_list_final

  type(TOPO_info), public, allocatable, target :: TOPO_info_list(:,:)
  type(TOPO_nodemap_info), public, allocatable, target :: node_map_info_list(:)

  integer :: tile_num_lon
  integer :: tile_num_lat
  logical, allocatable :: TOPO_info_flag(:,:)
  
contains
  subroutine TOPO_info_list_setup( &
    tile_num_lon_, tile_nlon, tile_num_lat_, tile_nlat, &
    in_topo_filebase, mesh2D )

    implicit none
    integer, intent(in) :: tile_num_lon_
    integer, intent(in) :: tile_nlon
    integer, intent(in) :: tile_num_lat_
    integer, intent(in) :: tile_nlat
    character(len=*), intent(in) :: in_topo_filebase
    class(MeshCubedSphereDom2D), target, intent(in) :: mesh2D
    
    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase2D), pointer :: elem2D

    integer :: domid
    integer :: ke, p

    integer :: file_i, file_j
    integer :: file_ii, file_jj
    integer :: i, j
    
    type(TOPO_nodemap_info), pointer :: node_map
    type(TOPO_info), pointer :: tinfo

    real(RP) :: fac_rad2deg
    real(RP) :: neigh_topo(0:TILE_nlon+1,0:TILE_nlat+1)

    real(RP) :: tmp_pole(2)
    !-----------------------------------------

    tile_num_lon = tile_num_lon_
    tile_num_lat = tile_num_lat_

    fac_rad2deg = 180.0E0_RP / PI

    allocate( TOPO_info_flag(tile_num_lon,tile_num_lat) )
    allocate( TOPO_info_list(tile_num_lon,tile_num_lat) )
    allocate( node_map_info_list(mesh2D%LOCAL_MESH_NUM) )

    TOPO_info_flag(:,:) = .false.

    do domid=1, mesh2D%LOCAL_MESH_NUM
      lcmesh2D => mesh2D%lcmesh_list(domid)
      elem2D => lcmesh2D%refElem2D

      node_map => node_map_info_list(domid)      
      call node_map%Init(lcmesh2D, elem2D)
      call node_map%Set_file_ijindx( &
        lcmesh2D, elem2D, TILE_num_lon, TILE_num_lat, &
        TOPO_info_flag )
    end do

    do file_j=1, TILE_num_lat
    do file_i=1, TILE_num_lon
      if (topo_info_flag(file_i,file_j)) then
        LOG_INFO("MKTOPO_DATA",*) "Initialize topo_info i,j=", file_i, file_j
        call topo_info_list(file_i,file_j)%Init(TILE_nlon, TILE_nlat)

        call topo_info_list(file_i,file_j)%Read_data(in_topo_filebase, file_i, file_j)
      end if
    end do
    end do

    do domid=1, mesh2D%LOCAL_MESH_NUM
      lcmesh2D => mesh2D%lcmesh_list(domid)
      elem2D => lcmesh2D%refElem2D
      node_map => node_map_info_list(domid)      

      do ke=lcmesh2D%NeS, lcmesh2D%NeE
      do p=1, elem2D%Np
        file_i = node_map%file_i_ind(p,ke)
        file_j = node_map%file_j_ind(p,ke)
        tinfo => topo_info_list(file_i,file_j)

        loop_ne: do j=-1, TILE_nlat
        do i=-1, TILE_nlon
          if (  tinfo%lon(i) <= lcmesh2D%lon(p,ke) .and. tinfo%lon(i+1) >= lcmesh2D%lon(p,ke) &
          .and. tinfo%lat(j) <= lcmesh2D%lat(p,ke) .and. tinfo%lat(j+1) >= lcmesh2D%lat(p,ke) ) then
            node_map%i_ind(p,ke) = i
            node_map%j_ind(p,ke) = j
            exit loop_ne
          end if
        end do    
        end do loop_ne
      end do
      end do

      do ke=lcmesh2D%NeS, lcmesh2D%NeE
      do p=1, elem2D%Np
        if (node_map%i_ind(p,ke) < 0 .or. node_map%j_ind(p,ke) < 0) then
          file_i = node_map%file_i_ind(p,ke)
          file_j = node_map%file_j_ind(p,ke)
          LOG_INFO("TOPO_info_list_setup",*) "Error", p, ke, "lon, lat=", lcmesh2D%lon(p,ke)*fac_rad2deg, lcmesh2D%lat(p,ke)*fac_rad2deg
          LOG_INFO("TOPO_info_list_setup",*) "file_i, file_j", file_i, file_j
          LOG_INFO("TOPO_info_list_setup",*) "lon range:", TOPO_info_list(file_i,file_j)%lon(1)*fac_rad2deg, TOPO_info_list(file_i,file_j)%lon(tile_nlon)*fac_rad2deg
          LOG_INFO("TOPO_info_list_setup",*) "lat range:", TOPO_info_list(file_i,file_j)%lat(1)*fac_rad2deg, TOPO_info_list(file_i,file_j)%lat(tile_nlat)*fac_rad2deg
          call PRC_abort
        end if
      end do
      end do
    end do
    
    !---

    do file_j=1, TILE_num_lat
    do file_i=1, TILE_num_lon
      if ( .not. TOPO_info_flag(file_i,file_j) ) cycle

      LOG_INFO("MKTOPO_DATA",*) "Copy boundary for tile i,j=", file_i, file_j

      ! WEST
      file_ii = file_i - 1; file_jj = file_j
      if ( file_i == 1 ) file_ii = TILE_num_lon 

      if ( .not. TOPO_info_flag(file_ii,file_jj) ) then
        call get_neighbor_topo_data( TILE_nlon, TILE_nlat, in_topo_filebase, file_ii, file_jj, neigh_topo)
        TOPO_info_list(file_i,file_j)%topo(0,1:TILE_nlat) = neigh_topo(TILE_nlon,1:TILE_nlat)
      else
        TOPO_info_list(file_i,file_j)%topo(0,1:TILE_nlat) = TOPO_info_list(file_ii,file_jj)%topo(TILE_nlon,1:TILE_nlat)
      end if

      ! EAST
      file_ii = file_i + 1; file_jj = file_j
      if ( file_i == TILE_num_lon ) file_ii = 1

      if ( .not. TOPO_info_flag(file_ii,file_jj) ) then
        call get_neighbor_topo_data( TILE_nlon, TILE_nlat, in_topo_filebase, file_ii, file_jj, neigh_topo)
        TOPO_info_list(file_i,file_j)%topo(TILE_nlon+1,1:TILE_nlat) = neigh_topo(1,1:TILE_nlat)
      else
        TOPO_info_list(file_i,file_j)%topo(TILE_nlon+1,1:TILE_nlat) = TOPO_info_list(file_ii,file_jj)%topo(1,1:TILE_nlat)
      end if
      
      ! SOUTH
      file_ii = file_i; file_jj = file_j - 1
      if ( file_j > 1 ) then
        if ( .not. TOPO_info_flag(file_ii,file_jj) ) then
          call get_neighbor_topo_data( TILE_nlon, TILE_nlat, in_topo_filebase, file_ii, file_jj, neigh_topo)
          TOPO_info_list(file_i,file_j)%topo(1:TILE_nlon,0) = neigh_topo(1:TILE_nlon,TILE_nlat)
        else
          TOPO_info_list(file_i,file_j)%topo(1:TILE_nlon,0) = TOPO_info_list(file_ii,file_jj)%topo(1:TILE_nlon,TILE_nlat)
        end if
      end if
      
      ! NORTH
      file_ii = file_i; file_jj = file_j + 1
      if ( file_j < TILE_num_lat ) then
        if ( .not. TOPO_info_flag(file_ii,file_jj) ) then
          call get_neighbor_topo_data( TILE_nlon, TILE_nlat, in_topo_filebase, file_ii, file_jj, neigh_topo)
          TOPO_info_list(file_i,file_j)%topo(1:TILE_nlon,TILE_nlat+1) = neigh_topo(1:TILE_nlon,1)
        else          
          TOPO_info_list(file_i,file_j)%topo(1:TILE_nlon,TILE_nlat+1) = TOPO_info_list(file_ii,file_jj)%topo(1:TILE_nlon,1)
        end if
      end if

    end do
    end do

    ! Pole
    tmp_pole(:) = 0.0_RP
    do file_i=1, TILE_num_lon

      file_ii = file_i; file_jj = TILE_num_lat
      if ( .not. TOPO_info_flag(file_ii,file_jj) ) then
        call get_neighbor_topo_data( TILE_nlon, TILE_nlat, in_topo_filebase, file_ii, file_jj, neigh_topo)
        tmp_pole(1) = tmp_pole(1) + sum(neigh_topo(1:TILE_nlon,TILE_nlat))
      else          
        tmp_pole(1) = tmp_pole(1) + sum( TOPO_info_list(file_ii,file_jj)%topo(1:TILE_nlon,TILE_nlat) )
      end if

      file_ii = file_i; file_jj = 1
      if ( .not. TOPO_info_flag(file_ii,file_jj) ) then
        call get_neighbor_topo_data( TILE_nlon, TILE_nlat, in_topo_filebase, file_ii, file_jj, neigh_topo)
        tmp_pole(2) = tmp_pole(2) + sum(neigh_topo(1:TILE_nlon,1))
      else          
        tmp_pole(2) = tmp_pole(2) + sum( TOPO_info_list(file_ii,file_jj)%topo(1:TILE_nlon,1) )
      end if
    end do      

    do file_i=1, TILE_num_lon
      file_j = TILE_num_lat
      if ( TOPO_info_flag(file_i,file_j) ) then
        TOPO_info_list(file_i,file_j)%topo(1:TILE_nlon,TILE_nlat+1) = tmp_pole(1) / dble( TILE_nlon * TILE_num_lon )
      end if

      file_j = 1
      if ( TOPO_info_flag(file_i,file_j) ) then
        TOPO_info_list(file_i,file_j)%topo(1:TILE_nlon,0) = tmp_pole(2) / dble( TILE_nlon * TILE_num_lon )
      end if  
    end do

    return
  end subroutine TOPO_info_list_setup

  subroutine get_neighbor_topo_data( &
    TILE_nlon, TILE_nlat, in_topo_filebase, file_i, file_j, &
    topo )
    implicit none
    integer, intent(in) :: TILE_nlon
    integer, intent(in) :: TILE_nlat
    character(len=*), intent(in) :: in_topo_filebase
    integer, intent(in) :: file_i
    integer, intent(in) :: file_j
    real(RP), intent(out) :: topo(0:TILE_nlon+1,0:TILE_nlat+1)

    type(TOPO_info) :: dummy
    !---------------------------

    call dummy%Init(TILE_nlon, TILE_nlat)
    call dummy%Read_data(in_topo_filebase, file_i, file_j)
    topo(:,:) = dummy%topo(:,:)
    call dummy%Final()

    return
  end subroutine get_neighbor_topo_data

  subroutine TOPO_info_list_final()
    implicit none

    integer :: file_i, file_j
    integer :: domid
    !-----------------------------------------

    do domid=1, size(node_map_info_list)
      call node_map_info_list(domid)%Final()
    end do

    do file_j=1, tile_num_lat
    do file_i=1, tile_num_lon
      if ( TOPO_info_flag(file_i,file_j) ) then
        call TOPO_info_list(file_i,file_j)%Final()
      end if
    end do
    end do

    deallocate( TOPO_info_list, TOPO_info_flag )
    deallocate( node_map_info_list )

    return
  end subroutine TOPO_info_list_final

  subroutine TOPO_info_Init(this, nlon, nlat)
    implicit none
    class(TOPO_info), intent(inout) :: this
    integer, intent(in) :: nlon, nlat
    !-----------------------------------------

    this%nlon = nlon
    this%nlat = nlat
    allocate( this%lon(0:nlon+1) )
    allocate( this%lat(0:nlat+1) )
    allocate( this%topo(0:nlon+1,0:nlat+1) )

    return
  end subroutine TOPO_info_Init

  subroutine TOPO_info_Final(this)
    implicit none
    class(TOPO_info), intent(inout) :: this
    !-----------------------------------------

    deallocate( this%lon, this%lat, this%topo )

    return
  end subroutine TOPO_info_Final

  subroutine TOPO_info_read_data(this, in_topo_filebase, file_i, file_j)
    implicit none
    class(TOPO_info), intent(inout) :: this
    character(len=*), intent(in) :: in_topo_filebase    
    integer, intent(in) :: file_i, file_j

    integer :: ncid

    integer :: dimid
    integer :: varid

    real(RP) :: lon(this%nlon)
    real(RP) :: lat(this%nlat)
    real(RP) :: topo(this%nlon,this%nlat)

    real(RP) :: fac_deg2rad

    character(len=H_MID) :: nc_path
    !-----------------------------------------
    
    fac_deg2rad = PI / 180.E0_RP

    write(nc_path,'(a,i0,a,i0,a)') trim(in_topo_filebase)//"_", file_i-1, "_", file_j-1, ".nc"
    LOG_INFO("TOPO_info_read_data",*) "Read TOPO data :", trim(nc_path)

    call check( nf90_open(nc_path, nf90_nowrite, ncid) )

    !-
    call check ( nf90_inq_varid( ncid, "lon", varid) )
    call check( nf90_get_var(ncid, varid, lon) )

    this%lon(1:this%nlon) = lon(:) * fac_deg2rad + PI
    this%lon(0) = 2.0_RP * this%lon(1) - this%lon(2)
    this%lon(this%nlon+1) = 2.0_RP * this%lon(this%nlon) - this%lon(this%nlon-1) 

    !-
    call check ( nf90_inq_varid( ncid, "lat", varid) )
    call check( nf90_get_var(ncid, varid, lat) )

    this%lat(1:this%nlat) = lat(:) * fac_deg2rad
    this%lat(0) = max(2.0_RP * this%lat(1) - this%lat(2), - 0.5_RP * PI)
    this%lat(this%nlat+1) = max(2.0_RP * this%lat(this%nlat) - this%lat(this%nlat-1), + 0.5_RP * PI)

    !-
    call check ( nf90_inq_varid( ncid, "topo", varid) )
    call check( nf90_get_var(ncid, varid, topo) )
    this%topo(:,:) = 0.0_RP
    this%topo(1:this%nlon,1:this%nlat) = topo(:,:)

    LOG_INFO("TOPO_info_read_data",*) "lon range:", this%lon(1) / fac_deg2rad, ":", this%lon(this%nlon) / fac_deg2rad
    LOG_INFO("TOPO_info_read_data",*) "lat range:", this%lat(1) / fac_deg2rad, ":", this%lat(this%nlat) / fac_deg2rad

    call check( nf90_close(ncid) )

    return
  end subroutine TOPO_info_read_data

!--

  subroutine check(status)
    implicit none
    integer, intent (in) :: status

    !-----------------------------------------    
    if(status /= nf90_noerr) then 
      LOG_INFO("TOPO_INFO",*) "NC CHECK ERROR:", trim(nf90_strerror(status))
      call PRC_abort
    end if

    return
   end subroutine check  
end module mod_topo_info