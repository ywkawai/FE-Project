#include "scalelib.h"
module mod_topo_nodemap
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
    PRC_abort
  use scale_const, only: &
    PI => CONST_PI

  use scale_element_base, only: ElementBase2D
  use scale_localmesh_2d, only: LocalMesh2D

  implicit none
  private
  
  type, public :: TOPO_nodemap_info
    integer, allocatable :: file_i_ind(:,:)
    integer, allocatable :: file_j_ind(:,:)
    integer, allocatable :: i_ind(:,:)
    integer, allocatable :: j_ind(:,:)
  contains
    procedure :: Init => TOPO_nodemap_info_init
    procedure :: Final => TOPO_nodemap_info_final
    procedure :: Set_file_ijindx => TOPO_nodemap_set_file_ijindx
  end type

contains

  subroutine TOPO_nodemap_info_init(this, lcmesh, elem)
    implicit none
    class(TOPO_nodemap_info), intent(inout) :: this
    class(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem

    integer :: ke
    !-----------------------------------------    

    allocate( this%file_i_ind(elem%Np,lcmesh%Ne) )
    allocate( this%file_j_ind(elem%Np,lcmesh%Ne) )

    allocate( this%i_ind(elem%Np,lcmesh%Ne) )
    allocate( this%j_ind(elem%Np,lcmesh%Ne) )

    do ke=lcmesh%NeS, lcmesh%NeE
      this%file_i_ind(:,ke) = -1
      this%file_j_ind(:,ke) = -1
      this%i_ind(:,ke) = -1
      this%j_ind(:,ke) = -1
    end do

    return
  end subroutine TOPO_nodemap_info_init

  subroutine TOPO_nodemap_info_final(this)
    implicit none
    class(TOPO_nodemap_info), intent(inout) :: this
    !-----------------------------------------    

    deallocate(this%file_i_ind, this%file_j_ind)
    deallocate(this%i_ind, this%j_ind)

    return
  end subroutine TOPO_nodemap_info_final

  subroutine TOPO_nodemap_set_file_ijindx( this, &
    lcmesh2D, elem2D, tile_num_lon, tile_num_lat, &
    topo_info_flag )

    implicit none
    class(TOPO_nodemap_info), intent(inout) :: this
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    integer, intent(in) :: tile_num_lon
    integer, intent(in) :: tile_num_lat
    logical, intent(inout) :: topo_info_flag(tile_num_lon,tile_num_lat)

    integer :: ke, p
    integer :: file_i, file_j

    real(RP) :: DLON_per_TILE, DLAT_per_TILE
    real(RP) :: fac_rad2deg
    
    !-----------------------------

    fac_rad2deg = 180.0E0_RP / PI

    DLON_per_TILE = 360.0_RP / float(tile_num_lon)    
    DLAT_per_TILE = 180.0_RP / float(tile_num_lat)

    do ke=lcmesh2D%NeS, lcmesh2D%NeE
    do p=1, elem2D%Np
        file_i = min( int( lcmesh2D%lon(p,ke) * fac_rad2deg / DLON_per_TILE ) + 1, tile_num_lon )
        file_j = min( int( ( lcmesh2D%lat(p,ke) + 0.5_RP * PI ) * fac_rad2deg / DLAT_per_TILE ) + 1, tile_num_lat )
        topo_info_flag(file_i,file_j) = .true.  
        
        this%file_i_ind(p,ke) = file_i
        this%file_j_ind(p,ke) = file_j
    end do
    end do

    return
  end subroutine TOPO_nodemap_set_file_ijindx
end module mod_topo_nodemap
