!-------------------------------------------------------------------------------
!> module MONITOR
!!
!! @par Description
!!          Monitor output module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_file_monitor_meshfield
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use scale_meshfield_base, only: MeshField3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase3D

  use scale_monitor, only: &
    FILE_monitor_meshfield_setup => MONITOR_setup,   &
    FILE_monitor_meshfield_reg   => MONITOR_reg,     &
    FILE_monitor_meshfield_write => MONITOR_write,   &
    FILE_monitor_meshfield_final => MONITOR_finalize


  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  !-----------------------------------------------------------------------------
  public :: FILE_monitor_meshfield_setup
  public :: FILE_monitor_meshfield_set_dim
  public :: FILE_monitor_meshfield_reg
  public :: FILE_monitor_meshfield_put
  public :: FILE_monitor_meshfield_write
  public :: FILE_monitor_meshfield_final


  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !

contains
  subroutine FILE_monitor_meshfield_set_dim( mesh3D, dim_type )
    use scale_monitor, only: MONITOR_set_dim
    implicit none

    class(MeshBase3D), intent(in) :: mesh3D
    character(*), intent(in) :: dim_type
    
    real(RP) :: area(1,1)
    real(RP) :: total_area, total_area_lc

    integer :: n
    integer :: ke
    type(LocalMesh3D), pointer :: lcmesh
    type(ElementBase3D), pointer :: elem3D
    !---------------------------------------------------------------------------

    total_area = cal_total_lc( mesh3D )
    area(1,1) = 1.0_RP

    call MONITOR_set_dim( 1, 1, 1, 1, 1, 1, 1, 1, 1,     &
      dim_type, 2, area=area(:,:), total_area=total_area )
    
    return
  end subroutine FILE_monitor_meshfield_set_dim

  subroutine FILE_monitor_meshfield_put( itemid, field )
    use scale_monitor, only: MONITOR_put
    implicit none

    class(MeshField3D), intent(in) :: field
    integer, intent(in) :: itemid

    real(RP) :: total_lc(1,1)
    !---------------------------------------------------------------------------

    if ( itemid <= 0 ) return

    total_lc(1,1) = cal_total_lc( field%mesh, field )
    call MONITOR_put( itemid, total_lc(:,:) )

    return
  end subroutine FILE_monitor_meshfield_put

!-- private----------------------

  function cal_total_lc( mesh3D, field ) result(total)
    implicit none

    class(MeshBase3D), intent(in), target :: mesh3D
    class(MeshField3D), intent(in), optional :: field
    real(RP) :: total

    integer :: n
    integer :: ke
    class(LocalMesh3D), pointer :: lcmesh
    class(ElementBase3D), pointer :: elem3D

    real(RP) :: total_lc
    !---------------------------------------------------------------------------

    total = 0.0_RP
    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      elem3D => lcmesh%refElem3D

      total_lc = 0.0_RP
      if ( present(field) ) then
        !$omp parallel do reduction(+:total_lc)
        do ke=lcmesh%NeS, lcmesh%NeE
          total_lc = total_lc &
            + sum( elem3D%IntWeight_lgl(:) * lcmesh%J(:,ke) * field%local(n)%val(:,ke) )
        end do
      else
        !$omp parallel do reduction(+:total_lc)
        do ke=lcmesh%NeS, lcmesh%NeE
          total_lc = total_lc + sum( elem3D%IntWeight_lgl(:) * lcmesh%J(:,ke) )
        end do
      end if
      total = total + total_lc
    end do

    return
  end function cal_total_lc

end module scale_file_monitor_meshfield
