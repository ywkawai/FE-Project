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
  use scale_monitor, only: MONITOR_set_dim
  use scale_monitor, only: MONITOR_put

  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_meshfield_base, only: MeshField2D
  use scale_meshfield_base, only: MeshField3D

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

  interface FILE_monitor_meshfield_set_dim
    module procedure FILE_monitor_meshfield_set_dim2D
    module procedure FILE_monitor_meshfield_set_dim3D
  end interface
  public :: FILE_monitor_meshfield_set_dim

  public :: FILE_monitor_meshfield_reg

  interface FILE_monitor_meshfield_put
    module procedure FILE_monitor_meshfield_put2D
    module procedure FILE_monitor_meshfield_put3D
  end interface
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
  subroutine FILE_monitor_meshfield_set_dim2D( mesh2D, dim_type )
    implicit none

    class(MeshBase2D), intent(in) :: mesh2D
    character(*), intent(in) :: dim_type
    
    real(RP) :: area(1,1)
    real(RP) :: total_area, total_area_lc

    integer :: n
    integer :: ke
    !---------------------------------------------------------------------------

    total_area = 0.0_RP
    do n=1, mesh2D%LOCAL_MESH_NUM
      total_area = total_area + cal_total_lc( mesh2D%lcmesh_list(n) )
    end do
    area(1,1) = 1.0_RP

    call MONITOR_set_dim( 1, 1, 1, 1, 1, 1, 1, 1, 1,     &
      dim_type, 2, area=area(:,:), total_area=total_area )
    
    return
  end subroutine FILE_monitor_meshfield_set_dim2D

  subroutine FILE_monitor_meshfield_set_dim3D( mesh3D, dim_type )
    implicit none

    class(MeshBase3D), intent(in) :: mesh3D
    character(*), intent(in) :: dim_type
    
    real(RP) :: area(1,1)
    real(RP) :: total_area, total_area_lc

    integer :: n
    integer :: ke
    !---------------------------------------------------------------------------

    total_area = 0.0_RP
    do n=1, mesh3D%LOCAL_MESH_NUM
      total_area = total_area + cal_total_lc( mesh3D%lcmesh_list(n) )
    end do
    area(1,1) = 1.0_RP

    call MONITOR_set_dim( 1, 1, 1, 1, 1, 1, 1, 1, 1,     &
      dim_type, 2, area=area(:,:), total_area=total_area )
    
    return
  end subroutine FILE_monitor_meshfield_set_dim3D

  subroutine FILE_monitor_meshfield_put2D( itemid, field )
    implicit none

    class(MeshField2D), intent(in) :: field
    integer, intent(in) :: itemid

    real(RP) :: total_lc(1,1)
    integer :: n
    !---------------------------------------------------------------------------

    if ( itemid <= 0 ) return

    total_lc(1,1) = 0.0_RP
    do n=1, field%mesh%LOCAL_MESH_NUM
      total_lc(1,1) = total_lc(1,1) &
                    + cal_total_lc( field%mesh%lcmesh_list(n), field%local(n)%val )
    end do
    call MONITOR_put( itemid, total_lc(:,:) )

    return
  end subroutine FILE_monitor_meshfield_put2D

  subroutine FILE_monitor_meshfield_put3D( itemid, field )
    implicit none

    class(MeshField3D), intent(in) :: field
    integer, intent(in) :: itemid

    real(RP) :: total_lc(1,1)
    integer :: n
    !---------------------------------------------------------------------------

    if ( itemid <= 0 ) return

    total_lc(1,1) = 0.0_RP
    do n=1, field%mesh%LOCAL_MESH_NUM
      total_lc(1,1) = total_lc(1,1) &
                    + cal_total_lc( field%mesh%lcmesh_list(n), field%local(n)%val )
    end do
    call MONITOR_put( itemid, total_lc(:,:) )

    return
  end subroutine FILE_monitor_meshfield_put3D

!-- private----------------------

  function cal_total_lc( lcmesh, field_val ) result(total)
    use scale_localmesh_base, only: LocalMeshBase
    use scale_element_base, only: ElementBase
    implicit none

    class(LocalMeshBase), intent(in) :: lcmesh
    real(RP), intent(in), optional :: field_val(lcmesh%refElem%Np,lcmesh%NeA)
    real(RP) :: total

    class(ElementBase), pointer :: elem
    integer :: ke
    !---------------------------------------------------------------------------

    total = 0.0_RP
    elem => lcmesh%refElem

    if ( present(field_val) ) then
      !$omp parallel do reduction(+:total)
      do ke=lcmesh%NeS, lcmesh%NeE
        total = total &
              + sum( elem%IntWeight_lgl(:) * lcmesh%J(:,ke) * lcmesh%Gsqrt(:,ke) * field_val(:,ke) )
      end do
    else
      !$omp parallel do reduction(+:total)
      do ke=lcmesh%NeS, lcmesh%NeE
        total = total + sum( elem%IntWeight_lgl(:) * lcmesh%J(:,ke) * lcmesh%Gsqrt(:,ke) )
      end do
    end if

    return
  end function cal_total_lc

end module scale_file_monitor_meshfield
