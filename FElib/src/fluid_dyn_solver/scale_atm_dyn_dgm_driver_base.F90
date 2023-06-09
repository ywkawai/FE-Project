!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / driver (base)
!!
!! @par Description
!!      Driver module for dynamical core based on DGM 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_driver_base
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_timeint_rk, only: TimeInt_RK
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public :: AtmDynDGMDriver_base
    integer :: PROGVAR_NUM
    type(TimeInt_RK), allocatable :: tint(:)
  end type AtmDynDGMDriver_base

  type, extends(AtmDynDGMDriver_base), public :: AtmDynDGMDriver_base3D
  end type AtmDynDGMDriver_base3D

  public :: AtmDynDGMDriver_base3D_Init
  public :: AtmDynDGMDriver_base3D_Final

  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  private :: AtmDynDGMDriver_base_Init
  private :: AtmDynDGMDriver_base_Final

contains
  subroutine AtmDynDGMDriver_base_Init( this, &
    PRGVAR_NUM, LOCAL_MESH_NUM )
    implicit none

    class(AtmDynDGMDriver_base), intent(inout) :: this
    integer, intent(in) :: PRGVAR_NUM
    integer, intent(in) :: LOCAL_MESH_NUM
    !-----------------------------------------------------------------------------

    this%PROGVAR_NUM = PRGVAR_NUM
    allocate( this%tint(LOCAL_MESH_NUM) )

    return
  end subroutine AtmDynDGMDriver_base_Init

  subroutine AtmDynDGMDriver_base_Final( this  )
    implicit none

    class(AtmDynDGMDriver_base), intent(inout) :: this
    !-----------------------------------------------------------------------------

    deallocate( this%tint )

    return
  end subroutine AtmDynDGMDriver_base_Final

!-- 3D

  subroutine AtmDynDGMDriver_base3D_Init( this, &
    prgvar_num,                                 &
    tint_type, dtsec,                           &
    mesh )
    implicit none

    class(AtmDynDGMDriver_base3D), intent(inout) :: this
    integer, intent(in) :: prgvar_num 
    character(len=*), intent(in) :: tint_type
    real(DP), intent(in) :: dtsec
    class(MeshBase3D), intent(in), target :: mesh

    class(LocalMesh3D), pointer :: lcmesh
    integer :: n
    !-----------------------------------------------------------------------------

    call AtmDynDGMDriver_base_Init( this, prgvar_num, mesh%LOCAL_MESH_NUM )

    do n = 1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)

      call this%tint(n)%Init( tint_type, dtsec, prgvar_num, 2, &
        (/ lcmesh%refElem%Np, lcmesh%NeA /) )
    end do

    return
  end subroutine AtmDynDGMDriver_base3D_Init

  subroutine AtmDynDGMDriver_base3D_Final( this )
    implicit none

    class(AtmDynDGMDriver_base3D), intent(inout) :: this

    integer :: n
    !-----------------------------------------------------------------------------

    do n = 1, size(this%tint)
      call this%tint(n)%Final()
    end do
    call AtmDynDGMDriver_base_Final( this )   

    return
  end subroutine AtmDynDGMDriver_base3D_Final

end module scale_atm_dyn_dgm_driver_base