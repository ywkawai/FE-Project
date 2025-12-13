!-------------------------------------------------------------------------------
!> module FElib / Mesh / 3D domain
!!
!! @par Description
!!      Manage mesh hierarchy of 3D domain for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_mesh_hierarchy_3d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement

  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D

  use scale_mesh_hierarchy_base, only: &
    MeshHierarchyBase_Init, &
    MeshHierarchyBase_Final, &
    MeshHierarchyBase,       &
    MeshHierarchyLevelBase,  &
    MeshHierarchyLocalMGDataBase, &
    MESH_HIERARCHY_pMG_FINEST_LEVEL, &
    MESH_HIERARCHY_hMG_FINEST_LEVEL

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type :: MeshPtr3D
    class(MeshBase3D), pointer :: ptr => null()
  end type MeshPtr3D

  type, extends(MeshHierarchyLocalMGDataBase), public :: MeshHierarchyLocalMGData3D
  contains
    procedure :: Init => MeshHierarchyLocalMGData3D_Init
    procedure :: Final => MeshHierarchyLocalMGData3D_Final
  end type MeshHierarchyLocalMGData3D

  type, extends(MeshHierarchyLevelBase), public :: MeshHierarchyLevel3D
    logical :: is_h_level = .true.
    type(MeshPtr3D), pointer :: fine_mesh => null()
    type(MeshPtr3D), pointer :: coarse_mesh => null()
    type(MeshHierarchyLocalMGData3D), allocatable :: mg_local(:)
  contains
    procedure :: Init => MeshHierarchyLevel3D_Init
    procedure :: Final => MeshHierarchyLevel3D_Final
  end type MeshHierarchyLevel3D

  type, extends(MeshHierarchyBase), public :: MeshHierarchy3D
    type(MeshHierarchyLevel3D), allocatable :: level(:)

    type(MeshPtr3D), allocatable :: p_mesh_list(:)
    type(MeshPtr3D), allocatable :: h_mesh_list(:)

    type(HexahedralElement), allocatable :: elem3D_list(:)
  contains
    procedure :: Init => MeshHierarchy3D_Init
    procedure :: Final => MeshHierarchy3D_Final
  end type MeshHierarchy3D

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
!OCL SERIAL
  subroutine MeshHierarchy3D_Init( this, &
    parent_mesh,                       &
    porder_list, p_LEVEL_NUM,          &
    NeGX_list, NeGY_list, h_LEVEL_NUM  )
    implicit none
    class(MeshHierarchy3D), intent(inout), target :: this
    class(MeshBase3D), intent(in), target :: parent_mesh
    integer, intent(in) :: p_LEVEL_NUM
    integer, intent(in) :: porder_list(p_LEVEL_NUM)
    integer, intent(in) :: h_LEVEL_NUM
    integer, intent(in) :: NeGX_list(h_LEVEL_NUM)
    integer, intent(in) :: NeGY_list(h_LEVEL_NUM)

    integer :: poly_lev
    integer :: h_lev

    integer :: ldom_id

    type(MeshHierarchyLevel3D), pointer :: level_ptr
    class(MeshBase3D), pointer :: mesh3D_ptr

    class(MeshPtr3D), pointer :: fine_mesh_ptr
    class(MeshPtr3D), pointer :: coarse_mesh_ptr
    !-------------------------------------------------------------

    call MeshHierarchyBase_Init( this, p_LEVEL_NUM, h_LEVEL_NUM )

    allocate( this%p_mesh_list(this%NUM_pMG_LEVEL) )
    allocate( this%h_mesh_list(this%NUM_hMG_LEVEL) )
    allocate( this%level(this%NUM_hMG_LEVEL) )

    !--
    if ( porder_list(p_LEVEL_NUM) /= 1 ) then
      LOG_INFO('MeshHierarchy3D_Init',*) 'Currently, only p=1 is supported for the coarsest p-mesh level. Check!'
      call PRC_abort
    end if

    allocate( this%elem3D_list(p_LEVEL_NUM) )
    do poly_lev=1, p_LEVEL_NUM
      call this%elem3D_list(poly_lev)%Init( porder_list(poly_lev), porder_list(poly_lev), .false. )
    end do

    !- Setup p-mesh hierarchy

    this%p_mesh_list(MESH_HIERARCHY_pMG_FINEST_LEVEL)%ptr => parent_mesh
    do poly_lev=MESH_HIERARCHY_pMG_FINEST_LEVEL + 1, this%NUM_pMG_LEVEL
      select type(parent_mesh)
      type is (MeshCubeDom3D)
        call construct_cubedom3D_mesh( this, this%p_mesh_list(poly_lev)%ptr, &
          parent_mesh%NeGX, parent_mesh%NeGY,     &
          parent_mesh, this%elem3D_list(poly_lev) )
      end select
    end do

    !- Setup h-mesh hierarchy

    this%h_mesh_list(MESH_HIERARCHY_hMG_FINEST_LEVEL)%ptr => parent_mesh
    do h_lev=MESH_HIERARCHY_hMG_FINEST_LEVEL + 1, this%NUM_hMG_LEVEL
      select type(parent_mesh)
      type is (MeshCubeDom3D)
        call construct_cubedom3D_mesh( this, this%h_mesh_list(h_lev)%ptr, &
          NeGX_list(h_lev), NeGY_list(h_lev),        &
          parent_mesh, this%elem3D_list(p_LEVEL_NUM) )
      end select
    end do

    do h_lev=MESH_HIERARCHY_hMG_FINEST_LEVEL, this%NUM_hMG_LEVEL
      call this%level(h_lev)%Init( h_lev, this%h_mesh_list(h_lev)%ptr, &
        this%h_mesh_list, this%NUM_hMG_LEVEL )
    end do

    return
  end subroutine MeshHierarchy3D_Init

!OCL SERIAL
  subroutine MeshHierarchy3D_Final(this)
    implicit none
    class(MeshHierarchy3D), intent(inout) :: this

    class(MeshBase3D), pointer :: mesh3D_ptr
    integer :: poly_lev
    integer :: h_lev
    integer :: ldom_id
    !-------------------------------------------------------------

    ! p-mesh hierarchy finalization
    do poly_lev=MESH_HIERARCHY_pMG_FINEST_LEVEL + 1, this%NUM_pMG_LEVEL
      mesh3D_ptr => this%p_mesh_list(poly_lev)%ptr
      select type(mesh3D_ptr)
      type is (MeshCubeDom3D)
        call mesh3D_ptr%Final()
      end select
    end do
    deallocate( this%p_mesh_list )

    ! h-mesh hierarchy finalization
    do h_lev=MESH_HIERARCHY_hMG_FINEST_LEVEL, this%NUM_hMG_LEVEL
      call this%level(h_lev)%Final()
    end do
    deallocate( this%level )

    do h_lev=MESH_HIERARCHY_hMG_FINEST_LEVEL + 1, this%NUM_hMG_LEVEL
      mesh3D_ptr => this%h_mesh_list(h_lev)%ptr
      select type(mesh3D_ptr)
      type is (MeshCubeDom3D)
        call mesh3D_ptr%Final()
      end select
    end do
    deallocate( this%h_mesh_list )
    
    !-
    call MeshHierarchyBase_Final( this )
    return
  end subroutine MeshHierarchy3D_Final

!-- private --------------------------------------------------------------

!OCL SERIAL
  subroutine MeshHierarchyLevel3D_Init( this, &
    hlevel_id, mesh3D, h_mesh_list, h_LEVEL_NUM )
    implicit none
    class(MeshHierarchyLevel3D), intent(inout) :: this
    integer, intent(in) :: hlevel_id
    class(MeshBase3D), intent(in) :: mesh3D
    integer, intent(in) :: h_LEVEL_NUM
    class(MeshPtr3D), intent(in), target :: h_mesh_list(h_LEVEL_NUM)

    integer :: ldom_id
    !-------------------------------------------------------------

    this%level_id = hlevel_id

    !- set mesh pointers with finer/coarser meshes

    if ( hlevel_id > MESH_HIERARCHY_hMG_FINEST_LEVEL ) then
      this%fine_mesh => h_mesh_list(hlevel_id-1)
    else
      this%fine_mesh => null()
    end if

    if ( hlevel_id < h_LEVEL_NUM ) then
      this%coarse_mesh => h_mesh_list(hlevel_id+1)
    else
      this%coarse_mesh => null()
    end if

    !-
    if ( hlevel_id < h_LEVEL_NUM ) then
      allocate( this%mg_local( mesh3D%LOCAL_MESH_NUM ) )

      do ldom_id=1, mesh3D%LOCAL_MESH_NUM
        call this%mg_local(ldom_id)%Init( mesh3D%lcmesh_list(ldom_id), &
          this%coarse_mesh%ptr%lcmesh_list, mesh3D%refElem3D           )
      end do
    end if

    return
  end subroutine MeshHierarchyLevel3D_Init

!OCL SERIAL
  subroutine MeshHierarchyLevel3D_Final( this )
    implicit none
    class(MeshHierarchyLevel3D), intent(inout) :: this

    integer :: ldom_id
    !-------------------------------------------------------------

    if ( allocated(this%mg_local) ) then
      do ldom_id=1, size(this%mg_local)
        call this%mg_local(ldom_id)%Final()
      end do
      deallocate( this%mg_local )
    end if
    return
  end subroutine MeshHierarchyLevel3D_Final

!OCL SERIAL
  subroutine MeshHierarchyLocalMGData3D_Init( this, &
    lcmesh3D, coarse_lcmesh_list, elem3D )
    use scale_mesh_hierarchy_base, only: MeshHierarchyLocalMGDataBase_Init
    implicit none
    class(MeshHierarchyLocalMGData3D), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(LocalMesh3D), intent(in), target :: coarse_lcmesh_list(:)
    class(ElementBase3D), intent(in) :: elem3D

    integer :: i, j
    integer :: ke

    integer :: i_c, j_c
    integer :: ke_c

    integer ke2i(lcmesh3D%NeX), ke2j(lcmesh3D%NeY)

    integer :: lcTileID_c
    type(LocalMesh3D), pointer :: lcmesh_c

    real(RP) :: vx_c(elem3D%Nv)
    real(RP) :: vy_c(elem3D%Nv)
    real(RP) :: vz_c(elem3D%Nv)
    integer :: i_EtoV(elem3D%Nv)

    integer :: p
    real(RP) :: r_c, s_c

    integer :: l
    !-------------------------------------------------------------  

    call MeshHierarchyLocalMGDataBase_Init( this, &
      lcmesh3D )

    this%CoarseLocalMesh_tileIDlist(1) = lcmesh3D%tileID

    ! Preparetion
    do j=1, lcmesh3D%NeY
    do i=1, lcmesh3D%NeX
      ke = i + (j-1)*lcmesh3D%NeX
      ke2i(ke) = i; ke2j(ke) = j
    end do
    end do

    ! Set the relation between coarse and fine grid,
    ! and construct interpolation operator

    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      lcTileID_c = this%CoarseLocalMesh_tileIDlist(1)
      lcmesh_c => coarse_lcmesh_list(lcTileID_c)

      i_c = (ke2i(ke)+1)/2; j_c = (ke2j(ke)+1)/2
      ke_c = i_c + (j_c-1)*lcmesh3D%NeX/2
      this%Ic2f_emap(ke) = ke_c

      i_EtoV(:) = lcmesh_c%EToV(ke_c,:)
      vx_c(:) = lcmesh_c%pos_ev(i_EtoV(:),1)
      vy_c(:) = lcmesh_c%pos_ev(i_EtoV(:),2)
      vz_c(:) = lcmesh_c%pos_ev(i_EtoV(:),3)
      do p=1, elem3D%Np
        r_c = -1.0_RP + 2.0_RP * ( lcmesh3D%pos_en(p,ke,1) - vx_c(1) ) / ( vx_c(2) - vx_c(1) )
        s_c = -1.0_RP + 2.0_RP * ( lcmesh3D%pos_en(p,ke,2) - vy_c(1) ) / ( vy_c(4) - vy_c(1) )

        this%Ic2f(p,1:4,ke) = 0.25_RP * &
          (/ ( 1.0_RP - r_c ) * ( 1.0_RP - s_c ), ( 1.0_RP + r_c ) * ( 1.0_RP - s_c ),  &
             ( 1.0_RP - r_c ) * ( 1.0_RP + s_c ), ( 1.0_RP + r_c ) * ( 1.0_RP + s_c )  /)
      end do

    end do    
    return
  end subroutine MeshHierarchyLocalMGData3D_Init

!OCL SERIAL
  subroutine MeshHierarchyLocalMGData3D_Final( this )
    use scale_mesh_hierarchy_base, only: MeshHierarchyLocalMGDataBase_Final
    implicit none
    class(MeshHierarchyLocalMGData3D), intent(inout) :: this
    !-------------------------------------------------------------  
    call MeshHierarchyLocalMGDataBase_Final( this )
    return
  end subroutine MeshHierarchyLocalMGData3D_Final

!-
!OCL SERIAL
  subroutine construct_cubedom3D_mesh( this, child_mesh_base_ptr, &
    NeGX, NeGY,         &
    parent_mesh, elem3D )
    implicit none
    class(MeshHierarchy3D), intent(inout) :: this
    class(MeshBase3D), intent(out), pointer :: child_mesh_base_ptr
    integer, intent(in) :: NeGX
    integer, intent(in) :: NeGY
    type(MeshCubeDom3D), intent(in) :: parent_mesh
    type(HexahedralElement), intent(in) :: elem3D

    type(MeshCubeDom3D), pointer :: child_mesh_ptr
    !-------------------------------------------------------------

    allocate( child_mesh_ptr )
    call child_mesh_ptr%Init( NeGX, NeGY, parent_mesh%NeGZ,                                &
      parent_mesh%xmin_gl, parent_mesh%xmax_gl, parent_mesh%ymin_gl, parent_mesh%ymax_gl,  &
      parent_mesh%zmin_gl, parent_mesh%zmax_gl,                                            &
      parent_mesh%isPeriodicX, parent_mesh%isPeriodicY, parent_mesh%isPeriodicZ,           &
      elem3D, 1, parent_mesh%NprcX, parent_mesh%NprcY,                                     &
      FZ=parent_mesh%FZ )

    call child_mesh_ptr%Generate()

    child_mesh_base_ptr => child_mesh_ptr
    return
  end subroutine construct_cubedom3D_mesh
end module scale_mesh_hierarchy_3d
