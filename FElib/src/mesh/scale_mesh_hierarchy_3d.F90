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
    MESH_HIERARCHY_hMG_FINEST_LEVEL, &
    MESH_HIERARCHY_TYPE_pMG,         &
    MESH_HIERARCHY_TYPE_hMG

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
    type(MeshPtr3D), pointer :: fine_mesh => null()
    type(MeshPtr3D), pointer :: coarse_mesh => null()
    type(MeshHierarchyLocalMGData3D), allocatable :: mg_local(:)

    real(RP), allocatable :: pMat1D_f2c(:,:)
    real(RP), allocatable :: pMat1D_c2f(:,:)
  contains
    procedure :: Init => MeshHierarchyLevel3D_Init
    procedure :: Final => MeshHierarchyLevel3D_Final
  end type MeshHierarchyLevel3D

  type, extends(MeshHierarchyBase), public :: MeshHierarchy3D
    type(MeshHierarchyLevel3D), allocatable :: p_level(:)
    type(MeshHierarchyLevel3D), allocatable :: h_level(:)

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
    parent_mesh,                                  &
    porder_list, p_LEVEL_NUM,                     &
    NeGX_list, NeGY_list, NeGZ_list, h_LEVEL_NUM  )
    implicit none
    class(MeshHierarchy3D), intent(inout), target :: this
    class(MeshBase3D), intent(in), target :: parent_mesh
    integer, intent(in) :: p_LEVEL_NUM
    integer, intent(in) :: porder_list(p_LEVEL_NUM)
    integer, intent(in) :: h_LEVEL_NUM
    integer, intent(in) :: NeGX_list(h_LEVEL_NUM)
    integer, intent(in) :: NeGY_list(h_LEVEL_NUM)
    integer, intent(in) :: NeGZ_list(h_LEVEL_NUM)

    integer :: poly_lev
    integer :: h_lev

    integer :: ldom_id

    type(MeshHierarchyLevel3D), pointer :: level_ptr
    class(MeshBase3D), pointer :: mesh3D_ptr

    class(MeshPtr3D), pointer :: fine_mesh_ptr
    class(MeshPtr3D), pointer :: coarse_mesh_ptr

    integer :: NeGZ
    !-------------------------------------------------------------

    call MeshHierarchyBase_Init( this, p_LEVEL_NUM, h_LEVEL_NUM )

    !--
    if ( parent_mesh%refElem3D%PolyOrder_h /= parent_mesh%refElem3D%PolyOrder_v ) then
      LOG_INFO('MeshHierarchy3D_Init',*) 'Currently, PolyOrder_h should equal PolyOrder_v. Check!'
      call PRC_abort
    end if
    if ( porder_list(p_LEVEL_NUM) /= 1 .and. h_LEVEL_NUM > 0 ) then
      LOG_INFO('MeshHierarchy3D_Init',*) 'Currently, only p=1 is supported for the coarsest p-mesh level. Check!'
      call PRC_abort
    end if
    if ( h_LEVEL_NUM > 0 ) then
      NeGZ = NeGZ_list(1)
      do h_lev=2, h_LEVEL_NUM
        if ( NeGZ_list(h_lev) /= NeGZ ) then
          LOG_INFO('MeshHierarchy3D_Init',*) 'Currently, the number of elements should be same in the vertical direction. Check!'
          call PRC_abort
        end if
      end do
    end if

    !- Setup p-mesh hierarchy

    allocate( this%elem3D_list(p_LEVEL_NUM) )
    do poly_lev=1, p_LEVEL_NUM
      call this%elem3D_list(poly_lev)%Init( porder_list(poly_lev), porder_list(poly_lev), .false. )
    end do

    allocate( this%p_mesh_list(this%NUM_pMG_LEVEL) )

    this%p_mesh_list(MESH_HIERARCHY_pMG_FINEST_LEVEL)%ptr => parent_mesh
    do poly_lev=MESH_HIERARCHY_pMG_FINEST_LEVEL + 1, this%NUM_pMG_LEVEL
      select type(parent_mesh)
      type is (MeshCubeDom3D)
        call construct_cubedom3D_mesh( this, this%p_mesh_list(poly_lev)%ptr,    &
          parent_mesh%NeGX, parent_mesh%NeGY, parent_mesh%NeGZ, parent_mesh%FZ, &
          parent_mesh, this%elem3D_list(poly_lev) )
      end select
    end do

    allocate( this%p_level(this%NUM_pMG_LEVEL) )

    do poly_lev=MESH_HIERARCHY_pMG_FINEST_LEVEL, this%NUM_pMG_LEVEL
      call this%p_level(poly_lev)%Init( poly_lev, this%p_mesh_list(poly_lev)%ptr, &
        this%p_mesh_list, this%NUM_pMG_LEVEL, MESH_HIERARCHY_TYPE_pMG )
    end do

    !- Setup h-mesh hierarchy

    allocate( this%h_mesh_list(this%NUM_hMG_LEVEL) )

    this%h_mesh_list(MESH_HIERARCHY_hMG_FINEST_LEVEL)%ptr => parent_mesh
    do h_lev=MESH_HIERARCHY_hMG_FINEST_LEVEL + 1, this%NUM_hMG_LEVEL
      select type(parent_mesh)
      type is (MeshCubeDom3D)
        call construct_cubedom3D_mesh( this, this%h_mesh_list(h_lev)%ptr,       &
          NeGX_list(h_lev), NeGY_list(h_lev), NeGZ_list(h_lev), parent_mesh%FZ, &
          parent_mesh, this%elem3D_list(p_LEVEL_NUM) )
      end select
    end do

    allocate( this%h_level(this%NUM_hMG_LEVEL) )

    do h_lev=MESH_HIERARCHY_hMG_FINEST_LEVEL, this%NUM_hMG_LEVEL
      call this%h_level(h_lev)%Init( h_lev, this%h_mesh_list(h_lev)%ptr, &
        this%h_mesh_list, this%NUM_hMG_LEVEL, MESH_HIERARCHY_TYPE_hMG )
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
    do poly_lev=MESH_HIERARCHY_pMG_FINEST_LEVEL, this%NUM_pMG_LEVEL
      call this%p_level(poly_lev)%Final()
    end do
    deallocate( this%p_level )

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
      call this%h_level(h_lev)%Final()
    end do
    deallocate( this%h_level )

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
    level_id, mesh3D, mesh_list, LEVEL_NUM,   &
    hierarchy_type )
    use scale_mesh_hierarchy_base, only: &
      MeshHierarchy_construct_pMG_mat1D
    implicit none
    class(MeshHierarchyLevel3D), intent(inout) :: this
    integer, intent(in) :: level_id
    class(MeshBase3D), intent(in) :: mesh3D
    integer, intent(in) :: LEVEL_NUM
    class(MeshPtr3D), intent(in), target :: mesh_list(LEVEL_NUM)
    integer, intent(in) :: hierarchy_type

    integer :: ldom_id
    !-------------------------------------------------------------

    this%hierarchy_type = hierarchy_type
    this%level_id = level_id

    !- set mesh pointers with finer/coarser meshes

    if ( level_id > 1 ) then
      this%fine_mesh => mesh_list(level_id-1)
    else
      this%fine_mesh => null()
    end if

    if ( level_id < LEVEL_NUM ) then
      this%coarse_mesh => mesh_list(level_id+1)
    else
      this%coarse_mesh => null()
    end if

    !- p-hierarchy
    if ( hierarchy_type == MESH_HIERARCHY_TYPE_pMG &
         .and. level_id < LEVEL_NUM                ) then
      
      allocate( this%pMat1D_c2f(mesh3D%refElem3D%Nnode_h1D,this%coarse_mesh%ptr%refElem3D%Nnode_h1D) )
      call MeshHierarchy_construct_pMG_mat1D( this%pMat1D_c2f,   &
        this%coarse_mesh%ptr%refElem3D%Nnode_h1D, mesh3D%refElem3D%Nnode_h1D )

      allocate( this%pMat1D_f2c(this%coarse_mesh%ptr%refElem3D%Nnode_h1D,mesh3D%refElem3D%Nnode_h1D) )
      call MeshHierarchy_construct_pMG_mat1D( this%pMat1D_f2c,   &
        mesh3D%refElem3D%Nnode_h1D, this%coarse_mesh%ptr%refElem3D%Nnode_h1D )
    end if

    !- h-hierarchy
    if ( hierarchy_type == MESH_HIERARCHY_TYPE_hMG &
         .and. level_id < LEVEL_NUM                ) then

      allocate( this%mg_local( mesh3D%LOCAL_MESH_NUM ) )
      
      do ldom_id=1, mesh3D%LOCAL_MESH_NUM
        call this%mg_local(ldom_id)%Init( mesh_list(level_id)%ptr%lcmesh_list(ldom_id), &
          this%coarse_mesh%ptr%lcmesh_list, mesh_list(level_id)%ptr%refElem3D           )
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

    integer :: i, j, k
    integer :: ke

    integer :: i_c, j_c, k_c
    integer :: ke_c

    integer :: ke2i(lcmesh3D%Ne)
    integer :: ke2j(lcmesh3D%Ne)
    integer :: ke2k(lcmesh3D%Ne)

    integer :: lcdomID_c
    integer :: lcTileID_c
    class(LocalMesh3D), pointer :: lcmesh_c

    real(RP) :: vx_c(elem3D%Nv)
    real(RP) :: vy_c(elem3D%Nv)
    real(RP) :: vz_c(elem3D%Nv)
    integer :: i_EtoV(elem3D%Nv)

    integer :: p
    real(RP) :: r_c, s_c, t_c

    integer :: l
    !-------------------------------------------------------------  

    call MeshHierarchyLocalMGDataBase_Init( this, &
      lcmesh3D )

    ! Current implementation assumes identical tileID & lcdomID in both coarsened meshes          
    this%CoarseLocalMesh_tileIDlist(1) = lcmesh3D%tileID
    this%CoarseLocalMesh_lcdomIDlist(1) = 1

    ! Preparation
    do k=1, lcmesh3D%NeZ
    do j=1, lcmesh3D%NeY
    do i=1, lcmesh3D%NeX
      ke = i + (j-1)*lcmesh3D%NeX + (k-1)*lcmesh3D%NeX*lcmesh3D%NeY
      ke2i(ke) = i; ke2j(ke) = j; ke2k(ke) = k
    end do
    end do
    end do

    ! Set the relation between coarse and fine grid,
    ! and construct interpolation operator

    allocate( this%If2c_emap(4,lcmesh3D%Ne/4) )

    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      lcTileID_c = this%CoarseLocalMesh_tileIDlist(1)
      lcdomID_c = this%CoarseLocalMesh_lcdomIDlist(1)

      lcmesh_c => coarse_lcmesh_list(lcdomID_c)

      i_c = (ke2i(ke)+1)/2; j_c = (ke2j(ke)+1)/2; k_c = ke2k(ke)
      ke_c = i_c + (j_c-1)*lcmesh_c%NeX + (k_c-1)*lcmesh_c%NeX*lcmesh_c%NeY
      this%Ic2f_emap(ke) = ke_c

      i = ke2i(ke) - 2*(i_c-1)
      j = ke2j(ke) - 2*(j_c-1)
!      k = ke2k(ke) - 2*(k_c-1)
      this%If2c_emap(i+(j-1)*2,ke_c) = ke

      i_EtoV(:) = lcmesh_c%EToV(ke_c,:)
      vx_c(:) = lcmesh_c%pos_ev(i_EtoV(:),1)
      vy_c(:) = lcmesh_c%pos_ev(i_EtoV(:),2)
      vz_c(:) = lcmesh_c%pos_ev(i_EtoV(:),3)
      do p=1, elem3D%Np
        r_c = -1.0_RP + 2.0_RP * ( lcmesh3D%pos_en(p,ke,1) - vx_c(1) ) / ( vx_c(2) - vx_c(1) )
        s_c = -1.0_RP + 2.0_RP * ( lcmesh3D%pos_en(p,ke,2) - vy_c(1) ) / ( vy_c(4) - vy_c(1) )
        t_c = -1.0_RP + 2.0_RP * ( lcmesh3D%pos_en(p,ke,3) - vz_c(1) ) / ( vz_c(5) - vz_c(1) )

        this%Ic2f(p,:,ke) = 0.125_RP * &
            (/ ( 1.0_RP - r_c ) * ( 1.0_RP - s_c ) * ( 1.0_RP - t_c ), ( 1.0_RP + r_c ) * ( 1.0_RP - s_c ) * ( 1.0_RP - t_c ),  &
               ( 1.0_RP - r_c ) * ( 1.0_RP + s_c ) * ( 1.0_RP - t_c ), ( 1.0_RP + r_c ) * ( 1.0_RP + s_c ) * ( 1.0_RP - t_c ),  &
               ( 1.0_RP - r_c ) * ( 1.0_RP - s_c ) * ( 1.0_RP + t_c ), ( 1.0_RP + r_c ) * ( 1.0_RP - s_c ) * ( 1.0_RP + t_c ),  &
               ( 1.0_RP - r_c ) * ( 1.0_RP + s_c ) * ( 1.0_RP + t_c ), ( 1.0_RP + r_c ) * ( 1.0_RP + s_c ) * ( 1.0_RP + t_c )  /)
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
    NeGX, NeGY, NeGZ, FZ, &
    parent_mesh, elem3D   )
    implicit none
    class(MeshHierarchy3D), intent(inout) :: this
    class(MeshBase3D), intent(out), pointer :: child_mesh_base_ptr
    integer, intent(in) :: NeGX
    integer, intent(in) :: NeGY
    integer, intent(in) :: NeGZ
    real(RP), intent(in) :: FZ(NeGZ+1)
    type(MeshCubeDom3D), intent(in) :: parent_mesh
    type(HexahedralElement), intent(in) :: elem3D

    type(MeshCubeDom3D), pointer :: child_mesh_ptr
    !-------------------------------------------------------------

    write(*,*) "Fz=", FZ(:)
    write(*,*) parent_mesh%isPeriodicX, parent_mesh%isPeriodicY, parent_mesh%isPeriodicZ

    allocate( child_mesh_ptr )
    call child_mesh_ptr%Init( NeGX, NeGY, NeGZ,                                            &
      parent_mesh%xmin_gl, parent_mesh%xmax_gl, parent_mesh%ymin_gl, parent_mesh%ymax_gl,  &
      parent_mesh%zmin_gl, parent_mesh%zmax_gl,                                            &
      parent_mesh%isPeriodicX, parent_mesh%isPeriodicY, parent_mesh%isPeriodicZ,           &
      elem3D, 1, parent_mesh%NprcX, parent_mesh%NprcY,                                     &
      FZ=FZ )

    call child_mesh_ptr%Generate()

    child_mesh_base_ptr => child_mesh_ptr
    return
  end subroutine construct_cubedom3D_mesh
end module scale_mesh_hierarchy_3d
