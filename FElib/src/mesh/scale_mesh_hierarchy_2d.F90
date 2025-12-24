!-------------------------------------------------------------------------------
!> module FElib / Mesh / 2D domain
!!
!! @par Description
!!      Manage mesh hierarchy of 2D domain for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_mesh_hierarchy_2d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement

  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D

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
  type :: MeshPtr2D
    class(MeshBase2D), pointer :: ptr => null()
  end type MeshPtr2D

  type, extends(MeshHierarchyLocalMGDataBase), public :: MeshHierarchyLocalMGData2D
  contains
    procedure :: Init => MeshHierarchyLocalMGData2D_Init
    procedure :: Final => MeshHierarchyLocalMGData2D_Final
  end type MeshHierarchyLocalMGData2D

  type, extends(MeshHierarchyLevelBase), public :: MeshHierarchyLevel2D
    type(MeshPtr2D), pointer :: fine_mesh => null()
    type(MeshPtr2D), pointer :: coarse_mesh => null()
    type(MeshHierarchyLocalMGData2D), allocatable :: mg_local(:)

    real(RP), allocatable :: pMat1D_f2c(:,:)
    real(RP), allocatable :: pMat1D_c2f(:,:)
  contains
    procedure :: Init => MeshHierarchyLevel2D_Init
    procedure :: Final => MeshHierarchyLevel2D_Final
  end type MeshHierarchyLevel2D

  type, extends(MeshHierarchyBase), public :: MeshHierarchy2D
    type(MeshHierarchyLevel2D), allocatable :: p_level(:)
    type(MeshHierarchyLevel2D), allocatable :: h_level(:)

    type(MeshPtr2D), allocatable :: p_mesh_list(:)
    type(MeshPtr2D), allocatable :: h_mesh_list(:)

    type(QuadrilateralElement), allocatable :: elem2D_list(:)
  contains
    procedure :: Init => MeshHierarchy2D_Init
    procedure :: Final => MeshHierarchy2D_Final
  end type MeshHierarchy2D

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
  subroutine MeshHierarchy2D_Init( this, &
    parent_mesh,                       &
    porder_list, p_LEVEL_NUM,          &
    NeGX_list, NeGY_list, h_LEVEL_NUM  )
    implicit none
    class(MeshHierarchy2D), intent(inout), target :: this
    class(MeshBase2D), intent(in), target :: parent_mesh
    integer, intent(in) :: p_LEVEL_NUM
    integer, intent(in) :: porder_list(p_LEVEL_NUM)
    integer, intent(in) :: h_LEVEL_NUM
    integer, intent(in) :: NeGX_list(h_LEVEL_NUM)
    integer, intent(in) :: NeGY_list(h_LEVEL_NUM)

    integer :: poly_lev
    integer :: h_lev

    integer :: ldom_id

    type(MeshHierarchyLevel2D), pointer :: level_ptr
    class(MeshBase2D), pointer :: mesh2D_ptr

    class(MeshPtr2D), pointer :: fine_mesh_ptr
    class(MeshPtr2D), pointer :: coarse_mesh_ptr
    !-------------------------------------------------------------

    call MeshHierarchyBase_Init( this, p_LEVEL_NUM, h_LEVEL_NUM )

    !- Setup p-mesh hierarchy

    if ( porder_list(p_LEVEL_NUM) /= 1 .and. h_LEVEL_NUM > 0 ) then
      LOG_INFO('MeshHierarchy2D_Init',*) 'Currently, only p=1 is supported for the coarsest p-mesh level when h_LEVEL_NUM > 0. Check!'
      call PRC_abort
    end if

    !- Setup p-mesh hierarchy

    allocate( this%elem2D_list(p_LEVEL_NUM) )
    do poly_lev=1, p_LEVEL_NUM
      call this%elem2D_list(poly_lev)%Init( porder_list(poly_lev), .false. )
    end do

    allocate( this%p_mesh_list(this%NUM_pMG_LEVEL) )

    this%p_mesh_list(MESH_HIERARCHY_pMG_FINEST_LEVEL)%ptr => parent_mesh
    do poly_lev=MESH_HIERARCHY_pMG_FINEST_LEVEL + 1, this%NUM_pMG_LEVEL
      select type(parent_mesh)
      type is (MeshRectDom2D)
        call construct_rectdom2D_mesh( this, this%p_mesh_list(poly_lev)%ptr, &
          parent_mesh%NeGX, parent_mesh%NeGY,     &
          parent_mesh, this%elem2D_list(poly_lev) )
      end select
    end do

    allocate( this%p_level(this%NUM_pMG_LEVEL) )

    do poly_lev=MESH_HIERARCHY_pMG_FINEST_LEVEL, this%NUM_pMG_LEVEL
      call this%p_level(poly_lev)%Init( poly_lev, this%p_mesh_list(poly_lev)%ptr, &
        this%p_mesh_list, this%NUM_pMG_LEVEL, MESH_HIERARCHY_TYPE_pMG )
    end do

    !- Setup h-mesh hierarchy

    allocate( this%h_mesh_list(this%NUM_hMG_LEVEL) )

    this%h_mesh_list(MESH_HIERARCHY_hMG_FINEST_LEVEL)%ptr => this%p_mesh_list(this%NUM_pMG_LEVEL)%ptr
    do h_lev=MESH_HIERARCHY_hMG_FINEST_LEVEL + 1, this%NUM_hMG_LEVEL
      select type(parent_mesh)
      type is (MeshRectDom2D)
        call construct_rectdom2D_mesh( this, this%h_mesh_list(h_lev)%ptr, &
          NeGX_list(h_lev), NeGY_list(h_lev),                             &
          parent_mesh, this%elem2D_list(p_LEVEL_NUM) )
      end select
    end do

    allocate( this%h_level(this%NUM_hMG_LEVEL) )

    do h_lev=MESH_HIERARCHY_hMG_FINEST_LEVEL, this%NUM_hMG_LEVEL
      call this%h_level(h_lev)%Init( h_lev, this%h_mesh_list(h_lev)%ptr, &
        this%h_mesh_list, this%NUM_hMG_LEVEL, MESH_HIERARCHY_TYPE_hMG )
    end do

    return
  end subroutine MeshHierarchy2D_Init

!OCL SERIAL
  subroutine MeshHierarchy2D_Final(this)
    implicit none
    class(MeshHierarchy2D), intent(inout) :: this

    class(MeshBase2D), pointer :: mesh2D_ptr
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
      mesh2D_ptr => this%p_mesh_list(poly_lev)%ptr
      select type(mesh2D_ptr)
      type is (MeshRectDom2D)
        call mesh2D_ptr%Final()
      end select
    end do
    deallocate( this%p_mesh_list )

    ! h-mesh hierarchy finalization
    do h_lev=MESH_HIERARCHY_hMG_FINEST_LEVEL, this%NUM_hMG_LEVEL
      call this%h_level(h_lev)%Final()
    end do
    deallocate( this%h_level )

    do h_lev=MESH_HIERARCHY_hMG_FINEST_LEVEL + 1, this%NUM_hMG_LEVEL
      mesh2D_ptr => this%h_mesh_list(h_lev)%ptr
      select type(mesh2D_ptr)
      type is (MeshRectDom2D)
        call mesh2D_ptr%Final()
      end select
    end do
    deallocate( this%h_mesh_list )
    
    !-
    call MeshHierarchyBase_Final( this )
    return
  end subroutine MeshHierarchy2D_Final

!-- private --------------------------------------------------------------

!OCL SERIAL
  subroutine MeshHierarchyLevel2D_Init( this, &
    level_id, mesh2D, mesh_list, LEVEL_NUM,   &
    hierarchy_type )
    use scale_mesh_hierarchy_base, only: &
      MeshHierarchy_construct_pMG_mat1D
    implicit none
    class(MeshHierarchyLevel2D), intent(inout) :: this
    integer, intent(in) :: level_id
    class(MeshBase2D), intent(in) :: mesh2D
    integer, intent(in) :: LEVEL_NUM
    class(MeshPtr2D), intent(in), target :: mesh_list(LEVEL_NUM)
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
      
      allocate( this%pMat1D_c2f(mesh2D%refElem2D%Nfp,this%coarse_mesh%ptr%refElem2D%Nfp) )
      call MeshHierarchy_construct_pMG_mat1D( this%pMat1D_c2f,   &
        this%coarse_mesh%ptr%refElem2D%Nfp, mesh2D%refElem2D%Nfp )

      allocate( this%pMat1D_f2c(this%coarse_mesh%ptr%refElem2D%Nfp,mesh2D%refElem2D%Nfp) )
      call MeshHierarchy_construct_pMG_mat1D( this%pMat1D_f2c,   &
        mesh2D%refElem2D%Nfp, this%coarse_mesh%ptr%refElem2D%Nfp )
    end if

    !- h-hierarchy
    if ( hierarchy_type == MESH_HIERARCHY_TYPE_hMG &
         .and. level_id < LEVEL_NUM                ) then

      allocate( this%mg_local( mesh2D%LOCAL_MESH_NUM ) )

      do ldom_id=1, mesh2D%LOCAL_MESH_NUM
        call this%mg_local(ldom_id)%Init( mesh_list(level_id)%ptr%lcmesh_list(ldom_id), &
          this%coarse_mesh%ptr%lcmesh_list, mesh_list(level_id)%ptr%refElem2D           )
      end do
    end if

    return
  end subroutine MeshHierarchyLevel2D_Init

!OCL SERIAL
  subroutine MeshHierarchyLevel2D_Final( this )
    implicit none
    class(MeshHierarchyLevel2D), intent(inout) :: this

    integer :: ldom_id
    !-------------------------------------------------------------

    !-
    if ( allocated(this%pMat1D_c2f) ) deallocate(this%pMat1D_c2f)
    if ( allocated(this%pMat1D_f2c) ) deallocate(this%pMat1D_f2c)

    !-
    if ( allocated(this%mg_local) ) then
      do ldom_id=1, size(this%mg_local)
        call this%mg_local(ldom_id)%Final()
      end do
      deallocate( this%mg_local )
    end if

    return
  end subroutine MeshHierarchyLevel2D_Final

!OCL SERIAL
  subroutine MeshHierarchyLocalMGData2D_Init( this, &
    lcmesh2D, coarse_lcmesh_list, elem2D )
    use scale_mesh_hierarchy_base, only: MeshHierarchyLocalMGDataBase_Init
    implicit none
    class(MeshHierarchyLocalMGData2D), intent(inout) :: this
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(LocalMesh2D), intent(in), target :: coarse_lcmesh_list(:)
    class(ElementBase2D), intent(in) :: elem2D

    integer :: i, j
    integer :: ke

    integer :: i_c, j_c
    integer :: ke_c

    integer :: ke2i(lcmesh2D%Ne)
    integer :: ke2j(lcmesh2D%Ne)

    integer :: lcdomID_c
    integer :: lcTileID_c
    class(LocalMesh2D), pointer :: lcmesh_c

    real(RP) :: vx_c(elem2D%Nv)
    real(RP) :: vy_c(elem2D%Nv)
    integer :: i_EtoV(elem2D%Nv)

    integer :: p
    real(RP) :: r_c, s_c
    !-------------------------------------------------------------  

    call MeshHierarchyLocalMGDataBase_Init( this, &
      lcmesh2D )

    
    ! Current implementation assumes identical tileID & lcdomID in both coarsened meshes          
    this%CoarseLocalMesh_tileIDlist(1) = lcmesh2D%tileID
    this%CoarseLocalMesh_lcdomIDlist(1) = 1

    ! Preparation
    do j=1, lcmesh2D%NeY
    do i=1, lcmesh2D%NeX
      ke = i + (j-1)*lcmesh2D%NeX
      ke2i(ke) = i; ke2j(ke) = j
    end do
    end do

    ! Set the relation between coarse and fine grid,
    ! and construct interpolation operator

!    write(*,*) "  Constructing local MG data for h-mesh..."

    allocate( this%If2c_emap(4,lcmesh2D%Ne/4) )

    do ke=lcmesh2D%NeS, lcmesh2D%NeE
      lcTileID_c = this%CoarseLocalMesh_tileIDlist(1)
      lcdomID_c = this%CoarseLocalMesh_lcdomIDlist(1)

      lcmesh_c => coarse_lcmesh_list(lcdomID_c)

      i_c = (ke2i(ke)+1)/2; j_c = (ke2j(ke)+1)/2
      ke_c = i_c + (j_c-1)*lcmesh2D%NeX/2
      this%Ic2f_emap(ke) = ke_c

      i = ke2i(ke) - 2*(i_c-1)
      j = ke2j(ke) - 2*(j_c-1)
      this%If2c_emap(i+(j-1)*2,ke_c) = ke

      i_EtoV(:) = lcmesh_c%EToV(ke_c,:)
      vx_c(:) = lcmesh_c%pos_ev(i_EtoV(:),1)
      vy_c(:) = lcmesh_c%pos_ev(i_EtoV(:),2)

      do p=1, elem2D%Np
        r_c = -1.0_RP + 2.0_RP * ( lcmesh2D%pos_en(p,ke,1) - vx_c(1) ) / ( vx_c(2) - vx_c(1) )
        s_c = -1.0_RP + 2.0_RP * ( lcmesh2D%pos_en(p,ke,2) - vy_c(1) ) / ( vy_c(3) - vy_c(1) )

        this%Ic2f(p,1:4,ke) = 0.25_RP * &
          (/ ( 1.0_RP - r_c ) * ( 1.0_RP - s_c ), ( 1.0_RP + r_c ) * ( 1.0_RP - s_c ),  &
             ( 1.0_RP - r_c ) * ( 1.0_RP + s_c ), ( 1.0_RP + r_c ) * ( 1.0_RP + s_c )  /)
      end do
    end do  
    return
  end subroutine MeshHierarchyLocalMGData2D_Init

!OCL SERIAL
  subroutine MeshHierarchyLocalMGData2D_Final( this )
    use scale_mesh_hierarchy_base, only: MeshHierarchyLocalMGDataBase_Final
    implicit none
    class(MeshHierarchyLocalMGData2D), intent(inout) :: this
    !-------------------------------------------------------------  
    call MeshHierarchyLocalMGDataBase_Final( this )
    return
  end subroutine MeshHierarchyLocalMGData2D_Final

!-
!OCL SERIAL
  subroutine construct_rectdom2D_mesh( this, child_mesh_base_ptr, &
    NeGX, NeGY,         &
    parent_mesh, elem2D )
    implicit none
    class(MeshHierarchy2D), intent(inout) :: this
    class(MeshBase2D), intent(out), pointer :: child_mesh_base_ptr
    integer, intent(in) :: NeGX
    integer, intent(in) :: NeGY
    type(MeshRectDom2D), intent(in) :: parent_mesh
    type(QuadrilateralElement), intent(in) :: elem2D

    type(MeshRectDom2D), pointer :: child_mesh_ptr
    !-------------------------------------------------------------

    allocate( child_mesh_ptr )
    call child_mesh_ptr%Init( NeGX, NeGY,                                                  &
      parent_mesh%xmin_gl, parent_mesh%xmax_gl, parent_mesh%ymin_gl, parent_mesh%ymax_gl,  &
      parent_mesh%isPeriodicX, parent_mesh%isPeriodicY,                                    &
      elem2D, 1, parent_mesh%NprcX, parent_mesh%NprcY )

    call child_mesh_ptr%Generate()

    child_mesh_base_ptr => child_mesh_ptr
    return
  end subroutine construct_rectdom2D_mesh
end module scale_mesh_hierarchy_2d
