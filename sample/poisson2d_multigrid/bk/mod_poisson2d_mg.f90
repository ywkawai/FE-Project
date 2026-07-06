#include "scaleFElib.h"
module mod_poisson2d_mg
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale_prc
  use scale_io

  use scale_sparsemat, only: SparseMat
  
  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D
  
  use scale_mesh_hierarchy_base, only: &
    pMG_FINEST_LEVEL => MESH_HIERARCHY_pMG_FINEST_LEVEL, &
    hMG_FINEST_LEVEL => MESH_HIERARCHY_hMG_FINEST_LEVEL

  use scale_mesh_hierarchy_2d, only: MeshHierarchy2D
  use scale_multigrid_solver_2d, only: MultiGridSolver2D

  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  use mod_poisson2d_smoother, only: MGSmoother_Poisson2D
  
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  public :: Poisson2d_mg_Init
  public :: Poisson2d_mg_Final
  public :: Poisson2d_MG_solve
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure & variables
  !
  type(MeshHierarchy2D) :: mesh_hierarchy
  type(MGSmoother_Poisson2D) :: smoother
  type(MultiGridSolver2D) :: mg_solver

contains
  !> Initialization
  subroutine Poisson2d_mg_Init( mesh )
    implicit none
    class(MeshBase2D), intent(in) :: mesh
    
    integer :: lev_h
    integer :: lev_p

    integer, parameter :: pMG_LV_LIST_MAX = 16
    integer :: P_LEVEL_NUM
    integer :: P_LEVEL_LIST(pMG_LV_LIST_MAX)
    integer :: P_ITR_NUM_LIST(pMG_LV_LIST_MAX)


    integer, parameter :: hMG_LV_LIST_MAX = 16
    integer :: H_LEVEL_NUM
    integer :: NeGX_list(hMG_LV_LIST_MAX)
    integer :: NeGY_list(hMG_LV_LIST_MAX)
    integer :: H_ITR_NUM_LIST(hMG_LV_LIST_MAX)

    namelist /PARAM_Poisson2D_MG/ &
      H_LEVEL_NUM, NeGX_list, NeGY_list, H_ITR_NUM_LIST, &
      P_LEVEL_NUM, P_LEVEL_LIST, P_ITR_NUM_LIST

    integer :: ierr
    !---------------------------------------------------------------------------


    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_Poisson2D_MG,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("Poisson2D_MG_Init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("Poisson2D_MG_Init",*) 'Not appropriate names in namelist PARAM_Poisson2D_MG. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_Poisson2D_MG)

    !- Setup an object to manage the mesh hierarchy

    LOG_INFO("Poisson2D_MG_Init",*) "Setup mesh hierarchy and MG solver"
    call mesh_hierarchy%Init( mesh, &
      P_LEVEL_LIST(1:P_LEVEL_NUM), P_LEVEL_NUM,                       &
      NeGX_list(1:H_LEVEL_NUM), NeGY_list(1:H_LEVEL_NUM), H_LEVEL_NUM )

    !- Setup an object to provide a smoother
    LOG_INFO("Poisson2D_MG_Init",*) "Setup smoother"
    call smoother%Init( mesh )

    !- Setup an object to provide a MG solver   
    LOG_INFO("Poisson2D_MG_Init",*) "Setup MG solver"    
    call mg_solver%Init( mesh_hierarchy, smoother, 0, 2, &
      P_ITR_NUM_LIST(1:P_LEVEL_NUM), H_ITR_NUM_LIST(1:H_LEVEL_NUM) )
    return
  end subroutine Poisson2d_mg_Init

  !> Finalization
  subroutine Poisson2d_mg_Final()
    use mod_poisson2d_smoother, only: Poisson2d_smoother_Final
    implicit none
    !---------------------------------------------------------------------------
    call smoother%Final()
    call mg_solver%Final()
    call mesh_hierarchy%Final()
    return
  end subroutine Poisson2d_mg_Final

!OCL SERIAL
  subroutine Poisson2d_MG_solve( q, &
    f )
    implicit none
    class(MeshField2D), intent(inout), target :: q
    class(MeshField2D), intent(inout) :: f
    !-------------------------------------
    call mg_solver%Solve( q, &
      f )
    return
  end subroutine Poisson2d_MG_solve
end module mod_poisson2d_mg
