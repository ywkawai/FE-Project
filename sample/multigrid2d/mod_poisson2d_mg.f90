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
  use scale_mesh_hierarchy_2d, only: MeshHierarchy2D

  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  public :: Poisson2d_mg_Init
  public :: Poisson2d_mg_Final
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !

  type, public :: FieldSet
    type(MeshField2D) :: f
    type(MeshField2D) :: res
    type(MeshField2D) :: qx
    type(MeshField2D) :: qy
    integer :: level_id

    type(MeshFieldCommRectDom2D) :: var_comm
    type(MeshFieldCommRectDom2D) :: aux_comm
  contains
    procedure :: Init => FieldSet_Init
    procedure :: Final => FieldSet_Final
  end type FieldSet

  type(FieldSet), allocatable :: fields_h(:)
  type(FieldSet), allocatable :: fields_p(:)

contains
  !> Initialization
  subroutine Poisson2d_mg_Init( mesh_hierarchy )
    implicit none
    class(MeshHierarchy2D), intent(inout) :: mesh_hierarchy

    integer :: lev_h
    integer :: lev_p
    !---------------------------------------------------------------------------

    allocate( fields_h(mesh_hierarchy%NUM_hMG_LEVEL) )
    do lev_h=1, mesh_hierarchy%NUM_hMG_LEVEL
      call fields_h(lev_h)%Init( mesh_hierarchy%h_mesh_list(lev_h)%ptr, lev_h )
    end do

    allocate( fields_p(mesh_hierarchy%NUM_pMG_LEVEL) )
    do lev_p=1, mesh_hierarchy%NUM_pMG_LEVEL
      call fields_p(lev_p)%Init( mesh_hierarchy%p_mesh_list(lev_p)%ptr, lev_p )
    end do

    return
  end subroutine Poisson2d_mg_Init

  !> Finalization
  subroutine Poisson2d_mg_Final()
    implicit none

    integer :: lev_h
    integer :: lev_p
    !---------------------------------------------------------------------------

    do lev_h=1, size(fields_h)
      call fields_h(lev_h)%Final()
    end do

    do lev_p=1, size(fields_p)
      call fields_p(lev_p)%Final()
    end do    
    return
  end subroutine Poisson2d_mg_Final


!OCL SERIAL
  recursive subroutine Poisson2d_mg_Vcycle( mg_level, f_in, &
    mesh_hierarchy )
    implicit none
    integer, intent(in) :: mg_level
    type(MeshField2D), intent(inout) :: f_in
    class(MeshHierarchy2D), intent(in) :: mesh_hierarchy

    real(RP) :: itr_res_eps
    integer :: m

    logical :: invoke_hMG
    !----------------------------------------------

    invoke_hMG = ( mg_level+1 > mesh_hierarchy%NUM_pMG_LEVEL .and. mesh_hierarchy%NUM_hMG_LEVEL > 0 )

    !- Pre-relaxation
    !do m=1, ITR_NUM
      ! call 
    !end do

    !-
    if ( invoke_hMG ) then
    else
      call mesh_hierarchy%Operate_pMG_restriction( fields_p(mg_level+1)%f, &
        fields_p(mg_level)%res, mg_level )

      call Poisson2d_mg_Vcycle( mg_level+1, fields_p(mg_level+1)%f, mesh_hierarchy)
    end if

    !- Post-relaxation
    !do m=1, ITR_NUM
      ! call 
    !end do

    return
  end subroutine Poisson2d_mg_Vcycle

!OCL SERIAL
  subroutine Poisson2d_mg_h_correction( lcmesh, elem2D )
    implicit none
    class(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem2D
    !-------------------------------------------
    return
  end subroutine Poisson2d_mg_h_correction

!OCL SERIAL
  subroutine Poisson2d_mg_h_restriction( lcmesh, elem2D )
    implicit none
    class(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem2D
    !-------------------------------------------
    return
  end subroutine Poisson2d_mg_h_restriction

!--------------

  subroutine FieldSet_Init(this, mesh2D, level_id)
    use scale_mesh_rectdom2d, only: MeshRectDom2D
    implicit none
    class(FieldSet), intent(inout) :: this
    class(MeshBase2D), intent(in) :: mesh2D
    integer, intent(in) :: level_id
    !---------------------------------

    this%level_id = level_id
    call this%f%Init( "f", "1", mesh2D )
    call this%res%Init( "res", "1", mesh2D )
    call this%qx%Init( "qx", "1", mesh2D )
    call this%qy%Init( "qy", "1", mesh2D )

    select type(mesh2D)
    type is (MeshRectDom2D)
      call this%var_comm%Init( 1, 0, 0, mesh2D )
      call this%aux_comm%Init( 0, 1, 0, mesh2D )
    end select
    return
  end subroutine FieldSet_Init

  subroutine FieldSet_Final(this)
    implicit none
    class(FieldSet), intent(inout) :: this
    !---------------------------------
    call this%var_comm%Final()
    call this%aux_comm%Final()

    call this%f%Final()    
    call this%res%Final()
    call this%qx%Final()
    call this%qy%Final()
    return
  end subroutine FieldSet_Final

end module mod_poisson2d_mg
