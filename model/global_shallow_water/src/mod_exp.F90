!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_exp
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_const, only: &
    PI => CONST_PI,       &
    GRAV => CONST_GRAV
    
  use scale_meshfield_base, only: MeshField2D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D

  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_meshfieldcomm_cubedspheredom2d, only: MeshFieldCommCubedSphereDom2D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, abstract, public :: experiment
    character(len=H_SHORT) :: label
  contains
    procedure, public :: Init => exp_Init
    procedure, public :: Final => exp_Final
    procedure, public :: SetInitCond => exp_SetInitCond
    procedure(exp_SetInitCond_lc), deferred :: setInitCond_lc
  end type

  abstract interface
    subroutine exp_SetInitCond_lc( &
      this, h, U, V, hs, u1, u2,   &
      x, y, lcmesh, elem )

      import experiment
      import LocalMesh2D 
      import ElementBase2D
      import RP

      class(experiment), intent(inout) :: this
      type(LocalMesh2D), intent(in) :: lcmesh
      class(ElementBase2D), intent(in) :: elem
      real(RP), intent(out) :: h(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: U(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: V(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: hs(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: u1(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: u2(elem%Np,lcmesh%NeA)
      real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
      real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
    end subroutine exp_SetInitCond_lc
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
  subroutine exp_Init( this, exp_name )
    implicit none
    class(experiment), intent(inout) :: this
    character(len=*), intent(in) :: exp_name
    !----------------------------------------------------------------------

    this%label = exp_name
    return
  end subroutine exp_Init

  subroutine exp_Final( this )
    implicit none
    class(experiment), intent(inout) :: this
    !----------------------------------------------------------------------

    return
  end subroutine exp_Final
  
  subroutine exp_SetInitCond( this, &
    model_mesh, sw_prgvars_manager, sw_auxvars_manager )
    
    use scale_meshfield_base, only: MeshFieldBase
    use scale_model_var_manager, only: ModelVarManager
    use scale_meshfieldcomm_base, only: MeshFieldContainer  
    use mod_sw_vars, only: &
      SWVars_GetLocalMeshPrgVars
    use mod_sw_mesh, only: SWMesh

    implicit none

    class(experiment), intent(inout) :: this
    class(SWMesh), target, intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: sw_prgvars_manager
    class(ModelVarManager), intent(inout) :: sw_auxvars_manager

    class(LocalMeshFieldBase), pointer :: h, U, V
    class(LocalMeshFieldBase), pointer :: hs, u1, u2

    integer :: n
    class(LocalMesh2D), pointer :: lcmesh
    class(MeshCubedSphereDom2D), pointer :: mesh
    !----------------------------------------------------------------------
    
    mesh => model_mesh%mesh

    do n=1, mesh%LOCAL_MESH_NUM
      call SWVars_GetLocalMeshPrgVars( &
        n, mesh, sw_prgvars_manager, sw_auxvars_manager, &
        h, U, V, hs, u1, u2, lcmesh                      )

      call this%setInitCond_lc( &
        h%val, U%val, V%val, hs%val, u1%val, u2%val,   & ! (out)
        lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),    & ! (in)
        lcmesh, lcmesh%refElem2D                       ) ! (in) 
    end do

    return
  end subroutine exp_SetInitCond
  
end module mod_exp