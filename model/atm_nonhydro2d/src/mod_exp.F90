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
    GRAV => CONST_GRAV,   &
    Rdry => CONST_Rdry,   &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00
    
  use scale_meshfield_base, only: MeshField2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D

  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

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
    procedure(exp_SetInitCond_lc), private, deferred :: setInitCond_lc
  end type

  abstract interface
    subroutine exp_SetInitCond_lc( &
      this, DENS_hyd, PRES_hyd, DDENS, MOMX, MOMZ, DRHOT, &
      x, z, dom_xmin, dom_xmax, dom_zmin, dom_zmax, lcmesh, elem )

      import experiment
      import LocalMesh2D 
      import QuadrilateralElement
      import RP

      class(experiment), intent(inout) :: this
      type(LocalMesh2D), intent(in) :: lcmesh
      class(QuadrilateralElement), intent(in) :: elem
      real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
      real(RP), intent(in) :: x(elem%Np,lcmesh%NeA)
      real(RP), intent(in) :: z(elem%Np,lcmesh%NeA)
      real(RP), intent(in) :: dom_xmin, dom_xmax
      real(RP), intent(in) :: dom_zmin, dom_zmax
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
  
  subroutine exp_SetInitCond( this )
    use mod_atmos_vars, only: &
      DENS_hydro, PRES_hydro,   &
      DDENS, MOMX, MOMZ, DRHOT
    use mod_atmos_mesh, only: &
      mesh, refElem
    implicit none

    class(experiment), intent(inout) :: this

    integer :: n
    type(LocalMesh2D), pointer :: lcmesh
    !----------------------------------------------------------------------
    
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call this%setInitCond_lc( &
        DENS_hydro%local(n)%val, PRES_hydro%local(n)%val,                                   & ! (out)
        DDENS%local(n)%val, MOMX%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,       & ! (out)
        lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),                                         & ! (in)
        mesh%xmin_gl, mesh%xmax_gl, mesh%ymin_gl, mesh%ymax_gl, lcmesh, refElem )             ! (in) 
    end do

    return
  end subroutine exp_SetInitCond
  
end module mod_exp