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
    
  use scale_meshfield_base, only: MeshField3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D

  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D

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
    procedure(exp_geostrophic_balance_correction_lc), deferred :: geostrophic_balance_correction_lc
  end type

  abstract interface
    subroutine exp_SetInitCond_lc( &
      this, DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,  &
      x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, &
      dom_zmax, lcmesh, elem )

      import experiment
      import LocalMesh3D 
      import ElementBase3D
      import RP

      class(experiment), intent(inout) :: this
      type(LocalMesh3D), intent(in) :: lcmesh
      class(ElementBase3D), intent(in) :: elem
      real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: MOMY(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
      real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
      real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
      real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
      real(RP), intent(in) :: dom_xmin, dom_xmax
      real(RP), intent(in) :: dom_ymin, dom_ymax      
      real(RP), intent(in) :: dom_zmin, dom_zmax
    end subroutine exp_SetInitCond_lc

    subroutine exp_geostrophic_balance_correction_lc( this,                              &
      DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
      lcmesh, elem )
    
      import experiment
      import LocalMesh3D 
      import ElementBase3D
      import RP

      class(experiment), intent(inout) :: this
      type(LocalMesh3D), intent(in) :: lcmesh
      class(ElementBase3D), intent(in) :: elem
      real(RP), intent(inout) :: DENS_hyd(elem%Np,lcmesh%NeA)
      real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
      real(RP), intent(inout) :: DDENS(elem%Np,lcmesh%NeA)
      real(RP), intent(inout) :: MOMX(elem%Np,lcmesh%NeA)
      real(RP), intent(inout) :: MOMY(elem%Np,lcmesh%NeA)    
      real(RP), intent(inout) :: MOMZ(elem%Np,lcmesh%NeA)
      real(RP), intent(inout) :: DRHOT(elem%Np,lcmesh%NeA)
    end subroutine exp_geostrophic_balance_correction_lc

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
    model_mesh, atm_prgvars_manager, atm_auxvars_manager )
    
    use scale_meshfield_base, only: MeshFieldBase
    use scale_model_var_manager, only: ModelVarManager
    use scale_meshfieldcomm_base, only: MeshFieldContainer  
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshFields, &
      ATMOS_AUXVARS_PRESHYDRO_ID, &
      ATMOS_AUXVARS_DENSHYDRO_ID
    use mod_atmos_mesh, only: AtmosMesh
    
    implicit none

    class(experiment), intent(inout) :: this
    class(AtmosMesh), target, intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: atm_prgvars_manager
    class(ModelVarManager), intent(inout) :: atm_auxvars_manager

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: &
      GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW, GxPT, GyPT, GzPT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd

    integer :: n
    class(LocalMesh3D), pointer :: lcmesh3D
    class(MeshCubeDom3D), pointer :: mesh

    type(MeshFieldCommCubeDom3D) :: hydvars_comm
    type(MeshFieldContainer) :: hydvars_comm_list(1)
    class(MeshFieldBase), pointer :: field_ptr
    !----------------------------------------------------------------------
    
    mesh => model_mesh%mesh
    
    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshFields( n, mesh, atm_prgvars_manager, atm_auxvars_manager, &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                                                     &
        GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW, GxPT, GyPT, GzPT,                      &
        DENS_hyd, PRES_hyd, lcmesh3D                                                        )

      call this%setInitCond_lc( &
        DENS_hyd%val, PRES_hyd%val,                                                         & ! (out)
        DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                                 & ! (out)
        lcmesh3D%pos_en(:,:,1), lcmesh3D%pos_en(:,:,2), lcmesh3D%pos_en(:,:,3),             & ! (in)
        mesh%xmin_gl, mesh%xmax_gl, mesh%ymin_gl, mesh%ymax_gl, mesh%zmin_gl, mesh%zmax_gl, & ! (in)
        lcmesh3D, lcmesh3D%refElem3D )                                                        ! (in) 
    end do

    !------------------------------------------------------
    call hydvars_comm%Init(1, 0, mesh)

    call atm_auxvars_manager%Get(ATMOS_AUXVARS_PRESHYDRO_ID, field_ptr)
    select type(field_ptr) 
    type is (MeshField3D)
      hydvars_comm_list(1)%field3d => field_ptr
    end select
    call hydvars_comm%Put(hydvars_comm_list, 1)
    call hydvars_comm%Exchange()
    call hydvars_comm%Get(hydvars_comm_list, 1)

    !--------------------------------
    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshFields( n, mesh, atm_prgvars_manager, atm_auxvars_manager, &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                                                     &
        GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW, GxPT, GyPT, GzPT,                      &
        DENS_hyd, PRES_hyd, lcmesh3D                                                        )

      call this%geostrophic_balance_correction_lc( &
        DENS_hyd%val, PRES_hyd%val,                                                         & ! (out)
        DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                                 & ! (out)
        lcmesh3D, lcmesh3D%refElem3D )                                                        ! (in) 
    end do

    call atm_auxvars_manager%Get(ATMOS_AUXVARS_DENSHYDRO_ID, field_ptr)
    select type(field_ptr) 
    type is (MeshField3D)
      hydvars_comm_list(1)%field3d => field_ptr
    end select
    call hydvars_comm%Put(hydvars_comm_list, 1)
    call hydvars_comm%Exchange()
    call hydvars_comm%Get(hydvars_comm_list, 1)
    call hydvars_comm%Final()

    return
  end subroutine exp_SetInitCond
  
end module mod_exp