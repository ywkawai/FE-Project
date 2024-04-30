!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_experiment
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

  use scale_mesh_base3d, only: MeshBase3D    
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_meshfieldcomm_base, only: MeshFieldCommBase
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D  
  use scale_meshfieldcomm_cubedspheredom3d, only: MeshFieldCommCubedSphereDom3D
  use scale_meshfield_base, only: MeshField3D
  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  type, public :: Experiment
    character(len=H_SHORT) :: label
    procedure(exp_SetInitCond_lc), pointer :: setInitCond_lc => null()
    procedure(exp_geostrophic_balance_correction_lc), pointer :: geostrophic_balance_correction_lc => null()
  contains
    procedure, public :: Init_Base => experiment_Init
    generic :: Init => Init_Base
    procedure, public :: Final_Base => experiment_Final
    generic :: Final => Final_Base
    procedure, public :: SetInitCond => experiment_SetInitCond
    procedure, public :: Regist_SetInitCond => experiment_regist_set_initcond
    procedure, public :: Regist_geostrophic_balance_correction => experiment_regist_geostrophic_balance_correction
  end type Experiment

  type, public :: TracerLocalMeshField_ptr
    class(LocalMeshFieldBase), pointer :: ptr
  end type

  interface
    subroutine exp_SetInitCond_lc( &
      this, DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,  &
      tracer_field_list,                                         &
      x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, &
      dom_zmax, lcmesh, elem )

      import Experiment
      import LocalMesh3D 
      import ElementBase3D
      import TracerLocalMeshField_ptr
      import RP

      class(Experiment), intent(inout) :: this
      type(LocalMesh3D), intent(in) :: lcmesh
      class(ElementBase3D), intent(in) :: elem
      real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: MOMY(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
      type(TracerLocalMeshField_ptr), intent(inout) :: tracer_field_list(:)
      real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
      real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
      real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
      real(RP), intent(in) :: dom_xmin, dom_xmax
      real(RP), intent(in) :: dom_ymin, dom_ymax      
      real(RP), intent(in) :: dom_zmin, dom_zmax
    end subroutine exp_SetInitCond_lc

    subroutine exp_geostrophic_balance_correction_lc( this,                &
      DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
      lcmesh, elem )
    
      import Experiment
      import LocalMesh3D 
      import ElementBase3D
      import RP

      class(Experiment), intent(inout) :: this
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
  subroutine experiment_Init( this, exp_name )
    implicit none
    class(Experiment), intent(inout) :: this
    character(len=*), intent(in) :: exp_name
    !----------------------------------------------------------------------

    this%label = exp_name

    this%setInitCond_lc => experiment_SetInitCond_lc_dummy
    this%geostrophic_balance_correction_lc => experiment_geostrophic_balance_correction_lc_dummy

    return
  end subroutine experiment_Init

  subroutine experiment_Final( this )
    implicit none
    class(Experiment), intent(inout) :: this
    !----------------------------------------------------------------------

    return
  end subroutine experiment_Final

  subroutine experiment_regist_set_initcond( this, exp_SetInitCond_lc )
    implicit none
    class(Experiment), intent(inout) :: this
    interface
      subroutine exp_SetInitCond_lc( &
        this, &
        DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,  &
        tracer_field_list,                                         &
        x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, &
        dom_zmax, lcmesh, elem )

        import Experiment
        import LocalMesh3D 
        import ElementBase3D
        import TracerLocalMeshField_ptr
        import RP

        class(Experiment), intent(inout) :: this
        type(LocalMesh3D), intent(in) :: lcmesh
        class(ElementBase3D), intent(in) :: elem
        real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
        real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
        real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
        real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
        real(RP), intent(out) :: MOMY(elem%Np,lcmesh%NeA)
        real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
        real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
        type(TracerLocalMeshField_ptr), intent(inout) :: tracer_field_list(:)
        real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
        real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
        real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
        real(RP), intent(in) :: dom_xmin, dom_xmax
        real(RP), intent(in) :: dom_ymin, dom_ymax      
        real(RP), intent(in) :: dom_zmin, dom_zmax
      end subroutine exp_SetInitCond_lc
    end interface    
    !----------------------------------------------------------------------

    this%setInitCond_lc => exp_SetInitCond_lc
    return
  end subroutine experiment_regist_set_initcond

  subroutine experiment_regist_geostrophic_balance_correction( this, exp_geostrophic_balance_correction_lc )
    implicit none
    class(Experiment), intent(inout) :: this
    interface
      subroutine exp_geostrophic_balance_correction_lc( this,                &
        DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
        lcmesh, elem )
      
        import Experiment
        import LocalMesh3D 
        import ElementBase3D
        import RP

        class(Experiment), intent(inout) :: this
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
    !----------------------------------------------------------------------

    this%geostrophic_balance_correction_lc => exp_geostrophic_balance_correction_lc
    return
  end subroutine experiment_regist_geostrophic_balance_correction

  subroutine experiment_SetInitCond( this, &
    model_mesh, atm_prgvars_manager, atm_auxvars_manager, atm_trcvars_manager )
    
    use scale_tracer, only: QA

    use scale_meshfield_base, only: MeshFieldBase
    use scale_model_var_manager, only: ModelVarManager
    use scale_meshfieldcomm_base, only: MeshFieldContainer 

    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      AUXVAR_PRESHYDRO_ID, &
      AUXVAR_DENSHYDRO_ID

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars, &
      AtmosVars_GetLocalMeshQTRCVar
    use mod_atmos_mesh, only: AtmosMesh
    
    implicit none

    class(Experiment), intent(inout) :: this
    class(AtmosMesh), target, intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: atm_prgvars_manager
    class(ModelVarManager), intent(inout) :: atm_auxvars_manager
    class(ModelVarManager), intent(inout) :: atm_trcvars_manager

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot

    integer :: n
    class(LocalMesh3D), pointer :: lcmesh3D
    class(MeshBase3D), pointer :: mesh

    type(MeshFieldCommCubeDom3D), target :: hydvars_comm_rm
    type(MeshFieldCommCubedSphereDom3D), target :: hydvars_comm_gm
    class(MeshFieldCommBase), pointer :: hydvars_comm
    type(MeshFieldContainer) :: hydvars_comm_list(1)
    class(MeshFieldBase), pointer :: field_ptr

    type(TracerLocalMeshField_ptr) :: tracer_field_list(max(1,QA))
    integer :: iq
    !----------------------------------------------------------------------
    
    mesh => model_mesh%ptr_mesh
    
    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, &
        mesh, atm_prgvars_manager, atm_auxvars_manager, &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                 &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot,         &
        lcmesh3D                                        )
      
      do iq=1, QA
        call AtmosVars_GetLocalMeshQTRCVar( n, mesh, atm_trcvars_manager, &
          iq, tracer_field_list(iq)%ptr )
      end do

      select type (mesh)
      type is (MeshCubeDom3D)
        call this%setInitCond_lc( &
          DENS_hyd%val, PRES_hyd%val,                                                         & ! (out)
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                                 & ! (out)
          tracer_field_list,                                                                  & ! (inout)
          lcmesh3D%pos_en(:,:,1), lcmesh3D%pos_en(:,:,2), lcmesh3D%pos_en(:,:,3),             & ! (in)
          mesh%xmin_gl, mesh%xmax_gl, mesh%ymin_gl, mesh%ymax_gl, mesh%zmin_gl, mesh%zmax_gl, & ! (in)
          lcmesh3D, lcmesh3D%refElem3D )                                                        ! (in) 
      type is (MeshCubedSphereDom3D)
        call this%setInitCond_lc( &
          DENS_hyd%val, PRES_hyd%val,                                                         & ! (out)
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                                 & ! (out)
          tracer_field_list,                                                                  & ! (inout)          
          lcmesh3D%pos_en(:,:,1), lcmesh3D%pos_en(:,:,2), lcmesh3D%pos_en(:,:,3),             & ! (in)
          mesh%xmin_gl, mesh%xmax_gl, mesh%ymin_gl, mesh%ymax_gl, mesh%zmin_gl, mesh%zmax_gl, & ! (in)
          lcmesh3D, lcmesh3D%refElem3D )                                                        ! (in)   
      end select
    end do

    !------------------------------------------------------

    call atm_auxvars_manager%Get(AUXVAR_PRESHYDRO_ID, field_ptr)
    select type(field_ptr) 
    type is (MeshField3D)
      hydvars_comm_list(1)%field3d => field_ptr
    end select
    select type (mesh)
    type is (MeshCubeDom3D)
      call hydvars_comm_rm%Init(1, 0, 0, mesh)
      hydvars_comm => hydvars_comm_rm
    type is (MeshCubedSphereDom3D)
      call hydvars_comm_gm%Init(1, 0, 0, mesh)
      hydvars_comm => hydvars_comm_gm
    end select

    call hydvars_comm%Put(hydvars_comm_list, 1)
    call hydvars_comm%Exchange()
    call hydvars_comm%Get(hydvars_comm_list, 1)

    !--------------------------------
    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, &
        mesh, atm_prgvars_manager, atm_auxvars_manager, &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                 &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot,         &
        lcmesh3D                                        )

      call this%geostrophic_balance_correction_lc( &
        DENS_hyd%val, PRES_hyd%val,                          & ! (out)
        DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,  & ! (out)
        lcmesh3D, lcmesh3D%refElem3D )                         ! (in) 
    end do

    call atm_auxvars_manager%Get(AUXVAR_DENSHYDRO_ID, field_ptr)
    select type(field_ptr) 
    type is (MeshField3D)
      hydvars_comm_list(1)%field3d => field_ptr
    end select
    call hydvars_comm%Put(hydvars_comm_list, 1)
    call hydvars_comm%Exchange()
    call hydvars_comm%Get(hydvars_comm_list, 1)

    select type (mesh)
    type is (MeshCubeDom3D)
      call hydvars_comm_rm%Final()
    type is (MeshCubedSphereDom3D)
      call hydvars_comm_gm%Final()
    end select    

    return
  end subroutine experiment_SetInitCond
  
  !------
!OCL SERIAL  
  subroutine experiment_SetInitCond_lc_dummy( this,                   &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    implicit none

    class(Experiment), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
    type(TracerLocalMeshField_ptr), intent(inout) :: tracer_field_list(:)    
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_ymin, dom_ymax
    real(RP), intent(in) :: dom_zmin, dom_zmax
    !---------------------------------------------------

    return
  end subroutine experiment_SetInitCond_lc_dummy

  subroutine experiment_geostrophic_balance_correction_lc_dummy( this,    &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,   &
    lcmesh, elem )
    
    implicit none
    class(Experiment), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(inout) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DRHOT(elem%Np,lcmesh%NeA)
    !---------------------------------------------------
    return
  end subroutine experiment_geostrophic_balance_correction_lc_dummy

end module mod_experiment