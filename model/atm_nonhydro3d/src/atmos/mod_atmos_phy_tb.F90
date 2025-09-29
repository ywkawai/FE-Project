!-------------------------------------------------------------------------------
!> module ATMOSPHERE physics / sub-grid scale turbulence process
!!
!! @par Description
!!          Module for sub-grid scale turbulence process
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_tb
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof
  use scale_const, only: &
    UNDEF8 => CONST_UNDEF8
  use scale_tracer, only: QA

  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
    
  use scale_meshfield_base, only: MeshFieldBase
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_component_proc, only:  ModelComponentProc

  use scale_atm_dyn_dgm_bnd, only: AtmDynBnd
  
  use scale_atm_phy_tb_dgm_driver, only: &
    AtmPhyTbDGMDriver
  use mod_atmos_phy_tb_vars, only: AtmosPhyTbVars

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage a component of sub-grid scale turbulent process
  !!
  type, extends(ModelComponentProc), public :: AtmosPhyTb
    integer :: TB_TYPEID                  !< Type ID of SGS turbulent scheme
    type(AtmPhyTbDGMDriver) :: tb_driver  !< Object to represent a driver for SGS turbulent schemes

    type(AtmosPhyTbVars) :: vars          !< Object to manage variables with SGS turbulent model
    type(AtmDynBnd), pointer :: dyn_bnd   !< Pointer to object for treating boundary conditions with atmospheric dynamics  
  contains
    procedure, public :: setup => AtmosPhyTb_setup 
    procedure, public :: calc_tendency => AtmosPhyTb_calc_tendency
    procedure, public :: update => AtmosPhyTb_update
    procedure, public :: finalize => AtmosPhyTb_finalize
    procedure, public :: SetDynBC => AtmosPhyTB_SetDynBC
  end type AtmosPhyTb

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
    
contains
!> Setup a component of SGS turbulence process
!!
!! @param model_mesh a object to manage computational mesh of atmospheric model 
!! @param tm_parent_comp a object to mange a temporal scheme in a parent component
!!
!OCL SERIAL
  subroutine AtmosPhyTb_setup( this, model_mesh, tm_parent_comp )
    use mod_atmos_mesh, only: AtmosMesh
    use scale_time_manager, only: TIME_manager_component

    implicit none
    class(AtmosPhyTb), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8  !< Timestep for sub-grid scale turbulent process
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'   !< Unit of timestep

    character(len=H_MID) :: TB_TYPE = 'SMAGORINSKY'           !< Type of sub-grid scale turbulent scheme
    namelist /PARAM_ATMOS_PHY_TB/ &
      TIME_DT,          & 
      TIME_DT_UNIT,     &
      TB_TYPE
    
    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh
    real(DP) :: dtsec

    integer :: ierr
    !--------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_TB_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_PHY_TB_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_TB. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_TB)

    !- get mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select

    !--- Register this component in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_TB', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                       ! (out)

    dtsec = tm_parent_comp%process_list(this%tm_process_id)%dtsec

    !- initialize the variables 
    call this%vars%Init( model_mesh )      

    !- Initialize a module for turbulence model
    call this%tb_driver%Init( TB_TYPE, dtsec, atm_mesh )

    !--
    this%dyn_bnd => null()

    return
  end subroutine AtmosPhyTb_setup

!> Calculate tendencies associated with a turbulent model
!!
!!
!! @param model_mesh a object to manage computational mesh of atmospheric model 
!! @param prgvars_list a object to mange prognostic variables with atmospheric dynamical core
!! @param trcvars_list a object to mange auxiliary variables 
!! @param forcing_list a object to mange forcing terms
!! @param is_update Flag to specify whether the tendencies are updated in this call
!!
!OCL SERIAL
  subroutine AtmosPhyTb_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )

    use scale_tracer, only: &
      TRACER_ADVC   
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRGVAR_DDENS_ID
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPrgVar,     &
      AtmosVars_GetLocalMeshQTRCVar,    &
      AtmosVars_GetLocalMeshPhyAuxVars, &
      AtmosVars_GetLocalMeshPhyTends
    use mod_atmos_phy_tb_vars, only:          &
      AtmosPhyTbVars_GetLocalMeshFields_tend
     
    implicit none
    
    class(AtmosPhyTb), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    class(MeshBase), pointer :: mesh
    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lcmesh
    integer :: n
    integer :: ke
    integer :: iq

    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_P
    type(LocalMeshFieldBaseList) :: RHOQ_tp(QA)
    class(LocalMeshFieldBase), pointer :: tb_MOMX_t, tb_MOMY_t, tb_MOMZ_t, tb_RHOT_t
    class(LocalMeshFieldBase), pointer :: tb_RHOQ_t
    type(LocalMeshFieldBaseList) :: tb_RHOQ_t_list(QA)
    !--------------------------------------------------

    if (.not. this%IsActivated()) return

    !LOG_INFO('AtmosDyn_tendency',*)

    LOG_PROGRESS(*) 'atmosphere / physics / turbulence'

    call model_mesh%GetModelMesh( mesh )
    select type(mesh)
    class is (MeshBase3D)
      mesh3D => mesh
    end select

    if (is_update) then
      call this%tb_driver%Tendency( this%vars%tends_manager, &
        prgvars_list, trcvars_list, auxvars_list,                                  &
        this%vars%auxvars_manager, this%vars%diagvars_manager,                     &
        this%dyn_bnd,                                                              &
        model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),    &
        model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3),    &
        model_mesh%LiftOptrMat, mesh3D                                             )
    end if
      
    call PROF_rapstart('ATM_PHY_TB_add_tend', 2)
    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPhyTends( n,        &
        mesh, forcing_list,                          &
        DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, &
        RHOH_p, RHOQ_tp  )

      call AtmosPhyTbVars_GetLocalMeshFields_tend( n,               &
        mesh, this%vars%tends_manager,                              &
        tb_MOMX_t, tb_MOMY_t, tb_MOMZ_t, tb_RHOT_t, tb_RHOQ_t_list, &
        lcmesh                                                      )

      !$omp parallel private(ke, iq)
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
        MOMX_tp%val(:,ke) = MOMX_tp%val(:,ke) + tb_MOMX_t%val(:,ke)
        MOMY_tp%val(:,ke) = MOMY_tp%val(:,ke) + tb_MOMY_t%val(:,ke)
        MOMZ_tp%val(:,ke) = MOMZ_tp%val(:,ke) + tb_MOMZ_t%val(:,ke)
        RHOT_tp%val(:,ke) = RHOT_tp%val(:,ke) + tb_RHOT_t%val(:,ke)
      end do
      !$omp end do    
      do iq = 1, QA
        if ( .not. TRACER_ADVC(iq) ) cycle
        !$omp do
        do ke = lcmesh%NeS, lcmesh%NeE
          RHOQ_tp(iq)%ptr%val(:,ke) = RHOQ_tp(iq)%ptr%val(:,ke)        &
                                    + tb_RHOQ_t_list(iq)%ptr%val(:,ke)
        end do
        !$omp end do
      end do
      !$omp end parallel
    end do
    call PROF_rapend('ATM_PHY_TB_add_tend', 2)

    return  
  end subroutine AtmosPhyTb_calc_tendency

!> Update variables in a turbulent model
!!
!! @param model_mesh a object to manage computational mesh of atmospheric model 
!! @param prgvars_list a object to mange prognostic variables with atmospheric dynamical core
!! @param trcvars_list a object to mange auxiliary variables 
!! @param forcing_list a object to mange forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL
  subroutine AtmosPhyTb_update( this, model_mesh, prgvars_list, trcvars_list, auxvars_list, forcing_list, is_update )
    implicit none

    class(AtmosPhyTb), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------

    return
  end subroutine AtmosPhyTb_update

!> Finalize a component of sub-grid scale turbulent process
!!
!OCL SERIAL
  subroutine AtmosPhyTb_finalize( this )
    implicit none
    class(AtmosPhyTb), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    call this%tb_driver%Final()
    call this%vars%Final()

    return
  end subroutine AtmosPhyTb_finalize

!> Set boundary conditions to a turbulent model
!!
!! @param dyn_bnd A object to manage boundary conditions of dynamical core
!!
!OCL SERIAL
  subroutine AtmosPhyTb_setDynBC( this, dyn_bnd )
    implicit none
    class(AtmosPhyTb), intent(inout) :: this
    type(AtmDynBnd), intent(in), target :: dyn_bnd
    !--------------------------------------------------

    this%dyn_bnd => dyn_bnd

    return
  end subroutine AtmosPhyTb_setDynBC
  
!- private ------------------------------------------------

end module mod_atmos_phy_tb