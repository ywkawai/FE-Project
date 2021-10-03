!> module ATMOSPHERE phyics / cloud microphysics
!!
!! @par Description
!!          Module for cloud microphysics
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_mp
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

  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
    
  use scale_meshfield_base, only: MeshFieldBase
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_component_proc, only:  ModelComponentProc
  
  use mod_atmos_phy_mp_vars, only: AtmosPhyMpVars

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, extends(ModelComponentProc), public :: AtmosPhyMp
    integer :: MP_TYPEID
    type(AtmosPhyMpVars) :: vars
  contains
    procedure, public :: setup => AtmosPhyMp_setup 
    procedure, public :: calc_tendency => AtmosPhyMp_calc_tendency
    procedure, public :: update => AtmosPhyMp_update
    procedure, public :: finalize => AtmosPhyMp_finalize
  end type AtmosPhyMp

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------

  integer, parameter :: MP_TYPEID_KESSLER = 1

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
  subroutine AtmosPhyMp_setup( this, model_mesh, tm_parent_comp )
    !-----------------------------
    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_regist
    use scale_atmos_phy_mp_kessler, only: &
      ATMOS_PHY_MP_KESSLER_ntracers,            &
      ATMOS_PHY_MP_KESSLER_nwaters,             &
      ATMOS_PHY_MP_KESSLER_nices,               &
      ATMOS_PHY_MP_KESSLER_tracer_names,        &
      ATMOS_PHY_MP_KESSLER_tracer_descriptions, &
      ATMOS_PHY_MP_KESSLER_tracer_units
    !------------------------------
    use scale_time_manager, only: TIME_manager_component    
    use mod_atmos_mesh, only: AtmosMesh
    use mod_atmos_vars, only: ATMOS_PROGVARS_NUM

    implicit none
    class(AtmosPhyMp), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  

    character(len=H_MID) :: MP_TYPE = ''
    namelist /PARAM_ATMOS_PHY_MP/ &
      TIME_DT,          &
      TIME_DT_UNIT,     &
      MP_TYPE
    
    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh

    integer :: ierr

    integer :: QS_MP, QE_MP, QA_MP
    !--------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_MP_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_PHY_MP_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP)

    !- get mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select

    !--- Regist this compoent in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_MP', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                       ! (out) 

    !--- Set the type of turbulence scheme

    select case( MP_TYPE )
    case ('KESSLER')
      this%MP_TYPEID = MP_TYPEID_KESSLER
      call ATMOS_HYDROMETEOR_regist( &
           ATMOS_PHY_MP_KESSLER_nwaters,                & ! (in)
           ATMOS_PHY_MP_KESSLER_nices,                  & ! (in)
           ATMOS_PHY_MP_KESSLER_tracer_names(:),        & ! (in)
           ATMOS_PHY_MP_KESSLER_tracer_descriptions(:), & ! (in)
           ATMOS_PHY_MP_KESSLER_tracer_units(:),        & ! (in)
           QS_MP                                        ) ! (out)
      QA_MP = ATMOS_PHY_MP_KESSLER_ntracers
  
    case default
      LOG_ERROR("ATMOS_PHY_MP_setup",*) 'Not appropriate names of MP_TYPE in namelist PARAM_ATMOS_PHY_MP. Check!'
      call PRC_abort
    end select

    QE_MP = QS_MP + QA_MP - 1

    !- initialize the variables 
    call this%vars%Init( model_mesh, QS_MP, QE_MP, QA_MP )     

    return
  end subroutine AtmosPhyMp_setup

!OCL SERIAL
  subroutine AtmosPhyMp_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars, &
      AtmosVars_GetLocalMeshPhyTends
    use mod_atmos_phy_mp_vars, only:          &
      AtmosPhyMPVars_GetLocalMeshFields_tend,     &
      AtmosPhyMpVars_GetLocalMeshFields_tend_qtrc
    
    implicit none
    
    class(AtmosPhyMp), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    class(MeshBase), pointer :: mesh
    class(LocalMesh3D), pointer :: lcmesh

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_P
    class(LocalMeshFieldBase), pointer :: mp_DENS_t, mp_MOMX_t, mp_MOMY_t, mp_MOMZ_t, mp_RHOT_t, mp_RHOQ_t

    integer :: n
    integer :: ke
    !--------------------------------------------------

    if (.not. this%IsActivated()) return

    !LOG_INFO('AtmosDyn_tendency',*)

    LOG_PROGRESS(*) 'atmosphere / physics / cloud microphysics'

    call model_mesh%GetModelMesh( mesh )

    do n=1, mesh%LOCAL_MESH_NUM
      call PROF_rapstart('ATM_PHY_MP_get_localmesh_ptr', 2)         
      call AtmosVars_GetLocalMeshPrgVars( n,  &
        mesh, prgvars_list, auxvars_list,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,       &
        DENS_hyd, PRES_hyd, lcmesh            )
      call AtmosVars_GetLocalMeshPhyAuxVars( n,      &
        mesh, auxvars_list,                          &
        PRES, PT )
      
      call AtmosVars_GetLocalMeshPhyTends( n,        &
        mesh, forcing_list,                          &
        DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, &
        RHOH_p  )

      call AtmosPhyMpVars_GetLocalMeshFields_tend( n, &
        mesh, this%vars%tends_manager,                                    &
        mp_DENS_t, mp_MOMX_t, mp_MOMY_t, mp_MOMZ_t, mp_RHOT_t, mp_RHOQ_t  )

      call PROF_rapend('ATM_PHY_MP_get_localmesh_ptr', 2)   

      call PROF_rapstart('ATM_PHY_MP_cal_tend', 2)
      ! if (is_update) then

      !   select case( this%MP_TYPEID )
      !   case ( MP_TYPEID_KESSLER )
      !   end select

      ! end if
      
      ! !$omp parallel do
      ! do ke=lcmesh%NeS, lcmesh%NeE
      !   MOMX_tp%val(:,ke) = MOMX_tp%val(:,ke) + mp_MOMX_t%val(:,ke)
      !   MOMY_tp%val(:,ke) = MOMY_tp%val(:,ke) + mp_MOMY_t%val(:,ke)
      !   MOMZ_tp%val(:,ke) = MOMZ_tp%val(:,ke) + mp_MOMZ_t%val(:,ke)
      !   RHOT_tp%val(:,ke) = RHOT_tp%val(:,ke) + mp_RHOT_t%val(:,ke)
      ! end do
      call PROF_rapend('ATM_PHY_MP_cal_tend', 2)
    end do

    return  
  end subroutine AtmosPhyMp_calc_tendency

!OCL SERIAL  
  subroutine AtmosPhyMp_update( this, model_mesh, &
    prgvars_list, trcvars_list,                   &
    auxvars_list, forcing_list, is_update         )  
    
    implicit none
    class(AtmosPhyMp), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------

    return
  end subroutine AtmosPhyMp_update

!OCL SERIAL  
  subroutine AtmosPhyMp_finalize( this )
    implicit none
    class(AtmosPhyMp), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    ! select case ( this%MP_TYPEID )
    ! case( MP_TYPEID_KESSLER )
    ! end select

    call this%vars%Final()

    return
  end subroutine AtmosPhyMp_finalize
  
!- private ------------------------------------------------

end module mod_atmos_phy_mp