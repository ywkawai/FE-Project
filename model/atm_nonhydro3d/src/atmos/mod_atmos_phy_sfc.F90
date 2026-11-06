!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / surface process
!!
!! @par Description
!!          Module for surface process
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_sfc
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
  use scale_atmos_hydrometeor, only: &
    ATMOS_HYDROMETEOR_dry
      
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
  
  use mod_atmos_mesh, only: AtmosMesh
  use mod_atmos_phy_sfc_vars, only: AtmosPhySfcVars

  use mod_atmos_vars_container, only: AtmosVarsContainer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage a component of surface process
  !!
  type, extends(ModelComponentProc), public :: AtmosPhySfc
    class(AtmosMesh), pointer :: mesh    !< Pointer to a object to manage the mesh with atmospheric component

    integer :: atm_var_container_typeid !< Type ID of variable container for surface process

    integer :: SFCFLX_TYPEID            !< Type id of surface scheme
    type(AtmosPhySfcVars) :: vars       !< An object to manage variables with surface component

    logical :: CPL_sw  !< Switch to indicate if the coupler is activated

    type(AtmosVarsContainer), pointer :: primary_atmvars_container !< Pointer to the primary atmospheric variable container
  contains
    procedure, public :: setup => AtmosPhySfc_setup
    procedure, public :: set_coupler_flag => AtmosPhySfc_set_coupler_flag
    procedure, public :: calc_tendency => AtmosPhySfc_calc_tendency
    procedure, public :: update => AtmosPhySfc_update
    procedure, public :: finalize => AtmosPhySfc_finalize
    procedure, public :: set_primary_atmvars_container => AtmosPhySfc_set_primary_atmvars_container
  end type AtmosPhySfc

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------

  integer, parameter :: SFCFLX_TYPEID_CONST           = 1 !< Type ID of constant surface flux scheme
  integer, parameter :: SFCFLX_TYPEID_SIMPLE          = 2 !< Type ID of simple bulk surface flux scheme

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  private :: cal_tend_from_sfcflx
  private :: cal_del_flux
  private :: convert_UV2LocalOrthVec
  private :: convert_LocalOrth2UVVec

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
    
contains

!> Setup an object to manage a component of surface process
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param tm_parent_comp Object to mange a temporal scheme in a parent component
!!
!OCL SERIAL
  subroutine AtmosPhySfc_setup( this, model_mesh, tm_parent_comp )
    use scale_time_manager, only: TIME_manager_component

    use scale_atmos_phy_sf_const, only: &
       ATMOS_PHY_SF_const_setup
    use scale_atm_phy_sf_bulk_simple, only: &
       ATMOS_PHY_SF_simple_setup

    use mod_atmos_vars, only: ATM_VARS_CONTAINER_PRIMARY_ID

    implicit none
    class(AtmosPhySfc), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                       = UNDEF8   !< Timestep for surface process
    character(len=H_SHORT) :: TIME_DT_UNIT    = 'SEC'    !< Unit of timestep
    character(len=H_MID) :: SFCFLX_TYPE       = "CONST"  !< Type of surface flux scheme
    integer :: atm_var_container_typeid                  !< Type ID of variable container for surface process
    real(RP) :: DEFAULT_SFC_TEMP              = 300.0_RP !< Default value of surface temperature
    namelist /PARAM_ATMOS_PHY_SFC/ &
      TIME_DT,                  &
      TIME_DT_UNIT,             &
      SFCFLX_TYPE,              &
      atm_var_container_typeid, &
      DEFAULT_SFC_TEMP
    
    integer :: ierr

    integer :: ldomID
    class(LocalMesh2D), pointer :: lcmesh2D
    integer :: ke
    !--------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SFC_setup",*) 'Setup'

    atm_var_container_typeid = ATM_VARS_CONTAINER_PRIMARY_ID

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SFC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_SFC_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_PHY_SFC_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_SFC. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_SFC)

    this%atm_var_container_typeid = atm_var_container_typeid

    !--- Regist this compoent in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_SFC', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                        ! (out)

    !- initialize the variables 
    call this%vars%Init( model_mesh, DEFAULT_SFC_TEMP )      

    !--- Set the type of surface flux scheme

    select case( SFCFLX_TYPE )
    case ('CONST')
      call ATMOS_PHY_SF_const_setup()
      this%SFCFLX_TYPEID = SFCFLX_TYPEID_CONST
    case ('SIMPLE')
      call ATMOS_PHY_SF_simple_setup()
      this%SFCFLX_TYPEID = SFCFLX_TYPEID_SIMPLE
    case default
      LOG_ERROR("ATMOS_PHY_SFC_setup",*) 'Not appropriate names of SFCFLX_TYPE in namelist PARAM_PHY_SFC. Check!'
      call PRC_abort
    end select

    !--
    select type(model_mesh)
    class is (AtmosMesh)
      this%mesh => model_mesh
    end select

    return
  end subroutine AtmosPhySfc_setup

  !> Set switch to indicate if the coupler is activated
  !!
  subroutine AtmosPhySfc_set_coupler_flag( this, cpl_sw )
    implicit none
    class(AtmosPhySfc), intent(inout) :: this
    logical, intent(in) :: cpl_sw
    !--------------------------------------------------
    this%CPL_sw = cpl_sw
    return
  end subroutine AtmosPhySfc_set_coupler_flag

  !> Set a pointer to the primary atmospheric variable container
  !!
  subroutine AtmosPhySfc_set_primary_atmvars_container( this, primary_atmvars_container )
    implicit none
    class(AtmosPhySfc), intent(inout) :: this
    type(AtmosVarsContainer), intent(in), target :: primary_atmvars_container
    !--------------------------------------------------

    this%primary_atmvars_container => primary_atmvars_container
    return
  end subroutine AtmosPhySfc_set_primary_atmvars_container

  !> Calculate tendencies associated with a surface model
  !!
  !! @param model_mesh Object to manage computational mesh of atmospheric model 
  !! @param prgvars_list Object to mange prognostic variables with atmospheric dynamical core
  !! @param trcvars_list Object to mange auxiliary variables 
  !! @param forcing_list Object to mange forcing terms
  !! @param is_update Flag to specify whether the tendencies are updated in this call
  !!
!OCL SERIAL
  subroutine AtmosPhySfc_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )

    use scale_tracer, only: &
      TRACER_inq_id

    use scale_meshfield_base, only: MeshField3D

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyTends,   &
      AtmosVars_GetLocalMeshPhyAuxVars, &
      AtmosVars_GetLocalMeshQTRC_Qv,    &
      AtmosVars_GetLocalMeshQTRCPhyTend
    use mod_atmos_phy_sfc_vars, only: &
      AtmosPhySfcVars_GetLocalMeshFields

    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      DENS_VID => PRGVAR_DDENS_ID, MOMX_VID => PRGVAR_MOMX_ID, &
      MOMY_VID => PRGVAR_MOMY_ID,  MOMZ_VID => PRGVAR_MOMZ_ID, &
      PRES_VID => AUXVAR_PRES_ID   
    implicit none
    
    class(AtmosPhySfc), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    class(MeshBase), pointer :: mesh
    class(LocalMesh3D), pointer :: lcmesh

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT, QV
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p, RHOQv_tp
    class(LocalMeshFieldBase), pointer :: SFLX_MU, SFLX_MV, SFLX_MW, SFLX_SH, SFLX_LH, SFLX_QV
    class(LocalMeshFieldBase), pointer :: SFC_TEMP

    class(MeshField3D), pointer :: DDENS_c, MOMX_c, MOMY_c, MOMZ_c, PRES_c

    integer :: n
    integer :: iq
    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    call model_mesh%GetModelMesh( mesh )

    if (is_update) then

      call prgvars_list%Get3D( MOMZ_VID, MOMZ_c  )
      call prgvars_list%Get3D( MOMX_VID, MOMX_c  )
      call prgvars_list%Get3D( MOMY_VID, MOMY_c  )
      call prgvars_list%Get3D( DENS_VID, DDENS_c )
      call auxvars_list%Get3D( PRES_VID, PRES_c )

      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart( 'ATM_PHY_SFC_get_localmesh_ptr', 2)         
        call AtmosVars_GetLocalMeshPrgVars( n,  &
          mesh, this%primary_atmvars_container%PROGVARS_manager, &
          this%primary_atmvars_container%AUXVARS_manager,        &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                        &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot,                &
          lcmesh                                                 )
        call AtmosVars_GetLocalMeshPhyTends( n,        &
          mesh, forcing_list,                          &
          DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, &
          RHOH_p                                       )
        call AtmosVars_GetLocalMeshQTRC_Qv( n, mesh,                     &
          this%primary_atmvars_container%QTRCVARS_manager, forcing_list, &
          QV, RHOQv_tp                                                   ) 
        call AtmosVars_GetLocalMeshPhyAuxVars( n, &
          mesh, auxvars_list,                     &
          PRES, PT                                )

        call AtmosPhySfcVars_GetLocalMeshFields( n,                      &
          mesh, this%vars%SFCVARS_manager, this%vars%SFCFLX_manager,     &
          SFC_TEMP, SFLX_MU, SFLX_MV, SFLX_MW, SFLX_SH, SFLX_LH, SFLX_QV )
        call PROF_rapend( 'ATM_PHY_SFC_get_localmesh_ptr', 2)   

        call PROF_rapstart( 'ATM_PHY_SFC_cal_tend', 2)
        call cal_tend_from_sfcflx( this, is_update,                                     &
          DENS_tp%val, MOMX_tp%val, MOMY_tp%val, MOMZ_tp%val, RHOH_p%val, RHOQv_tp%val, &
          SFLX_MU%val, SFLX_MV%val, SFLX_MW%val, SFLX_SH%val, SFLX_LH%val, SFLX_QV%val, &
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val, QV%val,                   &
          DDENS_c%local(n)%val, MOMX_c%local(n)%val, MOMY_c%local(n)%val,               &
          MOMZ_c%local(n)%val, PRES_c%local(n)%val,                                     &          
          DENS_hyd%val, PRES_hyd%val,                                                   &
          PRES%val, PT%val, Rtot%val, CVtot%val, CPtot%val, SFC_TEMP%val,               & 
          model_mesh%DOptrMat(3), model_mesh%SOptrMat(3), model_mesh%LiftOptrMat,       &
          lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D          )
        call PROF_rapend( 'ATM_PHY_SFC_cal_tend', 2)
      end do

    end if

    return  
  end subroutine AtmosPhySfc_calc_tendency

!> Update variables in a surface model
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param prgvars_list Object to mange prognostic variables with atmospheric dynamical core
!! @param trcvars_list Object to mange auxiliary variables 
!! @param forcing_list Object to mange forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL
  subroutine AtmosPhySfc_update( this, model_mesh,           &
    prgvars_list, trcvars_list, auxvars_list,  forcing_list, &
    is_update )
    implicit none

    class(AtmosPhySfc), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------

    return
  end subroutine AtmosPhySfc_update

!> Finalize a component of surface process
!!
!OCL SERIAL
  subroutine AtmosPhySfc_finalize( this )
    implicit none
    class(AtmosPhySfc), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    call this%vars%Final()

    return
  end subroutine AtmosPhySfc_finalize

!- private ------------------------------------------------

  !> Calculate tendencies associated with a surface model
!OCL SERIAL  
  subroutine cal_tend_from_sfcflx( this, is_update_sflx, &
    DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOH_p, RHOQV_tp, &
    SFLX_MU, SFLX_MV, SFLX_MW, SFLX_SH, SFLX_LH, SFLX_QV, &
    DDENS, MOMX, MOMY, MOMZ, DRHOT, QV,                   &
    DDENS_c, MOMX_c, MOMY_c, MOMZ_c, PRES_c,              &
    DENS_hyd, PRES_hyd, PRES, PT, Rtot, CVtot, CPtot,     &
    SFC_TEMP,                                             & 
    Dz, Sz, Lift,                                         &
    lcmesh, elem, lcmesh2D, elem2D )

    use scale_const, only: &
      Rdry => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00, &
      RPlanet => CONST_RADIUS    
    use scale_atmos_hydrometeor, only: &
      LHV, &
      I_QV
    use scale_tracer, only: &
      TRACER_CV, TRACER_CP, &
      TRACER_MASS,          &
      TRACER_ENGI0,         &
      QA
    use scale_atmos_phy_sf_const, only: &
       ATMOS_PHY_SF_const_flux
    use scale_atm_phy_sf_bulk_simple, only: &
       ATMOS_PHY_SF_simple_flux
       
    use scale_sparsemat, only: &
      SparseMat, &
      sparsemat_matmul    
    implicit none

    class(AtmosPhySfc), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    logical, intent(in) :: is_update_sflx
    real(RP), intent(inout) :: DENS_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMZ_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOH_p (elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOQV_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: SFLX_MU(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(inout) :: SFLX_MV(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(inout) :: SFLX_MW(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(inout) :: SFLX_SH(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(inout) :: SFLX_LH(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(inout) :: SFLX_QV(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: QV(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS_c(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX_c(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY_c(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ_c(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_c(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PT(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem%Np,lcmesh%NeA)    
    real(RP), intent(in) :: SFC_TEMP(elem2D%Np,lcmesh2D%NeA)
    type(SparseMat), intent(in) :: Dz, Sz, Lift

    real(RP) :: ATM_W   (elem2D%Np,lcmesh2D%Ne)
    real(RP) :: ATM_U   (elem2D%Np,lcmesh2D%Ne)
    real(RP) :: ATM_V   (elem2D%Np,lcmesh2D%Ne)
    real(RP) :: ATM_TEMP(elem2D%Np,lcmesh2D%Ne)
    real(RP) :: ATM_PRES(elem2D%Np,lcmesh2D%Ne)
    real(RP) :: ATM_QV  (elem2D%Np,lcmesh2D%Ne)
    real(RP) :: DZ1     (elem2D%Np,lcmesh2D%Ne)
    real(RP) :: Z1      (elem2D%Np,lcmesh2D%Ne)
    real(RP) :: SFC_PRES(elem2D%Np,lcmesh2D%Ne)
    real(RP) :: SFC_DENS(elem2D%Np,lcmesh2D%Ne)

    real(RP) :: SFLX_QTRC(elem2D%Np,lcmesh2D%Ne,max(1,QA))

    ! dummy
    real(RP) :: U10(elem2D%Np,lcmesh2D%Ne)
    real(RP) :: V10(elem2D%Np,lcmesh2D%Ne)

    integer :: ke
    integer :: ke2D
    integer :: hsliceZ0, hsliceZ1
    integer :: ij
    integer :: iq
    integer :: p

    real(RP) :: dens
    real(RP) :: LiftDelFlx(elem%Np), LiftDelFlx_c(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lcmesh%Ne,6)

    real(RP) :: TEMP(elem%Np)
    real(RP) :: RHOCV_t(elem%Np), RHOCP_t(elem%Np)
    real(RP) :: CP_t_x_dens, CV_t_x_dens
    real(RP) :: DENS_t_SF(elem%Np)
    real(RP) :: work

    real(RP) :: SFLX_RHOH_trc(elem%Np,lcmesh2D%Ne)
    integer :: ps, pe    

    real(RP) :: rPRE00
    !-------------------------------------------------

    if ( is_update_sflx .and. (.not. this%CPL_sw) ) then

      !$omp parallel do collapse(2) private( &
      !$omp ke, hSliceZ0, hsliceZ1,          &
      !$omp dens                             )
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
      do ij=1, elem2D%Np

        ke = ke2D
        hsliceZ0 = elem%Hslice(ij,1)
        hsliceZ1 = elem%Hslice(ij,2)

        dens = DENS_hyd(hsliceZ1,ke) + DDENS_c(hsliceZ1,ke)
        ATM_U(ij,ke2D) = MOMX_c(hsliceZ1,ke) / dens
        ATM_V(ij,ke2D) = MOMY_c(hsliceZ1,ke) / dens
        ATM_W(ij,ke2D) = 0.0_RP
        ATM_PRES(ij,ke2D) = PRES(hsliceZ1,ke)
!        ATM_TEMP(ij,ke2D) = ATM_PRES(ij,ke2D) / ( Rtot(hsliceZ1,ke) * dens )

        SFC_DENS(ij,ke2D) = DENS_hyd(hsliceZ0,ke) + DDENS(hsliceZ0,ke)

        ! Tentative: 
        ATM_TEMP(ij,ke2D) = PRES(hsliceZ0,ke) / ( Rtot(hsliceZ0,ke) * SFC_DENS(ij,ke2D) )
        ATM_QV(ij,ke2D) = QV(hsliceZ0,ke)
        SFC_DENS(ij,ke2D) = DENS_hyd(hsliceZ0,ke) + DDENS_c(hsliceZ0,ke)
        SFC_PRES(ij,ke2D) = PRES_c(hsliceZ0,ke)

        Z1(ij,ke2D) = lcmesh%zlev(hsliceZ1,ke)
        DZ1(ij,ke2D) = Z1(ij,ke2D) - lcmesh%zlev(hsliceZ0,ke)
        Z1(ij,ke2D) = Z1(ij,ke2D) + RPlanet
      end do
      end do

      call convert_UV2LocalOrthVec( &
        this%mesh%ptr_mesh, lcmesh2D%pos_en(:,:,1), lcmesh2D%pos_en(:,:,2), Z1(:,:), elem2D%Np*lcmesh2D%Ne, &
        ATM_U, ATM_V ) ! (inout)
      
      !---
      select case ( this%SFCFLX_TYPEID )
      case ( SFCFLX_TYPEID_CONST )
        call ATMOS_PHY_SF_const_flux( &
          elem2D%Np, 1, elem2D%Np, lcmesh2D%NeA, 1, lcmesh2D%Ne, & ! (in)
          ATM_W(:,:), ATM_U(:,:), ATM_V(:,:), SFC_TEMP(:,:),     & ! (in)
          DZ1(:,:), SFC_DENS(:,:),                               & ! (in)
          SFLX_MW(:,:), SFLX_MU(:,:), SFLX_MV(:,:),              & ! (out)
          SFLX_SH(:,:), SFLX_LH(:,:), SFLX_QV(:,:),              & ! (out)
          U10(:,:), V10(:,:)                                     ) ! (out)
      case ( SFCFLX_TYPEID_SIMPLE )
        call ATMOS_PHY_SF_simple_flux( &
          elem2D%Np, 1, elem2D%Np, lcmesh2D%NeA, 1, lcmesh2D%Ne,            & ! (in)
          ATM_W(:,:), ATM_U(:,:), ATM_V(:,:), ATM_TEMP(:,:), ATM_PRES(:,:), & ! (in) Note: ATM_PRES is not used
          ATM_QV(:,:),                                                      & ! (in)
          SFC_DENS(:,:), SFC_TEMP(:,:), SFC_PRES(:,:),                      & ! (in)
          DZ1(:,:),                                                         & ! (in)
          SFLX_MW(:,:), SFLX_MU(:,:), SFLX_MV(:,:),                         & ! (out)
          SFLX_SH(:,:), SFLX_LH(:,:), SFLX_QV(:,:),                         & ! (out)
          U10(:,:), V10(:,:)                                                ) ! (out)
      end select

      call convert_LocalOrth2UVVec( &
        this%mesh%ptr_mesh, lcmesh2D%pos_en(:,:,1), lcmesh2D%pos_en(:,:,2), Z1(:,:), elem2D%Np*lcmesh2D%Ne, &
        SFLX_MU(:,lcmesh2D%NeS:lcmesh2D%NeE), SFLX_MV(:,lcmesh2D%NeS:lcmesh2D%NeE) ) ! (inout)

    end if

    !- Add the tendency due to the surface flux

    call cal_del_flux( del_flux,            &
      SFLX_MU, SFLX_MV, SFLX_MW, SFLX_SH,   &
      SFLX_QV,                              &
      SFC_TEMP,                             &
      lcmesh%normal_fn(:,:,3),              &
      lcmesh, elem, lcmesh2D, elem2D        )
    
    if ( .not. ATMOS_HYDROMETEOR_dry ) then
    
      !$omp parallel do private(ke, iq, LiftDelFlx,    &
      !$omp RHOCP_t, RHOCV_t, DENS_t_SF, work, &
      !$omp CP_t_x_dens, CV_t_x_dens )
      do ke2D=1, lcmesh2D%Ne
        ke = ke2D
        RHOCP_t(:) = 0.0_RP; RHOCV_t(:) = 0.0_RP
        DENS_t_SF(:) = 0.0_RP

        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,6), LiftDelFlx )
        SFLX_RHOH_trc(:,ke) = - LiftDelFlx(:)

!       do iq=1, QA
          iq = I_QV
          SFLX_QTRC(:,ke2D,iq) = SFLX_QV(:,ke2D)

          call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,5), LiftDelFlx )
          RHOQv_tp(:,ke) = RHOQv_tp(:,ke) - LiftDelFlx(:)

          do p=1, elem%Np
            work = - LiftDelFlx(p)
            DENS_t_SF(p) = DENS_t_SF(p) + work * TRACER_MASS(iq)
            RHOCP_t(p) = RHOCP_t(p) + work * TRACER_CP(iq)
            RHOCV_t(p) = RHOCV_t(p) + work * TRACER_CV(iq)
          end do
!       end do

        DENS_tp(:,ke) = DENS_tp(:,ke) + DENS_t_SF(:)
        do p=1, elem%Np
          CP_t_x_dens = RHOCP_t(p) - CPtot(p,ke) * DENS_t_SF(p)  ! DENS d(CP) = d(RHOCP) - CP d(DENS)
          CV_t_x_dens = RHOCV_t(p) - CVtot(p,ke) * DENS_t_SF(p)  ! DENS d(CV) = d(RHOCV) - CV d(DENS)

          work = PRES(p,ke) / ( Rtot(p,ke) * ( DENS_hyd(p,ke) + DDENS(p,ke) ) )
          SFLX_RHOH_trc(p,ke) = SFLX_RHOH_trc(p,ke) &
            - ( CP_t_x_dens + log(PRES(p,ke) * RPRE00 ) * ( CVtot(p,ke) / CPtot(p,ke) * CP_t_x_dens - CV_t_x_dens ) ) * work
        end do
      end do
    end if

    !$omp parallel do private(ke, iq, LiftDelFlx)
    do ke2D=1, lcmesh2D%Ne
      ke = ke2D

      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,1), LiftDelFlx )
      MOMX_tp(:,ke) = MOMX_tp(:,ke) - LiftDelFlx(:) 

      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,2), LiftDelFlx )
      MOMY_tp(:,ke) = MOMY_tp(:,ke) - LiftDelFlx(:) 

      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,3), LiftDelFlx )
      MOMZ_tp(:,ke) = MOMZ_tp(:,ke) - LiftDelFlx(:) 
      
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,4), LiftDelFlx )
      RHOH_p (:,ke) = RHOH_p (:,ke) - LiftDelFlx(:) 

      if ( .not. ATMOS_HYDROMETEOR_dry ) then
        RHOH_p(:,ke) = RHOH_p(:,ke) + SFLX_RHOH_trc(:,ke)
      end if 
    end do

    return
  end subroutine cal_tend_from_sfcflx

!OCL SERIAL  
  subroutine cal_del_flux( del_flux,      &
    sflx_mu, sflx_mv, sflx_mw, sflx_sh,   &
    sflx_qv,                              &
    SFC_TEMP,                             &
    nz,                                   &
    lmesh, elem, lmesh2D, elem2D          )

    use scale_atmos_hydrometeor, only: &
      LHV, &
      I_QV
    use scale_tracer, only: &
      QA, TRACER_CV, TRACER_ENGI0  
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D   
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne,6)
    real(RP), intent(in) :: sflx_mu(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: sflx_mv(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: sflx_mw(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: sflx_sh(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: sflx_qv(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: SFC_TEMP(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)

    integer :: ke2D, p
    integer :: i
    integer :: iq
    real(RP) :: sflx_dens    
    !------------------------------------------------------------------------

    !$omp parallel private(i, iq, sflx_dens)
    !$omp workshare
    del_flux(:,:) = 0.0_RP
    !$omp end workshare

    !$omp do collapse(2)
    do ke2D=1, lmesh2D%Ne
    do p=1, elem2D%Np
      i = elem%Nfaces_h*elem%Nfp_h + p + (ke2D-1)*elem%NfpTot
      del_flux(i,1) = sflx_mu(p,ke2D) * nz(i)
      del_flux(i,2) = sflx_mv(p,ke2D) * nz(i)
      del_flux(i,3) = sflx_mw(p,ke2D) * nz(i)
      del_flux(i,4) = sflx_sh(p,ke2D) * nz(i)
      del_flux(i,5) = sflx_qv(p,ke2D) * nz(i)
    end do
    end do
    !$omp end do

    if ( .not. ATMOS_HYDROMETEOR_dry ) then
      !$omp do collapse(2)
      do ke2D=1, lmesh2D%Ne
      do p=1, elem2D%Np
        i = elem%Nfaces_h*elem%Nfp_h + p + (ke2D-1)*elem%NfpTot
        del_flux(i,6) = sflx_qv(p,ke2D) * ( TRACER_CV(I_QV) * SFC_TEMP(p,ke2D) + LHV ) * nz(i)

        ! do iq=1, QA
          iq = I_QV
          sflx_dens = sflx_qv(p,ke2D)

          del_flux(i,6) = del_flux(i,6) &
            - TRACER_ENGI0(iq) * sflx_dens * nz(i)
        ! end do          
      end do
      end do
    end if     
    !$omp end parallel

    return
  end subroutine cal_del_flux

!OCL SERIAL
  subroutine convert_UV2LocalOrthVec( mesh, x1, x2, r, N, & ! (in)
    U, V ) ! (inout)
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LocalOrthVec_alpha
    use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D

    implicit none
    integer, intent(in) :: N
    class(MeshBase3D), intent(in) :: mesh
    real(RP), intent(in) :: x1(N)
    real(RP), intent(in) :: x2(N)
    real(RP), intent(in) :: r(N)
    real(RP), intent(inout) :: U(N)
    real(RP), intent(inout) :: V(N)

    select type(mesh)
    type is (MeshCubedSphereDom3D)
      call CubedSphereCoordCnv_CS2LocalOrthVec_alpha( &
        x1, x2, r, N, & ! (in)
        U, V          ) ! (inout)
    end select
    return
  end subroutine convert_UV2LocalOrthVec

!OCL SERIAL
  subroutine convert_LocalOrth2UVVec( mesh, x1, x2, r, N, & ! (in)
    U, V ) ! (inout)
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LocalOrth2CSVec_alpha
    use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D

    implicit none
    integer, intent(in) :: N
    class(MeshBase3D), intent(in) :: mesh
    real(RP), intent(in) :: x1(N)
    real(RP), intent(in) :: x2(N)
    real(RP), intent(in) :: r(N)
    real(RP), intent(inout) :: U(N)
    real(RP), intent(inout) :: V(N)

    select type(mesh)
    type is (MeshCubedSphereDom3D)
      call CubedSphereCoordCnv_LocalOrth2CSVec_alpha( &
        x1, x2, r, N, & ! (in)
        U, V          ) ! (inout)
    end select
    return
  end subroutine convert_LocalOrth2UVVec 
end module mod_atmos_phy_sfc