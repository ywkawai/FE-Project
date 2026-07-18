!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Radiation
!!
!! @par Description
!!          Module for radiation component in atmospheric model
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_rd
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
  use scale_element_base, only: ElementBase, &
    ElementBase1D, ElementBase2D, ElementBase3D

  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField2D, MeshField3D
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_component_proc, only:  ModelComponentProc

  use mod_atmos_phy_rd_vars, only: AtmosPhyRdVars
  use mod_atmos_vars_container, only: &
    AtmosVarsContainer

  use scale_atm_phy_rd_dgm_simple, only: AtmPhyRadSimple

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage a component of radiation in atmospheric model
  !!
  type, extends(ModelComponentProc), public :: AtmosPhyRd
    integer :: RD_TYPEID          !< Type id of radiation scheme
    integer :: RD_SOLARINS_TYPEID !< Type id of solar insolation scheme

    type(AtmosPhyRdVars) :: vars               !< Object to manage variables with radiation
    type(MeshField2D), pointer :: SFC_TEMP_ptr !< Pointer to an object with surface temperature field
    type(MeshField2D), pointer :: SFC_ALB_ptr  !< Pointer to an object with surface albedo field

    integer :: atm_var_container_typeid        !< Type ID of variable container for radiation

    !-
    type(AtmPhyRadSimple) :: simple_rad !< Object to manage a simplified-radiation scheme
  contains
    procedure, public :: setup => AtmosPhyRd_setup
    procedure, public :: calc_tendency => AtmosPhyRd_calc_tendency
    procedure, public :: update => AtmosPhyRd_update
    procedure, public :: finalize => AtmosPhyRd_finalize
    procedure, public :: SetSfcVars => AtmosPhyRd_set_SfcVars
    procedure, private :: calc_tendency_core => AtmosPhyRd_calc_tendency_core
  end type AtmosPhyRd

  !-----------------------------------------------------------------------------
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
  integer, parameter :: RD_TYPEID_SIMPLE  = 1 !< Type ID of a simple radiation scheme

  integer, parameter :: RD_INSOLATION_TYPEID_NONE   = 0 !< Type ID of no solar insolation
  integer, parameter :: RD_INSOLATION_TYPEID_SIMPLE = 1 !< Type ID of a simple solar insolation scheme
contains

!> Setup a component of radiation in atmospheric model
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param tm_parent_comp Object to mange a temporal scheme in a parent component
!!
  subroutine AtmosPhyRd_setup( this, model_mesh, tm_parent_comp )
    use mod_atmos_mesh, only: AtmosMesh
    use scale_time_manager, only: TIME_manager_component
    use mod_atmos_vars, only: ATM_VARS_CONTAINER_PRIMARY_ID
    use scale_atm_phy_rd_solarins_simple, only: &
      atm_phy_rd_solarins_simple_setup
    implicit none
    class(AtmosPhyRd), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8 !< Timestep for radiation
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  !< Unit of timestep

    character(len=H_MID) :: RD_TYPE          = 'NONE'        !< Type of a radiation scheme [NONE, GRAYRAD]
    character(len=H_MID) :: RD_SOLARINS_TYPE = 'SIMPLE'      !< Type of solar insolation scheme [NONE, SIMPLE, REAL]
    integer :: atm_var_container_typeid

    namelist /PARAM_ATMOS_PHY_RD/ &
      TIME_DT,             &
      TIME_DT_UNIT,        &
      RD_TYPE,             &
      RD_SOLARINS_TYPE,    &
      atm_var_container_typeid

    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D

    integer :: ierr
    !-----------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_setup",*) 'Setup'

    atm_var_container_typeid = ATM_VARS_CONTAINER_PRIMARY_ID

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_RD_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_PHY_RD_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_RD. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_RD)

    this%atm_var_container_typeid = atm_var_container_typeid

    !--- Register this component in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_RD', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                       ! (out) 

    !--- Set the type of radiation scheme
    
    select case ( RD_TYPE )
    case ("SIMPLE")
      this%RD_TYPEID = RD_TYPEID_SIMPLE
      call this%simple_rad%Init()
    case default
      LOG_ERROR("ATMOS_PHY_RD_setup",*) 'Not appropriate RD_TYPE. Check!'
      call PRC_abort
    end select

    select case ( RD_SOLARINS_TYPE )
    case ("SIMPLE")
      this%RD_SOLARINS_TYPEID = RD_INSOLATION_TYPEID_SIMPLE
      call atm_phy_rd_solarins_simple_setup()
    case ("NONE")
      this%RD_SOLARINS_TYPEID = RD_INSOLATION_TYPEID_NONE
    case default
      LOG_ERROR("ATMOS_PHY_RD_setup",*) 'Not appropriate RD_SOLARINS_TYPE. Check!'
      call PRC_abort
    end select


    !- Initialize the variables 
    call this%vars%Init( model_mesh )

    return
  end subroutine AtmosPhyRd_setup

!OCL SERIAL
  subroutine AtmosPhyRd_set_SfcVars( this, sfc_temp, sfc_alb )
    implicit none
    class(AtmosPhyRd), intent(inout) :: this
    type(MeshField2D), target, intent(in) :: sfc_temp
    type(MeshField2D), target, intent(in) :: sfc_alb
    !-----------------------------------------------------
    this%SFC_TEMP_ptr => sfc_temp
    this%SFC_ALB_ptr => sfc_alb
    return
  end subroutine AtmosPhyRd_set_SfcVars

!> Calculate tendencies associated with radiation in atmospheric model
!!
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param prgvars_list Object to manage prognostic variables with atmospheric dynamical core
!! @param trcvars_list Object to manage auxiliary variables 
!! @param forcing_list Object to manage forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL
  subroutine AtmosPhyRd_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )
    use scale_atm_phy_rd_dgm_common, only: ATM_PHY_RD_DGM_calc_heating
    use scale_atm_phy_rd_solarins_simple, only: atm_phy_rd_solarins_simple_get
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,     &
      AtmosVars_GetLocalMeshPhyAuxVars,  &
      AtmosVars_GetLocalMeshQTRC_Qv,     &
      AtmosVars_GetLocalMeshPhyTends
    use mod_atmos_phy_rd_vars, only: &
      RD_RHOH_ID => ATMOS_PHY_RD_RHOH_ID, &
      SOLINS_ID => ATMOS_PHY_RD_AUX2D_SOLINS_ID,         &
      COS_SZA_ID => ATMOS_PHY_RD_AUX2D_COSSZA_ID,        &
      SFLX_SW_up_ID => ATMOS_PHY_RD_AUX2D_SFLX_SW_up_ID, &
      SFLX_SW_dn_ID => ATMOS_PHY_RD_AUX2D_SFLX_SW_dn_ID, &
      SFLX_LW_up_ID => ATMOS_PHY_RD_AUX2D_SFLX_LW_up_ID, &
      SFLX_LW_dn_ID => ATMOS_PHY_RD_AUX2D_SFLX_LW_dn_ID
    implicit none
    class(AtmosPhyRd), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    class(MeshBase), pointer :: mesh
    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMesh2D), pointer :: lcmesh2D

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: QV, QV_tp
    class(LocalMeshFieldBase), pointer :: rd_RHOH
    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p

    integer :: n
    integer :: ke
    !------------------------------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_PROGRESS(*) 'atmosphere / physics / radiation'

    call model_mesh%GetModelMesh( mesh )

    if (is_update) then
    
      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart('ATM_PHY_RD_get_localmesh_ptr', 2)    

        !- Get pointers to the fields in the variable container for cloud microphysics
        call AtmosVars_GetLocalMeshPrgVars( n,  &
          mesh, prgvars_list, auxvars_list,       &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,         &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
          lcmesh                                  )

        call AtmosVars_GetLocalMeshPhyAuxVars( n,  &
          mesh, auxvars_list,                      &
          PRES, PT )

        call AtmosVars_GetLocalMeshQTRC_Qv( n, &
          mesh, trcvars_list, forcing_list,    &
          QV, QV_tp )
        call PROF_rapend('ATM_PHY_RD_get_localmesh_ptr', 2)

        !- Calculate tendencies associated with radiation

        call PROF_rapstart('ATM_PHY_RD_cal_tend', 2)
        
        lcmesh2D => lcmesh%lcmesh2D

        select case ( this%RD_SOLARINS_TYPEID )
        case ( RD_INSOLATION_TYPEID_SIMPLE )
          call atm_phy_rd_solarins_simple_get( &
            this%vars%auxvars2D(SOLINS_ID)%local(n)%val(:,lcmesh2D%NeS:lcmesh2D%NeE),  & ! (out)
            this%vars%auxvars2D(COS_SZA_ID)%local(n)%val(:,lcmesh2D%NeS:lcmesh2D%NeE), & ! (out)
            lcmesh2D%lat, lcmesh2D%Ne * lcmesh2D%refElem2D%Np )
        end select

        call this%vars%tends(RD_RHOH_ID)%GetLocalMeshField( n, rd_RHOH )
        call this%calc_tendency_core( rd_RHOH%val,                & ! (out)
          this%vars%auxvars2D(SFLX_SW_up_ID)%local(n)%val,        & ! (out)
          this%vars%auxvars2D(SFLX_SW_dn_ID)%local(n)%val,        & ! (out)
          this%vars%auxvars2D(SFLX_LW_up_ID)%local(n)%val,        & ! (out)
          this%vars%auxvars2D(SFLX_LW_dn_ID)%local(n)%val,        & ! (out)
          this%vars%auxvars2D(SOLINS_ID)%local(n)%val,            & ! (in)
          this%vars%auxvars2D(COS_SZA_ID)%local(n)%val,           & ! (in)
          DDENS%val, PRES%val, QV%val,                            & ! (in)
          this%SFC_TEMP_ptr%local(n)%val,                         & ! (in)
          this%SFC_ALB_ptr%local(n)%val,                          & ! (in)
          DENS_hyd%val, Rtot%val, CVtot%val,                      & ! (in)
          lcmesh, lcmesh%refElem3D, lcmesh2D, lcmesh2D%refElem2D, & ! (in)
          model_mesh%element3D_operation )                          ! (in)
        call PROF_rapend('ATM_PHY_RD_cal_tend', 2)
      end do
    
    end if

    !- Add tendencies calculated in this component to the total tendencies

    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPhyTends( n,        &
        mesh, forcing_list,                          &
        DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, &
        RHOH_p, lcmesh3D=lcmesh )

      call this%vars%tends(RD_RHOH_ID)%GetLocalMeshField( n, rd_RHOH )
      !$omp parallel private(ke)
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
        RHOH_p %val(:,ke) = RHOH_p %val(:,ke) + rd_RHOH%val(:,ke)
      end do
      !$omp end parallel
    end do

    return
  end subroutine AtmosPhyRd_calc_tendency

!> Update variables in a component of radiation in atmospheric model
!!
!! @param model_mesh Object to manage computational mesh of atmospheric model 
!! @param prgvars_list Object to manage prognostic variables with atmospheric dynamical core
!! @param trcvars_list Object to manage auxiliary variables 
!! @param forcing_list Object to manage forcing terms
!! @param is_update Flag to speicfy whether the tendencies are updated in this call
!!
!OCL SERIAL  
  subroutine AtmosPhyRd_update( this, model_mesh, &
    prgvars_list, trcvars_list,                   &
    auxvars_list, forcing_list, is_update         )  
    
    implicit none
    class(AtmosPhyRd), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------
    return
  end subroutine AtmosPhyRd_update

!> Finalize a component of radiation in atmospheric model
!!
!OCL SERIAL  
  subroutine AtmosPhyRd_finalize( this )
    implicit none
    class(AtmosPhyRd), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    select case ( this%RD_TYPEID )
    case( RD_TYPEID_SIMPLE )
      call this%simple_rad%Final()
    end select

    call this%vars%Final()
    return
  end subroutine AtmosPhyRd_finalize

!- private ------------------------------------------------
!OCL SERIAL
  subroutine AtmosPhyRd_calc_tendency_core( this, &
    RHOH,                                           & ! (out)
    SFLX_SW_up, SFLX_SW_dn, SFLX_LW_up, SFLX_LW_dn, & ! (out)
    SOLINS, COS_SZA,                                & ! (in)
    DDENS, PRES, QV, SFC_TEMP, SFC_ALB,             & ! (in)
    DENS_hyd, Rtot, CVtot,                          & ! (in)
    lcmesh, elem3D, lcmesh2D, elem2D,               & ! (in)
    elem3D_operation )                                ! (in)
    use scale_atmos_phy_rd_common, only: &
      I_up, I_dn, I_LW, I_SW
    use scale_element_operation_base, only: ElementOperationBase3D
    use scale_atm_phy_rd_dgm_common, only: ATM_PHY_RD_DGM_calc_heating
    implicit none
    class(AtmosPhyRd), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: RHOH(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: SFLX_SW_up(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(out) :: SFLX_SW_dn(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(out) :: SFLX_LW_up(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(out) :: SFLX_LW_dn(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: SOLINS(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: COS_SZA(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: QV(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: SFC_TEMP(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: SFC_ALB(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA)
    class(ElementOperationBase3D), intent(in) :: elem3D_operation

    real(RP) :: flux_rad(elem3D%Np,lcmesh%Ne,2,2,2)
    real(RP) :: flux_rad_top(elem3D%Nnode_h1D**2,lcmesh%Ne2D,2,2,2)
    real(RP) :: sflux_rad_up(elem3D%Nnode_h1D**2,lcmesh%Ne2D,2,2)
    real(RP) :: sflux_rad_dn(elem3D%Nnode_h1D**2,lcmesh%Ne2D,2,2)
    real(RP) :: TEMP_(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: DENS_(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: PRES_(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP) :: QV_(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    
    ! real(RP) :: flux_up(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    ! real(RP) :: flux_dn(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    ! real(RP) :: flux_net(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    ! real(RP) :: flux_net_sfc(elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    ! real(RP) :: flux_net_toa(elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    ! real(RP) :: flux_net_tom(elem3D%Nnode_h1D**2,lcmesh%Ne2D)

    integer :: ke, ke_h, ke_z
    integer :: p, ph, pz

    real(RP) :: dens_tmp
    !--------------------------------------------------

    !$omp parallel do private(ke,p, dens_tmp) collapse(2)
    do ke_z=1, lcmesh%NeZ
    do ke_h=1, lcmesh%Ne2D
      ke = ke_h + (ke_z-1)*lcmesh%Ne2D
      do pz=1, elem3D%Nnode_v
      do ph=1, elem3D%Nnode_h1D**2
        p = ph + (pz-1)*elem3D%Nnode_h1D**2
        dens_tmp = DENS_hyd(p,ke) + DDENS(p,ke)
        
        PRES_(pz,ke_z,ph,ke_h) = PRES(p,ke)
        TEMP_(pz,ke_z,ph,ke_h) = PRES(p,ke) / ( Rtot(p,ke) * dens_tmp )
        DENS_(pz,ke_z,ph,ke_h) = dens_tmp
        QV_(pz,ke_z,ph,ke_h)   = QV(p,ke)
      end do
      end do
    end do
    end do

    select case( this%RD_TYPEID )
    case ( RD_TYPEID_SIMPLE )
      call this%simple_rad%calculate_rad_flux( flux_rad(:,:,:,:,1),            & ! (out)
        flux_rad_top(:,:,:,:,1), sflux_rad_up(:,:,:,1), sflux_rad_dn(:,:,:,1), & ! (out)
        SOLINS, PRES_, TEMP_, DENS_, QV_, SFC_TEMP, SFC_ALB,                   & ! (in)
        lcmesh, elem3D, lcmesh2D, elem2D )                                       ! (in)
      
      !$omp parallel do
      do ke=lcmesh2D%NeS, lcmesh2D%NeE
        SFLX_SW_dn(:,ke) = sflux_rad_dn(:,ke,I_SW,1)
        SFLX_LW_dn(:,ke) = sflux_rad_dn(:,ke,I_LW,1)
        SFLX_SW_up(:,ke) = sflux_rad_up(:,ke,I_SW,1)
        SFLX_LW_up(:,ke) = sflux_rad_up(:,ke,I_LW,1)
      end do
    end select

    call ATM_PHY_RD_DGM_calc_heating( RHOH, &
      flux_rad(:,:,:,:,1), DDENS, DENS_hyd, CVtot, &
      lcmesh, elem3D, elem2D, &
      elem3D_operation )

    return
  end subroutine AtmosPhyRd_calc_tendency_core
end module mod_atmos_phy_rd
