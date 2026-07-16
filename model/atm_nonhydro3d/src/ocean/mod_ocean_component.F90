!-------------------------------------------------------------------------------
!> module Oceanic component
!!
!! @par Description
!!          Oceanic component module
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_ocean_component
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc

  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList
  use scale_model_component, only: ModelComponent
  use scale_model_mesh_manager, only: ModelMesh3D

  use mod_ocean_mesh, only: OceanMesh
  use mod_ocean_mesh_rm, only: OceanMeshRM
  use mod_ocean_mesh_gm, only: OceanMeshGM

  use mod_cpl_component, only: CouplerComponent

  use mod_ocean_vars, only: &
    OceanVars,                                       &
    PRGVAR_THERM_ID => PRGVAR_THERM_ID,              &
    SFC_TEMP_ID => AUXVAR2D_SFC_TEMP_ID,             &
    ATM_SFC_DENS_ID => ATMVAR2D_SFC_DENS_ID,         &
    ATM_SFC_PRES_ID => ATMVAR2D_SFC_PRES_ID,         &
    ATM_TEMP_ID => ATMVAR2D_ATM_TEMP_ID,             &
    ATM_DENS_ID => ATMVAR2D_ATM_DENS_ID,             &
    ATM_PRES_ID => ATMVAR2D_ATM_PRES_ID,             &
    ATM_W_ID    => ATMVAR2D_ATM_W_ID,                &
    ATM_U_ID    => ATMVAR2D_ATM_U_ID,                &
    ATM_V_ID    => ATMVAR2D_ATM_V_ID,                &
    ATM_QV_ID   => ATMVAR2D_ATM_QV_ID,               &
    SFLX_RD_SW_DIR_ID => ATMVAR2D_SFLX_RD_SW_DIR_ID, &
    SFLX_RD_LW_DIF_ID => ATMVAR2D_SFLX_RD_LW_DIF_ID, &
    ZLEV_A_ID => ATMVAR2D_LOWEST_LAYER_ZLEV_ID


  use mod_ocean_dyn, only: OceanDyn

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage oceanic component
  type, extends(ModelComponent), public :: OceanComponent
    type(OceanVars) :: vars             !< Object to mange variables with oceanic component

    character(len=H_SHORT) :: mesh_type !< Type name of mesh oceanic component
    class(OceanMesh), pointer :: mesh   !< Pointer of mesh oceanic component
    type(OceanMeshRM) :: mesh_rm        !< Object to manage mesh for the case of regional mode
    type(OceanMeshGM) :: mesh_gm        !< Object to manage mesh for the case of global mode

    type(OceanDyn) :: dyn_proc          !< Object to manage dynamical process

    type(CouplerComponent), pointer :: coupler_ptr !< Pointer of coupler component
  contains
    procedure, public :: setup => Ocean_setup 
    procedure, public :: setup_vars => Ocean_setup_vars   
    procedure, public :: set_coupler => Ocean_set_coupler 
    procedure, public :: calc_tendency => Ocean_calc_tendency
    procedure, public :: update => Ocean_update
    procedure, public :: set_surface => Ocean_set_surface
    procedure, public :: get_surface => Ocean_get_surface
    procedure, public :: finalize => Ocean_finalize
  end type OceanComponent

  !-----------------------------------------------------------------------------
  !
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
  !-----------------------------------------------------------------------------
contains

  !> Setup an object to mange oceanic component
!OCL SERIAL
  subroutine Ocean_setup( this )
    use scale_const, only: &
      UNDEF8 => CONST_UNDEF8
    
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_setup  
    use scale_time_manager, only: &
      TIME_manager_Regist_component

    implicit none
    
    class(OceanComponent), intent(inout), target :: this

    logical :: ACTIVATE_FLAG = .false. !< Flag whether oceanic component is activated

    real(DP) :: TIME_DT                             = UNDEF8  !< Timestep value of oceanic component
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'   !< Timestep unit of oceanic component
    real(DP) :: TIME_DT_RESTART                     = UNDEF8  !< Timestep value when outputting restart file for oceanic component
    character(len=H_SHORT) :: TIME_DT_RESTART_UNIT  = 'SEC'   !< Timestep unit when outputting restart file for oceanic component

    logical :: OCEAN_DYN_DO    = .true.  !< Flag whether dynamics process is considered
    character(len=H_SHORT) :: OCEAN_MESH_TYPE = 'REGIONAL'  !< Name of mesh type for oceanic component ('REGIONAL' or 'GLOBAL')


    namelist / PARAM_OCEAN / &
      ACTIVATE_FLAG,         &
      TIME_DT,               &
      TIME_DT_UNIT,          &
      TIME_DT_RESTART,       &
      TIME_DT_RESTART_UNIT,  &
      OCEAN_MESH_TYPE,       &
      OCEAN_DYN_DO
    
    integer :: ierr
    !--------------------------------------------------

    LOG_INFO('OceanComponent_setup',*) 'Oceanic model components '
    
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("Ocean_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("Ocean_setup",*) 'Not appropriate names in namelist PARAM_OCEAN. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN)

    !************************************************
    call this%ModelComponent_Init('OCEAN', ACTIVATE_FLAG )

    if ( .not. ACTIVATE_FLAG ) return
    
    call PROF_rapstart( 'Ocean_setup', 1)

    !- Setup time manager

    call this%time_manager%Init( this%GetComponentName(), &
      TIME_DT, TIME_DT_UNIT,                              &
      TIME_DT_RESTART, TIME_DT_RESTART_UNIT               ) 
    
    call TIME_manager_Regist_component( this%time_manager )
    
    !- Setup mesh & file I/O for oceanic component
    
    this%mesh_type = OCEAN_MESH_TYPE
    select case( this%mesh_type )
    case('REGIONAL')
      call this%mesh_rm%Init()
      call FILE_HISTORY_meshfield_setup( mesh3d_=this%mesh_rm%mesh, & ! (in)
        dim_name_postfix_='_O',                                     & ! (in)
        registered_comp_id=this%vars%hist_comp_id                   ) ! (out)
      this%mesh => this%mesh_rm
    case('GLOBAL')
      call this%mesh_gm%Init()
      call FILE_HISTORY_meshfield_setup( meshCubedSphere3D_=this%mesh_gm%mesh, & ! (in)
        dim_name_postfix_='_O',                                                & ! (in)      
        registered_comp_id=this%vars%hist_comp_id )                              ! (out)
      this%mesh => this%mesh_gm
    case default
      LOG_ERROR("Ocean_setup",*) 'Unsupported type of mesh is specified. Check!', this%mesh_type
      call PRC_abort    
    end select
    
    !- setup common tools for oceanic model

    !- Setup each processes in ocean model ------------------------------------

    !- Setup the module for ocean / dynamics 
    call this%dyn_proc%ModelComponentProc_Init( 'OceanDyn', OCEAN_DYN_DO )
    call this%dyn_proc%setup( this%mesh, this%time_manager )
    
    !- Setup

    LOG_NEWLINE
    LOG_INFO('OceanComponent_setup',*) 'Finish setup of each oceanic component.'

    call PROF_rapend( 'Ocean_setup', 1)

    return
  end subroutine Ocean_setup

  !> Setup variables with the oceanic component
!OCL SERIAL
  subroutine Ocean_setup_vars( this )
    implicit none
    class(OceanComponent), intent(inout) :: this
    !----------------------------------------------------------

    call PROF_rapstart( 'Ocean_setup_vars', 1)
    call this%vars%Init( this%mesh )
    call PROF_rapend( 'Ocean_setup_vars', 1)
    return
  end subroutine Ocean_setup_vars

!> Set coupler component to the oceanic component
!OCL SERIAL
  subroutine Ocean_set_coupler( this, coupler )
    implicit none
    class(OceanComponent), intent(inout) :: this
    class(CouplerComponent), target, intent(inout) :: coupler
    !----------------------------------------------------------
    this%coupler_ptr => coupler
    return
  end subroutine Ocean_set_coupler  

!> Calculate tendencies with the oceanic component
!OCL SERIAL
  subroutine Ocean_calc_tendency( this, force )
    use mod_ocean_vars, only: &
      ALB_VIS_DIR_ID => AUXVAR2D_SFC_ALB_VIS_dir_ID, &
      SFLX_MW_ID => OCN_SFLX_MW_ID, &
      SFLX_MU_ID => OCN_SFLX_MU_ID, &
      SFLX_MV_ID => OCN_SFLX_MV_ID, &
      SFLX_SH_ID => OCN_SFLX_SH_ID, &
      SFLX_LH_ID => OCN_SFLX_LH_ID, &
      SFLX_RD_SW_DIR_ID => ATMVAR2D_SFLX_RD_SW_DIR_ID, &
      SFLX_RD_LW_DIF_ID => ATMVAR2D_SFLX_RD_LW_DIF_ID, &
      RHOH_ID => PHYTEND_RHOH_ID
    implicit none
    class(OceanComponent), intent(inout) :: this
    logical, intent(in) :: force

    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh2D), pointer :: lmesh2D
    class(LocalMesh3D), pointer :: lmesh3D
    integer :: idom

    real(RP), allocatable :: SFLX_GH(:,:)
    !------------------------------------------------------------------
    
    !- Get surface data from coupler 
    call this%get_surface()

    call PROF_rapstart( 'OCN_tendency', 1)

    mesh3D => this%mesh%ptr_mesh

    do idom=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(idom)
      lmesh2D => lmesh3D%lcmesh2D

      allocate( SFLX_GH(lmesh2D%refElem2D%Np,lmesh2D%Ne) )

      !- Calculate surface fluxes

      call calculate_surface_flux( &
        this%vars%OCN_SFLX(SFLX_MW_ID)%local(idom)%val, &
        this%vars%OCN_SFLX(SFLX_MU_ID)%local(idom)%val, &
        this%vars%OCN_SFLX(SFLX_MV_ID)%local(idom)%val, &
        this%vars%OCN_SFLX(SFLX_SH_ID)%local(idom)%val, &
        this%vars%OCN_SFLX(SFLX_LH_ID)%local(idom)%val, &
        !-
        this%vars%AUX_VARS2D(SFC_TEMP_ID)%local(idom)%val, &
        this%vars%ATM_VARS2D(ATM_SFC_DENS_ID)%local(idom)%val, &
        this%vars%ATM_VARS2D(ATM_SFC_PRES_ID)%local(idom)%val, &
        this%vars%ATM_VARS2D(ATM_TEMP_ID)%local(idom)%val, &
        this%vars%ATM_VARS2D(ATM_DENS_ID)%local(idom)%val, &
        this%vars%ATM_VARS2D(ATM_PRES_ID)%local(idom)%val, &
        this%vars%ATM_VARS2D(ATM_W_ID)%local(idom)%val, &
        this%vars%ATM_VARS2D(ATM_U_ID)%local(idom)%val, &
        this%vars%ATM_VARS2D(ATM_V_ID)%local(idom)%val, &
        this%vars%ATM_VARS2D(ATM_QV_ID)%local(idom)%val, &
        this%vars%ATM_VARS2D(ZLEV_A_ID)%local(idom)%val, &
        lmesh2D, lmesh2D%refElem2D )

      !- Calculate ground heat flux
      
      call calculate_ground_heat_flux( SFLX_GH,                  & ! (out)
        this%vars%AUX_VARS2D(SFC_TEMP_ID)%local(idom)%val,       & ! (in)
        this%vars%ATM_VARS2D(SFLX_RD_SW_DIR_ID)%local(idom)%val, & ! (in)
        this%vars%ATM_VARS2D(SFLX_RD_LW_DIF_ID)%local(idom)%val, & ! (in)
        this%vars%AUX_VARS2D(ALB_VIS_DIR_ID)%local(idom)%val,    & ! (in)
        lmesh2D, lmesh2D%refElem2D )                               ! (in)

      !- Calculate tendencies with physics

      call calculate_phys_tendency( this%vars%PHY_TEND(RHOH_ID)%local(idom)%val, & ! (out)
        SFLX_GH, lmesh3D%zlev,                                                   & ! (in)
        lmesh3D, lmesh3D%refElem3D, lmesh2D, lmesh2D%refElem2D )                   ! (in)
      
      ! write(*,*) "------"
      ! write(*,*) "SFC_TEMP = ", this%vars%AUX_VARS2D(SFC_TEMP_ID)%local(idom)%val(1,1)
      ! write(*,*) "SFC_DENS = ", this%vars%ATM_VARS2D(ATM_SFC_DENS_ID)%local(idom)%val(1,1)
      ! write(*,*) "SFC_PRES = ", this%vars%ATM_VARS2D(ATM_SFC_PRES_ID)%local(idom)%val(1,1)
      ! write(*,*) "ATM_TEMP = ", this%vars%ATM_VARS2D(ATM_TEMP_ID)%local(idom)%val(1,1)
      ! write(*,*) "ATM_PRES = ", this%vars%ATM_VARS2D(ATM_PRES_ID)%local(idom)%val(1,1)
      ! write(*,*) "ATM_QV= ", this%vars%ATM_VARS2D(ATM_QV_ID)%local(idom)%val(1,1)
      ! write(*,*) "RD_SW_DIR = ", this%vars%ATM_VARS2D(SFLX_RD_SW_DIR_ID)%local(idom)%val(1,1)
      ! write(*,*) "RD_LW_DIF = ", this%vars%ATM_VARS2D(SFLX_RD_LW_DIF_ID)%local(idom)%val(1,1)
      ! write(*,*) "ZLEV_A = ", this%vars%ATM_VARS2D(ZLEV_A_ID)%local(idom)%val(1,1)
      ! write(*,*) "SLFX_MW=", this%vars%OCN_SFLX(SFLX_MW_ID)%local(idom)%val(1,1)
      ! write(*,*) "SLFX_MU=", this%vars%OCN_SFLX(SFLX_MU_ID)%local(idom)%val(1,1)
      ! write(*,*) "SLFX_MV=", this%vars%OCN_SFLX(SFLX_MV_ID)%local(idom)%val(1,1)
      ! write(*,*) "SH=", this%vars%OCN_SFLX(SFLX_SH_ID)%local(idom)%val(1,1)
      ! write(*,*) "LH=", this%vars%OCN_SFLX(SFLX_LH_ID)%local(idom)%val(1,1)
      ! write(*,*) 'SFLX_GH = ', SFLX_GH(1,1)

      deallocate( SFLX_GH )
    end do

    call PROF_rapend( 'OCN_tendency', 1)

    !-
    call this%set_surface( countup=.true. )

    return  
  end subroutine Ocean_calc_tendency

!> Update variables with the oceanic component
!OCL SERIAL
  subroutine Ocean_update( this )
    implicit none
    class(OceanComponent), intent(inout) :: this
    
    integer :: tm_process_id
    logical :: is_update
    integer :: inner_itr

    integer :: idom
    class(LocalMesh2D), pointer :: lmesh2D
    class(LocalMesh3D), pointer :: lmesh
    !--------------------------------------------------
    call PROF_rapstart( 'OCN_update', 1)

    if ( this%dyn_proc%IsActivated() ) then
      call PROF_rapstart('OCN_Dynamics', 1)
      tm_process_id = this%dyn_proc%tm_process_id
      is_update = this%time_manager%Do_process( tm_process_id )

      LOG_PROGRESS(*) 'ocean / dynamics'
      do inner_itr=1, this%time_manager%Get_process_inner_itr_num( tm_process_id )
        call this%dyn_proc%update( &
          this%mesh, this%vars%PROGVARS_manager, this%vars%QTRCVARS_manager, &
          this%vars%AUXVARS_manager, this%vars%PHYTENDS_manager, is_update   )
      end do
      call PROF_rapend('OCN_Dynamics', 1)
    end if

    !- Set surface temperature
    do idom=1, this%mesh%ptr_mesh%LOCAL_MESH_NUM
      lmesh => this%mesh%ptr_mesh%lcmesh_list(idom)
      lmesh2D => lmesh%lcmesh2D
      call set_sfctemp_lc( this%vars%AUX_VARS2D(SFC_TEMP_ID)%local(idom)%val, & ! (out)
        this%vars%PROG_VARS(PRGVAR_THERM_ID)%local(idom)%val,                 & ! (in)
        lmesh, lmesh%refElem3D, lmesh2D, lmesh2D%refElem2D )                    ! (in)
      ! write(*,*) "OCN_update: SFC_TEMP = ", this%vars%AUX_VARS2D(SFC_TEMP_ID)%local(idom)%val(1,1)
    end do

    call PROF_rapend('OCN_update', 1)
    return  
  end subroutine Ocean_update
!OCL SERIAL
  subroutine set_sfctemp_lc( SFC_TEMP, &
    THERM, lmesh, elem, lmesh2D, elem2D )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: SFC_TEMP(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: THERM(elem%Np,lmesh%NeA)

    integer :: ke, ke2D
    integer :: hSlice(elem2D%Np)
    !------------------------------------------------------------------

    hSlice(:) = elem%Hslice(:,elem%Nnode_v)
    !$omp parallel do private(ke)
    do ke2D=lmesh2D%NeS, lmesh2D%NeE
      ke = ke2D + (lmesh%NeZ-1)*lmesh%Ne2D
      SFC_TEMP(:,ke2D) = THERM(hSlice(:),ke)
    end do
    return
  end subroutine set_sfctemp_lc

  !> Set ocean quantities to coupler component
!OCL SERIAL
  subroutine Ocean_set_surface( this, countup )
    implicit none
    class(OceanComponent), intent(inout) :: this
    logical, intent(in) :: countup
    !--------------------------------------------------

    call PROF_rapstart( 'OCN_sfc_exch', 1)

    call this%coupler_ptr%vars%PutOCN( this%vars, countup )

    call PROF_rapend( 'OCN_sfc_exch', 1)
    return
  end subroutine Ocean_set_surface

  !> Get atmospheric quantities from coupler component
!OCL SERIAL
  subroutine Ocean_get_surface( this )
    implicit none
    class(OceanComponent), intent(inout) :: this
    !--------------------------------------------------

    call PROF_rapstart( 'OCN_sfc_exch', 1)

    call this%coupler_ptr%vars%Get_ATM_OCN( &
      this%vars%ATM_VARS2D(ATM_SFC_DENS_ID),   &
      this%vars%ATM_VARS2D(ATM_SFC_PRES_ID),   &
      this%vars%ATM_VARS2D(ATM_TEMP_ID),       &
      this%vars%ATM_VARS2D(ATM_PRES_ID),       &
      this%vars%ATM_VARS2D(ATM_W_ID),          &
      this%vars%ATM_VARS2D(ATM_U_ID),          &
      this%vars%ATM_VARS2D(ATM_V_ID),          &
      this%vars%ATM_VARS2D(ATM_QV_ID),         &
      this%vars%ATM_VARS2D(SFLX_RD_SW_DIR_ID), &
      this%vars%ATM_VARS2D(SFLX_RD_LW_DIF_ID), &
      this%vars%ATM_VARS2D(ZLEV_A_ID) )
    
    call PROF_rapend( 'OCN_sfc_exch', 1)
    return
  end subroutine Ocean_get_surface

!> Finalize an object to manage the ocean component
!OCL SERIAL
  subroutine Ocean_finalize( this )
    implicit none
    class(OceanComponent), intent(inout) :: this
    !--------------------------------------------------

    LOG_INFO('OceanComponent_finalize',*)

    if ( .not. this%IsActivated() ) return

    ! call this%dyn_proc%finalize()
    call this%vars%Final()

    select case( this%mesh_type )
    case('REGIONAL')
      call this%mesh_rm%Final()
    case('GLOBAL')
      call this%mesh_gm%Final()
    end select  
    this%mesh => null()

    call this%time_manager%Final()

    return  
  end subroutine Ocean_finalize

!- Private subroutines -------------------------------------------------------------

  !> Calculate momentum and heat flux at the surface (tentative)
!OCL SERIAL
  subroutine calculate_surface_flux( SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_LH, &
    SFC_TEMP, SFC_DENS, SFC_PRES,                                       &
    ATM_TEMP, ATM_DENS, ATM_PRES, ATM_W, ATM_U, ATM_V, ATM_QV, zlev_a, &
    lmesh, elem )
    use scale_const, only: &
      RPlanet => CONST_RADIUS    
    use scale_atm_phy_sf_bulk_simple, only: &
       ATMOS_PHY_SF_simple_flux
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: SFLX_MW(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: SFLX_MU(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: SFLX_MV(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: SFLX_SH(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: SFLX_LH(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: SFC_TEMP(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: SFC_DENS(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: SFC_PRES(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ATM_TEMP(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ATM_DENS(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ATM_PRES(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ATM_W(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ATM_U(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ATM_V(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ATM_QV(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: zlev_a(elem%Np,lmesh%NeA)
    
    real(RP) :: SFLX_QV(elem%Np,lmesh%NeA)

    ! dummy
    real(RP) :: U10(elem%Np,lmesh%NeA)
    real(RP) :: V10(elem%Np,lmesh%NeA)

    integer :: ke
    integer :: ke2D, ij

    real(RP) :: DZ1     (elem%Np,lmesh%NeA)
    real(RP) :: Z1      (elem%Np,lmesh%NeA)
    !--------------------------------------------------

    !$omp parallel do collapse(2)
    do ke2D=lmesh%NeS, lmesh%NeE
    do ij=1, elem%Np

      Z1(ij,ke2D) = zlev_a(ij,ke2D)
      DZ1(ij,ke2D) = Z1(ij,ke2D) - zlev_a(ij,ke2D)
      Z1(ij,ke2D) = Z1(ij,ke2D) + RPlanet
    end do
    end do

    call ATMOS_PHY_SF_simple_flux( &
      elem%Np, 1, elem%Np, lmesh%NeA, 1, lmesh%Ne,                      & ! (in)
      ATM_W(:,:), ATM_U(:,:), ATM_V(:,:), ATM_TEMP(:,:), ATM_PRES(:,:), & ! (in) Note: ATM_PRES is not used
      ATM_QV(:,:),                                                      & ! (in)
      SFC_DENS(:,:), SFC_TEMP(:,:), SFC_PRES(:,:),                      & ! (in)
      DZ1(:,:),                                                         & ! (in)
      SFLX_MW(:,:), SFLX_MU(:,:), SFLX_MV(:,:),                         & ! (out)
      SFLX_SH(:,:), SFLX_LH(:,:), SFLX_QV(:,:),                         & ! (out)
      U10(:,:), V10(:,:)                                                ) ! (out)

    return
  end subroutine calculate_surface_flux

  !> Calculate ground heat flux (tentative)
  !!
!OCL SERIAL
  subroutine calculate_ground_heat_flux( sflx_GH, &
    SFC_TEMP, & 
    SFLX_RD_SW_dn_dir, SFLX_RD_LW_dn_dif, &
    SFC_ALB_dir_vis, &
    lmesh, elem )
    use scale_const, only: &
       STB   => CONST_STB
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: sflx_GH(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: SFC_TEMP(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: SFLX_RD_SW_dn_dir(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: SFLX_RD_LW_dn_dif(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: SFC_ALB_dir_vis(elem%Np,lmesh%NeA)

    integer :: ke
    real(RP) :: emis(elem%Np)
    real(RP) :: SWU(elem%Np), SWD(elem%Np), LWU(elem%Np), LWD(elem%Np)
    !--------------------------------------------------

    !$omp parallel do private(SWU,SWD,LWU,LWD,emis)
    do ke=lmesh%NeS, lmesh%NeE
      emis(:) = STB * SFC_TEMP(:,ke)**4

      LWD(:) = SFLX_RD_LW_dn_dif(:,ke)
      LWU(:) = emis(:)
      SWD(:) = SFLX_RD_SW_dn_dir(:,ke)
      SWU(:) = SFC_ALB_dir_vis(:,ke) * SWD(:)

      sflx_GH(:,ke) = SWD(:) - SWU(:) + LWD(:) - LWU(:)
    end do
    return
  end subroutine calculate_ground_heat_flux

!OCL SERIAL
  subroutine calculate_phys_tendency( RHOH, &
    sflx_GH, zlev, lmesh, elem, lmesh2D, elem2D )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: RHOH(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: sflx_GH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: zlev(elem%Np,lmesh%Ne2D,lmesh%NeZ)

    integer :: ke, ke_z, ke2D
    integer :: ph, pz, p
    real(RP) :: r_ocn_depth(elem2D%Np,lmesh2D%Ne)
    integer :: hSlice_t(elem2D%Np), hSlice_b(elem2D%Np)
    !---------------------------------------

    hSlice_t(:) = elem%Hslice(:,elem%Nnode_v)
    hSlice_b(:) = elem%Hslice(:,1)

    !$omp parallel private(ke)
    !$omp do
    do ke2D=1, lmesh%Ne2D
      r_ocn_depth(:,ke2D) = 1.0_RP / ( zlev(hSlice_t(:),ke2D,lmesh%NeZ) - zlev(hSlice_b(:),ke2D,1) )
    end do
    !$omp do collapse(2)
    do ke_z=1, lmesh%NeZ
    do ke2D=1, lmesh%Ne2D
      ke = ke2D + (ke_z-1)*lmesh%Ne2D
      do pz=1, elem%Nnode_v
      do ph=1, elem%Nnode_h1D**2
        p = ph + (pz-1)*elem%Nnode_h1D**2
        RHOH(p,ke) = sflx_GH(ph,ke2D) * r_ocn_depth(ph,ke2D)
      end do
      end do
    end do
    end do
    !$omp end parallel
    return
  end subroutine calculate_phys_tendency
end module mod_ocean_component