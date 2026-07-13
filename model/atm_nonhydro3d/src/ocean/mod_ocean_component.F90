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

  use mod_ocean_vars, only: &
    OceanVars

!   use mod_ocean_dyn, only: OceanDyn

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

        ! type(OceanDyn) :: dyn_proc          !< Object to manage dynamical process
  contains
    procedure, public :: setup => Ocean_setup 
    procedure, public :: setup_vars => Ocean_setup_vars    
    procedure, public :: calc_tendency => Ocean_calc_tendency
    procedure, public :: update => Ocean_update
    procedure, public :: set_surface => Ocean_set_surface
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

    logical :: OCEAN_DYN_DO    = .false.  !< Flag whether dynamics process is considered
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

    !- Setup the module for Oceanic / dynamics 
    ! call this%dyn_proc%ModelComponentProc_Init( 'OceanDyn', OCEAN_DYN_DO )
    ! call this%dyn_proc%setup( this%mesh, this%time_manager )
    
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

!> Calculate tendencies with the oceanic component
!OCL SERIAL
  subroutine Ocean_calc_tendency( this, force )
    ! use mod_ocean_vars, only: &
    !   OceanVars_GetLocalMeshPhyTends

    implicit none
    class(OceanComponent), intent(inout) :: this
    logical, intent(in) :: force
    !------------------------------------------------------------------
    
    call PROF_rapstart( 'Ocean_calc_tendency', 1)
    call PROF_rapend( 'Ocean_calc_tendency', 1)
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
    !--------------------------------------------------
    call PROF_rapstart( 'OCN_update', 1)
    call PROF_rapend('OCN_update', 1)
    return  
  end subroutine Ocean_update

!OCL SERIAL
  subroutine Ocean_set_surface( this )
    implicit none
    class(OceanComponent), intent(inout) :: this
    !--------------------------------------------------
    call PROF_rapstart( 'OCN_sfc_exch', 1)
    call PROF_rapend( 'OCN_sfc_exch', 1)
    return
  end subroutine Ocean_set_surface

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

end module mod_ocean_component