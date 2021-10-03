!-------------------------------------------------------------------------------
!> module INITIAL
!!
!! @par Description
!!          subroutines for preparing initial data
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_mkinit
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc

  use scale_const, only: &
    PI => CONST_PI,          &
    GRAV => CONST_GRAV,      &
    RPlanet => CONST_RADIUS, &
    OHM => CONST_OHM,        &
    Rdry => CONST_Rdry,      &
    CPdry => CONST_CPdry,    &
    CVdry => CONST_CVdry,    &
    PRES00 => CONST_PRE00,   &
    Pstd   => CONST_Pstd  


  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_atm_dyn_dgm_hydrostatic, only: &
    hydrostaic_build_rho_XYZ

  use mod_atmos_component, only: &
    AtmosComponent
      
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKINIT_setup
  public :: MKINIT

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  integer, public :: MKINIT_TYPE   = -1
  integer, public :: I_IGNORE      = 0

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKINIT_setup
    implicit none

    character(len=H_SHORT) :: initname = 'NONE'

    namelist / PARAM_MKINIT / &
      initname

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("MKINIT_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("MKINIT_setup",*) 'Not appropriate names in namelist PARAM_MKINIT. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT)

    select case(trim(initname))
    case('NONE')
      MKINIT_TYPE = I_IGNORE
    end select
    
    return
  end subroutine MKINIT_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine MKINIT( output,                                      & ! (out)
    model_mesh,                                                   & ! (in)
    atm_prgvars_manager, atm_auxvars_manager, atm_trcvars_manager )
  
    use scale_model_var_manager, only: ModelVarManager
    use mod_atmos_mesh, only: AtmosMesh
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars

    implicit none

    logical, intent(out) :: output
    class(AtmosMesh), target, intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: atm_prgvars_manager
    class(ModelVarManager), intent(inout) :: atm_auxvars_manager
    class(ModelVarManager), intent(inout) :: atm_trcvars_manager

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd

    integer :: n
    integer :: ke
    class(LocalMesh3D), pointer :: lcmesh3D
    class(MeshBase3D), pointer :: mesh
    !---------------------------------------------------------------------------

    mesh => model_mesh%ptr_mesh

    if ( MKINIT_TYPE == I_IGNORE ) then
      LOG_NEWLINE
      LOG_PROGRESS(*) 'skip  making initial data'
      output = .false.
    else
      LOG_NEWLINE
      LOG_PROGRESS(*) 'start making initial data'

      ! call PROF_rapstart('_MkInit_main',3)   
      
      do n=1, mesh%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshPrgVars( n, mesh, atm_prgvars_manager, atm_auxvars_manager, &
           DDENS, MOMX, MOMY, MOMZ, DRHOT,                                                     &
           DENS_hyd, PRES_hyd, lcmesh3D                                                        )

        !$omp parallel do
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          MOMX %val(:,ke) = 0.0_RP
          MOMY %val(:,ke) = 0.0_RP
          MOMZ %val(:,ke) = 0.0_RP
          DDENS%val(:,ke) = 0.0_RP
          DRHOT%val(:,ke) = 0.0_RP
        end do
      end do
      ! call PROF_rapend  ('_MkInit_main',3)
        
      ! LOG_PROGRESS(*) 'end   making initial data'
      output = .true.
    end if

    return
  end subroutine MKINIT

  !-- private---------------------------------------------------------------------
end module mod_mkinit