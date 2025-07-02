!-------------------------------------------------------------------------------
!> module SCALE-DG prep
!!
!! @par Description
!!          This program is driver of preprocess tools
!!          1) boundary data (e.g. topography, land use index)
!!          2) initial data for ideal/real test cases
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_dg_prep
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc

  use mod_atmos_component, only: &
    AtmosComponent
  
  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_write
  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate

  use mod_user, only: &
    USER
  use mod_mktopo, only: &
    MKTOPO, MKTOPO_write
  use mod_mkinit, only: &
    MKINIT
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
#include "scale-dg.h" 

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: dg_prep
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
  character(len=H_MID), private, parameter :: MODELNAME = "SCALE-DG ver. "//VERSION

  type(AtmosComponent) :: atmos
  type(User) :: user_

contains
  subroutine dg_prep(                     &
    comm_world, cnf_fname, path, add_path )

    use scale_time_manager, only: &
      TIME_manager_checkstate, TIME_manager_advance,          &
      TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP, TIME_NSTEP, &
      TIME_DOresume, TIME_DOend

    implicit none

    integer,          intent(in) :: comm_world
    character(len=*), intent(in) :: cnf_fname
    character(len=*), intent(in) :: path
    logical,          intent(in) :: add_path


    integer :: myrank
    logical :: ismaster

    logical :: output
    logical :: output_topo
    !---------------------------------------------------------------------------

    !########## Initial setup ##########

#ifdef SCALE_DEVELOP
    ! setup standard I/O
    if ( add_path .and. path /= "" ) then
      call IO_setup( MODELNAME, trim(path)//cnf_fname, prefix=path )
    else
#endif
      call IO_setup( MODELNAME, trim(path)//cnf_fname )
#ifdef SCALE_DEVELOP
    end if
#endif

    ! setup MPI
    call PRC_LOCAL_setup( comm_world, & ! [IN]
                          myrank,     & ! [OUT]
                          ismaster    ) ! [OUT]

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
 
    call initialize()

    !###########################################################################
 
    !########## main ##########
    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_prep', 0)

    !- Execute preprocess

    !- Execute mktopo
    call PROF_rapstart('MkTopo',1)
    call MKTOPO( output_topo,           &
      atmos%mesh, atmos%mesh%topography )
    call PROF_rapend  ('MkTopo',1)

    !- Re-setup
    call atmos%mesh%Setup_vcoordinate()

    !- Execute mkinit
    call PROF_rapstart('MkInit',1)
    call MKINIT( output, &
      atmos%mesh,                  &
      atmos%vars%PROGVARS_manager, &
      atmos%vars%AUXVARS_manager,  &
      atmos%vars%QTRCVARS_manager  )
    
!    call USER_mkinit( atmos )
    call user_%mkinit( atmos )
    if ( atmos%dyn_proc%dyncore_driver%ENTOT_CONSERVE_SCHEME_FLAG ) then
      call set_total_energy( atmos%vars%PROGVARS_manager, &
        atmos%vars%AUXVARS_manager, atmos%mesh )      
    end if
 
    call PROF_rapend  ('MkInit',1)
    call PROF_rapend('Main_prep', 0)

    !- Output

    if ( output_topo ) call MKTOPO_write( atmos%mesh, atmos%mesh%topography )

    if ( output ) then
      call PROF_rapstart('MkInit_restart',1)
      call restart_write()
      call PROF_rapend  ('MkInit_restart',1)
    end if

    !########## Finalize ##########
    call user_%mkfinal()
    call finalize()    

    return
  end subroutine dg_prep

  !----------------------------
  
  subroutine initialize()

    use scale_const, only: CONST_setup
    use scale_calendar, only: CALENDAR_setup
    use scale_random, only: RANDOM_setup
    use scale_atmos_hydrometeor, only: ATMOS_HYDROMETEOR_setup    

    use scale_file_restart_meshfield, only: &
      FILE_restart_meshfield_setup
    use scale_time_manager, only: TIME_manager_Init

    use mod_mktopo, only: MKTOPO_setup
    use mod_mkinit, only: MKINIT_setup

    implicit none

    !----------------------------------------------

    ! namelist compatibility check
    !call ADMIN_versioncheck

    ! setup PROF
    call PROF_setup

    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize', 0)
    
    ! setup constants
    call CONST_setup

    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init( .false. )

    ! setup random number
    call RANDOM_setup

    ! setup tracer index
    call ATMOS_HYDROMETEOR_setup

    ! setup a module for restart file
    call FILE_restart_meshfield_setup

    ! setup submodels
    call  atmos%setup()
    call user_%setup( atmos )

    call atmos%setup_vars()

    ! setup mktopo
    call MKTOPO_setup

    ! setup mkinit
    call MKINIT_setup()

    call PROF_rapend('Initialize', 0)

    return
  end subroutine initialize

  subroutine finalize()
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_finalize
    use scale_time_manager, only: TIME_manager_Final   
    implicit none
    
    !----------------------------------------------
    call PROF_setprefx('FIN')
    call PROF_rapstart('All', 1)

    !-
    call FILE_HISTORY_meshfield_finalize()

    ! finalization submodels
    call  atmos%finalize()

    !-
    call TIME_manager_Final()

    call PROF_rapend  ('All', 1)
    call PROF_rapreport()

    return
  end subroutine finalize

  subroutine set_total_energy( atm_prgvars_manager, & ! (inout)
      atm_auxvars_manager, model_mesh               ) 
  
    use scale_mesh_base3d, only: MeshBase3D
    use scale_localmesh_3d, only: LocalMesh3D
    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use scale_model_var_manager, only: ModelVarManager
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      atm_dyn_dgm_nonhydro3d_common_DRHOT2EnTot, &
      atm_dyn_dgm_nonhydro3d_common_calc_RHOT_hyd

    use mod_atmos_mesh, only: AtmosMesh
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars

    implicit none

    class(ModelVarManager), intent(inout) :: atm_prgvars_manager
    class(ModelVarManager), intent(inout) :: atm_auxvars_manager
    class(AtmosMesh), target, intent(in) :: model_mesh

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, THERM
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: Rtot, CPtot, CVtot

    integer :: n
    integer :: ke
    class(LocalMesh3D), pointer :: lcmesh3D
    class(MeshBase3D), pointer :: mesh

    real(RP), allocatable :: DRHOT_save(:,:)
    real(RP), allocatable :: RHOT_hyd_save(:,:)
    !---------------------------------------------------------------------------

    mesh => model_mesh%ptr_mesh
    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, &
        mesh, atm_prgvars_manager, atm_auxvars_manager, &
        DDENS, MOMX, MOMY, MOMZ, THERM,                 &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot,         &
        lcmesh3D                                        )
      
      allocate( DRHOT_save(lcmesh3D%refElem3D%Np,lcmesh3D%NeA) )
      allocate( RHOT_hyd_save(lcmesh3D%refElem3D%Np,lcmesh3D%NeA) )

      !$omp parallel do
      do ke=lcmesh3D%NeS, lcmesh3D%NeE
        DRHOT_save(:,ke) = THERM%val(:,ke)
      end do

      call atm_dyn_dgm_nonhydro3d_common_calc_RHOT_hyd( RHOT_hyd_save, &
        PRES_hyd%val, lcmesh3D, lcmesh3D%refElem3D )
      
      call atm_dyn_dgm_nonhydro3d_common_DRHOT2EnTot( THERM%val, &
        DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT_save,     &
        DENS_hyd%val, PRES_hyd%val, RHOT_hyd_save,               &
        Rtot%val, CVtot%val, CPtot%val,                          &
        lcmesh3D, lcmesh3D%refElem3D                             )
  
      deallocate( DRHOT_save )
    end do

    return
  end subroutine set_total_energy

  subroutine restart_write
    implicit none    
    !----------------------------------------

    if ( atmos%isActivated() ) call atmos%vars%Write_restart_file()

    return
  end subroutine restart_write

end module mod_dg_prep