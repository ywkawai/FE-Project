#include "scaleFElib.h"
program test_multigrid2d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_const, only: &
    PI  => CONST_PI
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof

  use scale_element_base, only: ElementBase1D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_hierarchy_2d, only: MeshHierarchy2D

  use scale_meshfield_base, only: MeshField2D

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_put,           &
    FILE_HISTORY_meshfield_write
  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate


  implicit none
  !-----------------------------------------------------------------------------

  type(QuadrilateralElement)  :: refElem
  type(MeshRectDom2D), target :: mesh
  type(MeshHierarchy2D) :: mesh_hierarchy

  type(MeshField2D) :: q
  type(MeshField2D) :: qexact
  integer, save :: HST_ID(2)

  !-----------------------------------------------------------------------------

  call init()
  call final()

contains
!OCL SERIAL

  !> Initialization
  subroutine init() 
    use scale_time_manager, only:           &
      TIME_manager_Init,                    &
      TIME_manager_report_timeintervals
    use scale_calendar, only: CALENDAR_setup
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_setup  
    use scale_file_history, only: FILE_HISTORY_reg  
    use scale_polynominal, only: Polynominal_GenLagrangePoly

    use mod_poisson2d_mg, only: Poisson2d_mg_Init
    implicit none

    real(RP), parameter :: dom_xmin =   0.0_RP
    real(RP), parameter :: dom_xmax = + 1.0_RP
    real(RP), parameter :: dom_ymin =   0.0_RP
    real(RP), parameter :: dom_ymax = + 1.0_RP

    integer            :: NprcX               = 2
    integer            :: NprcY               = 2
    integer            :: NeGX                = 2
    integer            :: NeGY                = 2
    integer            :: PolyOrder           = 1
    logical, parameter :: LumpedMassMatFlag   = .false.
    character(len=H_SHORT) :: TINTEG_SCHEME_TYPE

    integer, parameter :: pMG_LV_LIST_MAX = 16
    integer :: P_LEVEL_NUM
    integer :: P_LEVEL_LIST(pMG_LV_LIST_MAX)

    integer, parameter :: hMG_LV_LIST_MAX = 16
    integer :: H_LEVEL_NUM
    integer :: NeGX_list(hMG_LV_LIST_MAX)
    integer :: NeGY_list(hMG_LV_LIST_MAX)

    namelist /PARAM_TEST/ &
      NprcX, NprcY,                      &
      NeGX, NeGY, PolyOrder,             &
      H_LEVEL_NUM, NeGX_list, NeGY_list, &
      P_LEVEL_NUM, P_LEVEL_LIST

    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr
    
    class(LocalMesh2D), pointer :: lmesh    
    integer :: idom

    character(len=H_MID) :: conf_name
    !-----------------------------------------------------------------------------

    !-- setup MPI
    call PRC_MPIstart( comm )

    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    !-- setup scale_io
    conf_name = IO_ARG_getfname( ismaster )
    call IO_setup( "test_multigrid2d", conf_name )
    
    !-- setup log
    call IO_LOG_setup( myrank, ismaster )   

    !-- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TEST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_TEST. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_TEST)
    
    !-- setup profiler
    call PROF_setup
    call PROF_rapstart( "total", 0 )
    call PROF_rapstart( "init", 1 )

    !-- setup calendar & time manager

    call CALENDAR_setup
    call TIME_manager_Init

    !-- setup reference element and spatial operators

    call refElem%Init(PolyOrder, LumpedMassMatFlag)

    !-- setup mesh

    call mesh%Init( &
      NeGX, NeGY,                             &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
      .true., .true.,                         &
      refElem, 1, NprcX, NprcY )
    call mesh%Generate()

    call mesh_hierarchy%Init( mesh, &
      P_LEVEL_LIST, P_LEVEL_NUM,                                      &
      NeGX_list(1:H_LEVEL_NUM), NeGY_list(1:H_LEVEL_NUM), H_LEVEL_NUM )
      
    !-- setup fields

    call q%Init( "q", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )


    !-- setup history files    

    call FILE_HISTORY_meshfield_setup( mesh2D_=mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='XY')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='XY')

    !-- report information of time intervals
    call TIME_manager_report_timeintervals

    !- setup poisson solver
    call Poisson2d_mg_Init( mesh_hierarchy )

    return
  end subroutine init

  !> Finalization
  subroutine final()
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    ! use scale_time_manager, only: TIME_manager_Final 
    
    use mod_poisson2d_mg, only: Poisson2d_mg_Final
    implicit none
    integer :: idom
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )

    call FILE_HISTORY_meshfield_finalize()

    call Poisson2d_mg_Final()
    
    call q%Final()
    call qexact%Final()

    call mesh_hierarchy%Final()
    call mesh%Final()    
    call refElem%Final()

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()
    return
  end subroutine final  

end program test_multigrid2d