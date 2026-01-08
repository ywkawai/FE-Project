!-------------------------------------------------------------------------------
!> Program A sample program: 3-dimensional poisson equation
!! This equation is solved with a multigrid method
!! 
!! 
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
program test_multigrid3d
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

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_hierarchy_3d, only: MeshHierarchy3D

  use scale_meshfield_base, only: MeshField3d

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_put,           &
    FILE_HISTORY_meshfield_write
  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate

  use scale_mesh_hierarchy_3d, only: MeshHierarchy3D
  use scale_multigrid_solver_3d, only: MultiGridSolver3D
  use mod_poisson3d_smoother, only: MGSmoother_Poisson3D

  implicit none
  !-----------------------------------------------------------------------------

  type(HexahedralElement)  :: refElem
  type(MeshCubeDom3D), target :: mesh

  type(MeshField3D) :: q
  type(MeshField3D) :: qexact

  type(MeshField3D) :: f

  integer, save :: HST_ID(3)

  !-
  type(MeshHierarchy3D) :: mesh_hierarchy
  type(MGSmoother_Poisson3D) :: smoother
  type(MultiGridSolver3D) :: mg_solver

  !-----------------------------------------------------------------------------

  call init()
  !-
  call FILE_HISTORY_meshfield_put( HST_ID(2), qexact )
  call FILE_HISTORY_meshfield_put( HST_ID(3), f )

  call mg_solver%Solve( q, &
    f )

  call FILE_HISTORY_meshfield_put( HST_ID(1), q )
  call FILE_HISTORY_meshfield_write()
  !-
  call final()

contains
!OCL SERIAL
  subroutine set_rhs_lc( f_, qexact_, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: f_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qexact_(elem%Np,lmesh%NeA)

    integer :: ke
    real(RP) :: x(elem%Np)
    real(RP) :: y(elem%Np)
    real(RP) :: z(elem%Np)
    !-----------------------------------------

    do ke=lmesh%NeS, lmesh%NeE
      x(:) = lmesh%pos_en(:,ke,1)
      y(:) = lmesh%pos_en(:,ke,2)
      z(:) = lmesh%pos_en(:,ke,3)
      
      qexact_(:,ke) = - sin( 2.0_RP * PI * ( x(:) - 0.25_RP) ) &
                      * sin( 2.0_RP * PI * ( y(:) - 0.25_RP) )
      f_(:,ke) = - qexact_(:,ke) * (2.0_RP*PI)**2 * 2.0_RP

    end do
    return
  end subroutine set_rhs_lc

  !> Initialization
!OCL SERIAL
  subroutine init() 
    use scale_time_manager, only:           &
      TIME_manager_Init,                    &
      TIME_manager_report_timeintervals
    use scale_const, only: CONST_setup
    use scale_calendar, only: CALENDAR_setup
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_setup  
    use scale_file_history, only: FILE_HISTORY_reg  

    use mod_poisson3d_smoother, only: &
      MGSmoother_Possion3D_AUX_SCALAR_NUM, &
      MGSmoother_Possion3D_AUX_HVEC_NUM
    implicit none

    !-
    real(RP), parameter :: dom_xmin =   0.0_RP
    real(RP), parameter :: dom_xmax = + 1.0_RP
    real(RP), parameter :: dom_ymin =   0.0_RP
    real(RP), parameter :: dom_ymax = + 1.0_RP
    real(RP), parameter :: dom_zmin =   0.0_RP
    real(RP), parameter :: dom_zmax = + 1.0_RP

    integer            :: NprcX               = 2
    integer            :: NprcY               = 2
    integer            :: NeGX                = 2
    integer            :: NeGY                = 2
    integer            :: NeGZ                = 2
    integer            :: PolyOrder           = 1
    logical, parameter :: LumpedMassMatFlag   = .false.
    character(len=H_SHORT) :: TINTEG_SCHEME_TYPE

    namelist /PARAM_TEST/ &
      NprcX, NprcY,                &
      NeGX, NeGY, NeGZ, PolyOrder

    !-
    integer, parameter :: pMG_LV_LIST_MAX = 16
    integer :: P_LEVEL_NUM
    integer :: P_LEVEL_LIST(pMG_LV_LIST_MAX)


    integer, parameter :: hMG_LV_LIST_MAX = 16
    integer :: H_LEVEL_NUM
    integer :: NeGX_list(hMG_LV_LIST_MAX)
    integer :: NeGY_list(hMG_LV_LIST_MAX)
    integer :: NeGZ_list(hMG_LV_LIST_MAX)

    integer :: MG_Vcyc_Num_Max = 10
    real(RP) :: MG_Threshold_Ratio_Residual_L2 = 1.0e-6_RP
    real(RP) :: MG_Threshold_Residual_L2  = 1.0E-12_RP
    real(RP) :: MG_Threshold_Residual_Max = 1.0E-12_RP

    namelist /PARAM_Poisson3d_MG/ &
      H_LEVEL_NUM, NeGX_list, NeGY_list, NeGZ_list,       &
      P_LEVEL_NUM, P_LEVEL_LIST,                          &
      MG_Vcyc_Num_Max, MG_Threshold_Ratio_Residual_L2,    &
      MG_Threshold_Residual_L2, MG_Threshold_Residual_Max

    !-
    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr
    
    class(LocalMesh3D), pointer :: lmesh    
    integer :: ldom

    character(len=H_MID) :: conf_name
    !-----------------------------------------------------------------------------

    !-- setup MPI
    call PRC_MPIstart( comm )

    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    !-- setup scale_io
    conf_name = IO_ARG_getfname( ismaster )
    call IO_setup( "test_multigrid3d", conf_name )
    
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
    
    !-- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_Poisson3d_MG,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("Poisson3D_MG_Init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("Poisson3D_MG_Init",*) 'Not appropriate names in namelist PARAM_Poisson3d_MG. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_Poisson3D_MG)

    !-- setup profiler
    call PROF_setup
    call PROF_rapstart( "total", 0 )
    call PROF_rapstart( "init", 1 )

    ! setup constants
    call CONST_setup

    !-- setup calendar & time manager

    call CALENDAR_setup
    call TIME_manager_Init

    !-- setup reference element and spatial operators

    call refElem%Init( PolyOrder, PolyOrder, LumpedMassMatFlag )

    !-- setup mesh

    call mesh%Init( &
      NeGX, NeGY, NeGZ,                                           &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
      .true., .true., .true.,                                     &
      refElem, 1, NprcX, NprcY )
    call mesh%Generate()
    
    !-- setup fields

    call q%Init( "q", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )

    call f%Init( "f", "1", mesh )


    !-- setup history files

    call FILE_HISTORY_meshfield_setup( mesh3D_=mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='XYZ')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='XYZ')
    call FILE_HISTORY_reg( f%varname, "f", f%unit, HST_ID(3), dim_type='XYZ')

    !-- report information of time intervals
    
    call TIME_manager_report_timeintervals

    !- Setup an object to manage the mesh hierarchy

    LOG_INFO("Poisson3d_MG_Init",*) "Setup mesh hierarchy and MG solver"
    call mesh_hierarchy%Init( mesh, &
      P_LEVEL_LIST(1:P_LEVEL_NUM), P_LEVEL_NUM,                                                 &
      NeGX_list(1:H_LEVEL_NUM), NeGY_list(1:H_LEVEL_NUM), NeGZ_list(1:H_LEVEL_NUM), H_LEVEL_NUM )

    !- Setup an object to provide a smoother
    LOG_INFO("Poisson3d_MG_Init",*) "Setup smoother"
    call smoother%Init( mesh, P_LEVEL_LIST(1:P_LEVEL_NUM) )
    
    !- Setup an object to provide a MG solver   
    LOG_INFO("Poisson3d_MG_Init",*) "Setup MG solver"
    call mg_solver%Init( mesh_hierarchy, smoother, &
      MGSmoother_Possion3D_AUX_SCALAR_NUM, MGSmoother_Possion3D_AUX_HVEC_NUM )

    call mg_solver%Set_vcycle_parameter( MG_Vcyc_Num_Max, &
      MG_Threshold_Ratio_Residual_L2,  &
      MG_Threshold_Residual_L2,        &
      MG_Threshold_Residual_Max )

    !- Set the right-hand side and initial guess
    do ldom=1, mesh%LOCAL_MESH_NUM
      call set_rhs_lc( f%local(ldom)%val, qexact%local(ldom)%val, mesh%lcmesh_list(ldom), refElem )
      q%local(ldom)%val(:,:) = 0.0_RP
    end do

    call PROF_rapend( "init", 1 )
    LOG_INFO("Poisson3d_MG_Init",*) "The setup has succeeded."    
    return
  end subroutine init

!> Finalization
  subroutine final()
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    ! use scale_time_manager, only: TIME_manager_Final 
    implicit none
    integer :: idom
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )

    call FILE_HISTORY_meshfield_finalize()

    call smoother%Final()
    call mg_solver%Final()
    call mesh_hierarchy%Final()
    
    call q%Final()
    call qexact%Final()

    call f%Final()

    call mesh%Final()    
    call refElem%Final()

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()
    return
  end subroutine final
end program test_multigrid3d