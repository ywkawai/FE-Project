!-------------------------------------------------------------------------------
!> Program A sample program: 2-dimensional poisson equation on the sphere
!! This equation is solved with a multigrid method
!!
!! The Poisson equation, lap_S q = f, on a cubed-sphere mesh is multiplied by Gsqrt,
!! and the right-hand side, Gsqrt * f, is passed to the MG solver.
!! Since the sphere has no boundary, the (Gsqrt-weighted) mean of f is removed before solving,
!! and that of q is removed after solving.
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
program test_poisson2dGlobal_multigrid
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

  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D

  use scale_meshfield_base, only: MeshField2D

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_put,           &
    FILE_HISTORY_meshfield_write
  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate

  use scale_mesh_hierarchy_2d, only: MeshHierarchy2D
  use scale_multigrid_solver_2d, only: MultiGridSolver2D
  use mod_poisson2dGlobal_smoother, only: MGSmoother_Poisson2D
  use mod_poisson2dGlobal_numerror, only: Poisson2DGlobalNumErrorAnalysis

  implicit none
  !-----------------------------------------------------------------------------

  type(QuadrilateralElement)  :: refElem
  type(MeshCubedSphereDom2D), target :: mesh

  type(MeshField2D) :: q
  type(MeshField2D) :: qexact
  type(MeshField2D) :: f
  type(MeshField2D) :: rhs     !< Right-hand side multiplied by Gsqrt (passed to the MG solver)

  integer, save :: HST_ID(3)

  !-
  type(MeshHierarchy2D) :: mesh_hierarchy
  type(MGSmoother_Poisson2D) :: smoother
  type(MultiGridSolver2D) :: mg_solver

  !-
  logical :: Do_NumErrorAnalysis            !< Flag whether analysis of numerical error is performed
  type(Poisson2DGlobalNumErrorAnalysis) :: numerror_analysis

  !-----------------------------------------------------------------------------

  call init()
  !-
  call FILE_HISTORY_meshfield_put( HST_ID(2), qexact )
  call FILE_HISTORY_meshfield_put( HST_ID(3), f )

  call PROF_rapstart( "MG_solve", 0 )
  call mg_solver%Solve( q, &
    rhs )
  call PROF_rapend( "MG_solve", 0 )

  ! The solution is determined up to a constant: remove its mean
  call remove_mean( q, "q" )

  if ( Do_NumErrorAnalysis ) &
    call numerror_analysis%Eval( qexact, q, 1, 0.0_RP )

  call FILE_HISTORY_meshfield_put( HST_ID(1), q )
  call FILE_HISTORY_meshfield_write()
  !-
  call final()

contains
  !> Set the exact solution and the right-hand side
!OCL SERIAL
  subroutine set_rhs_lc( f_, qexact_, lmesh, elem, test_case_name, radius )
    use mod_poisson2dGlobal_numerror, only: poisson2dGlobal_get_exact_sol
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: f_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qexact_(elem%Np,lmesh%NeA)
    character(len=*), intent(in) :: test_case_name
    real(RP), intent(in) :: radius

    integer :: ke
    real(RP) :: lap_fac
    !-----------------------------------------

    do ke=lmesh%NeS, lmesh%NeE
      call poisson2dGlobal_get_exact_sol( qexact_(:,ke), lap_fac, & ! (out)
        test_case_name, lmesh%lon(:,ke), lmesh%lat(:,ke), radius, elem%Np ) ! (in)
      f_(:,ke) = lap_fac * qexact_(:,ke)
    end do

    return
  end subroutine set_rhs_lc

  !> Calculate the surface integral of a field over the local mesh: sum_K 1^T M ( J Gsqrt val )
!OCL SERIAL
  function integrate_lc( val, lmesh, elem ) result( val_int )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(in) :: val(elem%Np,lmesh%NeA)
    real(RP) :: val_int

    integer :: ke
    real(RP) :: int_w(elem%Np)
    !-----------------------------------------

    ! 1^T M (integration weights consistent with the mass matrix of the reference element)
    int_w(:) = sum( elem%M(:,:), dim=1 )

    val_int = 0.0_RP
    !$omp parallel do reduction(+:val_int)
    do ke=lmesh%NeS, lmesh%NeE
      val_int = val_int &
        + sum( int_w(:) * lmesh%J(:,ke) * lmesh%Gsqrt(:,ke) * val(:,ke) )
    end do
    return
  end function integrate_lc

  !> Calculate the global surface integral of a field
!OCL SERIAL
  function integrate_global( field ) result( val_int )
    use mpi
    implicit none
    type(MeshField2D), intent(in) :: field
    real(RP) :: val_int

    integer :: ldom
    real(RP) :: val_int_lc
    integer :: ierr
    !-----------------------------------------

    val_int_lc = 0.0_RP
    do ldom=1, mesh%LOCAL_MESH_NUM
      val_int_lc = val_int_lc + integrate_lc( field%local(ldom)%val, mesh%lcmesh_list(ldom), refElem )
    end do
    call MPI_Allreduce( val_int_lc, val_int, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr )
    return
  end function integrate_global

  !> Remove the (Gsqrt-weighted) mean of a field over the sphere
!OCL SERIAL
  subroutine remove_mean( field, label )
    implicit none
    type(MeshField2D), intent(inout) :: field
    character(len=*), intent(in) :: label

    type(MeshField2D) :: one
    real(RP) :: area, field_mean
    integer :: ldom
    class(LocalMesh2D), pointer :: lmesh
    integer :: ke
    !-----------------------------------------

    call one%Init( "one", "1", mesh )
    do ldom=1, mesh%LOCAL_MESH_NUM
      one%local(ldom)%val(:,:) = 1.0_RP
    end do
    area = integrate_global( one )
    call one%Final()

    field_mean = integrate_global( field ) / area

    do ldom=1, mesh%LOCAL_MESH_NUM
      lmesh => mesh%lcmesh_list(ldom)
      !$omp parallel do
      do ke=lmesh%NeS, lmesh%NeE
        field%local(ldom)%val(:,ke) = field%local(ldom)%val(:,ke) - field_mean
      end do
    end do

    LOG_INFO("remove_mean",'(a,a,a,ES15.7,a,ES15.7,a)') "Mean of ", trim(label), " = ", field_mean, &
      " (area=", area / ( 4.0_RP * PI * mesh%RPlanet**2 ), " x 4 pi a^2) is removed."
    return
  end subroutine remove_mean

  !> Initialization
!OCL SERIAL
  subroutine init()
    use scale_const, only: &
      CONST_setup,         &
      RPlanet => CONST_RADIUS
    use scale_time_manager, only:           &
      TIME_manager_Init,                    &
      TIME_manager_report_timeintervals
    use scale_calendar, only: CALENDAR_setup
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_setup
    use scale_file_history, only: FILE_HISTORY_reg

    use mod_poisson2dGlobal_smoother, only: &
      MGSmoother_Poisson2D_AUX_SCALAR_NUM, &
      MGSmoother_Poisson2D_AUX_VEC_NUM
    implicit none

    !-
    integer            :: NeGX                = 2
    integer            :: NeGY                = 2
    integer            :: PolyOrder           = 1
    integer            :: NLocalMeshPerPrc    = 1
    character(len=H_SHORT) :: TestCaseName    = 'Y22'  !< 'Y22': q=cos^2(lat)cos(2lon), 'Y10': q=sin(lat)
    logical, parameter :: LumpedMassMatFlag   = .false.

    namelist /PARAM_TEST/ &
      NeGX, NeGY, PolyOrder, NLocalMeshPerPrc, &
      TestCaseName, Do_NumErrorAnalysis

    !-
    integer, parameter :: pMG_LV_LIST_MAX = 16
    integer :: P_LEVEL_NUM
    integer :: P_LEVEL_LIST(pMG_LV_LIST_MAX)

    integer, parameter :: hMG_LV_LIST_MAX = 16
    integer :: H_LEVEL_NUM
    integer :: NeGX_list(hMG_LV_LIST_MAX)
    integer :: NeGY_list(hMG_LV_LIST_MAX)

    integer :: MG_Vcyc_Num_Max = 10
    real(RP) :: MG_Threshold_Ratio_Residual_L2 = 1.0e-6_RP
    real(RP) :: MG_Threshold_Residual_L2  = 1.0E-12_RP
    real(RP) :: MG_Threshold_Residual_Max = 1.0E-12_RP

    namelist /PARAM_Poisson2D_MG/ &
      H_LEVEL_NUM, NeGX_list, NeGY_list,                  &
      P_LEVEL_NUM, P_LEVEL_LIST,                          &
      MG_Vcyc_Num_Max, MG_Threshold_Ratio_Residual_L2,    &
      MG_Threshold_Residual_L2, MG_Threshold_Residual_Max

    !-
    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr

    class(LocalMesh2D), pointer :: lmesh
    integer :: ldom
    integer :: ke

    character(len=H_MID) :: conf_name
    !-----------------------------------------------------------------------------

    !-- setup MPI
    call PRC_MPIstart( comm )

    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]

    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]

    !-- setup scale_io
    conf_name = IO_ARG_getfname( ismaster )
    call IO_setup( "test_poisson2dGlobal_multigrid", conf_name )

    !-- setup log
    call IO_LOG_setup( myrank, ismaster )

    !-- read namelist
    Do_NumErrorAnalysis = .true.

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
    read(IO_FID_CONF,nml=PARAM_Poisson2D_MG,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("Poisson2DGlobal_MG_Init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("Poisson2DGlobal_MG_Init",*) 'Not appropriate names in namelist PARAM_Poisson2D_MG. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_Poisson2D_MG)

    !-- setup profiler
    call PROF_setup
    call PROF_rapstart( "total", 0 )
    call PROF_rapstart( "init", 1 )

    !-- setup constants (CONST_RADIUS can be changed by PARAM_CONST)
    call CONST_setup

    !-- setup calendar & time manager

    call CALENDAR_setup
    call TIME_manager_Init

    !-- setup reference element and spatial operators

    call refElem%Init(PolyOrder, LumpedMassMatFlag)

    !-- setup mesh

    call mesh%Init( &
      NeGX, NeGY, RPlanet,                    &
      refElem, NLocalMeshPerPrc, nproc=nprocs )

    call mesh%Generate()

    !-- setup fields

    call q%Init( "q", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )
    call f%Init( "f", "1", mesh )
    call rhs%Init( "rhs", "1", mesh )

    !-- setup history files

    call FILE_HISTORY_meshfield_setup( meshcubedsphere2D_=mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='XY')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='XY')
    call FILE_HISTORY_reg( f%varname, "f", f%unit, HST_ID(3), dim_type='XY')

    !-- report information of time intervals
    call TIME_manager_report_timeintervals

    !- Setup an object to manage the mesh hierarchy

    LOG_INFO("Poisson2DGlobal_MG_Init",*) "Setup mesh hierarchy and MG solver"
    call mesh_hierarchy%Init( mesh, &
      P_LEVEL_LIST(1:P_LEVEL_NUM), P_LEVEL_NUM,                       &
      NeGX_list(1:H_LEVEL_NUM), NeGY_list(1:H_LEVEL_NUM), H_LEVEL_NUM )

    !- Setup an object to provide a smoother
    LOG_INFO("Poisson2DGlobal_MG_Init",*) "Setup smoother"
    call smoother%Init( mesh, P_LEVEL_LIST(1:P_LEVEL_NUM) )

    !- Setup an object to provide a MG solver
    LOG_INFO("Poisson2DGlobal_MG_Init",*) "Setup MG solver"
    call mg_solver%Init( mesh_hierarchy, smoother, &
      MGSmoother_Poisson2D_AUX_SCALAR_NUM, MGSmoother_Poisson2D_AUX_VEC_NUM )

    call mg_solver%Set_vcycle_parameter( MG_Vcyc_Num_Max, &
      MG_Threshold_Ratio_Residual_L2,  &
      MG_Threshold_Residual_L2,        &
      MG_Threshold_Residual_Max )

    !- Set the right-hand side and initial guess
    do ldom=1, mesh%LOCAL_MESH_NUM
      call set_rhs_lc( f%local(ldom)%val, qexact%local(ldom)%val, mesh%lcmesh_list(ldom), refElem, &
        TestCaseName, mesh%RPlanet )
      q%local(ldom)%val(:,:) = 0.0_RP
    end do

    ! Discrete compatibility condition: sum_K 1^T M ( J Gsqrt f ) = 0
    call remove_mean( f, "f" )
    call remove_mean( qexact, "qexact" )

    ! Right-hand side multiplied by Gsqrt
    do ldom=1, mesh%LOCAL_MESH_NUM
      lmesh => mesh%lcmesh_list(ldom)
      rhs%local(ldom)%val(:,:) = 0.0_RP
      do ke=lmesh%NeS, lmesh%NeE
        rhs%local(ldom)%val(:,ke) = lmesh%Gsqrt(:,ke) * f%local(ldom)%val(:,ke)
      end do
    end do

    !- Setup an object to evaluate numerical errors
    if ( Do_NumErrorAnalysis ) &
      call numerror_analysis%Init( TestCaseName, mesh, refElem )

    call PROF_rapend( "init", 1 )
    LOG_INFO("Poisson2DGlobal_MG_Init",*) "The setup has succeeded."
    return
  end subroutine init

!> Finalization
  subroutine final()
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    implicit none
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )

    call FILE_HISTORY_meshfield_finalize()

    if ( Do_NumErrorAnalysis ) &
      call numerror_analysis%Final()

    call smoother%Final()
    call mg_solver%Final()
    call mesh_hierarchy%Final()

    call q%Final()
    call qexact%Final()
    call f%Final()
    call rhs%Final()

    call mesh%Final()
    call refElem%Final()

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()
    return
  end subroutine final
end program test_poisson2dGlobal_multigrid
