!-------------------------------------------------------------------------------
!> Program spectral analysis tool (SCALE-DG)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          spectral analysis tool for models based DG methods
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
program prg_spectral_analysis
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !  
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort
  use scale_const, only: &
    EPS => CONST_EPS, &
    PI => CONST_PI,   &
    RPlanet => CONST_RADIUS
  !-
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D

  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  !-
  use mod_spectral_analysis, only: &
    spectral_analysis_Init, &
    spectral_analysis_do
  use mod_spectral_analysis_mesh, only: &
    MeshList
  use mod_spectral_analysis_vars, only: &
    vars_list, vars2D_list, var_num_step, &
    vars_read, vars_write,   &
    g_var3D_list, s_var3D,   &
    g_var2D_list, s_var2D
  implicit none
  !-----------------------------------------------------------------------------
  !

  ! MPI parameters

  integer                 :: nprocs                      !< Number of MPI processes for analysis
  integer                 :: myrank                      !< My MPI rank                             
  logical                 :: ismaster                    !< Flag whether this is the master process?
  integer                 :: target_proc_s               !< Start rank of target MPI processes for analysis
  integer                 :: target_proc_e               !< End rank of target MPI processes for analysis
  integer                 :: target_proc_num             !< Number of target MPI processes for analysis

  ! Mesh variables

  type(MeshList), target :: mesh_list
  integer :: LevelNum                  !< Number of target levels for spectral analysis of 3D variables
  real(RP), allocatable :: levels(:)   !< Array of target levels for spectral analysis of 3D variables

  ! Spectral analysis parameters

  integer :: ST_NintGLPt              !< Number of Gauss-Legendre quadrature points used when calculating spectral coefficients (for L2 projection method)
  integer :: ST_NsamplePtPerElem1D    !< Number of equidistant sample points per element when calculating spectral coefficients (for uniform sampling method)
  
  integer :: ST_ks  !< Starting index of k wave number
  integer :: ST_ke  !< Ending index of k wave number
  integer :: ST_ls  !< Starting index of l wave number
  integer :: ST_le  !< Ending index of l wave number

  !-
  integer :: istep !< Time step index for spectral analysis

  ! integer :: m, l

  !-------------------------------------------------------------------


  call initialize()

  do istep=1, var_num_step
    LOG_INFO('Spectral analysis',*) "Read vars: istep=", istep
    call vars_read( istep, mesh_list%mesh3D_list, levels )

    LOG_INFO('Spectral analysis',*) "Spectral transformation"
    call spectral_analysis_do( g_var3D_list, g_var2D_list, &
      size(vars_list), size(vars2D_list), mesh_list%mesh_num_x, mesh_list%mesh_num_y, LevelNum, &
      ST_ks, ST_ke, ST_ls, ST_le, s_var3D, s_var2D )

    LOG_INFO('Spectral analysis',*) "Write vars: istep=", istep
    call vars_write( istep, LevelNum )

    call flush(IO_FID_LOG)  
  end do

  ! do l=0, Mt
  ! do m=0, Mt
  !   if ( m > l ) cycle
  !   LOG_INFO('Check:',*) m, l, ":", s_var2D(1,1,l,m,1), s_var2D(1,1,l,m,2)
  ! end do
  ! end do
  call finalize()

contains
!--
  !> Initialize the spectral analysis tool
  subroutine initialize()
    use scale_prc, only: &
      PRC_MPIstart,         &
      PRC_SINGLECOM_setup,  &
      PRC_ERRHANDLER_setup, &
      PRC_myrank
    use scale_const, only: &
      CONST_setup, &
      RPlanet => CONST_RADIUS
    use scale_file, only: &
      FILE_setup    

    use scale_mesh_base2d, only: &
      MFTYPE2D_XY => MeshBase2D_DIMTYPEID_XY
    use mod_spectral_analysis_vars, only: &
      vars_init
    use mod_spectral_analysis_mesh, only: &
      MeshList
    implicit none

    integer :: ierr
    integer :: comm                     !< MPI communicator   (execution)
    character(len=H_LONG) :: cnf_fname  !< Name of the configuration file

    character(len=H_LONG) :: in_bs_filebase !< Base name of the input file saving basic state
    character(len=H_LONG) :: in_filebase    !< Base name of the input file saving 3D and 2D variables
    character(len=H_LONG) :: out_filebase   !< Base name of the output file saving spectral analysis results

    integer :: target_proc_num_tot          !< Total number of MPI processes which were used for the original simulation

    integer, parameter :: Layer_nmax = 1000
    real(RP) :: TARGET_LEVELS(Layer_nmax)              !< Target levels for spectral analysis of 3D variables
    character(len=H_SHORT) :: level_units              !< Units of target levels (e.g., Pa or m)

    integer, parameter :: VarNum_nmax = 20
    character(len=H_SHORT)  :: vars(VarNum_nmax) = ''  !< Name of variables
    logical :: vars_full_output_flag = .false.         !< Flag whether to output full 2D spectral data for each variable at certain level

    integer, parameter :: Var2DNum_nmax = 20
    character(len=H_SHORT)  :: vars2D(Var2DNum_nmax) = ''  !< Name of 2D variables
    logical :: vars2D_full_output_flag = .false.           !< Flag whether to output full 2D spectral data for each variable

    logical :: KinEnergyAnalysisFlag = .false.                   !< Flag whether to output kinetic energy spectrum
    character(len=H_SHORT) :: KinEnergyAnalysisVarPostfix = ''   !< Postfix of variable name (U, V, and W) used for kinetic energy analysis

    character(len=H_MID) :: SpectralAnalysisType = 'Uniform_Sampling' ! Type of spectral analysis (How to calculate spectral coefficients: 'Uniform_Sampling' or 'L2_Projection')
    character(len=H_MID) :: ZInterpType_name  = 'SamplingUniPt_LinearIntrp' !< Type of vertical interpolation ('SamplingUniPt_LinearIntrp' or 'Nearest')

    !-
    namelist / PARAM_SPECTRAL_ANALYSIS / &
        in_bs_filebase,              &
        in_filebase,                 &
        out_filebase,                &
        target_proc_num_tot,         &
        LevelNum,                    &
        TARGET_LEVELS,               &
        level_units,                 &
        ZInterpType_name,            &
        VARS,                        &
        vars_full_output_flag,       &
        VARS2D,                      &
        vars2D_full_output_flag,     &
        KinEnergyAnalysisFlag,       &
        KinEnergyAnalysisVarPostfix, &
        SpectralAnalysisType,        &
        ST_ks, ST_ke,                &
        ST_ls, ST_le,                &
        ST_NsamplePtPerElem1D,       &
        ST_NintGLPt
    
    integer :: Ne2D
    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase2D), pointer :: elem2D

    integer :: m
    integer :: ke2D
    integer :: k
    !-----------------------------------------------------------------------

    ! Start MPI
    call PRC_MPIstart( comm ) ! [OUT]

    ! Setup MPI communicator
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
                              nprocs,  & ! [OUT]
                              myrank,  & ! [OUT]
                              ismaster ) ! [OUT]
  
    call PRC_ERRHANDLER_setup( use_fpm = .false., & ! [IN]
                               master  = .false.  ) ! [IN]

    cnf_fname = IO_ARG_getfname( ismaster )

    ! Setup standard I/O
    call IO_setup( "Spectral Analysis", cnf_fname )

    ! Setup log
    call IO_LOG_setup( myrank, ismaster )
    ! Setup constants
    call CONST_setup
    ! Setup file I/O
    call FILE_setup( myrank )
    
    LOG_NEWLINE
    LOG_INFO("Spectral Analysis",*) 'Setup'

    !- Read configuration file

    in_filebase    = "./in_data/topo"
    in_bs_filebase = ""
    out_filebase   = "./out_data/topo"
    LevelNum    =  1
    level_units = "Pa"
   
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SPECTRAL_ANALYSIS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("Spectral analysis setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("Spectral analysis setup",*) 'Not appropriate names in namelist PARAM_SPECTRAL_ANALYSIS. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_SPECTRAL_ANALYSIS)

    !- Setup target levels for 3D variables
    
    if ( LevelNum > 0 ) then
      allocate( levels(LevelNum) )
      do k=1, LevelNum
        levels(k) = TARGET_LEVELS(k)
      end do
    end if

    !- Calculate target MPI processes for spectral analysis

    target_proc_num = target_proc_num_tot / nprocs
    target_proc_s = myrank * target_proc_num
    target_proc_e = ( myrank + 1 ) * target_proc_num - 1
    LOG_INFO("Spectral analysis setup",*) 'target_proc_num_tot:', target_proc_num_tot
    LOG_INFO("Spectral analysis setup",*) 'target_proc_num:', target_proc_num
    LOG_INFO("Spectral analysis setup",*) 'target_proc:', target_proc_s, target_proc_e

    !- Initialize mesh

    call mesh_list%Init( target_proc_num, target_proc_s )

    elem2D => mesh_list%mesh3D_list(1)%refElem2D
    lcmesh2D => mesh_list%mesh3D_list(1)%mesh2D%lcmesh_list(1)
    Ne2D = lcmesh2D%Ne

    !- Inialize variables and spectral analysis tool

    call vars_init( &
      vars2D, vars2D_full_output_flag,                       &
      vars, vars_full_output_flag,                           &
      KinEnergyAnalysisFlag,                                 &
      KinEnergyAnalysisVarPostfix,                           &
      levels, level_units, ZInterpType_name,                 &
      elem2D%Np, Ne2D, mesh_list%mesh3D_list, target_proc_s, &
      in_filebase, in_bs_filebase, out_filebase, myrank,     &
      ST_ks, ST_ke, ST_ls, ST_le )

    call spectral_analysis_Init( SpectralAnalysisType, &
      ST_ks, ST_ke, ST_ls, ST_le, size(vars_list), LevelNum, size(vars2D_list), &
      mesh_list%mesh3D_list(1), ST_NintGLPt, ST_NsamplePtPerElem1D )
    
    !- Prepare test data
    
    do m=1, target_proc_num
      elem2D => mesh_list%mesh3D_list(m)%refElem2D
      lcmesh2D => mesh_list%mesh3D_list(m)%mesh2D%lcmesh_list(1)
  
      !$omp parallel do
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
        !** Y^0_0
        ! g_var2D(1,:,ke2D,m) = sqrt(1.0_RP / (4.0_RP*PI))
        !* Y^0_1
        ! g_var2D(1,1,:,ke2D,m) = sqrt(3.0_RP / (4.0_RP*PI)) &
        !  * sin(lcmesh2D%lat(:,ke2D)) 
        !** 1/2 * ( Y^1_2 + Y*^1_2 ) 
        ! g_var2D(1,:,ke2D,m) = - sqrt(15.0_RP / (8.0_RP*PI)) &
        !   * sin(lcmesh2D%lat(:,ke2D)) * cos(lcmesh2D%lat(:,ke2D)) &
        !   * cos(lcmesh2D%lon(:,ke2D))
        !** - i 1/2 * ( Y^1_2 - Y*^1_2 )
        ! g_var2D(1,:,ke2D,m) = - sqrt(15.0_RP / (8.0_RP*PI)) &
        !   * sin(lcmesh2D%lat(:,ke2D)) * cos(lcmesh2D%lat(:,ke2D)) &
        !   * sin(lcmesh2D%lon(:,ke2D))
        !** 1/2 * ( Y^2_3 + Y*^2_3 ) 
        ! g_var2D(1,:,ke2D,m) = sqrt(105.0_RP / (32.0_RP*PI)) &
        !   * cos(lcmesh2D%lat(:,ke2D)) * cos(lcmesh2D%lat(:,ke2D)) &
        !   * sin(lcmesh2D%lat(:,ke2D))                             &
        !   * cos(2.0_RP * lcmesh2D%lon(:,ke2D))
        !** 1/2 * ( Y^5_6 + Y*^5_6 ) 
        ! g_var2D(1,:,ke2D,m) = - 3.0_RP / 32.0_RP * sqrt(10001_RP / PI) &
        !   * cos(lcmesh2D%lat(:,ke2D)) **5                         &
        !   * sin(lcmesh2D%lat(:,ke2D))                             &
        !   * cos(5.0_RP * lcmesh2D%lon(:,ke2D))
      end do
    end do

    !-
    LOG_INFO("Spectral analysis setup",*) "Successfully initialized spectral analysis tool"
    call flush(IO_FID_LOG)    
    return
  end subroutine initialize

  !> Finalize the spectral analysis tool
!OCL SERIAL  
  subroutine finalize()
    use scale_prc, only: &
      PRC_mpibarrier, &
      PRC_MPIfinish
    use scale_file, only: &
      FILE_close_all
    !-
    use mod_spectral_analysis_vars, only: vars_final
    implicit none
    !-----------------------------------------------------------------------

    LOG_INFO("final",*) "Final"

    ! Finalize variable module
    call vars_final()
    ! Finalize mesh module
    call mesh_list%Final()

    ! Stop MPI
    call PRC_mpibarrier
    call PRC_MPIfinish 

    return
  end subroutine finalize

end program prg_spectral_analysis