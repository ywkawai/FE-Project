#include "scalelib.h"
program prg_sh_transform
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort
  use scale_const, only: &
    EPS => CONST_EPS, &
    PI => CONST_PI,   &
    RPlanet => CONST_RADIUS

  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D

  use mod_spectral_transform, only: &
    spectral_tranform
  
  use mod_grid, only: &
    lon2D_intrp, lat2D_intrp, Gsqrt2D, J2D, mesh3D_list, &
    Ne2D, refElem2D,                         &
    IntrpMat2D, Intrp_intw2D
  use mod_vars, only: &
    vars_list, vars2D_list, var_num_step, &
    vars_read, vars_write,   &
    g_var3D, s_var3D,        &
    g_var2D, s_var2D
  implicit none

  ! MPI parameters
  integer                 :: nprocs                      ! number of processes               (execution)
  integer                 :: myrank                      ! my rank                           (execution)
  logical                 :: ismaster                    ! master process?                   (execution)
  integer                 :: target_proc_s
  integer                 :: target_proc_e
  integer                 :: target_proc_num_tot
  integer                 :: target_proc_num

  integer :: Mt  
  integer :: LevelNum
  real(RP), allocatable :: levels(:)

  integer :: num_step
  integer :: istep

  integer :: m, l
  !--------

  call initialize()

  do istep=1, var_num_step
    LOG_INFO('SH_Transform',*) "Read vars: istep=", istep
    call vars_read( istep, mesh3D_list, levels )

    LOG_INFO('SH_Transform',*) "Spectral transformation"
    call spectral_tranform( g_var3D, g_var2D, &
      lon2D_intrp, lat2D_intrp, Gsqrt2D, J2D, mesh3D_list,    &
      size(vars_list), size(vars2D_list), LevelNum, refElem2D, Ne2D, Mt, target_proc_num, &
      IntrpMat2D, intrp_intw2D, size(intrp_intw2D),                                       &
      s_var3D, s_var2D )

    call vars_write( istep, Mt, LevelNum )

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
    use mod_vars, only: &
      vars_init
    use mod_grid, only: &
      grid_init
    implicit none

    integer :: ierr
    integer :: comm                     ! communicator   (execution)
    logical :: fileexist
    character(len=H_LONG) :: cnf_fname  ! config file for launcher

    character(len=H_LONG) :: in_bs_filebase
    character(len=H_LONG) :: in_filebase
    character(len=H_LONG) :: out_filebase

    integer :: m
    integer :: ke2D

    integer :: k
    integer, parameter :: Layer_nmax = 1000
    real(RP) :: TARGET_LEVELS(Layer_nmax)
    character(len=H_SHORT) :: level_units

    integer, parameter :: VarNum_nmax = 20
    character(len=H_SHORT)  :: vars(VarNum_nmax) = ''       ! name of variables

    integer, parameter :: Var2DNum_nmax = 20
    character(len=H_SHORT)  :: vars2D(VarNum_nmax) = ''       ! name of variables

    logical :: KinEnergyAnalysisFlag = .false. 

    !-
    namelist / PARAM_SH_TRANSFORM / &
        in_bs_filebase,           &
        in_filebase,              &
        out_filebase,             &
        target_proc_num_tot,      &
        Mt,                       &
        LevelNum,                 &
        TARGET_LEVELS,            &
        level_units,              &
        VARS, VARS2D,             &
        KinEnergyAnalysisFlag
    
    integer :: Ne2D
    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase2D), pointer :: elem2D
    class(ElementBase3D), pointer :: elem3D
    !-----------------------------------------------------------------------

    ! start MPI
    call PRC_MPIstart( comm ) ! [OUT]

    ! setup MPI communicator
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
                              nprocs,  & ! [OUT]
                              myrank,  & ! [OUT]
                              ismaster ) ! [OUT]
  
    call PRC_ERRHANDLER_setup( use_fpm = .false., & ! [IN]
                               master  = .false.  ) ! [IN]

    cnf_fname = IO_ARG_getfname( ismaster )

    ! setup standard I/O
    call IO_setup( "SH_TRANSFORM", cnf_fname )

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
    ! setup constants
    call CONST_setup
    ! setup fie I/O
    call FILE_setup( myrank )
    
    LOG_NEWLINE
    LOG_INFO("SH_TRANSFORM",*) 'Setup'

    !-
    in_filebase    = "./in_data/topo"
    in_bs_filebase = ""
    out_filebase   = "./out_data/topo"
    Mt       = 2
    LevelNum = 1
    level_units = "Pa"
   
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SH_TRANSFORM,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("SH_TRANSFORM_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("SH_TRANSFORM_setup",*) 'Not appropriate names in namelist PARAM_SH_TRANSFORM. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_SH_TRANSFORM)

    !-
    if ( LevelNum > 0 ) then
      allocate( levels(LevelNum) )
      do k=1, LevelNum
        levels(k) = TARGET_LEVELS(k)
      end do
    end if

    !--
    target_proc_num = target_proc_num_tot / nprocs
    target_proc_s = myrank * target_proc_num
    target_proc_e = ( myrank + 1 ) * target_proc_num - 1
    LOG_INFO("SH_TRANSFORM_setup",*) 'target_proc_num_tot:', target_proc_num_tot
    LOG_INFO("SH_TRANSFORM_setup",*) 'target_proc_num:', target_proc_num
    LOG_INFO("SH_TRANSFORM_setup",*) 'target_proc:', target_proc_s, target_proc_e

    call grid_init( target_proc_num, target_proc_s )

    elem2D => mesh3D_list(1)%refElem2D
    lcmesh2D => mesh3D_list(1)%mesh2D%lcmesh_list(1)
    Ne2D = lcmesh2D%Ne
    call vars_init( vars2D, vars, KinEnergyAnalysisFlag, &
      levels, level_units, elem2D%Np, Ne2D, Mt, mesh3D_list, target_proc_s, &
      in_filebase, in_bs_filebase, out_filebase, myrank )

    !--- Test data
    
    do m=1, target_proc_num
      elem2D => mesh3D_list(m)%refElem2D
      lcmesh2D => mesh3D_list(m)%mesh2D%lcmesh_list(1)
  
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
    
    call flush(IO_FID_LOG)  
    return
  end subroutine initialize

!OCL SERIAL  
  subroutine finalize()
    use scale_prc, only: &
      PRC_mpibarrier, &
      PRC_MPIfinish
    use scale_file, only: &
      FILE_close_all

    use mod_grid, only: grid_final
    use mod_vars, only: vars_final
    implicit none
    !-----------------------------------------------------------------------

    LOG_INFO("final",*) "FInal"

    !
    call vars_final()
    call grid_final()

    ! stop MPI
    call PRC_mpibarrier
    call PRC_MPIfinish 

    return
  end subroutine finalize

end program prg_sh_transform
