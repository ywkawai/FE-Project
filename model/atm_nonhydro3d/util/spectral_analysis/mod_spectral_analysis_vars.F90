!-------------------------------------------------------------------------------
!> module Variable containers for spectral analysis
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_spectral_analysis_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort

  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  
  use scale_meshfield_base, only: MeshField2D
  use scale_file_base_meshfield, only: FILE_base_meshfield

  use mod_spectral_analysis_vintrp, only: SpectralAnalysisVIntrp

  use netcdf
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: vars_init
  public :: vars_final
  public :: vars_read
  public :: vars_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public, allocatable :: vars2D_list(:) !< Array of 2D variable names
  character(len=H_SHORT), public, allocatable :: vars_list(:)   !< Array of variable names
  integer, public :: var_num_step !< Number of time steps for spectral analysis

  integer, public :: var2DNum   !< Number of 2D variables for spectral analysis
  integer, public :: varNum0    !< Number of 3D variables for spectral analysis (excluding kinetic energy analysis variables)
  integer, public :: varNum     !< Total number of 3D variables for spectral analysis (including kinetic energy analysis variables)

  type(MeshField2D), public, allocatable :: g_var2D_list(:,:)
  real(RP), public, allocatable :: s_var2D(:,:,:,:)

  type(MeshField2D), public, allocatable :: g_var3D_list(:,:,:)
  real(RP), public, allocatable :: s_var3D(:,:,:,:,:)

  type(MeshField2D), allocatable :: g_bs_var3D_list(:,:,:)
  integer, public, parameter :: BS_PRESHYD_ID = 1
  integer, public, parameter :: BS_DENSHYD_ID = 2

  integer, public, parameter :: KE_RsqrtUmet_ID = 1
  integer, public, parameter :: KE_RsqrtVmet_ID = 2
  integer, public, parameter :: KE_RsqrtWmet_ID = 3
  integer, public, parameter :: KE_U_ID         = 4
  integer, public, parameter :: KE_V_ID         = 5
  integer, public, parameter :: KE_W_ID         = 6

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
    
  integer :: ncido                     !< NetCDF file ID for output
  integer :: k_wave_number_dim_id      !< Dimension ID for k wave number in output NetCDF file
  integer :: l_wave_number_dim_id      !< Dimension ID for l wave number in output NetCDF file
  integer :: oneD_wave_number_dim_id   !< Dimension ID for 1D wave number in output NetCDF file
  integer :: level_dim_id              !< Dimension ID for vertical level in output NetCDF file
  integer :: time_dim_id               !< Dimension ID for time in output NetCDF file
  integer :: k_wave_number_id          !< Variable ID for k wave number in output NetCDF file
  integer :: l_wave_number_id          !< Variable ID for l wave number in output NetCDF file
  integer :: oneD_wave_number_id       !< Variable ID for 1D wave number in output NetCDF file
  integer :: level_id                  !< Variable ID for vertical level in output NetCDF file
  integer :: time_id                   !< Variable ID for time in output NetCDF file

  integer, allocatable :: var2D_id_list_r(:)     !< 2D Variable IDs for output NetCDF file (real part) 
  integer, allocatable :: var2D_id_list_i(:)     !< 2D Variable IDs for output NetCDF file (imaginary part)
  integer, allocatable :: var2D_id_spectra_1d(:) !< 2D Variable IDs for 1D spectra output NetCDF file
  logical :: vars2D_full_output_flag             !< Flag whether to output full 2D spectral data for each 2D variable

  integer, allocatable :: var_id_list_r(:)       !< 3D Variable IDs for output NetCDF file (real part)
  integer, allocatable :: var_id_list_i(:)       !< 3D Variable IDs for output NetCDF file (imaginary part)
  integer, allocatable :: var_id_spectra_1d(:)   !< 3D Variable IDs for 1D spectra output NetCDF file
  logical :: vars_full_output_flag               !< Flag whether to output full 2D spectral data for each 3D variable

  integer :: kin_energy_spectra_h_id             !< ID for horizontal kinetic energy spectra output NetCDF file
  integer :: kin_energy_spectra_v_id             !< ID for vertical kinetic energy spectra output NetCDF file
  integer :: kin_energy_spectra_1d_h_id          !< ID for horizontal kinetic energy 1D spectra output NetCDF file
  integer :: kin_energy_spectra_1d_v_id          !< ID for vertical kinetic energy 1D spectra output NetCDF file

  type(FILE_base_meshfield), allocatable :: in_files(:)    !< Array of objects to read 3D variables from input files
  type(FILE_base_meshfield), allocatable :: in_bs_files(:) !< Array of objects to read 3D variables from input files (basic state)
  type(FILE_base_meshfield), allocatable :: in_files_2d(:) !< Array of objects to read 2D variables from input files

  real(DP) :: output_dt
  real(DP) :: start_sec

  logical :: output_flag

  logical :: KinEnergyAnalysisFlag
  character(len=H_SHORT) :: KEAnalysisVarPostfix

  integer :: s_ks  !< Starting index of k wave number
  integer :: s_ke  !< Ending index of k wave number
  integer :: s_ls  !< Starting index of l wave number
  integer :: s_le  !< Ending index of l wave number

  type(SpectralAnalysisVIntrp) :: VIntrpTool !< Object to manage vertical interpolation for spectral analysis

contains
  !> Setup variables for spectral analysis
!OCL SERIAL
  subroutine vars_init( &
    vars2D, vars2D_full_output_flag_,                     &
    vars, vars_full_output_flag_,                         &
    KinEnergyAnalysisFlag_, KinEnergyAnalysisVarPostfix_, &
    levels, level_units, zintrp_type_name,                &
    Np2D, Ne2D, mesh3D_list, target_proc_s,               &
    in_fbasename, in_bs_fbasename, out_fbasename, myrank, &
    s_ks_, s_ke_, s_ls_, s_le_ )

    use scale_const, only: EPS => CONST_EPS
    implicit none
    character(len=H_SHORT), intent(in) :: vars2D(:)  !< Array of 2D variable names
    logical, intent(in) :: vars2D_full_output_flag_  !< Flag whether to output full 2D spectral data for each 2D variable
    character(len=H_SHORT), intent(in) :: vars(:)    !< Array of 3D variable names
    logical, intent(in) :: vars_full_output_flag_    !< Flag whether to output full 2D spectral data for each 3D variable
    logical, intent(in) :: KinEnergyAnalysisFlag_                       !< Flag whether to perform kinetic energy spectrum analysis
    character(len=H_SHORT), intent(in) :: KinEnergyAnalysisVarPostfix_  !< Postfix for variables used in kinetic energy analysis
    real(RP), intent(in) :: levels(:)                                   !< Array of vertical levels
    character(len=*), intent(in) :: level_units                         !< Units of vertical levels
    character(len=*), intent(in) :: zintrp_type_name                    !< Name of vertical interpolation type
    integer, intent(in) :: Np2D                                         !< Number of nodes in 2D element
    integer, intent(in) :: Ne2D                                         !< Number of elements in 2D mesh
    type(MeshCubeDom3D), intent(in), target :: mesh3D_list(:)
    integer, intent(in) :: target_proc_s
    character(len=*), intent(in) :: in_fbasename
    character(len=*), intent(in) :: in_bs_fbasename
    character(len=*), intent(in) :: out_fbasename
    integer, intent(in) :: myrank
    integer, intent(in) :: s_ks_, s_ke_
    integer, intent(in) :: s_ls_, s_le_
    
    integer :: k_wave_number(s_ks_:s_ke_)
    integer :: l_wave_number(s_ls_:s_le_)
    integer :: oneD_wave_number(0:s_ke_)
    integer :: k, l, m

    integer :: nn
    integer :: target_myrank

    integer :: LevelNum
    integer :: mesh_num

    real(DP) :: time_endsec

    integer :: v

    class(LocalMesh3D), pointer :: lmesh3D
    !------------------------------------------------------------


    s_ks = s_ks_; s_ke = s_ke_
    s_ls = s_ls_; s_le = s_le_

    vars2D_full_output_flag = vars2D_full_output_flag_
    vars_full_output_flag = vars_full_output_flag_

    KinEnergyAnalysisFlag = KinEnergyAnalysisFlag_

    lmesh3D => mesh3D_list(1)%lcmesh_list(1)

    !- Setup variable names for spectral analysis

    var2DNum = count_named_entries( vars2D )
    allocate( vars2D_list(var2DNum) )
    do nn=1, var2DNum
      vars2D_list(nn) = vars2D(nn)
    end do

    varNum0 = count_named_entries( vars )
    if ( KinEnergyAnalysisFlag ) then
      varNum = varNum0 + 6
    else
      varNum = varNum0
    end if

    allocate( vars_list(varNum) )
    do nn=1, varNum0
      vars_list(nn) = vars(nn)
    end do
    if ( KinEnergyAnalysisFlag ) then
      KEAnalysisVarPostfix = KinEnergyAnalysisVarPostfix_
      vars_list(varNum0+KE_RsqrtUmet_ID) = 'RsqrtUmet'
      vars_list(varNum0+KE_RsqrtVmet_ID) = 'RsqrtVmet'
      vars_list(varNum0+KE_RsqrtWmet_ID) = 'RsqrtW'
      vars_list(varNum0+KE_U_ID) = 'U'//trim(KEAnalysisVarPostfix)
      vars_list(varNum0+KE_V_ID) = 'V'//trim(KEAnalysisVarPostfix)
      vars_list(varNum0+KE_W_ID) = 'W'//trim(KEAnalysisVarPostfix)
    end if

    !- Setup variables for spectral analysis
    
    mesh_num = size(mesh3D_list)
    LevelNum = size(levels)

    ! 2D variables
    if (var2DNum > 0) then
      allocate( g_var2D_list(var2DNum,mesh_num) )
      do v=1, var2DNum
      do m=1, mesh_num
        call g_var2D_list(v,m)%Init(vars2D(v), "", mesh3D_list(m)%mesh2D )
      end do
      end do
      
      allocate( s_var2D(s_ks:s_ke,s_ls:s_le,2,var2DNum) )
      s_var2D(:,:,:,:) = 0.0_RP  
    end if
    ! 3D variables
    if (varNum > 0) then
      allocate( g_var3D_list(LevelNum,varNum,mesh_num) )
      do m=1, mesh_num
      do l=1, LevelNum
      do v=1, varNum
        call g_var3D_list(l,v,m)%Init( trim(vars(v)), "", mesh3D_list(m)%mesh2D )
      end do
      end do
      end do

      allocate( s_var3D(s_ks:s_ke,s_ls:s_le,2,LevelNum,varNum) )
      s_var3D(:,:,:,:,:) = 0.0_RP  
    end if

    ! Basic state variables
    if ( in_bs_fbasename /= "" ) then
      allocate( g_bs_var3D_list(LevelNum,2,mesh_num) )
      do m=1, mesh_num
      do l=1, LevelNum
        call g_bs_var3D_list(l,BS_DENSHYD_ID,m)%Init("DENS_hyd", "", mesh3D_list(m)%mesh2D )
        call g_bs_var3D_list(l,BS_PRESHYD_ID,m)%Init("PRES_hyd", "", mesh3D_list(m)%mesh2D )
      end do
      end do
    end if
    
    !- Setup vertical interpolation strategy used in spectral analysis

    call VIntrpTool%Init( zintrp_type_name, levels, lmesh3D )

    !- Setup objects to read data files -----------

    if ( in_fbasename /= "" ) then
      LOG_INFO("VARS_INIT",*) "Data file:", trim(in_fbasename)
      allocate( in_files(mesh_num) )
      allocate( in_files_2d(mesh_num) )
    end if
    if ( in_bs_fbasename /= "" ) then
      LOG_INFO("VARS_INIT",*) "Data file with basic state:", trim(in_bs_fbasename)
      allocate( in_bs_files(mesh_num) )
    end if

    do m=1, mesh_num
      target_myrank = target_proc_s + m - 1      
      !-
      if ( allocated(in_files) ) then
        call in_files(m)%Init( varNum, mesh3D=mesh3D_list(m) )
        call in_files(m)%Open( in_fbasename, myrank=target_myrank )
      end if
      if ( allocated(in_files_2d) ) then
        call in_files_2d(m)%Init( var2DNum, mesh2D=mesh3D_list(m)%mesh2D )
        call in_files_2d(m)%Open( in_fbasename, myrank=target_myrank )
      end if

      if ( allocated(in_bs_files) ) then
        call in_bs_files(m)%Init( 2, mesh3D=mesh3D_list(m) )
        call in_bs_files(m)%Open( in_bs_fbasename, myrank=target_myrank )  
      end if      
    end do

    ! Get the number of time steps and output interval from the input files
    if ( allocated(in_files) .and. varNum > 0) then
      call in_files(1)%Get_VarStepSize( trim(vars_list(1)),   & ! (in)
        var_num_step                                          ) ! (out)
      call in_files(1)%Get_dataInfo( trim(vars_list(1)), istep=1,  & ! (in)
        time_start=start_sec,                       & ! (out)
        time_end=time_endsec                        ) ! (out)
      output_dt = time_endsec - start_sec 
    end if 
    if ( allocated(in_files_2d) .and. var2DNum > 0) then
      call in_files_2d(1)%Get_VarStepSize( trim(vars2D_list(1)),   & ! (in)
        var_num_step                                               ) ! (out)
      call in_files_2d(1)%Get_dataInfo( trim(vars2D_list(1)), istep=1,  & ! (in)
        time_start=start_sec,                       & ! (out)
        time_end=time_endsec                        ) ! (out)
      output_dt = time_endsec - start_sec 
    end if
    if ( ( .not. allocated(in_files_2d) ) .and. ( .not. allocated(in_files) ) ) then
      output_dt = 0.0_RP; var_num_step = 0
    end if

    if ( output_dt < EPS .and. var_num_step == 0 ) then
      var_num_step   = 1
    end if

    !- Setup wave number arrays
    
    do k=s_ks, s_ke
      k_wave_number(k) = k
    end do
    do k=0, s_ke
      oneD_wave_number(k) = k
    end do
    do l=s_ls, s_le
      l_wave_number(l) = l
    end do

    !- Setup output variables with spectral analysis

    if (myrank==0) then
      output_flag = .true.
    else
      output_flag = .false.
    end if
    if (output_flag) then

      ! Define dimensions in the output NetCDF file
    
      call nc_check( nf90_create( trim(out_fbasename)//".nc", NF90_clobber, ncido) )
      call nc_check( nf90_def_dim( ncido, 'k', s_ke-s_ks+1, k_wave_number_dim_id ) )
      call nc_check( nf90_def_dim( ncido, 'l', s_le-s_ls+1, l_wave_number_dim_id ) )
      call nc_check( nf90_def_dim( ncido, 'K', s_ke+1, oneD_wave_number_dim_id ) )
      if (varNum > 0)  call nc_check( nf90_def_dim( ncido, 'level', LevelNum, level_dim_id ) )
      call nc_check( nf90_def_dim( ncido, 'time', NF90_UNLIMITED, time_dim_id ) )

      call nc_check( nf90_def_var( ncido, 'k', NF90_INT, k_wave_number_dim_id, k_wave_number_id ) )
      call nc_check( nf90_def_var( ncido, 'l', NF90_INT, l_wave_number_dim_id, l_wave_number_id ) )
      call nc_check( nf90_def_var( ncido, 'K', NF90_INT, oneD_wave_number_dim_id, oneD_wave_number_id ) )
      if (varNum > 0)  call nc_check( nf90_def_var( ncido, 'level', NF90_DOUBLE, level_dim_id, level_id ) )
      call nc_check( nf90_def_var( ncido, 'time', NF90_DOUBLE, time_dim_id, time_id ) )

      call nc_check( nf90_put_att( ncido, k_wave_number_id, 'units', 'm-1' ) )
      call nc_check( nf90_put_att( ncido, l_wave_number_id, 'units', 'm-1' ) )
      call nc_check( nf90_put_att( ncido, oneD_wave_number_id, 'units', 'm-1' ) )
      if (varNum > 0)  call nc_check( nf90_put_att( ncido, level_id, 'units', trim(level_units) ) )
      call nc_check( nf90_put_att( ncido, time_id, 'units', 'sec' ) )

      ! Define variables in the output NetCDF file

      if ( varNum0 > 0 ) then
        allocate( var_id_list_r(varNum0), var_id_list_i(varNum0) )
        allocate( var_id_spectra_1d(varNum0) )
        var_id_list_r(:) = -1
        var_id_list_i(:) = -1
      end if
      do nn=1, varNum0
        call define_spectral_output_var( ncido, trim(vars_list(nn)), vars_full_output_flag, & ! (in)
          (/ k_wave_number_dim_id, l_wave_number_dim_id, level_dim_id, time_dim_id /),      & ! (in)
          (/ oneD_wave_number_dim_id, level_dim_id, time_dim_id /),                         & ! (in)
          var_id_list_r(nn), var_id_list_i(nn), var_id_spectra_1d(nn) ) ! (out)
      end do

      if ( var2DNum > 0 ) then
        allocate( var2D_id_list_r(var2DNum), var2D_id_list_i(var2DNum) )
        allocate( var2D_id_spectra_1d(var2DNum) )
        var2D_id_list_r(:) = -1
        var2D_id_list_i(:) = -1
      end if
      do nn=1, var2DNum
        call define_spectral_output_var( ncido, trim(vars2D_list(nn)), vars2D_full_output_flag, & ! (in)
          (/ k_wave_number_dim_id, l_wave_number_dim_id, time_dim_id /),                        & ! (in)
          (/ oneD_wave_number_dim_id, time_dim_id /),                                           & ! (in)
          var2D_id_list_r(nn), var2D_id_list_i(nn), var2D_id_spectra_1d(nn) ) ! (out)
      end do

      ! Define variables for kinetic energy spectra

      if ( KinEnergyAnalysisFlag ) then
        if ( vars_full_output_flag ) then
          call nc_check( nf90_def_var( ncido, "kin_energy_spectra_h", &
            NF90_DOUBLE, (/ k_wave_number_dim_id, l_wave_number_dim_id, level_dim_id, time_dim_id /), kin_energy_spectra_h_id ) )
          call nc_check( nf90_def_var( ncido, "kin_energy_spectra_v", &
            NF90_DOUBLE, (/ k_wave_number_dim_id, l_wave_number_dim_id, level_dim_id, time_dim_id /), kin_energy_spectra_v_id ) )
        else
          kin_energy_spectra_h_id = -1; kin_energy_spectra_v_id = -1
        end if
        call nc_check( nf90_def_var( ncido, "kin_energy_spectra_1d_h", &
          NF90_DOUBLE, (/ oneD_wave_number_dim_id, level_dim_id, time_dim_id /), kin_energy_spectra_1d_h_id ) )
        call nc_check( nf90_def_var( ncido, "kin_energy_spectra_1d_v", &
          NF90_DOUBLE, (/ oneD_wave_number_dim_id, level_dim_id, time_dim_id /), kin_energy_spectra_1d_v_id ) )
      end if
      call nc_check( nf90_enddef(ncido) )
      
      ! Put wave number and level data into the output NetCDF file

      call nc_check( nf90_put_var(ncido, k_wave_number_id, k_wave_number ) )
      call nc_check( nf90_put_var(ncido, l_wave_number_id, l_wave_number ) )
      call nc_check( nf90_put_var(ncido, oneD_wave_number_id, oneD_wave_number ) )
      if (varNum > 0) call nc_check( nf90_put_var(ncido, level_id, levels ) )
    end if

    return
  end subroutine vars_init

  !> Finalize variables for spectral analysis
!OCL SERIAL
  subroutine vars_final()
    implicit none

    integer :: m
    !----------------------------

    deallocate( vars_list, vars2D_list )

    call VIntrpTool%Final()

    if ( allocated( in_files) ) then
      do m=1, size(in_files)
        call in_files(m)%Close()
        call in_files(m)%Final()
      end do
      deallocate( in_files )
    end if
    if ( allocated( in_files_2d) ) then
      do m=1, size(in_files_2d)
        call in_files_2d(m)%Close()
        call in_files_2d(m)%Final()
      end do
      deallocate( in_files_2d )
    end if
    if ( allocated(in_bs_files) ) then
      do m=1, size(in_bs_files)
        call in_bs_files(m)%Close()
        call in_bs_files(m)%Final()
      end do
      deallocate( in_bs_files )
    end if

    if (output_flag) then
      ! Close the output NetCDF file
      call nc_check( nf90_close(ncido) )

      if ( varNum0 > 0 ) then
        deallocate( var_id_list_r, var_id_list_i, var_id_spectra_1d )
      end if
      if ( var2DNum > 0 ) then
        deallocate( var2D_id_list_r, var2D_id_list_i, var2D_id_spectra_1d )
      end if
    end if
    
    return
  end subroutine vars_final

  !> Output variables for spectral analysis
!OCL SERIAL
  subroutine vars_write( istep, LevelNum )
    implicit none
    integer, intent(in) :: istep
    integer, intent(in) :: LevelNum

    integer :: vid
    real(RP) :: s_var3D_r_buf(s_ks:s_ke,s_ls:s_le,LevelNum)
    real(RP) :: s_var3D_i_buf(s_ks:s_ke,s_ls:s_le,LevelNum)

    real(RP) :: s_kin1D_buf(0:s_ke,LevelNum)

    real(RP) :: s_power_buf(s_ks:s_ke,s_ls:s_le,LevelNum)

    integer :: k
    integer :: s_kall, s_lall
    integer :: s_K1Dall
    !----------------------------------

    LOG_INFO("vars_write",*) "istep=", istep !, ":",  start_sec, output_dt, start_sec + istep * output_dt
    if (output_flag) then

      call nc_check( nf90_put_var( ncido, time_id, (/ start_sec + istep * output_dt /), &
        start=[istep], count=[1]) )
      
      s_kall = s_ke - s_ks + 1
      s_lall = s_le - s_ls + 1
      s_K1Dall = s_ke + 1

      !- Output spectral data with 3D variables

      do vid=1, varNum0
        !$omp parallel do
        do k=1, LevelNum
          s_var3D_r_buf(:,:,k) = s_var3D(:,:,1,k,vid)
          s_var3D_i_buf(:,:,k) = s_var3D(:,:,2,k,vid)
          s_power_buf(:,:,k) = s_var3D(:,:,1,k,vid)**2 + s_var3D(:,:,2,k,vid)**2
          call calc_oneD_spectra( s_power_buf(:,:,k), s_kin1D_buf(:,k) )
        end do
        if ( vars_full_output_flag ) then
          call nc_check( nf90_put_var( ncido, var_id_list_r(vid), s_var3D_r_buf, &
            start=[1,1,1,istep], count=[s_kall,s_lall,LevelNum,1]) )
          call nc_check( nf90_put_var( ncido, var_id_list_i(vid), s_var3D_i_buf, &
            start=[1,1,1,istep], count=[s_kall,s_lall,LevelNum,1]) )
        end if
        call nc_check( nf90_put_var( ncido, var_id_spectra_1d(vid), s_kin1D_buf, &
          start=[1,1,istep], count=[s_K1Dall,LevelNum,1]) )
      end do

      !- Output spectral data with 2D variables

      do vid=1, var2DNum
        s_power_buf(:,:,1) = s_var2D(:,:,1,vid)**2 + s_var2D(:,:,2,vid)**2
        call calc_oneD_spectra( s_power_buf(:,:,1), s_kin1D_buf(:,1) )
        if ( vars2D_full_output_flag ) then
          call nc_check( nf90_put_var( ncido, var2D_id_list_r(vid), s_var2D(:,:,1,vid), &
            start=[1,1,istep], count=[s_kall,s_lall,1]) )
          call nc_check( nf90_put_var( ncido, var2D_id_list_i(vid), s_var2D(:,:,2,vid), &
            start=[1,1,istep], count=[s_kall,s_lall,1]) )
        end if
        call nc_check( nf90_put_var( ncido, var2D_id_spectra_1d(vid), s_kin1D_buf(:,1), &
          start=[1,istep], count=[s_K1Dall,1]) )
      end do

      !- Output kinetic energy spectra

      if ( KinEnergyAnalysisFlag ) then
        ! Output horizontal kinetic energy spectra
        !$omp parallel do
        do k=1, LevelNum
          s_power_buf(:,:,k) = 0.5_RP * ( s_var3D(:,:,1,k,varNum0+KE_RsqrtUmet_ID)**2 + s_var3D(:,:,2,k,varNum0+KE_RsqrtUmet_ID)**2 &
                                        + s_var3D(:,:,1,k,varNum0+KE_RsqrtVmet_ID)**2 + s_var3D(:,:,2,k,varNum0+KE_RsqrtVmet_ID)**2 )
          call calc_oneD_spectra( s_power_buf(:,:,k), s_kin1D_buf(:,k) )
        end do
        if ( vars_full_output_flag ) then
          call nc_check( nf90_put_var( ncido, kin_energy_spectra_h_id, s_power_buf, &
            start=[1,1,1,istep], count=[s_kall,s_lall,LevelNum,1]) )
        end if
        call nc_check( nf90_put_var( ncido, kin_energy_spectra_1d_h_id, s_kin1D_buf, &
          start=[1,1,istep], count=[s_K1Dall,LevelNum,1]) )
        
        ! Output vertical kinetic energy spectra
        !$omp parallel do
        do k=1, LevelNum
          s_power_buf(:,:,k) = 0.5_RP * ( s_var3D(:,:,1,k,varNum0+KE_RsqrtWmet_ID)**2 + s_var3D(:,:,2,k,varNum0+KE_RsqrtWmet_ID)**2 )
          call calc_oneD_spectra( s_power_buf(:,:,k), s_kin1D_buf(:,k) )
        end do
        if ( vars_full_output_flag ) then
          call nc_check( nf90_put_var( ncido, kin_energy_spectra_v_id, s_power_buf, &
            start=[1,1,1,istep], count=[s_kall,s_lall,LevelNum,1]) )
        end if
        call nc_check( nf90_put_var( ncido, kin_energy_spectra_1d_v_id, s_kin1D_buf, &
          start=[1,1,istep], count=[s_K1Dall,LevelNum,1]) )
      end if
    end if

    return
  contains
    !> Calculate ring-averaged 1D spectra from 2D spectra
    subroutine calc_oneD_spectra( s_2D, s_1D )
      implicit none
      real(RP), intent(in) :: s_2D(s_ks:s_ke,s_ls:s_le)
      real(RP), intent(out) :: s_1D(0:s_ke)
      
      integer :: k1D, kk, ll
      real(RP) :: dist
      !---------------------------------------------

      s_1D(:) = 0.0_RP
      do k1D=0, s_ke
        do ll=s_ls, s_le
        do kk=s_ks, s_ke
          dist = sqrt( real(kk**2 + ll**2, kind=RP) )
          if ( k1D - 0.5_RP < dist .and. dist <= k1D + 0.5_RP ) then
            s_1D(k1D) = s_1D(k1D) + s_2D(kk,ll)
          end if
        end do
        end do
      end do
      return
    end subroutine calc_oneD_spectra
  end subroutine vars_write

  !> Read variables from input files for spectral analysis
!OCL SERIAL
  subroutine vars_read( istep, mesh3D_list, levels )
    use scale_meshfield_base, only: MeshField2D, MeshField3D
    use scale_mesh_base2d, only: &
      MF2D_XYT => MeshBase2D_DIMTYPEID_XYT
    use scale_mesh_base3d, only: &
      MF3D_XYZT => MeshBase3D_DIMTYPEID_XYZT
    implicit none
    integer, intent(in) :: istep
    type(MeshCubeDom3D), intent(in), target :: mesh3D_list(:)
    real(RP), intent(in) :: levels(:)

    integer :: m
    integer :: vid
    character(len=H_SHORT) :: varname
    integer :: ke2D
    integer :: ph
    integer :: k
    integer :: ldomID

    class(LocalMesh3D), pointer :: lcmesh3D
    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase3D), pointer :: elem3D

    ! real(RP), allocatable :: tmp_U(:,:)
    ! real(RP), allocatable :: tmp_V(:,:)
    ! real(RP), allocatable :: tmp_Umet(:,:)
    ! real(RP), allocatable :: tmp_Vmet(:,:)
    ! real(RP), allocatable :: gam(:,:)
    integer :: LayerNum
    integer :: vid_U, vid_V, vid_W

    type(MeshField2D) :: g_tmp_ddens(size(levels))
    real(RP), allocatable :: RsqrtDENS(:,:,:)
    !----------------------------

    LayerNum = size(levels)
    vid_U = -1; vid_V = -1; vid_W = -1

    do m=1, size(mesh3D_list)
      lcmesh2D => mesh3D_list(m)%mesh2D%lcmesh_list(1)
      do vid=1, var2DNum
        varname = trim(vars2D_list(vid))
        call in_files_2d(m)%Read_Var( &
          MF2D_XYT, varname, g_var2D_list(vid,m), step=istep )
      end do
    end do

    do m=1, size(mesh3D_list)
      lcmesh3D => mesh3D_list(m)%lcmesh_list(1)
      elem3D => lcmesh3D%refElem3D
      lcmesh2D => mesh3D_list(m)%mesh2D%lcmesh_list(1)

      if ( allocated(in_bs_files) .and. istep == 1 ) then
        call get_gvar3D( istep, in_bs_files(m), "PRES_hyd", mesh3D_list(m), levels, &
          g_bs_var3D_list(:,BS_PRESHYD_ID,m) )
        call get_gvar3D( istep, in_bs_files(m), "DENS_hyd", mesh3D_list(m), levels, &
          g_bs_var3D_list(:,BS_DENSHYD_ID,m) )
      end if

      do vid=1, varNum0
        varname = trim(vars_list(vid))
        call get_gvar3D( istep, in_files(m), varname, mesh3D_list(m), levels, &
          g_var3D_list(:,vid,m) )

        if ( trim(varname) == "U"//trim(KEAnalysisVarPostfix) ) vid_U = vid
        if ( trim(varname) == "V"//trim(KEAnalysisVarPostfix) ) vid_V = vid
        if ( trim(varname) == "W"//trim(KEAnalysisVarPostfix) ) vid_W = vid
      end do

      if ( KinEnergyAnalysisFlag ) then
        if (vid_U < 0 ) then
          vid_U = varNum0 + KE_U_ID
          call get_gvar3D( istep, in_files(m), "U"//trim(KEAnalysisVarPostfix), mesh3D_list(m), levels, g_var3D_list(:,vid_U,m) )
        end if
        if (vid_V < 0 ) then
          vid_V = varNum0 + KE_V_ID
          call get_gvar3D( istep, in_files(m), "V"//trim(KEAnalysisVarPostfix), mesh3D_list(m), levels, g_var3D_list(:,vid_V,m) )
        end if
        if (vid_W < 0 ) then
          vid_W = varNum0 + KE_W_ID
          call get_gvar3D( istep, in_files(m), "W"//trim(KEAnalysisVarPostfix), mesh3D_list(m), levels, g_var3D_list(:,vid_W,m) )
        end if

        !- Tentative implementation of coordinate conversion for global model data (Cubed-sphere to lon-lat) for kinetic energy analysis
        ! allocate( tmp_U(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        ! allocate( tmp_V(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        ! allocate( tmp_Umet(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        ! allocate( tmp_Vmet(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        ! allocate( gam(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        !
        ! gam(:,:) = 1.0_RP
        ! do k = 1, LayerNum
        !   tmp_U(:,:) = g_var3D(vid_U,k,:,:,m)
        !   tmp_V(:,:) = g_var3D(vid_V,k,:,:,m)
        !
        !   call CubedSphereCoordCnv_CS2LonLatVec( lcmesh2D%panelID, lcmesh2D%pos_en(:,:,1), lcmesh2D%pos_en(:,:,2), gam(:,:), &
        !     lcmesh2D%Ne*lcmesh2D%refElem2D%Np, tmp_U, tmp_V, &
        !     tmp_Umet, tmp_Vmet )
        !  
        !   g_var3D(vid_U,k,:,:,m) = tmp_Umet(:,:)
        !   g_var3D(vid_V,k,:,:,m) = tmp_Vmet(:,:)
        ! end do
        !---------------------------------------------------------------------------------
        
        do k=1, LayerNum
          call g_tmp_ddens(k)%Init("DDENS"//trim(KEAnalysisVarPostfix), "", mesh3D_list(m)%mesh2D )
        end do
        allocate( RsqrtDENS(lcmesh2D%refElem2D%Np,lcmesh2D%Ne,LayerNum) )

        call get_gvar3D( istep, in_files(m), "DDENS"//trim(KEAnalysisVarPostfix), mesh3D_list(m), levels, &
          g_tmp_ddens )

        !$omp parallel private(k,ke2D,ph)
        !$omp do collapse(2)
        do k=1, LayerNum
        do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
          do ph=1, lcmesh2D%refElem2D%Np
            RsqrtDENS(ph,ke2D,k) = sqrt( g_bs_var3D_list(k,BS_DENSHYD_ID,m)%local(1)%val(ph,ke2D) + g_tmp_ddens(k)%local(1)%val(ph,ke2D) )
          end do
        end do
        end do
        !$omp do collapse(2)
        do k=1, LayerNum
          do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
          do ph=1, lcmesh2D%refElem2D%Np  
            g_var3D_list(k,varNum0+KE_RsqrtUmet_ID,m)%local(1)%val(ph,ke2D) = RsqrtDENS(ph,ke2D,k) * g_var3D_list(k,vid_U,m)%local(1)%val(ph,ke2D)
            g_var3D_list(k,varNum0+KE_RsqrtVmet_ID,m)%local(1)%val(ph,ke2D) = RsqrtDENS(ph,ke2D,k) * g_var3D_list(k,vid_V,m)%local(1)%val(ph,ke2D)
            g_var3D_list(k,varNum0+KE_RsqrtWmet_ID,m)%local(1)%val(ph,ke2D) = RsqrtDENS(ph,ke2D,k) * g_var3D_list(k,vid_W,m)%local(1)%val(ph,ke2D)
          end do
          end do
        end do
        !$omp end parallel

        do k=1, LayerNum
          call g_tmp_ddens(k)%Final()
        end do
        deallocate( RsqrtDENS )

        ! deallocate( tmp_U, tmp_V, tmp_Umet, tmp_Vmet, gam )
      end if
    end do
    return
  end subroutine vars_read

!-- private subroutines ---------------------------------------------------

!OCL SERIAL
  subroutine get_gvar3D( istep, in_file, varname, mesh3D, levels, &
    gvar_out )
    use scale_meshfield_base, only: MeshField3D
    use scale_mesh_base3d, only: &
      MeshBase3D, &
      MF3D_XYZT => MeshBase3D_DIMTYPEID_XYZT
    implicit none
    integer, intent(in) :: istep
    class(FILE_base_meshfield), intent(inout) :: in_file
    character(len=*), intent(in) :: varname
    class(MeshBase3D), intent(in), target :: mesh3D
    real(RP), intent(in) :: levels(:)
    type(MeshField2D), intent(inout) :: gvar_out(size(levels))

    type(MeshField3D) :: tmp_field3D
    class(LocalMesh3D), pointer :: lcmesh3D

    integer :: k
    !---------------------------------------------------

    call tmp_field3D%Init( varname, "", mesh3D )
    call in_file%Read_Var( MF3D_XYZT, varname, tmp_field3D, step=istep )
    
    lcmesh3D => mesh3D%lcmesh_list(1)
    do k = 1, size(levels)
      call VIntrpTool%GetLevelValue( k, levels(k), lcmesh3D, lcmesh3D%refElem3D, tmp_field3D, & ! (in)
        gvar_out(k) ) ! (inout)
    end do

    call tmp_field3D%Final()
    return
  end subroutine get_gvar3D

  !> Count the number of non-empty entries at the head of a namelist name array
!OCL SERIAL
  pure function count_named_entries( names ) result( n )
    implicit none
    character(len=*), intent(in) :: names(:)
    integer :: n
    integer :: nn
    !---------------------------------------------
    n = 0
    do nn = 1, size(names)
      if ( names(nn) == '' ) exit
      n = n + 1
    end do
    return
  end function count_named_entries

  !-
  ! Define spectral output variables in the output NetCDF file
!OCL SERIAL
  subroutine define_spectral_output_var( ncid, varname, full_flag, full_dimids, oneD_dimids, &
      id_r, id_i, id_1d )
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varname
    logical, intent(in) :: full_flag
    integer, intent(in) :: full_dimids(:)
    integer, intent(in) :: oneD_dimids(:)

    integer, intent(out) :: id_r, id_i, id_1d
    !---------------------------------------------

    id_r = -1; id_i = -1

    if ( full_flag ) then
      call nc_check( nf90_def_var( ncid, trim(varname)//"_r", NF90_DOUBLE, full_dimids, id_r ) )
      call nc_check( nf90_def_var( ncid, trim(varname)//"_i", NF90_DOUBLE, full_dimids, id_i ) )
    end if
    call nc_check( nf90_def_var( ncid, trim(varname)//"_1D", NF90_DOUBLE, oneD_dimids, id_1d ) )

    return
  end subroutine define_spectral_output_var

!OCL SERIAL
  subroutine nc_check( status )
    implicit none
    integer, intent (in) :: status
    !---------------------------------------------
    if(status /= nf90_noerr) then 
      LOG_INFO('NetCDF_check',*) trim(nf90_strerror(status))
      call flush(IO_FID_LOG)
      call PRC_abort
    end if

    return
  end subroutine nc_check

end module mod_spectral_analysis_vars