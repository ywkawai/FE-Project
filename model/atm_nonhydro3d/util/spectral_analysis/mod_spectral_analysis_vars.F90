#include "scaleFElib.h"
module mod_spectral_analysis_vars
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort

  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  
  use scale_meshfield_base, only: &
    MeshField2D

  use scale_file_base_meshfield, only: FILE_base_meshfield

  use netcdf
  implicit none
  private

  public :: vars_init
  public :: vars_final
  public :: vars_read
  public :: vars_write

  character(len=H_SHORT), public, allocatable :: vars2D_list(:)
  character(len=H_SHORT), public, allocatable :: vars_list(:)
  integer, public :: var_num_step

  integer, public :: var2DNum
  integer, public :: varNum0
  integer, public :: varNum

  real(RP), public, allocatable :: g_var2D(:,:,:,:)
  type(MeshField2D), public, allocatable :: g_var2D_list(:,:)
  real(RP), public, allocatable :: s_var2D(:,:,:,:)

  real(RP), public, allocatable :: g_var3D(:,:,:,:,:)
  type(MeshField2D), public, allocatable :: g_var3D_list(:,:,:)
  real(RP), public, allocatable :: s_var3D(:,:,:,:,:)

  real(RP), public, allocatable :: g_bs_var3D(:,:,:,:,:)
  type(MeshField2D), allocatable :: g_bs_var3D_list(:,:,:)
  integer, public, parameter :: BS_PRESHYD_ID = 1
  integer, public, parameter :: BS_DENSHYD_ID = 2

  real(RP), public, allocatable :: g_ke_var3D(:,:,:,:,:)
  integer, public, parameter :: KE_RsqrtUmet_ID = 1
  integer, public, parameter :: KE_RsqrtVmet_ID = 2
  integer, public, parameter :: KE_RsqrtWmet_ID = 3

  integer :: ncido
  integer :: k_wave_number_dim_id
  integer :: l_wave_number_dim_id
  integer :: oneD_wave_number_dim_id
  integer :: level_dim_id
  integer :: time_dim_id
  integer :: k_wave_number_id
  integer :: l_wave_number_id
  integer :: oneD_wave_number_id
  integer :: level_id
  integer :: time_id
  integer, allocatable :: var2D_id_list_r(:)
  integer, allocatable :: var2D_id_list_i(:)
  integer, allocatable :: var_id_list_r(:)
  integer, allocatable :: var_id_list_i(:)
  integer :: kin_energy_spectra_h_id
  integer :: kin_energy_spectra_v_id
  integer :: kin_energy_spectra_1d_h_id
  integer :: kin_energy_spectra_1d_v_id

  type(FILE_base_meshfield), allocatable :: in_files(:)
  type(FILE_base_meshfield), allocatable :: in_bs_files(:)
  type(FILE_base_meshfield), allocatable :: in_files_2d(:)

  real(DP) :: output_dt
  real(DP) :: start_sec

  logical :: output_flag
  logical :: KinEnergyAnalysisFlag

  integer :: s_ks, s_ke
  integer :: s_ls, s_le

  real(RP), allocatable :: VIntrpMat(:,:,:)
  real(RP), allocatable :: VIntrp_coef1(:)
  integer, allocatable :: VIntrp_keZ1(:), VIntrp_keZ2(:)

  integer, parameter :: ZintrpType_Nearest_ID = 1
  integer, parameter :: ZintrpType_SamplingUniPt_LinInterp_ID = 2
  integer :: ZintrpType_ID

contains
!OCL SERIAL
  subroutine vars_init( &
    vars2D, vars, KinEnergyAnalysisFlag_,                            &
    levels, level_units, zintrp_type_name,                           &
    Np2D, Ne2D, mesh3D_list, target_proc_s,                          &
    in_fbasename, in_bs_fbasename, out_fbasename, myrank,            &
    s_ks_, s_ke_, s_ls_, s_le_ )

    use scale_const, only: EPS => CONST_EPS
    implicit none
    character(len=H_SHORT), intent(in) :: vars2D(:)
    character(len=H_SHORT), intent(in) :: vars(:)
    logical, intent(in) :: KinEnergyAnalysisFlag_
    real(RP), intent(in) :: levels(:)
    character(len=*), intent(in) :: level_units
    character(len=*), intent(in) :: zintrp_type_name
    integer, intent(in) :: Np2D
    integer, intent(in) :: Ne2D
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
    KinEnergyAnalysisFlag = KinEnergyAnalysisFlag_
    !-

    var2DNum = 0
    do nn= 1, size(vars2D)
      if ( vars2D(nn) == '' ) then
        exit
      else
        var2DNum = var2DNum + 1
      end if
    end do
    allocate( vars2D_list(var2DNum) )
    do nn=1, var2DNum
      vars2D_list(nn) = vars2D(nn)
    end do

    varNum0 = 0
    do nn= 1, size(vars)
      if ( vars(nn) == '' ) then
        exit
      else
        varNum0 = varNum0 + 1
      end if
    end do
    varNum = varNum0

    if ( KinEnergyAnalysisFlag ) varNum = varNum0 + 3

    allocate( vars_list(varNum) )
    do nn=1, varNum0
      vars_list(nn) = vars(nn)
    end do
    if ( KinEnergyAnalysisFlag ) then
      vars_list(varNum0+1) = 'RsqrtUmet'
      vars_list(varNum0+2) = 'RsqrtVmet'
      vars_list(varNum0+3) = 'RsqrtW'
    end if

    !--
    mesh_num = size(mesh3D_list)
    LevelNum = size(levels)
    if (var2DNum > 0) then
      allocate( g_var2D(var2DNum,Np2D,Ne2D,mesh_num) )
      allocate( g_var2D_list(var2DNum,mesh_num) )
      do v=1, var2DNum
      do m=1, mesh_num
        call g_var2D_list(v,m)%Init(vars2D(v), "", mesh3D_list(m)%mesh2D )
      end do
      end do
      
      allocate( s_var2D(s_ks:s_ke,s_ls:s_le,2,var2DNum) )
      s_var2D(:,:,:,:) = 0.0_RP

      allocate( var2D_id_list_r(var2DNum) )
      allocate( var2D_id_list_i(var2DNum) )  
    end if
    if (varNum > 0) then
      allocate( g_var3D(varNum,LevelNum,Np2D,Ne2D,mesh_num) )
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

      allocate( var_id_list_r(varNum) )
      allocate( var_id_list_i(varNum) )  
    end if

    !---
    if ( in_fbasename /= "" ) then
      allocate( in_files(mesh_num) )
      allocate( in_files_2d(mesh_num) )
    end if
    if ( in_bs_fbasename /= "" ) then
      LOG_INFO("VARS_INIT",*) "Data file with basic state:", trim(in_bs_fbasename)
      allocate( in_bs_files(mesh_num) )
      allocate( g_bs_var3D(2,LevelNum,Np2D,Ne2D,mesh_num) )

      allocate( g_bs_var3D_list(LevelNum,2,mesh_num) )
      do m=1, mesh_num
      do l=1, LevelNum
        call g_bs_var3D_list(l,BS_DENSHYD_ID,m)%Init("DENS_hyd", "", mesh3D_list(m)%mesh2D )
        call g_bs_var3D_list(l,BS_PRESHYD_ID,m)%Init("PRES_hyd", "", mesh3D_list(m)%mesh2D )
      end do
      end do
    end if
    
    !-------------
    select case( trim(zintrp_type_name) )
    case( "Nearest" )
      ZintrpType_ID = ZintrpType_Nearest_ID
    case( "SamplingUniPt_LinearIntrp" )
      ZintrpType_ID = ZintrpType_SamplingUniPt_LinInterp_ID
    case default
      LOG_INFO("SpectralAnalysis_vars_Init",*) "The specified zintrp_type_name is invalid. Check!", trim(zintrp_type_name)
      call PRC_abort
    end select

    lmesh3D => mesh3D_list(1)%lcmesh_list(1)

    allocate( VIntrpMat(lmesh3D%refElem3D%Nnode_v,2,LevelNum))
    allocate( VIntrp_keZ1(LevelNum), VIntrp_keZ2(LevelNum) )
    allocate( VIntrp_coef1(LevelNum) )

    do k=1, LevelNum
      call prep_VIntrp( VIntrpMat(:,:,k), VIntrp_keZ1(k), VIntrp_keZ2(k), VIntrp_coef1(k), &
        levels(k), lmesh3D, lmesh3D%refElem3D )
    end do

    !-------------
    do m=1, mesh_num
      target_myrank = target_proc_s + m - 1      
      !-
      if ( allocated(in_files) ) then
        call in_files(m)%Init( varNum, mesh3D=mesh3D_list(m) )
        call in_files(m)%Open( in_fbasename, myrank=target_myrank )
      end if

      if ( allocated(in_bs_files) ) then
        call in_bs_files(m)%Init( 2, mesh3D=mesh3D_list(m) )
        call in_bs_files(m)%Open( in_bs_fbasename, myrank=target_myrank )  
      end if

      if ( allocated(in_files_2d) ) then
        call in_files_2d(m)%Init( var2DNum, mesh2D=mesh3D_list(m)%mesh2D )
        call in_files_2d(m)%Open( in_fbasename, myrank=target_myrank )
      end if
    end do

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
    if ( ( .not. allocated(in_files_2d) ) .and. ( .not. allocated(in_files_2d) ) ) then
      output_dt = 0.0_RP; var_num_step = 0
    end if

    if ( output_dt < EPS .and. var_num_step == 0 ) then
      var_num_step   = 1
    end if

    !--
    do k=s_ks, s_ke
      k_wave_number(k) = k
    end do
    do k=0, s_ke
      oneD_wave_number(k) = k
    end do
    do l=s_ls, s_le
      l_wave_number(l) = l
    end do

    !--
    if (myrank==0) then
      output_flag = .true.
    else
      output_flag = .false.
    end if
    if (output_flag) then
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

      do nn=1, varNum
        call nc_check( nf90_def_var( ncido, trim(vars_list(nn))//"_r", &
          NF90_DOUBLE, (/ k_wave_number_dim_id, l_wave_number_dim_id, level_dim_id, time_dim_id /), var_id_list_r(nn) ) )
        call nc_check( nf90_def_var( ncido, trim(vars_list(nn))//"_i", &
          NF90_DOUBLE, (/ k_wave_number_dim_id, l_wave_number_dim_id, level_dim_id, time_dim_id /), var_id_list_i(nn) ) )
      end do
      do nn=1, var2DNum  
        call nc_check( nf90_def_var( ncido, trim(vars2D_list(nn))//"_r", &
          NF90_DOUBLE, (/ k_wave_number_dim_id, l_wave_number_dim_id, time_dim_id /), var2D_id_list_r(nn) ) )
        call nc_check( nf90_def_var( ncido, trim(vars2D_list(nn))//"_i", &
          NF90_DOUBLE, (/ k_wave_number_dim_id, l_wave_number_dim_id, time_dim_id /), var2D_id_list_i(nn) ) )
      end do

      if ( KinEnergyAnalysisFlag ) then
        call nc_check( nf90_def_var( ncido, "kin_energy_spectra_h", &
          NF90_DOUBLE, (/ k_wave_number_dim_id, l_wave_number_dim_id, level_dim_id, time_dim_id /), kin_energy_spectra_h_id ) )
        call nc_check( nf90_def_var( ncido, "kin_energy_spectra_v", &
          NF90_DOUBLE, (/ k_wave_number_dim_id, l_wave_number_dim_id, level_dim_id, time_dim_id /), kin_energy_spectra_v_id ) )
        call nc_check( nf90_def_var( ncido, "kin_energy_spectra_1d_h", &
          NF90_DOUBLE, (/ oneD_wave_number_dim_id, level_dim_id, time_dim_id /), kin_energy_spectra_1d_h_id ) )
        call nc_check( nf90_def_var( ncido, "kin_energy_spectra_1d_v", &
          NF90_DOUBLE, (/ oneD_wave_number_dim_id, level_dim_id, time_dim_id /), kin_energy_spectra_1d_v_id ) )
      end if
      call nc_check( nf90_enddef(ncido) )
      
      call nc_check( nf90_put_var(ncido, k_wave_number_id, k_wave_number ) )
      call nc_check( nf90_put_var(ncido, l_wave_number_id, l_wave_number ) )
      call nc_check( nf90_put_var(ncido, oneD_wave_number_id, oneD_wave_number ) )
      if (varNum > 0) call nc_check( nf90_put_var(ncido, level_id, levels ) )
    end if

    return
  end subroutine vars_init

!OCL SERIAL
  subroutine vars_final()
    implicit none

    integer :: m
    !----------------------------

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
      call nc_check( nf90_close(ncido) )
    end if

    return
  end subroutine vars_final

  subroutine vars_write( istep, LevelNum )
    implicit none
    integer, intent(in) :: istep
    integer, intent(in) :: LevelNum

    integer :: vid
    real(RP) :: s_var3D_r_buf(s_ks:s_ke,s_ls:s_le,LevelNum)
    real(RP) :: s_var3D_i_buf(s_ks:s_ke,s_ls:s_le,LevelNum)

    real(RP) :: s_kin1D_buf(0:s_ke,LevelNum)

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
      do vid=1, varNum
        do k=1, LevelNum
          s_var3D_r_buf(:,:,k) = s_var3D(:,:,1,k,vid)
          s_var3D_i_buf(:,:,k) = s_var3D(:,:,2,k,vid)
        end do
        call nc_check( nf90_put_var( ncido, var_id_list_r(vid), s_var3D_r_buf, &
          start=[1,1,1,istep], count=[s_kall,s_lall,LevelNum,1]) )
        call nc_check( nf90_put_var( ncido, var_id_list_i(vid), s_var3D_i_buf, &
          start=[1,1,1,istep], count=[s_kall,s_lall,LevelNum,1]) )
      end do
      do vid=1, var2DNum
        call nc_check( nf90_put_var( ncido, var2D_id_list_r(vid), s_var2D(:,:,1,vid), &
          start=[1,1,istep], count=[s_kall,s_lall,1]) )
        call nc_check( nf90_put_var( ncido, var2D_id_list_i(vid), s_var2D(:,:,2,vid), &
          start=[1,1,istep], count=[s_kall,s_lall,1]) )
      end do
      if ( KinEnergyAnalysisFlag ) then
        do k=1, LevelNum
          s_var3D_r_buf(:,:,k) = 0.5_RP * ( s_var3D(:,:,1,k,varNum0+1)**2 + s_var3D(:,:,2,k,varNum0+1)**2 &
                                          + s_var3D(:,:,1,k,varNum0+2)**2 + s_var3D(:,:,2,k,varNum0+2)**2 )
          call calc_oneD_spectra( s_var3D_r_buf(:,:,k), s_kin1D_buf(:,k) )
        end do
        call nc_check( nf90_put_var( ncido, kin_energy_spectra_h_id, s_var3D_r_buf, &
          start=[1,1,1,istep], count=[s_kall,s_lall,LevelNum,1]) )        
        call nc_check( nf90_put_var( ncido, kin_energy_spectra_1d_h_id, s_kin1D_buf, &
          start=[1,1,istep], count=[s_K1Dall,LevelNum,1]) )
        
        do k=1, LevelNum
          s_var3D_r_buf(:,:,k) = 0.5_RP * ( s_var3D(:,:,1,k,varNum0+3)**2 + s_var3D(:,:,2,k,varNum0+3)**2 )
          call calc_oneD_spectra( s_var3D_r_buf(:,:,k), s_kin1D_buf(:,k) )
        end do
        call nc_check( nf90_put_var( ncido, kin_energy_spectra_v_id, s_var3D_r_buf, &
          start=[1,1,1,istep], count=[s_kall,s_lall,LevelNum,1]) )
        call nc_check( nf90_put_var( ncido, kin_energy_spectra_1d_v_id, s_kin1D_buf, &
          start=[1,1,istep], count=[s_K1Dall,LevelNum,1]) )
      end if
    end if
    return
  contains
!OCL SERAIL
    subroutine calc_oneD_spectra( s_2D, s_1D )
      implicit none
      real(RP), intent(in) :: s_2D(s_ks:s_ke,s_ls:s_le)
      real(RP), intent(out) :: s_1D(0:s_ke)
      
      integer :: k1D, kk, ll
      real(RP) :: dk, dist
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

    type(MeshField3D) :: tmp_field3D
    type(MeshField2D) :: tmp_field2D

    integer :: m
    integer :: vid
    character(len=H_SHORT) :: varname
    integer :: ke2D, ke
    integer :: p, pp, ph, pz
    integer :: k

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

        if ( trim(varname) == "U" ) vid_U = vid
        if ( trim(varname) == "V" ) vid_V = vid
        if ( trim(varname) == "W" ) vid_W = vid
      end do

      if ( vid_U > 0 .and. vid_V > 0 .and. vid_W > 0 ) then
        ! allocate( tmp_U(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        ! allocate( tmp_V(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        ! allocate( tmp_Umet(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        ! allocate( tmp_Vmet(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        ! allocate( gam(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
  
        ! gam(:,:) = 1.0_RP
        ! do k = 1, LayerNum
        !   tmp_U(:,:) = g_var3D(vid_U,k,:,:,m)
        !   tmp_V(:,:) = g_var3D(vid_V,k,:,:,m)

        !   call CubedSphereCoordCnv_CS2LonLatVec( lcmesh2D%panelID, lcmesh2D%pos_en(:,:,1), lcmesh2D%pos_en(:,:,2), gam(:,:), &
        !     lcmesh2D%Ne*lcmesh2D%refElem2D%Np, tmp_U, tmp_V, &
        !     tmp_Umet, tmp_Vmet )
          
        !   g_var3D(vid_U,k,:,:,m) = tmp_Umet(:,:)
        !   g_var3D(vid_V,k,:,:,m) = tmp_Vmet(:,:)
        ! end do
        
        if ( KinEnergyAnalysisFlag ) then
          do k=1, LayerNum
            call g_tmp_ddens(k)%Init("DDENS", "", mesh3D_list(m)%mesh2D )
          end do
          allocate( RsqrtDENS(lcmesh2D%refElem2D%Np,lcmesh2D%Ne,LayerNum) )

          call get_gvar3D( istep, in_files(m), "DDENS", mesh3D_list(m), levels, &
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
              g_var3D_list(k,varNum0+1,m)%local(1)%val(ph,ke2D) = RsqrtDENS(ph,ke2D,k) * g_var3D_list(k,vid_U,m)%local(1)%val(ph,ke2D)
              g_var3D_list(k,varNum0+2,m)%local(1)%val(ph,ke2D) = RsqrtDENS(ph,ke2D,k) * g_var3D_list(k,vid_V,m)%local(1)%val(ph,ke2D)
              g_var3D_list(k,varNum0+3,m)%local(1)%val(ph,ke2D) = RsqrtDENS(ph,ke2D,k) * g_var3D_list(k,vid_W,m)%local(1)%val(ph,ke2D)
            end do
            end do
          end do
          !$omp end parallel

          do k=1, LayerNum
            call g_tmp_ddens(k)%Final()
          end do
          deallocate( RsqrtDENS )
        end if

        ! deallocate( tmp_U, tmp_V, tmp_Umet, tmp_Vmet, gam )
      end if
    end do
    return
  end subroutine vars_read

!-------
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

    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D

    type(MeshField3D) :: tmp_field3D

    integer :: LayerNum

    integer :: ke2D, ke, ke2
    integer :: ke_z1, ke_z2
    integer :: p, pp, ph, pz
    integer :: k
    real(RP), allocatable :: tmp3D_zxy(:,:,:)
    real(RP) :: coef1, coef2
    !---------------------------------------------------

    LayerNum = size(levels)

    call tmp_field3D%Init( varname, "", mesh3D )
    call in_file%Read_Var( &
      MF3D_XYZT, varname, tmp_field3D, step=istep )
    
    lcmesh3D => mesh3D%lcmesh_list(1)
    elem3D => lcmesh3D%refElem3D

    allocate( tmp3D_zxy(elem3D%Nnode_v,elem3D%Nnode_h1D**2,2) )

    do k = 1, LayerNum
      if ( ZintrpType_ID == ZintrpType_Nearest_ID ) then
        do ke = lcmesh3D%NeS, lcmesh3D%NeE
            ke2D = lcmesh3D%EMap3Dto2D(ke)
            do pz=1, elem3D%Nnode_v
            do ph=1, elem3D%Nnode_h1D**2
              p = ph + (pz-1)*elem3D%Nnode_h1D**2
              pp = ph + min(pz,elem3D%Nnode_v-1)*elem3D%Nnode_h1D**2
              if ( lcmesh3D%pos_en(p,ke,3) <= levels(k)   &
                .and. lcmesh3D%pos_en(pp,ke,3) >= levels(k) ) then
    !            LOG_INFO("get_gvar3D",*) trim(varname), ": level=", levels(k), "k,pz=", k, pz, "val=", tmp_field3D%local(1)%val(p,ke)
                gvar_out(k)%local(1)%val(ph,ke2D) = tmp_field3D%local(1)%val(p,ke)
              end if
            end do
            end do
          end do
      else if ( ZintrpType_ID == ZintrpType_SamplingUniPt_LinInterp_ID ) then
        coef1 = VIntrp_coef1(k); coef2 = 1.0_RP - coef1
        !$omp parallel do private(ke2D, ke_z1, ke_z2, ke, ke2, ph, tmp3D_zxy)
        do ke2D=1, lcmesh3D%Ne2D
          ke_z1 = VIntrp_keZ1(k); ke_z2 = VIntrp_keZ2(k)
          ke = ke2D + (ke_z1-1)*lcmesh3D%Ne2D
          ke2 = ke2D + (ke_z2-1)*lcmesh3D%Ne2D
          do ph=1, elem3D%Nnode_h1D**2
            tmp3D_zxy(:,ph,1) = tmp_field3D%local(1)%val(elem3D%Colmask(:,ph),ke)
            tmp3D_zxy(:,ph,2) = tmp_field3D%local(1)%val(elem3D%Colmask(:,ph),ke2)
          end do
          do ph=1, elem3D%Nnode_h1D**2
            gvar_out(k)%local(1)%val(ph,ke2D) = &
                coef1 * sum( VIntrpMat(:,1,k) * tmp3D_zxy(:,ph,1) ) &
              + coef2 * sum( VIntrpMat(:,2,k) * tmp3D_zxy(:,ph,2) ) 
          end do
        end do
      end if
    end do

    call tmp_field3D%Final()

    return
  end subroutine get_gvar3D

!OCL SERIAL
  subroutine prep_VIntrp( IntrpMat, VIntrp_k1, VIntrp_k2, VIntrp_coef1, target_lev, lmesh3D, elem3D )
    use scale_polynominal, only: &
      Polynominal_GenLagrangePoly
    implicit none
    class(ElementBase3D), intent(in) :: elem3D
    class(LocalMesh3D), intent(in) :: lmesh3D
    real(RP), intent(out) :: IntrpMat(elem3D%Nnode_v,2)
    integer, intent(out) :: VIntrp_k1, VIntrp_k2
    real(RP), intent(out) :: VIntrp_coef1
    real(RP), intent(in) :: target_lev

    integer :: ke, ke_z, ke_z2
    integer :: pz, pz1, pz2
    integer :: k, kk
    real(RP) :: vz(8)
    real(RP) :: delz
    real(RP) :: lev_uniform(elem3D%Nnode_v*lmesh3D%NeZ)
    real(RP) :: lev_uniform_xi(elem3D%Nnode_v*lmesh3D%NeZ)

    real(RP) :: lag_poly(2,elem3D%Nnode_v)
    !---------------------------------

    !$omp parallel do private(ke,pz,k,vz,delz)
    do ke_z=1, lmesh3D%NeZ
      ke = 1 + (ke_z-1)*lmesh3D%Ne2D
      vz(:) = lmesh3D%pos_ev(lmesh3D%EToV(ke,:),3)
      delz = ( vz(5) - vz(1) ) / real(elem3D%Nnode_v, kind=RP)
      do pz=1, elem3D%Nnode_v
        k = pz + (ke_z-1)*elem3D%Nnode_v
        lev_uniform(k) = vz(1) + (pz - 0.5_RP)*delz
        lev_uniform_xi(k) = - 1.0_RP + 2.0_RP * ( lev_uniform(k) - vz(1) ) / ( vz(5) - vz(1) )
      end do
    end do
    VIntrp_k1 = -1
    VIntrp_k2 = -1
    do ke_z=1, lmesh3D%NeZ
    do pz=1, elem3D%Nnode_v
      k = pz+(ke_z-1)*elem3D%Nnode_v
      kk = min(k+1, size(lev_uniform))
      if ( lev_uniform(k) <= target_lev .and. target_lev <= lev_uniform(kk) ) then
        if ( pz+1 > elem3D%Nnode_v ) then
          ke_z2 = ke_z + 1; pz2 = 1
        else
          ke_z2 = ke_z; pz2 = pz + 1
        end if
        VIntrp_k1 = ke_z; VIntrp_k2 = ke_z2; pz1 = pz
        exit
      end if
    end do
    end do

    if ( VIntrp_k1 < 0 .or. VIntrp_k2 < 0 ) then
      LOG_INFO("prep_VIntrp",*) "The specified level is invalid. Check!", target_lev
      call PRC_abort
    end if

    !-
    k = pz1 + (VIntrp_k1-1)*elem3D%Nnode_v
    kk = pz2 + (VIntrp_k2-1)*elem3D%Nnode_v

    lag_poly(:,:) = Polynominal_GenLagrangePoly( elem3D%PolyOrder_v, elem3D%x3(elem3D%Colmask(:,1)), (/ lev_uniform_xi(k), lev_uniform_xi(kk) /) )
    IntrpMat(:,:) = transpose(lag_poly)
    VIntrp_coef1 = ( lev_uniform(kk) - target_lev ) / ( lev_uniform(kk) - lev_uniform(k) )

    LOG_INFO("prep_VIntrp",*) "k1,k2=", VIntrp_k1, VIntrp_k2, ": pz1, pz2=", pz1, pz2
    LOG_INFO("prep_VIntrp",*) "lev=", lev_uniform(k), lev_uniform(kk)
    LOG_INFO("prep_VIntrp",*) "lev_xi=", lev_uniform_xi(k), lev_uniform_xi(kk)
    LOG_INFO("prep_VIntrp",*) "VIntrp_coef1=", VIntrp_coef1


    return
  end subroutine prep_VIntrp

!OCL SERIAL
  subroutine nc_check( status )
    integer, intent (in) :: status
    implicit none
    if(status /= nf90_noerr) then 
      LOG_INFO('NetCDF_check',*) trim(nf90_strerror(status))
      call flush(IO_FID_CONF)
      call PRC_abort
    end if

    return
  end subroutine nc_check

end module mod_spectral_analysis_vars