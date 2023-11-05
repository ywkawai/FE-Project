#include "scalelib.h"
module mod_vars
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort

  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
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
  real(RP), public, allocatable :: s_var2D(:,:,:,:)
  real(RP), public, allocatable :: g_var3D(:,:,:,:,:)
  real(RP), public, allocatable :: s_var3D(:,:,:,:,:)

  real(RP), public, allocatable :: g_bs_var3D(:,:,:,:,:)
  integer, public, parameter :: BS_PRESHYD_ID = 1
  integer, public, parameter :: BS_DENSHYD_ID = 2

  real(RP), public, allocatable :: g_ke_var3D(:,:,:,:,:)
  integer, public, parameter :: KE_RsqrtUmet_ID = 1
  integer, public, parameter :: KE_RsqrtVmet_ID = 2
  integer, public, parameter :: KE_RsqrtWmet_ID = 3

  integer :: ncido
  integer :: zonal_wave_number_dim_id
  integer :: total_wave_number_dim_id
  integer :: level_dim_id
  integer :: time_dim_id
  integer :: zonal_wave_number_id
  integer :: total_wave_number_id
  integer :: level_id
  integer :: time_id
  integer, allocatable :: var2D_id_list_r(:)
  integer, allocatable :: var2D_id_list_i(:)
  integer, allocatable :: var_id_list_r(:)
  integer, allocatable :: var_id_list_i(:)

  type(FILE_base_meshfield), allocatable :: in_files(:)
  type(FILE_base_meshfield), allocatable :: in_bs_files(:)
  type(FILE_base_meshfield), allocatable :: in_files_2d(:)

  real(DP) :: output_dt
  real(DP) :: start_sec

  logical :: output_flag
  logical :: KinEnergyAnalysisFlag

contains
!OCL SERIAL
  subroutine vars_init( &
    vars2D, vars, KinEnergyAnalysisFlag_,                            &
    levels, level_units, Np2D, Ne2D, Mt, mesh3D_list, target_proc_s, &
    in_fbasename, in_bs_fbasename, out_fbasename, myrank )

    use scale_const, only: EPS => CONST_EPS
    implicit none
    character(len=H_SHORT), intent(in) :: vars2D(:)
    character(len=H_SHORT), intent(in) :: vars(:)
    logical, intent(in) :: KinEnergyAnalysisFlag_
    real(RP), intent(in) :: levels(:)
    character(*), intent(in) :: level_units
    integer, intent(in) :: Np2D
    integer, intent(in) :: Ne2D
    integer, intent(in) :: Mt
    type(MeshCubedSphereDom3D), intent(in), target :: mesh3D_list(:)
    integer, intent(in) :: target_proc_s
    character(len=*), intent(in) :: in_fbasename
    character(len=*), intent(in) :: in_bs_fbasename
    character(len=*), intent(in) :: out_fbasename
    integer, intent(in) :: myrank
    
    integer :: zonal_wave_number(0:Mt)
    integer :: total_wave_number(0:Mt)
    integer :: l, m

    integer :: nn
    integer :: target_myrank

    integer :: LevelNum
    integer :: mesh_num

    real(DP) :: time_endsec
    !------------------------------------------------------------

    KinEnergyAnalysisFlag = KinEnergyAnalysisFlag_

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
      allocate( s_var2D(var2DNum,0:Mt,0:Mt,2) )
      s_var2D(:,:,:,:) = 0.0_RP

      allocate( var2D_id_list_r(varNum) )
      allocate( var2D_id_list_i(varNum) )  
    end if
    if (varNum > 0) then
      allocate( g_var3D(varNum,LevelNum,Np2D,Ne2D,mesh_num) )
      allocate( s_var3D(varNum,LevelNum,0:Mt,0:Mt,2) )
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
    end if
    
    !-------------

    !-------------
    do m=1, mesh_num
      target_myrank = target_proc_s + m - 1      
      !-
      if ( allocated(in_files) ) then
        call in_files(m)%Init( varNum, meshCubedSphere3D=mesh3D_list(m) )
        call in_files(m)%Open( in_fbasename, myrank=target_myrank )
      end if

      if ( allocated(in_bs_files) ) then
        call in_bs_files(m)%Init( 2, meshCubedSphere3D=mesh3D_list(m) )
        call in_bs_files(m)%Open( in_bs_fbasename, myrank=target_myrank )  
      end if

      if ( allocated(in_files_2d) ) then
        call in_files_2d(m)%Init( var2DNum, meshCubedSphere2D=mesh3D_list(m)%mesh2D )
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
    do l=0, Mt
      zonal_wave_number(l) = l
      total_wave_number(l) = l
    end do

    !--
    if (myrank==0) then
      output_flag = .true.
    else
      output_flag = .false.
    end if
    if (output_flag) then
      call nc_check( nf90_create( trim(out_fbasename)//".nc", NF90_clobber, ncido) )
      call nc_check( nf90_def_dim( ncido, 'm', Mt+1, zonal_wave_number_dim_id ) )
      call nc_check( nf90_def_dim( ncido, 'n', Mt+1, total_wave_number_dim_id ) )
      if (varNum > 0)  call nc_check( nf90_def_dim( ncido, 'level', LevelNum, level_dim_id ) )
      call nc_check( nf90_def_dim( ncido, 'time', NF90_UNLIMITED, time_dim_id ) )

      call nc_check( nf90_def_var( ncido, 'm', NF90_INT, zonal_wave_number_dim_id, zonal_wave_number_id ) )
      call nc_check( nf90_def_var( ncido, 'n', NF90_INT, total_wave_number_dim_id, total_wave_number_id ) )
      if (varNum > 0)  call nc_check( nf90_def_var( ncido, 'level', NF90_DOUBLE, level_dim_id, level_id ) )
      call nc_check( nf90_def_var( ncido, 'time', NF90_DOUBLE, time_dim_id, time_id ) )

      call nc_check( nf90_put_att( ncido, zonal_wave_number_id, 'units', '1' ) )
      call nc_check( nf90_put_att( ncido, total_wave_number_id, 'units', '1' ) )
      if (varNum > 0)  call nc_check( nf90_put_att( ncido, level_id, 'units', trim(level_units) ) )
      call nc_check( nf90_put_att( ncido, time_id, 'units', 'sec' ) )

      do nn=1, varNum  
        call nc_check( nf90_def_var( ncido, trim(vars_list(nn))//"_r", &
          NF90_DOUBLE, (/ total_wave_number_dim_id, zonal_wave_number_dim_id, level_dim_id, time_dim_id /), var_id_list_r(nn) ) )
        call nc_check( nf90_def_var( ncido, trim(vars_list(nn))//"_i", &
          NF90_DOUBLE, (/ total_wave_number_dim_id, zonal_wave_number_dim_id, level_dim_id, time_dim_id /), var_id_list_i(nn) ) )
      end do
      do nn=1, var2DNum  
        call nc_check( nf90_def_var( ncido, trim(vars2D_list(nn))//"_r", &
          NF90_DOUBLE, (/ total_wave_number_dim_id, zonal_wave_number_dim_id, time_dim_id /), var2D_id_list_r(nn) ) )
        call nc_check( nf90_def_var( ncido, trim(vars2D_list(nn))//"_i", &
          NF90_DOUBLE, (/ total_wave_number_dim_id, zonal_wave_number_dim_id, time_dim_id /), var2D_id_list_i(nn) ) )
      end do

      call nc_check( nf90_enddef(ncido) )
      
      call nc_check( nf90_put_var(ncido, zonal_wave_number_id, zonal_wave_number ) )
      call nc_check( nf90_put_var(ncido, total_wave_number_id, total_wave_number ) )
      if (varNum > 0) call nc_check( nf90_put_var(ncido, level_id, levels ) )
    end if

    return
  end subroutine vars_init

!OCL SERIAL
  subroutine vars_final()
    implicit none

    integer :: m
    !----------------------------

    if (output_flag) then
      call nc_check( nf90_close(ncido) )
    end if
    
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

    return
  end subroutine vars_final

  subroutine vars_write( istep, Mt, LevelNum )
    implicit none
    integer, intent(in) :: istep
    integer, intent(in) :: Mt
    integer, intent(in) :: LevelNum

    integer :: vid
    real(RP), allocatable :: s_var3D_r_buf(:,:,:)
    real(RP), allocatable :: s_var3D_i_buf(:,:,:)
    real(RP) :: s_var2D_r_buf(Mt+1,Mt+1)
    real(RP) :: s_var2D_i_buf(Mt+1,Mt+1)

    integer :: k
    !----------------------------------

    LOG_INFO("vars_write",*) "istep=", istep !, ":",  start_sec, output_dt, start_sec + istep * output_dt
    if (output_flag) then
      call nc_check( nf90_put_var( ncido, time_id, (/ start_sec + istep * output_dt /), &
        start=[istep], count=[1]) )

      allocate( s_var3D_r_buf(Mt+1,Mt+1,LevelNum) )        
      allocate( s_var3D_i_buf(Mt+1,Mt+1,LevelNum) )
      do vid=1, varNum
        do k=1, LevelNum
          s_var3D_r_buf(1:Mt+1,1:Mt+1,k) = s_var3D(vid,k,0:Mt,0:Mt,1)
          s_var3D_i_buf(1:Mt+1,1:Mt+1,k) = s_var3D(vid,k,0:Mt,0:Mt,2)
        end do
        call nc_check( nf90_put_var( ncido, var_id_list_r(vid), s_var3D_r_buf, &
          start=[1,1,1,istep], count=[Mt+1,Mt+1,LevelNum,1]) )
        call nc_check( nf90_put_var( ncido, var_id_list_i(vid), s_var3D_i_buf, &
          start=[1,1,1,istep], count=[Mt+1,Mt+1,LevelNum,1]) )
      end do

      do vid=1, var2DNum
        s_var2D_r_buf(1:Mt+1,1:Mt+1) = s_var2D(vid,0:Mt,0:Mt,1)
        s_var2D_i_buf(1:Mt+1,1:Mt+1) = s_var2D(vid,0:Mt,0:Mt,2)
        call nc_check( nf90_put_var( ncido, var2D_id_list_r(vid), s_var2D_r_buf, &
          start=[1,1,istep], count=[Mt+1,Mt+1,1]) )
        call nc_check( nf90_put_var( ncido, var2D_id_list_i(vid), s_var2D_i_buf, &
          start=[1,1,istep], count=[Mt+1,Mt+1,1]) )
      end do
    end if

    return
  end subroutine vars_write

!OCL SERIAL
  subroutine vars_read( istep, mesh3D_list, levels )
    use scale_meshfield_base, only: MeshField2D, MeshField3D
    use scale_mesh_base2d, only: &
      MF2D_XYT => MeshBase2D_DIMTYPEID_XYT
    use scale_mesh_base3d, only: &
      MF3D_XYZT => MeshBase3D_DIMTYPEID_XYZT
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatVec
    implicit none
    integer, intent(in) :: istep
    type(MeshCubedSphereDom3D), intent(in), target :: mesh3D_list(:)
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

    real(RP), allocatable :: tmp_U(:,:)
    real(RP), allocatable :: tmp_V(:,:)
    real(RP), allocatable :: tmp_Umet(:,:)
    real(RP), allocatable :: tmp_Vmet(:,:)
    real(RP), allocatable :: gam(:,:)
    integer :: LayerNum
    integer :: vid_U, vid_V, vid_W

    real(RP), allocatable :: g_tmp_ddens(:,:,:,:)
    real(RP), allocatable :: RsqrtDENS(:,:,:)
    !----------------------------

    LayerNum = size(levels)
    vid_U = -1; vid_V = -1; vid_W = -1

    do m=1, size(mesh3D_list)
      lcmesh2D => mesh3D_list(m)%mesh2D%lcmesh_list(1)
      do vid=1, var2DNum
        varname = trim(vars2D_list(vid))
        call tmp_field2D%Init( varname, "", mesh3D_list(m)%mesh2D )
        call in_files_2d(m)%Read_Var( &
          MF2D_XYT, varname, tmp_field2D, step=istep )

        do ke2D = lcmesh2D%NeS, lcmesh2D%NeE
          g_var2D(vid,:,ke2D,m) = tmp_field2D%local(1)%val(:,ke2D)
        end do
        call tmp_field2D%Final()
      end do
    end do

    do m=1, size(mesh3D_list)
      lcmesh3D => mesh3D_list(m)%lcmesh_list(1)
      elem3D => lcmesh3D%refElem3D
      lcmesh2D => mesh3D_list(m)%mesh2D%lcmesh_list(1)

      if ( allocated(in_bs_files) .and. istep == 1 ) then
        call get_gvar3D( istep, BS_PRESHYD_ID, in_bs_files(m), "PRES_hyd", mesh3D_list(m), levels, &
          g_bs_var3D(:,:,:,:,m) )
        call get_gvar3D( istep, BS_DENSHYD_ID, in_bs_files(m), "DENS_hyd", mesh3D_list(m), levels, &
          g_bs_var3D(:,:,:,:,m) )
      end if

      do vid=1, varNum0
        varname = trim(vars_list(vid))
        call get_gvar3D( istep, vid, in_files(m), varname, mesh3D_list(m), levels, &
          g_var3D(:,:,:,:,m) )

        if ( trim(varname) == "U" ) vid_U = vid
        if ( trim(varname) == "V" ) vid_V = vid
        if ( trim(varname) == "W" ) vid_W = vid
      end do

      if ( vid_U > 0 .and. vid_V > 0 .and. vid_W > 0 ) then
        allocate( tmp_U(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        allocate( tmp_V(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        allocate( tmp_Umet(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        allocate( tmp_Vmet(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
        allocate( gam(lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
  
        gam(:,:) = 1.0_RP
        do k = 1, LayerNum
          tmp_U(:,:) = g_var3D(vid_U,k,:,:,m)
          tmp_V(:,:) = g_var3D(vid_V,k,:,:,m)

          call CubedSphereCoordCnv_CS2LonLatVec( lcmesh2D%panelID, lcmesh2D%pos_en(:,:,1), lcmesh2D%pos_en(:,:,2), gam(:,:), &
            lcmesh2D%Ne*lcmesh2D%refElem2D%Np, tmp_U, tmp_V, &
            tmp_Umet, tmp_Vmet )
          
          g_var3D(vid_U,k,:,:,m) = tmp_Umet(:,:)
          g_var3D(vid_V,k,:,:,m) = tmp_Vmet(:,:)
        end do
        
        if ( KinEnergyAnalysisFlag ) then
          allocate( g_tmp_ddens(1,LayerNum,lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )
          allocate( RsqrtDENS(LayerNum,lcmesh2D%refElem2D%Np,lcmesh2D%Ne) )

          call get_gvar3D( istep, 1, in_files(m), "DDENS", mesh3D_list(m), levels, &
            g_tmp_ddens )

          !$omp parallel 
          !$omp do collapse(2)
          do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
          do ph=1, lcmesh2D%refElem2D%Np
            RsqrtDENS(:,ph,ke2D) = sqrt( g_bs_var3D(BS_DENSHYD_ID,:,ph,ke2D,m) + g_tmp_ddens(1,:,ph,ke2D) )
          end do
          end do
          !$omp do collapse(2)
          do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
          do ph=1, lcmesh2D%refElem2D%Np  
            g_var3D(varNum0+1,:,ph,ke2D,m) = RsqrtDENS(:,ph,ke2D) * g_var3D(vid_U,:,ph,ke2D,m)
            g_var3D(varNum0+2,:,ph,ke2D,m) = RsqrtDENS(:,ph,ke2D) * g_var3D(vid_V,:,ph,ke2D,m)
            g_var3D(varNum0+3,:,ph,ke2D,m) = RsqrtDENS(:,ph,ke2D) * g_var3D(vid_W,:,ph,ke2D,m)
          end do
          end do
          !$omp end parallel

          deallocate( g_tmp_ddens, RsqrtDENS )
        end if

        deallocate( tmp_U, tmp_V, tmp_Umet, tmp_Vmet, gam )
      end if

    end do

    return
  end subroutine vars_read

!-------
!OCL SERIAL
  subroutine get_gvar3D( istep, vid, in_file, varname, mesh3D, levels, &
    gvar3D )
    use scale_meshfield_base, only: MeshField3D
    use scale_mesh_base3d, only: &
      MF3D_XYZT => MeshBase3D_DIMTYPEID_XYZT
    implicit none
    integer, intent(in) :: istep
    integer, intent(in) :: vid
    class(FILE_base_meshfield), intent(inout) :: in_file
    character(len=*), intent(in) :: varname
    class(MeshCubedSphereDom3D), intent(in), target :: mesh3D
    real(RP), intent(in) :: levels(:)
    real(RP), intent(inout) :: gvar3D(:,:,:,:)

    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D

    type(MeshField3D) :: tmp_field3D

    integer :: LayerNum

    integer :: ke2D, ke
    integer :: p, pp, ph, pz
    integer :: k
    !---------------------------------------------------

    LayerNum = size(levels)

    call tmp_field3D%Init( varname, "", mesh3D )
    call in_file%Read_Var( &
      MF3D_XYZT, varname, tmp_field3D, step=istep )
    
    lcmesh3D => mesh3D%lcmesh_list(1)
    elem3D => lcmesh3D%refElem3D
    do k = 1, LayerNum
    do ke = lcmesh3D%NeS, lcmesh3D%NeE
        ke2D = lcmesh3D%EMap3Dto2D(ke)
        do pz=1, elem3D%Nnode_v
        do ph=1, elem3D%Nnode_h1D**2
          p = ph + (pz-1)*elem3D%Nnode_h1D**2
          pp = ph + min(pz,elem3D%Nnode_v-1)*elem3D%Nnode_h1D**2
          if ( lcmesh3D%pos_en(p,ke,3) <= levels(k)   &
            .and. lcmesh3D%pos_en(pp,ke,3) >= levels(k) ) then
            gvar3D(vid,k,ph,ke2D) = tmp_field3D%local(1)%val(p,ke)
          end if
        end do
        end do
      end do
    end do

    call tmp_field3D%Final()

    return
  end subroutine get_gvar3D

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

end module mod_vars