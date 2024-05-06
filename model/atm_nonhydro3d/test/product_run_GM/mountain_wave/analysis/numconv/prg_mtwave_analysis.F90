#include "scalelib.h"
program prg_mtwave_analysis
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort

  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D

  use scale_element_line, only: LineElement
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_mesh_linedom1d, only: MeshLineDom1D

  use scale_localmesh_2d, only: LocalMesh2D

  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D

  use scale_meshutil_vcoord, only: MESH_VCOORD_TERRAIN_FOLLOWING_ID
  use scale_mesh_topography, only: MeshTopography

  use scale_sparsemat, only: SparseMat
  
  use scale_meshfield_base, only: &
    MeshField1D, MeshField3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldContainer
    use scale_meshfieldcomm_rectdom2d, only: &
    MeshFieldCommRectDom2D
  use scale_meshfieldcomm_cubedom3d, only: &
    MeshFieldCommCubeDom3D
  use scale_meshfieldcomm_cubedspheredom2d, only: &
    MeshFieldCommCubedSphereDom2D
  use scale_meshfieldcomm_cubedspheredom3d, only: &
    MeshFieldCommCubedSphereDom3D
  
  use scale_file_base_meshfield, only: &
    FILE_base_meshfield

  use scale_mesh_base3d, only: &
    DIMTYPEID_XYZT => MeshBase3D_DIMTYPEID_XYZT
  
  use scale_meshfield_analysis_numerror, only: MeshFieldAnalysisNumerror3D

  implicit none

   ! MPI parameters
  integer                 :: nprocs                      ! number of processes               (execution)
  integer                 :: myrank                      ! my rank                           (execution)
  logical                 :: ismaster                    ! master process?                   (execution)

  type(HexahedralElement) :: refElem3D
  type(MeshCubedSphereDom3D), target :: mesh3D
  type(MeshTopography), target :: topography

  type(LineElement) :: refElemV1D
  type(MeshLineDom1D) :: meshV1D

  integer :: istep
  integer :: istep_ref_offset
  integer :: istep_ref

  type(MeshField3D) :: DDENS_ref
  type(MeshField3D) :: U
  type(MeshField3D) :: U_ref
  type(MeshField3D) :: V
  type(MeshField3D) :: V_ref
  type(MeshField3D) :: W
  type(MeshField3D) :: W_ref
  type(MeshField3D) :: DRHOT_ref
  type(MeshField3D), target :: PT
  type(MeshField3D), target :: PRES

  type(MeshField3D), target :: DDENS
  type(MeshField3D), target :: MOMX
  type(MeshField3D), target :: MOMY
  type(MeshField3D), target :: MOMZ
  type(MeshField3D), target :: DRHOT
  type(MeshField3D), target :: PRES_hyd
  type(MeshField3D), target :: DENS_hyd
  type(MeshField3D), target :: PRES_hyd_ref
  type(MeshField3D), target :: DENS_hyd_ref

  type(MeshFieldCommCubedSphereDom3D) :: vars_comm
  type(MeshFieldContainer) :: comm_vars_ptr(9)

  type(FILE_base_meshfield) :: in_file
  type(FILE_base_meshfield) :: in_file_ref
  type(FILE_base_meshfield) :: in_bs_file
  type(FILE_base_meshfield) :: in_bs_file_ref

  integer :: num_step       = 1
  real(RP) :: start_time0   = 12600.0_RP 
  real(RP) :: output_tintrv = 60.0_RP

  integer :: n
  class(LocalMesh3D), pointer :: lcmesh3D
  real(RP) :: harea

  type(FILE_base_meshfield) :: out_file_V1D
  character(*), parameter :: dtype = 'REAL8'

  real(RP) :: start_time, end_time

  type(MeshFieldAnalysisNumerror3D) :: numerror_analysis
  integer :: numerror_vid(5)

  logical :: apply_window = .false.
  real(RP) :: window_xs = 0.0_RP
  real(RP) :: window_xe = 144.E3_RP
  real(RP) :: window_ys = 0.0_RP
  real(RP) :: window_ye = 144.E3_RP
  real(RP) :: window_zs = 0.0_RP
  real(RP) :: window_ze = 12.0E3_RP

  !---------------------------------------

  call initialize()
  !-
  do istep=1, num_step
    LOG_PROGRESS(*) "istep=", istep
    start_time = start_time0 + (istep-2)*output_tintrv
    end_time   = start_time0 + (istep-1)*output_tintrv
    istep_ref = istep_ref_offset + istep

    if (istep==1) then
      ! call in_bs_file%Read_Var(DIMTYPEID_XYZT, "DENS_hyd", DENS_hyd, step=1)
      ! call in_bs_file%Read_Var(DIMTYPEID_XYZT, "PRES_hyd", PRES_hyd, step=1)
      ! call in_bs_file_ref%Read_Var(DIMTYPEID_XYZT, "DENS_hyd", DENS_hyd_ref, step=1)
      ! call in_bs_file_ref%Read_Var(DIMTYPEID_XYZT, "PRES_hyd", PRES_hyd_ref, step=1)
    end if
    call in_file%Read_Var(DIMTYPEID_XYZT, "DDENS", DDENS, step=istep)
    call in_file_ref%Read_Var(DIMTYPEID_XYZT, "DDENS", DDENS_ref, step=istep_ref)
    call in_file%Read_Var(DIMTYPEID_XYZT, "U", U, step=istep)
    call in_file_ref%Read_Var(DIMTYPEID_XYZT, "U", U_ref, step=istep_ref)
    call in_file%Read_Var(DIMTYPEID_XYZT, "V", V, step=istep)
    call in_file_ref%Read_Var(DIMTYPEID_XYZT, "V", V_ref, step=istep_ref)
    call in_file%Read_Var(DIMTYPEID_XYZT, "W", W, step=istep)
    call in_file_ref%Read_Var(DIMTYPEID_XYZT, "W", W_ref, step=istep_ref)
    call in_file%Read_Var(DIMTYPEID_XYZT, "THERM", DRHOT, step=istep)
    call in_file_ref%Read_Var(DIMTYPEID_XYZT, "THERM", DRHOT_ref, step=istep_ref)
!    call in_file%Read_Var(DIMTYPEID_XYZT, "PT", PT, step=istep)

    call numerror_analysis%Evaluate( &
      istep, end_time, mesh3D, evaluate_error_lc )

    ! do n=1, mesh3D%LOCAL_MESH_NUM
    !   lcmesh3D => mesh3D%lcmesh_list(n)
    !   call set_prgvars( &
    !     MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,         & ! (out)
    !     PRES%local(n)%val,                                                                   & ! (out)
    !     DDENS%local(n)%val, U%local(n)%val, V%local(n)%val, W%local(n)%val, PT%local(n)%val, & ! (in)
    !     DENS_hyd%local(n)%val, PRES_hyd%local(n)%val, lcmesh3D, lcmesh3D%refElem3D           ) ! (in)
    ! end do
    ! call vars_comm%Put( comm_vars_ptr, 1 )
    ! call vars_comm%Exchange()
    ! call vars_comm%Get( comm_vars_ptr, 1 )
    
    ! do n=1, mesh3D%LOCAL_MESH_NUM
    !   lcmesh3D => mesh3D%lcmesh_list(n)
    !   call analyze_lc_1( &
    !     DENS_V1D%local(1)%val, RHOT_V1D%local(1)%val, MOMZ_V1D%local(1)%val,    &
    !     DDENS%local(n)%val, PT%local(n)%val, W%local(n)%val,                    &
    !     DENS_hyd%local(n)%val,                                                  &
    !     n, lcmesh3D, refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D, &
    !     meshV1D%lcmesh_list(1), refElemV1D )
    ! end do

    ! call global_horizontal_mean( DENS_V1D )
    ! call global_horizontal_mean( RHOT_V1D )
    ! call global_horizontal_mean( MOMZ_V1D )

    ! do n=1, mesh3D%LOCAL_MESH_NUM
    !   lcmesh3D => mesh3D%lcmesh_list(n)
    !   call analyze_lc_2( &
    !     PT_V1D%local(1)%val, W_V1D%local(1)%val,                                           &
    !     W_PRIM3_V1D%local(1)%val,                                                          &
    !     HEAT_TOTFLX_V1D%local(1)%val, HEAT_MEANFLX_V1D%local(1)%val, HEAT_EDDYFLX_V1D%local(1)%val, &
    !     MOMZ_EDDYFLX_V1D%local(1)%val,                      &
    !     DDENS%local(n)%val, PT%local(n)%val, W%local(n)%val,                               &
    !     DENS_V1D%local(1)%val, RHOT_V1D%local(1)%val, MOMZ_V1D%local(1)%val,               &
    !     DENS_hyd%local(n)%val,                                                             &
    !     n, lcmesh3D, refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D,            &
    !     meshV1D%lcmesh_list(1), refElemV1D )
    ! end do

    ! call global_horizontal_mean( W_PRIM3_V1D )
    ! call global_horizontal_mean( HEAT_TOTFLX_V1D )
    ! call global_horizontal_mean( HEAT_MEANFLX_V1D )
    ! call global_horizontal_mean( HEAT_EDDYFLX_V1D )
    ! call global_horizontal_mean( MOMZ_EDDYFLX_V1D )

    ! call diag_tb_process( istep, start_time, end_time, &
    !   DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    !   DENS_hyd, PRES_hyd, PRES, PT, DENS_V1D,          &
    !   mesh3D, meshV1D                                  )

    ! if (ismaster) then
    !   call out_file_V1D%Write_var1D( DIAGVID_DENS, DENS_V1D, start_time, end_time )
    !   call out_file_V1D%Write_var1D( DIAGVID_RHOT, RHOT_V1D, start_time, end_time )
    !   call out_file_V1D%Write_var1D( DIAGVID_PT, PT_V1D, start_time, end_time )
    !   call out_file_V1D%Write_var1D( DIAGVID_W , W_V1D , start_time, end_time )
    !   call out_file_V1D%Write_var1D( DIAGVID_W_PRIM3 , W_PRIM3_V1D , start_time, end_time )
    !   call out_file_V1D%Write_var1D( DIAGVID_HEAT_TOTFLX, HEAT_TOTFLX_V1D, start_time, end_time )
    !   call out_file_V1D%Write_var1D( DIAGVID_HEAT_MEANFLX, HEAT_MEANFLX_V1D, start_time, end_time )
    !   call out_file_V1D%Write_var1D( DIAGVID_HEAT_EDDYFLX, HEAT_EDDYFLX_V1D, start_time, end_time )
    !   call out_file_V1D%Write_var1D( DIAGVID_MOMZ_EDDYFLX, MOMZ_EDDYFLX_V1D, start_time, end_time )
    ! end if
  end do
  !-
  call Finalize()

contains
!OCL SERIAL
    subroutine evaluate_error_lc( this, q, qexact, qexact_intrp, lcmesh, elem3D, intrp_epos )
      use scale_localmeshfield_base, only: LocalMeshFieldBase
      implicit none
      class(MeshFieldAnalysisNumerror3D), intent(in) :: this
      class(LocalMesh3D), intent(in) :: lcmesh
      class(ElementBase3D) :: elem3D
      real(RP), intent(out) :: q(elem3D%Np,lcmesh%Ne,this%var_num)
      real(RP), intent(out) :: qexact(elem3D%Np,lcmesh%Ne,this%var_num)
      real(RP), intent(out) :: qexact_intrp(this%intrp_np,lcmesh%Ne,this%var_num)
      real(RP), intent(in) :: intrp_epos(this%intrp_np,3)

      integer :: n
      integer :: ke, ke2D
      integer :: ke_x, ke_y, ke_z
      real(RP) :: DENS(elem3D%Np)
      real(RP) :: lon(elem3D%Np), lat(elem3D%Np)
      integer :: vid
      !---------------------------------------------

      n = lcmesh%lcdomID
        
      !$omp parallel private(vid, ke, ke2D, lon, lat, ke_x, ke_y, ke_z)
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE

        vid = numerror_vid(1)
        q(:,ke,vid) = DDENS%local(n)%val(:,ke)
        qexact(:,ke,vid) = DDENS_ref%local(n)%val(:,ke)

        vid = numerror_vid(2)
        q(:,ke,vid) = U%local(n)%val(:,ke)
        qexact(:,ke,vid) = U_ref%local(n)%val(:,ke)

        vid = numerror_vid(3)
        q(:,ke,vid) = V%local(n)%val(:,ke)
        qexact(:,ke,vid) = V_ref%local(n)%val(:,ke)        

        vid = numerror_vid(4)
        q(:,ke,vid) = W%local(n)%val(:,ke)
        qexact(:,ke,vid) = W_ref%local(n)%val(:,ke) 
        
        vid = numerror_vid(5)
        q(:,ke,vid) = DRHOT%local(n)%val(:,ke)
        qexact(:,ke,vid) = DRHOT_ref%local(n)%val(:,ke)        
      end do
      !$omp end do

      if ( apply_window ) then
        !$omp do collapse(2)
        do vid=1, size(numerror_vid)
          do ke_z=1, lcmesh%NeZ
          do ke_y=1, lcmesh%NeY
          do ke_x=1, lcmesh%NeX
            ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1)*lcmesh%NeX*lcmesh%NeY
            ke2D = lcmesh%EMap3Dto2D(ke)
            lon(:) = lcmesh%lon2D(elem3D%IndexH2Dto3D,ke2D)
            lat(:) = lcmesh%lat2D(elem3D%IndexH2Dto3D,ke2D)

            where ( &
              lon(:) < window_xs .or. lon(:) > window_xe .or. &
              lat(:) < window_ys .or. lat(:) > window_ye .or. &
              lcmesh%pos_en(:,ke,3) < window_zs .or. lcmesh%pos_en(:,ke,3) > window_ze )
              q(:,ke,vid) = 0.0_RP
              qexact(:,ke,vid) = 0.0_RP
            end where
          end do
          end do
          end do
        end do
        !$omp end do
      end if

      !$omp do collapse(2)
      do vid=1, this%var_num
      do ke=lcmesh%NeS, lcmesh%NeE
        qexact_intrp(:,ke,vid) = matmul( this%IntrpMat, qexact(:,ke,vid) )
      end do
      end do
      !$omp end do      

      !$omp end parallel

      return
    end subroutine evaluate_error_lc

!OCL SERIAL
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

    integer :: ierr
    integer :: comm                     ! communicator   (execution)
    logical :: fileexist
    character(len=H_LONG) :: cnf_fname  ! config file for launcher

    character(len=H_LONG) :: in_filebase
    character(len=H_LONG) :: in_bs_filebase
    character(len=H_LONG) :: in_ref_filebase
    character(len=H_LONG) :: in_bs_ref_filebase
    character(len=H_LONG) :: in_topo_ref_filebase
  
    character(len=H_LONG) :: out_filebase_V1D

    integer :: PolyOrderErrorCheck
    character(len=H_MID) :: NUMERROR_LOG_OUT_BASENAME = 'LOG_NUMERROR'

    !-
    namelist / PARAM_MTWAVE_ANALYSIS / &
      in_filebase,                    &
      in_bs_filebase,                 &
      in_ref_filebase,                &
      in_bs_ref_filebase,             &
      in_topo_ref_filebase,           &
      out_filebase_V1D,               &
      start_time0,                    &
      output_tintrv,                  &
      num_step,                       &
      istep_ref_offset,               &
      PolyOrderErrorCheck,            &
      NUMERROR_LOG_OUT_BASENAME,      &
      apply_window,                   &
      window_xs, window_xe, &
      window_ys, window_ye, &
      window_zs, window_ze

    !-
    integer :: k
    logical :: is_spec_FZ    
    integer, parameter :: FZ_nmax = 1000
    real(RP) :: FZ(FZ_nmax)
  
    integer :: Nprc = 6
    integer :: NeGX = 2
    integer :: NeGY = 2
    integer :: NeZ  = 12
    integer, parameter :: NLocalMeshPerPrc = 1
    
    integer :: PolyOrder_h = 7
    integer :: PolyOrder_v = 7
    logical  :: LumpedMassMatFlag = .false.

    real(RP) :: dom_zmin = 0.0_RP
    real(RP) :: dom_zmax = 3.0E3_RP    

    namelist / PARAM_ATMOS_MESH / &
      dom_zmin, dom_zmax, FZ,                      &
      NeGX, NeGY, NeZ,                             &
      PolyOrder_h, PolyOrder_v, LumpedMassMatFlag, &
      Nprc
        
    type(FILE_base_meshfield) :: file_topo

    type(MeshFieldCommCubedSphereDom3D) :: comm3D
    type(MeshFieldCommCubedSphereDom2D) :: comm2D
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
    call IO_setup( "MTWAVE_ANALYSIS", cnf_fname )

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
    ! setup constants
    call CONST_setup
    ! setup fie I/O
    call FILE_setup( myrank )
    
    LOG_NEWLINE
    LOG_INFO("MTWAVE_ANALYSIS",*) 'Setup'

    !-
    in_filebase = "./history"
    out_filebase_V1D = "./analysis/history_v1D"
    PolyOrderErrorCheck = PolyOrder_h
    FZ(:) = - 1.0_RP
    istep_ref_offset = 0

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MTWAVE_ANALYSIS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("MTWAVE_ANALYSIS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("MTWAVE_ANALYSIS_setup",*) 'Not appropriate names in namelist PARAM_MTWAVE_ANALYSIS. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_MTWAVE_ANALYSIS)

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_MESH,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("MTWAVE_ANALYSIS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("MTWAVE_ANALYSIS_setup",*) 'Not appropriate names in namelist PARAM_ATM_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_MESH)

    !- Mesh
    LOG_INFO("MTWAVE_ANALYSIS",*) 'Setup mesh..'

    call refElem3D%Init( PolyOrder_h, PolyOrder_v, LumpedMassMatFlag )
    call refElemV1D%Init( PolyOrder_v, LumpedMassMatFlag )

    is_spec_FZ = .true.
    do k=1, NeZ+1
      if (FZ(k) < 0.0_RP) then
        is_spec_FZ = .false.
      end if
    end do

    if (is_spec_FZ) then
      call mesh3D%Init( NeGX, NeGY, NeZ, RPlanet, dom_zmin, dom_zmax, refElem3D, NLocalMeshPerPrc, &
        Nprc, FZ=FZ(1:NeZ+1) )
    else
      call mesh3D%Init( NeGX, NeGY, NeZ, RPlanet, dom_zmin, dom_zmax, refElem3D, NLocalMeshPerPrc, &
        Nprc )
    end if
    call mesh3D%Generate()

    call meshV1D%Init( NeZ, dom_zmin, dom_zmax, refElemV1D, 1, &
                       nproc=1, myrank=0, FX=FZ(1:NeZ+1)       )
    call meshV1D%Generate()

    ! Topography
    call topography%Init( "topo", mesh3D%mesh2D )

    call file_topo%Init(1, meshCubedSphere2D=mesh3D%mesh2D )
    call file_topo%Open( in_topo_ref_filebase, myrank=PRC_myrank )
    call file_topo%Read_Var( MFTYPE2D_XY, "topo", topography%topo )
    call file_topo%Close()
    call file_topo%Final()

    call comm2D%Init( 1, 0, 0, mesh3D%mesh2D )
    call comm3D%Init( 2, 1, 0, mesh3D )

    call topography%SetVCoordinate( mesh3D,   &
      MESH_VCOORD_TERRAIN_FOLLOWING_ID, mesh3D%zmax_gl, comm3D, comm2D )

    call comm2D%Final()
    call comm3D%Final()

    ! Data
    LOG_INFO("MTWAVE_ANALYSIS",*) 'Setup variables..'

    call DDENS_ref%Init("DDENS", "kg/m3", mesh3D)
    call U%Init("U", "m/s", mesh3D)
    call U_ref%Init("U", "m/s", mesh3D)
    call V%Init("V", "m/s", mesh3D)
    call V_ref%Init("V", "m/s", mesh3D)
    call W%Init("W", "m/s", mesh3D)
    call W_ref%Init("W", "m/s", mesh3D)
    call DRHOT_ref%Init("DRHOT_ref", "kg/m3.K", mesh3D)
    call PT%Init("PT", "K", mesh3D)
    call PRES%Init("PRES", "Pa", mesh3D)

    call DDENS%Init("DDENS", "kg/m3", mesh3D)
    call MOMX%Init("MOMX", "kg/m3.m/s", mesh3D)
    call MOMY%Init("MOMY", "kg/m3.m/s", mesh3D)
    call MOMZ%Init("MOMZ", "kg/m3.m/s", mesh3D)
    call DRHOT%Init("DRHOT", "kg/m3.K", mesh3D)

    call PRES_hyd%Init("PRES_hyd", "Pa", mesh3D)
    call DENS_hyd%Init("DENS_hyd", "kg/m3", mesh3D)

    call PRES_hyd_ref%Init("PRES_hyd", "Pa", mesh3D)
    call DENS_hyd_ref%Init("DENS_hyd", "kg/m3", mesh3D)

    ! call DENS_V1D%Init("DENS", "kg/m3", meshV1D)
    ! call RHOT_V1D%Init("RHOT", "kg/m3.K", meshV1D)
    ! call MOMZ_V1D%Init("MOMZ", "kg/m3.m/s", meshV1D)
    ! call PT_V1D%Init("PT", "K", meshV1D)
    ! call W_V1D %Init("W", "m/s", meshV1D)
    ! call W_PRIM3_V1D%Init("W_PRIM3", "m/s", meshV1D)
    ! call HEAT_TOTFLX_V1D%Init("HEAT_TOTFLX", "J.m-2.s-1", meshV1D)
    ! call HEAT_MEANFLX_V1D%Init("HEAT_MEANFLX", "J.m-2.s-1", meshV1D)
    ! call HEAT_EDDYFLX_V1D%Init("HEAT_EDDYFLX", "J.m-2.s-1", meshV1D)
    ! call MOMZ_EDDYFLX_V1D%Init("MOMZ_EDDYFLX", "m/s.kg.m-2.s-1", meshV1D)

    ! Input & Output   

    LOG_INFO("MTWV_ANALYSIS",*) 'Open input file ..'

    call in_file%Init(5, meshCubedSphere3D=mesh3D)
    call in_file%Open(in_filebase)

    call in_file_ref%Init(5, meshCubedSphere3D=mesh3D)
    call in_file_ref%Open(in_ref_filebase)

    ! call in_bs_file%Init(2, meshCubedSphere3D=mesh3D)
    ! call in_bs_file%Open(in_bs_filebase)

    ! call in_bs_file_ref%Init(2, meshCubedSphere3D=mesh3D)
    ! call in_bs_file_ref%Open(in_bs_ref_filebase)

    if (ismaster) then
      LOG_INFO("MTWV_ANALYSIS",*) 'Setup output file for vertical 1D data..'

      ! call out_file_V1D%Init(10, mesh1D=meshV1D)
      ! call out_file_V1D%Create( out_filebase_V1D, 'Mountain wave analysis', &
!        dtype, fileexist, myrank=myrank )
    !   call out_file_V1D%Def_Var(DENS_V1D, "horizontal averaged density", DIAGVID_DENS, DIMTYPEID_ZT, dtype, &
    !     timeinv=output_tintrv )
    !   call out_file_V1D%Def_Var(RHOT_V1D, "horizontal averaged density*potential temperature", DIAGVID_RHOT, DIMTYPEID_ZT, dtype, &
    !     timeinv=output_tintrv )
    !   call out_file_V1D%Def_Var(PT_V1D, "horizontal averaged PT", DIAGVID_PT, DIMTYPEID_ZT, dtype, &
    !     timeinv=output_tintrv )
    !   call out_file_V1D%Def_Var( W_V1D,  "horizontal averaged W",  DIAGVID_W, DIMTYPEID_ZT, dtype, &
    !     timeinv=output_tintrv )
    !   call out_file_V1D%Def_Var( W_PRIM3_V1D,  "horizontal averaged W_PRIM3",  DIAGVID_W_PRIM3, DIMTYPEID_ZT, dtype, &
    !     timeinv=output_tintrv )
    !   call out_file_V1D%Def_Var( HEAT_TOTFLX_V1D, "horizontal averaged heat flux", DIAGVID_HEAT_TOTFLX, DIMTYPEID_ZT, dtype, &
    !     timeinv=output_tintrv )
    !   call out_file_V1D%Def_Var( HEAT_MEANFLX_V1D, "horizontal averaged mean heat flux", DIAGVID_HEAT_MEANFLX, DIMTYPEID_ZT, dtype, &
    !     timeinv=output_tintrv )
    !   call out_file_V1D%Def_Var( HEAT_EDDYFLX_V1D, "horizontal averaged eddy heat flux", DIAGVID_HEAT_EDDYFLX, DIMTYPEID_ZT, dtype, &
    !     timeinv=output_tintrv )
    !   call out_file_V1D%Def_Var( MOMZ_EDDYFLX_V1D, "horizontal averaged eddy momentum flux in z-direction", &
    !     DIAGVID_MOMZ_EDDYFLX, DIMTYPEID_ZT, dtype, timeinv=output_tintrv )      
      ! call out_file_V1D%End_def()
    end if

    call numerror_analysis%Init( &
      PolyOrderErrorCheck, NUMERROR_LOG_OUT_BASENAME, 1, refElem3D )
    
    call numerror_analysis%Regist( "DDENS", "m/s", numerror_vid(1) )
    call numerror_analysis%Regist( "U", "m/s", numerror_vid(2) )
    call numerror_analysis%Regist( "V", "m/s", numerror_vid(3) )
    call numerror_analysis%Regist( "W", "m/s", numerror_vid(4) )
    call numerror_analysis%Regist( "DRHOT", "m/s", numerror_vid(5) )

    ! Communication
    ! call vars_comm%Init( 9, 0, mesh3D )
    ! comm_vars_ptr(1)%field3d => DDENS
    ! comm_vars_ptr(2)%field3d => MOMX
    ! comm_vars_ptr(3)%field3d => MOMY
    ! comm_vars_ptr(4)%field3d => MOMZ
    ! comm_vars_ptr(5)%field3d => DRHOT
    ! comm_vars_ptr(6)%field3d => DENS_hyd
    ! comm_vars_ptr(7)%field3d => PRES_hyd
    ! comm_vars_ptr(8)%field3d => PT
    ! comm_vars_ptr(9)%field3d => PRES

    ! Diagnose
    ! call common_Init( mesh3D )
    ! call diag_tb_Init( mesh3D, meshV1D, &
    !   out_filebase_tb, dtype, output_tintrv, &
    !   myrank, ismaster, NLocalMeshPerPrc )


    LOG_INFO("MTWAVE_ANALYSIS",*) 'Setup has been finished.'

    return
  end subroutine initialize

!OCL SERIAL
  subroutine finalize()
    use scale_prc, only: &
      PRC_mpibarrier, &
      PRC_MPIfinish
    use scale_file, only: &
      FILE_close_all
    implicit none
    !---------------------------------------------

    ! IO
    call in_file%Close()
    call in_file%Final()

    ! call in_bs_file%Close()
    ! call in_bs_file%Final()

    ! Communication
!    call vars_comm%Final()

    ! Analysis
    call numerror_analysis%Final()

    ! Field
    call U%Final()
    call V%Final()
    call W%Final()
    call W_ref%Final()
    call PT%Final()
    call PRES%Final()

    call DDENS%Final()
    call MOMX%Final()
    call MOMY%Final()
    call MOMZ%Final()
    call DRHOT%Final()

    call PRES_hyd%Final()
    call DENS_hyd%Final()
    
    ! Mesh
    call mesh3D%Final()
    call meshV1D%Final()
    ! call topography%Final()
    ! Element
    call refElem3D%Final()
    call refElemV1D%Final()

    ! stop MPI
    call PRC_mpibarrier
    call PRC_MPIfinish 

    return
  end subroutine finalize

end program prg_mtwave_analysis
