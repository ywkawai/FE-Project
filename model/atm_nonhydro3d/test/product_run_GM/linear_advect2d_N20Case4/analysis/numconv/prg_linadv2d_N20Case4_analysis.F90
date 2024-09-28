#include "scalelib.h"
program prg_linadv2d_N20Case4_analysis
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

  type(LineElement) :: refElem1D
  type(MeshLineDom1D) :: mesh1D

  integer :: istep
  integer :: istep_ref_offset
  integer :: istep_ref

  type(MeshField3D) :: PTracer_ref
  type(MeshField3D) :: PTracer

  type(FILE_base_meshfield) :: in_file
  type(FILE_base_meshfield) :: in_file_ref

  integer :: num_step       = 1
  real(RP) :: start_time0   = 0.0_RP 
  real(RP) :: output_tintrv = 43200.0_RP

  integer :: n
  class(LocalMesh3D), pointer :: lcmesh3D
  real(RP) :: harea
  character(*), parameter :: dtype = 'REAL8'

  real(RP) :: start_time, end_time

  type(MeshFieldAnalysisNumerror3D) :: numerror_analysis
  integer :: numerror_vid(1)
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
    call in_file%Read_Var(DIMTYPEID_XYZT, "PTracer", PTracer, step=istep)
    call in_file_ref%Read_Var(DIMTYPEID_XYZT, "PTracer", PTracer_ref, step=istep_ref)

    call numerror_analysis%Evaluate( &
      istep, end_time, mesh3D, evaluate_error_lc )

    if (istep==1) call print_normalized_factor(PTracer_ref)
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
      integer :: ke
      real(RP) :: DENS(elem3D%Np)

      integer :: vid
      !---------------------------------------------

      n = lcmesh%lcdomID

      !$omp parallel private(vid, ke)
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE

        vid = numerror_vid(1)
        q(:,ke,vid) = PTracer%local(n)%val(:,ke)
        qexact(:,ke,vid) = PTracer_ref%local(n)%val(:,ke)
      end do
      !$omp end do

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
    character(len=H_LONG) :: in_ref_filebase
  

    integer :: PolyOrderErrorCheck
    character(len=H_MID) :: NUMERROR_LOG_OUT_BASENAME = 'LOG_NUMERROR'

    !-
    namelist / PARAM_LINADV_ANALYSIS / &
      in_filebase,                    &
      in_ref_filebase,                &
      start_time0,                    &
      output_tintrv,                  &
      num_step,                       &
      istep_ref_offset,               &
      PolyOrderErrorCheck,            &
      NUMERROR_LOG_OUT_BASENAME
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
    call IO_setup( "LINADV_ANALYSIS", cnf_fname )

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
    ! setup constants
    call CONST_setup
    ! setup fie I/O
    call FILE_setup( myrank )
    
    LOG_NEWLINE
    LOG_INFO("LINADV_ANALYSIS",*) 'Setup'

    !-
    in_filebase = "./history"
    PolyOrderErrorCheck = PolyOrder_h
    FZ(:) = - 1.0_RP
    istep_ref_offset = 0

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LINADV_ANALYSIS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("LINADV_ANALYSIS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("LINADV_ANALYSIS_setup",*) 'Not appropriate names in namelist PARAM_LINADV_ANALYSIS. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_LINADV_ANALYSIS)

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_MESH,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("LINADV_ANALYSIS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("LINADV_ANALYSIS_setup",*) 'Not appropriate names in namelist PARAM_ATM_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_MESH)

    !- Mesh
    LOG_INFO("LINADV_ANALYSIS",*) 'Setup mesh..'

    call refElem3D%Init( PolyOrder_h, PolyOrder_v, LumpedMassMatFlag )

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

    ! 1D mesh
    call refElem1D%Init( PolyOrder_v, .false. )
    call mesh1D%Init( NeZ, dom_zmin, dom_zmax, refElem1D, NLocalMeshPerPrc, &
        1, FX=mesh3D%FZ )
    
    ! Data
    LOG_INFO("LINADV_ANALYSIS",*) 'Setup variables..'

    call PTracer_ref%Init("PTracer", "1", mesh3D)
    call PTracer%Init("PTracer", "1", mesh3D)

    ! Input & Output   

    LOG_INFO("LINADV_ANALYSIS",*) 'Open input file ..'

    call in_file%Init(1, meshCubedSphere3D=mesh3D)
    call in_file%Open(in_filebase)

    call in_file_ref%Init(1, meshCubedSphere3D=mesh3D)
    call in_file_ref%Open(in_ref_filebase)

    call numerror_analysis%Init( &
      PolyOrderErrorCheck, NUMERROR_LOG_OUT_BASENAME, 1, refElem3D )
    
    call numerror_analysis%Regist( "PTracer", "1", numerror_vid(1) )

    LOG_INFO("LINADV_ANALYSIS",*) 'Setup has been finished.'

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

    ! Analysis
    call numerror_analysis%Final()

    ! Field
    call PTracer%Final()
    call PTracer_ref%Final()
    
    ! Mesh
    call mesh3D%Final()
    call mesh1D%Final()

    ! Element
    call refElem3D%Final()
    call refElem1D%Final()

    ! stop MPI
    call PRC_mpibarrier
    call PRC_MPIfinish 

    return
  end subroutine finalize

  subroutine print_normalized_factor( q )
    type(MeshField3D), intent(in) :: q
    integer :: lcdomid
    integer :: ke

    real(RP) :: lcval1, lcval2, lcvol
    real(RP) :: vol, glval1, glval2
    real(RP) :: int_weight(refElem3D%Np)
    !-------------------------------

    lcvol = 0.0_RP
    lcval1 = 0.0_RP; lcval2 = 0.0_RP
    do lcdomid=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(lcdomid)

      !$omp parallel do private(int_weight) reduction(+:lcvol, lcval1, lcval2)
      do ke=lcmesh3D%NeS, lcmesh3D%NeE
        int_weight(:) = lcmesh3D%Gsqrt(:,ke) * lcmesh3D%J(:,ke) * refElem3D%IntWeight_lgl(:)
        lcvol = lcvol + sum( int_weight(:) )
        lcval1 = lcval1 + sum( int_weight(:) * q%local(lcdomid)%val(:,ke) )
        lcval2 = lcval2 + sum( int_weight(:) * q%local(lcdomid)%val(:,ke)**2 )
      end do
    end do
    
    vol = global_global_int(lcvol)
    glval1 = global_global_int(lcval1)
    glval2 = global_global_int(lcval2)
    LOG_INFO("Normalized Factor",*) "vol:", vol, "L1_fac:", vol/glval1, "L2_fac:", sqrt(vol/glval2)

    return
  end subroutine print_normalized_factor

!OCL SERIAL
  function global_global_int(lc_val) result(global_val)
    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD    
    implicit none

    real(RP), intent(in), target :: lc_val
    real(RP) :: global_val

    integer :: nproc
    integer :: ierr
    !-----------------------------------------

    ! global sum
    call MPI_AllReduce( lc_val, global_val, 1, &
      MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr       )
    
    call MPI_Comm_size( PRC_LOCAL_COMM_WORLD, nproc, ierr )
    return
  end function global_global_int

end program prg_linadv2d_N20Case4_analysis
