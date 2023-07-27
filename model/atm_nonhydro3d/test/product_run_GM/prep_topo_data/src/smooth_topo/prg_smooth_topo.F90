#include "scalelib.h"
program prg_smooth_topo
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort
  use scale_const, only: &
    PI => CONST_PI

  use scale_localmesh_2d, only: LocalMesh2D

  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  
  use scale_meshfield_base, only: &
    MeshField2D

  use scale_file_base_meshfield, only: &
    FILE_base_meshfield

  use scale_mesh_base2d, only: &
    DIMTYPEID_XY => MeshBase2D_DIMTYPEID_XY
    
  use mod_smooth_topo_lap, only: &
    smooth_topo_lap_apply
  
  implicit none

   ! MPI parameters
  integer                 :: nprocs                      ! number of processes               (execution)
  integer                 :: myrank                      ! my rank                           (execution)
  logical                 :: ismaster                    ! master process?                   (execution)

  type(QuadrilateralElement) :: refElem2D
  type(HexahedralElement) :: refElem3D

  type(MeshCubedSphereDom2D), target :: mesh2D

  type(MeshField2D) :: in_topography
  type(MeshField2D) :: out_topography

  type(FILE_base_meshfield) :: out_file
  integer, parameter :: vid_topo = 1
  integer, parameter :: vid_topo_ori = 2

  call initialize()

  call smooth_topo_lap_apply( out_topography, &
    in_topography )

  call out_file%Write_var2D( vid_topo_ori, in_topography, 0.0_RP, 0.0_RP )
  call out_file%Write_var2D( vid_topo, out_topography, 0.0_RP, 0.0_RP )

  call finalize()

contains
!--
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

    use mod_smooth_topo_lap, only: &
      smooth_topo_lap_Init
    
    implicit none

    integer :: ierr
    integer :: comm                     ! communicator   (execution)
    logical :: fileexist
    character(len=H_LONG) :: cnf_fname  ! config file for launcher

    character(len=H_LONG) :: in_topo_filebase
    character(len=H_LONG) :: out_topo_filebase

    integer :: IDEALIZED_TOPO_ZONAL_WAVE_NUM
      
    integer :: Nprc = 6
    integer :: NeGX = 2
    integer :: NeGY = 2

    integer, parameter :: NLocalMeshPerPrc = 1
    
    integer :: PolyOrder_h = 7

    !-
    namelist / PARAM_SMOOTH_TOPO / &
        in_topo_filebase,  &
        out_topo_filebase, &
        Nprc, NeGX, NeGY,  &
        PolyOrder_h,       &
        IDEALIZED_TOPO_ZONAL_WAVE_NUM
    
    type(FILE_base_meshfield) :: in_file
    class(LocalMesh2D), pointer :: lcmesh2D
    integer :: domid, ke
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
    call IO_setup( "SMOOTH_TOPO", cnf_fname )

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
    ! setup constants
    call CONST_setup
    ! setup fie I/O
    call FILE_setup( myrank )

    LOG_NEWLINE
    LOG_INFO("SMOOTH_TOPO",*) 'Setup'

    !-
    in_topo_filebase = ""
    out_topo_filebase = "./out_data/topo"

    IDEALIZED_TOPO_ZONAL_WAVE_NUM = 360
   
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SMOOTH_TOPO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("SMOOTH_TOPO_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("SMOOTH_TOPO_setup",*) 'Not appropriate names in namelist PARAM_SMOOTH_TOPO. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_SMOOTH_TOPO)

    call refElem2D%Init(PolyOrder_h, .false.)
    call mesh2D%Init( NeGX, NeGY, RPlanet, refElem2D, NLocalMeshPerPrc, &
      Nprc, myrank )
    call mesh2D%Generate()

    call in_topography%Init("topo", "m", mesh2D)
    call out_topography%Init("topo", "m", mesh2D)

    !-
    call out_file%Init( 1, meshCubedSphere2D=mesh2D )
    call out_file%Create( out_topo_filebase, "Interp topo data", 'REAL8',     & ! (in)
                          fileexist,                                          & ! (out)
                          myrank=PRC_myrank ) ! (in)
    if (.not. fileexist) then
      call out_file%Put_GlobalAttribute_time( (/ 0, 1, 1, 0, 0, 0 /), 0.0_RP )
    end if
    
    call out_file%Def_Var( out_topography, "topo", vid_topo, MFTYPE2D_XY, "REAL8")
    call out_file%Def_Var( "topo_ori", "m", "topo", vid_topo_ori, MFTYPE2D_XY, "REAL8")
    call out_file%End_def()

    if ( trim(in_topo_filebase) /= "" ) then
      !-- Read topography data

      LOG_INFO("SMOOTH_TOPO_setup",*) 'Read topography data:', trim(in_topo_filebase)

      call in_file%Init(1, meshCubedSphere2D=mesh2D)
      call in_file%Open(in_topo_filebase, myrank=PRC_myrank)
      call in_file%Read_Var(MFTYPE2D_XY, "topo", in_topography)

      call in_file%Close()

      ! Mask
      do domid=1, mesh2D%LOCAL_MESH_NUM
        lcmesh2D => mesh2D%lcmesh_list(domid)
        do ke = lcmesh2D%NeS, lcmesh2D%NeE
          in_topography%local(domid)%val(:,ke) = max(-1.0E-16_RP, in_topography%local(domid)%val(:,ke))
        end do
      end do      
    else
      do domid=1, mesh2D%LOCAL_MESH_NUM
        lcmesh2D => mesh2D%lcmesh_list(domid)
        do ke = lcmesh2D%NeS, lcmesh2D%NeE
          in_topography%local(domid)%val(:,ke) = cos( dble(IDEALIZED_TOPO_ZONAL_WAVE_NUM) * lcmesh2D%lon(:,ke) ) &
                                               * cos(lcmesh2D%lat(:,ke))
        end do
      end do
    end if

    !-- Initialize module for smoothing topography
    call smooth_topo_lap_Init( mesh2D )

    return
  end subroutine initialize

!OCL SERAIL
  subroutine finalize()
    use scale_prc, only: &
      PRC_mpibarrier, &
      PRC_MPIfinish
    use scale_file, only: &
      FILE_close_all

    use mod_smooth_topo_lap, only: &
      smooth_topo_lap_Final
    implicit none
    !-----------------------------------------------------------------------

    ! 
    call smooth_topo_lap_Final()

    !
    call out_file%Close()
    call out_file%Final()

    !
    call in_topography%Final()
    call out_topography%Final()
    
    call mesh2D%Final()
    call refElem2D%Final()

    ! stop MPI
    call PRC_mpibarrier
    call PRC_MPIfinish 

    return
  end subroutine finalize

end program prg_smooth_topo
