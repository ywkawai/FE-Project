#include "scalelib.h"
program prg_mktopo_data
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort
  use scale_const, only: &
    EPS => CONST_EPS, &
    PI => CONST_PI

  use scale_localmesh_2d, only: LocalMesh2D

  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D

  use scale_sparsemat, only: &
    SparseMat, sparsemat_matmul
  
  use scale_meshfield_base, only: &
    MeshField2D
  use scale_meshfieldcomm_base, only: &
    MeshFieldContainer
    use scale_meshfieldcomm_rectdom2d, only: &
    MeshFieldCommRectDom2D
  use scale_meshfieldcomm_cubedspheredom2d, only: &
    MeshFieldCommCubedSphereDom2D
  
  use scale_file_base_meshfield, only: &
    FILE_base_meshfield

  use scale_mesh_base2d, only: &
    DIMTYPEID_XY => MeshBase2D_DIMTYPEID_XY
  
  use mod_topo_nodemap, only: &
    TOPO_nodemap_info
  use mod_topo_info, only: &
    TOPO_info, &
    TOPO_info_list, node_map_info_list, &
    TOPO_info_list_setup, TOPO_info_list_final
  
  implicit none

   ! MPI parameters
  integer                 :: nprocs                      ! number of processes               (execution)
  integer                 :: myrank                      ! my rank                           (execution)
  logical                 :: ismaster                    ! master process?                   (execution)

  type(QuadrilateralElement) :: refElem2D
  type(HexahedralElement) :: refElem3D
  class(LocalMesh2D), pointer :: lcmesh2D

  type(MeshCubedSphereDom2D), target :: mesh2D
  type(MeshField2D) :: out_topography
  type(MeshField2D) :: out_i_indx
  type(MeshField2D) :: out_j_indx
  type(MeshField2D) :: out_file_i_indx
  type(MeshField2D) :: out_file_j_indx

  integer, parameter :: vid_topo = 1
  integer, parameter :: vid_file_i_indx = 2
  integer, parameter :: vid_file_j_indx = 3
  integer, parameter :: vid_i_indx = 4
  integer, parameter :: vid_j_indx = 5


  type(FILE_base_meshfield) :: out_file

  integer, parameter :: TILE_num_lon = 36
  integer, parameter :: TILE_num_lat = 18
  integer, parameter :: TILE_nlon    = 600
  integer, parameter :: TILE_nlat    = 600
  !-------------------------------------------

  call initialize()
  call interpolate_topo()
  call finalize()

contains
  subroutine interpolate_topo()
    implicit none

    integer :: domid
    !-------------------------------------------

    do domid=1, mesh2D%LOCAL_MESH_NUM
      call interpolate_topo_lcdata( out_topography%local(domid)%val, &
        node_map_info_list(domid), mesh2D%lcmesh_list(domid), refElem2D )
    end do
    
    call out_file%Write_var2D(vid_topo, out_topography, 0.0_RP, 0.0_RP)
    call out_file%Write_var2D(vid_file_i_indx, out_file_i_indx, 0.0_RP, 0.0_RP)
    call out_file%Write_var2D(vid_file_j_indx, out_file_j_indx, 0.0_RP, 0.0_RP)
    call out_file%Write_var2D(vid_i_indx, out_i_indx, 0.0_RP, 0.0_RP)
    call out_file%Write_var2D(vid_j_indx, out_j_indx, 0.0_RP, 0.0_RP)

    return
  end subroutine interpolate_topo

  subroutine interpolate_topo_lcdata( out_topo, &
    nodemap, lcmesh, elem )

    implicit none
    class(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem
    class(TOPO_nodemap_info), intent(in) :: nodemap
    real(RP), intent(out) :: out_topo(elem%Np,lcmesh%NeA)

    integer :: p, ke
    integer :: file_i, file_j, i, j

    type(TOPO_info), pointer :: tinfo

    real(RP) :: lon_, lat_
    real(RP) :: lon0, lon1, lat0, lat1

    real(RP) :: weight_lon, weight_lat
    !-------------------------------------------

    do ke=lcmesh%NeS, lcmesh%NeE
    do p=1, elem%Np
      file_i = nodemap%file_i_ind(p,ke)
      file_j = nodemap%file_j_ind(p,ke)
      i = nodemap%i_ind(p,ke)
      j = nodemap%j_ind(p,ke)

      if (file_i < 0 .or. file_j < 0) then
        LOG_INFO("interpolate_topo_lcdata",*) p,ke,"ERROR", file_i, file_j
        call PRC_abort
      end if
      if (i < 0 .or. j < 0) then
        LOG_INFO("interpolate_topo_lcdata",*) p,ke,"ERROR", file_i, file_j
        call PRC_abort
      end if

      tinfo => TOPO_info_list(file_i,file_j)
      out_topo(p,ke) = tinfo%topo(i,j)

      lon_ = lcmesh%lon(p,ke)
      lat_ = lcmesh%lat(p,ke)
      lon0 = tinfo%lon(i); lon1 = tinfo%lon(i+1) 
      lat0 = tinfo%lat(j); lat1 = tinfo%lat(j+1)

      weight_lon = ( lon1 - lon_ ) / ( lon1 - lon0 )
      weight_lat = ( lat1 - lat_ ) / ( lat1 - lat0 )

      if ( abs(lat1 - lat0) < EPS ) then
        out_topo(p,ke) = weight_lon * tinfo%topo(i,j) + ( 1.0_RP - weight_lon ) * tinfo%topo(i+1,j)
      else
        out_topo(p,ke) = &
                         weight_lat * ( weight_lon * tinfo%topo(i,j)   + ( 1.0_RP - weight_lon ) * tinfo%topo(i+1,j)   ) &
          + ( 1.0_RP - weight_lat ) * ( weight_lon * tinfo%topo(i,j+1) + ( 1.0_RP - weight_lon ) * tinfo%topo(i+1,j+1) )
      end if


      out_file_i_indx%local(1)%val(p,ke) = file_i
      out_file_j_indx%local(1)%val(p,ke) = file_i
      out_i_indx%local(1)%val(p,ke) = i
      out_j_indx%local(1)%val(p,ke) = j

    end do
    end do

    return
  end subroutine interpolate_topo_lcdata

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
    implicit none

    integer :: ierr
    integer :: comm                     ! communicator   (execution)
    logical :: fileexist
    character(len=H_LONG) :: cnf_fname  ! config file for launcher

    character(len=H_LONG) :: in_topo_filebase
    character(len=H_LONG) :: out_topo_filebase
      
    integer :: Nprc = 6
    integer :: NeGX = 2
    integer :: NeGY = 2

    integer, parameter :: NLocalMeshPerPrc = 1
    
    integer :: PolyOrder_h = 7
    class(ElementBase2D), pointer :: elem2D

    !-
    namelist / PARAM_MKTOPO_DATA / &
        in_topo_filebase,  &
        out_topo_filebase, &
        Nprc, NeGX, NeGY,  &
        PolyOrder_h
    
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
    call IO_setup( "MKTOPO_DATA", cnf_fname )

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
    ! setup constants
    call CONST_setup
    ! setup fie I/O
    call FILE_setup( myrank )
    
    LOG_NEWLINE
    LOG_INFO("MKTOPO_DATA",*) 'Setup'

    !-
    in_topo_filebase = "./in_data/topo"
    out_topo_filebase = "./out_data/topo"
   
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_DATA,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("MKTOPO_DATA_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("MKTOPO_DATA_setup",*) 'Not appropriate names in namelist PARAM_MKTOPO_DATA. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_DATA)

    call refElem2D%Init(PolyOrder_h, .false.)
    call mesh2D%Init( NeGX, NeGY, RPlanet, refElem2D, NLocalMeshPerPrc, &
      Nprc, myrank )
    call mesh2D%Generate()

    call out_topography%Init("topo", "m", mesh2D)
    call out_i_indx%Init("i_indx", "m", mesh2D)
    call out_j_indx%Init("j_indx", "m", mesh2D)
    call out_file_i_indx%Init("file_i_indx", "m", mesh2D)
    call out_file_j_indx%Init("file_j_indx", "m", mesh2D)


    !-
    call TOPO_info_list_setup( &
      TILE_num_lon, TILE_nlon, TILE_num_lat, TILE_nlat, &
      in_topo_filebase, mesh2D )
    
    !-
    call out_file%Init( 5, meshCubedSphere2D=mesh2D )
    call out_file%Create( out_topo_filebase, "Interp topo data", 'REAL8',     & ! (in)
                          fileexist,                                          & ! (out)
                          myrank=PRC_myrank ) ! (in)
    if (.not. fileexist) then
      call out_file%Put_GlobalAttribute_time( (/ 0, 1, 1, 0, 0, 0 /), 0.0_RP )
    end if
    
    call out_file%Def_Var( out_topography, "topo", vid_topo, MFTYPE2D_XY, "REAL8")
    call out_file%Def_Var( out_file_i_indx, "file_i_indx", vid_file_i_indx, MFTYPE2D_XY, "REAL8")
    call out_file%Def_Var( out_file_j_indx, "file_j_indx", vid_file_j_indx, MFTYPE2D_XY, "REAL8")
    call out_file%Def_Var( out_i_indx, "i_indx", vid_i_indx, MFTYPE2D_XY, "REAL8")
    call out_file%Def_Var( out_j_indx, "j_indx", vid_j_indx, MFTYPE2D_XY, "REAL8")

    call out_file%End_def()

    return
  end subroutine initialize

  subroutine finalize()
    use scale_prc, only: &
      PRC_mpibarrier, &
      PRC_MPIfinish
    use scale_file, only: &
      FILE_close_all
    implicit none
    !-----------------------------------------------------------------------

    !
    call TOPO_info_list_final()

    call out_file%Close()
    call out_file%Final()

    !
    call out_topography%Final()
    
    call mesh2D%Final()
    call refElem2D%Final()

    ! stop MPI
    call PRC_mpibarrier
    call PRC_MPIfinish 

    return
  end subroutine finalize

end program prg_mktopo_data