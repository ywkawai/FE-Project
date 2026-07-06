!-------------------------------------------------------------------------------
!> module to manage meshes for spectral analysis
!!
!! @author Yuta Kawai, Team SCALE
!<-
#include "scaleFElib.h"
module mod_spectral_analysis_mesh
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort
  use scale_const, only: &
    PI => CONST_PI,   &
    RPlanet => CONST_RADIUS

  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D

  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage multiple meshes for spectral analysis
  type, public :: MeshList
    integer :: mesh_num     !< Number of meshes
    integer :: mesh_num_x   !< Number of meshes in x-direction
    integer :: mesh_num_y   !< Number of meshes in y-direction

    type(MeshCubeDom3D), allocatable :: mesh3D_list(:)  !< Array of 3D meshes for spectral analysis
    type(HexahedralElement) :: refElem3D                !< Object of reference element for 3D mesh
  contains
    procedure :: Init => MeshList_Init
    procedure :: Final => MeshList_Final
  end type

contains

  !> Initialize an object to manage multiple meshes for spectral analysis
  !! @param[inout] this Object of MeshList
  !! @param[in] target_proc_num Number of processes to be used for spectral analysis
  !! @param[in] target_proc_s Starting rank of processes to be used for spectral analysis
!OCL SERIAL
  subroutine MeshList_init( this, target_proc_num, target_proc_s )   
    use scale_element_quadrilateral, only: &
      QuadrilateralElement 
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatPos
    implicit none
    class(MeshList), intent(inout) :: this
    integer, intent(in) :: target_proc_num
    integer, intent(in) :: target_proc_s

    logical :: SHALLOW_ATM_APPROX_FLAG = .true.

    integer :: NprcX = 2   !< Number of MPI processes in x-direction
    integer :: NprcY = 2   !< Number of MPI processes in y-direction
    integer :: NeX   = 2   !< Number of elements in x-direction
    integer :: NeY   = 2   !< Number of elements in y-direction
    integer :: NeZ   = 4   !< Number of elements in z-direction
    integer, parameter :: NLocalMeshPerPrc = 1  !< Number of local meshes per process
    
    integer :: PolyOrder_h = 7   !< Polynomial order in horizontal direction
    integer :: PolyOrder_v = 7   !< Polynomial order in vertical direction

    real(RP) :: dom_xmin   !< Minimum x-coordinate of the domain
    real(RP) :: dom_xmax   !< Maximum x-coordinate of the domain
    real(RP) :: dom_ymin   !< Minimum y-coordinate of the domain
    real(RP) :: dom_ymax   !< Maximum y-coordinate of the domain
    real(RP) :: dom_zmin   !< Minimum z-coordinate of the domain
    real(RP) :: dom_zmax   !< Maximum z-coordinate of the domain

    integer, parameter :: FZ_nmax = 1000
    real(RP) :: FZ(FZ_nmax)  !< Array of vertical coordinates

    namelist / PARAM_SPECTRAL_ANALYSIS_MESH / &
        dom_xmin, dom_xmax,          &
        dom_ymin, dom_ymax,          &
        dom_zmin, dom_zmax,          &
        NprcX, NprcY, NeX, NeY, NeZ, &
        PolyOrder_h, PolyOrder_v,    &
        FZ
        
    integer :: ierr
    
    integer :: m
    integer :: target_myrank

    integer :: k
    logical :: is_spec_FZ 
    !-------------------------------------------

    FZ(:) = - 1.0_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SPECTRAL_ANALYSIS_MESH,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("SPECTRAL_ANALYSIS_MESH_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("SPECTRAL_ANALYSIS_MESH_setup",*) 'Not appropriate names in namelist PARAM_SPECTRAL_ANALYSIS_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_SPECTRAL_ANALYSIS_MESH)
    
    !- Check if FZ is specified

    is_spec_FZ = .true.
    do k=1, NeZ+1
      if (FZ(k) < 0.0_RP) then
        is_spec_FZ = .false.
      end if
    end do
    
    !- Setup meshes
    
    call this%refElem3D%Init(PolyOrder_h, PolyOrder_v, .false.)

    this%mesh_num = target_proc_num
    allocate( this%mesh3D_list(this%mesh_num) )
    do m=1, this%mesh_num
      target_myrank = target_proc_s + m - 1
      if (is_spec_FZ) then
        call this%mesh3D_list(m)%Init( NprcX*NeX, NprcY*NeY, NeZ, &
          dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,              &
          .true., .true., .false., this%refElem3D, NLocalMeshPerPrc, NprcX, NprcY, &
          nproc=NprcX*NprcY, FZ=FZ(1:NeZ+1), myrank=target_myrank )
      else
        call this%mesh3D_list(m)%Init( NprcX*NeX, NprcY*NeY, NeZ, &
          dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,              &
          .true., .true., .false., this%refElem3D, NLocalMeshPerPrc, NprcX, NprcY, &
          nproc=NprcX*NprcY, myrank=target_myrank )
      end if
      call this%mesh3D_list(m)%Generate()
    end do

    !-
    this%mesh_num_x = target_proc_num
    this%mesh_num_y = 1
    return
  end subroutine MeshList_init

  !> Finalize an object to manage multiple meshes for spectral analysis
!OCL SERIAL
  subroutine MeshList_final( this )
    use scale_prc
    implicit none
    class(MeshList), intent(inout) :: this
    integer :: m
    !------------------------------------------

    do m=1, this%mesh_num
      call this%mesh3D_list(m)%Final()
    end do
    deallocate( this%mesh3D_list )

    call this%refElem3D%Final()
    return
  end subroutine MeshList_final

end module mod_spectral_analysis_mesh