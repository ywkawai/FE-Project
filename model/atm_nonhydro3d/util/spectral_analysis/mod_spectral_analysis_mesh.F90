#include "scaleFElib.h"
module mod_spectral_analysis_mesh
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

  implicit none
  private

  type, public :: MeshList
    integer :: mesh_num
    integer :: mesh_num_x
    integer :: mesh_num_y
    type(MeshCubeDom3D), allocatable :: mesh3D_list(:)
    type(HexahedralElement) :: refElem3D
  contains
    procedure :: Init => MeshList_Init
    procedure :: Final => MeshList_Final
  end type

contains
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

    integer :: NprcX = 2
    integer :: NprcY = 2
    integer :: NeX   = 2
    integer :: NeY   = 2
    integer :: NeZ   = 4
    integer, parameter :: NLocalMeshPerPrc = 1
    
    integer :: PolyOrder_h = 7
    integer :: PolyOrder_v = 7

    real(RP) :: dom_xmin, dom_xmax
    real(RP) :: dom_ymin, dom_ymax
    real(RP) :: dom_zmin, dom_zmax

    integer :: k
    logical :: is_spec_FZ    
    integer, parameter :: FZ_nmax = 1000
    real(RP) :: FZ(FZ_nmax)

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
    
    !--
    is_spec_FZ = .true.
    do k=1, NeZ+1
      if (FZ(k) < 0.0_RP) then
        is_spec_FZ = .false.
      end if
    end do
    
    !--
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