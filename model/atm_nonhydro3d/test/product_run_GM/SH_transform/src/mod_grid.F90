#include "scalelib.h"
module mod_grid
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
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D

  implicit none
  private


  type(HexahedralElement) :: refElem3D
  real(RP), public, allocatable :: lon2D(:,:,:)
  real(RP), public, allocatable :: lat2D(:,:,:)
  real(RP), public, allocatable :: lon2D_intrp(:,:,:)
  real(RP), public, allocatable :: lat2D_intrp(:,:,:)
  real(RP), public, allocatable :: Gsqrt2D(:,:,:)
  real(RP), public, allocatable :: J2D(:,:,:)

  real(RP), public, allocatable :: intrp_intw2D(:)
  real(RP), public, allocatable :: IntrpMat2D(:,:)

  integer, public :: Ne2D
  class(ElementBase2D), public, pointer :: refElem2D

  integer, public :: mesh_num
  type(MeshCubedSphereDom3D), public, allocatable, target :: mesh3D_list(:)

  public :: grid_init
  public :: grid_final

contains
!OCL SERIAL
  subroutine grid_init( target_proc_num, target_proc_s )   
    use scale_element_quadrilateral, only: &
      QuadrilateralElement 
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatPos
    implicit none
    integer, intent(in) :: target_proc_num
    integer, intent(in) :: target_proc_s

    logical :: SHALLOW_ATM_APPROX_FLAG = .true.

    integer :: Nprc = 6
    integer :: NeGX = 2
    integer :: NeGY = 2
    integer :: NeZ  = 12
    integer, parameter :: NLocalMeshPerPrc = 1
    
    integer :: PolyOrder_h = 7
    integer :: PolyOrder_v = 7
    integer :: IntrpPolyOrder_h = 7
    real(RP) :: dom_zmin, dom_zmax

    integer :: k
    logical :: is_spec_FZ    
    integer, parameter :: FZ_nmax = 1000
    real(RP) :: FZ(FZ_nmax)

    namelist / PARAM_SH_MESH / &
        SHALLOW_ATM_APPROX_FLAG,  &
        dom_zmin, dom_zmax,       &
        Nprc, NeGX, NeGY, NeZ,    &
        PolyOrder_h, PolyOrder_v, &
        IntrpPolyOrder_h,         &
        FZ
    integer :: ierr
    
    integer :: m
    integer :: target_myrank

    class(LocalMesh2D), pointer :: lcmesh2D
    integer :: ke2D
    integer :: Np2D

    real(RP), allocatable :: intrp_x1(:)
    real(RP), allocatable :: intrp_x2(:)
    real(RP), allocatable :: vx(:), vy(:)
    real(RP), allocatable :: alph2D(:,:), beta2D(:,:)
    real(RP), allocatable :: gam(:,:)
    !-------------------------------------------

    FZ(:) = - 1.0_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SH_MESH,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("SH_MESH_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("SH_MESH_setup",*) 'Not appropriate names in namelist PARAM_SH_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_SH_MESH)
    
    !--
    is_spec_FZ = .true.
    do k=1, NeZ+1
      if (FZ(k) < 0.0_RP) then
        is_spec_FZ = .false.
      end if
    end do
    
    !--
    call refElem3D%Init(PolyOrder_h, PolyOrder_v, .false.)

    mesh_num = target_proc_num
    allocate( mesh3D_list(mesh_num) )
    do m=1, mesh_num
      target_myrank = target_proc_s + m - 1
      if (is_spec_FZ) then
        call mesh3D_list(m)%Init( NeGX, NeGY, NeZ, RPlanet, dom_zmin, dom_zmax, refElem3D, NLocalMeshPerPrc, &
          nproc=Nprc, FZ=FZ(1:NeZ+1), myrank=target_myrank, shallow_approx=SHALLOW_ATM_APPROX_FLAG )
      else
        call mesh3D_list(m)%Init( NeGX, NeGY, NeZ, RPlanet, dom_zmin, dom_zmax, refElem3D, NLocalMeshPerPrc, &
          nproc=Nprc, myrank=target_myrank, shallow_approx=SHALLOW_ATM_APPROX_FLAG )
      end if
      call mesh3D_list(m)%Generate()
    end do

    !---
    Ne2D = mesh3D_list(1)%mesh2D%lcmesh_list(1)%Ne
    refElem2D => mesh3D_list(1)%mesh2D%refElem2D
    Np2D = refElem2D%Np

    allocate( lon2D(Np2D,Ne2D,mesh_num) )
    allocate( lat2D(Np2D,Ne2D,mesh_num) )
    allocate( Gsqrt2D(Np2D,Ne2D,mesh_num) )
    allocate( J2D(Np2D,Ne2D,mesh_num) )

    do m=1, target_proc_num
      lcmesh2D => mesh3D_list(m)%mesh2D%lcmesh_list(1)
  
      !$omp parallel do
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
        lon2D(:,ke2D,m) = lcmesh2D%lon(:,ke2D)
        lat2D(:,ke2D,m) = lcmesh2D%lat(:,ke2D)
        Gsqrt2D(:,ke2D,m) = lcmesh2D%Gsqrt(:,ke2D)
        J2D(:,ke2D,m) = lcmesh2D%J(:,ke2D)
      end do
    end do

    !--
    allocate( intrp_intw2D(IntrpPolyOrder_h**2) )
    allocate( intrpMat2D(IntrpPolyOrder_h**2,refElem2D%Np) )
    allocate( lon2D_intrp(IntrpPolyOrder_h**2,Ne2D,mesh_num) )
    allocate( lat2D_intrp(IntrpPolyOrder_h**2,Ne2D,mesh_num) )

    allocate( intrp_x1(IntrpPolyOrder_h**2), intrp_x2(IntrpPolyOrder_h**2) )
    IntrpMat2D(:,:) = mesh3D_list(1)%refElem2D%GenIntGaussLegendreIntrpMat( IntrpPolyOrder_h, &
      intw_intrp=intrp_intw2D, x_intrp=intrp_x1, y_intrp=intrp_x2 )

    allocate( vx(refElem2D%Nv), vy(refElem2D%Nv) )
    allocate( alph2D(IntrpPolyOrder_h**2,lcmesh2D%Ne), beta2D(IntrpPolyOrder_h**2,lcmesh2D%Ne) )    
    allocate( gam(IntrpPolyOrder_h**2,lcmesh2D%Ne) )
    do m=1, target_proc_num
      lcmesh2D => mesh3D_list(m)%mesh2D%lcmesh_list(1)

      !$omp parallel do private(vx, vy)
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
        vx(:) = lcmesh2D%pos_ev(lcmesh2D%EToV(ke2D,:),1)
        vy(:) = lcmesh2D%pos_ev(lcmesh2D%EToV(ke2D,:),2)
        alph2D(:,ke2D) = vx(1) + 0.5_RP * ( intrp_x1(:) + 1.0_RP ) * ( vx(2) - vx(1) ) 
        beta2D(:,ke2D) = vy(1) + 0.5_RP * ( intrp_x2(:) + 1.0_RP ) * ( vy(4) - vy(1) )
        gam(:,ke2D) = 1.0_RP
      end do
      call CubedSphereCoordCnv_CS2LonLatPos( &
        lcmesh2D%panelID, alph2D, beta2D, gam, IntrpPolyOrder_h**2 * lcmesh2D%Ne, & ! (in)
        lon2D_intrp(:,:,m), lat2D_intrp(:,:,m)                                    ) ! (out)
    end do

    return
  end subroutine grid_init

!OCL SERIAL
  subroutine grid_final()
    use scale_prc
    implicit none
    integer :: m
    !------------------------------------------

    do m=1, mesh_num
      call mesh3D_list(m)%Final()
    end do
    deallocate( mesh3D_list )

    call refElem3D%Final()

    deallocate( lon2D_intrp, lat2D_intrp )
    deallocate( intrp_intw2D, IntrpMat2D )

    return
  end subroutine grid_final

end module mod_grid