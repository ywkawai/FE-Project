!-------------------------------------------------------------------------------
!> module Numerical error analysis for the Poisson equation in a spherical shell
!!
!! @par Description
!!      Evaluate numerical errors of the solution of the 3D Poisson equation
!!      on a cubed-sphere mesh, using MeshFieldAnalysisNumerror3D in FElib.
!!      The exact solutions are also provided by this module (shallow approximation):
!!        q = R(z) Y(lon,lat),  R = cos(k z) (Neumann) or sin(k z) (Dirichlet), k = pi / H,
!!        lap q = - ( k^2 + l(l+1)/a^2 ) q.
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_poisson3dGlobal_numerror
  use scale_precision
  use scale_prc
  use scale_io

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_meshfield_base, only: MeshField3D
  use scale_meshfield_analysis_numerror_base, only: &
    MeshFieldAnalysisNumerrorInfoBase
  use scale_meshfield_analysis_numerror, only: &
    MeshFieldAnalysisNumerror3D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public types & variables
  !
  type, extends(MeshFieldAnalysisNumerrorInfoBase) :: Poisson3DGlobal_Numerror_Info
    character(len=H_SHORT) :: TestCaseName
    character(len=H_SHORT) :: BCType
    real(RP) :: radius                       !< Radius of the sphere (bottom of the shell)
    real(RP) :: dom_zmin                     !< Bottom of the domain (height)
    real(RP) :: dom_zmax                     !< Top of the domain (height)
    integer ::  numerror_vid(1)
    !-
    type(MeshField3D), pointer :: qexact
    type(MeshField3D), pointer :: q
  end type Poisson3DGlobal_Numerror_Info

  type, public :: Poisson3DGlobalNumErrorAnalysis
    type(Poisson3DGlobal_Numerror_Info) :: info
    type(MeshFieldAnalysisNumerror3D) :: numerror_analysis
  contains
    procedure :: Init => poisson3dGlobal_numerror_Init
    procedure :: Eval => poisson3dGlobal_numerror_eval
    procedure :: Final => poisson3dGlobal_numerror_Final
  end type Poisson3DGlobalNumErrorAnalysis

  public :: poisson3dGlobal_get_exact_sol

contains
  !> Initialization
  subroutine poisson3dGlobal_numerror_Init( this, &
    TestCaseName, BCType, mesh, elem )
    implicit none
    class(Poisson3DGlobalNumErrorAnalysis), intent(inout) :: this
    character(*), intent(in) :: TestCaseName
    character(*), intent(in) :: BCType
    class(MeshCubedSphereDom3D), intent(in) :: mesh
    type(HexahedralElement), intent(in) :: elem

    integer :: ierr

    integer :: polyOrderErrorCheck
    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    integer :: LOG_STEP_INTERVAL

    namelist / PARAM_POISSON3DGLOBAL_NUMERROR / &
      polyOrderErrorCheck, &
      LOG_OUT_BASENAME,    &
      LOG_STEP_INTERVAL
    !----------------------------------------------

    LOG_NEWLINE
    LOG_INFO("poisson3dGlobal_numerror_Init",*) 'Setup'

    polyOrderErrorCheck = elem%PolyOrder_h + 3
    LOG_STEP_INTERVAL   = 1

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_POISSON3DGLOBAL_NUMERROR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO('poisson3dGlobal_numerror_Init',*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR('poisson3dGlobal_numerror_Init',*) 'Not appropriate names in namelist PARAM_POISSON3DGLOBAL_NUMERROR. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_POISSON3DGLOBAL_NUMERROR)

    !--
    call this%numerror_analysis%Init( &
      polyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, &
      mesh, elem, set_data_lc, this%info                        )
    call this%numerror_analysis%Regist( "q", "1", this%info%numerror_vid(1) )

    this%info%TestCaseName = TestCaseName
    this%info%BCType = BCType
    this%info%radius = mesh%RPlanet
    this%info%dom_zmin = mesh%zmin_gl
    this%info%dom_zmax = mesh%zmax_gl
    return
  end subroutine poisson3dGlobal_numerror_Init

  !> Evaluate numerical errors
  subroutine poisson3dGlobal_numerror_eval( this, qexact, & ! (in)
      q, istep, tsec                                      ) ! (in)
    implicit none
    class(Poisson3DGlobalNumErrorAnalysis), intent(inout) :: this
    type(MeshField3D), intent(in), target :: qexact
    type(MeshField3D), intent(in), target :: q
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    !------------------------------------------------------------------------
    this%info%qexact => qexact
    this%info%q => q
    call this%numerror_analysis%Evaluate( istep, tsec )
    return
  end subroutine poisson3dGlobal_numerror_eval

  !> Get the exact solution and the factor of the Laplacian: lap q = lap_fac * q (shallow approximation)
  !! q = R(z) Y(lon,lat), where Y is a spherical harmonic (lap_S Y = - l(l+1) Y / a^2) and
  !! R = cos(k (z-zmin)) for the Neumann condition or sin(k (z-zmin)) for the Dirichlet condition (k = pi / H).
!OCL SERIAL
  subroutine poisson3dGlobal_get_exact_sol( qexact, lap_fac, &
    TestCaseName, BCType, lon, lat, z, radius, dom_zmin, dom_zmax, Np )
    use scale_const, only: PI => CONST_PI
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: qexact(Np)
    real(RP), intent(out) :: lap_fac
    character(*), intent(in) :: TestCaseName
    character(*), intent(in) :: BCType
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: z(Np)
    real(RP), intent(in) :: radius
    real(RP), intent(in) :: dom_zmin
    real(RP), intent(in) :: dom_zmax

    real(RP) :: Y(Np), R(Np)
    real(RP) :: kz, ll1
    !---------------------------------------------

    select case( trim(TestCaseName) )
    case ( 'Y22' ) ! Y = cos^2(lat) cos(2 lon)
      Y(:) = cos( lat(:) )**2 * cos( 2.0_RP * lon(:) )
      ll1 = 6.0_RP
    case ( 'Y10' ) ! Y = sin(lat)
      Y(:) = sin( lat(:) )
      ll1 = 2.0_RP
    case default
      LOG_ERROR("poisson3dGlobal_get_exact_sol",*) 'Unknown TestCaseName: '//trim(TestCaseName)//'. Check!'
      call PRC_abort
    end select

    kz = PI / ( dom_zmax - dom_zmin )
    select case( trim(BCType) )
    case ( 'NEUMANN' )
      R(:) = cos( kz * ( z(:) - dom_zmin ) )
    case ( 'DIRICHLET' )
      R(:) = sin( kz * ( z(:) - dom_zmin ) )
    case default
      LOG_ERROR("poisson3dGlobal_get_exact_sol",*) 'Unknown BCType: '//trim(BCType)//'. Check!'
      call PRC_abort
    end select

    qexact(:) = R(:) * Y(:)
    lap_fac = - ( kz**2 + ll1 / radius**2 )
    return
  end subroutine poisson3dGlobal_get_exact_sol

!- private -----------------------------------------

!OCL SERIAL
  subroutine set_data_lc( analysis,  &
    q, qexact, qexact_intrp,         &
    lcmesh, elem3D, intrp_epos, tsec )
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatPos
    implicit none
    class(MeshFieldAnalysisNumerror3D), intent(in) :: analysis
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(out) :: q(elem3D%Np,lcmesh%Ne,analysis%var_num)
    real(RP), intent(out) :: qexact(elem3D%Np,lcmesh%Ne,analysis%var_num)
    real(RP), intent(out) :: qexact_intrp(analysis%intrp_np,lcmesh%Ne,analysis%var_num)
    real(RP), intent(in) :: intrp_epos(analysis%intrp_np,analysis%ndim)
    real(RP), intent(in) :: tsec

    integer :: n
    integer :: ke
    integer :: vid

    real(RP) :: vx(elem3D%Nv), vy(elem3D%Nv), vz(elem3D%Nv)
    real(RP) :: alpha_intrp(analysis%intrp_np)
    real(RP) :: beta_intrp(analysis%intrp_np)
    real(RP) :: z_intrp(analysis%intrp_np)
    real(RP) :: gam_intrp(analysis%intrp_np)
    real(RP) :: lon_intrp(analysis%intrp_np)
    real(RP) :: lat_intrp(analysis%intrp_np)
    real(RP) :: lap_fac

    class(Poisson3DGlobal_Numerror_Info), pointer :: info
    !---------------------------------------------

    select type(info_base => analysis%info)
    class is (Poisson3DGlobal_Numerror_Info)
      info => info_base
    end select

    n = lcmesh%lcdomID
    vid = info%numerror_vid(1)
    gam_intrp(:) = 1.0_RP

    do ke=lcmesh%NeS, lcmesh%NeE
      !- Positions of interpolation points in the (alpha, beta, z) coordinates
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
      vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)
      alpha_intrp(:) = minval(vx) + 0.5_RP * ( intrp_epos(:,1) + 1.0_RP ) * ( maxval(vx) - minval(vx) )
      beta_intrp (:) = minval(vy) + 0.5_RP * ( intrp_epos(:,2) + 1.0_RP ) * ( maxval(vy) - minval(vy) )
      z_intrp    (:) = minval(vz) + 0.5_RP * ( intrp_epos(:,3) + 1.0_RP ) * ( maxval(vz) - minval(vz) )

      call CubedSphereCoordCnv_CS2LonLatPos( lcmesh%panelID,   & ! (in)
        alpha_intrp, beta_intrp, gam_intrp, analysis%intrp_np, & ! (in)
        lon_intrp, lat_intrp                                   ) ! (out)

      call poisson3dGlobal_get_exact_sol( qexact_intrp(:,ke,vid), lap_fac,   & ! (out)
        info%TestCaseName, info%BCType, lon_intrp, lat_intrp, z_intrp,       & ! (in)
        info%radius, info%dom_zmin, info%dom_zmax, analysis%intrp_np         ) ! (in)

      q(:,ke,vid) = info%q%local(n)%val(:,ke)
      qexact(:,ke,vid) = info%qexact%local(n)%val(:,ke)
    end do

    return
  end subroutine set_data_lc

  !> Finalization
  subroutine poisson3dGlobal_numerror_Final( this )
    implicit none
    class(Poisson3DGlobalNumErrorAnalysis), intent(inout) :: this
    !---------------------------------
    call this%numerror_analysis%Final()
    return
  end subroutine poisson3dGlobal_numerror_Final

end module mod_poisson3dGlobal_numerror
