!-------------------------------------------------------------------------------
!> module Numerical error analysis for the Poisson equation on the sphere
!!
!! @par Description
!!      Evaluate numerical errors of the solution of the 2D Poisson equation
!!      on a cubed-sphere mesh, using MeshFieldAnalysisNumerror2D in FElib.
!!      The exact solutions (spherical harmonics) are also provided by this module.
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_poisson2dGlobal_numerror
  use scale_precision
  use scale_prc
  use scale_io

  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  use scale_meshfield_base, only: MeshField2D
  use scale_meshfield_analysis_numerror_base, only: &
    MeshFieldAnalysisNumerrorInfoBase
  use scale_meshfield_analysis_numerror, only: &
    MeshFieldAnalysisNumerror2D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public types & variables
  !
  type, extends(MeshFieldAnalysisNumerrorInfoBase) :: Poisson2DGlobal_Numerror_Info
    character(len=H_SHORT) :: TestCaseName
    real(RP) :: radius                       !< Radius of the sphere
    integer ::  numerror_vid(1)
    !-
    type(MeshField2D), pointer :: qexact
    type(MeshField2D), pointer :: q
  end type Poisson2DGlobal_Numerror_Info

  type, public :: Poisson2DGlobalNumErrorAnalysis
    type(Poisson2DGlobal_Numerror_Info) :: info
    type(MeshFieldAnalysisNumerror2D) :: numerror_analysis
  contains
    procedure :: Init => poisson2dGlobal_numerror_Init
    procedure :: Eval => poisson2dGlobal_numerror_eval
    procedure :: Final => poisson2dGlobal_numerror_Final
  end type Poisson2DGlobalNumErrorAnalysis

  public :: poisson2dGlobal_get_exact_sol

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !

contains
  !> Initialization
  subroutine poisson2dGlobal_numerror_Init( this, &
    TestCaseName, mesh, elem )
    implicit none
    class(Poisson2DGlobalNumErrorAnalysis), intent(inout) :: this
    character(*), intent(in) :: TestCaseName
    class(MeshCubedSphereDom2D), intent(in) :: mesh
    type(QuadrilateralElement), intent(in) :: elem

    integer :: ierr

    integer :: polyOrderErrorCheck
    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    integer :: LOG_STEP_INTERVAL

    namelist / PARAM_POISSON2DGLOBAL_NUMERROR / &
      polyOrderErrorCheck, &
      LOG_OUT_BASENAME,    &
      LOG_STEP_INTERVAL
    !----------------------------------------------

    LOG_NEWLINE
    LOG_INFO("poisson2dGlobal_numerror_Init",*) 'Setup'

    polyOrderErrorCheck = elem%PolyOrder + 3
    LOG_STEP_INTERVAL   = 1

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_POISSON2DGLOBAL_NUMERROR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO('poisson2dGlobal_numerror_Init',*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR('poisson2dGlobal_numerror_Init',*) 'Not appropriate names in namelist PARAM_POISSON2DGLOBAL_NUMERROR. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_POISSON2DGLOBAL_NUMERROR)

    !--
    call this%numerror_analysis%Init( &
      polyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, &
      mesh, elem, set_data_lc, this%info                        )
    call this%numerror_analysis%Regist( "q", "1", this%info%numerror_vid(1) )

    this%info%TestCaseName = TestCaseName
    this%info%radius = mesh%RPlanet
    return
  end subroutine poisson2dGlobal_numerror_Init

  !> Evaluate numerical errors
  subroutine poisson2dGlobal_numerror_eval( this, qexact, & ! (in)
      q, istep, tsec                                      ) ! (in)
    implicit none
    class(Poisson2DGlobalNumErrorAnalysis), intent(inout) :: this
    type(MeshField2D), intent(in), target :: qexact
    type(MeshField2D), intent(in), target :: q
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    !------------------------------------------------------------------------
    this%info%qexact => qexact
    this%info%q => q
    call this%numerror_analysis%Evaluate( istep, tsec )
    return
  end subroutine poisson2dGlobal_numerror_eval

  !> Get the exact solution and the factor of the Laplacian: lap_S q = lap_fac * q
  !! Spherical harmonics satisfy lap_S Y_l^m = - l(l+1) Y_l^m / a^2
!OCL SERIAL
  subroutine poisson2dGlobal_get_exact_sol( qexact, lap_fac, &
    TestCaseName, lon, lat, radius, Np                       )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: qexact(Np)
    real(RP), intent(out) :: lap_fac
    character(*), intent(in) :: TestCaseName
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: radius
    !---------------------------------------------

    select case( trim(TestCaseName) )
    case ( 'Y22' ) ! q = cos^2(lat) cos(2 lon)
      qexact(:) = cos( lat(:) )**2 * cos( 2.0_RP * lon(:) )
      lap_fac = - 6.0_RP / radius**2
    case ( 'Y10' ) ! q = sin(lat)
      qexact(:) = sin( lat(:) )
      lap_fac = - 2.0_RP / radius**2
    case default
      LOG_ERROR("poisson2dGlobal_get_exact_sol",*) 'Unknown TestCaseName: '//trim(TestCaseName)//'. Check!'
      call PRC_abort
    end select
    return
  end subroutine poisson2dGlobal_get_exact_sol

!- private -----------------------------------------

!OCL SERIAL
  subroutine set_data_lc( analysis,  &
    q, qexact, qexact_intrp,         &
    lcmesh, elem2D, intrp_epos, tsec )
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatPos
    implicit none
    class(MeshFieldAnalysisNumerror2D), intent(in) :: analysis
    class(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: q(elem2D%Np,lcmesh%Ne,analysis%var_num)
    real(RP), intent(out) :: qexact(elem2D%Np,lcmesh%Ne,analysis%var_num)
    real(RP), intent(out) :: qexact_intrp(analysis%intrp_np,lcmesh%Ne,analysis%var_num)
    real(RP), intent(in) :: intrp_epos(analysis%intrp_np,analysis%ndim)
    real(RP), intent(in) :: tsec

    integer :: n
    integer :: ke
    integer :: vid

    real(RP) :: vx(elem2D%Nv), vy(elem2D%Nv)
    real(RP) :: alpha_intrp(analysis%intrp_np)
    real(RP) :: beta_intrp(analysis%intrp_np)
    real(RP) :: gam_intrp(analysis%intrp_np)
    real(RP) :: lon_intrp(analysis%intrp_np)
    real(RP) :: lat_intrp(analysis%intrp_np)
    real(RP) :: lap_fac

    class(Poisson2DGlobal_Numerror_Info), pointer :: info
    !---------------------------------------------

    select type(info_base => analysis%info)
    class is (Poisson2DGlobal_Numerror_Info)
      info => info_base
    end select

    n = lcmesh%lcdomID
    vid = info%numerror_vid(1)
    gam_intrp(:) = 1.0_RP

    do ke=lcmesh%NeS, lcmesh%NeE
      !- Positions of interpolation points in the cubed-sphere (alpha, beta) coordinates
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
      alpha_intrp(:) = vx(1) + 0.5_RP * ( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
      beta_intrp (:) = vy(1) + 0.5_RP * ( intrp_epos(:,2) + 1.0_RP ) * ( vy(3) - vy(1) )

      call CubedSphereCoordCnv_CS2LonLatPos( lcmesh%panelID, & ! (in)
        alpha_intrp, beta_intrp, gam_intrp, analysis%intrp_np, & ! (in)
        lon_intrp, lat_intrp                                   ) ! (out)

      call poisson2dGlobal_get_exact_sol( qexact_intrp(:,ke,vid), lap_fac, & ! (out)
        info%TestCaseName, lon_intrp, lat_intrp, info%radius, analysis%intrp_np ) ! (in)

      q(:,ke,vid) = info%q%local(n)%val(:,ke)
      qexact(:,ke,vid) = info%qexact%local(n)%val(:,ke)
    end do

    return
  end subroutine set_data_lc

  !> Finalization
  subroutine poisson2dGlobal_numerror_Final( this )
    implicit none
    class(Poisson2DGlobalNumErrorAnalysis), intent(inout) :: this
    !---------------------------------
    call this%numerror_analysis%Final()
    return
  end subroutine poisson2dGlobal_numerror_Final

end module mod_poisson2dGlobal_numerror
