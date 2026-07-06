#include "scalelib.h"
module mod_advect3d_numerror
  use scale_precision
  use scale_prc
  use scale_io

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
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
  type, extends(MeshFieldAnalysisNumerrorInfoBase) :: Advect3D_Numerror_Info
    character(len=H_SHORT) :: VelTypeName
    real(RP) :: VelTypeParams(6)
    character(len=H_SHORT) :: InitShapeName
    real(RP) :: InitShapeParams(6)

    integer ::  numerror_vid(1)

    real(RP) :: dom_xmin, dom_xmax
    real(RP) :: dom_ymin, dom_ymax
    real(RP) :: dom_zmin, dom_zmax
    !-
    type(MeshField3D), pointer :: qtrc_exact
    type(MeshField3D), pointer :: qtrc
  end type Advect3D_Numerror_Info

  type, public :: Advect3DNumErrorAnalysis
    type(Advect3D_Numerror_Info) :: info
    type(MeshFieldAnalysisNumerror3D) :: numerror_analysis
  contains
    procedure :: Init => advect3d_numerror_Init
    procedure :: Eval => advect3d_numerror_eval
    procedure :: Final => advect3d_numerror_Final
  end type Advect3DNumErrorAnalysis

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !

contains
  !> Initialization
  subroutine advect3d_numerror_Init( this,                      & 
    VelTypeName, VelTypeParams, InitShapeName, InitShapeParams, &
    mesh, elem ) 
    use scale_prc, only: &
      PRC_abort    
    implicit none
    class(Advect3DNumErrorAnalysis), intent(inout) :: this
    character(len=*), intent(in) :: VelTypeName
    real(RP), intent(in) :: VelTypeParams(6)
    character(len=*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(6)
    type(MeshCubeDom3D), intent(in) :: mesh
    type(HexahedralElement), intent(in) :: elem

    integer :: ierr

    integer :: PolyOrderErrorCheck
    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    integer :: LOG_STEP_INTERVAL

    namelist / PARAM_ADVECT3D_NUMERROR / &
      PolyOrderErrorCheck, &
      LOG_OUT_BASENAME,    &
      LOG_STEP_INTERVAL    
    !----------------------------------------------

    LOG_NEWLINE
    LOG_INFO("advect3d_numerror_init",*) 'Setup'

    PolyOrderErrorCheck = 6
    LOG_STEP_INTERVAL   = 5

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ADVECT3D_NUMERROR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO('advect3d_numerror_init',*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR('advect3d_numerror_init',*) 'Not appropriate names in namelist PARAM_ADVECT3D_NUMERROR. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ADVECT3D_NUMERROR)

    !--
    call this%numerror_analysis%Init( &
      polyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, &
      mesh, elem, set_data_lc, this%info                        )
    call this%numerror_analysis%Regist( "q", "1", this%info%numerror_vid(1) )

    this%info%VelTypeName = VelTypeName
    this%info%VelTypeParams = VelTypeParams
    this%info%InitShapeName = InitShapeName
    this%info%InitShapeParams = InitShapeParams
    this%info%dom_xmin = mesh%xmin_gl; this%info%dom_xmax = mesh%xmax_gl
    this%info%dom_ymin = mesh%ymin_gl; this%info%dom_ymax = mesh%ymax_gl
    this%info%dom_zmin = mesh%zmin_gl; this%info%dom_zmax = mesh%zmax_gl
    return
  end subroutine advect3d_numerror_Init

  !> Evaluate numerical errors
  subroutine advect3d_numerror_eval( this, qtrc_exact, & ! (inout)
      qtrc, istep, tsec                                ) ! (in)
    implicit none
    class(Advect3DNumErrorAnalysis), intent(inout) :: this
    class(MeshField3D), intent(inout), target :: qtrc_exact
    class(MeshField3D), intent(in), target :: qtrc
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    !------------------------------------------------------------------------
    this%info%qtrc_exact => qtrc_exact
    this%info%qtrc => qtrc
    call this%numerror_analysis%Evaluate( istep, tsec ) 
    return
  end subroutine advect3d_numerror_eval

!OCL SERIAL
  subroutine set_data_lc( analysis, q, qexact, qexact_intrp, lcmesh, elem3D, intrp_epos, tsec )
    use mod_fieldutil, only: &
      get_upwind_pos3d => fieldutil_get_upwind_pos3d,         &
      get_profile3d_tracer => fieldutil_get_profile3d_tracer, &
      get_profile3d_flow => fieldutil_get_profile3d_flow
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

    real(RP) :: vx(8), vy(8), vz(8)
    real(RP) :: x_uwind(elem3D%Np), y_vwind(elem3D%Np), z_wwind(elem3D%Np)
    real(RP) :: x_uwind_intrp(analysis%intrp_np), y_vwind_intrp(analysis%intrp_np), z_wwind_intrp(analysis%intrp_np)
    real(RP) :: pos_intrp(analysis%intrp_np,3)  

    class(Advect3D_Numerror_Info), pointer :: info 
    class(MeshFieldAnalysisNumerrorInfoBase), pointer :: info_base
    !---------------------------------------------

    select type(info_base => analysis%info)
    class is (Advect3D_Numerror_Info)
      info => info_base
    end select

    n = lcmesh%lcdomID
    !$acc update host( info%qtrc%local(n)%val )

    !$omp parallel do private(vid, ke, vx, vy, vz, x_uwind, y_vwind, z_wwind, x_uwind_intrp, y_vwind_intrp, z_wwind_intrp, pos_intrp)
    do ke=lcmesh%NeS, lcmesh%NeE

      vid = info%numerror_vid(1)

      call get_upwind_pos3d( x_uwind, y_vwind, z_wwind, & !(out) 
        lcmesh%pos_en(:,ke,1), lcmesh%pos_en(:,ke,2), lcmesh%pos_en(:,ke,3), info%VelTypeName, info%VelTypeParams, tsec, & ! (in)
        info%dom_xmin, info%dom_xmax, info%dom_ymin, info%dom_ymax, info%dom_zmin, info%dom_zmax                         ) ! (in)

      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
      vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)
      pos_intrp(:,1) = vx(1) + 0.5_RP*( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
      pos_intrp(:,2) = vy(1) + 0.5_RP*( intrp_epos(:,2) + 1.0_RP ) * ( vy(3) - vy(1) )
      pos_intrp(:,3) = vz(1) + 0.5_RP*( intrp_epos(:,3) + 1.0_RP ) * ( vz(5) - vz(1) )
      call get_upwind_pos3d( x_uwind_intrp, y_vwind_intrp, z_wwind_intrp,                           & !(out) 
        pos_intrp(:,1), pos_intrp(:,2), pos_intrp(:,3), info%VelTypeName, info%VelTypeParams, tsec, & ! (in)
        info%dom_xmin, info%dom_xmax, info%dom_ymin, info%dom_ymax, info%dom_zmin, info%dom_zmax    ) ! (in)

      call get_profile3d_tracer( info%qtrc_exact%local(n)%val(:,ke),                    & ! (out)
        info%InitShapeName, x_uwind, y_vwind, z_wwind, info%InitShapeParams, elem3D%Np  ) ! (in)

      call get_profile3d_tracer( qexact_intrp(:,ke,vid),                                                         & ! (out)
        info%InitShapeName, x_uwind_intrp, y_vwind_intrp, z_wwind_intrp, info%InitShapeParams, analysis%intrp_np ) ! (in)

      q(:,ke,vid) = info%qtrc%local(n)%val(:,ke)
      qexact(:,ke,vid) = info%qtrc_exact%local(n)%val(:,ke)
    end do

    !$acc update device( info%qtrc_exact%local(n)%val )    
    return
  end subroutine set_data_lc

  !> Finalization
  subroutine advect3d_numerror_Final( this )
    implicit none
    class(Advect3DNumErrorAnalysis), intent(inout) :: this
    !---------------------------------
    call this%numerror_analysis%Final()
    return
  end subroutine advect3d_numerror_Final

end module mod_advect3d_numerror