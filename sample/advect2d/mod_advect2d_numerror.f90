#include "scalelib.h"
module mod_advect2d_numerror
  use scale_precision
  use scale_prc
  use scale_io

  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
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
  type, extends(MeshFieldAnalysisNumerrorInfoBase) :: Advect2D_Numerror_Info
    character(len=H_SHORT) :: VelTypeName
    real(RP) :: VelTypeParams(4)
    character(len=H_SHORT) :: InitShapeName
    real(RP) :: InitShapeParams(4)

    integer ::  numerror_vid(1)

    real(RP) :: dom_xmin, dom_xmax
    real(RP) :: dom_ymin, dom_ymax
    !-
    type(MeshField2D), pointer :: qtrc_exact
    type(MeshField2D), pointer :: qtrc
  end type Advect2D_Numerror_Info

  type, public :: Advect2DNumErrorAnalysis
    type(Advect2D_Numerror_Info) :: info
    type(MeshFieldAnalysisNumerror2D) :: numerror_analysis
  contains
    procedure :: Init => advect2d_numerror_Init
    procedure :: Eval => advect2d_numerror_eval
    procedure :: Final => advect2d_numerror_Final
  end type Advect2DNumErrorAnalysis

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !

contains
  !> Initialization
  subroutine advect2d_numerror_Init( this,                      & 
    VelTypeName, VelTypeParams, InitShapeName, InitShapeParams, &
    mesh, elem ) 
    use scale_prc, only: &
      PRC_abort    
    implicit none
    class(Advect2DNumErrorAnalysis), intent(inout) :: this
    character(len=*), intent(in) :: VelTypeName
    real(RP), intent(in) :: VelTypeParams(4)
    character(len=*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(4)
    type(MeshRectDom2D), intent(in) :: mesh
    type(QuadrilateralElement), intent(in) :: elem

    integer :: ierr

    integer :: PolyOrderErrorCheck
    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    integer :: LOG_STEP_INTERVAL

    namelist / PARAM_ADVECT2D_NUMERROR / &
      PolyOrderErrorCheck, &
      LOG_OUT_BASENAME,    &
      LOG_STEP_INTERVAL    
    !----------------------------------------------

    LOG_NEWLINE
    LOG_INFO("advect2d_numerror_init",*) 'Setup'

    PolyOrderErrorCheck = 6
    LOG_STEP_INTERVAL   = 5

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ADVECT2D_NUMERROR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO('advect2d_numerror_init',*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR('advect2d_numerror_init',*) 'Not appropriate names in namelist PARAM_ADVECT2D_NUMERROR. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ADVECT2D_NUMERROR)

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
    return
  end subroutine advect2d_numerror_Init

  !> Evaluate numerical errors
  subroutine advect2d_numerror_eval( this, qtrc_exact, & ! (inout)
      qtrc, istep, tsec                                ) ! (in)
    implicit none
    class(Advect2DNumErrorAnalysis), intent(inout) :: this
    class(MeshField2D), intent(inout), target :: qtrc_exact
    class(MeshField2D), intent(in), target :: qtrc
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    !------------------------------------------------------------------------
    this%info%qtrc_exact => qtrc_exact
    this%info%qtrc => qtrc
    call this%numerror_analysis%Evaluate( istep, tsec ) 
    return
  end subroutine advect2d_numerror_eval

!OCL SERIAL
  subroutine set_data_lc( analysis, q, qexact, qexact_intrp, lcmesh, elem2D, intrp_epos, tsec )
    use mod_fieldutil, only: &
      get_upwind_pos2d => fieldutil_get_upwind_pos2d,        &
      get_profile2d_tracer => fieldutil_get_profile2d_tracer 
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

    real(RP) :: vx(4), vy(4)
    real(RP) :: x_uwind(elem2D%Np), y_vwind(elem2D%Np)
    real(RP) :: x_uwind_intrp(analysis%intrp_np), y_vwind_intrp(analysis%intrp_np)
    real(RP) :: pos_intrp(analysis%intrp_np,2)  

    class(Advect2D_Numerror_Info), pointer :: info 
    class(MeshFieldAnalysisNumerrorInfoBase), pointer :: info_base
    !---------------------------------------------

    select type(info_base => analysis%info)
    class is (Advect2D_Numerror_Info)
      info => info_base
    end select

    n = lcmesh%lcdomID

    !$omp parallel do private(vid, ke, vx, vy, x_uwind, y_vwind, x_uwind_intrp, y_vwind_intrp, pos_intrp)
    do ke=lcmesh%NeS, lcmesh%NeE

      vid = info%numerror_vid(1)

      call get_upwind_pos2d( x_uwind, y_vwind, & !(out) 
        lcmesh%pos_en(:,ke,1), lcmesh%pos_en(:,ke,2), info%VelTypeName, info%VelTypeParams, tsec, & ! (in)
        info%dom_xmin, info%dom_xmax, info%dom_ymin, info%dom_ymax                                ) ! (in)

      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
      pos_intrp(:,1) = vx(1) + 0.5_RP*( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
      pos_intrp(:,2) = vy(1) + 0.5_RP*( intrp_epos(:,2) + 1.0_RP ) * ( vy(3) - vy(1) )
      call get_upwind_pos2d( x_uwind_intrp, y_vwind_intrp, & !(out) 
        pos_intrp(:,1), pos_intrp(:,2), info%VelTypeName, info%VelTypeParams, tsec, & ! (in)
        info%dom_xmin, info%dom_xmax, info%dom_ymin, info%dom_ymax                  ) ! (in)

      call get_profile2d_tracer( info%qtrc_exact%local(n)%val(:,ke),           & ! (out)
        info%InitShapeName, x_uwind, y_vwind, info%InitShapeParams, elem2D%Np  ) ! (in)

      call get_profile2d_tracer( qexact_intrp(:,ke,vid),                                          & ! (out)
        info%InitShapeName, x_uwind_intrp, y_vwind_intrp, info%InitShapeParams, analysis%intrp_np ) ! (in)

      q(:,ke,vid) = info%qtrc%local(n)%val(:,ke)
      qexact(:,ke,vid) = info%qtrc_exact%local(n)%val(:,ke)
    end do
    return
  end subroutine set_data_lc

  !> Finalization
  subroutine advect2d_numerror_Final( this )
    implicit none
    class(Advect2DNumErrorAnalysis), intent(inout) :: this
    !---------------------------------
    call this%numerror_analysis%Final()
    return
  end subroutine advect2d_numerror_Final

end module mod_advect2d_numerror