#include "scaleFElib.h"
module mod_advect1d_numerror
  use scale_precision
  use scale_prc
  use scale_io

  use scale_element_base, only: ElementBase1D
  use scale_element_line, only: LineElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_mesh_linedom1d, only: MeshLineDom1D
  use scale_meshfield_base, only: MeshField1D
  use scale_meshfield_analysis_numerror_base, only: &
    MeshFieldAnalysisNumerrorInfoBase
  use scale_meshfield_analysis_numerror, only: &
    MeshFieldAnalysisNumerror1D
  
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public types & variables
  !
  type, extends(MeshFieldAnalysisNumerrorInfoBase) :: Advect1D_Numerror_Info
    real(RP) :: ADV_VEL
    character(len=H_MID) :: InitShapeName
    real(RP) :: InitShapeParams(2)
    real(RP) :: dom_xmin
    real(RP) :: dom_xmax
    integer ::  numerror_vid(1)
    !-
    type(MeshField1D), pointer :: qtrc_exact
    type(MeshField1D), pointer :: qtrc
  end type Advect1D_Numerror_Info

  type, public :: Advect1DNumErrorAnalysis
    type(Advect1D_Numerror_Info) :: info
    type(MeshFieldAnalysisNumerror1D) :: numerror_analysis
  contains
    procedure :: Init => advect1d_numerror_Init
    procedure :: Eval => advect1d_numerror_eval
    procedure :: Final => advect1d_numerror_Final
  end type Advect1DNumErrorAnalysis
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !

contains
  !> Initialization
  subroutine advect1d_numerror_Init( this,   &
    ADV_VEL, InitShapeName, InitShapeParams, &
    mesh, elem )
    use scale_prc, only: &
       PRC_abort    
    implicit none
    class(Advect1DNumErrorAnalysis), intent(inout) :: this
    real(RP), intent(in) :: ADV_VEL
    character(*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(2)
    type(MeshLineDom1D), intent(in) :: mesh
    type(LineElement), intent(in) :: elem

    integer :: ierr

    integer :: polyOrderErrorCheck
    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    integer :: LOG_STEP_INTERVAL

    namelist / PARAM_ADVECT1D_NUMERROR / &
      polyOrderErrorCheck, &
      LOG_OUT_BASENAME,    &
      LOG_STEP_INTERVAL    
    !----------------------------------------------

    LOG_NEWLINE
    LOG_INFO("advect1d_numerror_init",*) 'Setup'

    polyOrderErrorCheck = 6
    LOG_STEP_INTERVAL   = 5

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ADVECT1D_NUMERROR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO('advect1d_numerror_init',*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR('advect1d_numerror_init',*) 'Not appropriate names in namelist PARAM_ADVECT1D_NUMERROR. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ADVECT1D_NUMERROR)

    !--
    call this%numerror_analysis%Init( &
      polyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, &
      mesh, elem, set_data_lc, this%info                        )
    call this%numerror_analysis%Regist( "q", "1", this%info%numerror_vid(1) )

    this%info%ADV_VEL = ADV_VEL
    this%info%InitShapeName = InitShapeName
    this%info%InitShapeParams = InitShapeParams
    this%info%dom_xmin = mesh%xmin_gl
    this%info%dom_xmax = mesh%xmax_gl
    return
  end subroutine advect1d_numerror_Init

  !> Evaluate numerical errors
  subroutine advect1d_numerror_eval( this, qtrc_exact, & ! (inout)
      qtrc, istep, tsec                                ) ! (in)
    implicit none
    class(Advect1DNumErrorAnalysis), intent(inout) :: this
    class(MeshField1D), intent(inout), target :: qtrc_exact
    class(MeshField1D), intent(in), target :: qtrc
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    !------------------------------------------------------------------------
    this%info%qtrc_exact => qtrc_exact
    this%info%qtrc => qtrc
    call this%numerror_analysis%Evaluate( istep, tsec ) 
  end subroutine advect1d_numerror_eval

!OCL SERIAL
  subroutine set_data_lc( analysis,  &
    q, qexact, qexact_intrp,         &
    lcmesh, elem1D, intrp_epos, tsec )
    use mod_fieldutil, only: &
      get_upwind_pos1d => fieldutil_get_upwind_pos1d,        &
      get_profile1d_tracer => fieldutil_get_profile1d_tracer 
    implicit none
    class(MeshFieldAnalysisNumerror1D), intent(in) :: analysis
    class(LocalMesh1D), intent(in) :: lcmesh
    class(ElementBase1D), intent(in) :: elem1D
    real(RP), intent(out) :: q(elem1D%Np,lcmesh%Ne,analysis%var_num)
    real(RP), intent(out) :: qexact(elem1D%Np,lcmesh%Ne,analysis%var_num)
    real(RP), intent(out) :: qexact_intrp(analysis%intrp_np,lcmesh%Ne,analysis%var_num)
    real(RP), intent(in) :: intrp_epos(analysis%intrp_np,analysis%ndim)
    real(RP), intent(in) :: tsec

    integer :: n
    integer :: ke
    integer :: vid

    real(RP) :: vx(2)
    real(RP) :: x_uwind(elem1D%Np)
    real(RP) :: x_uwind_intrp(analysis%intrp_np)
    real(RP) :: pos_intrp(analysis%intrp_np)  

    class(Advect1D_Numerror_Info), pointer :: info 
    class(MeshFieldAnalysisNumerrorInfoBase), pointer :: info_base
    !---------------------------------------------

    select type(info_base => analysis%info)
    class is (Advect1D_Numerror_Info)
      info => info_base
    end select

    n = lcmesh%lcdomID

    !$omp parallel do private(vid, ke, vx, x_uwind, x_uwind_intrp, pos_intrp)
    do ke=lcmesh%NeS, lcmesh%NeE
      vid = info%numerror_vid(1)

      x_uwind(:) = get_upwind_pos1d( lcmesh%pos_en(:,ke,1), info%ADV_VEL, tsec, info%dom_xmin, info%dom_xmax )
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      pos_intrp(:) = vx(1) + 0.5_RP * ( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
      x_uwind_intrp(:) = get_upwind_pos1d(pos_intrp(:), info%ADV_VEL, tsec, info%dom_xmin, info%dom_xmax)

      call get_profile1d_tracer( info%qtrc_exact%local(n)%val(:,ke),  & ! (out)
        info%InitShapeName, x_uwind, info%InitShapeParams, elem1D%Np  ) ! (in)

      call get_profile1d_tracer( qexact_intrp(:,ke,vid),                           & ! (out)
        info%InitShapeName, x_uwind_intrp, info%InitShapeParams, analysis%intrp_np ) ! (in)

      q(:,ke,vid) = info%qtrc%local(n)%val(:,ke)
      qexact(:,ke,vid) = info%qtrc_exact%local(n)%val(:,ke)
    end do
    return
  end subroutine set_data_lc

  !> Finalization
  subroutine advect1d_numerror_Final( this )
    implicit none
    class(Advect1DNumErrorAnalysis), intent(inout) :: this
    !---------------------------------
    call this%numerror_analysis%Final()
    return
  end subroutine advect1d_numerror_Final

end module mod_advect1d_numerror