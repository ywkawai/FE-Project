#include "scaleFElib.h"
module mod_advect1d_fvm_numerror
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
  !++ Public parameters & variables
  !
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
    real(RP), pointer :: qtrc_exact(:)
    real(RP), pointer :: qtrc(:)

    !-
    type(MeshLineDom1D) :: mesh_dg
    type(LineElement) :: elem

    !-
    integer :: FV_IS, FV_IE, FV_IA, FV_IHALO
  end type Advect1D_Numerror_Info

  type, public :: Advect1DNumErrorAnalysis
    type(Advect1D_Numerror_Info) :: info
    type(MeshFieldAnalysisNumerror1D) :: numerror_analysis
  contains
    procedure :: Init => advect1d_fvm_numerror_Init
    procedure :: Eval => advect1d_fvm_numerror_eval
    procedure :: Final => advect1d_fvm_numerror_Final
  end type Advect1DNumErrorAnalysis

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !

contains
  !> Initialization
  subroutine advect1d_fvm_numerror_Init( this, &
    ADV_VEL, InitShapeName, InitShapeParams,   &
    FX, IS, IE, IA, IHALO )
    use scale_prc, only: &
      PRC_abort    
    implicit none
    class(Advect1DNumErrorAnalysis), intent(inout) :: this
    real(RP), intent(in) :: ADV_VEL
    character(*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(2)
    integer, intent(in) :: IS, IE, IA, IHALO
    real(RP), intent(in) :: FX(0:IA)

    integer :: polyOrderErrorCheck
    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    integer :: LOG_STEP_INTERVAL

    namelist / PARAM_ADVECT1D_FVM_NUMERROR / &
      polyOrderErrorCheck, &
      LOG_OUT_BASENAME,    &
      LOG_STEP_INTERVAL    
    integer :: ierr
    !----------------------------------------------

    LOG_NEWLINE
    LOG_INFO("advect1d_fvm_numerror_init",*) 'Setup'

    polyOrderErrorCheck = 2
    LOG_STEP_INTERVAL   = 5

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ADVECT1D_FVM_NUMERROR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO('advect1d_fvm_numerror_init',*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR('advect1d_fvm_numerror_init',*) 'Not appropriate names in namelist PARAM_ADVECT1D_FVM_NUMERROR. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ADVECT1D_FVM_NUMERROR)

    !--
    call this%info%elem%Init( polyOrderErrorCheck, .false. )
    call this%info%mesh_dg%Init( IE-IS+1, FX(IS-1), FX(IE), this%info%elem, 1 )
    call this%info%mesh_dg%Generate()

    !--
    call this%numerror_analysis%Init( &
      polyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, &
      this%info%mesh_dg, this%info%elem, set_data_lc, this%info )
    call this%numerror_analysis%Regist( "q", "1", this%info%numerror_vid(1) )

    !-
    this%info%ADV_VEL = ADV_VEL
    this%info%InitShapeName = InitShapeName
    this%info%InitShapeParams = InitShapeParams
    this%info%dom_xmin = this%info%mesh_dg%xmin_gl
    this%info%dom_xmax = this%info%mesh_dg%xmax_gl
    this%info%FV_IS = IS
    this%info%FV_IE = IE
    this%info%FV_IA = IA
    this%info%FV_IHALO = IHALO
    return
  end subroutine advect1d_fvm_numerror_Init

  !> Evaluate numerical errors
  subroutine advect1d_fvm_numerror_eval( this, &
    qtrc_exact,            & ! (inout)
    qtrc, istep, tsec, IA  ) ! (in)

    implicit none
    class(Advect1DNumErrorAnalysis), intent(inout) :: this
    integer, intent(in) :: IA
    real(RP), intent(out), target :: qtrc_exact(IA)
    real(RP), intent(in), target :: qtrc(IA)
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    !------------------------------------------------------------------------

    this%info%qtrc_exact => qtrc_exact
    this%info%qtrc => qtrc
    call this%numerror_analysis%Evaluate( istep, tsec ) 
    return
  end subroutine advect1d_fvm_numerror_eval

!OCL SERIAL
  subroutine set_data_lc( analysis, q, qexact, qexact_intrp, lcmesh, elem1D, intrp_epos, tsec )
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
    integer :: kelem

    integer :: vid

    real(RP) :: vx(2)
    real(RP) :: x_uwind(elem1D%Np)
    real(RP) :: x_uwind_intrp(analysis%intrp_np)
    real(RP) :: pos_intrp(analysis%intrp_np)  
    real(RP) :: int_w(elem1D%Np)
    real(RP) :: qtrc_exact_dg(elem1D%Np)
    real(RP) :: qtrc_dg(elem1D%Np,lcmesh%NeA)

    class(Advect1D_Numerror_Info), pointer :: info 
    class(MeshFieldAnalysisNumerrorInfoBase), pointer :: info_base
    !---------------------------------------------

    select type(info_base => analysis%info)
    class is (Advect1D_Numerror_Info)
      info => info_base
    end select

    n = lcmesh%lcdomID

    call interp_FVtoDG_1D( qtrc_dg, &
      info%qtrc, lcmesh, lcmesh%refElem1D, info%FV_IS, info%FV_IE, info%FV_IA, info%FV_IHALO, 2 )

    !$omp parallel do private(vid, kelem, vx, x_uwind, x_uwind_intrp, pos_intrp, int_w, &
    !$omp qtrc_exact_dg )
    do kelem=lcmesh%NeS, lcmesh%NeE

      vid = info%numerror_vid(1)

      x_uwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,kelem,1), info%ADV_VEL, tsec, info%dom_xmin, info%dom_xmax)
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(kelem,:),1)
      pos_intrp(:) = vx(1) + 0.5_RP*( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
      x_uwind_intrp(:) = get_upwind_pos1d(pos_intrp(:), info%ADV_VEL, tsec, info%dom_xmin, info%dom_xmax)

      call get_profile1d_tracer( qtrc_exact_dg(:),                   & ! (out)
        info%InitShapeName, x_uwind, info%InitShapeParams, elem1D%Np ) ! (in)

      call get_profile1d_tracer( qexact_intrp(:,kelem,vid),                       & ! (out)
        info%InitShapeName, x_uwind_intrp, info%InitShapeParams(:), analysis%intrp_np ) ! (in)

      int_w(:) = lcmesh%Gsqrt(:,kelem) * lcmesh%J(:,kelem) * elem1D%IntWeight_lgl(:)
      int_w(:) = int_w(:) / sum(int_w(:))
      info%qtrc_exact(info%FV_IHALO+kelem) = sum( int_w(:) * qtrc_exact_dg(:) )

      q(:,kelem,vid) = qtrc_dg(:,kelem)
      qexact(:,kelem,vid) = qtrc_exact_dg(:)
    end do
    
    return
  end subroutine set_data_lc

  !> Finalization
  subroutine advect1d_fvm_numerror_Final( this )
    implicit none
    class(Advect1DNumErrorAnalysis), intent(inout) :: this
    !---------------------------------
    call this%info%mesh_dg%Final()
    call this%info%elem%Final()
    call this%numerror_analysis%Final()
    return
  end subroutine advect1d_fvm_numerror_Final

!- private
  subroutine interp_FVtoDG_1D( out_dg, &
    in_fv, lcmesh, elem1D_, IS, IE, IA, IHALO, interp_ord )
    implicit none
    class(LocalMesh1D), intent(in) :: lcmesh
    class(ElementBase1D), intent(in) :: elem1D_
    integer, intent(in) :: IS, IE, IA, IHALO
    real(RP), intent(inout) :: out_dg(elem1D_%Np,lcmesh%NeA)
    real(RP), intent(in) :: in_fv(IA)
    integer, intent(in) :: interp_ord

    integer :: i
    integer :: kelem

    real(RP) :: var(IA)
    real(RP) :: dqdx(IA), dqdxx(IA)
    !---------------------------------------
    
    var(:) = in_fv(:)
    do i=1, IHALO
      var(IS-i) = var(IE-i+1)
      var(IE+i) = var(IS+i-1)
    end do

    if ( interp_ord > 0 ) then
      do i=IS, IE
        dqdx(i) = 0.25_RP * ( var(i+1) - var(i-1) )
      end do
    end if
    if ( interp_ord > 1 ) then
      do i=IS, IE
        dqdxx(i) = 0.25_RP * ( var(i+1) - 2.0_RP * var(i) + var(i-1) )
      end do
    end if

    select case(interp_ord)
    case (0)
      do kelem=lcmesh%NeS, lcmesh%NeE
        i = IHALO + kelem
        out_dg(:,kelem) = in_fv(i)
      end do
    case (1)
      do kelem=lcmesh%NeS, lcmesh%NeE
        i = IHALO + kelem
        out_dg(:,kelem) = in_fv(i) &
          + dqdx(i) * elem1D_%x1(:)
      end do
    case (2)
      do kelem=lcmesh%NeS, lcmesh%NeE
        i = IHALO + kelem
        out_dg(:,kelem) = in_fv(i) &
          + dqdx(i) * elem1D_%x1(:) &
          + 0.5_RP * dqdxx(i) * ( elem1D_%x1(:)**2 - 1.0_RP / 3.0_RP )
      end do
    end select
    return
  end subroutine interp_FVtoDG_1D

end module mod_advect1d_fvm_numerror