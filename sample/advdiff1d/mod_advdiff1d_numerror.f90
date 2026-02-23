#include "scalelib.h"
module mod_advdiff1d_numerror
  use scale_const, only: &
    PI  => CONST_PI  
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
  
  use mod_dft, only: &
    DFT_dft1d, DFT_idft1d, DFT_interp1d
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public types & variables
  !
  type, extends(MeshFieldAnalysisNumerrorInfoBase) :: AdvDiff1D_Numerror_Info
    real(RP) :: ADV_VEL
    real(RP) :: DIFF_COEF
    character(len=H_MID) :: InitShapeName
    real(RP) :: InitShapeParams(2)
    real(RP) :: dom_xmin
    real(RP) :: dom_xmax
    integer ::  numerror_vid(1)
    !-
    type(MeshField1D), pointer :: qtrc_exact
    type(MeshField1D), pointer :: qtrc

    integer :: qexact_nmax
    integer :: qexact_mmax
    real(RP), allocatable :: x_qexact(:)
    complex(RP), allocatable :: s_qini(:)
    integer, allocatable :: s_mm(:)
  end type AdvDiff1D_Numerror_Info

  type, public :: AdvDiff1DNumErrorAnalysis
    type(AdvDiff1D_Numerror_Info) :: info
    type(MeshFieldAnalysisNumerror1D) :: numerror_analysis
  contains
    procedure :: Init => advdiff1d_numerror_Init
    procedure :: Eval => advdiff1d_numerror_eval
    procedure :: Final => advdiff1d_numerror_Final
  end type AdvDiff1DNumErrorAnalysis

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !

contains
  !> Initialization
  subroutine advdiff1d_numerror_Init(this,              &
    ADV_VEL, DIFF_COEF, InitShapeName, InitShapeParams, &
    mesh, elem )
    use scale_prc, only: &
       PRC_abort    
    implicit none
    class(AdvDiff1DNumErrorAnalysis), intent(inout) :: this
    real(RP), intent(in) :: ADV_VEL
    real(RP), intent(in) :: DIFF_COEF
    character(*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(2)
    type(MeshLineDom1D), intent(in) :: mesh
    type(LineElement), intent(in) :: elem

    integer :: ierr

    integer :: PolyOrderErrorCheck
    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    integer :: LOG_STEP_INTERVAL

    integer :: qexact_nmax
    integer :: qexact_mmax
    namelist / PARAM_ADVDIFF1D_NUMERROR / &
      PolyOrderErrorCheck,    &
      LOG_OUT_BASENAME,       &
      LOG_STEP_INTERVAL,      &
      qexact_nmax, qexact_mmax

    integer :: domid
    !----------------------------------------------

    LOG_NEWLINE
    LOG_INFO("advdiff1d_numerror_init",*) 'Setup'

    PolyOrderErrorCheck = 6
    LOG_STEP_INTERVAL   = 5
    qexact_nmax         = 64
    qexact_mmax         = 8

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ADVDIFF1D_NUMERROR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO('advdiff1d_numerror_init',*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR('advdiff1d_numerror_init',*) 'Not appropriate names in namelist PARAM_ADVDIFF1D_NUMERROR. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ADVDIFF1D_NUMERROR)

    !--
    call this%numerror_analysis%Init( &
      PolyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, &
      mesh, elem, set_data_lc, this%info                        )
    call this%numerror_analysis%Regist( "q", "1", this%info%numerror_vid(1) )

    this%info%ADV_VEL = ADV_VEL
    this%info%DIFF_COEF = DIFF_COEF
    this%info%InitShapeName = InitShapeName
    this%info%InitShapeParams = InitShapeParams
    this%info%dom_xmin = mesh%xmin_gl
    this%info%dom_xmax = mesh%xmax_gl
    this%info%qexact_nmax = qexact_nmax
    this%info%qexact_mmax = qexact_mmax

    !--

    allocate( this%info%s_qini(-qexact_mmax/2:qexact_mmax/2) )
    allocate( this%info%s_mm(-qexact_mmax/2:qexact_mmax/2) )
    allocate( this%info%x_qexact(0:qexact_nmax) )
    do domid=0, qexact_nmax
      this%info%x_qexact(domid) = real(domid,kind=RP) / real(qexact_nmax+1,kind=RP)
    end do
    do domid=-qexact_mmax/2, qexact_mmax/2
      this%info%s_mm(domid) = domid
    end do
    return
  end subroutine advdiff1d_numerror_Init

  !> Evaluate numerical errors
  subroutine advdiff1d_numerror_eval( this, qtrc_exact, & ! (inout)
      qtrc, istep, tsec                                 ) ! (in)

    use mod_fieldutil, only: &
      get_profile1d_tracer => fieldutil_get_profile1d_tracer
    implicit none
    class(AdvDiff1DNumErrorAnalysis), intent(inout) :: this
    class(MeshField1D), intent(inout), target :: qtrc_exact
    class(MeshField1D), intent(in), target :: qtrc
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec

    real(RP) :: qini_exact(0:this%info%qexact_nmax)
    !------------------------------------------------------------------------

    if (tsec == 0.0_RP) then
      call get_profile1d_tracer( qini_exact,                       & ! (out)
      this%info%InitShapeName, this%info%x_qexact(:), this%info%InitShapeParams, this%info%qexact_nmax+1 )   ! (in)   

      call DFT_dft1d( this%info%s_qini, qini_exact, this%info%qexact_nmax, this%info%qexact_mmax )
      LOG_INFO("CHECK_spectral real",*) real(this%info%s_qini, kind=RP)
      LOG_INFO("CHECK_spectral imag",*) imag(this%info%s_qini)  
    end if

    this%info%qtrc_exact => qtrc_exact
    this%info%qtrc => qtrc
    call this%numerror_analysis%Evaluate( istep, tsec ) 
    return
  end subroutine advdiff1d_numerror_eval

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
    integer :: ke

    integer :: vid

    real(RP) :: vx(2)
    real(RP) :: x_uwind(elem1D%Np)
    real(RP) :: x_uwind_intrp(analysis%intrp_np)
    real(RP) :: pos_intrp(analysis%intrp_np)
    complex(RP), allocatable :: s_qtrc(:)


    class(AdvDiff1D_Numerror_Info), pointer :: info 
    class(MeshFieldAnalysisNumerrorInfoBase), pointer :: info_base
    !---------------------------------------------

    select type(info_base => analysis%info)
    class is (AdvDiff1D_Numerror_Info)
      info => info_base
    end select  

    n = lcmesh%lcdomID

    allocate( s_qtrc(-info%qexact_mmax/2:info%qexact_mmax/2) )

    !$omp parallel do private(vid, ke, vx, x_uwind, x_uwind_intrp, pos_intrp, s_qtrc)
    do ke=lcmesh%NeS, lcmesh%NeE

      vid = info%numerror_vid(1)

      x_uwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,ke,1), info%ADV_VEL, tsec, info%dom_xmin, info%dom_xmax)
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      pos_intrp(:) = vx(1) + 0.5_RP*( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
      x_uwind_intrp(:) = get_upwind_pos1d(pos_intrp(:), info%ADV_VEL, tsec, info%dom_xmin, info%dom_xmax)

      s_qtrc(:) = info%s_qini(:) * exp(- info%DIFF_COEF * (2.0_RP * PI * info%s_mm(:))**2 * tsec)

      call DFT_interp1d( info%qtrc_exact%local(n)%val(:,ke), &
        s_qtrc(:), x_uwind(:), info%qexact_mmax )

      call DFT_interp1d( qexact_intrp(:,ke,vid), &
        s_qtrc(:), x_uwind_intrp(:), info%qexact_mmax )

      q(:,ke,vid) = info%qtrc%local(n)%val(:,ke)
      qexact(:,ke,vid) = info%qtrc_exact%local(n)%val(:,ke)
    end do
    return
  end subroutine set_data_lc

  !> Finalization
  subroutine advdiff1d_numerror_Final( this )
    implicit none
    class(AdvDiff1DNumErrorAnalysis), intent(inout) :: this
    !---------------------------------
    deallocate( this%info%x_qexact, this%info%s_qini, this%info%s_mm )
    call this%numerror_analysis%Final()
    return
  end subroutine advdiff1d_numerror_Final
end module mod_advdiff1d_numerror