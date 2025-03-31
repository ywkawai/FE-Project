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
  use scale_meshfield_analysis_numerror, only: &
    MeshFieldAnalysisNumerror1D

  use mod_dft, only: &
    DFT_dft1d, DFT_idft1d, DFT_interp1d
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  public :: advdiff1d_numerror_Init
  public :: advdiff1d_numerror_eval
  public :: advdiff1d_numerror_Final

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  type(MeshFieldAnalysisNumerror1D) :: numerror_analysis
  integer ::  numerror_vid(1)

  integer :: qexact_nmax
  integer :: qexact_mmax
  real(RP), allocatable :: x_qexact(:)
  complex(RP), allocatable :: s_qini(:)
  integer, allocatable :: s_mm(:)

contains
  !> Initialization
  subroutine advdiff1d_numerror_Init( elem )
    use scale_prc, only: &
       PRC_abort    
    implicit none
    class(LineElement), intent(in) :: elem

    integer :: ierr

    integer :: PolyOrderErrorCheck
    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    integer :: LOG_STEP_INTERVAL

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
    call numerror_analysis%Init( PolyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, elem )
    call numerror_analysis%Regist( "q", "1", numerror_vid(1) )

    !--

    allocate(s_qini(-qexact_mmax/2:qexact_mmax/2), s_mm(-qexact_mmax/2:qexact_mmax/2))
    allocate(x_qexact(0:qexact_nmax))
    do domid=0, qexact_nmax
      x_qexact(domid) = real(domid,kind=RP) / real(qexact_nmax+1,kind=RP)
    end do
    do domid=-qexact_mmax/2, qexact_mmax/2
      s_mm(domid) = domid
    end do

    return
  end subroutine advdiff1d_numerror_Init

  !> Evaluate numerical errors
  subroutine advdiff1d_numerror_eval( qtrc_exact,                  & ! (inout)
      qtrc, istep, tsec, ADV_VEL, DIFF_COEF, InitShapeName, InitShapeParams, & ! (in)
      mesh, elem                                                             ) ! (in)

    use scale_polynominal, only: Polynominal_genLegendrePoly
    use mod_fieldutil, only: &
      get_upwind_pos1d => fieldutil_get_upwind_pos1d,        &
      get_profile1d_tracer => fieldutil_get_profile1d_tracer 

    implicit none
    class(MeshField1D), intent(inout) :: qtrc_exact
    class(MeshField1D), intent(in) :: qtrc
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    real(RP), intent(in) :: ADV_VEL
    real(RP), intent(in) :: DIFF_COEF
    character(*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(2)
    class(MeshLineDom1D), intent(in) :: mesh
    class(ElementBase1D), intent(in) :: elem

    real(RP) :: dom_xmin, dom_xmax
    real(RP) :: qini_exact(0:qexact_nmax)
    !------------------------------------------------------------------------

    dom_xmin = mesh%xmin_gl
    dom_xmax = mesh%xmax_gl

    if (tsec == 0.0_RP) then
      call get_profile1d_tracer( qini_exact,                       & ! (out)
      InitShapeName, x_qexact(:), InitShapeParams, qexact_nmax+1 )   ! (in)   

      call DFT_dft1d( s_qini, qini_exact, qexact_nmax, qexact_mmax )
      LOG_INFO("CHECK_spectral real",*) real(s_qini, kind=RP)
      LOG_INFO("CHECK_spectral imag",*) imag(s_qini)  
    end if

    call numerror_analysis%Evaluate( istep, tsec, mesh, set_data_lc ) 

    return

  contains
!OCL SERIAL
    subroutine set_data_lc( this, q, qexact, qexact_intrp, lcmesh, elem1D, intrp_epos )
      use scale_localmeshfield_base, only: LocalMeshFieldBase
      implicit none
      class(MeshFieldAnalysisNumerror1D), intent(in) :: this
      class(LocalMesh1D), intent(in) :: lcmesh
      class(ElementBase1D) :: elem1D
      real(RP), intent(out) :: q(elem1D%Np,lcmesh%Ne,this%var_num)
      real(RP), intent(out) :: qexact(elem1D%Np,lcmesh%Ne,this%var_num)
      real(RP), intent(out) :: qexact_intrp(this%intrp_np,lcmesh%Ne,this%var_num)
      real(RP), intent(in) :: intrp_epos(this%intrp_np,this%ndim)

      integer :: n
      integer :: ke

      integer :: vid

      real(RP) :: vx(2)
      real(RP) :: x_uwind(elem%Np)
      real(RP) :: x_uwind_intrp(this%intrp_np)
      real(RP) :: pos_intrp(this%intrp_np)
      complex(RP) :: s_qtrc(-qexact_mmax/2:qexact_mmax/2)
      !---------------------------------------------

      n = lcmesh%lcdomID

      !$omp parallel do private(vid, ke, vx, x_uwind, x_uwind_intrp, pos_intrp, s_qtrc)
      do ke=lcmesh%NeS, lcmesh%NeE

        vid = numerror_vid(1)

        x_uwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,ke,1), ADV_VEL, tsec, dom_xmin, dom_xmax)
        vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
        pos_intrp(:) = vx(1) + 0.5_RP*( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
        x_uwind_intrp(:) = get_upwind_pos1d(pos_intrp(:), ADV_VEL, tsec, dom_xmin, dom_xmax)

        s_qtrc(:) = s_qini(:) * exp(- DIFF_COEF * (2.0_RP * PI * s_mm(:))**2 * tsec)

        ! call get_profile1d_tracer( qtrc_exact%local(n)%val(:,ke),    & ! (out)
        !   InitShapeName, x_uwind, InitShapeParams, elem%Np           ) ! (in)
        call DFT_interp1d( qtrc_exact%local(n)%val(:,ke), &
          s_qtrc(:), x_uwind(:), qexact_mmax )

        ! call get_profile1d_tracer( qexact_intrp(:,ke,vid),                  & ! (out)
        !   InitShapeName, x_uwind_intrp, InitShapeParams(:), this%intrp_np )   ! (in)
        call DFT_interp1d( qexact_intrp(:,ke,vid), &
          s_qtrc(:), x_uwind_intrp(:), qexact_mmax )

        q(:,ke,vid) = qtrc%local(n)%val(:,ke)
        qexact(:,ke,vid) = qtrc_exact%local(n)%val(:,ke)
      end do
      
      return
    end subroutine set_data_lc
  end subroutine advdiff1d_numerror_eval

  !> Finalization
  subroutine advdiff1d_numerror_Final()
    implicit none
    !---------------------------------

    deallocate( x_qexact, s_qini, s_mm )
    call numerror_analysis%Final()

    return
  end subroutine advdiff1d_numerror_Final

end module mod_advdiff1d_numerror