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
  use scale_meshfield_analysis_numerror, only: &
    MeshFieldAnalysisNumerror1D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  public :: advect1d_fvm_numerror_Init
  public :: advect1d_fvm_numerror_eval
  public :: advect1d_fvm_numerror_Final

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  type(MeshLineDom1D) :: mesh_dg
  type(LineElement) :: elem

  type(MeshFieldAnalysisNumerror1D) :: numerror_analysis
  integer ::  numerror_vid(1)

contains
  !> Initialization
  subroutine advect1d_fvm_numerror_Init( FX, IS, IE, IA )
    use scale_prc, only: &
      PRC_abort    
    implicit none
    integer, intent(in) :: IS, IE, IA
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
    call elem%Init( polyOrderErrorCheck, .false. )
    call mesh_dg%Init( IE-IS+1, FX(IS-1), FX(IE), elem, 1 )
    call mesh_dg%Generate()

    !--
    call numerror_analysis%Init( polyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, elem )
    call numerror_analysis%Regist( "q", "1", numerror_vid(1) )

    return
  end subroutine advect1d_fvm_numerror_Init

  !> Evaluate numerical errors
  subroutine advect1d_fvm_numerror_eval( qtrc_exact,              & ! (inout)
      qtrc, istep, tsec, ADV_VEL, InitShapeName, InitShapeParams, &
      CX, FX, IS, IE, IA, IHALO ) ! (in)

    use scale_polynominal, only: Polynominal_genLegendrePoly
    use mod_fieldutil, only: &
      get_upwind_pos1d => fieldutil_get_upwind_pos1d,        &
      get_profile1d_tracer => fieldutil_get_profile1d_tracer 

    implicit none
    integer, intent(in) :: IS, IE, IA, IHALO
    real(RP), intent(out) :: qtrc_exact(IA)
    real(RP), intent(in) :: qtrc(IA)
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    real(RP), intent(in) :: ADV_VEL
    character(*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(2)
    real(RP), intent(in) :: CX(IA)
    real(RP), intent(in) :: FX(IA)

    real(RP) :: dom_xmin, dom_xmax

    !------------------------------------------------------------------------

    dom_xmin = mesh_dg%xmin_gl
    dom_xmax = mesh_dg%xmax_gl
    call numerror_analysis%Evaluate( istep, tsec, mesh_dg, set_data_lc ) 

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
      integer :: kelem

      integer :: vid

      real(RP) :: vx(2)
      real(RP) :: x_uwind(elem%Np)
      real(RP) :: x_uwind_intrp(this%intrp_np)
      real(RP) :: pos_intrp(this%intrp_np)  
      real(RP) :: int_w(elem%Np)
      real(RP) :: qtrc_exact_dg(elem%Np)
      real(RP) :: qtrc_dg(elem%Np,lcmesh%NeA)
      !---------------------------------------------

      n = lcmesh%lcdomID

      call interp_FVtoDG_1D( qtrc_dg, &
        qtrc, lcmesh, lcmesh%refElem1D, IS, IE, IA, IHALO, 2 )

      !$omp parallel do private(vid, kelem, vx, x_uwind, x_uwind_intrp, pos_intrp, int_w, &
      !$omp qtrc_exact_dg )
      do kelem=lcmesh%NeS, lcmesh%NeE

        vid = numerror_vid(1)

        x_uwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,kelem,1), ADV_VEL, tsec, dom_xmin, dom_xmax)
        vx(:) = lcmesh%pos_ev(lcmesh%EToV(kelem,:),1)
        pos_intrp(:) = vx(1) + 0.5_RP*( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
        x_uwind_intrp(:) = get_upwind_pos1d(pos_intrp(:), ADV_VEL, tsec, dom_xmin, dom_xmax)

        call get_profile1d_tracer( qtrc_exact_dg(:),       & ! (out)
          InitShapeName, x_uwind, InitShapeParams, elem%Np ) ! (in)

        call get_profile1d_tracer( qexact_intrp(:,kelem,vid),                & ! (out)
          InitShapeName, x_uwind_intrp, InitShapeParams(:), this%intrp_np )   ! (in)

        int_w(:) = lcmesh%Gsqrt(:,kelem) * lcmesh%J(:,kelem) * elem%IntWeight_lgl(:)
        int_w(:) = int_w(:) / sum(int_w(:))
        qtrc_exact(IHALO+kelem) = sum( int_w(:) * qtrc_exact_dg(:) )

        q(:,kelem,vid) = qtrc_dg(:,kelem)
        qexact(:,kelem,vid) = qtrc_exact_dg(:)
      end do
      
      return
    end subroutine set_data_lc
  end subroutine advect1d_fvm_numerror_eval

  !> Finalization
  subroutine advect1d_fvm_numerror_Final()
    implicit none
    !---------------------------------
    call mesh_dg%Final()
    call elem%Final()
    call numerror_analysis%Final()
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