#include "scalelib.h"
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
  use scale_meshfield_analysis_numerror, only: &
    MeshFieldAnalysisNumerror1D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  public :: advect1d_numerror_Init
  public :: advect1d_numerror_eval
  public :: advect1d_numerror_Final

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  type(MeshFieldAnalysisNumerror1D) :: numerror_analysis
  integer ::  numerror_vid(1)

contains
  !> Initialization
  subroutine advect1d_numerror_Init( elem )
    use scale_prc, only: &
       PRC_abort    
    implicit none
    class(LineElement), intent(in) :: elem

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
    call numerror_analysis%Init( polyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, elem )
    call numerror_analysis%Regist( "q", "1", numerror_vid(1) )

    return
  end subroutine advect1d_numerror_Init

  !> Evaluate numerical errors
  subroutine advect1d_numerror_eval( qtrc_exact,                  & ! (inout)
      qtrc, istep, tsec, ADV_VEL, InitShapeName, InitShapeParams, & ! (in)
      mesh, elem                                                  ) ! (in)

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
    character(*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(2)
    class(MeshLineDom1D), intent(in) :: mesh
    class(ElementBase1D), intent(in) :: elem

    real(RP) :: dom_xmin, dom_xmax

    !------------------------------------------------------------------------

    dom_xmin = mesh%xmin_gl
    dom_xmax = mesh%xmax_gl

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
      !---------------------------------------------

      n = lcmesh%lcdomID

      !$omp parallel do private(vid, ke, vx, x_uwind, x_uwind_intrp, pos_intrp)
      do ke=lcmesh%NeS, lcmesh%NeE

        vid = numerror_vid(1)

        x_uwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,ke,1), ADV_VEL, tsec, dom_xmin, dom_xmax)
        vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
        pos_intrp(:) = vx(1) + 0.5_RP*( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
        x_uwind_intrp(:) = get_upwind_pos1d(pos_intrp(:), ADV_VEL, tsec, dom_xmin, dom_xmax)

        call get_profile1d_tracer( qtrc_exact%local(n)%val(:,ke),    & ! (out)
          InitShapeName, x_uwind, InitShapeParams, elem%Np           ) ! (in)

        call get_profile1d_tracer( qexact_intrp(:,ke,vid),                & ! (out)
          InitShapeName, x_uwind_intrp, InitShapeParams(:), this%intrp_np )   ! (in)

        q(:,ke,vid) = qtrc%local(n)%val(:,ke)
        qexact(:,ke,vid) = qtrc_exact%local(n)%val(:,ke)
      end do
      
      return
    end subroutine set_data_lc
  end subroutine advect1d_numerror_eval

  !> Finalization
  subroutine advect1d_numerror_Final()
    implicit none
    !---------------------------------
    
    call numerror_analysis%Final()
    return
  end subroutine advect1d_numerror_Final

end module mod_advect1d_numerror