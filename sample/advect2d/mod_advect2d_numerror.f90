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
  use scale_meshfield_analysis_numerror, only: &
    MeshFieldAnalysisNumerror2D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  public :: advect2d_numerror_Init
  public :: advect2d_numerror_eval
  public :: advect2d_numerror_Final

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  type(MeshFieldAnalysisNumerror2D) :: numerror_analysis
  integer ::  numerror_vid(1)

contains
  !> Initialization
  subroutine advect2d_numerror_Init( elem )
    use scale_prc, only: &
       PRC_abort    
    implicit none
    class(QuadrilateralElement), intent(in) :: elem

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
    call numerror_analysis%Init( polyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, elem )
    call numerror_analysis%Regist( "q", "1", numerror_vid(1) )

    return
  end subroutine advect2d_numerror_Init

  !> Evaluate numerical errors
  subroutine advect2d_numerror_eval( qtrc_exact,                                     & ! (inout)
      qtrc, istep, tsec, VelTypeName, VelTypeParams, InitShapeName, InitShapeParams, & ! (in)
      mesh, elem                                                                     ) ! (in)

    use scale_polynominal, only: Polynominal_genLegendrePoly
    use mod_fieldutil, only: &
      get_upwind_pos2d => fieldutil_get_upwind_pos2d,        &
      get_profile2d_tracer => fieldutil_get_profile2d_tracer 

    implicit none
    class(MeshField2D), intent(inout) :: qtrc_exact
    class(MeshField2D), intent(in) :: qtrc
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    character(*), intent(in) :: VelTypeName
    real(RP), intent(in) :: VelTypeParams(4)
    character(*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(2)
    class(MeshRectDom2D), intent(in) :: mesh
    class(ElementBase2D), intent(in) :: elem

    real(RP) :: dom_xmin, dom_xmax
    real(RP) :: dom_ymin, dom_ymax
    !------------------------------------------------------------------------

    dom_xmin = mesh%xmin_gl
    dom_xmax = mesh%xmax_gl
    dom_ymin = mesh%ymin_gl
    dom_ymax = mesh%ymax_gl

    call numerror_analysis%Evaluate( istep, tsec, mesh, set_data_lc ) 

    return

  contains
!OCL SERIAL
    subroutine set_data_lc( this, q, qexact, qexact_intrp, lcmesh, elem2D, intrp_epos )
      use scale_localmeshfield_base, only: LocalMeshFieldBase
      implicit none
      class(MeshFieldAnalysisNumerror2D), intent(in) :: this
      class(LocalMesh2D), intent(in) :: lcmesh
      class(ElementBase2D) :: elem2D
      real(RP), intent(out) :: q(elem2D%Np,lcmesh%Ne,this%var_num)
      real(RP), intent(out) :: qexact(elem2D%Np,lcmesh%Ne,this%var_num)
      real(RP), intent(out) :: qexact_intrp(this%intrp_np,lcmesh%Ne,this%var_num)
      real(RP), intent(in) :: intrp_epos(this%intrp_np,this%ndim)

      integer :: n
      integer :: ke

      integer :: vid

      real(RP) :: vx(4), vy(4)
      real(RP) :: x_uwind(elem%Np), y_vwind(elem%Np)
      real(RP) :: x_uwind_intrp(this%intrp_np), y_vwind_intrp(this%intrp_np)
      real(RP) :: pos_intrp(this%intrp_np,2)  
      !---------------------------------------------

      n = lcmesh%lcdomID

      !$omp parallel do private(vid, ke, vx, vy, x_uwind, y_vwind, x_uwind_intrp, y_vwind_intrp, pos_intrp)
      do ke=lcmesh%NeS, lcmesh%NeE

        vid = numerror_vid(1)

        call get_upwind_pos2d( x_uwind, y_vwind, & !(out) 
          lcmesh%pos_en(:,ke,1), lcmesh%pos_en(:,ke,2), VelTypeName, VelTypeParams, tsec, & ! (in)
          dom_xmin, dom_xmax, dom_ymin, dom_ymax                                          ) ! (in)

        vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
        vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
        pos_intrp(:,1) = vx(1) + 0.5_RP*( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
        pos_intrp(:,2) = vy(1) + 0.5_RP*( intrp_epos(:,2) + 1.0_RP ) * ( vy(3) - vy(1) )
        call get_upwind_pos2d( x_uwind_intrp, y_vwind_intrp, & !(out) 
          pos_intrp(:,1), pos_intrp(:,2), VelTypeName, VelTypeParams, tsec, & ! (in)
          dom_xmin, dom_xmax, dom_ymin, dom_ymax                            ) ! (in)

        call get_profile2d_tracer( qtrc_exact%local(n)%val(:,ke),    & ! (out)
          InitShapeName, x_uwind, y_vwind, InitShapeParams, elem%Np  ) ! (in)

        call get_profile2d_tracer( qexact_intrp(:,ke,vid),                               & ! (out)
          InitShapeName, x_uwind_intrp, y_vwind_intrp, InitShapeParams, this%intrp_np ) ! (in)

        q(:,ke,vid) = qtrc%local(n)%val(:,ke)
        qexact(:,ke,vid) = qtrc_exact%local(n)%val(:,ke)
      end do
      
      return
    end subroutine set_data_lc
  end subroutine advect2d_numerror_eval

  !> Finalization
  subroutine advect2d_numerror_Final()
    implicit none
    !---------------------------------
    
    call numerror_analysis%Final()
    return
  end subroutine advect2d_numerror_Final

end module mod_advect2d_numerror