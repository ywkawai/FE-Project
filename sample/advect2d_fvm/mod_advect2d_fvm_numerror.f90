#include "scaleFElib.h"
module mod_advect2d_fvm_numerror
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
  public :: advect2d_fvm_numerror_Init
  public :: advect2d_fvm_numerror_eval
  public :: advect2d_fvm_numerror_Final

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  type(MeshRectDom2D) :: mesh_dg
  type(QuadrilateralElement) :: elem

  type(MeshFieldAnalysisNumerror2D) :: numerror_analysis
  integer ::  numerror_vid(1)

contains
  !> Initialization
  subroutine advect2d_fvm_numerror_Init( FX, FY, IS, IE, IA, JS, JE, JA )
    use scale_prc, only: &
      PRC_abort
    use scale_prc_cartesC, only: &
      PRC_NUM_X, PRC_NUM_Y, &
      PRC_PERIODIC_X, PRC_PERIODIC_Y
    implicit none
    integer, intent(in) :: IS, IE, IA
    integer, intent(in) :: JS, JE, JA
    real(RP), intent(in) :: FX(0:IA)
    real(RP), intent(in) :: FY(0:JA)

    integer :: polyOrderErrorCheck
    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    integer :: LOG_STEP_INTERVAL

    namelist / PARAM_ADVECT2D_FVM_NUMERROR / &
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
    read(IO_FID_CONF,nml=PARAM_ADVECT2D_FVM_NUMERROR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO('advect2d_fvm_numerror_init',*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR('advect2d_fvm_numerror_init',*) 'Not appropriate names in namelist PARAM_ADVECT2D_FVM_NUMERROR. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ADVECT2D_FVM_NUMERROR)

    !--
    call elem%Init( polyOrderErrorCheck, .false. )
    call mesh_dg%Init( IE-IS+1, JE-JS+1, FX(IS-1), FX(IE), FY(JS-1), FY(JE), PRC_PERIODIC_X, PRC_PERIODIC_Y, &
      elem, 1, PRC_NUM_X, PRC_NUM_Y )
    call mesh_dg%Generate()

    !--
    call numerror_analysis%Init( polyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, elem )
    call numerror_analysis%Regist( "q", "1", numerror_vid(1) )

    return
  end subroutine advect2d_fvm_numerror_Init

  !> Evaluate numerical errors
  subroutine advect2d_fvm_numerror_eval( qtrc_exact,                               & ! (inout)
    qtrc, istep, tsec, VelTypeName, VelTypeParams, InitShapeName, InitShapeParams, & ! (in)
    CX, FX, CY, FY, IS, IE, IA, IHALO, JS, JE, JA, JHALO  ) ! (in)

    use scale_polynominal, only: Polynominal_genLegendrePoly
    use mod_fieldutil, only: &
      get_upwind_pos2d => fieldutil_get_upwind_pos2d,        &
      get_profile2d_tracer => fieldutil_get_profile2d_tracer 

    implicit none
    integer, intent(in) :: IS, IE, IA, IHALO
    integer, intent(in) :: JS, JE, JA, JHALO
    real(RP), intent(out) :: qtrc_exact(IA,JA)
    real(RP), intent(in) :: qtrc(IA,JA)
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    character(*), intent(in) :: VelTypeName
    real(RP), intent(in) :: VelTypeParams(4)
    character(*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(2)
    real(RP), intent(in) :: CX(IA)
    real(RP), intent(in) :: FX(IA)
    real(RP), intent(in) :: CY(JA)
    real(RP), intent(in) :: FY(JA)

    real(RP) :: dom_xmin, dom_xmax
    real(RP) :: dom_ymin, dom_ymax
    !------------------------------------------------------------------------

    dom_xmin = mesh_dg%xmin_gl
    dom_xmax = mesh_dg%xmax_gl
    dom_ymin = mesh_dg%ymin_gl
    dom_ymax = mesh_dg%ymax_gl
    call numerror_analysis%Evaluate( istep, tsec, mesh_dg, set_data_lc ) 

    return

  contains
!OCL SERIAL
    subroutine set_data_lc( this, q, qexact, qexact_intrp, lcmesh, elem2D, intrp_epos )
      use scale_localmeshfield_base, only: LocalMeshFieldBase
      use scale_meshfield_fvm_util, only: MeshFieldFVMUtil_interp_FVtoDG
      implicit none
      class(MeshFieldAnalysisNumerror2D), intent(in) :: this
      class(LocalMesh2D), intent(in) :: lcmesh
      class(ElementBase2D) :: elem2D
      real(RP), intent(out) :: q(elem2D%Np,lcmesh%Ne,this%var_num)
      real(RP), intent(out) :: qexact(elem2D%Np,lcmesh%Ne,this%var_num)
      real(RP), intent(out) :: qexact_intrp(this%intrp_np,lcmesh%Ne,this%var_num)
      real(RP), intent(in) :: intrp_epos(this%intrp_np,this%ndim)

      integer :: n
      integer :: kelem, ke_x, ke_y

      integer :: vid

      real(RP) :: vx(4), vy(4)
      real(RP) :: x_uwind(elem%Np), y_vwind(elem%Np)
      real(RP) :: x_uwind_intrp(this%intrp_np), y_vwind_intrp(this%intrp_np)
      real(RP) :: pos_intrp(this%intrp_np,2)  
      real(RP) :: int_w(elem%Np)
      real(RP) :: qtrc_exact_dg(elem%Np)
      real(RP) :: qtrc_dg(elem%Np,lcmesh%NeA)
      !---------------------------------------------

      n = lcmesh%lcdomID

      call MeshFieldFVMUtil_interp_FVtoDG( qtrc_dg, &
        qtrc, lcmesh, elem2D, IS, IE, IA, IHALO, JS, JE, JA, JHALO, &
        1 )

      !$omp parallel do private(vid, ke_x, ke_y, kelem, vx, vy, x_uwind, y_vwind, x_uwind_intrp, y_vwind_intrp, pos_intrp, int_w, &
      !$omp qtrc_exact_dg ) collapse(2)
      do ke_y=1, lcmesh%NeY
      do ke_x=1, lcmesh%NeX
        kelem = ke_x + (ke_y-1)*lcmesh%NeX

        vid = numerror_vid(1)

        call get_upwind_pos2d( x_uwind, y_vwind, & !(out) 
          lcmesh%pos_en(:,kelem,1), lcmesh%pos_en(:,kelem,2), VelTypeName, VelTypeParams, tsec, & ! (in)
          dom_xmin, dom_xmax, dom_ymin, dom_ymax                                                ) ! (in)

        vx(:) = lcmesh%pos_ev(lcmesh%EToV(kelem,:),1)
        vy(:) = lcmesh%pos_ev(lcmesh%EToV(kelem,:),2)
        pos_intrp(:,1) = vx(1) + 0.5_RP*( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
        pos_intrp(:,2) = vy(1) + 0.5_RP*( intrp_epos(:,2) + 1.0_RP ) * ( vy(3) - vy(1) )
        call get_upwind_pos2d( x_uwind_intrp, y_vwind_intrp, & !(out) 
          pos_intrp(:,1), pos_intrp(:,2), VelTypeName, VelTypeParams, tsec, & ! (in)
          dom_xmin, dom_xmax, dom_ymin, dom_ymax                            ) ! (in)

        call get_profile2d_tracer( qtrc_exact_dg(:),                 & ! (out)
          InitShapeName, x_uwind, y_vwind, InitShapeParams, elem%Np  ) ! (in)

        call get_profile2d_tracer( qexact_intrp(:,kelem,vid),                         & ! (out)
          InitShapeName, x_uwind_intrp, y_vwind_intrp, InitShapeParams, this%intrp_np ) ! (in)

        int_w(:) = lcmesh%Gsqrt(:,kelem) * lcmesh%J(:,kelem) * elem%IntWeight_lgl(:)
        int_w(:) = int_w(:) / sum(int_w(:))
        qtrc_exact(IHALO+ke_x,JHALO+ke_y) = sum( int_w(:) * qtrc_exact_dg(:) )

        q(:,kelem,vid) = qtrc_dg(:,kelem)
        qexact(:,kelem,vid) = qtrc_exact_dg(:)
      end do
      end do
      
      return
    end subroutine set_data_lc
  end subroutine advect2d_fvm_numerror_eval

  !> Finalization
  subroutine advect2d_fvm_numerror_Final()
    implicit none
    !---------------------------------
    call mesh_dg%Final()
    call elem%Final()
    call numerror_analysis%Final()
    return
  end subroutine advect2d_fvm_numerror_Final
  
end module mod_advect2d_fvm_numerror