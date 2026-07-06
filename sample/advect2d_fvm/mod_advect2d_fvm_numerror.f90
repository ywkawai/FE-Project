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
    type(MeshRectDom2D) :: mesh_dg
    type(QuadrilateralElement) :: elem

    !-
    real(RP), pointer :: qtrc_exact(:,:,:)
    real(RP), pointer :: qtrc(:,:,:)

    !-
    integer :: FV_IS, FV_IE, FV_IA, FV_IHALO
    integer :: FV_JS, FV_JE, FV_JA, FV_JHALO
    integer :: FV_KS
  contains
    procedure :: Init => advect2d_numerror_info_Init
    procedure :: Final => advect2d_numerror_info_Final
  end type Advect2D_Numerror_Info

  type, public :: Advect2DNumErrorAnalysis
    type(Advect2D_Numerror_Info) :: info
    type(MeshFieldAnalysisNumerror2D) :: numerror_analysis
  contains
    procedure :: Init => advect2d_fvm_numerror_Init
    procedure :: Eval => advect2d_fvm_numerror_Eval
    procedure :: Final => advect2d_fvm_numerror_Final
  end type Advect2DNumErrorAnalysis

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
contains
  !> Initialization
  subroutine advect2d_fvm_numerror_Init( this, &
    VelTypeName, VelTypeParams, InitShapeName, InitShapeParams, &
    FX, FY, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS )
    use scale_prc, only: &
      PRC_abort
    implicit none
    class(Advect2DNumErrorAnalysis), intent(inout) :: this
    character(len=*), intent(in) :: VelTypeName
    real(RP), intent(in) :: VelTypeParams(4)
    character(len=*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(4)
    integer, intent(in) :: IS, IE, IA, IHALO
    integer, intent(in) :: JS, JE, JA, JHALO
    integer, intent(in) :: KS
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
    LOG_INFO("advect2d_fvm_numerror_init",*) 'Setup'

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
    call this%info%Init( polyOrderErrorCheck, IS, IE, IA, IHALO, FX, JS, JE, JA, JHALO, FY, KS, &
      VelTypeName, VelTypeParams, InitShapeName, InitShapeParams )

    !--
    call this%numerror_analysis%Init( &
      polyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, &
      this%info%mesh_dg, this%info%elem, set_data_lc, this%info )
    call this%numerror_analysis%Regist( "q", "1", this%info%numerror_vid(1) )
    return
  end subroutine advect2d_fvm_numerror_Init

  !> Evaluate numerical errors
  subroutine advect2d_fvm_numerror_eval( this, &
    qtrc_exact,                   & ! (inout)
    qtrc, istep, tsec, KA, IA, JA ) ! (in)
    implicit none
    class(Advect2DNumErrorAnalysis), intent(inout) :: this
    integer, intent(in) :: KA, IA, JA
    real(RP), intent(out), target :: qtrc_exact(KA,IA,JA)
    real(RP), intent(in), target :: qtrc(KA,IA,JA)
    integer, intent(in) :: istep
    real(RP), intent(in) :: tsec
    !------------------------------------------------------------------------
    this%info%qtrc_exact => qtrc_exact
    this%info%qtrc => qtrc
    call this%numerror_analysis%Evaluate( istep, tsec ) 
    return
  end subroutine advect2d_fvm_numerror_eval

!OCL SERIAL
  subroutine set_data_lc( analysis, q, qexact, qexact_intrp, lcmesh, elem2D, intrp_epos, tsec )
    use mod_fieldutil, only: &
      get_upwind_pos2d => fieldutil_get_upwind_pos2d,        &
      get_profile2d_tracer => fieldutil_get_profile2d_tracer 
    use scale_meshfield_fvm_util, only: MeshFieldFVMUtil_interp_FVtoDG
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
    integer :: kelem, ke_x, ke_y

    integer :: vid

    real(RP) :: vx(4), vy(4)
    real(RP) :: x_uwind(elem2D%Np), y_vwind(elem2D%Np)
    real(RP) :: x_uwind_intrp(analysis%intrp_np), y_vwind_intrp(analysis%intrp_np)
    real(RP) :: pos_intrp(analysis%intrp_np,2)  
    real(RP) :: int_w(elem2D%Np)
    real(RP) :: qtrc_exact_dg(elem2D%Np)
    real(RP) :: qtrc_dg(elem2D%Np,lcmesh%NeA)

    class(Advect2D_Numerror_Info), pointer :: info 
    class(MeshFieldAnalysisNumerrorInfoBase), pointer :: info_base

    integer :: KS
    !---------------------------------------------

    select type(info_base => analysis%info)
    class is (Advect2D_Numerror_Info)
      info => info_base
    end select    

    n = lcmesh%lcdomID
    !$acc update host( info%qtrc )

    KS = info%FV_KS
    
    call MeshFieldFVMUtil_interp_FVtoDG( qtrc_dg, &
      info%qtrc(KS,:,:), lcmesh, elem2D,                 &
      info%FV_IS, info%FV_IE, info%FV_IA, info%FV_IHALO, &
      info%FV_JS, info%FV_JE, info%FV_JA, info%FV_JHALO, &
      1 )

    !$omp parallel do private(vid, ke_x, ke_y, kelem, vx, vy, x_uwind, y_vwind, x_uwind_intrp, y_vwind_intrp, pos_intrp, int_w, &
    !$omp qtrc_exact_dg ) collapse(2)
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
      kelem = ke_x + (ke_y-1)*lcmesh%NeX

      vid = info%numerror_vid(1)

      call get_upwind_pos2d( x_uwind, y_vwind, & !(out) 
        lcmesh%pos_en(:,kelem,1), lcmesh%pos_en(:,kelem,2), info%VelTypeName, info%VelTypeParams, tsec, & ! (in)
        info%dom_xmin, info%dom_xmax, info%dom_ymin, info%dom_ymax                                      ) ! (in)

      vx(:) = lcmesh%pos_ev(lcmesh%EToV(kelem,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(kelem,:),2)
      pos_intrp(:,1) = vx(1) + 0.5_RP*( intrp_epos(:,1) + 1.0_RP ) * ( vx(2) - vx(1) )
      pos_intrp(:,2) = vy(1) + 0.5_RP*( intrp_epos(:,2) + 1.0_RP ) * ( vy(3) - vy(1) )
      call get_upwind_pos2d( x_uwind_intrp, y_vwind_intrp, & !(out) 
        pos_intrp(:,1), pos_intrp(:,2), info%VelTypeName, info%VelTypeParams, tsec, & ! (in)
        info%dom_xmin, info%dom_xmax, info%dom_ymin, info%dom_ymax                  ) ! (in)

      call get_profile2d_tracer( qtrc_exact_dg(:),                             & ! (out)
        info%InitShapeName, x_uwind, y_vwind, info%InitShapeParams, elem2D%Np  ) ! (in)

      call get_profile2d_tracer( qexact_intrp(:,kelem,vid),                                       & ! (out)
        info%InitShapeName, x_uwind_intrp, y_vwind_intrp, info%InitShapeParams, analysis%intrp_np ) ! (in)

      int_w(:) = lcmesh%Gsqrt(:,kelem) * lcmesh%J(:,kelem) * elem2D%IntWeight_lgl(:)
      int_w(:) = int_w(:) / sum(int_w(:))
      info%qtrc_exact(KS,info%FV_IHALO+ke_x,info%FV_JHALO+ke_y) = sum( int_w(:) * qtrc_exact_dg(:) )

      q(:,kelem,vid) = qtrc_dg(:,kelem)
      qexact(:,kelem,vid) = qtrc_exact_dg(:)
    end do
    end do
    
    !$acc update device( info%qtrc_exact )
    return
  end subroutine set_data_lc

  !> Finalization
  subroutine advect2d_fvm_numerror_Final( this )
    implicit none
    class(Advect2DNumErrorAnalysis), intent(inout) :: this
    !---------------------------------
    call this%info%Final()
    call this%numerror_analysis%Final()
    return
  end subroutine advect2d_fvm_numerror_Final
  
!- private --------------------

  subroutine advect2d_numerror_info_Init( this, &
    polyOrder, IS, IE, IA, IHALO, FX,  JS, JE, JA, JHALO, FY, KS, &
    VelTypeName, VelTypeParams, InitShapeName, InitShapeParams )
    use scale_prc_cartesC, only: &
      PRC_NUM_X, PRC_NUM_Y, &
      PRC_PERIODIC_X, PRC_PERIODIC_Y    
    implicit none
    class(Advect2D_Numerror_Info), intent(inout) :: this
    integer, intent(in) :: polyOrder
    integer, intent(in) :: IS, IE, IA, IHALO
    real(RP), intent(in) :: FX(0:IA)
    integer, intent(in) :: JS, JE, JA, JHALO
    real(RP), intent(in) :: FY(0:JA)
    integer, intent(in) :: KS
    character(len=*), intent(in) :: VelTypeName
    real(RP), intent(in) :: VelTypeParams(4)
    character(len=*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(4)
    !----------------------------------------------

    call this%elem%Init( polyOrder, .false. )
    call this%mesh_dg%Init( IE-IS+1, JE-JS+1, FX(IS-1), FX(IE), FY(JS-1), FY(JE), PRC_PERIODIC_X, PRC_PERIODIC_Y, &
      this%elem, 1, PRC_NUM_X, PRC_NUM_Y )
    call this%mesh_dg%Generate()

    this%VelTypeName = VelTypeName
    this%VelTypeParams = VelTypeParams
    this%InitShapeName = InitShapeName
    this%InitShapeParams = InitShapeParams
    this%dom_xmin = this%mesh_dg%xmin_gl
    this%dom_xmax = this%mesh_dg%xmax_gl
    this%dom_ymin = this%mesh_dg%ymin_gl
    this%dom_ymax = this%mesh_dg%ymax_gl
    this%FV_IS = IS
    this%FV_IE = IE
    this%FV_IA = IA
    this%FV_IHALO = IHALO
    this%FV_JS = JS
    this%FV_JE = JE
    this%FV_JA = JA
    this%FV_JHALO = JHALO
    this%FV_KS = KS
    return
  end subroutine advect2d_numerror_info_Init

  subroutine advect2d_numerror_info_Final( this )
    implicit none
    class(Advect2D_Numerror_Info), intent(inout) :: this
    !---------------------------------
    call this%mesh_dg%Final()
    call this%elem%Final()
    return
  end subroutine advect2d_numerror_info_Final

end module mod_advect2d_fvm_numerror