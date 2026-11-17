!-------------------------------------------------------------------------------
!> module Multigrid smoother for the Poisson equation in a spherical shell
!!
!! @par Description
!!      A smoother of the multigrid method for the 3D Poisson equation
!!      on a cubed-sphere mesh with the shallow approximation (see formulation3d.pdf).
!!      The equation multiplied by Gsqrt,
!!
!!        d/dxi^i ( Gsqrt G^{ij} dq/dxi^j ) = Gsqrt f,
!!
!!      is discretized with a symmetric interior penalty (SIP) type nodal DG method,
!!      using the mass matrix of the reference element (without Gsqrt).
!!      The right-hand side given to the smoother must be already multiplied by Gsqrt.
!!
!!      Boundary conditions at the bottom and top boundaries (Neumann or Dirichlet) are imposed
!!      through ghost values in the numerical fluxes (halo values are not used there).
!!
!!      Smoother: element-block Gauss-Seidel with four colors
!!      (2 * panel group + red-black color in each tile).
!!      The element-local problem is solved by a preconditioned conjugate gradient method
!!      applied to the system multiplied by the mass matrix.
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_poisson3dGlobal_smoother
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof

  use scale_sparsemat, only: &
    SparseMat, SparseMat_matmul

  use scale_element_base, only: ElementBase3D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfield_base, only: MeshField3D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_base, only: MeshFieldCommBase

  use scale_multigrid_smoother_base, only: &
    MGSmootherBase3D, MGSmoother_PRE_ID, MGSmoother_POST_ID
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type, public, extends(MGSmootherBase3D) :: MGSmoother_Poisson3D
    integer :: num_internal_itr_max                !< Maximum number of internal iterations associated with local element solve
    real(RP) :: residual_threshold_l2_internal_itr !< Threshold for each internal iteration associated with local element solve
    real(RP) :: penalty_coef                       !< Coefficient C_tau in the penalty parameter (p-MG levels)
    real(RP) :: penalty_coef_hMG                   !< Coefficient C_tau in the penalty parameter (h-MG levels)
    integer :: bc_type_id                          !< Type of boundary conditions at the bottom and top boundaries
    real(RP) :: dom_zmin                           !< Bottom of the computational domain
    real(RP) :: dom_zmax                           !< Top of the computational domain
  contains
    procedure :: Init => MGSmoother_Poisson3D_Init
    procedure :: Final => MGSmoother_Poisson3D_Final
    procedure :: Advance_itr_1step => MGSmoother_Poisson3D_advance_itr_1step
    procedure, private :: advance_itr_1step_color => Poisson3dGlobal_smoother_advance_itr_1step_color
    procedure, private :: cal_q_lc_elem_PCG => cal_q_lc_elem_PCG
  end type MGSmoother_Poisson3D

  !> Number of auxiliary variables:
  !! one scalar (vertical component of the element-wise gradient) and
  !! one horizontal vector (contravariant components of the element-wise gradient)
  integer, parameter, public :: MGSmoother_Poisson3D_AUX_SCALAR_NUM  = 1
  integer, parameter, public :: MGSmoother_Poisson3D_AUX_HVEC_NUM    = 1

  integer, parameter, public :: BCTYPE_NEUMANN_ID   = 1  !< dq/dz = 0 at the bottom and top boundaries
  integer, parameter, public :: BCTYPE_DIRICHLET_ID = 2  !< q = 0 at the bottom and top boundaries

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: cal_grad_lc

  integer, parameter :: AUXVAR_QZ_ID = 1  !< Vertical component of grad(q)
  integer, parameter :: AUXVAR_QX_ID = 2  !< Contravariant component of grad(q) in the alpha direction
  integer, parameter :: AUXVAR_QY_ID = 3  !< Contravariant component of grad(q) in the beta direction

  !> Default value of the coefficient C_tau in the penalty parameter (see formulation.pdf of the 2D version)
  real(RP), parameter :: PENALTY_COEF_DEFAULT = 2.5_RP

  !> Number of colors in the multicolor (block) Gauss-Seidel method:
  !! color = 2 * panel_group + mod(i+j+k,2), where panel_group = 1 for panels 5 and 6 (otherwise 0).
  integer, parameter :: NUM_COLOR = 4

  !> Face IDs of the bottom and top faces of the hexahedral element
  integer, parameter :: FACE_BOTTOM = 5
  integer, parameter :: FACE_TOP    = 6

contains
  !> Initialization
  subroutine MGSmoother_Poisson3D_Init( this, mesh, porder_list, bc_type )
    use scale_multigrid_smoother_base, only: MGSmootherBase3D_Init
    use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
    implicit none
    class(MGSmoother_Poisson3D), intent(inout) :: this
    class(MeshBase3D), intent(in), target :: mesh
    integer, intent(in) :: porder_list(:)
    character(len=*), intent(in) :: bc_type   !< Boundary condition at the bottom and top boundaries ('NEUMANN' or 'DIRICHLET')

    integer :: num_smooth_itr_max            = 10
    real(RP) :: residual_threshold_ratio_l2  = 1E-2_RP
    real(RP) :: residual_threshold_l2        = 1E-3_RP
    real(RP) :: residual_threshold_max       = 1E-3_RP
    integer :: history_residual_step   = 1

    integer :: num_smooth_itr_max_internal = 50                     !< Maximum number of internal iterations associated with local element solve
    real(RP) :: residual_threshold_l2_internal = 1E-3_RP !< Threshold for each internal iteration associated with local element solve

    character(len=H_SHORT) :: elem_solver = "PCG" !< Element solver type (only "PCG" is supported in this sample)

    real(RP) :: penalty_coef = PENALTY_COEF_DEFAULT !< Coefficient C_tau in the penalty parameter (p-MG levels)
    real(RP) :: penalty_coef_hMG = -1.0_RP           !< Coefficient C_tau in the penalty parameter (h-MG levels; if not positive, penalty_coef is used)

    namelist /PARAM_MGSmoother/ &
      num_smooth_itr_max,             &
      residual_threshold_ratio_l2,    &
      residual_threshold_l2,          &
      residual_threshold_max,         &
      history_residual_step,          &
      num_smooth_itr_max_internal,    &
      residual_threshold_l2_internal, &
      elem_solver,                    &
      penalty_coef,                   &
      penalty_coef_hMG

    integer :: ierr
    !---------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF, nml=PARAM_MGSmoother, iostat=ierr)
    if ( ierr < 0 ) then
      LOG_INFO("MGSmoother_Poisson3D_Init",*) "Not found namelist PARAM_MGSmoother in the configuration file."
    else if( ierr > 0 ) then
      LOG_ERROR("MGSmoother_Poisson3D_Init",*) "Error reading namelist PARAM_MGSmoother in the configuration file."
      call PRC_abort
    end if
    LOG_NML(PARAM_MGSmoother)

    !-
    call MGSmootherBase3D_Init( this )

    call this%Set_smoothing_param( num_smooth_itr_max, &
      residual_threshold_ratio_l2,residual_threshold_l2, residual_threshold_max, &
      history_residual_step )

    this%num_internal_itr_max = num_smooth_itr_max_internal
    this%residual_threshold_l2_internal_itr = residual_threshold_l2_internal

    if ( penalty_coef <= 0.0_RP ) then
      LOG_ERROR("MGSmoother_Poisson3D_Init",*) "penalty_coef must be positive. Check!", penalty_coef
      call PRC_abort
    end if
    this%penalty_coef = penalty_coef
    if ( penalty_coef_hMG > 0.0_RP ) then
      this%penalty_coef_hMG = penalty_coef_hMG
    else
      this%penalty_coef_hMG = penalty_coef
    end if
    LOG_INFO("MGSmoother_Poisson3D_Init",*) "C_tau (p-MG levels, h-MG levels) =", this%penalty_coef, this%penalty_coef_hMG

    !-
    select case ( trim(adjustl(elem_solver)) )
    case ( "PCG" )
    case ( "DIRECT" )
      LOG_ERROR("MGSmoother_Poisson3D_Init",*) "elem_solver='DIRECT' is not supported for cubed-sphere meshes. Use 'PCG'."
      call PRC_abort
    case default
      LOG_ERROR("MGSmoother_Poisson3D_Init",*) "Unknown elem_solver type: "//trim(elem_solver)
      call PRC_abort
    end select

    !- Boundary conditions at the bottom and top boundaries
    select case ( trim(adjustl(bc_type)) )
    case ( "NEUMANN" )
      this%bc_type_id = BCTYPE_NEUMANN_ID
    case ( "DIRICHLET" )
      this%bc_type_id = BCTYPE_DIRICHLET_ID
    case default
      LOG_ERROR("MGSmoother_Poisson3D_Init",*) "Unknown bc_type: "//trim(bc_type)//". Use 'NEUMANN' or 'DIRICHLET'."
      call PRC_abort
    end select

    !- Domain information (only the shallow approximation is supported)
    select type(mesh)
    type is (MeshCubedSphereDom3D)
      if ( .not. mesh%shallow_approx ) then
        LOG_ERROR("MGSmoother_Poisson3D_Init",*) "Currently, only shallow_approx=.true. is supported. Check!"
        call PRC_abort
      end if
      this%dom_zmin = mesh%zmin_gl
      this%dom_zmax = mesh%zmax_gl
    class default
      LOG_ERROR("MGSmoother_Poisson3D_Init",*) "This smoother supports only MeshCubedSphereDom3D. Check!"
      call PRC_abort
    end select

    if ( mesh%refElem3D%PolyOrder_h /= mesh%refElem3D%PolyOrder_v ) then
      LOG_ERROR("MGSmoother_Poisson3D_Init",*) "PolyOrder_h must equal PolyOrder_v. Check!"
      call PRC_abort
    end if

    return
  end subroutine MGSmoother_Poisson3D_Init

  !> Finalization
  subroutine MGSmoother_Poisson3D_Final(this)
    use scale_multigrid_smoother_base, only: MGSmootherBase3D_Final
    implicit none
    class(MGSmoother_Poisson3D), intent(inout) :: this
    !---------------------------------------------------------------------------

    call MGSmootherBase3D_Final( this )
    return
  end subroutine MGSmoother_Poisson3D_Final

  !> Smoother
  !! Note: f must be the right-hand side multiplied by Gsqrt.
!OCL SERIAL
  subroutine MGSmoother_Poisson3D_advance_itr_1step( this, &
    q, res,                                            &
    f, aux_var,                                        &
    itr, var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D, &
    cal_res_flag, zero_initial_guess,                  &
    mg_p_level, mg_h_level, pre_or_post_smooth )
    implicit none
    class(MGSmoother_Poisson3D), intent(inout) :: this
    type(MeshField3D), intent(inout), target :: q
    type(MeshField3D), intent(inout) :: res
    type(MeshField3D), intent(in) :: f
    type(MeshField3D), intent(inout), target :: aux_var(:)
    integer, intent(in) :: itr
    class(MeshFieldCommBase), intent(inout) :: var_comm
    class(MeshFieldCommBase), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    class(MeshBase3D), intent(in), target :: mesh3D
    logical, intent(in) :: cal_res_flag
    logical, intent(in) :: zero_initial_guess
    integer, intent(in) :: mg_p_level
    integer, intent(in) :: mg_h_level
    integer, intent(in) :: pre_or_post_smooth

    integer :: ldomID
    class(LocalMesh3D), pointer :: lmesh3D
    integer :: ke

    type(MeshFieldContainer) :: var_comm_list(1)
    integer :: ic
    integer :: color_seq
    real(RP) :: penalty_coef
    !---------------------------------------------------------------------------

    !- The coefficient of the penalty parameter: mg_h_level >= 1 on the h-MG levels (otherwise 0)
    if ( mg_h_level >= 1 ) then
      penalty_coef = this%penalty_coef_hMG
    else
      penalty_coef = this%penalty_coef
    end if

    if ( zero_initial_guess ) then
      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)
        !$omp parallel do
        do ke=lmesh3D%NeS, lmesh3D%NeE
          q%local(ldomID)%val(:,ke) = 0.0_RP
        end do
      end do
    end if

    if ( itr == 1 ) then
      call PROF_rapstart( "Poisson3dGlobal_smoother_res", 1 )
      call this%advance_itr_1step_color( &
        q, res,                                                                 &
        f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID), &
        var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,                           &
        .true., 0, penalty_coef                                                 )

      call this%Set_initial_residual( res )
      call flush(IO_FID_LOG)
      call PROF_rapend( "Poisson3dGlobal_smoother_res", 1 )
    end if

    call PROF_rapstart( "Poisson3dGlobal_smoother_mc", 1 )
    ! Multicolor ordering: 0,1,2,3 for pre-smoothing and 3,2,1,0 for post-smoothing
    do ic=1, NUM_COLOR
      if ( pre_or_post_smooth == MGSmoother_PRE_ID ) then
        color_seq = ic - 1
      else
        color_seq = NUM_COLOR - ic
      end if
      call this%advance_itr_1step_color( &
        q, res,                                                                 &
        f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID), &
        var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,                           &
        .false., color_seq, penalty_coef                                        )
    end do
    call PROF_rapend( "Poisson3dGlobal_smoother_mc", 1 )

    if ( cal_res_flag ) then
      call this%advance_itr_1step_color( &
        q, res,                                                                 &
        f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID), &
        var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,                           &
        cal_res_flag, 0, penalty_coef                                           )
    end if

    var_comm_list(1)%field3d => q
    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )
    return
  end subroutine MGSmoother_Poisson3D_advance_itr_1step

!-- private -----------------------

  !> Update q for elements with a specified color (or evaluate the residual if cal_res_flag = .true.)
  !! color_seq (0, ..., NUM_COLOR-1) = 2 * panel_group + color_id,
  !! where panel_group = 1 for panels 5 and 6 (otherwise 0) and color_id is the red-black color in each tile.
!OCL SERIAL
  subroutine Poisson3dGlobal_smoother_advance_itr_1step_color( this, &
    q, res,                                                    &
    f, qx, qy, qz,                                             &
    var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,              &
    cal_res_flag, color_seq, penalty_coef )
    implicit none
    class(MGSmoother_Poisson3D), intent(inout) :: this
    type(MeshField3D), intent(inout), target :: q
    type(MeshField3D), intent(inout) :: res
    type(MeshField3D), intent(in) :: f
    type(MeshField3D), intent(inout), target :: qx
    type(MeshField3D), intent(inout), target :: qy
    type(MeshField3D), intent(inout), target :: qz
    class(MeshFieldCommBase), intent(inout) :: var_comm
    class(MeshFieldCommBase), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    class(MeshBase3D), intent(in), target :: mesh3D
    logical, intent(in) :: cal_res_flag
    integer, intent(in) :: color_seq
    real(RP), intent(in) :: penalty_coef  !< Coefficient C_tau in the penalty parameter used at the current MG level

    integer :: ldomID
    class(LocalMesh3D), pointer :: lmesh3D
    integer :: color_id
    integer :: panel_group

    type(MeshFieldContainer) :: var_comm_list(1)
    type(MeshFieldContainer) :: aux_comm_list(3)
    !---------------------------------------------------------------------------

    !- Exchange q
    var_comm_list(1)%field3d => q

    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    !- Element-wise gradient (contravariant components)
    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(ldomID)
      call cal_grad_lc( qx%local(ldomID)%val, qy%local(ldomID)%val, qz%local(ldomID)%val, & ! (out)
        q%local(ldomID)%val, Dx, Dy, Dz, lmesh3D, lmesh3D%refElem3D                       ) ! (in)
    end do

    !- Exchange the gradient: the vertical component as a scalar and the horizontal components as a vector.
    !  The conversion of horizontal vector components between different panels is performed in the communicator.
    aux_comm_list(AUXVAR_QZ_ID)%field3d => qz
    aux_comm_list(AUXVAR_QX_ID)%field3d => qx
    aux_comm_list(AUXVAR_QY_ID)%field3d => qy

    call aux_comm%Put( aux_comm_list, 1 )
    call aux_comm%Exchange()
    call aux_comm%Get( aux_comm_list, 1 )

    color_id = mod( color_seq, 2 )
    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(ldomID)

      ! Skip local meshes which do not belong to the panel group of the current color
      if ( lmesh3D%panelID == 5 .or. lmesh3D%panelID == 6 ) then
        panel_group = 1
      else
        panel_group = 0
      end if
      if ( ( .not. cal_res_flag ) .and. panel_group /= color_seq / 2 ) cycle

      call this%cal_q_lc_elem_PCG( &
        q%local(ldomID)%val, res%local(ldomID)%val,             &
        color_id, cal_res_flag, f%local(ldomID)%val,            &
        qx%local(ldomID)%val, qy%local(ldomID)%val, qz%local(ldomID)%val, &
        Dx, Dy, Dz, Lift, lmesh3D, lmesh3D%refElem3D, penalty_coef )
    end do
    return
  end subroutine Poisson3dGlobal_smoother_advance_itr_1step_color

  !> Calculate contravariant components of the element-wise gradient (without lifting terms)
  !! sigma~^I = G^{IJ} d q / d xi^J (horizontal), sigma~^3 = d q / d z (vertical; shallow approximation)
!OCL SERIAL
  subroutine cal_grad_lc( qx, qy, qz, &
    q, Dx, Dy, Dz, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qy(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qz(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz

    integer :: ke
    real(RP) :: q1(elem%Np), q2(elem%Np), q3(elem%Np)
    real(RP) :: GIJ(elem%Np,2,2)
    !----------------------------
    !$omp parallel do private(q1, q2, q3, GIJ)
    do ke=lmesh%NeS, lmesh%NeE
      call Sparsemat_matmul(Dx, q(:,ke), q1)
      call Sparsemat_matmul(Dy, q(:,ke), q2)
      call Sparsemat_matmul(Dz, q(:,ke), q3)

      ! Covariant components
      q1(:) = lmesh%Escale(:,ke,1,1) * q1(:)
      q2(:) = lmesh%Escale(:,ke,2,2) * q2(:)
      q3(:) = lmesh%Escale(:,ke,3,3) * q3(:)

      ! Contravariant components
      call get_GIJ_elem( GIJ, ke, lmesh, elem )
      qx(:,ke) = GIJ(:,1,1) * q1(:) + GIJ(:,1,2) * q2(:)
      qy(:,ke) = GIJ(:,2,1) * q1(:) + GIJ(:,2,2) * q2(:)
      qz(:,ke) = q3(:)
    end do
    return
  end subroutine cal_grad_lc

  !> Get the horizontal contravariant metric tensor at the 3D nodes of an element
!OCL SERIAL
  subroutine get_GIJ_elem( GIJ, ke, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: GIJ(elem%Np,2,2)
    integer, intent(in) :: ke

    integer :: ke2D
    !----------------------------
    ke2D = lmesh%EMap3Dto2D(ke)
    GIJ(:,1,1) = lmesh%GIJ(elem%IndexH2Dto3D(:),ke2D,1,1)
    GIJ(:,1,2) = lmesh%GIJ(elem%IndexH2Dto3D(:),ke2D,1,2)
    GIJ(:,2,1) = lmesh%GIJ(elem%IndexH2Dto3D(:),ke2D,2,1)
    GIJ(:,2,2) = lmesh%GIJ(elem%IndexH2Dto3D(:),ke2D,2,2)
    return
  end subroutine get_GIJ_elem

  !> Element-block Gauss-Seidel update with the PCG method for element-local problems
  !! (or evaluation of the residual, res = rhs - A q, if cal_res = .true.)
  !! The element-local problem A~_KK dq = r is solved as (M A~_KK) dq = M r,
  !! because M A~_KK (M: mass matrix of the reference element) is (approximately) symmetric,
  !! while A~_KK itself is not symmetric in the Euclidean inner product.
!OCL SERIAL
  subroutine cal_q_lc_elem_PCG( this, q, res,  &
    color_id, cal_res, rhs, qx, qy, qz,        &
    Dx, Dy, Dz, Lift,                          &
    lmesh, elem, penalty_coef                  )
    implicit none
    class(MGSmoother_Poisson3D), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: q(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: res(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: rhs(elem%Np,lmesh%NeA)   !< Right-hand side multiplied by Gsqrt
    real(RP), intent(in) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qy(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qz(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: color_id
    logical, intent(in) :: cal_res
    real(RP), intent(in) :: penalty_coef  !< Coefficient C_tau in the penalty parameter

    integer :: i, j, k, ke
    integer :: i0, ii, ncol

    real(RP) :: gtau(elem%NfpTot)
    real(RP) :: GIJ(elem%Np,2,2)
    logical :: is_bnd_bottom, is_bnd_top

    real(RP) :: Ax(elem%Np)
    real(RP) :: x(elem%Np)
    real(RP) :: r(elem%Np)
    real(RP) :: z(elem%Np)
    real(RP) :: q_(elem%Np)
    real(RP) :: Mr(elem%Np)
    real(RP) :: rnorm0, rnorm
    real(RP) :: rz_old, rz_new
    real(RP) :: alpha, beta

    integer :: itr

    real(RP) :: PmatDiag(elem%Np)
    real(RP) :: hk
    real(RP), parameter :: c1 = 1.0_RP
    real(RP), parameter :: c2 = 1.0_RP
    real(RP) :: G11(elem%Np), G22(elem%Np)
    integer :: ke2D

    real(RP) :: int_w(elem%Np)            !< Row sums of the mass matrix (used to scale the preconditioner)
    real(RP) :: M1D_h(elem%Nnode_h1D,elem%Nnode_h1D) !< 1D mass matrix (horizontal)
    real(RP) :: M1D_v(elem%Nnode_v,elem%Nnode_v)     !< 1D mass matrix (vertical)
    real(RP) :: ones(elem%Np)

    integer :: iM(elem%NfpTot)
    integer :: p, f
    real(RP), parameter :: EPS_Z = 1.0E-10_RP
    !---------------------------

    call PROF_rapstart( "Poisson3dGlobal_smoother_cal_q_lc", 1 )

    !- Local node indices on the faces (horizontal faces 1-4, bottom face 5, top face 6)
    do f=1, elem%Nfaces_h
    do p=1, elem%Nfp_h
      iM(p + (f-1)*elem%Nfp_h) = elem%Fmask_h(p,f)
    end do
    end do
    do f=1, elem%Nfaces_v
    do p=1, elem%Nfp_v
      iM(p + (f-1)*elem%Nfp_v + elem%Nfp_h*elem%Nfaces_h) = elem%Fmask_v(p,f)
    end do
    end do

    !- 1D mass matrices: M = M1D_v x M1D_h x M1D_h (PolyOrder_h = PolyOrder_v is assumed)
    !  M(i1,i2) for j=k=1 equals M1D_h(i1,i2) * M1D_h(1,1) * M1D_v(1,1) = M1D_h(i1,i2) * m11^2 with m11 = M(1,1)^(1/3)
    call get_mass1D( M1D_h, M1D_v, elem )
    ones(:) = 1.0_RP
    call apply_mass_tensor3D( int_w, ones, M1D_h, M1D_v, elem%Nnode_h1D, elem%Nnode_v )

    !$omp parallel do collapse(2) private(i, ke, i0, ii, ncol, gtau, GIJ, is_bnd_bottom, is_bnd_top, &
    !$omp Ax, x, r, z, q_, Mr, rnorm0, rnorm, rz_old, rz_new, alpha, beta, itr,                    &
    !$omp PmatDiag, hk, G11, G22, ke2D )
    do k=1, lmesh%NeZ
    do j=1, lmesh%NeY

      if ( cal_res ) then
        i0 = 1
        ncol = lmesh%NeX
      else
        i0 = 1 + modulo(color_id - modulo(j + k, 2), 2)  ! 1 or 2
        ncol = (lmesh%NeX - i0 + 2)/2   ! number of i in {i0, i0+2, ...} <= NeX (0 if NeX < i0)
      end if

      do ii=1, ncol
        if (cal_res) then
          i = ii
        else
          i = i0 + (ii-1)*2
        end if
        ke = i + (j-1)*lmesh%NeX + (k-1)*lmesh%NeX*lmesh%NeY

        !- Bottom and top boundaries of the computational domain
        is_bnd_bottom = ( k == 1         .and. abs( lmesh%zmin - this%dom_zmin ) < EPS_Z * max( 1.0_RP, abs(this%dom_zmax - this%dom_zmin) ) )
        is_bnd_top    = ( k == lmesh%NeZ .and. abs( lmesh%zmax - this%dom_zmax ) < EPS_Z * max( 1.0_RP, abs(this%dom_zmax - this%dom_zmin) ) )

        call get_GIJ_elem( GIJ, ke, lmesh, elem )
        call cal_penalty( gtau, &
          ke, iM, GIJ, penalty_coef, lmesh, elem )

        call eval_Ax( Ax, &
          ke, q, qx, qy, qz, gtau, GIJ,                         &
          is_bnd_bottom, is_bnd_top, this%bc_type_id,           &
          Dx, Dy, Dz, Lift, iM, lmesh, elem                     )

        r(:) = rhs(:,ke) - Ax(:)

        if ( cal_res ) then
          res(:,ke) = r(:)
        else
          ! Mass-matrix-weighted residual for the system (M A~_KK) dq = M r
          call apply_mass_tensor3D( Mr, r, M1D_h, M1D_v, elem%Nnode_h1D, elem%Nnode_v )
          r(:) = Mr(:)
          rnorm0 = sqrt( sum( r(:) * r(:) ) )

          if ( rnorm0 > 1.0E-16_RP ) then
            ! Precond: physical element width is estimated using the covariant metric
            ke2D = lmesh%EMap3Dto2D(ke)
            G11(:) = lmesh%G_ij(elem%IndexH2Dto3D(:),ke2D,1,1)
            G22(:) = lmesh%G_ij(elem%IndexH2Dto3D(:),ke2D,2,2)
            hk = 6.0_RP / maxval( lmesh%Escale(:,ke,1,1) / sqrt( G11(:) ) &
                                + lmesh%Escale(:,ke,2,2) / sqrt( G22(:) ) &
                                + lmesh%Escale(:,ke,3,3)                  )
            PmatDiag(:) = int_w(:) * lmesh%Gsqrt(:,ke) * ( c1 * elem%PolyOrder_h**2 / hk + c2 / hk**2 )

            x(:) = q(:,ke)
            z(:) = r(:) / PmatDiag(:)

            !-
            q_(:) = z(:)
            rz_old = sum( r(:) * z(:) )

            do itr=1, this%num_internal_itr_max
              !-
              call eval_Ax_itr( Ax, &
                ke, q_, gtau, GIJ,                                  &
                is_bnd_bottom, is_bnd_top, this%bc_type_id,         &
                Dx, Dy, Dz, Lift, iM, lmesh, elem                   )
              call apply_mass_tensor3D( Mr, Ax, M1D_h, M1D_v, elem%Nnode_h1D, elem%Nnode_v )
              Ax(:) = Mr(:)

              alpha = rz_old / sum(q_(:) * Ax(:))
              x(:) = x(:) + alpha * q_(:)
              r(:) = r(:) - alpha * Ax(:)

              rnorm = sqrt( sum( r(:) * r(:) ) )
              if ( rnorm / rnorm0 < this%residual_threshold_l2_internal_itr ) then
                exit
              else if ( itr == this%num_internal_itr_max ) then
                LOG_INFO("Poisson3dGlobal_smoother_cal_q_lc",*) "ke=", ke, "rnorm/rnorm0=", rnorm / rnorm0, rnorm, rnorm0
                exit
              end if

              ! Precond
              z(:) = r(:) / PmatDiag(:)

              rz_new = sum( r(:) * z(:) )

              !-
              beta = rz_new / rz_old
              q_(:) = z(:) + beta * q_(:)

              rz_old = rz_new
            end do
            q(:,ke) = x(:)
          end if
        end if
      end do
    end do
    end do

    call PROF_rapend( "Poisson3dGlobal_smoother_cal_q_lc", 1 )
    return
  end subroutine cal_q_lc_elem_PCG

  !> Calculate the penalty parameter on the face nodes:
  !! tau~ = C_tau p(p+1)/2 Fscale * Gsqrt * G^{ij} n_i n_j  (G^{33} = 1 with the shallow approximation)
!OCL SERIAL
  subroutine cal_penalty( gtau, &
    ke, iM, GIJ, penalty_coef, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: gtau(elem%NfpTot)
    integer, intent(in) :: ke
    integer, intent(in) :: iM(elem%NfpTot)
    real(RP), intent(in) :: GIJ(elem%Np,2,2)
    real(RP), intent(in) :: penalty_coef  !< Coefficient C_tau

    real(RP) :: nx(elem%NfpTot), ny(elem%NfpTot), nz(elem%NfpTot)
    !------------------------------------------

    nx(:) = lmesh%normal_fn(:,ke,1)
    ny(:) = lmesh%normal_fn(:,ke,2)
    nz(:) = lmesh%normal_fn(:,ke,3)

    gtau(:) = penalty_coef * 0.5_RP * real( elem%PolyOrder_h * ( elem%PolyOrder_h + 1 ), kind=RP ) &
            * lmesh%Fscale(:,ke) * lmesh%Gsqrt(iM,ke)                                       &
            * (            GIJ(iM,1,1) * nx(:) * nx(:)                                       &
                + 2.0_RP * GIJ(iM,1,2) * nx(:) * ny(:)                                       &
                +          GIJ(iM,2,2) * ny(:) * ny(:)                                       &
                +                        nz(:) * nz(:)                                       )
    return
  end subroutine cal_penalty

!- Evaluation of A q ------------------------------------------------------------

  !> Evaluate (A~ q)_K with the values of neighbor elements
  !!  A~ q = D_i( Gsqrt sigma^i ) + Lift[ Fscale ( F^ - Gsqrt sigma^{i,-} n_i ) ]
!OCL SERIAL
  subroutine eval_Ax( Ax, &
    ke, q, qx, qy, qz, gtau, GIJ,          &
    is_bnd_bottom, is_bnd_top, bc_type_id, &
    Dx, Dy, Dz, Lift, iM, lmesh, elem      )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np)
    integer, intent(in) :: ke
    real(RP), intent(in) :: q(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: qx(elem%Np*lmesh%NeA)    !< Contravariant element-wise gradient (alpha)
    real(RP), intent(in) :: qy(elem%Np*lmesh%NeA)    !< Contravariant element-wise gradient (beta)
    real(RP), intent(in) :: qz(elem%Np*lmesh%NeA)    !< Element-wise gradient (vertical)
    real(RP), intent(in) :: gtau(elem%NfpTot)
    real(RP), intent(in) :: GIJ(elem%Np,2,2)
    logical, intent(in) :: is_bnd_bottom, is_bnd_top
    integer, intent(in) :: bc_type_id
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: iM(elem%NfpTot)

    real(RP) :: q_e(elem%Np), sig_e(elem%Np,3)
    real(RP) :: qP(elem%NfpTot), sigP(elem%NfpTot,3)
    integer :: vP(elem%NfpTot)
    !---------------------------------

    vP(:) = lmesh%vmapP(:,ke)
    q_e(:) = q((ke-1)*elem%Np+1:ke*elem%Np)
    sig_e(:,1) = qx((ke-1)*elem%Np+1:ke*elem%Np)
    sig_e(:,2) = qy((ke-1)*elem%Np+1:ke*elem%Np)
    sig_e(:,3) = qz((ke-1)*elem%Np+1:ke*elem%Np)

    ! Values of neighbor elements on the faces (halo values are used for inter-tile faces)
    qP(:) = q(vP); sigP(:,1) = qx(vP); sigP(:,2) = qy(vP); sigP(:,3) = qz(vP)

    call eval_Ax_core( Ax, &
      ke, q_e, sig_e, qP, sigP, gtau, GIJ,   &
      is_bnd_bottom, is_bnd_top, bc_type_id, &
      Dx, Dy, Dz, Lift, iM, lmesh, elem      )
    return
  end subroutine eval_Ax

  !> Evaluate (A~_KK q)_K for the element-local problem (values of neighbor elements are set to zero)
!OCL SERIAL
  subroutine eval_Ax_itr( Ax, &
    ke, q, gtau, GIJ,                      &
    is_bnd_bottom, is_bnd_top, bc_type_id, &
    Dx, Dy, Dz, Lift, iM, lmesh, elem      )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np)
    integer, intent(in) :: ke
    real(RP), intent(in) :: q(elem%Np)
    real(RP), intent(in) :: gtau(elem%NfpTot)
    real(RP), intent(in) :: GIJ(elem%Np,2,2)
    logical, intent(in) :: is_bnd_bottom, is_bnd_top
    integer, intent(in) :: bc_type_id
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: iM(elem%NfpTot)

    real(RP) :: sig_e(elem%Np,3)
    real(RP) :: qP(elem%NfpTot), sigP(elem%NfpTot,3)
    !---------------------------------

    call cal_sig_e( sig_e, q, GIJ, ke, Dx, Dy, Dz, lmesh, elem )
    qP(:) = 0.0_RP; sigP(:,:) = 0.0_RP

    call eval_Ax_core( Ax, &
      ke, q, sig_e, qP, sigP, gtau, GIJ,     &
      is_bnd_bottom, is_bnd_top, bc_type_id, &
      Dx, Dy, Dz, Lift, iM, lmesh, elem      )
    return
  end subroutine eval_Ax_itr

  !> Contravariant components of the element-wise gradient for a single element
!OCL SERIAL
  subroutine cal_sig_e( sig_e, q, GIJ, ke, Dx, Dy, Dz, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: sig_e(elem%Np,3)
    real(RP), intent(in) :: q(elem%Np)
    real(RP), intent(in) :: GIJ(elem%Np,2,2)
    integer, intent(in) :: ke
    type(SparseMat), intent(in) :: Dx, Dy, Dz

    real(RP) :: g1(elem%Np), g2(elem%Np), g3(elem%Np)
    !---------------------------------
    call sparsemat_matmul( Dx, q, g1 )
    call sparsemat_matmul( Dy, q, g2 )
    call sparsemat_matmul( Dz, q, g3 )
    g1(:) = lmesh%Escale(:,ke,1,1) * g1(:)
    g2(:) = lmesh%Escale(:,ke,2,2) * g2(:)
    g3(:) = lmesh%Escale(:,ke,3,3) * g3(:)
    sig_e(:,1) = GIJ(:,1,1) * g1(:) + GIJ(:,1,2) * g2(:)
    sig_e(:,2) = GIJ(:,2,1) * g1(:) + GIJ(:,2,2) * g2(:)
    sig_e(:,3) = g3(:)
    return
  end subroutine cal_sig_e

  !> Core of the evaluation of A~ q for an element with given neighbor values on the faces.
  !! At the bottom and top boundaries, the neighbor values are replaced by ghost values:
  !!  Neumann  : q+ = q-,  sigma~+ . n = - sigma~- . n (the tangential components are unchanged)
  !!  Dirichlet: q+ = -q-, sigma~+ = sigma~-
!OCL SERIAL
  subroutine eval_Ax_core( Ax, &
    ke, q_e, sig_e, qP_in, sigP_in, gtau, GIJ, &
    is_bnd_bottom, is_bnd_top, bc_type_id,     &
    Dx, Dy, Dz, Lift, iM, lmesh, elem          )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np)
    integer, intent(in) :: ke
    real(RP), intent(in) :: q_e(elem%Np)            !< q of the element
    real(RP), intent(in) :: sig_e(elem%Np,3)        !< Element-wise contravariant gradient of the element
    real(RP), intent(in) :: qP_in(elem%NfpTot)      !< q of neighbor elements on the faces
    real(RP), intent(in) :: sigP_in(elem%NfpTot,3)  !< Element-wise contravariant gradient of neighbor elements on the faces
    real(RP), intent(in) :: gtau(elem%NfpTot)
    real(RP), intent(in) :: GIJ(elem%Np,2,2)
    logical, intent(in) :: is_bnd_bottom, is_bnd_top
    integer, intent(in) :: bc_type_id
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: iM(elem%NfpTot)

    real(RP) :: qM(elem%NfpTot), qP(elem%NfpTot)
    real(RP) :: sigM_e(elem%NfpTot,3), sigP(elem%NfpTot,3)
    real(RP) :: nx(elem%NfpTot), ny(elem%NfpTot), nz(elem%NfpTot)
    real(RP) :: dq(elem%NfpTot)
    real(RP) :: L1(elem%Np), L2(elem%Np), L3(elem%Np)
    real(RP) :: Sig(elem%Np,3)
    real(RP) :: GsqrtM(elem%NfpTot)
    real(RP) :: BndFlx(elem%NfpTot)
    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDFlx(elem%Np)
    integer :: is, ie
    !---------------------------------

    nx(:) = lmesh%normal_fn(:,ke,1)
    ny(:) = lmesh%normal_fn(:,ke,2)
    nz(:) = lmesh%normal_fn(:,ke,3)
    qM(:) = q_e(iM)
    sigM_e(:,1) = sig_e(iM,1); sigM_e(:,2) = sig_e(iM,2); sigM_e(:,3) = sig_e(iM,3)
    qP(:) = qP_in(:); sigP(:,:) = sigP_in(:,:)

    !- Ghost values at the bottom and top boundaries
    if ( is_bnd_bottom ) then
      is = elem%Nfp_h*elem%Nfaces_h + (FACE_BOTTOM-elem%Nfaces_h-1)*elem%Nfp_v + 1; ie = is + elem%Nfp_v - 1
      call set_ghost( qP(is:ie), sigP(is:ie,:), qM(is:ie), sigM_e(is:ie,:), bc_type_id, elem%Nfp_v )
    end if
    if ( is_bnd_top ) then
      is = elem%Nfp_h*elem%Nfaces_h + (FACE_TOP-elem%Nfaces_h-1)*elem%Nfp_v + 1; ie = is + elem%Nfp_v - 1
      call set_ghost( qP(is:ie), sigP(is:ie,:), qM(is:ie), sigM_e(is:ie,:), bc_type_id, elem%Nfp_v )
    end if

    !- (a) Lifting terms of the gradient (covariant components): 1/2 (q+ - q-) n_i
    dq(:) = qP(:) - qM(:)
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * 0.5_RP * dq(:) * nx(:), L1 )
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * 0.5_RP * dq(:) * ny(:), L2 )
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * 0.5_RP * dq(:) * nz(:), L3 )

    !- (b) Contravariant components: sigma = sigma~ + G^{ij} L_j
    Sig(:,1) = sig_e(:,1) + GIJ(:,1,1) * L1(:) + GIJ(:,1,2) * L2(:)
    Sig(:,2) = sig_e(:,2) + GIJ(:,2,1) * L1(:) + GIJ(:,2,2) * L2(:)
    Sig(:,3) = sig_e(:,3) + L3(:)

    !- (c) Divergence: F^ - Gsqrt sigma^{i,-} n_i, where
    !      F^ = 1/2 Gsqrt ( sigma~^{i,-} + sigma~^{i,+} ) n_i + 1/2 tau~ ( q+ - q- )
    GsqrtM(:) = lmesh%Gsqrt(iM,ke)
    BndFlx(:) = GsqrtM(:) * (                                         &
          ( 0.5_RP * ( sigM_e(:,1) + sigP(:,1) ) - Sig(iM,1) ) * nx(:) &
        + ( 0.5_RP * ( sigM_e(:,2) + sigP(:,2) ) - Sig(iM,2) ) * ny(:) &
        + ( 0.5_RP * ( sigM_e(:,3) + sigP(:,3) ) - Sig(iM,3) ) * nz(:) ) &
      + 0.5_RP * gtau(:) * dq(:)

    call sparsemat_matmul( Dx, lmesh%Gsqrt(:,ke) * Sig(:,1), Fx )
    call sparsemat_matmul( Dy, lmesh%Gsqrt(:,ke) * Sig(:,2), Fy )
    call sparsemat_matmul( Dz, lmesh%Gsqrt(:,ke) * Sig(:,3), Fz )
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * BndFlx(:), LiftDFlx )
    Ax(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDFlx(:)
    return
  end subroutine eval_Ax_core

  !> Set ghost values on a bottom/top boundary face (the normal is +/- e_z)
!OCL SERIAL
  subroutine set_ghost( qP, sigP, qM, sigM, bc_type_id, Nfp )
    implicit none
    integer, intent(in) :: Nfp
    real(RP), intent(out) :: qP(Nfp)
    real(RP), intent(out) :: sigP(Nfp,3)
    real(RP), intent(in) :: qM(Nfp)
    real(RP), intent(in) :: sigM(Nfp,3)
    integer, intent(in) :: bc_type_id
    !---------------------------------
    select case( bc_type_id )
    case( BCTYPE_NEUMANN_ID )
      qP(:) = qM(:)
      sigP(:,1) = sigM(:,1); sigP(:,2) = sigM(:,2)
      sigP(:,3) = - sigM(:,3)
    case( BCTYPE_DIRICHLET_ID )
      qP(:) = - qM(:)
      sigP(:,:) = sigM(:,:)
    end select
    return
  end subroutine set_ghost

!- Mass matrix with the tensor-product structure ---------------------------------

  !> Extract 1D mass matrices from the 3D mass matrix of a hexahedral element
  !! (PolyOrder_h = PolyOrder_v is assumed, so that M1D_h = M1D_v).
!OCL SERIAL
  subroutine get_mass1D( M1D_h, M1D_v, elem )
    implicit none
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: M1D_h(elem%Nnode_h1D,elem%Nnode_h1D)
    real(RP), intent(out) :: M1D_v(elem%Nnode_v,elem%Nnode_v)

    real(RP) :: m11
    integer :: k1, k2, Nh2
    !---------------------------------
    m11 = elem%M(1,1)**(1.0_RP/3.0_RP)
    M1D_h(:,:) = elem%M(1:elem%Nnode_h1D,1:elem%Nnode_h1D) / m11**2
    Nh2 = elem%Nnode_h1D**2
    do k2=1, elem%Nnode_v
    do k1=1, elem%Nnode_v
      M1D_v(k1,k2) = elem%M(1+(k1-1)*Nh2,1+(k2-1)*Nh2) / m11**2
    end do
    end do
    return
  end subroutine get_mass1D

  !> Apply the mass matrix of a hexahedral element using its tensor-product structure:
  !! y(i,j,k) = sum_{a,b,c} M1D_h(i,a) M1D_h(j,b) M1D_v(k,c) x(a,b,c)
!OCL SERIAL
  subroutine apply_mass_tensor3D( y, x, M1D_h, M1D_v, Nh, Nv )
    implicit none
    integer, intent(in) :: Nh, Nv
    real(RP), intent(out) :: y(Nh,Nh,Nv)
    real(RP), intent(in) :: x(Nh,Nh,Nv)
    real(RP), intent(in) :: M1D_h(Nh,Nh)
    real(RP), intent(in) :: M1D_v(Nv,Nv)

    real(RP) :: t1(Nh,Nh,Nv), t2(Nh,Nh,Nv)
    integer :: i, j, k, a
    !-------------------------------------------
    do k=1, Nv
    do j=1, Nh
    do i=1, Nh
      t1(i,j,k) = 0.0_RP
      do a=1, Nh
        t1(i,j,k) = t1(i,j,k) + M1D_h(i,a) * x(a,j,k)
      end do
    end do
    end do
    end do
    do k=1, Nv
    do j=1, Nh
    do i=1, Nh
      t2(i,j,k) = 0.0_RP
      do a=1, Nh
        t2(i,j,k) = t2(i,j,k) + M1D_h(j,a) * t1(i,a,k)
      end do
    end do
    end do
    end do
    do k=1, Nv
    do j=1, Nh
    do i=1, Nh
      y(i,j,k) = 0.0_RP
      do a=1, Nv
        y(i,j,k) = y(i,j,k) + M1D_v(k,a) * t2(i,j,a)
      end do
    end do
    end do
    end do
    return
  end subroutine apply_mass_tensor3D

end module mod_poisson3dGlobal_smoother
