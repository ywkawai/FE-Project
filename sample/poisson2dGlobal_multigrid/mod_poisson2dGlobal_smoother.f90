!-------------------------------------------------------------------------------
!> module Multigrid smoother for the Poisson equation on the sphere
!!
!! @par Description
!!      A smoother of the multigrid method for the 2D Poisson equation on a cubed-sphere mesh. 
!!      The discretization follows the equation multiplied by Gsqrt,
!!
!!        d/dxi^i ( Gsqrt G^{ij} dq/dxi^j ) = Gsqrt f,
!!
!!      is discretized with a symmetric interior penalty (SIP) type nodal DG method,
!!      using the mass matrix of the reference element (without Gsqrt).
!!      The right-hand side given to the smoother must be already multiplied by Gsqrt.
!!      Note that the restricted residual on coarse levels also includes Gsqrt.
!!
!!      Smoother: element-block Gauss-Seidel with red-black ordering in each tile.
!!      The element-local problem is solved by a preconditioned conjugate gradient method.
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_poisson2dGlobal_smoother
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

  use scale_element_base, only: ElementBase2D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D

  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_base, only: MeshFieldCommBase

  use scale_multigrid_smoother_base, only: &
    MGSmootherBase2D, MGSmoother_PRE_ID, MGSmoother_POST_ID
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type, public, extends(MGSmootherBase2D) :: MGSmoother_Poisson2D
    integer :: num_internal_itr_max                !< Maximum number of internal iterations associated with local element solve
    real(RP) :: residual_threshold_l2_internal_itr !< Threshold for each internal iteration associated with local element solve
  contains
    procedure :: Init => MGSmoother_Poisson2D_Init
    procedure :: Final => MGSmoother_Poisson2D_Final
    procedure :: Advance_itr_1step => MGSmoother_Poisson2D_advance_itr_1step
    procedure, private :: advance_itr_1step_color => Poisson2dGlobal_smoother_advance_itr_1step_color
    procedure, private :: cal_q_lc_elem_PCG => cal_q_lc_elem_PCG
  end type MGSmoother_Poisson2D

  !> Number of auxiliary variables: one horizontal vector (contravariant components of the element-wise gradient)
  integer, parameter, public :: MGSmoother_Poisson2D_AUX_SCALAR_NUM  = 0
  integer, parameter, public :: MGSmoother_Poisson2D_AUX_VEC_NUM     = 1

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: cal_grad_lc

  integer, parameter :: AUXVAR_QX_ID = 1  !< Contravariant component of grad(q) in the alpha direction
  integer, parameter :: AUXVAR_QY_ID = 2  !< Contravariant component of grad(q) in the beta direction

  real(RP), parameter :: PENALTY_COEF = 1.2_RP  !< Coefficient C_tau in the penalty parameter

contains
  !> Initialization
  subroutine MGSmoother_Poisson2D_Init( this, mesh, porder_list )
    use scale_multigrid_smoother_base, only: MGSmootherBase2D_Init
    implicit none
    class(MGSmoother_Poisson2D), intent(inout) :: this
    class(MeshBase2D), intent(in), target :: mesh
    integer, intent(in) :: porder_list(:)

    integer :: num_smooth_itr_max            = 10
    real(RP) :: residual_threshold_ratio_l2  = 1E-2_RP
    real(RP) :: residual_threshold_l2        = 1E-3_RP
    real(RP) :: residual_threshold_max       = 1E-3_RP
    integer :: history_residual_step   = 1

    integer :: num_smooth_itr_max_internal = 20                     !< Maximum number of internal iterations associated with local element solve
    real(RP) :: residual_threshold_l2_internal = 1E-3_RP !< Threshold for each internal iteration associated with local element solve

    character(len=H_SHORT) :: elem_solver = "PCG" !< Element solver type (only "PCG" is supported in this sample)

    namelist /PARAM_MGSmoother/ &
      num_smooth_itr_max,             &
      residual_threshold_ratio_l2,    &
      residual_threshold_l2,          &
      residual_threshold_max,         &
      history_residual_step,          &
      num_smooth_itr_max_internal,    &
      residual_threshold_l2_internal, &
      elem_solver

    integer :: ierr
    !---------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF, nml=PARAM_MGSmoother, iostat=ierr)
    if ( ierr < 0 ) then
      LOG_INFO("MGSmoother_Poisson2D_Init",*) "Not found namelist PARAM_MGSmoother in the configuration file."
    else if( ierr > 0 ) then
      LOG_ERROR("MGSmoother_Poisson2D_Init",*) "Error reading namelist PARAM_MGSmoother in the configuration file."
      call PRC_abort
    end if
    LOG_NML(PARAM_MGSmoother)

    !-
    call MGSmootherBase2D_Init( this )

    call this%Set_smoothing_param( num_smooth_itr_max, &
      residual_threshold_ratio_l2,residual_threshold_l2, residual_threshold_max, &
      history_residual_step )

    this%num_internal_itr_max = num_smooth_itr_max_internal
    this%residual_threshold_l2_internal_itr = residual_threshold_l2_internal

    !-
    select case ( trim(adjustl(elem_solver)) )
    case ( "PCG" )
    case ( "DIRECT" )
      LOG_ERROR("MGSmoother_Poisson2D_Init",*) "elem_solver='DIRECT' is not supported for cubed-sphere meshes. Use 'PCG'."
      call PRC_abort
    case default
      LOG_ERROR("MGSmoother_Poisson2D_Init",*) "Unknown elem_solver type: "//trim(elem_solver)
      call PRC_abort
    end select

    return
  end subroutine MGSmoother_Poisson2D_Init

  !> Finalization
  subroutine MGSmoother_Poisson2D_Final(this)
    use scale_multigrid_smoother_base, only: MGSmootherBase2D_Final
    implicit none
    class(MGSmoother_Poisson2D), intent(inout) :: this
    !---------------------------------------------------------------------------

    call MGSmootherBase2D_Final( this )
    return
  end subroutine MGSmoother_Poisson2D_Final

  !> Smoother
  !! Note: f must be the right-hand side multiplied by Gsqrt.
!OCL SERIAL
  subroutine MGSmoother_Poisson2D_advance_itr_1step( this, &
    q, res,                                         &
    f, aux_var,                                     &
    itr, var_comm, aux_comm, Dx, Dy, Lift, mesh2D,  &
    cal_res_flag, zero_initial_guess,               &
    mg_p_level, mg_h_level, pre_or_post_smooth      )
    implicit none
    class(MGSmoother_Poisson2D), intent(inout) :: this
    type(MeshField2D), intent(inout), target :: q
    type(MeshField2D), intent(inout) :: res
    type(MeshField2D), intent(in) :: f
    type(MeshField2D), intent(inout), target :: aux_var(:)
    integer, intent(in) :: itr
    class(MeshFieldCommBase), intent(inout) :: var_comm
    class(MeshFieldCommBase), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift
    class(MeshBase2D), intent(in), target :: mesh2D
    logical, intent(in) :: cal_res_flag
    logical, intent(in) :: zero_initial_guess
    integer, intent(in) :: mg_p_level
    integer, intent(in) :: mg_h_level
    integer, intent(in) :: pre_or_post_smooth

    integer :: ldomID
    class(LocalMesh2D), pointer :: lmesh2D
    integer :: ke

    type(MeshFieldContainer) :: var_comm_list(1)
    integer :: color_id_offset
    !---------------------------------------------------------------------------

    if ( pre_or_post_smooth == MGSmoother_PRE_ID ) then
      color_id_offset = 2
    else
      color_id_offset = 1
    end if

    if ( zero_initial_guess ) then
      do ldomID=1, mesh2D%LOCAL_MESH_NUM
        lmesh2D => mesh2D%lcmesh_list(ldomID)
        !$omp parallel do
        do ke=lmesh2D%NeS, lmesh2D%NeE
          q%local(ldomID)%val(:,ke) = 0.0_RP
        end do
      end do
    end if

    if ( itr == 1 ) then
      call PROF_rapstart( "Poisson2dGlobal_smoother_res", 1 )
      call this%advance_itr_1step_color( &
        q, res,                                          &
        f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), &
        var_comm, aux_comm, Dx, Dy, Lift, mesh2D,        &
        .true., 0                                        )

      call this%Set_initial_residual( res )
      call flush(IO_FID_LOG)
      call PROF_rapend( "Poisson2dGlobal_smoother_res", 1 )
    end if

    call PROF_rapstart( "Poisson2dGlobal_smoother_rb", 1 )
    ! Red
    call this%advance_itr_1step_color( &
      q, res,                                          &
      f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), &
      var_comm, aux_comm, Dx, Dy, Lift, mesh2D,        &
      .false., mod(color_id_offset,2)                  )

    ! Black
    call this%advance_itr_1step_color( &
      q, res,                                          &
      f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), &
      var_comm, aux_comm, Dx, Dy, Lift, mesh2D,        &
      .false., mod(color_id_offset+1,2)                )
    call PROF_rapend( "Poisson2dGlobal_smoother_rb", 1 )

    if ( cal_res_flag ) then
      call this%advance_itr_1step_color( &
        q, res,                                          &
        f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), &
        var_comm, aux_comm, Dx, Dy, Lift, mesh2D,        &
        cal_res_flag, 0                                  )
    end if

    var_comm_list(1)%field2d => q
    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )
    return
  end subroutine MGSmoother_Poisson2D_advance_itr_1step

!-- private -----------------------

  !> Update q for elements with a specified color (or evaluate the residual if cal_res_flag = .true.)
!OCL SERIAL
  subroutine Poisson2dGlobal_smoother_advance_itr_1step_color( this, &
    q, res,                                                    &
    f, qx, qy,                                                 &
    var_comm, aux_comm, Dx, Dy, Lift, mesh2D,                  &
    cal_res_flag, color_id )
    implicit none
    class(MGSmoother_Poisson2D), intent(inout) :: this
    type(MeshField2D), intent(inout), target :: q
    type(MeshField2D), intent(inout) :: res
    type(MeshField2D), intent(in) :: f
    type(MeshField2D), intent(inout), target :: qx
    type(MeshField2D), intent(inout), target :: qy
    class(MeshFieldCommBase), intent(inout) :: var_comm
    class(MeshFieldCommBase), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift
    class(MeshBase2D), intent(in), target :: mesh2D
    logical, intent(in) :: cal_res_flag
    integer, intent(in) :: color_id

    integer :: ldomID
    class(LocalMesh2D), pointer :: lmesh2D

    type(MeshFieldContainer) :: var_comm_list(1)
    type(MeshFieldContainer) :: aux_comm_list(2)
    !---------------------------------------------------------------------------

    !- Exchange q
    var_comm_list(1)%field2d => q

    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    !- Element-wise gradient (contravariant components)
    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(ldomID)
      call cal_grad_lc( qx%local(ldomID)%val, qy%local(ldomID)%val, & ! (out)
        q%local(ldomID)%val, Dx, Dy, lmesh2D, lmesh2D%refElem2D     ) ! (in)
    end do

    !- Exchange the gradient as a horizontal vector.
    !  The conversion of vector components between different panels is performed in the communicator.
    aux_comm_list(1)%field2d => qx
    aux_comm_list(2)%field2d => qy

    call aux_comm%Put( aux_comm_list, 1 )
    call aux_comm%Exchange()
    call aux_comm%Get( aux_comm_list, 1 )

    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(ldomID)
      call this%cal_q_lc_elem_PCG( &
        q%local(ldomID)%val, res%local(ldomID)%val,               &
        color_id, cal_res_flag, f%local(ldomID)%val,              &
        qx%local(ldomID)%val, qy%local(ldomID)%val, Dx, Dy, Lift, &
        lmesh2D, lmesh2D%refElem2D                                )
    end do
    return
  end subroutine Poisson2dGlobal_smoother_advance_itr_1step_color

  !> Calculate contravariant components of the element-wise gradient (without lifting terms)
  !! sigma~^i = G^{ij} d q / d xi^j
!OCL SERIAL
  subroutine cal_grad_lc( qx, qy, &
    q, Dx, Dy, lmesh, elem )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qy(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy

    integer :: ke
    real(RP) :: q1(elem%Np)
    real(RP) :: q2(elem%Np)
    !----------------------------
    !$omp parallel do private(q1, q2)
    do ke=lmesh%NeS, lmesh%NeE
      call Sparsemat_matmul(Dx, q(:,ke), q1)
      call Sparsemat_matmul(Dy, q(:,ke), q2)

      ! Covariant components: dq/dalpha, dq/dbeta
      q1(:) = lmesh%Escale(:,ke,1,1) * q1(:)
      q2(:) = lmesh%Escale(:,ke,2,2) * q2(:)

      ! Contravariant components
      qx(:,ke) = lmesh%GIJ(:,ke,1,1) * q1(:) &
               + lmesh%GIJ(:,ke,1,2) * q2(:)
      qy(:,ke) = lmesh%GIJ(:,ke,2,1) * q1(:) &
               + lmesh%GIJ(:,ke,2,2) * q2(:)
    end do
    return
  end subroutine cal_grad_lc

  !> Element-block Gauss-Seidel update with the PCG method for element-local problems
  !! (or evaluation of the residual, res = rhs - A q, if cal_res = .true.)
!OCL SERIAL
  subroutine cal_q_lc_elem_PCG( this, q, res,  &
    color_id, cal_res, rhs, qx, qy,            &
    Dx, Dy, Lift,                              &
    lmesh, elem                                )
    implicit none
    class(MGSmoother_Poisson2D), intent(in) :: this
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(inout) :: q(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: res(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: rhs(elem%Np,lmesh%NeA)   !< Right-hand side multiplied by Gsqrt
    real(RP), intent(in) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qy(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: color_id
    logical, intent(in) :: cal_res

    integer :: i, j, ke
    integer :: parity, parity_max
    integer :: j0
    integer :: j_int
    integer :: ii, ncol

    real(RP) :: Fx(elem%Np), Fy(elem%Np)
    real(RP) :: LiftDFlx(elem%Np)
    real(RP) :: BndFlx(elem%NfpTot,3)
    real(RP) :: gtau(elem%NfpTot)

    real(RP) :: Sig(elem%Np,2)
    real(RP) :: Sig_e(elem%Np,2)
    real(RP) :: Ax(elem%Np)
    real(RP) :: x(elem%Np)
    real(RP) :: r(elem%Np)
    real(RP) :: z(elem%Np)
    real(RP) :: q_(elem%Np)
    real(RP) :: rnorm0, rnorm
    real(RP) :: rz_old, rz_new
    real(RP) :: alpha, beta

    integer :: itr

    real(RP) :: PmatDiag(elem%Np)
    real(RP) :: hk
    real(RP), parameter :: c1 = 1.0_RP
    real(RP), parameter :: c2 = 1.0_RP

    integer :: iM(elem%NfpTot)
    integer :: p,f
    !---------------------------

    call PROF_rapstart( "Poisson2dGlobal_smoother_cal_q_lc", 1 )

    do f=1, elem%Nfaces
    do p=1, elem%Nfp
      i = p + (f-1)*elem%Nfp
      iM(i) = elem%Fmask(p,f)
    end do
    end do

    parity_max = merge(0, 1, cal_res)

    !$omp parallel private( j0, j_int, ncol, i, ke, &
    !$omp gtau, Fx,Fy,BndFlx,LiftDFlx,              &
    !$omp Sig,Sig_e, Ax,x,r,z,q_,                   &
    !$omp rnorm0,rnorm, rz_old,rz_new, alpha,beta,  &
    !$omp itr, PmatDiag, hk )
    do parity=0, parity_max

      if ( cal_res ) then
        j0 = 1; j_int = 1; ncol = lmesh%NeX
      else
        j0 = 1 + modulo(parity-color_id, 2)
        j_int = 2
        ncol = (lmesh%NeX+1-parity)/2
      end if

      !$omp do collapse(2)
      do j=j0, lmesh%NeY, j_int
      do ii=1, ncol
        if (cal_res) then
          i = ii
        else
          i = 2*ii - 1 + parity
        end if
        ke = i + (j-1)*lmesh%NeX

        !- Penalty parameter (tau~ = tau * Gsqrt * G^{ij} n_i n_j on the face nodes)
        call cal_penalty( gtau, &
          ke, iM, lmesh, elem   )

        call eval_Ax( Ax, &
          ke, q, qx, qy, gtau,            &
          Sig, BndFlx, Fx, Fy, LiftDFlx,  &
          Dx, Dy, Lift, iM, lmesh, elem   )

        r(:) = rhs(:,ke) - Ax(:)

        if ( cal_res ) then
          res(:,ke) = r(:)
        else
          rnorm0 = sqrt( sum( r(:) * r(:) ) )

          if ( rnorm0 > 1.0E-16_RP ) then
            ! Precond: physical element width is estimated using the covariant metric
            hk = 4.0_RP / maxval( lmesh%Escale(:,ke,1,1) / sqrt( lmesh%G_ij(:,ke,1,1) ) &
                                + lmesh%Escale(:,ke,2,2) / sqrt( lmesh%G_ij(:,ke,2,2) ) )
            PmatDiag(:) = lmesh%Gsqrt(:,ke) * ( c1 * elem%PolyOrder**2 / hk + c2 / hk**2 )

            x(:) = q(:,ke)
            z(:) = r(:) / PmatDiag(:)

            !-
            q_(:) = z(:)
            rz_old = sum( r(:) * z(:) )

            do itr=1, this%num_internal_itr_max
              !-
              call eval_Ax_itr( Ax, &
                ke, q_, gtau, Sig, Sig_e,  &
                BndFlx, Fx, Fy, LiftDFlx,  &
                Dx, Dy, Lift, iM, lmesh, elem )

              alpha = rz_old / sum(q_(:) * Ax(:))
              x(:) = x(:) + alpha * q_(:)
              r(:) = r(:) - alpha * Ax(:)

              rnorm = sqrt( sum( r(:) * r(:) ) )
              if ( rnorm / rnorm0 < this%residual_threshold_l2_internal_itr ) then
                exit
              else if ( itr == this%num_internal_itr_max ) then
                LOG_INFO("Poisson2dGlobal_smoother_cal_q_lc",*) "ke=", ke, "rnorm/rnorm0=", rnorm / rnorm0, rnorm, rnorm0
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
   !$omp end parallel

    call PROF_rapend( "Poisson2dGlobal_smoother_cal_q_lc", 1 )
    return
  end subroutine cal_q_lc_elem_PCG

  !> Calculate the penalty parameter on the face nodes:
  !! tau~ = C_tau p(p+1)/2 Fscale * Gsqrt * G^{ij} n_i n_j
  !! Since Gsqrt and G^{ij} n_i n_j coincide on both sides of the face (including panel boundaries),
  !! the metric of the own element is used.
!OCL SERIAL
  subroutine cal_penalty( gtau, &
    ke, iM, lmesh, elem )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: gtau(elem%NfpTot)
    integer, intent(in) :: ke
    integer, intent(in) :: iM(elem%NfpTot)

    real(RP) :: nx(elem%NfpTot), ny(elem%NfpTot)
    !------------------------------------------

    nx(:) = lmesh%normal_fn(:,ke,1)
    ny(:) = lmesh%normal_fn(:,ke,2)

    gtau(:) = PENALTY_COEF * 0.5_RP * real( elem%PolyOrder * ( elem%PolyOrder + 1 ), kind=RP ) &
            * lmesh%Fscale(:,ke) * lmesh%Gsqrt(iM,ke)                                  &
            * (            lmesh%GIJ(iM,ke,1,1) * nx(:) * nx(:)                        &
                + 2.0_RP * lmesh%GIJ(iM,ke,1,2) * nx(:) * ny(:)                        &
                +          lmesh%GIJ(iM,ke,2,2) * ny(:) * ny(:)                        )
    return
  end subroutine cal_penalty

!- Evaluation of A q ------------------------------------------------------------

  !> Evaluate (A~ q)_K with the values of neighbor elements
  !!  A~ q = D_i( Gsqrt sigma^i ) + Lift[ Fscale ( F^ - Gsqrt sigma^{i,-} n_i ) ]
!OCL SERIAL
  subroutine eval_Ax( Ax, &
    ke, q, qx, qy, gtau,            &
    Sig, BndFlx, Fx, Fy, LiftDFlx,  &
    Dx, Dy, Lift, iM, lmesh, elem   )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np)
    integer, intent(in) :: ke
    real(RP), intent(in) :: q(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qx(elem%Np,lmesh%NeA)    !< Contravariant element-wise gradient (alpha)
    real(RP), intent(in) :: qy(elem%Np,lmesh%NeA)    !< Contravariant element-wise gradient (beta)
    real(RP), intent(in) :: gtau(elem%NfpTot)
    real(RP), intent(inout) :: Sig(elem%Np,2) ! [work]
    real(RP), intent(inout) :: BndFlx(elem%NfpTot,3) ! [work]
    real(RP), intent(inout) :: Fx(elem%Np) ! [work]
    real(RP), intent(inout) :: Fy(elem%Np) ! [work]
    real(RP), intent(inout) :: LiftDFlx(elem%Np) ! [work]
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: iM(elem%NfpTot)

    real(RP) :: L1(elem%Np), L2(elem%Np)
    !---------------------------------

    !- (a) Lifting terms of the gradient (covariant components)
    call cal_bndflx_1( BndFlx(:,1:2), &
      q, lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2), &
      lmesh%vmapM(:,ke), lmesh%vmapP(:,ke), lmesh, elem    )

    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,1), L1 )
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,2), L2 )

    !- (b) Contravariant components: sigma = sigma~ + G^{ij} L_j
    Sig(:,1) = qx(:,ke) + lmesh%GIJ(:,ke,1,1) * L1(:) + lmesh%GIJ(:,ke,1,2) * L2(:)
    Sig(:,2) = qy(:,ke) + lmesh%GIJ(:,ke,2,1) * L1(:) + lmesh%GIJ(:,ke,2,2) * L2(:)

    !- (c) Divergence
    call cal_bndflx_2( BndFlx(:,3), &
      q, Sig(iM,1), Sig(iM,2), qx, qy, lmesh%Gsqrt(iM,ke),    &
      lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2),       &
      gtau, lmesh%vmapM(:,ke), lmesh%vmapP(:,ke), lmesh, elem )

    call sparsemat_matmul( Dx, lmesh%Gsqrt(:,ke) * Sig(:,1), Fx )
    call sparsemat_matmul( Dy, lmesh%Gsqrt(:,ke) * Sig(:,2), Fy )
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,3), LiftDFlx )
    Ax(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + LiftDFlx(:)
    return
  end subroutine eval_Ax

  !> Numerical flux for the gradient: 1/2 (q^+ - q^-) n_i
!OCL SERIAL
  subroutine cal_bndflx_1( BndFlx, &
    q, nx, ny,          &
    iM, iP, lmesh, elem )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: BndFlx(elem%NfpTot,2)
    real(RP), intent(in) :: q(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot)
    real(RP), intent(in) :: ny(elem%NfpTot)
    integer, intent(in) :: iM(elem%NfpTot)
    integer, intent(in) :: iP(elem%NfpTot)

    real(RP) :: dq(elem%NfpTot)
    !----------------------------------------------------

    dq(:) = q(iP) - q(iM)
    BndFlx(:,1) = 0.5_RP * dq(:) * nx(:)
    BndFlx(:,2) = 0.5_RP * dq(:) * ny(:)
    return
  end subroutine cal_bndflx_1

  !> Numerical flux for the divergence: F^ - Gsqrt sigma^{i,-} n_i, where
  !!  F^ = 1/2 Gsqrt ( sigma~^{i,-} + sigma~^{i,+} ) n_i + 1/2 tau~ ( q^+ - q^- )
!OCL SERIAL
  subroutine cal_bndflx_2( BndFlx, &
    q, qx_M, qy_M, qx_e, qy_e, Gsqrt_M, &
    nx, ny, gtau,                       &
    iM, iP, lmesh, elem )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: BndFlx(elem%NfpTot)
    real(RP), intent(in) :: q(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: qx_M(elem%NfpTot)
    real(RP), intent(in) :: qy_M(elem%NfpTot)
    real(RP), intent(in) :: qx_e(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: qy_e(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: Gsqrt_M(elem%NfpTot)
    real(RP), intent(in) :: nx(elem%NfpTot)
    real(RP), intent(in) :: ny(elem%NfpTot)
    real(RP), intent(in) :: gtau(elem%NfpTot)
    integer, intent(in) :: iM(elem%NfpTot)
    integer, intent(in) :: iP(elem%NfpTot)

    real(RP) :: dq(elem%NfpTot)
    !----------------------------------------------------

    dq(:) = q(iP) - q(iM)
    BndFlx(:) = Gsqrt_M(:) * (                                     &
          ( 0.5_RP * ( qx_e(iP) + qx_e(iM) ) - qx_M(:) ) * nx(:)   &
        + ( 0.5_RP * ( qy_e(iP) + qy_e(iM) ) - qy_M(:) ) * ny(:) ) &
      + 0.5_RP * gtau(:) * dq(:)
    return
  end subroutine cal_bndflx_2

  !> Evaluate (A~_KK q)_K for the element-local problem (values of neighbor elements are set to zero)
!OCL SERIAL
  subroutine eval_Ax_itr( Ax, &
    ke, q, gtau,                         &
    Sig, Sig_e, BndFlx, Fx, Fy, LiftDFlx, &
    Dx, Dy, Lift, iM, lmesh, elem        )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np)
    integer, intent(in) :: ke
    real(RP), intent(in) :: q(elem%Np)
    real(RP), intent(in) :: gtau(elem%NfpTot)
    real(RP), intent(inout) :: Sig(elem%Np,2) ! [work]
    real(RP), intent(inout) :: Sig_e(elem%Np,2) ! [work]
    real(RP), intent(inout) :: BndFlx(elem%NfpTot,3) ! [work]
    real(RP), intent(inout) :: Fx(elem%Np) ! [work]
    real(RP), intent(inout) :: Fy(elem%Np) ! [work]
    real(RP), intent(inout) :: LiftDFlx(elem%Np) ! [work]
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: iM(elem%NfpTot)

    real(RP) :: L1(elem%Np), L2(elem%Np)
    !---------------------------------

    call cal_bndflx_itr_1( BndFlx(:,1:2), &
      q(iM),                                                        &
      lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2), lmesh, elem )

    !- Element-wise gradient (covariant -> contravariant)
    call sparsemat_matmul( Dx, q(:), Fx )
    call sparsemat_matmul( Dy, q(:), Fy )
    Fx(:) = lmesh%Escale(:,ke,1,1) * Fx(:)
    Fy(:) = lmesh%Escale(:,ke,2,2) * Fy(:)

    Sig_e(:,1) = lmesh%GIJ(:,ke,1,1) * Fx(:) + lmesh%GIJ(:,ke,1,2) * Fy(:)
    Sig_e(:,2) = lmesh%GIJ(:,ke,2,1) * Fx(:) + lmesh%GIJ(:,ke,2,2) * Fy(:)

    !- Lifting terms
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,1), L1 )
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,2), L2 )

    Sig(:,1) = Sig_e(:,1) + lmesh%GIJ(:,ke,1,1) * L1(:) + lmesh%GIJ(:,ke,1,2) * L2(:)
    Sig(:,2) = Sig_e(:,2) + lmesh%GIJ(:,ke,2,1) * L1(:) + lmesh%GIJ(:,ke,2,2) * L2(:)

    !- Divergence
    call cal_bndflx_itr_2( BndFlx(:,3), &
      q(iM), Sig(iM,1), Sig(iM,2), Sig_e(iM,1), Sig_e(iM,2), lmesh%Gsqrt(iM,ke), &
      lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2), gtau, lmesh, elem        )

    call sparsemat_matmul( Dx, lmesh%Gsqrt(:,ke) * Sig(:,1), Fx )
    call sparsemat_matmul( Dy, lmesh%Gsqrt(:,ke) * Sig(:,2), Fy )
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,3), LiftDFlx )
    Ax(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + LiftDFlx(:)
    return
  end subroutine eval_Ax_itr

!OCL SERIAL
  subroutine cal_bndflx_itr_1( BndFlx, &
    q_M, nx, ny, lmesh, elem )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: BndFlx(elem%NfpTot,2)
    real(RP), intent(in) :: q_M(elem%NfpTot)
    real(RP), intent(in) :: nx(elem%NfpTot)
    real(RP), intent(in) :: ny(elem%NfpTot)

    real(RP) :: dq(elem%NfpTot)
    !----------------------------------------------------

    dq(:) =      - q_M(:)
    BndFlx(:,1) = 0.5_RP * dq(:) * nx(:)
    BndFlx(:,2) = 0.5_RP * dq(:) * ny(:)
    return
  end subroutine cal_bndflx_itr_1

!OCL SERIAL
  subroutine cal_bndflx_itr_2( BndFlx, &
    q, qx, qy, qx_e, qy_e, Gsqrt_M, &
    nx, ny, gtau, lmesh, elem )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: BndFlx(elem%NfpTot)
    real(RP), intent(in) :: q(elem%NfpTot)
    real(RP), intent(in) :: qx(elem%NfpTot)
    real(RP), intent(in) :: qy(elem%NfpTot)
    real(RP), intent(in) :: qx_e(elem%NfpTot)
    real(RP), intent(in) :: qy_e(elem%NfpTot)
    real(RP), intent(in) :: Gsqrt_M(elem%NfpTot)
    real(RP), intent(in) :: nx(elem%NfpTot)
    real(RP), intent(in) :: ny(elem%NfpTot)
    real(RP), intent(in) :: gtau(elem%NfpTot)

    real(RP) :: dq(elem%NfpTot)
    !----------------------------------------------------

    dq(:) =      - q(:)
    BndFlx(:) = Gsqrt_M(:) * (                                  &
          ( 0.5_RP * (          + qx_e(:) ) - qx(:) ) * nx(:)   &
        + ( 0.5_RP * (          + qy_e(:) ) - qy(:) ) * ny(:) ) &
      + 0.5_RP * gtau(:) * dq(:)
    return
  end subroutine cal_bndflx_itr_2

end module mod_poisson2dGlobal_smoother
