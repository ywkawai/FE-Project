#include "scaleFElib.h"
module mod_poisson2d_smoother
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
  
  use scale_element_base, only: ElementBase1D, ElementBase2D
  use scale_element_line, only: LineElement

  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D

  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_base, only: MeshFieldCommBase
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

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

    integer :: elem_solver_type_id  !< Element solver type ID

    ! For direct solver
    type(LineElement), allocatable :: elem1D_list(:)
  contains
    procedure :: Init => MGSmoother_Poisson2D_Init
    procedure :: Final => MGSmoother_Poisson2D_Final
    procedure :: Advance_itr_1step => MGSmoother_Poisson2D_advance_itr_1step
    procedure, private :: advance_itr_1step_color => Poisson2d_smoother_advance_itr_1step_color
    procedure, private :: cal_q_lc_elem_Direct => cal_q_lc_elem_direct
    procedure, private :: cal_q_lc_elem_PCG => cal_q_lc_elem_PCG
  end type MGSmoother_Poisson2D

  ! public :: poisson_smoother_evaluate_error_norm
  integer, parameter, public :: MGSmoother_Possion2D_AUX_SCALAR_NUM  = 0
  integer, parameter, public :: MGSmoother_Possion2D_AUX_VEC_NUM     = 1

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: cal_grad_lc

  integer, parameter :: AUXVAR_QX_ID = 1
  integer, parameter :: AUXVAR_QY_ID = 2

  integer, parameter :: ELEM_SOLVER_DIRECT_ID = 1
  integer, parameter :: ELEM_SOLVER_PCG_ID    = 2
  
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

    character(len=H_SHORT) :: elem_solver = "PCG" !< Element solver type ("DIRECT" or "PCG")

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

    integer :: plev
    integer :: porder_list_size
    !---------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF, nml=PARAM_MGSmoother, iostat=ierr)
    if ( ierr < 0 ) then
      LOG_INFO("MGSmoother_Poisson2D_Init",*) "Not found namelist PARAM_MGSmoother in the configuration file."
    else if( ierr > 0 ) then
      LOG_INFO("MGSmoother_Poisson2D_Init",*) "Error reading namelist PARAM_MGSmoother in the configuration file."
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
    case ( "DIRECT" )
      this%elem_solver_type_id = ELEM_SOLVER_DIRECT_ID
      ! Allocate 1D element for direct solver
      porder_list_size = size( porder_list )
      allocate( this%elem1D_list( porder_list_size ) )
      do plev=1, porder_list_size
        call this%elem1D_list(plev)%Init( porder_list(plev), .false. )
      end do
    case ( "PCG" )
      this%elem_solver_type_id = ELEM_SOLVER_PCG_ID
    case default
      LOG_INFO("MGSmoother_Poisson2D_Init",*) "Unknown elem_solver type: "//trim(elem_solver)
      call PRC_abort
    end select

    return
  end subroutine MGSmoother_Poisson2D_Init

  !> Finalization
  subroutine MGSmoother_Poisson2D_Final(this)
    use scale_multigrid_smoother_base, only: MGSmootherBase2D_Final
    implicit none
    class(MGSmoother_Poisson2D), intent(inout) :: this

    integer :: plev
    !---------------------------------------------------------------------------

    if ( this%elem_solver_type_id == ELEM_SOLVER_DIRECT_ID ) then
      if ( allocated(this%elem1D_list) ) then
        do plev=1, size(this%elem1D_list)
          call this%elem1D_list(plev)%Final()
        end do
        deallocate( this%elem1D_list )
      end if
    end if

    call MGSmootherBase2D_Final( this )
    return
  end subroutine MGSmoother_Poisson2D_Final

  !> Smoother
!OCL SERIAL
  subroutine MGSmoother_Poisson2D_advance_itr_1step( this, &
    q, res,                                    & 
    f, aux_var,                                &
    itr, var_comm, aux_comm, Dx, Dy, Lift, mesh2D,   &
    cal_res_flag, zero_initial_guess,            &
    mg_p_level, mg_h_level, pre_or_post_smooth )

    use scale_mesh_hierarchy_base, only: &
      pMG_FINEST_LEVEL => MESH_HIERARCHY_pMG_FINEST_LEVEL
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
    logical :: is_mg_top_level
    integer :: color_id_offset
    !---------------------------------------------------------------------------

    is_mg_top_level = ( mg_p_level == pMG_FINEST_LEVEL )

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
      call PROF_rapstart( "Poisson2d_smoother_res", 1 )
      call this%advance_itr_1step_color( &
        q, res,                                          & 
        f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), &
        itr, var_comm, aux_comm, Dx, Dy, Lift, mesh2D,   &
        .true., 0, is_mg_top_level, mg_p_level           )
      
      call this%Set_initial_residual( res )
      call flush(IO_FID_LOG)
      call PROF_rapend( "Poisson2d_smoother_res", 1 )
    end if

    call PROF_rapstart( "Poisson2d_smoother_rb", 1 )
    ! Red
    call this%advance_itr_1step_color( &
      q, res,                                                      & 
      f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID),             &
      itr, var_comm, aux_comm, Dx, Dy, Lift, mesh2D,               &
      .false., mod(color_id_offset,2), is_mg_top_level, mg_p_level )

    ! Black
    call this%advance_itr_1step_color( &
      q, res,                                                        & 
      f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID),               &
      itr, var_comm, aux_comm, Dx, Dy, Lift, mesh2D,                 &
      .false., mod(color_id_offset+1,2), is_mg_top_level, mg_p_level )      
    call PROF_rapend( "Poisson2d_smoother_rb", 1 )

    if ( cal_res_flag ) then
      call this%advance_itr_1step_color( &
        q, res,                                          & 
        f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), &
        itr, var_comm, aux_comm, Dx, Dy, Lift, mesh2D,   &
        cal_res_flag, 0, is_mg_top_level, mg_p_level     )
    end if
    
    var_comm_list(1)%field2d => q
    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )
    return
  end subroutine MGSmoother_Poisson2D_advance_itr_1step

!-- private -----------------------

!OCL SERIAL
  subroutine Poisson2d_smoother_advance_itr_1step_color( this, &
    q, res,                                                    & 
    f, qx, qy,                                                 &
    itr, var_comm, aux_comm, Dx, Dy, Lift, mesh2D,             &
    cal_res_flag, color_id, is_mg_top_level, mg_p_level )
    implicit none
    class(MGSmoother_Poisson2D), intent(inout) :: this
    type(MeshField2D), intent(inout), target :: q
    type(MeshField2D), intent(inout) :: res
    type(MeshField2D), intent(in) :: f
    type(MeshField2D), intent(inout), target :: qx
    type(MeshField2D), intent(inout), target :: qy
    integer, intent(in) :: itr
    class(MeshFieldCommBase), intent(inout) :: var_comm
    class(MeshFieldCommBase), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift
    class(MeshBase2D), intent(in), target :: mesh2D
    logical, intent(in) :: cal_res_flag
    integer, intent(in) :: color_id
    logical, intent(in) :: is_mg_top_level
    integer, intent(in) :: mg_p_level

    integer :: ldomID
    class(LocalMesh2D), pointer :: lmesh2D

    type(MeshFieldContainer) :: var_comm_list(1)
    type(MeshFieldContainer) :: aux_comm_list(2)
    !---------------------------------------------------------------------------

    !-
    var_comm_list(1)%field2d => q

    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(ldomID)
      call cal_grad_lc( qx%local(ldomID)%val, qy%local(ldomID)%val,       & ! (out)
        color_id, q%local(ldomID)%val, Dx, Dy, lmesh2D, lmesh2D%refElem2D ) ! (in)
    end do

    !- 
    aux_comm_list(1)%field2d => qx
    aux_comm_list(2)%field2d => qy

    call aux_comm%Put( aux_comm_list, 1 )
    call aux_comm%Exchange()
    call aux_comm%Get( aux_comm_list, 1 )

    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(ldomID)

      select case ( this%elem_solver_type_id )
      case ( ELEM_SOLVER_DIRECT_ID )
        call this%cal_q_lc_elem_Direct( &
          q%local(ldomID)%val, res%local(ldomID)%val,                   &
          color_id, cal_res_flag, is_mg_top_level, f%local(ldomID)%val, &
          qx%local(ldomID)%val, qy%local(ldomID)%val,                   &
          lmesh2D%VMapM, lmesh2D%VMapP, lmesh2D, lmesh2D%refElem2D,     &
          this%elem1D_list(mg_p_level) )
      case ( ELEM_SOLVER_PCG_ID )
        call this%cal_q_lc_elem_PCG( &
          q%local(ldomID)%val, res%local(ldomID)%val,                   &
          color_id, cal_res_flag, is_mg_top_level, f%local(ldomID)%val, &
          qx%local(ldomID)%val, qy%local(ldomID)%val, Dx, Dy, Lift,     &
          lmesh2D%VMapM, lmesh2D%VMapP, lmesh2D, lmesh2D%refElem2D )
      end select
    end do
    return
  end subroutine Poisson2d_smoother_advance_itr_1step_color

!- private
  subroutine cal_grad_lc( qx, qy, &
    color_id, q, Dx, Dy, lmesh, elem )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qy(elem%Np,lmesh%NeA)
    integer, intent(in) :: color_id
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

      q1(:) = lmesh%Escale(:,ke,1,1) * q1(:)
      q2(:) = lmesh%Escale(:,ke,2,2) * q2(:)

      qx(:,ke) = lmesh%G_ij(:,ke,1,1) * q1(:) &
              +  lmesh%G_ij(:,ke,2,1) * q2(:)
      qy(:,ke) = lmesh%G_ij(:,ke,1,2) * q1(:) &
              +  lmesh%G_ij(:,ke,2,2) * q2(:)
    end do
    return
  end subroutine cal_grad_lc

!OCL SERIAL
  subroutine cal_q_lc_elem_PCG( this, q, res,        &
    color_id, cal_res, is_mg_top_level, rhs, qx, qy, &
    Dx, Dy, Lift,                                    &
    vmapM, vmapP, lmesh, elem                        )
    use scale_prc
    use scale_linalgebra, only: Linalgebra_SolveLinEq
    implicit none
    class(MGSmoother_Poisson2D), intent(in) :: this
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(inout) :: q(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: res(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: rhs(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qy(elem%Np,lmesh%NeA)    
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: color_id
    logical, intent(in) :: cal_res
    logical, intent(in) :: is_mg_top_level
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)

    integer :: i, j, ke
    integer :: parity, parity_max
    integer :: j0
    integer :: j_int
    integer :: ii, ncol

    real(RP) :: Fx(elem%Np), Fy(elem%Np)
    real(RP) :: LiftDFlx(elem%Np)
    real(RP) :: BndFlx(elem%NfpTot,3)
    real(RP) :: gtau(elem%NfpTot)

    real(RP) :: Grad_q(elem%Np,2)
    real(RP) :: Grad_q_e(elem%Np,2)
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

    call PROF_rapstart( "Poisson3d_smoother_cal_q_lc", 1 )

    do f=1, elem%Nfaces
    do p=1, elem%Nfp
      i = p + (f-1)*elem%Nfp
      iM(i) = elem%Fmask(p,f)
    end do
    end do

    parity_max = merge(0, 1, cal_res)

    !$omp parallel private( j0, j_int, ncol, i, ke, &
    !$omp gtau, Fx,Fy,BndFlx,LiftDFlx,              &
    !$omp Grad_q,Grad_q_e, Ax,x,r,z,q_,             &
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
        
        !-
        gtau(:) = lmesh%Fscale(:,ke) * elem%PolyOrder * ( elem%PolyOrder + 1 ) * 0.5_RP * 1.2_RP

        call eval_Ax( Ax, &
          ke, q, qx, qy, gtau,              &
          Grad_q, BndFlx, Fx, Fy, LiftDFlx, &
          Dx, Dy, Lift, iM, lmesh, elem     )

        r(:) = rhs(:,ke) - Ax(:)

        if ( cal_res ) then
          res(:,ke) = r(:)
        else
          rnorm0 = sqrt( sum( r(:) * r(:) ) )

          if ( rnorm0 > 1.0E-16_RP ) then
            ! Precond
            hk = 4.0_RP / maxval(lmesh%Escale(:,ke,1,1) + lmesh%Escale(:,ke,2,2))
            PmatDiag(:) = c1 * elem%PolyOrder**2 / hk + c2 / hk**2

            x(:) = q(:,ke)
            z(:) = r(:) / PmatDiag(:)

            !-
            q_(:) = z(:)
            rz_old = sum( r(:) * z(:) )

            do itr=1, this%num_internal_itr_max
              !-
              call eval_Ax_itr( Ax, &
                ke, q_, gtau, Grad_q, Grad_q_e, &
                BndFlx, Fx, Fy, LiftDFlx,       &
                Dx, Dy, Lift, iM, lmesh, elem   )

              alpha = rz_old / sum(q_(:) * Ax(:))
              x(:) = x(:) + alpha * q_(:)
              r(:) = r(:) - alpha * Ax(:)

              rnorm = sqrt( sum( r(:) * r(:) ) )
              if ( rnorm / rnorm0 < this%residual_threshold_l2_internal_itr ) then
                exit
              else if ( itr == this%num_internal_itr_max ) then
                LOG_INFO("Poisson2d_smoother_cal_q_lc",*) "ke=", ke, "rnorm/rnorm0=", rnorm / rnorm0, rnorm, rnorm0
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

    call PROF_rapend( "Poisson3d_smoother_cal_q_lc", 1 )
    return
  end subroutine cal_q_lc_elem_PCG

!- Preconitioned Conjugate Gradient method in internal iteration --------------
!OCL SERIAL
  subroutine eval_Ax( Ax, &
    ke, q, qx, qy, gtau,              &
    Grad_q, BndFlx, Fx, Fy, LiftDFlx, &
    Dx, Dy, Lift, iM, lmesh, elem     )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np)
    integer, intent(in) :: ke
    real(RP), intent(in) :: q(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qy(elem%Np,lmesh%NeA)    
    real(RP), intent(in) :: gtau(elem%NfpTot)
    real(RP), intent(inout) :: Grad_q(elem%Np,2) ! [work]
    real(RP), intent(inout) :: BndFlx(elem%NfpTot,3) ! [work]
    real(RP), intent(inout) :: Fx(elem%Np) ! [work]
    real(RP), intent(inout) :: Fy(elem%Np) ! [work]
    real(RP), intent(inout) :: LiftDFlx(elem%Np) ! [work]
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: iM(elem%NfpTot)
    !---------------------------------

    call cal_bndflx_1( BndFlx(:,1:2), &
      q, lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2), &
      lmesh%vmapM(:,ke), lmesh%vmapP(:,ke), lmesh, elem    )
    
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,1), LiftDFlx )
    Grad_q(:,1) = qx(:,ke) + LiftDFlx(:)

    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,2), LiftDFlx )
    Grad_q(:,2) = qy(:,ke) + LiftDFlx(:)
    
    !-
    call cal_bndflx_2( BndFlx(:,3), &
      q, Grad_q(iM,1), Grad_q(iM,2), qx, qy,                  &
      lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2),       &
      gtau, lmesh%vmapM(:,ke), lmesh%vmapP(:,ke), lmesh, elem )

    call sparsemat_matmul( Dx, Grad_q(:,1), Fx )
    call sparsemat_matmul( Dy, Grad_q(:,2), Fy )
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,3), LiftDFlx )
    Ax(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + LiftDFlx(:)
    return
  end subroutine eval_Ax

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

!OCL SERIAL
  subroutine cal_bndflx_2( BndFlx, &
    q, qx_M, qy_M, qx_e, qy_e,     &
    nx, ny, gtau,                  &
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
    real(RP), intent(in) :: nx(elem%NfpTot)
    real(RP), intent(in) :: ny(elem%NfpTot)
    real(RP), intent(in) :: gtau(elem%NfpTot)
    integer, intent(in) :: iM(elem%NfpTot)
    integer, intent(in) :: iP(elem%NfpTot)

    real(RP) :: dq(elem%NfpTot)
    !----------------------------------------------------

    dq(:) = q(iP) - q(iM)
    BndFlx(:) = &
        ( 0.5_RP * ( qx_e(iP) + qx_e(iM) ) - qx_M(:) ) * nx(:) &
      + ( 0.5_RP * ( qy_e(iP) + qy_e(iM) ) - qy_M(:) ) * ny(:) &
      + 0.5_RP * gtau(:) * dq(:)
    return
  end subroutine cal_bndflx_2  

!OCL SERIAL
  subroutine eval_Ax_itr( Ax, &
    ke, q, gtau,                                &
    Grad_q, Grad_q_e, BndFlx, Fx, Fy, LiftDFlx, &
    Dx, Dy, Lift, iM, lmesh, elem               )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np)
    integer, intent(in) :: ke
    real(RP), intent(in) :: q(elem%Np)
    real(RP), intent(in) :: gtau(elem%NfpTot)
    real(RP), intent(inout) :: Grad_q(elem%Np,2) ! [work]
    real(RP), intent(inout) :: Grad_q_e(elem%Np,2) ! [work]
    real(RP), intent(inout) :: BndFlx(elem%NfpTot,3) ! [work]
    real(RP), intent(inout) :: Fx(elem%Np) ! [work]
    real(RP), intent(inout) :: Fy(elem%Np) ! [work]
    real(RP), intent(inout) :: LiftDFlx(elem%Np) ! [work]
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: iM(elem%NfpTot)
    !---------------------------------

    call cal_bndflx_itr_1( BndFlx(:,1:3), &
      q(iM),                                                        &
      lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2), lmesh, elem )

    call sparsemat_matmul( Dx, q(:), Fx )
    call sparsemat_matmul( Dy, q(:), Fy )

    Grad_q_e(:,1) = lmesh%Escale(:,ke,1,1) * Fx(:) 
    Grad_q_e(:,2) = lmesh%Escale(:,ke,2,2) * Fy(:) 

    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,1), LiftDFlx )
    Grad_q(:,1) = Grad_q_e(:,1) + LiftDFlx(:)

    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,2), LiftDFlx )
    Grad_q(:,2) = Grad_q_e(:,2) + LiftDFlx(:)

    !-
    call cal_bndflx_itr_2( BndFlx(:,3), &
      q(iM), Grad_q(iM,1), Grad_q(iM,2), Grad_q_e(iM,1), Grad_q_e(iM,2),  &
      lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2), gtau, lmesh, elem )
    
    call sparsemat_matmul( Dx, Grad_q(:,1), Fx )
    call sparsemat_matmul( Dy, Grad_q(:,2), Fy )
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
    q, qx, qy, qx_e, qy_e,    &
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
    real(RP), intent(in) :: nx(elem%NfpTot)
    real(RP), intent(in) :: ny(elem%NfpTot)
    real(RP), intent(in) :: gtau(elem%NfpTot)

    real(RP) :: dq(elem%NfpTot)
    !----------------------------------------------------

    dq(:) =      - q(:)
    BndFlx(:) = &
        ( 0.5_RP * (          + qx_e(:) ) - qx(:) ) * nx(:) &
      + ( 0.5_RP * (          + qy_e(:) ) - qy(:) ) * ny(:) &
      + 0.5_RP * gtau(:) * dq(:)
    return
  end subroutine cal_bndflx_itr_2  
  
!- Direct solver within an elemnt -----------------------------

!OCL SERIAL
  subroutine cal_q_lc_elem_Direct( this, q, res,     &
    color_id, cal_res, is_mg_top_level, rhs, qx, qy, &
    vmapM, vmapP, lmesh, elem, elem1D                )
    use scale_linalgebra, only: Linalgebra_SolveLinEq
    implicit none
    class(MGSmoother_Poisson2D), intent(inout) :: this
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    class(ElementBase1D), intent(in) :: elem1D
    real(RP), intent(inout) :: q(elem%Np*lmesh%NeA)
    real(RP), intent(out) :: res(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: rhs(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qx(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: qy(elem%Np*lmesh%NeA)    
    integer, intent(in) :: color_id
    logical, intent(in) :: cal_res
    logical, intent(in) :: is_mg_top_level
    integer, intent(in) :: vmapM(elem%Nfp,elem%Nfaces,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%Nfp,elem%Nfaces,lmesh%Ne)

    integer :: ke, ke2
    integer :: i, j

    integer :: p

    integer :: f1, f2

    real(RP) :: E11, E22

    real(RP) :: Dx(elem%Np,elem%Np)
    real(RP) :: Dy(elem%Np,elem%Np)
    real(RP) :: Sx_tr0(elem%Np,elem%Np)
    real(RP) :: Sy_tr0(elem%Np,elem%Np)
    real(RP) :: Sx_tr(elem%Np,elem%Np)
    real(RP) :: Sy_tr(elem%Np,elem%Np)

    real(RP) :: GsqrtDx(elem%Np,elem%Np)
    real(RP) :: GsqrtDy(elem%Np,elem%Np)
    real(RP) :: GsqrtI(elem%Np,elem%Np)

    real(RP) :: Dn1(elem%Np,elem%Np)    
    real(RP) :: GsqrtDn1(elem%Np,elem%Np)

    real(RP) :: OPT11(elem%Np,elem%Np)
    real(RP) :: mmE(elem%Np,elem%Np)   
    real(RP) :: Lift_(elem%Np,elem%Np)
    real(RP) :: Lift_x(elem%Np,elem%Np)       
    real(RP) :: Lift_y(elem%Np,elem%Np)
    real(RP) :: Lift_x_q(elem%Np)       
    real(RP) :: Lift_y_q(elem%Np)
    real(RP) :: Msrc(elem%Np)

    real(RP) :: massEdge(elem%Np,elem%Np,elem%Nfaces)
    real(RP) :: Emat(elem%Np,elem%Nfp*elem%Nfaces)
    integer :: Fm1(elem%Nfp)
    integer :: id
    real(RP) :: lnx, lny  

    real(RP) :: q_(elem%Np)
    real(RP) :: q_P(elem%Nfp)
    real(RP) :: qx_P(elem%Nfp)
    real(RP) :: qy_P(elem%Nfp)
    real(RP) :: GsqrtDn2(elem%Nfp)

    real(RP) :: gtau_

    integer :: parity, parity_max
    integer :: j0
    integer :: j_int
    integer :: ii, ncol

    real(RP) :: q_old
    real(RP), parameter :: omg = 1.0_RP

    real(RP) :: MassMat(elem%Np,elem%Np)
    real(RP) :: S1Dtr_D1D(elem1D%Np,elem1D%Np)
    !---------------------------

    call PROF_rapstart( "Poisson2d_smoother_cal_q_lc", 1 )

    call PROF_rapstart( "Poisson2d_smoother_cal_q_lc_0", 1 )

    Emat(:,:) = matmul(elem%M, elem%Lift)
    massEdge(:,:,:)  = 0.0_RP
    do f1=1, elem%Nfaces
      Fm1(:) = elem%Fmask(:,f1)
      massEdge(Fm1,Fm1,f1) = Emat(Fm1, (f1-1)*elem%Nfp+1:f1*elem%Nfp)
    end do

    Sx_tr0(:,:) = matmul(transpose(elem%Dx1),elem%M)
    Sy_tr0(:,:) = matmul(transpose(elem%Dx2),elem%M)

    MassMat(:,:) = elem%M(:,:)

    S1Dtr_D1D(:,:) = matmul(transpose(elem1D%Dx1), elem1D%M)
    S1Dtr_D1D(:,:) = matmul(S1Dtr_D1D,elem1D%Dx1)

    parity_max = merge(0, 1, cal_res)
    call PROF_rapend( "Poisson2d_smoother_cal_q_lc_0", 1 )

    !$omp parallel private( j0, j_int, ncol, i, ke, &
    !$omp E11, E22, Dx, Dy, Sx_tr, Sy_tr, GsqrtI, GsqrtDx, GsqrtDy,                       &
    !$omp Dn1, GsqrtDn1, GsqrtDn2, OPT11, mmE, Lift_, Lift_x, Lift_y, Lift_x_q, Lift_y_q, &
    !$omp Msrc, f1, ke2, f2, id, lnx, lny, Fm1, q_, q_P, qx_P, qy_P, gtau_, p, q_old      )

    do parity=0, parity_max

      if (cal_res) then
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
        
        ! Assume E11 and E22 are constant in an element
        E11 = lmesh%Escale(1,ke,1,1)
        E22 = lmesh%Escale(1,ke,2,2)

        Dx(:,:) = E11 * elem%Dx1
        Dy(:,:) = E22 * elem%Dx2
        Sx_tr(:,:) = E11 * Sx_tr0(:,:)
        Sy_tr(:,:) = E22 * Sy_tr0(:,:)

        do p=1, elem%Np
          GsqrtI(:,p) = lmesh%Gsqrt(:,ke)
          GsqrtDx(:,p) = GsqrtI(:,p) * ( lmesh%GIJ(:,ke,1,1) * Dx(:,p) + lmesh%GIJ(:,ke,1,2) * Dy(:,p) )
          GsqrtDy(:,p) = GsqrtI(:,p) * ( lmesh%GIJ(:,ke,2,1) * Dx(:,p) + lmesh%GIJ(:,ke,2,2) * Dy(:,p) )
        end do

        Msrc(:) = matmul(elem%M, rhs(:,ke) * lmesh%Gsqrt(:,ke) )

        call construct_OPT11_elem( OPT11, &
          S1Dtr_D1D, elem1D%M, E11, E22, lmesh%Gsqrt(1,ke), elem1D%Np )

        do f1=1, 4        
          ke2 = lmesh%EToE(ke,f1); f2 = lmesh%EToF(ke,f1)
          
          id = 1 + (f1-1)*elem%Nfp
          lnx = lmesh%normal_fn(id,ke,1)
          lny = lmesh%normal_fn(id,ke,2)
          Fm1(:) = elem%Fmask(:,f1)
                    
          mmE(:,:) = lmesh%Fscale(id,ke) * massEdge(:,:,f1)  
                          
          Lift_(:,:) = matmul(elem%invM, mmE)
          do p=1, elem%Np
            Lift_x(:,p) = lmesh%Gsqrt(:,ke) * ( lnx * lmesh%GIJ(:,ke,1,1) + lny * lmesh%GIJ(:,ke,1,2) ) * Lift_(:,p)
            Lift_y(:,p) = lmesh%Gsqrt(:,ke) * ( lnx * lmesh%GIJ(:,ke,2,1) + lny * lmesh%GIJ(:,ke,2,2) ) * Lift_(:,p)
          end do
          
          !-
          gtau_ = elem%PolyOrder * (elem%PolyOrder + 1) * lmesh%Fscale(id,ke) * 0.5_RP * 1.2E0_RP

          !-
          Dn1(:,:) = lnx * Dx(:,:) + lny * Dy(:,:)
          GsqrtDn1(:,:) = lnx * GsqrtDx(:,:) + lny * GsqrtDy(:,:)

          OPT11(:,:) = OPT11(:,:) + 0.5_RP * ( &
            - gtau_ * mmE(:,:) + matmul(mmE, GsqrtDn1)        &
            + ( matmul(Sx_tr,Lift_x) + matmul(Sy_tr,Lift_y) ) )
          
          q_P(:) = q(VMapP(:,f1,ke))
          Lift_x_q(:) = matmul(Lift_x(:,Fm1), q_P(:))
          Lift_y_q(:) = matmul(Lift_y(:,Fm1), q_P(:))
          Msrc(:) = Msrc(:) + 0.5_RP * ( &
            + matmul(Sx_tr,Lift_x_q) + matmul(Sy_tr,Lift_y_q) )

          qx_P(:) = qx(VMapP(:,f1,ke))
          qy_P(:) = qy(VMapP(:,f1,ke))
          GsqrtDn2(:) = lmesh%Gsqrt(Fm1(:),ke) * ( lnx * qx_P(:) + lny * qy_P(:) )
          Msrc(Fm1) = Msrc(Fm1) - 0.5_RP * (  &
            matmul(mmE(Fm1,Fm1), gtau_ * q_P(:) + GsqrtDn2(:)) )
        end do 
        
        if ( cal_res ) then
          do p=1, elem%Np
            q_(p) = q(p+(ke-1)*elem%Np)
          end do
          res(:,ke) = Msrc(:) - matmul(OPT11(:,:), q_(:))
          res(:,ke) = matmul(elem%invM, res(:,ke))
        else
          call LinAlgebra_SolveLinEq( OPT11, Msrc, q_ )
          do p=1, elem%Np
            q_old = q(p+(ke-1)*elem%Np)
            q(p+(ke-1)*elem%Np) = q_old + omg * ( q_(p) - q_old )
          end do
        end if        
      end do
      end do
    end do
    !$omp end parallel

    call PROF_rapend( "Poisson2d_smoother_cal_q_lc", 1 )
    return
  end subroutine cal_q_lc_elem_Direct

!OCL SERIAL
  subroutine construct_OPT11_elem( OPT11, &
    S1Dtr_D1D, M1D,       &
    E11, E22, Gsqrt, Np1D )
    implicit none
    integer, intent(in) :: Np1D
    real(RP), intent(out) :: OPT11(Np1D,Np1D,Np1D,Np1D)
    real(RP), intent(in) :: S1Dtr_D1D(Np1D,Np1D)
    real(RP), intent(in) :: M1D(Np1D,Np1D)
    real(RP), intent(in) :: E11
    real(RP), intent(in) :: E22
    real(RP), intent(in) :: Gsqrt

    integer :: i, j
    integer :: l, m

    real(RP) :: coef
    real(RP) :: tmp
    !-----------------------------------------------

    coef = - E11**2 * Gsqrt
    do m=1, Np1D; do l=1, Np1D
      do j=1, Np1D
        tmp = coef * M1D(j,m)
        do i=1, Np1D
          OPT11(i,j,l,m) = S1Dtr_D1D(i,l) * tmp
        end do
      end do
    end do; end do

    coef = - E22**2 * Gsqrt
    do m=1, Np1D; do l=1, Np1D
      do j=1, Np1D
        tmp = coef * S1Dtr_D1D(j,m)
        do i=1, Np1D
          OPT11(i,j,l,m) = OPT11(i,j,l,m) + M1D(i,l) * tmp
        end do
      end do
    end do; end do
    return
  end subroutine construct_OPT11_elem

end module mod_poisson2d_smoother
