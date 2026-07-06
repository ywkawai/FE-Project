
#include "scaleFElib.h"
module mod_poisson3d_smoother
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
  
  use scale_element_base, only: ElementBase1D, ElementBase3D
  use scale_element_line, only: LineElement

  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfield_base, only: MeshField3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase, MeshFieldContainer
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D

  use scale_multigrid_smoother_base, only: &
    MGSmootherBase3D, MGSmoother_PRE_ID
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type, public, extends(MGSmootherBase3D) :: MGSmoother_Poisson3D
    integer :: num_internal_itr_max                !< Maximum number of internal iterations associated with local element solve
    real(RP) :: residual_threshold_l2_internal_itr !< Threshold for each internal iteration associated with local element solve

    integer :: elem_solver_type_id  !< Element solver type ID

    ! For direct solver
    type(LineElement), allocatable :: elem1D_list(:)
  contains
    procedure :: Init => MGSmoother_Poisson3D_Init
    procedure :: Final => MGSmoother_Poisson3D_Final
    procedure :: Advance_itr_1step => MGSmoother_Poisson3D_advance_itr_1step
    procedure, private :: advance_itr_1step_color => Poisson3d_smoother_advance_itr_1step_color
    procedure, private :: cal_q_lc_elem_Direct => cal_q_lc_elem_direct
    procedure, private :: cal_q_lc_elem_PCG => cal_q_lc_elem_PCG
  end type MGSmoother_Poisson3D

  ! public :: poisson_smoother_evaluate_error_norm
  integer, parameter, public :: MGSmoother_Possion3D_AUX_SCALAR_NUM  = 1
  integer, parameter, public :: MGSmoother_Possion3D_AUX_HVEC_NUM    = 1

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: cal_grad_lc

  integer, parameter :: AUXVAR_QZ_ID = 1
  integer, parameter :: AUXVAR_QX_ID = 2
  integer, parameter :: AUXVAR_QY_ID = 3

  integer, parameter :: ELEM_SOLVER_DIRECT_ID = 1
  integer, parameter :: ELEM_SOLVER_PCG_ID    = 2

contains
  !> Initialization
  subroutine MGSmoother_Poisson3D_Init( this, mesh, porder_list )
    use scale_prc, only: PRC_abort
    use scale_multigrid_smoother_base, only: MGSmootherBase3D_Init
    implicit none
    class(MGSmoother_Poisson3D), intent(inout) :: this
    class(MeshBase3D), intent(in), target :: mesh
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
      LOG_INFO("MGSmoother_Poisson3D_Init",*) "Not found namelist PARAM_MGSmoother in the configuration file."
    else if( ierr > 0 ) then
      LOG_INFO("MGSmoother_Poisson3D_Init",*) "Error reading namelist PARAM_MGSmoother in the configuration file."
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
      LOG_INFO("MGSmoother_Poisson3D_Init",*) "Unknown elem_solver type: "//trim(elem_solver)
      call PRC_abort
    end select

    return
  end subroutine MGSmoother_Poisson3D_Init

  !> Finalization
  subroutine MGSmoother_Poisson3D_Final(this)
    use scale_multigrid_smoother_base, only: MGSmootherBase3D_Final
    implicit none
    class(MGSmoother_Poisson3D), intent(inout) :: this

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

    call MGSmootherBase3D_Final( this )
    return
  end subroutine MGSmoother_Poisson3D_Final

  !> Smoother
!OCL SERIAL
  subroutine MGSmoother_Poisson3D_advance_itr_1step( this, &
    q, res,                                      & 
    f, aux_var,                                  &
    itr, var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D, &
    cal_res_flag, zero_initial_guess,            &
    mg_p_level, mg_h_level, pre_or_post_smooth )

    use scale_mesh_hierarchy_base, only: &
      pMG_FINEST_LEVEL => MESH_HIERARCHY_pMG_FINEST_LEVEL
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
    type(SparseMat), intent(in) :: DZ
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
      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)
        !$omp parallel do
        do ke=lmesh3D%NeS, lmesh3D%NeE
          q%local(ldomID)%val(:,ke) = 0.0_RP
        end do
      end do
    end if

    if ( itr == 1 ) then
      call PROF_rapstart( "Poisson3d_smoother_res", 1 )
      call this%advance_itr_1step_color( &
        q, res,                                                                 & 
        f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID), &
        itr, var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,                      &
        .true., 0, is_mg_top_level, mg_p_level )
      
      call this%Set_initial_residual( res )
      call flush(IO_FID_LOG)
      call PROF_rapend( "Poisson3d_smoother_res", 1 )
    end if


    call PROF_rapstart( "Poisson3d_smoother_rb", 1 )
    ! Red
    call this%advance_itr_1step_color( &
      q, res,                                                                 & 
      f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID), &
      itr, var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,                      &
      .false., mod(color_id_offset,2), is_mg_top_level, mg_p_level )
    
    ! Black
    call this%advance_itr_1step_color( &
      q, res,                                                                 & 
      f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID), &
      itr, var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,                      &
      .false., mod(color_id_offset+1,2), is_mg_top_level, mg_p_level )
    call PROF_rapend( "Poisson3d_smoother_rb", 1 )

    if ( cal_res_flag ) then
      call this%advance_itr_1step_color( &
        q, res,                                                                  & 
        f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID),  &
        itr, var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,                       &
        cal_res_flag, 0, is_mg_top_level, mg_p_level )
    end if
    
    var_comm_list(1)%field3d => q
    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    return
  end subroutine MGSmoother_Poisson3D_advance_itr_1step

!-- private -----------------------

!OCL SERIAL
  subroutine Poisson3d_smoother_advance_itr_1step_color( this, &
    q, res,                                                    & 
    f, qx, qy, qz,                                             &
    itr, var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,         &
    cal_res_flag, color_id, is_mg_top_level, mg_p_level )

    use scale_const, only: &
      UNDEF => CONST_UNDEF
    implicit none
    class(MGSmoother_Poisson3D), intent(inout) :: this
    type(MeshField3D), intent(inout), target :: q
    type(MeshField3D), intent(inout) :: res
    type(MeshField3D), intent(in) :: f
    type(MeshField3D), intent(inout), target :: qx
    type(MeshField3D), intent(inout), target :: qy
    type(MeshField3D), intent(inout), target :: qz
    integer, intent(in) :: itr
    class(MeshFieldCommBase), intent(inout) :: var_comm
    class(MeshFieldCommBase), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    class(MeshBase3D), intent(in), target :: mesh3D
    logical, intent(in) :: cal_res_flag
    integer, intent(in) :: color_id
    logical, intent(in) :: is_mg_top_level
    integer, intent(in) :: mg_p_level

    integer :: ldomID
    class(LocalMesh3D), pointer :: lmesh3D

    type(MeshFieldContainer) :: var_comm_list(1)
    type(MeshFieldContainer) :: aux_comm_list(3)
    !---------------------------------------------------------------------------

    !-
    var_comm_list(1)%field3d => q

    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(ldomID)
      call cal_grad_lc( qx%local(ldomID)%val, qy%local(ldomID)%val,  qz%local(ldomID)%val, & ! (out)
        color_id, q%local(ldomID)%val, Dx, Dy, Dz, lmesh3D, lmesh3D%refElem3D ) ! (in)
    end do

    !- 
    aux_comm_list(AUXVAR_QZ_ID)%field3d => qz
    aux_comm_list(AUXVAR_QX_ID)%field3d => qx
    aux_comm_list(AUXVAR_QY_ID)%field3d => qy

    call aux_comm%Put( aux_comm_list, 1 )
    call aux_comm%Exchange()
    call aux_comm%Get( aux_comm_list, 1 )

    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(ldomID)

      select case ( this%elem_solver_type_id )
      case ( ELEM_SOLVER_DIRECT_ID )
        call this%cal_q_lc_elem_Direct( &
          q%local(ldomID)%val, res%local(ldomID)%val,                       &
          color_id, cal_res_flag, is_mg_top_level,                          &
          f%local(ldomID)%val,                                              &
          qx%local(ldomID)%val, qy%local(ldomID)%val, qz%local(ldomID)%val, &
          lmesh3D%VMapM, lmesh3D%VMapP, lmesh3D, lmesh3D%refElem3D,         &
          this%elem1D_list(mg_p_level) )
      case ( ELEM_SOLVER_PCG_ID )
        call this%cal_q_lc_elem_PCG( &
          q%local(ldomID)%val, res%local(ldomID)%val,                       &
          color_id, cal_res_flag, is_mg_top_level,                          &
          f%local(ldomID)%val,                                              &
          qx%local(ldomID)%val, qy%local(ldomID)%val, qz%local(ldomID)%val, &
          Dx, Dy, Dz, Lift,                                                 &
          lmesh3D%VMapM, lmesh3D%VMapP, lmesh3D, lmesh3D%refElem3D )
      end select
    end do
    
    return
  end subroutine Poisson3d_smoother_advance_itr_1step_color

!- private
  subroutine cal_grad_lc( qx, qy, qz, &
    color_id, q, Dx, Dy, Dz, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qy(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qz(elem%Np,lmesh%NeA)
    integer, intent(in) :: color_id
    real(RP), intent(in) :: q(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz

    integer :: ke
    real(RP) :: q1(elem%Np)
    real(RP) :: q2(elem%Np)
    real(RP) :: q3(elem%Np)
    !----------------------------
    
    call PROF_rapstart( "Poisson3d_smoother_cal_grad_lc", 1 )

    !$omp parallel do private(q1, q2, q3)
    do ke=lmesh%NeS, lmesh%NeE
      call Sparsemat_matmul(Dx, q(:,ke), q1)
      call Sparsemat_matmul(Dy, q(:,ke), q2)
      call Sparsemat_matmul(Dz, q(:,ke), q3)

      qx(:,ke) = lmesh%Escale(:,ke,1,1) * q1(:)
      qy(:,ke) = lmesh%Escale(:,ke,2,2) * q2(:)
      qz(:,ke) = lmesh%Escale(:,ke,3,3) * q3(:)
    end do

    call PROF_rapend( "Poisson3d_smoother_cal_grad_lc", 1 )
    return
  end subroutine cal_grad_lc

!OCL SERIAL
  subroutine cal_q_lc_elem_PCG( this, q, res,       &
    color_id, cal_res, is_mg_top_level, rhs, qx, qy,qz, &
    Dx, Dy, Dz, Lift,                                   &
    vmapM, vmapP, lmesh, elem                           )
    use scale_prc
    use scale_linalgebra, only: Linalgebra_SolveLinEq
    implicit none
    class(MGSmoother_Poisson3D), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: q(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: res(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: rhs(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qy(elem%Np,lmesh%NeA)    
    real(RP), intent(in) :: qz(elem%Np,lmesh%NeA)   
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: color_id
    logical, intent(in) :: cal_res
    logical, intent(in) :: is_mg_top_level
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)

    integer :: i, j, k, ke
    integer :: i0, ii
    integer :: ncol

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np)
    real(RP) :: LiftDFlx(elem%Np)
    real(RP) :: BndFlx(elem%NfpTot,4)
    real(RP) :: gtau(elem%NfpTot)

    real(RP) :: Grad_q(elem%Np,3)
    real(RP) :: Grad_q_e(elem%Np,3)
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

    do f=1, elem%Nfaces_h
    do p=1, elem%Nfp_h
      i = p + (f-1)*elem%Nfp_h
      iM(i) = elem%Fmask_h(p,f)
    end do
    end do
    do f=1, elem%Nfaces_v
    do p=1, elem%Nfp_v
      i = p + (f-1)*elem%Nfp_v + elem%Nfp_h * elem%Nfaces_h
      iM(i) = elem%Fmask_v(p,f)
    end do
    end do

  !$omp parallel private(i,j,k,ke,ii,i0,ncol,    &
  !$omp gtau, Fx,Fy,Fz,BndFlx,LiftDFlx,          &
  !$omp Grad_q,Grad_q_e, Ax,x,r,z,q_,            &
  !$omp rnorm0,rnorm, rz_old,rz_new, alpha,beta, &
  !$omp itr, PmatDiag, hk ) 
  !$omp do collapse(2)
    do k=1, lmesh%NeZ
    do j=1, lmesh%NeY

      if ( cal_res ) then
        i0 = 1
        ncol = lmesh%NeX
      else
        i0 = 1 + modulo(color_id - modulo(j + k, 2), 2)  ! 1 or 2
        ncol = (lmesh%NeX - i0)/2 + 1
        if ( ncol < 1 ) cycle
      end if

      do ii=1, ncol
        if (cal_res) then
          i = ii
        else
          i = i0 + (ii-1)*2
        end if
        ke = i + (j-1)*lmesh%NeX + (k-1)*lmesh%NeX*lmesh%NeY
        
        !-
        gtau(:) = lmesh%Fscale(:,ke) * elem%PolyOrder_h * ( elem%PolyOrder_h + 1 ) * 0.5_RP * 1.2_RP

        call eval_Ax( Ax, &
          ke, q, qx, qy, qz, gtau,              &
          Grad_q, BndFlx, Fx, Fy, Fz, LiftDFlx, &
          Dx, Dy, Dz, Lift, iM, lmesh, elem     )

        r(:) = rhs(:,ke) - Ax(:)

        if ( cal_res ) then
          res(:,ke) = r(:)
        else
          rnorm0 = sqrt( sum( r(:) * r(:) ) )

          if ( rnorm0 > 1.0E-16_RP ) then
            ! Precond
            hk = 6.0_RP / maxval(lmesh%Escale(:,ke,1,1) + lmesh%Escale(:,ke,2,2) + lmesh%Escale(:,ke,3,3))
            PmatDiag(:) = c1 * elem%PolyOrder_h**2 / hk + c2 / hk**2

            x(:) = q(:,ke)
            z(:) = r(:) / PmatDiag(:)

            !-
            q_(:) = z(:)
            rz_old = sum( r(:) * z(:) )

            do itr=1, this%num_internal_itr_max
              !-
              call eval_Ax_itr( Ax, &
                ke, q_, gtau, Grad_q, Grad_q_e, &
                BndFlx, Fx, Fy, Fz, LiftDFlx,   &
                Dx, Dy, Dz, Lift, iM, lmesh, elem )

              alpha = rz_old / sum(q_(:) * Ax(:))
              x(:) = x(:) + alpha * q_(:)
              r(:) = r(:) - alpha * Ax(:)

              rnorm = sqrt( sum( r(:) * r(:) ) )
              if ( rnorm / rnorm0 < this%residual_threshold_l2_internal_itr ) then
                exit
              else if ( itr == this%num_internal_itr_max ) then
                LOG_INFO("Poisson3d_smoother_cal_q_lc",*) "ke=", ke, "rnorm/rnorm0=", rnorm / rnorm0, rnorm, rnorm0
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
    ke, q, qx, qy, qz, gtau, &
    Grad_q, BndFlx, Fx, Fy, Fz, LiftDFlx, &
    Dx, Dy, Dz, Lift, iM, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np)
    integer, intent(in) :: ke
    real(RP), intent(in) :: q(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qy(elem%Np,lmesh%NeA)    
    real(RP), intent(in) :: qz(elem%Np,lmesh%NeA)   
    real(RP), intent(in) :: gtau(elem%NfpTot)
    real(RP), intent(inout) :: Grad_q(elem%Np,3) ! [work]
    real(RP), intent(inout) :: BndFlx(elem%NfpTot,4) ! [work]
    real(RP), intent(inout) :: Fx(elem%Np) ! [work]
    real(RP), intent(inout) :: Fy(elem%Np) ! [work]
    real(RP), intent(inout) :: Fz(elem%Np) ! [work]
    real(RP), intent(inout) :: LiftDFlx(elem%Np) ! [work]
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: iM(elem%NfpTot)
    !---------------------------------

    call cal_bndflx_1( BndFlx(:,1:3), &
      q, lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2), lmesh%normal_fn(:,ke,3), &
      lmesh%vmapM(:,ke), lmesh%vmapP(:,ke), lmesh, elem )
    
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,1), LiftDFlx )
    Grad_q(:,1) = qx(:,ke) + LiftDFlx(:)

    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,2), LiftDFlx )
    Grad_q(:,2) = qy(:,ke) + LiftDFlx(:)

    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,3), LiftDFlx )
    Grad_q(:,3) = qz(:,ke) + LiftDFlx(:)
    
    !-
    call cal_bndflx_2( BndFlx(:,4), &
      q, Grad_q(iM,1), Grad_q(iM,2), Grad_q(iM,3), qx, qy, qz,                   &
      lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2), lmesh%normal_fn(:,ke,3), &
      gtau, lmesh%vmapM(:,ke), lmesh%vmapP(:,ke), lmesh, elem )

    call sparsemat_matmul( Dx, Grad_q(:,1), Fx )
    call sparsemat_matmul( Dy, Grad_q(:,2), Fy )
    call sparsemat_matmul( Dz, Grad_q(:,3), Fz )
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,4), LiftDFlx )
    Ax(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDFlx(:)
    return
  end subroutine eval_Ax

!OCL SERIAL
  subroutine cal_bndflx_1( BndFlx, &
    q, nx, ny, nz,        &
    iM, iP, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: BndFlx(elem%NfpTot,3)
    real(RP), intent(in) :: q(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot)
    real(RP), intent(in) :: ny(elem%NfpTot)
    real(RP), intent(in) :: nz(elem%NfpTot)
    integer, intent(in) :: iM(elem%NfpTot)
    integer, intent(in) :: iP(elem%NfpTot)

    real(RP) :: dq(elem%NfpTot)
    !----------------------------------------------------

    dq(:) = q(iP) - q(iM)
    BndFlx(:,1) = 0.5_RP * dq(:) * nx(:)
    BndFlx(:,2) = 0.5_RP * dq(:) * ny(:)
    BndFlx(:,3) = 0.5_RP * dq(:) * nz(:)
    return
  end subroutine cal_bndflx_1

!OCL SERIAL
  subroutine cal_bndflx_2( BndFlx,   &
    q, qx_M, qy_M, qz_M, qx_e, qy_e, qz_e, &
    nx, ny, nz, gtau,                &
    iM, iP, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: BndFlx(elem%NfpTot)
    real(RP), intent(in) :: q(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: qx_M(elem%NfpTot)
    real(RP), intent(in) :: qy_M(elem%NfpTot)    
    real(RP), intent(in) :: qz_M(elem%NfpTot)   
    real(RP), intent(in) :: qx_e(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: qy_e(elem%Np*lmesh%NeA)    
    real(RP), intent(in) :: qz_e(elem%Np*lmesh%NeA)   
    real(RP), intent(in) :: nx(elem%NfpTot)
    real(RP), intent(in) :: ny(elem%NfpTot)
    real(RP), intent(in) :: nz(elem%NfpTot)
    real(RP), intent(in) :: gtau(elem%NfpTot)
    integer, intent(in) :: iM(elem%NfpTot)
    integer, intent(in) :: iP(elem%NfpTot)

    real(RP) :: dq(elem%NfpTot)
    !----------------------------------------------------

    dq(:) = q(iP) - q(iM)
    BndFlx(:) = &
        ( 0.5_RP * ( qx_e(iP) + qx_e(iM) ) - qx_M(:) ) * nx(:) &
      + ( 0.5_RP * ( qy_e(iP) + qy_e(iM) ) - qy_M(:) ) * ny(:) &
      + ( 0.5_RP * ( qz_e(iP) + qz_e(iM) ) - qz_M(:) ) * nz(:) &
      + 0.5_RP * gtau(:) * dq(:)
    return
  end subroutine cal_bndflx_2  

!OCL SERIAL
  subroutine eval_Ax_itr( Ax, &
    ke, q, gtau, &
    Grad_q, Grad_q_e, BndFlx, Fx, Fy, Fz, LiftDFlx, &
    Dx, Dy, Dz, Lift, iM, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np)
    integer, intent(in) :: ke
    real(RP), intent(in) :: q(elem%Np)
    real(RP), intent(in) :: gtau(elem%NfpTot)
    real(RP), intent(inout) :: Grad_q(elem%Np,3) ! [work]
    real(RP), intent(inout) :: Grad_q_e(elem%Np,3) ! [work]
    real(RP), intent(inout) :: BndFlx(elem%NfpTot,4) ! [work]
    real(RP), intent(inout) :: Fx(elem%Np) ! [work]
    real(RP), intent(inout) :: Fy(elem%Np) ! [work]
    real(RP), intent(inout) :: Fz(elem%Np) ! [work]
    real(RP), intent(inout) :: LiftDFlx(elem%Np) ! [work]
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    integer, intent(in) :: iM(elem%NfpTot)
    !---------------------------------

    call cal_bndflx_itr_1( BndFlx(:,1:3), &
      q(iM),                                                                     &
      lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2), lmesh%normal_fn(:,ke,3), &
      lmesh, elem )

    call sparsemat_matmul( Dx, q(:), Fx )
    call sparsemat_matmul( Dy, q(:), Fy )
    call sparsemat_matmul( Dz, q(:), Fz )

    Grad_q_e(:,1) = lmesh%Escale(:,ke,1,1) * Fx(:) 
    Grad_q_e(:,2) = lmesh%Escale(:,ke,2,2) * Fy(:) 
    Grad_q_e(:,3) = lmesh%Escale(:,ke,3,3) * Fz(:) 

    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,1), LiftDFlx )
    Grad_q(:,1) = Grad_q_e(:,1) + LiftDFlx(:)

    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,2), LiftDFlx )
    Grad_q(:,2) = Grad_q_e(:,2) + LiftDFlx(:)

    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,3), LiftDFlx )
    Grad_q(:,3) = Grad_q_e(:,3) !+ LiftDFlx(:)

    !-
    call cal_bndflx_itr_2( BndFlx(:,4), &
      q(iM), Grad_q(iM,1), Grad_q(iM,2), Grad_q(iM,3), Grad_q_e(iM,1), Grad_q_e(iM,2), Grad_q_e(iM,3), &
      lmesh%normal_fn(:,ke,1), lmesh%normal_fn(:,ke,2), lmesh%normal_fn(:,ke,3), gtau, lmesh, elem     )

    call sparsemat_matmul( Dx, Grad_q(:,1), Fx )
    call sparsemat_matmul( Dy, Grad_q(:,2), Fy )
    call sparsemat_matmul( Dz, Grad_q(:,3), Fz )
    call sparsemat_matmul( Lift, lmesh%Fscale(:,ke)*BndFlx(:,4), LiftDFlx )
    Ax(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + 0.0_RP * lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDFlx(:)
    return
  end subroutine eval_Ax_itr

!OCL SERIAL
  subroutine cal_bndflx_itr_1( BndFlx, &
    q_M, nx, ny, nz,      &
    lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: BndFlx(elem%NfpTot,3)
    real(RP), intent(in) :: q_M(elem%NfpTot)
    real(RP), intent(in) :: nx(elem%NfpTot)
    real(RP), intent(in) :: ny(elem%NfpTot)
    real(RP), intent(in) :: nz(elem%NfpTot)

    real(RP) :: dq(elem%NfpTot)
    !----------------------------------------------------

    dq(:) =      - q_M(:)
    BndFlx(:,1) = 0.5_RP * dq(:) * nx(:)
    BndFlx(:,2) = 0.5_RP * dq(:) * ny(:)
    BndFlx(:,3) = 0.5_RP * dq(:) * nz(:)
    return
  end subroutine cal_bndflx_itr_1

!OCL SERIAL
  subroutine cal_bndflx_itr_2( BndFlx, &
    q, qx, qy, qz, qx_e, qy_e, qz_e,   &
    nx, ny, nz, gtau,                  &
    lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: BndFlx(elem%NfpTot)
    real(RP), intent(in) :: q(elem%NfpTot)
    real(RP), intent(in) :: qx(elem%NfpTot)
    real(RP), intent(in) :: qy(elem%NfpTot)    
    real(RP), intent(in) :: qz(elem%NfpTot)   
    real(RP), intent(in) :: qx_e(elem%NfpTot)
    real(RP), intent(in) :: qy_e(elem%NfpTot)    
    real(RP), intent(in) :: qz_e(elem%NfpTot)   
    real(RP), intent(in) :: nx(elem%NfpTot)
    real(RP), intent(in) :: ny(elem%NfpTot)
    real(RP), intent(in) :: nz(elem%NfpTot)
    real(RP), intent(in) :: gtau(elem%NfpTot)

    real(RP) :: dq(elem%NfpTot)
    !----------------------------------------------------

    dq(:) =      - q(:)
    BndFlx(:) = &
        ( 0.5_RP * (          + qx_e(:) ) - qx(:) ) * nx(:) &
      + ( 0.5_RP * (          + qy_e(:) ) - qy(:) ) * ny(:) &
      + ( 0.5_RP * (          + qz_e(:) ) - qz(:) ) * nz(:) &
      + 0.5_RP * gtau(:) * dq(:) * (1.0_RP-nz(:)**2)
    return
  end subroutine cal_bndflx_itr_2  

!- Direct solver within an elemnt -----------------------------

!OCL SERIAL
  subroutine cal_q_lc_elem_direct( this, q, res,        &
    color_id, cal_res, is_mg_top_level, rhs, qx, qy,qz, &
    vmapM, vmapP, lmesh, elem, elem1D                   )
    use scale_prc
    use scale_linalgebra, only: Linalgebra_SolveLinEq
    implicit none
    class(MGSmoother_Poisson3D), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(ElementBase1D), intent(in) :: elem1D
    real(RP), intent(inout) :: q(elem%Np*lmesh%NeA)
    real(RP), intent(out) :: res(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: rhs(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qx(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: qy(elem%Np*lmesh%NeA)    
    real(RP), intent(in) :: qz(elem%Np*lmesh%NeA)    
    integer, intent(in) :: color_id
    logical, intent(in) :: cal_res
    logical, intent(in) :: is_mg_top_level
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)

    integer :: ke, f1
    integer :: i, j, k

    integer :: p

    real(RP) :: E11, E22, E33

    real(RP) :: Dx(elem%Np,elem%Np)
    real(RP) :: Dy(elem%Np,elem%Np)
    real(RP) :: Dz(elem%Np,elem%Np)
    real(RP) :: Sx_tr0(elem%Np,elem%Np)
    real(RP) :: Sy_tr0(elem%Np,elem%Np)
    real(RP) :: Sz_tr0(elem%Np,elem%Np)
    real(RP) :: Sx_tr(elem%Np,elem%Np)
    real(RP) :: Sy_tr(elem%Np,elem%Np)
    real(RP) :: Sz_tr(elem%Np,elem%Np)

    real(RP) :: GsqrtDx(elem%Np,elem%Np)
    real(RP) :: GsqrtDy(elem%Np,elem%Np)
    real(RP) :: GsqrtDz(elem%Np,elem%Np)
    real(RP) :: GsqrtI(elem%Np,elem%Np)
  
    real(RP) :: GsqrtDn1(elem%Np,elem%Np)

    real(RP) :: OPT11(elem%Np,elem%Np)
    real(RP) :: mmE(elem%Np,elem%Np)   
    real(RP) :: Lift_(elem%Np,elem%Np)
    real(RP) :: Lift_x(elem%Np,elem%Np)       
    real(RP) :: Lift_y(elem%Np,elem%Np)
    real(RP) :: Lift_h_tmp(elem%Np,elem%Nfp_h)       
    real(RP) :: Lift_z(elem%Np,elem%Np)
    real(RP) :: Lift_x_q(elem%Np)       
    real(RP) :: Lift_y_q(elem%Np)
    real(RP) :: Lift_z_q(elem%Np)
    real(RP) :: Msrc(elem%Np)

    real(RP) :: Emat(elem%Np,elem%NfpTot)
    real(RP) :: massEdge_h(elem%Np,elem%Np,elem%Nfaces_h)
    real(RP) :: massEdge_v(elem%Np,elem%Np,elem%Nfaces_v)
    integer :: Fm1_h(elem%Nfp_h)
    integer :: Fm1_v(elem%Nfp_v)
    integer :: vid_s, vid_e
    real(RP) :: lnx, lny, lnz

    real(RP) :: q_(elem%Np)
    real(RP) :: q_P_h(elem%Nfp_h)
    real(RP) :: q_P_v(elem%Nfp_v)
    real(RP) :: gradQ_P_h(elem%Nfp_h,2)
    real(RP) :: gradQ_P_v(elem%Nfp_v)
    real(RP) :: GsqrtDn2_h(elem%Nfp_h)
    real(RP) :: GsqrtDn2_v(elem%Nfp_v)

    real(RP) :: tmp, tmp1, tmp2
    integer :: m, n

    real(RP) :: gtau_

    integer :: i0, i_int
    integer :: ii, ncol

    integer :: nfp_os

    real(RP) :: q_old
    real(RP), parameter :: omg = 1.0_RP

    real(RP) :: MassMat(elem%Np,elem%Np)
    real(RP) :: S1Dtr_D1D(elem1D%Np,elem1D%Np)
    !---------------------------

    call PROF_rapstart( "Poisson3d_smoother_cal_q_lc", 1 )

    call PROF_rapstart( "Poisson3d_smoother_cal_q_lc_0", 1 )

    Emat(:,:) = matmul(elem%M, elem%Lift)

    massEdge_h(:,:,:)  = 0.0_RP
    do f1=1, elem%Nfaces_h
      Fm1_h(:) = elem%Fmask_h(:,f1)
      massEdge_h(Fm1_h,Fm1_h,f1) = Emat(Fm1_h, (f1-1)*elem%Nfp_h+1:f1*elem%Nfp_h)
    end do

    massEdge_v(:,:,:)  = 0.0_RP
    nfp_os = elem%Nfaces_h * elem%Nfp_h
    do f1=1, elem%Nfaces_v
      Fm1_v(:) = elem%Fmask_v(:,f1)
      massEdge_v(Fm1_v,Fm1_v,f1) = Emat(Fm1_v, nfp_os+(f1-1)*elem%Nfp_v+1:nfp_os+f1*elem%Nfp_v)
    end do

    Sx_tr0(:,:) = matmul(transpose(elem%Dx1),elem%M)
    Sy_tr0(:,:) = matmul(transpose(elem%Dx2),elem%M)
    Sz_tr0(:,:) = matmul(transpose(elem%Dx3),elem%M)

    MassMat(:,:) = elem%M(:,:)

    S1Dtr_D1D(:,:) = matmul(transpose(elem1D%Dx1), elem1D%M)
    S1Dtr_D1D(:,:) = matmul(S1Dtr_D1D,elem1D%Dx1)
    call PROF_rapend( "Poisson3d_smoother_cal_q_lc_0", 1 )

    !$omp parallel private( i0, i_int, ncol, i, ii, j, k, ke, &
    !$omp E11, E22, E33, Dx, Dy, Dz, Sx_tr, Sy_tr, Sz_tr, GsqrtI, GsqrtDx, GsqrtDy, GsqrtDz,   &
    !$omp GsqrtDn1, GsqrtDn2_h, GsqrtDn2_v,                                          &
    !$omp OPT11, mmE, Lift_, Lift_x, Lift_y, Lift_h_tmp, Lift_z, Lift_x_q, Lift_y_q,  Lift_z_q,       &
    !$omp Msrc, f1, vid_s, vid_e, lnx, lny, lnz, Fm1_h, Fm1_v,                   &
    !$omp q_, q_P_h, q_P_v, gradQ_P_h, gradQ_P_v, gtau_, p, q_old, & 
    !$omp tmp, tmp1, tmp2     )
    !$omp do collapse(2)
    do k=1, lmesh%NeZ
    do j=1, lmesh%NeY

      if ( cal_res ) then
        i0 = 1
        ncol = lmesh%NeX
      else
        i0 = 1 + modulo(color_id - modulo(j + k, 2), 2)  ! 1 or 2
        ncol = (lmesh%NeX - i0)/2 + 1
        if ( ncol < 1 ) cycle
      end if

      do ii=1, ncol
        if (cal_res) then
          i = ii
        else
          i = i0 + (ii-1)*2
        end if
        ke = i + (j-1)*lmesh%NeX + (k-1)*lmesh%NeX*lmesh%NeY
        
        ! Assume E11, E22, and E33 are constant in an element
        E11 = lmesh%Escale(1,ke,1,1)
        E22 = lmesh%Escale(1,ke,2,2)
        E33 = lmesh%Escale(1,ke,3,3)

        Dx(:,:) = E11 * elem%Dx1
        Dy(:,:) = E22 * elem%Dx2
        Dz(:,:) = E33 * elem%Dx3
        Sx_tr(:,:) = E11 * Sx_tr0(:,:)
        Sy_tr(:,:) = E22 * Sy_tr0(:,:)
        Sz_tr(:,:) = E33 * Sz_tr0(:,:)

        do p=1, elem%Np
          GsqrtI(:,p) = lmesh%Gsqrt(:,ke)
          GsqrtDx(:,p) = GsqrtI(:,p) * Dx(:,p)
          GsqrtDy(:,p) = GsqrtI(:,p) * Dy(:,p)
          GsqrtDz(:,p) = GsqrtI(:,p) * Dz(:,p)
        end do

        Msrc(:) = matmul(MassMat, rhs(:,ke) * lmesh%Gsqrt(:,ke) )

        call construct_OPT11_elem( OPT11, &
          S1Dtr_D1D, elem1D%M, E11, E22, lmesh%Gsqrt(1,ke), elem1D%Np )

        do f1=1, elem%Nfaces_h       
          vid_s = 1 + (f1-1)*elem%Nfp_h
          vid_e = vid_s + elem%Nfp_h - 1

          lnx = lmesh%normal_fn(vid_s,ke,1)
          lny = lmesh%normal_fn(vid_s,ke,2)
          Fm1_h(:) = elem%Fmask_h(:,f1)
                      
          mmE(:,:) = lmesh%Fscale(vid_s,ke) * massEdge_h(:,:,f1)  
                          
          Lift_(:,:) = matmul(elem%invM, mmE)
          do p=1, elem%Np
            Lift_x(:,p) = lmesh%Gsqrt(:,ke) * lnx * Lift_(:,p)
            Lift_y(:,p) = lmesh%Gsqrt(:,ke) * lny * Lift_(:,p)
          end do
          
          gtau_ = elem%PolyOrder_h * ( elem%PolyOrder_h + 1 ) * lmesh%Fscale(vid_s,ke) * 0.5_RP * 1.2_RP
          
          GsqrtDn1(:,:) = lnx * GsqrtDx(:,:) + lny * GsqrtDy(:,:)

          OPT11(:,:) = OPT11(:,:) + 0.5_RP * ( &
            - gtau_ * mmE(:,:) + matmul(mmE, GsqrtDn1)         &
            + ( matmul(Sx_tr,Lift_x) + matmul(Sy_tr,Lift_y) ) )
          
          q_P_h(:) = q(VMapP(vid_s:vid_e,ke))
          
          Lift_h_tmp(:,:) = Lift_x(:,Fm1_h)
          Lift_x_q(:) = matmul(Lift_h_tmp, q_P_h(:))

          Lift_h_tmp(:,:) = Lift_y(:,Fm1_h)
          Lift_y_q(:) = matmul(Lift_h_tmp, q_P_h(:))

          Msrc(:) = Msrc(:) + 0.5_RP * ( &
            + matmul(Sx_tr,Lift_x_q) + matmul(Sy_tr,Lift_y_q) )

          gradQ_P_h(:,1) = qx(VMapP(vid_s:vid_e,ke))
          gradQ_P_h(:,2) = qy(VMapP(vid_s:vid_e,ke))
          GsqrtDn2_h(:) = lmesh%Gsqrt(Fm1_h(:),ke) * ( lnx * gradQ_P_h(:,1) + lny * gradQ_P_h(:,2) )
          Msrc(Fm1_h) = Msrc(Fm1_h) - 0.5_RP * (  &
            matmul(mmE(Fm1_h,Fm1_h), gtau_ * q_P_h(:) + GsqrtDn2_h(:)) )
        end do 

        do f1=1, -1!elem%Nfaces_v       
          vid_s = 1 + elem%Nfaces_h*elem%Nfp_h + (f1-1)*elem%Nfp_v
          vid_e = vid_s + elem%Nfp_v - 1

          lnz = lmesh%normal_fn(vid_s,ke,3)
          Fm1_v(:) = elem%Fmask_v(:,f1)
          
          GsqrtDn1(:,:) = lnz * GsqrtDz(:,:)
          
          mmE(:,:) = lmesh%Fscale(vid_s,ke) * massEdge_v(:,:,f1)  
                          
          Lift_(:,:) = matmul(elem%invM, mmE)
          do p=1, elem%Np
            Lift_z(:,p) = lmesh%Gsqrt(:,ke) * Lift_(:,p)
          end do
          
          gtau_ = elem%PolyOrder_v * (elem%PolyOrder_v + 1 ) * lmesh%Fscale(vid_s,ke) * 0.5_RP * 1.2_RP
          
          OPT11(:,:) = OPT11(:,:) + 0.5_RP * ( &
            - gtau_ * mmE(:,:) + matmul(mmE, GsqrtDn1) &
            + matmul(Sz_tr,Lift_z)                     )
          
          q_P_v(:) = q(VMapP(vid_s:vid_e,ke))
          Lift_z_q(:) = matmul(Lift_z(:,Fm1_v), q_P_v(:))
          Msrc(:) = Msrc(:) + 0.5_RP * ( &
            + matmul(Sz_tr,Lift_z_q)     )

          gradQ_P_v(:) = qz(VMapP(vid_s:vid_e,ke))
          GsqrtDn2_v(:) = lmesh%Gsqrt(Fm1_v(:),ke) * ( lnz * gradQ_P_v(:) )
          Msrc(Fm1_v) = Msrc(Fm1_v) - 0.5_RP * (  &
            matmul(mmE(Fm1_v,Fm1_v), gtau_ * q_P_v(:) + GsqrtDn2_v(:)) )
        end do 
        
        if ( cal_res ) then
          do p=1, elem%Np
            q_(p) = q(p+(ke-1)*elem%Np)
          end do
          res(:,ke) = Msrc(:) - matmul(OPT11(:,:), q_(:))
          res(:,ke) = matmul(elem%invM, res(:,ke))
        else
          ! call PROF_rapstart( "Poisson3d_smoother_cal_q_lc_linsol", 1 )
          call LinAlgebra_SolveLinEq( OPT11, Msrc, q_ )
          do p=1, elem%Np
            q_old = q(p+(ke-1)*elem%Np)
            q(p+(ke-1)*elem%Np) = q_old + omg * ( q_(p) - q_old )
          end do
          ! call PROF_rapend( "Poisson3d_smoother_cal_q_lc_linsol", 1 )
        end if        
      end do
    end do
    end do
    !$omp end parallel

    call PROF_rapend( "Poisson3d_smoother_cal_q_lc", 1 )
    return
  end subroutine cal_q_lc_elem_direct

!OCL SERIAL
  subroutine construct_OPT11_elem( OPT11, &
    S1Dtr_D1D, M1D,       &
    E11, E22, Gsqrt, Np1D )
    implicit none
    integer, intent(in) :: Np1D
    real(RP), intent(out) :: OPT11(Np1D,Np1D,Np1D,Np1D,Np1D,Np1D)
    real(RP), intent(in) :: S1Dtr_D1D(Np1D,Np1D)
    real(RP), intent(in) :: M1D(Np1D,Np1D)
    real(RP), intent(in) :: E11
    real(RP), intent(in) :: E22
    real(RP), intent(in) :: Gsqrt

    integer :: i, j, k
    integer :: l, m, n

    real(RP) :: coef
    real(RP) :: tmp
    !-----------------------------------------------

    coef = - E11**2 * Gsqrt
    do n=1, Np1D; do m=1, Np1D; do l=1, Np1D
      do k=1, Np1D; do j=1, Np1D
        tmp = coef * M1D(j,m) * M1D(k,n)
        do i=1, Np1D
          OPT11(i,j,k,l,m,n) = S1Dtr_D1D(i,l) * tmp
        end do
      end do; end do
    end do; end do; end do

    coef = - E22**2 * Gsqrt
    do n=1, Np1D; do m=1, Np1D; do l=1, Np1D
      do k=1, Np1D; do j=1, Np1D
        tmp = coef * S1Dtr_D1D(j,m) * M1D(k,n)
        do i=1, Np1D
          OPT11(i,j,k,l,m,n) = OPT11(i,j,k,l,m,n) + M1D(i,l) * tmp
        end do
      end do; end do
    end do; end do; end do
    return
  end subroutine construct_OPT11_elem
  
end module mod_poisson3d_smoother

