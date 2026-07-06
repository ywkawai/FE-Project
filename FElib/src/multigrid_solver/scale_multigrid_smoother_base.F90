!-------------------------------------------------------------------------------
!> module FElib / Multigrid / Smoother base
!!
!! @par Description
!!      A module to provide derived type to manage multigrid smoother
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_multigrid_smoother_base
   !-----------------------------------------------------------------------------
   !
   !++ used modules
   !
   !
   use scale_precision
   use scale_io
   use scale_prc, only: PRC_abort

   use scale_sparsemat, only: SparseMat

   use scale_element_quadrilateral, only: QuadrilateralElement
   use scale_element_hexahedral, only: HexahedralElement

   use scale_mesh_base2d, only: MeshBase2D
   use scale_localmesh_2d, only: LocalMesh2D
   use scale_mesh_base3d, only: MeshBase3D
   use scale_localmesh_3d, only: LocalMesh3D
   use scale_mesh_rectdom2d, only: MeshRectDom2D
   use scale_mesh_cubedom3d, only: MeshCubeDom3D

   use scale_meshfield_base, only: &
     MeshField2D, MeshField3D
   use scale_meshfieldcomm_base, only: MeshFieldCommBase
   use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D
   use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
   !-----------------------------------------------------------------------------
   implicit none
   private

   !-----------------------------------------------------------------------------
   !
   !++ Public type & procedure
   ! 
   integer, public, parameter :: MGSmoother_PRE_ID  = 1  !< ID to represent pre-smoothing
   integer, public, parameter :: MGSmoother_POST_ID = 2  !< ID to represent post-smoothing

   !> Base type for multigrid smoother
   type, public :: MGSmootherBase
      integer :: itration           !< Current smoothing iteration
      integer :: num_smooth_itr_max !< Maximum number of smoothing iterations

      real(RP) :: threshold_ratio_residual_l2  !< Threshold of decreasing ratio of residual for smoothing (L2 norm)
      real(RP) :: threshold_residual_l2  !< Residual threshold for smoothing (L2 norm)
      real(RP) :: threshold_residual_max !< Residual threshold for smoothing (max norm)

      integer :: history_residual_step                 !< Step recording residual history
      integer :: history_residual_rstep                !< Remaining step to record residual history
      integer :: history_residual_count                !< Count of recorded residuals
      real(RP), allocatable :: history_residual_l2(:)  !< Residual history of L2 norm
      real(RP) :: history_residual_l2_initial          !< Initial residual L2 norm
      real(RP), allocatable :: history_residual_max(:) !< Residual history of max norm
      real(RP) :: history_residual_max_initial         !< Initial residual max norm
   contains
      procedure :: Set_smoothing_param => MGSmootherBase_set_smoothing_param
      procedure :: Output_residual_history => MGSmootherBase_output_residual_history
      procedure :: Is_converged => MGSmootherBase_is_converged
      procedure :: Get_initial_residual_statistics => MGSmootherBase_get_initial_residual_statistics
      procedure :: Get_current_residual_statistics => MGSmootherBase_get_current_residual_statistics
      procedure, private :: Is_residual_check_step => MGSmootherBase_is_residual_check_step
      procedure, private :: Eval_global_statistic => MGSmootherBase_eval_global_statistic
      
   end type MGSmootherBase
   public :: MGSmootherBase_Init
   public :: MGSmootherBase_Final
   public :: MGSmootherBase_put_residual_history
   public :: MGSmootherBase_get_statistic_lc

   integer, public, parameter :: MGSmoother_DEFAULT_NUM_SMOOTH_ITE_MAX = 10
   real(RP), public, parameter :: MGSmoother_DEFAULT_THRESHOLD_RATIO_RESIDUAL_L2 = 1.0E-3_RP
   real(RP), public, parameter :: MGSmoother_DEFAULT_THRESHOLD_RESIDUAL_L2  = 1.0E-3_RP
   real(RP), public, parameter :: MGSmoother_DEFAULT_THRESHOLD_RESIDUAL_MAX = 1.0E-3_RP


   !> Derived type for 2D multigrid smoother
   type, abstract, extends(MGSmootherBase), public :: MGSmootherBase2D
   contains
     procedure :: Do_smoothing => MGSmootherBase2D_Do_smoothing
     procedure(MGSmootherBase2D_Advance_itr_1step), deferred :: Advance_itr_1step
     procedure :: Eval_residual_statistic => MGSmootherBase2D_eval_residual_statistic
     procedure :: Set_initial_residual => MGSmootherBase2D_set_initial_residual
     procedure, private :: Put_residual => MGSmootherBase2D_put_residual
   end type MGSmootherBase2D
   public :: MGSmootherBase2D_Init
   public :: MGSmootherBase2D_Final

   abstract interface
    subroutine MGSmootherBase2D_Advance_itr_1step( this, q, res, &
      f, aux_var, itr, var_comm, aux_comm, Dx, Dy, Lift, mesh2D, &
      cal_res_flag, zero_initial_guess, mg_p_level, mg_h_level,  &
      pre_or_post_smooth )
      import MeshBase2D
      import :: MeshField2D
      import MeshFieldCommBase
      import SparseMat
      import :: MGSmootherBase2D
      class(MGSmootherBase2D), intent(inout) :: this
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
    end subroutine MGSmootherBase2D_Advance_itr_1step
   end interface  

   !> Derived type for 3D multigrid smoother
   type, abstract, extends(MGSmootherBase), public :: MGSmootherBase3D
   contains
    procedure, public :: Do_smoothing => MGSmootherBase3D_Do_smoothing
    procedure(MGSmootherBase3D_Advance_itr_1step), deferred :: Advance_itr_1step
    procedure :: Eval_residual_statistic => MGSmootherBase3D_eval_residual_statistic
    procedure :: Set_initial_residual => MGSmootherBase3D_set_initial_residual
    procedure, private :: Put_residual => MGSmootherBase3D_put_residual
   end type MGSmootherBase3D
   public :: MGSmootherBase3D_Init
   public :: MGSmootherBase3D_Final

  abstract interface
    subroutine MGSmootherBase3D_Advance_itr_1step( this, q, res,      &
      f, aux_var, itr, var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,  &
      cal_res_flag, zero_initial_guess, mg_p_level, mg_h_level,       &
      pre_or_post_smooth )
      import MeshBase3D
      import :: MeshField3D
      import MeshFieldCommBase
      import SparseMat
      import :: MGSmootherBase3D
      class(MGSmootherBase3D), intent(inout) :: this
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
    end subroutine MGSmootherBase3D_Advance_itr_1step
  end interface  

contains
!-- Base
!OCL SERIAL
  subroutine MGSmootherBase_Init(this)
    implicit none
    class(MGSmootherBase), intent(inout) :: this
    !---------------------------------
    this%num_smooth_itr_max  = MGSmoother_DEFAULT_NUM_SMOOTH_ITE_MAX
    this%threshold_ratio_residual_l2 = MGSmoother_DEFAULT_THRESHOLD_RATIO_RESIDUAL_L2
    this%threshold_residual_l2  = MGSmoother_DEFAULT_THRESHOLD_RESIDUAL_L2
    this%threshold_residual_max = MGSmoother_DEFAULT_THRESHOLD_RESIDUAL_MAX
    return
  end subroutine MGSmootherBase_Init

!OCL SERIAL
  subroutine MGSmootherBase_Final(this)
    implicit none
    class(MGSmootherBase), intent(inout) :: this
    !---------------------------------
    return
  end subroutine MGSmootherBase_Final

  !> Set smoothing parameters
!OCL SERIAL
  subroutine MGSmootherBase_set_smoothing_param( this, &
    num_smooth_itr_max,       &
    threshold_ratio_residual_l2, &
    threshold_residual_l2,    &
    threshold_residual_max,   &
    history_residual_step     )
    implicit none
    class(MGSmootherBase), intent(inout) :: this
    integer, intent(in) :: num_smooth_itr_max      !< Maximum number of smoothing iterations
    real(RP), intent(in) :: threshold_ratio_residual_l2  !< Threshold of decreasing ratio of residual for smoothing (L2 norm)
    real(RP), intent(in) :: threshold_residual_l2  !< Residual threshold for smoothing (L2 norm)
    real(RP), intent(in) :: threshold_residual_max !< Residual threshold for smoothing (max norm)
    integer, intent(in) :: history_residual_step   !< Step recording residual history

    integer :: num_history_residual
    !---------------------------------

    this%num_smooth_itr_max = num_smooth_itr_max
    this%threshold_ratio_residual_l2 = threshold_ratio_residual_l2
    this%threshold_residual_l2  = threshold_residual_l2
    this%threshold_residual_max = threshold_residual_max 

    !-
    this%history_residual_step = history_residual_step
    num_history_residual = this%num_smooth_itr_max / this%history_residual_step + 1

    allocate( this%history_residual_l2(num_history_residual) )
    allocate( this%history_residual_max(num_history_residual) )

    LOG_INFO( 'MGSmootherBase_set_smoothing_param',* ) "  Smoothing parameters set:"
    LOG_INFO( 'MGSmootherBase_set_smoothing_param',* ) "    num_smooth_itr_max          =", this%num_smooth_itr_max
    LOG_INFO( 'MGSmootherBase_set_smoothing_param',* ) "    threshold_ratio_residual_l2 =", this%threshold_ratio_residual_l2
    LOG_INFO( 'MGSmootherBase_set_smoothing_param',* ) "    threshold_residual_l2       =", this%threshold_residual_l2
    LOG_INFO( 'MGSmootherBase_set_smoothing_param',* ) "    threshold_residual_max      =", this%threshold_residual_max
    LOG_INFO( 'MGSmootherBase_set_smoothing_param',* ) "    history_residual_step       =", this%history_residual_step
    LOG_INFO( 'MGSmootherBase_set_smoothing_param',* ) "    num of history residual     =", num_history_residual
    LOG_INFO( 'MGSmootherBase_set_smoothing_param',* ) "----------------------------------------"
    return
  end subroutine MGSmootherBase_set_smoothing_param

!OCL SERIAL
  function MGSmootherBase_is_residual_check_step( this ) result( is_check_step )
    implicit none
    class(MGSmootherBase), intent(in) :: this
    logical :: is_check_step
    !---------------------------------

    if ( this%history_residual_rstep == 1 ) then
      is_check_step = .true.
    else
      is_check_step = .false.
    end if
    return
  end function MGSmootherBase_is_residual_check_step

!OCL SERIAL
  subroutine MGSmootherBase_get_initial_residual_statistics( this, res_l2_init, res_max_init )
    implicit none
    class(MGSmootherBase), intent(in) :: this
    real(RP), intent(out) :: res_l2_init
    real(RP), intent(out) :: res_max_init
    !---------------------------------

    res_l2_init = this%history_residual_l2_initial
    res_max_init = this%history_residual_max_initial
    return
  end subroutine MGSmootherBase_get_initial_residual_statistics

!OCL SERIAL
  subroutine MGSmootherBase_get_current_residual_statistics( this, res_l2, res_max )
    implicit none
    class(MGSmootherBase), intent(in) :: this
    real(RP), intent(out) :: res_l2
    real(RP), intent(out) :: res_max
    !---------------------------------

    res_l2 = this%history_residual_l2( this%history_residual_count )
    res_max = this%history_residual_max( this%history_residual_count )
    return
  end subroutine MGSmootherBase_get_current_residual_statistics

!OCL SERIAL
  function MGSmootherBase_is_converged( this, is_Vcycle_end ) result( conv_flag )
    implicit none
    class(MGSmootherBase), intent(in) :: this
    logical, intent(in) :: is_Vcycle_end
    logical :: conv_flag

    real(RP) :: res_l2, res_max
    !---------------------------------

    call this%Get_current_residual_statistics( res_l2, res_max )

    conv_flag = .false.
    if ( is_Vcycle_end ) then
      if ( (      res_l2 < this%threshold_residual_l2         &
                .and. res_max < this%threshold_residual_max ) ) then
        conv_flag = .true.
      end if
    else
      if ( ( res_l2 / this%history_residual_l2_initial < this%threshold_ratio_residual_l2 ) &
          .or. (      res_l2 < this%threshold_residual_l2                                   &
                .and. res_max < this%threshold_residual_max ) ) then
        conv_flag = .true.
      end if
    end if
    return
  end function MGSmootherBase_is_converged

!OCL SERIAL
  subroutine MGSmootherBase_output_residual_history( this )
    implicit none
    class(MGSmootherBase), intent(in) :: this
    
    integer :: i
    !---------------------------------

    LOG_INFO( 'MGSmootherBase_output_residual_history','(A,2ES12.4)') "Smoothing residual history (ini, L2 norm, Max norm) : ", this%history_residual_l2_initial, this%history_residual_max_initial
    LOG_INFO( 'MGSmootherBase_output_residual_history',*) "Smoothing residual history (itr, L2 norm, Max norm) : "
    do i=1, this%history_residual_count
      LOG_INFO( 'MGSmootherBase_output_residual_history', '(I8, 2ES12.4)' ) i*this%history_residual_step, &
        this%history_residual_l2(i), this%history_residual_max(i) 
    end do
    return
  end subroutine MGSmootherBase_output_residual_history

!-- 2D -----------------------------------------

!OCL SERIAL
  subroutine MGSmootherBase2D_Init(this)
    implicit none
    class(MGSmootherBase2D), intent(inout) :: this
    !---------------------------------
    call MGSmootherBase_Init( this )
    return
  end subroutine MGSmootherBase2D_Init

!OCL SERIAL
  subroutine MGSmootherBase2D_Final(this)
    implicit none
    class(MGSmootherBase2D), intent(inout) :: this
    !---------------------------------
    call MGSmootherBase_Final( this )
    return
  end subroutine MGSmootherBase2D_Final

!OCL SERIAL
  subroutine MGSmootherBase2D_Do_smoothing( this, q, res,  &
    f, aux_var, var_comm, aux_comm, Dx, Dy, Lift, mesh2D,  &
    mg_p_level, mg_h_level,                                &
    pre_or_post_smooth )
    implicit none
    class(MGSmootherBase2D), intent(inout) :: this
    type(MeshField2D), intent(inout), target :: q
    type(MeshField2D), intent(inout) :: res
    type(MeshField2D), intent(in) :: f
    type(MeshField2D), intent(inout), target :: aux_var(:)
    class(MeshFieldCommBase), intent(inout) :: var_comm
    class(MeshFieldCommBase), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift
    class(MeshBase2D), intent(in), target :: mesh2D
    integer, intent(in) :: mg_p_level
    integer, intent(in) :: mg_h_level
    integer, intent(in) :: pre_or_post_smooth

    integer :: itr
    logical :: zero_initial_guess
    logical :: cal_res_flag
    logical :: conv_flag
    !--------------------------------------

    this%history_residual_count = 0
    this%history_residual_rstep = this%history_residual_step
    
    do itr=1, this%num_smooth_itr_max
      zero_initial_guess = ( pre_or_post_smooth == MGSmoother_PRE_ID .and. mg_p_level /= 1 .and. itr == 1 )
      cal_res_flag = this%Is_residual_check_step()

      call this%Advance_itr_1step( q, res,                          &
        f, aux_var, itr, var_comm, aux_comm, Dx, Dy, Lift, mesh2D,  &
        cal_res_flag, zero_initial_guess, mg_p_level, mg_h_level,   &
        pre_or_post_smooth )

      call this%put_residual( res, conv_flag )
      if ( cal_res_flag ) then
        if ( conv_flag ) then
          exit
        else if ( itr == this%num_smooth_itr_max ) then
          LOG_INFO( 'MGSmootherBase2D_Do_smoothing',*) "Smoothing did not converge in the maximum number of iterations.."
          LOG_INFO( 'MGSmootherBase2D_Do_smoothing',*) "Residual L2, max norm   : ", this%history_residual_l2(this%history_residual_count), this%history_residual_max(this%history_residual_count)
          exit
        end if
      end if
    end do
    return
  end subroutine MGSmootherBase2D_Do_smoothing

!OCL SERIAL
  subroutine MGSmootherBase2D_eval_residual_statistic( this, res, &
    res_l2, res_max )
    implicit none
    class(MGSmootherBase2D), intent(in) :: this
    type(MeshField2D), intent(in), target :: res
    real(RP), intent(out) :: res_l2
    real(RP), intent(out) :: res_max

    integer :: ldom
    class(MeshBase2D), pointer :: mesh2D
    class(LocalMesh2D), pointer :: lmesh2D

    real(RP) :: res_l2_lc, res_max_lc
    !---------------------------------

    mesh2D => res%mesh
    res_l2_lc = 0.0_RP; res_max_lc = 0.0_RP

    do ldom=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(ldom)
      call MGSmootherBase_get_statistic_lc( res_l2_lc, res_max_lc, &
        res%local(ldom)%val, lmesh2D, mesh2D%refElem2D )
    end do
    call this%Eval_global_statistic( res_l2, res_max, &
      res_l2_lc, res_max_lc )
    return
  end subroutine MGSmootherBase2D_eval_residual_statistic

!OCL SERIAL
  subroutine MGSmootherBase2D_set_initial_residual( this, res )
    implicit none
    class(MGSmootherBase2D), intent(inout) :: this
    type(MeshField2D), intent(in) :: res
    !----------------------------------------
    call this%Eval_residual_statistic( res, & ! (in)
      this%history_residual_l2_initial, & ! (out)
      this%history_residual_max_initial ) ! (out)
    return
  end subroutine MGSmootherBase2D_set_initial_residual
  
  !> Put residual history for 2D smoother.
  !! To update internal counters, this subroutine should be called after each smoothing iteration 
  !! even when the evaluation of residual is not performed. 
!OCL SERIAL
  subroutine MGSmootherBase2D_put_residual( this, res, conv_flag)
    implicit none
    class(MGSmootherBase2D), intent(inout) :: this
    type(MeshField2D), intent(in), target :: res !< Residual field
    logical, intent(out) :: conv_flag

    real(RP) :: res_l2, res_max
    logical :: is_residual_check_step
    !---------------------------------

    is_residual_check_step = this%Is_residual_check_step()
    if ( is_residual_check_step ) then
      call this%Eval_residual_statistic( res, & ! (in)
        res_l2, res_max ) ! (out)
    end if
    call MGSmootherBase_put_residual_history( this, &
      conv_flag,                              & ! (out)
      is_residual_check_step, res_l2, res_max ) ! (in)
    return
  end subroutine MGSmootherBase2D_put_residual

!-- 3D -----------------------------------------

!OCL SERIAL
  subroutine MGSmootherBase3D_Init(this)
    implicit none
    class(MGSmootherBase3D), intent(inout) :: this
    !---------------------------------
    call MGSmootherBase_Init( this )
    return
  end subroutine MGSmootherBase3D_Init

!OCL SERIAL
  subroutine MGSmootherBase3D_Final(this)
    implicit none
    class(MGSmootherBase3D), intent(inout) :: this
    !---------------------------------
    call MGSmootherBase_Final( this )
    return
  end subroutine MGSmootherBase3D_Final

!OCL SERIAL
  subroutine MGSmootherBase3D_Do_smoothing( this, q, res,      &
    f, aux_var, var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,  &
    mg_p_level, mg_h_level,                                    &
    pre_or_post_smooth )
    implicit none
    class(MGSmootherBase3D), intent(inout) :: this
    type(MeshField3D), intent(inout), target :: q
    type(MeshField3D), intent(inout) :: res
    type(MeshField3D), intent(in) :: f
    type(MeshField3D), intent(inout), target :: aux_var(:)
    class(MeshFieldCommBase), intent(inout) :: var_comm
    class(MeshFieldCommBase), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    class(MeshBase3D), intent(in), target :: mesh3D
    integer, intent(in) :: mg_p_level
    integer, intent(in) :: mg_h_level
    integer, intent(in) :: pre_or_post_smooth

    integer :: itr
    logical :: zero_initial_guess
    logical :: cal_res_flag
    logical :: conv_flag
    !--------------------------------------

    this%history_residual_count = 0
    this%history_residual_rstep = this%history_residual_step
    
    do itr=1, this%num_smooth_itr_max
      cal_res_flag = this%Is_residual_check_step()
      zero_initial_guess = ( pre_or_post_smooth == MGSmoother_PRE_ID .and. mg_p_level /= 1 .and. itr == 1 )

      call this%Advance_itr_1step( q, res,                              &
        f, aux_var, itr, var_comm, aux_comm, Dx, Dy, Dz, Lift, mesh3D,  &
        cal_res_flag, zero_initial_guess, mg_p_level, mg_h_level,       &
        pre_or_post_smooth )

      call this%put_residual( res, conv_flag )
      if ( cal_res_flag ) then
        if ( conv_flag ) then
          exit
        else if ( itr == this%num_smooth_itr_max ) then
          LOG_INFO( 'MGSmootherBase3D_Do_smoothing',*) "Smoothing did not converge in the maximum number of iterations.."
          LOG_INFO( 'MGSmootherBase3D_Do_smoothing',*) "Residual L2, max norm   : ", this%history_residual_l2(this%history_residual_count), this%history_residual_max(this%history_residual_count)
          exit
        end if
      end if
    end do
    return
  end subroutine MGSmootherBase3D_Do_smoothing

!OCL SERIAL
  subroutine MGSmootherBase3D_eval_residual_statistic( this, res, &
    res_l2, res_max )
    implicit none
    class(MGSmootherBase3D), intent(in) :: this
    type(MeshField3D), intent(in) :: res
    real(RP), intent(out) :: res_l2
    real(RP), intent(out) :: res_max

    integer :: ldom
    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lmesh3D

    real(RP) :: res_l2_lc, res_max_lc
    !---------------------------------

    mesh3D => res%mesh
    res_l2_lc = 0.0_RP; res_max_lc = 0.0_RP

    do ldom=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(ldom)
      call MGSmootherBase_get_statistic_lc( res_l2_lc, res_max_lc, &
        res%local(ldom)%val, lmesh3D, mesh3D%refElem3D )
    end do
    call this%Eval_global_statistic( res_l2, res_max, &
      res_l2_lc, res_max_lc )
    return
  end subroutine MGSmootherBase3D_eval_residual_statistic

!OCL SERIAL
  subroutine MGSmootherBase3D_set_initial_residual( this, res )
    implicit none
    class(MGSmootherBase3D), intent(inout) :: this
    type(MeshField3D), intent(in) :: res
    !----------------------------------------
    call this%Eval_residual_statistic( res, & ! (in)
      this%history_residual_l2_initial, & ! (out)
      this%history_residual_max_initial ) ! (out)
    return
  end subroutine MGSmootherBase3D_set_initial_residual

!OCL SERIAL
  subroutine MGSmootherBase3D_put_residual( this, res, conv_flag )
    implicit none
    class(MGSmootherBase3D), intent(inout) :: this
    type(MeshField3D), intent(in), target :: res
    logical, intent(out) :: conv_flag

    real(RP) :: res_l2, res_max
    logical :: is_residual_check_step
    !---------------------------------

    is_residual_check_step = this%Is_residual_check_step()
    if ( is_residual_check_step ) then
      call this%Eval_residual_statistic( res, & ! (in)
        res_l2, res_max ) ! (out)
    end if
    call MGSmootherBase_put_residual_history( this, &
      conv_flag,                              & ! (out)
      is_residual_check_step, res_l2, res_max ) ! (in)
    return
  end subroutine MGSmootherBase3D_put_residual

!-- private subroutines ----------------------------------------------------

  !> Put residual history and evaluate whether convergence is achieved.
  !! This subroutine also updates internal counters.
!OCL SERIAL
  subroutine MGSmootherBase_put_residual_history( this, conv_flag, &
    is_residual_check_step, res_l2, res_max )
    implicit none
    class(MGSmootherBase), intent(inout) :: this
    logical, intent(out) :: conv_flag
    logical, intent(in) :: is_residual_check_step
    real(RP), intent(in) :: res_l2
    real(RP), intent(in) :: res_max
    !---------------------------------

    if ( is_residual_check_step ) then
      !- Store the statistics of residual to history
      this%history_residual_count = this%history_residual_count + 1
      this%history_residual_l2 (this%history_residual_count) = res_l2
      this%history_residual_max(this%history_residual_count) = res_max

      !- Reset remaining step
      this%history_residual_rstep = this%history_residual_step

      !- Evaluate convergence
      conv_flag = this%Is_converged( .false. )
    else
      this%history_residual_rstep = this%history_residual_rstep - 1
      conv_flag = .false.
    end if
    return
  end subroutine MGSmootherBase_put_residual_history

  !> Get residual statistics on local mesh
!OCL SERIAL
  subroutine MGSmootherBase_get_statistic_lc( res_l2_lc, res_max_lc, &
    res_lc, lmesh, elem )
    use scale_localmesh_base, only: LocalMeshBase
    use scale_element_base, only: ElementBase
    implicit none
    class(LocalMeshBase), intent(in), target :: lmesh
    class(ElementBase), intent(in), target :: elem
    real(RP), intent(inout) :: res_l2_lc
    real(RP), intent(inout) :: res_max_lc
    real(RP), intent(in) :: res_lc(elem%Np,lmesh%NeA)

    integer :: ke
    real(RP) :: res_(elem%Np)
    !---------------------------------

    !$omp parallel do private(ke,res_) reduction(+:res_l2_lc) reduction(max:res_max_lc)
    do ke=lmesh%NeS, lmesh%NeE
      res_(:) = res_lc(:,ke)
      res_l2_lc = res_l2_lc + sum( res_(:)**2 )
      res_max_lc = max(res_max_lc, maxval( abs(res_(:)) ))
    end do
    return
  end subroutine MGSmootherBase_get_statistic_lc

!OCL SERIAL
  subroutine MGSmootherBase_eval_global_statistic( this, res_l2_global, res_max_global, &
    res_l2_lc, res_max_lc)
    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD
    implicit none
    class(MGSmootherBase), intent(in) :: this
    real(RP), intent(out) :: res_l2_global
    real(RP), intent(out) :: res_max_global
    real(RP), intent(in) :: res_l2_lc
    real(RP), intent(in) :: res_max_lc

    integer :: ierr
    !---------------------------------

    !- Gather L2 residual from all processes
    call MPI_Allreduce( res_l2_lc, res_l2_global, &
      1, MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr )
      
    res_l2_global = sqrt( res_l2_global )

    !- Gather max residual from all processes
    call MPI_Allreduce( res_max_lc, res_max_global, &
      1, MPI_DOUBLE_PRECISION, MPI_MAX, PRC_LOCAL_COMM_WORLD, ierr )
    return
  end subroutine MGSmootherBase_eval_global_statistic

end module scale_multigrid_smoother_base
