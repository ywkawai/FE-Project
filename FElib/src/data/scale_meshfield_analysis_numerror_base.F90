!-------------------------------------------------------------------------------
!> module FElib / Data / Statistics / numerical error
!!
!! @par Description
!!           This module provides a base class useful for evaluating numerical errrors
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_meshfield_analysis_numerror_base
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: &
    PRC_ismaster, PRC_abort
  use scale_time_manager, only: &
    TIME_NSTEP

  use scale_element_base, only: ElementBase
  use scale_element_line, only: LineElement
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_mesh_base, only: MeshBase
  use scale_meshfield_base, only: MeshFieldBase
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, public :: MeshFieldAnalysisNumerrorBase
    integer :: var_num
    integer :: log_fid
    integer :: log_step_interval
    integer :: log_rstep
    logical :: output_error_first

    integer :: PolyOrderErrorCheck
    integer :: intrp_np
    integer :: ndim
    real(RP), allocatable :: IntrpMat(:,:)
    real(RP), allocatable :: intw_intrp(:)
    real(RP), allocatable :: epos_intrp(:,:) 
  contains
    procedure :: Init_base => meshfield_analysis_numerror_base_Init 
    generic :: Init => Init_base
    procedure :: Final => meshfield_analysis_numerror_base_Final
    procedure :: Regist => meshfield_analysis_numerror_base_regist
    procedure :: Evaluate_base => meshfield_analysis_numerror_base_evaluate
    !
    procedure :: Evaluate_error_lc => meshfield_analysis_numerror_base_evaluate_error_lc
    procedure :: Evaluate_covariance_lc => meshfield_analysis_numerror_base_evaluate_covariance_lc
  end type MeshFieldAnalysisNumerrorBase

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  private :: numerror_do_step

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains

  !-----------------------------------------------------------------------------
  !> Setup
!OCL SERIAL
  subroutine meshfield_analysis_numerror_base_Init( this,                 &
    porder_error_check, ndim, np, intrp_np, log_fname_base, log_step_interval )
    implicit none
    class(MeshFieldAnalysisNumerrorBase), intent(inout) :: this
    integer, intent(in) :: porder_error_check
    integer, intent(in) :: ndim
    integer, intent(in) :: np
    integer, intent(in) :: intrp_np
    character(len=*), intent(in) :: log_fname_base
    integer, intent(in) :: log_step_interval

    character(len=H_MID) :: fname
    integer :: ierr
    !---------------------------------------------------------------------------

    this%var_num = 0

    this%log_step_interval = log_step_interval
    this%log_rstep         = log_step_interval

    this%output_error_first = .true. 

    !--
    if ( PRC_ismaster ) then
      this%log_fid = IO_get_available_fid()
      fname = trim(log_fname_base)//".peall"
      open( unit   = this%log_fid,                 &
        file   = fname,                            &
        form   = 'formatted',                      &
        iostat = ierr                              )
      if ( ierr /= 0 ) then
        LOG_ERROR('MeshField_NumErrorAnalysis_Init',*) 'File open error! :', trim(fname)
        call PRC_abort
      endif 
    end if

    !--
    this%PolyOrderErrorCheck = porder_error_check
    this%intrp_np = intrp_np
    this%ndim = ndim
    allocate( this%IntrpMat(this%intrp_np,np) )
    allocate( this%intw_intrp(this%intrp_np) )
    allocate( this%epos_intrp(this%intrp_np,this%ndim) )

    return
  end subroutine meshfield_analysis_numerror_base_Init

  !-----------------------------------------------------------------------------
  !> Finalization
!OCL SERIAL
  subroutine meshfield_analysis_numerror_base_Final( this )
    implicit none
    class(MeshFieldAnalysisNumerrorBase), intent(inout) :: this
    !---------------------------------------------------------------------------
    
    deallocate( this%IntrpMat, this%intw_intrp, this%epos_intrp )

    close( this%log_fid )
    this%log_fid = -1
    
    return
  end subroutine meshfield_analysis_numerror_base_Final

!OCL SERIAL
  subroutine meshfield_analysis_numerror_base_regist( this, varname, unit, varid )
    implicit none
    class(MeshFieldAnalysisNumerrorBase), intent(inout) :: this
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: unit
    integer, intent(out) :: varid
    !---------------------------------------------------------------------------

    this%var_num = this%var_num + 1
    varid = this%var_num

    if ( PRC_ismaster ) then
      write(this%log_fid,'(A25)',advance='no') 'L1_error   ('//trim(varname)//')'
      write(this%log_fid,'(A25)',advance='no') 'L2_error   ('//trim(varname)//')'
      write(this%log_fid,'(A25)',advance='no') 'Linf_error ('//trim(varname)//')'
      write(this%log_fid,'(A25)',advance='no') 'Ediss      ('//trim(varname)//')'
      write(this%log_fid,'(A25)',advance='no') 'Edisp      ('//trim(varname)//')'
    end if

    return
  end subroutine meshfield_analysis_numerror_base_regist

!OCL SERIAL
  subroutine meshfield_analysis_numerror_base_evaluate( &
    this, tstep, tsec, dom_vol, evaluate_error, calc_covariance )
    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD
    implicit none
    class(MeshFieldAnalysisNumerrorBase), intent(inout) :: this
    integer, intent(in) :: tstep
    real(RP) :: tsec
    real(RP), intent(in) :: dom_vol
    interface 
      subroutine evaluate_error( this_, &
        num_error_l1_lc_, num_error_l2_lc_, num_error_linf_lc_, &
        numsol_mean, exactsol_mean )
        import RP
        import MeshFieldAnalysisNumerrorBase
        class(MeshFieldAnalysisNumerrorBase), intent(in) :: this_
        real(RP), intent(inout) :: num_error_l1_lc_(this_%var_num)
        real(RP), intent(inout) :: num_error_l2_lc_(this_%var_num)
        real(RP), intent(inout) :: num_error_linf_lc_(this_%var_num)    
        real(RP), intent(inout) :: numsol_mean(this_%var_num) 
        real(RP), intent(inout) :: exactsol_mean(this_%var_num) 
      end subroutine evaluate_error

      subroutine calc_covariance( this_, &
        cov_numsol_numsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc, &
        numsol_mean, exactsol_mean                      )
        import RP
        import MeshFieldAnalysisNumerrorBase
        class(MeshFieldAnalysisNumerrorBase), intent(in) :: this_
        real(RP), intent(inout) :: cov_numsol_numsol_lc(this_%var_num)
        real(RP), intent(inout) :: cov_numsol_exactsol_lc(this_%var_num)
        real(RP), intent(inout) :: cov_exactsol_exactsol_lc(this_%var_num)
        real(RP), intent(in) :: numsol_mean(this_%var_num) 
        real(RP), intent(in) :: exactsol_mean(this_%var_num)   
      end subroutine calc_covariance   
    end interface

    logical :: do_check_numerror

    integer, parameter :: SUM_BUF_L1_ID = 1
    integer, parameter :: SUM_BUF_L2_ID = 2
    integer, parameter :: SUM_BUF_MEAN_NUMSOL_ID   = 3
    integer, parameter :: SUM_BUF_MEAN_EXACTSOL_ID = 4
    integer, parameter :: SUM_BUF_NUM   = 4

    real(RP) :: num_error_linf_lc(this%var_num)
    real(RP) :: sum_buf_lc(this%var_num,SUM_BUF_NUM)
    real(RP) :: num_error_linf(this%var_num)
    real(RP) :: sum_buf(this%var_num,SUM_BUF_NUM)

    real(RP) :: covariance_lc(this%var_num,3)
    real(RP) :: covariance(this%var_num,3)
    real(RP) :: Ediss(this%var_num), Edisp(this%var_num)

    integer :: iv
    integer :: ierr
    !---------------------------------------------------------------------------

    if ( this%var_num == 0 ) return
    
    if ( this%output_error_first ) then
      if ( PRC_ismaster ) write(this%log_fid,*)
      this%output_error_first = .false.
    end if

    do_check_numerror = numerror_do_step( this, tstep )
    if ( .not. do_check_numerror  ) then
      !if ( PRC_ismaster .and. TIME_NOWSTEP == TIME_NSTEP ) close( NUMERROR_LOG_FID )
      return
    end if

    !---
    num_error_linf_lc(:) = 0.0_RP
    sum_buf_lc(:,:)      = 0.0_RP

    call evaluate_error( this, &
      sum_buf_lc(:,SUM_BUF_L1_ID), sum_buf_lc(:,SUM_BUF_L2_ID), num_error_linf_lc(:), &
      sum_buf_lc(:,SUM_BUF_MEAN_NUMSOL_ID), sum_buf_lc(:,SUM_BUF_MEAN_EXACTSOL_ID)    )

    call MPI_Allreduce( sum_buf_lc(:,:), sum_buf(:,:), &
      this%var_num * SUM_BUF_NUM,                      &
      MPI_DOUBLE_PRECISION,   &
      MPI_SUM,                &
      PRC_LOCAL_COMM_WORLD,   &
      ierr                    )
    call MPI_Allreduce( num_error_linf_lc(:), num_error_linf(:), &
      this%var_num,           &
      MPI_DOUBLE_PRECISION,   &
      MPI_MAX,                &
      PRC_LOCAL_COMM_WORLD,   &
      ierr                    )

    sum_buf(:,SUM_BUF_MEAN_NUMSOL_ID) = sum_buf(:,SUM_BUF_MEAN_NUMSOL_ID) / dom_vol
    sum_buf(:,SUM_BUF_MEAN_EXACTSOL_ID) = sum_buf(:,SUM_BUF_MEAN_EXACTSOL_ID) / dom_vol
  
    !--
    covariance_lc(:,:) = 0.0_RP

    call calc_covariance( this, &
      covariance_lc(:,1), covariance_lc(:,2), covariance_lc(:,3),            &
      sum_buf(:,SUM_BUF_MEAN_NUMSOL_ID), sum_buf(:,SUM_BUF_MEAN_EXACTSOL_ID) )
    
    call MPI_Allreduce( covariance_lc(:,:), covariance(:,:), &
      this%var_num * 3,                                      &
      MPI_DOUBLE_PRECISION,   &
      MPI_SUM,                &
      PRC_LOCAL_COMM_WORLD,   &
      ierr                    )

    if ( PRC_ismaster ) then
      
      Ediss(:) = ( sqrt(covariance(:,1) / dom_vol) -  sqrt(covariance(:,3) / dom_vol) )**2      &
               + ( sum_buf(:,SUM_BUF_MEAN_NUMSOL_ID) - sum_buf(:,SUM_BUF_MEAN_EXACTSOL_ID) )**2
      
      Edisp(:) = 2.0_RP * (   sqrt( covariance(:,1)  *  covariance(:,3) ) &
                            - covariance(:,2)                             ) / dom_vol

      write(this%log_fid,'(A,ES18.8)',advance='no') 'tsec=', tsec
      do iv=1, this%var_num
       write(this%log_fid,'(A,ES18.8)',advance='no') ' ', sum_buf(iv,SUM_BUF_L1_ID) / dom_vol
       write(this%log_fid,'(A,ES18.8)',advance='no') ' ', sqrt(sum_buf(iv,SUM_BUF_L2_ID) / dom_vol)
       write(this%log_fid,'(A,ES18.8)',advance='no') ' ', num_error_linf(iv)
       write(this%log_fid,'(A,ES18.8)',advance='no') ' ', Ediss(iv)
       write(this%log_fid,'(A,ES18.8)',advance='no') ' ', Edisp(iv)
      end do
      write(this%log_fid,*)

!      if ( tstep > TIME_NSTEP ) close( NUMERROR_LOG_FID )
    end if

    return
  end subroutine meshfield_analysis_numerror_base_evaluate

!OCL SERIAL
  subroutine meshfield_analysis_numerror_base_evaluate_error_lc( base, &
    num_error_l1_lc, num_error_l2_lc, num_error_linf_lc, &
    numsol_mean_lc, exactsol_mean_lc,                    &
    q, qexact, qexact_intrp, lcmesh, elem                )
    implicit none
    class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
    class(LocalMeshBase), intent(in) :: lcmesh
    class(ElementBase), intent(in) :: elem
    real(RP), intent(inout) :: num_error_l1_lc(base%var_num)
    real(RP), intent(inout) :: num_error_l2_lc(base%var_num)
    real(RP), intent(inout) :: num_error_linf_lc(base%var_num)
    real(RP), intent(inout) :: numsol_mean_lc(base%var_num) 
    real(RP), intent(inout) :: exactsol_mean_lc(base%var_num) 
    real(RP), intent(in) :: q(elem%Np,lcmesh%Ne,base%var_num)
    real(RP), intent(in) :: qexact(elem%Np,lcmesh%Ne,base%var_num)
    real(RP), intent(in) :: qexact_intrp(base%intrp_np,lcmesh%Ne,base%var_num)
    
    integer :: iv
    integer :: ke
    
    real(RP) :: linf_max_tmp(lcmesh%Ne,base%var_num)
    real(RP) :: q_intrp(base%intrp_np)
    real(RP) :: JGsqrtxIntw_intrp(base%intrp_np)
    !---------------------------------------------------------------------------

    do iv=1, base%var_num
      !$omp parallel do private(ke, q_intrp, JGsqrtxIntw_intrp) reduction(+: num_error_l1_lc, num_error_l2_lc, numsol_mean_lc, exactsol_mean_lc )
      do ke=lcmesh%NeS, lcmesh%NeE
        q_intrp(:) = matmul( base%IntrpMat, q(:,ke,iv) )
        JGsqrtxIntw_intrp(:) = matmul( base%IntrpMat, lcmesh%J(:,ke) * lcmesh%Gsqrt(:,ke) )
        JGsqrtxIntw_intrp(:) = JGsqrtxIntw_intrp(:) * base%intw_intrp(:)

        num_error_l1_lc(iv) = num_error_l1_lc(iv) + sum( JGsqrtxIntw_intrp(:) * abs( q_intrp(:) - qexact_intrp(:,ke,iv) ) )
        num_error_l2_lc(iv) = num_error_l2_lc(iv) + sum( JGsqrtxIntw_intrp(:) * ( q_intrp(:) - qexact_intrp(:,ke,iv) )**2 )
        linf_max_tmp(ke,iv) = maxval(abs(q(:,ke,iv) - qexact(:,ke,iv)))

        numsol_mean_lc(iv) = numsol_mean_lc(iv) + sum( JGsqrtxIntw_intrp(:) * q_intrp(:) )
        exactsol_mean_lc(iv) = exactsol_mean_lc(iv) + sum( JGsqrtxIntw_intrp(:) * qexact_intrp(:,ke,iv) )
      end do
      
      num_error_linf_lc(iv) = max( num_error_linf_lc(iv), maxval( linf_max_tmp(:,iv) ) )
    end do

    return
  end subroutine meshfield_analysis_numerror_base_evaluate_error_lc

!OCL SERIAL
  subroutine meshfield_analysis_numerror_base_evaluate_covariance_lc( base, &
    cov_numsol_numsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc, &
    q, qexact, qexact_intrp, numsol_mean, exactsol_mean,                    &
    lcmesh, elem                                                            )
  
    implicit none
    class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
    class(LocalMeshBase), intent(in) :: lcmesh
    class(ElementBase), intent(in) :: elem
    real(RP), intent(inout) :: cov_numsol_numsol_lc(base%var_num)
    real(RP), intent(inout) :: cov_numsol_exactsol_lc(base%var_num)
    real(RP), intent(inout) :: cov_exactsol_exactsol_lc(base%var_num)
    real(RP), intent(in) :: q(elem%Np,lcmesh%Ne,base%var_num)
    real(RP), intent(in) :: qexact(elem%Np,lcmesh%Ne,base%var_num)
    real(RP), intent(in) :: qexact_intrp(base%intrp_np,lcmesh%Ne,base%var_num)
    real(RP), intent(in) :: numsol_mean(base%var_num) 
    real(RP), intent(in) :: exactsol_mean(base%var_num) 

    integer :: iv
    integer :: ke
    
    real(RP) :: dq_intrp(base%intrp_np)
    real(RP) :: JGsqrtxIntw_intrp(base%intrp_np)
    real(RP) :: dq_exact_intrp(base%intrp_np)
    !---------------------------------------------------------------------------

    do iv=1, base%var_num
      !$omp parallel do private(ke, dq_intrp, dq_exact_intrp, JGsqrtxIntw_intrp) reduction(+: cov_numsol_numsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc )
      do ke=lcmesh%NeS, lcmesh%NeE
        dq_intrp(:) = matmul( base%IntrpMat, q(:,ke,iv) ) - numsol_mean(iv)
        dq_exact_intrp(:) = qexact_intrp(:,ke,iv)- exactsol_mean(iv)
        JGsqrtxIntw_intrp(:) = matmul( base%IntrpMat, lcmesh%J(:,ke) * lcmesh%Gsqrt(:,ke) )
        JGsqrtxIntw_intrp(:) = JGsqrtxIntw_intrp(:) * base%intw_intrp(:)

        cov_numsol_numsol_lc(iv)   = cov_numsol_numsol_lc(iv) + sum( JGsqrtxIntw_intrp(:) * dq_intrp(:) * dq_intrp(:) )
        cov_numsol_exactsol_lc(iv) = cov_numsol_exactsol_lc(iv) + sum( JGsqrtxIntw_intrp(:) * dq_intrp(:) * dq_exact_intrp(:) )
        cov_exactsol_exactsol_lc(iv) = cov_exactsol_exactsol_lc(iv) + sum( JGsqrtxIntw_intrp(:) * dq_exact_intrp(:) * dq_exact_intrp(:) )
      end do
    end do

    return
  end subroutine meshfield_analysis_numerror_base_evaluate_covariance_lc


!--- private

!OCL SERIAL
  function numerror_do_step( this, step ) result(do_flag)
    implicit none

    type(MeshFieldAnalysisNumerrorBase), intent(inout) :: this
    integer, intent(in) :: step
    logical :: do_flag    
    !------------------------------------------------------------------------

    if ( step == 1 .or. TIME_NSTEP < step ) then
      do_flag = .true.
      return
    end if

    this%log_rstep = this%log_rstep - 1
    if ( this%log_rstep == 0 ) then
      do_flag = .true.
      this%log_rstep  = this%log_step_interval
    else
      do_flag = .false. 
    end if

    return
  end function numerror_do_step

end module scale_meshfield_analysis_numerror_base