#include "scaleFElib.h"
module scale_meshfield_analysis_numerror
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

  use scale_element_base, only: &
    ElementBase, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField1D, MeshField2D, MeshField3D
    
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
  end type MeshFieldAnalysisNumerrorBase

  type, public, extends(MeshFieldAnalysisNumerrorBase) :: MeshFieldAnalysisNumerror3D
  contains
    procedure :: Init_3D => meshfield_analysis_numerror_3D_Init
    generic :: Init => Init_3D
    procedure :: Evaluate => meshfield_analysis_numerror_3D_evaluate
  end type MeshFieldAnalysisNumerror3D

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
      write(this%log_fid,'(A25)',advance='no') 'L1_error ('//trim(varname)//')'
      write(this%log_fid,'(A25)',advance='no') 'L2_error ('//trim(varname)//')'
      write(this%log_fid,'(A25)',advance='no') 'Linf_error ('//trim(varname)//')'
    end if

    return
  end subroutine meshfield_analysis_numerror_base_regist

!OCL SERIAL
  subroutine meshfield_analysis_numerror_base_evaluate( &
    this, tstep, tsec, dom_vol, evaluate_error          )
    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD
    implicit none
    class(MeshFieldAnalysisNumerrorBase), intent(inout) :: this
    integer, intent(in) :: tstep
    real(RP) :: tsec
    real(RP), intent(in) :: dom_vol
    interface 
      subroutine evaluate_error( this_, num_error_l1_lc_, num_error_l2_lc_, num_error_linf_lc_ )
        import RP
        import MeshFieldAnalysisNumerrorBase
        class(MeshFieldAnalysisNumerrorBase), intent(in) :: this_
        real(RP), intent(inout) :: num_error_l1_lc_(this_%var_num)
        real(RP), intent(inout) :: num_error_l2_lc_(this_%var_num)
        real(RP), intent(inout) :: num_error_linf_lc_(this_%var_num)    
      end subroutine evaluate_error
    end interface

    logical :: do_check_numerror

    real(RP) :: num_error_l1_lc(this%var_num)
    real(RP) :: num_error_l2_lc(this%var_num)
    real(RP) :: num_error_linf_lc(this%var_num)
    real(RP) :: num_error_l1(this%var_num)
    real(RP) :: num_error_l2(this%var_num)
    real(RP) :: num_error_linf(this%var_num)

    integer :: iv
    integer :: ierr
    !---------------------------------------------------------------------------

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
    num_error_l1_lc(:)   = 0.0_RP
    num_error_l2_lc(:)   = 0.0_RP
    num_error_linf_lc(:) = 0.0_RP
    call evaluate_error( this, num_error_l1_lc, num_error_l2_lc, num_error_linf_lc )

    call MPI_Allreduce( num_error_l1_lc(:), num_error_l1(:), &
      this%var_num,           &
      MPI_DOUBLE_PRECISION,   &
      MPI_SUM,                &
      PRC_LOCAL_COMM_WORLD,   &
      ierr                    )
   call MPI_Allreduce( num_error_l2_lc(:), num_error_l2(:), &
      this%var_num,           &
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
        
    if ( PRC_ismaster ) then
      write(this%log_fid,'(A,ES18.8)',advance='no') 'tsec=', tsec
      do iv=1, this%var_num
       write(this%log_fid,'(A,ES18.8)',advance='no') ' ', num_error_l1(iv) / dom_vol
       write(this%log_fid,'(A,ES18.8)',advance='no') ' ', sqrt(num_error_l2(iv) / dom_vol)
       write(this%log_fid,'(A,ES18.8)',advance='no') ' ', num_error_linf(iv)
      end do
      write(this%log_fid,*)

!      if ( tstep > TIME_NSTEP ) close( NUMERROR_LOG_FID )
    end if

    return
  end subroutine meshfield_analysis_numerror_base_evaluate

!-- 3D --

  !-----------------------------------------------------------------------------
  !> Setup
!OCL SERIAL
  subroutine meshfield_analysis_numerror_3D_Init( this, &
    porder_error_check, log_fname_base, log_step_interval, &
    refElem3D )
    implicit none
    class(MeshFieldAnalysisNumerror3D), intent(inout) :: this
    integer, intent(in) :: porder_error_check
    character(len=*), intent(in) :: log_fname_base
    integer, intent(in) :: log_step_interval
    type(HexahedralElement), intent(in) :: refElem3D
    !---------------------------------------------------------------------------

    call this%MeshFieldAnalysisNumerrorBase%Init( &
      porder_error_check, 3, refElem3D%Np, porder_error_check**3, &
      log_fname_base, log_step_interval                           )

    this%IntrpMat(:,:) = refElem3D%GenIntGaussLegendreIntrpMat( &
      this%PolyOrderErrorCheck,                                                         &  ! (in)
      this%intw_intrp, this%epos_intrp(:,1), this%epos_intrp(:,2), this%epos_intrp(:,3) )  ! (out)

    return
  end subroutine meshfield_analysis_numerror_3D_Init

!OCL SERIAL
  subroutine meshfield_analysis_numerror_3D_evaluate( &
    this, tstep, tsec, mesh, evaluate_error_lc        )
    implicit none
    class(MeshFieldAnalysisNumerror3D), intent(inout) :: this
    integer, intent(in) :: tstep
    real(RP) :: tsec
    class(MeshBase3D), target, intent(in) :: mesh
    interface 
      subroutine evaluate_error_lc( this_, q, qexact, qexact_intrp, lcmesh, elem3D, intrp_epos )
        import RP
        import ElementBase3D
        import LocalMesh3D
        import MeshFieldAnalysisNumerror3D
        class(MeshFieldAnalysisNumerror3D), intent(in) :: this_
        class(LocalMesh3D), intent(in) :: lcmesh
        class(ElementBase3D) :: elem3D
        real(RP), intent(out) :: q(elem3D%Np,lcmesh%Ne,this_%var_num)
        real(RP), intent(out) :: qexact(elem3D%Np,lcmesh%Ne,this_%var_num)
        real(RP), intent(out) :: qexact_intrp(this_%intrp_np,lcmesh%Ne,this_%var_num)
        real(RP), intent(in) :: intrp_epos(this_%intrp_np,this_%ndim)
      end subroutine evaluate_error_lc
    end interface
    !---------------------------------------------------------------------------

    call this%Evaluate_base( tstep, tsec, mesh%dom_vol, evaluate_error_core )

    return
  contains 
    subroutine evaluate_error_core( base, num_error_l1_lc, num_error_l2_lc, num_error_linf_lc )
      implicit none
      class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
      real(RP), intent(inout) :: num_error_l1_lc(base%var_num)
      real(RP), intent(inout) :: num_error_l2_lc(base%var_num)
      real(RP), intent(inout) :: num_error_linf_lc(base%var_num)

      integer :: lcdomid
      integer :: iv
      integer :: ke
      
      class(LocalMesh3D), pointer :: lcmesh
      real(RP), allocatable :: q(:,:,:)
      real(RP), allocatable :: qexact(:,:,:)
      real(RP), allocatable :: qexact_intrp(:,:,:)
      real(RP), allocatable :: linf_max_tmp(:,:)
      real(RP) :: q_intrp(base%intrp_np)
      real(RP) :: JGsqrtxIntw_intrp(base%intrp_np)
     !---------------------------------------------------------------------------

      do lcdomid=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(lcdomid)
        allocate( q(lcmesh%refElem3D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact(lcmesh%refElem3D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact_intrp(base%intrp_np,lcmesh%Ne,base%var_num) )
        allocate( linf_max_tmp(lcmesh%Ne,base%var_num) )
  
        call evaluate_error_lc( this, q, qexact, qexact_intrp, & ! (in)
          lcmesh, lcmesh%refElem3D, base%epos_intrp            ) ! (out)
  
        do iv=1, this%var_num
          !$omp parallel do private(ke, q_intrp, JGsqrtxIntw_intrp) reduction(+: num_error_l1_lc, num_error_l2_lc )
          do ke=lcmesh%NeS, lcmesh%NeE
            q_intrp(:) = matmul( this%IntrpMat, q(:,ke,iv) )
            JGsqrtxIntw_intrp(:) = matmul( this%IntrpMat, lcmesh%J(:,ke) * lcmesh%Gsqrt(:,ke) )
            JGsqrtxIntw_intrp(:) = JGsqrtxIntw_intrp(:) * this%intw_intrp(:)

            num_error_l1_lc(iv) = num_error_l1_lc(iv) + sum( JGsqrtxIntw_intrp(:) * abs( q_intrp(:) - qexact_intrp(:,ke,iv) ) )
            num_error_l2_lc(iv) = num_error_l2_lc(iv) + sum( JGsqrtxIntw_intrp(:) * ( q_intrp(:) - qexact_intrp(:,ke,iv) )**2 )
            linf_max_tmp(ke,iv) = maxval(abs(q(:,ke,iv) - qexact(:,ke,iv)))
          end do
          
          num_error_linf_lc(iv) = max( num_error_linf_lc(iv), maxval( linf_max_tmp(:,iv) ) )
        end do
  
        deallocate( q, qexact, qexact_intrp, linf_max_tmp )
      end do    
    end subroutine evaluate_error_core
  end subroutine meshfield_analysis_numerror_3D_evaluate

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

end module scale_meshfield_analysis_numerror