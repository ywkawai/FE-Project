#include "scaleFElib.h"
program test_advect1d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC, only: &
    CZ             => ATMOS_GRID_CARTESC_CZ,             &  
    CDZ             => ATMOS_GRID_CARTESC_CDZ,           &
    RCDZ            => ATMOS_GRID_CARTESC_RCDZ

  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate, FILE_HISTORY_put, FILE_HISTORY_write         

  use scale_time_manager, only: &
    TIME_manager_advance,                              &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,        &
    TIME_DTSEC, TIME_NSTEP 
  
  use scale_timeint_rk, only: &
    timeint_rk 

  use mod_fieldutil, only: &
    get_upwind_pos1d => fieldutil_get_upwind_pos1d,         &
    get_profile1d_tracer => fieldutil_get_profile1d_tracer
  
  use mod_operator_fvm, only: &
    operator_fvm    
   
  !-----------------------------------------------------------------------------
  implicit none

  integer :: NeGX, GXHALO
  integer, parameter :: NLocalMeshPerPrc = 1

  ! The type of initial q (sin, gaussian-hill, cosine-bell, top-hat)
  character(len=H_SHORT) :: InitShapeName
  real(RP) :: InitShapeParams(2)
  ! The type of specified velocify field (constant)
  real(RP) :: ADV_VEL

  real(RP), parameter :: dom_xmin =  0.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP

  character(len=H_SHORT) :: FLUX_SCHEME_TYPE
  type(operator_fvm) :: optr_fvm

  character(len=H_SHORT) :: TINTEG_SCHEME_TYPE
  type(timeint_rk) :: tinteg 
  integer :: nowstep
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1

  real(RP) :: tsec_
  character(len=H_MID), parameter :: APPNAME = "advect1d with FVM"

  integer :: HST_ID(2)
  real(RP), allocatable :: q(:,:,:)
  real(RP), allocatable :: qexact(:,:,:)
  real(RP), allocatable :: u(:,:,:)

  integer :: nstep_eval_error  
  !-------------------------------------------------------

  call init()
  call PROF_rapstart( 'set_initcond', 1)  
  call set_initcond()
  call PROF_rapend( 'set_initcond', 1)  

  call PROF_rapstart( 'TimeLoop', 1)
  do nowstep=1, TIME_NSTEP
    do rkstage=1, tinteg%nstage
      !* Exchange halo data
      call PROF_rapstart( 'exchange_halo', 1)
      call excahge_halo( q(:,IS,JS) )
      call PROF_rapend( 'exchange_halo', 1)

      !* Update prognostic variables
      call PROF_rapstart( 'cal_dyn_tend', 1)
      tintbuf_ind = tinteg%tend_buf_indmap(rkstage)
      call cal_dyn_tend( tinteg%tend_buf1D_ex(:,RKVAR_Q,tintbuf_ind), q, u )
      call PROF_rapend( 'cal_dyn_tend', 1) 

      call PROF_rapstart( 'update_var', 1)
      call tinteg%Advance(rkstage, q(:,IS,JS), RKVAR_Q, KS, KE)
      call PROF_rapend('update_var', 1)
    end do
    
    !* Advance time
    call TIME_manager_advance()

    tsec_ = TIME_NOWDATE(6) + TIME_NOWSUBSEC
    if (mod(nowstep,nstep_eval_error) == 0) then 
      LOG_PROGRESS('(A,F13.5,A)') "t=", real(tsec_), "[s]"
      call evaluate_error(tsec_)
    end if
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )

    !* Output
    call FILE_HISTORY_put(HST_ID(1), q(KS:KE,IS,JS))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS:KE,IS,JS))
    call FILE_HISTORY_write()
  end do
  call PROF_rapend( 'TimeLoop', 1)    

  call final()

contains
subroutine cal_dyn_tend( dqdt, q_, u_ )

  implicit none

  real(RP), intent(out) :: dqdt(KA,IS:IS,JS:JS)
  real(RP), intent(in) :: q_(KA,IA,JA)
  real(RP), intent(in)  :: u_(KA,IA,JA)

  integer :: k, i, j
  real(RP) :: qflux(KA,IA,JA)

    !------------------------------------------------------------------------
    call PROF_rapstart( 'update_dyn_cal_flux', 2)
    call optr_fvm%C_flux_XYW( qflux, u_, q_ )
    call PROF_rapend( 'update_dyn_cal_flux', 2)

    call PROF_rapstart( 'cal_dyn_tend_dqdt', 2)
    do j=JS, JS
    do i=IS, IS
    do k=KS, KE
      dqdt(k,i,j) = - (qflux(k,i,j) - qflux(k-1,i,j)) * RCDZ(k)
    end do
    end do
    end do
    call PROF_rapend( 'cal_dyn_tend_dqdt', 2)

    return
  end subroutine cal_dyn_tend
 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine excahge_halo(var)
    real(RP), intent(inout) :: var(KA)
    integer :: k 

    do k=1, KHALO
      var(KS-k) = var(KE-k+1)
      var(KE+k) = var(KS+k-1)
    end do
  end subroutine excahge_halo

  subroutine evaluate_error(tsec)
    
    implicit none
    real(DP), intent(in) :: tsec

    real(RP) :: diff(KA)
    real(RP) :: exact_sol(KA)
    real(RP) :: x_uwind_1d(KA)
    integer :: k
    real(RP) :: l2error  
    real(RP) :: linferror
    !------------------------------------------------------------------------
    
    x_uwind_1d(KS:KE) =  get_upwind_pos1d(CZ(KS:KE), ADV_VEL, tsec, dom_xmin, dom_xmax)
    call get_profile1d_tracer( qexact(KS:KE,IS,JS),                & ! (out)
      initShapeName, x_uwind_1d(KS:KE), InitShapeParams, KE-KS+1 )   ! (in)    
    
    l2error   = 0.0_RP   
    linferror = 0.0_RP      
    do k=KS, KE
      l2error = l2error + CDZ(k) * ( q(k,IS,JS) - qexact(k,IS,JS) )**2
      linferror = max(linferror, abs(q(k,IS,JS) - qexact(k,IS,JS)))
    end do
    LOG_INFO("evaluate_error_l2",*) sqrt(l2error)/(dom_xmax - dom_xmin)
    LOG_INFO("evaluate_error_linf",*) linferror

  end subroutine evaluate_error

  subroutine set_initcond() 
    implicit none

    integer :: i, j

    !-----------------------------------------
    do j=JS, JE
    do i=IS, IE
      call get_profile1d_tracer( q(KS:KE,i,j),               & ! (out)
        initShapeName, CZ(KS:KE), InitShapeParams, KE-KS+1 )   ! (in)
      qexact(:,i,j) = q(:,i,j)
    end do
    end do
    u(:,:,:)      = ADV_VEL

    call FILE_HISTORY_put(HST_ID(1), q(KS:KE,IS,JS))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS:KE,IS,JS))
    call FILE_HISTORY_write()   
    
    return
  end subroutine set_initcond

  subroutine init()
    use scale_prc_cartesC, only: PRC_CARTESC_setup     
    use scale_comm_cartesC, only: COMM_setup
    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only: TIME_manager_Init 
    use scale_atmos_grid_cartesC, only: &
      ATMOS_GRID_CARTESC_allocate, &
      ATMOS_GRID_CARTESC_generate
    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_regist
    use scale_file_cartesC, only: &
      FILE_CARTESC_setup
    use mod_output_fvm, only: output_fvm_setup
    use scale_file_history, only: FILE_HISTORY_reg

    implicit none

    namelist /PARAM_TEST/ &
      NeGX, GXHALO,                         &
      FLUX_SCHEME_TYPE, TINTEG_SCHEME_TYPE, &
      InitShapeName, InitShapeParams,       &
      ADV_VEL,                              &
      nstep_eval_error
        
    integer :: comm, myrank, nprocs
    logical :: ismaster
    real(RP) :: del

    integer :: ierr
    !----------------------------------------------
    
    ! scale setup
    call SCALE_init( APPNAME )  
    call PROF_rapstart( "init", 1 )
    
    !--- read namelist

    NeGX = 2; GXHALO = 2
    FLUX_SCHEME_TYPE = 'CD2'
    TINTEG_SCHEME_TYPE = 'ERK_SSP_3s3o'
    InitShapeName    = 'sin'
    InitShapeParams  = (/ 1.0_RP, 0.0_RP /)
    ADV_VEL          = 1.0_RP
    nstep_eval_error = 5

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TEST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_TEST. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_TEST)

    ! setup process
    call PRC_CARTESC_setup

    call ATMOS_GRID_CARTESC_INDEX_setup( &
      KMAX=NeGX, KHALO=GXHALO, IMAX=1, IHALO=1, JMAX=1, JHALO=1, &
      IBLOCK=1, JBLOCK=1 ) 

    ! setup horizontal/veritical grid system
    call ATMOS_GRID_CARTESC_allocate
    del = (dom_xmax - dom_xmin)/real(NeGX,kind=RP)
    call ATMOS_GRID_CARTESC_generate( DZ=del, DX=del, DY=del )

    ! setup file I/O
    call FILE_CARTESC_setup

    ! setup mpi communication
    call COMM_setup

    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init
    
    ! setup a module for FVM operator
    call optr_fvm%Init( FLUX_SCHEME_TYPE, KS, KE, KA, IS, IS, IA, JS, JS, JA )
    
    ! setup a module for time integrator
    call tinteg%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1, &
                      1, (/ KA /) )

    ! setup variables and history files
    allocate( q(KA,IA,JA), qexact(KA,IA,JA), u(KA,IA,JA) )

    call output_fvm_setup('linedom1d')
    call FILE_HISTORY_reg( "q", "q", "1", HST_ID(1), dim_type='X')
    call FILE_HISTORY_reg( "qexact", "qexact", "1", HST_ID(2), dim_type='X')

    call PROF_rapend( "init", 1 )    
    return
  end subroutine init

  subroutine final()
    use mod_output_fvm, only: output_fvm_finalize
    use scale_time_manager, only: TIME_manager_Final    
    implicit none

    call PROF_rapstart( "final", 1 )

    call optr_fvm%Final()
    call output_fvm_finalize
    call TIME_manager_Final
    
    call PROF_rapend( "final", 1 )  
    call SCALE_finalize()

    return
  end subroutine final

end program test_advect1d
