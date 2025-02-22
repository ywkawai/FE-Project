!-------------------------------------------------------------------------------
!> Program A sample program: 1-dimensional linear advection test using FVM
!! 
!! 
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
program test_advect1d_fvm
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
    CZ             => ATMOS_GRID_CARTESC_CZ,   &  
    FZ             => ATMOS_GRID_CARTESC_FZ,   &
    CDZ            => ATMOS_GRID_CARTESC_CDZ,  &
    RCDZ           => ATMOS_GRID_CARTESC_RCDZ

  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate, &
    FILE_HISTORY_put, FILE_HISTORY_write         

  use scale_time_manager, only: &
    TIME_manager_checkstate, TIME_manager_advance,     &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,        &
    TIME_DTSEC, TIME_NSTEP, TIME_DOresume, TIME_DOend
  
  use scale_timeint_rk, only: timeint_rk 
  use mod_operator_fvm, only: operator_fvm    
   
  use mod_advect1d_fvm_numerror, only: advect1d_fvm_numerror_eval
  !-----------------------------------------------------------------------------
  implicit none

  integer :: NeGX, GXHALO

  character(len=H_SHORT) :: InitShapeName   !< The type of initial profile (sin, gaussian-hill, cosine-bell, top-hat)
  real(RP), save :: InitShapeParams(2)
  real(RP) :: ADV_VEL                       !< The constant speed of advection

  type(operator_fvm) :: optr_fvm

  real(RP), allocatable :: q(:,:,:)
  real(RP), allocatable :: qexact(:,:,:)
  real(RP), allocatable :: u(:,:,:)
  integer, save :: HST_ID(2)

  real(RP) :: tsec_
  type(timeint_rk) :: tinteg 
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1

  character(len=H_MID), parameter :: APPNAME = "advect1d with FVM"
  !-------------------------------------------------------

  call init()

  do 
    !* Report current time
    call TIME_manager_checkstate

    if (TIME_DOresume) call set_initcond()

    !* Advance time
    call TIME_manager_advance()
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )

    do rkstage=1, tinteg%nstage

      !* Exchange halo data

      call PROF_rapstart( 'exchange_halo', 1)
      call excahge_halo( q(:,IS,JS) )
      call PROF_rapend( 'exchange_halo', 1)

      !* Update prognostic variables

      call PROF_rapstart( 'cal_tend', 1)
      tintbuf_ind = tinteg%tend_buf_indmap(rkstage)
      call cal_dyn_tend( tinteg%tend_buf1D_ex(:,RKVAR_Q,tintbuf_ind), q, u )
      call PROF_rapend( 'cal_tend', 1) 

      call PROF_rapstart( 'update_var', 1)
      call tinteg%Advance(rkstage, q(:,IS,JS), RKVAR_Q, KS, KE)
      call PROF_rapend('update_var', 1)
    end do

    tsec_ = TIME_DTSEC * real(TIME_NOWSTEP-1, kind=RP)
    call advect1d_fvm_numerror_eval( qexact(:,IS,JS), & ! (out)
      q(:,IS,JS), TIME_NOWSTEP, tsec_, ADV_VEL, InitShapeName, InitShapeParams, & ! (in)
      CZ, FZ, KS, KE, KA, KHALO                                                 ) ! (in)


    !* Output history file

    call FILE_HISTORY_put(HST_ID(1), q(KS:KE,IS,JS))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS:KE,IS,JS))
    call FILE_HISTORY_write()
  
    if (TIME_DOend) exit
  end do

  call final()

contains
  !> Calculate the tendency
  !! dqdt = - ( <u q>_k - <uq>_k-1 ) / DZ
  !!
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
    implicit none
    real(RP), intent(inout) :: var(KA)

    integer :: k 
    !-------------------------------
    do k=1, KHALO
      var(KS-k) = var(KE-k+1)
      var(KE+k) = var(KS+k-1)
    end do
    return
  end subroutine excahge_halo

  !> Set inital data
  subroutine set_initcond() 
    use mod_fieldutil, only: fieldutil_get_profile1d_tracer 
    implicit none

    integer :: i, j
    !-----------------------------------------

    do j=JS, JE
    do i=IS, IE
      call fieldutil_get_profile1d_tracer( q(KS:KE,i,j),   & ! (out)
        InitShapeName, CZ(KS:KE), InitShapeParams, KE-KS+1 ) ! (in)
    end do
    end do
    u(:,:,:) = ADV_VEL

    call advect1d_fvm_numerror_eval( qexact(:,IS,JS),                 & ! (out)
      q(:,IS,JS), 1, 0.0_RP, ADV_VEL, InitShapeName, InitShapeParams, & ! (in)
      CZ, FZ, KS, KE, KA, KHALO                                       ) ! (in)

    call FILE_HISTORY_put(HST_ID(1), q(KS:KE,IS,JS))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS:KE,IS,JS))
    call FILE_HISTORY_write()
    
    return
  end subroutine set_initcond

  !> Initialization
  subroutine init()
    use scale_prc_cartesC, only: PRC_CARTESC_setup     
    use scale_comm_cartesC, only: COMM_setup
    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only:           &
      TIME_manager_Init,                    &
      TIME_manager_report_timeintervals
    use scale_atmos_grid_cartesC, only: &
      ATMOS_GRID_CARTESC_allocate, &
      ATMOS_GRID_CARTESC_generate
    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_regist
    use scale_file_cartesC, only: &
      FILE_CARTESC_setup
    use mod_output_fvm, only: output_fvm_setup
    use scale_file_history, only: FILE_HISTORY_reg

    use mod_advect1d_fvm_numerror, only: advect1d_fvm_numerror_Init
    implicit none

    real(RP), parameter :: dom_xmin =  0.0_RP
    real(RP), parameter :: dom_xmax = +1.0_RP
    character(len=H_SHORT) :: FLUX_SCHEME_TYPE
    character(len=H_SHORT) :: TINTEG_SCHEME_TYPE

    namelist /PARAM_TEST/ &
      NeGX, GXHALO,                         &
      FLUX_SCHEME_TYPE, TINTEG_SCHEME_TYPE, &
      InitShapeName, InitShapeParams,       &
      ADV_VEL
    
    integer :: ierr

    real(RP) :: del
    !----------------------------------------------
    
    !-- setup SCALE modules

    call SCALE_init( APPNAME )  
    call PROF_rapstart( "init", 1 )
    
    !-- read namelist

    NeGX = 2; GXHALO = 2
    FLUX_SCHEME_TYPE = 'CD2'
    TINTEG_SCHEME_TYPE = 'ERK_SSP_3s3o'
    InitShapeName    = 'sin'
    InitShapeParams  = (/ 1.0_RP, 0.0_RP /)
    ADV_VEL          = 1.0_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TEST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_TEST. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_TEST)

    !-- setup process
    call PRC_CARTESC_setup

    !-- setup horizontal/veritical grid system

    call ATMOS_GRID_CARTESC_INDEX_setup( &
      KMAX=NeGX, KHALO=GXHALO, IMAX=1, IHALO=1, JMAX=1, JHALO=1, &
      IBLOCK=1, JBLOCK=1 ) 

    call ATMOS_GRID_CARTESC_allocate
    del = (dom_xmax - dom_xmin)/real(NeGX,kind=RP)
    call ATMOS_GRID_CARTESC_generate( DZ=del, DX=del, DY=del )

    !-- setup file I/O
    call FILE_CARTESC_setup

    !-- setup mpi communication
    call COMM_setup

    !-- setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init
    
    !-- setup a module for FVM operator
    call optr_fvm%Init( FLUX_SCHEME_TYPE, KS, KE, KA, IS, IS, IA, JS, JS, JA )
    
    !-- setup a module for time integrator
    call tinteg%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1, &
                      1, (/ KA /) )

    !-- setup variables and history files

    allocate( q(KA,IA,JA), qexact(KA,IA,JA), u(KA,IA,JA) )

    call output_fvm_setup('linedom1d')
    call FILE_HISTORY_reg( "q", "q", "1", HST_ID(1), dim_type='X')
    call FILE_HISTORY_reg( "qexact", "qexact", "1", HST_ID(2), dim_type='X')

    !-- setup a module for evaluating numerical errors 
    call advect1d_fvm_numerror_Init( FZ, KS, KE, KA )

    !-- report information of time intervals
    call TIME_manager_report_timeintervals

    call PROF_rapend( "init", 1 )    
    return
  end subroutine init

  subroutine final()
    use mod_output_fvm, only: output_fvm_finalize
    use scale_time_manager, only: TIME_manager_Final 
    use mod_advect1d_fvm_numerror, only: advect1d_fvm_numerror_Final   
    implicit none
    !-----------------------------------------
    call PROF_rapstart( "final", 1 )

    call advect1d_fvm_numerror_Final()
    
    call optr_fvm%Final()
    call output_fvm_finalize
    call tinteg%Final()
    call TIME_manager_Final
    
    call PROF_rapend( "final", 1 )  
    call SCALE_finalize()

    return
  end subroutine final

end program test_advect1d_fvm
