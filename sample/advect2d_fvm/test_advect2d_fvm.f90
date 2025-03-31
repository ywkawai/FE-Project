!-------------------------------------------------------------------------------
!> Program A sample program: 2-dimensional linear advection test using FVM
!! 
!! 
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
program test_advect2d_fvm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_prof
  use scale_io

  use scale_atmos_grid_cartesC, only: &
     CX              => ATMOS_GRID_CARTESC_CX,              &
     FX              => ATMOS_GRID_CARTESC_FX,              &
     CY              => ATMOS_GRID_CARTESC_CY,              &
     FY              => ATMOS_GRID_CARTESC_FY,              &
     CDX             => ATMOS_GRID_CARTESC_CDX,             &
     CDY             => ATMOS_GRID_CARTESC_CDY,             &
     RCDX            => ATMOS_GRID_CARTESC_RCDX,            &
     RCDY            => ATMOS_GRID_CARTESC_RCDY

  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate, &
    FILE_HISTORY_put, FILE_HISTORY_write         

  use scale_time_manager, only: &
    TIME_manager_checkstate, TIME_manager_advance,     &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,        &
    TIME_DTSEC, TIME_NSTEP, TIME_DOresume, TIME_DOend

  use scale_timeint_rk, only: timeint_rk 
  use mod_operator_fvm, only: operator_fvm    
  
  use mod_advect2d_fvm_numerror, only: advect2d_fvm_numerror_eval
  !-----------------------------------------------------------------------------
  implicit none

  integer :: NeGX, GXHALO
  integer :: NeGY, GYHALO

  character(len=H_SHORT) :: InitShapeName !< The type of initial q (sin, gaussian-hill, cosine-bell, top-hat)
  real(RP), save :: InitShapeParams(4)
  character(len=H_SHORT) :: VelTypeName   !< The type of specified velocify field (constant, rigid-body-rot)
  real(RP), save :: VelTypeParams(4)
  logical :: Do_NumErrorAnalysis          !< Flag wheter analysis of numerical error is performed

  type(operator_fvm) :: optr_fvm

  real(RP), allocatable :: q(:,:,:)
  real(RP), allocatable :: qexact(:,:,:)
  real(RP), allocatable :: u(:,:,:)
  real(RP), allocatable :: v(:,:,:)
  integer, save :: HST_ID(2)

  real(RP) :: tsec_
  type(timeint_rk) :: tinteg 
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1

  character(len=H_MID), parameter :: APPNAME = "advect2d with FVM"
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
      call excahge_halo( q(KS,:,:) )
      call PROF_rapend( 'exchange_halo', 1)

      call PROF_rapstart( 'set_velocity', 1)
      call set_velocity( u, v, tsec_ )
      call PROF_rapend( 'set_velocity', 1)  

      !* Update prognostic variables

      call PROF_rapstart( 'cal_tend', 1)
      tintbuf_ind = tinteg%tend_buf_indmap(rkstage)
      call cal_tend( tinteg%tend_buf2D_ex(:,:,RKVAR_Q,tintbuf_ind), q, u, v )
      call PROF_rapend( 'cal_tend', 1) 

      call PROF_rapstart( 'update_var', 1)
      call tinteg%Advance(rkstage, q(KS,:,:), RKVAR_Q, IS, IE, JS, JE)
      call PROF_rapend('update_var', 1)
    end do

    tsec_ = TIME_DTSEC * real(TIME_NOWSTEP-1, kind=RP)
    if ( Do_NumErrorAnalysis ) then
      call advect2d_fvm_numerror_eval( qexact(KS,:,:),                                              & ! (out)
        q(KS,:,:), TIME_NOWSTEP, tsec_, VelTypeName, VelTypeParams, InitShapeName, InitShapeParams, & ! (in)
        CX, FX, CY, FY, IS, IE, IA, IHALO, JS, JE, JA, IHALO                                        ) ! (in)
    end if

    !* Output history file

    call FILE_HISTORY_put(HST_ID(1), q(KS,IS:IE,JS:JE))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS,IS:IE,JS:JE))
    call FILE_HISTORY_write()

    if (TIME_DOend) exit
  end do

  call final()

contains
  !> Calculate the tendency
  !! dqdt = - ( <u q>_i - <u q>_i-1 ) / DX - ( <v q>_j - <v q>_j-1 ) / DY
  !!
  subroutine cal_tend( dqdt, q_, u_, v_ )
    implicit none
    real(RP), intent(out) :: dqdt(KS:KS,IA,JA)
    real(RP), intent(in) :: q_(KA,IA,JA)
    real(RP), intent(in)  :: u_(KA,IA,JA)
    real(RP), intent(in)  :: v_(KA,IA,JA)

    integer :: i, j
    real(RP) :: qflux(KA,IA,JA,2)
    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_tend_flux', 2)
    call optr_fvm%C_flux_UYZ( qflux(:,:,:,1), u_, q_ )
    call optr_fvm%C_flux_XVZ( qflux(:,:,:,2), v_, q_ )
    call PROF_rapend( 'cal_tend_flux', 2)

    !----

    call PROF_rapstart( 'cal_tend_dqdt', 2)
    do j=JS, JE
    do i=IS, IE
      dqdt(KS,i,j) = - (qflux(KS,i,j,1) - qflux(KS,i-1,j,1)) * RCDX(i) &
                     - (qflux(KS,i,j,2) - qflux(KS,i,j-1,2)) * RCDY(j)
    end do
    end do
    call PROF_rapend( 'cal_tend_dqdt', 2)

    return
  end subroutine cal_tend
 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine excahge_halo(var)
    real(RP), intent(inout) :: var(IA,JA)
    integer :: i, j

    do j=IS, IE
    do i=1, IHALO
      var(IS-i,j) = var(IE-i+1,j)
      var(IE+i,j) = var(IS+i-1,j)
    end do
    end do

    do j=1, JHALO
    do i=IS, IE
      var(i,JS-j) = var(i,JE-j+1)
      var(i,JE+j) = var(i,JS+j-1)
    end do
    end do
  end subroutine excahge_halo

  !> Set velocity data at tsec
  subroutine set_velocity( u_, v_, tsec )
    use mod_fieldutil, only: fieldutil_get_profile2d_flow
    implicit none
    real(RP), intent(inout) :: u_(KA,IA,JA)
    real(RP), intent(inout) :: v_(KA,IA,JA)
    real(RP), intent(in) :: tsec

    real(DP) :: y_tmp(IA), v_uyz(IA), u_xvz(IA)
    integer :: j
    !-----------------------------------------------

    VelTypeParams(4) = tsec

    do j=JS, JE
      y_tmp(:) = CY(j)      
      call fieldutil_get_profile2d_flow( u_(KS,IS-1:IE,j), v_uyz(IS-1:IE), & ! (out)
        VelTypeName, FX(IS-1:IE), y_tmp(IS-1:IE), VelTypeParams, IE-IS+2   ) ! (in)
    end do
    do j=JS-1, JE
      y_tmp(:) = FY(j)      
      call fieldutil_get_profile2d_flow( u_xvz(IS:IE), v_(KS,IS:IE,j), & ! (out)
        VelTypeName, CX(IS:IE), y_tmp(IS:IE), VelTypeParams, IE-IS+1   ) ! (in)
    end do

    return
  end subroutine set_velocity

  !> Set inital data
  subroutine set_initcond() 
    use mod_fieldutil, only: fieldutil_get_profile2d_tracer 
    implicit none

    real(RP) :: x(KA,IA), y(KA,IA)
    real(DP) :: y_tmp(IA)
    integer :: j  

    !----------------------------------------- 
    do j=JS, JE
      y_tmp(:) = CY(j)
      call fieldutil_get_profile2d_tracer ( q(KS,IS:IE,j),               & ! (out)
        InitShapeName, CX(IS:IE), y_tmp(KS:KE), InitShapeParams, IE-IS+1 )  ! (in)
    end do
    call set_velocity( u, v, 0.0_RP )

    if ( Do_NumErrorAnalysis ) then
      call advect2d_fvm_numerror_eval( qexact(KS,:,:),                                    & ! (out)
        q(KS,:,:), 1, 0.0_RP, VelTypeName, VelTypeParams, InitShapeName, InitShapeParams, & ! (in)
        CX, FX, CY, FY, IS, IE, IA, IHALO, JS, JE, JA, JHALO                              ) ! (in)
    end if
    call FILE_HISTORY_put(HST_ID(1), q(KS,IS:IE,JS:JE))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS,IS:IE,JS:JE))
    call FILE_HISTORY_write()   
    
    LOG_PROGRESS('(A,F13.5,A)') "t=", real(0.0_RP), "[s]"
    return
  end subroutine set_initcond

  !> Initialization
  subroutine init()
    use scale_prc_cartesC, only: PRC_CARTESC_setup     
    use scale_comm_cartesC, only: &
      COMM_setup, COMM_regist
    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only:           &
      TIME_manager_Init,                    &
      TIME_manager_report_timeintervals
    use scale_atmos_grid_cartesC, only: &
      ATMOS_GRID_CARTESC_allocate, &
      ATMOS_GRID_CARTESC_generate
    use scale_file_cartesC, only: FILE_CARTESC_setup
    use mod_output_fvm, only: output_fvm_setup
    use scale_file_history, only: FILE_HISTORY_reg

    use mod_advect2d_fvm_numerror, only: advect2d_fvm_numerror_Init
    implicit none

    real(RP), parameter :: dom_xmin =  0.0_RP
    real(RP), parameter :: dom_xmax = +1.0_RP
    real(RP), parameter :: dom_ymin =  0.0_RP
    real(RP), parameter :: dom_ymax = +1.0_RP
    character(len=H_SHORT) :: FLUX_SCHEME_TYPE
    character(len=H_SHORT) :: TINTEG_SCHEME_TYPE

    namelist /PARAM_TEST/ &
      NeGX, GXHALO, NeGY, GYHALO,            &
      FLUX_SCHEME_TYPE, TINTEG_SCHEME_TYPE,  &
      InitShapeName, InitShapeParams,        &
      VelTypeName, VelTypeParams,            &
      Do_NumErrorAnalysis

    integer :: ierr      
    real(RP) :: delx, dely    
    integer :: gid
    !----------------------------------------------
    
    !-- setup SCALE modules

    call SCALE_init( APPNAME ) 
    call PROF_rapstart( "init", 1 )

    !-- read namelist

    NeGX = 2; GXHALO =2
    NeGY = 2; GYHALO =2
    FLUX_SCHEME_TYPE    = 'CD2'
    TINTEG_SCHEME_TYPE  = 'RK_SSP_3s3o'
    InitShapeName       = 'sin'
    InitShapeParams(:)  = (/ 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP /)
    VelTypeName         = 'const'
    VelTypeParams(:)    = (/ 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP /)
    Do_NumErrorAnalysis = .false.

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
      KMAX=2, KHALO=2, IMAX=NeGX, IHALO=GXHALO, JMAX=NeGY, JHALO=GYHALO, &
      IBLOCK=1, JBLOCK=1 )

    call ATMOS_GRID_CARTESC_allocate
    delx = (dom_xmax - dom_xmin)/real(NeGX,kind=RP)
    dely = (dom_ymax - dom_ymin)/real(NeGY,kind=RP)
    call ATMOS_GRID_CARTESC_generate( &
      DZ=delx, DX=delx, DY=dely )

    !-- setup file I/O
    call FILE_CARTESC_setup

    !-- setup mpi communication
    call COMM_setup
    call COMM_regist( KA, IA, JA, IHALO, JHALO, gid )

    !-- setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init
    
    !-- setup a module for FVM operator
    call optr_fvm%Init( FLUX_SCHEME_TYPE, KS, KE, KA, IS, IE, IA, JS, JE, JA )

    !-- setup a module for time integrator
    call tinteg%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1, &
                      2, (/ IA, JA /) )
    
    !-- setup variables and history files

    allocate( q(KA,IA,JA), qexact(KA,IA,JA) )
    allocate( u(KA,IA,JA), v(KA,IA,JA) )

    call output_fvm_setup( 'rectdom2d' )
    call FILE_HISTORY_reg( "q", "q", "1", HST_ID(1), dim_type='XY')
    call FILE_HISTORY_reg( "qexact", "qexact", "1", HST_ID(2), dim_type='XY')

    !-- setup a module for evaluating numerical errors 
    if ( Do_NumErrorAnalysis ) &
      call advect2d_fvm_numerror_Init( FX, FY, IS, IE, IA, JS, JE, JA )

    !-- report information of time intervals
    call TIME_manager_report_timeintervals

    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()
    use mod_output_fvm, only: output_fvm_finalize
    use scale_time_manager, only: TIME_manager_Final    
    use mod_advect2d_fvm_numerror, only: advect2d_fvm_numerror_Final   
    implicit none
    !----------------------------------------------

    call PROF_rapstart( "final", 1 )
    if ( Do_NumErrorAnalysis ) &
      call advect2d_fvm_numerror_Final()

    call optr_fvm%Final()
    call output_fvm_finalize
    call tinteg%Final()
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call SCALE_finalize()

    return
  end subroutine final

end program test_advect2d_fvm
