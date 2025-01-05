#include "scalelib.h"
program test_advect2d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer
  use scale_prof
  use scale_io
  use scale_atmos_grid_cartesC, only: &
     CX              => ATMOS_GRID_CARTESC_CX,              &
     FX              => ATMOS_GRID_CARTESC_FX,              &
     CZ              => ATMOS_GRID_CARTESC_CZ,              &
     FZ              => ATMOS_GRID_CARTESC_FZ,              &
     CDZ             => ATMOS_GRID_CARTESC_CDZ,             &
     CDX             => ATMOS_GRID_CARTESC_CDX,             &
     RCDZ            => ATMOS_GRID_CARTESC_RCDZ,            &
     RCDX            => ATMOS_GRID_CARTESC_RCDX,            &
     RFDZ            => ATMOS_GRID_CARTESC_RFDZ

  use scale_atmos_hydrometeor, only: &
     I_QV

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
    get_profile2d_tracer => fieldutil_get_profile2d_tracer, &
    get_profile2d_flow => fieldutil_get_profile2d_flow    

  use mod_operator_fvm, only: &
    operator_fvm

  
  !-----------------------------------------------------------------------------
  implicit none

  integer :: NeGX, GXHALO
  integer :: NeGY, GYHALO

  ! The type of initial q (sin, gaussian-hill, cosine-bell, top-hat)
  character(len=H_SHORT) :: InitShapeName
  real(RP), save :: InitShapeParams(4)
  ! The type of specified velocify field (constant, rigid-body-rot)
  character(len=H_SHORT) :: VelTypeName 
  real(RP), save :: VelTypeParams(4)

  real(RP), parameter :: dom_xmin =  0.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: dom_ymin =  0.0_RP
  real(RP), parameter :: dom_ymax = +1.0_RP

  character(len=H_SHORT) :: FLUX_SCHEME_TYPE
  type(operator_fvm) :: optr_fvm

  character(len=H_SHORT) :: TINTEG_SCHEME_TYPE
  type(timeint_rk) :: tinteg
  integer :: nowstep
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1

  real(RP) :: tsec_
  character(len=H_MID), parameter :: APPNAME = "advect2d with FVM"

  integer, save :: HST_ID(2)
  real(RP), allocatable :: q(:,:,:)
  real(RP), allocatable :: qexact(:,:,:)
  real(RP), allocatable :: u(:,:,:)
  real(RP), allocatable :: v(:,:,:)

  integer :: j
  integer :: nstep_eval_error
  !-------------------------------------------------------

  call init()
  call PROF_rapstart( 'set_initcond', 1)  
  call set_initcond()
  call PROF_rapend( 'set_initcond', 1)  

  call PROF_rapstart( 'TimeLoop', 1)
  do nowstep=1, TIME_NSTEP
    do rkstage=1, tinteg%nstage
      tsec_ = TIME_NOWDATE(6) + TIME_NOWSUBSEC

      !* Exchange halo data
      call PROF_rapstart( 'exchange_halo', 1)
      call excahge_halo( q(:,:,JS) )
      call PROF_rapend( 'exchange_halo', 1)

      call PROF_rapstart( 'set_velocity', 1)
      call set_velocity( u, v, tsec_ )
      call PROF_rapend( 'set_velocity', 1)  

      !* Update prognostic variables
      call PROF_rapstart( 'cal_dyn_tend', 1)
      tintbuf_ind = tinteg%tend_buf_indmap(rkstage)
      call cal_dyn_tend( tinteg%tend_buf2D_ex(:,:,RKVAR_Q,tintbuf_ind), q, u, v )
      call PROF_rapend( 'cal_dyn_tend', 1) 

      call PROF_rapstart( 'update_var', 1)
      call tinteg%Advance(rkstage, q(:,:,JS), RKVAR_Q, KS, KE, IS, IE)
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
    call FILE_HISTORY_put(HST_ID(1), q(KS:KE,IS:IE,JS))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS:KE,IS:IE,JS))
    call FILE_HISTORY_write()
  end do
  call PROF_rapend( 'TimeLoop', 1)    

  call final()

contains
  subroutine cal_dyn_tend( dqdt, q_, u_, v_ )

    implicit none

    real(RP), intent(out) :: dqdt(KA,IA,JS:JS)
    real(RP), intent(in) :: q_(KA,IA,JA)
    real(RP), intent(in)  :: u_(KA,IA,JA)
    real(RP), intent(in)  :: v_(KA,IA,JA)

    integer :: k, i, j
    real(RP) :: qflux(KA,IA,JA,2)

    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_flux', 2)
    call optr_fvm%C_flux_XYW( qflux(:,:,:,1), u_, q_ )
    call optr_fvm%C_flux_UYZ( qflux(:,:,:,2), v_, q_ )
    call PROF_rapend( 'cal_dyn_tend_flux', 2)

    !----

    call PROF_rapstart( 'cal_dyn_tend_dqdt', 2)
    do j=JS, JS
    do i=IS, IE
    do k=KS, KE
      dqdt(k,i,j) = - (qflux(k,i,j,1) - qflux(k-1,i,j,1)) * RCDZ(k) &
                    - (qflux(k,i,j,2) - qflux(k,i-1,j,2)) * RCDX(i)
    end do
    end do
    end do
    call PROF_rapend( 'cal_dyn_tend_dqdt', 2)

    return
  end subroutine cal_dyn_tend
 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine excahge_halo(var)
    real(RP), intent(inout) :: var(KA,IA)
    integer :: k, i 

    do i=IS, IE
    do k=1, KHALO
      var(KS-k,i) = var(KE-k+1,i)
      var(KE+k,i) = var(KS+k-1,i)
    end do
    end do

    do i=1, IHALO
    do k=KS, KE
      var(k,IS-i) = var(k,IE-i+1)
      var(k,IE+i) = var(k,IS+i-1)
    end do
    end do
  end subroutine excahge_halo

  subroutine evaluate_error(tsec)
    
    implicit none
    real(DP), intent(in) :: tsec

    real(RP) :: diff(KA)
    real(RP) :: exact_sol(KA)
    real(RP) :: x_uwind_1d(KA), y_vwind_1d(IA)
    real(RP) :: y_vwind(KA)
    integer :: k, i

    real(RP) :: l2error
    real(RP) :: linferror
    real(RP) :: ADV_VELX, ADV_VELY    
    !------------------------------------------------------------------------

    l2error   = 0.0_RP   
    linferror = 0.0_RP
    ADV_VELX = VelTypeParams(1); ADV_VELY = VelTypeParams(2)

    x_uwind_1d(KS:KE) =  get_upwind_pos1d(CZ(KS:KE), ADV_VELX, tsec, dom_xmin, dom_xmax)
    y_vwind_1d(IS:IE) =  get_upwind_pos1d(CX(IS:IE), ADV_VELY, tsec, dom_ymin, dom_ymax)
    do i=IS, IE
      y_vwind(:) = y_vwind_1d(i)
      call get_profile2d_tracer( qexact(KS:KE,i,j),                                 & ! (out)
        InitShapeName, x_uwind_1d(KS:KE), y_vwind(KS:KE), InitShapeParams, KE-KS+1 )  ! (in)      
    end do

    do i=IS, IE
    do k=KS, KE
      l2error = l2error &
               + CDZ(k)*CDX(i) * ( q(k,i,JS) - qexact(k,i,JS) )**2
      linferror = max(linferror, abs(q(k,i,JS) - qexact(k,i,JS)))
    end do
    end do
    LOG_INFO("evaluate_error_l2",*) sqrt(l2error)/((dom_xmax - dom_xmin)*(dom_ymax - dom_ymin))
    LOG_INFO("evaluate_error_linf",*) linferror

  end subroutine evaluate_error

  subroutine set_velocity( u_, v_, tsec )
    real(RP), intent(inout) :: u_(KA,IA,JA)
    real(RP), intent(inout) :: v_(KA,IA,JA)
    real(RP), intent(in) :: tsec

    real(DP) :: y_tmp(KA), v_uyz(KA), u_xvz(KA)
    integer :: k, i    
    !-----------------------------------------------

    VelTypeParams(4) = tsec

    do j=JS, JE
    do i=IS, IE
      y_tmp(:) = CX(i)      
      call get_profile2d_flow( u_(KS-1:KE,i,j), v_uyz(KS-1:KE),              & ! (out)
        VelTypeName, FZ(KS-1:KE), y_tmp(KS-1:KE), VelTypeParams, KE-KS+2 )    ! (in)
    end do
    end do
    do j=JS, JE
    do i=IS-1, IE
        y_tmp(:) = FX(i)      
        call get_profile2d_flow( u_xvz(KS:KE), v_(KS:KE,i,j),                & ! (out)
          VelTypeName, CZ(KS:KE), y_tmp(KS:KE), VelTypeParams, KE-KS+1 )      ! (in)
    end do
    end do

  end subroutine set_velocity

  subroutine set_initcond() 
    implicit none

    real(RP) :: x(KA,IA), y(KA,IA)
    real(DP) :: y_tmp(KA), v_uyz(KA), u_xvz(KA)
    integer :: k, i    

    !----------------------------------------- 
    do j=JS, JE
    do i=IS, IE
      y_tmp(:) = CX(i)
      call get_profile2d_tracer( q(KS:KE,i,j),                            & ! (out)
        InitShapeName, CZ(KS:KE), y_tmp(KS:KE), InitShapeParams, KE-KS+1 )  ! (in)

      qexact(:,i,j) = q(:,i,j)
    end do
    end do
    call set_velocity( u, v, 0.0_RP )

    call FILE_HISTORY_put(HST_ID(1), q(KS:KE,IS:IE,JS))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS:KE,IS:IE,JS))
    call FILE_HISTORY_write()   
    
    LOG_PROGRESS('(A,F13.5,A)') "t=", real(0.0_RP), "[s]"
    call evaluate_error(0.0_RP)
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
      NeGX, GXHALO, NeGY, GYHALO,            &
      FLUX_SCHEME_TYPE, TINTEG_SCHEME_TYPE,  &
      InitShapeName, InitShapeParams,        &
      VelTypeName, VelTypeParams,            &
      nstep_eval_error
        
    integer :: comm, myrank, nprocs
    logical :: ismaster
    real(RP) :: delx, dely
    
    integer :: ierr
    !----------------------------------------------
    
    ! scale setup
    call SCALE_init( APPNAME ) 
    call PROF_rapstart( "init", 1 )

    !--- read namelist

    NeGX = 2; GXHALO =2
    NeGY = 2; GYHALO =2
    FLUX_SCHEME_TYPE   = 'CD2'
    TINTEG_SCHEME_TYPE = 'RK_SSP_3s3o'
    InitShapeName      = 'sin'
    InitShapeParams(:) = (/ 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP /)
    VelTypeName        = 'const'
    VelTypeParams(:)   = (/ 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP /)
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
      KMAX=NeGX, KHALO=GXHALO, IMAX=NeGY, IHALO=GYHALO, JMAX=1, JHALO=1, &
      IBLOCK=1, JBLOCK=1 )

    ! setup horizontal/veritical grid system
    call ATMOS_GRID_CARTESC_allocate
    delx = (dom_xmax - dom_xmin)/real(NeGX,kind=RP)
    dely = (dom_ymax - dom_ymin)/real(NeGY,kind=RP)
    call ATMOS_GRID_CARTESC_generate( &
      DZ=delx, DX=dely, DY=delx )

    ! setup file I/O
    call FILE_CARTESC_setup

    ! setup mpi communication
    call COMM_setup

    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init
    
    ! setup a module for FVM operator
    call optr_fvm%Init( FLUX_SCHEME_TYPE, KS, KE, KA, IS, IE, IA, JS, JS, JA )

    ! setup a module for time integrator
    call tinteg%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1, &
                      2, (/ KA, IA /) )
    
    ! setup variables and history files
    allocate( q(KA,IA,JA), qexact(KA,IA,JA) )
    allocate( u(KA,IA,JA), v(KA,IA,JA) )

    call output_fvm_setup( 'rectdom2d' )
    call FILE_HISTORY_reg( "q", "q", "1", HST_ID(1), dim_type='XY')
    call FILE_HISTORY_reg( "qexact", "qexact", "1", HST_ID(2), dim_type='XY')

    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()

    use mod_output_fvm, only: output_fvm_finalize
    use scale_time_manager, only: TIME_manager_Final    
    implicit none
    !----------------------------------------------

    call PROF_rapstart( "final", 1 )

    call optr_fvm%Final()
    call tinteg%Final()

    call output_fvm_finalize
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call SCALE_finalize()

    return
  end subroutine final

end program test_advect2d
