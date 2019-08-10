#include "scalelib.h"
program test_advect1d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer
  use scale_io, only: &
     H_SHORT
  use scale_prof
  use scale_prc_cartesC, only: &
     PRC_CARTESC_setup     
  use scale_comm_cartesC, only: &
     COMM_setup, &
     COMM_vars8, &
     COMM_wait
  use scale_atmos_grid_cartesC, only: &
     DOMAIN_CENTER_Y => ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
     CY              => ATMOS_GRID_CARTESC_CY,              &
     CZ              => ATMOS_GRID_CARTESC_CZ,              &
     FZ              => ATMOS_GRID_CARTESC_FZ,              &
     CDZ             => ATMOS_GRID_CARTESC_CDZ,             &
     RCDZ            => ATMOS_GRID_CARTESC_RCDZ,            &
     RFDZ            => ATMOS_GRID_CARTESC_RFDZ
  use scale_atmos_hydrometeor, only: &
     I_QV

  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate, FILE_HISTORY_put, FILE_HISTORY_write         

  use scale_time_manager, only: &
    TIME_manager_advance,                              &
    TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP,            &
    TIME_DTSEC, TIME_NSTEP 

  !-----------------------------------------------------------------------------
  implicit none

  integer, parameter :: NeGX = 512
  integer, parameter :: NLocalMeshPerPrc = 1

  ! sin, cosbell, hat
  character(len=H_SHORT) :: InitShapeName = 'cosbell'

  real(RP), parameter :: dom_xmin = 0.0_RP
  real(RP), parameter :: dom_xmax = +2.0_RP
  real(RP), parameter :: ADV_VEL  = 1.0_RP

  integer :: nowstep
  integer :: rkstage
  integer, parameter :: nrkstage = 3
  real(RP) :: rkcoef1(nrkstage) = (/ 0.0_RP, 3.0_RP/4.0_RP, 1.0_RP/3.0_RP /)
  real(RP) :: rkcoef2(nrkstage) = (/ 1.0_RP, 1.0_RP/4.0_RP, 2.0_RP/3.0_RP /)

  real(RP) :: tsec_
  character(len=H_MID), parameter :: APPNAME = "advect1d with FVM"

  integer :: HST_ID(2)
  real(RP), allocatable :: q(:,:,:)
  real(RP), allocatable :: q_A(:,:,:)
  real(RP), allocatable :: q0(:,:,:)
  real(RP), allocatable :: qexact(:,:,:)
  real(RP), allocatable :: u(:,:,:)

  !-------------------------------------------------------

  call init()
  call set_initcond()

  do nowstep=1, TIME_NSTEP
  
    q0(:,:,:) = q(:,:,:)

    do rkstage=1, nrkstage
      !* Exchange halo data
      call PROF_rapstart( 'exchange_halo', 1)
      call excahge_halo( q(:,IS,JS) )
      call PROF_rapend( 'exchange_halo', 1)

      !* Update prognostic variables
      call PROF_rapstart( 'update_dyn', 1)
      call update_dyn( &
          TIME_DTSEC, rkcoef1(rkstage), rkcoef2(rkstage), &
          q(:,IS,JS), q0(:,IS,JS), u(:,IS,JS) )
      call PROF_rapend( 'update_dyn', 1)
    end do
    
    !* Advance time
    call TIME_manager_advance()

    tsec_ = TIME_NOWDATE(6) + TIME_NOWMS
    if (mod(tsec_,0.25_RP) == 0.0_RP) then 
      write(*,*) "t=", real(tsec_), "[s]"
      call evaluate_error(tsec_)
    end if
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP )

    !* Output
    call FILE_HISTORY_put(HST_ID(1), q(KS:KE,IS,JS))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS:KE,IS,JS))
    call FILE_HISTORY_write()
  end do

  call final()

contains
  subroutine update_dyn( dt, rkcoef_1, rkcoef_2, q_, q0_, u_ )

    implicit none

    real(RP), intent(in) :: dt
    real(RP), intent(in) :: rkcoef_1, rkcoef_2
    real(RP), intent(out) :: q_(KA)
    real(RP), intent(in)    :: q0_(KA)
    real(RP), intent(in)    :: u_(KA)

    integer :: k
    real(RP) :: qflux(KA)
    real(RP) :: dqdt
    
    real(RP) :: F2  = 0.5_RP
    real(RP) :: F41 =   7.0_RP/12.0_RP
    real(RP) :: F42 = - 1.0_RP/12.0_RP

    !------------------------------------------------------------------------
    call PROF_rapstart( 'update_dyn_cal_flux', 2)
    do k=KS-1, KE
      !qflux(k) = u_(k)*( F2*(q_(k) + q_(k+1)) )
      qflux(k) = u_(k)*( F41*(q_(k) + q_(k+1)) + F42*(q_(k+2) + q_(k-1)) )
    end do
    call PROF_rapend( 'update_dyn_cal_flux', 2)

    call PROF_rapstart( 'update_dyn_cal_tend', 2)
    do k=KS, KE
      dqdt = - (qflux(k) - qflux(k-1)) * RCDZ(k)
      q_(k) = rkcoef_1*q0_(k) + rkcoef_2*(q_(k) + dt * dqdt)
    end do
    call PROF_rapend( 'update_dyn_cal_tend', 2)

    return
  end subroutine update_dyn
 
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

    real(RP) :: l2error
    real(RP) :: diff(KA)
    real(RP) :: exact_sol(KA)
    real(RP) :: x_uwind(KA)
    integer :: period

    real(RP) :: center_x

    l2error = 0.0_RP
    center_x = 0.5_RP*(dom_xmax + dom_xmin)
    period = ADV_VEL*tsec/(dom_xmax - dom_xmin)


    x_uwind(KS:KE) =  CZ(KS:KE) - center_x &
                    - (ADV_VEL*tsec - dble(period)*(dom_xmax - dom_xmin))
    where (x_uwind < -1.0_RP)
      x_uwind = dom_xmax + (x_uwind - dom_xmin)
    end where
    qexact(KS:KE,IS,JS) = get_profile(InitShapeName, x_uwind(KS:KE))
    
    l2error = l2error &
          + sum(  CDZ(KS:KE)                                     &
                * ( q(KS:KE,IS,JS) - qexact(KS:KE,IS,JS) )**2 ) 

    write(*,*) "L2 error:", sqrt(l2error)/(dom_xmax - dom_xmin)
  end subroutine evaluate_error

  subroutine set_initcond() 
    implicit none

    integer :: i, j

    !-----------------------------------------
    do j=JS, JE
    do i=IS, IE
      q(:,i,j)      = get_profile(initShapeName, CZ(:) - 0.5_RP*(dom_xmin + dom_xmax))
      qexact(:,i,j) = q(:,i,j)
    end do
    end do
    u(:,:,:)      = ADV_VEL

    call FILE_HISTORY_put(HST_ID(1), q(KS:KE,IS,JS))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS:KE,IS,JS))
    call FILE_HISTORY_write()   
    
    return
  end subroutine set_initcond

  function get_profile(profile_name, x) result(profile)
    use scale_const, only: PI => CONST_PI 
    implicit none

    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: x(:)
    real(RP) :: profile(size(x))

    real(RP) :: half_width = 0.3_RP

    profile(:) = 0.0_RP

    select case(InitShapeName)
    case ('sin')
      profile(:) = sin( PI*x(:) )
    case ('cosbell')
      where( abs(x) <= half_width )
        profile(:) = (1.0_RP + cos(PI*x/half_width))*0.5_RP
      end where
    case ('hat')
      where( abs(x) <= half_width )
        profile(:) = 1.0_RP
      end where
    end select

    return
  end function get_profile

  subroutine init()

    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only: TIME_manager_Init 
    use scale_atmos_grid_cartesC, only: &
      ATMOS_GRID_CARTESC_allocate, &
      ATMOS_GRID_CARTESC_generate
    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_regist
    use scale_file_cartesC, only: &
      FILE_CARTESC_setup
    use mod_output, only: output_setup
    use scale_file_history, only: FILE_HISTORY_reg

    implicit none

    integer :: comm, myrank, nprocs
    logical :: ismaster
    real(RP) :: del

    !----------------------------------------------
    
    ! scale setup
    call SCALE_init( APPNAME )  
    
    ! setup profiler
    call PROF_setup
    call PROF_rapstart( "total", 0 )
    call PROF_rapstart( "init", 1 )

    ! setup process
    call PRC_CARTESC_setup

    call ATMOS_GRID_CARTESC_INDEX_setup( KMAX=NeGX, IMAX=2, JMAX=2, IBLOCK=1, JBLOCK=1 )


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
    
    ! setup variables and history files
    allocate( q(KA,IA,JA), q_A(KA,IA,JA), q0(KA,IA,JA), qexact(KA,IA,JA), u(KA,IA,JA) )

    call output_setup
    call FILE_HISTORY_reg( "q", "q", "1", HST_ID(1), dim_type='X')
    call FILE_HISTORY_reg( "qexact", "qexact", "1", HST_ID(2), dim_type='X')

    call PROF_rapend( "init", 1 )    
    return
  end subroutine init

  subroutine final()
    use mod_output, only: output_finalize
    use scale_time_manager, only: TIME_manager_Final    
    implicit none

    call PROF_rapstart( "final", 1 )

    call output_finalize
    call TIME_manager_Final
    
    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport    
    call SCALE_finalize()

    return
  end subroutine final

end program test_advect1d
