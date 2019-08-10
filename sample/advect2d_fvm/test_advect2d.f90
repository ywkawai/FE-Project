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
  use scale_io, only: &
     H_SHORT
  use scale_prc_cartesC, only: &
     PRC_CARTESC_setup     
  use scale_comm_cartesC, only: &
     COMM_setup, &
     COMM_vars8, &
     COMM_wait
  use scale_atmos_grid_cartesC, only: &
     CX              => ATMOS_GRID_CARTESC_CX,              &
     CZ              => ATMOS_GRID_CARTESC_CZ,              &
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
    TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP,            &
    TIME_DTSEC, TIME_NSTEP 

  !-----------------------------------------------------------------------------
  implicit none

  integer, parameter :: NeGX = 128*2
  integer, parameter :: NeGY = 128*2
  integer, parameter :: NLocalMeshPerPrc = 1

  ! sin, cosbell, hat
  character(len=H_SHORT) :: InitShapeName = 'cosbell'

  real(RP), parameter :: dom_xmin = 0.0_RP
  real(RP), parameter :: dom_xmax = +2.0_RP
  real(RP), parameter :: dom_ymin = 0.0_RP
  real(RP), parameter :: dom_ymax = +2.0_RP

  real(RP), parameter :: ADV_VELX  = sqrt(1.0_RP)
  real(RP), parameter :: ADV_VELY  = sqrt(1.0_RP)
  
  integer, parameter :: FLUX_CD2_SCHEME = 1
  integer, parameter :: FLUX_CD4_SCHEME = 2
  integer :: FLUX_SCHEME_TYPE = FLUX_CD4_SCHEME

  integer :: nowstep
  integer :: rkstage
  integer, parameter :: nrkstage = 3
  real(RP) :: rkcoef1(nrkstage) = (/ 0.0_RP, 3.0_RP/4.0_RP, 1.0_RP/3.0_RP /)
  real(RP) :: rkcoef2(nrkstage) = (/ 1.0_RP, 1.0_RP/4.0_RP, 2.0_RP/3.0_RP /)

  real(RP) :: tsec_
  character(len=H_MID), parameter :: APPNAME = "advect2d with FVM"

  integer :: HST_ID(2)
  real(RP), allocatable :: q(:,:,:)
  real(RP), allocatable :: q0(:,:,:)
  real(RP), allocatable :: qexact(:,:,:)
  real(RP), allocatable :: u(:,:,:)
  real(RP), allocatable :: v(:,:,:)

  integer :: j

  !-------------------------------------------------------

  call init()
  call PROF_rapstart( 'set_initcond', 1)  
  call set_initcond()
  call PROF_rapend( 'set_initcond', 1)  

  call PROF_rapstart( 'TimeLoop', 1)
  do nowstep=1, TIME_NSTEP

    q0(:,:,JS) = q(:,:,JS)
    do rkstage=1, nrkstage

      !* Exchange halo data
      call PROF_rapstart( 'exchange_halo', 1)
      call excahge_halo( q(:,:,JS) )
      call PROF_rapend( 'exchange_halo', 1)

      !* Update prognostic variables
      call PROF_rapstart( 'update_dyn', 1)
      call update_dyn( &
          TIME_DTSEC, rkcoef1(rkstage), rkcoef2(rkstage), &
          q(:,:,JS), q0(:,:,JS), u(:,:,JS), v(:,:,JS) )
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
    call FILE_HISTORY_put(HST_ID(1), q(KS:KE,IS:IE,JS))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS:KE,IS:IE,JS))
    call FILE_HISTORY_write()
  end do
  call PROF_rapend( 'TimeLoop', 1)    

  call final()

contains
  subroutine update_dyn( dt, rkcoef_1, rkcoef_2, q_, q0_, u_, v_ )

    implicit none

    real(RP), intent(in) :: dt
    real(RP), intent(in) :: rkcoef_1, rkcoef_2
    real(RP), intent(out) :: q_(KA,IA)
    real(RP), intent(in)  :: q0_(KA,IA)
    real(RP), intent(in)  :: u_(KA,IA)
    real(RP), intent(in)  :: v_(KA,IA)

    integer :: k, i
    real(RP) :: qfluxUY(KA,IA)
    real(RP) :: qfluxXV(KA,IA)
    real(RP) :: dqdt
    
    real(RP) :: F2  = 0.5_RP
    real(RP) :: F41 =   7.0_RP/12.0_RP
    real(RP) :: F42 = - 1.0_RP/12.0_RP

    !------------------------------------------------------------------------

    call PROF_rapstart( 'update_dyn_cal_flux', 2)

    select case(FLUX_SCHEME_TYPE)
    case (FLUX_CD2_SCHEME) !----------------------------------------
      do i=IS, IE
      do k=KS-1, KE
        qfluxUY(k,i) = u_(k,i)*( F2*(q_(k,i) + q_(k+1,i)) )
      end do
      end do
      do i=IS-1, IE
      do k=KS, KE
        qfluxXV(k,i) = v_(k,i)*( F2*(q_(k,i) + q_(k,i+1)) )
      end do
      end do
    case (FLUX_CD4_SCHEME) !------------------------------------------------------------
      do i=IS, IE
      do k=KS-1, KE
        qfluxUY(k,i) = u_(k,i)*( F41*(q_(k,i) + q_(k+1,i)) + F42*(q_(k+2,i) + q_(k-1,i)) )
      end do
      end do
      do i=IS-1, IE
      do k=KS, KE
        qfluxXV(k,i) = v_(k,i)*( F41*(q_(k,i) + q_(k,i+1)) + F42*(q_(k,i+2) + q_(k,i-1)) )
      end do
      end do      
    end select

    call PROF_rapend( 'update_dyn_cal_flux', 2)

    !----

    call PROF_rapstart( 'update_dyn_cal_tend', 2)
    do i=IS, IE
    do k=KS, KE
      dqdt = - (qfluxUY(k,i) - qfluxUY(k-1,i  )) * RCDZ(k) &
             - (qfluxXV(k,i) - qfluxXV(k  ,i-1)) * RCDX(i)
      q_(k,i) = rkcoef_1*q0_(k,i) + rkcoef_2*(q_(k,i) + dt * dqdt)
    end do
    end do
    call PROF_rapend( 'update_dyn_cal_tend', 2)

    return
  end subroutine update_dyn
 
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

    real(RP) :: l2error
    real(RP) :: diff(KA)
    real(RP) :: exact_sol(KA)
    real(RP) :: x_uwind(KA), y_vwind(IA)
    integer :: periodx, periody

    real(RP) :: center_x, center_y
    integer :: k, i

    l2error = 0.0_RP
    center_x = 0.5_RP*(dom_xmax + dom_xmin)
    center_y = 0.5_RP*(dom_ymax + dom_ymin)
    periodx = ADV_VELX*tsec/(dom_xmax - dom_xmin)
    periody = ADV_VELY*tsec/(dom_ymax - dom_ymin)


    x_uwind(KS:KE) =  CZ(KS:KE) - center_x &
                    - (ADV_VELX*tsec - dble(periodx)*(dom_xmax - dom_xmin))
    where (x_uwind < -1.0_RP)
      x_uwind = dom_xmax + (x_uwind - dom_xmin)
    end where
    
    y_vwind(IS:IE) =  CX(IS:IE) - center_y &
                   - (ADV_VELY*tsec - dble(periody)*(dom_ymax - dom_ymin))
    where (y_vwind < -1.0_RP)
      y_vwind = dom_ymax + (y_vwind - dom_ymin)
    end where
    
    qexact(KS:KE,IS:IE,JS) = get_profile(InitShapeName, x_uwind(KS:KE), y_vwind(IS:IE))
    
    l2error = 0.0_RP
    do i=IS, IE
    do k=KS, KE
      l2error = l2error &
               + CDZ(k)*CDX(i) * ( q(k,i,JS) - qexact(k,i,JS) )**2
    end do
    end do
    write(*,*) "L2 error:", sqrt(l2error)/((dom_xmax - dom_xmin)*(dom_ymax - dom_ymin))

  end subroutine evaluate_error

  subroutine set_initcond() 
    implicit none

    integer :: i, j

    !-----------------------------------------
    do j=JS, JE
      q(:,:,j)      = get_profile( initShapeName,                        &
                                   CZ(:) - 0.5_RP*(dom_xmin + dom_xmax), &
                                   CX(:) - 0.5_RP*(dom_ymin + dom_ymax)   )
      qexact(:,:,j) = q(:,:,j)
    end do
    u(:,:,:)      = ADV_VELX
    v(:,:,:)      = ADV_VELY

    call FILE_HISTORY_put(HST_ID(1), q(KS:KE,IS:IE,JS))
    call FILE_HISTORY_put(HST_ID(2), qexact(KS:KE,IS:IE,JS))
    call FILE_HISTORY_write()   
    
    return
  end subroutine set_initcond

  function get_profile(profile_name, x, y) result(profile)
    use scale_const, only: &
      PI => CONST_PI
    implicit none

    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: x(:)
    real(RP), intent(in) :: y(:)
    real(RP) :: profile(size(x),size(y))

    real(RP) :: half_width = 0.3_RP
    real(RP) :: dist(size(x),size(y))

    integer :: j
    !------------------------------------------------------------------------

    profile(:,:) = 0.0_RP

    do j=1,size(y)
      dist(:,j) = sqrt(x(:)**2 + y(j)**2)
    end do

    select case(InitShapeName)
    case ('sin')
      do j=1, size(y)
        profile(:,j) = sin( PI*x(:) )
      end do
    case ('cosbell')
      where( dist <= half_width )
        profile(:,:) = (1.0_RP + cos(PI*dist(:,:)/half_width))*0.5_RP
      end where
    case ('hat')
      where( dist <= half_width )
        profile(:,:) = 1.0_RP
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
    real(RP) :: delx, dely
    !----------------------------------------------
    
    ! scale setup
    call SCALE_init( APPNAME ) 
    call PROF_rapstart( "init", 1 )

    ! setup process
    call PRC_CARTESC_setup

    call ATMOS_GRID_CARTESC_INDEX_setup( KMAX=NeGX, IMAX=NeGY, JMAX=1, JHALO=1, IBLOCK=1, JBLOCK=1 )


    ! setup horizontal/veritical grid system
    call ATMOS_GRID_CARTESC_allocate
    delx = (dom_xmax - dom_xmin)/real(NeGX,kind=RP)
    dely = (dom_ymax - dom_ymin)/real(NeGY,kind=RP)
    call ATMOS_GRID_CARTESC_generate( DZ=delx, DX=dely, DY=delx )

    ! setup file I/O
    call FILE_CARTESC_setup

    ! setup mpi communication
    call COMM_setup

    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init
    
    ! setup variables and history files
    allocate( q(KA,IA,JA), q0(KA,IA,JA), qexact(KA,IA,JA) )
    allocate( u(KA,IA,JA), v(KA,IA,JA) )

    call output_setup
    call FILE_HISTORY_reg( "q", "q", "1", HST_ID(1), dim_type='XY')
    call FILE_HISTORY_reg( "qexact", "qexact", "1", HST_ID(2), dim_type='XY')

    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()

    use mod_output, only: output_finalize
    use scale_time_manager, only: TIME_manager_Final    
    implicit none
    !----------------------------------------------

    call PROF_rapstart( "final", 1 )

    call output_finalize
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call SCALE_finalize()

    return
  end subroutine final

end program test_advect2d
