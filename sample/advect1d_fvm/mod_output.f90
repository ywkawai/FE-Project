!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_output
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_atmos_grid_cartesC_index
  use scale_file_history, only: &
    FILE_HISTORY_Setup, &
    FILE_HISTORY_Set_NowDate, &
    FILE_HISTORY_truncate_1D, &
    FILE_HISTORY_truncate_2D, &
    FILE_HISTORY_truncate_3D, &
    FILE_HISTORY_put,         &
    FILE_HISTORY_write,       &
    FILE_HISTORY_Set_Dim,     &
    FILE_HISTORY_Set_Axis,    &
    FILE_HISTORY_finalize

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: output_setup
  public :: output_finalize
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  integer,          parameter :: nzs = 1
  character(len=8), parameter :: zs(nzs) = (/  "model   " /)

  integer  :: FILE_HISTORY_MESHFIELD_STARTDATE(6) !< start time [YYYY MM DD HH MM SS]
  real(DP) :: FILE_HISTORY_MESHFIELD_STARTMS      !< subsecond part of start time [millisec]

contains

!----------------

subroutine output_setup()

  use scale_file_h, only: &
    FILE_HSHORT
  use scale_prc, only: &
    PRC_masterrank, &
    PRC_myrank
  use scale_time, only: &
    TIME_NOWDATE,       &
    TIME_NOWMS,         &
    TIME_STARTDAYSEC,   &
    TIME_DTSEC,         &
    TIME_NOWSTEP
  use scale_calendar, only: &
    CALENDAR_get_name
  implicit none

  character(len=H_MID) :: FILE_HISTORY_MESHFILED_H_TITLE = 'SCALE-FEM FILE_HISTORY_MESHFIELD' !< title of the output file
  character(len=H_MID) :: FILE_HISTORY_MESHFIELD_T_SINCE

  
  character(len=FILE_HSHORT) :: calendar
  real(DP) :: start_daysec
  integer  :: ierr
  integer  :: k

  !---------------------------------------------------------------------------


  FILE_HISTORY_MESHFIELD_STARTDATE(:) = TIME_NOWDATE
  FILE_HISTORY_MESHFIELD_STARTMS      = TIME_NOWMS

  start_daysec = TIME_STARTDAYSEC
  if ( TIME_NOWDATE(1) > 0 ) then
     write(FILE_HISTORY_MESHFIELD_T_SINCE,'(I4.4,5(A1,I2.2))') TIME_NOWDATE(1), &
                                                        '-', TIME_NOWDATE(2), &
                                                        '-', TIME_NOWDATE(3), &
                                                        ' ', TIME_NOWDATE(4), &
                                                        ':', TIME_NOWDATE(5), &
                                                        ':', TIME_NOWDATE(6)
     start_daysec = TIME_NOWMS
  else
     FILE_HISTORY_MESHFIELD_T_SINCE = ''
  endif

  ! get calendar name
  call CALENDAR_get_name( calendar )

  call FILE_HISTORY_Setup( FILE_HISTORY_MESHFILED_H_TITLE,        & ! [IN]
    H_SOURCE, H_INSTITUTE,                                        & ! [IN]
    start_daysec, TIME_DTSEC,                                     & ! [IN]
    time_since = FILE_HISTORY_MESHFIELD_T_SINCE,                  & ! [IN]
    calendar = calendar,                                          & ! [IN]
    default_zcoord = 'model',                                     & ! [IN]
    myrank = PRC_myrank                        )                    ! [IN]

  call FILE_HISTORY_Set_NowDate( TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP )

  call set_dims()
  call set_axis()

end subroutine output_setup
subroutine output_finalize
  call FILE_HISTORY_finalize
end subroutine output_finalize
!----------------

subroutine set_dims()
  character(len=H_SHORT) :: dims(3,3)
  integer :: start(3,3), count(3,3)  

  start(1,1) = 1
  dims(1,1)  = "x"
  count(1,1) = KMAX
  call FILE_HISTORY_Set_Dim( "X", 1, 1, dims(:,:), zs(:), start(:,:), count(:,:))
end subroutine set_dims

subroutine set_axis()
  use scale_atmos_grid_cartesC, only: &
     DOMAIN_CENTER_Y => ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
     CY              => ATMOS_GRID_CARTESC_CY,              &
     CZ              => ATMOS_GRID_CARTESC_CZ,              &
     FZ              => ATMOS_GRID_CARTESC_FZ,              &
     CDZ             => ATMOS_GRID_CARTESC_CDZ,             &
     RCDZ            => ATMOS_GRID_CARTESC_RCDZ,            &
     RFDZ            => ATMOS_GRID_CARTESC_RFDZ
    

  integer :: k
  integer :: startz
  real(RP) :: z_bnds(2,KA)

  !--------------------
    !
    startz = 1

    ! bounds
    do k = KS, KE
      z_bnds(1,k) = FZ(k-1)
      z_bnds(2,k) = FZ(k  )
    end do
  call FILE_HISTORY_Set_Axis( 'x', 'X-coordinate', '1', 'x', CZ(KS:KE))!,               &
                              !bounds=z_bnds (:,KS  :KE), gsize=KMAX, start=startZ )
end subroutine set_axis

!----------------

end module mod_output
