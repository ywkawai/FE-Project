!-------------------------------------------------------------------------------
!> Program regrid tool (SCALE-DG)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          regrid tool for models based DG methods
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
program regrid_tool
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort

  use mod_regrid_mesh, only: &
    out_mesh, nodemap
  use mod_regrid_interp_field, only: &
    OutVarInfo,                      &
    out_vinfo, out_var_num,          &
    out_var3D, out_var2D,            &
    regrid_interp_field_Interpolate
  use mod_regrid_interp_vcoord, only: &
    regrid_interp_vcoord,             &
    REGRID_VCOORD_MODEL_ID
  use mod_regrid_file, only: &
    regrid_file_write_var

  
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

  ! namelist parameters
  logical                 :: debug            = .false.

  ! MPI parameters
  integer                 :: nprocs                      ! number of processes               (execution)
  integer                 :: myrank                      ! my rank                           (execution)
  logical                 :: ismaster                    ! master process?                   (execution)

  logical :: do_output
  integer :: vid
  integer :: istep
  real(DP) :: start_sec

  type(regrid_interp_vcoord) :: vintrp
  type(OutVarInfo), pointer :: vinfo

  !-----------------------------------------------------------------------------

  call initialize()
   
  !-
  !########## main ##########
  call PROF_rapstart('Main', 0)

  LOG_NEWLINE
  LOG_PROGRESS(*) 'START LOOP'

  do vid =1, out_var_num
    vinfo => out_vinfo(vid)
    do istep=1, vinfo%num_step, vinfo%out_tintrv
      start_sec = vinfo%start_sec + (istep - 1) * vinfo%dt

      LOG_INFO("regrid_tool",'(a,i4)') ' interpolate :' // trim(out_vinfo(vid)%varname) // " step=", istep
      if ( associated( out_mesh%ptr_mesh3D ) ) then

        call regrid_interp_field_Interpolate( istep, vinfo%varname, &
            out_mesh, out_var3D, nodemap                            )

        if ( vintrp%vintrp_typeid == REGRID_VCOORD_MODEL_ID ) then
            call regrid_file_write_var( vid, out_var3D,          &
              start_sec, start_sec + vinfo%dt )
        else
            call vintrp%Update_weight( istep, out_mesh, nodemap )               
            call vintrp%Interpolate( istep, out_mesh%ptr_mesh3D, out_var3D )

            call regrid_file_write_var( vid, vintrp%vintrp_var3D, &
              start_sec, start_sec + vinfo%dt )               
        end if

      else

        call regrid_interp_field_Interpolate( istep, vinfo%varname, &
            out_mesh, out_var2D, nodemap                            )
        call regrid_file_write_var( vid, out_var2D,          &
            start_sec, start_sec + vinfo%dt )
        
      end if

      if( IO_L ) call flush(IO_FID_LOG)      
    end do
  end do

  LOG_PROGRESS(*) 'END LOOP'
  LOG_NEWLINE
  call PROF_rapend  ('Main', 0)

  !-
  call finalize()

contains

!OCL SERIAL
  subroutine initialize()
    use scale_prc, only: &
        PRC_MPIstart,        &
        PRC_SINGLECOM_setup, &
        PRC_ERRHANDLER_setup
    use scale_const, only: &
        CONST_setup
    use scale_calendar, only: &
        CALENDAR_setup
    use scale_file, only: &
        FILE_setup
    
    use mod_regrid_mesh, only: &
        regrid_mesh_Init
    use mod_regrid_interp_field, only: &
        regrid_interp_field_Init, &
        in_basename
    use mod_regrid_file, only: &
        regrid_file_Init
    
    implicit none

    namelist / PARAM_REGRID_TOOL / &
        debug      

    integer :: ierr
    integer :: comm     ! communicator                      (execution)

    !-----------------------------------------------------------------------

    ! start MPI
    call PRC_MPIstart( comm ) ! [OUT]

    ! setup MPI communicator
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
                              nprocs,  & ! [OUT]
                              myrank,  & ! [OUT]
                              ismaster ) ! [OUT]

    call PRC_ERRHANDLER_setup( use_fpm = .false., & ! [IN]
                                master  = .false.  ) ! [IN]

    ! setup standard I/O
    call IO_setup( "regrid_tool" )

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )

    ! setup profiler
    call PROF_setup
    call PROF_rapstart ('Initialize', 0)

    ! setup constants
    call CONST_setup

    ! setup calendar
    call CALENDAR_setup

    ! setup fie I/O
    call FILE_setup( myrank )

    LOG_NEWLINE
    LOG_INFO("regrid_tool",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_REGRID_TOOL,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_tool",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_tool",*) 'Not appropriate names in namelist PARAM_REGRID_TOOL. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_REGRID_TOOL)

    !
    call regrid_mesh_Init()
    call regrid_interp_field_Init( out_mesh )
    
    if ( associated(out_mesh%ptr_mesh3D) ) then
      call vintrp%Init( out_mesh, nodemap )
      call regrid_file_Init( in_basename, out_vinfo, vintrp%out_mesh_ptr )
    else
      call regrid_file_Init( in_basename, out_vinfo, out_mesh )      
    end if

    !-
    do_output = .true.

    LOG_INFO("regrid_tool",*) 'Setup has been finished.'
    if( IO_L ) call flush(IO_FID_LOG)

    call PROF_rapend ('Initialize', 0)

    return
  end subroutine initialize

!OCL SERIAL
  subroutine finalize
    use scale_prc, only: &
        PRC_mpibarrier, &
        PRC_MPIfinish
    use scale_file, only: &
        FILE_close_all

    use mod_regrid_interp_field, only: regrid_interp_field_Final
    use mod_regrid_mesh, only: regrid_mesh_Final
    use mod_regrid_file, only: regrid_file_Final

    implicit none
    !--------------------------------------

    call PROF_rapstart ('Finalize', 0)

    LOG_NEWLINE
    LOG_INFO("regrid_tool",*) 'Shutdown'

    call regrid_file_Final
    if ( associated(out_mesh%ptr_mesh3D) ) call vintrp%Final()
    call regrid_interp_field_Final( out_mesh )
    call regrid_mesh_Final()

    !
    call FILE_close_all
    call PROF_rapend ('Finalize', 0)
    
    call PROF_rapreport

    ! stop MPI
    call PRC_mpibarrier
    call PRC_MPIfinish

    return
  end subroutine finalize
  
end program regrid_tool