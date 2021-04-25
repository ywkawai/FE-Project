!-------------------------------------------------------------------------------
!> Program convert_cs2lonlat (SCALE-DG)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Interpolation tool (convert_cs2lonlat)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
program convert_cs2lonlat
   !-----------------------------------------------------------------------------
   !
   !++ used modules
   !
   use mpi
   use scale_precision
   use scale_io
   use scale_prof
   use scale_prc, only: PRC_abort

   use mod_cs2lonlat_interp_mesh, only: &
      out_mesh2D,               &
      nodeMap_list  
   use mod_cs2lonlat_interp_field, only: &
      OutVarInfo,              &
      out_var2D_num,           &
      out_var2D,               &
      out_vinfo,               &
      interp_field_Interpolate
   use mod_cs2lonlat_interp_file, only: &
      interp_file_write_var
   
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

   type(OutVarInfo), pointer :: vinfo
   !-----------------------------------------------------------------------------

   call initialize()
   
   !-
   !########## main ##########
   call PROF_rapstart('Main', 0)

   LOG_NEWLINE
   LOG_PROGRESS(*) 'START LOOP'

   do vid =1, out_var2D_num
      vinfo => out_vinfo(vid)
      do istep=1, vinfo%num_step, vinfo%out_tintrv
         start_sec = vinfo%start_sec + (istep - 1) * vinfo%dt

         LOG_INFO("INTERP",'(a,i4)') 'Interpolate :' // trim(out_vinfo(vid)%varname) // " step=", istep
         call interp_field_Interpolate( istep, vinfo%varname, &
            out_mesh2D, out_var2D, nodeMap_list               )
         call interp_file_write_var( vid, out_var2D,          &
            start_sec, start_sec + vinfo%dt )
         
         if( IO_L ) call flush(IO_FID_LOG)      
      end do
   end do

   LOG_PROGRESS(*) 'END LOOP'
   LOG_NEWLINE
   call PROF_rapend  ('Main', 0)

   !-
   call finalize()

contains
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
      
      use mod_cs2lonlat_interp_mesh, only: &
         interp_mesh_init
      use mod_cs2lonlat_interp_field, only: &
         interp_field_Init, &
         in_basename
      use mod_cs2lonlat_interp_file, only: &
         interp_file_Init
      
      implicit none

      namelist / PARAM_INTERP / &
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
      call IO_setup( "convert_cs2lonlat" )
    
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
      LOG_INFO("INTERP",*) 'Setup'
    
      !--- read namelist
      rewind(IO_FID_CONF)
      read(IO_FID_CONF,nml=PARAM_INTERP,iostat=ierr)
      if ( ierr < 0 ) then !--- missing
         LOG_INFO("convert_cs2lonlat",*) 'Not found namelist. Default used.'
      elseif( ierr > 0 ) then !--- fatal error
         LOG_ERROR("convert_cs2lonlat",*) 'Not appropriate names in namelist PARAM_INTERP. Check!'
         call PRC_abort
      endif
      LOG_NML(PARAM_INTERP)
    
      !
      call interp_mesh_Init
      call interp_field_Init( out_mesh2D, nodeMap_list )
      call interp_file_Init( in_basename, out_vinfo, out_mesh2D )

      !-
      do_output = .true.
    
      LOG_INFO("convert_cs2lonlat",*) 'Setup has been finished.'

      call PROF_rapend ('Initialize', 0)

      return
   end subroutine initialize

   subroutine finalize
      use scale_prc, only: &
         PRC_mpibarrier, &
         PRC_MPIfinish
      use scale_file, only: &
         FILE_close_all

      use mod_cs2lonlat_interp_field, only: interp_field_Final
      use mod_cs2lonlat_interp_mesh, only: interp_mesh_Final
      use mod_cs2lonlat_interp_file, only: interp_file_Final

      implicit none
      !--------------------------------------

      call PROF_rapstart ('Finalize', 0)

      call interp_file_Final
      call interp_field_Final
      call interp_mesh_Final

      !
      call FILE_close_all
      call PROF_rapend ('Finalize', 0)
      
      call PROF_rapreport
    
      ! stop MPI
      call PRC_mpibarrier
      call PRC_MPIfinish

      return
   end subroutine finalize
end program convert_cs2lonlat
