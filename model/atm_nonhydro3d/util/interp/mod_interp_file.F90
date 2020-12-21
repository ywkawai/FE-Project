!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_interp_file
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
    PRC_myrank, PRC_abort

  use scale_file_base_meshfield, only: &
    FILE_base_meshfield
  use scale_mesh_cubedom3d, only: &
    MeshCubeDom3D
  use scale_meshfield_base, only: &
    MeshField3D
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  public :: interp_file_Init
  public :: interp_file_write_var
  public :: interp_file_Final

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
 
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  type(FILE_base_meshfield) :: out_file
  logical, private  :: out_UniformGrid = .false. 

contains
  subroutine interp_file_Init( in_basename, out_vinfo, mesh3D )
    use scale_file_h
    use scale_file_common_meshfield, only: &
      DIMTYPE_XYZT => FILE_COMMON_MESHFILED3D_DIMTYPEID_XYZT
    use mod_interp_field, only: OutVarInfo
    implicit none

    character(*), intent(in) :: in_basename
    type(OutVarInfo), intent(in) :: out_vinfo(:)
    class(MeshCubeDom3D), intent(in) :: mesh3D

    character(len=H_LONG )   :: out_basename     = ''       ! Basename of the output file
    character(len=H_MID)     :: out_title        = ''        !< Title    of the output file
    character(len=H_SHORT)   :: out_dtype        = 'DEFAULT' !< REAL4 or REAL8

    namelist / PARAM_INTERP_FILE / &
      out_basename, &
      out_title,    &
      out_dtype,    &
      out_UniformGrid
    
    integer :: ierr
    logical :: fileexisted

    type(FILE_base_meshfield) :: in_file
    integer :: nn
    integer :: var_num

    real(DP) :: tstart, tend
    real(DP) :: tintv
    integer :: num_step
    character(len=FILE_HMID) :: tunits
    character(len=FILE_HSHORT) :: calendar
    character(len=FILE_HMID) :: desc
    character(len=FILE_HMID) :: standard_name
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("interp_file",*) 'Setup'
  
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_INTERP_FILE,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
        LOG_INFO("interp_file",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("interp_file",*) 'Not appropriate names in namelist PARAM_INTERP_FILE. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_INTERP_FILE)


    !--
    var_num = size(out_vinfo)

    call in_file%Init( var_num, mesh3D=mesh3D )
    call in_file%Open( in_basename, myrank=0 )
    if (out_title=='') then
      call in_file%Get_commonInfo( title=out_title ) ! (out)
    end if

    call in_file%Get_dataInfo( out_vinfo(1)%varname, 1,     & ! (in)
      time_units=tunits, calendar=calendar                  ) ! (out)

    call out_file%Init( var_num, mesh3D=mesh3D, force_uniform_grid=out_UniformGrid )
    call out_file%Create( out_basename, out_title, out_dtype,                 & ! (in)
                          fileexisted,                                        & ! (out)
                          myrank=PRC_myrank, calendar=calendar, tunits=tunits ) ! (in)

    do nn=1, var_num
      call in_file%Get_dataInfo( out_vinfo(nn)%varname, 1,     & ! (in)
        description=desc, standard_name=standard_name,         & ! (in)
        time_start=tstart, time_end=tend                       ) ! (out)
      call in_file%Get_VarStepSize( out_vinfo(nn)%varname,     & ! (in)
        num_step                                               )
      
      call out_file%Def_Var( out_vinfo(nn)%varname, out_vinfo(nn)%units, &
        desc, nn, DIMTYPE_XYZT,  out_dtype,                              &
        standard_name=standard_name, timeinv=out_vinfo(nn)%dt  )
    end do
    call out_file%End_def()

    call in_file%Close()
    call in_file%Final()
                                                
    return
  end subroutine interp_file_Init

  subroutine interp_file_write_var( vid, field, start_sec, end_sec )
    implicit none
    integer, intent(in) :: vid
    class(MeshField3D), intent(in) :: field
    real(DP), intent(in) :: start_sec
    real(DP), intent(in) :: end_sec
    !-------------------------------------------

    call out_file%Write_var3D( vid, field, start_sec, end_sec )
    return
  end subroutine interp_file_write_var

  subroutine interp_file_Final()
    implicit none
    !-------------------------------------------

    call out_file%Close()
    call out_file%Final()
    return
  end subroutine interp_file_Final

!-- private -----------------------------------


end module mod_interp_file