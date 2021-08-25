!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_cs2lonlat_interp_file
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
  use scale_mesh_rectdom2d, only: &
    MeshRectDom2D
  use scale_mesh_cubedom3d, only: &
    MeshCubeDom3D
  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  public :: interp_file_Init
  interface interp_file_write_var
    module procedure interp_file_write_var2D
    module procedure interp_file_write_var3D
  end interface
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
  subroutine interp_file_Init( in_basename, out_vinfo, mesh2D, mesh3D, is_mesh3D )
    use scale_file_h
    use scale_mesh_base2d, only: &
      DIMTYPE2D_XYT  => MeshBase2D_DIMTYPEID_XYT
    use scale_mesh_base3d, only: &
      DIMTYPE3D_XYZT  => MeshBase3D_DIMTYPEID_XYZT
    use mod_cs2lonlat_interp_field, only: OutVarInfo
    implicit none

    character(*), intent(in) :: in_basename
    type(OutVarInfo), intent(in) :: out_vinfo(:)
    class(MeshRectDom2D), intent(in) :: mesh2D
    class(MeshCubeDom3D), intent(in) :: mesh3D
    logical, intent(in) :: is_mesh3D

    character(len=H_LONG )   :: out_basename     = ''       ! Basename of the output file
    character(len=H_MID)     :: out_title        = ''        !< Title    of the output file
    character(len=H_SHORT)   :: out_dtype        = 'DEFAULT' !< REAL4 or REAL8

    namelist / PARAM_INTERP_FILE / &
      out_basename,   &
      out_title,      &
      out_dtype,      &
      out_UniformGrid
    
    integer :: ierr
    logical :: fileexisted

    type(FILE_base_meshfield) :: in_file
    integer :: nn
    integer :: var_num

    integer :: dimtype
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

    if (is_mesh3D) then
      call in_file%Init( var_num, mesh3D=mesh3D )
      dimtype = DIMTYPE3D_XYZT
    else
      call in_file%Init( var_num, mesh2D=mesh2D )
      dimtype = DIMTYPE2D_XYT
    end if

    call in_file%Open( in_basename, myrank=0 )
    if (out_title=='') then
      call in_file%Get_commonInfo( title=out_title ) ! (out)
    end if

    call in_file%Get_dataInfo( out_vinfo(1)%varname, 1,     & ! (in)
      time_units=tunits, calendar=calendar                  ) ! (out)

    if (is_mesh3D) then
      call out_file%Init( var_num, mesh3D=mesh3D, force_uniform_grid=out_UniformGrid )
    else
      call out_file%Init( var_num, mesh2D=mesh2D, force_uniform_grid=out_UniformGrid )
    end if
    call out_file%Create( out_basename, out_title, out_dtype,                 & ! (in)
                          fileexisted,                                        & ! (out)
                          myrank=PRC_myrank, calendar=calendar, tunits=tunits ) ! (in)

    do nn=1, var_num      
      call out_file%Def_Var( out_vinfo(nn)%varname, out_vinfo(nn)%units, &
        desc, nn, dimtype, out_dtype,                                    &
        standard_name=standard_name,                                     &
        timeinv=out_vinfo(nn)%dt * dble(out_vinfo(nn)%out_tintrv)        )
    end do
    call out_file%End_def()

    call in_file%Close()
    call in_file%Final()
                                                
    return
  end subroutine interp_file_Init

  subroutine interp_file_write_var2D( vid, field, start_sec, end_sec )
    implicit none
    integer, intent(in) :: vid
    class(MeshField2D), intent(in) :: field
    real(DP), intent(in) :: start_sec
    real(DP), intent(in) :: end_sec
    !-------------------------------------------

    call PROF_rapstart('INTERP_file_write_var2D', 0)
    
    call out_file%Write_var2D( vid, field, start_sec, end_sec )

    call PROF_rapend('INTERP_file_write_var2D', 0)
    return
  end subroutine interp_file_write_var2D

  subroutine interp_file_write_var3D( vid, field, start_sec, end_sec )
    implicit none
    integer, intent(in) :: vid
    class(MeshField3D), intent(in) :: field
    real(DP), intent(in) :: start_sec
    real(DP), intent(in) :: end_sec
    !-------------------------------------------

    call PROF_rapstart('INTERP_file_write_var3D', 0)
    
    call out_file%Write_var3D( vid, field, start_sec, end_sec )

    call PROF_rapend('INTERP_file_write_var3D', 0)
    return
  end subroutine interp_file_write_var3D  

  subroutine interp_file_Final()
    implicit none
    !-------------------------------------------

    call out_file%Close()
    call out_file%Final()
    return
  end subroutine interp_file_Final

!-- private -----------------------------------

end module mod_cs2lonlat_interp_file