!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_regrid_file
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
  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  
  use mod_regrid_mesh_base, only: &
    regrid_mesh_base
  use mod_regrid_outvar_info, only: &
    OutVarInfoList, OutVarInfo
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  public :: regrid_file_Init
  interface regrid_file_write_var
    module procedure regrid_file_write_var2D
    module procedure regrid_file_write_var3D
  end interface
  public :: regrid_file_write_var
  public :: regrid_file_Final

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
  subroutine regrid_file_Init( in_basename, out_vinfo, out_vinfo_oper, out_mesh )
    use scale_file_h
    use mod_regrid_outvar_info, only: OutVarInfoList
    implicit none

    character(*), intent(in) :: in_basename
    type(OutVarInfoList), intent(inout) :: out_vinfo
    type(OutVarInfoList), intent(inout) :: out_vinfo_oper
    class(regrid_mesh_base), intent(in) :: out_mesh

    character(len=H_LONG )   :: out_basename     = ''       ! Basename of the output file
    character(len=H_MID)     :: out_title        = ''        !< Title    of the output file
    character(len=H_SHORT)   :: out_dtype        = 'DEFAULT' !< REAL4 or REAL8

    namelist / PARAM_REGRID_FILE / &
      out_basename,   &
      out_title,      &
      out_dtype,      &
      out_UniformGrid
    
    integer :: ierr
    logical :: fileexisted

    type(FILE_base_meshfield) :: in_file
    integer :: nn
    integer :: var_num

    character(len=FILE_HMID) :: tunits
    character(len=FILE_HSHORT) :: calendar
    character(len=FILE_HMID) :: desc
    character(len=FILE_HMID) :: standard_name

    integer :: dimtype_id
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("regrid_file",*) 'Setup'
  
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_REGRID_FILE,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_file",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_file",*) 'Not appropriate names in namelist PARAM_REGRID_FILE. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_REGRID_FILE)


    !--
    var_num = out_vinfo%item_num + out_vinfo_oper%item_num

    call file_init( in_file, dimtype_id, 1, out_mesh )
    call in_file%Open( in_basename, myrank=0 )
    if (out_title=='') then
      call in_file%Get_commonInfo( title=out_title ) ! (out)
    end if

    call in_file%Get_dataInfo( out_vinfo%items(1)%varname, 1,     & ! (in)
      time_units=tunits, calendar=calendar                        ) ! (out)

    call file_init( out_file, dimtype_id, var_num, out_mesh, out_UniformGrid )

    call out_file%Create( out_basename, out_title, out_dtype,                 & ! (in)
                          fileexisted,                                        & ! (out)
                          myrank=PRC_myrank, calendar=calendar, tunits=tunits ) ! (in)

    call out_vinfo%DefVarForNetCDF( out_file, dimtype_id, out_dtype, 0 )
    call out_vinfo_oper%DefVarForNetCDF( out_file, dimtype_id, out_dtype, out_vinfo%item_num )

    call out_file%End_def()

    call in_file%Close()
    call in_file%Final()
                                                
    return
  end subroutine regrid_file_Init

  subroutine regrid_file_write_var2D( vinfo, field, istep )
    implicit none
    type(OutVarInfo), intent(in) :: vinfo
    class(MeshField2D), intent(in) :: field
    integer, intent(in) :: istep

    real(DP) :: start_sec
    !-------------------------------------------

    call PROF_rapstart('regrid_file_write_var2D', 0)
    
    start_sec = vinfo%start_sec + dble(istep - 1) * vinfo%dt

    call out_file%Write_var2D( vinfo%vidForOutput, field, &
      start_sec, start_sec + vinfo%dt )

    call PROF_rapend('regrid_file_write_var2D', 0)

    return
  end subroutine regrid_file_write_var2D

  subroutine regrid_file_write_var3D(  vinfo, field, istep )
    implicit none
    type(OutVarInfo), intent(in) :: vinfo
    class(MeshField3D), intent(in) :: field
    integer, intent(in) :: istep

    real(DP) :: start_sec
    !-------------------------------------------

    call PROF_rapstart('regrid_file_write_var3D', 0)
    
    start_sec = vinfo%start_sec + dble(istep - 1) * vinfo%dt

    call out_file%Write_var3D( vinfo%vidForOutput, field, &
      start_sec, start_sec + vinfo%dt )

    call PROF_rapend('regrid_file_write_var3D', 0)
    
    return
  end subroutine regrid_file_write_var3D  

  subroutine regrid_file_Final()
    implicit none
    !-------------------------------------------

    call out_file%Close()
    call out_file%Final()
    return
  end subroutine regrid_file_Final

!-- private -----------------------------------

  subroutine file_init( file, dimtype_id, &
    var_num, mesh, force_uniform_grid     )

    use scale_mesh_base2d, only: MeshBase2D
    use scale_mesh_base3d, only: MeshBase3D
    use scale_mesh_rectdom2d, only: MeshRectDom2D
    use scale_mesh_cubedom3d, only: MeshCubeDom3D
    use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
    use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D

    use scale_mesh_base2d, only: &
      DIMTYPE2D_XYT  => MeshBase2D_DIMTYPEID_XYT
    use scale_mesh_base3d, only: &
      DIMTYPE3D_XYZT  => MeshBase3D_DIMTYPEID_XYZT  
      
    implicit none

    type(FILE_base_meshfield), intent(inout) :: file
    integer, intent(out) :: dimtype_id
    integer, intent(in) :: var_num
    class(regrid_mesh_base), intent(in) :: mesh
    logical, intent(in), optional :: force_uniform_grid

    class(MeshBase2D), pointer :: ptr_mesh2D
    class(MeshBase3D), pointer :: ptr_mesh3D

    logical :: force_uniform_grid_
    !-------------------------------------------

    if ( present(force_uniform_grid) ) then
      force_uniform_grid_ = force_uniform_grid
    else
      force_uniform_grid_ = .false.
    end if

    if ( associated( mesh%ptr_mesh2D ) ) then

      select type( ptr_mesh2D => mesh%ptr_mesh2D  )
      class is ( MeshRectDom2D )
        call file%Init( var_num, mesh2D=ptr_mesh2D, force_uniform_grid=force_uniform_grid_ )
      class is ( MeshCubedSphereDom2D )
        call file%Init( var_num, meshCubedSphere2D=ptr_mesh2D, force_uniform_grid=force_uniform_grid_ )
      end select
      dimtype_id = DIMTYPE2D_XYT

    else if ( associated( mesh%ptr_mesh3D ) ) then

      select type( ptr_mesh3D => mesh%ptr_mesh3D  )
      class is ( MeshCubeDom3D )
        call file%Init( var_num, mesh3D=ptr_mesh3D, force_uniform_grid=force_uniform_grid_ )
      class is ( MeshCubedSphereDom3D )
        call file%Init( var_num, meshCubedSphere3D=ptr_mesh3D, force_uniform_grid=force_uniform_grid_ )
      end select
      dimtype_id = DIMTYPE3D_XYZT

    end if

    return
  end subroutine file_init

end module mod_regrid_file