!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_regrid_mesh
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort
  use scale_const, only: &
    PI => CONST_PI,         &
    RPlanet => CONST_radius

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D    

  use mod_regrid_mesh_base, only: &
    regrid_mesh_base,                    &
    MESHIN_ID => REGRID_MESH_BASE_IN_ID, &
    MESHOUT_ID => REGRID_MESH_BASE_OUT_ID

  use mod_regrid_nodemap, only: &
    regrid_nodemap
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  public :: regrid_mesh_Init
  public :: regrid_mesh_Final

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  type(regrid_mesh_base), public :: out_mesh
  type(regrid_nodemap), public, allocatable :: nodemap(:)

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !



contains
  subroutine regrid_mesh_Init()
    implicit none

    character(len=H_MID) :: in_meshType
    character(len=H_MID) :: out_meshType

    namelist / PARAM_REGRID_MESH / &
      in_MeshType, &
      out_MeshType
        
    integer :: ierr
    integer :: n
    !-------------------------------

    LOG_NEWLINE
    LOG_INFO("regrid_mesh",*) 'Setup .. '
    
    !--- read namelist (PARAM_REGRID_MESH)
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_REGRID_MESH,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_mesh",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_mesh",*) 'Not appropriate names in namelist PARAM_REGRID_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_REGRID_MESH)

    !--
    LOG_INFO("regrid_mesh",*) 'Initialize mesh for output data..'
    if( IO_L ) call flush(IO_FID_LOG)

    call out_mesh%Init( MESHOUT_ID, out_meshType )
    call out_mesh%Generate()
    
    !--
    LOG_INFO("regrid_mesh",*) 'Initialize node mapping ..'        
    if( IO_L ) call flush(IO_FID_LOG)

    allocate( nodemap( out_mesh%NLocalMeshPerPRC ) )
    do n=1, size(nodemap)
      call nodemap(n)%Init( in_meshType, out_mesh )
    end do

    return
  end subroutine regrid_mesh_Init

  subroutine regrid_mesh_Final()
    implicit none

    integer :: n
    !-------------------------------

    call out_mesh%Final()

    do n=1, size(nodemap)
      call nodemap(n)%Final()
    end do
    deallocate( nodemap )

    return
  end subroutine regrid_mesh_Final
  
end module mod_regrid_mesh
