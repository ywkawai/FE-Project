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

  use scale_element_base, only: ElementBase2D, ElementBase3D
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

  type(regrid_mesh_base), public, target :: out_mesh

  logical, public :: out_GarlerkinProjection
  type(regrid_mesh_base), public, target :: out_mesh_GP
  real(RP), public, allocatable :: GPMat(:,:)

  type(regrid_nodemap), public, allocatable :: nodemap(:)

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
!OCL SERIAL
  subroutine regrid_mesh_Init()
    implicit none

    character(len=H_MID) :: in_meshType
    character(len=H_MID) :: out_meshType

    namelist / PARAM_REGRID_MESH / &
      in_MeshType,             &
      out_MeshType,            &
      out_GarlerkinProjection
        
    integer :: ierr
    integer :: n

    type(regrid_mesh_base), pointer :: ptr_out_mesh

    !-------------------------------

    LOG_NEWLINE
    LOG_INFO("regrid_mesh",*) 'Setup .. '
    
    !--- read namelist (PARAM_REGRID_MESH)

    out_GarlerkinProjection = .false.

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
    call out_mesh%Init( MESHOUT_ID, out_meshType, .false. )
    call out_mesh%Generate()
    ptr_out_mesh => out_mesh

    if ( out_GarlerkinProjection ) then
      call out_mesh_GP%Init( MESHOUT_ID, out_meshType, .true. )
      call out_mesh_GP%Generate()
      ptr_out_mesh => out_mesh

      if ( out_mesh%polyorder_v < 0 ) then
        allocate( GPMat(out_mesh%elem2D%Np,out_mesh_GP%elem2D%Np) )
        call out_mesh%elem2D%Generate_L2ProjMat( out_mesh_GP%elem2D, & ! (in)
          GPMat ) ! (out)
      else
        allocate( GPMat(out_mesh%elem3D%Np,out_mesh_GP%elem3D%Np) )
        call out_mesh%elem3D%Generate_L2ProjMat( out_mesh_GP%elem3D, & ! (in)
          GPMat ) ! (out)
      end if
    end if

    !--
    LOG_INFO("regrid_mesh",*) 'Initialize node mapping ..'        
    allocate( nodemap( ptr_out_mesh%NLocalMeshPerPRC ) )
    do n=1, size(nodemap)
      call nodemap(n)%Init( in_meshType, ptr_out_mesh )
    end do

    return
  end subroutine regrid_mesh_Init

!OCL SERIAL
  subroutine regrid_mesh_Final()
    implicit none

    integer :: n
    !-------------------------------

    call out_mesh%Final()

    do n=1, size(nodemap)
      call nodemap(n)%Final()
    end do
    deallocate( nodemap )

    if ( allocated(GPMat) ) deallocate( GPMat )

    return
  end subroutine regrid_mesh_Final
end module mod_regrid_mesh
