!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_mesh
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_meshfield_base, only: MeshField3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
  use scale_sparsemat, only: sparsemat

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: ATMOS_MESH_setup
  public :: ATMOS_MESH_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  type(MeshCubeDom3D),        public, save, target :: mesh
  type(HexahedralElement), public, save :: refElem
  integer, public, parameter :: NLocalMeshPerPrc = 1
  
  type(sparsemat), public :: Dx
  type(sparsemat), public :: Sx
  type(sparsemat), public :: Dy
  type(sparsemat), public :: Sy  
  type(sparsemat), public :: Dz
  type(sparsemat), public :: Sz
  type(sparsemat), public :: Lift

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: dom_xmin         = 0.0_RP 
  real(RP), private :: dom_xmax         = 100.0E3_RP
  real(RP), private :: dom_ymin         = 0.0_RP 
  real(RP), private :: dom_ymax         = 100.0E3_RP
  real(RP), private :: dom_zmin         = 0.0_RP
  real(RP), private :: dom_zmax         = 10.0E3_RP
  logical, private :: isPeriodicX       = .true.
  logical, private :: isPeriodicY       = .true.
  logical, private :: isPeriodicZ       = .false.

  integer, private :: NeGX              = 2
  integer, private :: NeGY              = 2
  integer, private :: NeGZ              = 2
  integer, public  :: PolyOrder_h       = 2
  integer, public  :: PolyOrder_v       = 2
  logical, private :: LumpedMassMatFlag = .false.

contains
  subroutine ATMOS_MESH_setup()
    implicit none

    namelist / PARAM_ATMOS_MESH / &
      dom_xmin, dom_xmax,                        &
      dom_ymin, dom_ymax,                        &
      dom_zmin, dom_zmax,                        &
      isPeriodicX, isPeriodicY, isPeriodicZ,     &
      NeGX, NeGY, NeGZ, PolyOrder_h, PolyOrder_v
    
    integer :: ierr
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_MESH_driver_setup",*) 'Setup'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_MESH,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("ATMOS_MESH_driver_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("ATMOS_MESH_driver_setup",*) 'Not appropriate names in namelist PARAM_ATM_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_MESH)
    !----

    call refElem%Init(PolyOrder_h, PolyOrder_v, LumpedMassMatFlag)
    call Dx%Init(refElem%Dx1)
    call Sx%Init(refElem%Sx1)
    call Dy%Init(refElem%Dx2)
    call Sy%Init(refElem%Sx2)
    call Dz%Init(refElem%Dx3)
    call Sz%Init(refElem%Sx3)
    call Lift%Init(refElem%Lift)

    call mesh%Init( &
      NeGX, NeGY, NeGZ,                                          &
      dom_xmin, dom_xmax,dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
      isPeriodicX, isPeriodicY, isPeriodicZ,                     &
      refElem, NLocalMeshPerPrc )
    
    call mesh%Generate()

    return
  end subroutine ATMOS_MESH_setup

  subroutine ATMOS_MESH_finalize
    implicit none
    !-------------------------------------------

    call mesh%Final()
    
    call Dx%Final()
    call Sx%Final()
    call Dy%Final()
    call Sy%Final()
    call Dz%Final()
    call Sz%Final()
    call Lift%Final()
    call refElem%Final()

    return
  end subroutine ATMOS_MESH_finalize

end module mod_atmos_mesh