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

  use scale_model_mesh_manager, only: &
    ModelMesh3D
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, extends(ModelMesh3D), public :: AtmosMesh
    type(MeshCubeDom3D) :: mesh
    type(HexahedralElement) :: element
  contains
    procedure :: Init => AtmosMesh_Init
    procedure :: Final => AtmosMesh_Final
  end type AtmosMesh
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, parameter :: ATMOS_MESH_NLocalMeshPerPrc = 1
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
  subroutine AtmosMesh_Init( this )
    implicit none
    class(AtmosMesh), target, intent(inout) :: this

    real(RP) :: dom_xmin         = 0.0_RP 
    real(RP) :: dom_xmax         = 100.0E3_RP
    real(RP) :: dom_ymin         = 0.0_RP 
    real(RP) :: dom_ymax         = 100.0E3_RP
    real(RP) :: dom_zmin         = 0.0_RP
    real(RP) :: dom_zmax         = 10.0E3_RP
    logical  :: isPeriodicX       = .true.
    logical  :: isPeriodicY       = .true.
    logical  :: isPeriodicZ       = .false.
  
    integer  :: NeGX              = 2
    integer  :: NeGY              = 2
    integer  :: NeGZ              = 2
    integer  :: PolyOrder_h       = 2
    integer  :: PolyOrder_v       = 2
    logical  :: LumpedMassMatFlag = .false.

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

    call this%element%Init( PolyOrder_h, PolyOrder_v, LumpedMassMatFlag )
    
    call this%mesh%Init( &
      NeGX, NeGY, NeGZ,                                          &
      dom_xmin, dom_xmax,dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
      isPeriodicX, isPeriodicY, isPeriodicZ,                     &
      this%element, ATMOS_MESH_NLocalMeshPerPrc )
    
    call this%mesh%Generate()
    
    !-
    call this%ModelMesh3D_Init( this%mesh )

    call this%DOptrMat(1)%Init( this%element%Dx1 )
    call this%DOptrMat(2)%Init( this%element%Dx2 )
    call this%DOptrMat(3)%Init( this%element%Dx3 )

    call this%SOptrMat(1)%Init( this%element%Sx1 )
    call this%SOptrMat(2)%Init( this%element%Sx2 )
    call this%SOptrMat(3)%Init( this%element%Sx3 )

    call this%LiftOptrMat%Init( this%element%Lift )

    return
  end subroutine AtmosMesh_Init

  subroutine AtmosMesh_Final(this)
    implicit none

    class(AtmosMesh), intent(inout) :: this
    !-------------------------------------------

    call this%mesh%Final()
    call this%ModelMesh3D_Final()

    return
  end subroutine AtmosMesh_Final

end module mod_atmos_mesh