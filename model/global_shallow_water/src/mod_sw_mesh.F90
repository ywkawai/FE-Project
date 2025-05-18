!-------------------------------------------------------------------------------
!> module  Global SW / Mesh
!!
!! @par Description
!!          Module to manage the mesh for global shallow water equation
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_sw_mesh
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_meshfield_base, only: MeshField3D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_meshfieldcomm_cubedspheredom2d, only: MeshFieldCommCubedSphereDom2D
  use scale_sparsemat, only: sparsemat
  use scale_model_mesh_manager, only: &
    ModelMesh2D
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, extends(ModelMesh2D), public :: SWMesh
    type(MeshCubedSphereDom2D) :: mesh
    type(QuadrilateralElement) :: element
  contains
    procedure :: Init => SWMesh_Init
    procedure :: Final => SWMesh_Final
    procedure :: Construct_ModalFilter2D => SWMesh_construct_ModalFilter2D
  end type SWMesh
  
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

contains
  subroutine SWMesh_Init( this )
    use scale_const, only: &
      RPlanet => CONST_RADIUS
    use scale_FILE_monitor_meshfield, only: &
      FILE_monitor_meshfield_set_dim
    
    implicit none
    class(SWMesh), target, intent(inout) :: this
  
    integer  :: NeGX               = 2
    integer  :: NeGY               = 2
    integer  :: Nprc               = 1
    integer  :: NLocalMeshPerPrc   = 6
    logical  :: LumpedMassMatFlag = .false.
    integer :: PolyOrder           = 1

    namelist / PARAM_SW_MESH / &
      NeGX, NeGY, NLocalMeshPerPrc, &
      Nprc,                         &
      PolyOrder, LumpedMassMatFlag
    
    integer :: n
    character(len=H_SHORT) :: dim_type
    class(LocalMesh2D), pointer :: lcmesh 

    character(len=H_SHORT) :: SpMV_storage_format = 'ELL' ! CSR or ELL

    integer :: ierr
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SW_MESH_setup",*) 'Setup'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SW_MESH,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("SW_MESH_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("SW_MESH_setup",*) 'Not appropriate names in namelist PARAM_SW_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_SW_MESH)

    !----

    ! Setup the element

    call this%element%Init( PolyOrder, LumpedMassMatFlag )

    ! Setup the mesh
    call this%mesh%Init( &
      NeGX, NeGY, RPlanet,                    &
      this%element, NLocalMeshPerPrc, Nprc    )
    
    call this%mesh%Generate()
    
    !-
    call this%ModelMesh2D_Init( this%mesh )

    call this%DOptrMat(1)%Init( this%element%Dx1, storage_format=SpMV_storage_format )
    call this%DOptrMat(2)%Init( this%element%Dx2, storage_format=SpMV_storage_format )

    call this%SOptrMat(1)%Init( this%element%Sx1, storage_format=SpMV_storage_format )
    call this%SOptrMat(2)%Init( this%element%Sx2, storage_format=SpMV_storage_format )

    call this%LiftOptrMat%Init( this%element%Lift, storage_format=SpMV_storage_format )

    !-
    call FILE_monitor_meshfield_set_dim( this%mesh, 'ATM2D' )

    return
  end subroutine SWMesh_Init

  subroutine SWMesh_Final(this)
    implicit none

    class(SWMesh), intent(inout) :: this
    !-------------------------------------------

    call this%mesh%Final()
    call this%ModelMesh2D_Final()

    return
  end subroutine SWMesh_Final

  subroutine SWMesh_construct_ModalFilter2D( this, &
    filter,                                        &
    etac, alpha, ord                               )
    
    use scale_element_modalfilter, only: ModalFilter
    implicit none
    class(SWMesh), intent(in) :: this
    class(ModalFilter), intent(inout) :: filter
    real(RP), intent(in) :: etac
    real(RP), intent(in) :: alpha
    integer, intent(in) :: ord
    !-------------------------------------------

    call filter%Init( this%element, &
      etac, alpha, ord              )

    return
  end subroutine SWMesh_construct_ModalFilter2D

end module mod_sw_mesh