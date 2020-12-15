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
    procedure :: Construct_ModalFilter3D => AtmosMesh_construct_ModalFilter3D
    procedure :: Construct_ModalFilterHV => AtmosMesh_construct_ModalFilterHV
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

    use scale_FILE_monitor_meshfield, only: &
      FILE_monitor_meshfield_set_dim
    
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
  
    integer  :: NeX               = 2
    integer  :: NeY               = 2
    integer  :: NeZ               = 2
    integer  :: NprcX             = 1
    integer  :: NprcY             = 1 
    integer  :: PolyOrder_h       = 2
    integer  :: PolyOrder_v       = 2
    logical  :: LumpedMassMatFlag = .false.

    integer, parameter :: FZ_nmax = 1000
    real(RP) :: FZ(FZ_nmax)

    namelist / PARAM_ATMOS_MESH / &
      dom_xmin, dom_xmax,                          &
      dom_ymin, dom_ymax,                          &
      dom_zmin, dom_zmax,                          &
      FZ,                                          &
      isPeriodicX, isPeriodicY, isPeriodicZ,       &
      NeX, NeY, NeZ,                               &
      PolyOrder_h, PolyOrder_v, LumpedMassMatFlag, &
      NprcX, NprcY
    
    integer :: n
    character(len=H_SHORT) :: dim_type
    class(LocalMesh3D), pointer :: lcmesh 

    integer :: k
    logical :: is_spec_FZ
    
    integer :: ierr
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_MESH_setup",*) 'Setup'

    FZ(:) = -1.0_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_MESH,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("ATMOS_MESH_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("ATMOS_MESH_setup",*) 'Not appropriate names in namelist PARAM_ATM_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_MESH)

    !----

    ! Setup the element

    call this%element%Init( PolyOrder_h, PolyOrder_v, LumpedMassMatFlag )

    ! Setup the mesh
    
    is_spec_FZ = .true.
    do k=1, NeZ+1
      if (FZ(k) < 0.0_RP) then
        is_spec_FZ = .false.
      end if
    end do
    if (is_spec_FZ) then
      call this%mesh%Init( &
        NprcX*NeX, NprcY*NeY, NeZ,                                 &
        dom_xmin, dom_xmax,dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
        isPeriodicX, isPeriodicY, isPeriodicZ,                     &
        this%element, ATMOS_MESH_NLocalMeshPerPrc, NprcX, NprcY,   &
        FZ=FZ(1:NeZ+1)    )
    else
      call this%mesh%Init( &
        NprcX*NeX, NprcY*NeY, NeZ,                                 &
        dom_xmin, dom_xmax,dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
        isPeriodicX, isPeriodicY, isPeriodicZ,                     &
        this%element, ATMOS_MESH_NLocalMeshPerPrc, NprcX, NprcY    )
    end if
    
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

    !-
    call FILE_monitor_meshfield_set_dim( this%mesh, 'ATM3D' )

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

  subroutine AtmosMesh_construct_ModalFilter3D( this, &
    filter,                                           &
    etac_h, alpha_h, ord_h,                           & 
    etac_v, alpha_v, ord_v                            )
    
    use scale_element_modalfilter, only: ModalFilter
    implicit none
    class(AtmosMesh), intent(in) :: this
    class(ModalFilter), intent(inout) :: filter
    real(RP), intent(in) :: etac_h
    real(RP), intent(in) :: alpha_h
    integer, intent(in) :: ord_h
    real(RP), intent(in) :: etac_v
    real(RP), intent(in) :: alpha_v
    integer, intent(in) :: ord_v
    !-------------------------------------------

    call filter%Init( this%element, &
      etac_h, alpha_h, ord_h,       & 
      etac_v, alpha_v, ord_v        )

    return
  end subroutine AtmosMesh_construct_ModalFilter3D

  subroutine AtmosMesh_construct_ModalFilterHV( this, &
    filterH3D, filterV1D,                             &
    etac_h, alpha_h, ord_h,                           & 
    etac_v, alpha_v, ord_v                            )
    
    use scale_element_modalfilter, only: ModalFilter
    use scale_element_line, only: LineElement
    implicit none
    class(AtmosMesh), intent(in) :: this
    class(ModalFilter), intent(inout) :: filterH3D
    class(ModalFilter), intent(inout) :: filterV1D
    real(RP), intent(in) :: etac_h
    real(RP), intent(in) :: alpha_h
    integer, intent(in) :: ord_h
    real(RP), intent(in) :: etac_v
    real(RP), intent(in) :: alpha_v
    integer, intent(in) :: ord_v

    type(LineElement) :: elemV1D
    !-------------------------------------------

    call filterH3D%Init( this%element, &
      etac_h, alpha_h, ord_h,          & 
      1.0_RP,  0.0_RP, ord_v           )

    call elemV1D%Init( this%element%PolyOrder_v, this%element%IsLumpedMatrix() )      
    call filterV1D%Init( elemV1D, &
      etac_v, alpha_v, ord_v,     &
      tend_flag = .true.          )
    call elemV1D%Final()

    return
  end subroutine AtmosMesh_construct_ModalFilterHV

end module mod_atmos_mesh