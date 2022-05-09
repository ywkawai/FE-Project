!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_mesh_rm
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_meshfield_base, only: MeshField3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
  use scale_sparsemat, only: sparsemat
  use scale_file_restart_meshfield, only: FILE_restart_meshfield_component

  use scale_model_var_manager, only: ModelVarManager
  use scale_model_mesh_manager, only: ModelMesh3D
  
  use mod_atmos_mesh, only: &
    AtmosMesh, ATM_MESH_MAX_COMMNUICATOR_NUM

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, extends(AtmosMesh), public :: AtmosMeshRM
    type(MeshCubeDom3D) :: mesh
    type(MeshFieldCommCubeDom3D) :: comm_list(ATM_MESH_MAX_COMMNUICATOR_NUM)
  contains
    procedure :: Init => AtmosMeshRM_Init
    procedure :: Final => AtmosMeshRM_Final
    procedure :: Create_communicator => AtmosMeshRM_Create_communicator
    procedure :: Setup_restartfile1 => AtmosMeshRM_setup_restartfile1
    procedure :: Setup_restartfile2 => AtmosMeshRM_setup_restartfile2
    procedure :: Calc_UVmet => AtmosMeshRM_calc_UVmet
    procedure :: Setup_vcoordinate => AtmosMeshRM_setup_vcoordinate
  end type AtmosMeshRM

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

  !- AtmosMesh RM -----------------------------------------

  subroutine AtmosMeshRM_Init( this )
    use scale_file_base_meshfield, only: FILE_base_meshfield
    use scale_mesh_base2d, only: &
      MFTYPE2D_XY => MeshBase2D_DIMTYPEID_XY
    use scale_meshutil_vcoord, only: &
      MeshUtil_get_VCoord_TypeID

    implicit none
    class(AtmosMeshRM), target, intent(inout) :: this

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
    character(len=H_LONG) :: TOPO_IN_BASENAME    = ''     !< basename of the input file
    character(len=H_MID)  :: TOPO_IN_VARNAME     = 'topo' !< variable name of topo in the input file
    character(len=H_MID)  :: VERTICAL_COORD_NAME = "TERRAIN_FOLLOWING"

    integer, parameter :: ATMOS_MESH_NLocalMeshPerPrc = 1

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
      NprcX, NprcY,                                &
      TOPO_IN_BASENAME, TOPO_IN_VARNAME,           &
      VERTICAL_COORD_NAME
    
    integer :: k
    logical :: is_spec_FZ

    integer :: ierr

    type(FILE_base_meshfield) :: file_topo
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
    call this%element_v1D%Init( PolyOrder_v, LumpedMassMatFlag )

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
    call this%AtmosMesh_Init( this%mesh )

    !- Set topography & vertical coordinate

    if ( TOPO_IN_BASENAME /= '' ) then
      LOG_INFO("ATMOS_MESH_setup",*) 'Read topography data'

      call file_topo%Init(1, mesh2D=this%mesh%mesh2D )
      call file_topo%Open( TOPO_IN_BASENAME, myrank=PRC_myrank )
      call file_topo%Read_Var( MFTYPE2D_XY, TOPO_IN_VARNAME, this%topography%topo )
      call file_topo%Close()
      call file_topo%Final()
    end if

    this%vcoord_type_id = MeshUtil_get_VCoord_TypeID( VERTICAL_COORD_NAME )
    call this%Setup_vcoordinate()

    return
  end subroutine AtmosMeshRM_Init

  subroutine AtmosMeshRM_Final(this)
    implicit none

    class(AtmosMeshRM), intent(inout) :: this
    integer :: commid
    !-------------------------------------------

    do commid=1, this%communicator_num
      call this%comm_list(commid)%Final()
    end do

    call this%mesh%Final()
    call this%AtmosMesh_Final()

    return
  end subroutine AtmosMeshRM_Final

  subroutine AtmosMeshRM_create_communicator( this, sfield_num, hvfield_num, var_manager, field_list, commid )
    implicit none
    class(AtmosMeshRM), target, intent(inout) :: this
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    class(ModelVarManager), intent(inout) :: var_manager
    class(MeshField3D), intent(in) :: field_list(:)
    integer, intent(out) :: commid
    !-----------------------------------------------------

    commid = this%Get_communicatorID( ATM_MESH_MAX_COMMNUICATOR_NUM )
    call this%comm_list(commid)%Init(sfield_num,  hvfield_num, this%mesh )
    call var_manager%MeshFieldComm_Prepair( this%comm_list(commid), field_list )

    return
  end subroutine AtmosMeshRM_create_communicator 
  
  subroutine AtmosMeshRM_setup_restartfile1( this, restart_file, var_num )
    implicit none
    class(AtmosMeshRM), target, intent(inout) :: this
    class(FILE_restart_meshfield_component), intent(inout) :: restart_file
    integer, intent(in) :: var_num  
    !------------------------------------------------

    call restart_file%Init('ATMOS', var_num, mesh3D=this%mesh )
    return
  end subroutine AtmosMeshRM_setup_restartfile1

  subroutine AtmosMeshRM_setup_restartfile2( this, restart_file, &
    in_basename, in_postfix_timelabel,                         &
    out_basename, out_postfix_timelabel,                       &
    out_dtype, out_title, var_num                              )
    implicit none
    class(AtmosMeshRM), target, intent(inout) :: this
    class(FILE_restart_meshfield_component), intent(inout) :: restart_file
    character(*), intent(in) :: in_basename
    logical, intent(in) :: in_postfix_timelabel
    character(*), intent(in) :: out_basename
    logical, intent(in) :: out_postfix_timelabel
    character(*), intent(in) :: out_title
    character(*), intent(in) :: out_dtype  
    integer, intent(in) :: var_num  
    !-----------------------------------------------------------

    call restart_file%Init('ATMOS', in_basename, in_postfix_timelabel,    &
      out_basename, out_postfix_timelabel, out_dtype, out_title, var_num, &
      mesh3D=this%mesh )

    return
  end subroutine AtmosMeshRM_setup_restartfile2

  subroutine AtmosMeshRM_calc_UVMet( this, U, V, &
    Umet, Vmet )
    implicit none
    class(AtmosMeshRM), target, intent(in) :: this
    type(MeshField3D), intent(in) :: U
    type(MeshField3D), intent(in) :: V
    type(MeshField3D), intent(inout) :: Umet
    type(MeshField3D), intent(inout) :: Vmet

    integer :: n
    integer :: ke
    type(LocalMesh3D), pointer :: lcmesh
    !------------------------------------------

    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh%lcmesh_list(n)
      !$omp parallel do
      do ke=lcmesh%NeS, lcmesh%NeE
        Umet%local(n)%val(:,ke) = U%local(n)%val(:,ke)
        Vmet%local(n)%val(:,ke) = V%local(n)%val(:,ke)
      end do
    end do

    return
  end subroutine AtmosMeshRM_calc_UVMet

  subroutine AtmosMeshRM_setup_vcoordinate( this )
    use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D
    implicit none
    class(AtmosMeshRM), TARGET, INTENT(INOUT) :: this

    type(MeshFieldCommCubeDom3D) :: comm3D
    type(MeshFieldCommRectDom2D) :: comm2D
    !-------------------------------------------------

    call comm2D%Init( 1, 0, this%mesh%mesh2D )
    call comm3D%Init( 1, 1, this%mesh )

    call this%topography%SetVCoordinate( this%ptr_mesh,   &
      this%vcoord_type_id, this%mesh%zmax_gl, comm3D, comm2D )

    call comm2D%Final()
    call comm3D%Final()

    return
  end SUBROUTINE AtmosMeshRM_setup_vcoordinate

end module mod_atmos_mesh_rm