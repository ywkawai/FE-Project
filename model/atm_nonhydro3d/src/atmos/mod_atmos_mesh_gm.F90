!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_mesh_gm
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_meshfield_base, only: MeshField3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_meshfieldcomm_cubedspheredom3d, only: MeshFieldCommCubedSphereDom3D
  use scale_sparsemat, only: sparsemat
  use scale_file_base_meshfield, only: FILE_base_meshfield
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

  type, extends(AtmosMesh), public :: AtmosMeshGM
    type(MeshCubedSphereDom3D) :: mesh
    type(MeshFieldCommCubedSphereDom3D) :: comm_list(ATM_MESH_MAX_COMMNUICATOR_NUM)
  contains
    procedure :: Init => AtmosMeshGM_Init
    procedure :: Final => AtmosMeshGM_Final
    procedure :: Create_communicator => AtmosMeshGM_Create_communicator
    procedure :: Setup_restartfile1 => AtmosMeshGM_setup_restartfile1
    procedure :: Setup_restartfile2 => AtmosMeshGM_setup_restartfile2
    procedure :: Calc_UVmet => AtmosMeshGM_calc_UVMet
    procedure :: Setup_vcoordinate => AtmosMeshGM_setup_vcoordinate
  end type AtmosMeshGM

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
!OCL SERIAL
  subroutine AtmosMeshGM_Init( this )    
    use scale_const, only: &
      RPlanet => CONST_RADIUS
    use scale_mesh_base2d, only: &
      MFTYPE2D_XY => MeshBase2D_DIMTYPEID_XY
    use scale_meshutil_vcoord, only: &
      MeshUtil_get_VCoord_TypeID
    
    implicit none
    class(AtmosMeshGM), target, intent(inout) :: this

    real(RP) :: dom_zmin          = 0.0_RP
    real(RP) :: dom_zmax          = 10.0E3_RP
    logical  :: isPeriodicZ       = .false.
  
    integer, parameter :: FZ_nmax = 1000
    real(RP) :: FZ(FZ_nmax)

    !* Global
    logical :: SHALLOW_ATM_APPROX_FLAG = .true.
    integer  :: NeGX               = 2
    integer  :: NeGY               = 2
    integer  :: NeZ                = 2
    integer  :: NLocalMeshPerPrc   = 6
    integer  :: Nprc               = 1
    integer  :: PolyOrder_h        = 2
    integer  :: PolyOrder_v        = 2
    logical  :: LumpedMassMatFlag  = .false.
    character(len=H_MID)  :: VERTICAL_COORD_NAME = "TERRAIN_FOLLOWING"
    character(len=H_LONG) :: TOPO_IN_BASENAME    = ''       !< basename of the input file
    character(len=H_MID)  :: TOPO_IN_VARNAME     = 'topo'   !< variable name of topo in the input file
    logical :: COMM_USE_MPI_PC     = .false.
    
    namelist / PARAM_ATMOS_MESH / &
      SHALLOW_ATM_APPROX_FLAG,                     &
      dom_zmin, dom_zmax,                          &
      FZ, isPeriodicZ,                             &
      NeGX, NeGY, NeZ, NLocalMeshPerPrc, Nprc,     &
      PolyOrder_h, PolyOrder_v, LumpedMassMatFlag, &
      VERTICAL_COORD_NAME,                         &
      TOPO_IN_BASENAME, TOPO_IN_VARNAME,           &
      COMM_USE_MPI_PC
    
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
        NeGX, NeGY, NeZ, RPlanet, dom_zmin, dom_zmax,          &
        this%element, NLocalMeshPerPrc, nproc=Nprc,            &
        FZ=FZ(1:NeZ+1), shallow_approx=SHALLOW_ATM_APPROX_FLAG )
    else
      call this%mesh%Init( &
        NeGX, NeGY, NeZ, RPlanet, dom_zmin, dom_zmax, &
        this%element, NLocalMeshPerPrc, nproc=Nprc,   &
        shallow_approx=SHALLOW_ATM_APPROX_FLAG        )
    end if
    
    call this%mesh%Generate()
    
    !-

    call this%AtmosMesh_Init( this%mesh )

    !- Set topography & vertical coordinate
    
    if ( TOPO_IN_BASENAME /= '' ) then
      LOG_INFO("ATMOS_MESH_setup",*) 'Read topography data'

      call file_topo%Init(1, meshcubedsphere2D=this%mesh%mesh2D )
      call file_topo%Open( TOPO_IN_BASENAME, myrank=PRC_myrank )
      call file_topo%Read_Var( MFTYPE2D_XY, TOPO_IN_VARNAME, this%topography%topo )
      call file_topo%Close()
      call file_topo%Final()
    end if

    this%vcoord_type_id = MeshUtil_get_VCoord_TypeID( VERTICAL_COORD_NAME )
    call this%Setup_vcoordinate()

    !-
    this%comm_use_mpi_pc = COMM_USE_MPI_PC

    return
  end subroutine AtmosMeshGM_Init

!OCL SERIAL
  subroutine AtmosMeshGM_Final(this)
    implicit none

    class(AtmosMeshGM), intent(inout) :: this
    integer :: commid
    !-------------------------------------------

    do commid=1, this%communicator_num
      call this%comm_list(commid)%Final()
    end do

    call this%mesh%Final()
    call this%AtmosMesh_Final()
    
    return
  end subroutine AtmosMeshGM_Final

!OCL SERIAL
  subroutine AtmosMeshGM_create_communicator( this, sfield_num, hvfield_num, htensorfield_num, &
    var_manager, field_list, commid )
    implicit none
    class(AtmosMeshGM), target, intent(inout) :: this
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    integer, intent(in) :: htensorfield_num
    class(ModelVarManager), intent(inout) :: var_manager
    class(MeshField3D), intent(in) :: field_list(:)
    integer, intent(out) :: commid
    !-----------------------------------------------------

    commid = this%Get_communicatorID( ATM_MESH_MAX_COMMNUICATOR_NUM )
    call this%comm_list(commid)%Init( sfield_num, hvfield_num, htensorfield_num, this%mesh )
    if ( this%comm_use_mpi_pc ) call this%comm_list(commid)%Prepare_PC()
    call var_manager%MeshFieldComm_Prepair( this%comm_list(commid), field_list )

    return
  end subroutine AtmosMeshGM_create_communicator  

!OCL SERIAL
  subroutine AtmosMeshGM_setup_restartfile1( this, restart_file, var_num )
    implicit none
    class(AtmosMeshGM), target, intent(inout) :: this
    class(FILE_restart_meshfield_component), intent(inout) :: restart_file
    integer, intent(in) :: var_num  
    !------------------------------------------------

    call restart_file%Init('ATMOS', var_num, meshcubedsphere3D=this%mesh )
    return
  end subroutine AtmosMeshGM_setup_restartfile1

!OCL SERIAL
  subroutine AtmosMeshGM_setup_restartfile2( this, restart_file, &
    in_basename, in_postfix_timelabel,                         &
    out_basename, out_postfix_timelabel,                       &
    out_dtype, out_title, var_num                              )
    implicit none
    class(AtmosMeshGM), target, intent(inout) :: this
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
      meshcubedsphere3D=this%mesh )

  end subroutine AtmosMeshGM_setup_restartfile2

!OCL SERIAL
  subroutine AtmosMeshGM_calc_UVMet( this, U, V, &
    Umet, Vmet )

    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatVec
    implicit none
    class(AtmosMeshGM), target, intent(in) :: this
    type(MeshField3D), intent(in) :: U
    type(MeshField3D), intent(in) :: V
    type(MeshField3D), intent(inout) :: Umet
    type(MeshField3D), intent(inout) :: Vmet

    integer :: n
    integer :: ke, ke2D
    type(LocalMesh3D), pointer :: lcmesh
    class(ElementBase3D), pointer :: elem
    !------------------------------------------

    do n=1, this%mesh%LOCAL_MESH_NUM
      lcmesh => this%mesh%lcmesh_list(n)
      elem => lcmesh%refElem3D
      call CubedSphereCoordCnv_CS2LonLatVec( &
        lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), & ! (in)
        lcmesh%gam, elem%Np * lcmesh%Ne,                            & ! (in)
        U%local(n)%val(:,lcmesh%NeS:lcmesh%NeE),                    & ! (in)
        V%local(n)%val(:,lcmesh%NeS:lcmesh%NeE),                    & ! (in)
        Umet%local(n)%val(:,lcmesh%NeS:lcmesh%NeE),                 & ! (out)
        Vmet%local(n)%val(:,lcmesh%NeS:lcmesh%NeE)                  ) ! (out)
    end do

    return
  end subroutine AtmosMeshGM_calc_UVMet

!OCL SERIAL
  subroutine AtmosMeshGM_setup_vcoordinate( this )
    use scale_meshfieldcomm_cubedspheredom2d, only: MeshFieldCommCubedSphereDom2D
    implicit none
    class(AtmosMeshGM), TARGET, INTENT(INOUT) :: this

    type(MeshFieldCommCubedSphereDom3D) :: comm3D
    type(MeshFieldCommCubedSphereDom2D) :: comm2D
    !-------------------------------------------------

    call comm2D%Init( 1, 0, 0, this%mesh%mesh2D )
    call comm3D%Init( 2, 1, 0, this%mesh )

    call this%topography%SetVCoordinate( this%ptr_mesh,   &
      this%vcoord_type_id, this%mesh%zmax_gl, comm3D, comm2D )

    call comm2D%Final()
    call comm3D%Final()

    return
  end SUBROUTINE AtmosMeshGM_setup_vcoordinate

end module mod_atmos_mesh_gm