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

  subroutine AtmosMeshGM_Init( this )    
    use scale_const, only: &
      RPlanet => CONST_RADIUS
    
    implicit none
    class(AtmosMeshGM), target, intent(inout) :: this

    real(RP) :: dom_zmin          = 0.0_RP
    real(RP) :: dom_zmax          = 10.0E3_RP
    logical  :: isPeriodicZ       = .false.
  
    integer, parameter :: FZ_nmax = 1000
    real(RP) :: FZ(FZ_nmax)

    !* Global
    integer  :: NeGX               = 2
    integer  :: NeGY               = 2
    integer  :: NeZ                = 2
    integer  :: NLocalMeshPerPrc   = 6
    integer  :: Nprc               = 1
    integer  :: PolyOrder_h        = 2
    integer  :: PolyOrder_v        = 2
    logical  :: LumpedMassMatFlag  = .false.

    namelist / PARAM_ATMOS_MESH / &
      dom_zmin, dom_zmax,                          &
      FZ, isPeriodicZ,                             &
      NeGX, NeGY, NeZ, NLocalMeshPerPrc, Nprc,     &
      PolyOrder_h, PolyOrder_v, LumpedMassMatFlag    

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
        NeGX, NeGY, NeZ, RPlanet, dom_zmin, dom_zmax, &
        this%element, NLocalMeshPerPrc, nproc=Nprc,   &
        FZ=FZ(1:NeZ+1)    )
    else
      call this%mesh%Init( &
        NeGX, NeGY, NeZ, RPlanet, dom_zmin, dom_zmax, &
        this%element, NLocalMeshPerPrc, nproc=Nprc    )
    end if
    
    call this%mesh%Generate()
    
    !-
    call this%AtmosMesh_Init( this%mesh )

    return
  end subroutine AtmosMeshGM_Init

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

  subroutine AtmosMeshGM_create_communicator( this, sfield_num, hvfield_num, var_manager, field_list, commid )
    implicit none
    class(AtmosMeshGM), target, intent(inout) :: this
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    class(ModelVarManager), intent(inout) :: var_manager
    class(MeshField3D), intent(in) :: field_list(:)
    integer, intent(out) :: commid
    !-----------------------------------------------------

    commid = this%Get_communicatorID( ATM_MESH_MAX_COMMNUICATOR_NUM )
    call this%comm_list(commid)%Init( sfield_num, hvfield_num, this%mesh )
    call var_manager%MeshFieldComm_Prepair( this%comm_list(commid), field_list )

    return
  end subroutine AtmosMeshGM_create_communicator  

  subroutine AtmosMeshGM_setup_restartfile1( this, restart_file, var_num )
    implicit none
    class(AtmosMeshGM), target, intent(inout) :: this
    class(FILE_restart_meshfield_component), intent(inout) :: restart_file
    integer, intent(in) :: var_num  
    !------------------------------------------------

    call restart_file%Init('ATMOS', var_num, meshcubedsphere3D=this%mesh )
    return
  end subroutine AtmosMeshGM_setup_restartfile1

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

end module mod_atmos_mesh_gm