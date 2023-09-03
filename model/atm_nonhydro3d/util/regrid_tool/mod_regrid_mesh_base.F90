!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_regrid_mesh_base
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort
  use scale_const, only: &
    UNDEF => CONST_UNDEF,   &
    PI => CONST_PI,         &
    RPlanet => CONST_radius

  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D

  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom2d, only: &
    MeshCubedSphereDom2D,                &
    MeshCubedSphereDom2D_check_division_params
  use scale_mesh_cubedspheredom3d, only: &
    MeshCubedSphereDom3D
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: regrid_mesh_base
    integer :: NprcX
    integer :: NprcY
    integer :: NprcZ
    integer :: Nprc
    real(RP) :: dom_xmin
    real(RP) :: dom_xmax
    real(RP) :: dom_ymin
    real(RP) :: dom_ymax
    real(RP) :: dom_zmin
    real(RP) :: dom_zmax
    logical  :: isPeriodicX  = .false.
    logical  :: isPeriodicY  = .false.
    logical  :: isPeriodicZ  = .false.
    integer :: NLocalMeshPerPRC

    integer :: NeX
    integer :: NeY
    integer :: NeZ
    integer :: NeGX
    integer :: NeGY
    integer :: NeGZ

    integer :: polyorder_h
    integer :: polyorder_v

    type(QuadrilateralElement) :: elem2D
    type(HexahedralElement) :: elem3D
    
    type(MeshRectDom2D) :: mesh2D
    type(MeshCubeDom3D) :: mesh3D
    type(MeshCubedSphereDom2D) :: csmesh2D
    type(MeshCubedSphereDom3D) :: csmesh3D

    class(MeshBase2D), pointer :: ptr_mesh2D
    class(MeshBase3D), pointer :: ptr_mesh3D

    integer :: inout_id
    integer :: mesh_type_id

    real(RP), allocatable :: FZ(:)
  contains
    procedure :: Init1 => regrid_mesh_base_init_1
    procedure :: Init2 => regrid_mesh_base_init_2
    generic :: Init => Init1, Init2
    procedure :: Final => regrid_mesh_base_final
    procedure :: Generate => regrid_mesh_base_generate
    procedure :: Get_inmesh_hmapinfo => regrid_mesh_base_get_inmesh_hmapinfo
  end type regrid_mesh_base

  public :: regrid_mesh_base_meshtype_name2id

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
    
  integer, public, parameter :: REGRID_MESH_BASE_IN_ID  = 1
  integer, public, parameter :: REGRID_MESH_BASE_OUT_ID = 2

  integer, public, parameter :: REGRID_MESHTYPE_STRUCTURED2D_ID  = 1
  integer, public, parameter :: REGRID_MESHTYPE_CUBEDSPHERE2D_ID = 2
  integer, public, parameter :: REGRID_MESHTYPE_LONLAT2D_ID      = 3
  integer, public, parameter :: REGRID_MESHTYPE_STRUCTURED3D_ID  = 4
  integer, public, parameter :: REGRID_MESHTYPE_CUBEDSPHERE3D_ID = 5
  integer, public, parameter :: REGRID_MESHTYPE_LONLAT3D_ID      = 6

contains

!OCL SERIAL
  subroutine regrid_mesh_base_init_1( this, & ! (inout)
    mesh_inout_id, mesh_type_id, GP_flag    ) ! (in)

    implicit none
    class(regrid_mesh_base), intent(inout) :: this
    integer, intent(in) :: mesh_inout_id
    integer, intent(in) :: mesh_type_id
    logical, intent(in) :: GP_flag
    !----------------------------------------

    select case( mesh_inout_id )
    case ( REGRID_MESH_BASE_IN_ID, REGRID_MESH_BASE_OUT_ID )
      this%inout_id = mesh_inout_id
    case default
      LOG_ERROR("regrid_mesh_base",*) 'Not appropriate inout id of mesh. Check!'
      call PRC_abort
    end select

    this%mesh_type_id = mesh_type_id

    select case (this%mesh_type_id)
    case( REGRID_MESHTYPE_STRUCTURED2D_ID )
      call regrid_mesh_base_init_mesh2D( this, GP_flag )
    case( REGRID_MESHTYPE_STRUCTURED3D_ID )
      call regrid_mesh_base_init_mesh3D( this, GP_flag )
    case( REGRID_MESHTYPE_LONLAT2D_ID )
      call regrid_mesh_base_init_mesh2D( this, GP_flag )
    case( REGRID_MESHTYPE_LONLAT3D_ID )
      call regrid_mesh_base_init_mesh3D( this, GP_flag )
    case( REGRID_MESHTYPE_CUBEDSPHERE2D_ID )
      call regrid_mesh_base_init_csmesh2D( this, GP_flag )
    case( REGRID_MESHTYPE_CUBEDSPHERE3D_ID )
      call regrid_mesh_base_init_csmesh3D( this, GP_flag )
    case default
      LOG_ERROR("regrid_mesh_base_init",*) 'Not supported type id of mesh. Check! ', mesh_type_id
      call PRC_abort
    end select

    nullify( this%ptr_mesh2D, this%ptr_mesh3D )

    return
  end subroutine regrid_mesh_base_init_1

!OCL SERIAL
  subroutine regrid_mesh_base_init_2( this, & ! (inout)
    mesh_inout_id, mesh_type_name, GP_flag  ) ! (in)

    implicit none
    class(regrid_mesh_base), intent(inout) :: this
    integer, intent(in) :: mesh_inout_id
    character(len=*), intent(in) :: mesh_type_name
    logical, intent(in) :: GP_flag

    integer :: mesh_type_id
    !----------------------------------------
    
    mesh_type_id = regrid_mesh_base_meshtype_name2id( mesh_type_name )
    call this%Init1( mesh_inout_id, mesh_type_id, GP_flag )

    return
  end subroutine regrid_mesh_base_init_2

!OCL SERIAL
  function regrid_mesh_base_meshtype_name2id( mesh_type_name ) result(mesh_type_id)
    implicit none
    character(len=*), intent(in) :: mesh_type_name
    integer :: mesh_type_id
    !----------------------------------------
    
    select case(trim(mesh_type_name))
    case( "STRUCTURED2D" )
      mesh_type_id = REGRID_MESHTYPE_STRUCTURED2D_ID
    case( "STRUCTURED3D" )
      mesh_type_id = REGRID_MESHTYPE_STRUCTURED3D_ID
    case( "LONLAT2D" )
      mesh_type_id = REGRID_MESHTYPE_LONLAT2D_ID
    case( "LONLAT3D" )
      mesh_type_id = REGRID_MESHTYPE_LONLAT3D_ID
    case( "CUBEDSPHERE2D" )
      mesh_type_id = REGRID_MESHTYPE_CUBEDSPHERE2D_ID
    case( "CUBEDSPHERE3D" )
      mesh_type_id = REGRID_MESHTYPE_CUBEDSPHERE3D_ID
    case default
      LOG_ERROR("regrid_mesh_meshtype_name2id",*) 'Not supported type of mesh. Check! ', trim(mesh_type_name)
      call PRC_abort
    end select

    return    
  end function regrid_mesh_base_meshtype_name2id

!OCL SERIAL  
  subroutine regrid_mesh_base_final( this )
    implicit none
    class(regrid_mesh_base), intent(inout) :: this
    !----------------------------------------

    select case( this%mesh_type_id )
    case( REGRID_MESHTYPE_STRUCTURED2D_ID, REGRID_MESHTYPE_LONLAT2D_ID )
      call this%mesh2D%Final()
    case( REGRID_MESHTYPE_STRUCTURED3D_ID, REGRID_MESHTYPE_LONLAT3D_ID )
      call this%mesh3D%Final()
    case( REGRID_MESHTYPE_CUBEDSPHERE2D_ID )
      call this%csmesh2D%Final()
    case( REGRID_MESHTYPE_CUBEDSPHERE3D_ID )
      call this%csmesh3D%Final()
    end select
    
    if ( allocated(this%FZ) ) then
      deallocate( this%FZ )
    end if

    return
  end subroutine regrid_mesh_base_final

!OCL SERIAL  
  subroutine regrid_mesh_base_generate( this, myrank )
    use scale_prc, only: &
      PRC_myrank 
    
    use scale_mesh_base2d, only: &
      MeshBase2D_DIMTYPEID_X, MeshBase2D_DIMTYPEID_Y,   &
      MeshBase2D_DIMTYPEID_XY, MeshBase2D_DIMTYPEID_XYT
    use scale_mesh_base3d, only: &
      MeshBase3D_DIMTYPEID_X, MeshBase3D_DIMTYPEID_Y, MeshBase3D_DIMTYPEID_Z,  &
      MeshBase3D_DIMTYPEID_XYZ, MeshBase3D_DIMTYPEID_XYZT

    implicit none
    class(regrid_mesh_base), intent(inout), target :: this
    integer, intent(in), optional :: myrank

    integer :: myrank_
    !-----------------------------------------------

    if ( present(myrank) ) then
      myrank_ = myrank
    else
      myrank_ = PRC_myrank
    end if

    select case( this%mesh_type_id )
    case( REGRID_MESHTYPE_STRUCTURED2D_ID )

      call this%mesh2D%Init( this%NeGX, this%NeGY, &
        this%dom_xmin, this%dom_xmax, this%dom_ymin, this%dom_ymax,             &
        this%isPeriodicX, this%isPeriodicY, this%elem2D, this%NLocalMeshPerPrc, &
        this%NprcX, this%NprcY, myrank=myrank_  )
      
      this%ptr_mesh2D => this%mesh2D

    case( REGRID_MESHTYPE_LONLAT2D_ID )

      call this%mesh2D%Init( this%NeGX, this%NeGY, &
        this%dom_xmin, this%dom_xmax, this%dom_ymin, this%dom_ymax,             &
        this%isPeriodicX, .false., this%elem2D, this%NLocalMeshPerPrc,          &
        this%NprcX, this%NprcY, myrank=myrank_ )

      call this%mesh2D%SetDimInfo( MeshBase2D_DIMTYPEID_X, 'lon', 'degree_east', 'longitude' )
      call this%mesh2D%SetDimInfo( MeshBase2D_DIMTYPEID_Y, 'lat', 'degree_north', 'latitude' )
      call this%mesh2D%SetDimInfo( MeshBase2D_DIMTYPEID_XY, 'lonlat', 'degree', 'longitude,latitude' )
      call this%mesh2D%SetDimInfo( MeshBase2D_DIMTYPEID_XYT, 'lonlatt', 'degree', 'longitude,latitude' )

      this%ptr_mesh2D => this%mesh2D

    case( REGRID_MESHTYPE_STRUCTURED3D_ID )

      if ( allocated( this%FZ ) ) then
        call this%mesh3D%Init( this%NeGX, this%NeGY, this%NeGZ,                                     &
          this%dom_xmin, this%dom_xmax, this%dom_ymin, this%dom_ymax, this%dom_zmin, this%dom_zmax, &
          this%isPeriodicX, this%isPeriodicY, .false., this%elem3D, this%NLocalMeshPerPrc,          &
          this%NprcX, this%NprcY, nproc=this%Nprc, myrank=myrank_, FZ=this%FZ )
      else
        call this%mesh3D%Init( this%NeGX, this%NeGY, this%NeGZ,                                     &
          this%dom_xmin, this%dom_xmax, this%dom_ymin, this%dom_ymax, this%dom_zmin, this%dom_zmax, &
          this%isPeriodicX, this%isPeriodicY, .false., this%elem3D, this%NLocalMeshPerPrc,          &
          this%NprcX, this%NprcY, nproc=this%Nprc, myrank=myrank_ )
      end if

      this%ptr_mesh3D => this%mesh3D

    case( REGRID_MESHTYPE_LONLAT3D_ID )

      if ( allocated( this%FZ ) ) then
        call this%mesh3D%Init( this%NeGX, this%NeGY, this%NeGZ,                                     &
          this%dom_xmin, this%dom_xmax, this%dom_ymin, this%dom_ymax, this%dom_zmin, this%dom_zmax, &
          this%isPeriodicX, .false., .false., this%elem3D, this%NLocalMeshPerPrc,                   &
          this%NprcX, this%NprcY, nproc=this%Nprc, myrank=myrank_, FZ=this%FZ )
      else
        call this%mesh3D%Init( this%NeGX, this%NeGY, this%NeGZ,                                     &
          this%dom_xmin, this%dom_xmax, this%dom_ymin, this%dom_ymax, this%dom_zmin, this%dom_zmax, &
          this%isPeriodicX, .false., .false., this%elem3D, this%NLocalMeshPerPrc,                   &
          this%NprcX, this%NprcY, nproc=this%Nprc, myrank=myrank_ )
      end if
      call this%mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_X, 'lon', 'degree_east', 'longitude' )
      call this%mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_Y, 'lat', 'degree_north', 'latitude' )
      call this%mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_XYZ, 'lonlatz', '', 'longitude,latitude,altitude' )
      call this%mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_XYZT, 'lonlatzt', '', 'longitude,latitude,altitude' )
      
      this%ptr_mesh3D => this%mesh3D

    case( REGRID_MESHTYPE_CUBEDSPHERE2D_ID )
      
      call this%csmesh2D%Init( this%NeGX, this%NeGY, RPlanet,               &
        this%elem2D, this%NLocalMeshPerPrc, nproc=this%Nprc, myrank=myrank_ )
      
      this%ptr_mesh2D => this%csmesh2D

    case( REGRID_MESHTYPE_CUBEDSPHERE3D_ID )

      if ( allocated( this%FZ ) ) then
        call this%csmesh3D%Init( this%NeGX, this%NeGY, this%NeGZ, RPlanet,  &
          this%dom_zmin, this%dom_zmax, this%elem3D, this%NLocalMeshPerPrc, &
          nproc=this%Nprc, myrank=myrank_, FZ=this%FZ )
      else
        call this%csmesh3D%Init( this%NeGX, this%NeGY, this%NeGZ, RPlanet,  &
          this%dom_zmin, this%dom_zmax, this%elem3D, this%NLocalMeshPerPrc, &
          nproc=this%Nprc, myrank=myrank_ )        
      end if
      
      this%ptr_mesh3D => this%csmesh3D

    end select

    if ( associated( this%ptr_mesh2D ) ) call this%ptr_mesh2D%Generate()
    if ( associated( this%ptr_mesh3D ) ) call this%ptr_mesh3D%Generate()

    return
  end subroutine regrid_mesh_base_generate

  subroutine regrid_mesh_base_get_inmesh_hmapinfo( this, &
    tileID_table, panelID_table, pi_table, pj_table       )

    implicit none
    class(regrid_mesh_base), intent(inout) :: this
    integer, intent(out)  :: tileID_table(this%NLocalMeshPerPrc,this%Nprc)
    integer, intent(out)  :: panelID_table(this%NLocalMeshPerPrc,this%Nprc) 
    integer, intent(out) :: pi_table(this%NLocalMeshPerPrc*this%Nprc)
    integer, intent(out) :: pj_table(this%NLocalMeshPerPrc*this%Nprc)
    
    type(MeshCubedSphereDom2D) :: csmesh_dummy
    type(MeshRectDom2D) :: mesh_dummy
    !----------------------------------------

    select case( this%mesh_type_id )
    case ( REGRID_MESHTYPE_CUBEDSPHERE2D_ID, REGRID_MESHTYPE_CUBEDSPHERE3D_ID )

      call csmesh_dummy%Init( this%NeGX, this%NeGY, RPlanet,  & ! (in)
        this%elem2D, this%NLocalMeshPerPRC,                   & ! (in)
        nproc=this%Nprc, myrank=0                             ) ! (in)
          
      call csmesh_dummy%AssignDomID( this%NprcX, this%NprcY,  & ! (in)
        tileID_table, panelID_table, pi_table, pj_table       ) ! (out)
      
      call csmesh_dummy%Final()

    case ( REGRID_MESHTYPE_STRUCTURED2D_ID, REGRID_MESHTYPE_STRUCTURED3D_ID, &
           REGRID_MESHTYPE_LONLAT2D_ID, REGRID_MESHTYPE_LONLAT3D_ID          )

      call mesh_dummy%Init( this%NeGX, this%NeGY,                   & ! (in)
        this%dom_xmin, this%dom_xmax, this%dom_ymin, this%dom_ymax, & ! (in)
        this%isPeriodicX, this%isPeriodicY, this%elem2D,            & ! (in)
        this%NLocalMeshPerPRC, this%NprcX, this%NprcY,              & ! (in)
        nproc=this%Nprc, myrank=0                                   ) ! (in)

      call mesh_dummy%AssignDomID( &
        tileID_table, panelID_table, pi_table, pj_table )
      
      call mesh_dummy%Final()

    end select

    return
  end subroutine regrid_mesh_base_get_inmesh_hmapinfo

!-- private --------------------------------

!OCL SERIAL  
  subroutine regrid_mesh_base_init_mesh2D( this, GP_flag )

    implicit none
    class(regrid_mesh_base), intent(inout) :: this
    logical, intent(in) :: GP_flag

    ! Structured mesh  
    integer :: NprcX         = 1    
    integer :: NprcY         = 1       
    integer :: NeX           = 1
    integer :: NeY           = 1
    real(RP) :: dom_xmin     = 0.0_RP
    real(RP) :: dom_xmax     = 0.0_RP    
    real(RP) :: dom_ymin     = 0.0_RP
    real(RP) :: dom_ymax     = 0.0_RP    
    logical  :: isPeriodicX  = .false.
    logical  :: isPeriodicY  = .false.
    integer :: NLocalMeshPerPrc = 1

    integer :: PolyOrder_h     = 1
    integer :: PolyOrder_h_GP  = 1

    namelist / PARAM_REGRID_INMESH2D_STRUCTURED / &
      NprcX, NprcY, NeX, NeY, NLocalMeshPerPrc,   &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax,     &
      PolyOrder_h,                                &
      isPeriodicX, isPeriodicY
    
    namelist / PARAM_REGRID_OUTMESH2D_STRUCTURED / &
      NprcX, NprcY, NeX, NeY, NLocalMeshPerPrc,    &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax,      &
      PolyOrder_h, PolyOrder_h_GP,                 &
      isPeriodicX, isPeriodicY

    integer :: ierr
    !----------------------------------------

    rewind(IO_FID_CONF)

    if (this%inout_id == REGRID_MESH_BASE_IN_ID) then

      read(IO_FID_CONF,nml=PARAM_REGRID_INMESH2D_STRUCTURED,iostat=ierr)
      if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_mesh_base_init2D",*) 'Not found namelist. Default used.'
      elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_mesh_base_init2D",*) 'Not appropriate names in namelist PARAM_REGRID_INMESH2D_STRUCTURED. Check!'
        call PRC_abort
      endif
      LOG_NML(PARAM_REGRID_INMESH2D_STRUCTURED)

    else if (this%inout_id == REGRID_MESH_BASE_OUT_ID ) then

      read(IO_FID_CONF,nml=PARAM_REGRID_OUTMESH2D_STRUCTURED,iostat=ierr)
      if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_mesh_base_init2D",*) 'Not found namelist. Default used.'
      elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_mesh_base_init2D",*) 'Not appropriate names in namelist PARAM_REGRID_OUTMESH2D_STRUCTURED. Check!'
        call PRC_abort
      endif
      LOG_NML(PARAM_REGRID_OUTMESH2D_STRUCTURED)

    end if

    this%NprcX = NprcX
    this%NprcY = NprcY
    this%NprcZ = 1
    this%Nprc = NprcX * NprcY

    this%NeX = NeX
    this%NeY = NeY
    this%NeZ = 1
    this%NeGX = NeX * NprcX
    this%NeGY = NeY * NprcY
    this%NeGZ = 1
    this%NLocalMeshPerPRC = NLocalMeshPerPrc
    
    this%dom_xmin = dom_xmin
    this%dom_xmax = dom_xmax
    this%dom_ymin = dom_ymin
    this%dom_ymax = dom_ymax

    this%polyorder_h = PolyOrder_h
    this%polyorder_v = -1
    if ( this%inout_id == REGRID_MESH_BASE_OUT_ID .and. GP_flag ) then
      this%polyorder_h = PolyOrder_h_GP
    end if

    this%isPeriodicX = isPeriodicX
    this%isPeriodicY = isPeriodicY

    call this%elem2D%Init( PolyOrder_h, .true. )

    return
  end subroutine regrid_mesh_base_init_mesh2D

!OCL SERIAL  
  subroutine regrid_mesh_base_init_csmesh2D( this, GP_flag )
    implicit none

    class(regrid_mesh_base), intent(inout) :: this
    logical, intent(in) :: GP_flag

    ! Structured mesh  
    integer :: Nprc             = 1    
    integer :: NeGX             = 1
    integer :: NeGY             = 1
    integer :: NLocalMeshPerPrc = 1

    integer :: PolyOrder_h      = 1
    integer :: PolyOrder_h_GP   = 1

    namelist / PARAM_REGRID_INMESH2D_CUBEDSPHERE / &
      Nprc, NeGX, NeGY, NLocalMeshPerPrc,          &
      PolyOrder_h
   
    namelist / PARAM_REGRID_OUTMESH2D_CUBEDSPHERE / &
      Nprc, NeGX, NeGY, NLocalMeshPerPrc,          &
      PolyOrder_h, PolyOrder_h_GP

    integer :: ierr

    type(MeshCubedSphereDom2D) :: csmesh2D_dummy
    !----------------------------------------

    rewind(IO_FID_CONF)

    if (this%inout_id == REGRID_MESH_BASE_IN_ID) then

      read(IO_FID_CONF,nml=PARAM_REGRID_INMESH2D_CUBEDSPHERE,iostat=ierr)
      if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_mesh_base_init2D",*) 'Not found namelist. Default used.'
      elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_mesh_base_init2D",*) 'Not appropriate names in namelist PARAM_REGRID_INMESH2D_CUBEDSPHERE. Check!'
        call PRC_abort
      endif
      LOG_NML(PARAM_REGRID_INMESH2D_CUBEDSPHERE)

    else if (this%inout_id == REGRID_MESH_BASE_OUT_ID ) then

      read(IO_FID_CONF,nml=PARAM_REGRID_OUTMESH2D_CUBEDSPHERE,iostat=ierr)
      if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_mesh_base_init2D",*) 'Not found namelist. Default used.'
      elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_mesh_base_init2D",*) 'Not appropriate names in namelist PARAM_REGRID_OUTMESH2D_CUBEDSPHERE. Check!'
        call PRC_abort
      endif
      LOG_NML(PARAM_REGRID_OUTMESH2D_CUBEDSPHERE)

    end if

    this%Nprc = Nprc

    this%NeGX = NeGX
    this%NeGY = NeGY
    this%NeGZ = 1    
    this%NLocalMeshPerPRC = NLocalMeshPerPrc

    this%polyorder_h = PolyOrder_h
    this%polyorder_v = -1
    if ( this%inout_id == REGRID_MESH_BASE_OUT_ID .and. GP_flag ) then
      this%polyorder_h = PolyOrder_h_GP
    end if

    call this%elem2D%Init( PolyOrder_h, .true. )

    call csmesh2D_dummy%Init( NeGX, NeGY, RPlanet,        &
      this%elem2D, NLocalMeshPerPrc, nproc=Nprc, myrank=0 )
   
    call MeshCubedSphereDom2D_check_division_params(  &
      this%NprcX, this%NprcY,                         & ! (out)
      Nprc, NLocalMeshPerPrc * Nprc                   ) ! (in)
    
    this%NeX = NeGX / this%NprcX
    this%NeY = NeGY / this%NprcY
    this%NeZ = 1
    this%dom_xmin = csmesh2D_dummy%xmin_gl
    this%dom_xmax = csmesh2D_dummy%xmax_gl
    this%dom_ymin = csmesh2D_dummy%ymin_gl
    this%dom_ymax = csmesh2D_dummy%ymax_gl
  
   call csmesh2D_dummy%Final()

   return
  end subroutine regrid_mesh_base_init_csmesh2D

!OCL SERIAL  
  subroutine regrid_mesh_base_init_mesh3D( this, GP_flag )
    
    implicit none

    class(regrid_mesh_base), intent(inout) :: this
    logical, intent(in) :: GP_flag

    ! Structured mesh  
    integer :: NprcX         = 1    
    integer :: NprcY         = 1       
    integer :: NprcZ         = 1      
    integer :: NeX           = 1
    integer :: NeY           = 1
    integer :: NeGZ          = 1
    real(RP) :: dom_xmin     = 0.0_RP
    real(RP) :: dom_xmax     = 0.0_RP    
    real(RP) :: dom_ymin     = 0.0_RP
    real(RP) :: dom_ymax     = 0.0_RP    
    real(RP) :: dom_zmin     = 0.0_RP
    real(RP) :: dom_zmax     = 0.0_RP    
    logical  :: isPeriodicX  = .false.
    logical  :: isPeriodicY  = .false.
    logical  :: isPeriodicZ  = .false.
    integer :: NLocalMeshPerPrc = 1

    integer :: PolyOrder_h     = 1
    integer :: PolyOrder_h_GP  = 1
    integer :: PolyOrder_v     = 1
    integer :: PolyOrder_v_GP  = 1
    
    logical :: is_spec_FZ         
    integer, parameter :: FZ_nmax = 1000
    real(RP) :: FZ(FZ_nmax)

    namelist / PARAM_REGRID_INMESH3D_STRUCTURED / &
      NprcX, NprcY, NprcZ, NeX, NeY, NeGZ, NLocalMeshPerPrc,           &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,      &
      PolyOrder_h, PolyOrder_v,                                        &
      isPeriodicX, isPeriodicY, isPeriodicZ, &
      FZ
   
    namelist / PARAM_REGRID_OUTMESH3D_STRUCTURED / &
      NprcX, NprcY, NprcZ, NeX, NeY, NeGZ, NLocalMeshPerPrc,           &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,      &
      PolyOrder_h, PolyOrder_v,                                        &
      PolyOrder_h_GP, PolyOrder_v_GP,                                  &
      isPeriodicX, isPeriodicY, isPeriodicZ, &
      FZ

    integer :: ierr
    integer :: k
    !----------------------------------------

    NeGZ       = -1
    FZ(:)      = UNDEF

    rewind(IO_FID_CONF)

    if (this%inout_id == REGRID_MESH_BASE_IN_ID) then
 
      read(IO_FID_CONF,nml=PARAM_REGRID_INMESH3D_STRUCTURED,iostat=ierr)
      if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_mesh_base_init3D",*) 'Not found namelist. Default used.'
      elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_mesh_base_init3D",*) 'Not appropriate names in namelist PARAM_REGRID_INMESH3D_STRUCTURED. Check!'
        call PRC_abort
      endif
      LOG_NML(PARAM_REGRID_INMESH3D_STRUCTURED)

    else if (this%inout_id == REGRID_MESH_BASE_OUT_ID ) then

      read(IO_FID_CONF,nml=PARAM_REGRID_OUTMESH3D_STRUCTURED,iostat=ierr)
      if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_mesh_base_init3D",*) 'Not found namelist. Default used.'
      elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_mesh_base_init3D",*) 'Not appropriate names in namelist PARAM_REGRID_OUTMESH3D_STRUCTURED. Check!'
        call PRC_abort
      endif
      LOG_NML(PARAM_REGRID_OUTMESH3D_STRUCTURED)

    end if

    this%NprcX = NprcX
    this%NprcY = NprcY
    this%NprcZ = NprcZ
    this%Nprc = NprcX * NprcY * NprcZ

    this%NeX = NeX
    this%NeY = NeY
    this%NeZ = NeGZ / NprcZ
    this%NeGX = NeX * NprcX
    this%NeGY = NeY * NprcY
    this%NeGZ = NeGZ
    this%NLocalMeshPerPRC = NLocalMeshPerPrc
   
    this%dom_xmin = dom_xmin
    this%dom_xmax = dom_xmax
    this%dom_ymin = dom_ymin
    this%dom_ymax = dom_ymax
    this%dom_zmin = dom_zmin
    this%dom_zmax = dom_zmax

    this%polyorder_h = PolyOrder_h
    this%polyorder_v = PolyOrder_v
    if ( this%inout_id == REGRID_MESH_BASE_OUT_ID .and. GP_flag ) then
      this%polyorder_h = PolyOrder_h_GP
      this%polyorder_v = PolyOrder_v_GP
    end if

    this%isPeriodicX = isPeriodicX
    this%isPeriodicY = isPeriodicY
    this%isPeriodicZ = isPeriodicZ

    call this%elem2D%Init( PolyOrder_h, .true. )
    call this%elem3D%Init( PolyOrder_h, PolyOrder_v, .true. )
   
    is_spec_FZ = .true.    
    do k=1, NeGZ+1
      if ( FZ(k) == UNDEF ) then
        is_spec_FZ = .false.; exit
      end if
    end do
    if ( is_spec_FZ ) then
      allocate( this%FZ(1:NeGZ+1) )
      this%FZ(:) = FZ(1:NeGZ+1)
    end if
    
    return
  end subroutine regrid_mesh_base_init_mesh3D

!OCL SERIAL  
  subroutine regrid_mesh_base_init_csmesh3D( this, GP_flag )
    implicit none

    class(regrid_mesh_base), intent(inout) :: this
    logical, intent(in) :: GP_flag

    ! Structured mesh
    logical :: SHALLOW_ATM_APPROX_FLAG = .true.
    integer :: Nprc             = 1    
    integer :: NeGX             = 1
    integer :: NeGY             = 1
    integer :: NeGZ             = 1
    integer :: NLocalMeshPerPrc = 1
    real(RP) :: dom_zmin        = 0.0_RP
    real(RP) :: dom_zmax        = 0.0_RP    

    integer :: PolyOrder_h      = 1
    integer :: PolyOrder_h_GP   = 1
    integer :: PolyOrder_v      = 1
    integer :: PolyOrder_v_GP   = 1

    logical :: is_spec_FZ
    integer, parameter :: FZ_nmax = 1000
    real(RP) :: FZ(FZ_nmax)

    namelist / PARAM_REGRID_INMESH3D_CUBEDSPHERE / &
      SHALLOW_ATM_APPROX_FLAG,                     &
      Nprc, NeGX, NeGY, NeGZ, NLocalMeshPerPrc,    &
      PolyOrder_h, PolyOrder_v,                    &
      dom_zmin, dom_zmax,                          &
      FZ
   
    namelist / PARAM_REGRID_OUTMESH3D_CUBEDSPHERE / &
      SHALLOW_ATM_APPROX_FLAG,                      &
      Nprc, NeGX, NeGY, NeGZ, NLocalMeshPerPrc,     &
      PolyOrder_h, PolyOrder_v,                     &
      PolyOrder_h_GP, PolyOrder_v_GP,               &
      dom_zmin, dom_zmax,                           &
      FZ

    integer :: ierr
    integer :: k

    type(MeshCubedSphereDom2D) :: csmesh2D_dummy
    !----------------------------------------

    NeGZ       = -1
    FZ(:)      = UNDEF

    rewind(IO_FID_CONF)

    if (this%inout_id == REGRID_MESH_BASE_IN_ID) then

      read(IO_FID_CONF,nml=PARAM_REGRID_INMESH3D_CUBEDSPHERE,iostat=ierr)
      if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_mesh_base_init3D",*) 'Not found namelist. Default used.'
      elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_mesh_base_init3D",*) 'Not appropriate names in namelist PARAM_REGRID_INMESH3D_CUBEDSPHERE. Check!'
        call PRC_abort
      endif
      LOG_NML(PARAM_REGRID_INMESH3D_CUBEDSPHERE)

    else if (this%inout_id == REGRID_MESH_BASE_OUT_ID ) then

      read(IO_FID_CONF,nml=PARAM_REGRID_OUTMESH3D_CUBEDSPHERE,iostat=ierr)
      if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_mesh_base_init3D",*) 'Not found namelist. Default used.'
      elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_mesh_base_init3D",*) 'Not appropriate names in namelist PARAM_REGRID_OUTMESH3D_CUBEDSPHERE. Check!'
        call PRC_abort
      endif
      LOG_NML(PARAM_REGRID_OUTMESH3D_CUBEDSPHERE)

    end if

    this%Nprc  = Nprc
    this%NprcZ = 1

    this%NeGX = NeGX
    this%NeGY = NeGY
    this%NeGZ = NeGZ        
    this%NLocalMeshPerPRC = NLocalMeshPerPrc

    this%polyorder_h = PolyOrder_h
    this%polyorder_v = PolyOrder_v
    if ( this%inout_id == REGRID_MESH_BASE_OUT_ID .and. GP_flag ) then
      this%polyorder_h = PolyOrder_h_GP
      this%polyorder_v = PolyOrder_v_GP
    end if

    call this%elem2D%Init( PolyOrder_h, .true. )
    call this%elem3D%Init( PolyOrder_h, PolyOrder_v, .true. )

    call csmesh2D_dummy%Init( NeGX, NeGY, RPlanet,        &
      this%elem2D, NLocalMeshPerPrc, nproc=Nprc, myrank=0 )
   
    call MeshCubedSphereDom2D_check_division_params(  &
      this%NprcX, this%NprcY,                         & ! (out)
      Nprc, NLocalMeshPerPrc * Nprc                   ) ! (in)
    
    this%NeX = NeGX / this%NprcX
    this%NeY = NeGY / this%NprcY
    this%NeZ = NeGZ / this%NprcZ

    this%dom_xmin = csmesh2D_dummy%xmin_gl
    this%dom_xmax = csmesh2D_dummy%xmax_gl
    this%dom_ymin = csmesh2D_dummy%ymin_gl
    this%dom_ymax = csmesh2D_dummy%ymax_gl
    this%dom_zmin = dom_zmin
    this%dom_zmax = dom_zmax
         
    call csmesh2D_dummy%Final()

    is_spec_FZ = .true.    
    do k=1, NeGZ+1
      if ( FZ(k) == UNDEF ) then
        is_spec_FZ = .false.; exit
      end if
    end do
    if ( is_spec_FZ ) then
      allocate( this%FZ(1:NeGZ+1) )
      this%FZ(:) = FZ(1:NeGZ+1)
    end if

    return
  end subroutine regrid_mesh_base_init_csmesh3D

end module mod_regrid_mesh_base
