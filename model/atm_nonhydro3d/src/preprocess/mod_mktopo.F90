!-------------------------------------------------------------------------------
!> module INITIAL
!!
!! @par Description
!!          subroutines for preparing topography data 
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_mktopo
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc

  use scale_const, only: &
    PI => CONST_PI,          &
    RPlanet => CONST_RADIUS

  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_meshfield_base, only: MeshField2D
  use scale_mesh_topography, only: MeshTopography

  use mod_atmos_component, only: &
    AtmosComponent
  use mod_atmos_mesh, only: AtmosMesh

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKTOPO_setup
  public :: MKTOPO
  public :: MKTOPO_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  integer, public :: MKTOPO_TYPE                   = -1
  integer, parameter, public :: I_IGNORE           = 0
  integer, parameter, public :: I_FLAT             = 1
  integer, parameter, public :: I_BELLSHAPE        = 2
  integer, parameter, public :: I_BELLSHAPE_GLOBAL = 3

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  private :: MKTOPO_flat
  private :: MKTOPO_bellshape_global

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  character(len=H_LONG)  :: OUT_BASENAME = ''     !< basename of the output file
  character(len=H_MID)   :: OUT_VARNAME  = 'topo' !< variable name of topo in the output file
  character(len=H_MID)   :: OUT_TITLE    = 'SCALE-DG TOPOGRAPHY'  !< title    of the output file
  character(len=H_SHORT) :: OUT_DTYPE    = 'DEFAULT'              !< REAL4 or REAL8
  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKTOPO_setup
    implicit none

    character(len=H_SHORT) :: toponame = 'NONE'

    namelist / PARAM_MKTOPO /  &
      toponame,                  &
      OUT_BASENAME, OUT_VARNAME, &
      OUT_TITLE, OUT_DTYPE


    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKTOPO_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("MKTOPO_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("MKTOPO_setup",*) 'Not appropriate names in namelist PARAM_MKTOPO. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO)

    select case(trim(toponame))
    case('NONE')
      MKTOPO_TYPE = I_IGNORE
    case('FLAT')
      MKTOPO_TYPE = I_FLAT  
    case('BELLSHAPE')
      MKTOPO_TYPE = I_BELLSHAPE
    case('BELLSHAPE_GLOBAL')
      MKTOPO_TYPE = I_BELLSHAPE_GLOBAL 
    case default
      LOG_ERROR("MKTOPO_setup",*) 'Not appropriate toponame. Check!', toponame
      call PRC_abort      
    end select
    
    return
  end subroutine MKTOPO_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine MKTOPO( output, model_mesh, topography )
  
    use scale_model_var_manager, only: ModelVarManager

    implicit none
    logical, intent(out) :: output
    class(AtmosMesh), target, intent(in) :: model_mesh
    class(MeshTopography), intent(inout) :: topography

    integer :: n
    integer :: ke
    class(LocalMesh3D), pointer :: lcmesh3D
    class(MeshBase3D), pointer :: mesh
    class(LocalMesh2D), pointer :: lcmesh2D
    class(MeshBase2D), pointer :: mesh2D
    !---------------------------------------------------------------------------

    mesh => model_mesh%ptr_mesh

    if ( MKTOPO_TYPE == I_IGNORE ) then
      LOG_NEWLINE
      LOG_PROGRESS(*) 'skip  making topography data'
      output = .false.
    else
      LOG_NEWLINE
      LOG_PROGRESS(*) 'start making topography data'

      ! call PROF_rapstart('_MkTOPO_main',3)   
      
      call mesh%GetMesh2D( mesh2D )

      select case( MKTOPO_TYPE )
      case ( I_FLAT )
        call MKTOPO_flat( mesh2D, topography%topo )
      case ( I_BELLSHAPE )
        call MKTOPO_bellshape( mesh2D, topography%topo )
      case ( I_BELLSHAPE_GLOBAL )
        call MKTOPO_bellshape_global( mesh2D, topography%topo )
      end select

      ! call PROF_rapend  ('_MkTOPO_main',3)
      ! LOG_PROGRESS(*) 'end   making topography data'

      output = .true.
    end if

    return
  end subroutine MKTOPO

  subroutine MKTOPO_write( model_mesh, topography )

    use scale_time, only: &
      NOWDATE   => TIME_NOWDATE,   &
      NOWSUBSEC => TIME_NOWSUBSEC, &
      NOWDAYSEC => TIME_NOWDAYSEC
    use scale_prc, only: PRC_myrank

    use scale_file_base_meshfield, only: FILE_base_meshfield
    use scale_mesh_rectdom2d, only: MeshRectDom2D
    use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
    use scale_mesh_base2d, only: &
      MFTYPE2D_XY => MeshBase2D_DIMTYPEID_XY

    implicit none
    class(AtmosMesh), target, intent(in) :: model_mesh
    class(MeshTopography), intent(inout) :: topography

    type(FILE_base_meshfield) :: file
    class(MeshBase2D), pointer :: mesh

    logical :: file_existed
    integer, parameter :: vid_topo = 1
    !--------------------------------------------------
   
    nullify( mesh )
    call model_mesh%ptr_mesh%GetMesh2D( mesh )
    select type( mesh )
    class is (MeshCubedSphereDom2D)
      call file%Init( 1, meshcubedsphere2D=mesh )
    class is (MeshRectDom2D)      
      call file%Init( 1, mesh2D=mesh )
    end select
       
    call file%Create( OUT_BASENAME, OUT_TITLE, OUT_DTYPE, &
                      file_existed,                       &
                      myrank=PRC_myrank                   )

    if ( .not. file_existed ) then
      call file%Put_GlobalAttribute_time( NOWDATE, NOWSUBSEC )
    end if

    call file%Def_Var( topography%topo, "TOPOGRAPHY", vid_topo, MFTYPE2D_XY, 'XY' )
    call file%End_def()
    call file%Write_var2D( vid_topo, topography%topo, NOWDAYSEC, NOWDAYSEC )
    call file%Close()
    call file%Final()

    return
  end subroutine MKTOPO_write

  !-- private---------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> Make flat mountain  
  subroutine MKTOPO_flat( mesh, topo )
    implicit none

    class(MeshBase2D), intent(in), target :: mesh
    class(MeshField2D), intent(inout) :: topo

    ! flat mountain parameter
    real(RP) :: FLAT_HEIGHT   =  100.0_RP ! height of mountain [m]

    namelist / PARAM_MKTOPO_FLAT / &
       FLAT_HEIGHT

    integer :: ierr    
    integer :: n
    integer :: ke2d
    type(LocalMesh2D), pointer :: lmesh2D
    !--------------------------------

    LOG_INFO("MKTOPO_flat",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_FLAT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_flat",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_flat",*) 'Not appropriate names in namelist PARAM_MKTOPO_FLAT. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_FLAT)

    do n=1, mesh%LOCAL_MESH_NUM
      lmesh2D => mesh%lcmesh_list(n)
      !$omp parallel do
      do ke2D=lmesh2D%NeS, lmesh2D%NeE
        topo%local(n)%val(:,ke2d) = FLAT_HEIGHT
      end do
    end do
    
    return
  end subroutine MKTOPO_flat


  !-----------------------------------------------------------------------------
  !> Make mountain with bell shape (global)
  subroutine MKTOPO_bellshape( mesh, topo )
    implicit none

    class(MeshBase2D), intent(in), target :: mesh
    class(MeshField2D), intent(inout) :: topo

    ! bell-shaped mountain parameter
    real(RP) :: BELL_cx     =   0.0E0_RP ! center location [m]: x-coordinate
    real(RP) :: BELL_cy     =   0.0E3_RP ! center location [m]: y-coordinate
    real(RP) :: BELL_R      =   2.E3_RP  ! half radius        [m]
    real(RP) :: BELL_HEIGHT =  100.0_RP  ! height of mountain [m]
    logical  :: BELL_QUASI_2D = .false. 
    namelist / PARAM_MKTOPO_BELLSHAPE / &
      BELL_cx,     &
      BELL_cy,     &
      BELL_R,      &    
      BELL_HEIGHT, &
      BELL_QUASI_2D

    integer :: ierr    
    integer :: n
    integer :: ke2d
    type(LocalMesh2D), pointer :: lmesh2D
    class(ElementBase2D), pointer :: elem

    real(RP), allocatable :: r2(:)
    real(RP) :: fac_y_quasi2D
    !--------------------------------

    LOG_INFO("MKTOPO_bellshape",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_BELLSHAPE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_bellshape",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_bellshape",*) 'Not appropriate names in namelist PARAM_MKTOPO_BELLSHAPE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_BELLSHAPE)

    if ( BELL_QUASI_2D ) then
      fac_y_quasi2D = 0.0_RP
    else
      fac_y_quasi2D = 1.0_RP
    end if

    do n=1, mesh%LOCAL_MESH_NUM
      lmesh2D => mesh%lcmesh_list(n)
      elem => lmesh2D%refElem2D

      allocate( r2(elem%Np) )
      !$omp parallel do private(r2)
      do ke2D=lmesh2D%NeS, lmesh2D%NeE
        r2(:) = (                   ( lmesh2D%pos_en(:,ke2d,1) - BELL_cx )**2   & 
                  + fac_y_quasi2D * ( lmesh2D%pos_en(:,ke2d,2) - BELL_cy )**2 ) &
                / BELL_R**2
        topo%local(n)%val(:,ke2d) = BELL_HEIGHT / ( 1.0_RP + r2(:) )
      end do
      deallocate( r2 )
    end do
    
    return
  end subroutine MKTOPO_bellshape

  !-----------------------------------------------------------------------------
  !> Make mountain with bell shape (global)
  subroutine MKTOPO_bellshape_global( mesh, topo )
    use scale_const, only: &
      RPlanet => CONST_RADIUS
    implicit none

    class(MeshBase2D), intent(in), target :: mesh
    class(MeshField2D), intent(inout) :: topo

    ! bell-shaped mountain parameter
    real(RP) :: BELL_Clon   =   0.0E0_RP ! center location [rad]: longitude
    real(RP) :: BELL_Clat   =   0.0E3_RP ! center location [rad]: latitude
    real(RP) :: BELL_R      =   2.E3_RP  ! half radius        [m]
    real(RP) :: BELL_HEIGHT =  100.0_RP  ! height of mountain [m]
    namelist / PARAM_MKTOPO_BELLSHAPE_GLOBAL / &
      BELL_Clon,     &
      BELL_Clat,     &
      BELL_R,        &    
      BELL_HEIGHT

    integer :: ierr    
    integer :: n
    integer :: ke2d
    type(LocalMesh2D), pointer :: lmesh2D
    class(ElementBase2D), pointer :: elem

    real(RP), allocatable :: r(:)
    !--------------------------------

    LOG_INFO("MKTOPO_bellshape",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_BELLSHAPE_GLOBAL,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_bellshape",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_bellshape",*) 'Not appropriate names in namelist PARAM_MKTOPO_BELLSHAPE_GLOBAL. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_BELLSHAPE_GLOBAL)

    do n=1, mesh%LOCAL_MESH_NUM
      lmesh2D => mesh%lcmesh_list(n)
      elem => lmesh2D%refElem2D

      allocate( r(elem%Np) )
      !$omp parallel do private(r)
      do ke2D=lmesh2D%NeS, lmesh2D%NeE
        r(:) = RPlanet / BELL_R &
             * acos( sin(BELL_Clat) * sin(lmesh2D%lat(:,ke2d))                                        &
                   + cos(BELL_Clat) * cos(lmesh2D%lat(:,ke2d)) * cos(lmesh2D%lon(:,ke2d) - BELL_Clon) )
        topo%local(n)%val(:,ke2d) = BELL_HEIGHT / ( 1.0_RP + r(:)**2 )
      end do
      deallocate( r )
    end do
    
    return
  end subroutine MKTOPO_bellshape_global

end module mod_mktopo