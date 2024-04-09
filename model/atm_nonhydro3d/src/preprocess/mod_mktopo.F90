!-------------------------------------------------------------------------------
!> module INITIAL
!!
!! @par Description
!!          A module for preparing topography data 
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
  
  integer, public :: MKTOPO_TYPE                          = -1
  integer, parameter, public :: I_IGNORE                  = 0
  integer, parameter, public :: I_INPUT_FILE              = 1
  integer, parameter, public :: I_FLAT                    = 2
  integer, parameter, public :: I_BELLSHAPE               = 3
  integer, parameter, public :: I_SCAER                   = 4
  integer, parameter, public :: I_BELLSHAPE_GLOBAL        = 5
  integer, parameter, public :: I_SCAHER_GLOBAL           = 6
  integer, parameter, public :: I_BAROCWAVE_GLOBAL_JW2006 = 7  

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  private :: MKTOPO_flat
  private :: MKTOPO_bellshape
  private :: MKTOPO_Schaer_type_mountain
  private :: MKTOPO_bellshape_global
  private :: MKTOPO_barocwave_global_JW2006

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
!OCL SERIAL
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
    case('INPUT_FILE')
      MKTOPO_TYPE = I_INPUT_FILE
    case('FLAT')
      MKTOPO_TYPE = I_FLAT  
    case('BELLSHAPE')
      MKTOPO_TYPE = I_BELLSHAPE
    case('SCHAER')
      MKTOPO_TYPE = I_SCAER
    case('BELLSHAPE_GLOBAL')
      MKTOPO_TYPE = I_BELLSHAPE_GLOBAL
    case('SCHAER_GLOBAL')
      MKTOPO_TYPE = I_SCAHER_GLOBAL
    case('BAROCWAVE_GLOBAL_JW2006')
      MKTOPO_TYPE = I_BAROCWAVE_GLOBAL_JW2006      
    case default
      LOG_ERROR("MKTOPO_setup",*) 'Not appropriate toponame. Check!', toponame
      call PRC_abort      
    end select
    
    return
  end subroutine MKTOPO_setup

  !-----------------------------------------------------------------------------
  !> Driver
!OCL SERIAL
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

      call PROF_rapstart('_MkTOPO_main', 3)   
      
      call mesh%GetMesh2D( mesh2D )

      select case( MKTOPO_TYPE )
      case ( I_INPUT_FILE )
        call MKTOPO_input_file( mesh2D, topography%topo )
      case ( I_FLAT )
        call MKTOPO_flat( mesh2D, topography%topo )
      case ( I_BELLSHAPE )
        call MKTOPO_bellshape( mesh2D, topography%topo )
      case ( I_SCAER )
        call MKTOPO_Schaer_type_mountain( mesh2D, topography%topo )
      case ( I_BELLSHAPE_GLOBAL )
        call MKTOPO_bellshape_global( mesh2D, topography%topo )
      case ( I_SCAHER_GLOBAL )
        call MKTOPO_Schaer_type_global( mesh2D, topography%topo )        
      case ( I_BAROCWAVE_GLOBAL_JW2006 )
        call MKTOPO_barocwave_global_JW2006( mesh2D, topography%topo )
      end select

      call PROF_rapend  ('_MkTOPO_main', 3)
      LOG_PROGRESS(*) 'end making topography data'

      output = .true.
    end if

    return
  end subroutine MKTOPO

  !-----------------------------------------------------------------------------
  !> Output topography data
!OCL SERIAL
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
  !> Make topography from input file
  !!
!OCL SERIAL  
  subroutine MKTOPO_input_file( mesh, topo )
    use scale_file_base_meshfield, only: FILE_base_meshfield
    use scale_mesh_rectdom2d, only: MeshRectDom2D
    use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
    use scale_mesh_base2d, only: &
      MFTYPE2D_XY => MeshBase2D_DIMTYPEID_XY

    implicit none

    class(MeshBase2D), intent(in), target :: mesh
    class(MeshField2D), intent(inout) :: topo

    character(len=H_LONG) :: IN_BASENAME
    character(len=H_MID) :: VARNAME

    namelist / PARAM_MKTOPO_INPUT_FILE / &
      IN_BASENAME, &
      VARNAME

    integer :: ierr    
    integer :: n
    integer :: ke2d
    type(LocalMesh2D), pointer :: lmesh2D

    type(FILE_base_meshfield) :: file

    logical :: file_existed
    !--------------------------------

    LOG_INFO("MKTOPO_input_file",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_INPUT_FILE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_input_file",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_input_file",*) 'Not appropriate names in namelist PARAM_MKTOPO_INPUT_FILE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_INPUT_FILE)

    !-- Read topography data

    select type( mesh )
    class is (MeshCubedSphereDom2D)
      call file%Init( 1, meshcubedsphere2D=mesh )
    class is (MeshRectDom2D)      
      call file%Init( 1, mesh2D=mesh )
    end select

    call file%Open( IN_BASENAME, myrank=PRC_myrank)
    call file%Read_Var( MFTYPE2D_XY, VARNAME, topo )

    call file%Close()
    call file%Final()

    return
  end subroutine MKTOPO_input_file


  !-----------------------------------------------------------------------------
  !> Make flat mountain  
  !!
!OCL SERIAL  
  subroutine MKTOPO_flat( mesh, topo )
    implicit none

    class(MeshBase2D), intent(in), target :: mesh
    class(MeshField2D), intent(inout) :: topo

    ! flat mountain parameter
    real(RP) :: FLAT_HEIGHT   =  100.0_RP !< height of mountain [m]

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
  !> Make mountain with bell shape
  !!
!OCL SERIAL
  subroutine MKTOPO_bellshape( mesh, topo )
    implicit none

    class(MeshBase2D), intent(in), target :: mesh
    class(MeshField2D), intent(inout) :: topo

    ! bell-shaped mountain parameter
    real(RP) :: BELL_cx     =   0.0E0_RP !< center location [m]: x-coordinate
    real(RP) :: BELL_cy     =   0.0E3_RP !< center location [m]: y-coordinate
    real(RP) :: BELL_R      =   2.E3_RP  !< half radius        [m]
    real(RP) :: BELL_HEIGHT =  100.0_RP  !< height of mountain [m]
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
  !> Make Schaer-type mountain
  !!
  !! @par Reference
  !! - Schaer et al., 2002, MWR, Vol.130, 2459-2480
  !! - Klemp et al., 2003,  MWR, Vol.131, 1229-1239
!OCL SERIAL
  subroutine MKTOPO_Schaer_type_mountain( mesh, topo )
    implicit none

    class(MeshBase2D), intent(in), target :: mesh
    class(MeshField2D), intent(inout) :: topo

    ! Schaer-type mountain parameter
    real(RP) :: SCHAER_CX       =  25.E3_RP !< center location [m]: x
    real(RP) :: SCHAER_RX       =   5.E3_RP !< bubble radius   [m]: x
    real(RP) :: SCHAER_LAMBDA   =   4.E3_RP !< wavelength of wavelike perturbation [m]: x
    real(RP) :: SCHAER_HEIGHT   =  250.0_RP !< height of mountain [m]
    namelist / PARAM_MKTOPO_SCHAER / &
       SCHAER_CX,     &
       SCHAER_RX,     &
       SCHAER_LAMBDA, &
       SCHAER_HEIGHT

    integer :: ierr    
    integer :: n
    integer :: ke2d
    type(LocalMesh2D), pointer :: lmesh2D
    class(ElementBase2D), pointer :: elem

    real(RP), allocatable :: dist(:)
    real(RP) :: fac_y_quasi2D
    !--------------------------------

    LOG_INFO("MKTOPO_Schaer",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_SCHAER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_Schaer",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_Schaer",*) 'Not appropriate names in namelist PARAM_MKTOPO_SCAHER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_SCHAER)

    do n=1, mesh%LOCAL_MESH_NUM
      lmesh2D => mesh%lcmesh_list(n)
      elem => lmesh2D%refElem2D

      allocate( dist(elem%Np) )
      !$omp parallel do private(dist)
      do ke2D=lmesh2D%NeS, lmesh2D%NeE
        dist(:) = exp( - ( lmesh2D%pos_en(:,ke2d,1) - SCHAER_CX )**2 / SCHAER_RX**2 )
        topo%local(n)%val(:,ke2d) = SCHAER_HEIGHT * dist(:) * ( cos( PI * ( lmesh2D%pos_en(:,ke2d,1) - SCHAER_CX ) / SCHAER_LAMBDA ) )**2
      end do
      deallocate( dist )
    end do
    
    return
  end subroutine MKTOPO_Schaer_type_mountain

  !-----------------------------------------------------------------------------
  !> Make mountain with bell shape (global)
  !!
  subroutine MKTOPO_bellshape_global( mesh, topo )
    implicit none

    class(MeshBase2D), intent(in), target :: mesh
    class(MeshField2D), intent(inout) :: topo

    ! bell-shaped mountain parameter
    real(RP) :: BELL_Clon   =   0.0E0_RP !< center location [rad]: longitude
    real(RP) :: BELL_Clat   =   0.0E0_RP !< center location [rad]: latitude
    real(RP) :: BELL_R      =   2.E3_RP  !< half radius        [m]
    real(RP) :: BELL_HEIGHT =  100.0_RP  !< height of mountain [m]
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

    LOG_INFO("MKTOPO_bellshape_global",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_BELLSHAPE_GLOBAL,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_bellshape_global",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_bellshape_global",*) 'Not appropriate names in namelist PARAM_MKTOPO_BELLSHAPE_GLOBAL. Check!'
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

  !-----------------------------------------------------------------------------
  !> Make Schaer-type mountain (global)
  !!
  !! @par Reference
  !! - Klemp, J. B., Skamarock, W. C., & Park, S. H., 2015: Idealized global nonhydrostatic atmospheric test cases on a reduced‐radius sphere. Journal of Advances in Modeling Earth Systems, 7(3), 1155-1177.
  !!
  subroutine MKTOPO_Schaer_type_global( mesh, topo )
    implicit none

    class(MeshBase2D), intent(in), target :: mesh
    class(MeshField2D), intent(inout) :: topo

    ! bell-shaped mountain parameter
    real(RP) :: SCHAER_Clon     !< center location [rad]: longitude
    real(RP) :: SCHAER_Clat     !< center location [rad]: latitude
    real(RP) :: SCHAER_R        !< half-width of mountain   [m]
    real(RP) :: SCHAER_LAMBDA   !< wavelength of mountain oscillation [m]
    integer  :: SCHAER_SHAPE_ID !< 1: gaussian, 2: cosine bell
    logical  :: SCHAER_MERI_TAPER_FLAG
    real(RP)  :: SCHAER_MERI_TAPER_TANH_Clat
    real(RP)  :: SCHAER_MERI_TAPER_TANH_LatWidth
    real(RP) :: SCHAER_HEIGHT   !< height of mountain [m]
    logical :: quasi_2D_flag

    namelist / PARAM_MKTOPO_SCHAER_GLOBAL / &
      SCHAER_Clon,     &
      SCHAER_Clat,     &
      SCHAER_R,        &
      SCHAER_LAMBDA,   &
      SCHAER_SHAPE_ID, &
      SCHAER_HEIGHT,   &
      SCHAER_MERI_TAPER_FLAG,          &
      SCHAER_MERI_TAPER_TANH_Clat,     &
      SCHAER_MERI_TAPER_TANH_LatWidth, &
      quasi_2D_flag

    integer :: ierr    
    integer :: n
    integer :: ke2d
    type(LocalMesh2D), pointer :: lmesh2D
    class(ElementBase2D), pointer :: elem

    real(RP), allocatable :: r(:)
    real(RP), allocatable :: mwt_shape(:)

    !--------------------------------

    LOG_INFO("MKTOPO_Schaer_global",*) 'Setup'

    SCHAER_Clon     = 1.0_RP / 4.0_RP * PI
    SCHAER_Clat     = 0.0_RP
    SCHAER_R        = RPlanet * 3.0_RP / 4.0_RP * PI
    SCHAER_LAMBDA   = RPlanet * PI / 16.0_RP
    SCHAER_SHAPE_ID = 1
    SCHAER_HEIGHT   = 2000.0_RP
    quasi_2D_flag   = .false.

    SCHAER_MERI_TAPER_FLAG = .false.

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_SCHAER_GLOBAL,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_Schaer_global",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_Schaer_global",*) 'Not appropriate names in namelist PARAM_MKTOPO_SCHAER_GLOBAL. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_SCHAER_GLOBAL)

    do n=1, mesh%LOCAL_MESH_NUM
      lmesh2D => mesh%lcmesh_list(n)
      elem => lmesh2D%refElem2D

      allocate( r(elem%Np), mwt_shape(elem%Np) )
      !$omp parallel do private(r, mwt_shape)
      do ke2D=lmesh2D%NeS, lmesh2D%NeE
        if ( quasi_2D_flag ) then
!          r(:) = RPlanet * ( min( abs(lmesh2D%lon(:,ke2d) - SCHAER_Clon), abs(lmesh2D%lon(:,ke2d) - (SCHAER_Clon + 2.0_RP * PI ) ) ) )
          r(:) = RPlanet * abs(lmesh2D%lon(:,ke2d) - SCHAER_Clon)
        else
          r(:) = RPlanet &
              * acos( sin(SCHAER_Clat) * sin(lmesh2D%lat(:,ke2d))                                          &
                    + cos(SCHAER_Clat) * cos(lmesh2D%lat(:,ke2d)) * cos(lmesh2D%lon(:,ke2d) - SCHAER_Clon) )
        end if

        if ( SCHAER_SHAPE_ID == 1 ) then
          mwt_shape(:) = exp( - ( r(:) / SCHAER_R )**2 ) 
        else if ( SCHAER_SHAPE_ID == 2 ) then
          where ( r(:) < SCHAER_R )
            mwt_shape(:) = 0.5_RP * ( 1.0_RP + cos( PI * r(:) / SCHAER_R ) )
          elsewhere
            mwt_shape(:) = 0.0_RP
          end where          
        end if

        if ( quasi_2D_flag ) then
          mwt_shape(:) = mwt_shape(:) * cos(lmesh2D%lat(:,ke2d))
          if ( SCHAER_MERI_TAPER_FLAG ) then
            mwt_shape(:) = mwt_shape(:) * & 
              0.5_RP * ( 1.0_RP - tanh( ( abs(lmesh2D%lat(:,ke2d)) - SCHAER_MERI_TAPER_TANH_Clat ) / SCHAER_MERI_TAPER_TANH_LatWidth ) )
          end if
        end if

        topo%local(n)%val(:,ke2d) = SCHAER_HEIGHT *  mwt_shape(:) * cos( PI * r(:) / SCHAER_LAMBDA )**2
      end do
      deallocate( r, mwt_shape )
    end do
    
    return
  end subroutine MKTOPO_Schaer_type_global

  !-----------------------------------------------------------------------------
  !> Make a topography for a global test case of idealized baroclinic wave in Jablonowski and Williamson (2006)
  !!
  !! @par Reference
  !! - Jablonowski, C., & Williamson, D. L., 2006: A baroclinic instability test case for atmospheric model dynamical cores. Quarterly Journal of the Royal Meteorological Society: A journal of the atmospheric sciences, applied meteorology and physical oceanography, 132(621C), 2943-2975.
  !!
!OCL SERIAL
  subroutine MKTOPO_barocwave_global_JW2006( mesh, topo )
    use mod_mktopo_util, only: &
      calc_topo => mktopoutil_barocwave_global_JW2006_calc_topo
    implicit none

    class(MeshBase2D), intent(in), target :: mesh
    class(MeshField2D), intent(inout) :: topo

    ! parameters of baroclinic wave test case in JW2006
    real(RP) :: ETA0 = 0.252_RP  !< The value of η at a reference level (position of the jet)
    real(RP) :: U0   = 35.E0_RP  !< The parameter associated with zonal jet maximum amplitude  [m/s]

    namelist / PARAM_MKTOPO_BAROCWAVE_GLOBAL_JW2006 / &
      ETA0, U0
    
    integer :: ierr    
    integer :: n
    integer :: ke2d
    type(LocalMesh2D), pointer :: lmesh2D
    class(ElementBase2D), pointer :: elem
    !--------------------------------

    LOG_INFO("MKTOPO_barocwave_JW2006",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF, nml=PARAM_MKTOPO_BAROCWAVE_GLOBAL_JW2006, iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_barocwave_JW2006",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_barocwave_JW2006",*) 'Not appropriate names in namelist PARAM_MKTOPO_BAROCWAVE_GLOBAL_JW2006. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_BAROCWAVE_GLOBAL_JW2006)

    do n=1, mesh%LOCAL_MESH_NUM
      lmesh2D => mesh%lcmesh_list(n)
      elem => lmesh2D%refElem2D
      
      !$omp parallel do
      do ke2D=lmesh2D%NeS, lmesh2D%NeE
        call calc_topo( topo%local(n)%val(:,ke2d), & ! (out)
          U0, ETA0, lmesh2D%lat(:,ke2d), elem%Np )   ! (in)
      end do
    end do
    
    return
  end subroutine MKTOPO_barocwave_global_JW2006

end module mod_mktopo