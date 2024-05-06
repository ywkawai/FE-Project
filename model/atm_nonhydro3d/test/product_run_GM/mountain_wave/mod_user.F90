!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of mountain wave in global model
!!          Note that we evaluate the forcing terms with inital spinup and sponge layer in RK loop of dynamical process. 
!!          The corresponding program is implemented in FElib/fluid_dyn_solver/scale_atm_dyn_dgm_globanonhydro3d_rhot_heve.F90. 
!!          When compiling it, we need to specify -DSCALE_PRODUCT_RUN_GM_MOUNTAIN_WAVE -DSCALE_PRODUCT_RUN_GM_MOUNTAIN_WAVE_BC for FFLAG_TESTCASE (see Makefile).
!!
!!
!! @p
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort  

  use scale_const, only: &
    PI => CONST_PI,        &
    GRAV => CONST_GRAV,    &
    OHM => CONST_OHM,      &
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  &
    CVdry => CONST_CVdry,  &
    PRES00 => CONST_PRE00, &
    RPlanet => CONST_RADIUS

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: MeshField3D

  use mod_user_base, only: UserBase
  use mod_experiment, only: Experiment

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  type, public, extends(UserBase) :: User
  contains
    procedure :: mkinit_ => USER_mkinit
    generic :: mkinit => mkinit_
    procedure :: setup_ => USER_setup
    generic :: setup => setup_
    procedure :: calc_tendency => USER_calc_tendency
    procedure :: update => USER_update
  end type User


  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT) :: DCMIP_case = '2-0-0' !< 2-0-0, 2-0-1, 2-1, 2-2
  real(RP) :: TEMP0  = 300.0_RP
  real(RP) :: TLAPS  = 6.5E-3_RP  
  real(RP) :: Teq    = 300_RP    !< Reference surface temperature at the equator [K]
  real(RP) :: Ueq    =  0.0_RP   !< Reference zonal wind velocity [m/s]
  real(RP) :: Ueq0   =  0.0_RP   !< Reference zonal wind velocity [m/s]
  real(RP) :: Cs     =  0.0_RP   !< Equatorial surface wind shear (for sheared flow) [m-1]

  integer :: IniIntrpPolyOrder_h = 8
  integer :: IniIntrpPolyOrder_v = 8

  type(MeshField3D), private :: PT_diff
  logical :: is_PREShyd_ref_set

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit ( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('mountain_wave_global')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_mountain_wave )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL
  subroutine USER_setup( this, atm )
    implicit none
    class(User), intent(inout) :: this    
    class(AtmosComponent), intent(inout) :: atm

    logical :: USER_do                   = .false. !< do user step?
    character(len=H_SHORT) :: SPONGE_LAYER_FUNC_NAME

    namelist / PARAM_USER / &
       USER_do

    integer :: ierr    
    !------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    call this%UserBase%Setup( atm, USER_do )

    !-
    call read_exp_params()

    is_PREShyd_ref_set = .false. 

    if ( this%USER_do ) call PT_diff%Init( 'PT_diff', 'K', atm%mesh%ptr_mesh )

    return
  end subroutine USER_setup

!OCL SERIAL  
  subroutine USER_calc_tendency( this, atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in

    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRESHYD_VID => AUXVAR_PRESHYDRO_ID,         &
      PRESHYD_REF_VID => AUXVAR_PRESHYDRO_REF_ID

    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(MeshField3D), pointer :: PRES_hyd
    class(MeshField3D), pointer :: PRES_hyd_ref

    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem
    integer :: domid
    integer :: ke
    !------------------------------------------

    ! Set reference hydrostatic pressure
    if ( .not. is_PREShyd_ref_set ) then
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_VID, PRES_hyd )
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_REF_VID, PRES_hyd_ref )

      do domid=1, PRES_hyd_ref%mesh%LOCAL_MESH_NUM
        lcmesh3D => PRES_hyd_ref%mesh%lcmesh_list(domid)
        elem => lcmesh3D%refElem3D

        select case( trim(DCMIP_case) )
        case ('2-1', '2-2')
          !$omp parallel do
          do ke=lcmesh3D%NeS, lcmesh3D%NeE
            PRES_hyd_ref%local(domid)%val(:,ke) = PRES00 * exp(- Grav * lcmesh3D%zlev(:,ke) / ( Rdry * Teq ) )
          end do
        end select
      end do
      is_PREShyd_ref_set = .true.
    end if

    return
  end subroutine USER_calc_tendency

!OCL SERIAL
  subroutine USER_update( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm
    !----------------------------------------------------------
    return
  end subroutine USER_update

  !------

!OCL SERIAL  
  subroutine exp_SetInitCond_mountain_wave( this,                 &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
     
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constTLAPS
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSVec
  
    use mod_experiment, only: &
      TracerLocalMeshField_ptr   
       
    implicit none

    class(Experiment), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
    type(TracerLocalMeshField_ptr), intent(inout) :: tracer_field_list(:)    
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_ymin, dom_ymax
    real(RP), intent(in) :: dom_zmin, dom_zmax    

    real(RP) :: T(elem%Np), T0(elem%Np), DENS0(elem%Np)
    real(RP) :: PRES_hyd0(elem%Np)
    real(RP) :: sin_lat(elem%Np), cos_lat(elem%Np)
    real(RP) :: MOMX_met(elem%Np,lcmesh%Ne)
    real(RP) :: MOMY_met(elem%Np,lcmesh%Ne)

    integer :: ke, ke2d

    type(LocalMesh2D), pointer :: lmesh2D
    real(RP) :: H0_pres

    real(RP) :: rgamm
    !-----------------------------------------------------------------------------

    call read_exp_params()

    rgamm = CVDry / CPDry

    select case( trim(DCMIP_case) )
    case( '2-0-0', '2-0-1' )
      call hydrostatic_calc_basicstate_constTLAPS( DENS_hyd, PRES_hyd, & ! (out)
      TLAPS, TEMP0, PRES00,                                            & ! (in)
      x, y, lcmesh%zlev, lcmesh, elem                                  ) ! (in)
    case ('2-1', '2-2')
      !$omp parallel do private( ke2d, sin_lat, cos_lat, T, T0, DENS0, PRES_hyd0 )
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2d = lcmesh%EMap3Dto2D(ke)
        sin_lat(:) = sin(lcmesh%lat2D(elem%IndexH2Dto3D(:),ke2D))
        cos_lat(:) = cos(lcmesh%lat2D(elem%IndexH2Dto3D(:),ke2D))

        T (:) = Teq * ( 1.0_RP - Cs * Ueq **2 * sin_lat(:)**2 / Grav )
        T0(:) = Teq * ( 1.0_RP - Cs * Ueq0**2 * sin_lat(:)**2 / Grav )

        PRES_hyd(:,ke) = PRES00 * exp( - 0.5_RP * Ueq**2 / ( Rdry * Teq ) * sin_lat(:)**2 - Grav * lcmesh%zlev(:,ke) / ( Rdry * T(:) ) )
        DENS_hyd(:,ke) = PRES_hyd(:,ke) / ( Rdry * T(:) )

        PRES_hyd0(:) = PRES00 * exp( - 0.5_RP * Ueq0**2 / ( Rdry * Teq ) * sin_lat(:)**2 - Grav * lcmesh%zlev(:,ke) / ( Rdry * T0(:) ) )
        DRHOT(:,ke) = PRES00 / Rdry * ( &
          ( PRES_hyd0(:) / PRES00 )**rgamm - ( PRES_hyd(:,ke) / PRES00 )**rgamm )
        
        DENS0(:) = PRES_hyd0(:) / ( Rdry * T0(:) ) 
        DDENS(:,ke) = DENS0(:) - PRES_hyd(:,ke) / ( Rdry * T(:) )

        MOMX_met(:,ke) = DENS0(:) * Ueq0 * cos_lat(:) * sqrt( 2.0_RP * Teq / T0(:) * Cs * lcmesh%zlev(:,ke) + T0(:) / Teq )
        MOMY_met(:,ke) = 0.0_RP
      end do

      call CubedSphereCoordCnv_LonLat2CSVec( &
        lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),    & ! (in)
        lcmesh%gam(:,lcmesh%NeS:lcmesh%NeE), elem%Np * lcmesh%Ne,      & ! (in)
        MOMX_met(:,:), MOMY_met(:,:),                                  & ! (in)
        MOMX(:,lcmesh%NeS:lcmesh%NeE), MOMY(:,lcmesh%NeS:lcmesh%NeE)   ) ! (out)

    case default
      LOG_ERROR("MOUNTAIN_WAVE_setup",*) 'Not appropriate DCMIP case. Check!'
      call PRC_abort
    end select
    

    return
  end subroutine exp_SetInitCond_mountain_wave

  subroutine read_exp_params()
    implicit none

    namelist /PARAM_EXP/ &
      DCMIP_case,          &
      TEMP0, TLAPS,        &
      Ueq, Ueq0, Cs,       &
      TLAPS,               &
      IniIntrpPolyOrder_h, &
      IniIntrpPolyOrder_v!, &
!      SPONGE_LAYER_tau

    integer :: ierr
    !-----------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MOUNTAIN_WAVE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MOUNTAIN_WAVE_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    return
  end subroutine read_exp_params

end module mod_user
