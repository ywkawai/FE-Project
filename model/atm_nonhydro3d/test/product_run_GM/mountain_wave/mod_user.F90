!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of mountain wave in global model
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
  real(RP) :: Cs     =  0.0_RP   !< Equatorial surface wind shear (for sheared flow) [m-1]
!  real(RP) :: SPONGE_LAYER_tau = 1800.0_RP

  logical :: sponge_layer_flag = .false.
  real(RP), private :: zTop
  real(RP), private :: SPONGE_HEIGHT
  real(RP), private :: SPONGE_EFOLD_SEC         = 600.0_RP
  logical :: lateral_sponge_layer_flag = .false.
  real(RP), private :: LATERAL_SPONGE_EFOLD_SEC = 600.0_RP

  integer :: IniIntrpPolyOrder_h = 8
  integer :: IniIntrpPolyOrder_v = 8

  type(MeshField3D), private :: PT_diff

  type(MeshField3D) :: U_bg, V_bg, T_bg

  logical :: ini_bg_force_flag = .false.
  real(RP) :: ini_bg_force_tend = - 60.0_RP
  real(RP) :: ini_bg_force_tscale = 10.0_RP

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
    namelist / PARAM_USER / &
       USER_do,             &
       ini_bg_force_flag,   &
       ini_bg_force_tend,   &
       ini_bg_force_tscale, &
       sponge_layer_flag,   &
       zTop,                &
       SPONGE_HEIGHT,       &
       SPONGE_EFOLD_SEC,    &
       lateral_sponge_layer_flag, &
       LATERAL_SPONGE_EFOLD_SEC


    integer :: ierr    
    !------------------------------------------

    zTop                 = 20.E3_RP
    SPONGE_HEIGHT        = 15.E3_RP
    SPONGE_EFOLD_SEC     = 900.0_RP

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

    call U_bg%Init( "U_bg", "s-1", atm%mesh%ptr_mesh )
    call V_bg%Init( "V_bg", "s-1", atm%mesh%ptr_mesh )
    call T_bg%Init( "T_bg", "K", atm%mesh%ptr_mesh )

    return
  end subroutine USER_setup

!OCL SERIAL  
  subroutine USER_calc_tendency( this, atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in

    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use scale_cubedsphere_coord_cnv, only: CubedSphereCoordCnv_LonLat2CSVec
  
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRESHYD_VID => AUXVAR_PRESHYDRO_ID,         &
      PRESHYD_REF_VID => AUXVAR_PRESHYDRO_REF_ID, &
      MOMX_p  => PHYTEND_MOMX_ID, &
      MOMY_p  => PHYTEND_MOMY_ID, &
      MOMZ_p  => PHYTEND_MOMZ_ID, &
      RHOH_p  => PHYTEND_RHOH_ID

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars

    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(MeshField3D), pointer :: PRES_hyd
    class(MeshField3D), pointer :: PRES_hyd_ref

    real(RP), allocatable :: sin_lat(:), cos_lat(:)
    real(RP), allocatable :: T(:)
    real(RP), allocatable :: Umet(:,:), Vmet(:,:)

    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem
    integer :: domid
    integer :: ke, ke2D
    !------------------------------------------

    if ( this%USER_do ) then
      call atm%vars%Calc_diagVar( 'PT_diff', PT_diff )
      call FILE_HISTORY_meshfield_in( PT_diff, "perturbation of potential temperature" )
    end if

    ! Set reference hydrostatic pressure
    if ( .not. is_PREShyd_ref_set ) then
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_VID, PRES_hyd )
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_REF_VID, PRES_hyd_ref )

      do domid=1, PRES_hyd_ref%mesh%LOCAL_MESH_NUM
        lcmesh3D => PRES_hyd_ref%mesh%lcmesh_list(domid)
        elem => lcmesh3D%refElem3D

        allocate( Umet(elem%Np,lcmesh3D%Ne), Vmet(elem%Np,lcmesh3D%Ne) )
        allocate( T(elem%Np) )
        allocate( sin_lat(elem%Np), cos_lat(elem%Np) )

        select case( trim(DCMIP_case) )
        case ('2-1', '2-2')
          !$omp parallel do private(ke2D,sin_lat,cos_lat,T)
          do ke=lcmesh3D%NeS, lcmesh3D%NeE
            ke2d = lcmesh3D%EMap3Dto2D(ke)
            sin_lat(:) = sin(lcmesh3D%lat2D(elem%IndexH2Dto3D(:),ke2D))
            cos_lat(:) = cos(lcmesh3D%lat2D(elem%IndexH2Dto3D(:),ke2D))

            T(:) = Teq * ( 1.0_RP - Cs * Ueq**2 * sin_lat(:)**2 / Grav )
            Umet(:,ke) = Ueq * cos_lat(:) * sqrt( 2.0_RP * Teq / T(:) * Cs * lcmesh3D%zlev(:,ke) + T(:) / Teq )
            Vmet(:,ke) = 0.0_RP
            
            PRES_hyd_ref%local(domid)%val(:,ke) = PRES00 * exp(- Grav * lcmesh3D%zlev(:,ke) / ( Rdry * Teq ) )
            T_bg%local(domid)%val(:,ke) = T(:)
          end do

          call CubedSphereCoordCnv_LonLat2CSVec( &
            lcmesh3D%panelID, lcmesh3D%pos_en(:,:,1), lcmesh3D%pos_en(:,:,2),  & ! (in)
            lcmesh3D%gam(:,lcmesh3D%NeS:lcmesh3D%NeE), elem%Np * lcmesh3D%Ne,  & ! (in)
            Umet(:,:), Vmet(:,:),                                              & ! (in)
            U_bg%local(domid)%val(:,lcmesh3D%NeS:lcmesh3D%NeE),                & ! (out)
            V_bg%local(domid)%val(:,lcmesh3D%NeS:lcmesh3D%NeE)                 ) ! (out)
        end select

        deallocate( Umet, Vmet )
        deallocate( sin_lat, cos_lat )
      end do
      is_PREShyd_ref_set = .true.
    end if

    return
  end subroutine USER_calc_tendency

!OCL SERIAL
  subroutine USER_update( this, atm )
    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      MOMX_p  => PHYTEND_MOMX_ID, &
      MOMY_p  => PHYTEND_MOMY_ID, &
      MOMZ_p  => PHYTEND_MOMZ_ID, &
      RHOH_p  => PHYTEND_RHOH_ID
    use scale_time_manager, only:  TIME_NOWSTEP
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars

    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    type(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke2D

    real(RP), allocatable :: DENS(:), T(:)
    real(RP), allocatable :: rtau(:), sfac(:), rsfac(:)
    real(RP), allocatable :: lon(:), lat(:)

    real(RP) :: rtau_ini_bg
    real(RP) :: rtau_sponge
    real(RP) :: rtau_lateral_sponge
 
    real(DP) :: dt
    real(RP) :: tsec
    real(RP) :: Gamm
    !----------------------------------------------------------

    dt = atm%time_manager%dtsec
    tsec = dt * real( TIME_NOWSTEP, kind=RP )
    gamm = CpDry / CvDry 

    rtau_ini_bg = 0.0_RP
    if ( ini_bg_force_flag .and. ( tsec < ini_bg_force_tend ) ) then
      rtau_ini_bg = 1.0_RP / ini_bg_force_tscale
    end if

    rtau_sponge = 0.0_RP
    if ( sponge_layer_flag ) then
      rtau_sponge = 1.0_RP / SPONGE_EFOLD_SEC
    end if

    rtau_lateral_sponge = 0.0_RP
    if ( lateral_sponge_layer_flag ) then
      rtau_lateral_sponge = 1.0_RP / LATERAL_SPONGE_EFOLD_SEC
    end if

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      

      call AtmosVars_GetLocalMeshPhyAuxVars( n, atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                         )
      
      elem3D => lcmesh%refElem3D

      allocate( DENS(elem3D%Np), T(elem3D%Np) )
      allocate( rtau(elem3D%Np), sfac(elem3D%Np), rsfac(elem3D%Np) )
      allocate( lon(elem3D%Np), lat(elem3D%Np) )

      !$omp parallel do private(DENS, T, lon, lat, rtau, sfac, rsfac, ke2D)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)

        lon(:) = lcmesh%lon2D(elem3D%IndexH2Dto3D,ke2D)
        lat(:) = lcmesh%lat2D(elem3D%IndexH2Dto3D,ke2D)

        rtau(:) = rtau_ini_bg
        where ( lcmesh%zlev(:,ke) > SPONGE_HEIGHT )
          rtau(:) = rtau(:) + rtau_sponge * 0.5_RP * ( 1.0_RP - cos( PI * ( lcmesh%zlev(:,ke) - SPONGE_HEIGHT ) / ( zTop - SPONGE_HEIGHT ) ) )
        end where
        where ( lon(:) < PI * 0.5_RP )
          rtau(:) = rtau(:) + rtau_lateral_sponge * 0.5_RP * ( 1.0_RP - cos( PI * ( lon(:) - PI * 0.5_RP ) / ( PI * 0.5_RP ) ) )
        end where
        where ( lon(:) > 1.5_RP * PI )
          rtau(:) = rtau(:) + rtau_lateral_sponge * 0.5_RP * ( 1.0_RP - cos( PI * ( lon(:) - PI * 1.5_RP ) / ( PI * 0.5_RP ) ) )
        end where

        sfac(:) = dt * rtau(:)
        rsfac(:) = 1.0_RP / ( 1.0_RP + sfac(:) )

        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)
        T(:) =  PRES%val(:,ke) / ( Rdry * DENS(:) )

        MOMX%val(:,ke) = ( MOMX%val(:,ke) + sfac(:) * DENS(:) * U_bg%local(n)%val(:,ke) ) * rsfac(:)
        MOMY%val(:,ke) = ( MOMY%val(:,ke) + sfac(:) * DENS(:) * V_bg%local(n)%val(:,ke) ) * rsfac(:)
        MOMZ%val(:,ke) = MOMZ%val(:,ke) * rsfac(:)

        !- For the case of d DRHOT /dt = dens * Cp * ( Teq - T ) / tauT     
        !  <- It is based on the forcing form in Held and Surez in which the temperature evolution equation is assumed to be dT/dt = R/Cp * T/p * dp/dt + (Teq - T) / tauT
        DRHOT%val(:,ke) = DRHOT%val(:,ke) &
                        - dt * rtau(:) * ( 1.0_RP - T_bg%local(n)%val(:,ke) / T(:) ) * DENS(:) * PT%val(:,ke)     &
                        / ( 1.0_RP + dt * rtau(:) * ( 1.0_RP + (gamm - 1.0_RP) * T_bg%local(n)%val(:,ke) / T(:) ) )  
      end do

      deallocate( DENS, T )
      deallocate( rtau, sfac, rsfac )
      deallocate( lat )
    end do
    
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

    real(RP) :: T(elem%Np)
    real(RP) :: sin_lat(elem%Np), cos_lat(elem%Np)
    real(RP) :: MOMX_met(elem%Np,lcmesh%Ne)
    real(RP) :: MOMY_met          (elem%Np,lcmesh%Ne)

    integer :: ke, ke2d

    type(LocalMesh2D), pointer :: lmesh2D
    real(RP) :: H0_pres
    real(RP) :: Ueq_
    !-----------------------------------------------------------------------------

    call read_exp_params()

    select case( trim(DCMIP_case) )
    case( '2-0-0', '2-0-1' )
      call hydrostatic_calc_basicstate_constTLAPS( DENS_hyd, PRES_hyd, & ! (out)
      TLAPS, TEMP0, PRES00,                                            & ! (in)
      x, y, lcmesh%zlev, lcmesh, elem                                  ) ! (in)
    case ('2-1', '2-2')
      !$omp parallel do private( ke2d, sin_lat, cos_lat, T, Ueq_ )
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2d = lcmesh%EMap3Dto2D(ke)
        sin_lat(:) = sin(lcmesh%lat2D(elem%IndexH2Dto3D(:),ke2D))
        cos_lat(:) = cos(lcmesh%lat2D(elem%IndexH2Dto3D(:),ke2D))

        T(:) = Teq * ( 1.0_RP - Cs * Ueq**2 * sin_lat(:)**2 / Grav )

        PRES_hyd(:,ke) = PRES00 * exp( - 0.5_RP * Ueq**2 / ( Rdry * Teq ) * sin_lat(:)**2 - Grav * lcmesh%zlev(:,ke) / ( Rdry * T(:) ) )
        DENS_hyd(:,ke) = PRES_hyd(:,ke) / ( Rdry * T(:) )
  
        MOMX_met(:,ke) = DENS_hyd(:,ke) * Ueq * cos_lat(:) * sqrt( 2.0_RP * Teq / T(:) * Cs * lcmesh%zlev(:,ke) + T(:) / Teq )
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
      Ueq, Cs,             &
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
