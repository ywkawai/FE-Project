!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of gravity wave
!!
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
    PI => CONST_PI
  
  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: MeshField3D

  use mod_user_base, only: UserBase
  use mod_experiment, only: Experiment
  use mod_atmos_component, only: &
    AtmosComponent

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
  type(MeshField3D), private :: PT_diff

  real(RP) :: lonc   = 0.0_RP
  real(RP) :: latc   = 0.0_RP
  integer  :: nv     = 1
  real(RP) :: rh
  real(RP) :: Zt

  real(RP) :: THETA0 = 300.0_RP
  real(RP) :: TEMP0  = 300.0_RP
  real(RP) :: BruntVaisalaFreq = 0.01_RP  
  real(RP) :: DTHETA = 0.01_RP
!  real(RP) :: DTEMP  = 0.01_RP

  integer :: Ini_GP_PolyOrder_h = 8
  integer :: Ini_GP_PolyOrder_v = 8

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

    call exp_manager%Init( 'inertia_gravity_wave_global' )
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_inertia_gravity_wave )
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

    !-
    call this%UserBase%Setup( atm, USER_do )

    if ( this%USER_do ) then
      call PT_diff%Init( 'PT_diff', 'K', atm%mesh%ptr_mesh )
    end if

    is_PREShyd_ref_set = .false. 

    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRESHYD_VID => AUXVAR_PRESHYDRO_ID,        &
      PRESHYD_REF_VID => AUXVAR_PRESHYDRO_REF_ID    
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(MeshField3D), pointer :: PRES_hyd
    class(MeshField3D), pointer :: PRES_hyd_ref
    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: domid, ke
    !------------------------------------------


    ! Set reference hydrostatic pressure
    if ( .not. is_PREShyd_ref_set ) then
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_VID, PRES_hyd )
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_REF_VID, PRES_hyd_ref )

      do domid=1, PRES_hyd_ref%mesh%LOCAL_MESH_NUM
        lcmesh3D => PRES_hyd_ref%mesh%lcmesh_list(domid)
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          PRES_hyd_ref%local(domid)%val(:,ke) = PRES_hyd%local(domid)%val(:,ke)
        end do
      end do
      is_PREShyd_ref_set = .true.
    end if

    if ( this%USER_do ) then
      call atm%vars%Calc_diagVar( 'PT_diff', PT_diff )
      call FILE_HISTORY_meshfield_in( PT_diff, "perturbation of potential temperature" )
    end if

    return
  end subroutine USER_calc_tendency

  !------
!OCL SERIAL  
  subroutine exp_SetInitCond_inertia_gravity_wave( this,                 &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
        
    use scale_const, only: &
      GRAV => CONST_GRAV,    &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      PRES00 => CONST_PRE00, &
      RPlanet => CONST_RADIUS    
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constT,     &
      hydrostatic_calc_basicstate_constBVFreq
    use mod_mkinit_util, only: &
      mkinitutil_calc_cosinebell_global, &
      mkinitutil_GalerkinProjection_global

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
    
    character(H_MID) :: InitAtmType     = 'ConstBVFreq' ! or 'Isothermal'
    character(H_MID) :: InitPerturbType = 'GAUSSIAN'    ! or 'COSBELL'
  
    namelist /PARAM_EXP/ &
      InitAtmType,       &
      InitPerturbType,   &
      THETA0, DTHETA,    &
      BruntVaisalaFreq,  &
      TEMP0,             &
!      DTEMP,            &
      lonc, latc, rh,    &
      nv,                &
      Ini_GP_PolyOrder_h, &
      Ini_GP_PolyOrder_v
    
    real(RP) :: EXNER(elem%Np)
  
    real(RP) :: RovCp

    real(RP), allocatable :: PT_purtub(:,:)
    real(RP), allocatable :: TEMP_purtub(:,:) 

    integer :: ke
    integer :: ierr
    !-----------------------------------------------------------------------------

    rh   = RPlanet / 3.0_RP
    lonc = PI

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("INERTIA_GRAVITY_WAVE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("INERTIA_GRAVITY_WAVE_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    
    Zt = dom_zmax - dom_zmin

    !---
    select case(InitAtmType)
    case( 'ISOTHERMAL' )
      BruntVaisalaFreq = Grav / sqrt( CpDry * TEMP0 )
    end select

    call hydrostatic_calc_basicstate_constBVFreq( DENS_hyd, PRES_hyd, & ! (out)
      BruntVaisalaFreq, THETA0, PRES00,                               & ! (in)
      x, y, z, lcmesh, elem                                           ) ! (in)

    allocate( PT_purtub(elem%Np,lcmesh%NeA) )
    !$omp parallel do
    do ke=lcmesh%NeS, lcmesh%NeE
      PT_purtub(:,ke) = 0.0_RP
    end do
    
    select case(InitPerturbType)
    case( 'GAUSSIAN' )
    
      call mkinitutil_GalerkinProjection_global( &
        DDENS,                                                           &  ! (out)
!        calc_ddens_perturbation, Ini_GP_PolyOrder_h, Ini_GP_PolyOrder_v, &  ! (in)
        calc_ddens_perturbation_tdash, Ini_GP_PolyOrder_h, Ini_GP_PolyOrder_v, &  ! (in)
        lcmesh, elem, RPlanet                                            )  ! (in)   

    case( 'COSBELL' )

      call mkinitutil_calc_cosinebell_global( PT_purtub,    & ! (out)
        DTHETA, rh, lonc, latc, RPlanet,                    & ! (in)
        x, y, z, lcmesh, elem,                              & ! (in)
        Ini_GP_PolyOrder_h, Ini_GP_PolyOrder_v,             & ! (in)
        'sin', (/ real(nv,kind=RP), dom_zmax - dom_zmin /)  ) ! (in)      

    end select
    
    !---
    RovCp = Rdry / CPdry

    !$omp parallel do private(EXNER)
    do ke=lcmesh%NeS, lcmesh%NeE
      EXNER(:) = ( PRES_hyd(:,ke) / PRES00 )**(RovCp)
      ! DDENS(:,ke) = PRES_hyd(:,ke) / ( Rdry * ( PRES_hyd(:,ke) / ( Rdry * DENS_hyd(:,ke) ) +  PT_purtub(:,ke) * EXNER(:) ) ) &
      !             - DENS_hyd(:,ke)
    end do

    return
  end subroutine exp_SetInitCond_inertia_gravity_wave

!OCL SERIAL
  subroutine calc_ddens_perturbation_tdash(  DDENS_intrp,   &
    lon_intrp, lat_intrp, z_intrp, elem_intrp, RPlanet )

    use scale_const, only: &
      GRAV => CONST_GRAV,   &
      Rdry => CONST_Rdry,   &
      CPdry => CONST_CPdry, &
      PRES00 => CONST_PRE00    
    implicit none
    class(ElementBase3D), intent(in) :: elem_intrp
    real(RP), intent(out) :: DDENS_intrp(elem_intrp%Np)
    real(RP), intent(in) :: lon_intrp(elem_intrp%Np)
    real(RP), intent(in) :: lat_intrp(elem_intrp%Np)
    real(RP), intent(in) :: z_intrp(elem_intrp%Np)
    real(RP), intent(in) :: RPlanet

    real(RP) :: r_intrp(elem_intrp%Np)

    real(RP) :: DT_intrp(elem_intrp%Np)
    real(RP) :: PT_hyd(elem_intrp%Np)
    real(RP) :: PRES_hyd(elem_intrp%Np)
    real(RP) :: DENS_hyd(elem_intrp%Np)
    real(RP) :: EXNER(elem_intrp%Np)
    !---------------------------------------------


    PT_hyd(:) = THETA0 * exp( BruntVaisalaFreq**2 / Grav * z_intrp(:) )
    
    EXNER(:) = 1.0_RP + Grav**2 / ( CpDry * BruntVaisalaFreq**2 ) * ( 1.0_RP / PT_hyd(:) - 1.0_RP / THETA0 )
    PRES_hyd(:) = PRES00 * EXNER(:)**(CpDry/Rdry)
    DENS_hyd(:) =  PRES_hyd(:) / ( Rdry * exner(:) * PT_hyd(:) )

    r_intrp(:) = RPlanet / rh * acos( sin(latc) * sin(lat_intrp(:)) + cos(latc) * cos(lat_intrp(:)) * cos(lon_intrp(:) - lonc) )

    DT_intrp(:) = DTHETA * exp( - r_intrp(:)**2 )           &
                   * sin( dble(nv) * PI * z_intrp(:) / Zt ) &
                   * exp( 0.5_RP * z_intrp(:) * Grav / (Rdry * THETA0 ) )

    DDENS_intrp(:) = PRES_hyd(:) / ( Rdry * THETA0 ) &
                  * ( - DT_intrp(:) / THETA0  ) / ( 1.0_RP +  DT_intrp(:) / THETA0 )
    return
  end subroutine calc_ddens_perturbation_tdash

end module mod_user
