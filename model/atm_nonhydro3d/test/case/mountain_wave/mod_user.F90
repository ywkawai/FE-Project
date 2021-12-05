!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
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
  use mod_exp, only: experiment

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: MeshField3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_mkinit
  public :: USER_setup
  public :: USER_calc_tendency
  public :: USER_update

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

  type, private, extends(experiment) :: Exp_mountain_wave
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_mountain_wave
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_mountain_wave), private :: exp_manager

  logical, private :: USER_do        = .false. !< do user step?
  real(RP), private :: U0            
  real(RP), private :: zTop
  real(RP), private :: SPONGE_HEIGHT

  type(MeshField3D), private :: PT_diff

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit ( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    call exp_manager%Init('mountain_wave')

    call exp_manager%SetInitCond( atm%mesh,                &
      atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager, &
      atm%vars%QTRCVARS_manager                            )
    
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL  
  subroutine USER_setup( atm )
    implicit none
    
    class(AtmosComponent), intent(inout) :: atm

    namelist / PARAM_USER / &
       USER_do,             &
       U0,                  &
       zTop, Sponge_Height

    integer :: ierr    
    !------------------------------------------


    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    U0             = 10.0_RP
    zTop           = 30.E3_RP
    SPONGE_HEIGHT  = 15.E3_RP

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
    if ( USER_do ) call PT_diff%Init( 'PT_diff', 'K', atm%mesh%ptr_mesh )

    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_calc_tendency( atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in

    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CpDry => CONST_CPdry,    &
      RPlanet => CONST_RADIUS, &
      PI => CONST_PI

    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars, &
      MOMX_p  => ATMOS_PHYTEND_MOMX_ID, &
      MOMY_p  => ATMOS_PHYTEND_MOMY_ID, &
      MOMZ_p  => ATMOS_PHYTEND_MOMZ_ID, &
      RHOH_p  => ATMOS_PHYTEND_RHOH_ID
    use scale_localmeshfield_base, only: LocalMeshFieldBase

    implicit none

    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(ElementBase3D), pointer :: elem
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT

    real(RP), parameter :: rtau = 1.0_RP / 600.0_RP

    integer :: n
    integer :: ke

    real(RP), allocatable :: DENS(:)
    real(RP), allocatable :: sfac(:)
    !------------------------------------------

    if ( USER_do ) then
      call atm%vars%Calc_diagVar( 'PT_diff', PT_diff )
      call FILE_HISTORY_meshfield_in( PT_diff, "perturbation of potential temperature" )
    end if
    
    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, lcmesh                               )      
      
      call AtmosVars_GetLocalMeshPhyAuxVars( n,  atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                          )
      
      elem => lcmesh%refElem3D
      allocate( DENS(elem%Np), sfac(elem%Np) )

      !$omp parallel do private(DENS, sfac)
      do ke=lcmesh%NeS, lcmesh%NeE
        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)
        sfac(:) = rtau * 0.5_RP * ( 1.0_RP - cos( PI * ( lcmesh%pos_en(:,ke,3) - SPONGE_HEIGHT ) / ( zTop - SPONGE_HEIGHT ) ) ) 

        where ( lcmesh%pos_en(:,ke,3) > SPONGE_HEIGHT ) 
          atm%vars%PHY_TEND(MOMX_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMX_p)%local(n)%val(:,ke)   &
            - sfac(:) * ( MOMX%val(:,ke) - DENS(:) * U0 )
          atm%vars%PHY_TEND(MOMY_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMY_p)%local(n)%val(:,ke)   &
            - sfac(:) * MOMY%val(:,ke)
          atm%vars%PHY_TEND(MOMZ_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMZ_p)%local(n)%val(:,ke)   &
            - sfac(:) * MOMZ%val(:,ke)

          atm%vars%PHY_TEND(RHOH_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(RHOH_p)%local(n)%val(:,ke)   &
            - DENS(:) * sfac(:) * CpDry * ( PRES%val(:,ke) / DENS(:) - PRES_hyd%val(:,ke) / DENS_hyd%val(:,ke) ) / Rdry
        end where
      end do
      deallocate( DENS, sfac )
    end do

    return
  end subroutine USER_calc_tendency

  subroutine USER_update( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_update

  !------

!OCL SERIAL  
  subroutine exp_SetInitCond_mountain_wave( this,                 &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,        &
      GRAV => CONST_GRAV,    &
      OHM => CONST_OHM,      &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      PRES00 => CONST_PRE00, &
      RPlanet => CONST_RADIUS
    
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constBVFreq
    
    use mod_exp, only: &
      TracerLocalMeshField_ptr
    
    implicit none

    class(Exp_mountain_wave), intent(inout) :: this
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
    
    real(RP) :: THETA0 = 280.0_RP
    real(RP) :: BruntVaisalaFreq = 0.01_RP  
    namelist /PARAM_EXP/ &
      U0,                &
      BruntVaisalaFreq,  &
      THETA0
    integer, parameter :: IntrpPolyOrder_h = 8
    integer, parameter :: IntrpPolyOrder_v = 8

    integer :: ke
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
    !---

    call hydrostatic_calc_basicstate_constBVFreq( DENS_hyd, PRES_hyd, & ! (out)
      BruntVaisalaFreq, THETA0, PRES00,                               & ! (in)
      x, y, lcmesh%zlev, lcmesh, elem                                 ) ! (in)
    !---

    !$omp parallel do
    do ke=lcmesh%NeS, lcmesh%NeE
      MOMX(:,ke) = DENS_hyd(:,ke) * U0
      MOMY(:,ke) = 0.0_RP
    end do

    return
  end subroutine exp_SetInitCond_mountain_wave

  subroutine exp_geostrophic_balance_correction( this,  &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, &
    lcmesh, elem )
    
    implicit none

    class(Exp_mountain_wave), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(inout) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DRHOT(elem%Np,lcmesh%NeA)

    !---------------------------------------------------
    return
  end subroutine exp_geostrophic_balance_correction 

end module mod_user
