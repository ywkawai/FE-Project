!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a Kelvin-Helmholtz wave experiment
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

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: MeshField3D
  
  use mod_user_base, only: UserBase
  use mod_experiment, only: Experiment
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(UserBase) :: User
  contains
    procedure :: mkinit_ => USER_mkinit
    generic :: mkinit => mkinit_
    procedure :: setup_ => USER_setup
    generic :: setup => setup_
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

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('Kelvin_Helmholtz_wave')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_Kelvin_Helmholtz_wave )
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

    call this%UserBase%Setup( atm, USER_do )

    return
  end subroutine USER_setup

  !------

!OCL SERIAL  
  subroutine exp_SetInitCond_Kelvin_Helmholtz_wave( this,                  &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,        &
      GRAV => CONST_GRAV,    &
      Rdry => CONST_Rdry,    &
      Rvap => CONST_Rvap,    &      
      CPdry => CONST_CPdry,  &
      CVdry => CONST_CVdry,  &
      PRES00 => CONST_PRE00, &
      Pstd   => CONST_Pstd
    use scale_random, only: &
      RANDOM_uniform      
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constPT, &
      hydrostaic_build_rho_XYZ
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
    
    ! Surface state
    real(RP) :: SFC_THETA
    real(RP) :: SFC_PRES
    ! Environment state
    real(RP) :: ENV_L1_ZTOP    = 1.9E3_RP ! top height of the layer1 (low THETA)  [m]
    real(RP) :: ENV_L3_ZBOTTOM = 2.1E3_RP ! bottom height of the layer3 (high THETA) [m]
    real(RP) :: ENV_L1_THETA   = 300.0_RP ! THETA in the layer1 (small THETA)     [K]
    real(RP) :: ENV_L3_THETA   = 301.0_RP ! THETA in the layer3 (large THETA)     [K]
    real(RP) :: ENV_L1_U       =   0.0_RP ! velocity u in the layer1 (low  THETA) [K]
    real(RP) :: ENV_L3_U       =  20.0_RP ! velocity u in the layer3 (high THETA) [K]    
    ! Disturbance
    real(RP) :: RANDOM_U     = 0.0_RP   ! amplitude of random disturbance of U [m]

    namelist /PARAM_EXP/ &
      SFC_PRES,          &      
      ENV_L1_ZTOP,       &
      ENV_L3_ZBOTTOM,    &
      ENV_L1_THETA,      &
      ENV_L3_THETA,      &
      ENV_L1_U,          &
      ENV_L3_U,          &
      RANDOM_U


    real(RP) :: fact(elem%Np)
    real(RP) :: rndm(elem%Np)   
    real(RP) :: PT_tmp(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: DENS(elem%Np)

    integer :: ke, p
    integer :: ke_x, ke_y, ke_z
    integer :: ierr
    !-----------------------------------------------------------------------------

    SFC_PRES  = Pstd

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("Kelvin_Helmholtz_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("Kelvin_Helmholtz_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    SFC_THETA = ENV_L1_THETA
    
    call hydrostatic_calc_basicstate_constPT( DENS_hyd, PRES_hyd, &
      SFC_THETA, SFC_PRES, x, y, z,  lcmesh, elem                 )
    
    !$omp parallel do collapse(3) private(ke,ke_x,ke_y,ke_z, fact)
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
    do ke_z=1, lcmesh%NeZ

      ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1)*lcmesh%Ne2D
      fact(:) = ( lcmesh%zlev(:,ke) - ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM - ENV_L1_ZTOP )
      fact(:) = max( min( fact(:), 1.0_RP ),  0.0_RP )

      PT_tmp(:,ke_z,ke_x,ke_y) = ENV_L1_THETA * ( 1.0_RP - fact(:) ) &
                               + ENV_L3_THETA * (          fact(:) )
    end do
    end do
    end do

    call hydrostaic_build_rho_XYZ( DDENS, & ! (out)
      DENS_hyd, PRES_hyd, PT_tmp,         & ! (in)
      x, y, z, lcmesh, elem               ) ! (in)
      
    !$omp parallel do collapse(3) private(DENS, ke,ke_x,ke_y,ke_z)
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
    do ke_z=1, lcmesh%NeZ
      ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1)*lcmesh%Ne2D

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      PRES_hyd(:,ke) = PRES00 * ( Rdry / PRES00 * DENS(:) * PT_tmp(:,ke_z,ke_x,ke_y) )**(CpDry/CvDry)
      DENS_hyd(:,ke) = DENS(:)
      DDENS(:,ke) = 0.0_RP
    end do
    end do
    end do

    !$omp parallel do private(fact, rndm, DENS)
    do ke=lcmesh%NeS, lcmesh%NeE
      fact(:) = ( lcmesh%zlev(:,ke) - ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM - ENV_L1_ZTOP )
      fact(:) = max( min(fact(:), 1.0_RP ),  0.0_RP )

      ! 
      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)

      call RANDOM_uniform( rndm )
      MOMX(:,ke) = DENS(:) * ( &
               ENV_L1_U * ( 1.0_RP - fact(:) )       &
             + ENV_L3_U * (          fact(:) )       &
             + ( rndm(:) * 2.0_RP - 1.0_RP ) * RANDOM_U  )
    end do

    return
  end subroutine exp_SetInitCond_Kelvin_Helmholtz_wave

end module mod_user
