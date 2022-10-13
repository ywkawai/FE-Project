!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a benchmark test of dynamical core described in Menou & Rauscher (2009) for Earth-like atmospheric simulation. 
!!          This is a simplified variant of the Held & Suarez (1994) test. Different from the HS test, relaxation timescales in the radiative forcing and 
!!          is constant values and the level of tropopause is explicity specified. This experimental setup is considered to be useful to explore qualitative features 
!!          of atmospheric flows on generalized Earth-like planet. 
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
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  & 
    CVdry => CONST_CVdry,  &      
    PRES00 => CONST_PRE00, &
    Grav => CONST_GRAV,    &
    PI =>   CONST_PI
  
  use mod_exp, only: experiment

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_line, only: LineElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_mesh_base3d, only: MeshBase3D
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

  type, private, extends(experiment) :: Exp_MR2009_EarthLike
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_MR2009_EarthLike
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_MR2009_EarthLike), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  real(RP), private, parameter :: kf              = 1.0_RP / ( 86400.0_RP * 1.0_RP  )
  real(RP), private, parameter :: krad            = 1.0_RP / ( 86400.0_RP * 15.0_RP )
  real(RP), private, parameter :: LAPSE_RATE_trop = 6.5E-3_RP !< Lapse rate of troposphere
  real(RP), private, parameter :: Zstrato         = 1.2E4_RP  !< Specified height of tropopause 
  real(RP), private, parameter :: TEMP_strato     = 212.0_RP  !< Surface temperature at equator associated with radiative forcing
  real(RP), private, parameter :: DTEMP_strato    =   2.0_RP  !< Tropopause temperature increment
  real(RP), private, parameter :: SFCTEMP_eq      = 288.0_RP  !< Surface temperature at equator associated with radiative forcing
  real(RP), private, parameter :: DelT_y          = 60.0_RP   !< Equator-to-pole temperature difference associated with radiative forcing
  real(RP), private, parameter :: sigb            = 0.7_RP

  !-----------------------------------------------------------------------------

contains

!OCL SERIAL
  subroutine USER_mkinit ( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    call exp_manager%Init('Earth-like atmosphere')

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

    return
  end subroutine USER_setup

!OCL SERIAL  
  subroutine USER_calc_tendency( atm )  
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_calc_tendency

!OCL SERIAL
  subroutine USER_update( atm )
    use scale_localmeshfield_base, only: LocalMeshFieldBase

    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      MOMX_p  => PHYTEND_MOMX_ID, &
      MOMY_p  => PHYTEND_MOMY_ID, &
      MOMZ_p  => PHYTEND_MOMZ_ID, &
      RHOH_p  => PHYTEND_RHOH_ID

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars

    implicit none

    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    type(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke1, ke2
    integer :: ke2D, ke_z
    integer :: p2D, p_z

    real(RP), allocatable :: DENS(:), T(:), Teq(:), sig(:), PRES_sfc(:)
    real(RP), allocatable :: rtauV(:)
    real(RP), allocatable :: lat(:), z(:)
    real(RP), allocatable :: Tvert(:), Beta_trop(:)
    real(RP), allocatable :: sig_strat(:)

    integer :: ke_z2, p_z2
    real(RP) :: w1_zstrat
    integer :: ke1_zstrat, p1_zstrat 
    integer :: ke2_zstrat, p2_zstrat
    real(RP), allocatable :: sig_strat2D(:,:)
    real(RP), allocatable :: zlev_col(:,:)

    real(DP) :: dt
    real(RP) :: Gamm
    real(RP) :: SFCTEMP0
    !----------------------------------------------------------

    dt = atm%time_manager%dtsec
    gamm = CpDry / CvDry 
    SFCTEMP0 = TEMP_strato + Zstrato * LAPSE_RATE_trop - DTEMP_strato

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      

      call AtmosVars_GetLocalMeshPhyAuxVars( n, atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                         )
      
      elem3D => lcmesh%refElem3D

      allocate( DENS(elem3D%Np), T(elem3D%Np), Teq(elem3D%Np), sig(elem3D%Np) )
      allocate( PRES_sfc(elem3D%Nnode_h1D**2) )
      allocate( rtauV(elem3D%Np) )
      allocate( lat(elem3D%Np), z(elem3D%Np) )
      allocate( Tvert(elem3D%Np), Beta_trop(elem3D%Np) )
      allocate( sig_strat(elem3D%Np) )
      allocate( sig_strat2D(elem3D%Nnode_h1D**2,lcmesh%Ne2D) )
      allocate( zlev_col(elem3D%Nnode_v,lcmesh%NeZ) )
      
      !--

      do ke2D=1, 1
      do p2D=1, 1
        do ke_z=1, lcmesh%NeZ
          ke = ke2D + (ke_z-1)*lcmesh%Ne2D
          zlev_col(:,ke_z) = lcmesh%zlev(elem3D%Colmask(:,p2D),ke)
        end do
        zloop: do ke_z=1, lcmesh%NeZ
          p_z2 = 2; ke_z2 = ke_z
          do p_z=1, elem3D%Nnode_v
            if ( ( zlev_col(p_z,ke_z) - Zstrato ) * ( zlev_col(p_z2,ke_z2) - Zstrato ) <= 0.0_RP )  then
              ke1_zstrat = ke_z; ke2_zstrat = ke_z2
              p1_zstrat = p_z; p2_zstrat = p_z2
              w1_zstrat = ( zlev_col(p_z2,ke_z2) - Zstrato ) / ( zlev_col(p_z2,ke_z2) - zlev_col(p_z,ke_z) )
              exit zloop
            end if
          end do

          p_z2 = p_z2 + 1
          if ( p_z2 == elem3D%Nnode_v ) then
            p_z2 = 2; ke_z2 = ke_z + 1
          end if
          if ( ke_z2 > lcmesh%NeZ ) then
            LOG_ERROR("USER_update",*) 'Altitidue of Zstrato is outside the vertical computational domain. Check!'
            call PRC_abort
          end if
        end do zloop
      end do
      end do

      !$omp parallel do private(ke2D, ke1, ke2)
      do ke2D=1, lcmesh%Ne2D
        ke1 = ke2D +(ke1_zstrat - 1)*lcmesh%Ne2D
        ke2 = ke2D +(ke2_zstrat - 1)*lcmesh%Ne2D
        sig_strat2D(:,ke2D) = &
        (              w1_zstrat   * PRES%val(elem3D%Hslice(:,p1_zstrat),ke1)   &
          + ( 1.0_RP - w1_zstrat ) * PRES%val(elem3D%Hslice(:,p2_zstrat),ke2) ) &
        / PRES%val(elem3D%Hslice(:,1),ke2D)
      end do

      !$omp parallel do private( &
      !$omp DENS, T, Teq, PRES_sfc, sig, lat, z, rtauV, ke2D, &
      !$omp Tvert, Beta_trop, sig_strat )
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)

        PRES_sfc(:) = PRES%val(elem3D%Hslice(:,1),ke2D)
        sig(:) = PRES%val(:,ke) / PRES_sfc(elem3D%IndexH2Dto3D)
        lat(:) = lcmesh%lat2D(elem3D%IndexH2Dto3D,ke2D)
        z(:) = lcmesh%zlev(:,ke)

        ! Calculate radiative equilirbium temperature

        where ( z(:) <= Zstrato )
          Tvert(:) = &
                SFCTEMP0 - LAPSE_RATE_trop * ( Zstrato + 0.5_RP * (z(:) - Zstrato) )       &
              + sqrt( ( 0.5_RP * LAPSE_RATE_trop * (z(:) - Zstrato) )**2 + DTEMP_strato**2 )
        elsewhere
          Tvert(:) = TEMP_strato 
        end where

        sig_strat(:) = sig_strat2D(elem3D%IndexH2Dto3D(:),ke2D)
        Beta_trop(:) = max( 0.0_RP, sin( 0.5_RP * PI * (sig(:) - sig_strat(:))/(1.0_RP - sig_strat(:)) ) )

        Teq(:) = Tvert(:) + Beta_trop(:) * DelT_y * ( 1.0_RP / 3.0_RP - sin(lat(:))**2 )


        ! Adjust momentum and temperature

        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)
        T(:) =  PRES%val(:,ke) / ( Rdry * DENS(:) )

        rtauV(:) = kf * max( 0.0_RP, (sig(:) - sigb)/(1.0_RP - sigb) )
        MOMX%val(:,ke) = MOMX%val(:,ke) / ( 1.0_RP + dt * rtauV )
        MOMY%val(:,ke) = MOMY%val(:,ke) / ( 1.0_RP + dt * rtauV )

        !- For the case of d DRHOT /dt = dens * Cp * ( Teq - T ) / tauT     
        !  <- It is based on the forcing form in Held and Surez in which the temperature evolution equation is assumed to be dT/dt = R/Cp * T/p * dp/dt + (Teq - T) / tauT
        DRHOT%val(:,ke) = DRHOT%val(:,ke) &
                        - dt * krad * ( 1.0_RP - Teq(:) / T(:) ) * DENS(:) * PT%val(:,ke)     &
                        / ( 1.0_RP + dt * krad * ( 1.0_RP + (gamm - 1.0_RP) * Teq(:) / T(:) ) )  
      end do

      deallocate( DENS, T, Teq, sig, PRES_sfc )
      deallocate( rtauV )
      deallocate( lat, z )
      deallocate( Tvert, Beta_trop )
      deallocate( sig_strat, sig_strat2D )
      deallocate( zlev_col )
    end do
    
    return
  end subroutine USER_update

  !------

!OCL SERIAL
  subroutine exp_SetInitCond_MR2009_EarthLike( this,                            &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use mod_exp, only: &
      TracerLocalMeshField_ptr

    use scale_const, only: &
      PRES00 => CONST_PRE00, &
      RPlanet => CONST_Radius, &
      CPdry => CONST_CPdry, &     
      CVdry => CONST_CVdry

    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constT
  
    implicit none

    class(Exp_MR2009_EarthLike), intent(inout) :: this
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
    
    real(RP) :: TEMP0           = 264.0_RP !< Temperature of initial isothermal atmosphere based on Table 1 of Heng et al. (2011,  Monthly Notices of the Royal Astronomical Society)
    namelist /PARAM_EXP/ &
      TEMP0

    integer :: ierr
    !-----------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("Held_Suarez_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("Held_Suarez_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    call hydrostatic_calc_basicstate_constT( DENS_hyd, PRES_hyd, & ! (out)
      TEMP0, PRES00,                                             & ! (in)
      x, y, z, lcmesh, elem                                      ) ! (in)

    return
  end subroutine exp_SetInitCond_MR2009_EarthLike

!OCL SERIAL
  subroutine exp_geostrophic_balance_correction( this,  &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, &
    lcmesh, elem )
    
    implicit none

    class(Exp_MR2009_EarthLike), intent(inout) :: this
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
