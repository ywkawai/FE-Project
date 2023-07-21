!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!          a test case proposed by Galewsky et al. (2004)
!!          barotropic instability test
!!
!! @par Reference
!!  - J. Galewsky, R.K. Scott, L.M. Polvani, 2004:
!!    An initial-value problem for testing numerical models of the global shallow-water equations. 
!!    Tellus, Series A: Dynamic Meteorology and Oceanography, 56(5), 429-440.
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
    PI => CONST_PI,          &
    GRAV => CONST_GRAV,      &
    RPlanet => CONST_RADIUS, &
    OMG => CONST_OHM 
      
  use mod_exp, only: experiment
  use mod_globalsw_component, only: &
    GlobalSWComponent

  use scale_element_line, only: LineElement
  use scale_element_base, only: ElementBase2D
  use scale_localmesh_2d, only: LocalMesh2D   

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

  type, private, extends(experiment) :: Exp_G04_BarotropicInstability
  contains 
   procedure :: setInitCond_lc => exp_SetInitCond_G04_BarotropicInstability
  end type
  type(Exp_G04_BarotropicInstability), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  subroutine USER_mkinit ( swmodel )
    implicit none

    class(GlobalSWComponent), intent(inout) :: swmodel
    !------------------------------------------

    call exp_manager%Init('DB89_CrossPolarFlow')
    call exp_manager%SetInitCond( &
      swmodel%mesh, swmodel%vars%PROGVARS_manager, swmodel%vars%AUXVARS_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

  subroutine USER_setup( atm )
    implicit none
    
    class(GlobalSWComponent), intent(inout) :: atm

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
    return
  end subroutine USER_setup

  subroutine USER_calc_tendency( atm )
    implicit none
    class(GlobalSWComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_calc_tendency

  subroutine USER_update( atm )
    implicit none
    class(GlobalSWComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_update

  !------
  subroutine exp_SetInitCond_G04_BarotropicInstability( this,  &
    h, U, V, hs, u1, u2,                                 &
    x, y, lcmesh, elem                                   )
    
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSVec
    implicit none

    class(Exp_G04_BarotropicInstability), intent(inout) :: this
    type(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: h(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: U(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: V(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: hs(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: u1(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: u2(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)  

    real(RP) :: lat0, lat1
    real(RP) :: umax
    real(RP) :: H0

    real(RP) :: hhat
    real(RP) :: phi2
    real(RP) :: alph
    real(RP) :: beta

    namelist /PARAM_EXP/ &
      lat0, lat1, umax,      &
      hhat, phi2, alph, beta
    
  
    real(RP) :: VelLon(elem%Np,lcmesh%Ne)
    real(RP) :: VelLat(elem%Np,lcmesh%Ne)
    integer :: ke
    integer :: ierr

    type(LineElement) :: elem1D
    integer :: p

    real(RP) :: lon_(elem%Np)
    !-----------------------------------------------------------------------------

    lat0 = PI / 7.0_RP
    lat1 = 0.5_RP * PI - lat0
    umax = 80.0_RP
    H0 = 10.E3_RP

    hhat = 120.0_RP
    phi2 = 0.25_RP * PI
    alph = 1.0_RP / 3.0_RP
    beta = 1.0_RP / 15.0_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("DB89_CrossPolarFlow_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("DB89_CrossPolarFlow_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    call elem1D%Init( 12, .false. )

    !$omp parallel do private (p, lon_)
    do ke=lcmesh%NeS, lcmesh%NeE

      call cal_zonal_vel( VelLon(:,ke),             &
        lcmesh%lat(:,ke), elem%Np, lat0, lat1, umax )
      VelLat(:,ke) = 0.0_RP    
                  
      do p=1, elem%Np
        call cal_height_variation( h(p,ke), lcmesh%lat(p,ke), elem1D, lat0, lat1, umax )
        h(p,ke) = H0 + h(p,ke)
      end do

      lon_(:) = lcmesh%lon(:,ke)
      where ( lon_(:) > PI ) 
        lon_(:) = lon_(:) - 2.0_RP * PI
      end where
      where ( abs(lon_(:)) < PI )
        h(:,ke) = h(:,ke) &
          + hhat * cos(lcmesh%lat(:,ke)) * exp(- (lon_(:)/alph)**2) &
            * exp(- ((phi2 - lcmesh%lat(:,ke))/beta)**2 )
      end where
      hs(:,ke) = 0.0_RP
    end do

    call CubedSphereCoordCnv_LonLat2CSVec( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), elem%Np * lcmesh%Ne, RPlanet, & ! (in)
      VelLon(:,:), VelLat(:,:),                                                                 & ! (in)
      U(:,lcmesh%NeS:lcmesh%NeE), V(:,lcmesh%NeS:lcmesh%NeE)                                    ) ! (out)
    !$omp parallel do
    do ke=lcmesh%NeS, lcmesh%NeE
      u1(:,ke) = lcmesh%G_ij(:,ke,1,1) * U(:,ke) + lcmesh%G_ij(:,ke,1,2) * V(:,ke)
      u2(:,ke) = lcmesh%G_ij(:,ke,2,1) * U(:,ke) + lcmesh%G_ij(:,ke,2,2) * V(:,ke)
    end do

    !--
    call elem1D%Final()

    return
  end subroutine exp_SetInitCond_G04_BarotropicInstability

  subroutine cal_height_variation( height, &
      lat, elem1D, lat0, lat1, umax )
    implicit none
    real(RP), intent(out) :: height
    real(RP), intent(in) :: lat
    type(LineElement), intent(in) :: elem1D
    real(RP), intent(in) :: lat0, lat1
    real(RP), intent(in) :: umax

    integer :: ke
    integer :: ke_max
    real(RP) :: dlat_e_default
    real(RP), allocatable :: lat_e(:)
    real(RP) :: lat_lc(elem1D%Np)
    real(RP) :: u_lc(elem1D%Np)
    !--------------------------------------

    dlat_e_default = 2.5_RP / 180.0_RP * PI
    ke_max = int( abs(lat) / dlat_e_default )
    if ( ke_max * dlat_e_default < abs(lat) ) ke_max = ke_max + 1
    
    allocate( lat_e(ke_max+1) )
    lat_e(1) = 0.0_RP  
    do ke=2, ke_max
      lat_e(ke) = lat_e(ke-1) + sign(dlat_e_default, lat)
    end do
    lat_e(ke_max+1) = lat

    !
    height = 0.0_RP
    do ke=1, ke_max
      lat_lc(:) =  lat_e(ke) + 0.5_RP * ( lat_e(ke+1) - lat_e(ke) ) * ( 1.0_RP + elem1D%x1(:) ) 
      call cal_zonal_vel( u_lc(:), lat_lc(:), elem1D%Np, lat0, lat1, umax )

      height = height &
         - 0.5_RP * ( lat_e(ke+1) - lat_e(ke) ) * RPlanet / Grav                        &
          * sum(  elem1D%IntWeight_lgl(:) *  u_lc(:) * (                                &
                   2.0_RP * OMG * sin(lat_lc(:)) + u_lc(:) * tan(lat_lc(:)) / RPlanet ) )
    end do

    return
  end subroutine cal_height_variation

  subroutine cal_zonal_vel( u, lat, Np, lat0, lat1, umax )
    integer, intent(in) :: Np
    real(RP), intent(out) :: u(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: lat0, lat1
    real(RP), intent(in) :: umax

    real(RP) :: en
    !---------------------------------------------------

    en = exp( - 4.0_RP / (lat1 - lat0)**2 )

    where( lat0 < lat(:) .and. lat(:) < lat1 )
      u(:) = umax / en * exp( 1.0_RP / ( ( lat(:) - lat0 ) * ( lat(:) - lat1 ) ) )
    elsewhere
      u(:) = 0.0_RP
    end where

    return
  end subroutine cal_zonal_vel

end module mod_user
