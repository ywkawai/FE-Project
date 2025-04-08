!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!    Set the initial data for baroclinic wave in global model based on Jablonowski and Williamson (2006).
!!
!! @author Yuta Kawai, Team SCALE
!!
!! @par Reference
!!  - Jablonowski, C. and Williamson, D.L., 2006:
!!    A baroclinic instability test case for atmospheric model dynamical cores. 
!!    Q.J.R. Meteorol. Soc., 132, 2943-2975.
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
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
    EPS => CONST_EPS,        &  
    PI => CONST_PI,          &
    GRAV => CONST_GRAV,      &
    Rdry => CONST_Rdry,      &
    CPdry => CONST_CPdry,    &
    CVdry => CONST_CVdry,    &
    PRES00 => CONST_PRE00,   &
    RPlanet => CONST_RADIUS, &
    OHM     => CONST_OHM
 
  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D  

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

  ! Parameters for inital stratification
  real(RP), private :: REF_TEMP    = 288.E0_RP !< The reference temperature [K]
  real(RP), private :: REF_PRES    = 1.E5_RP   !< The reference pressure [Pa]
  real(RP), private :: LAPSE_RATE  = 5.E-3_RP  !< The lapse rate [K/m]

  ! Parameters for background zonal jet
  real(RP), private :: U0 = 35.E0_RP           !< The parameter associated with zonal jet maximum amplitude  [m/s]

  ! Parameters for inital perturbation of zonal wind with a Gaussian profile
  !
  real(RP), private :: Up     = 1.E0_RP         !< The maximum amplitude of zonal wind perturbation [m/s]
  real(RP), private :: Lp                       !< The width of Gaussian profile
  real(RP), private :: Lon_c                    !< The center point (x) of inital perturbation
  real(RP), private :: Lat_c                    !< The center point (y) of inital perturbation

  real(RP), private :: ETA0 = 0.252_RP
  real(RP), private :: ETAt = 0.2_RP
  real(RP), private :: DELTAT = 4.8E5_RP

  ! Polynomial order used in Galerkin projection for initial data
  !
  integer, private :: IniIntrpPolyOrder_h = 8
  integer, private :: IniIntrpPolyOrder_v = 8


  ! SPONGE Layer
  ! logical :: SPONGE_FLAG = .false.
  ! real(RP) :: SPONGE_sigR = 0.7_RP
  ! real(RP) :: SPONGE_TauR = 600.0_RP !86400.0_RP

  ! Flag with topography
  logical :: skip_topo = .false.

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init( 'baroclinic_wave_global' )
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_baroclinicwave )
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

    return
  end subroutine USER_setup

  !------

!OCL SERIAL  
  subroutine exp_SetInitCond_baroclinicwave( this, &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use mod_mkinit_util, only: &
      mkinitutil_gen_GPMat
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
    
    
    namelist /PARAM_EXP/ &
      REF_TEMP, REF_PRES, LAPSE_RATE, &
      U0,                             &
      Up, Lp, Lon_c, Lat_c,           &
      IniIntrpPolyOrder_h,            &
      IniIntrpPolyOrder_v,            &
      skip_topo

    integer :: ke

    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: IntrpMat(:,:)

    real(RP), allocatable :: PRES_balance(:,:)
    real(RP), allocatable :: TEMP_balance(:,:)
    real(RP), allocatable :: U_balance(:,:)
    real(RP), allocatable :: V_balance(:,:)
    real(RP), allocatable :: U_dash(:,:)
    real(RP), allocatable :: V_dash(:,:)
    real(RP), allocatable :: DENS_ip(:)
    integer :: ierr
    !-----------------------------------------------------------------------------

    Lp = RPlanet / 10.0_RP
    Lon_c = PI / 9.0_RP
    Lat_c = 2.0_RP / 9.0_RP * PI

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("BAROCLINIC_WAVE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("BAROCLINIC_WAVE_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    
    call elem_intrp%Init( IniIntrpPolyOrder_h, IniIntrpPolyOrder_v, elem%IsLumpedMatrix() )

    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    call mkinitutil_gen_GPMat( IntrpMat, elem_intrp, elem )

    !---

    allocate( PRES_balance(elem_intrp%Np,lcmesh%Ne), TEMP_balance(elem_intrp%Np,lcmesh%Ne) ) 
    allocate( U_balance(elem_intrp%Np,lcmesh%Ne), V_balance(elem_intrp%Np,lcmesh%Ne) ) 
    allocate( U_dash(elem_intrp%Np,lcmesh%Ne), V_dash(elem_intrp%Np,lcmesh%Ne) ) 
    allocate( DENS_ip(elem_intrp%Np) )

    call calc_balanced_field( &
      PRES_balance, TEMP_balance, U_balance, V_balance, U_dash, V_dash,                  &
      lcmesh, elem_intrp%x1, elem_intrp%x2, elem_intrp%x3, elem_intrp%Np, elem_intrp%Nv, &
      dom_zmax )
    
    !$omp parallel do private( ke, DENS_ip )
    do ke = lcmesh%NeS, lcmesh%NeE
      DENS_ip(:) = PRES_balance(:,ke) / ( Rdry * TEMP_balance(:,ke) )

      PRES_hyd(:,ke) = matmul( IntrpMat, PRES_balance(:,ke) )
      DENS_hyd(:,ke) = matmul( IntrpMat, DENS_ip(:) )

      MOMX(:,ke) = matmul( IntrpMat, DENS_ip(:) * ( U_balance(:,ke) + U_dash(:,ke) ) )
      MOMY(:,ke) = matmul( IntrpMat, DENS_ip(:) * ( V_balance(:,ke) + V_dash(:,ke) ) )
    end do

    !------
    call elem_intrp%Final()

    return
  end subroutine exp_SetInitCond_baroclinicwave

!OCL SERIAL    
  subroutine calc_balanced_field( &
    PRES, TEMP, U, V, U_dash, V_dash,                   & ! (out)
    lcmesh3D, elem_x1, elem_x2, elem_x3, Np, Nv,        & ! (in)
    dom_zmax,                                           & ! (in)
    lat_         ) ! (in)
    
    use scale_const, only: &
      RPlanet => CONST_RADIUS, &
      EPS => CONST_EPS
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatPos, &
      CubedSphereCoordCnv_LonLat2CSVec
    use mod_mktopo_util, only: &
      calc_topo => mktopoutil_barocwave_global_JW2006_calc_topo
    
    implicit none
    class(LocalMesh3D), intent(in), target :: lcmesh3D
    integer, intent(in) :: Np, Nv
    real(RP), intent(out) :: PRES(Np,lcmesh3D%NeA)
    real(RP), intent(out) :: TEMP(Np,lcmesh3D%NeA)
    real(RP), intent(out) :: U(Np,lcmesh3D%NeA)
    real(RP), intent(out) :: V(Np,lcmesh3D%NeA)
    real(RP), intent(out) :: U_dash(Np,lcmesh3D%NeA)
    real(RP), intent(out) :: V_dash(Np,lcmesh3D%NeA)
    real(RP), intent(in) :: elem_x1(Np)
    real(RP), intent(in) :: elem_x2(Np)
    real(RP), intent(in) :: elem_x3(Np)
    real(RP), intent(in) :: dom_zmax
    real(RP), intent(out), optional :: lat_(Np,lcmesh3D%Ne)

    real(RP) :: vx(Nv), vy(Nv), vz(Nv)
    real(RP) :: alpha(Np,lcmesh3D%Ne), beta(Np,lcmesh3D%Ne), eta(Np,lcmesh3D%Ne)
    real(RP) :: zlev(Np,lcmesh3D%Ne)
    real(RP) :: lon(Np,lcmesh3D%Ne), lat(Np,lcmesh3D%Ne)
    real(RP) :: gam(Np,lcmesh3D%Ne)
    real(RP) :: VelLon(Np,lcmesh3D%Ne), VelLat(Np,lcmesh3D%Ne)
    real(RP) :: VelLon_dash(Np,lcmesh3D%Ne), VelLat_dash(Np,lcmesh3D%Ne)

    real(RP) :: geopot_hvari1(Np,lcmesh3D%Ne2D), geopot_hvari2(Np,lcmesh3D%Ne2D)
    real(RP) :: r_intrp(Np) 

    real(RP) :: topo(Np,lcmesh3D%Ne2D)

    integer :: ke, ke2D
    integer :: i, j, k
    integer :: p

    class(LocalMesh2D), pointer :: lcmesh2D
    !---------------------------------------------

    lcmesh2D => lcmesh3D%lcmesh2D

    !$omp parallel do private(vx, vy, vz)
    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      vx(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),1)
      vy(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),2)
      vz(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),3)
      alpha(:,ke) = vx(1) + 0.5_RP * ( elem_x1(:) + 1.0_RP ) * ( vx(2) - vx(1) ) 
      beta (:,ke) = vy(1) + 0.5_RP * ( elem_x2(:) + 1.0_RP ) * ( vy(4) - vy(1) )
      eta  (:,ke) = vz(1) + 0.5_RP * ( elem_x3(:) + 1.0_RP ) * ( vz(5) - vz(1) )
      gam  (:,ke) = 1.0_RP
    end do

    call CubedSphereCoordCnv_CS2LonLatPos( lcmesh3D%panelID, alpha, beta, gam, Np * lcmesh3D%Ne, &
      lon(:,:), lat(:,:) )
    
    if ( present(lat_) ) then
      do ke=lcmesh3D%NeS, lcmesh3D%NeE
        lat_(:,ke) = lat(:,ke)
      end do
    end if
    
    !$omp parallel do private(ke2D)
    do ke2D = lcmesh2D%NeS, lcmesh2D%NeE
      call get_thermal_wind_balance_1point_geopot_hvari( geopot_hvari1(:,ke2D), geopot_hvari2(:,ke2D), &
        lat(:,ke2D), Np )

        if (skip_topo) then
          topo(:,ke2D) = 0.0_RP
        else
          call calc_topo( topo(:,ke2D), &
            U0, ETA0, lat(:,ke2D), Np )
        end if  
    end do

    !$omp parallel do private(ke2D)
    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      ke2D = lcmesh3D%EMap3Dto2D(ke)

      zlev(:,ke) = eta(:,ke) &
                 + ( 1.0_RP - eta(:,ke) / dom_zmax ) * topo(:,ke2D)
    end do

    !$omp parallel do collapse(3) private(i,j,k,p,ke,ke2D, r_intrp)
    do k = 1, lcmesh3D%NeZ
    do j = 1, lcmesh3D%NeY
    do i = 1, lcmesh3D%NeX
      ke = i + (j-1)*lcmesh3D%NeX + (k-1)*lcmesh3D%NeX*lcmesh3D%NeY
      ke2D = lcmesh3D%EMap3Dto2D(ke)

      do p=1, Np
        call get_thermal_wind_balance_1point_itr( &
          pres(p,ke), temp(p,ke), VelLon(p,ke),                                        &
          lon(p,ke), lat(p,ke), zlev(p,ke), geopot_hvari1(p,ke2D), geopot_hvari2(p,ke2D), &
          REF_TEMP, REF_PRES, LAPSE_RATE,                                              &
          p, ke )
      end do
      VelLat(:,ke) = 0.0_RP

      r_intrp(:) = RPlanet / Lp * acos( sin(Lat_c) * sin(lat(:,ke2D)) + cos(Lat_c) * cos(lat(:,ke2D)) * cos(lon(:,ke2D) - Lon_c) )
      VelLon_dash(:,ke) = Up * exp( - r_intrp(:)**2 )
      VelLat_dash(:,ke) = 0.0_RP
      do p=1, Np
        if (k==1 .and. cos(lat(p,ke2D)) < EPS  ) then
          LOG_INFO("Replace VelLon with zero: (lat,VelLon)=", *) lat(p,ke2D), ":", VelLon(p,ke)
          VelLon(p,ke) = 0.0_RP
        end if
      end do
    end do
    end do
    end do

    call CubedSphereCoordCnv_LonLat2CSVec( &
      lcmesh3D%panelID, alpha, beta, gam, Np * lcmesh3D%Ne,             & ! (in)
      VelLon(:,:), VelLat(:,:),                                         & ! (in)
      U(:,lcmesh3D%NeS:lcmesh3D%NeE), V(:,lcmesh3D%NeS:lcmesh3D%NeE)    ) ! (out)

    call CubedSphereCoordCnv_LonLat2CSVec( &
      lcmesh3D%panelID, alpha, beta, gam, Np * lcmesh3D%Ne,                    & ! (in)
      VelLon_dash(:,:), VelLat_dash(:,:),                                      & ! (in)
      U_dash(:,lcmesh3D%NeS:lcmesh3D%NeE), V_dash(:,lcmesh3D%NeS:lcmesh3D%NeE) ) ! (out)

    return
  end subroutine calc_balanced_field 

!OCL SERIAL
  subroutine get_thermal_wind_balance_1point_geopot_hvari( geopot_hvari1, geopot_hvari2, &
    lat, Np )
    implicit none

    integer, intent(in) :: Np
    real(RP), intent(out) :: geopot_hvari1(Np)
    real(RP), intent(out) :: geopot_hvari2(Np)
    real(RP), intent(in) :: lat(Np)

    real(RP) :: sin_lat(Np)
    real(RP) :: cos_lat(Np)
    !-------------------------------------------

    sin_lat(:) = sin(lat(:))
    cos_lat(:) = cos(lat(:))

    ! Calc horizontal variation of geopotential height
    geopot_hvari1(:) = U0 * ( - 2.0_RP * sin_lat(:)**6 * ( cos_lat(:)**2 + 1.0_RP / 3.0_RP ) + 10.0_RP / 63.0_RP )
    geopot_hvari2(:) = RPlanet * OHM * ( 8.0_RP / 5.0_RP * cos_lat(:)**3 * ( sin_lat(:)**2 + 2.0_RP / 3.0_RP ) - 0.25_RP * PI )

    return
  end subroutine get_thermal_wind_balance_1point_geopot_hvari

  subroutine get_thermal_wind_balance_1point_itr( &
    pres, temp, vel_lon,                          &
    lon, lat, z, geopot_hvari1, geopot_hvari2,    &
    REF_TEMP, REF_PRES, LAPSE_RATE,               &
    p_, ke )

    implicit none

    real(RP), intent(out) :: pres
    real(RP), intent(out) :: temp
    real(RP), intent(out) :: vel_lon
    real(RP), intent(in) :: lon
    real(RP), intent(in) :: lat
    real(RP), intent(in) :: z
    real(RP), intent(in) :: geopot_hvari1
    real(RP), intent(in) :: geopot_hvari2
    real(RP), intent(in) :: REF_TEMP, REF_PRES, LAPSE_RATE
    integer, intent(in) :: p_, ke

    integer :: itr
    real(RP) :: eta, eta_save, etav, del_eta
    real(RP) :: ln_eta
    real(RP) :: temp_vfunc
    real(RP) :: temp_
    real(RP) :: geopot
    real(RP) :: cos_etav
    real(RP) :: cos_1ov2_etav
    real(RP) :: cos_3ov2_etav

    integer,  parameter :: ITRMAX = 1000
    real(RP), parameter :: CONV_EPS = 5E-15_RP    
    !------------------------------------------------

    del_eta = 1.0_RP

    !-- The loop for iteration
    itr = 0
    eta = 1.0E-8_RP ! Set first guess of eta
    do while( abs(del_eta) > CONV_EPS ) 
      etav = 0.5_RP * PI * ( eta - ETA0 ) 
      cos_etav = cos(etav)
      cos_1ov2_etav = sqrt(cos_etav)
      cos_3ov2_etav = cos_1ov2_etav**3

      temp_ = REF_TEMP * eta**( Rdry * LAPSE_RATE / Grav )
      geopot = GRAV / LAPSE_RATE * (  REF_TEMP  - temp_ )

      if ( ETAt > eta ) then
        temp_ = temp_ + DELTAT * ( ETAt - eta )**5
        geopot = geopot &
               - Rdry * DELTAT * (   ( log(eta/ETAt) + 137.0_RP / 60.0_RP ) * ETAt**5 - 5.0_RP * ETAt**4 * eta &
                                   + 5.0_RP * ETAt**3 * eta**2 - 10.0_RP / 3.0_RP * ETAt**2 * eta**3           &
                                   + 1.25_RP * ETAt * eta**4 - 0.2_RP * eta**5                                 ) 
      end if

      temp_  = temp_ &
             + 0.75_RP * eta * PI * U0 / Rdry * sin(etav) * cos_1ov2_etav &
                * ( geopot_hvari1 * 2.0_RP * cos_3ov2_etav + geopot_hvari2 )
      
      geopot = geopot  &
             + U0 * cos_3ov2_etav * ( geopot_hvari1 * cos_3ov2_etav + geopot_hvari2 )
      
      del_eta = -  ( - Grav * z + geopot )         & ! <- F
                 * ( - eta / ( Rdry * temp_ ) )      ! <- (dF/deta)^-1

      eta_save = eta
      eta = eta + del_eta
      itr = itr + 1

      if ( itr > ITRMAX ) then
          LOG_ERROR("BAROCLINIC_WAVE_setup",*) "Fail the convergence of iteration. Check!"
          LOG_ERROR_CONT(*) "* (lon,lat,z)=", lon, lat, z
          LOG_ERROR_CONT(*) "itr=", itr, "del_eta=", del_eta, "eta=", eta, "temp=", temp_
          call PRC_abort
      end if                                   
    end do  !- End of loop for iteration ----------------------------

    pres = eta_save * REF_PRES
    temp = temp_
    vel_lon = U0 * cos_3ov2_etav * sin( 2.0_RP * lat )**2

    return
  end subroutine get_thermal_wind_balance_1point_itr


end module mod_user
