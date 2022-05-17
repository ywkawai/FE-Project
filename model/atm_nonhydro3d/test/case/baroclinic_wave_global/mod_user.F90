!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!    Set the initial data for baroclinic wave in global model based on Jablonowski and Williamson (2006).
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
    EPS => CONST_EPS,        &  
    PI => CONST_PI,          &
    GRAV => CONST_GRAV,      &
    Rdry => CONST_Rdry,      &
    CPdry => CONST_CPdry,    &
    CVdry => CONST_CVdry,    &
    PRES00 => CONST_PRE00,   &
    RPlanet => CONST_RADIUS, &
    OHM     => CONST_OHM
 
  use mod_exp, only: experiment

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D  

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

  real(RP), private :: f0, beta0
  real(RP), private :: y0

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  type, private, extends(experiment) :: Exp_baroclinic_wave_global
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_baroclinicwave
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_baroclinic_wave_global), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm

    !------------------------------------------

    call exp_manager%Init('baroclinic_wave')

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

  subroutine USER_calc_tendency( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

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
  subroutine exp_SetInitCond_baroclinicwave( this,                       &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_cubedsphere_cnv, only: &
      CubedSphereCnv_CS2LonLatCoord, &
      CubedSphereCnv_LonLat2CSVec

    use mod_mkinit_util, only: &
      mkinitutil_gen_GPMat,    &
      mkinitutil_gen_Vm1Mat
    use mod_exp, only: &
      TracerLocalMeshField_ptr
    
    implicit none

    class(Exp_baroclinic_wave_global), intent(inout) :: this
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
    
    ! Parameters for inital stratification
    real(RP) :: REF_TEMP    = 288.E0_RP ! The reference temperature [K]
    real(RP) :: REF_PRES    = 1.E5_RP   ! The reference pressure [Pa]
    real(RP) :: LAPSE_RATE  = 5.E-3_RP  ! The lapse rate [K/m]

    real(RP) :: ETA_tropo   = 0.2_RP    ! The value oftropopause level
    real(RP) :: ETA0        = 0.252_RP  ! The value of η at a reference level (position of the jet)
    real(RP) :: DT_strat    = 4.8E5_RP  ! The empirical temperature difference

    ! Parameters for background zonal jet
    real(RP) :: U0 = 35.E0_RP          ! The parameter associated with zonal jet maximum amplitude  [m/s]

    ! Parameters for inital perturbation of zonal wind with a Gaussian profile
    !
    real(RP) :: Up   = 1.E0_RP         ! The maximum amplitude of zonal wind perturbation [m/s]
    real(RP) :: Lp                     ! The width of Gaussian profile
    real(RP) :: latc                   ! The center point (lat) of inital perturbation
    real(RP) :: lonc                   ! The center point (lon) of inital perturbation

    namelist /PARAM_EXP/ &
      REF_TEMP, REF_PRES, LAPSE_RATE, &
      ETA_tropo, ETA0,  U0,           &
      Up, Lp, latc, lonc

    real(RP) :: geopot_hvari1, geopot_hvari2
    real(RP), allocatable :: pres_intrp(:,:,:,:,:,:)
    real(RP), allocatable :: temp_intrp(:,:,:,:,:,:)

    integer :: ke, ke2D
    integer :: i, j, k
    integer :: p1, p2, p3, p_, p2D_
      
    integer, parameter :: IntrpPolyOrder_h = 8
    integer, parameter :: IntrpPolyOrder_v = 8
    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: x_intrp(:), y_intrp(:), z_intrp(:)
    real(RP), allocatable :: x_intrp2D(:), y_intrp2D(:)
    real(RP), allocatable :: gauss_intrp2D(:)
    real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)
    real(RP), allocatable :: IntrpMat(:,:)
    real(RP), allocatable :: IntrpVM1Mat(:,:)
    integer :: p_intrp

    real(RP), allocatable :: dens_intrp(:)
    real(RP), allocatable :: etav_intrp(:)
    real(RP), allocatable :: pres_hyd_intrp(:)
    real(RP), allocatable :: temp_hyd_intrp(:)
    real(RP), allocatable :: dens_hyd_intrp(:)
    real(RP), allocatable :: cos_3ov2_etav_intrp(:)

    real(RP), allocatable :: lon(:,:,:)
    real(RP), allocatable :: lat(:,:,:)
    real(RP) :: MOMX_met_ov_coslat(elem%Np,lcmesh%Ne)
    real(RP) :: MOMY_met          (elem%Np,lcmesh%Ne)

    real(RP) :: Gamm, RGamma
    real(RP) :: Fy(elem%Np)

    integer :: ierr
    !-----------------------------------------------------------------------------

    Lp = RPlanet / 10.0_RP
    lonc = PI / 9.0_RP
    latc = PI * 2.0_RP / 9.0_RP

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
    
    call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, elem%IsLumpedMatrix() )

    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    call mkinitutil_gen_GPMat( IntrpMat, elem_intrp, elem )

    allocate( IntrpVm1Mat(elem%Np,elem_intrp%Np) )
    call mkinitutil_gen_Vm1Mat( IntrpVM1Mat, elem_intrp, elem )

    !---

    allocate( pres_intrp(elem_intrp%Nnode_h1D,elem_intrp%Nnode_h1D,elem_intrp%Nnode_v,lcmesh%NeX,lcmesh%NeY,lcmesh%NeZ) )
    allocate( temp_intrp(elem_intrp%Nnode_h1D,elem_intrp%Nnode_h1D,elem_intrp%Nnode_v,lcmesh%NeX,lcmesh%NeY,lcmesh%NeZ) )

    allocate( x_intrp(elem_intrp%Np), y_intrp(elem_intrp%Np), z_intrp(elem_intrp%Np) )
    allocate( x_intrp2D(elem_intrp%Nnode_h1D**2), y_intrp2D(elem_intrp%Nnode_h1D**2) )
    allocate( pres_hyd_intrp(elem_intrp%Np), temp_hyd_intrp(elem_intrp%Np) )
    allocate( dens_intrp(elem_intrp%Np), dens_hyd_intrp(elem_intrp%Np) )
    allocate( etav_intrp(elem_intrp%Np), cos_3ov2_etav_intrp(elem_intrp%Np) )
    allocate( gauss_intrp2D(elem_intrp%Nnode_h1D**2) )

    allocate( lon(elem_intrp%Nnode_h1D**2,lcmesh%NeX,lcmesh%NeY) )
    allocate( lat(elem_intrp%Nnode_h1D**2,lcmesh%NeX,lcmesh%NeY) )

    ! Calculate eta(=p/p_s) level corresponding to z level of each (y,z) grid point
    ! using Newton's iteration method
    
    do j=1, lcmesh%NeY
    do i=1, lcmesh%NeX
      ke2D = i + (j-1)*lcmesh%NeX
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke2D,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke2D,:),2)
      x_intrp2D(:) = vx(1) + 0.5_RP * (elem_intrp%x1(elem_intrp%Hslice(:,1)) + 1.0_RP) * ( vx(2) - vx(1) )
      y_intrp2D(:) = vy(1) + 0.5_RP * (elem_intrp%x2(elem_intrp%Hslice(:,1)) + 1.0_RP) * ( vy(3) - vy(1) )

      call CubedSphereCnv_CS2LonLatCoord( lcmesh%panelID, x_intrp2D(:), y_intrp2D(:), elem_intrp%Nnode_h1D**2, RPlanet, &
        lon(:,i,j), lat(:,i,j) )
      where( abs( lat(:,i,j) - 0.5_RP * PI ) < EPS .or.  abs( lat(:,i,j) + 0.5_RP * PI ) < EPS )
        lat(:,i,j) =lat(:,i,j) - sign(1.0E-16_RP, lat(:,i,j))
      end where
    end do
    end do

    !$omp parallel do collapse(2) private(    &
    !$omp ke2D, p2D_, geopot_hvari1,  geopot_hvari2, &
    !$omp k, ke, vz, z_intrp, p1, p2, p3, p_         )
    do j=1, lcmesh%NeY
    do i=1, lcmesh%NeX
    do p2=1, elem_intrp%Nnode_h1D
    do p1=1, elem_intrp%Nnode_h1D

      ke2D = i + (j-1)*lcmesh%NeX
      p2D_ = p1 + (p2-1)*elem_intrp%Nnode_h1D

      call get_thermal_wind_balance_1point_geopot_hvari( &
        geopot_hvari1, geopot_hvari2,                       & ! (out)
        lat(p2D_,i,j), U0                                   ) ! (in)

      do k=1, lcmesh%NeZ
        ke = ke2D + (k-1)*lcmesh%NeX*lcmesh%NeY

        vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)    
        z_intrp(:) = vz(1) + 0.5_RP*(elem_intrp%x3(:) + 1.0_RP)*(vz(5) - vz(1))
  
        do p3=1, elem_intrp%Nnode_v
          p_ = p2D_ + (p3-1)*elem_intrp%Nnode_h1D**2
          call get_thermal_wind_balance_1point_itr( &
            pres_intrp(p1,p2,p3,i,j,k), temp_intrp(p1,p2,p3,i,j,k),       & ! (out)
            geopot_hvari1, geopot_hvari2, REF_TEMP, REF_PRES, LAPSE_RATE, & ! (in)
            ETA0, ETA_tropo, DT_strat,                                    & ! (in)
            p_, ke, lat(p2D_,i,j), z_intrp(p_), .false. )                   ! (in)

        end do ! for p3
      end do ! for k
    end do ! for p1
    end do ! for p2
    end do ! for i
    end do ! for j

    Gamm = CpDry/CvDry
    RGamma = CvDry/CpDry

   !$omp parallel do collapse(3) private( ke, ke2D, p3, p2, p1, p_,       &
   !$omp  etav_intrp, cos_3ov2_etav_intrp,                                &
   !$omp  pres_hyd_intrp, temp_hyd_intrp, dens_hyd_intrp,                 &
   !$omp  gauss_intrp2D                                                   )
    do k = 1, lcmesh%NeZ
    do j = 1, lcmesh%NeY
    do i = 1, lcmesh%NeX
      ke   = i + (j-1)*lcmesh%NeX + (k-1)*lcmesh%NeX*lcmesh%NeY
      ke2D = i + (j-1)*lcmesh%NeX

      do p3 = 1, elem_intrp%Nnode_v
      do p2 = 1, elem_intrp%Nnode_h1D
      do p1 = 1, elem_intrp%Nnode_h1D
        p_ = p1 + (p2-1)*elem_intrp%Nnode_h1D + (p3-1)*elem_intrp%Nnode_h1D**2
        pres_hyd_intrp(p_) = pres_intrp(p1,p2,p3,i,j,k)
        temp_hyd_intrp(p_) = temp_intrp(p1,p2,p3,i,j,k)
      end do
      end do  
      end do

      PRES_hyd(:,ke) = matmul(IntrpMat, pres_hyd_intrp(:))    

      dens_hyd_intrp(:) = pres_hyd_intrp(:) / ( Rdry * temp_hyd_intrp(:) )
      DENS_hyd(:,ke) = matmul(IntrpVM1Mat, dens_hyd_intrp) 

      etav_intrp(:) = 0.5_RP * PI * ( pres_hyd_intrp(:) / REF_PRES - ETA0 )
      cos_3ov2_etav_intrp(:) = cos( etav_intrp(:) ) * sqrt( cos(etav_intrp(:)) )

      gauss_intrp2D(:) = exp ( - ( RPlanet / Lp &
        * acos( sin(latc) * sin(lat(:,i,j)) + cos(latc) * cos(lat(:,i,j)) * cos(lon(:,i,j) - lonc) ) )**2 )

      MOMX_met_ov_coslat(:,ke) =  & 
        + matmul( IntrpMat,                                                                       &
          dens_hyd_intrp(:) * (                                                                   &
            U0 * sin(2.0_RP * lat(elem_intrp%IndexH2Dto3D,i,j))**2 * cos_3ov2_etav_intrp(:)       &
          + Up * gauss_intrp2D(elem_intrp%IndexH2Dto3D(:)) )                                      )
              
      MOMX_met_ov_coslat(:,ke) = MOMX_met_ov_coslat(:,ke) / cos( lcmesh%lat2D(elem%IndexH2Dto3D,ke2D) )
      MOMY_met(:,ke) = 0.0_RP
    end do    
    end do  
    end do

    call CubedSphereCnv_LonLat2CSVec( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), elem%Np * lcmesh%Ne, RPlanet, &
      MOMX_met_ov_coslat(:,:), MOMY_met(:,:),                                                   &
      MOMX(:,lcmesh%NeS:lcmesh%NeE), MOMY(:,lcmesh%NeS:lcmesh%NeE)                              )

    !------
    call elem_intrp%Final()

    return
  end subroutine exp_SetInitCond_baroclinicwave

  subroutine get_thermal_wind_balance_1point_geopot_hvari( &
    geopot_hvari1, geopot_hvari2, &
    lat, U0 )
    implicit none

    real(RP), intent(out) :: geopot_hvari1
    real(RP), intent(out) :: geopot_hvari2
    real(RP), intent(in) :: lat
    real(RP), intent(in) :: U0

    real(RP) :: cos_lat, sin_lat
    !-------------------------------------------

    cos_lat = cos(lat)
    sin_lat = sin(lat)

    ! Calc horizontal variation of geopotential height
    geopot_hvari1 = ( - 2.0_RP * sin_lat**6 * ( cos_lat**2 + 1.0_RP/3.0_RP ) &
                      + 10.0_RP / 63.0_RP ) * U0**2
    geopot_hvari2 = ( 8.0_RP / 5.0_RP * cos_lat**3 * (sin_lat**2 + 2.0_RP/3.0_RP) &
                      - PI / 4.0_RP ) * RPlanet * OHM * U0           

    return
  end subroutine get_thermal_wind_balance_1point_geopot_hvari

  subroutine get_thermal_wind_balance_1point_itr( &
    pres_yz, temp_yz,                                             &
    geopot_hvari1, geopot_hvari2, REF_TEMP, REF_PRES, LAPSE_RATE, &
    ETA0, ETA_t, DT_strat,                                        &
    p_, ke, lat, z, debug )

    implicit none

    real(RP), intent(out) :: pres_yz
    real(RP), intent(out) :: temp_yz
    real(RP), intent(in) :: geopot_hvari1, geopot_hvari2 
    real(RP), intent(in) :: REF_TEMP, REF_PRES, LAPSE_RATE
    real(RP), intent(in) :: ETA0, ETA_t, DT_strat
    integer, intent(in) :: p_, ke
    real(RP), intent(in) :: lat, z
    logical, intent(in) :: debug

    integer :: itr
    real(RP) :: eta, del_eta
    real(RP) :: ln_eta
    real(RP) :: temp_vfunc
    real(RP) :: geop_vfunc
    real(RP) :: temp_
    real(RP) :: geopot
    real(RP) :: eta_v
    real(RP) :: cos_etav
    real(RP) :: sqrt_cos_etav

    integer,  parameter :: ITRMAX = 1000
    real(RP), parameter :: CONV_EPS = 1E-14_RP    
    !------------------------------------------------

    del_eta = 1.0_RP

    !-- The loop for iteration
    itr = 0
    eta = 1.0E-7_RP ! Set first guess of eta

    do while( abs(del_eta) > CONV_EPS ) 
      ln_eta = log(eta)
      eta_v = ( eta - ETA0 ) * 0.5_RP * PI
      cos_etav = cos(eta_v)
      sqrt_cos_etav = sqrt( cos_etav )

      temp_vfunc = eta**( Rdry * LAPSE_RATE / Grav )
      geop_vfunc = GRAV / LAPSE_RATE * (1.0_RP - temp_vfunc)
      if ( eta < ETA_t ) then
        temp_vfunc = temp_vfunc &
                   + DT_strat / REF_TEMP * (ETA_t - eta)**5 
        geop_vfunc = geop_vfunc &
                   - Rdry * DT_strat / REF_TEMP * ( &
                      (log(eta/ETA_t) + 137.0_RP/60.0_RP) * ETA_t**5 - 5.0_RP * ETA_t**4 * eta &
                     + 5.0_RP * ETA_t**3 * eta**2 - 10.0_RP/3.0_RP * ETA_t**2 * eta**3         &
                     + 1.25_RP * ETA_t * eta**4 - 0.2_RP * eta**5                              )
      end if

      temp_  = REF_TEMP * temp_vfunc &
             + ( geopot_hvari1 * 2.0_RP * cos_etav * sqrt_cos_etav + geopot_hvari2 ) &
                * 0.75_RP * PI * eta / Rdry * sin(eta_v) * sqrt_cos_etav
      geopot = REF_TEMP * geop_vfunc &
              + ( geopot_hvari1 * cos_etav * sqrt_cos_etav + geopot_hvari2 ) &
              * cos_etav * sqrt_cos_etav
      del_eta = -  ( - Grav * z + geopot )         & ! <- F
                 * ( - eta / ( Rdry * temp_ ) )      ! <- (dF/deta)^-1

      eta = eta + del_eta
      itr = itr + 1

      if ( debug ) then
        write(*,*) "itr=", itr, "del_eta=", del_eta, "eta=", eta, "temp=", temp_
        write(*,*) "     -- geopot=", geopot, "old_eta=", eta - del_eta
      end if
      if ( itr > ITRMAX ) then
          LOG_ERROR("BAROCLINIC_WAVE_setup",*) "Fail the convergence of iteration. Check!"
          LOG_ERROR_CONT(*) "* (lat,Z)=", lat, z
          LOG_ERROR_CONT(*) "itr=", itr, "del_eta=", del_eta, "eta=", eta, "temp=", temp_
          call PRC_abort
      end if                                   
    end do  !- End of loop for iteration ----------------------------

    pres_yz = eta*REF_PRES
    temp_yz = temp_

    return
  end subroutine get_thermal_wind_balance_1point_itr

  subroutine exp_geostrophic_balance_correction( this,                   &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    lcmesh, elem )
    
    implicit none

    class(Exp_baroclinic_wave_global), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(inout) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DRHOT(elem%Np,lcmesh%NeA)

    integer :: ke
    real(RP) :: LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lcmesh%Ne,5)
    !---------------------------------------------------
 
    return
  end subroutine exp_geostrophic_balance_correction

end module mod_user
