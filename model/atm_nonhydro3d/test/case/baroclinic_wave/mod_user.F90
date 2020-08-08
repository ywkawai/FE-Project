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

  use scale_const, only: &
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

  type, private, extends(experiment) :: Exp_baroclinic_wave
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_baroclinicwave
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_baroclinic_wave), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  subroutine USER_mkinit
    implicit none
    !------------------------------------------

    return
  end subroutine USER_mkinit

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

    !-
    call exp_manager%Init('baroclinic_wave')
    call exp_manager%SetInitCond( &
      atm%mesh, atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager )
    call exp_manager%Final()

    return
  end subroutine USER_setup

  subroutine USER_calc_tendency
    implicit none
    !------------------------------------------

    return
  end subroutine USER_calc_tendency

  subroutine USER_update
    implicit none
    !------------------------------------------

    return
  end subroutine USER_update

  !------

  subroutine exp_SetInitCond_baroclinicwave( this,                      &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
    lcmesh, elem )
    
    implicit none

    class(Exp_baroclinic_wave), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
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

    ! Parameters associated with coriolis parameter on a beta-plane
    real(RP) :: Phi0Deg     = 45.E0_RP  ! The central latitude [degree_north]

    ! Parameters for background zonal jet
    real(RP) :: U0 = 35.E0_RP          ! The parameter associated with zonal jet maximum amplitude  [m/s]
    real(RP) :: b  = 2.E0_RP           ! The vertical half-width [1]

    ! Parameters for inital perturbation of zonal wind with a Gaussian profile
    !
    real(RP) :: Up  = 1.E0_RP         ! The maximum amplitude of zonal wind perturbation [m/s]
    real(RP) :: Lp  = 600.E3_RP       ! The width of Gaussian profile
    real(RP) :: Xc  = 2000.E3_RP      ! The center point (x) of inital perturbation
    real(RP) :: Yc  = 2500.E3_RP      ! The center point (y) of inital perturbation

    namelist /PARAM_EXP/ &
      REF_TEMP, REF_PRES, LAPSE_RATE, &
      phi0Deg,                        &
      U0, b,                          &
      Up, Lp, Xc, Yc

    real(RP) :: Ly
    real(RP) :: y_
    real(RP) :: geopot_hvari
    real(RP), allocatable :: pres_z(:,:)
    real(RP), allocatable :: temp_z(:,:)
    real(RP), allocatable :: pres_yz(:,:,:,:)
    real(RP), allocatable :: temp_yz(:,:,:,:)

    integer :: ke, ke2D
    integer :: i, j, k
    integer :: p1, p2, p3, p_, p2D_
      
    integer, parameter :: IntrpPolyOrder_h = 8
    integer, parameter :: IntrpPolyOrder_v = 8
    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: x_intrp(:), y_intrp(:), z_intrp(:)
    real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)
    real(RP), allocatable :: IntrpMat(:,:), InvV_intrp(:,:)
    real(RP), allocatable :: IntrpVm1Mat(:,:), InvV_intrpVm1(:,:)
    integer :: p_intrp

    real(RP), allocatable :: pres_intrp(:)
    real(RP), allocatable :: temp_intrp(:)
    real(RP), allocatable :: dens_intrp(:)
    real(RP), allocatable :: pres_hyd_intrp(:)
    real(RP), allocatable :: rhot_hyd_intrp(:)  
    real(RP), allocatable :: temp_hyd_intrp(:)
    real(RP), allocatable :: dens_hyd_intrp(:)
    real(RP), allocatable :: ln_eta_intrp(:)

    real(RP) :: Gamm, RGamma
    real(RP) :: Fy(elem%Np)
    integer :: ierr
    !-----------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("exp_SetInitCond_baroclinicwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("exp_SetInitCond_baroclinicwave",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !----------------------------

    call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, elem%IsLumpedMatrix() )
    allocate( IntrpMat(elem%Np,elem_intrp%Np), InvV_intrp(elem%Np,elem_intrp%Np) )
    allocate( IntrpVm1Mat(elem%Np,elem_intrp%Np), InvV_intrpVm1(elem%Np,elem_intrp%Np) )

    InvV_intrp(:,:) = 0.0_RP
    do p3=1, elem%PolyOrder_v+1
    do p2=1, elem%PolyOrder_h+1
    do p1=1, elem%PolyOrder_h+1
      p_ = p1 + (p2-1)*(elem%PolyOrder_h + 1) + (p3-1)*(elem%PolyOrder_h + 1)**2
      p_intrp = p1 + (p2-1)*(elem_intrp%PolyOrder_h + 1) + (p3-1)*(elem_intrp%PolyOrder_h + 1)**2
      InvV_intrp(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do
    end do
    IntrpMat(:,:) = matmul(elem%V, InvV_intrp)

    InvV_intrpVm1(:,:) = 0.0_RP
    do p3=1, elem%PolyOrder_v
    do p2=1, elem%PolyOrder_h+1
    do p1=1, elem%PolyOrder_h+1
      p_ = p1 + (p2-1)*(elem%PolyOrder_h + 1) + (p3-1)*(elem%PolyOrder_h + 1)**2
      p_intrp = p1 + (p2-1)*(elem_intrp%PolyOrder_h + 1) + (p3-1)*(elem_intrp%PolyOrder_h + 1)**2
      InvV_intrpVm1(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do
    end do
    IntrpVm1Mat(:,:) = matmul(elem%V, InvV_intrpVm1)
    
    !----------------------------

    allocate( pres_z(elem_intrp%Nnode_v,lcmesh%NeZ) )
    allocate( temp_z(elem_intrp%Nnode_v,lcmesh%NeZ) )
    allocate( pres_yz(elem_intrp%Nnode_h1D,elem_intrp%Nnode_v,lcmesh%NeY,lcmesh%NeZ) )
    allocate( temp_yz(elem_intrp%Nnode_h1D,elem_intrp%Nnode_v,lcmesh%NeY,lcmesh%NeZ) )

    allocate( x_intrp(elem_intrp%Np), y_intrp(elem_intrp%Np), z_intrp(elem_intrp%Np) )
    allocate( pres_intrp(elem_intrp%Np), pres_hyd_intrp(elem_intrp%Np) )
    allocate( rhot_hyd_intrp(elem_intrp%Np) )
    allocate( temp_intrp(elem_intrp%Np), temp_hyd_intrp(elem_intrp%Np) )
    allocate( dens_intrp(elem_intrp%Np), dens_hyd_intrp(elem_intrp%Np) )
    allocate( ln_eta_intrp(elem_intrp%Np) )

    Ly = dom_ymax - dom_ymin
    y0 = 0.5_RP*(dom_ymax + dom_ymin)

    ! Set coriolis parameters
    f0 = 2.0_RP*OHM*sin(phi0Deg*PI/180.0_RP)
    beta0 = (2.0_RP*OHM/RPlanet)*cos(phi0Deg*PI/180.0_RP)

    ! Calculate eta(=p/p_s) level corresponding to z level of each (y,z) grid point
    ! using Newton's iteration method
    
    !$omp parallel do collapse(2) private(    &
    !$omp   ke2D, p2D_, vy, y_, geopot_hvari, &
    !$omp   k, ke, vz, z_intrp, p3, p_        )
    do j=1, lcmesh%NeY
    do p2=1, elem_intrp%Nnode_h1D

      ke2D = 1 + (j-1)*lcmesh%NeX
      p2D_ = 1 + (p2-1)*elem_intrp%Nnode_h1D

      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke2D,:),2)
      y_ = vy(1) + 0.5_RP*(elem_intrp%x2(p2D_) + 1.0_RP)*(vy(4) - vy(1))
      
      call get_thermal_wind_balance_1point_geopot_hvari( geopot_hvari, & ! (out)
        y_, y0, Ly, f0, beta0, U0 )                                      ! (in)

      do k=1, lcmesh%NeZ
        ke = ke2D + (k-1)*lcmesh%NeX*lcmesh%NeY

        vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)    
        z_intrp(:) = vz(1) + 0.5_RP*(elem_intrp%x3(:) + 1.0_RP)*(vz(5) - vz(1))
  
        do p3=1, elem_intrp%Nnode_v
          p_ = p2D_ + (p3-1)*elem_intrp%Nnode_h1D**2

          call get_thermal_wind_balance_1point_itr( &
            pres_yz(p2,p3,j,k), temp_yz(p2,p3,j,k),          & ! (out)
            geopot_hvari, b, REF_TEMP, REF_PRES, LAPSE_RATE, & ! (in)
            p_, ke, y_, z_intrp(p_) )                          ! (in)

          if (j==1 .and. p2==1) then
            temp_z(p3,k) = REF_TEMP - LAPSE_RATE * z_intrp(p_)
            pres_z(p3,k) = REF_PRES * (temp_z(p3,k) / REF_TEMP)**(Grav / (Rdry * LAPSE_RATE))
          end if

        end do ! for p3
      end do ! for k
    end do ! for p2
    end do ! for j
    
    Gamm = CpDry/CvDry
    RGamma = CvDry/CpDry

   !$omp parallel do collapse(3) private(& 
   !$omp  ke, p3, p2, p1, p_,                              &
   !$omp  pres_intrp, temp_intrp, ln_eta_intrp,            &
   !$omp  pres_hyd_intrp, temp_hyd_intrp, dens_hyd_intrp,  &
   !$omp  rhot_hyd_intrp,                                  &
   !$omp  vx, vy, x_intrp, y_intrp                         )
    do k = 1, lcmesh%NeZ
    do j = 1, lcmesh%NeY
    do i = 1, lcmesh%NeX
      ke = i + (j-1)*lcmesh%NeX + (k-1)*lcmesh%NeX*lcmesh%NeY

      do p3 = 1, elem_intrp%Nnode_v
      do p2 = 1, elem_intrp%Nnode_h1D
      do p1 = 1, elem_intrp%Nnode_h1D
        p_ = p1 + (p2-1)*elem_intrp%Nnode_h1D + (p3-1)*elem_intrp%Nnode_h1D**2
        pres_intrp(p_) = pres_yz(p2,p3,j,k)
        temp_intrp(p_) = temp_yz(p2,p3,j,k)
      end do
      end do  
      end do
      pres_hyd_intrp(:) = pres_intrp(:)
      temp_hyd_intrp(:) = temp_intrp(:)

      rhot_hyd_intrp(:) = PRES00/Rdry * (pres_hyd_intrp(:)/PRES00)**RGamma

      PRES_hyd(:,ke) = matmul(IntrpMat, pres_hyd_intrp(:))    
      !PRES_hyd(:,ke) = PRES00 * matmul( IntrpMat, (Rdry/PRES00*rhot_hyd_intrp(:))**Gamm )
      DRHOT(:,ke) =  0.0_RP
      ln_eta_intrp(:) = log(pres_hyd_intrp(:)/REF_PRES)

      dens_hyd_intrp(:) = pres_hyd_intrp(:) / (Rdry * temp_hyd_intrp(:))
      DENS_hyd(:,ke) = matmul(IntrpVM1Mat, dens_hyd_intrp) !- lcmesh%Escale(:,ke,3,3)*matmul(elem%Dx3,PRES_hyd(:,ke)) / Grav
      !DENS_hyd(:,ke) = - lcmesh%Escale(:,ke,3,3)*matmul(elem%Dx3,PRES_hyd(:,ke)) / Grav
      DDENS(:,ke) = 0.0_RP 

      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
      x_intrp(:) = vx(1) + 0.5_RP*(elem_intrp%x1(:) + 1.0_RP)*(vx(2) - vx(1))
      y_intrp(:) = vy(1) + 0.5_RP*(elem_intrp%x2(:) + 1.0_RP)*(vy(4) - vy(1))

      MOMX(:,ke) = & !- lcmesh%Escale(:,ke,2,2)*matmul(elem%Dx2,PRES_hyd(:,ke)) / (f0 + beta0*(y(:,ke) - y0)) &
        + matmul( IntrpMat, &
          dens_hyd_intrp(:)*(                                                               &
          - U0*sin(PI*y_intrp(:)/Ly)**2 * ln_eta_intrp(:) * exp(-(ln_eta_intrp(:)/b)**2)   &
          + Up*exp(- ((x_intrp - Xc)**2 + (y_intrp - Yc)**2)/Lp**2)                      ) )
   
      MOMY(:,ke) = 0.0_RP
      MOMZ(:,ke) = 0.0_RP
    end do    
    end do  
    end do

    !------
    call elem_intrp%Final()

    return
  end subroutine exp_SetInitCond_baroclinicwave

  subroutine get_thermal_wind_balance_1point_geopot_hvari( geopot_hvari, &
    y_, y0, Ly, f0, beta0, U0 )
    implicit none

    real(RP), intent(out) :: geopot_hvari
    real(RP), intent(in) :: y_
    real(RP), intent(in) :: y0
    real(RP), intent(in) :: Ly
    real(RP), intent(in) :: f0
    real(RP), intent(in) :: beta0
    real(RP), intent(in) :: U0

    real(RP) :: yphase
    !-------------------------------------------

    yphase = 2.0_RP*PI*y_/Ly

    ! Calc horizontal variation of geopotential height
    geopot_hvari = 0.5_RP*U0*(                                                                     &
      (f0 - beta0*y0)*(y_ - 0.5_RP*Ly*(1.0_RP + sin(yphase)/PI))                                   &
      + 0.5_RP*beta0*( y_**2 - Ly*y_/PI*sin(yphase) - 0.5_RP*(Ly/PI)**2*(cos(yphase) + 1.0_RP)     &
                        - Ly**2/3.0_RP                                                          )   &
    )

    return
  end subroutine get_thermal_wind_balance_1point_geopot_hvari

  subroutine get_thermal_wind_balance_1point_itr( &
    pres_yz, temp_yz,                                &
    geopot_hvari, b, REF_TEMP, REF_PRES, LAPSE_RATE, &
    p_, ke, y, z )

    implicit none

    real(RP), intent(out) :: pres_yz
    real(RP), intent(out) :: temp_yz
    real(RP), intent(in) :: geopot_hvari
    real(RP), intent(in) :: b, REF_TEMP, REF_PRES, LAPSE_RATE
    integer, intent(in) :: p_, ke
    real(RP), intent(in) :: y, z

    integer :: itr
    real(RP) :: eta, del_eta
    real(RP) :: ln_eta
    real(RP) :: temp_vfunc
    real(RP) :: temp_
    real(RP) :: geopot

    integer,  parameter :: ITRMAX = 1000
    real(RP), parameter :: CONV_EPS = 1E-15_RP    
    !------------------------------------------------

    del_eta = 1.0_RP

    !-- The loop for iteration
    itr = 0
    eta = 1.0E-8_RP ! Set first guess of eta
    do while( abs(del_eta) > CONV_EPS ) 
      ln_eta = log(eta)
      temp_vfunc = eta**(Rdry*LAPSE_RATE/Grav)   
      temp_  = REF_TEMP * temp_vfunc &
              + geopot_hvari / Rdry * (2.0_RP*(ln_eta/b)**2 - 1.0_RP)*exp(-(ln_eta/b)**2) 
      geopot = REF_TEMP * GRAV / LAPSE_RATE*(1.0_RP - temp_vfunc)  &
              + geopot_hvari * ln_eta * exp(-(ln_eta/b)**2)
      del_eta = -  ( - Grav*z + geopot )         & ! <- F
                 * ( - eta/(Rdry*temp_) )          ! <- (dF/deta)^-1

      eta = eta + del_eta
      itr = itr + 1

      if ( itr > ITRMAX ) then
          LOG_ERROR("MKINIT_barocwave",*) "Fail the convergence of iteration. Check!"
          LOG_ERROR_CONT(*) "* (Y,Z)=", y, z
          LOG_ERROR_CONT(*) "itr=", itr, "del_eta=", del_eta, "eta=", eta, "temp=", temp_
          call PRC_abort
      end if                                   
    end do  !- End of loop for iteration ----------------------------

    pres_yz = eta*REF_PRES
    temp_yz = temp_

    return
  end subroutine get_thermal_wind_balance_1point_itr

  subroutine exp_geostrophic_balance_correction( this,                              &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    lcmesh, elem )
    
    implicit none

    class(Exp_baroclinic_wave), intent(inout) :: this
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
    call cal_del_flux_dyn( del_flux, &
      PRES_hyd, lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), &
      lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3), &
      lcmesh%VMapM, lcmesh%VMapP, lcmesh, elem )

    do ke=lcmesh%NeS, lcmesh%NeE
      LiftDelFlx(:) = matmul(elem%Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,2))    
      MOMX(:,ke) = MOMX(:,ke) - LiftDelFlx(:) / (f0 + beta0*(lcmesh%pos_en(:,ke,2) - y0))


      LiftDelFlx(:) = matmul(elem%Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,3))     
      DENS_hyd(:,ke) = DENS_hyd(:,ke) - LiftDelFlx(:)/Grav
    end do

    return
  end subroutine exp_geostrophic_balance_correction

  subroutine cal_del_flux_dyn( del_flux, &
    PRES_hyd,  ny, nz, x, y, z, vmapM, vmapP, lmesh, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne,5) 
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: x(elem%Np*lmesh%Ne), y(elem%Np*lmesh%Ne), z(elem%Np*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    !------------------------------------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i) 
      del_flux(i,2) = 0.5_RP*( PRES_hyd(iP) - PRES_hyd(iM) )*ny(i)
      del_flux(i,3) = 0.5_RP*( PRES_hyd(iP) - PRES_hyd(iM) )*nz(i)
    end do

    return
  end subroutine cal_del_flux_dyn


end module mod_user
