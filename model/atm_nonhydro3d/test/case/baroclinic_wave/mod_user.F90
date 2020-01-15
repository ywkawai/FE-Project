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

  type, private, extends(experiment) :: Exp_baroclinic_wave
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_baroclinicwave
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
    call exp_manager%Init('density_current')
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
    
    use scale_const, only: &
      PI => CONST_PI,          &
      GRAV => CONST_GRAV,      &
      Rdry => CONST_Rdry,      &
      CPdry => CONST_CPdry,    &
      CVdry => CONST_CVdry,    &
      PRES00 => CONST_PRE00,   &
      RPlanet => CONST_RADIUS, &
      OHM     => CONST_OHM
    
    use scale_element_base, only: ElementBase3D
    use scale_element_hexahedral, only: HexahedralElement
    use scale_localmesh_3d, only: LocalMesh3D    
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

    real(RP) :: f0, beta0
    real(RP) :: y0, Ly

    real(RP) :: eta
    real(RP) :: ln_eta
    real(RP) :: del_eta
    real(RP) :: y_, yphase
    real(RP) :: temp_, temp_yz(elem%Nnode_h1D,elem%Nnode_v,lcmesh%NeY,lcmesh%NeZ)
    real(RP) :: temp_vfunc
    real(RP) :: geopot, geopot_hvari
    real(RP) :: pres_yz(elem%Nnode_h1D,elem%Nnode_v,lcmesh%NeY,lcmesh%NeZ)

    integer :: itr
    integer,  parameter :: ITRMAX = 1000
    real(RP), parameter :: CONV_EPS = 1E-15_RP

    integer :: ke, ke2D
    integer :: i, j, k
    integer :: p1, p2, p3, p_, p2D_
    
    real(RP) :: ln_eta_e(elem%Np)
    real(RP) :: yphase_e(elem%Np)
    real(RP) :: DPRES(elem%Np)
  
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

    !----
    Ly = dom_ymax - dom_ymin
    y0 = 0.5_RP*(dom_ymax + dom_ymin)

    ! Set coriolis parameters
    f0 = 2.0_RP*OHM*sin(phi0Deg*PI/180.0_RP)
    beta0 = (2.0_RP*OHM/RPlanet)*cos(phi0Deg*PI/180.0_RP)

    ! Calculate eta(=p/p_s) level corresponding to z level of each (y,z) grid point
    ! using Newton's iteration method
    
    do j=1, lcmesh%NeY
      ke2D = 1 + (j-1)*lcmesh%NeX
      do p2=1, elem%Nnode_h1D
        p2D_ = 1 + (p2-1)*elem%Nnode_h1D
        
        y_ = y(p2D_,ke2D)
        yphase = 2.0_RP*PI*y_/Ly
        ! Calc horizontal variation of geopotential height
        geopot_hvari = 0.5_RP*U0*(                                                                     &
          (f0 - beta0*y0)*(y_ - 0.5_RP*Ly*(1.0_RP + sin(yphase)/PI))                                   &
          + 0.5_RP*beta0*( y_**2 - Ly*y_/PI*sin(yphase) - 0.5_RP*(Ly/PI)**2*(cos(yphase) + 1.0_RP)     &
                           - Ly**2/3.0_RP                                                          )   &
          )
        do k=1, lcmesh%NeZ
          ke = ke2D + (k-1)*lcmesh%NeX*lcmesh%NeY
          do p3=1, elem%Nnode_v
            p_ = p2D_ + (p3-1)*elem%Nnode_h1D**2
            del_eta = 1.0_RP

            !-- The loop for iteration
            itr = 0
            eta = 1.0E-8_RP ! Set first guess of eta
            do while( abs(del_eta) > CONV_EPS ) 
              ln_eta = log(eta)
              temp_vfunc = eta**(Rdry*LAPSE_RATE/Grav)   
              temp_  = REF_TEMP*temp_vfunc &
                     + geopot_hvari/Rdry*(2.0_RP*(ln_eta/b)**2 - 1.0_RP)*exp(-(ln_eta/b)**2) 
              geopot = REF_TEMP*GRAV/LAPSE_RATE*(1.0_RP - temp_vfunc)  &
                     + geopot_hvari*ln_eta*exp(-(ln_eta/b)**2)
              del_eta = -  ( - Grav*z(p_,ke) + geopot )    & ! <- F
                             *( - eta/(Rdry*temp_) )            ! <- (dF/deta)^-1
   
              eta = eta + del_eta
              itr = itr + 1

              if ( itr > ITRMAX ) then
                  LOG_ERROR("MKINIT_barocwave",*) "Fail the convergence of iteration. Check!"
                  LOG_ERROR_CONT(*) "* (X,Y,Z)=", x(p_,ke), y(p_,ke), z(p_,ke)
                  LOG_ERROR_CONT(*) "itr=", itr, "del_eta=", del_eta, "eta=", eta, "temp=", temp_
                  call PRC_abort
              end if                                   
            end do  !- End of loop for iteration ----------------------------
            pres_yz(p2,p3,j,k) = eta*REF_PRES
            temp_yz(p2,p3,j,k) = temp_
          end do
        end do       
      end do
    end do
    
    do k = 1, lcmesh%NeZ
    do j = 1, lcmesh%NeY
    do i = 1, lcmesh%NeX
      ke = i + (j-1)*lcmesh%NeX + (k-1)*lcmesh%NeX*lcmesh%NeY
      do p3 = 1, elem%Nnode_v
      do p2 = 1, elem%Nnode_h1D
      do p1 = 1, elem%Nnode_h1D
        p_ = p1 + (p2-1)*elem%Nnode_h1D + (p3-1)*elem%Nnode_h1D**2

        PRES_hyd(p_,ke) = REF_PRES*(1.0_RP - LAPSE_RATE * z(p_,ke) / REF_TEMP)**(Grav / (Rdry * LAPSE_RATE))
        !DPRES(p_) = pres_yz(p2,p3,j,k) - PRES_hyd(p_,ke)
        DENS_hyd(p_,ke) = PRES_hyd(p_,ke) / (Rdry * (REF_TEMP - LAPSE_RATE * z(p_,ke)))

        DDENS(p_,ke) = pres_yz(p2,p3,j,k) / (Rdry * temp_yz(p2,p3,j,k)) - DENS_hyd(p_,ke)
        DRHOT(p_,ke) = PRES00 / Rdry * ( (pres_yz(p2,p3,j,k)/PRES00)**(CvDry/CpDry) - (PRES_hyd(p_,ke)/PRES00)**(CvDry/CpDry) )
        ln_eta_e(p_) = log(pres_yz(p2,p3,j,k)/REF_PRES)
      end do
      end do  
      end do

      !DENS_hyd(:,ke) = - lcmesh%Escale(:,ke,3,3)*matmul(elem%Dx3, PRES_hyd(:,ke)) / Grav
      !DDENS(:,ke) =  - lcmesh%Escale(:,ke,3,3)*matmul(elem%Dx3, DPRES(:)) / Grav

      yphase_e(:) = 2.0_RP*PI*y(:,ke)/Ly
      MOMX(:,ke) = (DENS_hyd(:,ke) + DDENS(:,ke))*( &
        - U0*sin(0.5_RP*yphase_e(:))**2*ln_eta_e(:)*exp(-(ln_eta_e(:)/b)**2)        &
        + Up*exp(- ((x(:,ke) - Xc)**2 + (y(:,ke) - Yc)**2)/Lp**2)  )
      MOMY(:,ke) = 0.0_RP
      MOMZ(:,ke) = 0.0_RP
    end do    
    end do  
    end do

    return
  end subroutine exp_SetInitCond_baroclinicwave
end module mod_user
