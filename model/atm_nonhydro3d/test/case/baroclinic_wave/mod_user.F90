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
    real(RP) :: temp_
    real(RP) :: temp_vfunc
    real(RP) :: geopot, geopot_hvari
    real(RP), allocatable :: pres_z(:,:)
    real(RP), allocatable :: temp_z(:,:)
    real(RP), allocatable :: pres_yz(:,:,:,:)
    real(RP), allocatable :: temp_yz(:,:,:,:)

    integer :: itr
    integer,  parameter :: ITRMAX = 1000
    real(RP), parameter :: CONV_EPS = 1E-15_RP

    integer :: ke, ke2D
    integer :: i, j, k
    integer :: p1, p2, p3, p_, p2D_
    
    real(RP) :: ln_eta_e(elem%Np)
    real(RP) :: yphase_e(elem%Np)
    real(RP) :: DPRES(elem%Np)
  
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
    real(RP), allocatable :: temp_hyd_intrp(:)
    real(RP), allocatable :: dens_hyd_intrp(:)
    real(RP), allocatable :: ln_eta_intrp(:)


    real(RP) :: RGamma

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
    do p3=1, 2!elem%PolyOrder_v+1
    do p2=1, 2!elem%PolyOrder_h+1
    do p1=1, 2!elem%PolyOrder_h+1
      p_ = p1 + (p2-1)*(elem%PolyOrder_h + 1) + (p3-1)*(elem%PolyOrder_h + 1)**2
      p_intrp = p1 + (p2-1)*(elem_intrp%PolyOrder_h + 1) + (p3-1)*(elem_intrp%PolyOrder_h + 1)**2
      InvV_intrp(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do
    end do
    IntrpMat(:,:) = matmul(elem%V, InvV_intrp)

    InvV_intrpVm1(:,:) = 0.0_RP
    do p3=1, 1!elem%PolyOrder_v+1
    do p2=1, 2!elem%PolyOrder_h+1
    do p1=1, 2!elem%PolyOrder_h+1
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
    
    do j=1, lcmesh%NeY
      ke2D = 1 + (j-1)*lcmesh%NeX

      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke2D,:),2)
      y_intrp(:) = vy(1) + 0.5_RP*(elem_intrp%x2(:) + 1.0_RP)*(vy(4) - vy(1))

      do p2=1, elem_intrp%Nnode_h1D
        p2D_ = 1 + (p2-1)*elem_intrp%Nnode_h1D
        
        y_ = y_intrp(p2D_)
        yphase = 2.0_RP*PI*y_/Ly
        ! Calc horizontal variation of geopotential height
        geopot_hvari = 0.5_RP*U0*(                                                                     &
          (f0 - beta0*y0)*(y_ - 0.5_RP*Ly*(1.0_RP + sin(yphase)/PI))                                   &
          + 0.5_RP*beta0*( y_**2 - Ly*y_/PI*sin(yphase) - 0.5_RP*(Ly/PI)**2*(cos(yphase) + 1.0_RP)     &
                           - Ly**2/3.0_RP                                                          )   &
          )
        do k=1, lcmesh%NeZ
          ke = ke2D + (k-1)*lcmesh%NeX*lcmesh%NeY

          vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)    
          z_intrp(:) = vz(1) + 0.5_RP*(elem_intrp%x3(:) + 1.0_RP)*(vz(5) - vz(1))
    
          do p3=1, elem_intrp%Nnode_v
            p_ = p2D_ + (p3-1)*elem_intrp%Nnode_h1D**2
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
              del_eta = -  ( - Grav*z_intrp(p_) + geopot )     & ! <- F
                             * ( - eta/(Rdry*temp_) )            ! <- (dF/deta)^-1
   
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
            if (j==1 .and. p2==1) then
              temp_z(p3,k) = REF_TEMP - LAPSE_RATE * z_intrp(p_)
              pres_z(p3,k) = REF_PRES * (temp_z(p3,k) / REF_TEMP)**(Grav / (Rdry * LAPSE_RATE))
            end if
          end do ! for p3
        end do ! for k
      end do ! for p2
    end do ! for j
    
    RGamma = CvDry/CpDry

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
      pres_hyd_intrp(:) = pres_z(elem_intrp%IndexZ1Dto3D,k)
      temp_hyd_intrp(:) = temp_z(elem_intrp%IndexZ1Dto3D,k)

      dens_intrp(:) = pres_intrp(:) / (Rdry * temp_intrp(:))
      dens_hyd_intrp(:) = pres_hyd_intrp(:) / (Rdry * temp_hyd_intrp(:))
      DENS_hyd(:,ke) = matmul( IntrpVm1Mat, dens_hyd_intrp(:) )
      DDENS(:,ke) = matmul(IntrpVm1Mat, dens_intrp(:)) - DENS_hyd(:,ke)

      PRES_hyd(:,ke) = matmul(IntrpMat, pres_hyd_intrp(:) )
      DRHOT(:,ke) =  PRES00 / Rdry  * matmul( IntrpMat,                     &
        (pres_intrp(:)/PRES00)**RGamma - (pres_hyd_intrp(:)/PRES00)**RGamma )
      ln_eta_intrp(:) = log(pres_intrp(:)/REF_PRES)

      vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
      vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
      x_intrp(:) = vx(1) + 0.5_RP*(elem_intrp%x1(:) + 1.0_RP)*(vx(2) - vx(1))
      y_intrp(:) = vy(1) + 0.5_RP*(elem_intrp%x2(:) + 1.0_RP)*(vy(4) - vy(1))

      MOMX(:,ke) = matmul( IntrpMat, &
        dens_intrp(:)*(                                                                   &
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
end module mod_user
