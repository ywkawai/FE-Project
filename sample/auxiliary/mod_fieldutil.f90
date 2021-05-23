!-------------------------------------------------------------------------------
!> module sample / auxiliary
!!
!! @par Description
!!      A module to provide idealized & analytic profiles for 1D-3D test cases
!!
!! @par Reference
!! - Kent et al. 2014:
!!   Dynamical core model intercomparison project: Tracer transport test cases. 
!!   Quarterly Journal of the Royal Meteorological Society, 140(681), 1279-1293.
!!
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_fieldutil
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_const, only: &
    PI => CONST_PI

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !

  public :: fieldutil_get_profile1d_tracer
  public :: fieldutil_get_upwind_pos1d
  
  public :: fieldutil_get_profile2d_tracer
  public :: fieldutil_get_profile2d_flow

  public :: fieldutil_get_profile2dGlobal_tracer
  public :: fieldutil_get_profile2dGlobal_flow

  public :: fieldutil_get_profile3d_tracer
  public :: fieldutil_get_profile3d_flow

  public :: fieldutil_get_profile3dGlobal_tracer
  public :: fieldutil_get_profile3dGlobal_flow

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !


  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------
contains  

  !--- 1d -----------------------------------------------------

  function fieldutil_get_upwind_pos1d(pos, ADV_VEL, nowtime, dom_min, dom_max) result(upos)
    real(RP), intent(in) :: pos(:)
    real(RP), intent(in) :: ADV_VEL
    real(RP), intent(in) :: nowtime
    real(RP), intent(in) :: dom_min, dom_max
    real(RP) :: upos(size(pos))

    integer :: period
    !-------

    period = ADV_VEL*nowtime/(dom_max - dom_min)
    upos(:) = pos(:) - (ADV_VEL*nowtime - dble(period)*(dom_max - dom_min))
    where (upos < dom_min)
      upos = dom_max + (upos - dom_min)
    end where
  end function fieldutil_get_upwind_pos1d

  !> Get 1D data whose value is set based on specified pattern. 
  !
  ! Assume that range of domain is 0 <= x <= 1.
  ! If the profile_name is 'sin', param1 indicates the wavenumber. 
  ! If the profile_name is 'cosbell', param1 indicates the half of width. 
  ! If the profile_name is 'top-hat', param1 indicates the half of width. 
  !
  subroutine fieldutil_get_profile1d_tracer( profile,  &  
    profile_name, x, params, N ) 
    implicit none
    !--------------------------------
    integer, intent(in) :: N
    real(RP), intent(out) :: profile(N)
    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: x(N)
    real(RP), intent(in) :: params(:)

    real(RP) :: dist(N)
    !--------------------------------

    profile(:) = 0.0_RP

    select case(profile_name)
    case ('constant')
      profile(:) = params(1)    
    case ('sin')
      profile(:) = sin( params(1)*2.0_RP*PI*x(:) )
    case ('gaussian-hill')
      profile(:) = exp( - 0.5_RP * ((x(:) - params(1))**2)/params(2)**2 )
    case ('cosine-bell')
      dist(:) = abs(x(:) - params(1))
      where( dist <= params(2) )
        profile(:) = (1.0_RP + cos(PI*dist/params(2)))*0.5_RP
      end where
    case ('top-hat')
      where( abs(x - params(1)) <= params(2) )
        profile(:) = 1.0_RP
      end where
    case default
      LOG_ERROR('fieldutil_get_profile2d',*) trim(profile_name)//' is not supported. Check!'
      call PRC_abort    
    end select

    return
  end subroutine fieldutil_get_profile1d_tracer

  !-- 2d--------------------------------------------------------------

  !> Get 2D data whose value is set based on specified pattern. 
  !
  ! Assume that range of domain is 0 <= x,y <= 1.
  ! If the profile_name is 'sin', param1 and param2 are the wavenumbers in x- and y- directions, respectively. 
  ! If the profile_name is 'gaussian-hill', (param1,param2) is the coordinate of center position, param3 is the half of width. The shape is isotropic about the center of domain. 
  ! If the profile_name is 'cosine-bell', (param1,param2) is the coordinate of center position, param3 is the half of width. The shape is isotropic about the center of domain. 
  ! If the profile_name is 'top-hat', (param1,param2) is the coordinate of center position, param3 and param4 are the half of width in x- and y- directions, respectively. 
  !
  subroutine fieldutil_get_profile2d_tracer( profile, &
    profile_name, x, y, params, N )
    implicit none
    !--------------------------------
    integer, intent(in) :: N
    real(RP), intent(out) :: profile(N)
    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: x(N)
    real(RP), intent(in) :: y(N)
    real(RP), intent(in) :: params(:)

    real(RP) :: dist(N)
    !------------------------------------------------------------------------

    profile(:) = 0.0_RP

    select case(profile_name)
    case ('constant')
      profile(:) = params(1)
    case ('sin')
      profile(:) = sin( params(1)*2.0_RP*PI*x(:) ) * sin( params(2)*2.0_RP*PI*y(:) )
    case ('gaussian-hill')
      profile(:) = exp( - 0.5_RP * ((x(:) - params(1))**2 + (y(:) - params(2))**2)/params(3)**2 )
    case ('cosine-bell')
      dist(:) = sqrt( (x(:) - params(1))**2 + (y(:) - params(2))**2 )
      where( dist <= params(3) )
        profile(:) = (1.0_RP + cos(PI*dist(:)/params(3)))*0.5_RP
      end where
    case ('top-hat')
      where( abs(x-params(1)) <= params(3) .and. abs(y-params(2)) <= params(4) )
        profile(:) = 1.0_RP
      end where
    case default
      LOG_ERROR('fieldutil_get_profile2d',*) trim(profile_name)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine fieldutil_get_profile2d_tracer

  !-------------------------------

  !> Get 2D velocity data whose values are set based on specified pattern. 
  !
  ! Assume that range of domain is 0 <= x,y <= 1.
  ! If the profile_name is 'constant', param1 and param2 are the compoent of velocity in x- and y- directions, respectively. 
  ! If the profile_name is 'rigid-body-rot', (param1,param2) is the coordinate of center position, param3 is the period of rigidlid rotation. 
  ! If the profile_name is 'swirling', param3 is a parameter of period and the flow reaches the maximum deformation at n x param3/2 where n is a natural number. param4 is the current time. 
  !
  subroutine fieldutil_get_profile2d_flow( flow_x, flow_y, &
    profile_name, x, y, params, N)

    implicit none
    !--------------------------------
    integer, intent(in) :: N
    real(RP), intent(out) :: flow_x(N)
    real(RP), intent(out) :: flow_y(N)
    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: x(N)
    real(RP), intent(in) :: y(N)
    real(RP), intent(in) :: params(:)

    real(RP) :: dist(N)
    real(RP) :: fac
    !------------------------------------------------------------------------

    flow_x(:) = 0.0_RP
    flow_y(:) = 0.0_RP

    select case(profile_name)
    case ('constant')
      flow_x(:) = params(1)
      flow_y(:) = params(2)
    case ('rigid-body-rot')
      flow_x(:) = - 2.0_RP*PI/params(3)*(y(:) - params(1))
      flow_y(:) = + 2.0_RP*PI/params(3)*(x(:) - params(2))
    case ('swirling')
      fac = cos(PI*params(4)/params(3))
      flow_x(:) = + sin(PI*x(:))**2 * sin(2.0_RP*PI*y(:)) * fac
      flow_y(:) = - sin(PI*y(:))**2 * sin(2.0_RP*PI*x(:)) * fac
    case default
      LOG_ERROR('fieldutil_get_flow2d',*) trim(profile_name)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine fieldutil_get_profile2d_flow

  !-- 2d Global--------------------------------------------------------------

  !> Get 2D data whose value is set based on specified pattern. 
  !
  ! Assume that range of domain is 0 <= x,y <= 1.
  ! If the profile_name is 'sin', param1 and param2 are the wavenumbers in x- and y- directions, respectively. 
  ! If the profile_name is 'gaussian-hill', (param1,param2) is the coordinate of center position, param3 is the half of width. The shape is isotropic about the center of domain. 
  ! If the profile_name is 'cosine-bell', (param1,param2) is the coordinate of center position, param3 is the half of width. The shape is isotropic about the center of domain. 
  ! If the profile_name is 'top-hat', (param1,param2) is the coordinate of center position, param3 and param4 are the half of width in x- and y- directions, respectively. 
  !
  subroutine fieldutil_get_profile2dGlobal_tracer( profile, &
    profile_name, lon, lat, RPlanet, params, N )
    implicit none
    !--------------------------------
    integer, intent(in) :: N
    real(RP), intent(out) :: profile(N)
    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: lon(N)
    real(RP), intent(in) :: lat(N)
    real(RP), intent(in) :: RPlanet
    real(RP), intent(in) :: params(:)

    real(RP) :: dist(N)
    !------------------------------------------------------------------------

    profile(:) = 0.0_RP

    select case(profile_name)
    ! case ('constant')
    !   profile(:) = params(1)
    ! case ('sin')
    !   profile(:) = sin( params(1)*2.0_RP*PI*x(:) ) * sin( params(2)*2.0_RP*PI*y(:) )
    ! case ('gaussian-hill')
    !   profile(:) = exp( - 0.5_RP * ((x(:) - params(1))**2 + (y(:) - params(2))**2)/params(3)**2 )
    case ('cosine-bell')
      dist(:) = RPlanet * acos( &
        sin(params(2))*sin(lat(:)) + cos(params(2))*cos(lat(:))*cos(lon(:) - params(1)) ) 
      where( dist <= params(3) )
        profile(:) = (1.0_RP + cos(PI*dist(:)/params(3)))*0.5_RP
      end where
    ! case ('top-hat')
    !   where( abs(x-params(1)) <= params(3) .and. abs(y-params(2)) <= params(4) )
    !     profile(:) = 1.0_RP
    !   end where
    case default
      LOG_ERROR('fieldutil_get_profile2dGlobal',*) trim(profile_name)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine fieldutil_get_profile2dGlobal_tracer

  !-------------------------------

  subroutine fieldutil_get_profile2dGlobal_flow( flow_lon, flow_lat, &
    profile_name, lon, lat, params, N)

    implicit none
    !--------------------------------
    integer, intent(in) :: N
    real(RP), intent(out) :: flow_lon(N)
    real(RP), intent(out) :: flow_lat(N)
    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: lon(N)
    real(RP), intent(in) :: lat(N)
    real(RP), intent(in) :: params(:)

    real(RP) :: dist(N)
    real(RP) :: fac
    !------------------------------------------------------------------------

    flow_lon(:) = 0.0_RP
    flow_lat(:) = 0.0_RP

    select case(profile_name)
    ! case ('constant')
    !   flow_x(:) = params(1)
    !   flow_y(:) = params(2)
    case ('rigid-body-rot')
      flow_lon(:) =   params(1) * ( cos(params(2)) * cos(lat(:)) + sin(params(2)) * cos(lon(:)) * sin(lat(:)) )
      flow_lat(:) = - params(1) * sin(params(2)) * sin(lon(:)) 
    ! case ('swirling')
    !   fac = cos(PI*params(4)/params(3))
    !   flow_x(:) = + sin(PI*x(:))**2 * sin(2.0_RP*PI*y(:)) * fac
    !   flow_y(:) = - sin(PI*y(:))**2 * sin(2.0_RP*PI*x(:)) * fac
    case default
      LOG_ERROR('fieldutil_get_flow2dGlobal',*) trim(profile_name)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine fieldutil_get_profile2dGlobal_flow

  !-- 3d--------------------------------------------------------------

  !> Get 3D data whose value is set based on specified pattern. 
  !
  ! Assume that range of domain is 0 <= x,y,z <= 1.
  ! If the profile_name is 'sin', param1, param2 and param3 are the wavenumbers in x-, y- and z- directions, respectively. 
  ! If the profile_name is 'gaussian-hill', (param1,param2,param3) is the coordinate of center position, param4 is the half of width. The shape is isotropic about the center of domain. 
  ! If the profile_name is 'cosine-bell', (param1,param2,param3) is the coordinate of center position, param4 is the half of width. The shape is isotropic about the center of domain. 
  ! If the profile_name is 'top-hat', (param1,param2,param3) is the coordinate of center position, param4, param5 and param6 are the half of width in x-, y- and z- directions, respectively. 
  !
  subroutine fieldutil_get_profile3d_tracer( profile, &
    profile_name, x, y, z, params, N )
    implicit none
    !--------------------------------
    integer, intent(in) :: N
    real(RP), intent(out) :: profile(N)
    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: x(N)
    real(RP), intent(in) :: y(N)
    real(RP), intent(in) :: z(N)
    real(RP), intent(in) :: params(:)

    real(RP) :: dist(N)
    !------------------------------------------------------------------------

    profile(:) = 0.0_RP

    select case(profile_name)
    case ('constant')
      profile(:) = params(1)
    case ('sin')
      profile(:) = sin( params(1)*2.0_RP*PI*x(:) ) * sin( params(2)*2.0_RP*PI*y(:) ) * sin( params(3)*2.0_RP*PI*z(:) )
    case ('gaussian-hill')
      dist(:) = ((x(:) - params(1))**2 + (y(:) - params(2))**2 + (z(:) - params(3))**2)/params(4)**2
      profile(:) = exp( - 0.5_RP * dist(:) )
    case ('cosine-bell')
      dist(:) = sqrt( (x(:) - params(1))**2 + (y(:) - params(2))**2 + (z(:) - params(3))**2)
      where( dist <= params(4) )
        profile(:) = (1.0_RP + cos(PI*dist(:)/params(4)))*0.5_RP
      end where
    case ('top-hat')
      where( abs(x-params(1)) <= params(4) .and. abs(y-params(2)) <= params(5) .and. abs(y-params(3)) <= params(6))
        profile(:) = 1.0_RP
      end where
    case default
      LOG_ERROR('fieldutil_get_profile2d',*) trim(profile_name)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine fieldutil_get_profile3d_tracer

  !-------------------------------

  subroutine fieldutil_get_profile3d_flow( flow_x, flow_y, flow_z, &
    profile_name, x, y, z, params, N)

    implicit none
    !--------------------------------
    integer, intent(in) :: N
    real(RP), intent(out) :: flow_x(N)
    real(RP), intent(out) :: flow_y(N)
    real(RP), intent(out) :: flow_z(N)
    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: x(N)
    real(RP), intent(in) :: y(N)
    real(RP), intent(in) :: z(N)
    real(RP), intent(in) :: params(:)

    real(RP) :: dist(N)
    real(RP) :: fac
    !------------------------------------------------------------------------

    flow_x(:) = 0.0_RP
    flow_y(:) = 0.0_RP
    flow_z(:) = 0.0_RP

    select case(profile_name)
    case ('constant')
      flow_x(:) = params(1)
      flow_y(:) = params(2)
      flow_z(:) = params(3)
    ! case ('rigid-body-rot')
    !   flow_x(:) = - 2.0_RP*PI/params(4)*(y(:) - params(1))
    !   flow_y(:) = + 2.0_RP*PI/params(4)*(x(:) - params(2))
    !   flow_z(:) = 0.0_RP
    ! case ('swirling')
    !   fac = cos(PI*params(5)/params(4))
    !   flow_x(:) = + sin(PI*x(:))**2 * sin(2.0_RP*PI*y(:)) * fac
    !   flow_y(:) = - sin(PI*y(:))**2 * sin(2.0_RP*PI*x(:)) * fac
    !   flow_z(:) = 0.0_RP
    case default
      LOG_ERROR('fieldutil_get_flow3d',*) trim(profile_name)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine fieldutil_get_profile3d_flow

  !-- 3d Global--------------------------------------------------------------

  !> Get 3D data whose value is set based on specified pattern. 
  !
  ! Assume that range of domain is 0 <= x,y <= 1.
  ! If the profile_name is 'sin', param1 and param2 are the wavenumbers in x- and y- directions, respectively. 
  ! If the profile_name is 'gaussian-hill', (param1,param2) is the coordinate of center position, param3 is the half of width. The shape is isotropic about the center of domain. 
  ! If the profile_name is 'cosine-bell', (param1,param2) is the coordinate of center position, param3 is the half of width. The shape is isotropic about the center of domain. 
  ! If the profile_name is 'top-hat', (param1,param2) is the coordinate of center position, param3 and param4 are the half of width in x- and y- directions, respectively. 
  !
!OCL SERIAL  
  subroutine fieldutil_get_profile3dGlobal_tracer( profile, &
    profile_name, lon, lat, z, RPlanet, params, N )
    implicit none
    !--------------------------------
    integer, intent(in) :: N
    real(RP), intent(out) :: profile(N)
    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: lon(N)
    real(RP), intent(in) :: lat(N)
    real(RP), intent(in) :: z(N)
    real(RP), intent(in) :: RPlanet
    real(RP), intent(in) :: params(:)

    real(RP) :: dist(N)
    character(2) :: DCMIP_TRCNAME
    !------------------------------------------------------------------------

    profile(:) = 0.0_RP

    select case(profile_name)
    ! case ('constant')
    !   profile(:) = params(1)
    ! case ('sin')
    !   profile(:) = sin( params(1)*2.0_RP*PI*x(:) ) * sin( params(2)*2.0_RP*PI*y(:) )
    ! case ('gaussian-hill')
    !   profile(:) = exp( - 0.5_RP * ((x(:) - params(1))**2 + (y(:) - params(2))**2)/params(3)**2 )
    case ('cosine-bell')
      dist(:) = RPlanet * acos( &
        sin(params(2))*sin(lat(:)) + cos(params(2))*cos(lat(:))*cos(lon(:) - params(1)) ) 
      where( dist <= params(3) )
        profile(:) = (1.0_RP + cos(PI*dist(:)/params(3)))*0.5_RP
      end where
    ! case ('top-hat')
    !   where( abs(x-params(1)) <= params(3) .and. abs(y-params(2)) <= params(4) )
    !     profile(:) = 1.0_RP
    !   end where
    case ('DCMIP_Test1-1_q1', 'DCMIP_Test1-1_q2', 'DCMIP_Test1-1_q3', 'DCMIP_Test1-1_q4')
      DCMIP_TRCNAME = profile_name(15:16)
      call DCMIP2012_tracer( profile(:),       &
        DCMIP_TRCNAME, lon(:), lat(:), z(:), N )
    case default
      LOG_ERROR('fieldutil_get_profile3dGlobal',*) trim(profile_name)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine fieldutil_get_profile3dGlobal_tracer

  !-------------------------------

!OCL SERIAL  
  subroutine fieldutil_get_profile3dGlobal_flow( flow_lon, flow_lat, flow_w, &
    profile_name, lon, lat, z, params, N)

    implicit none
    !--------------------------------
    integer, intent(in) :: N
    real(RP), intent(out) :: flow_lon(N)
    real(RP), intent(out) :: flow_lat(N)
    real(RP), intent(out) :: flow_w(N)
    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: lon(N)
    real(RP), intent(in) :: lat(N)
    real(RP), intent(in) :: z(N)
    real(RP), intent(in) :: params(:)

    real(RP) :: dist(N)
    real(RP) :: fac
    !------------------------------------------------------------------------

    flow_lon(:) = 0.0_RP
    flow_lat(:) = 0.0_RP
    flow_w(:) = 0.0_RP

    select case(profile_name)
    case ('rigid-body-rot')
      flow_lon(:) =   params(1) * ( cos(params(2)) * cos(lat(:)) + sin(params(2)) * cos(lon(:)) * sin(lat(:)) )
      flow_lat(:) = - params(1) * sin(params(2)) * sin(lon(:))
      flow_w(:)   = 0.0_RP
    case ('DCMIP_Test1-1') ! 3D deformational flow
      call DCMIP2012_deformation_flow( flow_lon(:), flow_lat(:), flow_w(:), &
        lon(:), lat(:), z(:), params(4), N )
    case default
      LOG_ERROR('fieldutil_get_flow3dGlobal',*) trim(profile_name)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine fieldutil_get_profile3dGlobal_flow

!OCL SERIAL
  subroutine DCMIP2012_tracer( q, &
    qtrcname, lon, lat, z, Np )

    use scale_const, only: &
      RPlanet => CONST_RADIUS
    implicit none

    integer, intent(in) :: Np
    real(RP), intent(out) :: q(Np)
    character(*), intent(in) :: qtrcname
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: z(Np)
    
    integer :: i
    real(RP) :: d(Np,2)
    real(RP) :: r(Np)

    real(RP) :: Rt                           !< Horizontal half-width of tracers
    real(RP), parameter :: Zt = 1000.0_RP    !< Vertical half-width of tracers
    real(RP) :: lon_c(2)                     !< Initial longitude of first and second tracers
    real(RP), parameter :: lat_c = 0.0_RP    !< Initial latitude of tracers
    real(RP), parameter :: z_c   = 5000.0_RP !< Initial altitude of tracers
    !----------------------------------------------

    Rt = 0.5_RP * RPlanet 
    lon_c(:) = (/ 5.0_RP, 7.0_RP /) * PI / 6.0_RP

    do i=1, 2
      r(:) = RPlanet * acos(sin(lat_c)*sin(lat(:)) + cos(lat_c)*cos(lat(:))*cos(lon(:)-lon_c(i)))
      d(:,i) = min( 1.0_RP, (r(:) / Rt)**2 + ((z(:) - z_c) / Zt)**2 )
    end do

    select case(qtrcname)
    case('q1')
      q(:) = 1.0_RP + 0.5_RP * ( cos(PI * d(:,1)) + cos(PI * d(:,2)) )
    case default
      LOG_ERROR('DCMIP2012_tracer',*) trim(qtrcname)//' is not supported. Check!'
      call PRC_abort
    end select
    return
  end subroutine DCMIP2012_tracer

!OCL SERIAL
  subroutine DCMIP2012_deformation_flow( U, V, W, &
    lon, lat, z, time, Np )

    use scale_const, only: &
      P00 => CONST_PRE00,     &
      Grav => CONST_GRAV,     &
      Rdry => CONST_Rdry,     &
      RPlanet => CONST_RADIUS
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: U(Np)
    real(RP), intent(out) :: V(Np)
    real(RP), intent(out) :: W(Np)
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: z(Np)
    real(RP), intent(in) :: time

    real(RP) :: lon2(Np)
    real(RP) :: omg(Np)
    real(RP) :: pres(Np)
    real(RP) :: ud(Np), ua(Np)
    real(RP) :: vd(Np), va(Np)

    real(RP), parameter :: tau = 1036800.0_RP  !> Period of motion [sec]
    real(RP) :: OMG0                           !> Maximum of the vertical pressure velocity in units Pa/s
    real(RP), parameter :: b    = 0.2_RP       !> Normalized pressure depth of the divergent layer
    real(RP), parameter :: ptop = 254.944E2_RP
    real(RP), parameter :: T0   = 300.0_RP     !> Isothermal atmospheric temperature [K]
    !----------------------------------------------

    lon2(:) = lon(:) - 2.0_RP * PI * time / tau
    OMG0 = 23000.0_RP * PI / tau

    pres(:) = P00 * exp(- Grav * z(:) / ( Rdry * T0 ) )
    omg(:) = OMG0 * sin(lon2(:)) * cos(lat(:)) * cos(2.0_RP * PI * time / tau)             &
      * ( 1.0_RP + exp( (ptop - P00) / (b * ptop) ) - exp( (pres(:) - P00) / (b * ptop) )  &
                                                    - exp( (ptop - pres(:)) / (b * ptop) ) )

    ua(:) = 10.0_RP * RPlanet / tau * sin(lon2(:))**2 * sin(2.0_RP*lat(:)) * cos(PI * time / tau) &
          + 2.0_RP * PI * RPlanet / tau *  cos(lat(:))
    va(:) = 10.0_RP * RPlanet / tau * sin(2.0_RP * lon2(:)) * cos(lat(:)) * cos(PI * time / tau)

    ud(:) = OMG0 * RPlanet / (b * ptop) * cos(lon2(:)) * cos(lat(:))**2 * cos(2.0_RP * PI * time / tau) &
            * ( - exp( (pres(:) - P00) / (b * ptop) ) + exp( (ptop - pres(:)) / (b * ptop) )  )
    vd(:) = 0.0_RP

    U(:) = ua(:) + ud(:)
    V(:) = va(:) + vd(:)
    W(:) = - omg(:) / (Grav * pres(:) / ( Rdry * T0 ) )

    return
  end subroutine DCMIP2012_deformation_flow
end module mod_fieldutil
