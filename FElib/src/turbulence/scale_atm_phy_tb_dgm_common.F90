!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Physics turbulence / Common
!!
!! @par Description
!!      A coomon modules for atmospheric turbulent parameterization
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_tb_dgm_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D

  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_phy_tb_dgm_common_get_varinfo
  public :: atm_phy_tb_dgm_common_calc_lambda
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T13_ID     = 1
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T23_ID     = 2
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T33_ID     = 3
  integer, public, parameter :: ATMOS_PHY_TB_AUX_DPTDZ_ID   = 4
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T31_ID     = 5
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T32_ID     = 6
  integer, public, parameter :: ATMOS_PHY_TB_AUX_DPTDX_ID   = 7
  integer, public, parameter :: ATMOS_PHY_TB_AUX_DPTDY_ID   = 8
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T11_ID     = 9
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T12_ID     = 10
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T21_ID     = 11
  integer, public, parameter :: ATMOS_PHY_TB_AUX_T22_ID     = 12
  integer, public, parameter :: ATMOS_PHY_TB_AUX_SCALAR_NUM  = 4
  integer, public, parameter :: ATMOS_PHY_TB_AUX_HVEC_NUM    = 2
  integer, public, parameter :: ATMOS_PHY_TB_AUX_HTENSOR_NUM = 1
  integer, public, parameter :: ATMOS_PHY_TB_AUX_NUM        = 12
  
  integer, public, parameter :: ATMOS_PHY_TB_DIAG_TKE_ID     = 1
  integer, public, parameter :: ATMOS_PHY_TB_DIAG_NU_ID      = 2
  integer, public, parameter :: ATMOS_PHY_TB_DIAG_KH_ID      = 3
  integer, public, parameter :: ATMOS_PHY_TB_DIAG_NUM        = 3

  integer, public, parameter :: ATMOS_PHY_TB_MOMX_t_ID  = 1
  integer, public, parameter :: ATMOS_PHY_TB_MOMY_t_ID  = 2
  integer, public, parameter :: ATMOS_PHY_TB_MOMZ_t_ID  = 3
  integer, public, parameter :: ATMOS_PHY_TB_RHOT_t_ID  = 4
  integer, public, parameter :: ATMOS_PHY_TB_TENDS_NUM1 = 4
  !-

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------
  
  private :: fact, p1, p2, p3
  
contains
!OCL SERIAL
  subroutine atm_phy_tb_dgm_common_get_varinfo( &
    auxvar_info, diagvar_info, tend_info        )

    implicit none

    type(VariableInfo), intent(out) :: auxvar_info(ATMOS_PHY_TB_AUX_NUM)
    type(VariableInfo), intent(out) :: diagvar_info(ATMOS_PHY_TB_DIAG_NUM)
    type(VariableInfo), intent(out) :: tend_info(ATMOS_PHY_TB_TENDS_NUM1)

    type(VariableInfo) :: ATMOS_PHY_TB_AUX_VINFO(ATMOS_PHY_TB_AUX_NUM)
    DATA ATMOS_PHY_TB_AUX_VINFO / &
      VariableInfo( ATMOS_PHY_TB_AUX_T13_ID, 'TB_T13', 'stress tensor (T13)',          & ! rho x nu x S : [kg/m3] x [m2/s] x [s-1]
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                 ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T23_ID, 'TB_T23', 'stress tensor (T23)',          &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                 ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T33_ID, 'TB_T33', 'stress tensor (T33)',          &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                 ), &
      VariableInfo( ATMOS_PHY_TB_AUX_DPTDZ_ID, 'DPTDZ', 'gradient of PT (z)',          &
                    'K/m',  3, 'XYZ',  ''                                           ), &                                      
      VariableInfo( ATMOS_PHY_TB_AUX_T31_ID, 'TB_T31', 'stress tensor (T31)',          &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                 ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T32_ID, 'TB_T32', 'stress tensor (T32)',          &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                 ), &
      VariableInfo( ATMOS_PHY_TB_AUX_DPTDX_ID, 'DPTDX', 'gradient of PT (x)',          &
                    'K/m',  3, 'XYZ',  ''                                           ), &
      VariableInfo( ATMOS_PHY_TB_AUX_DPTDY_ID, 'DPTDY', 'gradient of PT (y)',          &
                    'K/m',  3, 'XYZ',  ''                                           ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T11_ID, 'TB_T11', 'stress tensor (T11)',          &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                 ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T12_ID, 'TB_T12', 'stress tensor (T12)',          &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                 ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T21_ID, 'TB_T21', 'stress tensor (T21)',          &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                 ), &
      VariableInfo( ATMOS_PHY_TB_AUX_T22_ID, 'TB_T22', 'stress tensor (T22)',          &
                    'kg.m-3.m2.s-2',  3, 'XYZ',  ''                                 )  /

    type(VariableInfo) :: ATMOS_PHY_TB_DIAG_VINFO(ATMOS_PHY_TB_DIAG_NUM)
    DATA ATMOS_PHY_TB_DIAG_VINFO / &
      VariableInfo( ATMOS_PHY_TB_DIAG_TKE_ID, 'TKE', 'SGS turbluence kinetic energy',   &
                    'm2/s2',  3, 'XYZ',  ''                                          ), &
      VariableInfo( ATMOS_PHY_TB_DIAG_NU_ID, 'NU', 'eddy viscosity',                    &
                    'm2/s',  3, 'XYZ',  ''                                           ), &
      VariableInfo( ATMOS_PHY_TB_DIAG_KH_ID, 'KH', 'eddy diffusion',                    &
                    'm2/s',  3, 'XYZ',  ''                                           )  /
                    
    type(VariableInfo) :: ATMOS_PHY_TB_TEND_VINFO(ATMOS_PHY_TB_TENDS_NUM1)
    DATA ATMOS_PHY_TB_TEND_VINFO / &
      VariableInfo( ATMOS_PHY_TB_MOMX_t_ID, 'TB_MOMX_t', 'tendency of x-momentum in TB process',    &
                    'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
      VariableInfo( ATMOS_PHY_TB_MOMY_t_ID, 'TB_MOMY_t', 'tendency of y-momentum in TB process',    &
                    'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
      VariableInfo( ATMOS_PHY_TB_MOMZ_t_ID, 'TB_MOMZ_t', 'tendency of z-momentum in TB process',    &
                    'kg/m2/s2',  3, 'XYZ',  ''                                                   ), &
      VariableInfo( ATMOS_PHY_TB_RHOT_t_ID, 'TB_RHOT_t', 'tendency of rho*PT in TB process',        &
                    'kg/m3.K/s', 3, 'XYZ',  ''                                                   )  / 

    !----------------------------------------------------------

    auxvar_info(:) = ATMOS_PHY_TB_AUX_VINFO
    diagvar_info(:) = ATMOS_PHY_TB_DIAG_VINFO
    tend_info(:) = ATMOS_PHY_TB_TEND_VINFO

    return
  end subroutine atm_phy_tb_dgm_common_get_varinfo

!OCL SERIAL  
  subroutine atm_phy_tb_dgm_common_calc_lambda( lambda, & ! (out)
    Cs, filter_fac, lmesh, elem, lmesh2D, elem2D        ) ! (in)

    use scale_const, only: &
      KARMAN  => CONST_KARMAN
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: lambda(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: Cs
    real(RP), intent(in) :: filter_fac

    integer :: ke

    real(RP) :: vol
    real(RP) :: he
    real(RP) :: lambda0
    real(RP) :: Zs(elem2D%Np)
    real(RP) :: Z1(elem2D%Np)
    real(RP) :: dz(elem  %Np)
    real(RP) :: FZ
    real(RP) :: elem_aspect_eff

    real(RP), parameter :: OneOverThree  = 1.0_RP / 3.0_RP
    !--------------------------------------------------------------------

    !$omp parallel do private( &
    !$omp vol, lambda0, Zs, Z1, dz, he, FZ, elem_aspect_eff )
    do ke=lmesh%NeS, lmesh%NeE
      vol = sum( elem%IntWeight_lgl(:) * lmesh%J(:,ke) * lmesh%Gsqrt(:,ke) )
      he = ( vol / ( dble(elem%PolyOrder_h+1)**2 * dble(elem%PolyOrder_v+1) ) )**OneOverThree

      FZ = lmesh%zlev(elem%Colmask(elem%Nnode_v,1),ke) - lmesh%zlev(elem%Colmask(1,1),ke)      
      elem_aspect_eff = fact( FZ/dble(elem%PolyOrder_v+1), sqrt(vol/FZ)/dble(elem%PolyOrder_h+1), sqrt(vol/FZ)/dble(elem%PolyOrder_h+1) )
      lambda0 = elem_aspect_eff * Cs * filter_fac * he

      Zs(:) = lmesh%zlev(elem%Hslice(:,1),lmesh%EMap3Dto2D(ke))
      Z1(:) = lmesh%zlev(elem%Hslice(:,2),lmesh%EMap3Dto2D(ke))
      dz(:) = max( lmesh%zlev(:,ke) - Zs(elem%IndexH2Dto3D(:)), he )

      !lambda(:,ke) = sqrt( 1.0_RP / (1.0_RP / lambda0**2 + 1.0_RP / ( KARMAN * max( dz(:), EPS ) )**2 ) )
      lambda(:,ke) = sqrt( 1.0_RP / (1.0_RP / lambda0**2 + 1.0_RP / ( KARMAN * ( dz(:) + 1.0E-4_RP ) )**2 ) )    
    end do

    return
  end subroutine atm_phy_tb_dgm_common_calc_lambda

!-- private ---------------------------------------------

!OCL SERIAL
  elemental function fact(dz, dx, dy)
    implicit none
    real(RP), intent(in) :: dz
    real(RP), intent(in) :: dx
    real(RP), intent(in) :: dy
    real(RP) :: fact ! (out)

    real(RP), parameter :: oot = -1.0_RP/3.0_RP
    real(RP), parameter :: fot =  5.0_RP/3.0_RP
    real(RP), parameter :: eot = 11.0_RP/3.0_RP
    real(RP), parameter :: tof = -3.0_RP/4.0_RP
    real(RP) :: a1, a2, b1, b2, dmax
    !--------------------------------------------------------------------

    dmax = max(dz, dx, dy)
    if ( dz == dmax ) then
       a1 = dx / dmax
       a2 = dy / dmax
    else if ( dx == dmax ) then
       a1 = dz / dmax
       a2 = dy / dmax
    else ! dy == dmax
       a1 = dz / dmax
       a2 = dx / dmax
    end if
    b1 = atan( a1/a2 )
    b2 = atan( a2/a1 )

   fact = 1.736_RP * (a1*a2)**oot &
         * ( 4.0_RP*p1(b1)*a1**oot + 0.222_RP*p2(b1)*a1**fot + 0.077*p3(b1)*a1**eot - 3.0_RP*b1 &
           + 4.0_RP*p1(b2)*a2**oot + 0.222_RP*p2(b2)*a2**fot + 0.077*p3(b2)*a2**eot - 3.0_RP*b2 &
           )**tof
   return
  end function fact
!OCL SERIAL
  elemental function p1(z)
    implicit none
    real(RP), intent(in) :: z
    real(RP) :: p1 ! (out)

    real(RP), parameter :: TwoOverThree  = 2.0_RP / 3.0_RP
    !--------------------------------------------------------------------

    p1 = 2.5_RP * p2(z) - 1.5_RP * sin(z) * cos(z)**TwoOverThree
    return
  end function p1
!OCL SERIAL
  elemental function p2(z)
    implicit none
    real(RP), intent(in) :: z
    real(RP) :: p2 ! (out)
    !--------------------------------------------------------------------

    p2 = 0.986_RP * z + 0.073_RP * z**2 - 0.418_RP * z**3 + 0.120_RP * z**4
    return
  end function p2
!OCL SERIAL
  elemental function p3(z)
    implicit none
    real(RP), intent(in) :: z
    real(RP) :: p3 ! (out)
    !--------------------------------------------------------------------

    p3 = 0.976_RP * z + 0.188_RP * z**2 - 1.169_RP * z**3 + 0.755_RP * z**4 - 0.151_RP * z**5
    return
  end function p3

end module scale_atm_phy_tb_dgm_common
