!> module FElib / Atmosphere / Physics radiation
!!
!! @par Description
!!      A module to provide simplified radiation schemes based on Vallis et al. (2018)
!!
!! @author Yuta Kawai, Team SCALE
!!
!! @par Reference
!!  - Vallis, G.K., Colyer, G., Geen, R., Gerber, E., Jucker, M., Maher, P.,
!!    Paterson, A., Pietschnig, M., Penn, J., and Thomson, S.I. 2018:
!!    Isca, v1.0: a framework for the global modelling of the atmospheres of Earth and other planets at varying levels of complexity.
!!    Geosci. Model Dev., 11, 843–859.
!!  - Geen, R., Czaja, A., and Haigh, J.D. 2016:
!!    The effects of increasing humidity on heat transport by extratropical waves.
!!    Geophys. Res. Lett., 43, 8314–8321.
!!
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_rd_dgm_simple
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_atmos_phy_rd_common, only: &
    I_up, I_dn, I_LW, I_SW
  
  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  !-----------------------------------------------------------------------------
  implicit none
  private  
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !


  integer, parameter :: N_LW_BND_MAX = 2

  type, public :: AtmPhyRadSimple
    integer :: optdep_type                  !< Type of optical depth calculation

    !- Parameters for optical depth calculation with a gray-radiation scheme
    real(RP) :: OPTDEP_VALLIS_EQ5_PARAMS(4) !< Parameters associated with optical depth calculation (Eq.5 in Vallis et al. (2018))
                                            !! A, mu, B, C

    !- Parameters for optical depth calculation with a two-band longwave and one-band shortwave radiation scheme
    real(RP) :: OPTDEP_VALLIS_EQ6a_PARAMS(2) !< Parameters associated with optical depth calculation (Eq.6a in Vallis et al. (2018))
                                             !! a_SW, c_SW
    real(RP) :: OPTDEP_VALLIS_EQ7a_PARAMS(4) !< Parameters associated with optical depth calculation (Eq.7a in Vallis et al. (2018))
                                             !! a_LW, b_LW, c_LW, d_LW
    real(RP) :: OPTDEP_VALLIS_EQ7b_PARAMS(4) !< Parameters associated with optical depth calculation (Eq.7b in Vallis et al. (2018))
                                             !! a_win, b_win, c_win, d_win

    real(RP) :: FRAC_LW(N_LW_BND_MAX)    !< The fraction of the longwave spectrum at each band

    real(RP) :: pCO2 !< CO2 concentration in ppmv
    real(RP) :: diffFactor !< Diffusivity factor for the two-stream approximation
  contains
    procedure, public :: Init  => atm_phy_rd_dgm_simple_Init
    procedure, public :: Final => atm_phy_rd_dgm_simple_Final
    procedure, public :: calculate_rad_flux  => atm_phy_rd_dgm_simple_flux
    procedure, private :: calc_optical_thick
  end type AtmPhyRadSimple
  
  !-----------------------------------------------------------------------------
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
  integer, parameter :: OPTDEP_TYPE_VALLIS2018_EQ5 = 1 !< An idealized radiation scheme based on Eq.5
                                                       !! This is gray in infrared so that a single optical thickness is defined for the entire longwave spectrum, 
                                                       !! which includes a parameterization of long-wave absorption by CO2.
  integer, parameter :: OPTDEP_TYPE_VALLIS2018_EQ6 = 2 !< An idealized radiation scheme based on Eq.6
                                                       !! This scheme has two infrared bands and one solar band, as described in Green et al. (2016),
                                                       !! which provides an intermediate complexity between gray radiation and more sophisticated radiation schemes.

contains
  !> Initialize an object to represent a gray-radiation scheme
!OCL SERIAL
  subroutine atm_phy_rd_dgm_simple_Init( this )
    implicit none
    class(AtmPhyRadSimple), intent(inout) :: this

    character(len=H_MID) :: SIMPLE_RD_TYPE                = 'Gray'     !< Type of radiation scheme. 'Gray', 'LWbnd2+SWbnd1'

    character(len=H_MID) :: OPTDEP_VALLIS_EQ5_PARAMS_TYPE = 'V2018'    !< How to specify the parameters for optical depth calculation (Eq.5 in Vallis et al. (2018)). 'BO2013', 'V2018', 'USER'
    real(RP) :: OPTDEP_VALLIS_EQ5_PARAMS(4)                            !< Array of parameters for optical depth calculation (Eq.5 in Vallis et al. (2018)). A, mu, B, C

    character(len=H_MID) :: OPTDEP_VALLIS_EQ6_PARAMS_TYPE = 'V2018'    !< How to specify the parameters for optical depth calculation (Eq.6 in Vallis et al. (2018)). 'V2018', 'USER'
    real(RP) :: OPTDEP_VALLIS_EQ6a_PARAMS(2)                           !< Array of parameters for optical depth calculation (Eq.6 in Vallis et al. (2018)). a_SW, c_SW
    real(RP) :: OPTDEP_VALLIS_EQ7a_PARAMS(4)                           !< Array of parameters for optical depth calculation (Eq.6 in Vallis et al. (2018)). a_LW, b_LW, c_LW, d_LW
    real(RP) :: OPTDEP_VALLIS_EQ7b_PARAMS(4)                           !< Array of parameters for optical depth calculation (Eq.6 in Vallis et al. (2018)). a_win, b_win, c_win, d_win
    real(RP) :: WINDOW_FRAC = 0.2_RP                                   !< Window fraction for LWbnd2+SWbnd1 scheme. This is the fraction of the longwave spectrum that is transparent to radiation, which is set to 0.2 by default.

    real(RP) :: DiffusivityFactor = 1.66_RP !< Diffusivity factor for the two-stream approximation
    real(RP) :: pCO2 = 360.0_RP             !< CO2 concentration in ppmv

    namelist / PARAM_ATMOS_PHY_RD_DGM_SIMPLE / &
      SIMPLE_RD_TYPE,                &
      !-
      OPTDEP_VALLIS_EQ5_PARAMS_TYPE, &
      OPTDEP_VALLIS_EQ5_PARAMS,      &
      !-
      OPTDEP_VALLIS_EQ6_PARAMS_TYPE, &
      OPTDEP_VALLIS_EQ6a_PARAMS,     &
      OPTDEP_VALLIS_EQ7a_PARAMS,     &
      OPTDEP_VALLIS_EQ7b_PARAMS,     &
      WINDOW_FRAC,                   &
      !-
      DiffusivityFactor,             &
      pCO2
        
    integer :: ierr
    !--------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_dgm_simple_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_DGM_SIMPLE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_RD_dgm_simple_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_RD_dgm_simple_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_RD_DGM_SIMPLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_RD_DGM_SIMPLE)

    select case(trim(SIMPLE_RD_TYPE))
    case('Gray')
      this%optdep_type = OPTDEP_TYPE_VALLIS2018_EQ5
      this%FRAC_LW = [ 1.0_RP, 0.0_RP ]
    case('LWbnd2+SWbnd1')
      this%optdep_type = OPTDEP_TYPE_VALLIS2018_EQ6
    case default
      LOG_ERROR("ATMOS_PHY_RD_dgm_simple_setup",*) 'SIMPLE_RD_TYPE is invalid. Check!', trim(SIMPLE_RD_TYPE)
      call PRC_abort
    end select

    this%pCO2 = pCO2
    this%diffFactor = DiffusivityFactor

    if ( this%optdep_type == OPTDEP_TYPE_VALLIS2018_EQ5 ) then
      select case(trim(OPTDEP_VALLIS_EQ5_PARAMS_TYPE))
      case('BO2013')
        this%OPTDEP_VALLIS_EQ5_PARAMS = [ 0.8678_RP, 1.0_RP, 1997.9_RP, 0.0_RP ]
      case('V2018')
        this%OPTDEP_VALLIS_EQ5_PARAMS = [ 0.1627_RP, 1.0_RP, 1997.9_RP, 0.17_RP ]
      case('USER')
        this%OPTDEP_VALLIS_EQ5_PARAMS(:) = OPTDEP_VALLIS_EQ5_PARAMS(:)
      case default
        LOG_ERROR("ATMOS_PHY_RD_dgm_simple_setup",*) 'OPTDEP_VALLIS_EQ5_PARAMS_TYPE is invalid. Check!', trim(OPTDEP_VALLIS_EQ5_PARAMS_TYPE)
        call PRC_abort
      end select
    end if

    if ( this%optdep_type == OPTDEP_TYPE_VALLIS2018_EQ6 ) then
      select case(trim(OPTDEP_VALLIS_EQ6_PARAMS_TYPE))
      case('V2018') 
        ! Based on Vallis et al. (2018), the default values for these coefficients were fitted to output from Santa Barbara DISORT Atmospheric Radiative Transfer 60 (SBDART).1
        this%OPTDEP_VALLIS_EQ6a_PARAMS = [ 5.96E-2_RP, 2.9E-3_RP ]
        this%OPTDEP_VALLIS_EQ7a_PARAMS = [ 0.1_RP, 23.8_RP, 254.0_RP, 0.0954_RP ]
        this%OPTDEP_VALLIS_EQ7b_PARAMS = [ 0.215_RP, 1.4711E2_RP, 1.0814E4_RP, 0.2023_RP ]
        this%FRAC_LW = [ 0.6268_RP, 0.3732_RP ]
      case('USER')
        this%OPTDEP_VALLIS_EQ6a_PARAMS(:) = OPTDEP_VALLIS_EQ6a_PARAMS(:)
        this%OPTDEP_VALLIS_EQ7a_PARAMS(:) = OPTDEP_VALLIS_EQ7a_PARAMS(:)
        this%OPTDEP_VALLIS_EQ7b_PARAMS(:) = OPTDEP_VALLIS_EQ7b_PARAMS(:)
        this%FRAC_LW = [ 1.0_RP - WINDOW_FRAC, WINDOW_FRAC ]
      case default
        LOG_ERROR("ATMOS_PHY_RD_dgm_simple_setup",*) 'OPTDEP_VALLIS_EQ6_PARAMS_TYPE is invalid. Check!', trim(OPTDEP_VALLIS_EQ6_PARAMS_TYPE)
        call PRC_abort
      end select
    end if
    return
  end subroutine atm_phy_rd_dgm_simple_Init

  !> Finalize an object to represent a gray-radiation scheme
!OCL SERIAL
  subroutine atm_phy_rd_dgm_simple_Final(this)
    implicit none
    class(AtmPhyRadSimple), intent(inout) :: this
    !--------------------------------------------------------------------
    return
  end subroutine atm_phy_rd_dgm_simple_Final

  !> Calculate radiative fluxes using a gray-radiation scheme
!OCL SERIAL
  subroutine atm_phy_rd_dgm_simple_flux( this, &
    flux, flux_top, sflx_dn,                         & ! (out)
    SOLINS, PRES, TEMP, DENS, QV, SFC_TEMP, SFC_ALB, & ! (in)
    lcmesh, elem3D, lcmesh2D, elem2D )                 ! (in)
    use scale_const, only: &
      PI => CONST_PI, &
      STB => CONST_STB
    implicit none
    class(AtmPhyRadSimple), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: flux(elem3D%Nnode_h1D**2,elem3D%Nnode_v,lcmesh%Ne2D,lcmesh%NeZ,2,2)
    real(RP), intent(out) :: flux_top(elem2D%Np,lcmesh2D%Ne,2,2)
    real(RP), intent(out) :: sflx_dn(elem2D%Np,lcmesh2D%Ne,2)
    real(RP), intent(in) :: PRES(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP), intent(in) :: TEMP(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP), intent(in) :: DENS(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP), intent(in) :: QV(elem3D%Nnode_v,lcmesh%NeZ,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP), intent(in) :: SOLINS(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: SFC_ALB(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: SFC_TEMP(elem2D%Np,lcmesh2D%NeA)

    integer :: ke, ke_z, ke_h
    integer :: p, p_z, p_h

    real(RP) :: dtau_lw(elem3D%Nnode_v-1,lcmesh%NeZ,N_LW_BND_MAX)
    real(RP) :: dtau_sw(elem3D%Nnode_v-1,lcmesh%NeZ)
    real(RP) :: trans_lw(elem3D%Nnode_v-1,lcmesh%NeZ,N_LW_BND_MAX) !< Transmission across the layer
    real(RP) :: trans_sw(elem3D%Nnode_v-1,lcmesh%NeZ)              !< Transmission across the layer
    real(RP) :: CO2(elem3D%Nnode_v,lcmesh%NeZ)
    real(RP) :: Src(elem3D%Nnode_v-1,lcmesh%NeZ)                   !< Source term across the layer
    real(RP) :: temp_
    real(RP) :: D

    real(RP) :: flux_up_lw(elem3D%Nnode_v,lcmesh%NeZ,N_LW_BND_MAX)
    real(RP) :: flux_dn_lw(elem3D%Nnode_v,lcmesh%NeZ,N_LW_BND_MAX)
    real(RP) :: flux_up_lw_tot(elem3D%Nnode_v,lcmesh%NeZ)
    real(RP) :: flux_dn_lw_tot(elem3D%Nnode_v,lcmesh%NeZ)
    real(RP) :: flux_up_sw(elem3D%Nnode_v,lcmesh%NeZ)
    real(RP) :: flux_dn_sw(elem3D%Nnode_v,lcmesh%NeZ)

    integer :: n_lw_bnd
    integer :: bnd_i

    real(RP) :: R_LW(N_LW_BND_MAX) !< Fraction of the longwave spectrum in each band
    !--------------------------------------------------------------------

    D = this%diffFactor
    R_LW(:) = this%FRAC_LW(:)

    select case(this%optdep_type)
    case(OPTDEP_TYPE_VALLIS2018_EQ5)
      n_lw_bnd = 1
    case(OPTDEP_TYPE_VALLIS2018_EQ6)
      n_lw_bnd = 2
    end select

    !$omp parallel do collapse(2) &
    !$omp private( dtau_lw, dtau_sw, trans_lw, trans_sw, flux_up_lw, flux_dn_lw, flux_up_sw, flux_dn_sw, &
    !$omp          flux_up_lw_tot, flux_dn_lw_tot, CO2, temp_, Src ) 
    do ke_h=1, lcmesh%Ne2D
    do p_h=1, elem3D%Nnode_h1D**2
      CO2(:,:) = this%pCO2
      call this%calc_optical_thick( dtau_lw, dtau_sw,  &
        PRES(:,:,p_h,ke_h), QV(:,:,p_h,ke_h), CO2(:,:), &
        elem3D%Nnode_v, lcmesh%NeZ                      )

      trans_lw(:,:,1:n_lw_bnd) = exp(- D * dtau_lw(:,:,1:n_lw_bnd))
      trans_sw(:,:) = exp(- D * dtau_sw(:,:))

      !- Calculate the source term for longwave radiation
      do ke_z=1, lcmesh%NeZ
      do p_z=1, elem3D%Nnode_v-1
        temp_ = 0.5_RP * ( TEMP(p_z,ke_z,p_h,ke_h) + TEMP(p_z+1,ke_z,p_h,ke_h) )
        Src(p_z,ke_z) = temp_**4 * STB
      end do
      end do

      !- Calculate downward radiative fluxes for longwave and shortwave radiation

      flux_dn_lw_tot(:,:) = 0.0_RP
      do bnd_i =1, n_lw_bnd
        flux_dn_lw(elem3D%Nnode_v,lcmesh%NeZ,bnd_i) = 0.0_RP
        do ke_z=lcmesh%NeZ, 1, -1
          do p_z=elem3D%Nnode_v-1, 1, -1
            flux_dn_lw(p_z,ke_z,bnd_i) = flux_dn_lw(p_z+1,ke_z,bnd_i) * trans_lw(p_z,ke_z,bnd_i) &
                                       + R_LW(bnd_i) * Src(p_z,ke_z) * ( 1.0_RP - trans_lw(p_z,ke_z,bnd_i) )
          end do
          if ( ke_z > 1 ) then
            flux_dn_lw(elem3D%Nnode_v,ke_z-1,bnd_i) = flux_dn_lw(1,ke_z,bnd_i)
          end if      
        end do

        flux_dn_lw_tot(:,:) = flux_dn_lw_tot(:,:) + flux_dn_lw(:,:,bnd_i)
      end do

      flux_dn_sw(elem3D%Nnode_v,lcmesh%NeZ) = SOLINS(p_h,ke_h)
      do ke_z=lcmesh%NeZ, 1, -1
      do p_z=elem3D%Nnode_v-1, 1, -1
        flux_dn_sw(p_z,ke_z) = flux_dn_sw(p_z+1,ke_z) * trans_sw(p_z,ke_z)
      end do
        if ( ke_z > 1 ) then
          flux_dn_sw(elem3D%Nnode_v,ke_z-1) = flux_dn_sw(1,ke_z)
        end if      
      end do
  
      !- Calculate upward radiative fluxes for longwave and shortwave radiation

      flux_up_lw_tot(:,:) = 0.0_RP
      do bnd_i =1, n_lw_bnd
        flux_up_lw(1,1,bnd_i) = R_LW(bnd_i) * SFC_TEMP(p_h,ke_h)**4 * STB
        do ke_z=1, lcmesh%NeZ
          do p_z=1, elem3D%Nnode_v-1
            flux_up_lw(p_z+1,ke_z,bnd_i) = flux_up_lw(p_z,ke_z,bnd_i) * trans_lw(p_z,ke_z,bnd_i) &
                                         + R_LW(bnd_i) * Src(p_z,ke_z) * ( 1.0_RP - trans_lw(p_z,ke_z,bnd_i) )
          end do
          if ( ke_z < lcmesh%NeZ ) then
            flux_up_lw(1,ke_z+1,bnd_i) = flux_up_lw(elem3D%Nnode_v,ke_z,bnd_i)
          end if
        end do

        flux_up_lw_tot(:,:) = flux_up_lw_tot(:,:) + flux_up_lw(:,:,bnd_i)
      end do

      flux_up_sw(1,1) = SFC_ALB(p_h,ke_h) * flux_dn_sw(1,1)
      do ke_z=1, lcmesh%NeZ
        do p_z=1, elem3D%Nnode_v-1
          flux_up_sw(p_z+1,ke_z) = flux_up_sw(p_z,ke_z) * trans_sw(p_z,ke_z)
        end do
        if ( ke_z < lcmesh%NeZ ) then
          flux_up_sw(1,ke_z+1) = flux_up_sw(elem3D%Nnode_v,ke_z)
        end if
      end do

      !- Store the calculated fluxes into the output arrays

      do ke_z=1, lcmesh%NeZ
      do p_z=1, elem3D%Nnode_v
        flux(p_h,p_z,ke_h,ke_z,I_LW,I_up) = flux_up_lw_tot(p_z,ke_z)
        flux(p_h,p_z,ke_h,ke_z,I_LW,I_dn) = flux_dn_lw_tot(p_z,ke_z)

        flux(p_h,p_z,ke_h,ke_z,I_SW,I_up) = flux_up_sw(p_z,ke_z)
        flux(p_h,p_z,ke_h,ke_z,I_SW,I_dn) = flux_dn_sw(p_z,ke_z)
      end do
      end do

      flux_top(p_h,ke_h,I_LW,I_up) = flux_up_lw_tot(elem3D%Nnode_v,lcmesh%NeZ)
      flux_top(p_h,ke_h,I_LW,I_dn) = flux_dn_lw_tot(elem3D%Nnode_v,lcmesh%NeZ)
      flux_top(p_h,ke_h,I_SW,I_up) = flux_up_sw(elem3D%Nnode_v,lcmesh%NeZ)
      flux_top(p_h,ke_h,I_SW,I_dn) = flux_dn_sw(elem3D%Nnode_v,lcmesh%NeZ)
      sflx_dn(p_h,ke_h,I_LW) = flux_dn_lw_tot(1,1)
      sflx_dn(p_h,ke_h,I_SW) = flux_dn_sw(1,1)

    end do
    end do
    return
  end subroutine atm_phy_rd_dgm_simple_flux

!-- private -----
!OCL SERIAL
  subroutine calc_optical_thick( this, dtau_lw, dtau_sw, &
    pres, qv, CO2, Nnode_v, NeZ )
    implicit none
    class(AtmPhyRadSimple), intent(in) :: this
    integer, intent(in) :: Nnode_v
    integer, intent(in) :: NeZ
    real(RP), intent(out) :: dtau_lw(Nnode_v-1,NeZ,N_LW_BND_MAX)
    real(RP), intent(out) :: dtau_sw(Nnode_v-1,NeZ)
    real(RP), intent(in) :: pres(Nnode_v,NeZ)
    real(RP), intent(in) :: qv(Nnode_v,NeZ)
    real(RP), intent(in) :: CO2(Nnode_v,NeZ)

    integer :: ke_z, p_z
    real(RP) :: A, B, C
    real(RP) :: dsig

    real(RP), parameter :: P0 = 1.0e5_RP
    real(RP) :: qv_tmp
    !---------------------------------------------------------

    if (this%optdep_type == OPTDEP_TYPE_VALLIS2018_EQ5) then
      A = this%OPTDEP_VALLIS_EQ5_PARAMS(1) * this%OPTDEP_VALLIS_EQ5_PARAMS(2)
      B = this%OPTDEP_VALLIS_EQ5_PARAMS(3)
      C = this%OPTDEP_VALLIS_EQ5_PARAMS(4)

      do ke_z=1, NeZ
      do p_z=1, Nnode_v-1
        qv_tmp = 0.5_RP * ( qv(p_z,ke_z) + qv(p_z+1,ke_z) )
        dsig = max( pres(p_z,ke_z) - pres(p_z+1,ke_z), 0.0_RP ) / P0
        dtau_lw(p_z,ke_z,1) = ( A + B * qv_tmp + C * log( CO2(p_z,ke_z) / 360.0_RP ) ) * dsig
        dtau_sw(p_z,ke_z) = 0.0_RP * dsig
      end do
      end do
    else if (this%optdep_type == OPTDEP_TYPE_VALLIS2018_EQ6) then
    end if

    return
  end subroutine calc_optical_thick

end  module scale_atm_phy_rd_dgm_simple
