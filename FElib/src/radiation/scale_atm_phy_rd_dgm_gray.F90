!> module FElib / Atmosphere / Physics radiation
!!
!! @par Description
!!      Radiation process with a gray-radiation scheme
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_rd_dgm_gray
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

  type, public :: AtmPhyRadGray
    integer :: optdep_type                  !< Type of optical depth calculation

    !
    real(RP) :: OPTDEP_VALLIS_EQ5_PARAMS(4) !< Parameters associated with optical depth calculation (Eq.5 in Vallis et al. (2018))
                                            !! A, mu, B, C

    real(RP) :: OPTDEP_VALLIS_EQ6a_PARAMS(3) !< Parameters associated with optical depth calculation (Eq.6a in Vallis et al. (2018))
                                             !! a_SW, b_SW, c_SW
    real(RP) :: OPTDEP_VALLIS_EQ6b_PARAMS(4) !< Parameters associated with optical depth calculation (Eq.6b in Vallis et al. (2018))
                                             !! a_LW, b_LW, c_LW, d_LW
    real(RP) :: OPTDEP_VALLIS_EQ6c_PARAMS(4) !< Parameters associated with optical depth calculation (Eq.6c in Vallis et al. (2018))
                                             !! a_win, b_win, c_win, d_win

    real(RP) :: pCO2 !< CO2 concentration in ppmv
    real(RP) :: diffFactor !< Diffusivity factor for the two-stream approximation
  contains
    procedure, public :: Init  => atm_phy_rd_dgm_gray_Init
    procedure, public :: Final => atm_phy_rd_dgm_gray_Final
    procedure, public :: calculate_rad_flux  => atm_phy_rd_dgm_gray_flux
    procedure, private :: calc_optical_thick
  end type AtmPhyRadGray
  
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

contains
  !> Initialize an object to represent a gray-radiation scheme
!OCL SERIAL
  subroutine atm_phy_rd_dgm_gray_Init( this )
    implicit none
    class(AtmPhyRadGray), intent(inout) :: this

    character(len=H_MID) :: OPTICAL_DEPTH_TYPE            = 'Vallis2018_Eq5' !< How to calculate optical depth. 'Vallis2018_Eq5'
    character(len=H_MID) :: OPTDEP_VALLIS_EQ5_PARAMS_TYPE = 'V2018'          !< How to specify the parameters for optical depth calculation (Eq.5 in Vallis et al. (2018)). 'BO2013', 'V2018', 'USER'
    real(RP) :: OPTDEP_VALLIS_EQ5_PARAMS(4) !< Array of parameters for optical depth calculation (Eq.5 in Vallis et al. (2018)). A, mu, B, C

    real(RP) :: DiffusivityFactor = 1.66_RP !< Diffusivity factor for the two-stream approximation
    real(RP) :: pCO2 = 360.0_RP             !< CO2 concentration in ppmv

    namelist / PARAM_ATMOS_PHY_RD_DGM_GRAY / &
      OPTICAL_DEPTH_TYPE,            &
      OPTDEP_VALLIS_EQ5_PARAMS_TYPE, &
      OPTDEP_VALLIS_EQ5_PARAMS,      &
      DiffusivityFactor,             &
      pCO2
        
    integer :: ierr
    !--------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_dgm_gray_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_DGM_GRAY,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_RD_dgm_gray_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_RD_dgm_gray_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_RD_DGM_GRAY. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_RD_DGM_GRAY)

    select case(trim(OPTICAL_DEPTH_TYPE))
    case('Vallis2018_Eq5')
      this%optdep_type = OPTDEP_TYPE_VALLIS2018_EQ5
    case default
      LOG_ERROR("ATMOS_PHY_RD_dgm_gray_setup",*) 'OPTICAL_DEPTH_TYPE is invalid. Check!', trim(OPTICAL_DEPTH_TYPE)
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
        LOG_ERROR("ATMOS_PHY_RD_dgm_gray_setup",*) 'OPTDEP_VALLIS_EQ5_PARAMS_TYPE is invalid. Check!', trim(OPTDEP_VALLIS_EQ5_PARAMS_TYPE)
        call PRC_abort
      end select
    end if
    return
  end subroutine atm_phy_rd_dgm_gray_Init

  !> Finalize an object to represent a gray-radiation scheme
!OCL SERIAL
  subroutine atm_phy_rd_dgm_gray_Final(this)
    implicit none
    class(AtmPhyRadGray), intent(inout) :: this
    !--------------------------------------------------------------------
    return
  end subroutine atm_phy_rd_dgm_gray_Final

  !> Calculate radiative fluxes using a gray-radiation scheme
!OCL SERIAL
  subroutine atm_phy_rd_dgm_gray_flux( this, &
    flux, flux_top, sflx_dn,                         & ! (out)
    SOLINS, PRES, TEMP, DENS, QV, SFC_TEMP, SFC_ALB, & ! (in)
    lcmesh, elem3D, lcmesh2D, elem2D )                 ! (in)
    use scale_const, only: &
      PI => CONST_PI, &
      STB => CONST_STB
    implicit none
    class(AtmPhyRadGray), intent(inout) :: this
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

    real(RP) :: dtau(elem3D%Nnode_v-1,lcmesh%NeZ,2)
    real(RP) :: CO2(elem3D%Nnode_v,lcmesh%NeZ)
    real(RP) :: trans(elem3D%Nnode_v-1,lcmesh%NeZ,2) !< Transmission across the layer
    real(RP) :: Src(elem3D%Nnode_v-1,lcmesh%NeZ)   !< Source term across the layer
    real(RP) :: temp_
    real(RP) :: D

    real(RP) :: flux_up(elem3D%Nnode_v,lcmesh%NeZ,2)
    real(RP) :: flux_dn(elem3D%Nnode_v,lcmesh%NeZ,2)
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------

    D = this%diffFactor

    !$omp parallel do private(CO2, trans, temp_, flux_up, flux_dn, Src, dtau) collapse(2)
    do ke_h=1, lcmesh%Ne2D
    do p_h=1, elem3D%Nnode_h1D**2
      CO2(:,:) = this%pCO2
      call this%calc_optical_thick( dtau,               &
        PRES(:,:,p_h,ke_h), QV(:,:,p_h,ke_h), CO2(:,:), &
        elem3D%Nnode_v, lcmesh%NeZ                      )

      trans(:,:,:) = exp(- D * dtau(:,:,:))

      !--
      do ke_z=1, lcmesh%NeZ
      do p_z=1, elem3D%Nnode_v-1
        temp_ = 0.5_RP * ( TEMP(p_z,ke_z,p_h,ke_h) + TEMP(p_z+1,ke_z,p_h,ke_h) )
        Src(p_z,ke_z) = temp_**4 * STB
      end do
      end do

      flux_dn(elem3D%Nnode_v,lcmesh%NeZ,I_LW) = 0.0_RP
      flux_dn(elem3D%Nnode_v,lcmesh%NeZ,I_SW) = SOLINS(p_h,ke_h)
      do ke_z=lcmesh%NeZ, 1, -1
        do p_z=elem3D%Nnode_v-1, 1, -1
          flux_dn(p_z,ke_z,I_LW) = flux_dn(p_z+1,ke_z,I_LW) * trans(p_z,ke_z,I_LW) + Src(p_z,ke_z) * ( 1.0_RP - trans(p_z,ke_z,I_LW) )
          flux_dn(p_z,ke_z,I_SW) = flux_dn(p_z+1,ke_z,I_SW) * trans(p_z,ke_z,I_SW)
        end do
        if ( ke_z > 1 ) then
          flux_dn(elem3D%Nnode_v,ke_z-1,I_LW) = flux_dn(1,ke_z,I_LW)
          flux_dn(elem3D%Nnode_v,ke_z-1,I_SW) = flux_dn(1,ke_z,I_SW)
        end if
      end do

      flux_up(1,1,I_LW) = SFC_TEMP(p_h,ke_h)**4 * STB
      flux_up(1,1,I_SW) = SFC_ALB(p_h,ke_h) * flux_dn(1,1,I_SW)
      do ke_z=1, lcmesh%NeZ
        do p_z=1, elem3D%Nnode_v-1
          flux_up(p_z+1,ke_z,I_LW) = flux_up(p_z,ke_z,I_LW) * trans(p_z,ke_z,I_LW) + Src(p_z,ke_z) * ( 1.0_RP - trans(p_z,ke_z,I_LW) )
          flux_up(p_z+1,ke_z,I_SW) = flux_up(p_z,ke_z,I_SW) * trans(p_z,ke_z,I_SW)
        end do
        if ( ke_z < lcmesh%NeZ ) then
          flux_up(1,ke_z+1,I_LW) = flux_up(elem3D%Nnode_v,ke_z,I_LW)
          flux_up(1,ke_z+1,I_SW) = flux_up(elem3D%Nnode_v,ke_z,I_SW)
        end if
      end do

      do ke_z=1, lcmesh%NeZ
      do p_z=1, elem3D%Nnode_v
        flux(p_h,p_z,ke_h,ke_z,I_LW,I_up) = flux_up(p_z,ke_z,I_LW)
        flux(p_h,p_z,ke_h,ke_z,I_LW,I_dn) = flux_dn(p_z,ke_z,I_LW)

        flux(p_h,p_z,ke_h,ke_z,I_SW,I_up) = flux_up(p_z,ke_z,I_SW)
        flux(p_h,p_z,ke_h,ke_z,I_SW,I_dn) = flux_dn(p_z,ke_z,I_SW)
      end do
      end do

      flux_top(p_h,ke_h,I_LW,I_up) = flux_up(elem3D%Nnode_v,lcmesh%NeZ,I_LW)
      flux_top(p_h,ke_h,I_LW,I_dn) = flux_dn(elem3D%Nnode_v,lcmesh%NeZ,I_LW)
      flux_top(p_h,ke_h,I_SW,I_up) = flux_up(elem3D%Nnode_v,lcmesh%NeZ,I_SW)
      flux_top(p_h,ke_h,I_SW,I_dn) = flux_dn(elem3D%Nnode_v,lcmesh%NeZ,I_SW)
      sflx_dn(p_h,ke_h,I_LW) = flux_dn(1,1,I_LW)
      sflx_dn(p_h,ke_h,I_SW) = flux_dn(1,1,I_SW)
    end do
    end do
    return
  end subroutine atm_phy_rd_dgm_gray_flux

!-- private -----
!OCL SERIAL
  subroutine calc_optical_thick( this, dtau, &
    pres, qv, CO2, Nnode_v, NeZ )
    implicit none
    class(AtmPhyRadGray), intent(in) :: this
    integer, intent(in) :: Nnode_v
    integer, intent(in) :: NeZ
    real(RP), intent(out) :: dtau(Nnode_v-1,NeZ,2)
    real(RP), intent(in) :: pres(Nnode_v,NeZ)
    real(RP), intent(in) :: qv(Nnode_v,NeZ)
    real(RP), intent(in) :: CO2(Nnode_v,NeZ)

    integer :: ke_z, p_z
    real(RP) :: A, B, C
    real(RP) :: dsig

    real(RP), parameter :: P0 = 1.0e5_RP
    real(RP) :: qv_tmp
    !---------------------------------------------------------

    A = this%OPTDEP_VALLIS_EQ5_PARAMS(1) * this%OPTDEP_VALLIS_EQ5_PARAMS(2)
    B = this%OPTDEP_VALLIS_EQ5_PARAMS(3)
    C = this%OPTDEP_VALLIS_EQ5_PARAMS(4)

    do ke_z=1, NeZ
    do p_z=1, Nnode_v-1
      qv_tmp = 0.5_RP * ( qv(p_z,ke_z) + qv(p_z+1,ke_z) )
      dsig = max( pres(p_z,ke_z) - pres(p_z+1,ke_z), 0.0_RP ) / P0
      dtau(p_z,ke_z,I_LW) = ( A + B * qv_tmp + C * log( CO2(p_z,ke_z) / 360.0_RP ) ) * dsig
      dtau(p_z,ke_z,I_SW) = 0.0_RP * dsig
    end do
    end do
    return
  end subroutine calc_optical_thick

end  module scale_atm_phy_rd_dgm_gray
