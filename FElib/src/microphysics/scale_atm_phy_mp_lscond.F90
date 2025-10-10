!> module FElib / Atmosphere / Physics cloud microphysics / Large-scale condensation
!!
!! @par Description
!!     A module to provide a large-scale condensation scheme
!!
!! @par Reference
!!
!! @author Yuta Kawai, Team SCALE
!!
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_mp_lscond
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort 
  use scale_const, only: &
    Rdry => CONST_Rdry,    &
    Rvap => CONST_Rvap,    &
    CPdry => CONST_CPdry,  & 
    CVdry => CONST_CVdry,  & 
    LHV0 => CONST_LHV0,    &
    PRES00 => CONST_PRE00, &
    Grav => CONST_GRAV
  
  use scale_element_base, only: &
    ElementBase1D, ElementBase3D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  public :: ATMOS_PHY_MP_lscond_setup
  public :: ATMOS_PHY_MP_lscond_adjustment
  public :: ATMOS_PHY_MP_lscond_precipitation
  public :: ATMOS_PHY_MP_lscond_precipitation_momentum

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !
  integer, private, parameter :: QA_MP  = 2

  integer,                parameter, public :: ATMOS_PHY_MP_lscond_ntracers = QA_MP
  integer,                parameter, public :: ATMOS_PHY_MP_lscond_nwaters = 1
  integer,                parameter, public :: ATMOS_PHY_MP_lscond_nices = 0
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_MP_lscond_tracer_names(QA_MP) = (/ &
       'QV', &
       'QC'  /)
  character(len=H_MID)  , parameter, public :: ATMOS_PHY_MP_lscond_tracer_descriptions(QA_MP) = (/ &
       'Ratio of Water Vapor mass to total mass (Specific humidity)', &
       'Ratio of Cloud Water mass to total mass                    '  /)
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_MP_lscond_tracer_units(QA_MP) = (/ &
       'kg/kg',  &
       'kg/kg'   /)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !


  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter   :: I_QV = 1
  integer,  private, parameter   :: I_QC = 2

  integer,  private, parameter   :: I_hyd_QC =  1

  logical,  private              :: flag_liquid = .true.     ! warm rain
  logical,  private              :: couple_aerosol = .false. ! Consider CCN effect ?

  ! real(RP), private, parameter   :: re_qc =  8.E-6_RP        ! effective radius for cloud water

contains

  !-----------------------------------------------------------------------------
  !> Setup a module for a large-scale condensation scheme
  !!
  subroutine ATMOS_PHY_MP_lscond_setup
    use scale_prc, only: &
       PRC_abort
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_lscond_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_MP_lscond_setup",*) 'large-scale condensation'

    if( couple_aerosol ) then
       LOG_ERROR("ATMOS_PHY_MP_lscond_setup",*) 'MP_aerosol_couple should be .false. for large-scale condensation type MP!'
       call PRC_abort
    endif

    return
  end subroutine ATMOS_PHY_MP_lscond_setup

!> Calculate a state after the saturation process
!!
!OCL SERIAL
  subroutine ATMOS_PHY_MP_lscond_adjustment( &
    KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! (in)
    DENS, PRES,                         & ! (in)
    TEMP, QTRC, CPtot, CVtot,           & ! (inout)
    RHOE_t, EVAPORATE                   ) ! (out)

    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_pres2qsat_liq
    use scale_atmos_hydrometeor, only: &
      CV_WATER      
    use scale_const, only: &
      CL => CONST_CL
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: PRES(KA,IA,JA)
    real(RP), intent(inout) :: TEMP(KA,IA,JA) 
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP)
    real(RP), intent(inout) :: CPtot(KA,IA,JA) 
    real(RP), intent(inout) :: CVtot(KA,IA,JA) 
    real(RP), intent(out) :: RHOE_t(KA,IA,JA)
    real(RP), intent(out) :: EVAPORATE(KA,IA,JA)

    real(RP) :: Qsat(KA,IA,JA)
    real(RP) :: tmp(KA)
    real(RP) :: coef1
    real(RP) :: coef2

    integer :: i, j
    !----------------------------------------------

    call ATMOS_SATURATION_pres2qsat_liq( &
      KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! (in)
      TEMP, PRES,                         & ! (in)
      Qsat                                ) ! (out)


    coef1 = LHV0**2 / ( CVDry * Rvap )
    coef2 = Rvap / LHV0

    !$omp parallel do private(i,j, tmp) collapse(2)
    do j=JS, JE
    do i=IS, IE
      tmp(:) = max( ( QTRC(:,i,j,I_QV) - Qsat(:,i,j) ) / ( 1.0_RP + coef1 / TEMP(:,i,j)**2 * ( 1.0_RP - coef2 * TEMP(:,i,j) ) * Qsat(:,i,j) ), &
                    0.0_RP )

      QTRC(:,i,j,I_QV) = QTRC(:,i,j,I_QV) - tmp(:)
      QTRC(:,i,j,I_QC) = tmp(:)

      RHOE_t(:,i,j) = LHV0 * DENS(:,i,j) * tmp(:)
    end do
    end do
    return
  end subroutine ATMOS_PHY_MP_lscond_adjustment

!OCL SERIAL
  subroutine ATMOS_PHY_MP_lscond_precipitation( &
    DENS, RHOQ, CPtot, CVtot, RHOE,            & ! (inout)
    sflx_rain, sflx_snow, esflx,               & ! (inout)
    TEMP,                                      & ! (in)
    QHA, QLA, QIA, lcmesh, elem, elem1D        ) ! (in)
    
    use scale_const, only: &
      UNDEF => CONST_UNDEF8
    use scale_atmos_hydrometeor, only: &
      CV_WATER, &
      CP_WATER, &
      CV_ICE,   &
      CP_ICE    
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: QHA                   !< hydrometeor (water + ice)
    real(RP), intent(inout) :: DENS (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: RHOQ (elem%Np,lcmesh%NeZ,lcmesh%Ne2D,QHA)
    real(RP), intent(inout) :: CPtot(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: CVtot(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: RHOE (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: sflx_rain(elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(inout) :: sflx_snow(elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(inout) :: esflx    (elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(in) :: TEMP (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    integer, intent(in) :: QLA, QIA
    class(ElementBase1D), intent(in) :: elem1D

    real(RP) :: RHOCP(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOCV(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)

    real(RP) :: dDENS(elem%Np)
    real(RP) :: dInternalEn(elem%Np)

    real(RP) :: vint_weight(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: condens_vint_lc(elem%Nnode_h1D**2)
    real(RP) :: ien_vint_lc(elem%Nnode_h1D**2)

    real(RP) :: eflx(elem%Np)
    real(RP) :: CP(QHA)
    real(RP) :: CV(QHA)

    integer :: ke2D
    integer :: p2D

    integer :: ke_z
    integer :: ke
    integer :: iq
    !-------------------------------------------------------

    do iq = 1, QHA
      if ( iq > QLA + QIA ) then
        CP(iq) = UNDEF 
        CV(iq) = UNDEF
      else if ( iq > QLA ) then ! ice water
        CP(iq) = CP_ICE
        CV(iq) = CV_ICE
      else                      ! liquid water
        CP(iq) = CP_WATER
        CV(iq) = CV_WATER
      end if
    end do

    !$omp parallel do collapse(2)
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ
      RHOCP(:,ke_z,ke2D) = CPtot(:,ke_z,ke2D) * DENS(:,ke_z,ke2D)
      RHOCV(:,ke_z,ke2D) = CVtot(:,ke_z,ke2D) * DENS(:,ke_z,ke2D)
    end do
    end do

    do iq = 1, QHA
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D

        dDENS(:) = - RHOQ(:,ke_z,ke2D,iq)
        RHOQ(:,ke_z,ke2D,iq) = 0.0_RP

        do p2D=1, elem%Nnode_h1D**2
          vint_weight(:,p2D) = 0.5_RP * elem1D%IntWeight_lgl(:) * ( lcmesh%zlev(elem%Colmask(elem%Nnode_v,p2D),ke) - lcmesh%zlev(elem%Colmask(1,p2D),ke) )
        end do

        do p2D=1, elem%Nnode_h1D**2
          condens_vint_lc(p2D) = sum( vint_weight(:,p2D) * dDENS(elem%Colmask(:,p2D)) )
        end do
      
        if ( ke_z == 1 ) then
          if ( iq > QLA ) then ! ice water
              sflx_snow(:,ke2D) = sflx_snow(:,ke2D)  &
                                + condens_vint_lc(:)
          else                 ! liquid water
              sflx_rain(:,ke2D) = sflx_rain(:,ke2D)  &
                                + condens_vint_lc(:)
          end if
        end if

        !--- update density

        RHOCP(:,ke_z,ke2D) = RHOCP(:,ke_z,ke2D) + CP(iq) * dDENS(:)
        RHOCV(:,ke_z,ke2D) = RHOCV(:,ke_z,ke2D) + CV(iq) * dDENS(:)
        DENS (:,ke_z,ke2D) = DENS(:,ke_z,ke2D) + dDENS(:)

        !--- update internal energy   

        dInternalEn(:) = CP(iq) * dDENS(:) * TEMP(:,ke_z,ke2D)

        do p2D=1, elem%Nnode_h1D**2
          ien_vint_lc(p2D) = sum( vint_weight(:,p2D) * dInternalEn(elem%Colmask(:,p2D))  )
        end do

        RHOE(:,ke_z,ke2D) = RHOE(:,ke_z,ke2D) + dInternalEn(:)

        if ( ke_z == 1 ) then
          esflx(:,ke2D) = esflx(:,ke2D) &
                        + ien_vint_lc(:)
        end if
      end do
      end do
    end do

    !$omp parallel do collapse(2)
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ  
      CPtot(:,ke_z,ke2D) = RHOCP(:,ke_z,ke2D) / DENS(:,ke_z,ke2D)
      CVtot(:,ke_z,ke2D) = RHOCV(:,ke_z,ke2D) / DENS(:,ke_z,ke2D)
    end do
    end do

    return
  end subroutine ATMOS_PHY_MP_lscond_precipitation

!OCL SERIAL
  subroutine ATMOS_PHY_MP_lscond_precipitation_momentum( &
    MOMU_t, MOMV_t, MOMZ_t,                & ! (out)
    DENS, MOMU, MOMV, MOMZ, DENS_new,      & ! (in)
    rdt_MP, lcmesh, elem                   ) ! (in)
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: MOMU_t(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMV_t(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMZ_t(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMU(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMV(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMZ(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: DENS_new(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: rdt_MP
    
    integer :: ke2D
    integer :: ke_z
    integer :: ke
    real(RP) :: coef(elem%Np)
    !----------------------------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke2D, ke_z, ke, coef )
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      coef(:) = ( DENS_new(:,ke_z,ke2D) / DENS(:,ke_z,ke2D) - 1.0_RP ) * rdt_MP

      MOMU_t(:,ke) = coef(:) * MOMU(:,ke_z,ke2D)
      MOMV_t(:,ke) = coef(:) * MOMV(:,ke_z,ke2D)
      MOMZ_t(:,ke) = coef(:) * MOMZ(:,ke_z,ke2D)
    end do
    end do
    return
  end subroutine ATMOS_PHY_MP_lscond_precipitation_momentum
  
end module scale_atm_phy_mp_lscond
