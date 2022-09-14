!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of rising warm bubble with moist process
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

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: MeshField3D
  
  
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

  type, private, extends(experiment) :: Exp_rising_therm_bubble_moist
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_rising_therm_bubble
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_rising_therm_bubble_moist), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    call exp_manager%Init('rising_therm_bubble')

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
  subroutine exp_SetInitCond_rising_therm_bubble( this,                    &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_tracer, only: &
      TRACER_inq_id
    use scale_const, only: &
      PI => CONST_PI,        &
      GRAV => CONST_GRAV,    &
      Rdry => CONST_Rdry,    &
      Rvap => CONST_Rvap,    &      
      CPdry => CONST_CPdry,  &
      CVdry => CONST_CVdry,  &
      PRES00 => CONST_PRE00, &
      Pstd   => CONST_Pstd,  &
      EPSvap => CONST_EPSvap
    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_psat_all,      &
      ATMOS_SATURATION_pres2qsat_all
    use scale_atmos_hydrometeor, only: &
      CV_VAPOR, &
      CV_WATER, &
      CP_VAPOR, &
      CP_WATER
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constPT, &
      hydrostaic_build_rho_XYZ
    use mod_mkinit_util, only: &
      mkinitutil_calc_cosinebell
    use mod_exp, only: &
      TracerLocalMeshField_ptr
    
    implicit none

    class(Exp_rising_therm_bubble_moist), intent(inout) :: this
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
    
    real(RP) :: ENV_RH       =  80.0_RP ! Relative Humidity of environment [%]
    real(RP) :: ENV_L1_ZTOP  =  1.E3_RP ! top height of the layer1 (constant THETA)       [m]
    real(RP) :: ENV_L2_ZTOP  = 14.E3_RP ! top height of the layer2 (small THETA gradient) [m]
    real(RP) :: ENV_L2_TLAPS = 4.E-3_RP ! Lapse rate of THETA in the layer2 (small THETA gradient) [K/m]
    real(RP) :: ENV_L3_TLAPS = 3.E-2_RP ! Lapse rate of THETA in the layer3 (large THETA gradient) [K/m]
    real(RP) :: SFC_THETA
    real(RP) :: SFC_PRES

    real(RP) :: DTHETA
    real(RP) :: x_c, y_c, z_c
    real(RP) :: r_x, r_y, r_z

    namelist /PARAM_EXP/ &
      ENV_RH,            &
      ENV_L1_ZTOP,       &
      ENV_L2_ZTOP,       &
      ENV_L2_TLAPS,      &
      ENV_L3_TLAPS,      &
      SFC_THETA,         &
      x_c, y_c, z_c,     &
      r_x, r_y, r_z,     &
      DTHETA

    integer, parameter :: IntrpPolyOrder_h = 7
    integer, parameter :: IntrpPolyOrder_v = 7
    real(RP), allocatable :: THETA_purtub(:,:)
    
    real(RP) :: PT  (elem%Np)
    real(RP) :: DENS (elem%Np)
    real(RP) :: DENS2(elem%Np)
    real(RP) :: PRES(elem%Np)
    real(RP) :: QV(elem%Np)

    real(RP) :: PT_tmp(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: Rtot  (elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: CPtot(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: CPtot_ov_CVtot(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: temp_z(elem%Nnode_v,lcmesh%NeZ)
    real(RP) :: pres_z(elem%Nnode_v,lcmesh%NeZ)
    real(RP) :: qdry_z(elem%Nnode_v,lcmesh%NeZ)
    real(RP) :: qv_z  (elem%Nnode_v,lcmesh%NeZ)
    real(RP) :: qsat_z(elem%Nnode_v,lcmesh%NeZ)
    real(RP) :: psat_sfc, qsat_sfc, qv_sfc

    integer :: ke, ke2D
    integer :: ke_x, ke_y, ke_z
    integer :: p, p3, p2D
    integer :: ierr

    integer :: iq_QV, iq_QC

    real(RP), allocatable :: bnd_SFC_PRES(:,:)
    real(RP) :: sfc_rhot(elem%Nfp_v)
    !-----------------------------------------------------------------------------

    SFC_THETA = 300.0_RP
    SFC_PRES  = Pstd

    x_c = 500.0_RP; y_c = 500.0_RP; z_c = 350.0_RP
    r_x = 250.0_RP; r_y = 250.0_RP; r_z = 250.0_RP;
    DTHETA    = 0.5_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("RISING_THERMAL_BUBBL_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("RISING_THERMAL_BUBBL_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    allocate( THETA_purtub(elem%Np,lcmesh%NeA) )
    call mkinitutil_calc_cosinebell( &
      THETA_purtub,                          &
      DTHETA, r_x, r_y, r_z, x_c, y_c, z_c,  &
      x, y, z, lcmesh, elem,                 &
      IntrpPolyOrder_h, IntrpPolyOrder_v     )  
    
    call hydrostatic_calc_basicstate_constPT( DENS_hyd, PRES_hyd,                            &
      SFC_THETA, SFC_PRES, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3), &
      lcmesh, elem )

    
    allocate( bnd_SFC_PRES(elem%Nfp_v,lcmesh%Ne2DA) )

    !$omp parallel do collapse(3) private( ke, ke2D, ke_x, ke_y, ke_z, sfc_rhot )
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
    do ke_z=1, lcmesh%NeZ
      ke2D = ke_x + (ke_y-1)*lcmesh%NeX
      ke = ke2D + (ke_z-1)*lcmesh%NeX*lcmesh%NeY

      where ( lcmesh%zlev(:,ke) <= ENV_L1_ZTOP )
        PT_tmp(:,ke_z,ke_x,ke_y) = SFC_THETA
      elsewhere ( lcmesh%zlev(:,ke) < ENV_L2_ZTOP  )
        PT_tmp(:,ke_z,ke_x,ke_y) = SFC_THETA + ( lcmesh%zlev(:,ke) - ENV_L1_ZTOP ) * ENV_L2_TLAPS
      elsewhere
        PT_tmp(:,ke_z,ke_x,ke_y) = SFC_THETA + ( ENV_L2_ZTOP - ENV_L1_ZTOP ) * ENV_L2_TLAPS       &
                                 + ( lcmesh%zlev(:,ke) - ENV_L2_ZTOP ) * ENV_L3_TLAPS
      end where

    end do
    end do
    end do

    call hydrostaic_build_rho_XYZ( DDENS, & ! (out)
      DENS_hyd, PRES_hyd, PT_tmp,         & ! (in)
      x, y, z, lcmesh, elem               ) ! (in)
      
    !$omp parallel do collapse(3) private(PT, DENS, PRES, ke,ke_x,ke_y,ke_z)
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
    do ke_z=1, lcmesh%NeZ
      ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1) * lcmesh%Ne2D

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      PRES_hyd(:,ke) = PRES00 * ( Rdry / PRES00 * DENS(:) * PT_tmp(:,ke_z,ke_x,ke_y) )**(CpDry/CvDry)
      DENS_hyd(:,ke) = DENS(:)
      DDENS(:,ke) = 0.0_RP
    end do
    end do
    end do

    !---
    call TRACER_inq_id( "QV", iq_QV )
    call TRACER_inq_id( "QC", iq_QC )

    ! calc QV from RH

    ! call ATMOS_SATURATION_psat_all( tempxxxxx, & ! [IN]
    !                                 psat_sfc   ) ! [OUT]
    ! qsat_sfc = EPSvap * psat_sfc / ( SFC_PRES - ( 1.0_RP - EPSvap ) * psat_sfc )
    ! qv_sfc = ENV_RH * 1.0E-2_RP * qsat_sfc
    do ke_z = 1, lcmesh%NeZ
    do p3=1, elem%Nnode_v
      ke = 1 + (ke_z - 1)*lcmesh%Ne2D
      p  = 1 + (p3 - 1)*elem%Nnode_h1D**2
      qdry_z(p3,ke_z) = 1.0_RP &
                      - tracer_field_list(iq_QV)%ptr%val(p,ke) &
                      - tracer_field_list(iq_QC)%ptr%val(p,ke)
      pres_z(p3,ke_z) = PRES_hyd(p,ke)
      temp_z(p3,ke_z) = PRES_hyd(p,ke) / ( Rdry * DENS_hyd(p,ke) )
    end do
    end do
    call ATMOS_SATURATION_pres2qsat_all( &
      elem%Nnode_v, 1, elem%Nnode_v, lcmesh%NeZ, 1, lcmesh%NeZ, &
      temp_z, pres_z, qdry_z,                                   & ! [IN]
      qsat_z                                                    ) ! [OUT]

    !$omp parallel do collapse(3) private(ke_z,ke_x,ke_y,ke,ke2D,p3,p2D,p, QV, sfc_rhot)
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
    do ke_z=1, lcmesh%NeZ
      ke2D = ke_x + (ke_y-1)*lcmesh%NeX
      ke = ke2D + (ke_z-1)*lcmesh%NeX*lcmesh%NeY

      do p3=1, elem%Nnode_v
      do p2D=1, elem%Nnode_h1D**2
        p = p2D + (p3 - 1)*elem%Nnode_h1D**2
        if ( lcmesh%zlev(p,ke) > ENV_L2_ZTOP ) then
          QV(p) = 0.0_RP
        else
          QV(p) = ENV_RH * 1.0E-2_RP * qsat_z(p3,ke_z)
        end if
      end do
      end do
      tracer_field_list(iq_QV)%ptr%val(:,ke) = QV(:)

      Rtot          (:,ke_z,ke_x,ke_y) = Rdry  * ( 1.0_RP - QV(:) ) + Rvap     * QV(:)
      CPtot         (:,ke_z,ke_x,ke_y) = CPdry * ( 1.0_RP - QV(:) ) + CP_VAPOR * QV(:)
      CPtot_ov_CVtot(:,ke_z,ke_x,ke_y) = CPtot(:,ke_z,ke_x,ke_y)                         &
                                       / ( CVdry * ( 1.0_RP - QV(:) ) + CV_VAPOR * QV(:) ) 
      if ( ke_z == 1 ) then
        sfc_rhot(:) = DENS_hyd(elem%Hslice(:,1),ke2D) * PT_tmp(elem%Hslice(:,1),ke_z,ke_x,ke_y)        
        bnd_SFC_PRES(:,ke2D) = PRES00 * ( Rtot(:,ke_z,ke_x,ke_y) * sfc_rhot(:) / PRES00 )**( CPtot_ov_CVtot(:,ke_z,ke_x,ke_y) )
      end if                              
    end do
    end do
    end do

    call hydrostaic_build_rho_XYZ( DDENS, & ! (out)
      DENS_hyd, PRES_hyd, PT_tmp,         & ! (in)
      Rtot, CPtot_ov_CVtot,               & ! (in)
      x, y, z, lcmesh, elem, bnd_SFC_PRES ) ! (in)


    !$omp parallel do collapse(3) private(PT, DENS, DENS2, PRES, ke,ke_x,ke_y,ke_z)
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
    do ke_z=1, lcmesh%NeZ
      ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1) * lcmesh%Ne2D

      ! 
      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      !PRES(:) = PRES00 * ( Rtot(:,ke_z,ke_x,ke_y) / PRES00 * DENS(:) * PT_tmp(:,ke_z,ke_x,ke_y) )**CPtot_ov_CVtot(:,ke_z,ke_x,ke_y)

      !
      PT(:) = PT_tmp(:,ke_z,ke_x,ke_y) + THETA_purtub(:,ke)
!      DENS2(:) = PRES(:) / ( Rtot(:,ke_z,ke_x,ke_y) * PT(:) * (PRES(:)/PRES00)**( Rtot(:,ke_z,ke_x,ke_y) / CPtot(:,ke_z,ke_x,ke_y) ) )
      DDENS(:,ke) = DENS(:) - DENS_hyd(:,ke)
      DRHOT(:,ke) = DENS(:) * PT(:) &
                  - PRES00 / Rdry * ( PRES_hyd(:,ke) / PRES00 )**(CVdry/CPdry)
    end do
    end do
    end do

    return
  end subroutine exp_SetInitCond_rising_therm_bubble

!OCL SERIAL
  subroutine exp_geostrophic_balance_correction( this,                   &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
    lcmesh, elem )
    
    implicit none

    class(Exp_rising_therm_bubble_moist), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(inout) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DRHOT(elem%Np,lcmesh%NeA)

    !---------------------------------------------------
    return
  end subroutine exp_geostrophic_balance_correction 

end module mod_user
