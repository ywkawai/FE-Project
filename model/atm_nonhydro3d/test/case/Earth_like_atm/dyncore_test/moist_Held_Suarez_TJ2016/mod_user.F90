!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for moist Held and Suarez test proposed by Thatcher and Jablonowski (2016).
!!          
!! @author Team SCALE
!!
!! @par Reference
!!  - Thatcher and  Jablonowski 2016:
!!    A moist aquaplanet variant of the Held–Suarez test for atmospheric model dynamical cores
!!    Geoscientific Model Development, 9, 1263–1292
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
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
    Rdry => CONST_Rdry,    &
    Rvap => CONST_Rvap,    &
    CPdry => CONST_CPdry,  & 
    CVdry => CONST_CVdry,  & 
    LHV0 => CONST_LHV0,    &
    PRES00 => CONST_PRE00, &
    Grav => CONST_GRAV, &
    OHM => CONST_OHM,   &
    RPlanet => CONST_RADIUS, &
    PI => CONST_PI
  use scale_tracer, only: &
    TRACER_inq_id
  
  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_line, only: LineElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_element_modalfilter, only: ModalFilter

  use mod_user_base, only: UserBase
  use mod_experiment, only: Experiment

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    DENS_p  => PHYTEND_DENS_ID, &
    MOMX_p  => PHYTEND_MOMX_ID, &
    MOMY_p  => PHYTEND_MOMY_ID, &
    MOMZ_p  => PHYTEND_MOMZ_ID, &
    RHOT_p  => PHYTEND_RHOT_ID, &
    RHOH_p  => PHYTEND_RHOH_ID

  use mod_atmos_vars, only: &
    AtmosVars_GetLocalMeshPrgVars,    &
    AtmosVars_GetLocalMeshPhyAuxVars, &
    AtmosVars_GetLocalMeshQTRC_Qv

  use mod_user_sub_LSC, only: &
    USER_sub_LSC_Init, USER_sub_LSC_calc_tendency
  use mod_user_sub_BLmixing, only: &
    USER_sub_BLmixing_Init, USER_sub_BLmixing_calc_tendency
  use mod_user_sub_Filter, only: &
    Filter
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(UserBase) :: User
  contains
    procedure :: mkinit_ => USER_mkinit
    generic :: mkinit => mkinit_
    procedure :: setup_ => USER_setup
    generic :: setup => setup_
    procedure :: calc_tendency => USER_calc_tendency
    procedure :: update => USER_update
  end type User

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
  real(RP), private, parameter :: kf          = 1.0_RP / ( 86400.0_RP * 1.0_RP  )
  real(RP), private, parameter :: ka          = 1.0_RP / ( 86400.0_RP * 40.0_RP )
  real(RP), private, parameter :: ks          = 1.0_RP / ( 86400.0_RP * 4.0_RP  )
  real(RP), private, parameter :: TEMP_strato = 200.0_RP
  real(RP), private, parameter :: SFCTEMP_eq  = 294.0_RP
  real(RP), private, parameter :: DelT_y      = 65.0_RP
  real(RP), private, parameter :: DelPT_z     = 10.0_RP
  real(RP), private, parameter :: sigb        = 0.7_RP
  real(RP), private, parameter :: Ts_DelT     = 29.0_RP
  real(RP), private, parameter :: Ts_DelLat   = 26.0_RP * PI / 180.0_RP
  real(RP), private, parameter :: Ts_Tmin     = 271.0_RP

  real(RP), parameter :: Temp0_E = 310.0_RP ! Surface equatorial temperature
  real(RP), parameter :: Temp0_P = 240.0_RP ! Surface polar temperature
  real(RP), parameter :: Temp0 = 0.5_RP * ( Temp0_E + Temp0_P )
  !-----------------------------------------------------------------------------

  logical :: APPLY_NewFilter
  type(Filter) :: newFilter

contains

!OCL SERIAL
  subroutine USER_mkinit ( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init( 'Held_Suarez' )
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_Held_Suarez )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL  
  subroutine USER_setup( this, atm )
    use scale_tracer, only: &
      TRACER_regist
    use mod_atmos_phy_sfc_vars, only: &
      ATMOS_PHY_SF_SVAR_TEMP_ID
    use scale_polynominal
    implicit none
    class(User), intent(inout) :: this    
    class(AtmosComponent), intent(inout), target :: atm

    logical :: USER_do                   = .false. !< do user step?

    character(len=H_SHORT) :: FilterShape
    real(RP) :: FilterWidthFac
    namelist / PARAM_USER / &
       USER_do, &
       APPLY_NewFilter, FilterShape, FilterWidthFac
    integer :: ierr    

    integer :: n, ke2D
    class(MeshBase2D), pointer :: mesh2D
    class(LocalMesh3D), pointer :: lcmesh3D
    class(LocalMesh2D), pointer :: lcmesh2D

    type(MeshField2D), pointer :: SfcTemp

    real(RP), allocatable :: gl_pts(:)
    integer :: p1,p2,p3,p_
    real(RP), allocatable :: invV_PordM1(:,:)
    class(ElementBase3D), pointer :: elem3D

    integer :: iv
    !------------------------------------------


    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    APPLY_NewFilter = .false.
    FilterShape     = "GAUSSIAN"
    FilterWidthFac  = 1.0_RP

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

    call this%UserBase%Setup( atm, USER_do )

    !--
    call atm%mesh%ptr_mesh%GetMesh2D( mesh2D )
    SfcTemp => atm%phy_sfc_proc%vars%SFC_VARS(ATMOS_PHY_SF_SVAR_TEMP_ID)

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      lcmesh3D => atm%mesh%ptr_mesh%lcmesh_list(n)
      lcmesh2D => lcmesh3D%lcmesh2D
      !$omp parallel do private(ke2D)
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
        SfcTemp%local(n)%val(:,ke2D) = &
            Ts_DelT * exp( - 0.5_RP * ( lcmesh3D%lat2D(:,ke2D) / Ts_DelLat )**2 ) &
          + Ts_Tmin
      end do
    end do

    !- Setup LSC 
    call USER_sub_LSC_Init( atm%mesh%ptr_mesh )
    !- Setup BL mixing
    call USER_sub_BLmixing_Init( atm%mesh%ptr_mesh )

    !-
    call newFilter%Init( FilterShape, FilterWidthFac, atm%mesh%ptr_mesh )
    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    use scale_time_manager, only:  TIME_NOWSTEP
    use scale_prc 
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in    
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      DENS_p  => PHYTEND_DENS_ID, &
      MOMX_p  => PHYTEND_MOMX_ID, &
      MOMY_p  => PHYTEND_MOMY_ID, &
      MOMZ_p  => PHYTEND_MOMZ_ID, &
      RHOT_p  => PHYTEND_RHOT_ID, &
      RHOH_p  => PHYTEND_RHOH_ID      
    implicit none
    class(User), intent(inout) :: this 
    class(AtmosComponent), intent(inout) :: atm

    real(RP) :: dt
    !----------------------------

    dt = atm%time_manager%dtsec

    if ( APPLY_NewFilter ) then
      call newFilter%Apply( atm%vars%PHY_TEND(DENS_p), atm%mesh%ptr_mesh )
      call newFilter%Apply( atm%vars%PHY_TEND(MOMX_p), atm%mesh%ptr_mesh )
      call newFilter%Apply( atm%vars%PHY_TEND(MOMY_p), atm%mesh%ptr_mesh )
      call newFilter%Apply( atm%vars%PHY_TEND(MOMZ_p), atm%mesh%ptr_mesh )
      call newFilter%Apply( atm%vars%PHY_TEND(RHOT_p), atm%mesh%ptr_mesh )
      call newFilter%Apply( atm%vars%PHY_TEND(RHOH_p), atm%mesh%ptr_mesh )
      call newFilter%Apply( atm%vars%PHY_TEND(RHOH_p+1), atm%mesh%ptr_mesh )
    end if

    !-- Large-scale condensation
    call USER_sub_LSC_calc_tendency( atm%vars, &
      atm%mesh%ptr_mesh, dt )
    
    !-- Boundary layer mixing
    call USER_sub_BLmixing_calc_tendency( atm%vars, &
      atm%mesh%DOptrMat(3), atm%mesh%LiftOptrMat, atm%mesh%ptr_mesh )

    return
  end subroutine USER_calc_tendency

!OCL SERIAL
  subroutine USER_update( this, atm )
    use scale_localmeshfield_base, only: LocalMeshFieldBase
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT, QV, RHOQv_tp
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    type(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke2D

    real(RP), allocatable :: DENS(:), T(:), Teq(:), sig(:), PRES_sfc(:)
    real(RP), allocatable :: rtauT(:), rtauV(:)
    real(RP), allocatable :: lat(:)

    real(DP) :: dt
    real(RP) :: Gamm
    !----------------------------------------------------------

    dt = atm%time_manager%dtsec
    gamm = CpDry / CvDry 

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      

      call AtmosVars_GetLocalMeshPhyAuxVars( n, atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                         )
      
      call AtmosVars_GetLocalMeshQTRC_Qv( n, atm%mesh%ptr_mesh, &
        atm%vars%QTRCVARS_manager, atm%vars%PHYTENDS_manager, QV, RHOQv_tp )

      elem3D => lcmesh%refElem3D

      allocate( DENS(elem3D%Np), T(elem3D%Np), Teq(elem3D%Np), sig(elem3D%Np) )
      allocate( PRES_sfc(elem3D%Nnode_h1D**2) )
      allocate( rtauT(elem3D%Np), rtauV(elem3D%Np) )
      allocate( lat(elem3D%Np) )

      !$omp parallel do private(DENS, T, Teq, PRES_sfc, sig, lat, rtauT, rtauV, ke2D)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)

        PRES_sfc(:) = PRES%val(elem3D%Hslice(:,1),ke2D)
        sig(:) = PRES%val(:,ke) / PRES_sfc(elem3D%IndexH2Dto3D)
        lat(:) = lcmesh%lat2D(elem3D%IndexH2Dto3D,ke2D)

        rtauT(:) = ka + (ks - ka) * max( 0.0_RP, (sig(:) - sigb)/(1.0_RP - sigb) ) * cos(lat(:))**4
        rtauV(:) = kf * max( 0.0_RP, (sig(:) - sigb)/(1.0_RP - sigb) )

        Teq(:) = max( TEMP_strato, &
          ( SFCTEMP_eq - DelT_y * sin(lat(:))**2 - DelPT_z * log(PRES%val(:,ke)/PRES00) * cos(lat(:))**2 ) &
          * (PRES%val(:,ke)/PRES00)**(Rdry/CPDry)                                                        )
        
        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)
        T(:) =  PRES%val(:,ke) / ( Rtot%val(:,ke) * DENS(:) )

        MOMX%val(:,ke) = MOMX%val(:,ke) / ( 1.0_RP + dt * rtauV )
        MOMY%val(:,ke) = MOMY%val(:,ke) / ( 1.0_RP + dt * rtauV )

        !--
        !- For the case of d DRHOT /dt = dens * Cp * ( Teq - T ) / tauT     
        !  <- It is based on the forcing form in Held and Surez in which the temperature evolution equation is assumed to be dT/dt = R/Cp * T/p * dp/dt + (Teq - T) / tauT
        DRHOT%val(:,ke) = DRHOT%val(:,ke) &
                        - dt * rtauT(:) * ( 1.0_RP - Teq(:) / T(:) ) * DENS(:) * PT%val(:,ke)     &
                        / ( 1.0_RP + dt * rtauT(:) * ( 1.0_RP + (gamm - 1.0_RP) * Teq(:) / T(:) ) )  
      end do

      deallocate( DENS, T, Teq, sig, PRES_sfc )
      deallocate( rtauT, rtauV )
      deallocate( lat )
    end do
    call atm%vars%Calc_diagnostics()

    return
  end subroutine USER_update

  !------

!OCL SERIAL
  subroutine exp_SetInitCond_Held_Suarez( this,                            &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use mod_experiment, only: &
      TracerLocalMeshField_ptr

    use scale_const, only: &
      PRES00 => CONST_PRE00, &
      RPlanet => CONST_Radius, &
      CPdry => CONST_CPdry, &     
      CVdry => CONST_CVdry

    use scale_tracer, only: &
      TRACER_CV, TRACER_CP
    use scale_atmos_hydrometeor, only: &
      LHV, &
      I_QV            
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSVec
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constT
    implicit none

    class(Experiment), intent(inout) :: this
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

    real(RP) :: Qv_max = 18E-3_RP ! Maximum specific humidity
    real(RP) :: Vp     = 1.0_RP 
    
    namelist /PARAM_EXP/ &
      Vp, Qv_max

    integer :: ke, ke2D
    integer :: p, p_h, p_v   
    integer :: ierr

    real(RP) :: DENS_ini(elem%Np)
    real(RP) :: PRES_ini(elem%Np)
    real(RP) :: QV_ini(elem%Np)
    real(RP) :: MOMX_met(elem%Np,lcmesh%Ne)
    real(RP) :: MOMY_met(elem%Np,lcmesh%Ne)
    real(RP) :: CPtot(elem%Np), CVtot(elem%Np)

    integer :: iq
    !-----------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("moist_Held_Suarez_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("moist_Held_Suarez_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    call hydrostatic_calc_basicstate_constT( DENS_hyd, PRES_hyd, & ! (out)
      TEMP0, PRES00,                                             & ! (in)
      x, y, z, lcmesh, elem                                      ) ! (in)

    !--
    call TRACER_inq_id( "QV", iq )
!    tracer_field_list(iq)%ptr%val(:,ke),         &

    !$omp parallel do private( p_h, p_v, p, ke2D, &
    !$omp DENS_ini, PRES_ini, QV_ini, CPtot, CVtot        )
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)
      do p_v=1, elem%Nnode_v
      do p_h=1, elem%Nnode_h1D**2
        p = p_h + (p_v-1)*elem%Nnode_h1D**2

        call calc_initial_state_1pt( &
          lcmesh%lon2D(p_h,ke2D), lcmesh%lat2D(p_h,ke2D), lcmesh%zlev(p,ke), Vp, Qv_max, & ! (in)
          DENS_ini(p), PRES_ini(p), MOMX_met(p,ke), MOMY_met(p,ke),              & ! (out)
          QV_ini(p) ) ! (out)
      end do
      end do

      tracer_field_list(iq)%ptr%val(:,ke) = QV_ini(:)

      CVtot(:) = CVdry + ( TRACER_CV(I_QV) - CVdry ) * QV_ini(:)
      CPtot(:) = CPdry + ( TRACER_CP(I_QV) - CPdry ) * QV_ini(:)
      DDENS(:,ke) = DENS_ini(:) - DENS_hyd(:,ke)
      DRHOT(:,ke) = PRES00 &
        * ( ( PRES_ini   (:) / PRES00 )**(CVtot(:)/CPtot(:)) / ( CPtot(:) - CVtot(:) ) &
          - ( PRES_hyd(:,ke) / PRES00 )**(CVdry/CPdry) / Rdry )
    end do

    call CubedSphereCoordCnv_LonLat2CSVec( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),    & ! (in)
      lcmesh%gam(:,lcmesh%NeS:lcmesh%NeE), elem%Np * lcmesh%Ne,      & ! (in)
      MOMX_met(:,:), MOMY_met(:,:),                                  & ! (in)
      MOMX(:,lcmesh%NeS:lcmesh%NeE), MOMY(:,lcmesh%NeS:lcmesh%NeE)   ) ! (out)

    return
  end subroutine exp_SetInitCond_Held_Suarez

  subroutine calc_initial_state_1pt( lon, lat, zlev, vp_perturb, Qv0, &
    dens, pres, momx_met, momy_met, qv )
    implicit none
    real(RP), intent(in) :: lon, lat, zlev
    real(RP), intent(in) :: vp_perturb, Qv0
    real(RP), intent(out) :: dens, pres, momx_met, momy_met
    real(RP), intent(out) :: qv

    integer, parameter :: b_hw = 2 ! Half-width parameter
    real(RP), parameter :: Gamma  = 5E-3_RP   ! Lapse rate [K/m]
    integer, parameter :: k  = 3         ! Power used for temperature field
    real(RP), parameter :: Lat_hw = 2.0_RP * PI / 9.0_RP !  Horizontal half-width of the specific humidity profile with latitude
    real(RP), parameter :: Pres_hw = 3E4_RP ! Vertical half-width of the specific humidity profile with pressure.

    real(RP), parameter :: lonc = PI / 9.0_RP
    real(RP), parameter :: latc = 2.0_RP * PI / 9.0_RP
    real(RP), parameter :: d0 = 1.0_RP / 6.0_RP
    real(RP), parameter :: zt = 1.5E4_RP ! Top of perturbation domain

    real(RP) :: H ! Scale height
    real(RP) :: A, B, C
    real(RP) :: fac1, fac2, fac3
    real(RP) :: tau1, tau2
    real(RP) :: int_tau1, int_tau2
    real(RP) :: U
    real(RP) :: umet
    real(RP) :: temp_v
    real(RP) :: umet_dash, vmet_dash
    real(RP) :: dist
    !------------------------------------------------------
    
    A = 1.0_RP / Gamma
    B = ( Temp0_E - Temp0_P ) / ( ( Temp0_E + Temp0_P ) * Temp0_P )
    C = 0.5_RP * dble( k + 2 ) * ( Temp0_E - Temp0_P ) / ( Temp0_E * Temp0_P )
    H = Rdry * Temp0 / Grav

    fac1 = exp( ( Gamma / TEMP0 ) * zlev )
    fac2 = exp( - ( zlev / ( b_hw * H ) )**2 )
    fac3 = 1.0_RP - 2.0_RP * ( zlev / ( b_hw * H ) )**2
    tau1 = 1.0_RP / Temp0 * fac1 &
         + B * fac3 * fac2
    tau2 = C * fac3 * fac2
    int_tau1 = A * ( fac1 - 1.0_RP ) &
             + B * zlev * fac2
    int_tau2 = C * zlev * fac2

    fac3 = cos(lat)**k - ( dble(k) / dble(k+2) ) * cos(lat)**(k+2)
    temp_v = 1.0_RP / ( tau1 - tau2 * fac3 )
    pres = PRES00 * exp( Grav / Rdry * ( - int_tau1 + int_tau2 * fac3 ) )
    
    U = Grav / RPlanet * k * int_tau2 * ( cos(lat)**(k-1) - cos(lat)**(k+1) ) * temp_v
    umet = - OHM * RPlanet * cos(lat) &
           + sqrt( ( OHM * RPlanet * cos(lat) )**2 + RPlanet * cos(lat) * U )

    dens = pres / ( Rdry * temp_v )

    !--
    if ( zlev/zt < 1.0_RP ) then
      fac1 = 16.0_RP * vp_perturb / ( 3.0_RP * sqrt(3.0_RP) ) &
        * ( 1.0_RP - 3.0_RP * ( zlev/zt )**2 + 2.0_RP * ( zlev/zt )**3 )
    else
      fac1 = 0.0_RP
    end if
    dist = acos( sin(latc) * sin(lat) + cos(latc) * cos(lat) * cos(lon-lonc) )
    fac2 = cos( 0.5_RP * PI * dist / d0 )**3 * sin( 0.5_RP * PI * dist / d0 )

    umet_dash = - fac1 * fac2 &
      * ( -sin(latc)*cos(lat) + cos(latc)*sin(lat)*cos(lon-lonc) ) / sin(dist)
    vmet_dash = fac1 * fac2 &
      * cos(lonc) * sin(lon-lonc) / sin(dist)
    
    !--
    momx_met = dens * ( umet + umet_dash )
    momy_met = dens * (        vmet_dash )

    !--
    qv = 0.0_RP
    if ( pres > 1E4_RP ) then
      qv = Qv0 * exp( - ( lat / Lat_hw )**4 ) * exp( - ( ( pres - PRES00 ) / Pres_hw )**2 ) 
    end if

    return
  end subroutine calc_initial_state_1pt

end module mod_user
