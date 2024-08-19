!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for moist Held and Suarez test proposed by Thatcher and  Jablonowski (2016).
!!          
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
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  & 
    CVdry => CONST_CVdry,  &      
    PRES00 => CONST_PRE00, &
    Grav => CONST_GRAV, &
    OHM => CONST_OHM,   &
    RPlanet => CONST_RADIUS, &
    PI => CONST_PI
  
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

  use mod_user_base, only: UserBase
  use mod_experiment, only: Experiment

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

  type(MeshField2D) :: SfcTemp

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
    implicit none
    class(User), intent(inout) :: this    
    class(AtmosComponent), intent(inout) :: atm

    logical :: USER_do                   = .false. !< do user step?
    namelist / PARAM_USER / &
       USER_do

    integer :: ierr    

    integer :: n, ke2D
    class(MeshBase2D), pointer :: mesh2D
    class(LocalMesh3D), pointer :: lcmesh3D
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

    call this%UserBase%Setup( atm, USER_do )

    !--
    call atm%mesh%ptr_mesh%GetMesh2D( mesh2D )
    call SfcTemp%Init( "sfc_temp", "K", mesh2D )

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      lcmesh3D => atm%mesh%ptr_mesh%lcmesh_list(n)

      SfcTemp%local(n)%val(:,ke2D) = &
          Ts_DelT * exp( - 0.5_RP * ( lcmesh3D%lat2D(:,ke2D) / Ts_DelLat )**2 ) &
        + Ts_Tmin
    end do

    !--
    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    use scale_time_manager, only:  TIME_NOWSTEP
    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars   
    use scale_prc 
    implicit none

    class(User), intent(inout) :: this 
    class(AtmosComponent), intent(inout) :: atm

  end subroutine USER_calc_tendency

!OCL SERIAL
  subroutine USER_update( this, atm )
    use scale_localmeshfield_base, only: LocalMeshFieldBase

    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      MOMX_p  => PHYTEND_MOMX_ID, &
      MOMY_p  => PHYTEND_MOMY_ID, &
      MOMZ_p  => PHYTEND_MOMZ_ID, &
      RHOH_p  => PHYTEND_RHOH_ID

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars

    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
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
        T(:) =  PRES%val(:,ke) / ( Rdry * DENS(:) )

        MOMX%val(:,ke) = MOMX%val(:,ke) / ( 1.0_RP + dt * rtauV )
        MOMY%val(:,ke) = MOMY%val(:,ke) / ( 1.0_RP + dt * rtauV )

        !-- For the case of d DRHOT / dt = dens * Cv * ( Teq - T ) / tauT <- ( Ullrich and Jablonowski (2012, JCP) )
        ! DRHOT%val(:,ke) = DRHOT%val(:,ke) &
        !                 - dt * rtauT(:) / gamm * ( 1.0_RP - Teq(:) / T(:) ) * DENS(:) * PT%val(:,ke)     &
        !                 / ( 1.0_RP + dt * rtauT(:) / gamm * ( 1.0_RP + (gamm - 1.0_RP) * Teq(:) / T(:) ) )

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
      TRACER_inq_id

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
    
    ! namelist /PARAM_EXP/

    integer :: ke, ke2D
    integer :: p, p_h, p_v   
    integer :: ierr

    real(RP) :: DENS_ini(elem%Np)
    real(RP) :: PRES_ini(elem%Np)
    real(RP) :: MOMX_met(elem%Np,lcmesh%Ne)
    real(RP) :: MOMY_met(elem%Np,lcmesh%Ne)

    integer :: iq
    !-----------------------------------------------------------------------------

    ! rewind(IO_FID_CONF)
    ! read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    ! if( ierr < 0 ) then !--- missing
    !    LOG_INFO("Held_Suarez_setup",*) 'Not found namelist. Default used.'
    ! elseif( ierr > 0 ) then !--- fatal error
    !    LOG_ERROR("Held_Suarez_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
    !    call PRC_abort
    ! endif
    ! LOG_NML(PARAM_EXP)

    !---
    call hydrostatic_calc_basicstate_constT( DENS_hyd, PRES_hyd, & ! (out)
      TEMP0, PRES00,                                             & ! (in)
      x, y, z, lcmesh, elem                                      ) ! (in)

    !--
    call TRACER_inq_id( "QV", iq )
!    tracer_field_list(iq)%ptr%val(:,ke),         &

    !$omp parallel do private( p_h, p_v, p, ke2D, &
    !$omp DENS_ini, PRES_ini                      )
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)
      do p_v=1, elem%Nnode_v
      do p_h=1, elem%Nnode_h1D**2
        p = p_h + (p_v-1)*elem%Nnode_h1D**2

        call calc_initial_state_1pt( &
          lcmesh%lon2D(p_h,ke2D), lcmesh%lat2D(p_h,ke2D), lcmesh%zlev(p,ke),            & ! (in)
          DENS_ini(p), PRES_ini(p), MOMX_met(p,ke), tracer_field_list(iq)%ptr%val(p,ke) ) ! (out)
      end do
      end do

      DDENS(:,ke) = DENS_ini(:) - DENS_hyd(:,ke)
      DRHOT(:,ke) = PRES00 / Rdry * (  ( PRES_ini   (:) / PRES00 )**(CVdry/CPdry) &
                                     - ( PRES_hyd(:,ke) / PRES00 )**(CVdry/CPdry) )
      MOMY_met(:,ke) = 0.0_RP
    end do

    call CubedSphereCoordCnv_LonLat2CSVec( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),    & ! (in)
      lcmesh%gam(:,lcmesh%NeS:lcmesh%NeE), elem%Np * lcmesh%Ne,      & ! (in)
      MOMX_met(:,:), MOMY_met(:,:),                                  & ! (in)
      MOMX(:,lcmesh%NeS:lcmesh%NeE), MOMY(:,lcmesh%NeS:lcmesh%NeE)   ) ! (out)

    return
  end subroutine exp_SetInitCond_Held_Suarez

  subroutine calc_initial_state_1pt( lon, lat, zlev, &
    dens, pres, momx_met, qv )
    implicit none
    real(RP), intent(in) :: lon, lat, zlev
    real(RP), intent(out) :: dens, pres, momx_met
    real(RP), intent(out) :: qv

    integer, parameter :: b_hw = 2 ! Half-width parameter
    real(RP), parameter :: Gamma  = 5E-3_RP   ! Lapse rate [K/m]
    integer, parameter :: k  = 3         ! Power used for temperature field
    real(RP), parameter :: Qv0 = 18E-3_RP ! Maximum specific humidity
    real(RP), parameter :: Lat_hw = 2.0_RP * PI / 9.0_RP !  Horizontal half-width of the specific humidity profile with latitude
    real(RP), parameter :: Pres_hw = 3E4_RP ! Vertical half-width of the specific humidity profile with pressure.
    real(RP) :: H ! Scale height
    real(RP) :: A, B, C
    real(RP) :: fac1, fac2, fac3
    real(RP) :: tau1, tau2
    real(RP) :: int_tau1, int_tau2
    real(RP) :: U
    real(RP) :: umet
    real(RP) :: temp_v
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
    temp_v = 1.0_RP / ( tau1 - tau2 * fac1 )
    pres = PRES00 * exp( Grav / Rdry * ( - int_tau1 + int_tau2 * fac3 ) )
    
    U = Grav / RPlanet * k * int_tau2 * ( cos(lat)**(k-1) - cos(lat)**(k+1) ) * temp_v
    umet = - OHM * RPlanet * cos(lat) &
           + sqrt( ( OHM * RPlanet * cos(lat) )**2 - RPlanet * cos(lat) * U )

    dens = pres / ( Rdry * temp_v )
    momx_met = dens * umet

    qv = 0.0_RP
    if ( pres > 1E4_RP ) then
      qv = Qv0 * exp( - ( lat / Lat_hw )**4 ) * exp( - ( ( pres - PRES00 ) / Pres_hw )**2 ) 
    end if

    return
  end subroutine calc_initial_state_1pt

end module mod_user
