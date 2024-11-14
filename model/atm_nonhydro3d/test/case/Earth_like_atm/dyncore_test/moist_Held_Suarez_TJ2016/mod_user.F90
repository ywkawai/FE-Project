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

  real(RP) :: BL_CE_PT = 0.0044_RP
  real(RP) :: BL_CE_QV = 0.0044_RP


  type(LineElement) :: elem1D
  type(MeshField2D) :: RAIN_LSC
  type(MeshField3D) :: CondensRate, CondensRate_ori

  integer :: GL1D_npts, LGL1D_npts
  real(RP), allocatable :: PhyIntrpMat1D(:,:) 
  real(RP), allocatable :: PhyIntrpMat1D_tr(:,:) 

  real(RP), allocatable :: PhyIntrpMat3D(:,:) 
  type(ModalFilter) :: PhyTendFilter

  integer :: LSC_nstep
  integer :: LSC_rnstep
  type(MeshField3D) :: LSC_TEND(6)
  !-----------------------------------------------------------------------------
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
    integer :: LSC_MF_ORDER_h
    real(RP) :: LSC_MF_ALPHA_h
    integer :: LSC_MF_ORDER_v
    real(RP) :: LSC_MF_ALPHA_v

    namelist / PARAM_USER / &
       USER_do, &
       LSC_nstep, &
       LSC_MF_ORDER_h, LSC_MF_ALPHA_h, &
       LSC_MF_ORDER_v, LSC_MF_ALPHA_v, &
       BL_CE_PT, BL_CE_QV

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

    LSC_nstep = 1
    LSC_MF_ORDER_h = 16; LSC_MF_ALPHA_h = 0.0_RP
    LSC_MF_ORDER_v = 16; LSC_MF_ALPHA_v = 0.0_RP

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

    !--
    call elem1D%Init( lcmesh3D%refElem3D%PolyOrder_v, .false. )
    call RAIN_LSC%Init( "RAIN_LSC", "kg/s", mesh2D )
    call CondensRate%Init( "CondensRate_LSC", "kg/m3.s-1", atm%mesh%ptr_mesh )
    call CondensRate_ori%Init( "CondensRate_LSC_ori", "kg/m3.s-1", atm%mesh%ptr_mesh )

    LGL1D_npts = elem1D%Np; GL1D_npts = LGL1D_npts
    allocate( PhyIntrpMat1D(GL1D_npts,LGL1D_npts), PhyIntrpMat1D_tr(LGL1D_npts,GL1D_npts) )
    allocate( gl_pts(GL1D_npts) )
    gl_pts(:) = Polynominal_GenGaussLegendrePt(GL1D_npts)
    PhyIntrpMat1D(:,:) = Polynominal_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, gl_pts )
    PhyIntrpMat1D_tr(:,:) = transpose(PhyIntrpMat1D)

    elem3D => lcmesh3D%refElem3D
    allocate( PhyIntrpMat3D(elem3D%Np,elem3D%Np) )
    allocate( invV_PordM1(elem3D%Np,elem3D%Np) )

    invV_PordM1(:,:) = elem3D%InvV(:,:)
    do p3=1, elem3D%Nnode_v
    do p2=1, elem3D%Nnode_h1D
    do p1=1, elem3D%Nnode_h1D
      if (p3 == elem3D%Nnode_v .or. p2 >= elem3D%Nnode_h1D-1 .or. p1 >= elem3D%Nnode_h1D-1  ) then
        p_ = p1 + (p2-1)*elem3D%Nnode_h1D + (p3-1)*elem3D%Nnode_h1D**2
        invV_PordM1(p_,:) = 0.0_RP
      end if
    end do
    end do
    end do
    PhyIntrpMat3D(:,:) = matmul(elem3D%V, invV_PordM1)

    call PhyTendFilter%Init(atm%mesh%element, 0.0_RP, LSC_MF_ALPHA_h, LSC_MF_ORDER_h, 0.0_RP, LSC_MF_ALPHA_v, LSC_MF_ORDER_v)

    LSC_rnstep = 0
    do iv=1, 6
      call LSC_TEND(iv)%Init( "LSC_TEND", "", atm%mesh%ptr_mesh )
    end do

    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    use scale_time_manager, only:  TIME_NOWSTEP
    use scale_prc 
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in    
    implicit none

    class(User), intent(inout) :: this 
    class(AtmosComponent), intent(inout) :: atm

    integer :: n
    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT, QV
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: RHOQv_tp
    type(ElementBase3D), pointer :: elem3D

    real(RP) :: dt
    integer :: ke
    ! real(RP) :: mass_check
    !----------------------------

    dt = atm%time_manager%dtsec

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      

      call AtmosVars_GetLocalMeshPhyAuxVars( n, atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                         )

      call AtmosVars_GetLocalMeshQTRC_Qv( n, atm%mesh%ptr_mesh, &
        atm%vars%QTRCVARS_manager, atm%vars%PHYTENDS_manager, QV, RHOQv_tp )

      !--
      if ( LSC_rnstep == 0 ) then
        call Large_scale_Precip( &
          LSC_TEND(1)%local(n)%val, LSC_TEND(2)%local(n)%val, LSC_TEND(3)%local(n)%val, LSC_TEND(4)%local(n)%val, &
!          atm%vars%PHY_TEND(DENS_p)%local(n)%val, atm%vars%PHY_TEND(MOMX_p)%local(n)%val, atm%vars%PHY_TEND(MOMY_p)%local(n)%val, atm%vars%PHY_TEND(MOMZ_p)%local(n)%val, &
          LSC_TEND(5)%local(n)%val, LSC_TEND(6)%local(n)%val, & 
!          RHOQv_tp%val, atm%vars%PHY_TEND(RHOH_p)%local(n)%val, &
          RAIN_LSC%local(n)%val, CondensRate%local(n)%val,                                  &
          QV%val, DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DENS_hyd%val, PRES%val, Rtot%val, CVtot%val, CPtot%val, &
          real(LSC_nstep, kind=RP) * dt, lcmesh, lcmesh%refElem3D )
        LSC_rnstep = LSC_nstep
      end if
      !$omp parallel do
      do ke=lcmesh%NeS, lcmesh%NeE
        atm%vars%PHY_TEND(DENS_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(DENS_p)%local(n)%val(:,ke) + LSC_TEND(1)%local(n)%val(:,ke)
        atm%vars%PHY_TEND(MOMX_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMX_p)%local(n)%val(:,ke) + LSC_TEND(2)%local(n)%val(:,ke)
        atm%vars%PHY_TEND(MOMY_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMY_p)%local(n)%val(:,ke) + LSC_TEND(3)%local(n)%val(:,ke)
        atm%vars%PHY_TEND(MOMZ_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMZ_p)%local(n)%val(:,ke) + LSC_TEND(4)%local(n)%val(:,ke)
        RHOQv_tp%val(:,ke) = RHOQv_tp%val(:,ke) + LSC_TEND(5)%local(n)%val(:,ke)
        atm%vars%PHY_TEND(RHOH_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(RHOH_p)%local(n)%val(:,ke) + LSC_TEND(6)%local(n)%val(:,ke)
      end do
      LSC_rnstep = LSC_rnstep - 1

      ! mass_check = 0.0_RP
      ! elem3D => lcmesh%refElem3D
      ! !$omp parallel do reduction(+: mass_check)
      ! do ke=lcmesh%NeS, lcmesh%NeE
      !   mass_check = mass_check &
      !     + sum( lcmesh%J(:,ke) * lcmesh%Gsqrt(:,ke) * elem3D%IntWeight_lgl(:) * LSC_TEND(5)%local(n)%val(:,ke) )
      ! end do
      ! LOG_INFO("LSC_mass_check: ",*) mass_check

      call BLmixing( atm%vars%PHY_TEND(RHOT_p)%local(n)%val, RHOQv_tp%val, &
        DDENS%val, MOMX%val, MOMY%val, MOMZ%val, QV%val,                            &
        DENS_hyd%val, PT%val, PRES%val, atm%mesh%DOptrMat(3), atm%mesh%LiftOptrMat, &
        lcmesh, lcmesh%refElem3D )
    end do

    call FILE_HISTORY_meshfield_in( RAIN_LSC, "RAIN with large scale condensation" )
    call FILE_HISTORY_meshfield_in( CondensRate, "Condensation rate (mass) with large scale condensation" )
    call FILE_HISTORY_meshfield_in( CondensRate_ori, "Condensation rate (mass) with large scale condensation" )

    return
  end subroutine USER_calc_tendency

  subroutine Large_scale_Precip( DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOQv_tp, RHOH_p, SFLX_RAIN, CondensRate, &
    QV, DDENS, MOMX, MOMY, MOMZ, DENS_hyd, PRES, Rtot, CVtot, CPtot, &
    DT, lcmesh, elem3D )
    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_pres2qsat_liq
    use scale_atmos_hydrometeor, only: &
      CV_WATER      
    use scale_const, only: &
      CL => CONST_CL    
    implicit none
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(inout) :: DENS_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMZ_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOQv_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOH_p(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: SFLX_RAIN(elem3D%Nnode_h1D**2,lcmesh%lcmesh2D%NeA)
    real(RP), intent(inout) :: CondensRate(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: QV(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DT

    integer :: ke, ke_xy, ke_z
    integer :: p2D
    real(RP) :: TEMP(elem3D%Np,lcmesh%NeA)
    real(RP) :: Qsat(elem3D%Np,lcmesh%NeA)
    real(RP) :: CondensRate_zxy(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D) !< Condensation rate
    real(RP) :: DENS(elem3D%Np)
    real(RP) :: rain_tmp(elem3D%Nnode_v)
    real(RP) :: coef1, coef2

    real(RP) :: work1(GL1D_npts,LGL1D_npts,LGL1D_npts)
    real(RP) :: work2(GL1D_npts,GL1D_npts,LGL1D_npts)
    real(RP) :: work(GL1D_npts**3)
    real(RP) :: work0(LGL1D_npts**3)
    real(RP) :: work_ip(GL1D_npts**3)

    real(RP) :: RHOE0(elem3D%Np), RHOE(elem3D%Np), CVtot1(elem3D%Np)
    !----------------------------------------------

    !$omp parallel do private(ke, work0,work,work1,work2,work_ip)
    do ke=lcmesh%NeS, lcmesh%NeE
      
      work0(:) = PRES(:,ke) / ( Rtot(:,ke) * ( DENS_hyd(:,ke) + DDENS(:,ke) ) )
      ! call interp_LGL2GL( work0, TEMP(:,ke), PhyIntrpMat1D, PhyIntrpMat1D_tr, work1, work2, LGL1D_npts, GL1D_npts )
      TEMP(:,ke) = work0(:)
      ! call interp_LGL2GL( PRES(:,ke), work_ip(:), PhyIntrpMat1D, PhyIntrpMat1D_tr, work1, work2, LGL1D_npts, GL1D_npts )
      work_ip(:) = PRES(:,ke)

      call ATMOS_SATURATION_pres2qsat_liq( elem3D%Np, 1, elem3D%Np, &
        TEMP(:,ke), work_ip(:), &
        Qsat(:,ke) )
    end do

    coef1 = LHV0**2 / ( CVDry * Rvap )
    coef2 = Rvap / LHV0

    !$omp parallel do private(ke, DENS, work_ip, RHOE, RHOE0, CVtot1) collapse(2)
    do ke_z=1, lcmesh%NeZ
    do ke_xy=1, lcmesh%Ne2D
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D
      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)

      ! call interp_LGL2GL( QV(:,ke), work_ip(:), PhyIntrpMat1D, PhyIntrpMat1D_tr, work1, work2, LGL1D_npts, GL1D_npts )
      work_ip(:) = QV(:,ke)

      CondensRate_ori%local(1)%val(:,ke) = DENS(:) &
        * max( 0.0_RP, ( work_ip(:) - Qsat(:,ke) ) / ( 1.0_RP + coef1 / TEMP(:,ke)**2 * ( 1.0_RP - coef2 * TEMP(:,ke) ) * Qsat(:,ke) ) / DT )
      ! CondensRate_zxy(:,ke_z,ke_xy) = matmul( PhyIntrpMat3D, CondensRate_zxy(:,ke_z,ke_xy) )
      ! CondensRate_zxy(:,ke_z,ke_xy) = matmul( PhyTendFilter%FilterMat, CondensRate_ori%local(1)%val(:,ke) )

      ! Qsat(:,ke) = matmul( PhyTendFilter, Qsat(:,ke) )
      ! work_ip(:) = matmul( PhyTendFilter, work_ip(:) )
      ! CondensRate_zxy(:,ke_z,ke_xy) = matmul( PhyTendFilter%FilterMat, DENS(:) * ( work_ip(:) - Qsat(:,ke) ) / ( 1.0_RP + coef1 / TEMP(:,ke)**2 * ( 1.0_RP - coef2 * TEMP(:,ke) ) * Qsat(:,ke) ) / DT )
      ! CondensRate_zxy(:,ke_z,ke_xy) = max( 0.0_RP, CondensRate_zxy(:,ke_z,ke_xy)  )
      CondensRate_zxy(:,ke_z,ke_xy) = CondensRate_ori%local(1)%val(:,ke)

      RHOH_p  (:,ke) = &
        LHV0 * CondensRate_zxy(:,ke_z,ke_xy) - CPtot(:,ke) * TEMP(:,ke) * CondensRate_zxy(:,ke_z,ke_xy)

      DENS_tp(:,ke) = - CondensRate_zxy(:,ke_z,ke_xy)
      RHOQv_tp(:,ke) = - CondensRate_zxy(:,ke_z,ke_xy)
      MOMX_tp(:,ke) = - CondensRate_zxy(:,ke_z,ke_xy) * MOMX(:,ke) / DENS(:)
      MOMY_tp(:,ke) = - CondensRate_zxy(:,ke_z,ke_xy) * MOMY(:,ke) / DENS(:)
      MOMZ_tp(:,ke) = - CondensRate_zxy(:,ke_z,ke_xy) * MOMZ(:,ke) / DENS(:)
      CondensRate(:,ke) = CondensRate_zxy(:,ke_z,ke_xy)
    end do
    end do

    SFLX_RAIN(:,:) = 0.0_RP
    !$omp parallel do private(ke_z,ke,p2D, rain_tmp, coef1)
    do ke_xy=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D
      do p2D=1, elem3D%Nnode_h1D**2
        rain_tmp(:) = CondensRate_zxy(elem3D%Colmask(:,p2D),ke_z,ke_xy)
        coef1 = 0.5_RP * ( lcmesh%zlev(elem3D%Colmask(elem3D%Nnode_v,p2D),ke) - lcmesh%zlev(elem3D%Colmask(1,p2D),ke) )
        SFLX_RAIN(p2D,ke_xy) = SFLX_RAIN(p2D,ke_xy) &
          + coef1 * sum(rain_tmp(:) * elem1D%IntWeight_lgl(:) )
      end do
    end do
    end do

    return
  end subroutine Large_scale_Precip

!OCL SERIAL
  subroutine BLmixing( RHOT_dt, RHOQV_dt, &
    DDENS, MOMX, MOMY, MOMZ, QV, &
    DENS_hyd, PT, PRES, &
    Dz, Lift, &
    lcmesh, elem3D )
    use scale_sparsemat, only: &
      sparsemat, sparsemat_matmul
    use scale_atm_dyn_dgm_nonhydro3d_rhot_hevi_common, only: &
      vi_gen_vmap => atm_dyn_dgm_nonhydro3d_rhot_hevi_common_gen_vmap      
    implicit none
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(inout) :: RHOT_dt(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOQV_dt(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: QV(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PT(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(elem3D%Np,lcmesh%NeA)
    class(SparseMat), intent(in) :: Dz, Lift

    real(RP) :: VDiffCoef(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: DiffFlux(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D,2)
    real(RP) :: PT_(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: QV_(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: nz(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: za_tmp(elem3D%Nnode_h1D**2)
    real(RP) :: za(elem3D%Np,lcmesh%Ne2D)
    real(RP) :: Vabs_tmp(elem3D%Nnode_h1D**2)
    real(RP) :: Vabs(elem3D%Np,lcmesh%Ne2D)

    integer :: ke, ke_xy, ke_z

    integer :: vmapM(elem3D%NfpTot,lcmesh%NeZ)
    integer :: vmapP(elem3D%NfpTot,lcmesh%NeZ)

    real(RP) :: del_flux(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,2)
    real(RP) :: Fz(elem3D%Np), LiftDelFlx(elem3D%Np)

    real(RP), parameter :: p_pbl = 850E2_RP
    real(RP), parameter :: p_strato = 100E2_RP

    real(RP) :: G_11(elem3D%Np), G_12(elem3D%Np), G_22(elem3D%Np)
    real(RP) :: DENS(elem3D%Np)
    integer :: hslice(elem3D%Nnode_h1D**2)

    ! real(RP) :: mass_check, mass_check2
    !-------------------------------------------------------------

    call vi_gen_vmap( vmapM, vmapP, lcmesh, elem3D )

    !$omp parallel private(ke, hslice, za_tmp, Vabs_tmp )

    !$omp do
    do ke_xy=1, lcmesh%Ne2D
      ! hslice(:) = elem3D%Hslice(:,elem3D%Nnode_v)
      ! za_tmp(:) = lcmesh%zlev(hslice(:),ke_xy) / dble(elem3D%Nnode_v)
      hslice(:) = elem3D%Hslice(:,2)
      za_tmp(:) = 60.0_RP !lcmesh%zlev(hslice(:),ke_xy)
      za(:,ke_xy) = za_tmp(elem3D%IndexH2Dto3D(:))

      hslice(:) = elem3D%Hslice(:,1)
      ke = ke_xy
      Vabs_tmp(:) = sqrt( &
          MOMX(hslice(:),ke) *  (lcmesh%G_ij(:,ke_xy,1,1) * MOMX(hslice(:),ke) + lcmesh%G_ij(:,ke_xy,2,1) * MOMY(hslice(:),ke)) &
        + MOMY(hslice(:),ke) *  (lcmesh%G_ij(:,ke_xy,2,1) * MOMX(hslice(:),ke) + lcmesh%G_ij(:,ke_xy,2,2) * MOMY(hslice(:),ke)) &
        + MOMZ(hslice(:),ke)**2 ) / ( DENS_hyd(hslice(:),ke) + DDENS(hslice(:),ke) )
      Vabs(:,ke_xy) = Vabs_tmp(elem3D%IndexH2Dto3D(:))
    end do

    !$omp do collapse(2)
    do ke_xy=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D

      where ( PRES(:,ke) > p_pbl )
        VDiffCoef(:,ke_z,ke_xy) = Vabs(:,ke_xy) * za(:,ke_xy)
      elsewhere
        VDiffCoef(:,ke_z,ke_xy) = Vabs(:,ke_xy) * za(:,ke_xy) * exp( - ((p_pbl-PRES(:,ke))/p_strato)**2 )
      endwhere
      
      PT_(:,ke_z,ke_xy) = PT(:,ke)
      QV_(:,ke_z,ke_xy) = QV(:,ke)
      nz (:,ke_z,ke_xy) = lcmesh%normal_fn(:,ke,3)
    end do
    end do

    !$omp end parallel

    call BLmixing_bnd_flux_grad( del_flux, &
      PT_, QV_, nz, vmapM, vmapP, lcmesh, elem3D )

    !$omp parallel do private(ke, Fz, LiftDelFlx, DENS) collapse(2)
    do ke_xy=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)

      call sparsemat_matmul(Dz, PT_(:,ke_z,ke_xy), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,1), LiftDelFlx)
      DiffFlux(:,ke_z,ke_xy,1) = &
        BL_CE_PT * DENS(:) * VDiffCoef(:,ke_z,ke_xy) * ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )

      call sparsemat_matmul(Dz, QV_(:,ke_z,ke_xy), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,2), LiftDelFlx)
      DiffFlux(:,ke_z,ke_xy,2) = &
        BL_CE_QV * DENS(:) * VDiffCoef(:,ke_z,ke_xy) * ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )        
    end do
    end do

    call BLmixing_bnd_flux( del_flux, &
      DiffFlux, nz, vmapM, vmapP, lcmesh, elem3D )

    !$omp parallel do private(ke, Fz, LiftDelFlx) collapse(2)
    do ke_xy=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D

      call sparsemat_matmul(Dz, DiffFlux(:,ke_z,ke_xy,1), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,1), LiftDelFlx)
      RHOT_dt(:,ke) = RHOT_dt(:,ke) &
        + ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )

      call sparsemat_matmul(Dz, DiffFlux(:,ke_z,ke_xy,2), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,2), LiftDelFlx)
      RHOQV_dt(:,ke) = RHOQV_dt(:,ke) &
        + ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )
    end do
    end do

    ! mass_check = 0.0_RP; mass_check2 = 0.0_RP
    ! !$omp parallel do private(ke, Fz, LiftDelFlx) reduction(+: mass_check, mass_check2)
    ! do ke_xy=1, lcmesh%Ne2D
    ! do ke_z=1, lcmesh%NeZ
    !   ke = ke_xy + (ke_z-1)*lcmesh%Ne2D
    !   call sparsemat_matmul(Dz, DiffFlux(:,ke_z,ke_xy,2), Fz)
    !   call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,2), LiftDelFlx)
    !   mass_check = mass_check &
    !     + sum( lcmesh%J(:,ke) * lcmesh%Gsqrt(:,ke) * elem3D%IntWeight_lgl(:) * ( lcmesh%Escale(:,ke,3,3) * ( Fz(:) ) ) )
    !   mass_check2 = mass_check2 &
    !     + sum( lcmesh%J(:,ke) * lcmesh%Gsqrt(:,ke) * elem3D%IntWeight_lgl(:) * ( LiftDelFlx(:) ) )
    ! end do
    ! end do
    ! LOG_INFO("BLM_mass_check: ++",*) mass_check + mass_check2      

    return
  end subroutine BLmixing

  subroutine BLmixing_bnd_flux_grad( bnd_flux, &
    PT, QV, nz, vmapM, vmapP, lcmesh, elem3D )
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: bnd_flux(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,2)
    real(RP), intent(in) :: PT(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: QV(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: nz(elem3D%NfpTot,lcmesh%NeZ,lcmesh%NeX*lcmesh%NeY)
    integer, intent(in) :: vmapM(elem3D%NfpTot,lcmesh%NeZ)
    integer, intent(in) :: vmapP(elem3D%NfpTot,lcmesh%NeZ)

    integer :: ke_z, ke2D, ke
    integer :: iM(elem3D%NfpTot), iP(elem3D%NfpTot)
    !--------------------------------------------------------

    !$omp parallel do private(iM, iP, ke2D, ke_z, ke) collapse(2)
    do ke2D=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)
      bnd_flux(:,ke_z,ke2D,1) = 0.5_RP * ( PT(iP,ke2D) - PT(iM,ke2D) ) * nz(:,ke_z,ke2D)
      bnd_flux(:,ke_z,ke2D,2) = 0.5_RP * ( QV(iP,ke2D) - QV(iM,ke2D) ) * nz(:,ke_z,ke2D)
    end do
    end do

    return
  end subroutine BLmixing_bnd_flux_grad

  subroutine BLmixing_bnd_flux( bnd_flux, &
    DiffFlux, nz, vmapM, vmapP, lcmesh, elem3D )
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: bnd_flux(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,2)
    real(RP), intent(in) :: DiffFlux(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D,2)
    real(RP), intent(in) :: nz(elem3D%NfpTot,lcmesh%NeZ,lcmesh%NeX*lcmesh%NeY)
    integer, intent(in) :: vmapM(elem3D%NfpTot,lcmesh%NeZ)
    integer, intent(in) :: vmapP(elem3D%NfpTot,lcmesh%NeZ)

    integer :: ke_z, ke2D, ke
    integer :: iM(elem3D%NfpTot), iP(elem3D%NfpTot)
    real(RP) :: DiffFlux_P(elem3D%NfpTot,2)
    !--------------------------------------------------------

    !$omp parallel do private(iM, iP, ke2D, ke_z, ke, DiffFlux_P) collapse(2)
    do ke2D=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      where ( abs(nz(:,ke_z,ke2D)) > 0.5_RP .and. iP(:) == iM(:) )
        DiffFlux_P(:,1) = - DiffFlux(iM,ke2D,1)
        DiffFlux_P(:,2) = - DiffFlux(iM,ke2D,2)
      elsewhere
        DiffFlux_P(:,1) = DiffFlux(iP,ke2D,1)
        DiffFlux_P(:,2) = DiffFlux(iP,ke2D,2)
      endwhere
      bnd_flux(:,ke_z,ke2D,1) = 0.5_RP * ( DiffFlux_P(:,1) - DiffFlux(iM,ke2D,1) ) * nz(:,ke_z,ke2D)
      bnd_flux(:,ke_z,ke2D,2) = 0.5_RP * ( DiffFlux_P(:,2) - DiffFlux(iM,ke2D,2) ) * nz(:,ke_z,ke2D)
    end do
    end do

    return
  end subroutine BLmixing_bnd_flux

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

    real(RP), allocatable :: DENS_tp_tmp(:,:), RHOT_tp_tmp(:,:), RHOH_p_tmp(:,:), RHOQv_tp_tmp(:,:)
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

      allocate( DENS_tp_tmp(elem3D%Np,lcmesh%NeA), RHOT_tp_tmp(elem3D%Np,lcmesh%NeA), RHOQv_tp_tmp(elem3D%Np,lcmesh%NeA)  )
      allocate( RHOH_p_tmp(elem3D%Np,lcmesh%NeA) )

      ! !--
      ! !$omp parallel do
      ! do ke=lcmesh%NeS, lcmesh%NeE
      !   DENS_tp_tmp(:,ke) = 0.0_RP
      !   RHOQv_tp_tmp(:,ke) = 0.0_RP
      !   RHOH_p_tmp(:,ke) = 0.0_RP
      ! end do
      ! call Large_scale_Precip( &
      !   DENS_tp_tmp, RHOQv_tp_tmp, RHOH_p_tmp, &
      !   RAIN_LSC%local(n)%val, CondensRate%local(n)%val,                                  &
      !   QV%val, DDENS%val, DENS_hyd%val, PRES%val, Rtot%val, &
      !   dt, lcmesh, lcmesh%refElem3D )
      ! !$omp parallel do
      ! do ke=lcmesh%NeS, lcmesh%NeE
      !   DDENS%val(:,ke) = DDENS%val(:,ke) + DENS_tp_tmp(:,ke) * dt
      !   QV%val(:,ke) = QV%val(:,ke) + RHOQv_tp_tmp(:,ke) * dt
      !   DRHOT%val(:,ke) = DRHOT%val(:,ke) &
      !     + RHOH_p_tmp(:,ke) / ( CPtot%val(:,ke) * (PRES%va(:,ke)/PRES00)**( Rtot%val(:,ke) / CPtot%val(:,ke)) ) * dt
      ! end do
      ! call atm%vars%Calc_diagnostics()
      

      ! !$omp parallel do
      ! do ke=lcmesh%NeS, lcmesh%NeE
      !   RHOT_tp_tmp(:,ke) = 0.0_RP
      !   RHOQv_tp_tmp(:,ke) = 0.0_RP
      ! end do
      ! call BLmixing( RHOT_tp_tmp, RHOQv_tp_tmp, &
      !   DDENS%val, MOMX%val, MOMY%val, MOMZ%val, QV%val,                            &
      !   DENS_hyd%val, PT%val, PRES%val, atm%mesh%DOptrMat(3), atm%mesh%LiftOptrMat, &
      !   lcmesh, lcmesh%refElem3D )
      ! !$omp parallel do
      ! do ke=lcmesh%NeS, lcmesh%NeE
      !   DDENS%val(:,ke) = DDENS%val(:,ke) + DENS_tp_tmp(:,ke) * dt
      !   QV%val(:,ke) = QV%val(:,ke) + RHOQv_tp_tmp(:,ke) * dt
      !   DRHOT%val(:,ke) = DRHOT%val(:,ke) + RHOT_tp_tmp(:,ke) * dt
      ! end do
      ! call atm%vars%Calc_diagnostics()

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
    
    ! call FILE_HISTORY_meshfield_in( RAIN_LSC, "RAIN with large scale condensation" )
    ! call FILE_HISTORY_meshfield_in( CondensRate, "Condensation rate (mass) with large scale condensation" )

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
