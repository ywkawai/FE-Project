#include "scalelib.h"
module mod_user_sub_LSC

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

  use scale_element_base, only: ElementBase1D, ElementBase3D
  use scale_element_line, only: LineElement
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_element_modalfilter, only: ModalFilter

  use mod_atmos_vars, only: &
    AtmosVars, &
    AtmosVars_GetLocalMeshPrgVars,    &
    AtmosVars_GetLocalMeshPhyAuxVars, &
    AtmosVars_GetLocalMeshQTRC_Qv

  use mod_user_sub_Filter, only: &
    Filter

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  public :: USER_sub_LSC_Init
  public :: USER_sub_LSC_calc_tendency
  
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
  integer :: LSC_nstep
  integer :: LSC_rnstep
  type(MeshField3D) :: LSC_TEND(6)

  type(MeshField2D) :: RAIN_LSC
  type(MeshField3D) :: CondensRate, CondensRate_ori
  type(MeshField3D) :: TEMP, TEMP_LSC

  type(LineElement) :: refElem1D
  type(ModalFilter) :: PhyTendFilter

  logical :: APPLY_ModalFilter

  logical :: APPLY_NewFilter
  type(Filter) :: newFilter

  real(RP), allocatable :: IntrpGL_mat(:,:)
  real(RP), allocatable :: IntPhytendMat(:,:)
  integer :: GL_PolyOrder_h

contains
  subroutine USER_sub_LSC_Init( mesh3D )
    use scale_polynominal, only: &
      Polynominal_GenLagrangePoly
    implicit none
    class(MeshBase3D), intent(in) :: mesh3D

    integer :: LSC_MF_ORDER_h
    real(RP) :: LSC_MF_ALPHA_h
    integer :: LSC_MF_ORDER_v
    real(RP) :: LSC_MF_ALPHA_v

    character(len=H_SHORT) :: FilterShape
    real(RP) :: FilterWidthFac

    namelist / PARAM_USER_LSC / &
       LSC_nstep, &
       APPLY_ModalFilter, &
       LSC_MF_ORDER_h, LSC_MF_ALPHA_h, &
       LSC_MF_ORDER_v, LSC_MF_ALPHA_v, &
       APPLY_NewFilter, FilterShape, FilterWidthFac, &
       GL_PolyOrder_h

    integer :: ierr   

    class(MeshBase2D), pointer :: mesh2D
    type(HexahedralElement) :: refElem3D_dummy

    type(LineElement) :: refElem1D_dummy
    type(QuadrilateralElement) :: refElem2D_dummy
    real(RP), allocatable :: IntrpGL_lag1D(:,:)
    real(RP), allocatable :: IntWGL(:)
    real(RP), allocatable :: x_gl(:), y_gl(:)
    integer :: p, p1, p2, p_, p1_, p2_
  
    integer :: iv
    !------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    LSC_nstep = 1
    APPLY_ModalFilter = .false.
    LSC_MF_ORDER_h = 16; LSC_MF_ALPHA_h = 0.0_RP
    LSC_MF_ORDER_v = 16; LSC_MF_ALPHA_v = 0.0_RP

    APPLY_NewFilter = .false.
    FilterShape     = "GAUSSIAN"
    FilterWidthFac  = 1.0_RP

    GL_PolyOrder_h = mesh3D%refElem3D%PolyOrder_h

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER_LSC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER_LSC. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER_LSC)

    !----
    LSC_rnstep = 0
    do iv=1, 6
      call LSC_TEND(iv)%Init( "LSC_TEND", "", mesh3D )
    end do

    !---
    call mesh3D%GetMesh2D(mesh2D)
    call RAIN_LSC%Init( "RAIN_LSC", "kg/s", mesh2D )
    call CondensRate%Init( "CondensRate_LSC", "kg/m3.s-1", mesh3D )
    call CondensRate_ori%Init( "CondensRate_LSC_ori", "kg/m3.s-1",mesh3D )

    call TEMP%Init( "TEMP", "K", mesh3D )
    call TEMP_LSC%Init( "TEMP_LSC", "K", mesh3D )

    !--
    call refElem1D%Init( mesh3D%refElem3D%PolyOrder_v, .false. )

    call refElem3D_dummy%Init( mesh3D%refElem3D%PolyOrder_h, mesh3D%refElem3D%PolyOrder_v, .false. )
    call PhyTendFilter%Init( refElem3D_dummy, 0.0_RP, LSC_MF_ALPHA_h, LSC_MF_ORDER_h, 0.0_RP, LSC_MF_ALPHA_v, LSC_MF_ORDER_v)
    call refElem3D_dummy%Final()

    !--
    call newFilter%Init( FilterShape, FilterWidthFac, mesh3D )

    !--
    ! call refElem1D_dummy%Init( mesh3D%refElem3D%PolyOrder_h, .false. )
    ! call refElem2D_dummy%Init( mesh3D%refElem3D%PolyOrder_h, .false. )
    
    ! allocate( IntrpGL_mat(GL_PolyOrder_h**2,refElem2D_dummy%Np) )
    ! allocate( IntPhytendMat(refElem2D_dummy%Np,GL_PolyOrder_h**2) )

    ! allocate( x_gl(GL_PolyOrder_h**2), y_gl(GL_PolyOrder_h**2) )
    ! allocate( IntWGL(GL_PolyOrder_h**2) )
    ! IntrpGL_mat(:,:) = refElem2D_dummy%GenIntGaussLegendreIntrpMat( GL_PolyOrder_h, &
    !                     IntWGL(:), x_gl(:), y_gl(:) )
    
    ! allocate( IntrpGL_lag1D(GL_PolyOrder_h,refElem1D_dummy%PolyOrder+1) )
    ! IntrpGL_lag1D(:,:) = Polynominal_GenLagrangePoly( refElem1D_dummy%PolyOrder, refElem1D_dummy%x1, x_gl(1:GL_PolyOrder_h) )
    ! do p2_=1, refElem1D_dummy%Np
    ! do p1_=1, refElem1D_dummy%Np
    !   p_ = p1_ + (p2_-1)*refElem1D_dummy%Np
    !   do p2=1, GL_PolyOrder_h
    !   do p1=1, GL_PolyOrder_h
    !     p = p1 + (p2-1)*refElem1D_dummy%Np
    !     IntPhytendMat(p_,p) = IntWGL(p) * IntrpGL_lag1D(p,p_)
    !   end do
    !   end do
    ! end do
    ! end do
    ! IntPhytendMat(:,:) = matmul(refElem2D_dummy%invM, IntPhytendMat)
    ! call refElem1D_dummy%Final()
    ! call refElem2D_dummy%Final()
    
    return
  end subroutine USER_sub_LSC_Init

!OCL SERIAL
  subroutine USER_sub_LSC_calc_tendency( vars, mesh3D, dt )
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
    class(AtmosVars), intent(inout) :: vars
    class(MeshBase3D), intent(in), target :: mesh3D
    real(RP), intent(in) :: dt

    class(LocalMesh3D), pointer :: lcmesh
    integer :: n
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT, QV
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: RHOQv_tp

    integer :: ke

    type(MeshField3D) :: tmp_field
    integer :: ke2D
    !----------------------------------------------------------

    if ( LSC_rnstep == 0 ) then
      call tmp_field%Init( "tmp_field", "K", mesh3D )

      do n=1, mesh3D%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshPrgVars( n, mesh3D,  &
          vars%PROGVARS_manager, vars%AUXVARS_manager,     &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh   )      

        call AtmosVars_GetLocalMeshPhyAuxVars( n, mesh3D, &
          vars%AUXVARS_manager, PRES, PT                  )
        !$omp parallel do private(ke2D)
        do ke=lcmesh%NeS, lcmesh%NeE
          ke2D = lcmesh%EMap3Dto2D(ke)
          TEMP%local(n)%val(:,ke) = PRES%val(:,ke) / ( Rtot%val(:,ke) * ( DENS_hyd%val(:,ke) + DDENS%val(:,ke) ) )
          tmp_field%local(n)%val(:,ke) = lcmesh%lat2D(lcmesh%refElem3D%IndexH2Dto3D(:),ke2D)
        end do
      end do
      ! call newFilter%Apply( TEMP_LSC, &
      !   tmp_field, mesh3D          )

      do n=1, mesh3D%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshPrgVars( n, mesh3D,  &
          vars%PROGVARS_manager, vars%AUXVARS_manager,     &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh   )      

        call AtmosVars_GetLocalMeshPhyAuxVars( n, mesh3D, &
          vars%AUXVARS_manager, PRES, PT                  )
        call AtmosVars_GetLocalMeshQTRC_Qv( n, mesh3D, &
          vars%QTRCVARS_manager, vars%PHYTENDS_manager, QV, RHOQv_tp )

        !--
        call Large_scale_Precip_core1( CondensRate_ori%local(n)%val, &
          QV%val, DDENS%val, DENS_hyd%val, PRES%val,                 &
!          TEMP_LSC%local(n)%val, &
          TEMP%local(n)%val, &
          real(LSC_nstep, kind=RP) * dt, lcmesh, lcmesh%refElem3D%Np )
      end do

      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh => mesh3D%lcmesh_list(n)
        !$omp parallel do
        do ke=lcmesh%NeS, lcmesh%NeE
          CondensRate%local(n)%val(:,ke) = CondensRate_ori%local(n)%val(:,ke)
        end do
      end do
      
      if ( APPLY_NewFilter ) call newFilter%Apply( CondensRate, mesh3D )

      do n=1, mesh3D%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshPrgVars( n, mesh3D,  &
          vars%PROGVARS_manager, vars%AUXVARS_manager,     &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh   )      
        call AtmosVars_GetLocalMeshQTRC_Qv( n, mesh3D, &
          vars%QTRCVARS_manager, vars%PHYTENDS_manager, QV, RHOQv_tp )

        !--
        call Large_scale_Precip_core2( &
          LSC_TEND(1)%local(n)%val, LSC_TEND(2)%local(n)%val, LSC_TEND(3)%local(n)%val, LSC_TEND(4)%local(n)%val,  &
          LSC_TEND(5)%local(n)%val, LSC_TEND(6)%local(n)%val, RAIN_LSC%local(n)%val,                               &
!          CondensRate%local(n)%val,                                                                                &
          CondensRate%local(n)%val,                                                                                &
          QV%val, DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DENS_hyd%val, PRES%val, Rtot%val, CVtot%val, CPtot%val, &
          TEMP%local(n)%val,                                                                                       &
          real(LSC_nstep, kind=RP) * dt, lcmesh, lcmesh%refElem3D, refElem1D )
      end do

      LSC_rnstep = LSC_nstep
    end if

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      call AtmosVars_GetLocalMeshQTRC_Qv( n, mesh3D, &
        vars%QTRCVARS_manager, vars%PHYTENDS_manager, QV, RHOQv_tp )

      !$omp parallel do
      do ke=lcmesh%NeS, lcmesh%NeE
        vars%PHY_TEND(DENS_p)%local(n)%val(:,ke) = vars%PHY_TEND(DENS_p)%local(n)%val(:,ke) + LSC_TEND(1)%local(n)%val(:,ke)
        vars%PHY_TEND(MOMX_p)%local(n)%val(:,ke) = vars%PHY_TEND(MOMX_p)%local(n)%val(:,ke) + LSC_TEND(2)%local(n)%val(:,ke)
        vars%PHY_TEND(MOMY_p)%local(n)%val(:,ke) = vars%PHY_TEND(MOMY_p)%local(n)%val(:,ke) + LSC_TEND(3)%local(n)%val(:,ke)
        vars%PHY_TEND(MOMZ_p)%local(n)%val(:,ke) = vars%PHY_TEND(MOMZ_p)%local(n)%val(:,ke) + LSC_TEND(4)%local(n)%val(:,ke)
        RHOQv_tp%val(:,ke) = RHOQv_tp%val(:,ke) + LSC_TEND(5)%local(n)%val(:,ke)
        vars%PHY_TEND(RHOH_p)%local(n)%val(:,ke) = vars%PHY_TEND(RHOH_p)%local(n)%val(:,ke) + LSC_TEND(6)%local(n)%val(:,ke)
      end do
      LSC_rnstep = LSC_rnstep - 1
    end do

    call FILE_HISTORY_meshfield_in( RAIN_LSC, "RAIN with large scale condensation" )
    call FILE_HISTORY_meshfield_in( CondensRate, "Condensation rate (mass) with large scale condensation" )
    call FILE_HISTORY_meshfield_in( CondensRate_ori, "Condensation rate (mass) with large scale condensation" )

    call FILE_HISTORY_meshfield_in( TEMP_LSC, "filtered TEMP" )
    call FILE_HISTORY_meshfield_in( TEMP, "Unfiltered TEMP" )

    return
  end subroutine USER_sub_LSC_calc_tendency

!OCL SERIAL
  subroutine Large_scale_Precip_core1( CondensRate_, &
    QV, DDENS, DENS_hyd, PRES, TEMP_, &
    DT, lcmesh, Np )
    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_pres2qsat_liq
    use scale_atmos_hydrometeor, only: &
      CV_WATER      
    use scale_const, only: &
      CL => CONST_CL    
    implicit none
    type(LocalMesh3D), intent(in) :: lcmesh
    integer, intent(in) :: Np
    real(RP), intent(out) :: CondensRate_(Np,lcmesh%NeA)
    real(RP), intent(in) :: QV(Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(Np,lcmesh%NeA)
    real(RP), intent(in) :: TEMP_(Np,lcmesh%NeA)
    real(RP), intent(in) :: DT

    integer :: ke, ke_xy, ke_z
    integer :: p2D
    real(RP) :: Qsat(Np,lcmesh%NeA)
    real(RP) :: DENS(Np)
    real(RP) :: coef1, coef2
    !----------------------------------------------

    !$omp parallel do private(ke)
    do ke=lcmesh%NeS, lcmesh%NeE      
      call ATMOS_SATURATION_pres2qsat_liq( Np, 1, Np, &
        TEMP_(:,ke), PRES(:,ke), &
        Qsat(:,ke) )
    end do

    coef1 = LHV0**2 / ( CVDry * Rvap )
    coef2 = Rvap / LHV0

    !$omp parallel do private(ke, DENS) collapse(2)
    do ke_z=1, lcmesh%NeZ
    do ke_xy=1, lcmesh%Ne2D
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D
      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)

      CondensRate_(:,ke) = DENS(:) &
        * max( 0.0_RP, ( QV(:,ke) - Qsat(:,ke) ) / ( 1.0_RP + coef1 / TEMP_(:,ke)**2 * ( 1.0_RP - coef2 * TEMP_(:,ke) ) * Qsat(:,ke) ) / DT )
    end do
    end do

    return
  end subroutine Large_scale_Precip_core1

!OCL SERIAL
  subroutine Large_scale_Precip_core2( DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOQv_tp, RHOH_p, SFLX_RAIN, &
    CondensRate_, QV, DDENS, MOMX, MOMY, MOMZ, DENS_hyd, PRES, Rtot, CVtot, CPtot, TEMP_, &
    DT, lcmesh, elem3D, elem1D )
    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_pres2qsat_liq
    use scale_atmos_hydrometeor, only: &
      CV_WATER      
    use scale_const, only: &
      CL => CONST_CL    
    implicit none
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    class(ElementBase1D), intent(in) :: elem1D
    real(RP), intent(inout) :: DENS_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMZ_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOQv_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOH_p(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: SFLX_RAIN(elem3D%Nnode_h1D**2,lcmesh%lcmesh2D%NeA)
    real(RP), intent(in) :: CondensRate_(elem3D%Np,lcmesh%NeA)
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
    real(RP), intent(in) :: TEMP_(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DT

    integer :: ke, ke_xy, ke_z
    integer :: p2D
    real(RP) :: Qsat(elem3D%Np,lcmesh%NeA)
    real(RP) :: CondensRate_zxy(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D) !< Condensation rate
    real(RP) :: DENS(elem3D%Np)
    real(RP) :: rain_tmp(elem3D%Nnode_v)

    real(RP) :: coef1
    !----------------------------------------------

    !$omp parallel do private(ke, DENS) collapse(2)
    do ke_z=1, lcmesh%NeZ
    do ke_xy=1, lcmesh%Ne2D
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D
      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      
      CondensRate_zxy(:,ke_z,ke_xy) = CondensRate_(:,ke)

      RHOH_p  (:,ke) = &
        LHV0 * CondensRate_zxy(:,ke_z,ke_xy) - CPtot(:,ke) * TEMP_(:,ke) * CondensRate_zxy(:,ke_z,ke_xy)

      DENS_tp(:,ke) = - CondensRate_zxy(:,ke_z,ke_xy)
      RHOQv_tp(:,ke) = - CondensRate_zxy(:,ke_z,ke_xy)
      MOMX_tp(:,ke) = - CondensRate_zxy(:,ke_z,ke_xy) * MOMX(:,ke) / DENS(:)
      MOMY_tp(:,ke) = - CondensRate_zxy(:,ke_z,ke_xy) * MOMY(:,ke) / DENS(:)
      MOMZ_tp(:,ke) = - CondensRate_zxy(:,ke_z,ke_xy) * MOMZ(:,ke) / DENS(:)
    end do
    end do

    !$omp parallel private(ke_z,ke,p2D, rain_tmp, coef1)
    !$omp workshare
    SFLX_RAIN(:,:) = 0.0_RP
    !$omp end workshare
    !$omp do 
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
    !$omp end parallel

    return
  end subroutine Large_scale_Precip_core2
end module mod_user_sub_LSC
