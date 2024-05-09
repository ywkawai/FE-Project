!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Global nonhydrostatic model / HEVE
!!
!! @par Description
!!      HEVE DGM scheme for Global Atmospheric Dynamical process. 
!!      The governing equations is a fully compressibile nonhydrostic equations, 
!!      which consist of mass, momentum, and thermodynamics (density * potential temperature conservation) equations. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_globalnonhydro3d_rhot_heve
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,    &
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  &
    CVdry => CONST_CVdry,  &
    PRES00 => CONST_PRE00, &
    RPlanet => CONST_RADIUS

  use scale_sparsemat  
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    atm_dyn_dgm_nonhydro3d_common_Init,                       &
    atm_dyn_dgm_nonhydro3d_common_Final,                      &
    DENS_VID => PRGVAR_DDENS_ID, RHOT_VID => PRGVAR_DRHOT_ID, &
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID,                               &
    PRGVAR_NUM, IntrpMat_VPOrdM1
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_Init
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_Final
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_shallow_atm
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_deep_atm
#ifdef SCALE_PRODUCT_RUN_GM_MOUNTAIN_WAVE
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_set_dampcoef
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

#ifdef SCALE_PRODUCT_RUN_GM_MOUNTAIN_WAVE
  type(MeshField3D), public :: forcing_U0
  type(MeshField3D), public :: forcing_V0
  type(MeshField3D), public :: forcing_W0
  real(RP) :: U0
  logical :: ini_bg_force_flag
  real(RP) :: ini_bg_force_tscale
  real(RP) :: ini_bg_force_turnoff_tstart
  real(RP) :: ini_bg_force_turnoff_tscale
  real(RP) :: ini_bg_force_tscale2
  real(RP) :: ini_bg_sfac

  real(RP), allocatable :: sfac(:,:)
  real(RP) :: sw
#endif

contains
!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_Init( mesh )
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSVec    
    use scale_const, only: &
      PI => CONST_PI
    implicit none
    class(MeshBase3D), target, intent(in) :: mesh

#ifdef SCALE_PRODUCT_RUN_GM_MOUNTAIN_WAVE
    integer :: n, ke, ke2d
    class(ElementBase3D), pointer :: elem
    class(LocalMesh3D), pointer :: lmesh3D
    real(RP), allocatable :: lon(:), lat(:)
    real(RP), allocatable :: Umet(:,:), Vmet(:,:)
    real(RP), allocatable :: GsqrtV(:)

    real(RP) :: zTop
    real(RP) :: SPONGE_HEIGHT
    real(RP) :: SPONGE_LATERAL_WIDTH
    real(RP) :: SPONGE_EFOLD_SEC
    logical :: lateral_sponge_layer_flag
    real(RP) :: LATERAL_SPONGE_EFOLD_SEC
    real(RP) :: SL_TANH_NONDIM_WIDTH 

    real(RP) :: rtau_sponge
    real(RP) :: rtau_lateral_sponge

    logical  :: SL_MERI_TAPER_FLAG
    real(RP)  :: SL_MERI_TAPER_TANH_Clat
    real(RP)  :: SL_MERI_TAPER_TANH_LatWidth

    namelist / PARAM_USER_MTWAVE / &
      U0,                  &
      zTop,                &
      SPONGE_HEIGHT,       &
      SPONGE_EFOLD_SEC,    &
      lateral_sponge_layer_flag, &
      LATERAL_SPONGE_EFOLD_SEC,  &
      SL_TANH_NONDIM_WIDTH,      &
      ini_bg_force_flag,         &
      ini_bg_force_tscale,       &
      ini_bg_force_turnoff_tstart, &       
      ini_bg_force_turnoff_tscale, &
      SL_MERI_TAPER_FLAG, &
      SL_MERI_TAPER_TANH_Clat,    &
      SL_MERI_TAPER_TANH_LatWidth    
    integer :: ierr
    
    real(RP), allocatable :: sfac_h(:)
#endif
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Init( mesh )
#ifdef SCALE_PRODUCT_RUN_GM_MOUNTAIN_WAVE
    !---
    U0 = 20.0_RP
    ini_bg_force_flag = .true.
    ini_bg_force_tscale = 60.0_RP
    ini_bg_force_turnoff_tstart = 120.0_RP
    ini_bg_force_turnoff_tscale = 1800.0_RP  
    ! ini_bg_force_tscale = 600.0_RP
    ! ini_bg_force_turnoff_tstart = 3000.0_RP!120.0_RP
    ! ini_bg_force_turnoff_tscale = 6000.0_RP  
    LOG_INFO("SCALE_PRODUCT_RUN_GM_MOUNTAIN_WAVE",*) ini_bg_force_tscale, ini_bg_force_turnoff_tstart, ini_bg_force_turnoff_tscale
    
    zTop = 30E3_RP
    SPONGE_HEIGHT = 15E3_RP
    SPONGE_EFOLD_SEC = 100E0_RP
    SPONGE_LATERAL_WIDTH = 120E3_RP
    lateral_sponge_layer_flag = .true. 
    LATERAL_SPONGE_EFOLD_SEC = 200E0_RP!100E0_RP
    SL_TANH_NONDIM_WIDTH = 0.16E0_RP

    SL_MERI_TAPER_FLAG = .true.
    SL_MERI_TAPER_TANH_Clat = 1.0471975511965976D0      ! 60 deg
    SL_MERI_TAPER_TANH_LatWidth = 0.13962634015954636D0 ! 8 deg          

    LOG_NEWLINE
    LOG_INFO("globalnonhydro3d_rhot_heve_Init",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER_MTWAVE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("globalnonhydro3d_rhot_heve_Init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("globalnonhydro3d_rhot_heve_Init",*) 'Not appropriate names in namelist PARAM_USER_MTWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER_MTWAVE)
    !-----------------------

    call forcing_U0%Init( "forcing_U0", "m/s", mesh )
    call forcing_V0%Init( "forcing_V0", "m/s", mesh )
    call forcing_W0%Init( "forcing_W0", "m/s", mesh )

    lmesh3D => mesh%lcmesh_list(1)
    elem => lmesh3D%refElem3D
    allocate( lon(elem%Np), lat(elem%Np) )
    allocate( Umet(elem%Np,lmesh3D%Ne), Vmet(elem%Np,lmesh3D%Ne) )
    allocate( GsqrtV(elem%Np) )
    allocate( sfac(elem%Np,lmesh3D%Ne) )
    allocate( sfac_h(elem%Np) )

    do n=1, mesh%LOCAL_MESH_NUM
      lmesh3D => mesh%lcmesh_list(n)
      !$omp parallel do private(ke2D,lon,lat,ke)
      do ke=lmesh3D%NeS, lmesh3D%NeE
        ke2D = lmesh3D%EMap3Dto2D(ke)
        lon(:) = lmesh3D%lon2D(elem%IndexH2Dto3D,ke2D)
        lat(:) = lmesh3D%lat2D(elem%IndexH2Dto3D,ke2D)

        Umet(:,ke) = U0 * cos(lat(:))
        Vmet(:,ke) = 0.0_RP
      end do
      call CubedSphereCoordCnv_LonLat2CSVec( &
        lmesh3D%panelID, lmesh3D%pos_en(:,:,1), lmesh3D%pos_en(:,:,2),     & ! (in)
        lmesh3D%gam(:,lmesh3D%NeS:lmesh3D%NeE), elem%Np * lmesh3D%Ne,      & ! (in)
        Umet(:,:), Vmet(:,:),                                              & ! (in)
        forcing_U0%local(n)%val(:,lmesh3D%NeS:lmesh3D%NeE),                & ! (out)
        forcing_V0%local(n)%val(:,lmesh3D%NeS:lmesh3D%NeE)                 ) ! (out)
      
      !$omp parallel do private( GsqrtV, ke2d )
      do ke=lmesh3D%NeS, lmesh3D%NeE
        ke2D = lmesh3D%EMap3Dto2D(ke)
        GsqrtV(:) = lmesh3D%Gsqrt(:,ke) / lmesh3D%GsqrtH(elem%IndexH2Dto3D,ke2d)        
        forcing_W0%local(n)%val(:,ke) = &
          - GsqrtV(:) * ( lmesh3D%GI3(:,ke,1) * forcing_U0%local(n)%val(:,ke) + lmesh3D%GI3(:,ke,2) * forcing_V0%local(n)%val(:,ke) ) &
            * exp(-lmesh3D%pos_en(:,ke,3)/2000.0_RP) * 0.0_RP
      end do
    end do

    rtau_sponge = 1.0_RP / SPONGE_EFOLD_SEC
    if ( lateral_sponge_layer_flag ) then
      rtau_lateral_sponge = 1.0_RP / LATERAL_SPONGE_EFOLD_SEC
    else
      rtau_lateral_sponge = 0.0_RP
    end if

    !$omp parallel do private(lon, lat, sfac_h, ke2D)
    do ke=lmesh3D%NeS, lmesh3D%NeE
      ke2D = lmesh3D%EMap3Dto2D(ke)
      lon(:) = lmesh3D%lon2D(elem%IndexH2Dto3D,ke2D)
      lat(:) = lmesh3D%lat2D(elem%IndexH2Dto3D,ke2D)

      sfac_h(:) = &
        rtau_lateral_sponge * 0.5_RP * ( 1.0_RP - tanh( ( lon(:) - PI * 0.25_RP ) / ( SL_TANH_NONDIM_WIDTH * PI * 0.5_RP ) ) ) &
      + rtau_lateral_sponge * 0.5_RP * ( 1.0_RP + tanh( ( lon(:) - PI * 1.75_RP ) / ( SL_TANH_NONDIM_WIDTH * PI * 0.5_RP ) ) )

      if ( SL_MERI_TAPER_FLAG ) then
        sfac_h(:) = sfac_h(:) &
          *  0.5_RP * ( 1.0_RP - tanh( ( abs(lat(:)) - SL_MERI_TAPER_TANH_Clat ) / SL_MERI_TAPER_TANH_LatWidth ) )
      end if
      sfac(:,ke) = &
          rtau_sponge * 0.5_RP * ( 1.0_RP + tanh( ( lmesh3D%zlev(:,ke) - 0.5_RP * ( zTop + SPONGE_HEIGHT ) ) / ( SL_TANH_NONDIM_WIDTH * ( zTop - SPONGE_HEIGHT ) ) ) ) &
        + sfac_h(:)
    end do    
#endif
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_Init

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_Final()
    implicit none
    !--------------------------------------------
    
    call atm_dyn_dgm_nonhydro3d_common_Final()    
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_Final  

#ifdef SCALE_PRODUCT_RUN_GM_MOUNTAIN_WAVE
!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_set_dampcoef( tsec )
    use scale_const, only: PI => CONST_PI
    implicit none
    real(RP), intent(in) :: tsec

    real(RP) :: ini_bg_off_tsec
    !-----------------------------------

    if ( ini_bg_force_flag ) then
      ini_bg_off_tsec = ini_bg_force_turnoff_tstart + ini_bg_force_turnoff_tscale
      if ( tsec < ini_bg_off_tsec ) then
        if ( tsec > ini_bg_force_turnoff_tstart ) then
          sw = 0.5_RP * ( 1.0_RP - cos( PI * ( ( tsec - ini_bg_force_turnoff_tstart ) / ini_bg_force_turnoff_tscale - 1.0_RP ) ) )
        else
          sw = 1.0_RP
        end if
      else
        sw = 0.0_RP
      end if
      ini_bg_sfac = sw * 1.0_RP / ini_bg_force_tscale
    else
      sw = 0.0_RP
      ini_bg_sfac = 0.0_RP
    end if

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_set_dampcoef
#endif

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_shallow_atm( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                   & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, PRES_hyd_ref, & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot,                                                  & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )                     ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc
    use scale_const, only: &
      OHM => CONST_OHM
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift
    real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DPRES_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd_ref(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot (elem%Np,lmesh%NeA)    
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DPRES_hyd(elem%Np), GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: RHOT_(elem%Np)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%np), drho(elem%Np)

    real(RP) :: G11(elem%Np), G12(elem%Np), G22(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np)
    real(RP) :: X2D(elem2D%Np,lmesh2D%Ne), Y2D(elem2D%Np,lmesh2D%Ne)
    real(RP) :: X(elem%Np), Y(elem%Np), twoOVdel2(elem%Np)
    real(RP) :: CORI(elem%Np,2)
    logical :: is_panel1to4
    real(RP) :: s

    integer :: ke, ke2d
    integer :: p, p12, p3

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR

    real(RP) :: damp_coef(elem%Np)
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux, del_flux_hyd,                                                  & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd,         & ! (in)
      Rtot, CVtot, CPtot,                                                      & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%GsqrtH, lmesh%gam, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),             & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%vmapM, lmesh%vmapP, elem%IndexH2Dto3D_bnd,                         & ! (in)
      lmesh, elem, lmesh2D, elem2D                                             ) ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    s = 1.0_RP 
    is_panel1to4 = .true.
    if ( lmesh%panelID == 5 ) then
      is_panel1to4 = .false.
    else if ( lmesh%panelID == 6 ) then
      is_panel1to4 = .false.
      s = - 1.0_RP
    end if

    !$omp parallel private(                        &
    !$omp RHOT_, rdens_, u_, v_, w_, wt_,          &
    !$omp Fx, Fy, Fz, LiftDelFlx,                  &
    !$omp drho, DPRES_hyd, GradPhyd_x, GradPhyd_y, &
    !$omp G11, G12, G22, GsqrtV, RGsqrtV,          &
    !$omp X, Y, twoOVdel2,                         &
    !$omp CORI, ke, ke2D, damp_coef                )

    !$omp do
    do ke2D = lmesh2D%NeS, lmesh2D%NeE
      X2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,1))
      Y2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,2))
    end do

    !$omp do
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      G11(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,1)
      G12(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,2)
      G22(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,2,2)
      GsqrtV(:)  = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2d)
      RGsqrtV(:) = 1.0_RP / GsqrtV(:)

      !--
      RHOT_(:) = P0ovR * ( PRES_hyd(:,ke) * rP0 )**rgamm + DRHOT_(:,ke)
      ! DPRES_(:) = PRES00 * ( Rtot(:,ke) * rP0 * RHOT_(:) )**( CPtot(:,ke) / CVtot(:,ke) ) &
      !           - PRES_hyd(:,ke)

      rdens_(:) = 1.0_RP / ( DDENS_(:,ke) + DENS_hyd(:,ke) )
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)
      wt_(:) = w_(:) * RGsqrtV(:) + lmesh%GI3(:,ke,1) * u_(:) + lmesh%GI3(:,ke,2) * v_(:) 

      X(:) = X2D(elem%IndexH2Dto3D,ke2d)
      Y(:) = Y2D(elem%IndexH2Dto3D,ke2d)
      twoOVdel2(:) = 2.0_RP / ( 1.0_RP + X(:)**2 + Y(:)**2 )

      CORI(:,1) = s * OHM * twoOVdel2(:) * ( - X(:) * Y(:)          * MOMX_(:,ke) + ( 1.0_RP + Y(:)**2 ) * MOMY_(:,ke) )
      CORI(:,2) = s * OHM * twoOVdel2(:) * ( - ( 1.0_RP + X(:)**2 ) * MOMX_(:,ke) +  X(:) * Y(:)         * MOMY_(:,ke) )
      if ( is_panel1to4 ) then
        CORI(:,1) = s * Y(:) * CORI(:,1)
        CORI(:,2) = s * Y(:) * CORI(:,2)
      end if

      drho(:) = matmul(IntrpMat_VPOrdM1, DDENS_(:,ke))

      !-- Gradient hydrostatic pressure
      
      DPRES_hyd(:) = PRES_hyd(:,ke) - PRES_hyd_ref(:,ke)

      call sparsemat_matmul(Dx, GsqrtV(:) * DPRES_hyd(:), Fx)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,1) * DPRES_hyd(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,1), LiftDelFlx)
      GradPhyd_x(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      call sparsemat_matmul(Dy, GsqrtV(:) * DPRES_hyd(:), Fy)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,2) * DPRES_hyd(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,2), LiftDelFlx)
      GradPhyd_y(:) = lmesh%Escale(:,ke,2,2) * Fy(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)
      
      !-- DENS
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * MOMY_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( DDENS_(:,ke) + DENS_hyd(:,ke) ) * wt_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DENS_VID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)     &
          + lmesh%Escale(:,ke,2,2) * Fy(:)     &
          + lmesh%Escale(:,ke,3,3) * Fz(:)     &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      !-- MOMX
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMX_(:,ke) + G11(:) * DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMX_(:,ke) + G12(:) * DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMX_(:,ke)                              &
                                                    + ( lmesh%GI3(:,ke,1) * G11(:) + lmesh%GI3(:,ke,2) * G12(:) ) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                         &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                         &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                         &
            + LiftDelFlx(:)                   ) / lmesh%Gsqrt(:,ke)                  &
          - twoOVdel2(:) * Y(:) *                                                    &
            ( X(:) * Y(:) * u_(:) - (1.0_RP + Y(:)**2) * v_(:) ) * MOMX_(:,ke)       &
          - ( G11(:) * GradPhyd_x(:) + G12(:) * GradPhyd_y(:) ) * RGsqrtV(:)         &
          + CORI(:,1)              

      !-- MOMY
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMY_(:,ke) + G12(:) * DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMY_(:,ke) + G22(:) * DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMY_(:,ke)                              &
                                                    + ( lmesh%GI3(:,ke,1) * G12(:) + lmesh%GI3(:,ke,2) * G22(:) ) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = &
            - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                         &
              + lmesh%Escale(:,ke,2,2) * Fy(:)                                         &
              + lmesh%Escale(:,ke,3,3) * Fz(:)                                         &
              + LiftDelFlx(:)                  ) / lmesh%Gsqrt(:,ke)                   &
            - twoOVdel2(:) * X(:) *                                                    &
              ( - (1.0_RP + X(:)**2) * u_(:) + X(:) * Y(:) * v_(:) ) * MOMY_(:,ke)     &
            - ( G12(:) * GradPhyd_x(:) + G22(:) * GradPhyd_y(:) ) * RGsqrtV(:)         &
            + CORI(:,2) 

      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *   u_ (:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *   v_ (:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMZ_(:,ke) + RGsqrtV(:) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)       &
            + lmesh%Escale(:,ke,2,2) * Fy(:)       &
            + lmesh%Escale(:,ke,3,3) * Fz(:)       &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)  &
          - Grav * drho(:)

      !-- RHOT
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * u_ (:) * RHOT_(:), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * v_ (:) * RHOT_(:), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * wt_(:) * RHOT_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,RHOT_VID), LiftDelFlx)
      
      RHOT_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 

#ifdef SCALE_PRODUCT_RUN_GM_MOUNTAIN_WAVE
      damp_coef(:) = ini_bg_sfac + ( 1.0_RP - sw ) * sfac(:,ke )
      DENS_dt(:,ke) = DENS_dt(:,ke) &
        - damp_coef(:) * DDENS_(:,ke)
      MOMX_dt(:,ke) = MOMX_dt(:,ke) &
        - damp_coef(:) * ( MOMX_(:,ke) - ( DENS_hyd(:,ke) + DDENS_(:,ke) ) * forcing_U0%local(lmesh%lcdomID)%val(:,ke) )! &
      !   + ( 1.0_RP + 0.0_RP*DDENS_(:,ke)/DENS_hyd(:,ke) ) * ( G11(:) * GradPhyd_x(:) + G12(:) * GradPhyd_y(:) ) * RGsqrtV(:)     
      MOMY_dt(:,ke) = MOMY_dt(:,ke) &
        - damp_coef(:) * ( MOMY_(:,ke) - ( DENS_hyd(:,ke) + DDENS_(:,ke) ) * forcing_V0%local(lmesh%lcdomID)%val(:,ke) ) !&
      !   + ( 1.0_RP + 0.0_RP*DDENS_(:,ke)/DENS_hyd(:,ke) ) * ( G12(:) * GradPhyd_x(:) + G22(:) * GradPhyd_y(:) ) * RGsqrtV(:)    
      MOMZ_dt(:,ke) = MOMZ_dt(:,ke) &
        - damp_coef(:) * ( MOMZ_(:,ke) - ( DENS_hyd(:,ke) + DDENS_(:,ke) ) * forcing_W0%local(lmesh%lcdomID)%val(:,ke) )
      RHOT_dt(:,ke) = RHOT_dt(:,ke) &
        - damp_coef(:) * DRHOT_(:,ke)
#endif    
    end do
    !$omp end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_shallow_atm


!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_deep_atm( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                   & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, PRES_hyd_ref, & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot,                                                  & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )                     ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc
    use scale_const, only: &
      OHM => CONST_OHM
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift
    real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DPRES_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd_ref(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot (elem%Np,lmesh%NeA)    
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DPRES_hyd(elem%Np), GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: RHOT_(elem%Np)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%np), drho(elem%Np)

    real(RP) :: G11(elem%Np), G12(elem%Np), G22(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np), Rgam2(elem%Np)
    real(RP) :: X2D(elem2D%Np,lmesh2D%Ne), Y2D(elem2D%Np,lmesh2D%Ne)
    real(RP) :: X(elem%Np), Y(elem%Np), twoOVdel2(elem%Np)
    real(RP) :: OM1(elem%Np), OM2(elem%Np), OM3(elem%Np), DEL(elem%Np), R(elem%Np)
    logical :: is_panel1to4
    real(RP) :: s

    integer :: ke, ke2d
    integer :: p
    
    real(RP) :: rgamm    
    real(RP) :: rP0
    real(RP) :: P0ovR
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux, del_flux_hyd,                                                  & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd,         & ! (in)
      Rtot, CVtot, CPtot,                                                      & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%GsqrtH, lmesh%gam, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),             & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%vmapM, lmesh%vmapP, elem%IndexH2Dto3D_bnd,                         & ! (in)
      lmesh, elem, lmesh2D, elem2D                                             ) ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    P0ovR = PRES00 / Rdry

    s = 2.0_RP * OHM
    is_panel1to4 = .true.
    if ( lmesh%panelID == 5 ) then
      is_panel1to4 = .false.
    else if ( lmesh%panelID == 6 ) then
      is_panel1to4 = .false.
      s = - s
    end if

    !$omp parallel private(                           &
    !$omp RHOT_, rdens_, u_, v_, w_, wt_,             &
    !$omp Fx, Fy, Fz, LiftDelFlx,                     &
    !$omp drho, DPRES_hyd, GradPhyd_x, GradPhyd_y,    &
    !$omp G11, G12, G22, Rgam2, GsqrtV, RGsqrtV,      &
    !$omp X, Y, twoOVdel2,                            &
    !$omp OM1, OM2, OM3, DEL, R, ke, ke2D             )

    !$omp do
    do ke2D = lmesh2D%NeS, lmesh2D%NeE
      X2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,1))
      Y2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,2))
    end do

    !$omp do
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      Rgam2(:) = 1.0_RP / lmesh%gam(:,ke)**2
      G11(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,1) * Rgam2(:)
      G12(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,2) * Rgam2(:)
      G22(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,2,2) * Rgam2(:)
      GsqrtV(:)  = lmesh%Gsqrt(:,ke) * Rgam2(:) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2d)
      RGsqrtV(:) = 1.0_RP / GsqrtV(:)

      !--
      RHOT_(:) = P0ovR * ( PRES_hyd(:,ke) * rP0 )**rgamm + DRHOT_(:,ke)
      ! DPRES_(:) = PRES00 * ( Rtot(:,ke) * rP0 * RHOT_(:) )**( CPtot(:,ke) / CVtot(:,ke) ) &
      !           - PRES_hyd(:,ke)

      rdens_(:) = 1.0_RP / ( DDENS_(:,ke) + DENS_hyd(:,ke) )
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)
      wt_(:) = w_(:) * RGsqrtV(:) + lmesh%GI3(:,ke,1) * u_(:) + lmesh%GI3(:,ke,2) * v_(:) 

      X(:) = X2D(elem%IndexH2Dto3D,ke2d)
      Y(:) = Y2D(elem%IndexH2Dto3D,ke2d)
      DEL(:) = sqrt( 1.0_RP + X(:)**2 + Y(:)**2 )
      twoOVdel2(:) = 2.0_RP / ( 1.0_RP + X(:)**2 + Y(:)**2 )

      R(:) = RPlanet * lmesh%gam(:,ke)

      !-  pnl=1~4: OM1:                     0, OM2: s del / (r (1+Y^2)) ,   s OM3 : Y /del
      !-  pnl=5,6: OM1: - s X del/(r (1+X^2)), OM2: -  s Y del/(r(1+Y^2)), OM3 : s/del
      if ( is_panel1to4 ) then
        OM1(:) = 0.0_RP
        OM2(:) = s * DEL(:) / ( R(:) * ( 1.0_RP + Y(:)**2 ) )
        OM3(:) = s * Y(:) / DEL(:)
      else
        OM1(:) = - s * X(:) * DEL(:) / ( R(:) * ( 1.0_RP + X(:)**2 ) )
        OM2(:) = - s * Y(:) * DEL(:) / ( R(:) * ( 1.0_RP + Y(:)**2 ) )
        OM3(:) =   s / DEL(:)
      end if

      drho(:) = matmul(IntrpMat_VPOrdM1, DDENS_(:,ke))

      !-- Gradient hydrostatic pressure
      
      DPRES_hyd(:) = PRES_hyd(:,ke) - PRES_hyd_ref(:,ke)

      call sparsemat_matmul(Dx, GsqrtV(:) * DPRES_hyd(:), Fx)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,1) * DPRES_hyd(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,1), LiftDelFlx)
      GradPhyd_x(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      call sparsemat_matmul(Dy, GsqrtV(:) * DPRES_hyd(:), Fy)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,2) * DPRES_hyd(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,2), LiftDelFlx)
      GradPhyd_y(:) = lmesh%Escale(:,ke,2,2) * Fy(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)
      
      !-- DENS
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * MOMY_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( DDENS_(:,ke) + DENS_hyd(:,ke) ) * wt_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DENS_VID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)     &
          + lmesh%Escale(:,ke,2,2) * Fy(:)     &
          + lmesh%Escale(:,ke,3,3) * Fz(:)     &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      !-- MOMX
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMX_(:,ke) + G11(:) * DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMX_(:,ke) + G12(:) * DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMX_(:,ke)                              &
                                                    + ( lmesh%GI3(:,ke,1) * G11(:) + lmesh%GI3(:,ke,2) * G12(:) ) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                                &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                                &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                                &
            + LiftDelFlx(:)                   ) / lmesh%Gsqrt(:,ke)                         &
          - twoOVdel2(:) * Y(:) *                                                           & !-> metric terms
            ( X(:) * Y(:) * u_(:) - ( 1.0_RP + Y(:)**2 ) * v_(:) ) * MOMX_(:,ke)            & !
          - 2.0_RP * u_(:) * MOMZ_(:,ke) / R(:)                                             & !<-
          - ( G11(:) * GradPhyd_x(:) + G12(:) * GradPhyd_y(:) ) * RGsqrtV(:)                & !-> gradient hydrostatic pressure
          - lmesh%Gsqrt(:,ke) * (  G11(:) * ( OM2(:) * MOMZ_(:,ke) - OM3(:) * MOMY_(:,ke) ) & !-> Coriolis term
                                 - G12(:) * ( OM1(:) * MOMZ_(:,ke) - OM3(:) * MOMX_(:,ke) ) )

      !-- MOMY
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMY_(:,ke) + G12(:) * DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMY_(:,ke) + G22(:) * DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMY_(:,ke)                              &
                                                    + ( lmesh%GI3(:,ke,1) * G12(:) + lmesh%GI3(:,ke,2) * G22(:) ) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = &
            - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                                &
              + lmesh%Escale(:,ke,2,2) * Fy(:)                                                &
              + lmesh%Escale(:,ke,3,3) * Fz(:)                                                &
              + LiftDelFlx(:)                  ) / lmesh%Gsqrt(:,ke)                          &
            - twoOVdel2(:) * X(:) *                                                           & !-> metric terms
              ( - (1.0_RP + X(:)**2) * u_(:) + X(:) * Y(:) * v_(:) ) * MOMY_(:,ke)            & !
            - 2.0_RP * v_(:) * MOMZ_(:,ke) / R(:)                                             & !<-
            - ( G12(:) * GradPhyd_x(:) + G22(:) * GradPhyd_y(:) ) * RGsqrtV(:)                & !-> gradient hydrostatic pressure
            - lmesh%Gsqrt(:,ke) * (  G12(:) * ( OM2(:) * MOMZ_(:,ke) - OM3(:) * MOMY_(:,ke) ) & !-> Coriolis term
                                   - G22(:) * ( OM1(:) * MOMZ_(:,ke) - OM3(:) * MOMX_(:,ke) ) )

      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *   u_ (:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *   v_ (:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMZ_(:,ke) + RGsqrtV(:) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                                                    &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                                                    &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                                                    &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)                                                               &
          - 0.25_RP * R(:) * twoOVdel2(:)**2 * ( 1.0_RP * X(:)**2 ) * ( 1.0_RP * Y(:)**2 )                      & !-> metric terms
            * ( - ( 1.0_RP + X(:)**2 ) * MOMX_(:,ke) * u_(:)                                                    & !
                + 2.0_RP * X(:) * Y(:) * MOMX_(:,ke) * v_(:)                                                    & !
                - ( 1.0_RP + Y(:)**2 ) * MOMY_(:,ke) * v_(:)  )                                                 & !<-
          + 2.0_RP * DPRES_(:,ke) / R(:)                                                                        & !-> metric term with gradient of pressure deviaition
          - lmesh%Gsqrt(:,ke) * ( OM1(:) * MOMY_(:,ke) - OM2(:) * MOMX_(:,ke) )                                 & !-> Coriolis term
          - Grav * Rgam2(:) * drho(:)                                                                             !-> buoyancy term

      !-- RHOT
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * u_ (:) * RHOT_(:), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * v_ (:) * RHOT_(:), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * wt_(:) * RHOT_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,RHOT_VID), LiftDelFlx)
      
      RHOT_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 
    end do
    !$omp end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend_deep_atm

end module scale_atm_dyn_dgm_globalnonhydro3d_rhot_heve
