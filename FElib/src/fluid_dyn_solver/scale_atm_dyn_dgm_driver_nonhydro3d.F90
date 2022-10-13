!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics / DGM driver (nonydro3d)
!!
!! @par Description
!!      Driver module for dynamical core based on DGM 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_driver_nonhydro3d
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,  &
    Rdry => CONST_Rdry,  &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00  

  use scale_timeint_rk, only: TimeInt_RK
  use scale_sparsemat, only: SparseMat

  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField2D, MeshField3D

  use scale_element_modalfilter, only: &
    ModalFilter

  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo

  use scale_atm_dyn_dgm_driver_base, only: &
    AtmDynDGMDriver_base3D, AtmDynDGMDriver_base3D_Init, AtmDynDGMDriver_base3D_Final

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    PRGVAR_NUM, &
    DENS_VID => PRGVAR_DDENS_ID, THERM_VID => PRGVAR_THERM_ID,&
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID,                               &
    AUXVAR_NUM, &
    PRESHYD_VID => AUXVAR_PRESHYDRO_ID, DENSHYD_VID => AUXVAR_DENSHYDRO_ID,                 &
    CPTOT_VID => AUXVAR_CPtot_ID, CVTOT_VID => AUXVAR_CVtot_ID, RTOT_VID => AUXVAR_Rtot_ID, &
    PRES_VID => AUXVAR_PRES_ID,                                                             &
    PHYTEND_DENS_ID, PHYTEND_MOMX_ID, PHYTEND_MOMY_ID, PHYTEND_MOMZ_ID, PHYTEND_RHOT_ID,    &
    PHYTEND_RHOH_ID    


  use scale_atm_dyn_dgm_nonhydro3d_rhot_heve, only: &
    atm_dyn_dgm_nonhydro3d_rhot_heve_Init,          &
    atm_dyn_dgm_nonhydro3d_rhot_heve_Final,         &
    atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend

  use scale_atm_dyn_dgm_nonhydro3d_etot_heve, only: &
    atm_dyn_dgm_nonhydro3d_etot_heve_Init,          &
    atm_dyn_dgm_nonhydro3d_etot_heve_Final,         &
    atm_dyn_dgm_nonhydro3d_etot_heve_cal_tend
    
  use scale_atm_dyn_dgm_nonhydro3d_rhot_hevi, only: &
    atm_dyn_dgm_nonhydro3d_rhot_hevi_Init,          &
    atm_dyn_dgm_nonhydro3d_rhot_hevi_Final,         &
    atm_dyn_dgm_nonhydro3d_rhot_hevi_cal_tend,      &
    atm_dyn_dgm_nonhydro3d_rhot_hevi_cal_vi

  use scale_atm_dyn_dgm_nonhydro3d_etot_hevi, only: &
    atm_dyn_dgm_nonhydro3d_etot_hevi_Init,          &
    atm_dyn_dgm_nonhydro3d_etot_hevi_Final,         &
    atm_dyn_dgm_nonhydro3d_etot_hevi_cal_tend,      &
    atm_dyn_dgm_nonhydro3d_etot_hevi_cal_vi

  use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_splitform, only: &
    atm_dyn_dgm_nonhydro3d_rhot_heve_splitform_Init,          &
    atm_dyn_dgm_nonhydro3d_rhot_heve_splitform_Final,         &
    atm_dyn_dgm_nonhydro3d_rhot_heve_splitform_cal_tend

  use scale_atm_dyn_dgm_nonhydro3d_rhot_hevi_splitform, only: &
    atm_dyn_dgm_nonhydro3d_rhot_hevi_splitform_Init,          &
    atm_dyn_dgm_nonhydro3d_rhot_hevi_splitform_Final,         &
    atm_dyn_dgm_nonhydro3d_rhot_hevi_splitform_cal_tend,      &
    atm_dyn_dgm_nonhydro3d_rhot_hevi_splitform_cal_vi    
  
  use scale_atm_dyn_dgm_globalnonhydro3d_rhot_heve, only: &
    atm_dyn_dgm_globalnonhydro3d_rhot_heve_Init,          &
    atm_dyn_dgm_globalnonhydro3d_rhot_heve_Final,         &
    atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend

  use scale_atm_dyn_dgm_globalnonhydro3d_etot_heve, only: &
    atm_dyn_dgm_globalnonhydro3d_etot_heve_Init,          &
    atm_dyn_dgm_globalnonhydro3d_etot_heve_Final,         &
    atm_dyn_dgm_globalnonhydro3d_etot_heve_cal_tend

  use scale_atm_dyn_dgm_globalnonhydro3d_rhot_hevi, only:  &
    atm_dyn_dgm_globalnonhydro3d_rhot_hevi_Init,           &
    atm_dyn_dgm_globalnonhydro3d_rhot_hevi_Final,          &
    atm_dyn_dgm_globalnonhydro3d_rhot_hevi_cal_tend,       &
    atm_dyn_dgm_globalnonhydro3d_rhot_hevi_cal_vi

  use scale_atm_dyn_dgm_globalnonhydro3d_etot_hevi, only:  &
    atm_dyn_dgm_globalnonhydro3d_etot_hevi_Init,           &
    atm_dyn_dgm_globalnonhydro3d_etot_hevi_Final,          &
    atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_tend,       &
    atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_vi
  
  use scale_atm_dyn_dgm_bnd, only: AtmDynBnd
  use scale_atm_dyn_dgm_spongelayer, only: AtmDynSpongeLayer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  abstract interface    
    subroutine atm_dyn_nonhydro3d_cal_tend_ex( &
      DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                        & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, THERM_, DENS_hyd, PRES_hyd, CORIOLIS,  & ! (in)
      Rtot, CVtot, CPtot,                                                 & ! (in)
      Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )          ! (in)

      import RP
      import LocalMesh3D
      import ElementBase3D
      import LocalMesh2D
      import elementbase2D
      import SparseMat
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
      real(RP), intent(in)  :: THERM_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
      real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    end subroutine atm_dyn_nonhydro3d_cal_tend_ex
  end interface

  abstract interface    
    subroutine atm_dyn_nonhydro3d_cal_vi( &
      DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,             & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, THERM_, DENS_hyd, PRES_hyd, & ! (in)
      DDENS0_, MOMX0_, MOMY0_, MOMZ0_, THERM0_,                & ! (in) 
      Rtot, CVtot, CPtot,                                      & ! (in)
      Dz, Lift,                                                & ! (in)
      modalFilterFlag, VModalFilter,                           & ! (in)
      impl_fac, dt,                                            & ! (in)
      lmesh, elem, lmesh2D, elem2D )
  
      import RP
      import LocalMesh3D
      import ElementBase3D
      import LocalMesh2D
      import ElementBase2D
      import ModalFilter
      import SparseMat
      implicit none
  
      class(LocalMesh3D), intent(in) :: lmesh
      class(elementbase3D), intent(in) :: elem
      class(LocalMesh2D), intent(in) :: lmesh2D
      class(ElementBase2D), intent(in) :: elem2D
      real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMY_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: RHOT_dt(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: THERM_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DDENS0_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMX0_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMY0_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMZ0_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: THERM0_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
      class(SparseMat), intent(in) :: Dz, Lift
      logical, intent(in) :: modalFilterFlag
      class(ModalFilter), intent(in) :: VModalFilter
      real(RP), intent(in) :: impl_fac
      real(RP), intent(in) :: dt
    end subroutine atm_dyn_nonhydro3d_cal_vi
  end interface

  abstract interface    
    subroutine atm_dyn_nonhydro3d_final()
    end subroutine atm_dyn_nonhydro3d_final
  end interface 

  type, extends(AtmDynDGMDriver_base3D), public :: AtmDynDGMDriver_nonhydro3d
    integer :: EQS_TYPEID
    logical :: ENTOT_CONSERVE_SCHEME_FLAG

    ! element-wise modal filter
    logical :: MODALFILTER_FLAG
    type(ModalFilter) :: modal_filter_3d
    type(ModalFilter) :: modal_filter_v1D

    ! sponge layer
    type(AtmDynSpongeLayer) :: sponge_layer
    logical :: SPONGELAYER_FLAG

    ! prognositc variables

    ! diagnostic variables
    type(MeshFIeld3D) :: DPRES

    ! boundary condition
    type(AtmDynBnd) :: boundary_cond

    !
    logical :: hevi_flag

    procedure (atm_dyn_nonhydro3d_cal_vi), pointer, nopass :: cal_vi => null()
    procedure (atm_dyn_nonhydro3d_cal_tend_ex), pointer, nopass :: cal_tend_ex => null()
    procedure (atm_dyn_nonhydro3d_final), pointer, nopass :: dynsolver_final => null()
  contains
    procedure :: Init => AtmDynDGMDriver_nonhydro3d_Init
    procedure :: Final => AtmDynDGMDriver_nonhydro3d_Final
    procedure :: Update => AtmDynDGMDriver_nonhydro3d_update
    procedure :: calc_pressure => AtmDynDGMDriver_nonhydro3d_calc_pressure
  end type AtmDynDGMDriver_nonhydro3d

  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_HEVE             = 1
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_HEVE_ENTOT       = 2
  integer, public, parameter :: EQS_TYPEID_GLOBALNONHYD3D_HEVE       = 3
  integer, public, parameter :: EQS_TYPEID_GLOBALNONHYD3D_HEVE_ENTOT = 4
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_SPLITFORM_HEVE   = 5
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_HEVI             = 6
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_HEVI_ENTOT       = 7
  integer, public, parameter :: EQS_TYPEID_GLOBALNONHYD3D_HEVI       = 8
  integer, public, parameter :: EQS_TYPEID_GLOBALNONHYD3D_HEVI_ENTOT = 9
  integer, public, parameter :: EQS_TYPEID_NONHYD3D_SPLITFORM_HEVI   = 10  

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  private :: setup_modalfilter

contains
!OCL SERIAL  
  subroutine AtmDynDGMDriver_nonhydro3d_Init( this, &
    eqs_type_name, tint_type_name, dtsec,           &
    sponge_layer_flag, modal_filter_flag,           &
    mesh3D )
    implicit none

    class(AtmDynDGMDriver_nonhydro3d), intent(inout) :: this
    character(len=*), intent(in) :: eqs_type_name
    character(len=*), intent(in) :: tint_type_name
    real(DP), intent(in) :: dtsec
    logical, intent(in) :: sponge_layer_flag
    logical, intent(in) :: modal_filter_flag
    class(MeshBase3D), intent(in), target :: mesh3D


    class(LocalMesh3D), pointer :: lcmesh
    integer :: domID
    class(HexahedralElement), pointer :: refElem3D
    class(ElementBase), pointer :: refElem

    integer :: iv
    logical :: reg_file_hist
    !-----------------------------------------------------------------------------

    call AtmDynDGMDriver_base3D_Init( this, &
      PRGVAR_NUM,                           &
      tint_type_name, dtsec,                &
      mesh3D                                )

    this%ENTOT_CONSERVE_SCHEME_FLAG = .false.
    this%hevi_flag                  = .false.

    select case(eqs_type_name)
    !-- HEVE ------------------
    case("NONHYDRO3D_HEVE", "NONHYDRO3D_RHOT_HEVE")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_HEVE
      call atm_dyn_dgm_nonhydro3d_rhot_heve_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend
      this%cal_vi => null()
      this%dynsolver_final => atm_dyn_dgm_nonhydro3d_rhot_heve_Final
    case("NONHYDRO3D_ETOT_HEVE")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_HEVE_ENTOT
      call atm_dyn_dgm_nonhydro3d_etot_heve_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_nonhydro3d_etot_heve_cal_tend
      this%cal_vi => null()
      this%dynsolver_final => atm_dyn_dgm_nonhydro3d_etot_heve_Final
      this%ENTOT_CONSERVE_SCHEME_FLAG = .true.
    case("GLOBALNONHYDRO3D_HEVE", "GLOBALNONHYDRO3D_RHOT_HEVE")
      this%EQS_TYPEID = EQS_TYPEID_GLOBALNONHYD3D_HEVE
      call atm_dyn_dgm_globalnonhydro3d_rhot_heve_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_globalnonhydro3d_rhot_heve_cal_tend
      this%cal_vi => null()
      this%dynsolver_final => atm_dyn_dgm_globalnonhydro3d_rhot_heve_Final
    case("GLOBALNONHYDRO3D_ETOT_HEVE")
      this%EQS_TYPEID = EQS_TYPEID_GLOBALNONHYD3D_HEVE_ENTOT
      call atm_dyn_dgm_globalnonhydro3d_etot_heve_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_globalnonhydro3d_etot_heve_cal_tend
      this%cal_vi => null()
      this%dynsolver_final => atm_dyn_dgm_globalnonhydro3d_etot_heve_Final
      this%ENTOT_CONSERVE_SCHEME_FLAG = .true.
    case("NONHYDRO3D_HEVE_SPLITFORM", "NONHYDRO3D_RHOT_HEVE_SPLITFORM")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_SPLITFORM_HEVE
      call atm_dyn_dgm_nonhydro3d_rhot_heve_splitform_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_nonhydro3d_rhot_heve_splitform_cal_tend
      this%cal_vi => null()
      this%dynsolver_final => atm_dyn_dgm_nonhydro3d_rhot_heve_splitform_Final
    !-- HEVI ------------------
    case("NONHYDRO3D_HEVI", "NONHYDRO3D_RHOT_HEVI") 
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_HEVI
      call atm_dyn_dgm_nonhydro3d_rhot_hevi_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_nonhydro3d_rhot_hevi_cal_tend
      this%cal_vi => atm_dyn_dgm_nonhydro3d_rhot_hevi_cal_vi
      this%dynsolver_final => atm_dyn_dgm_nonhydro3d_rhot_hevi_Final
      this%hevi_flag = .true.
    case("NONHYDRO3D_ETOT_HEVI")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_HEVI_ENTOT
      call atm_dyn_dgm_nonhydro3d_etot_hevi_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_nonhydro3d_etot_hevi_cal_tend
      this%cal_vi => atm_dyn_dgm_nonhydro3d_etot_hevi_cal_vi
      this%dynsolver_final => atm_dyn_dgm_nonhydro3d_etot_hevi_Final
      this%ENTOT_CONSERVE_SCHEME_FLAG = .true.
      this%hevi_flag = .true.
    case("GLOBALNONHYDRO3D_HEVI", "GLOBALNONHYDRO3D_RHOT_HEVI")
      this%EQS_TYPEID = EQS_TYPEID_GLOBALNONHYD3D_HEVI
      call atm_dyn_dgm_globalnonhydro3d_rhot_hevi_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_globalnonhydro3d_rhot_hevi_cal_tend
      this%cal_vi => atm_dyn_dgm_globalnonhydro3d_rhot_hevi_cal_vi
      this%dynsolver_final => atm_dyn_dgm_globalnonhydro3d_rhot_hevi_Final
      this%hevi_flag = .true.
    case("GLOBALNONHYDRO3D_ETOT_HEVI")
      this%EQS_TYPEID = EQS_TYPEID_GLOBALNONHYD3D_HEVI_ENTOT
      call atm_dyn_dgm_globalnonhydro3d_etot_hevi_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_tend
      this%cal_vi => atm_dyn_dgm_globalnonhydro3d_etot_hevi_cal_vi
      this%dynsolver_final => atm_dyn_dgm_globalnonhydro3d_etot_hevi_Final
      this%ENTOT_CONSERVE_SCHEME_FLAG = .true.
      this%hevi_flag = .true.
    case("NONHYDRO3D_HEVI_SPLITFORM", "NONHYDRO3D_RHOT_HEVI_SPLITFORM")
      this%EQS_TYPEID = EQS_TYPEID_NONHYD3D_SPLITFORM_HEVI
      call atm_dyn_dgm_nonhydro3d_rhot_hevi_splitform_Init( mesh3D )
      this%cal_tend_ex => atm_dyn_dgm_nonhydro3d_rhot_hevi_splitform_cal_tend
      this%cal_vi => atm_dyn_dgm_nonhydro3d_rhot_hevi_splitform_cal_vi
      this%dynsolver_final => atm_dyn_dgm_nonhydro3d_rhot_hevi_splitform_Final
      this%hevi_flag = .true.
    case default
      LOG_ERROR("AtmDynDGMDriver_nonhydro3d_Init",*) 'Invalid EQS_TYPE in namelist PARAM_ATMOS_DYN. Check!'
      call PRC_abort
    end select
    
    !- Initialize variables

    call this%DPRES%Init( 'DPRES', 'Pa', mesh3D )

    !- initialize an object to manage boundary conditions
    call this%boundary_cond%Init()
    call this%boundary_cond%SetBCInfo( mesh3D )

    !- initialize an object to manage sponge layer
    this%SPONGELAYER_FLAG = sponge_layer_flag
    if (this%SPONGELAYER_FLAG) then
      call this%sponge_layer%Init( mesh3D, dtsec )
    end if

    !- initialize an object to manage modal filter
    this%MODALFILTER_FLAG = modal_filter_flag
    if (this%MODALFILTER_FLAG) then
      refElem => mesh3D%refElem
      select type(refElem)
      class is (HexahedralElement)
        refElem3D => refElem
      end select

      call setup_modalfilter( this, refElem3D )
    end if

    return
  end subroutine AtmDynDGMDriver_nonhydro3d_Init

!OCL SERIAL
  subroutine AtmDynDGMDriver_nonhydro3d_update( this, &
    PROG_VARS, AUX_VARS, PHYTENDS,                       &
    DENS_TRC, DENS0_TRC,                                 &
    MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg,               &
    ALPH_DENS_M_tavg, ALPH_DENS_P_tavg,                  &
    Coriolis,                                            &
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, mesh3D                 )

    use scale_tracer, only: &
      QA, TRACER_ADVC, TRACER_NAME
    
    use scale_atm_dyn_dgm_trcadvect3d_heve, only: &
      atm_dyn_dgm_trcadvect3d_save_massflux
  
    use scale_model_var_manager, only: ModelVarManager

    use scale_atm_dyn_dgm_modalfilter, only: &
      atm_dyn_dgm_modalfilter_apply
    
    implicit none

    class(AtmDynDGMDriver_nonhydro3d), intent(inout) :: this
    class(MeshBase3D), intent(in), target :: mesh3D
    class(ModelVarManager), intent(inout) :: PROG_VARS
    class(ModelVarManager), intent(inout) :: AUX_VARS
    class(ModelVarManager), intent(inout) :: PHYTENDS
    type(MeshField3D), intent(inout) :: DENS_TRC
    type(MeshField3D), intent(inout) :: DENS0_TRC
    type(MeshField3D), intent(inout) :: MFLX_x_tavg
    type(MeshField3D), intent(inout) :: MFLX_y_tavg
    type(MeshField3D), intent(inout) :: MFLX_z_tavg
    type(MeshField3D), intent(inout) :: ALPH_DENS_M_tavg
    type(MeshField3D), intent(inout) :: ALPH_DENS_P_tavg
    class(MeshField2D), intent(in) :: Coriolis
    type(SparseMat), intent(in) :: Dx, Dy, Dz
    type(SparseMat), intent(in) :: Sx, Sy, Sz
    type(SparseMat), intent(in) :: Lift

    integer :: rkstage
    integer :: tintbuf_ind
    real(RP) :: dt
    real(RP) :: implicit_fac

    real(RP) :: tavg_coef_MFLXZ(this%tint(1)%nstage)

    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    integer :: ke

    class(MeshField3D), pointer :: DDENS, MOMX, MOMY, MOMZ, THERM
    class(MeshField3D), pointer :: PRES_hyd, DENS_hyd, Rtot, CVtot, CPtot, PRES
    class(MeshField3D), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p

    !-----------------------------------------------------------------------------
    
    !-
    call PROF_rapstart( 'ATM_DYN_update_pre', 2)
    call PROG_VARS%Get3D(DENS_VID , DDENS)
    call PROG_VARS%Get3D(THERM_VID, THERM)
    call PROG_VARS%Get3D(MOMZ_VID , MOMZ )
    call PROG_VARS%Get3D(MOMX_VID , MOMX )
    call PROG_VARS%Get3D(MOMY_VID , MOMY )

    call AUX_VARS%Get3D( PRESHYD_VID, PRES_hyd )
    call AUX_VARS%Get3D( DENSHYD_VID, DENS_hyd )
    call AUX_VARS%Get3D( PRES_VID, PRES )
    call AUX_VARS%Get3D( RTOT_VID, Rtot )
    call AUX_VARS%Get3D( CVTOT_VID, CVtot )
    call AUX_VARS%Get3D( CPTOT_VID, CPtot )

    call PHYTENDS%Get3D( PHYTEND_DENS_ID, DENS_tp )
    call PHYTENDS%Get3D( PHYTEND_MOMX_ID, MOMX_tp )
    call PHYTENDS%Get3D( PHYTEND_MOMY_ID, MOMY_tp )
    call PHYTENDS%Get3D( PHYTEND_MOMZ_ID, MOMZ_tp )
    call PHYTENDS%Get3D( PHYTEND_RHOT_ID, RHOT_tp )
    call PHYTENDS%Get3D( PHYTEND_RHOH_ID, RHOH_p )

    if (this%hevi_flag) then 
      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh3D => mesh3D%lcmesh_list(n)
        
        call this%tint(n)%StoreVar0( DDENS%local(n)%val,  DENS_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )
        call this%tint(n)%StoreVar0( THERM%local(n)%val, THERM_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )
        call this%tint(n)%StoreVar0(  MOMZ%local(n)%val,  MOMZ_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )
        call this%tint(n)%StoreVar0(  MOMX%local(n)%val,  MOMX_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )
        call this%tint(n)%StoreVar0(  MOMY%local(n)%val,  MOMY_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )        
      end do
    end if
    call PROF_rapend( 'ATM_DYN_update_pre', 2)

    do rkstage=1, this%tint(1)%nstage
    
      if (this%hevi_flag) then        
        do n=1, mesh3D%LOCAL_MESH_NUM
          lcmesh3D => mesh3D%lcmesh_list(n)

          call PROF_rapstart( 'ATM_DYN_cal_vi', 2)
          implicit_fac = this%tint(n)%Get_implicit_diagfac(rkstage)
          tintbuf_ind = this%tint(n)%tend_buf_indmap(rkstage)
          dt = this%tint(n)%Get_deltime()
          call this%cal_vi( &
            this%tint(n)%tend_buf2D_im(:,:,DENS_VID,tintbuf_ind),                           & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,MOMX_VID ,tintbuf_ind),                          & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,MOMY_VID ,tintbuf_ind),                          & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,MOMZ_VID ,tintbuf_ind),                          & ! (out)
            this%tint(n)%tend_buf2D_im(:,:,THERM_VID,tintbuf_ind),                          & ! (out)
            DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val,    & ! (in)
            THERM%local(n)%val,                                                             & ! (in)
            DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,                                   & ! (in)
            this%tint(n)%var0_2D(:,:,DENS_VID), this%tint(n)%var0_2D(:,:,MOMX_VID),         & ! (in)
            this%tint(n)%var0_2D(:,:,MOMY_VID ), this%tint(n)%var0_2D(:,:,MOMZ_VID),        & ! (in)
            this%tint(n)%var0_2D(:,:,THERM_VID ),                                           & ! (in)
            Rtot%local(n)%val, CVtot%local(n)%val, CPtot%local(n)%val,                      & ! (in)
            Dz, Lift,                                                                       & ! (in)
            this%MODALFILTER_FLAG, this%modal_filter_v1D,                                   & ! (in)
            implicit_fac, dt,                                                               & ! (in)
            lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D    ) ! (in)
          
          call PROF_rapend( 'ATM_DYN_cal_vi', 2)  
          
          call PROF_rapstart( 'ATM_DYN_store_impl', 2)      
          call this%tint(n)%StoreImplicit( rkstage, DDENS%local(n)%val,  DENS_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )
          call this%tint(n)%StoreImplicit( rkstage, THERM%local(n)%val, THERM_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )
          call this%tint(n)%StoreImplicit( rkstage,  MOMZ%local(n)%val,  MOMZ_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )
          call this%tint(n)%StoreImplicit( rkstage,  MOMX%local(n)%val,  MOMX_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )
          call this%tint(n)%StoreImplicit( rkstage,  MOMY%local(n)%val,  MOMY_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )
          call PROF_rapend( 'ATM_DYN_store_impl', 2) 
        end do
      
      end if

      !* Exchange halo data
      call PROF_rapstart( 'ATM_DYN_exchange_prgv', 2)
      call PROG_VARS%MeshFieldComm_Exchange()
      call PROF_rapend( 'ATM_DYN_exchange_prgv', 2)

      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh3D => mesh3D%lcmesh_list(n)

        !* Apply boundary conditions
        call PROF_rapstart( 'ATM_DYN_applyBC_prgv', 2)
        call this%boundary_cond%ApplyBC_PROGVARS_lc( n,                                        & ! (in)
          DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val,         & ! (inout)
          THERM%local(n)%val,                                                                  & ! (inout)
          DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,                                        & ! (in)
          lcmesh3D%Gsqrt(:,:), lcmesh3D%GsqrtH(:,:), lcmesh3D%GI3(:,:,1), lcmesh3D%GI3(:,:,2), & ! (in)
          lcmesh3D%normal_fn(:,:,1), lcmesh3D%normal_fn(:,:,2), lcmesh3D%normal_fn(:,:,3),     & ! (in)
          lcmesh3D%vmapM, lcmesh3D%vmapP, lcmesh3D%vmapB,                                      & ! (in)
          lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D         ) ! (in)
        call PROF_rapend( 'ATM_DYN_applyBC_prgv', 2)
      end do


      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh3D => mesh3D%lcmesh_list(n)

        call PROF_rapstart( 'ATM_DYN_update_caltend_ex', 2)
        tintbuf_ind = this%tint(n)%tend_buf_indmap(rkstage)
        call this%cal_tend_ex( &
          this%tint(n)%tend_buf2D_ex(:,:,DENS_VID,tintbuf_ind),                            & ! (out)
          this%tint(n)%tend_buf2D_ex(:,:,MOMX_VID ,tintbuf_ind),                           & ! (out)
          this%tint(n)%tend_buf2D_ex(:,:,MOMY_VID ,tintbuf_ind),                           & ! (out)
          this%tint(n)%tend_buf2D_ex(:,:,MOMZ_VID ,tintbuf_ind),                           & ! (out)
          this%tint(n)%tend_buf2D_ex(:,:,THERM_VID,tintbuf_ind),                           & ! (out)
          DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val,     & ! (in)
          THERM%local(n)%val,                                                              & ! (in)
          DENS_hyd%local(n)%val, PRES_hyd%local(n)%val, Coriolis%local(n)%val,             & ! (in)
          Rtot%local(n)%val, CVtot%local(n)%val, CPtot%local(n)%val,                       & ! (in)
          Dx, Dy, Dz, Sx, Sy, Sz, Lift,                                                    & ! (in)
          lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D     ) 
        call PROF_rapend( 'ATM_DYN_update_caltend_ex', 2)

        !- Sponge layer
        if (this%SPONGELAYER_FLAG) then
          call PROF_rapstart('ATM_DYN_caltend_sponge', 2)
          call this%sponge_layer%AddTend( &
            this%tint(n)%tend_buf2D_ex(:,:,MOMX_VID ,tintbuf_ind),   & ! (inout)
            this%tint(n)%tend_buf2D_ex(:,:,MOMY_VID ,tintbuf_ind),   & ! (inout)
            this%tint(n)%tend_buf2D_ex(:,:,MOMZ_VID ,tintbuf_ind),   & ! (inout)
            MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val, & ! (in)
            lcmesh3D, lcmesh3D%refElem3D                             ) ! (in)
          call PROF_rapend('ATM_DYN_caltend_sponge', 2)
        end if
      end do
      
      call PROF_rapstart( 'ATM_DYN_update_add_tp', 2)      
      call calc_pressure( this, PRES, this%DPRES,      &
        DDENS, MOMX, MOMY, MOMZ, THERM,                &
        PRES_hyd, DENS_hyd, Rtot, CVtot, CPtot, mesh3D )

      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh3D => mesh3D%lcmesh_list(n)
        
        tintbuf_ind = this%tint(n)%tend_buf_indmap(rkstage)
        call add_phy_tend( &
          this, this%tint(n)%tend_buf2D_ex(:,:,:,tintbuf_ind),              & ! (inout)
          DENS_tp%local(n)%val, MOMX_tp%local(n)%val, MOMY_tp%local(n)%val, & ! (in)
          MOMZ_tp%local(n)%val, RHOT_tp%local(n)%val, RHOH_p %local(n)%val, & ! (in)
          PRES%local(n)%val, Rtot%local(n)%val, CPtot%local(n)%val,         & ! (in)
          n, lcmesh3D, lcmesh3D%refElem3D                                   ) ! (in)
      end do
      call PROF_rapend( 'ATM_DYN_update_add_tp', 2)

      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh3D => mesh3D%lcmesh_list(n)

        if ( QA > 0 ) then
          call PROF_rapstart( 'ATM_DYN_tavg_mflx', 2)
          if (this%tint(1)%imex_flag) then
            tavg_coef_MFLXZ(:) = this%tint(n)%coef_b_im(:)
          else
            tavg_coef_MFLXZ(:) = this%tint(n)%coef_b_ex(:)
          end if
          call atm_dyn_dgm_trcadvect3d_save_massflux( &
            MFLX_x_tavg%local(n)%val, MFLX_y_tavg%local(n)%val, MFLX_z_tavg%local(n)%val,   & ! (inout)
            ALPH_DENS_M_tavg%local(n)%face_val, ALPH_DENS_P_tavg%local(n)%face_val,         & ! (inout)
            DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val,    & ! (in)
            DENS_hyd%local(n)%val, PRES_hyd%local(n)%val, Coriolis%local(n)%val,            & ! (in)
            Rtot%local(n)%val, CVtot%local(n)%val, CPtot%local(n)%val,                      & ! (in)
            lcmesh3D, lcmesh3D%refElem3D,                                                   & ! (in)
            rkstage, this%tint(n)%coef_b_ex(rkstage), tavg_coef_MFLXZ(rkstage)              ) ! (in)
          call PROF_rapend( 'ATM_DYN_tavg_mflx', 2)
        end if
        
        call PROF_rapstart( 'ATM_DYN_update_advance', 2)      
        call this%tint(n)%Advance( rkstage, DDENS%local(n)%val,  DENS_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )        
        call this%tint(n)%Advance( rkstage, THERM%local(n)%val, THERM_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )    
        call this%tint(n)%Advance( rkstage,  MOMZ%local(n)%val,  MOMZ_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )    
        call this%tint(n)%Advance( rkstage,  MOMX%local(n)%val,  MOMX_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )    
        call this%tint(n)%Advance( rkstage,  MOMY%local(n)%val,  MOMY_VID, 1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE )
        call PROF_rapend( 'ATM_DYN_update_advance', 2)
      end do
    
    end do ! end for RK loop

    call PROF_rapstart( 'ATM_DYN_update_post', 2)
    if ( QA > 0 ) then
      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh3D => mesh3D%lcmesh_list(n)

        !$omp parallel do
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          DENS0_TRC%local(n)%val(:,ke) = this%tint(n)%var0_2D(:,ke,DENS_VID)
          DENS_TRC %local(n)%val(:,ke) = DDENS%local(n)%val(:,ke)
        end do        
      end do
    end if
    call PROF_rapend( 'ATM_DYN_update_post', 2)      

    if ( this%MODALFILTER_FLAG ) then
      do n=1, mesh3D%LOCAL_MESH_NUM
        call PROF_rapstart( 'ATM_DYN_update_modalfilter', 2)
        call atm_dyn_dgm_modalfilter_apply(  & 
          DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val, THERM%local(n)%val, & ! (inout)
          lcmesh3D, lcmesh3D%refElem3D, this%modal_filter_3d,        & ! (in)
!          do_weight_Gsqrt = .false.                                 ) ! (in)          
          do_weight_Gsqrt = .true.                                   ) ! (in)        
        call PROF_rapend( 'ATM_DYN_update_modalfilter', 2)
      end do
    end if

    call PROF_rapstart( 'ATM_DYN_update_post', 2)
    call calc_pressure( this, PRES, this%DPRES,      &
      DDENS, MOMX, MOMY, MOMZ, THERM,                &
      PRES_hyd, DENS_hyd, Rtot, CVtot, CPtot, mesh3D )    
    call PROF_rapend( 'ATM_DYN_update_post', 2)      

    return
  end subroutine AtmDynDGMDriver_nonhydro3d_update

!OCL SERIAL
  subroutine AtmDynDGMDriver_nonhydro3d_Final( this )
    implicit none

    class(AtmDynDGMDriver_nonhydro3d), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call this%dynsolver_final()
    call this%boundary_cond%Final()
    
    if ( this%SPONGELAYER_FLAG ) call this%sponge_layer%Final()

    if ( this%MODALFILTER_FLAG ) then
      call this%modal_filter_3d%Final()
      if ( this%hevi_flag ) call this%modal_filter_v1D%Final()
    end if
    
    call this%DPRES%Final()
    call AtmDynDGMDriver_base3D_Final( this )   

    return
  end subroutine AtmDynDGMDriver_nonhydro3d_Final

!OCL SERIAL
  subroutine AtmDynDGMDriver_nonhydro3d_calc_pressure( this, &
    PRES, PROG_VARS, AUX_VARS                                )

    implicit none

    class(AtmDynDGMDriver_nonhydro3d), intent(inout) :: this
    class(MeshField3D), intent(inout) :: PRES
    class(ModelVarManager), intent(inout) :: PROG_VARS
    class(ModelVarManager), intent(inout) :: AUX_VARS

    class(MeshBase3D), pointer :: mesh3D
    class(MeshField3D), pointer :: DDENS, MOMX, MOMY, MOMZ, THERM
    class(MeshField3D), pointer :: PRES_hyd, DENS_hyd, Rtot, CVtot, CPtot
    
    type(MeshField3D) :: DPRES
    !-----------------------------------------------------------------------------

    call PROG_VARS%Get3D(DENS_VID , DDENS)
    call PROG_VARS%Get3D(THERM_VID, THERM)
    call PROG_VARS%Get3D(MOMZ_VID , MOMZ )
    call PROG_VARS%Get3D(MOMX_VID , MOMX )
    call PROG_VARS%Get3D(MOMY_VID , MOMY )

    call AUX_VARS%Get3D( PRESHYD_VID, PRES_hyd )
    call AUX_VARS%Get3D( DENSHYD_VID, DENS_hyd )
    call AUX_VARS%Get3D( RTOT_VID, Rtot )
    call AUX_VARS%Get3D( CVTOT_VID, CVtot )
    call AUX_VARS%Get3D( CPTOT_VID, CPtot )

    mesh3D => DDENS%mesh
    call DPRES%Init( "DPRES", "Pa", mesh3D )

    call calc_pressure( this, PRES, DPRES,      &
      DDENS, MOMX, MOMY, MOMZ, THERM,                &
      PRES_hyd, DENS_hyd, Rtot, CVtot, CPtot, mesh3D )

    call DPRES%Final()

    return
  end subroutine AtmDynDGMDriver_nonhydro3d_calc_pressure
    
!----- private --------------------------------------

!OCL SERIAL
  subroutine add_phy_tend( this, & ! (in)
    dyn_tends,                                           & ! (inout)
    DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p, & !(in)
    PRES, Rtot, CPtot,                                   & ! (in)
    domID, lcmesh, elem3D                                ) ! (in)


    implicit none

    class(AtmDynDGMDriver_nonhydro3d), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(inout) :: dyn_tends(elem3D%Np,lcmesh%NeA,PRGVAR_NUM)
    real(RP), intent(in) :: DENS_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: RHOT_tp(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: RHOH_p (elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA)
    integer, intent(in) :: domID

    integer :: ke
    real(RP) :: EXNER(elem3D%Np)
    real(RP) :: rP0
    !---------------------------------------------------------------------------------

    rP0   = 1.0_RP / PRES00

    !$omp parallel private( EXNER, ke )
    
    !$omp do
    do ke=lcmesh%NeS, lcmesh%NeE
      dyn_tends(:,ke,DENS_VID) = dyn_tends(:,ke,DENS_VID) + DENS_tp(:,ke)
      dyn_tends(:,ke,MOMZ_VID) = dyn_tends(:,ke,MOMZ_VID) + MOMZ_tp(:,ke)
      dyn_tends(:,ke,MOMX_VID) = dyn_tends(:,ke,MOMX_VID) + MOMX_tp(:,ke)
      dyn_tends(:,ke,MOMY_VID) = dyn_tends(:,ke,MOMY_VID) + MOMY_tp(:,ke)
    end do
    !$omp end do

    if ( this%ENTOT_CONSERVE_SCHEME_FLAG ) then
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
        EXNER(:) = ( PRES(:,ke) * rP0 )**( Rtot(:,ke) / CPtot(:,ke) )

        dyn_tends(:,ke,THERM_VID) = dyn_tends(:,ke,THERM_VID) &
                                        + RHOH_p(:,ke)                             &
                                        + ( CPtot(:,ke) * EXNER(:) ) * RHOT_tp(:,ke)
      end do
      !$omp end do
    else
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
        EXNER(:) = ( PRES(:,ke) * rP0 )**( Rtot(:,ke) / CPtot(:,ke) )

        dyn_tends(:,ke,THERM_VID) = dyn_tends(:,ke,THERM_VID) &
                                        + RHOT_tp(:,ke)                             &
                                       + RHOH_p  (:,ke) / ( CPtot(:,ke) * EXNER(:) )
      end do
      !$omp end do
    end if
    !$omp end parallel
    return
  end subroutine add_phy_tend

!OCL SERIAL
  subroutine calc_pressure( this, PRES, DPRES, & ! (inout)
    DDENS, MOMX, MOMY, MOMZ, THERM,            & ! (in)
    PRES_hyd, DENS_hyd, Rtot, CVtot, CPtot,    & ! (in)
    mesh3D                                     ) ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      atm_dyn_dgm_nonhydro3d_common_DRHOT2PRES, &
      atm_dyn_dgm_nonhydro3d_common_EnTot2PRES
    implicit none
    class(AtmDynDGMDriver_nonhydro3d), intent(inout) :: this
    class(MeshField3D), intent(inout) :: PRES
    class(MeshField3D), intent(inout) :: DPRES
    class(MeshField3D), intent(in) :: DDENS
    class(MeshField3D), intent(in) :: MOMX
    class(MeshField3D), intent(in) :: MOMY
    class(MeshField3D), intent(in) :: MOMZ
    class(MeshField3D), intent(in) :: THERM
    class(MeshField3D), intent(in) :: PRES_hyd
    class(MeshField3D), intent(in) :: DENS_hyd
    class(MeshField3D), intent(in) :: Rtot
    class(MeshField3D), intent(in) :: CVtot
    class(MeshField3D), intent(in) :: CPtot
    class(MeshBase3D), intent(in), target :: mesh3D
    
    integer :: n
    class(LocalMesh3D), pointer :: lcmesh3D
    !---------------------------

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(n)

      if ( this%ENTOT_CONSERVE_SCHEME_FLAG ) then
        call atm_dyn_dgm_nonhydro3d_common_EnTot2PRES( PRES%local(n)%val, DPRES%local(n)%val,                  &
          DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val, THERM%local(n)%val,     &
          PRES_hyd%local(n)%val, DENS_hyd%local(n)%val, Rtot%local(n)%val, CVtot%local(n)%val,                 &
          lcmesh3D, lcmesh3D%refElem3D )
      else
        call atm_dyn_dgm_nonhydro3d_common_DRHOT2PRES( PRES%local(n)%val, DPRES%local(n)%val,                   &
          THERM%local(n)%val, PRES_hyd%local(n)%val, Rtot%local(n)%val, CVtot%local(n)%val, CPtot%local(n)%val, &
          lcmesh3D, lcmesh3D%refElem3D )
      end if
    end do

    return
  end subroutine calc_pressure

  !-- Setup modal filter
!OCL SERIAL
  subroutine setup_modalfilter( this, refElem3D )
    use scale_element_line, only: LineElement
    implicit none

    class(AtmDynDGMDriver_nonhydro3d), target, intent(inout) :: this
    class(HexahedralElement), target, intent(in) :: refElem3D

    real(RP) :: MF_ETAC_h  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_h = 36.0_RP
    integer  :: MF_ORDER_h = 16
    real(RP) :: MF_ETAC_v  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_v = 36.0_RP
    integer  :: MF_ORDER_v = 16

    namelist /PARAM_ATMOS_DYN_MODALFILTER/ &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    

    integer :: ierr

    type(LineElement) :: elemV1D
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_MODALFILTER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_setup_modalfilter",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_setup_modalfilter",*) 'Invalid names in namelist PARAM_ATMOS_DYN_MODALFILTER. Check!'
      call PRC_abort
    endif

    LOG_NML(PARAM_ATMOS_DYN_MODALFILTER)

    if ( this%hevi_flag ) then

      call this%modal_filter_3d%Init( &
        refElem3D,                           & ! (in)
        MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   & ! (in)
        MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    ) ! (in)
    else
      call this%modal_filter_3d%Init( &
        refElem3D,                    & ! (in)
        MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   & ! (in)
        1.0_RP, 0.0_RP, MF_ORDER_v           ) ! (in)

      call elemV1D%Init( refElem3D%PolyOrder_v, refElem3D%IsLumpedMatrix() )      
      call this%modal_filter_v1D%Init( elemV1D, &
        MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v,      &
        tend_flag = .true.                      )
      
      call elemV1D%Final()
    end if  

    return
  end subroutine setup_modalfilter

end module scale_atm_dyn_dgm_driver_nonhydro3d