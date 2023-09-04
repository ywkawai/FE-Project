!-------------------------------------------------------------------------------
!> module FElib / Physics turbulence / Atmosphere / driver 
!!
!! @par Description
!!      Driver module for turbulent model based on DGM 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_tb_dgm_driver
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_tracer, only: QA

  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,  &
    Rdry => CONST_Rdry,  &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00  

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
  
  use scale_model_mesh_manager, only: ModelMesh3D
  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo

  use scale_atm_dyn_dgm_bnd, only: AtmDynBnd

  use scale_atm_phy_tb_dgm_smg, only: &
    atm_phy_tb_dgm_smg_Init,          &
    atm_phy_tb_dgm_smg_Final,         &
    atm_phy_tb_dgm_smg_cal_grad,      &
    atm_phy_tb_dgm_smg_cal_tend,      &
    atm_phy_tb_dgm_smg_cal_grad_qtrc, &
    atm_phy_tb_dgm_smg_cal_tend_qtrc

  use scale_atm_phy_tb_dgm_globalsmg, only: &
    atm_phy_tb_dgm_globalsmg_Init,          &
    atm_phy_tb_dgm_globalsmg_Final,         &
    atm_phy_tb_dgm_globalsmg_cal_grad,      &
    atm_phy_tb_dgm_globalsmg_cal_tend,      &
    atm_phy_tb_dgm_globalsmg_cal_grad_qtrc, &
    atm_phy_tb_dgm_globalsmg_cal_tend_qtrc

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    PRGVAR_NUM, &
    DENS_VID => PRGVAR_DDENS_ID, THERM_VID => PRGVAR_THERM_ID,&
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID,                               &
    AUXVAR_NUM, &
    PRESHYD_VID => AUXVAR_PRESHYDRO_ID, DENSHYD_VID => AUXVAR_DENSHYDRO_ID, &
    PRES_VID => AUXVAR_PRES_ID, PT_VID => AUXVAR_PT_ID

  use scale_atm_phy_tb_dgm_common, only: &
    ATMOS_PHY_TB_TENDS_NUM1, &
    TB_MOMX_t_VID => ATMOS_PHY_TB_MOMX_t_ID,TB_MOMY_t_VID => ATMOS_PHY_TB_MOMY_t_ID,  &
    TB_MOMZ_t_VID => ATMOS_PHY_TB_MOMZ_t_ID, TB_RHOT_t_VID => ATMOS_PHY_TB_RHOT_t_ID, &
    ATMOS_PHY_TB_AUX_NUM,                                                   &
    T11_VID => ATMOS_PHY_TB_AUX_T11_ID, T12_VID => ATMOS_PHY_TB_AUX_T12_ID, T13_VID => ATMOS_PHY_TB_AUX_T13_ID, &
    T21_VID => ATMOS_PHY_TB_AUX_T21_ID, T22_VID => ATMOS_PHY_TB_AUX_T22_ID, T23_VID => ATMOS_PHY_TB_AUX_T23_ID, &
    T31_VID => ATMOS_PHY_TB_AUX_T31_ID, T32_VID => ATMOS_PHY_TB_AUX_T32_ID, T33_VID => ATMOS_PHY_TB_AUX_T33_ID, &
    DPTDX_VID => ATMOS_PHY_TB_AUX_DPTDX_ID, DPTDY_VID => ATMOS_PHY_TB_AUX_DPTDY_ID, &
    DPTDZ_VID => ATMOS_PHY_TB_AUX_DPTDZ_ID,                                         &
    NU_VID => ATMOS_PHY_TB_DIAG_NU_ID, KH_VID => ATMOS_PHY_TB_DIAG_KH_ID,           &
    TKE_VID => ATMOS_PHY_TB_DIAG_TKE_ID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  abstract interface    
    subroutine atm_phy_tb_final()
      implicit none
    end subroutine atm_phy_tb_final
  end interface 

  abstract interface
    subroutine atm_phy_tb_cal_grad( &
      T11, T12, T13, T21, T22, T23, T31, T32, T33,                & ! (out)
      dPTdx, dPTdy, dPTdz,                                        & ! (out)
      TKE, Nu, Kh,                                                & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,    & ! (in)
      PRES, PT,                                                   & ! (in)
      Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
      is_bound                                                    ) ! (in)
      import RP
      import LocalMesh3D
      import LocalMesh2D
      import ElementBase3D
      import ElementBase2D
      import SparseMat
      implicit none

      class(LocalMesh3D), intent(in) :: lmesh
      class(ElementBase3D), intent(in) :: elem
      class(LocalMesh2D), intent(in) :: lmesh2D
      class(ElementBase2D), intent(in) :: elem2D
      real(RP), intent(out) :: T11(elem%Np,lmesh%NeA), T12(elem%Np,lmesh%NeA), T13(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: T21(elem%Np,lmesh%NeA), T22(elem%Np,lmesh%NeA), T23(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: T31(elem%Np,lmesh%NeA), T32(elem%Np,lmesh%NeA), T33(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: dPTdx(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: dPTdy(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: dPTdz(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: TKE(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: Nu(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: Kh(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PRES(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PT(elem%Np,lmesh%NeA)
      type(SparseMat), intent(in) :: Dx, Dy, Dz
      type(SparseMat), intent(in) :: Sx, Sy, Sz
      type(SparseMat), intent(in) :: Lift
      logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)
    end subroutine atm_phy_tb_cal_grad
  end interface

  abstract interface
    subroutine atm_phy_tb_cal_grad_qtrc( &
      dQTdx, dQTdy, dQTdz,                                        & ! (out)
      dRdx, dRdy, dRdz,                                           & ! (inout)
      QTRC, DDENS, DENS_hyd,                                      & ! (in)
      Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
      is_bound, cal_grad_dens                                     ) ! (in)
      import RP
      import LocalMesh3D
      import LocalMesh2D
      import ElementBase3D
      import ElementBase2D
      import SparseMat
      implicit none

      class(LocalMesh3D), intent(in) :: lmesh
      class(elementbase3D), intent(in) :: elem
      class(LocalMesh2D), intent(in) :: lmesh2D
      class(elementbase2D), intent(in) :: elem2D
      real(RP), intent(out) :: dQTdx(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: dQTdy(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: dQTdz(elem%Np,lmesh%NeA)
      real(RP), intent(inout) :: dRdx(elem%Np,lmesh%NeA)
      real(RP), intent(inout) :: dRdy(elem%Np,lmesh%NeA)
      real(RP), intent(inout) :: dRdz(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: QTRC(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DDENS(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)        
      type(SparseMat), intent(in) :: Dx, Dy, Dz
      type(SparseMat), intent(in) :: Sx, Sy, Sz
      type(SparseMat), intent(in) :: Lift
      logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)
      logical, intent(in) :: cal_grad_dens
    end subroutine atm_phy_tb_cal_grad_qtrc
  end interface

  abstract interface
    subroutine atm_phy_tb_cal_tend( &
      MOMX_t, MOMY_t, MOMZ_t, RHOT_t,                             & ! (out)
      T11, T12, T13, T21, T22, T23, T31, T32, T33,                & ! (in)
      dPTdx, dPTdy, dPTdz,                                        & ! (in)
      Nu, Kh,                                                     & ! (in)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,                        & ! (in)
      DENS_hyd, PRES_hyd,  PRES_, PT_,                            & ! (in)
      Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
      is_bound                                                    ) ! (in)

    import RP
    import LocalMesh3D
    import LocalMesh2D
    import ElementBase3D
    import ElementBase2D
    import SparseMat
    implicit none
      
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
    real(RP), intent(out) :: MOMX_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_t(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: T11(elem%Np,lmesh%NeA), T12(elem%Np,lmesh%NeA), T13(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: T21(elem%Np,lmesh%NeA), T22(elem%Np,lmesh%NeA), T23(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: T31(elem%Np,lmesh%NeA), T32(elem%Np,lmesh%NeA), T33(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dPTdx(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dPTdy(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dPTdz(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Nu   (elem%Np,lmesh%NeA) ! Eddy viscosity
    real(RP), intent(in)  :: Kh   (elem%Np,lmesh%NeA) ! Eddy diffusivity
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_ (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_ (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_ (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PT_  (elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx, Dy, Dz
    type(SparseMat), intent(in) :: Sx, Sy, Sz
    type(SparseMat), intent(in) :: Lift
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)

    end subroutine atm_phy_tb_cal_tend
  end interface

  abstract interface 
    subroutine atm_phy_tb_cal_tend_qtrc( &
      RHOQ_t,                                                     & ! (out)
      dQTdx, dQTdy, dQTdz,                                        & ! (in)
      Kh, DDENS_,DENS_hyd,                                        & ! (in)
      Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
      is_bound                                                    ) ! (in)

      import RP
      import LocalMesh3D
      import LocalMesh2D
      import ElementBase3D
      import ElementBase2D
      import SparseMat
      implicit none

      class(LocalMesh3D), intent(in) :: lmesh
      class(elementbase3D), intent(in) :: elem
      class(LocalMesh2D), intent(in) :: lmesh2D
      class(elementbase2D), intent(in) :: elem2D
      real(RP), intent(out) :: RHOQ_t(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: dQTdx(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: dQTdy(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: dQTdz(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: Kh   (elem%Np,lmesh%NeA) ! Eddy diffusivity
      real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
      type(SparseMat), intent(in) :: Dx, Dy, Dz
      type(SparseMat), intent(in) :: Sx, Sy, Sz
      type(SparseMat), intent(in) :: Lift
      logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)
    end subroutine atm_phy_tb_cal_tend_qtrc
  end interface 

  type, public :: AtmPhyTbDGMDriver
    !
    integer :: TB_TYPEID

    !
    type(MeshField3D), allocatable :: auxtrcvars(:)
    type(ModelVarManager) :: auxtrcvars_manager
    integer :: auxtrcvars_commid

    ! diagnostic variables
    type(MeshField3D), allocatable :: AUX_TBVARS3D(:)
    type(ModelVarManager) :: AUXTBVAR3D_manager
    integer :: AUXTBVAR3D_commid

    !
    type(MeshField3D) :: GRAD_DENS(3)

    procedure (atm_phy_tb_cal_grad), pointer, nopass :: tbsolver_cal_grad => null()
    procedure (atm_phy_tb_cal_grad_qtrc), pointer, nopass :: tbsolver_cal_grad_qtrc => null()
    procedure (atm_phy_tb_final), pointer, nopass :: tbsolver_final => null()
    procedure (atm_phy_tb_cal_tend), pointer, nopass :: tbsolver_cal_tend => null()
    procedure (atm_phy_tb_cal_tend_qtrc), pointer, nopass :: tbsolver_cal_tend_qtrc => null()
  contains
    procedure :: Init => AtmPhyTbDGMDriver_Init
    procedure :: Final => AtmPhyTbDGMDriver_Final
    procedure :: Tendency => AtmPhyTbDGMDriver_tendency
  end type AtmPhyTbDGMDriver

  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  integer, parameter :: TB_TYPEID_SMAGORINSKY         = 1
  integer, parameter :: TB_TYPEID_SMAGORINSKY_GLOBAL  = 2

  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_DQTDZ_ID   = 1  
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_DQTDX_ID   = 2
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_DQTDY_ID   = 3
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_SCALAR_NUM = 1 
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_HVEC_NUM   = 1
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_NUM        = 3
  type(VariableInfo), public :: ATMOS_PHY_TB_AUXTRC_VINFO(ATMOS_PHY_TB_AUXTRC_NUM)
  DATA ATMOS_PHY_TB_AUXTRC_VINFO / &
    VariableInfo( ATMOS_PHY_TB_AUXTRC_DQTDZ_ID, 'DQTDZ', 'gradient of QTRC (z)',      &
                  'K/m',  3, 'XYZ',  ''                                           ),  &
    VariableInfo( ATMOS_PHY_TB_AUXTRC_DQTDX_ID, 'DQTDX', 'gradient of QTRC (x)',      &
                  'K/m',  3, 'XYZ',  ''                                           ),  &
    VariableInfo( ATMOS_PHY_TB_AUXTRC_DQTDY_ID, 'DQTDY', 'gradient of QTRC (y)',      &
                  'K/m',  3, 'XYZ',  ''                                           )   /

contains
!OCL SERIAL  
  subroutine AtmPhyTbDGMDriver_Init( this, &
    tb_type_name, dtsec,                   &
    model_mesh3D )
    use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
    implicit none

    class(AtmPhyTbDGMDriver), intent(inout) :: this
    character(len=*), intent(in) :: tb_type_name
    real(DP), intent(in) :: dtsec
    class(ModelMesh3D), intent(inout), target :: model_mesh3D

    class(MeshBase3D), pointer :: mesh3D
    class(MeshCubedSphereDom3D), pointer :: gm_mesh3D

    integer :: iv
    logical :: reg_file_hist
                    
    integer :: idim
    !-----------------------------------------------------------------------------

    mesh3D => model_mesh3D%ptr_mesh

    !--- Set the type of turbulence scheme
    
    select case(tb_type_name)
    case ('SMAGORINSKY')
      this%TB_TYPEID = TB_TYPEID_SMAGORINSKY

      call atm_phy_tb_dgm_smg_Init( mesh3D )
      this%tbsolver_cal_grad => atm_phy_tb_dgm_smg_cal_grad
      this%tbsolver_cal_grad_qtrc => atm_phy_tb_dgm_smg_cal_grad_qtrc
      this%tbsolver_cal_tend => atm_phy_tb_dgm_smg_cal_tend
      this%tbsolver_cal_tend_qtrc => atm_phy_tb_dgm_smg_cal_tend_qtrc
      this%tbsolver_final => atm_phy_tb_dgm_smg_Final
    case ('SMAGORINSKY_GLOBAL')
      this%TB_TYPEID = TB_TYPEID_SMAGORINSKY_GLOBAL

      select type(mesh3D)
      class is (MeshCubedSphereDom3D)
        gm_mesh3D => mesh3D
      end select
      call atm_phy_tb_dgm_globalsmg_Init( mesh3D, gm_mesh3D%shallow_approx )
      this%tbsolver_cal_grad => atm_phy_tb_dgm_globalsmg_cal_grad
      this%tbsolver_cal_grad_qtrc => atm_phy_tb_dgm_globalsmg_cal_grad_qtrc
      this%tbsolver_cal_tend => atm_phy_tb_dgm_globalsmg_cal_tend
      this%tbsolver_cal_tend_qtrc => atm_phy_tb_dgm_globalsmg_cal_tend_qtrc
      this%tbsolver_final => atm_phy_tb_dgm_globalsmg_Final
    case default
      LOG_ERROR("AtmPhyTbDGMDriver_Init",*) 'Invalid TB_TYPE in namelist PARAM_ATMOS_TB. Check!'
      call PRC_abort
    end select
    
    !- Initialize variables

    !-
    if ( QA > 1 ) then
      call this%auxtrcvars_manager%Init()
      allocate( this%auxtrcvars(ATMOS_PHY_TB_AUXTRC_NUM) )

      reg_file_hist = .true.    
      do iv = 1, ATMOS_PHY_TB_AUXTRC_NUM
        call this%auxtrcvars_manager%Regist( &
          ATMOS_PHY_TB_AUXTRC_VINFO(iv), mesh3D,  & ! (in) 
          this%auxtrcvars(iv),                    & ! (inout)
          reg_file_hist, fill_zero=.true.            ) ! (in)
      end do

      call model_mesh3D%Create_communicator( &
        ATMOS_PHY_TB_AUXTRC_SCALAR_NUM, ATMOS_PHY_TB_AUXTRC_HVEC_NUM, 0, & ! (in) 
        this%auxtrcvars_manager,         & ! (inout)
        this%auxtrcvars(:),              & ! (in)
        this%auxtrcvars_commid           ) ! (out)

      do idim=1, 3
        call this%GRAD_DENS(idim)%Init( "Grad_DENS", "kg/m4", mesh3D )      
      end do
    end if

    return
  end subroutine AtmPhyTbDGMDriver_Init

!OCL SERIAL
  subroutine AtmPhyTbDGMDriver_tendency( this, &
    TB_TENDS, PROG_VARS, TRC_VARS, AUX_VARS,   &
    AUX_TB_VARS,  DIAG_TB_VARS,                            &
    boundary_cond,                                         &
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, mesh3D                   )

    use scale_tracer, only: &
      QA, TRACER_ADVC, TRACER_NAME
      
    use scale_model_var_manager, only: ModelVarManager
    implicit none

    class(AtmPhyTbDGMDriver), intent(inout) :: this
    class(MeshBase3D), intent(in), target :: mesh3D
    class(ModelVarManager), intent(inout) :: TB_TENDS
    class(ModelVarManager), intent(inout) :: PROG_VARS
    class(ModelVarManager), intent(inout) :: TRC_VARS
    class(ModelVarManager), intent(inout) :: AUX_VARS
    class(ModelVarManager), intent(inout) :: AUX_TB_VARS
    class(ModelVarManager), intent(inout) :: DIAG_TB_VARS
    type(AtmDynBnd), intent(in) :: boundary_cond
    type(SparseMat), intent(in) :: Dx, Dy, Dz
    type(SparseMat), intent(in) :: Sx, Sy, Sz
    type(SparseMat), intent(in) :: Lift

    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    integer :: ke

    class(MeshField3D), pointer :: DDENS, MOMX, MOMY, MOMZ, THERM
    class(MeshField3D), pointer :: PRES_hyd, DENS_hyd, PRES, PT
    class(MeshField3D), pointer :: T11, T12, T13, T21, T22, T23, T31, T32, T33
    class(MeshField3D), pointer :: dPTdx, dPTdy, dPTdz
    class(MeshField3D), pointer :: dQTdx, dQTdy, dQTdz
    class(MeshField3D), pointer :: Nu, Kh, TKE
    class(MeshField3D), pointer :: tb_MOMX_t, tb_MOMY_t, tb_MOMZ_t, tb_RHOT_t
    class(MeshField3D), pointer :: QTRC
    class(MeshField3D), pointer :: tb_RHOQ_t

    type DYN_BNDInfo
      logical, allocatable :: is_bound(:,:)
    end type
    type(DYN_BNDInfo), allocatable :: bnd_info(:)

    integer :: iq
    logical :: cal_grad_flag
    !-----------------------------------------------------------------------------
    
    !-
    call PROF_rapstart( 'ATM_TB_tendency_pre', 2)
    call PROG_VARS%Get3D(DENS_VID , DDENS)
    call PROG_VARS%Get3D(THERM_VID, THERM)
    call PROG_VARS%Get3D(MOMZ_VID , MOMZ )
    call PROG_VARS%Get3D(MOMX_VID , MOMX )
    call PROG_VARS%Get3D(MOMY_VID , MOMY )

    call AUX_VARS%Get3D( PRESHYD_VID, PRES_hyd )
    call AUX_VARS%Get3D( DENSHYD_VID, DENS_hyd )
    call AUX_VARS%Get3D( PRES_VID, PRES )
    call AUX_VARS%Get3D( PT_VID, PT )

    call AUX_TB_VARS%Get3D( T11_VID, T11 )
    call AUX_TB_VARS%Get3D( T12_VID, T12 )
    call AUX_TB_VARS%Get3D( T13_VID, T13 )
    call AUX_TB_VARS%Get3D( T21_VID, T21 )
    call AUX_TB_VARS%Get3D( T22_VID, T22 )
    call AUX_TB_VARS%Get3D( T23_VID, T23 )
    call AUX_TB_VARS%Get3D( T31_VID, T31 )
    call AUX_TB_VARS%Get3D( T32_VID, T32 )
    call AUX_TB_VARS%Get3D( T33_VID, T33 )
    call AUX_TB_VARS%Get3D( DPTDX_VID, dPTdx )
    call AUX_TB_VARS%Get3D( DPTDY_VID, dPTdy )
    call AUX_TB_VARS%Get3D( DPTDZ_VID, dPTdz )

    call DIAG_TB_VARS%Get3D( TKE_VID, TKE )
    call DIAG_TB_VARS%Get3D( NU_VID, Nu )
    call DIAG_TB_VARS%Get3D( KH_VID, Kh )

    call TB_TENDS%Get3D( TB_MOMX_t_VID, tb_MOMX_t )
    call TB_TENDS%Get3D( TB_MOMY_t_VID, tb_MOMY_t )
    call TB_TENDS%Get3D( TB_MOMZ_t_VID, tb_MOMZ_t )
    call TB_TENDS%Get3D( TB_RHOT_t_VID, tb_RHOT_t )

    allocate( bnd_info(mesh3D%LOCAL_MESH_NUM) )
    call PROF_rapend( 'ATM_TB_tendency_pre', 2)

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(n)

      !- Apply boundary conditions
      call PROF_rapstart('ATM_PHY_TB_inquire_bnd', 2)
      allocate( bnd_info(n)%is_bound(lcmesh3D%refElem3D%NfpTot,lcmesh3D%Ne) )
      call boundary_cond%Inquire_bound_flag(  bnd_info(n)%is_bound, &
        n, lcmesh3D%VMapM, lcmesh3D%VMapP, lcmesh3D%VMapB,          &
        lcmesh3D, lcmesh3D%refElem3D                                )
      call PROF_rapend('ATM_PHY_TB_inquire_bnd', 2)

      call PROF_rapstart('ATM_PHY_TB_cal_grad', 2)
      call this%tbsolver_cal_grad( &
          T11%local(n)%val, T12%local(n)%val, T13%local(n)%val,         & ! (out)
          T21%local(n)%val, T22%local(n)%val, T23%local(n)%val,         & ! (out)
          T31%local(n)%val, T32%local(n)%val, T33%local(n)%val,         & ! (out)
          dPTdx%local(n)%val, dPTdy%local(n)%val, dPTdz%local(n)%val,   & ! (out)
          TKE%local(n)%val, Nu%local(n)%val, Kh%local(n)%val,           & ! (out)
          DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val, THERM%local(n)%val,            & ! (in)
          DENS_hyd%local(n)%val, PRES_hyd%local(n)%val, PRES%local(n)%val, PT%local(n)%val,                           & ! (in)
          Dx, Dy, Dz, Sx, Sy, Sz, Lift, lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D, & ! (in)
          bnd_info(n)%is_bound                                                                                        ) ! (in)
      call PROF_rapend('ATM_PHY_TB_cal_grad', 2)
    end do

    !* Exchange halo data
    call PROF_rapstart('ATM_PHY_TB_exchange_prgv', 2)
    call AUX_TB_VARS%MeshFieldComm_Exchange()
    call PROF_rapend('ATM_PHY_TB_exchange_prgv', 2)

    do n=1, mesh3D%LOCAL_MESH_NUM
      call PROF_rapstart('ATM_PHY_TB_cal_tend', 2)
      call this%tbsolver_cal_tend( &
        tb_MOMX_t%local(n)%val, tb_MOMY_t%local(n)%val, tb_MOMZ_t%local(n)%val, tb_RHOT_t%local(n)%val,             & ! (out)
        T11%local(n)%val, T12%local(n)%val, T13%local(n)%val,                                                       & ! (in)
        T21%local(n)%val, T22%local(n)%val, T23%local(n)%val,                                                       & ! (in)
        T31%local(n)%val, T32%local(n)%val, T33%local(n)%val,                                                       & ! (in)
        dPTdx%local(n)%val, dPTdy%local(n)%val, dPTdz%local(n)%val, Nu%local(n)%val, Kh%local(n)%val,               &
        DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val, THERM%local(n)%val,            &
        DENS_hyd%local(n)%val, PRES_hyd%local(n)%val, PRES%local(n)%val, PT%local(n)%val,                           &
        Dx, Dy, Dz, Sx, Sy, Sz, Lift, lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D, &
        bnd_info(n)%is_bound                                                                                        )
      call PROF_rapend('ATM_PHY_TB_cal_tend', 2)
    end do

    cal_grad_flag = .true.
    do iq = 1, QA
      if ( .not. TRACER_ADVC(iq) ) cycle
      
      call TRC_VARS%Get3D( iq, QTRC )
      call TB_TENDS%Get3D( ATMOS_PHY_TB_TENDS_NUM1 + iq, tb_RHOT_t )
      
      do n=1, mesh3D%LOCAL_MESH_NUM
        call PROF_rapstart('ATM_PHY_TB_cal_grad_qtrc', 2)
        call this%tbsolver_cal_grad_qtrc( &
          dQTdx%local(n)%val, dQTdy%local(n)%val, dQTdz%local(n)%val,      & ! (out)
          this%GRAD_DENS(1)%local(n)%val, this%GRAD_DENS(2)%local(n)%val,  & ! (inout)
          this%GRAD_DENS(3)%local(n)%val,                                  & ! (inout)
          QTRC%local(n)%val, DDENS%local(n)%val, DENS_hyd%local(n)%val,    & ! (in) 
          Dx, Dy, Dx, Sx, Sy, Sz, Lift, lcmesh3D, lcmesh3D%refElem3D,      & ! (in)
          lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D,                  & ! (in)
          bnd_info(n)%is_bound, cal_grad_flag )                              ! (in)
        call PROF_rapend('ATM_PHY_TB_cal_grad_qtrc', 2)
      end do
      cal_grad_flag = .false.

      call this%auxtrcvars_manager%MeshFieldComm_Exchange()

      do n=1, mesh3D%LOCAL_MESH_NUM
        call PROF_rapstart('ATM_PHY_TB_cal_tend_qtrc', 2)
        call this%tbsolver_cal_tend_qtrc( tb_RHOQ_t%local(n)%val,          & ! (out)
          dQTdx%local(n)%val, dQTdy%local(n)%val, dQTdz%local(n)%val,      & ! (out)
          Kh%local(n)%val, DDENS%local(n)%val, DENS_hyd%local(n)%val,      & ! (in) 
          Dx, Dy, Dx, Sx, Sy, Sz, Lift, lcmesh3D, lcmesh3D%refElem3D,      & ! (in)
          lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D,                  & ! (in)
          bnd_info(n)%is_bound )                                             ! (in)
        call PROF_rapend('ATM_PHY_TB_cal_tend_qtrc', 2)
      end do        
    end do

    do n = 1, mesh3D%LOCAL_MESH_NUM
      deallocate( bnd_info(n)%is_bound )
    end do
  
    return
  end subroutine AtmPhyTbDGMDriver_tendency

!OCL SERIAL
  subroutine AtmPhyTbDGMDriver_Final( this )
    implicit none

    class(AtmPhyTbDGMDriver), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call this%tbsolver_final()

    if ( QA > 0 ) then
      call this%auxtrcvars_manager%Final()
    end if

    return
  end subroutine AtmPhyTbDGMDriver_Final

end module scale_atm_phy_tb_dgm_driver