!-------------------------------------------------------------------------------
!> module Coupler / Variables
!!
!! @par Description
!!          Module to manage variables with coupler component
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_cpl_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_debug

  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: &
    MeshBase2D
  use scale_mesh_base3d, only: &
    MeshBase3D,                              &
    DIMTYPE_XY  => MeshBase3D_DIMTYPEID_XY,  &
    DIMTYPE_XYZ  => MeshBase3D_DIMTYPEID_XYZ

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField2D, MeshField3D, MeshField3DList
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  !> Derived type to manage variables with coupler component
  type, public :: CouplerVars
    type(MeshField2D), allocatable :: OCN_vars(:)

    type(MeshField2D), allocatable :: OCN_atm_vars(:)

    real(RP) :: CNT_putATM_OCN  !< Counter for putting ATM variables to OCN component
    real(RP) :: CNT_putOCN_ATM  !< Counter for putting OCN variables to ATM component

    class(MeshBase3D), pointer :: mesh3D_atm => null()
    class(MeshBase2D), pointer :: mesh2D_atm => null()
    class(MeshBase2D), pointer :: mesh2D_ocn => null()
  contains
    procedure :: Init => CouplerVars_Init
    procedure :: Final => CouplerVars_Final
    procedure :: PutATM => CouplerVars_putATM
    procedure :: PutOCN => CouplerVars_putOCN
    procedure :: Get_SFC_ATM => CouplerVars_get_SFC_ATM
    procedure :: Get_ATM_OCN => CouplerVars_get_ATM_OCN
    !-
    procedure, private :: putATM_lc => putATM_local
    procedure, private :: putOCN_lc => putOCN_local
    procedure, private :: get_SFC_ATM_lc => get_SFC_ATM_local
    procedure, private :: get_ATM_OCN_lc => get_ATM_OCN_local
  end type CouplerVars

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  ! Input from ocean component
  integer, parameter, public :: OCN_SFC_TEMP_ID   = 1
  integer, parameter, public :: OCN_SFLX_MW_ID    = 2
  integer, parameter, public :: OCN_SFLX_MU_ID    = 3
  integer, parameter, public :: OCN_SFLX_MV_ID    = 4
  integer, parameter, public :: OCN_SFLX_SH_ID    = 5
  integer, parameter, public :: OCN_SFLX_LH_ID    = 6
  integer, parameter, public :: OCN_SFC_ALBEDO_ID = 7
  integer, parameter, public :: OCN_VAR_NUM       = 7

  ! Output to ocean component
  integer, parameter, public :: OCN_ATM_SFC_DENS_ID    = 1
  integer, parameter, public :: OCN_ATM_SFC_PRES_ID    = 2
  integer, parameter, public :: OCN_ATM_TEMP_ID        = 3
  integer, parameter, public :: OCN_ATM_PRES_ID        = 4
  integer, parameter, public :: OCN_ATM_W_ID           = 5
  integer, parameter, public :: OCN_ATM_U_ID           = 6
  integer, parameter, public :: OCN_ATM_V_ID           = 7
  integer, parameter, public :: OCN_ATM_QV_ID          = 8
  integer, parameter, public :: OCN_ATM_SFLX_RD_SW_DIR = 9
  integer, parameter, public :: OCN_ATM_SFLX_RD_LW_DIF = 10
  integer, parameter, public :: OCN_ATM_SFLX_ENGI      = 11
  integer, parameter, public :: OCN_ATM_VAR_NUM        = 11

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !

contains

!OCL SERIAL
  subroutine CouplerVars_Init(this, &
    mesh3D_atm, mesh2D_atm, mesh2D_ocn)
    implicit none
    class(CouplerVars), intent(inout) :: this
    class(MeshBase3D), intent(in), target :: mesh3D_atm
    class(MeshBase2D), intent(in), target :: mesh2D_atm
    class(MeshBase2D), intent(in), target :: mesh2D_ocn

    integer :: iv
    !--------------------------------------------------------------------------------

    allocate( this%OCN_vars(OCN_VAR_NUM) )
    do iv=1, OCN_VAR_NUM
      call this%OCN_vars(iv)%Init( "", "", mesh2D_ocn )
    end do

    allocate( this%OCN_atm_vars(OCN_ATM_VAR_NUM) )
    do iv=1, OCN_ATM_VAR_NUM
      call this%OCN_atm_vars(iv)%Init( "", "", mesh2D_atm )
    end do

    !-
    this%mesh3D_atm => mesh3D_atm
    this%mesh2D_atm => mesh2D_atm
    this%mesh2D_ocn => mesh2D_ocn
    
    !-
    this%CNT_putATM_OCN = 0.0_RP
    this%CNT_putOCN_ATM = 0.0_RP

    return
  end subroutine CouplerVars_Init

!OCL SERIAL
  subroutine CouplerVars_Final(this)
    implicit none
    class(CouplerVars), intent(inout) :: this

    integer :: iv
    !--------------------------------------------------------------------------------
    do iv=1, OCN_VAR_NUM
      call this%OCN_vars(iv)%Final()
    end do
    deallocate( this%OCN_vars )

    do iv=1, OCN_ATM_VAR_NUM
      call this%OCN_atm_vars(iv)%Final()
    end do
    deallocate( this%OCN_atm_vars )
    return
  end subroutine CouplerVars_Final

  !> Put ATM variables to CPL buffer for OCN component
  subroutine CouplerVars_putATM( this, &
    ATM_DDENS, ATM_MOMZ, ATM_MOMX, ATM_MOMY, ATM_QV, &
    ATM_PRES, ATM_Rtot, ATM_DENS_hyd, &
    SFLX_RAD_SW_dir_dn, SFLX_RAD_LW_dif_dn, SFLX_ENGI, &
    countup )
    use mod_atmos_vars_container, only: AtmosVarsContainer
    implicit none
    class(CouplerVars), intent(inout) :: this
    type(MeshField3D), intent(in), target :: ATM_DDENS
    type(MeshField3D), intent(in), target :: ATM_MOMZ
    type(MeshField3D), intent(in), target :: ATM_MOMX
    type(MeshField3D), intent(in), target :: ATM_MOMY
    type(MeshField3D), intent(in), target :: ATM_QV
    type(MeshField3D), intent(in), target :: ATM_PRES
    type(MeshField3D), intent(in), target :: ATM_Rtot
    type(MeshField3D), intent(in), target :: ATM_DENS_hyd
    type(MeshField2D), intent(in), target :: SFLX_RAD_SW_dir_dn
    type(MeshField2D), intent(in), target :: SFLX_RAD_LW_dif_dn
    type(MeshField2D), intent(in), target :: SFLX_ENGI
    logical, intent(in) :: countup

    class(LocalMesh3D), pointer :: lmesh
    class(LocalMesh2D), pointer :: lmesh2D
    integer :: ldomID
    real(RP), allocatable :: ATM_TEMP_lc(:,:)

    !----------------------------------------------------------------

    do ldomID=1, this%mesh2D_atm%LOCAL_MESH_NUM
      lmesh => ATM_DDENS%mesh%lcmesh_list(ldomID)
      lmesh2D => lmesh%lcmesh2D
      call this%putATM_lc( lmesh, lmesh%refElem3D, lmesh2D, lmesh2D%refElem2D, &
        ATM_DDENS%local(ldomID)%val, ATM_MOMZ%local(ldomID)%val,  ATM_MOMX%local(ldomID)%val, ATM_MOMY%local(ldomID)%val, &
        ATM_QV%local(ldomID)%val, ATM_PRES%local(ldomID)%val, ATM_Rtot%local(ldomID)%val, ATM_DENS_hyd%local(ldomID)%val, &
        SFLX_RAD_SW_dir_dn%local(ldomID)%val, SFLX_RAD_LW_dif_dn%local(ldomID)%val,                                       &
        SFLX_ENGI         %local(ldomID)%val,  &
        !-
        this%OCN_atm_vars(OCN_ATM_SFC_DENS_ID)%local(ldomID)%val,  &
        this%OCN_atm_vars(OCN_ATM_SFC_PRES_ID)%local(ldomID)%val,  &
        this%OCN_atm_vars(OCN_ATM_TEMP_ID)%local(ldomID)%val,      &
        this%OCN_atm_vars(OCN_ATM_PRES_ID)%local(ldomID)%val,      &
        this%OCN_atm_vars(OCN_ATM_W_ID)%local(ldomID)%val,         &
        this%OCN_atm_vars(OCN_ATM_U_ID)%local(ldomID)%val,         &
        this%OCN_atm_vars(OCN_ATM_V_ID)%local(ldomID)%val,         &
        this%OCN_atm_vars(OCN_ATM_QV_ID)%local(ldomID)%val,           &
        this%OCN_atm_vars(OCN_ATM_SFLX_RD_SW_DIR)%local(ldomID)%val,  &
        this%OCN_atm_vars(OCN_ATM_SFLX_RD_LW_DIF)%local(ldomID)%val,  &
        this%OCN_atm_vars(OCN_ATM_SFLX_ENGI)%local(ldomID)%val )
    end do
    
    ! Update counter
    if (countup) then
      this%CNT_putATM_OCN = this%CNT_putATM_OCN + 1.0_RP
    end if
    return
  end subroutine CouplerVars_putATM

  !> Put OCN variables to CPL buffer for ATM component
!OCL SERIAL
  subroutine CouplerVars_putOCN( this, ocn_vars, countup )
    use mod_ocean_vars, only: OceanVars, &
      SFC_TEMP_ID => AUXVAR2D_SFC_TEMP_ID, &
      SFLX_MW_ID => OCN_SFLX_MW_ID,        &
      SFLX_MU_ID => OCN_SFLX_MU_ID,        &
      SFLX_MV_ID => OCN_SFLX_MV_ID,        &
      SFLX_SH_ID => OCN_SFLX_SH_ID,        &
      SFLX_LH_ID => OCN_SFLX_LH_ID
    implicit none
    class(CouplerVars), intent(inout) :: this
    type(OceanVars), intent(in) :: ocn_vars
    logical, intent(in) :: countup

    class(LocalMesh2D), pointer :: lmesh
    integer :: ldomID    
    !----------------------------------------------------------------

    do ldomID=1, this%mesh2D_ocn%LOCAL_MESH_NUM
      lmesh => this%mesh2D_ocn%lcmesh_list(ldomID)

      call this%putOCN_lc( lmesh, lmesh%refElem2D, &
        ocn_vars%AUX_VARS2D(SFC_TEMP_ID)%local(ldomID)%val,  &
        ocn_vars%OCN_SFLX(SFLX_MW_ID)%local(ldomID)%val,     &
        ocn_vars%OCN_SFLX(SFLX_MU_ID)%local(ldomID)%val,     &
        ocn_vars%OCN_SFLX(SFLX_MV_ID)%local(ldomID)%val,     &
        ocn_vars%OCN_SFLX(SFLX_SH_ID)%local(ldomID)%val,     &
        ocn_vars%OCN_SFLX(SFLX_LH_ID)%local(ldomID)%val,     &
        !-
        this%OCN_vars(OCN_SFC_TEMP_ID)%local(ldomID)%val,    &
        this%OCN_vars(OCN_SFLX_MW_ID)%local(ldomID)%val,     &
        this%OCN_vars(OCN_SFLX_MU_ID)%local(ldomID)%val,     &
        this%OCN_vars(OCN_SFLX_MV_ID)%local(ldomID)%val,     &
        this%OCN_vars(OCN_SFLX_SH_ID)%local(ldomID)%val,     &
        this%OCN_vars(OCN_SFLX_LH_ID)%local(ldomID)%val      )
    end do

    ! Update counter
    if ( countup ) then
      this%CNT_putOCN_ATM = this%CNT_putOCN_ATM + 1.0_RP
    end if

    return
  end subroutine CouplerVars_putOCN

  !> Get SFC variables from CPL buffer for ATM component
!OCL SERIAL
  subroutine CouplerVars_get_SFC_ATM( this,     &
    SFC_TEMP, SFC_ALBEDO,                       &
    SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_LH )
    implicit none
    class(CouplerVars), intent(inout) :: this
    type(MeshField2D), intent(inout) :: SFC_TEMP
    type(MeshField2D), intent(inout) :: SFC_ALBEDO
    type(MeshField2D), intent(inout) :: SFLX_MW
    type(MeshField2D), intent(inout) :: SFLX_MU
    type(MeshField2D), intent(inout) :: SFLX_MV
    type(MeshField2D), intent(inout) :: SFLX_SH
    type(MeshField2D), intent(inout) :: SFLX_LH

    integer :: ldomID
    class(LocalMesh2D), pointer :: lmesh, lmesh_o
    !----------------------------------------------------------------

    do ldomID=1, this%mesh2D_atm%LOCAL_MESH_NUM
      lmesh => this%mesh2D_atm%lcmesh_list(ldomID)
      lmesh_o => this%mesh2D_ocn%lcmesh_list(ldomID)

      call this%get_SFC_ATM_lc( &
        lmesh, lmesh%refElem2D, lmesh_o, lmesh_o%refElem2D,       &
        this%OCN_vars(OCN_SFC_TEMP_ID)%local(ldomID)%val,         &
        this%OCN_vars(OCN_SFC_ALBEDO_ID)%local(ldomID)%val,       &
        this%OCN_vars(OCN_SFLX_MW_ID)%local(ldomID)%val,          &
        this%OCN_vars(OCN_SFLX_MU_ID)%local(ldomID)%val,          &
        this%OCN_vars(OCN_SFLX_MV_ID)%local(ldomID)%val,          &
        this%OCN_vars(OCN_SFLX_SH_ID)%local(ldomID)%val,          &
        this%OCN_vars(OCN_SFLX_LH_ID)%local(ldomID)%val,          &
        !-
        SFC_TEMP%local(ldomID)%val, SFC_ALBEDO%local(ldomID)%val, &
        SFLX_MW%local(ldomID)%val, SFLX_MU%local(ldomID)%val,     &
        SFLX_MV%local(ldomID)%val, SFLX_SH%local(ldomID)%val,     &
        SFLX_LH%local(ldomID)%val )
    end do
    return
  end subroutine CouplerVars_get_SFC_ATM

  !> Get ATM variables from CPL buffer for OCN component
!OCL SERIAL
  subroutine CouplerVars_get_ATM_OCN( this,          &
    ATM_SFC_DENS, ATM_SFC_PRES, ATM_TEMP, ATM_PRES, ATM_W, ATM_U, ATM_V, ATM_QV, &
    RD_SFLX_SW_DIR, RD_SFLX_LW_DIF, ZLEV_A )
    implicit none
    class(CouplerVars), intent(inout) :: this
    type(MeshField2D), intent(inout) :: ATM_SFC_DENS
    type(MeshField2D), intent(inout) :: ATM_SFC_PRES
    type(MeshField2D), intent(inout) :: ATM_TEMP
    type(MeshField2D), intent(inout) :: ATM_PRES
    type(MeshField2D), intent(inout) :: ATM_W
    type(MeshField2D), intent(inout) :: ATM_U
    type(MeshField2D), intent(inout) :: ATM_V
    type(MeshField2D), intent(inout) :: ATM_QV
    type(MeshField2D), intent(inout) :: RD_SFLX_SW_DIR
    type(MeshField2D), intent(inout) :: RD_SFLX_LW_DIF
    type(MeshField2D), intent(inout) :: ZLEV_A

    integer :: ldomID
    !----------------------------------------------------------------

    do ldomID=1, this%mesh2D_atm%LOCAL_MESH_NUM
      call this%get_ATM_OCN_lc( &
        this%mesh3D_atm%lcmesh_list(ldomID), this%mesh3D_atm%refElem3D, &
        this%mesh2D_atm%lcmesh_list(ldomID), this%mesh2D_atm%refElem2D, &
        this%mesh2D_ocn%lcmesh_list(ldomID), this%mesh2D_ocn%refElem2D, &
        !-
        this%OCN_atm_vars(OCN_ATM_SFC_DENS_ID)%local(ldomID)%val,     &
        this%OCN_atm_vars(OCN_ATM_SFC_PRES_ID)%local(ldomID)%val,     &
        this%OCN_atm_vars(OCN_ATM_TEMP_ID)%local(ldomID)%val,         &
        this%OCN_atm_vars(OCN_ATM_PRES_ID)%local(ldomID)%val,         &
        this%OCN_atm_vars(OCN_ATM_W_ID)%local(ldomID)%val,            &
        this%OCN_atm_vars(OCN_ATM_U_ID)%local(ldomID)%val,            &
        this%OCN_atm_vars(OCN_ATM_V_ID)%local(ldomID)%val,            &
        this%OCN_atm_vars(OCN_ATM_QV_ID)%local(ldomID)%val,           &
        this%OCN_atm_vars(OCN_ATM_SFLX_RD_SW_DIR)%local(ldomID)%val,  &
        this%OCN_atm_vars(OCN_ATM_SFLX_RD_LW_DIF)%local(ldomID)%val,  &
        !-
        ATM_SFC_DENS%local(ldomID)%val, ATM_SFC_PRES%local(ldomID)%val, &
        ATM_TEMP%local(ldomID)%val, ATM_PRES%local(ldomID)%val,                                              &
        ATM_W%local(ldomID)%val, ATM_U%local(ldomID)%val, ATM_V%local(ldomID)%val, ATM_QV%local(ldomID)%val, &
        RD_SFLX_SW_DIR%local(ldomID)%val,  &
        RD_SFLX_LW_DIF%local(ldomID)%val, &
        ZLEV_A%local(ldomID)%val )
    end do

    this%CNT_putATM_OCN = 0.0_RP
    return
  end subroutine CouplerVars_get_ATM_OCN

!-
!OCL SERIAL
  subroutine putATM_local( this, &
    lmesh_a, elem_a, lmesh2D_a, elem2D_a,&
    ATM_DDENS, ATM_MOMZ, ATM_MOMX, ATM_MOMY, ATM_QV,             &
    ATM_PRES, ATM_Rtot, ATM_DENS_hyd,                            &
    ATM_RD_SFLX_SW_DIR, ATM_RD_SFLX_LW_DIF, ATM_SFLX_ENGI,       &
    O_ATM_SFC_DENS, O_ATM_SFC_PRES,  &
    O_ATM_TEMP, O_ATM_PRES, O_ATM_W, O_ATM_U, O_ATM_V, O_ATM_QV, &
    O_ATM_RD_SFLX_SW_DIR, O_ATM_RD_SFLX_LW_DIF, O_ATM_SFLX_ENGI  )
    implicit none
    class(CouplerVars), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lmesh_a
    class(ElementBase3D), intent(in) :: elem_a
    class(LocalMesh2D), intent(in) :: lmesh2D_a
    class(ElementBase2D), intent(in) :: elem2D_a
    real(RP), intent(in) :: ATM_DDENS(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_MOMZ(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_MOMX(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_MOMY(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_QV(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_PRES(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_Rtot(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_DENS_hyd(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_RD_SFLX_SW_DIR(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(in) :: ATM_RD_SFLX_LW_DIF(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(in) :: ATM_SFLX_ENGI(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_SFC_DENS(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_SFC_PRES(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_TEMP(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_PRES(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_W(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_U(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_V(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_QV(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_RD_SFLX_SW_DIR(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_RD_SFLX_LW_DIF(elem2D_a%Np,lmesh2D_a%NeA)
    real(RP), intent(out) :: O_ATM_SFLX_ENGI(elem2D_a%Np,lmesh2D_a%NeA)

    integer :: ke2D
    integer :: coef1, coef2

    integer :: hSlice0(elem2D_a%Np), hSlice1(elem2D_a%Np)
    real(RP) :: dens_(elem2D_a%Np)
    real(RP) :: temp_(elem2D_a%Np)
    !----------------------------------------------------------------

    coef1 = this%CNT_putATM_OCN
    coef2 = 1.0_RP / ( 1.0_RP + coef1 )

    hSlice0(:) = elem_a%Hslice(:,1)
    hSlice1(:) = elem_a%Hslice(:,2)

    !$omp parallel private(temp_, dens_)
    !$omp do
    do ke2D=lmesh2D_a%NeS, lmesh2D_a%NeE
      dens_(:) = ATM_DENS_hyd(hSlice0(:),ke2D) + ATM_DDENS(hSlice0(:),ke2D)
      temp_(:) = ATM_PRES(hSlice0(:),ke2D) / (dens_(:) * ATM_Rtot(hSlice0(:),ke2D))

      O_ATM_SFC_DENS(:,ke2D) = O_ATM_SFC_DENS(:,ke2D) * coef1 + dens_(:)
      O_ATM_SFC_PRES(:,ke2D) = O_ATM_SFC_PRES(:,ke2D) * coef1 + ATM_PRES(hSlice0(:),ke2D)

      O_ATM_TEMP(:,ke2D) = O_ATM_TEMP(:,ke2D) * coef1 + temp_(:)
      O_ATM_PRES(:,ke2D) = O_ATM_PRES(:,ke2D) * coef1 + ATM_PRES(hSlice0(:),ke2D)
      O_ATM_W(:,ke2D)    = O_ATM_W(:,ke2D) * coef1 + ATM_MOMZ(hSlice0(:),ke2D) / dens_(:)
      O_ATM_U(:,ke2D)    = O_ATM_U(:,ke2D) * coef1 + ATM_MOMX(hSlice0(:),ke2D) / dens_(:)
      O_ATM_V(:,ke2D)    = O_ATM_V(:,ke2D) * coef1 + ATM_MOMY(hSlice0(:),ke2D) / dens_(:)
      O_ATM_QV(:,ke2D)   = O_ATM_QV(:,ke2D) * coef1 + ATM_QV(hSlice0(:),ke2D)

      O_ATM_RD_SFLX_SW_DIR(:,ke2D) = O_ATM_RD_SFLX_SW_DIR(:,ke2D) * coef1 + ATM_RD_SFLX_SW_DIR(:,ke2D)
      O_ATM_RD_SFLX_LW_DIF(:,ke2D) = O_ATM_RD_SFLX_LW_DIF(:,ke2D) * coef1 + ATM_RD_SFLX_LW_DIF(:,ke2D)
      O_ATM_SFLX_ENGI(:,ke2D) = O_ATM_SFLX_ENGI(:,ke2D) * coef1 + ATM_SFLX_ENGI(:,ke2D)
    end do
    !$omp do
    do ke2D=lmesh2D_a%NeS, lmesh2D_a%NeE
      O_ATM_SFC_DENS(:,ke2D) = O_ATM_SFC_DENS(:,ke2D) * coef2
      O_ATM_SFC_PRES(:,ke2D) = O_ATM_SFC_PRES(:,ke2D) * coef2

      O_ATM_TEMP(:,ke2D) = O_ATM_TEMP(:,ke2D) * coef2
      O_ATM_PRES(:,ke2D) = O_ATM_PRES(:,ke2D) * coef2
      O_ATM_W(:,ke2D)    = O_ATM_W(:,ke2D)    * coef2
      O_ATM_U(:,ke2D)    = O_ATM_U(:,ke2D)    * coef2
      O_ATM_V(:,ke2D)    = O_ATM_V(:,ke2D)    * coef2
      O_ATM_QV(:,ke2D)   = O_ATM_QV(:,ke2D)   * coef2

      O_ATM_RD_SFLX_SW_DIR(:,ke2D) = O_ATM_RD_SFLX_SW_DIR(:,ke2D) * coef2
      O_ATM_RD_SFLX_LW_DIF(:,ke2D) = O_ATM_RD_SFLX_LW_DIF(:,ke2D) * coef2
      O_ATM_SFLX_ENGI(:,ke2D) = O_ATM_SFLX_ENGI(:,ke2D) * coef2
    end do
    !$omp end parallel
    return
  end subroutine putATM_local

!OCL SERIAL
  subroutine putOCN_local( this, &
    lmesh_o, elem_o, &
    SFC_TEMP, SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_LH, &
    SFC_TEMP_out, SFLX_MW_out, SFLX_MU_out, SFLX_MV_out, SFLX_SH_out, SFLX_LH_out )
    implicit none
    class(CouplerVars), intent(inout) :: this
    class(LocalMesh2D), intent(in) :: lmesh_o
    class(ElementBase2D), intent(in) :: elem_o
    real(RP), intent(in) :: SFC_TEMP(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: SFLX_MW(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: SFLX_MU(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: SFLX_MV(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: SFLX_SH(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: SFLX_LH(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: SFC_TEMP_out(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: SFLX_MW_out(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: SFLX_MU_out(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: SFLX_MV_out(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: SFLX_SH_out(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: SFLX_LH_out(elem_o%Np,lmesh_o%NeA)

    integer :: ke
    real(RP) :: coef1, coef2
    !----------------------------------------------------------------

    coef1 = this%CNT_putOCN_ATM
    coef2 = 1.0_RP / ( 1.0_RP + coef1 )

    !$omp parallel
    !$omp do
    do ke=lmesh_o%NeS, lmesh_o%NeE
      SFC_TEMP_out(:,ke) = SFC_TEMP_out(:,ke) * coef1 + SFC_TEMP(:,ke)
      SFLX_MW_out(:,ke) = SFLX_MW_out(:,ke) * coef1 + SFLX_MW(:,ke)
      SFLX_MU_out(:,ke) = SFLX_MU_out(:,ke) * coef1 + SFLX_MU(:,ke)
      SFLX_MV_out(:,ke) = SFLX_MV_out(:,ke) * coef1 + SFLX_MV(:,ke)
      SFLX_SH_out(:,ke) = SFLX_SH_out(:,ke) * coef1 + SFLX_SH(:,ke)
      SFLX_LH_out(:,ke) = SFLX_LH_out(:,ke) * coef1 + SFLX_LH(:,ke)
    end do
    !$omp do
    do ke=lmesh_o%NeS, lmesh_o%NeE
      SFC_TEMP_out(:,ke) = SFC_TEMP_out(:,ke) * coef2
      SFLX_MW_out(:,ke) = SFLX_MW_out(:,ke) * coef2
      SFLX_MU_out(:,ke) = SFLX_MU_out(:,ke) * coef2
      SFLX_MV_out(:,ke) = SFLX_MV_out(:,ke) * coef2
      SFLX_SH_out(:,ke) = SFLX_SH_out(:,ke) * coef2
      SFLX_LH_out(:,ke) = SFLX_LH_out(:,ke) * coef2
    end do
    !$omp end parallel
    return
  end subroutine putOCN_local
  
!OCL SERIAL
  subroutine get_SFC_ATM_local( this, &
    lmesh_a, elem_a, lmesh_o, elem_o, &
    O_SFC_TEMP, O_SFC_ALBEDO, O_SFLX_MW, O_SFLX_MU, O_SFLX_MV, O_SFLX_SH, O_SFLX_LH, &
    SFC_TEMP, SFC_ALBEDO, SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_LH )
    implicit none
    class(CouplerVars), intent(inout) :: this
    class(LocalMesh2D), intent(in) :: lmesh_a
    class(ElementBase2D), intent(in) :: elem_a
    class(LocalMesh2D), intent(in) :: lmesh_o
    class(ElementBase2D), intent(in) :: elem_o
    real(RP), intent(in) :: O_SFC_TEMP(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: O_SFC_ALBEDO(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: O_SFLX_MW(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: O_SFLX_MU(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: O_SFLX_MV(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: O_SFLX_SH(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: O_SFLX_LH(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: SFC_TEMP(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: SFC_ALBEDO(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: SFLX_MW(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: SFLX_MU(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: SFLX_MV(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: SFLX_SH(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: SFLX_LH(elem_a%Np,lmesh_a%NeA)

    integer :: ke
    !----------------------------------------------------------------

    !$omp parallel do
    do ke=lmesh_a%NeS, lmesh_a%NeE
      SFC_TEMP(:,ke)  = O_SFC_TEMP(:,ke)
      SFC_ALBEDO(:,ke) = O_SFC_ALBEDO(:,ke)
      SFLX_MW(:,ke) = O_SFLX_MW(:,ke)
      SFLX_MU(:,ke) = O_SFLX_MU(:,ke)
      SFLX_MV(:,ke) = O_SFLX_MV(:,ke)
      SFLX_SH(:,ke) = O_SFLX_SH(:,ke)
      SFLX_LH(:,ke) = O_SFLX_LH(:,ke)
    end do
    return
  end subroutine get_SFC_ATM_local

!OCL SERIAL
  subroutine get_ATM_OCN_local( this, &
    lmesh3D_a, elem3D_a, lmesh_a, elem_a, lmesh_o, elem_o, &
    O_ATM_SFC_DENS, O_ATM_SFC_PRES, O_ATM_TEMP, O_ATM_PRES, O_ATM_W, O_ATM_U, O_ATM_V, O_ATM_QV, &
!    O_ATM_SFLX_ENGI, &
    O_ATM_RD_SFLX_SW_DIR, O_ATM_RD_SFLX_LW_DIF, &
    ATM_SFC_DENS, ATM_SFC_PRES, ATM_TEMP, ATM_PRES, ATM_W, ATM_U, ATM_V, ATM_QV, &
!    ATM_SFLX_ENGI, &
    ATM_RD_SFLX_SW_DIR, ATM_RD_SFLX_LW_DIF, &
    ZLEV_A )
    implicit none
    class(CouplerVars), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lmesh3D_a
    class(ElementBase3D), intent(in) :: elem3D_a
    class(LocalMesh2D), intent(in) :: lmesh_a
    class(ElementBase2D), intent(in) :: elem_a
    class(LocalMesh2D), intent(in) :: lmesh_o
    class(ElementBase2D), intent(in) :: elem_o
    real(RP), intent(in) :: O_ATM_SFC_DENS(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_SFC_PRES(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_TEMP(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_PRES(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_W(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_U(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_V(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_QV(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(in) :: O_ATM_SFLX_ENGI(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_RD_SFLX_SW_DIR(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_RD_SFLX_LW_DIF(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: ATM_SFC_DENS(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_SFC_PRES(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_TEMP(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_PRES(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_W(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_U(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_V(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_QV(elem_o%Np,lmesh_o%NeA)
    ! real(RP), intent(out) :: ATM_SFLX_ENGI(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_RD_SFLX_SW_DIR(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_RD_SFLX_LW_DIF(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ZLEV_A(elem_o%Np,lmesh_o%NeA)

    integer :: ke
    !----------------------------------------------------------------

    !$omp parallel do
    do ke=lmesh_a%NeS, lmesh_a%NeE
      ATM_SFC_DENS(:,ke) = O_ATM_SFC_DENS(:,ke)
      ATM_SFC_PRES(:,ke) = O_ATM_SFC_PRES(:,ke)

      ATM_TEMP(:,ke)  = O_ATM_TEMP(:,ke)
      ATM_PRES(:,ke)  = O_ATM_PRES(:,ke)
      ATM_W(:,ke)     = O_ATM_W(:,ke)
      ATM_U(:,ke)     = O_ATM_U(:,ke)
      ATM_V(:,ke)     = O_ATM_V(:,ke)
      ATM_QV(:,ke)    = O_ATM_QV(:,ke)
      ! ATM_SFLX_ENGI(:,ke) = O_ATM_SFLX_ENGI(:,ke)
      ATM_RD_SFLX_SW_DIR(:,ke) = O_ATM_RD_SFLX_SW_DIR(:,ke)
      ATM_RD_SFLX_LW_DIF(:,ke) = O_ATM_RD_SFLX_LW_DIF(:,ke)

      ZLEV_A(:,ke) = lmesh3D_a%zlev(elem3D_a%Hslice(:,2),ke)
    end do
    return
  end subroutine get_ATM_OCN_local
end module mod_cpl_vars