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
  integer, parameter, public :: OCN_SFC_ALBEDO_ID = 2
  integer, parameter, public :: OCN_VAR_NUM       = 2

  ! Output to ocean component
  integer, parameter, public :: OCN_ATM_TEMP_ID        = 1
  integer, parameter, public :: OCN_ATM_PRES_ID        = 2
  integer, parameter, public :: OCN_ATM_W_ID           = 3
  integer, parameter, public :: OCN_ATM_U_ID           = 4
  integer, parameter, public :: OCN_ATM_V_ID           = 5
  integer, parameter, public :: OCN_ATM_SFLX_RD_SW_DIR = 6
  integer, parameter, public :: OCN_ATM_SFLX_RD_LW_DIF = 7
  integer, parameter, public :: OCN_ATM_SFLX_ENGI      = 8
  integer, parameter, public :: OCN_ATM_VAR_NUM        = 8

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !

contains

!OCL SERIAL
  subroutine CouplerVars_Init(this, &
    mesh2D_atm, mesh2D_ocn)
    implicit none
    class(CouplerVars), intent(inout) :: this
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
    SFLX_RAD_SW_dir_dn, SFLX_RAD_LW_dif_dn, SFLX_ENGI, &
    countup )
    use mod_atmos_vars_container, only: AtmosVarsContainer
    implicit none
    class(CouplerVars), intent(inout) :: this
    type(MeshField2D), intent(in), target :: SFLX_RAD_SW_dir_dn
    type(MeshField2D), intent(in), target :: SFLX_RAD_LW_dif_dn
    type(MeshField2D), intent(in), target :: SFLX_ENGI
    logical, intent(in) :: countup

    class(LocalMesh2D), pointer :: lmesh
    integer :: ldomID
    !----------------------------------------------------------------

    do ldomID=1, this%mesh2D_atm%LOCAL_MESH_NUM
        lmesh => this%mesh2D_atm%lcmesh_list(ldomID)
        call this%putATM_lc( lmesh, lmesh%refElem2D, &
            ! atm_vars%AUX_VARS2D(ATM_TEMP_ID)%local(ldomID)%val,  &
            ! atm_vars%AUX_VARS2D(ATM_PRES_ID)%local(ldomID)%val,  &
            ! atm_vars%AUX_VARS2D(ATM_W_ID)%local(ldomID)%val,     &
            ! atm_vars%AUX_VARS2D(ATM_U_ID)%local(ldomID)%val,     &
            ! atm_vars%AUX_VARS2D(ATM_V_ID)%local(ldomID)%val,     &
            SFLX_RAD_SW_dir_dn%local(ldomID)%val,  &
            SFLX_RAD_LW_dif_dn%local(ldomID)%val,  &
            SFLX_ENGI         %local(ldomID)%val,  &
            !-
            ! this%OCN_atm_vars(OCN_ATM_TEMP_ID)%local(ldomID)%val,  &
            ! this%OCN_atm_vars(OCN_ATM_PRES_ID)%local(ldomID)%val,  &
            ! this%OCN_atm_vars(OCN_ATM_W_ID)%local(ldomID)%val,     &
            ! this%OCN_atm_vars(OCN_ATM_U_ID)%local(ldomID)%val,     &
            ! this%OCN_atm_vars(OCN_ATM_V_ID)%local(ldomID)%val,     &
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
        SFC_TEMP_ID => AUXVAR2D_SFC_TEMP_ID
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
        !-
        this%OCN_vars(OCN_SFC_TEMP_ID)%local(ldomID)%val )
    end do

    ! Update counter
    if ( countup ) then
      this%CNT_putOCN_ATM = this%CNT_putOCN_ATM + 1.0_RP
    end if

    return
  end subroutine CouplerVars_putOCN

  !> Get SFC variables from CPL buffer for ATM component
!OCL SERIAL
  subroutine CouplerVars_get_SFC_ATM( this, &
    SFC_TEMP, SFC_ALBEDO )
    implicit none
    class(CouplerVars), intent(inout) :: this
    type(MeshField2D), intent(inout) :: SFC_TEMP
    type(MeshField2D), intent(inout) :: SFC_ALBEDO

    integer :: ldomID
    class(LocalMesh2D), pointer :: lmesh, lmesh_o
    !----------------------------------------------------------------

    do ldomID=1, this%mesh2D_atm%LOCAL_MESH_NUM
      lmesh => this%mesh2D_atm%lcmesh_list(ldomID)
      lmesh_o => this%mesh2D_ocn%lcmesh_list(ldomID)

      call this%get_SFC_ATM_lc( &
        lmesh, lmesh%refElem2D, lmesh_o, lmesh_o%refElem2D,      &
        this%OCN_vars(OCN_SFC_TEMP_ID)%local(ldomID)%val,        &
        this%OCN_vars(OCN_SFC_ALBEDO_ID)%local(ldomID)%val,      &
        !-
        SFC_TEMP%local(ldomID)%val, SFC_ALBEDO%local(ldomID)%val )
    end do
    return
  end subroutine CouplerVars_get_SFC_ATM

  !> Get ATM variables from CPL buffer for OCN component
!OCL SERIAL
  subroutine CouplerVars_get_ATM_OCN( this, &
    RD_SFLX_SW_DIR, RD_SFLX_LW_DIF )
    implicit none
    class(CouplerVars), intent(inout) :: this
    type(MeshField2D), intent(inout) :: RD_SFLX_SW_DIR
    type(MeshField2D), intent(inout) :: RD_SFLX_LW_DIF

    integer :: ldomID
    !----------------------------------------------------------------

    do ldomID=1, this%mesh2D_atm%LOCAL_MESH_NUM
      call this%get_ATM_OCN_lc( &
        this%mesh2D_atm%lcmesh_list(ldomID), this%mesh2D_atm%refElem2D, &
        this%mesh2D_ocn%lcmesh_list(ldomID), this%mesh2D_ocn%refElem2D, &
        !-
        this%OCN_atm_vars(OCN_ATM_SFLX_RD_SW_DIR)%local(ldomID)%val,  &
        this%OCN_atm_vars(OCN_ATM_SFLX_RD_LW_DIF)%local(ldomID)%val,  &
        !-
        RD_SFLX_SW_DIR%local(ldomID)%val,  &
        RD_SFLX_LW_DIF%local(ldomID)%val   )
    end do


    this%CNT_putATM_OCN = 0.0_RP
    return
  end subroutine CouplerVars_get_ATM_OCN

!-
!OCL SERIAL
  subroutine putATM_local( this, &
    lmesh_a, elem_a, &
!    ATM_TEMP, ATM_PRES, ATM_W, ATM_U, ATM_V, &
    ATM_RD_SFLX_SW_DIR, ATM_RD_SFLX_LW_DIF, ATM_SFLX_ENGI, &
!    O_ATM_TEMP, O_ATM_PRES, O_ATM_W, O_ATM_U, O_ATM_V, &
    O_ATM_RD_SFLX_SW_DIR, O_ATM_RD_SFLX_LW_DIF, O_ATM_SFLX_ENGI )
    implicit none
    class(CouplerVars), intent(inout) :: this
    class(LocalMesh2D), intent(in) :: lmesh_a
    class(ElementBase2D), intent(in) :: elem_a
    ! real(RP), intent(in) :: ATM_TEMP(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(in) :: ATM_PRES(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(in) :: ATM_W(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(in) :: ATM_U(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(in) :: ATM_V(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_RD_SFLX_SW_DIR(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_RD_SFLX_LW_DIF(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: ATM_SFLX_ENGI(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(out) :: O_ATM_TEMP(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(out) :: O_ATM_PRES(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(out) :: O_ATM_W(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(out) :: O_ATM_U(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(out) :: O_ATM_V(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: O_ATM_RD_SFLX_SW_DIR(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: O_ATM_RD_SFLX_LW_DIF(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: O_ATM_SFLX_ENGI(elem_a%Np,lmesh_a%NeA)

    integer :: ke
    integer :: coef1, coef2
    !----------------------------------------------------------------

    coef1 = this%CNT_putATM_OCN
    coef2 = 1.0_RP / ( 1.0_RP + coef1 )

    !$omp parallel
    !$omp do
    do ke=lmesh_a%NeS, lmesh_a%NeE
      ! O_ATM_TEMP(:,ke) = O_ATM_TEMP(:,ke) + coef1 * ATM_TEMP(:,ke)
      ! O_ATM_PRES(:,ke) = O_ATM_PRES(:,ke) + coef1 * ATM_PRES(:,ke)
      ! O_ATM_W(:,ke)    = O_ATM_W(:,ke)    + coef1 * ATM_W(:,ke)
      ! O_ATM_U(:,ke)    = O_ATM_U(:,ke)    + coef1 * ATM_U(:,ke)
      ! O_ATM_V(:,ke)    = O_ATM_V(:,ke)    + coef1 * ATM_V(:,ke)

      O_ATM_RD_SFLX_SW_DIR(:,ke) = O_ATM_RD_SFLX_SW_DIR(:,ke) + coef1 * ATM_RD_SFLX_SW_DIR(:,ke)
      O_ATM_RD_SFLX_LW_DIF(:,ke) = O_ATM_RD_SFLX_LW_DIF(:,ke) + coef1 * ATM_RD_SFLX_LW_DIF(:,ke)
      O_ATM_SFLX_ENGI(:,ke) = O_ATM_SFLX_ENGI(:,ke) + coef1 * ATM_SFLX_ENGI(:,ke)
    end do
    !$omp do
    do ke=lmesh_a%NeS, lmesh_a%NeE
      ! O_ATM_TEMP(:,ke) = O_ATM_TEMP(:,ke) * coef2
      ! O_ATM_PRES(:,ke) = O_ATM_PRES(:,ke) * coef2
      ! O_ATM_W(:,ke)    = O_ATM_W(:,ke)    * coef2
      ! O_ATM_U(:,ke)    = O_ATM_U(:,ke)    * coef2
      ! O_ATM_V(:,ke)    = O_ATM_V(:,ke)    * coef2

      O_ATM_RD_SFLX_SW_DIR(:,ke) = O_ATM_RD_SFLX_SW_DIR(:,ke) * coef2
      O_ATM_RD_SFLX_LW_DIF(:,ke) = O_ATM_RD_SFLX_LW_DIF(:,ke) * coef2
      O_ATM_SFLX_ENGI(:,ke) = O_ATM_SFLX_ENGI(:,ke) * coef2
    end do
    !$omp end parallel
    return
  end subroutine putATM_local

!OCL SERIAL
  subroutine putOCN_local( this, &
    lmesh_o, elem_o, &
    SFC_TEMP,  &
    SFC_TEMP_out  )
    implicit none
    class(CouplerVars), intent(inout) :: this
    class(LocalMesh2D), intent(in) :: lmesh_o
    class(ElementBase2D), intent(in) :: elem_o
    real(RP), intent(in) :: SFC_TEMP(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: SFC_TEMP_out(elem_o%Np,lmesh_o%NeA)

    integer :: ke
    !--------------------------------------------------------------------
    !$omp parallel do
    do ke=lmesh_o%NeS, lmesh_o%NeE
      SFC_TEMP_out(:,ke) = SFC_TEMP(:,ke)
    end do

    return
  end subroutine putOCN_local
  
!OCL SERIAL
  subroutine get_SFC_ATM_local( this, &
    lmesh_a, elem_a, lmesh_o, elem_o, &
    O_SFC_TEMP, O_SFC_ALBEDO, &
    SFC_TEMP, SFC_ALBEDO )
    implicit none
    class(CouplerVars), intent(inout) :: this
    class(LocalMesh2D), intent(in) :: lmesh_a
    class(ElementBase2D), intent(in) :: elem_a
    class(LocalMesh2D), intent(in) :: lmesh_o
    class(ElementBase2D), intent(in) :: elem_o
    real(RP), intent(in) :: O_SFC_TEMP(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(in) :: O_SFC_ALBEDO(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: SFC_TEMP(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(out) :: SFC_ALBEDO(elem_a%Np,lmesh_a%NeA)

    integer :: ke
    !----------------------------------------------------------------

    !$omp parallel do
    do ke=lmesh_a%NeS, lmesh_a%NeE
      SFC_TEMP(:,ke)  = O_SFC_TEMP(:,ke)
      SFC_ALBEDO(:,ke) = O_SFC_ALBEDO(:,ke)
    end do
    return
  end subroutine get_SFC_ATM_local

!OCL SERIAL
  subroutine get_ATM_OCN_local( this, &
    lmesh_a, elem_a, lmesh_o, elem_o, &
!    O_ATM_TEMP, O_ATM_PRES, O_ATM_W, O_ATM_U, O_ATM_V, O_ATM_SFLX_ENGI, &
    O_ATM_RD_SFLX_SW_DIR, O_ATM_RD_SFLX_LW_DIF, &
!    ATM_TEMP, ATM_PRES, ATM_W, ATM_U, ATM_V, ATM_SFLX_ENGI, &
    ATM_RD_SFLX_SW_DIR, ATM_RD_SFLX_LW_DIF )
    implicit none
    class(CouplerVars), intent(inout) :: this
    class(LocalMesh2D), intent(in) :: lmesh_a
    class(ElementBase2D), intent(in) :: elem_a
    class(LocalMesh2D), intent(in) :: lmesh_o
    class(ElementBase2D), intent(in) :: elem_o
    ! real(RP), intent(in) :: O_ATM_TEMP(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(in) :: O_ATM_PRES(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(in) :: O_ATM_W(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(in) :: O_ATM_U(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(in) :: O_ATM_V(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(in) :: O_ATM_SFLX_ENGI(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_RD_SFLX_SW_DIR(elem_a%Np,lmesh_a%NeA)
    real(RP), intent(in) :: O_ATM_RD_SFLX_LW_DIF(elem_a%Np,lmesh_a%NeA)
    ! real(RP), intent(out) :: ATM_TEMP(elem_o%Np,lmesh_o%NeA)
    ! real(RP), intent(out) :: ATM_PRES(elem_o%Np,lmesh_o%NeA)
    ! real(RP), intent(out) :: ATM_W(elem_o%Np,lmesh_o%NeA)
    ! real(RP), intent(out) :: ATM_U(elem_o%Np,lmesh_o%NeA)
    ! real(RP), intent(out) :: ATM_V(elem_o%Np,lmesh_o%NeA)
    ! real(RP), intent(out) :: ATM_SFLX_ENGI(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_RD_SFLX_SW_DIR(elem_o%Np,lmesh_o%NeA)
    real(RP), intent(out) :: ATM_RD_SFLX_LW_DIF(elem_o%Np,lmesh_o%NeA)

    integer :: ke
    !----------------------------------------------------------------

    !$omp parallel do
    do ke=lmesh_a%NeS, lmesh_a%NeE
      ! ATM_TEMP(:,ke)  = O_ATM_TEMP(:,ke)
      ! ATM_PRES(:,ke)  = O_ATM_PRES(:,ke)
      ! ATM_W(:,ke)     = O_ATM_W(:,ke)
      ! ATM_U(:,ke)     = O_ATM_U(:,ke)
      ! ATM_V(:,ke)     = O_ATM_V(:,ke)
      ! ATM_SFLX_ENGI(:,ke) = O_ATM_SFLX_ENGI(:,ke)
      ATM_RD_SFLX_SW_DIR(:,ke) = O_ATM_RD_SFLX_SW_DIR(:,ke)
      ATM_RD_SFLX_LW_DIF(:,ke) = O_ATM_RD_SFLX_LW_DIF(:,ke)
    end do
    return
  end subroutine get_ATM_OCN_local
end module mod_cpl_vars