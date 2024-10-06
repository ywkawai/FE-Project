!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Boundary
!!
!! @par Description
!!          A module for setting halo data based on boundary conditions. 
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_bnd
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

  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
  use scale_mesh_base, only: MeshBase
  use scale_localmesh_base, only: LocalMeshBase

  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_base2d, only: MeshBase2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmesh_2d, only: LocalMesh2D  

  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_mesh_bndinfo, only: MeshBndInfo

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    DENS_VID => PRGVAR_DDENS_ID, RHOT_VID => PRGVAR_DRHOT_ID, &
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID,                               &
    PRGVAR_NUM

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameter & type & procedures
  !
  type, public :: AtmDynBnd
    type(MeshBndInfo), allocatable :: VelBC_list(:)
    type(MeshBndInfo), allocatable :: ThermalBC_list(:)

    integer, allocatable:: velBC_ids(:)
    integer, allocatable :: thermalBC_ids(:)
    real(RP), allocatable :: thermal_fixval(:)
  contains
    procedure :: Init => ATMOS_dyn_bnd_setup
    procedure :: Final => ATMOS_dyn_bnd_finalize
    procedure :: SetBCInfo => ATMOS_dyn_bnd_setBCInfo
    procedure :: ApplyBC_PROGVARS_lc => ATMOS_dyn_bnd_applyBC_prgvars_lc
    procedure :: ApplyBC_numdiff_odd_lc => ATMOS_dyn_bnd_applyBC_numdiff_odd_lc
    procedure :: ApplyBC_numdiff_even_lc => ATMOS_dyn_bnd_applyBC_numdiff_even_lc
    procedure :: ApplyBC_Grad_TBVARS_lc => ATMOS_dyn_bnd_applyBC_tbvars_lc
    procedure :: ApplyBC_Grad_TBStress_lc => ATMOS_dyn_bnd_applyBC_tbstress_lc
    procedure :: Inquire_bound_flag => ATMOS_dyn_bnd_inquire_bound_flag
  end type

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------


  integer, parameter :: domBnd_South_ID = 1
  integer, parameter :: domBnd_East_ID  = 2
  integer, parameter :: domBnd_North_ID = 3
  integer, parameter :: domBnd_West_ID  = 4
  integer, parameter :: domBnd_Btm_ID   = 5
  integer, parameter :: domBnd_Top_ID   = 6
  integer, parameter :: DOM_BND_NUM     = 6

contains

!OCL SERIAL
  subroutine ATMOS_dyn_bnd_setup( this )
    use scale_const, only: &
      UNDEF8 => CONST_UNDEF8    
    use scale_mesh_bndinfo, only: &
      BND_TYPE_NOSPEC_NAME, &
      BND_TYPE_PERIODIC_ID, &
      BndType_NameToID
    implicit none
    
    class(AtmDynBnd), intent(inout) :: this

    character(len=H_SHORT) :: btm_vel_bc, top_vel_bc
    character(len=H_SHORT) :: north_vel_bc, south_vel_bc
    character(len=H_SHORT) :: east_vel_bc, west_vel_bc
    character(len=H_SHORT) :: btm_thermal_bc, top_thermal_bc
    character(len=H_SHORT) :: north_thermal_bc, south_thermal_bc
    character(len=H_SHORT) :: east_thermal_bc, west_thermal_bc

    real(RP) :: btm_thermal_fixval, top_thermal_fixval
    real(RP) :: north_thermal_fixval, south_thermal_fixval
    real(RP) :: east_thermal_fixval, west_thermal_fixval

    namelist /PARAM_ATMOS_DYN_BND/ &
      btm_vel_bc, top_vel_bc, north_vel_bc, south_vel_bc, east_vel_bc, west_vel_bc,                         &
      btm_thermal_bc, top_thermal_bc, north_thermal_bc, south_thermal_bc, east_thermal_bc, west_thermal_bc, &
      btm_thermal_fixval, top_thermal_fixval, north_thermal_fixval, south_thermal_fixval, east_thermal_fixval, west_thermal_fixval
    
    integer :: ierr
    !-----------------------------------------------

    btm_vel_bc   = BND_TYPE_NOSPEC_NAME
    top_vel_bc   = BND_TYPE_NOSPEC_NAME
    north_vel_bc  = BND_TYPE_NOSPEC_NAME
    south_vel_bc = BND_TYPE_NOSPEC_NAME
    east_vel_bc  = BND_TYPE_NOSPEC_NAME
    west_vel_bc = BND_TYPE_NOSPEC_NAME

    btm_thermal_bc   = BND_TYPE_NOSPEC_NAME
    top_thermal_bc   = BND_TYPE_NOSPEC_NAME
    north_thermal_bc  = BND_TYPE_NOSPEC_NAME
    south_thermal_bc = BND_TYPE_NOSPEC_NAME
    east_thermal_bc  = BND_TYPE_NOSPEC_NAME
    west_thermal_bc = BND_TYPE_NOSPEC_NAME

    btm_thermal_fixval = UNDEF8
    top_thermal_fixval = UNDEF8    
    north_thermal_fixval = UNDEF8
    south_thermal_fixval = UNDEF8    
    east_thermal_fixval = UNDEF8
    west_thermal_fixval = UNDEF8    

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_BND,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_dyn_bnd_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_dyn_bnd_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_BND. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN_BND)
    
    !--
    allocate( this%velBC_ids(DOM_BND_NUM) )
    this%velBC_ids(domBnd_Btm_ID) = BndType_NameToID(btm_vel_bc)
    this%velBC_ids(domBnd_Top_ID) = BndType_NameToID(top_vel_bc)
    this%velBC_ids(domBnd_North_ID) = BndType_NameToID(north_vel_bc)
    this%velBC_ids(domBnd_South_ID) = BndType_NameToID(south_vel_bc)
    this%velBC_ids(domBnd_East_ID) = BndType_NameToID(east_vel_bc)
    this%velBC_ids(domBnd_West_ID) = BndType_NameToID(west_vel_bc)

    allocate( this%thermalBC_ids(DOM_BND_NUM) )
    this%thermalBC_ids(domBnd_Btm_ID) = BndType_NameToID(btm_thermal_bc)
    this%thermalBC_ids(domBnd_Top_ID) = BndType_NameToID(top_thermal_bc)
    this%thermalBC_ids(domBnd_North_ID) = BndType_NameToID(north_thermal_bc)
    this%thermalBC_ids(domBnd_South_ID) = BndType_NameToID(south_thermal_bc)    
    this%thermalBC_ids(domBnd_East_ID) = BndType_NameToID(east_thermal_bc)
    this%thermalBC_ids(domBnd_West_ID) = BndType_NameToID(west_thermal_bc)   
    
    allocate( this%thermal_fixval(DOM_BND_NUM) )
    this%thermal_fixval(domBnd_Btm_ID) = btm_thermal_fixval
    this%thermal_fixval(domBnd_Top_ID) = top_thermal_fixval
    this%thermal_fixval(domBnd_North_ID) = north_thermal_fixval
    this%thermal_fixval(domBnd_South_ID) = south_thermal_fixval
    this%thermal_fixval(domBnd_East_ID) = east_thermal_fixval
    this%thermal_fixval(domBnd_West_ID) = west_thermal_fixval

    !------

    return
  end subroutine ATMOS_dyn_bnd_setup

!OCL SERIAL
  subroutine ATMOS_dyn_bnd_finalize( this )
    implicit none
    class(AtmDynBnd), intent(inout) :: this

    integer :: n
    !--------------------------------------

    if ( allocated(this%VelBC_list) ) then
      do n=1, size(this%VelBC_list)
        call this%VelBC_list(n)%Final()
        call this%ThermalBC_list(n)%Final()
      end do
      deallocate( this%VelBC_list, this%ThermalBC_list )

      deallocate( this%velBC_ids, this%thermalBC_ids )
    end if

    return
  end subroutine ATMOS_dyn_bnd_finalize  

!OCL SERIAL
  subroutine ATMOS_dyn_bnd_setBCInfo( this, mesh )

    implicit none

    class(AtmDynBnd), intent(inout) :: this
    class(MeshBase), target, intent(in) :: mesh
    
    integer :: n
    integer :: b
    class(LocalMeshBase), pointer :: ptr_lcmesh
    class(LocalMesh3D), pointer :: lcmesh3D
    !--------------------------------------------------


    allocate( this%VelBC_list(mesh%LOCAL_MESH_NUM) )
    allocate( this%ThermalBC_list(mesh%LOCAL_MESH_NUM) )

    nullify( lcmesh3D )

    do n=1, mesh%LOCAL_MESH_NUM
      call mesh%GetLocalMesh( n, ptr_lcmesh )
      select type (ptr_lcmesh)
      type is (LocalMesh3D)
        lcmesh3D => ptr_lcmesh
      end select

      call bnd_Init_lc( &
        this%VelBC_list(n), this%ThermalBC_list(n),        & ! (inout)
        this%velBC_ids(:), this%thermalBC_ids(:),          & ! (in)
        this%thermal_fixval(:),                            & ! (in)
        lcmesh3D%VMapB, mesh, lcmesh3D, lcmesh3D%refElem3D ) ! (in)
    end do

    return
  end subroutine ATMOS_dyn_bnd_setBCInfo

!OCL SERIAL
  subroutine ATMOS_dyn_bnd_applyBC_prgvars_lc( this,  &
    domID,                                              & ! (in)
    DDENS, MOMX, MOMY, MOMZ, THERM,                     & ! (inout)
    DENS_hyd, PRES_hyd,                                 & ! (in)
    Gsqrt, GsqrtH, G11, G12, G22, G13, G23, nx, ny, nz, & ! (in)
    vmapM, vmapP, vmapB, lmesh, elem, lmesh2D, elem2D   ) ! (in)

    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID
        
    implicit none

    class(AtmDynBnd), intent(in) :: this    
    integer, intent(in) :: domID
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(inout) :: DDENS(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMZ(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: THERM(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: GsqrtH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G11(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G12(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G22(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)    
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapB(:)

    integer :: p, ke, ke2D
    integer :: i, i_, iM, iP
    real(RP) :: mom_normal

    real(RP) :: MOMW
    real(RP) :: GsqrtV, G11_, G12_, G22_
    real(RP) :: fac
    !-----------------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke, p, ke2D, i, i_, iM, iP,      &
    !$omp mom_normal, MOMW, GsqrtV, G11_, G12_, G22_, fac    )
    do ke=lmesh%NeS, lmesh%NeE
    do p=1, elem%NfpTot
      i = p + (ke-1)*elem%NfpTot
      iP = vmapP(i)
      i_ = iP - elem%Np*lmesh%NeE
      
      if (i_ > 0) then
        iM = vmapM(i)

        select case( this%VelBC_list(domID)%list(i_) )
        case ( BND_TYPE_SLIP_ID)
          ke2D = lmesh%EMap3Dto2D(ke)

          GsqrtV = Gsqrt(iM) / GsqrtH(elem%IndexH2Dto3D_bnd(p),ke2D) 
          G11_ = G11(elem%IndexH2Dto3D_bnd(p),ke2D)
          G12_ = G12(elem%IndexH2Dto3D_bnd(p),ke2D)
          G22_ = G22(elem%IndexH2Dto3D_bnd(p),ke2D)

          MOMW = MOMZ(iM) / GsqrtV &
             + G13(iM) * MOMX(iM) + G23(iM) * MOMY(iM)
          fac = nz(i) * GsqrtV**2 / ( 1.0_RP + G11_ * ( GsqrtV * G13(iM) )**2 + 2.0_RP * G12_ * ( GsqrtV**2 * G13(iM) * G23(iM) ) + G22_ * ( GsqrtV * G23(iM) )**2 )

          mom_normal = MOMX(iM) * nx(i) + MOMY(iM) * ny(i) + MOMW * nz(i)
          MOMX(iP) = MOMX(iM) - 2.0_RP * mom_normal * ( nx(i) + fac * ( G11_ * G13(iM) + G12_ * G23(iM) ) )
          MOMY(iP) = MOMY(iM) - 2.0_RP * mom_normal * ( ny(i) + fac * ( G12_ * G13(iM) + G22_ * G23(iM) ) )
          MOMZ(iP) = MOMZ(iM) - 2.0_RP * mom_normal * fac / GsqrtV

        case ( BND_TYPE_NOSLIP_ID )
          MOMX(iP) = - MOMX(iM)
          MOMY(iP) = - MOMY(iM)
          MOMZ(iP) = - MOMZ(iM)          
        end select
      end if

    end do
    end do
    
    return
  end  subroutine ATMOS_dyn_bnd_applyBC_prgvars_lc
  
!OCL SERIAL
  subroutine ATMOS_dyn_bnd_applyBC_numdiff_odd_lc(  this, & ! (in)
    GxVar, GyVar, GzVar,                                  & ! (inout)
    is_bound,                                             & ! (out)
    VarID, domID,                                         & ! (in)
    nx, ny, nz, vmapM, vmapP, vmapB, lmesh, elem )          ! (in)

    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID
    
    implicit none

    class(AtmDynBnd), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem    
    real(RP), intent(inout) :: GxVar(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GyVar(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GzVar(elem%Np*lmesh%NeA)
    logical, intent(out) :: is_bound(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: VarID
    integer, intent(in) :: domID
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapB(:)

    integer :: i, i_, iM, iP
    real(RP) :: grad_normal
    !-----------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iP = vmapP(i)
      i_ = iP - elem%Np*lmesh%NeE
      is_bound(i) = .false.
      
      if (i_ > 0) then  
        iM = vmapM(i)
        grad_normal = GxVar(iM) * nx(i) + GyVar(iM) * ny(i) + GzVar(iM) * nz(i)

        if ( this%VelBC_list(domID)%list(i_) == BND_TYPE_SLIP_ID ) then
          select case(VarID)
          case(MOMX_VID)
            GyVar(iP) = GyVar(iM) - 2.0_RP * grad_normal * ny(i)
            GzVar(iP) = GzVar(iM) - 2.0_RP * grad_normal * nz(i)
          case(MOMY_VID)          
            GxVar(iP) = GxVar(iM) - 2.0_RP * grad_normal * nx(i)
            GzVar(iP) = GzVar(iM) - 2.0_RP * grad_normal * nz(i)
          case(MOMZ_VID)          
            GxVar(iP) = GxVar(iM) - 2.0_RP * grad_normal * nx(i)
            GyVar(iP) = GyVar(iM) - 2.0_RP * grad_normal * ny(i)
          end select
          is_bound(i) = .true.

        end if
        if ( this%ThermalBC_list(domID)%list(i_) == BND_TYPE_ADIABAT_ID ) then
          select case(VarID)
          case(DENS_VID, RHOT_VID)
            GxVar(iP) = GxVar(iM) - 2.0_RP * grad_normal * nx(i)
            GyVar(iP) = GyVar(iM) - 2.0_RP * grad_normal * ny(i)
            GzVar(iP) = GzVar(iM) - 2.0_RP * grad_normal * nz(i)
          end select
          is_bound(i) = .true.

        end if
      
      end if

    end do

    return
  end subroutine ATMOS_dyn_bnd_applyBC_numdiff_odd_lc

!OCL SERIAL
  subroutine ATMOS_dyn_bnd_applyBC_numdiff_even_lc(  this, & ! (in)
    Var,                                                   & ! (inout)
    is_bound,                                              & ! (out)
    VarID, domID,                                          & ! (in)
    nx, ny, nz, vmapM, vmapP, vmapB, lmesh, elem )           ! (in)

    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID
    implicit none

    class(AtmDynBnd), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem    
    real(RP), intent(inout) :: Var(elem%Np*lmesh%NeA)
    logical, intent(out) :: is_bound(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: VarID
    integer, intent(in) :: domID
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapB(:)

    integer :: i, i_, iM, iP
    real(RP) :: grad_normal
    !-----------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iP = vmapP(i)
      i_ = iP - elem%Np*lmesh%NeE
      is_bound(i) = .false.

      if (i_ > 0) then  
        iM = vmapM(i)

        if ( this%VelBC_list(domID)%list(i_) == BND_TYPE_SLIP_ID ) then
          select case(VarID)
          case(MOMX_VID)
            grad_normal = Var(iM) * nx(i)         
            Var(iP) = Var(iM) - 2.0_RP * grad_normal * nx(i)
          case(MOMY_VID) 
            grad_normal = Var(iM) * ny(i)         
            Var(iP) = Var(iM) - 2.0_RP * grad_normal * ny(i)
          case(MOMZ_VID)          
            grad_normal = Var(iM) * nz(i)         
            Var(iP) = Var(iM) - 2.0_RP * grad_normal * nz(i)
          end select
          is_bound(i) = .true.

        else if ( this%VelBC_list(domID)%list(i_) == BND_TYPE_NOSLIP_ID ) then          
          select case(VarID)
          case( MOMX_VID, MOMY_VID, MOMZ_VID )
            Var(iP) = - Var(iM)
          end select
          is_bound(i) = .true.
        end if

      end if
    end do

    return
  end subroutine ATMOS_dyn_bnd_applyBC_numdiff_even_lc


!> Set exterior values at boundaries for the turbulent schemes 
!!
!OCL SERIAL
  subroutine ATMOS_dyn_bnd_applyBC_tbvars_lc( this,  &
    domID,                                              & ! (in)
    DDENS, MOMX, MOMY, MOMZ, PT, PRES,                  & ! (inout)
    DENS_hyd, PRES_hyd, Rtot, CPtot,                    & ! (in)
    Gsqrt, GsqrtH, G11, G12, G22, G13, G23, nx, ny, nz, & ! (in)
    vmapM, vmapP, vmapB, lmesh, elem, lmesh2D, elem2D   ) ! (in)

    use scale_const, only: &
      PRES00 => CONST_PRE00   
    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID, BND_TYPE_FIXVAL_ID
        
    implicit none

    class(AtmDynBnd), intent(in) :: this    
    integer, intent(in) :: domID
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(inout) :: DDENS(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMZ(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: PT(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: PRES(elem%Np*lmesh%NeA)    
    real(RP), intent(in) :: DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: Rtot(elem%Np*lmesh%NeA)    
    real(RP), intent(in) :: CPtot(elem%Np*lmesh%NeA)    
    real(RP), intent(in) :: Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: GsqrtH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G11(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G12(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G22(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)  
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapB(:)

    integer :: p, ke, ke2D
    integer :: i, i_, iM, iP
    real(RP) :: mom_normal

    real(RP) :: MOMW
    real(RP) :: TEMP_P, TEMP_B
    real(RP) :: GsqrtV, G11_, G12_, G22_
    real(RP) :: fac
    !-----------------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke, p, ke2D, i, i_, iM, iP,      &
    !$omp mom_normal, MOMW, GsqrtV, G11_, G12_, G22_, fac, &
    !$omp TEMP_P, TEMP_B )
    do ke=lmesh%NeS, lmesh%NeE
    do p=1, elem%NfpTot
      i = p + (ke-1)*elem%NfpTot
      iP = vmapP(i)
      i_ = iP - elem%Np*lmesh%NeE
      
      if (i_ > 0) then
        iM = vmapM(i)

        select case( this%VelBC_list(domID)%list(i_) )
        case ( BND_TYPE_SLIP_ID)
          ke2D = lmesh%EMap3Dto2D(ke)

          GsqrtV = Gsqrt(iM) / GsqrtH(elem%IndexH2Dto3D_bnd(p),ke2D) 
          G11_ = G11(elem%IndexH2Dto3D_bnd(p),ke2D)
          G12_ = G12(elem%IndexH2Dto3D_bnd(p),ke2D)
          G22_ = G22(elem%IndexH2Dto3D_bnd(p),ke2D)

          MOMW = MOMZ(iM) / GsqrtV &
             + G13(iM) * MOMX(iM) + G23(iM) * MOMY(iM)
          fac = nz(i) * GsqrtV**2 / ( 1.0_RP + G11_ * ( GsqrtV * G13(iM) )**2 + 2.0_RP * G12_ * ( GsqrtV**2 * G13(iM) * G23(iM) ) + G22_ * ( GsqrtV * G23(iM) )**2 )

          mom_normal = MOMX(iM) * nx(i) + MOMY(iM) * ny(i) + MOMW * nz(i)
          MOMX(iP) = MOMX(iM) - 2.0_RP * mom_normal * ( nx(i) + fac * ( G11_ * G13(iM) + G12_ * G23(iM) ) )
          MOMY(iP) = MOMY(iM) - 2.0_RP * mom_normal * ( ny(i) + fac * ( G12_ * G13(iM) + G22_ * G23(iM) ) )
          MOMZ(iP) = MOMZ(iM) - 2.0_RP * mom_normal * fac / GsqrtV
        
        case ( BND_TYPE_NOSLIP_ID )
          MOMX(iP) = - MOMX(iM)
          MOMY(iP) = - MOMY(iM)
          MOMZ(iP) = - MOMZ(iM)          
        end select

        select case( this%ThermalBC_list(domID)%list(i_) )
        case ( BND_TYPE_FIXVAL_ID )
          TEMP_B =  this%ThermalBC_list(domID)%val(i_)
          TEMP_P = 2.0_RP * TEMP_B - PRES(iM) / ( ( DENS_hyd(iM) + DDENS(iM) ) * Rtot(iM) )
          PT(iP) = TEMP_P * ( PRES00 / PRES(iM) )**(Rtot(iM)/CPtot(iM))
          DDENS(iP) = PRES(iM) / ( Rtot(iM) * TEMP_P ) - DENS_hyd(iM)
        end select
      end if

    end do
    end do
    
    return
  end  subroutine ATMOS_dyn_bnd_applyBC_tbvars_lc

!> Set exterior shear-stress tensor and heat flux at boundaries for the turbulent schemes 
!!
!! Note: This codes are tentatively implementated. 
!! We need to check whether the formulations is valid for the case when general vertical coordinate is introduced. 
!OCL SERIAL
  subroutine ATMOS_dyn_bnd_applyBC_tbstress_lc( this,  &
    domID,                                              & ! (in)
    T11, T12, T13, T21, T22, T23, T31, T32, T33,        & ! (inout)
    DF1, DF2, DF3,                                      & ! (in)
    Gsqrt, GsqrtH, G11, G12, G22, G13, G23, nx, ny, nz, & ! (in)
    vmapM, vmapP, vmapB, lmesh, elem, lmesh2D, elem2D   ) ! (in)
 
    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID, BND_TYPE_FIXVAL_ID
        
    implicit none

    class(AtmDynBnd), intent(in) :: this    
    integer, intent(in) :: domID
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(inout) :: T11(elem%Np*lmesh%NeA), T12(elem%Np*lmesh%NeA), T13(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: T21(elem%Np*lmesh%NeA), T22(elem%Np*lmesh%NeA), T23(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: T31(elem%Np*lmesh%NeA), T32(elem%Np*lmesh%NeA), T33(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: DF1(elem%Np*lmesh%NeA), DF2(elem%Np*lmesh%NeA), DF3(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: GsqrtH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G11(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G12(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G22(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)    
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapB(:)

    integer :: p, ke, ke2D
    integer :: i, i_, iM, iP

    real(RP) :: GsqrtV, G11_, G12_, G22_
    real(RP) :: fac, normal_vec(3)

    real(RP) :: stress_t_tmp, stress_t_tmpvec(3)
    real(RP) :: VECW
    real(RP) :: normal_flux
    !-----------------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke, p, ke2D, i, i_, iM, iP,                      &
    !$omp VECW, GsqrtV, G11_, G12_, G22_, fac, normal_vec, &
    !$omp stress_t_tmp, stress_t_tmpvec, normal_flux )
    do ke=lmesh%NeS, lmesh%NeE
    do p=1, elem%NfpTot
      i = p + (ke-1)*elem%NfpTot
      iP = vmapP(i)
      i_ = iP - elem%Np*lmesh%NeE
      
      if (i_ > 0) then
        ke2D = lmesh%EMap3Dto2D(ke)
        iM = vmapM(i)

        GsqrtV = Gsqrt(iM) / GsqrtH(elem%IndexH2Dto3D_bnd(p),ke2D) 
        G11_ = G11(elem%IndexH2Dto3D_bnd(p),ke2D)
        G12_ = G12(elem%IndexH2Dto3D_bnd(p),ke2D)
        G22_ = G22(elem%IndexH2Dto3D_bnd(p),ke2D)
        fac = nz(i) * GsqrtV**2 / ( 1.0_RP + G11_ * ( GsqrtV * G13(iM) )**2 + 2.0_RP * G12_ * ( GsqrtV**2 * G13(iM) * G23(iM) ) + G22_ * ( GsqrtV * G23(iM) )**2 )
        normal_vec(1) = nx(i) + fac * ( G11_ * G13(iM) + G12_ * G23(iM) )
        normal_vec(2) = ny(i) + fac * ( G12_ * G13(iM) + G22_ * G23(iM) )
        normal_vec(3) = fac / GsqrtV

        select case( this%VelBC_list(domID)%list(i_) )
        case ( BND_TYPE_SLIP_ID)
          VECW = T13(iM) / GsqrtV + G13(iM) * T11(iM) + G23(iM) * T12(iM)
          stress_t_tmpvec(1) = T11(iM) * nx(i) + T12(iM) * ny(i) + VECW * nz(i)

          VECW = T23(iM) / GsqrtV + G13(iM) * T21(iM) + G23(iM) * T22(iM)
          stress_t_tmpvec(2) = T21(iM) * nx(i) + T22(iM) * ny(i) + VECW * nz(i)

          VECW = T33(iM) / GsqrtV
          stress_t_tmpvec(3) = T31(iM) * nx(i) + T32(iM) * ny(i) + VECW * nz(i)

          stress_t_tmp = sum( stress_t_tmpvec(:) * normal_vec(:) )
          stress_t_tmpvec(:) = stress_t_tmpvec(:) - stress_t_tmp * normal_vec(:) 

          ! T1j
          T11(iP) = T11(iM) - 2.0_RP * stress_t_tmpvec(1) * normal_vec(1)
          T12(iP) = T12(iM) - 2.0_RP * stress_t_tmpvec(1) * normal_vec(2)
          T13(iP) = T13(iM) - 2.0_RP * stress_t_tmpvec(1) * normal_vec(3)

          ! T2j
          T21(iP) = T21(iM) - 2.0_RP * stress_t_tmpvec(2) * normal_vec(1)
          T22(iP) = T22(iM) - 2.0_RP * stress_t_tmpvec(2) * normal_vec(2)
          T23(iP) = T23(iM) - 2.0_RP * stress_t_tmpvec(2) * normal_vec(3)

          ! T3j
          T31(iP) = T31(iM) - 2.0_RP * stress_t_tmpvec(3) * normal_vec(1)
          T32(iP) = T32(iM) - 2.0_RP * stress_t_tmpvec(3) * normal_vec(2)
          T33(iP) = T33(iM) - 2.0_RP * stress_t_tmpvec(3) * normal_vec(3)

        case ( BND_TYPE_NOSLIP_ID )
        end select

        select case( this%ThermalBC_list(domID)%list(i_) )
        case ( BND_TYPE_ADIABAT_ID )
          normal_flux = DF1(iM) * nx(i) + DF2(iM) * ny(i) &
            + ( DF3(iM) / GsqrtV + G13(iM) * DF1(iM) + G23(iM) * DF2(iM) )
          DF1(iP) = DF1(iM) - 2.0_RP * normal_flux * normal_vec(1)
          DF2(iP) = DF2(iM) - 2.0_RP * normal_flux * normal_vec(2)
          DF3(iP) = DF3(iM) - 2.0_RP * normal_flux * normal_vec(3)
        end select
      end if
    end do
    end do
    
    return
  end  subroutine ATMOS_dyn_bnd_applyBC_tbstress_lc

!OCL SERIAL
  subroutine ATMOS_dyn_bnd_inquire_bound_flag(  this,      & ! (in)
    is_bound,                                              & ! (out)
    domID, vmapM, vmapP, vmapB, lmesh, elem                ) ! (in)

    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID
    implicit none

    class(AtmDynBnd), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem    
    logical, intent(out) :: is_bound(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: domID
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapB(:)

    integer :: i, i_, iM, iP
    !-----------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iP = vmapP(i)
      i_ = iP - elem%Np*lmesh%NeE
      is_bound(i) = .false.

      if (i_ > 0) then  
        iM = vmapM(i)
        if ( this%VelBC_list(domID)%list(i_) == BND_TYPE_SLIP_ID ) then
          is_bound(i) = .true.
        else if ( this%VelBC_list(domID)%list(i_) == BND_TYPE_NOSLIP_ID ) then          
          is_bound(i) = .true.
        end if

      end if
    end do

    return
  end subroutine ATMOS_dyn_bnd_inquire_bound_flag

  !-- private -------------------------------------------------------

!OCL SERIAL
  subroutine bnd_Init_lc(   &
    velBCInfo, thermalBCInfo,           &  ! (inout)
    velBC_ids, thermalBC_ids,           &  ! (in)
    thermal_fixval,                     &  ! (in)
    vmapB, mesh, lmesh, elem )             ! (in)

    use scale_mesh_bndinfo, only: BND_TYPE_NOSPEC_ID
    implicit none
    
    type(MeshBndInfo), intent(inout) :: velBCInfo
    type(MeshBndInfo), intent(inout) :: thermalBCInfo
    integer, intent(in) :: velBC_ids(DOM_BND_NUM)
    integer, intent(in) :: thermalBC_ids(DOM_BND_NUM)
    real(RP), intent(in) :: thermal_fixval(DOM_BND_NUM)
    class(MeshBase), intent(in) :: mesh
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: vmapB(:)

    integer :: tileID
    integer :: dom_bnd_sizes(DOM_BND_NUM)
    integer :: bnd_buf_size
    integer :: b, is_, ie_

    !-----------------------------------------------

    dom_bnd_sizes(:) = &
        elem%Nfp_h*lmesh%NeZ*(/ lmesh%NeX, lmesh%NeY, lmesh%NeX, lmesh%NeY, 0, 0 /) &
      + elem%Nfp_v*lmesh%NeX*lmesh%NeY*(/ 0, 0, 0, 0, 1, 1 /)
    bnd_buf_size = sum(dom_bnd_sizes)
    
    call velBCInfo%Init( bnd_buf_size )
    call velBCInfo%Set(1, bnd_buf_size, BND_TYPE_NOSPEC_ID)

    call thermalBCInfo%Init( bnd_buf_size )
    call thermalBCInfo%Set(1, bnd_buf_size, BND_TYPE_NOSPEC_ID)

    tileID = lmesh%tileID
    is_ = 1
    do b=1, DOM_BND_NUM
      ie_ = is_ + dom_bnd_sizes(b) - 1
      if ( mesh%tileID_globalMap(b,tileID) == tileID      &
           .and. mesh%tileFaceID_globalMap(b,tileID) == b ) then
        call velBCInfo%Set( is_, ie_, velbc_ids(b) )
        call thermalBCInfo%Set( is_, ie_, thermalbc_ids(b), thermal_fixval(b) )
      end if
      is_ = ie_ + 1
    end do

    return
  end  subroutine bnd_Init_lc

end module scale_atm_dyn_dgm_bnd