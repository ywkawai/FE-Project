!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_dyn_bnd
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
    ElementBase, ElementBase3D
  use scale_mesh_base, only: MeshBase
  use scale_localmesh_base, only: LocalMeshBase

  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_3d, only: LocalMesh3D

  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_mesh_bndinfo, only: MeshBndInfo

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameter & type & procedures
  !
  type, public :: AtmosDynBnd
    type(MeshBndInfo), allocatable :: VelBC_list(:)
    type(MeshBndInfo), allocatable :: ThermalBC_list(:)

    integer, allocatable:: velBC_ids(:)
    integer, allocatable :: thermalBC_ids(:)
  contains
    procedure :: Init => ATMOS_dyn_bnd_setup
    procedure :: Final => ATMOS_dyn_bnd_finalize
    procedure :: SetBCInfo => ATMOS_dyn_bnd_setBCInfo
    procedure :: ApplyBC_PROGVARS_lc => ATMOS_dyn_bnd_applyBC_prgvars_lc
    procedure :: ApplyBC_AUXVARS_lc => ATMOS_dyn_bnd_applyBC_auxvars_lc
  end type



  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------



  integer, parameter :: domBnd_North_ID = 1
  integer, parameter :: domBnd_South_ID = 2
  integer, parameter :: domBnd_East_ID  = 3
  integer, parameter :: domBnd_West_ID  = 4
  integer, parameter :: domBnd_Btm_ID   = 5
  integer, parameter :: domBnd_Top_ID   = 6
  integer, parameter :: DOM_BND_NUM     = 6


contains
  subroutine ATMOS_dyn_bnd_setup( this )
    use scale_mesh_bndinfo, only: &
      BND_TYPE_NOSPEC_NAME, &
      BND_TYPE_PERIODIC_ID, &
      BndType_NameToID
    implicit none
    
    class(AtmosDynBnd), intent(inout) :: this

    character(len=H_SHORT) :: btm_vel_bc, top_vel_bc
    character(len=H_SHORT) :: north_vel_bc, south_vel_bc
    character(len=H_SHORT) :: east_vel_bc, west_vel_bc
    character(len=H_SHORT) :: btm_thermal_bc, top_thermal_bc
    character(len=H_SHORT) :: north_thermal_bc, south_thermal_bc
    character(len=H_SHORT) :: east_thermal_bc, west_thermal_bc

    namelist /PARAM_ATMOS_DYN_BND/ &
      btm_vel_bc, top_vel_bc, north_vel_bc, south_vel_bc, east_vel_bc, west_vel_bc,                        &
      btm_thermal_bc, top_thermal_bc, north_thermal_bc, south_thermal_bc, east_thermal_bc, west_thermal_bc
    
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

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_BND,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_dyn_bnd_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_dyn_bnd_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_BND. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN_BND)
    
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
    !------

    return
  end subroutine ATMOS_dyn_bnd_setup

  subroutine ATMOS_dyn_bnd_finalize( this )
    implicit none
    class(AtmosDynBnd), intent(inout) :: this

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

  subroutine ATMOS_dyn_bnd_setBCInfo( this, mesh )

    implicit none

    class(AtmosDynBnd), intent(inout) :: this
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
        this%VelBC_list(n), this%ThermalBC_list(n),       & ! (inout)
        this%velBC_ids(:), this%thermalBC_ids(:),         & ! (in)
        lcmesh3D%VMapB, mesh, lcmesh3D, lcmesh3D%refElem3D )     ! (in)
    end do

    return
  end subroutine ATMOS_dyn_bnd_setBCInfo

  subroutine ATMOS_dyn_bnd_applyBC_prgvars_lc( this, domID,  &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,               & ! (inout)
    DENS_hyd, PRES_hyd,                           & ! (in)
    nx, ny, nz, vmapM, vmapP, vmapB, lmesh, elem )  ! (in)

    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID
    implicit none

    class(AtmosDynBnd), intent(in) :: this    
    integer, intent(in) :: domID
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem    
    real(RP), intent(inout) :: DDENS(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMZ(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: DRHOT(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapB(:)

    integer :: i, i_, iM, iP
    real(RP) :: mom_normal
    !-----------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iP = vmapP(i)
      i_ = iP - elem%Np*lmesh%NeE
      
      if (i_ > 0) then
        iM = vmapM(i)
        select case( this%VelBC_list(domID)%list(i_) )
        case ( BND_TYPE_SLIP_ID)
          mom_normal = MOMX(iM)*nx(i) + MOMY(iM)*ny(i) + MOMZ(iM)*nz(i)
          MOMX(iP) = MOMX(iM) - 2.0_RP*mom_normal*nx(i)
          MOMY(iP) = MOMY(iM) - 2.0_RP*mom_normal*ny(i)
          MOMZ(iP) = MOMZ(iM) - 2.0_RP*mom_normal*nz(i)
        case ( BND_TYPE_NOSLIP_ID )
          MOMX(iP) = - MOMX(iM)
          MOMY(iP) = - MOMY(iM)
          MOMZ(iP) = - MOMZ(iM)          
        end select
      end if
    end do
    
    return
  end subroutine ATMOS_dyn_bnd_applyBC_prgvars_lc

  subroutine ATMOS_dyn_bnd_applyBC_auxvars_lc(  this,  domID, &
    GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW,   & ! (inout)
    GxPT, GyPT, GzPT,                              & ! (inout)
    DENS_hyd, PRES_hyd,                            & ! (in)
    viscCoef_h, viscCoef_v,                        & ! (in)
    diffCoef_h, diffCoef_v,                        & ! (in)
    nx, ny, nz, vmapM, vmapP, vmapB, lmesh, elem )   ! (in)

    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID
    implicit none

    class(AtmosDynBnd), intent(in) :: this
    integer, intent(in) :: domID
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem    
    real(RP), intent(inout) :: GxU(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GyU(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GzU(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GxV(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GyV(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GzV(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GxW(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GyW(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GzW(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GxPT(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GyPT(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: GzPT(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: viscCoef_h
    real(RP), intent(in) :: viscCoef_v
    real(RP), intent(in) :: diffCoef_h
    real(RP), intent(in) :: diffCoef_v
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapB(:)

    integer :: i, i_, iM, iP
    real(RP) :: gradU_normal
    real(RP) :: gradV_normal
    real(RP) :: gradW_normal
    real(RP) :: gradPT_normal

    !-----------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iP = vmapP(i)
      i_ = iP - elem%Np*lmesh%NeE
      
      if (i_ > 0) then  
        iM = vmapM(i)

        if ( this%VelBC_list(domID)%list(i_) == BND_TYPE_SLIP_ID ) then
          gradU_normal = GxU(iM)*nx(i) + GyU(iM)*ny(i) + GzU(iM)*nz(i)
          !GxU(iP) = GxU(iM) - 2.0_RP*gradU_normal*nx(i)
          GyU(iP) = GyU(iM) - 2.0_RP*gradU_normal*ny(i)
          GzU(iP) = GzU(iM) - 2.0_RP*gradU_normal*nz(i)

          gradV_normal = GxV(iM)*nx(i) + GyV(iM)*ny(i) + GzV(iM)*nz(i)
          !GxU(iP) = GxU(iM) - 2.0_RP*gradU_normal*nx(i)
          GxV(iP) = GxV(iM) - 2.0_RP*gradV_normal*nx(i)
          GzV(iP) = GzV(iM) - 2.0_RP*gradV_normal*nz(i)

          gradW_normal = GxW(iM)*nx(i) + GyW(iM)*ny(i) + GzW(iM)*nz(i)
          GxW(iP) = GxW(iM) - 2.0_RP*gradW_normal*nx(i)
          GyW(iP) = GyW(iM) - 2.0_RP*gradW_normal*ny(i)
          !GzW(iP) = GzW(iM) - 2.0_RP*gradW_normal*nz(i)
        end if
        if ( this%ThermalBC_list(domID)%list(i_) == BND_TYPE_ADIABAT_ID ) then
          gradPT_normal = GxPT(iM)*nx(i) + GyPT(iM)*ny(i) + GzPT(iM)*nz(i)
          GxPT(iP) = GxPT(iM) - 2.0_RP*gradPT_normal*nx(i)
          GyPT(iP) = GyPT(iM) - 2.0_RP*gradPT_normal*ny(i)
          GzPT(iP) = GzPT(iM) - 2.0_RP*gradPT_normal*nz(i)
        end if

      end if
    end do

    return
  end subroutine ATMOS_dyn_bnd_applyBC_auxvars_lc

  !---------------------------------------------------------

  subroutine bnd_Init_lc(   &
    velBCInfo, thermalBCInfo,           &  ! (inout)
    velBC_ids, thermalBC_ids,           &  ! (in)
    vmapB, mesh, lmesh, elem )             ! (in)

    use scale_mesh_bndinfo, only: BND_TYPE_NOSPEC_ID
    implicit none
    
    type(MeshBndInfo), intent(inout) :: velBCInfo
    type(MeshBndInfo), intent(inout) :: thermalBCInfo
    integer, intent(in) :: velBC_ids(DOM_BND_NUM)
    integer, intent(in) :: thermalBC_ids(DOM_BND_NUM)
    class(MeshBase), intent(in) :: mesh
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
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
        call velBCInfo%Set(is_, ie_, velbc_ids(b))
        call thermalBCInfo%Set(is_, ie_, thermalbc_ids(b))
      end if
      is_ = ie_ + 1
    end do

    return
  end  subroutine bnd_Init_lc
end module mod_atmos_dyn_bnd