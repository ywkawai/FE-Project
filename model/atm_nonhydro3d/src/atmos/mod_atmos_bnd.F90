!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_bnd
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

  use scale_sparsemat  
  use scale_element_base
  use scale_element_hexahedral
  use scale_localmesh_3d
  use scale_mesh_cubedom3d

  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_mesh_bndinfo, only: &
    MeshBndInfo

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: ATMOS_bnd_setup
  public :: ATMOS_bnd_finalize
  public :: ATMOS_bnd_setBCInfo
  public :: ATMOS_bnd_applyBC_prgvars
  public :: ATMOS_bnd_applyBC_auxvars

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  type(MeshBndInfo), allocatable :: VelBC_list(:)
  type(MeshBndInfo), allocatable :: ThermalBC_list(:)

  integer, parameter :: domBnd_North_ID = 1
  integer, parameter :: domBnd_South_ID = 2
  integer, parameter :: domBnd_East_ID  = 3
  integer, parameter :: domBnd_West_ID  = 4
  integer, parameter :: domBnd_Btm_ID   = 5
  integer, parameter :: domBnd_Top_ID   = 6
  integer, parameter :: DOM_BND_NUM     = 6

  integer :: velBC_ids(DOM_BND_NUM)
  integer :: thermalBC_ids(DOM_BND_NUM)

  logical, save :: is_initialized = .false.

contains
  subroutine ATMOS_bnd_setup()
    use scale_mesh_bndinfo, only: &
      BND_TYPE_NOSPEC_NAME, &
      BND_TYPE_PERIODIC_ID, &
      BndType_NameToID
    implicit none
    
    character(len=H_SHORT) :: btm_vel_bc, top_vel_bc
    character(len=H_SHORT) :: north_vel_bc, south_vel_bc
    character(len=H_SHORT) :: east_vel_bc, west_vel_bc
    character(len=H_SHORT) :: btm_thermal_bc, top_thermal_bc
    character(len=H_SHORT) :: north_thermal_bc, south_thermal_bc
    character(len=H_SHORT) :: east_thermal_bc, west_thermal_bc

    namelist /PARAM_ATMOS_BND/ &
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
    read(IO_FID_CONF,nml=PARAM_ATMOS_BND,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_bnd_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_bnd_setup",*) 'Not appropriate names in namelist PARAM_BND. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_BND)
    
    velBC_ids(domBnd_Btm_ID) = BndType_NameToID(btm_vel_bc)
    velBC_ids(domBnd_Top_ID) = BndType_NameToID(top_vel_bc)
    velBC_ids(domBnd_North_ID) = BndType_NameToID(north_vel_bc)
    velBC_ids(domBnd_South_ID) = BndType_NameToID(south_vel_bc)
    velBC_ids(domBnd_East_ID) = BndType_NameToID(east_vel_bc)
    velBC_ids(domBnd_West_ID) = BndType_NameToID(west_vel_bc)

    thermalBC_ids(domBnd_Btm_ID) = BndType_NameToID(btm_thermal_bc)
    thermalBC_ids(domBnd_Top_ID) = BndType_NameToID(top_thermal_bc)
    thermalBC_ids(domBnd_North_ID) = BndType_NameToID(north_thermal_bc)
    thermalBC_ids(domBnd_South_ID) = BndType_NameToID(south_thermal_bc)    
    thermalBC_ids(domBnd_East_ID) = BndType_NameToID(east_thermal_bc)
    thermalBC_ids(domBnd_West_ID) = BndType_NameToID(west_thermal_bc)    
    !------

    is_initialized = .false.
    return
  end subroutine ATMOS_bnd_setup

  subroutine ATMOS_bnd_finalize()
    implicit none

    integer :: n
    !--------------------------------------

    if (is_initialized) then
      do n=1, size(VelBC_list)
        call VelBC_list(n)%Final()
        call ThermalBC_list(n)%Final()
      end do
      deallocate( VelBC_list, ThermalBC_list )
    end if

    is_initialized = .false.
    return
  end subroutine ATMOS_bnd_finalize  

  subroutine ATMOS_bnd_setBCInfo()
    use mod_atmos_mesh, only: mesh
    implicit none

    integer :: n
    integer :: b
    type(LocalMesh3D), pointer :: lcmesh  
    !--------------------------------------------------


    allocate( VelBC_list(mesh%LOCAL_MESH_NUM) )
    allocate( ThermalBC_list(mesh%LOCAL_MESH_NUM) )

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call bnd_Init_lc( &
        VelBC_list(n), ThermalBC_list(n),                       & ! (inout)
        lcmesh%VMapB, mesh,  lcmesh, lcmesh%refElem3D )           ! (in)
    end do

    is_initialized = .true.
    return
  end subroutine ATMOS_bnd_setBCInfo

  subroutine ATMOS_bnd_applyBC_prgvars(    &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,        &
    DENS_hydro, PRES_hydro,                &
    mesh )

    implicit none
    type(MeshCubeDom3D), intent(in), target :: mesh
    type(MeshField3D), intent(inout) :: DDENS
    type(MeshField3D), intent(inout) :: MOMX
    type(MeshField3D), intent(inout) :: MOMY
    type(MeshField3D), intent(inout) :: MOMZ
    type(MeshField3D), intent(inout) :: DRHOT
    type(MeshField3D), intent(in) :: DENS_hydro
    type(MeshField3D), intent(in) :: PRES_hydro

    integer :: n
    type(LocalMesh3D), pointer :: lcmesh
    !-----------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call applyBC_prgvars_lc( &
        DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val, & ! (inout)
        DENS_hydro%local(n)%val, PRES_hydro%local(n)%val,                                                & ! (in)
        VelBC_list(n), ThermalBC_list(n),                                                                & ! (in)
        lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3),                       & ! (in)
        lcmesh%VMapM, lcmesh%VMapP, lcmesh%VMapB,                                                        & ! (in)
        lcmesh, lcmesh%refElem3D )                                                                         ! (in)
    end do

    return
  end subroutine ATMOS_bnd_applyBC_prgvars

  subroutine ATMOS_bnd_applyBC_auxvars(  &
    GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW,    &
    GxPT, GyPT, GzPT,                               &
    DENS_hydro, PRES_hydro,                         &
    viscCoef_h, viscCoef_v, diffCoef_h, diffCoef_v, &
    mesh )

    implicit none
    type(MeshCubeDom3D), intent(in), target :: mesh
    type(MeshField3D), intent(inout) :: GxU, GyU, GzU
    type(MeshField3D), intent(inout) :: GxV, GyV, GzV
    type(MeshField3D), intent(inout) :: GxW, GyW, GzW
    type(MeshField3D), intent(inout) :: GxPT, GyPT, GzPT
    type(MeshField3D), intent(in) :: DENS_hydro
    type(MeshField3D), intent(in) :: PRES_hydro
    real(RP), intent(in) :: viscCoef_h
    real(RP), intent(in) :: viscCoef_v
    real(RP), intent(in) :: diffCoef_h
    real(RP), intent(in) :: diffCoef_v

    integer :: n
    type(LocalMesh3D), pointer :: lcmesh
    !-----------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call applyBC_auxvars_lc( &
        GxU%local(n)%val,  GyU%local(n)%val,  GzU%local(n)%val,                    & ! (inout)
        GxV%local(n)%val,  GyV%local(n)%val,  GzV%local(n)%val,                    & ! (inout)
        GxW%local(n)%val,  GyW%local(n)%val,  GzW%local(n)%val,                    & ! (inout)
        GxPT%local(n)%val, GyPT%local(n)%val, GzPT%local(n)%val,                   & ! (inour) 
        DENS_hydro%local(n)%val, PRES_hydro%local(n)%val,                          & ! (in)
        viscCoef_h, viscCoef_v, diffCoef_h, diffCoef_v,                            & ! (in)
        VelBC_list(n), ThermalBC_list(n),                                          & ! (in)
        lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), & ! (in)
        lcmesh%VMapM, lcmesh%VMapP, lcmesh%VMapB,                                  & ! (in)
        lcmesh, lcmesh%refElem3D )                                                   ! (in)
    end do

    return
  end subroutine ATMOS_bnd_applyBC_auxvars
    
  !---------------------------------------------------------

  subroutine bnd_Init_lc(   &
    velBCInfo, thermalBCInfo,           &  ! (inout)
    vmapB, mesh, lmesh, elem )             ! (in)

    use scale_mesh_bndinfo, only: BND_TYPE_NOSPEC_ID
    implicit none
    
    type(MeshBndInfo), intent(inout) :: velBCInfo
    type(MeshBndInfo), intent(inout) :: thermalBCInfo
    type(MeshCubeDom3D), intent(in) :: mesh
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(in) :: vmapB(:)

    integer :: tileID
    integer :: dom_bnd_sizes(DOM_BND_NUM)
    integer :: bnd_buf_size
    integer :: b, is_, ie_

    !-----------------------------------------------

    dom_bnd_sizes(:) = &
        elem%NfpTot*lmesh%NeZ*(/ lmesh%NeX, lmesh%NeY, lmesh%NeX, lmesh%NeY, 0, 0 /) &
      + elem%NfpTot*lmesh%NeX*lmesh%NeY*(/ 0, 0, 0, 0, 1, 1 /)
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
  
  subroutine applyBC_prgvars_lc(  &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,               & ! (inout)
    DENS_hyd, PRES_hyd,                           & ! (in)
    velBCInfo, thermalBCInfo,                     & ! (in)
    nx, ny, nz, vmapM, vmapP, vmapB, lmesh, elem )  ! (in)

    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem    
    real(RP), intent(inout) :: DDENS(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: MOMZ(elem%Np*lmesh%NeA)
    real(RP), intent(inout) :: DRHOT(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np*lmesh%NeA)
    type(MeshBndInfo), intent(in) :: velBCInfo
    type(MeshBndInfo), intent(in) :: thermalBCInfo
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
        if (velBCInfo%list(i_) == BND_TYPE_SLIP_ID) then
          mom_normal = MOMX(iM)*nx(i) + MOMY(iM)*ny(i) + MOMZ(iM)*nz(i)
          MOMX(iP) = MOMX(iM) - 2.0_RP*mom_normal*nx(i)
          MOMY(iP) = MOMY(iM) - 2.0_RP*mom_normal*ny(i)
          MOMZ(iP) = MOMZ(iM) - 2.0_RP*mom_normal*nz(i)
        end if
        if (velBCInfo%list(i_) == BND_TYPE_NOSLIP_ID) then
          MOMX(iP) = - MOMX(iM)
          MOMY(iP) = - MOMY(iM)
          MOMZ(iP) = - MOMZ(iM)
        end if
      end if
      
    end do
    
    return
  end subroutine applyBC_prgvars_lc

  subroutine applyBC_auxvars_lc(  &
    GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW,   & ! (inout)
    GxPT, GyPT, GzPT,                              & ! (inout)
    DENS_hyd, PRES_hyd,                            & ! (in)
    viscCoef_h, viscCoef_v,                        & ! (in)
    diffCoef_h, diffCoef_v,                        & ! (in)
    velBCInfo, thermalBCInfo,                      & ! (in)
    nx, ny, nz, vmapM, vmapP, vmapB, lmesh, elem )   ! (in)

    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID
    implicit none

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
    type(MeshBndInfo), intent(in) :: velBCInfo
    type(MeshBndInfo), intent(in) :: thermalBCInfo
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

        if (velBCInfo%list(i_) == BND_TYPE_SLIP_ID) then
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
        if (thermalBCInfo%list(i_) == BND_TYPE_ADIABAT_ID) then
          gradPT_normal = GxPT(iM)*nx(i) + GyPT(iM)*ny(i) + GzPT(iM)*nz(i)
          GxPT(iP) = GxPT(iM) - 2.0_RP*gradPT_normal*nx(i)
          GyPT(iP) = GyPT(iM) - 2.0_RP*gradPT_normal*ny(i)
          GzPT(iP) = GzPT(iM) - 2.0_RP*gradPT_normal*nz(i)
        end if

      end if
    end do

    return
  end subroutine applyBC_auxvars_lc

end module mod_atmos_bnd