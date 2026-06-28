!> module FElib / Atmosphere / Physics radiation
!!
!! @par Description
!!      Radiation process with a gray-radiation scheme
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_rd_dgm_gray
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  type, public :: AtmPhyRadGray
    integer :: optdep_type

    real(RP) :: OPTDEP_VALLIS_EQ5_PARAM_A
    real(RP) :: OPTDEP_VALLIS_EQ5_PARAM_MU
    real(RP) :: OPTDEP_VALLIS_EQ5_PARAM_C
  contains
  end type AtmPhyRadGray

  integer, parameter :: OPTDEP_TYPE_VALLIS2018_EQ5 = 1
  
contains
!OCL SERIAL
  subroutine atm_phy_rd_dgm_gray_Init( mesh )
    implicit none    
    class(MeshBase3D), intent(in) :: mesh

    character(len=H_MID) :: OPTICAL_DEPTH_TYPE
    namelist / PARAM_ATMOS_PHY_RD_DGM_GRAY / &
      OPTICAL_DEPTH_TYPE
        
    integer :: ierr
    !--------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_dgm_gray_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_DGM_GRAY,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_RD_dgm_gray_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_RD_dgm_gray_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_RD_DGM_GRAY. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_RD_DGM_GRAY)

    select case(trim(OPTICAL_DEPTH_TYPE))
    case('Vallis2018_Eq5')
      this%optdep_type = OPTDEP_TYPE_VALLIS2018_EQ5
    case default
      LOG_ERROR("ATMOS_PHY_RD_dgm_gray_setup",*) 'OPTICAL_DEPTH_TYPE is invalid. Check!', trim(OPTICAL_DEPTH_TYPE)
    end select
    return
  end subroutine atm_phy_rd_dgm_gray_Init

!OCL SERIAL
  subroutine atm_phy_rd_dgm_gray_Final()
    implicit none
    !--------------------------------------------------------------------
    return
  end subroutine atm_phy_rd_dgm_gray_Final

!OCL SERIAL
  subroutine atm_phy_rd_dgm_gray_flux( flux, &
    DENS, TEMP, &
    lcmesh, elem3D, lcmesh2D, elem2D, elem_v1D )
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    class(ElementBase1D), intent(in) :: elem_v1D
    implicit none
    real(RP), intent(in) :: flux(elem3D%Nnode_v,lcmesh%NeZ,,elem3D%Nnode_h1D**2,lcmesh%Ne2D,2,2,2)
    real(RP), intent(in) :: TEMP(elem3D%Nnode_v,lcmesh%NeZ,,elem3D%Nnode_h1D**2,lcmesh%Ne2D)
    real(RP), intent(in) :: DENS(elem3D%Nnode_v,lcmesh%NeZ,,elem3D%Nnode_h1D**2,lcmesh%Ne2D)

    integer :: ke_z, ke_h
    integer :: p_z, p_h

    real(RP) :: dtau(elem3D%Nnode_v,lcmesh%NeZ)
    real(RP) :: TransFunc(elem3D%Nnode_v,lcmesh%NeZ)

    !--------------------------------------------------------------------
    !$omp parallel do
    do ke_h=1, lcmesh%Ne2D
    do p_h=1, elem3D%Nnode_h1D**2
      do ke_z=1, lcmesh%NeZ
      do p_z=1, elem3D%Nnode_v
        dtau(p_z,ke_z) = 
      end do
      end do
    end do
    end do
    return
  end subroutine atm_phy_tb_dgm_gray_flux

!-- private -----
  subroutine calc_optical_thick( dtau, &
    pres, qv, Nnode_v, NeZ )
    implicit none
    integer, intent(in) :: Nnode_v
    integer, intent(in) :: NeZ
    real(RP), intent(out) :: dtau(Nnode_v-1,NeZ)
    real(RP), intent(in) :: pres(Nnode_v,NeZ)
    real(RP), intent(in) :: qv(Nnode_v,NeZ)

    integer :: ke_z, p_z
    !---------------------------------------------------------

    do ke_z=1, NeZ
    do p_z=1, Nnode_v-1
      dtau(p_z,ke_z) = 
    end do
    end do
    return
  end subroutine calc_optical_thick

end  module scale_atm_phy_rd_dgm_gray
