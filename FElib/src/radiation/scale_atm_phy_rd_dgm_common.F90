!> module FElib / Atmosphere / Physics radiation
!!
!! @par Description
!!      Common module for radiation process
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_rd_dgm_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_atmos_phy_rd_common, only: &
    I_up, I_dn, I_LW, I_SW
  
  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_element_operation_base, only: ElementOperationBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private  
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  public :: ATM_PHY_RD_DGM_calc_heating
  
  !-----------------------------------------------------------------------------
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
contains
!OCL SERIAL
  subroutine ATM_PHY_RD_DGM_calc_heating( RHOH, &
    flux_rad, DDENS, DENS_hyd, CVtot, &
    lcmesh, elem3D, elem2D, &
    element3D_operation,  &
    TEMP_t )
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: RHOH(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: flux_rad(elem3D%Np,lcmesh%Ne,2,2)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA)   
    class(ElementOperationBase3D), intent(in) :: element3D_operation
    real(RP), intent(out), optional :: TEMP_t(elem3D%Np,lcmesh%NeA,2)

    integer :: ke
    real(RP) :: Fz_LW(elem3D%Np)
    real(RP) :: Fz_SW(elem3D%Np)
    real(RP) :: DENS(elem3D%Np)

    logical :: cal_TEMP_t
    !------------------------------------------------------------------------------

    if (present(TEMP_t)) then
      cal_TEMP_t = .true.
    else
      cal_TEMP_t = .false.
    end if

    !$omp parallel do private(Fz_LW, Fz_SW, DENS)
    do ke=lcmesh%NeS, lcmesh%NeE
      call element3D_operation%Dz( flux_rad(:,ke,I_LW,I_up) - flux_rad(:,ke,I_LW,I_dn), Fz_LW )
      call element3D_operation%Dz( flux_rad(:,ke,I_SW,I_up) - flux_rad(:,ke,I_SW,I_dn), Fz_SW )

      RHOH(:,ke) = - lcmesh%Escale(:,ke,3,3) * ( Fz_LW(:) + Fz_SW(:) )

      if ( cal_TEMP_t ) then
        DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
        TEMP_t(:,ke,1) = - lcmesh%Escale(:,ke,3,3) * Fz_LW(:)  / ( DENS(:) * CVtot(:,ke) )
        TEMP_t(:,ke,2) = - lcmesh%Escale(:,ke,3,3) * Fz_SW(:)  / ( DENS(:) * CVtot(:,ke) )
      end if
    end do

    return
  end subroutine ATM_PHY_RD_DGM_calc_heating
end module scale_atm_phy_rd_dgm_common