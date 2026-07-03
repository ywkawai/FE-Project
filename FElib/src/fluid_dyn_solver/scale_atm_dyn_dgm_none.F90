!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Regional nonhydrostatic model / Skip dynamical process
!!
!! @par Description
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_none
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
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_element_operation_base, only: ElementOperationBase3D
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
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
  public :: atm_dyn_dgm_none_Init
  public :: atm_dyn_dgm_none_Final
  public :: atm_dyn_dgm_none_cal_tend

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

contains
  subroutine atm_dyn_dgm_none_Init( mesh )
    implicit none
    class(MeshBase3D), intent(in) :: mesh
    !--------------------------------------------
    return
  end subroutine atm_dyn_dgm_none_Init

  subroutine atm_dyn_dgm_none_Final()
    implicit none
    !--------------------------------------------
    return
  end subroutine atm_dyn_dgm_none_Final  

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_none_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,         & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_,         & ! (in)
    DENS_hyd, PRES_hyd, PRES_hyd_ref, THERM_hyd,         & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,      & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift,   & ! (in)
    lmesh, elem, lmesh2D, elem2D )                         ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    class(ElementOperationBase3D), intent(in) :: element3D_operation
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
    real(RP), intent(in)  :: THERM_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    integer :: ke
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_interior', 3)

    !$omp parallel do private( ke )
    do ke = lmesh%NeS, lmesh%NeE
      DENS_dt(:,ke) = 0.0_RP
      MOMX_dt(:,ke) = 0.0_RP
      MOMY_dt(:,ke) = 0.0_RP
      MOMZ_dt(:,ke) = 0.0_RP
      RHOT_dt(:,ke) = 0.0_RP 
    end do
    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_none_cal_tend
end module scale_atm_dyn_dgm_none
