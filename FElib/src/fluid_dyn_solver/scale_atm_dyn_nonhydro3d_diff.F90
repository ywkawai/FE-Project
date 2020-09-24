!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_nonhydro3d_diff
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
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D


  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_nonhydro3d_diff_Init
  public :: atm_dyn_nonhydro3d_diff_Final
  public :: atm_dyn_nonhydro3d_diff_cal_flx

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  integer, private, parameter :: VARS_DDENS_ID  = 1
  integer, private, parameter :: VARS_MOMX_ID   = 2
  integer, private, parameter :: VARS_MOMY_ID   = 3
  integer, private, parameter :: VARS_MOMZ_ID   = 4
  integer, private, parameter :: VARS_DRHOT_ID  = 5
  integer, private, parameter :: PROG_VARS_NUM  = 5
  
  integer, private, parameter :: VARS_GxU_ID      = 1
  integer, private, parameter :: VARS_GyU_ID      = 2
  integer, private, parameter :: VARS_GzU_ID      = 3
  integer, private, parameter :: VARS_GxV_ID      = 4
  integer, private, parameter :: VARS_GyV_ID      = 5
  integer, private, parameter :: VARS_GzV_ID      = 6
  integer, private, parameter :: VARS_GxW_ID      = 7
  integer, private, parameter :: VARS_GyW_ID      = 8
  integer, private, parameter :: VARS_GzW_ID      = 9
  integer, private, parameter :: VARS_GxPT_ID     = 10
  integer, private, parameter :: VARS_GyPT_ID     = 11
  integer, private, parameter :: VARS_GzPT_ID     = 12
  integer, private, parameter :: AUX_DIFFVARS_NUM = 12

  private :: cal_del_gradDiffVar

contains
  subroutine atm_dyn_nonhydro3d_diff_Init( mesh )

    implicit none
    class(MeshBase3D), intent(in) :: mesh

    integer :: p1, p2, p_
    type(ElementBase3D), pointer :: elem
    real(RP) :: invV_VPOrdM1(mesh%refElem3D%Np,mesh%refElem3D%Np)
    !--------------------------------------------

    elem => mesh%refElem3D

    return
  end subroutine atm_dyn_nonhydro3d_diff_Init


  subroutine atm_dyn_nonhydro3d_diff_Final()
    implicit none
    !--------------------------------------------
    
    return
  end subroutine atm_dyn_nonhydro3d_diff_Final  

  !-------------------------------

  subroutine atm_dyn_nonhydro3d_diff_cal_flx( &
    GxU_, GyU_, GzU_, GxV_, GyV_, GzV_, GxW_, GyW_, GzW_, GxPT_, GyPT_, GzPT_,   &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                     &
    Dx, Dy, Dz, Lift, lmesh, elem )
    
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
    real(RP), intent(out)  :: GxU_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GyU_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzU_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GxV_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GyV_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzV_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GxW_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GyW_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzW_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GxPT_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GyPT_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzPT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)

    real(RP) :: DENS_(elem%Np), U_(elem%Np), V_(elem%Np), W_(elem%Np)
    real(RP) :: DTHETA_(elem%Np), RHOT_(elem%Np), RHOT_hyd(elem%Np)
    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,AUX_DIFFVARS_NUM)

    integer :: ke
    !------------------------------------------------------------------------------

    call cal_del_gradDiffVar( del_flux,                                       & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem )                                                             ! (in)

    do ke=lmesh%NeS, lmesh%NeE
      DENS_(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) 
      RHOT_(:) = RHOT_hyd(:) + DRHOT_(:,ke)

      U_(:) = MOMX_(:,ke)/DENS_(:)
      V_(:) = MOMY_(:,ke)/DENS_(:)
      W_(:) = MOMZ_(:,ke)/DENS_(:)
      DTHETA_(:) = RHOT_(:)/DENS_(:) - RHOT_hyd(:)/DENS_hyd(:,ke)

      !- U

      call sparsemat_matmul(Dx, U_, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GxU_ID), LiftDelFlx)
      GxU_(:,ke) = lmesh%Escale(:,ke,1,1)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dy, U_, Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GyU_ID), LiftDelFlx)
      GyU_(:,ke) = lmesh%Escale(:,ke,2,2)*Fy(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, U_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GzU_ID), LiftDelFlx)
      GzU_(:,ke) = lmesh%Escale(:,ke,3,3)*Fz(:) + LiftDelFlx(:)

      !- V

      call sparsemat_matmul(Dx, V_, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GxV_ID), LiftDelFlx)
      GxV_(:,ke) = lmesh%Escale(:,ke,1,1)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dy, V_, Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GyV_ID), LiftDelFlx)
      GyV_(:,ke) = lmesh%Escale(:,ke,2,2)*Fy(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, V_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GzV_ID), LiftDelFlx)
      GzV_(:,ke) = lmesh%Escale(:,ke,3,3)*Fz(:) + LiftDelFlx(:)

      !- W

      call sparsemat_matmul(Dx, W_, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GxW_ID), LiftDelFlx)
      GxW_(:,ke) = lmesh%Escale(:,ke,1,1)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dy, W_, Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GyW_ID), LiftDelFlx)
      GyW_(:,ke) = lmesh%Escale(:,ke,2,2)*Fy(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, W_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GzW_ID), LiftDelFlx)
      GzW_(:,ke) = lmesh%Escale(:,ke,3,3)*Fz(:) + LiftDelFlx(:)

      !- PT

      call sparsemat_matmul(Dx, DTHETA_, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GxPT_ID), LiftDelFlx)
      GxPT_(:,ke) = lmesh%Escale(:,ke,1,1)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dy, DTHETA_, Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GyPT_ID), LiftDelFlx)
      GyPT_(:,ke) = lmesh%Escale(:,ke,2,2)*Fy(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, DTHETA_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GzPT_ID), LiftDelFlx)
      GzPT_(:,ke) = lmesh%Escale(:,ke,3,3)*Fz(:) + LiftDelFlx(:)

    end do

    return
  end subroutine atm_dyn_nonhydro3d_diff_cal_flx

  subroutine cal_del_gradDiffVar( del_flux,                  &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, &
    nx, ny, nz, vmapM, vmapP, lmesh, elem )
    
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne,AUX_DIFFVARS_NUM)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    real(RP) :: densM, densP, rhot_hyd, rhotM, rhotP
    real(RP) :: delU, delV, delW, delPT
    real(RP) :: MOMX_P, MOMY_P, MOMZ_P 
    real(RP) :: gamm, rgamm
    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      rhot_hyd = PRES00/Rdry * (PRES_hyd(iM)/PRES00)**rgamm
      rhotM = rhot_hyd + DRHOT_(iM)
      rhotP = rhot_hyd + DRHOT_(iP) 
      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iM)

      delU = 0.5_RP*(MOMX_(iP)/densP - MOMX_(iM)/densM)
      delV = 0.5_RP*(MOMY_(iP)/densP - MOMY_(iM)/densM)
      delW = 0.5_RP*(MOMZ_(iP)/densP - MOMZ_(iM)/densM)
      delPT = 0.5_RP*(rhotP/densP - rhotM/densM)

      !- U
      del_flux(i,VARS_GxU_ID) = delU * nx(i)
      del_flux(i,VARS_GyU_ID) = delU * ny(i)
      del_flux(i,VARS_GzU_ID) = delU * nz(i)

      !- V
      del_flux(i,VARS_GxV_ID) = delV * nx(i)
      del_flux(i,VARS_GyV_ID) = delV * ny(i)
      del_flux(i,VARS_GzV_ID) = delV * nz(i)

      !- W
      del_flux(i,VARS_GxW_ID) = delW * nx(i)
      del_flux(i,VARS_GyW_ID) = delW * ny(i)
      del_flux(i,VARS_GzW_ID) = delW * nz(i)

      !- PT
      del_flux(i,VARS_GxPT_ID) = delPT * nx(i)
      del_flux(i,VARS_GyPT_ID) = delPT * ny(i)
      del_flux(i,VARS_GzPT_ID) = delPT * nz(i)
    end do

    return
  end subroutine cal_del_gradDiffVar

  !------------------------------------------------
end module scale_atm_dyn_nonhydro3d_diff
