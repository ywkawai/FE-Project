!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_nonhydro3d_numdiff
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

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
  public :: atm_dyn_nonhydro3d_numdiff_Init
  public :: atm_dyn_nonhydro3d_numdiff_Final
  public :: atm_dyn_nonhydro3d_numdiff_tend 
  public :: atm_dyn_nonhydro3d_numdiff_cal_laplacian 
  public :: atm_dyn_nonhydro3d_numdiff_cal_flx
  
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
  
  private :: cal_del_gradDiffVar

contains
  subroutine atm_dyn_nonhydro3d_numdiff_Init( mesh )
    
    implicit none
    class(MeshBase3D), intent(in) :: mesh
    !--------------------------------------------


    return
  end subroutine atm_dyn_nonhydro3d_numdiff_Init


  subroutine atm_dyn_nonhydro3d_numdiff_Final()
    implicit none
    !--------------------------------------------
    
    return
  end subroutine atm_dyn_nonhydro3d_numdiff_Final  

  !-------------------------------

  subroutine atm_dyn_nonhydro3d_numdiff_tend( &
    tend_,                                                     & ! (out)
    GxV_, GyV_, GzV_,                                          & ! (in)
    DDENS_, DENS_hyd, diffcoef_h, diffcoef_v,                  & ! (in)
    Dx, Dy, Dz, Lift, lmesh, elem, mul_dens_flag               )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
    real(RP), intent(inout)  :: tend_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GxV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GyV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: diffcoef_h
    real(RP), intent(in)  :: diffcoef_v
    logical, intent(in) :: mul_dens_flag

    real(RP) ::  del_flux(elem%NfpTot,lmesh%Ne)
    real(RP) :: coef_h(elem%Np), coef_v(elem%Np)
    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    
    integer :: ke

    !--------------------------------------------
    
    call cal_del_flux_lap_with_coef( del_flux,                                & ! (out)
      GxV_, GyV_, GzV_,                                                       & ! (in)
      DDENS_, DENS_hyd, diffcoef_h, diffcoef_v,                               & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, mul_dens_flag )
        
    !$omp parallel do private( &
    !$omp coef_h, coef_v, Fx, Fy, Fz, LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE

      if (mul_dens_flag) then
        coef_h(:) = diffcoef_h * (DENS_hyd(:,ke) + DDENS_(:,ke))
        coef_v(:) = diffcoef_v * (DENS_hyd(:,ke) + DDENS_(:,ke))
      else
        coef_h(:) = diffcoef_h
        coef_v(:) = diffcoef_v
      end if
      call sparsemat_matmul(Dx, coef_h(:) * GxV_(:,ke), Fx)
      call sparsemat_matmul(Dy, coef_h(:) * GyV_(:,ke), Fy)
      call sparsemat_matmul(Dz, coef_v(:) * GzV_(:,ke) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx)

      tend_(:,ke) =   ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)      &
          + lmesh%Escale(:,ke,2,2) * Fy(:)      &
          + lmesh%Escale(:,ke,3,3) * Fz(:)      &
          + LiftDelFlx(:)                       )
    end do

    return
  end subroutine atm_dyn_nonhydro3d_numdiff_tend

  subroutine cal_del_flux_lap_with_coef( del_flux, & ! (out)
    GxV_, GyV_, GzV_,                              & ! (in)
    DDENS_, DENS_hyd, coef_h, coef_v,              & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, mul_dens_flag )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  GxV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GyV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzV_(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: coef_h
    real(RP), intent(in)  :: coef_v
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: mul_dens_flag
    
    integer :: i, iP, iM
    real(RP) :: densP, densM
    !------------------------------------------------------------------------

    !$omp parallel do                   &
    !$omp private( iM, iP, densP, densM  ) 
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      if (mul_dens_flag) then
        densM =  DENS_hyd(iM) + DDENS_(iM)      
        densP =  DENS_hyd(iP) + DDENS_(iP)
      else
        densM = 1.0_RP
        densP = 1.0_RP
      end if
      
      del_flux(i) = 0.5_RP * ( &
          ( 1.0_RP + sign(1.0_RP,nx(i)) ) * coef_h * ( densP * GxV_(iP) - densM * GxV_(iM) ) * nx(i) &
        + ( 1.0_RP + sign(1.0_RP,ny(i)) ) * coef_h * ( densP * GyV_(iP) - densM * GyV_(iM) ) * ny(i) &
        + ( 1.0_RP + sign(1.0_RP,nz(i)) ) * coef_v * ( densP * GzV_(iP) - densM * GzV_(iM) ) * nz(i) )
    end do

    return
  end subroutine cal_del_flux_lap_with_coef

  !--

  subroutine atm_dyn_nonhydro3d_numdiff_cal_laplacian( &
    lapla_h, lapla_v,                               & ! (out)
    GxV_, GyV_, GzV_,                               & ! (in)
    Dx, Dy, Dz, Lift, lmesh, elem                   ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
    real(RP), intent(out) :: lapla_h(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: lapla_v(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GxV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GyV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzV_(elem%Np,lmesh%NeA)

    real(RP) ::  del_flux_h(elem%NfpTot,lmesh%Ne)
    real(RP) ::  del_flux_v(elem%NfpTot,lmesh%Ne)
    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    
    integer :: ke
    !--------------------------------------------
    
    call cal_del_flux_lap( del_flux_h, del_flux_v,                            & ! (out)
      GxV_, GyV_, GzV_,                                                       & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem )
        
    !$omp parallel do private( &
    !$omp Fx, Fy, Fz, LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE

      call sparsemat_matmul(Dx, GxV_(:,ke), Fx)
      call sparsemat_matmul(Dy, GyV_(:,ke), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux_h(:,ke), LiftDelFlx)

      lapla_h(:,ke) = ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)      &
          + lmesh%Escale(:,ke,2,2) * Fy(:)      &
          + LiftDelFlx(:)                       &
          )
      
      call sparsemat_matmul(Dz, GzV_(:,ke) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux_v(:,ke), LiftDelFlx)
      lapla_v(:,ke) = (  &
            lmesh%Escale(:,ke,3,3) * Fz(:)      &
          + LiftDelFlx(:)                       &
          )
    end do

    return
  end subroutine atm_dyn_nonhydro3d_numdiff_cal_laplacian 


  subroutine cal_del_flux_lap( del_flux_h, del_flux_v, &
    GxV_, GyV_, GzV_,                     &
    nx, ny, nz, vmapM, vmapP, lmesh, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux_h(elem%NfpTot*lmesh%Ne)
    real(RP), intent(out) ::  del_flux_v(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  GxV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GyV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    !------------------------------------------------------------------------

    !$omp parallel do         &
    !$omp private( iM, iP ) 
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)
          
      del_flux_h(i) = 0.5_RP * ( &
          ( 1.0_RP + sign(1.0_RP,nx(i)) ) * ( GxV_(iP) - GxV_(iM) ) * nx(i) &
        + ( 1.0_RP + sign(1.0_RP,ny(i)) ) * ( GyV_(iP) - GyV_(iM) ) * ny(i) )
      
      del_flux_v(i) = 0.5_RP * ( &
          ( 1.0_RP + sign(1.0_RP,nz(i)) ) * ( GzV_(iP) - GzV_(iM) ) * nz(i) )
    end do

    return
  end subroutine cal_del_flux_lap

  !-------------------------------------------------------

  subroutine atm_dyn_nonhydro3d_numdiff_cal_flx( &
    GxV_, GyV_, GzV_,                                 & ! (out)
    Varh_, Varv_, DDENS_, DENS_hyd_,                  & ! (in)
    Dx, Dy, Dz, Lift, lmesh, elem, divide_dens_flag   ) ! (in)
    
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
    real(RP), intent(out)  :: GxV_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GyV_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Varh_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Varv_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd_(elem%Np,lmesh%NeA)
    logical, intent(in) :: divide_dens_flag

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: vh(elem%Np), vv(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,3)

    integer :: ke

    !------------------------------------------------------------------------------

    call cal_del_gradDiffVar( del_flux,                                       & ! (out)
      Varh_, Varv_, DDENS_, DENS_hyd_,                                        & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, divide_dens_flag )                                           ! (in)

    !$omp parallel do private(Fx, Fy, Fz, LiftDelFlx, vh, vv)
    do ke=lmesh%NeS, lmesh%NeE

      if (divide_dens_flag) then
        vh(:) = Varh_(:,ke) / (DDENS_(:,ke) + DENS_hyd_(:,ke))
        vv(:) = Varv_(:,ke) / (DDENS_(:,ke) + DENS_hyd_(:,ke))
      else
        vh(:) = Varh_(:,ke)
        vv(:) = Varv_(:,ke)      
      end if

      call sparsemat_matmul(Dx, vh, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,1), LiftDelFlx)
      GxV_(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dy, vh, Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,2), LiftDelFlx)
      GyV_(:,ke) = lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, vv, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,3), LiftDelFlx)
      GzV_(:,ke) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)
    end do

    return
  end subroutine atm_dyn_nonhydro3d_numdiff_cal_flx

  subroutine cal_del_gradDiffVar( del_flux, &
    VARh_, VARv_, DDENS_, DENS_hyd_,        &
    nx, ny, nz, vmapM, vmapP, lmesh, elem,  &
    divide_dens_flag )
    
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) :: del_flux(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(in) :: Varh_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: Varv_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: divide_dens_flag
    
    integer :: i, iP, iM
    real(RP) :: delVarh, delVarv
    real(RP) :: weight_P, weight_M
    !------------------------------------------------------------------------

    !$omp parallel do private(iM, iP, delVarh, delVarv,  weight_P, weight_M)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      if (divide_dens_flag) then
        weight_P = 1.0_RP / (DDENS_(iP) + DENS_hyd_(iP))
        weight_M = 1.0_RP / (DDENS_(iM) + DENS_hyd_(iM))
      else
        weight_M = 1.0_RP
        weight_P = 1.0_RP
      end if

      delVarh = 0.5_RP * (Varh_(iP) * weight_P - Varh_(iM) * weight_M)
      delVarv = 0.5_RP * (Varv_(iP) * weight_P - Varv_(iM) * weight_M)

      del_flux(i,1) = ( 1.0_RP - sign(1.0_RP,nx(i)) ) * delVarh * nx(i)
      del_flux(i,2) = ( 1.0_RP - sign(1.0_RP,ny(i)) ) * delVarh * ny(i)
      del_flux(i,3) = ( 1.0_RP - sign(1.0_RP,nz(i)) ) * delVarv * nz(i)
    end do

    return
  end subroutine cal_del_gradDiffVar

  !------------------------------------------------
end module scale_atm_dyn_nonhydro3d_numdiff
