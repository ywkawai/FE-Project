!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_nonhydro3d_hevi
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
  public :: atm_dyn_nonhydro3d_hevi_Init
  public :: atm_dyn_nonhydro3d_hevi_Final
  public :: atm_dyn_nonhydro3d_hevi_cal_tend
  public :: atm_dyn_nonhydro3d_hevi_cal_vi

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  integer, private, parameter :: DDENS_VID  = 1
  integer, private, parameter :: MOMX_VID   = 2
  integer, private, parameter :: MOMY_VID   = 3
  integer, private, parameter :: MOMZ_VID   = 4
  integer, private, parameter :: DRHOT_VID  = 5
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

  real(RP), private, allocatable :: IntrpMat_VPOrdM1(:,:)

  private :: cal_del_flux_dyn

contains
  subroutine atm_dyn_nonhydro3d_hevi_Init( mesh )

    implicit none
    class(MeshBase3D), intent(in) :: mesh

    integer :: p1, p2, p3, p_, p_intrp
    type(ElementBase3D), pointer :: elem
    real(RP) :: invV_POrdM1(mesh%refElem3D%Np,mesh%refElem3D%Np)

    real(RP), allocatable :: InvV_(:,:), InvV_ip(:,:)
    !--------------------------------------------

    elem => mesh%refElem3D
    allocate( IntrpMat_VPOrdM1(elem%Np,elem%Np) )
    
    invV_POrdM1(:,:) = elem%invV
    do p2=1, elem%Nnode_h1D
    do p1=1, elem%Nnode_h1D
      p_ = p1 + (p2-1)*elem%Nnode_h1D + (elem%Nnode_v-1)*elem%Nnode_h1D**2
      invV_POrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem%V, invV_POrdM1)

    !---
    return
  end subroutine atm_dyn_nonhydro3d_hevi_Init


  subroutine atm_dyn_nonhydro3d_hevi_Final()
    implicit none
    !--------------------------------------------

    deallocate( IntrpMat_VPOrdM1 )    
    return
  end subroutine atm_dyn_nonhydro3d_hevi_Final  

  !-------------------------------

  subroutine atm_dyn_nonhydro3d_hevi_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, CORIOLIS,          & ! (in)
    GxU_, GyU_, GzU_, GxV_, GyV_, GzV_, GxW_, GyW_, GzW_, GxPT_, GyPT_, GzPT_,  & ! (in)
    viscCoef_h, viscCoef_v, diffCoef_h, diffCoef_v,                             & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
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
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: GxU_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GyU_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzU_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GxV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GyV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GxW_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GyW_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzW_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GxPT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GyPT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzPT_(elem%Np,lmesh%NeA) 
    real(RP), intent(in) :: viscCoef_h
    real(RP), intent(in) :: viscCoef_v
    real(RP), intent(in) :: diffCoef_h
    real(RP), intent(in) :: diffCoef_v

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: dens_(elem%Np), RHOT_hyd(elem%Np), RHOT_(elem%Np), dpres_(elem%Np)
    real(RP) :: pres_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), pot_(elem%Np)
    real(RP) :: Cori(elem%Np)

    real(RP) :: tmp(elem%Np)
    integer :: ke, ke2d
    real(RP) :: gamm, rgamm
    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 3)
    call cal_del_flux_dyn( del_flux,                                          & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                & ! (in)
      GxU_, GyU_, GzU_, GxV_, GyV_, GzV_, GxW_, GyW_, GzW_,                   & ! (in)
      GxPT_, GyPT_, GzPT_,                                                    & ! (in)
      viscCoef_h, viscCoef_v,                                                 & ! (in)
      diffCoef_h, diffCoef_v,                                                 & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem )                                                             ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 3)
    gamm = CPDry / CvDry
    rgamm = CvDry / CpDry

    !$omp parallel do private( &
    !$omp RHOT_hyd, RHOT_, pres_, dpres_, dens_, u_, v_, w_, pot_, ke2d, Cori, &
    !$omp Fx, Fy, Fz,LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE
      !--
      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**rgamm
      RHOT_(:) = RHOT_hyd(:) + DRHOT_(:,ke)
      pres_(:) = PRES_hyd(:,ke) * (1.0_RP + DRHOT_(:,ke)/RHOT_hyd(:))**gamm
      dpres_(:) = pres_(:) - PRES_hyd(:,ke)
      dens_(:) = DDENS_(:,ke) + DENS_hyd(:,ke)

      u_(:) = MOMX_(:,ke)/dens_(:)
      v_(:) = MOMY_(:,ke)/dens_(:)
      w_(:) = MOMZ_(:,ke)/dens_(:)
      pot_(:) = RHOT_(:)/dens_(:)

      ke2d = lmesh%EMap3Dto2D(ke)
      Cori(:) = CORIOLIS(elem%IndexH2Dto3D(:),ke2d)

      !-- DENS
      call sparsemat_matmul(Dx, MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, MOMY_(:,ke), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,DDENS_VID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + LiftDelFlx(:) ) 
      
      !-- MOMX
      call sparsemat_matmul(Dx, u_(:)*MOMX_(:,ke) + pres_(:)  - viscCoef_h*dens_(:)*GxU_(:,ke), Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMX_(:,ke)              - viscCoef_h*dens_(:)*GyU_(:,ke), Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMX_(:,ke)              - viscCoef_v*dens_(:)*GzU_(:,ke) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:)      &
          + lmesh%Escale(:,ke,2,2) * Fy(:)      &
          + lmesh%Escale(:,ke,3,3) * Fz(:)      &
          + LiftDelFlx(:)                       &
          - Cori(:)*MOMY_(:,ke)                 &
          ) 

      !-- MOMY
      call sparsemat_matmul(Dx, u_(:)*MOMY_(:,ke)              - viscCoef_h*dens_(:)*GxV_(:,ke), Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMY_(:,ke) + pres_(:)   - viscCoef_h*dens_(:)*GyV_(:,ke), Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMY_(:,ke)              - viscCoef_v*dens_(:)*GzV_(:,ke) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDelFlx(:)                  &
          + Cori(:)*MOMX_(:,ke) &
          ) 

      !-- MOMZ
      call sparsemat_matmul(Dx, u_(:)*MOMZ_(:,ke) - viscCoef_h*dens_(:)*GxW_(:,ke), Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMZ_(:,ke) - viscCoef_h*dens_(:)*GyW_(:,ke), Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMZ_(:,ke) - viscCoef_v*dens_(:)*GzW_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)   &
          + lmesh%Escale(:,ke,2,2) * Fy(:)   &
          + lmesh%Escale(:,ke,3,3) * Fz(:)   &
          + LiftDelFlx(:)                )   

      !-- RHOT
      call sparsemat_matmul(Dx, pot_(:)*MOMX_(:,ke) - diffCoef_h*dens_(:)*GxPT_(:,ke), Fx)
      call sparsemat_matmul(Dy, pot_(:)*MOMY_(:,ke) - diffCoef_h*dens_(:)*GyPT_(:,ke), Fy)
      call sparsemat_matmul(Dz,                     - diffCoef_v*dens_(:)*GzPT_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,DRHOT_VID), LiftDelFlx)
      
      RHOT_dt(:,ke) =  - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDelFlx(:) ) 
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_nonhydro3d_hevi_cal_tend

  !------

  subroutine cal_del_flux_dyn( del_flux, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,   &
    GxU_, GyU_, GzU_, GxV_, GyV_, GzV_, GxW_, GyW_, GzW_,      &
    GxPT_, GyPT_, GzPT_,                                       &
    viscCoef_h, viscCoef_v,                                    &
    diffCoef_h, diffCoef_v,                                    &
    nx, ny, nz, vmapM, vmapP, lmesh, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne,PROG_VARS_NUM)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GxU_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GyU_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzU_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GxV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GyV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GxW_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GyW_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzW_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GxPT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GyPT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzPT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: viscCoef_h
    real(RP), intent(in) :: viscCoef_v
    real(RP), intent(in) :: diffCoef_h
    real(RP), intent(in) :: diffCoef_v
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, MOMZ_P, alpha, swV
    real(RP) :: presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P
    real(RP) :: dDiffFluxU, dDiffFluxV, dDiffFluxW, dDiffFluxPT
    real(RP) :: gamm, rgamm
    real(RP) :: mu
    logical :: visc_flag, diff_flag

    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    if (viscCoef_h > 0.0_RP .or. viscCoef_v > 0.0_RP) then
      visc_flag = .true.
    else 
      visc_flag = .false.
      dDiffFluxU  = 0.0_RP
      dDiffFluxV  = 0.0_RP
      dDiffFluxW  = 0.0_RP
    end if 
   
    if (diffCoef_h > 0.0_RP .or. diffCoef_v > 0.0_RP) then
      diff_flag = .true.
    else
      diff_flag = .false.
      dDiffFluxPT = 0.0_RP
    end if     

    !$omp parallel do  &
    !$omp private( iM, iP, VelP, VelM, MOMZ_P, alpha, swV, &
    !$omp presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P ) &
    !$omp firstprivate( dDiffFluxU, dDiffFluxV, dDiffFluxW, dDiffFluxPT, mu )
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      rhot_hyd_M = PRES00/Rdry * (PRES_hyd(iM)/PRES00)**rgamm
      rhot_hyd_P = PRES00/Rdry * (PRES_hyd(iP)/PRES00)**rgamm
      
      rhotM = rhot_hyd_M + DRHOT_(iM)
      presM = PRES_hyd(iM) * (1.0_RP + DRHOT_(iM)/rhot_hyd_M)**gamm
      dpresM = presM - PRES_hyd(iM)*abs(nz(i))

      rhotP = rhot_hyd_P + DRHOT_(iP) 
      presP = PRES_hyd(iP) * (1.0_RP + DRHOT_(iP)/rhot_hyd_P)**gamm
      dpresP = presP - PRES_hyd(iP)*abs(nz(i))

      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iP)

      swV = 1.0_RP - nz(i)**2
      VelM = (MOMX_(iM)*nx(i) + MOMY_(iM)*ny(i) + MOMZ_(iM)*nz(i))/densM
      VelP = (MOMX_(iP)*nx(i) + MOMY_(iP)*ny(i) + MOMZ_(iP)*nz(i))/densP
      MOMZ_P = MOMZ_(iP)
      
      alpha = swV*max( sqrt(gamm * presM / densM) + abs(VelM), &
                       sqrt(gamm * presP / densP) + abs(VelP)  )
          
      mu = (2.0_RP * dble((elem%PolyOrder_h+1)*(elem%PolyOrder_h+2)) / 2.0_RP / 600.0_RP) 

      
      if ( visc_flag ) then
        dDiffFluxU = ( &
            viscCoef_h*(densP*GxU_(iP) - densM*GxU_(iM))*nx(i)     &
          + viscCoef_h*(densP*GyU_(iP) - densM*GyU_(iM))*ny(i)     &
          + viscCoef_v*(densP*GzU_(iP) - densM*GzU_(iM))*nz(i)     &
          + mu*(densP + densM)*(MOMX_(iP)/densP - MOMX_(iM)/densM) )

        dDiffFluxV = ( &
            viscCoef_h*(densP*GxV_(iP) - densM*GxV_(iM))*nx(i)     &
          + viscCoef_h*(densP*GyV_(iP) - densM*GyV_(iM))*ny(i)     &
          + viscCoef_v*(densP*GzV_(iP) - densM*GzV_(iM))*nz(i)     &
          + mu*(densP + densM)*(MOMY_(iP)/densP - MOMY_(iM)/densM) )
        
        dDiffFluxW = ( &
            viscCoef_h*(densP*GxW_(iP) - densM*GxW_(iM))*nx(i)     &
          + viscCoef_h*(densP*GyW_(iP) - densM*GyW_(iM))*ny(i)     &
          + viscCoef_v*(densP*GzW_(iP) - densM*GzW_(iM))*nz(i)     &
          + mu*(densP + densM)*(MOMZ_(iP)/densP - MOMZ_(iM)/densM) )
      end if
      if ( diff_flag ) then
        dDiffFluxPT = ( &
            diffCoef_h*(densP*GxPT_(iP) - densM*GxPT_(iM))*nx(i)   &
          + diffCoef_h*(densP*GyPT_(iP) - densM*GyPT_(iM))*ny(i)   &            
          + diffCoef_v*(densP*GzPT_(iP) - densM*GzPT_(iM))*nz(i)   &
          + mu*(densP + densM)*(rhotP/densP - rhotM/densM)         )
      end if
      
      del_flux(i,DDENS_VID) = 0.5_RP*(                    &
                      (MOMX_(iP) - MOMX_(iM))*nx(i)       &
                    + (MOMY_(iP) - MOMY_(iM))*ny(i)       &
                    - alpha * (DDENS_(iP) - DDENS_(iM))   )
      
      del_flux(i,MOMX_VID) = 0.5_RP*(                     &
                      ( MOMX_(iP)*VelP - MOMX_(iM)*VelM ) &
                    + ( dpresP - dpresM )*nx(i)           &
                    - alpha * (MOMX_(iP) - MOMX_(iM))     &
                    - dDiffFluxU                   )

      del_flux(i,MOMY_VID) = 0.5_RP*(                      &
                      ( MOMY_(iP)*VelP - MOMY_(iM)*VelM )  &
                    + ( dpresP - dpresM )*ny(i)            &
                    - alpha * (MOMY_(iP) - MOMY_(iM))      &
                    - dDiffFluxV                   )               

      del_flux(i,MOMZ_VID) = 0.5_RP*(                      &
                      ( MOMZ_P*VelP - MOMZ_(iM)*VelM )     &
                    - alpha * (MOMZ_(iP) - MOMZ_(iM))      &
                    - dDiffFluxW                   )
      
      del_flux(i,DRHOT_VID) = 0.5_RP*(                    &
                      swV*( rhotP*VelP - rhotM*VelM )     &
                    - alpha *(DRHOT_(iP) - DRHOT_(iM))    &
                    - dDiffFluxPT                  )
    end do

    return
  end subroutine cal_del_flux_dyn

  subroutine atm_dyn_nonhydro3d_hevi_cal_vi( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,             & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, & ! (in)
    Dz, Lift, impl_fac, lmesh, elem, lmesh2D, elem2D )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
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
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: impl_fac

    real(RP) :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: PROG_VARS00(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: b(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: b1D(elem%Nnode_v,PROG_VARS_NUM,lmesh%NeZ,elem%Nnode_h1D**2)
    integer :: ipiv(elem%Nnode_v*PROG_VARS_NUM*lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: tend(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: DENS_hyd_z(elem%Np,lmesh%NeZ)
    real(RP) :: PRES_hyd_z(elem%Np,lmesh%NeZ)
    real(RP) :: nz(elem%NfpTot,lmesh%NeZ)
    integer :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer :: ke_x, ke_y, ke_z, ke, p, v
    integer :: itr_lin, itr_nlin
    integer :: f, vs, ve, kl, ku, nz_1D
    integer :: ij, info
    logical :: is_converged


    real(RP), allocatable :: PmatBnd(:,:,:)

    !------------------------------------------------------------------------

    
    call PROF_rapstart( 'hevi_cal_vi_prep', 3)

    nz_1D = elem%Nnode_v * PROG_VARS_NUM * lmesh%NeZ
    kl = 2 * elem%Nnode_v * PROG_VARS_NUM - 1
    ku = kl
    allocate( PmatBnd(2*kl+ku+1,nz_1D,elem%Nnode_h1D**2) )

    !$omp parallel private(f, vs, ve)
    !$omp do
    do ke_z=1, lmesh%NeZ
      do f=1, elem%Nfaces_h
        vs = 1 + (f-1)*elem%Nfp_h
        ve = vs + elem%Nfp_h - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_h(:,f) + (ke_z-1)*elem%Np
      end do
      do f=1, elem%Nfaces_v
        vs = elem%Nfp_h*elem%Nfaces_h + 1 + (f-1)*elem%Nfp_v
        ve = vs + elem%Nfp_v - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_v(:,f) + (ke_z-1)*elem%Np
      end do
      vmapP(:,ke_z) = vmapM(:,ke_z)
    end do
    !$omp do
    do ke_z=1, lmesh%NeZ
      vs = elem%Nfp_h*elem%Nfaces_h + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z > 1) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,2) + (ke_z-2)*elem%Np

      vs = elem%Nfp_h*elem%Nfaces_h + elem%Nfp_v + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z < lmesh%NeZ) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1) + ke_z*elem%Np
    end do
    !$omp end parallel
    call PROF_rapend( 'hevi_cal_vi_prep', 3)

    do ke_y=1, lmesh%NeY
    do ke_x=1, lmesh%NeX

      call PROF_rapstart( 'hevi_cal_vi_get_var', 3)

      !$omp parallel do private(ke)
      do ke_z=1, lmesh%NeZ
        ke = ke_x + (ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

        PROG_VARS(:,DDENS_VID,ke_z) = DDENS_(:,ke)
        PROG_VARS(:,MOMX_VID,ke_z) = MOMX_(:,ke)
        PROG_VARS(:,MOMY_VID,ke_z) = MOMY_(:,ke)
        PROG_VARS(:,MOMZ_VID,ke_z) = MOMZ_(:,ke)
        PROG_VARS(:,DRHOT_VID,ke_z) = DRHOT_(:,ke)
        DENS_hyd_z(:,ke_z) = DENS_hyd(:,ke)
        PRES_hyd_z(:,ke_z) = PRES_hyd(:,ke)

        PROG_VARS0(:,:,ke_z) = PROG_VARS(:,:,ke_z)
        PROG_VARS00(:,:,ke_z) = PROG_VARS(:,:,ke_z)
        nz(:,ke_z) = lmesh%normal_fn(:,ke,3)
      end do
      call PROF_rapend( 'hevi_cal_vi_get_var', 3)
      
      if ( abs(impl_fac) > 0.0_RP ) then
        call PROF_rapstart( 'hevi_cal_vi_itr', 3)

        do itr_nlin = 1, 1

          call vi_eval_Ax( Ax(:,:,:),                       & ! (out)
            PROG_VARS, PROG_VARS0, DENS_hyd_z, PRES_hyd_z,  & ! (in)
            Dz, Lift, impl_fac, lmesh, elem,                & ! (in)
            nz, vmapM, vmapP, ke_x, ke_y, .false. )

          do ke_z=1, lmesh%NeZ
            b(:,:,ke_z) = - Ax(:,:,ke_z) + PROG_VARS00(:,:,ke_z)
          end do

          call PROF_rapstart( 'hevi_cal_vi_matbnd', 3)

          call vi_construct_matbnd( PmatBnd,                & ! (out)
            kl, ku, nz_1D,                                  & ! (in)
            PROG_VARS0, DENS_hyd_z, PRES_hyd_z,             & ! (in)
            Dz, Lift, impl_fac, lmesh, elem,                & ! (in)
            nz, vmapM, vmapP, ke_x, ke_y  )

          call PROF_rapend( 'hevi_cal_vi_matbnd', 3)
          
          call PROF_rapstart( 'hevi_cal_vi_lin', 3)
          !$omp parallel private(ij, v, ke_z, info)
          !$omp do
          do ij=1, elem%Nnode_h1D**2
            do ke_z=1, lmesh%NeZ
            do v=1, PROG_VARS_NUM
              b1D(:,v,ke_z,ij) = b(elem%Colmask(:,ij),v,ke_z)
            end do
            end do

            call dgbsv( nz_1D, kl, ku, 1, PmatBnd(:,:,ij), 2*kl+ku+1, ipiv(:,ij), b1D(:,:,:,ij), nz_1D, info)

            do ke_z=1, lmesh%NeZ
            do v=1, PROG_VARS_NUM
              PROG_VARS(elem%Colmask(:,ij),v,ke_z) = PROG_VARS(elem%Colmask(:,ij),v,ke_z) + b1D(:,v,ke_z,ij)
            end do
            end do
          end do
          !$omp do 
          do ke_z=1, lmesh%NeZ
            PROG_VARS0(:,:,ke_z) = PROG_VARS(:,:,ke_z)
          end do
          !$omp end parallel
          call PROF_rapend( 'hevi_cal_vi_lin', 3)
        end do ! itr nlin

        call PROF_rapend( 'hevi_cal_vi_itr', 3)
      end if

      call PROF_rapstart( 'hevi_cal_vi_retrun_var', 3)
      if ( abs(impl_fac) > 0.0_RP) then
        !$omp parallel do 
        do ke_z=1, lmesh%NeZ
          tend(:,:,ke_z) = (- PROG_VARS(:,:,ke_z) + PROG_VARS00(:,:,ke_z))/impl_fac
        end do
      else
        call vi_eval_Ax( tend(:,:,:),                      & ! (out)
          PROG_VARS, PROG_VARS, DENS_hyd_z, PRES_hyd_z,    & ! (in)
          Dz, Lift, impl_fac, lmesh, elem,                 & ! (in)
          nz, vmapM, vmapP, ke_x, ke_y, .true. )
      end if

      !$omp parallel do private(ke)
      do ke_z=1, lmesh%NeZ
        ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
        DENS_dt(:,ke) = - tend(:,DDENS_VID,ke_z)
        MOMX_dt(:,ke) = - tend(:,MOMX_VID,ke_z)
        MOMY_dt(:,ke) = - tend(:,MOMY_VID,ke_z)
        MOMZ_dt(:,ke) = - tend(:,MOMZ_VID,ke_z)
        RHOT_dt(:,ke) = - tend(:,DRHOT_VID,ke_z)
      end do
      call PROF_rapend( 'hevi_cal_vi_retrun_var', 3)
    end do
    end do

    return
  end subroutine atm_dyn_nonhydro3d_hevi_cal_vi

  !------------------------------------------------

  !---
  subroutine vi_eval_Ax( Ax, &
    PROG_VARS, PROG_VARS0, DENS_hyd, PRES_hyd,        & ! (in)
    Dz, Lift, impl_fac, lmesh, elem,                  & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y, cal_tend_flag )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y
    logical, intent(in) :: cal_tend_flag

    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%NeZ,PROG_VARS_NUM)
    real(RP) :: RHOT_hyd(elem%Np), POT(elem%Np)
    real(RP) :: DPRES(elem%Np)
    integer :: ke_z
    integer :: ke
    integer :: v
    real(RP) :: gamm, rgamm

    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    call vi_cal_del_flux_dyn( del_flux,                    & ! (out)
      PROG_VARS(:,DDENS_VID,:), PROG_VARS(:,MOMX_VID,:),   & ! (in)
      PROG_VARS(:,MOMY_VID ,:), PROG_VARS(:,MOMZ_VID,:),   & ! (in)
      PROG_VARS(:,DRHOT_VID,:),                            & ! (in)
      PROG_VARS0(:,DDENS_VID,:), PROG_VARS0(:,MOMX_VID,:), & ! (in)
      PROG_VARS0(:,MOMY_VID ,:), PROG_VARS0(:,MOMZ_VID,:), & ! (in)
      PROG_VARS0(:,DRHOT_VID,:),                           & ! (in)
      DENS_hyd, PRES_hyd, nz, vmapM, vmapP,                & ! (in)
      lmesh, elem )                                          ! (in)

    !$omp parallel do private(ke, RHOT_hyd, DPRES, POT, Fz, LiftDelFlx)
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(:,ke_z)/PRES00)**rgamm

      DPRES(:) = PRES_hyd(:,ke_z) * ((1.0_RP + PROG_VARS(:,DRHOT_VID,ke_z)/RHOT_hyd(:))**gamm - 1.0_RP)
      POT(:) = (RHOT_hyd(:) + PROG_VARS(:,DRHOT_VID,ke_z))/(DENS_hyd(:,ke_z) + PROG_VARS(:,DDENS_VID,ke_z))

      !- DENS
      call sparsemat_matmul(Dz, PROG_VARS(:,MOMZ_VID,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,DDENS_VID), LiftDelFlx)
      Ax(:,DDENS_VID,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      !- MOMX
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,MOMX_VID), LiftDelFlx)
      Ax(:,MOMX_VID,ke_z) = LiftDelFlx(:)

      !-MOMY
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,MOMY_VID), LiftDelFlx)
      Ax(:,MOMY_VID,ke_z) = LiftDelFlx(:)

      !-MOMZ
      call sparsemat_matmul(Dz, DPRES(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,MOMZ_VID), LiftDelFlx)
      Ax(:,MOMZ_VID,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)             &
                          + Grav * matmul(IntrpMat_VPOrdM1, PROG_VARS(:,DDENS_VID,ke_z)) 


      !-RHOT
      call sparsemat_matmul(Dz, POT(:)*PROG_VARS(:,MOMZ_VID,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,DRHOT_VID), LiftDelFlx)
      Ax(:,DRHOT_VID,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      !--
      if ( .not. cal_tend_flag ) then
        Ax(:,:,ke_z) =  PROG_VARS(:,:,ke_z) + impl_fac * Ax(:,:,ke_z)
      end if 

    end do    

    return
  end subroutine vi_eval_Ax

  subroutine vi_cal_del_flux_dyn( del_flux, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,              &
    DDENS0_, MOMX0_, MOMY0_, MOMZ0_, DRHOT0_,              &
    DENS_hyd, PRES_hyd, nz, vmapM, vmapP, lmesh, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%NeZ,PROG_VARS_NUM)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  DDENS0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMX0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMY0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  MOMZ0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  DRHOT0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%NeZ)
    
    integer :: i, p, ke_z, iP, iM
    real(RP) :: alpha0, swV
    real(RP) ::  MOMZ_P
    real(RP) :: rhot_hyd_M, rhot_hyd_P
    real(RP) :: dpresM, dpresP, densM, densP,  pottM, pottP
    real(RP) :: pres0M, pres0P, dens0M, dens0P
    real(RP) :: gamm, rgamm
    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry
    
    !$omp parallel do private( p, i, iM, iP,                  &
    !$omp rhot_hyd_M, rhot_hyd_P, densM, densP, pottM, pottP, &
    !$omp dpresM, dpresP, MOMZ_P,                             &
    !$omp dens0M, dens0P, pres0M, pres0P,                     &
    !$omp swV, alpha0 )    
    do ke_z=1, lmesh%NeZ
    do p=1, elem%NfpTot
      i = p + (ke_z-1)*elem%NfpTot
      iM = vmapM(i); iP = vmapP(i)

      !-
      rhot_hyd_M = PRES00/Rdry * (PRES_hyd(iM)/PRES00)**rgamm
      rhot_hyd_P = PRES00/Rdry * (PRES_hyd(iP)/PRES00)**rgamm

      densM = DENS_hyd(iM) + DDENS_(iM)
      densP = DENS_hyd(iP) + DDENS_(iP)
      
      pottM = (rhot_hyd_M + DRHOT_(iM)) / densM
      pottP = (rhot_hyd_P + DRHOT_(iP)) / densP
   
      dpresM = PRES_hyd(iM) * ((1.0_RP + DRHOT_(iM)/rhot_hyd_M)**gamm - 1.0_RP) 
      dpresP = PRES_hyd(iP) * ((1.0_RP + DRHOT_(iP)/rhot_hyd_P)**gamm - 1.0_RP) 

      !-
      dens0M = DENS_hyd(iM) + DDENS0_(iM)
      dens0P = DENS_hyd(iP) + DDENS0_(iP)

      pres0M = PRES_hyd(iM) * (1.0_RP + DRHOT0_(iM)/rhot_hyd_M)**gamm
      pres0P = PRES_hyd(iP) * (1.0_RP + DRHOT0_(iP)/rhot_hyd_P)**gamm

      swV = nz(i)**2      
      alpha0 = swV * max( abs(MOMZ0_(iM)/dens0M) + sqrt(gamm * pres0M/dens0M), &
                          abs(MOMZ0_(iP)/dens0P) + sqrt(gamm * pres0P/dens0P)  )
      
      if (iM==iP .and. (ke_z == 1 .or. ke_z == lmesh%NeZ)) then
        MOMZ_P = - MOMZ_(iM)
      else
        MOMZ_P = MOMZ_(iP)
      end if

      del_flux(i,DDENS_VID) = 0.5_RP * (                   &
                    + ( MOMZ_P - MOMZ_(iM) ) * nz(i)       &
                    - alpha0 * ( DDENS_(iP) - DDENS_(iM) ) )
      
      del_flux(i,MOMX_VID) = 0.5_RP * (                   &
                    - alpha0 * (  MOMX_(iP) - MOMX_(iM) ) )
      
      del_flux(i,MOMY_VID) = 0.5_RP * (                   &  
                    - alpha0 * ( MOMY_(iP) - MOMY_(iM) )  )               
      
      del_flux(i,MOMZ_VID) = 0.5_RP * (                   &
                    + ( dpresP - dpresM ) * nz(i)         &                    
                    - alpha0 * ( MOMZ_P - MOMZ_(iM) )     )
      
      del_flux(i,DRHOT_VID) = 0.5_RP * (                               &
                    + ( pottP * MOMZ_P  - pottM * MOMZ_(iM) ) * nz(i)  &
                    - alpha0 * ( DRHOT_(iP) - DRHOT_(iM) )             )
    end do
    end do

    return
  end subroutine vi_cal_del_flux_dyn

  subroutine vi_construct_matbnd( PmatBnd,  & ! (out)
    kl, ku, nz_1D,                          & ! (in)
    PROG_VARS0, DENS_hyd, PRES_hyd,         & ! (in)
    Dz, Lift, impl_fac, lmesh, elem,        & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y )

    use scale_linalgebra, only: linalgebra_LU
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(in) :: kl, ku, nz_1D
    real(RP), intent(out) :: PmatBnd(2*kl+ku+1,elem%Nnode_v,PROG_VARS_NUM,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y

    real(RP) :: RHOT_hyd(elem%Nnode_v)
    real(RP) :: DENS0(elem%Nnode_v)
    real(RP) :: POT0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: Cs0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: W0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: DPDRHOT0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    integer :: ke_z, ke_z2
    integer :: v, ke, p, f1, fp, FmV
    real(RP) :: gamm, rgamm
    real(RP) :: fac_dz_p(elem%Nnode_v)
    real(RP) :: PmatD(elem%Nnode_v,elem%Nnode_v,PROG_VARS_NUM,PROG_VARS_NUM)
    real(RP) :: PmatL(elem%Nnode_v,elem%Nnode_v,PROG_VARS_NUM,PROG_VARS_NUM)
    real(RP) :: PmatU(elem%Nnode_v,elem%Nnode_v,PROG_VARS_NUM,PROG_VARS_NUM)
    integer :: Colmask(elem%Nnode_v)
    real(RP) :: Id(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: tmp
    real(RP) :: fac

    integer :: ij, v1, v2, pv1, pv2,  g_kj, g_kjp1, g_kjm1, pb, pb1
    logical :: bc_flag
    logical :: eval_flag(PROG_VARS_NUM,PROG_VARS_NUM)
    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    eval_flag(:,:) = .false.
    do v=1, PROG_VARS_NUM
      eval_flag(v,v) = .true.
    end do
    eval_flag(DDENS_VID,MOMZ_VID) = .true.
    eval_flag(MOMZ_VID,DDENS_VID) = .true.
    eval_flag(MOMZ_VID,DRHOT_VID) = .true.
    eval_flag(DRHOT_VID,MOMZ_VID) = .true.
    eval_flag(DRHOT_VID,DDENS_VID) = .true.

    Id(:,:) = 0.0_RP
    do p=1, elem%Nnode_v
      Id(p,p) = 1.0_RP
    end do

    !$omp parallel private(RHOT_hyd, DENS0, Colmask)
    !$omp do
    do v=1, PROG_VARS_NUM
      PmatD(:,:,:,v) = 0.0_RP
      PmatL(:,:,:,v) = 0.0_RP
      PmatU(:,:,:,v) = 0.0_RP  
    end do
    !$omp do
    do ij=1, elem%Nnode_h1D**2
      PmatBnd(:,:,:,:,ij) = 0.0_RP
    end do
    !$omp do collapse(2)
    do ij=1, elem%Nnode_h1D**2
    do ke_z=1, lmesh%NeZ
      Colmask(:) = elem%Colmask(:,ij)
      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(Colmask(:),ke_z)/PRES00)**rgamm

      DPDRHOT0(:,ke_z,ij) = gamm * PRES_hyd(Colmask(:),ke_z) / RHOT_hyd(:)                      &
                  * ( 1.0_RP + PROG_VARS0(Colmask(:),DRHOT_VID,ke_z) / RHOT_hyd(:) )**(gamm-1) 

      DENS0(:) = DENS_hyd(Colmask(:),ke_z) + PROG_VARS0(Colmask(:),DDENS_VID,ke_z)
      POT0(:,ke_z,ij) = ( RHOT_hyd(:) + PROG_VARS0(Colmask(:),DRHOT_VID,ke_z) ) / DENS0(:)
      W0(:,ke_z,ij) = PROG_VARS0(Colmask(:),MOMZ_VID,ke_z) / DENS0(:)
      Cs0(:,ke_z,ij) = sqrt( gamm * PRES_hyd(Colmask(:),ke_z) * ( 1.0_RP + PROG_VARS0(Colmask(:),DRHOT_VID,ke_z) / RHOT_hyd(:) )**gamm / DENS0(:) )            
    end do
    end do
    !$omp end parallel

    !$omp parallel do private(ke_z, ke, ColMask, p, fp, v, f1, ke_z2, fac_dz_p, &
    !$omp fac, tmp,  FmV,                               &
    !$omp ij, v1, v2, pv1, pv2, pb1, g_kj, g_kjp1, g_kjm1, bc_flag ) &
    !$omp firstprivate(PmatD, PmatL, PmatU)
    do ij=1, elem%Nnode_h1D**2
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
      Colmask(:) = elem%Colmask(:,ij)

      !-----  
      do p=1, elem%Nnode_v
        fac_dz_p(:) = impl_fac * lmesh%Escale(Colmask(:),ke,3,3) * elem%Dx3(Colmask(:),Colmask(p)) 

        ! DDENS
        PmatD(:,p,DDENS_VID,DDENS_VID) = Id(:,p)
        PmatD(:,p,DDENS_VID,MOMZ_VID) = fac_dz_p(:) 

        ! MOMX
        PmatD(:,p,MOMX_VID,MOMX_VID) = Id(:,p)

        ! MOMY
        PmatD(:,p,MOMY_VID,MOMY_VID) = Id(:,p)

        ! MOMZ
        PmatD(:,p,MOMZ_VID,MOMZ_VID) = Id(:,p)
        PmatD(:,p,MOMZ_VID,DDENS_VID) = impl_fac * Grav * IntrpMat_VPOrdM1(Colmask(:),Colmask(p))
        PmatD(:,p,MOMZ_VID,DRHOT_VID) = fac_dz_p(:) * DPDRHOT0(p,ke_z,ij)

        !DRHOT
        PmatD(:,p,DRHOT_VID,DDENS_VID) = - fac_dz_p(:) * POT0(p,ke_z,ij) * W0(p,ke_z,ij)
        PmatD(:,p,DRHOT_VID,MOMZ_VID ) =   fac_dz_p(:) * POT0(p,ke_z,ij)
        PmatD(:,p,DRHOT_VID,DRHOT_VID) = Id(:,p) + fac_dz_p(:) * W0(p,ke_z,ij)
      end do

      do f1=1, 2
        if (f1==1) then
          ke_z2 = max(ke_z-1,1)
          pv1 = 1; pv2 = elem%Nnode_v
        else
          ke_z2 = min(ke_z+1,lmesh%NeZ)
          pv1 = elem%Nnode_v; pv2 = 1
        end if
        fac  = 0.5_RP * impl_fac
        if ( (ke_z == 1 .and. f1==1) .or. (ke_z == lmesh%NeZ .and. f1==elem%Nfaces_v) ) then
          bc_flag = .true.
          pv2 = pv1
        else 
          bc_flag = .false.    
        end if

        FmV  = elem%Fmask_v(ij,f1)
        fp = elem%Nfp_h * elem%Nfaces_h + (f1-1)*elem%Nfp_v + ij
        
        !--
        tmp = fac * elem%Lift(FmV,fp) * lmesh%Fscale(fp,ke) * &
                        max( abs( W0(pv1,ke_z ,ij) ) + Cs0(pv1,ke_z ,ij), &
                             abs( W0(pv2,ke_z2,ij) ) + Cs0(pv2,ke_z2,ij)  )
                  
        if (bc_flag) then
          PmatD(pv1,pv1,MOMZ_VID,MOMZ_VID) = PmatD(pv1,pv1,MOMZ_VID,MOMZ_VID) + 2.0_RP * tmp
        else 
          do v=1, PROG_VARS_NUM
            PmatD(pv1,pv1,v,v) = PmatD(pv1,pv1,v,v) + tmp
            if (f1 == 1) then
              PmatL(pv1,pv2,v,v) = - tmp                                
            else
              PmatU(pv1,pv2,v,v) = - tmp
            end if
          end do
        end if 

        !--
        tmp = fac * elem%Lift(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)

        if (bc_flag) then
          PmatD(pv1,pv1,DDENS_VID,MOMZ_VID) = PmatD(pv1,pv1,DDENS_VID,MOMZ_VID) - 2.0_RP * tmp
          PmatD(pv1,pv1,DRHOT_VID,MOMZ_VID) = PmatD(pv1,pv1,DRHOT_VID,MOMZ_VID) - 2.0_RP * tmp * POT0(pv1,ke_z,ij)

        else 
          PmatD(pv1,pv1,DDENS_VID,MOMZ_VID ) = PmatD(pv1,pv1,DDENS_VID,MOMZ_VID ) - tmp
          PmatD(pv1,pv1,MOMZ_VID ,DRHOT_VID) = PmatD(pv1,pv1,MOMZ_VID ,DRHOT_VID) - tmp * DPDRHOT0(pv1,ke_z,ij)
          PmatD(pv1,pv1,DRHOT_VID,MOMZ_VID ) = PmatD(pv1,pv1,DRHOT_VID,MOMZ_VID ) - tmp * POT0(pv1,ke_z,ij)
          PmatD(pv1,pv1,DRHOT_VID,DDENS_VID) = PmatD(pv1,pv1,DRHOT_VID,DDENS_VID) + tmp * POT0(pv1,ke_z,ij) * W0(pv1,ke_z,ij)
          PmatD(pv1,pv1,DRHOT_VID,DRHOT_VID) = PmatD(pv1,pv1,DRHOT_VID,DRHOT_VID) - tmp * W0(pv1,ke_z,ij)

          if (f1 == 1) then
            PmatL(pv1,pv2,DDENS_VID,MOMZ_VID ) = + tmp
            PmatL(pv1,pv2,MOMZ_VID,DRHOT_VID ) = + tmp * DPDRHOT0(pv2,ke_z2,ij) 
            PmatL(pv1,pv2,DRHOT_VID,MOMZ_VID ) = + tmp * POT0(pv2,ke_z2,ij)
            PmatL(pv1,pv2,DRHOT_VID,DDENS_VID) = - tmp * POT0(pv2,ke_z2,ij) * W0(pv2,ke_z2,ij)
            PmatL(pv1,pv2,DRHOT_VID,DRHOT_VID) = PmatL(pv1,pv2,DRHOT_VID,DRHOT_VID) &
                                                 + tmp * W0(pv2,ke_z2,ij)
          else
            PmatU(pv1,pv2,DDENS_VID,MOMZ_VID ) = + tmp
            PmatU(pv1,pv2,MOMZ_VID,DRHOT_VID ) = + tmp * DPDRHOT0(pv2,ke_z2,ij)     
            PmatU(pv1,pv2,DRHOT_VID,MOMZ_VID ) = + tmp * POT0(pv2,ke_z2,ij)
            PmatU(pv1,pv2,DRHOT_VID,DDENS_VID) = - tmp * POT0(pv2,ke_z2,ij) * W0(pv2,ke_z2,ij)
            PmatU(pv1,pv2,DRHOT_VID,DRHOT_VID) = PmatU(pv1,pv2,DRHOT_VID,DRHOT_VID) &
                                                 + tmp * W0(pv2,ke_z2,ij)
          end if
        end if
      end do

      do v2=1, PROG_VARS_NUM
      do v1=1, PROG_VARS_NUM
        if ( eval_flag(v1,v2) ) then
          do pv2=1, elem%Nnode_v
            g_kj   = pv2 + (v2-1)*elem%Nnode_v + (ke_z-1)*elem%Nnode_v*PROG_VARS_NUM
            g_kjm1 = pv2 + (v2-1)*elem%Nnode_v + (ke_z-2)*elem%Nnode_v*PROG_VARS_NUM
            g_kjp1 = pv2 + (v2-1)*elem%Nnode_v + (ke_z  )*elem%Nnode_v*PROG_VARS_NUM

            do pv1=1, elem%Nnode_v            
              pb1 = pv1 + (v1-1)*elem%Nnode_v + (ke_z-1)*elem%Nnode_v*PROG_VARS_NUM
              if (ke_z > 1) then
                PmatBnd(kl+ku+1+pb1-g_kjm1, pv2,v2,ke_z-1, ij) = PmatL(pv1,pv2,v1,v2)
              end if
              PmatBnd(kl+ku+1+pb1-g_kj, pv2,v2,ke_z, ij) = PmatD(pv1,pv2,v1,v2)
              if (ke_z < lmesh%NeZ) then
                PmatBnd(kl+ku+1+pb1-g_kjp1, pv2,v2,ke_z+1, ij) = PmatU(pv1,pv2,v1,v2)
              end if
            end do
          end do            
        end if
      end do
      end do
      end do  
    end do  

    return
  end subroutine vi_construct_matbnd

end module scale_atm_dyn_nonhydro3d_hevi
