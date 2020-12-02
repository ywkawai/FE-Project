!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics HEVI 
!!
!! @par Description
!!      HEVE DGM scheme for Atmospheric dynamical process. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_nonhydro3d_heve
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
  public :: atm_dyn_nonhydro3d_heve_Init
  public :: atm_dyn_nonhydro3d_heve_Final
  public :: atm_dyn_nonhydro3d_heve_cal_tend

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
  
  real(RP), private, allocatable :: IntrpMat_VPOrdM1(:,:)

  private :: cal_del_flux_dyn

contains
  subroutine atm_dyn_nonhydro3d_heve_Init( mesh )

    implicit none
    class(MeshBase3D), intent(in) :: mesh

    integer :: p1, p2, p_
    type(ElementBase3D), pointer :: elem
    real(RP) :: invV_VPOrdM1(mesh%refElem3D%Np,mesh%refElem3D%Np)
    !--------------------------------------------

    elem => mesh%refElem3D
    allocate( IntrpMat_VPOrdM1(elem%Np,elem%Np) )
    
    InvV_VPOrdM1(:,:) = elem%invV
    do p2=1, elem%Nnode_h1D
    do p1=1, elem%Nnode_h1D
      p_ = p1 + (p2-1)*elem%Nnode_h1D + (elem%Nnode_v-1)*elem%Nnode_h1D**2
      InvV_VPOrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem%V, invV_VPOrdM1)

    return
  end subroutine atm_dyn_nonhydro3d_heve_Init


  subroutine atm_dyn_nonhydro3d_heve_Final()
    implicit none
    !--------------------------------------------
    
    deallocate( IntrpMat_VPOrdM1 )
    
    return
  end subroutine atm_dyn_nonhydro3d_heve_Final  

  !-------------------------------

  subroutine atm_dyn_nonhydro3d_heve_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, CORIOLIS,          & ! (in)
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

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: dens_(elem%Np), RHOT_(elem%Np), dpres_(elem%Np)
    real(RP) :: pres_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np)
    real(RP) :: Cori(elem%Np)

    integer :: ke, ke2d
    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 3)
    call cal_del_flux_dyn( del_flux,                                          & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem )                                                             ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 3)
    !$omp parallel do private(RHOT_,pres_,dpres_,dens_,u_,v_,w_,ke2d,Cori,Fx,Fy,Fz,LiftDelFlx)
    do ke = lmesh%NeS, lmesh%NeE
      !--

      RHOT_(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
      pres_(:) = PRES00 * (Rdry*RHOT_(:)/PRES00)**(CPdry/Cvdry)
      dpres_(:) = pres_(:) - PRES_hyd(:,ke)
      dens_(:) = DDENS_(:,ke) + DENS_hyd(:,ke)

      u_(:) = MOMX_(:,ke)/dens_(:)
      v_(:) = MOMY_(:,ke)/dens_(:)
      w_(:) = MOMZ_(:,ke)/dens_(:)

      ke2d = lmesh%EMap3Dto2D(ke)
      Cori(:) = CORIOLIS(elem%IndexH2Dto3D(:),ke2d)

      !-- DENS
      call sparsemat_matmul(Dx, MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, MOMY_(:,ke), Fy)
      call sparsemat_matmul(Dz, MOMZ_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_DDENS_ID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDelFlx(:) )
      
      !-- MOMX
      call sparsemat_matmul(Dx, u_(:)*MOMX_(:,ke) + pres_(:), Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMX_(:,ke)           , Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMX_(:,ke)           , Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_MOMX_ID), LiftDelFlx)

      MOMX_dt(:,ke) = - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:)      &
          + lmesh%Escale(:,ke,2,2) * Fy(:)      &
          + lmesh%Escale(:,ke,3,3) * Fz(:)      &
          + LiftDelFlx(:)                       &
          - Cori(:)*MOMY_(:,ke)                 &
          )

      !-- MOMY
      call sparsemat_matmul(Dx, u_(:)*MOMY_(:,ke)           , Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMY_(:,ke) + pres_(:), Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMY_(:,ke)           , Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_MOMY_ID), LiftDelFlx)

      MOMY_dt(:,ke) = - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDelFlx(:)                  &
          + Cori(:)*MOMX_(:,ke)            &
          )

      !-- MOMZ
      call sparsemat_matmul(Dx, u_(:)*MOMZ_(:,ke)            , Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMZ_(:,ke)            , Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMZ_(:,ke) + dpres_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_MOMZ_ID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)   &
          + lmesh%Escale(:,ke,2,2) * Fy(:)   &
          + lmesh%Escale(:,ke,3,3) * Fz(:)   &
          + LiftDelFlx(:)                )   &
          - matmul(IntrpMat_VPOrdM1, DDENS_(:,ke)) * Grav      

      !-- RHOT
      call sparsemat_matmul(Dx, u_(:)*RHOT_(:), Fx)
      call sparsemat_matmul(Dy, v_(:)*RHOT_(:), Fy)
      call sparsemat_matmul(Dz, w_(:)*RHOT_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_DRHOT_ID), LiftDelFlx)
      
      RHOT_dt(:,ke) =  - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDelFlx(:) )

    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_nonhydro3d_heve_cal_tend

  !------

  subroutine cal_del_flux_dyn( del_flux, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,   &
    nx, ny, nz, vmapM, vmapP, lmesh, elem                      )

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
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    real(RP) :: uM, uP, vM, vP, wM, wP, presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P
    real(RP) :: gamm, rgamm
    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    !$omp parallel do private( iM, iP, alpha, &
    !$omp uM, uP, vM, vP, wM, wP, VelP, VelM,                 &
    !$omp presM, presP, dpresM, dpresP,                       &
    !$omp densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P  )
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      rhot_hyd_M = PRES00/Rdry * (PRES_hyd(iM)/PRES00)**rgamm
      rhot_hyd_P = PRES00/Rdry * (PRES_hyd(iP)/PRES00)**rgamm
      
      rhotM = rhot_hyd_M + DRHOT_(iM)
      presM = PRES00 * (Rdry*rhotM/PRES00)**gamm
      dpresM = presM - PRES_hyd(iM)*abs(nz(i))

      rhotP = rhot_hyd_P + DRHOT_(iP) 
      presP = PRES00 * (Rdry*rhotP/PRES00)**gamm
      dpresP = presP - PRES_hyd(iP)*abs(nz(i))

      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iP)

      VelM = (MOMX_(iM)*nx(i) + MOMY_(iM)*ny(i) + MOMZ_(iM)*nz(i))/densM
      VelP = (MOMX_(iP)*nx(i) + MOMY_(iP)*ny(i) + MOMZ_(iP)*nz(i))/densP

      alpha = max( sqrt(gamm*presM/densM) + abs(VelM), sqrt(gamm*presP/densP) + abs(VelP)  )
      
      del_flux(i,VARS_DDENS_ID) = 0.5_RP*(               &
                    ( densP*VelP - densM*VelM )          &
                    - alpha*(DDENS_(iP) - DDENS_(iM))   )
      
      del_flux(i,VARS_MOMX_ID) = 0.5_RP*(                &
                    ( MOMX_(iP)*VelP - MOMX_(iM)*VelM )  &
                    + ( dpresP - dpresM )*nx(i)          &
                    - alpha*(MOMX_(iP) - MOMX_(iM))      )
      
      del_flux(i,VARS_MOMY_ID) = 0.5_RP*(                &
                    ( MOMY_(iP)*VelP - MOMY_(iM)*VelM )  &
                    + ( dpresP - dpresM )*ny(i)          &
                    - alpha*(MOMY_(iP) - MOMY_(iM))      )               
      
      del_flux(i,VARS_MOMZ_ID) = 0.5_RP*(                &
                    ( MOMZ_(iP)*VelP - MOMZ_(iM)*VelM)   &
                    + ( dpresP - dpresM )*nz(i)          &                    
                    - alpha*(MOMZ_(iP) - MOMZ_(iM))      )
      
      del_flux(i,VARS_DRHOT_ID) = 0.5_RP*(               &
                      ( rhotP*VelP - rhotM*VelM )        &
                    - alpha*(DRHOT_(iP) - DRHOT_(iM))    )
    end do

    return
  end subroutine cal_del_flux_dyn

  !------------------------------------------------
end module scale_atm_dyn_nonhydro3d_heve
