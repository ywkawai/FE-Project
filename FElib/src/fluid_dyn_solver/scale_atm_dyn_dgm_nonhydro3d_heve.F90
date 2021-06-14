!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics HEVE
!!
!! @par Description
!!      HEVE DGM scheme for Atmospheric dynamical process. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_heve
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
  public :: atm_dyn_dgm_nonhydro3d_heve_Init
  public :: atm_dyn_dgm_nonhydro3d_heve_Final
  public :: atm_dyn_dgm_nonhydro3d_heve_cal_tend

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
  subroutine atm_dyn_dgm_nonhydro3d_heve_Init( mesh )

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

    !--
    return
  end subroutine atm_dyn_dgm_nonhydro3d_heve_Init


  subroutine atm_dyn_dgm_nonhydro3d_heve_Final()
    implicit none
    !--------------------------------------------
    
    deallocate( IntrpMat_VPOrdM1 )
    
    return
  end subroutine atm_dyn_dgm_nonhydro3d_heve_Final  

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_heve_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, CORIOLIS,          & ! (in)
    SL_flag, wdamp_tau, wdamp_height, hveldamp_flag,                            & ! (in)
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )

    use scale_atm_dyn_dgm_spongelayer, only: &
      atm_dyn_dgm_spongelayer_add_tend

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
    logical, intent(in) :: SL_flag
    real(RP), intent(in) :: wdamp_tau
    real(RP), intent(in) :: wdamp_height
    logical, intent(in) :: hveldamp_flag

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: PRES_(elem%Np)
    real(RP) :: RHOT_(elem%Np)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np)
    real(RP) :: Cori(elem%Np)
    real(RP) :: drho(elem%Np)

    integer :: ke, ke2d

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR       
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call cal_del_flux_dyn( del_flux,                                          & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem )                                                             ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private(                                &
    !$omp RHOT_, PRES_, rdens_, u_, v_, w_, ke2d, Cori, drho, &
    !$omp Fx, Fy, Fz, LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE
      !--
      RHOT_(:) = P0ovR * (PRES_hyd(:,ke) * rP0)**rgamm + DRHOT_(:,ke)
      PRES_(:) = PRES00 * (RovP0 * RHOT_(:))**gamm

      rdens_(:) = 1.0_RP / (DDENS_(:,ke) + DENS_hyd(:,ke))
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)

      ke2d = lmesh%EMap3Dto2D(ke)
      Cori(:) = CORIOLIS(elem%IndexH2Dto3D(:),ke2d)

      drho(:) = matmul(IntrpMat_VPOrdM1, DDENS_(:,ke))

      !-- DENS
      call sparsemat_matmul(Dx, MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, MOMY_(:,ke), Fy)
      call sparsemat_matmul(Dz, MOMZ_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_DDENS_ID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDelFlx(:) )
      
      !-- MOMX
      call sparsemat_matmul(Dx, u_(:) * MOMX_(:,ke) + PRES_(:), Fx)
      call sparsemat_matmul(Dy, v_(:) * MOMX_(:,ke)           , Fy)
      call sparsemat_matmul(Dz, w_(:) * MOMX_(:,ke)           , Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_MOMX_ID), LiftDelFlx)

      MOMX_dt(:,ke) = - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:)      &
          + lmesh%Escale(:,ke,2,2) * Fy(:)      &
          + lmesh%Escale(:,ke,3,3) * Fz(:)      &
          + LiftDelFlx(:)                       &
          - Cori(:) * MOMY_(:,ke)               &
          )

      !-- MOMY
      call sparsemat_matmul(Dx, u_(:) * MOMY_(:,ke)           , Fx)
      call sparsemat_matmul(Dy, v_(:) * MOMY_(:,ke) + PRES_(:), Fy)
      call sparsemat_matmul(Dz, w_(:) * MOMY_(:,ke)           , Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_MOMY_ID), LiftDelFlx)

      MOMY_dt(:,ke) = - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDelFlx(:)                  &
          + Cori(:) * MOMX_(:,ke)          &
          ) 

      !-- MOMZ
      call sparsemat_matmul(Dx, u_(:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, v_(:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, w_(:) * MOMZ_(:,ke) + PRES_(:) - PRES_hyd(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_MOMZ_ID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)   &
          + lmesh%Escale(:,ke,2,2) * Fy(:)   &
          + lmesh%Escale(:,ke,3,3) * Fz(:)   &
          + LiftDelFlx(:)                )   &
          - Grav * drho(:)

      !-- RHOT
      call sparsemat_matmul(Dx, u_(:) * RHOT_(:), Fx)
      call sparsemat_matmul(Dy, v_(:) * RHOT_(:), Fy)
      call sparsemat_matmul(Dz, w_(:) * RHOT_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_DRHOT_ID), LiftDelFlx)
      
      RHOT_dt(:,ke) =  - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDelFlx(:) ) 
    end do

    call PROF_rapend('cal_dyn_tend_interior', 3)

    !- Sponge layer
    if (SL_flag) then
      call PROF_rapstart('cal_dyn_tend_sponge', 3)
      call atm_dyn_dgm_spongelayer_add_tend( &
        MOMX_dt, MOMY_dt, MOMZ_dt,                    & ! (out)
        MOMX_, MOMY_, MOMZ_, wdamp_tau, wdamp_height, & ! (in)
        hveldamp_flag, lmesh, elem                    ) ! (in)
      call PROF_rapend('cal_dyn_tend_sponge', 3)
    end if

    return
  end subroutine atm_dyn_dgm_nonhydro3d_heve_cal_tend

  !------

!OCL SERIAL
  subroutine cal_del_flux_dyn( del_flux, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,   &
    nx, ny, nz, vmapM, vmapP, lmesh, elem                      )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
     real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    
    integer :: ke, i, iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: VelP(elem%NfpTot), VelM(elem%NfpTot), alpha(elem%NfpTot)
    real(RP) :: dpres(elem%NfpTot), presM(elem%NfpTot), presP(elem%NfpTot)
    real(RP) :: densM(elem%NfpTot), densP(elem%NfpTot)
    real(RP) :: rhotM(elem%NfpTot), rhotP(elem%NfpTot)
    real(RP) :: DDENS_P(elem%NfpTot), DDENS_M(elem%NfpTot)
    real(RP) :: MOMX_P(elem%NfpTot), MOMX_M(elem%NfpTot)
    real(RP) :: MOMY_P(elem%NfpTot), MOMY_M(elem%NfpTot)
    real(RP) :: MOMZ_P(elem%NfpTot), MOMZ_M(elem%NfpTot)
    real(RP) :: DRHOT_P(elem%NfpTot), DRHOT_M(elem%NfpTot)
    real(RP) :: PRES_hyd_P(elem%NfpTot), PRES_hyd_M(elem%NfpTot)

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR     
    !------------------------------------------------------------------------

    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( i, iM, iP, &
    !$omp alpha, VelM, VelP,                               &
    !$omp dpres, presM, presP, densM, densP, rhotM, rhotP, &
    !$omp MOMX_M, MOMX_P, MOMY_M, MOMY_P, MOMZ_M, MOMZ_P,  &
    !$omp DDENS_M, DDENS_P, DRHOT_M, DRHOT_P,              &
    !$omp PRES_hyd_M, PRES_hyd_P   )
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)

      DDENS_M(:) = DDENS_(iM)
      DDENS_P(:) = DDENS_(iP)
      MOMX_M(:) = MOMX_(iM)
      MOMX_P(:) = MOMX_(iP)
      MOMY_M(:) = MOMY_(iM)
      MOMY_P(:) = MOMY_(iP)
      MOMZ_M(:) = MOMZ_(iM)
      MOMZ_P(:) = MOMZ_(iP)
      DRHOT_M(:) = DRHOT_(iM)
      DRHOT_P(:) = DRHOT_(iP)
      PRES_hyd_M(:) = PRES_hyd(iM)
      PRES_hyd_P(:) = PRES_hyd(iP)

      densM(:) = DDENS_M(:) + DENS_hyd(iM)
      densP(:) = DDENS_P(:) + DENS_hyd(iP)

      VelM(:) = ( MOMX_M(:) * nx(:,ke) + MOMY_M(:) * ny(:,ke) + MOMZ_M(:) * nz(:,ke) ) / densM(:)
      VelP(:) = ( MOMX_P(:) * nx(:,ke) + MOMY_P(:) * ny(:,ke) + MOMZ_P(:) * nz(:,ke) ) / densP(:)

      rhotM(:) = P0ovR * (PRES_hyd_M(:) * rP0)**rgamm + DRHOT_M(:)
      rhotP(:) = P0ovR * (PRES_hyd_P(:) * rP0)**rgamm + DRHOT_P(:)

      presM(:) = PRES00 * (RovP0 * rhotM(:))**gamm
      presP(:) = PRES00 * (RovP0 * rhotP(:))**gamm

      dpres(:)  =  presP(:) - presM(:)                         &
             - ( PRES_hyd_P(:) - PRES_hyd_M(:) ) * abs(nz(:,ke))

      alpha(:) = max( sqrt(gamm * presM(:) / densM(:)) + abs(VelM(:)), &
                      sqrt(gamm * presP(:) / densP(:)) + abs(VelP(:))  )

      del_flux(:,ke,VARS_DDENS_ID) = 0.5_RP * (                  &
                    ( densP(:) * VelP(:) - densM(:) * VelM(:) )  &
                    - alpha(:) * ( DDENS_P(:) - DDENS_M(:) )     )

      del_flux(:,ke,VARS_MOMX_ID) = 0.5_RP * (                     &
                    ( MOMX_P(:) * VelP(:) - MOMX_M(:) * VelM(:) )  &
                    + dpres(:) * nx(:,ke)                          &
                    - alpha(:) * ( MOMX_P(:) - MOMX_M(:) )         )
      
      del_flux(:,ke,VARS_MOMY_ID) = 0.5_RP * (                     &
                    ( MOMY_P(:) * VelP(:) - MOMY_M(:) * VelM(:) )  &
                    + dpres(:) * ny(:,ke)                          &
                    - alpha(:) * ( MOMY_P(:) - MOMY_M(:) )         )               
      
      del_flux(:,ke,VARS_MOMZ_ID) = 0.5_RP * (                     &
                    ( MOMZ_P(:) * VelP(:) - MOMZ_M(:) * VelM(:))   &
                    + dpres(:) * nz(:,ke)                          &                   
                    - alpha(:) * ( MOMZ_P(:) - MOMZ_M(:) )         )
      
      del_flux(:,ke,VARS_DRHOT_ID) = 0.5_RP * (                    &
                    ( rhotP(:) * VelP(:) - rhotM(:) * VelM(:) )    &
                    - alpha(:) * ( DRHOT_P(:) - DRHOT_M(:) )       )
    end do

    return
  end subroutine cal_del_flux_dyn

end module scale_atm_dyn_dgm_nonhydro3d_heve
