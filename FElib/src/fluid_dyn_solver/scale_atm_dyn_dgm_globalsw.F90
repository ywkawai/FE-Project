!-------------------------------------------------------------------------------
!> module Atmosphere / global shallow water
!!
!! @par Description
!!      DGM scheme for Atmospheric dynamical process. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_globalsw
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_const, only: &
    GRAV => CONST_GRAV
  
  use scale_sparsemat
  use scale_element_base, only: &
    ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_mesh_base2d, only: MeshBase2D
  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_globalsw_Init
  public :: atm_dyn_dgm_globalsw_Final
  public :: atm_dyn_dgm_globalsw_cal_tend

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  integer, private, parameter :: VARS_H_ID      = 1
  integer, private, parameter :: VARS_U1_ID     = 2
  integer, private, parameter :: VARS_u2_ID     = 3
  integer, private, parameter :: PROG_VARS_NUM  = 3
  
  private :: cal_del_flux_dyn

contains
  subroutine atm_dyn_dgm_globalsw_Init( mesh )
    implicit none
    class(MeshBase2D), intent(in) :: mesh
    !--------------------------------------------
    return
  end subroutine atm_dyn_dgm_globalsw_Init


  subroutine atm_dyn_dgm_globalsw_Final()
    implicit none
    !--------------------------------------------
    return
  end subroutine atm_dyn_dgm_globalsw_Final  

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_globalsw_cal_tend( &
    h_dt, U_dt, V_dt,                       & ! (out)
    h_, U_, V_, hs_, u1_, u2_, CORIOLIS,    & ! (in)
    Dx, Dy, Sx, Sy, Lift, lmesh, elem  )

    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dy, Sx, Sy, Lift
    real(RP), intent(out) :: h_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: U_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: V_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: h_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: U_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: V_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: hs_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u1_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u2_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: del_flux_aux(elem%NfpTot,lmesh%Ne,1)
    real(RP) :: u1_dt(elem%Np), u2_dt(elem%Np)
    real(RP) :: VOR(elem%Np), E(elem%Np)
    integer :: ke       
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call cal_del_flux_dyn( del_flux, del_flux_aux,         & ! (out)
      h_, U_, V_, hs_, u1_, u2_,                           & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2),      & ! (in)
      lmesh%vmapM, lmesh%vmapP,                            & ! (in)
      lmesh, elem )                                          ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)

    !$omp parallel do private( &
    !$omp E, VOR,                          &
    !$omp u1_dt, u2_dt, Fx, Fy, LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE
      !--
      call sparsemat_matmul(Dx, u2_(:,ke), Fx)
      call sparsemat_matmul(Dy, u1_(:,ke), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_aux(:,ke,1), LiftDelFlx)

      VOR(:) = ( &
          lmesh%Escale(:,ke,1,1) * Fx(:) &
        - lmesh%Escale(:,ke,2,2) * Fy(:) &
        + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      E(:) = Grav * ( hs_(:,ke) + h_(:,ke) )                        &
           + 0.5_RP * ( u1_(:,ke) * U_(:,ke) + u2_(:,ke) * V_(:,ke) )
      
      !-- h
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * h_(:,ke) * U_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * h_(:,ke) * V_(:,ke), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_H_ID), LiftDelFlx)

      h_dt(:,ke) = - ( &
              lmesh%Escale(:,ke,1,1) * Fx(:) &
            + lmesh%Escale(:,ke,2,2) * Fy(:) &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      !-- u1
      call sparsemat_matmul(Dx, E(:), Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_u1_ID), LiftDelFlx)

      u1_dt(:) = - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:)                                &
          + LiftDelFlx(:)                                                 &
          - lmesh%Gsqrt(:,ke) * ( CORIOLIS(:,ke) + VOR(:) ) * V_(:,ke)    )
      
      !-- u2
      call sparsemat_matmul(Dy, E(:), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,VARS_u2_ID), LiftDelFlx)

      u2_dt(:) = - (  &
            lmesh%Escale(:,ke,2,2) * Fy(:)                                &
          + LiftDelFlx(:)                                                 &
          + lmesh%Gsqrt(:,ke) * ( CORIOLIS(:,ke) + VOR(:) ) * U_(:,ke)    )
    
      !--
      U_dt(:,ke) = lmesh%GIJ(:,ke,1,1) * u1_dt(:) + lmesh%GIJ(:,ke,1,2) * u2_dt(:) 
      V_dt(:,ke) = lmesh%GIJ(:,ke,2,1) * u1_dt(:) + lmesh%GIJ(:,ke,2,2) * u2_dt(:) 
    end do

    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_globalsw_cal_tend

  !------

!OCL SERIAL
  subroutine cal_del_flux_dyn( del_flux, del_flux_aux,     &
    h_, U_, V_, hs_, u1_, u2_,                             &
    Gsqrt_, G11_, G22_, nx, ny, vmapM, vmapP, lmesh, elem  )

    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP), intent(out) ::  del_flux_aux(elem%NfpTot,lmesh%Ne,1)
    real(RP), intent(in) ::  h_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  U_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  V_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  hs_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  u1_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  u2_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt_(elem%Np*lmesh%Ne) 
    real(RP), intent(in) ::  G11_(elem%Np*lmesh%Ne) 
    real(RP), intent(in) ::  G22_(elem%Np*lmesh%Ne) 
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    
    integer :: ke, iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: VelP(elem%NfpTot), VelM(elem%NfpTot), alpha(elem%NfpTot)
    real(RP) :: h_P(elem%NfpTot), h_M(elem%NfpTot)
    real(RP) :: E_P(elem%NfpTot), E_M(elem%NfpTot)
    real(RP) :: Gii(elem%NfpTot)
    !------------------------------------------------------------------------


    !$omp parallel do private( iM, iP, &
    !$omp alpha, VelM, VelP,           &
    !$omp h_P, h_M, E_P, E_M, Gii )
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)

      h_M(:) = h_(iM)
      h_P(:) = h_(iP)
      VelM(:) = U_(iM) * nx(:,ke) + V_(iM) * ny(:,ke)
      VelP(:) = U_(iP) * nx(:,ke) + V_(iP) * ny(:,ke)
      Gii(:) = G11_(iM) * abs(nx(:,ke)) + G22_(iM) * abs(ny(:,ke))

      E_M(:) = 0.5_RP * ( U_(iM) * u1_(iM) + V_(iM) * u2_(iM) ) &
             + Grav * ( h_M(:) + hs_(iM) )
      E_P(:) = 0.5_RP * ( U_(iP) * u1_(iP) + V_(iP) * u2_(iP) ) &
             + Grav * ( h_P(:) + hs_(iP) )
      
      alpha(:) = max( sqrt( Gii(:) * Grav * (h_M(:) + hs_(iM)) ) + abs(VelM(:)), &
                      sqrt( Gii(:) * Grav * (h_P(:) + hs_(iP)) ) + abs(VelP(:))  )
      
      del_flux(:,ke,VARS_h_ID) = 0.5_RP * Gsqrt_(iM) * (     &
                      (h_P(:) * VelP(:) - h_M(:) * VelM(:) ) &
                    - alpha(:) * ( h_P(:) - h_M(:) )         )

      del_flux(:,ke,VARS_u1_ID) = 0.5_RP * (               &
                    ( E_P(:)  - E_M(:) ) * nx(:,ke)        &
                    - alpha(:) * ( u1_(iP) - u1_(iM) )     )

      del_flux(:,ke,VARS_u2_ID) = 0.5_RP * (               &
                    ( E_P(:)  - E_M(:) ) * ny(:,ke)        &
                    - alpha(:) * ( u2_(iP) - u2_(iM) )     )
      
      del_flux_aux(:,ke,1) = 0.5_RP * ( &
                     ( u2_(iP) - u2_(iM) ) * nx(:,ke) &
                   - ( u1_(iP) - u1_(iM) ) * ny(:,ke) )
    end do

    return
  end subroutine cal_del_flux_dyn

end module scale_atm_dyn_dgm_globalsw
