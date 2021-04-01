!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics HEVE 
!!
!! @par Description
!!      HEVE DGM scheme for Atmospheric dynamical process. 
!!      To improve the numerical instability due to the aliasing errors, 
!!      the split form based on Gassner et al. (2016, JCP) is used for advection terms. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_splitform_heve
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
  use scale_element_modalfilter, only: ModalFilter
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
  public :: atm_dyn_dgm_nonhydro3d_heve_splitform_Init
  public :: atm_dyn_dgm_nonhydro3d_heve_splitform_Final
  public :: atm_dyn_dgm_nonhydro3d_heve_splitform_cal_tend

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

  real(RP), private, allocatable :: IntrpMat_VPOrdM1(:,:)
  real(RP), private, allocatable :: DxT1D_(:,:)
  real(RP), private, allocatable :: DyT1D_(:,:)
  real(RP), private, allocatable :: DzT1D_(:,:)

  private :: cal_del_flux_dyn
  private :: dx_ab, dy_ab, dz_ab
  private :: dx_abc, dy_abc, dz_abc

contains
  subroutine atm_dyn_dgm_nonhydro3d_heve_splitform_Init( mesh )

    implicit none
    class(MeshBase3D), intent(in) :: mesh

    integer :: p1, p2, p3, p_, p_intrp
    type(ElementBase3D), pointer :: elem
    real(RP) :: invV_POrdM1(mesh%refElem3D%Np,mesh%refElem3D%Np)

    real(RP), allocatable :: InvV_(:,:)
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

    allocate( DxT1D_(elem%Nnode_h1D,elem%Nnode_h1D) )
    allocate( DyT1D_(elem%Nnode_h1D,elem%Nnode_h1D) )
    allocate( DzT1D_(elem%Nnode_v,elem%Nnode_v) )

    do p1=1, elem%Nnode_h1D
      DxT1D_(:,p1) = elem%Dx1(p1,1:elem%Nnode_h1D)
    end do

    do p2=1, elem%Nnode_h1D
    do p_=1, elem%Nnode_h1D
      DyT1D_(p_,p2) = elem%Dx2(1+(p2-1)*elem%Nnode_h1D,1+(p_-1)*elem%Nnode_h1D)
    end do
    end do

    do p3=1, elem%Nnode_v
      DzT1D_(:,p3) = elem%Dx3(elem%Colmask(p3,1),elem%Colmask(:,1))
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_heve_splitform_Init


  subroutine atm_dyn_dgm_nonhydro3d_heve_splitform_Final()
    implicit none
    !--------------------------------------------

    deallocate( IntrpMat_VPOrdM1 )
    deallocate( DxT1D_, DyT1D_, DzT1D_ )

    return
  end subroutine atm_dyn_dgm_nonhydro3d_heve_splitform_Final  

  !-------------------------------

  subroutine atm_dyn_dgm_nonhydro3d_heve_splitform_cal_tend( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, CORIOLIS,          & ! (in)
    SL_flag, wdamp_tau, wdamp_height,                                           & ! (in)
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

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: Fx_sp(elem%Np), Fy_sp(elem%Np), Fz_sp(elem%Np)    
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: dens_(elem%Np), RHOT_hyd(elem%Np), RHOT_(elem%Np)
    real(RP) :: pres_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), pot_(elem%Np)
    real(RP) :: Cori(elem%Np)
    real(RP) :: drho(elem%Np)

    real(RP) :: tmp(elem%Np)
    integer :: ke, ke2d
    real(RP) :: gamm, rgamm
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
    gamm = CPDry / CvDry
    rgamm = CvDry / CpDry

    !$omp parallel do private( &
    !$omp RHOT_hyd, RHOT_, pres_, dens_, drho, u_, v_, w_, pot_, ke2d, Cori, &
    !$omp Fx, Fy, Fz, Fx_sp, Fy_sp, Fz_sp, LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE
      !--
      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**rgamm

      RHOT_(:) = RHOT_hyd(:) + DRHOT_(:,ke)
      pres_(:) = PRES_hyd(:,ke) * (1.0_RP + DRHOT_(:,ke)/RHOT_hyd(:))**gamm
      dens_(:) = DDENS_(:,ke) + DENS_hyd(:,ke)

      u_  (:) = MOMX_(:,ke) / dens_(:)
      v_  (:) = MOMY_(:,ke) / dens_(:)
      w_  (:) = MOMZ_(:,ke) / dens_(:)
      pot_(:) = RHOT_(:)    / dens_(:)

      ke2d = lmesh%EMap3Dto2D(ke)
      Cori(:) = CORIOLIS(elem%IndexH2Dto3D(:),ke2d)

      drho(:) = matmul(IntrpMat_VPOrdM1, DDENS_(:,ke))

      !-- DENS
      call dx_ab( DxT1D_, dens_, u_, elem%Nnode_h1D, elem%Nnode_v, Fx_sp )
      call dy_ab( DyT1D_, dens_, v_, elem%Nnode_h1D, elem%Nnode_v, Fy_sp )
      call dz_ab( DzT1D_, dens_, w_, elem%Nnode_h1D, elem%Nnode_v, Fz_sp )
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DDENS_VID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx_sp(:) &
          + lmesh%Escale(:,ke,2,2) * Fy_sp(:) &
          + lmesh%Escale(:,ke,3,3) * Fz_sp(:) &
          + LiftDelFlx(:) )
      
      !-- MOMX
      call dx_abc( DxT1D_, dens_, u_, u_, elem%Nnode_h1D, elem%Nnode_v, Fx_sp )
      call dy_abc( DyT1D_, dens_, u_, v_, elem%Nnode_h1D, elem%Nnode_v, Fy_sp )
      call dz_abc( DzT1D_, dens_, u_, w_, elem%Nnode_h1D, elem%Nnode_v, Fz_sp )
      call sparsemat_matmul(Dx, pres_(:)           , Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = - (  &
            lmesh%Escale(:,ke,1,1) * (Fx_sp(:) + Fx(:)) &
          + lmesh%Escale(:,ke,2,2) *  Fy_sp(:)          &
          + lmesh%Escale(:,ke,3,3) *  Fz_sp(:)          &
          + LiftDelFlx(:)                               &
          - Cori(:) * MOMY_(:,ke)                       &
          )

      !-- MOMY
      call dx_abc( DxT1D_, dens_, v_, u_, elem%Nnode_h1D, elem%Nnode_v, Fx_sp )
      call dy_abc( DyT1D_, dens_, v_, v_, elem%Nnode_h1D, elem%Nnode_v, Fy_sp )  
      call dz_abc( DzT1D_, dens_, v_, w_, elem%Nnode_h1D, elem%Nnode_v, Fz_sp )
      call sparsemat_matmul(Dy, pres_(:)           , Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = - (  &
            lmesh%Escale(:,ke,1,1) *  Fx_sp(:)          &
          + lmesh%Escale(:,ke,2,2) * (Fy_sp(:) + Fy(:)) &
          + lmesh%Escale(:,ke,3,3) *  Fz_sp(:)          &
          + LiftDelFlx(:)                               &
          + Cori(:) * MOMX_(:,ke)                       &
          ) 

      !-- MOMZ
      call dx_abc( DxT1D_, dens_, w_, u_, elem%Nnode_h1D, elem%Nnode_v, Fx_sp )
      call dy_abc( DyT1D_, dens_, w_, v_, elem%Nnode_h1D, elem%Nnode_v, Fy_sp )
      call dz_abc( DzT1D_, dens_, w_, w_, elem%Nnode_h1D, elem%Nnode_v, Fz_sp )
      call sparsemat_matmul(Dz, pres_(:) - PRES_hyd(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx_sp(:)            &
          + lmesh%Escale(:,ke,2,2) * Fy_sp(:)            &
          + lmesh%Escale(:,ke,3,3) * (Fz_sp(:) + Fz(:))  &
          + LiftDelFlx(:)                              ) &  
          - Grav * drho(:)

      !-- RHOT
      call dx_abc( DxT1D_, dens_, pot_,  u_, elem%Nnode_h1D, elem%Nnode_v, Fx_sp )
      call dy_abc( DyT1D_, dens_, pot_,  v_, elem%Nnode_h1D, elem%Nnode_v, Fy_sp )
      call dz_abc( DzT1D_, dens_, pot_,  w_, elem%Nnode_h1D, elem%Nnode_v, Fz_sp )        
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DRHOT_VID), LiftDelFlx)
      
      RHOT_dt(:,ke) =  - (  &
            lmesh%Escale(:,ke,1,1) * Fx_sp(:) &
          + lmesh%Escale(:,ke,2,2) * Fy_sp(:) &
          + lmesh%Escale(:,ke,3,3) * Fz_sp(:) &
          + LiftDelFlx(:)                     ) 
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 3)

    !- Sponge layer
    if (SL_flag) then
      call PROF_rapstart( 'cal_dyn_tend_sponge', 3)
      call atm_dyn_dgm_spongelayer_add_tend( MOMZ_dt, &
        MOMZ_, wdamp_tau, wdamp_height, lmesh, elem   )
      call PROF_rapend( 'cal_dyn_tend_sponge', 3)
    end if

    return
  end subroutine atm_dyn_dgm_nonhydro3d_heve_splitform_cal_tend

  !------

  subroutine cal_del_flux_dyn( del_flux, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,   &
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
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    real(RP) :: presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P
    real(RP) :: gamm, rgamm

    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry 

    !$omp parallel do private( iM, iP, alpha,  &
    !$omp VelP, VelM, presM, presP, dpresM, dpresP,           &
    !$omp densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P  & 
    !$omp )
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

      VelM = (MOMX_(iM)*nx(i) + MOMY_(iM)*ny(i) + MOMZ_(iM)*nz(i))/densM
      VelP = (MOMX_(iP)*nx(i) + MOMY_(iP)*ny(i) + MOMZ_(iP)*nz(i))/densP

      alpha = max( sqrt(gamm * presM / densM) + abs(VelM), &
                   sqrt(gamm * presP / densP) + abs(VelP)  )
            
      del_flux(i,DDENS_VID) = 0.5_RP * ( &
                    ( MOMX_(iP) - MOMX_(iM) ) * nx(i)   &
                  + ( MOMY_(iP) - MOMY_(iM) ) * ny(i)   &
                  + ( MOMZ_(iP) - MOMZ_(iM) ) * nz(i)   &
                  - alpha * ( DDENS_(iP) - DDENS_(iM) ) )
    
      del_flux(i,MOMX_VID) = 0.5_RP * ( &
                    ( MOMX_(iP) * VelP - MOMX_(iM) * VelM ) &
                  + ( dpresP - dpresM ) * nx(i)             &
                  - alpha * ( MOMX_(iP) - MOMX_(iM) )       )   

      del_flux(i,MOMY_VID) = 0.5_RP* ( &
                    ( MOMY_(iP) * VelP - MOMY_(iM) * VelM ) &
                  + ( dpresP - dpresM ) * ny(i)             &
                  - alpha * ( MOMY_(iP) - MOMY_(iM) )       )

      del_flux(i,MOMZ_VID) = 0.5_RP * ( &
                    ( MOMZ_(iP) * VelP - MOMZ_(iM) * VelM ) &
                  + ( dpresP - dpresM ) * nz(i)             &
                  - alpha * ( MOMZ_(iP) - MOMZ_(iM) )       )      
      
      del_flux(i,DRHOT_VID) = 0.5_RP * ( &
                  ( rhotP * VelP - rhotM * VelM )     &
                - alpha * ( DRHOT_(iP) - DRHOT_(iM) ) )

    end do

    return
  end subroutine cal_del_flux_dyn

!-----------------------------------------

  subroutine dx_ab(DxT1D, a, b, Nnode_h1D, Nnode_v, fx)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DxT1D(Nnode_h1D,Nnode_h1D)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fx(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fx(i,j,k) = 0.5_RP * sum( DxT1D(:,i) * (a(i,j,k) + a(:,j,k)) * (b(i,j,k) + b(:,j,k)) )
    end do
    end do
    end do

    return
  end subroutine dx_ab

  subroutine dy_ab(DyT1D, a, b, Nnode_h1D, Nnode_v, fy)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DyT1D(Nnode_h1D,Nnode_h1D)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fy(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fy(i,j,k) = 0.5_RP * sum( DyT1D(:,j) * (a(i,j,k) + a(i,:,k)) * (b(i,j,k) + b(i,:,k)) )
    end do
    end do
    end do

    return
  end subroutine dy_ab
  
  subroutine dz_ab(DzT1D, a, b, Nnode_h1D, Nnode_v, fz)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DzT1D(Nnode_v,Nnode_v)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fz(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fz(i,j,k) = 0.5_RP * sum( DzT1D(:,k) * (a(i,j,k) + a(i,j,:)) * (b(i,j,k) + b(i,j,:)) )
    end do
    end do
    end do

    return
  end subroutine dz_ab

  subroutine dx_abc(DxT1D, a, b, c, Nnode_h1D, Nnode_v, fx)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DxT1D(Nnode_h1D,Nnode_h1D)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: c(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fx(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fx(i,j,k) = 0.25_RP * sum( DxT1D(:,i) * (a(i,j,k) + a(:,j,k)) * (b(i,j,k) + b(:,j,k)) * (c(i,j,k) + c(:,j,k)) )
    end do
    end do
    end do

    return
  end subroutine dx_abc

  subroutine dy_abc(DyT1D, a, b, c, Nnode_h1D, Nnode_v, fy)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DyT1D(Nnode_h1D,Nnode_h1D)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: c(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fy(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fy(i,j,k) = 0.25_RP * sum( DyT1D(:,j) * (a(i,j,k) + a(i,:,k)) * (b(i,j,k) + b(i,:,k)) * (c(i,j,k) + c(i,:,k)) )
    end do
    end do
    end do

    return
  end subroutine dy_abc

  subroutine dz_abc(DzT1D, a, b, c, Nnode_h1D, Nnode_v, fz)
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(in) :: DzT1D(Nnode_v,Nnode_v)
    real(RP), intent(in) :: a(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: b(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(in) :: c(Nnode_h1D,Nnode_h1D,Nnode_v)
    real(RP), intent(out) :: fz(Nnode_h1D,Nnode_h1D,Nnode_v)

    integer :: i, j, k, p
    integer :: i2, j2, k2, p2
    !-------------------------------------------------

    do k=1, Nnode_v
    do j=1, Nnode_h1D    
    do i=1, Nnode_h1D
      fz(i,j,k) = 0.25_RP * sum( DzT1D(:,k) * (a(i,j,k) + a(i,j,:)) * (b(i,j,k) + b(i,j,:)) * (c(i,j,k) + c(i,j,:)) )
    end do
    end do
    end do

    return
  end subroutine dz_abc

end module scale_atm_dyn_dgm_nonhydro3d_splitform_heve
