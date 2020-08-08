!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_nonhydro3d_hevi_gmres
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

  use scale_gmres, only: GMRES

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
  public :: atm_dyn_nonhydro3d_hevi_prepair_expfilter
  public :: atm_dyn_nonhydro3d_hevi_cal_tend
  public :: atm_dyn_nonhydro3d_hevi_cal_grad_diffVars
  public :: atm_dyn_nonhydro3d_hevi_filter_prgvar
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

  real(RP), private, allocatable :: FilterMat(:,:)
  real(RP), private, allocatable :: IntrpMat_VPOrdM1(:,:)
  real(RP), private, allocatable :: IntrpMat_YPOrdM1(:,:)

  private :: cal_del_flux_dyn
  private :: cal_del_gradDiffVar

contains
  subroutine atm_dyn_nonhydro3d_hevi_Init( mesh )

    implicit none
    class(MeshBase3D), intent(in) :: mesh

    integer :: p1, p2, p3, p_
    type(ElementBase3D), pointer :: elem
    real(RP) :: invV_POrdM1(mesh%refElem3D%Np,mesh%refElem3D%Np)
    !--------------------------------------------

    elem => mesh%refElem3D
    allocate( IntrpMat_VPOrdM1(elem%Np,elem%Np) )
    allocate( IntrpMat_YPOrdM1(elem%Np,elem%Np) )
    
    invV_POrdM1(:,:) = elem%invV
    do p2=1, elem%Nnode_h1D
    do p1=1, elem%Nnode_h1D
      p_ = p1 + (p2-1)*elem%Nnode_h1D + (elem%Nnode_v-1)*elem%Nnode_h1D**2
      invV_POrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem%V, invV_POrdM1)

    invV_POrdM1(:,:) = elem%invV
    do p3=1, elem%Nnode_v
    do p1=1, elem%Nnode_h1D
      p_ = p1 + (elem%Nnode_h1D-1)*elem%Nnode_h1D + (p3-1)*elem%Nnode_h1D**2
      invV_POrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_YPOrdM1(:,:) = matmul(elem%V, invV_POrdM1)

    return
  end subroutine atm_dyn_nonhydro3d_hevi_Init

  subroutine atm_dyn_nonhydro3d_hevi_prepair_expfilter(  &
    elem,                                           &
    etac, alpha, ord )

    implicit none
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(in) :: etac
    real(RP), intent(in) :: alpha
    integer, intent(in) :: ord

    real(RP) :: filter1D_h(elem%Nnode_h1D)
    real(RP) :: filter1D_v(elem%Nnode_v)
    real(RP) :: eta
    integer :: p1, p2, p3
    integer :: l
    !----------------------------------------------------

    filter1D_h(:) = 1.0_RP
    do p1=1, elem%Nnode_h1D
      eta = dble(p1-1)/dble(elem%PolyOrder_h)
      if ( eta >  etac .and. p1 /= 1) then
        filter1D_h(p1) = exp( -  alpha*( ((eta - etac)/(1.0_RP - etac))**ord ))
      end if
    end do

    filter1D_v(:) = 1.0_RP
    do p3=1, elem%Nnode_v
      eta = dble(p3-1)/dble(elem%PolyOrder_v)
      if ( eta >  etac .and. p3 /= 1) then
        filter1D_v(p3) = exp( -  alpha*( ((eta - etac)/(1.0_RP - etac))**ord ))
      end if
    end do

    allocate( FilterMat(elem%Np,elem%Np) )
    FilterMat(:,:) = 0.0_RP
    do p3=1, elem%Nnode_v
    do p2=1, elem%Nnode_h1D
    do p1=1, elem%Nnode_h1D
      l = p1 + (p2-1)*elem%Nnode_h1D + (p3-1)*elem%Nnode_h1D**2
      FilterMat(l,l) = filter1D_h(p1) * filter1D_h(p2) * filter1D_v(p3)
    end do  
    end do
    end do
    FilterMat(:,:) = matmul(FilterMat, elem%invV)
    FilterMat(:,:) = matmul(elem%V, FilterMat)
    
    return
  end subroutine atm_dyn_nonhydro3d_hevi_prepair_expfilter

  subroutine atm_dyn_nonhydro3d_hevi_Final()
    implicit none
    !--------------------------------------------

    deallocate( IntrpMat_YPOrdM1 )    
    deallocate( IntrpMat_VPOrdM1 )
    if( allocated(FilterMat) ) deallocate( FilterMat )
    
    return
  end subroutine atm_dyn_nonhydro3d_hevi_Final  

  !-------------------------------

  subroutine atm_dyn_nonhydro3d_hevi_filter_prgvar( &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, lmesh, elem  )
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem    
    real(RP), intent(inout)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: DRHOT_(elem%Np,lmesh%NeA)
    
    integer :: k
    !------------------------------------

    !$omp parallel do
    do k=1, lmesh%Ne
      DDENS_(:,k) = matmul(FilterMat,DDENS_(:,k))
      MOMX_(:,k) = matmul(FilterMat,MOMX_(:,k))
      MOMY_(:,k) = matmul(FilterMat,MOMY_(:,k))
      MOMZ_(:,k) = matmul(FilterMat,MOMZ_(:,k))
      DRHOT_(:,k) = matmul(FilterMat,DRHOT_(:,k))
    end do
    
    return
  end subroutine atm_dyn_nonhydro3d_hevi_filter_prgvar

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

    !$omp parallel do private(RHOT_hyd,RHOT_,pres_,dpres_,dens_,u_,v_,w_,pot_,ke2d,Cori,Fx,Fy,Fz,LiftDelFlx,tmp)
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
      
      alpha = swV*max( sqrt(gamm*presM/densM) + abs(VelM), &
                       sqrt(gamm*presP/densP) + abs(VelP)  )
          
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

      del_flux(i,MOMZ_VID) = 0.5_RP*(                     &
                      ( MOMZ_P*VelP - MOMZ_(iM)*VelM)     &
                    - alpha * (MOMZ_(iP) - MOMZ_(iM))     &
                    - dDiffFluxW                   )
      
      del_flux(i,DRHOT_VID) = 0.5_RP*(                    &
                      swV*( rhotP*VelP - rhotM*VelM )     &
                    - alpha *(DRHOT_(iP) - DRHOT_(iM))    &
                    - dDiffFluxPT                  )
    end do

    return
  end subroutine cal_del_flux_dyn

  subroutine atm_dyn_nonhydro3d_hevi_cal_grad_diffVars( &
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
      GyU_(:,ke) = lmesh%Escale(:,ke,2,2)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, U_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GzU_ID), LiftDelFlx)
      GzU_(:,ke) = lmesh%Escale(:,ke,3,3)*Fz(:) + LiftDelFlx(:)

      !- V

      call sparsemat_matmul(Dx, V_, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GxV_ID), LiftDelFlx)
      GxV_(:,ke) = lmesh%Escale(:,ke,1,1)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dy, V_, Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GyV_ID), LiftDelFlx)
      GyV_(:,ke) = lmesh%Escale(:,ke,2,2)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, V_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GzV_ID), LiftDelFlx)
      GzV_(:,ke) = lmesh%Escale(:,ke,3,3)*Fz(:) + LiftDelFlx(:)

      !- W

      call sparsemat_matmul(Dx, W_, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GxW_ID), LiftDelFlx)
      GxW_(:,ke) = lmesh%Escale(:,ke,1,1)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dy, W_, Fx)
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
      GyPT_(:,ke) = lmesh%Escale(:,ke,2,2)*Fz(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, DTHETA_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_GzPT_ID), LiftDelFlx)
      GzPT_(:,ke) = lmesh%Escale(:,ke,3,3)*Fz(:) + LiftDelFlx(:)

    end do

    return
  end subroutine atm_dyn_nonhydro3d_hevi_cal_grad_diffVars

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
    real(RP) :: PROG_DEL(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: b(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: tend(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: DENS_hyd_z(elem%Np,lmesh%NeZ)
    real(RP) :: PRES_hyd_z(elem%Np,lmesh%NeZ)
    real(RP) :: nz(elem%NfpTot,lmesh%NeZ)
    integer :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer :: ke_x, ke_y, ke_z, ke, p, v
    integer :: itr_lin, itr_nlin
    integer :: m, N
    integer :: f, vs, ve
    logical :: is_converged


    type(GMRES) :: gmres_hevi
    real(RP), allocatable :: wj(:), pinv_v(:)
    real(RP) :: PmatDlu(elem%Np*PROG_VARS_NUM,elem%Np*PROG_VARS_NUM,lmesh%NeZ)
    integer :: PmatDlu_ipiv(elem%Np*PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: PmatL(elem%Np,elem%Np,PROG_VARS_NUM,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: PmatU(elem%Np,elem%Np,PROG_VARS_NUM,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), parameter :: EPS0 = 1.0E-16_RP
    real(RP), parameter :: EPS = 1.0E-16_RP

    !------------------------------------------------------------------------

    call PROF_rapstart( 'hevi_cal_vi', 3)    
    call PROF_rapstart( 'hevi_cal_vi_prep', 3)

    N = elem%Np * PROG_VARS_NUM * lmesh%NeZ
    m = N / PROG_VARS_NUM / elem%Nnode_h1D**2

    call gmres_hevi%Init( N, m, EPS, EPS )
    allocate( wj(N), pinv_v(N) )

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

    do ke_z=1, lmesh%NeZ
      vs = elem%Nfp_h*elem%Nfaces_h + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z > 1) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,2) + (ke_z-2)*elem%Np
      ! if (ke_z > 1) then
      !   vmapP(vs:ve,ke_z) = elem%Fmask_v(:,2) + (ke_z-2)*elem%Np
      ! else
      !   vmapP(vs:ve,ke_z) = elem%Fmask_v(:,2) + (lmesh%NeZ-1)*elem%Np
      ! end if

      vs = elem%Nfp_h*elem%Nfaces_h + elem%Nfp_v + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z < lmesh%NeZ) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1) + ke_z*elem%Np
      ! if (ke_z < lmesh%NeZ) then
      !   vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1) + ke_z*elem%Np
      ! else
      !   vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1)
      ! end if      
    end do
    call PROF_rapend( 'hevi_cal_vi_prep', 3)

    do ke_y=1, lmesh%NeY
    do ke_x=1, lmesh%NeX

      call PROF_rapstart( 'hevi_cal_vi_get_var', 3)

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
      
      if ( impl_fac > 0.0_RP ) then
        call PROF_rapstart( 'hevi_cal_vi_itr', 3)

        do itr_nlin = 1, 1

          do ke_z=1, lmesh%NeZ
            PROG_DEL(:,:,ke_z) = 0.0_RP
          end do

          call vi_eval_Ax( Ax(:,:,:),                       & ! (out)
            PROG_VARS, PROG_VARS0, DENS_hyd_z, PRES_hyd_z,  & ! (in)
            Dz, Lift, impl_fac, lmesh, elem,                & ! (in)
            nz, vmapM, vmapP, ke_x, ke_y, .false. )

          do ke_z=1, lmesh%NeZ
            b(:,:,ke_z) = - Ax(:,:,ke_z) + PROG_VARS00(:,:,ke_z)
          end do

          call PROF_rapstart( 'hevi_cal_vi_pmatinv', 3)

          call vi_construct_pmatInv( PmatDlu, PmatDlu_ipiv, PmatL, PmatU,  & ! (out)
            PROG_VARS0, DENS_hyd_z, PRES_hyd_z,             & ! (in)
            Dz, Lift, impl_fac, lmesh, elem,                & ! (in)
            nz, vmapM, vmapP, ke_x, ke_y  )

          call PROF_rapend( 'hevi_cal_vi_pmatinv', 3)

          call PROF_rapstart( 'hevi_cal_vi_itr_lin', 3)            
          do itr_lin=1, 2*int(N/m)

            call vi_GMRES_core( gmres_hevi, PROG_DEL(:,:,:), wj, & ! (inout)
              is_converged,                                      & ! (out)
              PROG_VARS(:,:,:), b(:,:,:), N, m, PmatDlu,  PmatDlu_ipiv,  PmatL, PmatU, pinv_v, & ! (in)
              DENS_hyd_z, PRES_hyd_z,                            & ! (in)
              Dz, Lift, impl_fac, lmesh, elem,                   & ! (in)
              nz, vmapM, vmapP, ke_x, ke_y ) 
          
            if (is_converged) exit
          end do ! itr lin
!          write(*,*) lmesh%tileID, ke_x, ke_y, itr_lin, gmres_hevi%m_out
          call PROF_rapend( 'hevi_cal_vi_itr_lin', 3)

          do ke_z=1, lmesh%NeZ
            PROG_VARS(:,:,ke_z) = PROG_VARS(:,:,ke_z) + PROG_DEL(:,:,ke_z)
            PROG_VARS0(:,:,ke_z) = PROG_VARS(:,:,ke_z)
          end do  
        end do ! itr nlin

        call PROF_rapend( 'hevi_cal_vi_itr', 3)
      end if

      call PROF_rapstart( 'hevi_cal_vi_retrun_var', 3)
      call vi_eval_Ax( tend(:,:,:),                      & ! (out)
        PROG_VARS, PROG_VARS, DENS_hyd_z, PRES_hyd_z,    & ! (in)
        Dz, Lift, impl_fac, lmesh, elem,                 & ! (in)
        nz, vmapM, vmapP, ke_x, ke_y, .true. )
      
      !tend(:,:,:) = 0.0_RP
      do ke_z=1, lmesh%NeZ
        ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
        DENS_dt(:,ke) = - tend(:,DDENS_VID,ke_z)
        MOMX_dt(:,ke) = - tend(:,MOMX_VID,ke_z)
        MOMY_dt(:,ke) = - tend(:,MOMY_VID,ke_z)
        MOMZ_dt(:,ke) = - tend(:,MOMZ_VID,ke_z)
        RHOT_dt(:,ke) = - tend(:,DRHOT_VID,ke_z)
      end do
      !tend(:,:,:) = 0.0_RP
      call PROF_rapend( 'hevi_cal_vi_retrun_var', 3)
    end do
    end do

    call gmres_hevi%Final()

    call PROF_rapend( 'hevi_cal_vi', 3)    

    return
  end subroutine atm_dyn_nonhydro3d_hevi_cal_vi

  !------------------------------------------------

  subroutine vi_GMRES_core( gmres_hevi, x, wj, is_converged,  & ! (inout)
    x0, b, N, m, PmatDlu, PmatDlu_ipiv, PmatL, PmatU, pinv_v,   & ! (in)
    DENS_hyd, PRES_hyd,                       & ! (in)
    Dz, Lift, impl_fac, lmesh, elem,          & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y  )
    
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(in) :: N
    integer, intent(in) :: m    

    class(GMRES), intent(inout) :: gmres_hevi
    real(RP), intent(inout) :: x(N)
    real(RP), intent(inout) :: wj(N)
    logical, intent(out) :: is_converged
    real(RP), intent(in) :: x0(N)
    real(RP), intent(in) :: b(N)
    real(RP), intent(in) :: PmatDlu(elem%Np*PROG_VARS_NUM,elem%Np*PROG_VARS_NUM,lmesh%NeZ)
    integer, intent(in) :: PmatDlu_ipiv(elem%Np*PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in) :: PmatL(elem%Np,elem%Np,PROG_VARS_NUM,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in) :: PmatU(elem%Np,elem%Np,PROG_VARS_NUM,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(inout) :: pinv_v(N)
    !---
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y    

    integer :: j

    !--------------------------------------

!    call PROF_rapstart('vi_lin_core_eval_Ax', 3)
    call vi_eval_Ax_lin( wj(:),             & ! (out)
      x, x0, DENS_hyd, PRES_hyd,            & ! (in)
      Dz, Lift, impl_fac, lmesh, elem,      & ! (in)
      nz, vmapM, vmapP, ke_x, ke_y, .false. ) ! (in)
!    call PROF_rapend('vi_lin_core_eval_Ax', 3)

 !   call PROF_rapstart('vi_lin_core_pre', 3)
    call gmres_hevi%Iterate_pre( b, wj, is_converged)
    if (is_converged) return
 !   call PROF_rapend('vi_lin_core_pre', 3)

    do j=1, min(m, N)
      call PROF_rapstart('vi_lin_core_pinv_v', 3)
      call matmul_pinv_v( pinv_v, PmatDlu, PmatDlu_ipiv, PmatL, PmatU, gmres_hevi%v(:,j) )
      call PROF_rapend('vi_lin_core_pinv_v', 3)

      call PROF_rapstart('vi_lin_core_eval_Ax', 3)
      call vi_eval_Ax_lin( wj(:),                 & ! (out)
       pinv_v, x0, DENS_hyd, PRES_hyd, & ! (in)
       Dz, Lift, impl_fac, lmesh, elem,           & ! (in)
       nz, vmapM, vmapP, ke_x, ke_y, .false.      ) ! (in)
      call PROF_rapend('vi_lin_core_eval_Ax', 3)
      
      call PROF_rapstart('vi_lin_core_itr', 3)
      call gmres_hevi%Iterate_step_j( j, wj, is_converged)
      call PROF_rapend('vi_lin_core_itr', 3)
      if (is_converged) exit
    end do

    call PROF_rapstart('vi_lin_core_post', 3)
    do j=1, N
      wj(j) = 0.0_RP
    end do
    call gmres_hevi%Iterate_post( wj )
    call PROF_rapstart('vi_lin_core_pinv_v_px0', 3)
    call matmul_pinv_v_plus_x0( x, PmatDlu, PmatDlu_ipiv, PmatL, pmatU, wj)
    call PROF_rapend('vi_lin_core_pinv_v_px0', 3)
    call PROF_rapend('vi_lin_core_post', 3)
    
    return
  contains
    subroutine matmul_pinv_v( pinv_v_, pDlu_, PmatDlu_ipiv_, pL, pU, v)
      implicit none
      real(RP), intent(out) :: pinv_v_(elem%Np,PROG_VARS_NUM,lmesh%NeZ)      
      real(RP), intent(in) :: pDlu_(elem%Np*PROG_VARS_NUM,elem%Np*PROG_VARS_NUM,lmesh%NeZ)
      integer, intent(in) :: PmatDlu_ipiv_(elem%Np*PROG_VARS_NUM,lmesh%NeZ)
      real(RP), intent(in) :: pL(elem%Np,elem%Np,PROG_VARS_NUM,PROG_VARS_NUM,lmesh%NeZ)
      real(RP), intent(in) :: pU(elem%Np,elem%Np,PROG_VARS_NUM,PROG_VARS_NUM,lmesh%NeZ)
      real(RP), intent(in) :: v(elem%Np*PROG_VARS_NUM,lmesh%NeZ)      

      integer :: k, n
      real(RP) :: tmp(elem%Np*PROG_VARS_NUM)
      integer :: vid, vs, ve
      integer :: info
      !------------------------------------
      
      n = elem%Np * PROG_VARS_NUM

      do vid = 1, PROG_VARS_NUM
        vs = 1 + (vid-1) * elem%Np
        ve = vs + elem%Np - 1
        pinv_v_(:,vid,1) = v(vs:ve,1)
      end do
      call DGETRS('N', n, 1, pDlu_(:,:,1), n, PmatDlu_ipiv_(:,1), pinv_v_(:,:,1), n, info)

!       tmp(:) = matmul( pDlu_(:,:,1), v(:,1) )
!       do vid = 1, PROG_VARS_NUM
! !        pinv_v_(:,vid,1) = matmul( pDlu_(:,:,vid,1), v(:,1) )
!         vs = 1 + (vid-1) * elem%Np
!         ve = vs + elem%Np - 1
!         pinv_v_(:,vid,1) = tmp(vs:ve)
!       end do

      do k=2, lmesh%NeZ
        vs = 1; ve = elem%Np
        pinv_v_(:,DDENS_VID,k) = v(vs:ve,k) &
          - matmul( pL(:,:,DDENS_VID,DDENS_VID,k), pinv_v_(:,DDENS_VID,k-1) ) &
          - matmul( pL(:,:,MOMZ_VID ,DDENS_VID,k), pinv_v_(:,MOMZ_VID ,k-1) )

        vs = ve+1; ve = vs + elem%Np - 1
        pinv_v_(:,MOMX_VID,k) = v(vs:ve,k) &
          - matmul( pL(:,:,MOMX_VID,MOMX_VID,k), pinv_v_(:,MOMX_VID,k-1) )

        vs = ve+1; ve = vs + elem%Np - 1
        pinv_v_(:,MOMY_VID,k) = v(vs:ve,k) &
          - matmul( pL(:,:,MOMY_VID,MOMY_VID,k), pinv_v_(:,MOMY_VID,k-1) )

        vs = ve+1; ve = vs + elem%Np - 1          
        pinv_v_(:,MOMZ_VID,k) = v(vs:ve,k) &
          - matmul( pL(:,:,DDENS_VID,MOMZ_VID,k), pinv_v_(:,DDENS_VID,k-1) ) &
          - matmul( pL(:,:,MOMZ_VID ,MOMZ_VID,k), pinv_v_(:,MOMZ_VID ,k-1) ) &
          - matmul( pL(:,:,DRHOT_VID,MOMZ_VID,k), pinv_v_(:,DRHOT_VID,k-1) )

        vs = ve+1; ve = vs + elem%Np - 1          
        pinv_v_(:,DRHOT_VID,k) = v(vs:ve,k) &
          - matmul( pL(:,:,DDENS_VID ,DRHOT_VID,k), pinv_v_(:,DDENS_VID,k-1) ) &
          - matmul( pL(:,:,MOMZ_VID  ,DRHOT_VID,k), pinv_v_(:,MOMZ_VID ,k-1) ) &
          - matmul( pL(:,:,DRHOT_VID ,DRHOT_VID,k), pinv_v_(:,DRHOT_VID,k-1) )
        
        call DGETRS('N', n, 1, pDlu_(:,:,k), n, PmatDlu_ipiv_(:,k), pinv_v_(:,:,k), n, info)

        ! !$omp parallel do
        ! do vid = 1, PROG_VARS_NUM
        !   pinv_v_(:,vid,k) = matmul( pDlu_(:,:,vid,k), tmp(:) )
        ! end do

        ! do vid = 1, PROG_VARS_NUM
        !   vs = 1 + (vid-1) * elem%Np; ve = vs + elem%Np - 1
        !   tmp(vs:ve) = pinv_v_(:,vid,k)
        ! end do            
        ! tmp(:) = matmul( pDlu_(:,:,k), tmp(:) )
        ! do vid = 1, PROG_VARS_NUM
        !   vs = 1 + (vid-1) * elem%Np; ve = vs + elem%Np - 1
        !   pinv_v_(:,vid,k) = tmp(vs:ve)
        ! end do
      end do

      !
      do k=lmesh%NeZ-1, 1, -1
        vs = 1; ve = elem%Np
        tmp(vs:ve) = &
            matmul( pU(:,:,DDENS_VID,DDENS_VID,k), pinv_v_(:,DDENS_VID,k+1) ) &
          + matmul( pU(:,:,MOMZ_VID ,DDENS_VID,k), pinv_v_(:,MOMZ_VID ,k+1) )

        vs = ve+1; ve = vs + elem%Np - 1
        tmp(vs:ve) =  &
            matmul( pU(:,:,MOMX_VID,MOMX_VID,k), pinv_v_(:,MOMX_VID,k+1) )

        vs = ve+1; ve = vs + elem%Np - 1
        tmp(vs:ve) = &
            matmul( pU(:,:,MOMY_VID,MOMY_VID,k), pinv_v_(:,MOMY_VID,k+1) )

        vs = ve+1; ve = vs + elem%Np - 1          
        tmp(vs:ve) = &
            matmul( pU(:,:,DDENS_VID,MOMZ_VID,k), pinv_v_(:,DDENS_VID,k+1) ) &
          + matmul( pU(:,:,MOMZ_VID ,MOMZ_VID,k), pinv_v_(:,MOMZ_VID ,k+1) ) &
          + matmul( pU(:,:,DRHOT_VID,MOMZ_VID,k), pinv_v_(:,DRHOT_VID,k+1) )

        vs = ve+1; ve = vs + elem%Np - 1          
        tmp(vs:ve) = &
            matmul( pU(:,:,DDENS_VID ,DRHOT_VID,k), pinv_v_(:,DDENS_VID,k+1) ) &
          + matmul( pU(:,:,MOMZ_VID  ,DRHOT_VID,k), pinv_v_(:,MOMZ_VID ,k+1) ) &
          + matmul( pU(:,:,DRHOT_VID ,DRHOT_VID,k), pinv_v_(:,DRHOT_VID,k+1) )
  
        call DGETRS('N', n, 1, pDlu_(:,:,k), n, PmatDlu_ipiv_(:,k), tmp(:), n, info)

        !! $omp parallel do
        do vid = 1, PROG_VARS_NUM
          vs = 1 + (vid-1) * elem%Np
          ve = vs + elem%Np - 1  
          pinv_v_(:,vid,k) = pinv_v_(:,vid,k) - tmp(vs:ve)
        end do
        
        ! tmp(:) = matmul( pDlu_(:,:,k), tmp(:) )
        ! do vid = 1, PROG_VARS_NUM
        !   vs = 1 + (vid-1) * elem%Np; ve = vs + elem%Np - 1
        !   pinv_v_(:,vid,k) = pinv_v_(:,vid,k) - tmp(vs:ve)
        ! end do
  
      end do

      return
    end subroutine matmul_pinv_v

    subroutine matmul_pinv_v_plus_x0( x_, pDlu_, PmatDlu_ipiv_, pL, pU, v)
      implicit none
      real(RP), intent(inout) :: x_(elem%Np*PROG_VARS_NUM,lmesh%NeZ)      
      real(RP), intent(in) :: pDlu_(elem%Np*PROG_VARS_NUM,elem%Np*PROG_VARS_NUM,lmesh%NeZ)
      integer, intent(in) :: PmatDlu_ipiv_(elem%Np*PROG_VARS_NUM,lmesh%NeZ)
      real(RP), intent(in) :: pL(elem%Np*PROG_VARS_NUM,elem%Np*PROG_VARS_NUM,lmesh%NeZ)
      real(RP), intent(in) :: pU(elem%Np*PROG_VARS_NUM,elem%Np*PROG_VARS_NUM,lmesh%NeZ)
      real(RP), intent(in) :: v(elem%Np*PROG_VARS_NUM,lmesh%NeZ)      

      integer :: k
      real(RP) :: tmp(elem%Np*PROG_VARS_NUM,lmesh%NeZ)

      !------------------------------------
      
      call matmul_pinv_v( tmp, pDlu_, PmatDlu_ipiv_, pL, pU, v)
      !$omp parallel do
      do k=1, lmesh%NeZ
        x_(:,k) = x_(:,k) + tmp(:,k)
      end do

      return
    end subroutine matmul_pinv_v_plus_x0    
  end subroutine vi_GMRES_core

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
!-
  subroutine vi_eval_Ax_lin( Ax, &
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
    real(RP) :: RHOT_hyd(elem%Np)
    real(RP) :: POT0(elem%Np), DENS0(elem%Np)
    real(RP) :: DPRES(elem%Np)
    integer :: ke_z
    integer :: ke
    real(RP) :: gamm, rgamm

    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    call vi_cal_del_flux_dyn_lin( del_flux,                & ! (out)
      PROG_VARS(:,DDENS_VID,:), PROG_VARS(:,MOMX_VID,:),   & ! (in)
      PROG_VARS(:,MOMY_VID ,:), PROG_VARS(:,MOMZ_VID,:),   & ! (in)
      PROG_VARS(:,DRHOT_VID,:),                            & ! (in)
      PROG_VARS0(:,DDENS_VID,:), PROG_VARS0(:,MOMX_VID,:), & ! (in)
      PROG_VARS0(:,MOMY_VID ,:), PROG_VARS0(:,MOMZ_VID,:), & ! (in)
      PROG_VARS0(:,DRHOT_VID,:),                           & ! (in)
      DENS_hyd, PRES_hyd, nz, vmapM, vmapP,                & ! (in)
      lmesh, elem )                                          ! (in)


    !$omp parallel do private(ke, RHOT_hyd, DPRES, DENS0, POT0, Fz, LiftDelFlx)
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(:,ke_z)/PRES00)**rgamm

      DPRES(:) = gamm * PRES_hyd(:,ke_z) / RHOT_hyd(:)                               &
                 * ( 1.0_RP + PROG_VARS0(:,DRHOT_VID,ke_z) / RHOT_hyd(:) )**(gamm-1) &
                 * PROG_VARS(:,DRHOT_VID,ke_z)

      DENS0(:) = DENS_hyd(:,ke_z) + PROG_VARS0(:,DDENS_VID,ke_z)
      POT0(:) = ( RHOT_hyd(:) + PROG_VARS0(:,DRHOT_VID,ke_z) ) / DENS0(:)

      !- DENS
      call sparsemat_matmul(Dz, PROG_VARS(:,MOMZ_VID,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,DDENS_VID), LiftDelFlx)
      Ax(:,DDENS_VID,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      !- MOMX
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,MOMX_VID), LiftDelFlx)
      Ax(:,MOMX_VID,ke_z) = LiftDelFlx(:)

      !-MOMY
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,MOMY_VID), LiftDelFlx)
      Ax(:,MOMY_VID,ke_z) = LiftDelFlx(:)

      !-MOMZ
      call sparsemat_matmul(Dz, DPRES(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,MOMZ_VID), LiftDelFlx)
      Ax(:,MOMZ_VID,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)             &
                          + Grav * matmul(IntrpMat_VPOrdM1, PROG_VARS(:,DDENS_VID,ke_z))

      !-RHOT
      call sparsemat_matmul(Dz,   POT0(:) * PROG_VARS(:,MOMZ_VID,ke_z)                                            &
                                 + PROG_VARS0(:,MOMZ_VID,ke_z) / DENS0(:) * PROG_VARS(:,DRHOT_VID,ke_z)            &
                                 - POT0(:) * PROG_VARS0(:,MOMZ_VID,ke_z) / DENS0(:) * PROG_VARS(:,DDENS_VID,ke_z), & 
                           Fz )

      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,DRHOT_VID), LiftDelFlx)
      Ax(:,DRHOT_VID,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      !--
      if ( .not. cal_tend_flag ) then
        Ax(:,:,ke_z) =  PROG_VARS(:,:,ke_z) + impl_fac * Ax(:,:,ke_z)
      end if 

    end do    


    return
  end subroutine vi_eval_Ax_lin

  subroutine vi_cal_del_flux_dyn_lin( del_flux, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,         &
    DDENS0_, MOMX0_, MOMY0_, MOMZ0_, DRHOT0_,    &
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
    real(RP) :: rhot_hyd_M, rhot_hyd_P
    real(RP) :: dpresM, dpresP, densM, desnP, MOMZ_P
    real(RP) :: pres0M, pres0P, dens0M, dens0P, pott0M, pott0P, MOMZ0_P

    real(RP) :: gamm, rgamm
    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry
    
    !$omp parallel do private( p, i, iM, iP, &
    !$omp rhot_hyd_M, rhot_hyd_P, dpresM, dpresP, MOMZ_P,          &
    !$omp dens0M, dens0P, pott0M, pott0P, pres0M, pres0P, MOMZ0_P, &
    !$omp swV, alpha0 )
    do ke_z=1, lmesh%NeZ
    do p=1, elem%NfpTot
      i = p + (ke_z-1)*elem%NfpTot
      iM = vmapM(i); iP = vmapP(i)
      
      rhot_hyd_M = PRES00/Rdry * (PRES_hyd(iM)/PRES00)**rgamm
      rhot_hyd_P = PRES00/Rdry * (PRES_hyd(iP)/PRES00)**rgamm

      dpresM = gamm * PRES_hyd(iM) / rhot_hyd_M * ( 1.0_RP + DRHOT0_(iM) / rhot_hyd_M )**(gamm-1) * DRHOT_(iM)
      dpresP = gamm * PRES_hyd(iP) / rhot_hyd_P * ( 1.0_RP + DRHOT0_(iP) / rhot_hyd_P )**(gamm-1) * DRHOT_(iP)

      !-
      dens0M = DENS_hyd(iM) + DDENS0_(iM)
      dens0P = DENS_hyd(iP) + DDENS0_(iP)

      pott0M = ( rhot_hyd_M + DRHOT0_(iM) ) / dens0M
      pott0P = ( rhot_hyd_P + DRHOT0_(iP) ) / dens0P

      pres0M = PRES_hyd(iM) * ( 1.0_RP + DRHOT0_(iM) / rhot_hyd_M )**gamm
      pres0P = PRES_hyd(iP) * ( 1.0_RP + DRHOT0_(iP) / rhot_hyd_P )**gamm

      swV = nz(i)**2
      alpha0 = swV * max( abs(MOMZ0_(iM)/dens0M) + sqrt(gamm*pres0M/dens0M), &
                          abs(MOMZ0_(iP)/dens0P) + sqrt(gamm*pres0P/dens0P)  )

      if (iM==iP .and. (ke_z == 1 .or. ke_z == lmesh%NeZ)) then
        MOMZ_P = - MOMZ_(iM)
        MOMZ0_P = - MOMZ0_(iM)
      else
        MOMZ_P = MOMZ_(iP)
        MOMZ0_P = MOMZ0_(iP)
      end if

      del_flux(i,DDENS_VID) = 0.5_RP * (                    &
                    + ( MOMZ_P - MOMZ_(iM) ) * nz(i)        &
                    - alpha0 * ( DDENS_(iP) - DDENS_(iM) )  )
      
      del_flux(i,MOMX_VID) = 0.5_RP * (                     &
                    - alpha0 * ( MOMX_(iP) - MOMX_(iM) )    )
      
      del_flux(i,MOMY_VID) = 0.5_RP * (                     &  
                    - alpha0 * ( MOMY_(iP) - MOMY_(iM) )    )               
      
      del_flux(i,MOMZ_VID) = 0.5_RP * (                     &
                    + ( dpresP - dpresM ) * nz(i)           &                    
                    - alpha0 * ( MOMZ_P - MOMZ_(iM) )       )
      
      del_flux(i,DRHOT_VID) = 0.5_RP * (                             &
                       (  pott0P * MOMZ_P                            &
                        - pott0M * MOMZ_(iM)                         &
                        + MOMZ0_P    / dens0P * DRHOT_(iP)           &
                        - MOMZ0_(iM) / dens0M * DRHOT_(iM)           &
                        - pott0P * MOMZ0_P    / dens0P * DDENS_(iP)  &
                        + pott0M * MOMZ0_(iM) / dens0M * DDENS_(iM)  &
                       ) * nz(i)                                     &
                       - alpha0 * ( DRHOT_(iP) - DRHOT_(iM) )        )

    end do
    end do

    return
  end subroutine vi_cal_del_flux_dyn_lin

!-
  subroutine vi_construct_pmatInv( PmatDlu, PmatDlu_ipiv, PmatL, PmatU, &
    PROG_VARS0, DENS_hyd, PRES_hyd,                   & ! (in)
    Dz, Lift, impl_fac, lmesh, elem,                  & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y )

    use scale_linalgebra, only: linalgebra_LU
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: PmatDlu(elem%Np*PROG_VARS_NUM,elem%Np*PROG_VARS_NUM,lmesh%NeZ)
    integer, intent(out) :: PmatDlu_ipiv(elem%Np*PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(out) :: PmatL(elem%Np,elem%Np,PROG_VARS_NUM,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(out) :: PmatU(elem%Np,elem%Np,PROG_VARS_NUM,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y

    real(RP) :: RHOT_hyd(elem%Np)
    real(RP) :: POT0(elem%Np,lmesh%NeZ)
    real(RP) :: DENS0(elem%Np,lmesh%NeZ)
    real(RP) :: PRES0(elem%Np,lmesh%NeZ)
    real(RP) :: W0(elem%Np,lmesh%NeZ)
    real(RP) :: DPDRHOT0(elem%Np,lmesh%NeZ)
    integer :: ke_z, ke_z2
    integer :: ke, p, fp, v
    real(RP) :: gamm, rgamm
    real(RP) :: dz_p(elem%Np)
    real(RP) :: PmatD(elem%Np,PROG_VARS_NUM,elem%Np,PROG_VARS_NUM)

    integer :: f1, f2, fp_s, fp_e
    integer :: FmV(elem%Nfp_v)
    integer :: FmV2(elem%Nfp_v)
    real(RP) :: lift_op(elem%Np,elem%NfpTot)
    real(RP) :: Id(elem%Np,elem%Np)
    real(RP) :: lift_(elem%Np,elem%Np)
    real(RP) :: lift_2(elem%Np,elem%Np)
    real(RP) :: lift_s(elem%Np,elem%Np)
    real(RP) :: lift_s2(elem%Np,elem%Np)
    real(RP) :: lift_p(elem%Np,elem%Np)
    real(RP) :: lift_p2(elem%Np,elem%Np)
    real(RP) :: lift_pt(elem%Np,elem%Np)
    real(RP) :: lift_pt2(elem%Np,elem%Np)
    real(RP) :: lift_pt_rho(elem%Np,elem%Np)
    real(RP) :: lift_pt_rho2(elem%Np,elem%Np)
    real(RP) :: lift_pt_rhot(elem%Np,elem%Np)
    real(RP) :: lift_pt_rhot2(elem%Np,elem%Np)
    real(RP) :: tmp(elem%Nfp_v)
    real(RP) :: fac, facw

    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    lift_op(:,:) = elem%Lift

    Id(:,:) = 0.0_RP
    do p=1, elem%Np
      Id(p,p) = 1.0_RP
    end do

    !$omp parallel do private(RHOT_hyd)    
    do ke_z=1, lmesh%NeZ
      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(:,ke_z)/PRES00)**rgamm

      DPDRHOT0(:,ke_z) = gamm * PRES_hyd(:,ke_z) / RHOT_hyd(:)                      &
                  * ( 1.0_RP + PROG_VARS0(:,DRHOT_VID,ke_z) / RHOT_hyd(:) )**(gamm-1) 

      DENS0(:,ke_z) = DENS_hyd(:,ke_z) + PROG_VARS0(:,DDENS_VID,ke_z)
      POT0(:,ke_z) = ( RHOT_hyd(:) + PROG_VARS0(:,DRHOT_VID,ke_z) ) / DENS0(:,ke_z)

      W0(:,ke_z) = PROG_VARS0(:,MOMZ_VID,ke_z) / DENS0(:,ke_z)
      PRES0(:,ke_z) = PRES_hyd(:,ke_z) * ( 1.0_RP + PROG_VARS0(:,DRHOT_VID,ke_z) / RHOT_hyd(:) )**gamm
    end do

    !$omp parallel do private(ke, p, fp, v, f1, f2, ke_z2, dz_p, &
    !$omp fac, facw, tmp, lift_, lift_2, lift_s, lift_s2, lift_p, lift_p2, lift_pt, lift_pt2, lift_pt_rho, lift_pt_rho2, lift_pt_rhot, lift_pt_rhot2, &
    !$omp PmatD, FmV, FmV2, fp_s, fp_e)
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

      !-----
      PmatD(:,:,:,:) = 0.0_RP
      PmatL(:,:,:,:,ke_z) = 0.0_RP
      PmatU(:,:,:,:,ke_z) = 0.0_RP

      do v=1, PROG_VARS_NUM
        PmatD(:,v,:,v) = Id(:,:)
      end do

      do p=1, elem%Np
        dz_p(:) = lmesh%Escale(p,ke,3,3) * elem%Dx3(p,:) 

        ! DDENS
        PmatD(p,DDENS_VID,:,MOMZ_VID) = impl_fac * dz_p(:) 

        ! MOMZ
        PmatD(p,MOMZ_VID,:,DDENS_VID) = impl_fac * Grav * IntrpMat_VPOrdM1(p,:)
        PmatD(p,MOMZ_VID,:,DRHOT_VID) = impl_fac * dz_p(:) * DPDRHOT0(:,ke_z)

        !DRHOT
        PmatD(p,DRHOT_VID,:,DDENS_VID) = - impl_fac * POT0(:,ke_z) * W0(:,ke_z) * dz_p(:)
        PmatD(p,DRHOT_VID,:,MOMZ_VID) = impl_fac * POT0(:,ke_z) * dz_p(:) 
        PmatD(p,DRHOT_VID,:,DRHOT_VID) = PmatD(p,DRHOT_VID,:,DRHOT_VID) &
                                      + impl_fac * W0(:,ke_z) * dz_p(:) 
      end do

      do f1=1, 2
        if (f1==1) then
          f2 = 2; ke_z2 = max(ke_z-1, 1)
        else
          f2 = 1; ke_z2 = min(ke_z+1, lmesh%NeZ)
        end if
        fac  = 0.5_RP * impl_fac
        facw = fac
        if ( (ke_z == 1 .and. f1==1) .or. (ke_z == lmesh%NeZ .and. f1==elem%Nfaces_v) ) then
          facw = impl_fac
          f2 = f1
        end if

        FmV (:) = elem%Fmask_v(:,f1)
        FmV2(:) = elem%Fmask_v(:,f2)

        fp_s = elem%Nfp_h * elem%Nfaces_h + 1 + (f1-1)*elem%Nfp_v
        fp_e = fp_s + elem%Nfp_v - 1
        
        !
        lift_ (:,:) = 0.0_RP
        lift_2(:,:) = 0.0_RP
        !
        lift_s (:,:) = 0.0_RP
        lift_s2(:,:) = 0.0_RP
        !
        lift_p (:,:) = 0.0_RP
        lift_p2(:,:) = 0.0_RP
        !
        lift_pt (:,:) = 0.0_RP
        lift_pt2(:,:) = 0.0_RP
        !
        lift_pt_rho (:,:) = 0.0_RP
        lift_pt_rho2(:,:) = 0.0_RP
        !
        lift_pt_rhot (:,:) = 0.0_RP
        lift_pt_rhot2(:,:) = 0.0_RP   

        !--
        do fp=fp_s, fp_e
          p = fp-fp_s+1
          tmp(:) = lift_op(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)
          lift_ (FmV,FmV (p)) = tmp(:)
          lift_2(FmV,FmV2(p)) = tmp(:)
        end do

        !---
        do fp=fp_s, fp_e
          p = fp-fp_s+1
          tmp(:) = lift_op(FmV,fp) * lmesh%Fscale(fp,ke) * &
                         max( abs( W0(FmV (p),ke_z ) ) + sqrt( gamm * PRES0(FmV (p),ke_z ) / DENS0(FmV (p),ke_z ) ), &
                              abs( W0(FmV2(p),ke_z2) ) + sqrt( gamm * PRES0(FmV2(p),ke_z2) / DENS0(FmV2(p),ke_z2) )  )

          lift_s (FmV,FmV (p)) = tmp(:)
          lift_s2(FmV,FmV2(p)) = tmp(:)
        end do

        !--
        do fp=fp_s, fp_e
          p = fp-fp_s+1
          tmp(:) = lift_op(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)
          lift_p (FmV,FmV (p)) = tmp(:) * DPDRHOT0(FmV (p),ke_z )
          lift_p2(FmV,FmV2(p)) = tmp(:) * DPDRHOT0(FmV2(p),ke_z2)
        end do

        !--
        do fp=fp_s, fp_e
          p = fp-fp_s+1
          tmp(:) = lift_op(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)
          lift_pt (FmV,FmV (p)) = tmp(:) * POT0(FmV (p),ke_z )
          lift_pt2(FmV,FmV2(p)) = tmp(:) * POT0(FmV2(p),ke_z2)
        end do

        !--

        do fp=fp_s, fp_e
          p = fp-fp_s+1
          tmp(:) = lift_op(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)          
          lift_pt_rho (FmV,FmV (p)) = - tmp(:) * POT0(FmV (p),ke_z ) * W0(FmV (p),ke_z )
          lift_pt_rho2(FmV,FmV2(p)) = - tmp(:) * POT0(FmV2(p),ke_z2) * W0(FmV2(p),ke_z2)
        end do

        !--

        do fp=fp_s, fp_e
          p = fp-fp_s+1
          tmp(:) = lift_op(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)          
          lift_pt_rhot (FmV,FmV (p)) = tmp(:) * W0(FmV (p),ke_z )
          lift_pt_rhot2(FmV,FmV2(p)) = tmp(:) * W0(FmV2(p),ke_z2)
        end do

        !----        

        if ( (ke_z == 1 .and. f1==1) .or. (ke_z == lmesh%NeZ .and. f1==elem%Nfaces_v) ) then
          PmatD(:,DDENS_VID,:,MOMZ_VID) = PmatD(:,DDENS_VID,:,MOMZ_VID) &
            - facw * lift_(:,:)
                    
          PmatD(:,DRHOT_VID,:,MOMZ_VID) = PmatD(:,DRHOT_VID,:,MOMZ_VID) &
            - facw * lift_pt(:,:)

          PmatD(:,MOMZ_VID,:,MOMZ_VID) = PmatD(:,MOMZ_VID,:,MOMZ_VID) &
            + facw * lift_s(:,:)
        else
          
          do v=1, PROG_VARS_NUM
            PmatD(:,v,:,v) = PmatD(:,v,:,v) + fac * lift_s(:,:)
            if (f1 == 1) then
              PmatL(:,:,v,v,ke_z) = - fac * lift_s2(:,:)
            else
              PmatU(:,:,v,v,ke_z) = - fac * lift_s2(:,:)
            end if
          end do

          !-
          PmatD(:,DDENS_VID,:,MOMZ_VID) = PmatD(:,DDENS_VID,:,MOMZ_VID) &
            - fac * lift_(:,:)
          if (f1 == 1) then
            PmatL(:,:,MOMZ_VID,DDENS_VID,ke_z) = fac * lift_2(:,:)
          else
            PmatU(:,:,MOMZ_VID,DDENS_VID,ke_z) = fac * lift_2(:,:)
          end if
            
          !-
          PmatD(:,MOMZ_VID,:,DRHOT_VID) = PmatD(:,MOMZ_VID,:,DRHOT_VID) &
            - fac * lift_p(:,:)
          if (f1 == 1) then
            PmatL(:,:,DRHOT_VID,MOMZ_VID,ke_z) = fac * lift_p2(:,:)
          else
            PmatU(:,:,DRHOT_VID,MOMZ_VID,ke_z) = fac * lift_p2(:,:)
          end if

          PmatD(:,DRHOT_VID,:,MOMZ_VID) = PmatD(:,DRHOT_VID,:,MOMZ_VID) &
            - fac * lift_pt(:,:)
          PmatD(:,DRHOT_VID,:,DDENS_VID) = PmatD(:,DRHOT_VID,:,DDENS_VID) &
            - fac * lift_pt_rho(:,:)
          PmatD(:,DRHOT_VID,:,DRHOT_VID) = PmatD(:,DRHOT_VID,:,DRHOT_VID) &
            - fac * lift_pt_rhot(:,:)
          if (f1 == 1) then
            PmatL(:,:,MOMZ_VID ,DRHOT_VID,ke_z) = fac * lift_pt2(:,:)
            PmatL(:,:,DDENS_VID,DRHOT_VID,ke_z) = fac * lift_pt_rho2(:,:)
            PmatL(:,:,DRHOT_VID,DRHOT_VID,ke_z) = PmatL(:,:,DRHOT_VID,DRHOT_VID,ke_z) &
                                                + fac * lift_pt_rhot2(:,:)
          else
            PmatU(:,:,MOMZ_VID ,DRHOT_VID,ke_z) = fac * lift_pt2(:,:)
            PmatU(:,:,DDENS_VID,DRHOT_VID,ke_z) = fac * lift_pt_rho2(:,:)
            PmatU(:,:,DRHOT_VID,DRHOT_VID,ke_z) = PmatU(:,:,DRHOT_VID,DRHOT_VID,ke_z) &
                                                + fac * lift_pt_rhot2(:,:)
          end if
  
        end if
      end do

      call PROF_rapstart( 'hevi_cal_vi_cal_pinv', 3)
      call get_PmatD_LU( PmatDlu(:,:,ke_z), PmatD, PmatDlu_ipiv(:,ke_z), elem%Np * PROG_VARS_NUM )
      !call get_PmatD_inv( PmatDlu(:,:,ke_z), PmatD, PmatDlu_ipiv(:,ke_z), elem%Np * PROG_VARS_NUM )
      call PROF_rapend( 'hevi_cal_vi_cal_pinv', 3)
    end do    

    return

  contains 
    subroutine get_PmatD_LU( pmatDlu_, pmatD_, pmatDlu_ipiv_, N)
      use scale_linalgebra, only: linalgebra_LU
      implicit none
      integer, intent(in) :: N
      real(RP), intent(out) :: pmatDlu_(N,N)
      real(RP), intent(in) :: pmatD_(N,N)
      integer, intent(out) :: pmatDlu_ipiv_(N)
      integer :: info
      !------------------------------------------

      pmatDlu_(:,:) = pmatD_(:,:)
      call linalgebra_LU(pmatDlu_, pmatDlu_ipiv_)
      return
    end subroutine get_PmatD_LU
    
    subroutine get_PmatD_inv( pmatDinv_, pmatD_, pmatDlu_ipiv_, N)
      use scale_linalgebra, only: linalgebra_inv
      implicit none
      integer, intent(in) :: N
      real(RP), intent(out) :: pmatDinv_(N,N)
      real(RP), intent(in) :: pmatD_(N,N)
      integer, intent(out) :: pmatDlu_ipiv_(N)
      !------------------------------------------

      pmatDinv_(:,:) = linalgebra_inv(pmatD_)
      return
    end subroutine get_PmatD_inv    
  end subroutine vi_construct_pmatInv

    
end module scale_atm_dyn_nonhydro3d_hevi_gmres
