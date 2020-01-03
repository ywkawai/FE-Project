!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_nonhydro3d
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
  public :: atm_dyn_nonhydro3d_Init
  public :: atm_dyn_nonhydro3d_Final
  public :: atm_dyn_nonhydro3d_prepair_expfilter
  public :: atm_dyn_nonhydro3d_cal_tend
  public :: atm_dyn_nonhydro3d_cal_grad_diffVars
  public :: atm_dyn_nonhydro3d_filter_prgvar

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

  real(RP), private, allocatable :: FilterMat(:,:)
  real(RP), private, allocatable :: IntrpMat_VPOrdM1(:,:)

  private :: cal_del_flux_dyn
  private :: cal_del_gradDiffVar

contains
  subroutine atm_dyn_nonhydro3d_Init( mesh )

    implicit none
    class(MeshBase3D), intent(in) :: mesh

    integer :: p1, p2, p_
    type(ElementBase3D), pointer :: elem
    real(RP) :: invV_VPOrdM1(mesh%refElem3D%Np,mesh%refElem3D%Np)
    !--------------------------------------------

    elem => mesh%refElem3D
    allocate( IntrpMat_VPOrdM1(elem%Np,elem%Np) )
    
    InvV_VPOrdM1(:,:) = elem%invV
    do p2=1, elem%PolyOrder_h+1
    do p1=1, elem%PolyOrder_h+1
      p_ = p1 + (p2-1)*(elem%PolyOrder_h + 1) + elem%PolyOrder_v*(elem%PolyOrder_h + 1)**2
      InvV_VPOrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem%V, invV_VPOrdM1)

    return
  end subroutine atm_dyn_nonhydro3d_Init

  subroutine atm_dyn_nonhydro3d_prepair_expfilter(  &
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
    do p1=1, elem%Nnode_h1D
      eta = dble(p1-1)/dble(elem%PolyOrder_v)
      if ( eta >  etac .and. p1 /= 1) then
        filter1D_v(p1) = exp( -  alpha*( ((eta - etac)/(1.0_RP - etac))**ord ))
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
  end subroutine atm_dyn_nonhydro3d_prepair_expfilter

  subroutine atm_dyn_nonhydro3d_Final()
    implicit none
    !--------------------------------------------
    
    deallocate( IntrpMat_VPOrdM1 )
    if( allocated(FilterMat) ) deallocate( FilterMat )
    
    return
  end subroutine atm_dyn_nonhydro3d_Final  

  !-------------------------------

  subroutine atm_dyn_nonhydro3d_filter_prgvar( &
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

    do k=1, lmesh%Ne
      DDENS_(:,k) = matmul(FilterMat,DDENS_(:,k))
      MOMX_(:,k) = matmul(FilterMat,MOMX_(:,k))
      MOMY_(:,k) = matmul(FilterMat,MOMY_(:,k))
      MOMZ_(:,k) = matmul(FilterMat,MOMZ_(:,k))
      DRHOT_(:,k) = matmul(FilterMat,DRHOT_(:,k))
    end do
    
    return
  end subroutine atm_dyn_nonhydro3d_filter_prgvar

  subroutine atm_dyn_nonhydro3d_cal_tend( &
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
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%Ne)
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
    real(RP) :: dens_(elem%Np), RHOT_(elem%Np), dpres_(elem%Np)
    real(RP) :: pres_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np)
    real(RP) :: Cori(elem%Np)

    integer :: ke, ke2d
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
      call sparsemat_matmul(Dx, u_(:)*MOMX_(:,ke) + dpres_(:) - viscCoef_h*dens_(:)*GxU_(:,ke), Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMX_(:,ke)             - viscCoef_h*dens_(:)*GyU_(:,ke), Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMX_(:,ke)             - viscCoef_v*dens_(:)*GzU_(:,ke) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_MOMX_ID), LiftDelFlx)

      MOMX_dt(:,ke) = - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:)      &
          + lmesh%Escale(:,ke,2,2) * Fy(:)      &
          + lmesh%Escale(:,ke,3,3) * Fz(:)      &
          + LiftDelFlx(:)                       &
          - Cori(:)*MOMY_(:,ke)                 &
          )

      !-- MOMY
      call sparsemat_matmul(Dx, u_(:)*MOMY_(:,ke)             - viscCoef_h*dens_(:)*GxV_(:,ke), Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMY_(:,ke) + dpres_(:) - viscCoef_h*dens_(:)*GyV_(:,ke), Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMY_(:,ke)             - viscCoef_v*dens_(:)*GzV_(:,ke) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_MOMY_ID), LiftDelFlx)

      MOMY_dt(:,ke) = - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDelFlx(:)                  &
          + Cori(:)*MOMX_(:,ke)            &
          )
    
      !-- MOMZ
      call sparsemat_matmul(Dx, u_(:)*MOMZ_(:,ke)             - viscCoef_h*dens_(:)*GxW_(:,ke), Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMZ_(:,ke)             - viscCoef_h*dens_(:)*GyW_(:,ke), Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMZ_(:,ke) + dpres_(:) - viscCoef_v*dens_(:)*GzW_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_MOMZ_ID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)   &
          + lmesh%Escale(:,ke,2,2) * Fy(:)   &
          + lmesh%Escale(:,ke,3,3) * Fz(:)   &
          + LiftDelFlx(:)                )   &
          - matmul(IntrpMat_VPOrdM1, DDENS_(:,ke)) * Grav
        !- DDENS_(:,k)*Grav
        

      !-- RHOT
      call sparsemat_matmul(Dx, u_(:)*RHOT_(:) - diffCoef_h*dens_(:)*GxPT_(:,ke), Fx)
      call sparsemat_matmul(Dy, v_(:)*RHOT_(:) - diffCoef_h*dens_(:)*GyPT_(:,ke), Fy)
      call sparsemat_matmul(Dz, w_(:)*RHOT_(:) - diffCoef_v*dens_(:)*GzPT_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARS_DRHOT_ID), LiftDelFlx)
      
      RHOT_dt(:,ke) =  - (  &
            lmesh%Escale(:,ke,1,1) * Fx(:) &
          + lmesh%Escale(:,ke,2,2) * Fy(:) &
          + lmesh%Escale(:,ke,3,3) * Fz(:) &
          + LiftDelFlx )
    
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_nonhydro3d_cal_tend

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
    real(RP) :: VelP, VelM, alpha
    real(RP) :: uM, uP, vM, vP, wM, wP, presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd
    real(RP) :: dDiffFluxU, dDiffFluxV, dDiffFluxW, dDiffFluxPT
    real(RP) :: gamm, rgamm
    real(RP) :: mu, viscCoef, diffCoef
    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      rhot_hyd = PRES00/Rdry * (PRES_hyd(iM)/PRES00)**rgamm
      
      rhotM = rhot_hyd + DRHOT_(iM)
      presM = PRES00 * (Rdry*rhotM/PRES00)**gamm
      dpresM = presM - PRES_hyd(iM)

      rhotP = rhot_hyd + DRHOT_(iP) 
      presP = PRES00 * (Rdry*rhotP/PRES00)**gamm
      dpresP = presP - PRES_hyd(iM)

      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iM)

      VelM = (MOMX_(iM)*nx(i) + MOMY_(iM)*ny(i) + MOMZ_(iM)*nz(i))/densM
      VelP = (MOMX_(iP)*nx(i) + MOMY_(iP)*ny(i) + MOMZ_(iP)*nz(i))/densP

      alpha = max( sqrt(gamm*presM/densM) + abs(VelM), sqrt(gamm*presP/densP) + abs(VelP)  )
      mu = (2.0_RP * dble((elem%PolyOrder_h+1)*(elem%PolyOrder_h+2)) / 2.0_RP / 600.0_RP) 
      viscCoef = viscCoef_h*(abs(nx(i)) + abs(ny(i))) + viscCoef_v*abs(nz(i))
      diffCoef = diffCoef_h*(abs(nx(i)) + abs(ny(i))) + diffCoef_v*abs(nz(i))
      
      if (diffCoef > 0.0_RP) then
        dDiffFluxU = viscCoef*( &
            (densP*GxU_(iP) - densM*GxU_(iM))*nx(i)                    &
          + (densP*GyU_(iP) - densM*GyU_(iM))*ny(i)                    &
          + (densP*GzU_(iP) - densM*GzU_(iM))*nz(i)                    &
          + mu*(densP + densM)*(MOMX_(iP)/densP - MOMX_(iM)/densM) )

        dDiffFluxV = viscCoef*( &
          (densP*GxV_(iP) - densM*GxV_(iM))*nx(i)                    &
        + (densP*GyV_(iP) - densM*GyV_(iM))*ny(i)                    &
        + (densP*GzV_(iP) - densM*GzV_(iM))*nz(i)                    &
        + mu*(densP + densM)*(MOMY_(iP)/densP - MOMY_(iM)/densM) )
        
        dDiffFluxW = viscCoef*( &
            (densP*GxW_(iP) - densM*GxW_(iM))*nx(i)                    &
          + (densP*GyW_(iP) - densM*GyW_(iM))*ny(i)                    &
          + (densP*GzW_(iP) - densM*GzW_(iM))*nz(i)                    &
          + mu*(densP + densM)*(MOMZ_(iP)/densP - MOMZ_(iM)/densM) )

        dDiffFluxPT = diffCoef*( &
            (densP*GxPT_(iP) - densM*GxPT_(iM))*nx(i)     &
          + (densP*GzPT_(iP) - densM*GzPT_(iM))*nz(i)     &
          + mu*(densP + densM)*(rhotP/densP - rhotM/densM) )
      else
        dDiffFluxU  = 0.0_RP
        dDiffFluxV  = 0.0_RP
        dDiffFluxW  = 0.0_RP
        dDiffFluxPT = 0.0_RP
      end if
      
      del_flux(i,VARS_DDENS_ID) = 0.5_RP*(               &
                    ( densP*VelP - densM*VelM )          &
                    - alpha*(DDENS_(iP) - DDENS_(iM))   )
      
      del_flux(i,VARS_MOMX_ID) = 0.5_RP*(                &
                    ( MOMX_(iP)*VelP - MOMX_(iM)*VelM )  &
                    + ( dpresP - dpresM )*nx(i)          &
                    - alpha*(MOMX_(iP) - MOMX_(iM))      &
                    - dDiffFluxU                        )
      
      del_flux(i,VARS_MOMY_ID) = 0.5_RP*(                &
                    ( MOMY_(iP)*VelP - MOMY_(iM)*VelM )  &
                    + ( dpresP - dpresM )*ny(i)          &
                    - alpha*(MOMY_(iP) - MOMY_(iM))      &
                    - dDiffFluxV                        )               
      
      del_flux(i,VARS_MOMZ_ID) = 0.5_RP*(                &
                    ( MOMZ_(iP)*VelP - MOMZ_(iM)*VelM)   &
                    + ( dpresP - dpresM )*nz(i)          &                    
                    - alpha*(MOMZ_(iP) - MOMZ_(iM))      &
                    - dDiffFluxW                        )
      
      del_flux(i,VARS_DRHOT_ID) = 0.5_RP*(               &
                    ( rhotP*VelP - rhotM*VelM )          &
                    - alpha*(DRHOT_(iP) - DRHOT_(iM))    &
                    - dDiffFluxPT                       )

    end do

    return
  end subroutine cal_del_flux_dyn

  subroutine atm_dyn_nonhydro3d_cal_grad_diffVars( &
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
      GzW_(:,ke) = lmesh%Escale(:,ke,2,2)*Fz(:) + LiftDelFlx(:)

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
  end subroutine atm_dyn_nonhydro3d_cal_grad_diffVars

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

end module scale_atm_dyn_nonhydro3d
