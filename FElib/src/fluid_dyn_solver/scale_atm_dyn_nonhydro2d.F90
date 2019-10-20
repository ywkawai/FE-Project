!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_nonhydro2d
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
  use scale_element_base
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
  public :: atm_dyn_nonhydro2d_Init
  public :: atm_dyn_nonhydro2d_Final
  public :: atm_dyn_nonhydro2d_prepair_expfilter
  public :: atm_dyn_nonhydro2d_cal_tend
  public :: atm_dyn_nonhydro2d_cal_grad_diffVars
  public :: atm_dyn_nonhydro2d_filter_prgvar

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
  integer, private, parameter :: VARS_MOMZ_ID   = 3
  integer, private, parameter :: VARS_DRHOT_ID  = 4
  integer, private, parameter :: PROG_VARS_NUM  = 4
  
  integer, private, parameter :: VARS_GxU_ID      = 1
  integer, private, parameter :: VARS_GzU_ID      = 2
  integer, private, parameter :: VARS_GxW_ID      = 3
  integer, private, parameter :: VARS_GzW_ID      = 4
  integer, private, parameter :: VARS_GxPT_ID     = 5
  integer, private, parameter :: VARS_GzPT_ID     = 6
  integer, private, parameter :: AUX_DIFFVARS_NUM = 6

  real(RP), private, allocatable :: FilterMat(:,:)

  private :: cal_del_flux_dyn
  private :: cal_del_gradDiffVar

contains
  subroutine atm_dyn_nonhydro2d_Init( mesh )

    implicit none
    class(MeshBase2D), intent(in) :: mesh
    !--------------------------------------------

    return
  end subroutine atm_dyn_nonhydro2d_Init

  subroutine atm_dyn_nonhydro2d_prepair_expfilter(  &
    elem,                                           &
    etac, alpha, ord )

    implicit none
    class(elementbase2D), intent(in) :: elem
    real(RP), intent(in) :: etac
    real(RP), intent(in) :: alpha
    real(RP), intent(in) :: ord

    real(RP) :: filter1D(elem%Nfp)
    real(RP) :: eta
    integer :: p1, p2
    integer :: l
    !----------------------------------------------------

    filter1D(:) = 1.0_RP
    do p1=1, elem%Nfp
      eta = dble(p1-1)/dble(elem%PolyOrder)
      if ( eta >  etac .and. p1 /= 1) then
        filter1D(p1) = exp( -  alpha*( ((eta - etac)/(1.0_RP - etac))**ord ))
      end if
    end do

    allocate( FilterMat(elem%Np,elem%Np) )
    FilterMat(:,:) = 0.0_RP
    do p2=1, elem%Nfp
    do p1=1, elem%Nfp
      l = p1 + (p2-1)*elem%Nfp
      FilterMat(l,l) = filter1D(p1) * filter1D(p2)
    end do  
    end do
    FilterMat(:,:) = matmul(FilterMat, elem%invV)
    FilterMat(:,:) = matmul(elem%V, FilterMat)
    
    return
  end subroutine atm_dyn_nonhydro2d_prepair_expfilter

  subroutine atm_dyn_nonhydro2d_Final()
    implicit none
    !--------------------------------------------
    
    if( allocated(FilterMat) ) deallocate( FilterMat )
    
    return
  end subroutine atm_dyn_nonhydro2d_Final  

  !-------------------------------

  subroutine atm_dyn_nonhydro2d_filter_prgvar( &
    DDENS_, MOMX_, MOMZ_, DRHOT_, lmesh, elem  )
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem    
    real(RP), intent(inout)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: DRHOT_(elem%Np,lmesh%NeA)
    
    integer :: k
    !------------------------------------

    do k=1, lmesh%Ne
      DDENS_(:,k) = matmul(FilterMat,DDENS_(:,k))
      MOMX_(:,k) = matmul(FilterMat,MOMX_(:,k))
      MOMZ_(:,k) = matmul(FilterMat,MOMZ_(:,k))
      DRHOT_(:,k) = matmul(FilterMat,DRHOT_(:,k))
    end do
    
    return
  end subroutine atm_dyn_nonhydro2d_filter_prgvar

  subroutine atm_dyn_nonhydro2d_cal_tend( &
    DENS_dt, MOMX_dt, MOMZ_dt, RHOT_dt,               & ! (out)
    DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, & ! (in)
    GxU_, GzU_, GxW_, GzW_, GxPT_, GzPT_,             & ! (in)
    viscCoef_h, viscCoef_v, diffCoef_h, diffCoef_v,   & ! (in)
    Dx, Dz, Sx, Sz, Lift, lmesh, elem)

    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dz, Sx, Sz, Lift
    real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GxU_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzU_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GxW_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzW_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GxPT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzPT_(elem%Np,lmesh%NeA) 
    real(RP), intent(in) :: viscCoef_h
    real(RP), intent(in) :: viscCoef_v
    real(RP), intent(in) :: diffCoef_h
    real(RP), intent(in) :: diffCoef_v

    real(RP) :: Fx(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: dens_(elem%Np), RHOT_(elem%Np), dpres_(elem%Np)
    real(RP) :: pres_(elem%Np), u_(elem%Np), w_(elem%Np)

    integer :: k

    integer :: p1, p_
    real(RP) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)
    real(RP) :: invV_VPOrdM1(elem%Np,elem%Np)
    !------------------------------------------------------------------------

    InvV_VPOrdM1(:,:) = elem%invV
    do p1=1, elem%PolyOrder+1
      p_ = p1 + elem%PolyOrder*(elem%PolyOrder + 1)
      InvV_VPOrdM1(p_,:) = 0.0_RP
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem%V, invV_VPOrdM1)
    
    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 2)
    call cal_del_flux_dyn( del_flux,                              & ! (out)
      DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,           & ! (in)
      GxU_, GzU_, GxW_, GzW_, GxPT_, GzPT_,                       & ! (in)
      diffCoef_h, diffCoef_v,                                     & ! (in)
      viscCoef_h, viscCoef_v,                                     & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2),             & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                   & ! (in)
      lmesh, elem )                                                 ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 2)
    
    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 2)
    do k = lmesh%NeS, lmesh%NeE
      !--

      RHOT_(:) = PRES00/Rdry * (PRES_hyd(:,k)/PRES00)**(CVdry/CPdry) + DRHOT_(:,k)
      pres_(:) = PRES00 * (Rdry*RHOT_(:)/PRES00)**(CPdry/Cvdry)
      dpres_(:) = pres_(:) - PRES_hyd(:,k)
      dens_(:) = DDENS_(:,k) + DENS_hyd(:,k)

      u_(:) = MOMX_(:,k)/dens_(:)
      w_(:) = MOMZ_(:,k)/dens_(:)

      !-- DENS
      call sparsemat_matmul(Dx, MOMX_(:,k), Fx)
      call sparsemat_matmul(Dz, MOMZ_(:,k), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_DDENS_ID), LiftDelFlx)

      DENS_dt(:,k) = - ( &
            lmesh%Escale(:,k,1,1) * Fx(:) &
          + lmesh%Escale(:,k,2,2) * Fz(:) &
          + LiftDelFlx(:) )
      
      !-- MOMX
      call sparsemat_matmul(Dx, u_(:)*MOMX_(:,k) + dpres_(:) - viscCoef_h*dens_(:)*GxU_(:,k), Fx)
      call sparsemat_matmul(Dz, w_(:)*MOMX_(:,k)             - viscCoef_v*dens_(:)*GzU_(:,k) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_MOMX_ID), LiftDelFlx)

      MOMX_dt(:,k) = - (  &
            lmesh%Escale(:,k,1,1) * Fx(:) &
          + lmesh%Escale(:,k,2,2) * Fz(:) &
          + LiftDelFlx(:) )

      !-- MOMZ
      call sparsemat_matmul(Dx, u_(:)*MOMZ_(:,k)             - viscCoef_h*dens_(:)*GxW_(:,k), Fx)
      call sparsemat_matmul(Dz, w_(:)*MOMZ_(:,k) + dpres_(:) - viscCoef_v*dens_(:)*GzW_(:,k), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_MOMZ_ID), LiftDelFlx)
      
      MOMZ_dt(:,k) = - ( &
            lmesh%Escale(:,k,1,1) * Fx(:)   &
          + lmesh%Escale(:,k,2,2) * Fz(:)   &
          + LiftDelFlx(:)                )  &
          - matmul(IntrpMat_VPOrdM1, DDENS_(:,k)) * Grav
        !- DDENS_(:,k)*Grav
        

      !-- RHOT
      call sparsemat_matmul(Dx, u_(:)*RHOT_(:) - diffCoef_h*dens_(:)*GxPT_(:,k), Fx)
      call sparsemat_matmul(Dz, w_(:)*RHOT_(:) - diffCoef_v*dens_(:)*GzPT_(:,k), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_DRHOT_ID), LiftDelFlx)
      
      RHOT_dt(:,k) =  - (  &
            lmesh%Escale(:,k,1,1) * Fx(:) &
          + lmesh%Escale(:,k,2,2) * Fz(:) &
          + LiftDelFlx )
    
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine atm_dyn_nonhydro2d_cal_tend

  !------

  subroutine cal_del_flux_dyn( del_flux, &
    DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,     &
    GxU_, GzU_, GxW_, GzW_, GxPT_, GzPT_,                 &
    viscCoef_h, viscCoef_v,                               &
    diffCoef_h, diffCoef_v,                               &
    nx, nz, vmapM, vmapP, lmesh, elem )

    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne,PROG_VARS_NUM)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GxU_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzU_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GxW_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzW_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GxPT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzPT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: viscCoef_h
    real(RP), intent(in) :: viscCoef_v
    real(RP), intent(in) :: diffCoef_h
    real(RP), intent(in) :: diffCoef_v
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    real(RP) :: uM, uP, wM, wP, presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd
    real(RP) :: dDiffFluxU, dDiffFluxW, dDiffFluxPT
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

      VelM = (MOMX_(iM)*nx(i) + MOMZ_(iM)*nz(i))/densM
      VelP = (MOMX_(iP)*nx(i) + MOMZ_(iP)*nz(i))/densP

      alpha = max( sqrt(gamm*presM/densM) + abs(VelM), sqrt(gamm*presP/densP) + abs(VelP)  )
      mu = (2.0_RP * dble((elem%PolyOrder+1)*(elem%PolyOrder+2)) / 2.0_RP / 600.0_RP) 
      viscCoef = viscCoef_h*abs(nx(i)) + viscCoef_v*abs(nz(i))
      diffCoef = diffCoef_h*abs(nx(i)) + diffCoef_v*abs(nz(i))
      
      if (diffCoef > 0.0_RP) then
        dDiffFluxU = viscCoef*( &
            (densP*GxU_(iP) - densM*GxU_(iM))*nx(i)                    &
          + (densP*GzU_(iP) - densM*GzU_(iM))*nz(i)                    &
          + mu*(densP + densM)*(MOMX_(iP)/densP - MOMX_(iM)/densM) )
        
        dDiffFluxW = viscCoef*( &
            (densP*GxW_(iP) - densM*GxW_(iM))*nx(i)                    &
          + (densP*GzW_(iP) - densM*GzW_(iM))*nz(i)                    &
          + mu*(densP + densM)*(MOMZ_(iP)/densP - MOMZ_(iM)/densM) )

        dDiffFluxPT = diffCoef*( &
            (densP*GxPT_(iP) - densM*GxPT_(iM))*nx(i)     &
          + (densP*GzPT_(iP) - densM*GzPT_(iM))*nz(i)     &
          + mu*(densP + densM)*(rhotP/densP - rhotM/densM) )
      else
        dDiffFluxU  = 0.0_RP
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

  subroutine atm_dyn_nonhydro2d_cal_grad_diffVars( GxU_, GzU_, GxW_, GzW_, GxPT_, GzPT_,   &
    DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                       &
    Dx, Dz, Lift, lmesh, elem )
    
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dz, Lift
    real(RP), intent(out)  :: GxU_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzU_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GxW_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzW_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GxPT_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzPT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)

    integer :: k
    real(RP) :: DENS_(elem%Np), U_(elem%Np), W_(elem%Np), DTHETA_(elem%Np), RHOT_(elem%Np), RHOT_hyd(elem%Np)
    real(RP) :: Fx(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,AUX_DIFFVARS_NUM)
    !------------------------------------------------------------------------------

    call cal_del_gradDiffVar( del_flux,                           & ! (out)
      DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,           & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2),             & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                   & ! (in)
      lmesh, elem )                                                 ! (in)

    do k=lmesh%NeS, lmesh%NeE
      DENS_(:) = DDENS_(:,k) + DENS_hyd(:,k)
      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(:,k)/PRES00)**(CVdry/CPdry) 
      RHOT_(:) = RHOT_hyd(:) + DRHOT_(:,k)

      U_(:) = MOMX_(:,k)/DENS_(:)
      W_(:) = MOMZ_(:,k)/DENS_(:)
      DTHETA_(:) = RHOT_(:)/DENS_(:) - RHOT_hyd(:)/DENS_hyd(:,k)

      call sparsemat_matmul(Dx, U_, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_GxU_ID), LiftDelFlx)
      GxU_(:,k) = lmesh%Escale(:,k,1,1)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, U_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_GzU_ID), LiftDelFlx)
      GzU_(:,k) = lmesh%Escale(:,k,2,2)*Fz(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dx, W_, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_GxW_ID), LiftDelFlx)
      GxW_(:,k) = lmesh%Escale(:,k,1,1)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, W_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_GzW_ID), LiftDelFlx)
      GzW_(:,k) = lmesh%Escale(:,k,2,2)*Fz(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dx, DTHETA_, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_GxPT_ID), LiftDelFlx)
      GxPT_(:,k) = lmesh%Escale(:,k,1,1)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, DTHETA_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_GzPT_ID), LiftDelFlx)
      GzPT_(:,k) = lmesh%Escale(:,k,2,2)*Fz(:) + LiftDelFlx(:)
    end do

    return
  end subroutine atm_dyn_nonhydro2d_cal_grad_diffVars

  subroutine cal_del_gradDiffVar( del_flux, &
    DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, nx, nz, vmapM, vmapP, lmesh, elem )
    
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne,AUX_DIFFVARS_NUM)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    real(RP) :: densM, densP, rhot_hyd, rhotM, rhotP
    real(RP) :: delU, delW, delPT
    real(RP) :: MOMZ_P, MOMX_P
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
      delW = 0.5_RP*(MOMZ_(iP)/densP - MOMZ_(iM)/densM)
      delPT = 0.5_RP*(rhotP/densP - rhotM/densM)

      del_flux(i,VARS_GxU_ID) = delU * nx(i)
      del_flux(i,VARS_GzU_ID) = delU * nz(i)
      del_flux(i,VARS_GxW_ID) = delW * nx(i)
      del_flux(i,VARS_GzW_ID) = delW * nz(i)
      del_flux(i,VARS_GxPT_ID) = delPT * nx(i)
      del_flux(i,VARS_GzPT_ID) = delPT * nz(i)
    end do

    return
  end subroutine cal_del_gradDiffVar

end module scale_atm_dyn_nonhydro2d
