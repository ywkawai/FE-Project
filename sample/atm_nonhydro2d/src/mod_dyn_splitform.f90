!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_dyn
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
  use scale_element_quadrilateral
  use scale_localmesh_2d
  use scale_mesh_rectdom2d

  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D

  use mod_vars, only: &
    PROG_VARS_NUM, AUX_VARS_NUM,  &
    VARS_DDENS_ID, VARS_MOMX_ID, VARS_MOMZ_ID, VARS_DRHOT_ID, &
    VARS_GxU_ID, VARS_GzU_ID, VARS_GxW_ID, VARS_GzW_ID,       &
    VARS_GxPT_ID, VARS_GzPT_ID,                               &
    AUX_DIFFVARS_NUM

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: dyn_Init
  public :: dyn_Final
  public :: dyn_cal_tend
  public :: dyn_cal_tend2
  public :: dyn_cal_grad_diffVars
  public :: dyn_filtering

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains
  subroutine dyn_Init()
    return
  end subroutine dyn_Init

  subroutine dyn_Final()
    return
  end subroutine dyn_Final  

  !-------------------------------

  subroutine dyn_filtering( &
    DDENS_, MOMX_, MOMZ_, DRHOT_, lmesh, elem  )
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem    
    real(RP), intent(inout)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: DRHOT_(elem%Np,lmesh%NeA)
    
    integer :: k
    real(RP) :: Filter(elem%Np,elem%Np)

    !------------------------------------
    return

    call gen_filter(Filter, elem)
    do k=1, lmesh%Ne
      DDENS_(:,k) = matmul(Filter,DDENS_(:,k))
      MOMX_(:,k) = matmul(Filter,MOMX_(:,k))
      MOMZ_(:,k) = matmul(Filter,MOMZ_(:,k))
      DRHOT_(:,k) = matmul(Filter,DRHOT_(:,k))
    end do
    return
  end subroutine dyn_filtering

  subroutine dyn_cal_tend2( dQdt, &
    DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, &
    GxU_, GzU_, GxW_, GzW_, GxPT_, GzPT_,             &
    diffCoef_h, diffCoef_v,                           &
    Dx, Dz, Sx, Sz, Lift, lmesh, elem)
    
    use mod_vars, only: &
      DxMOMX, DzMOMZ, LiftDDENS
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dz, Sx, Sz, Lift
    real(RP), intent(out) :: dQdt(elem%Np,lmesh%NeA,PROG_VARS_NUM)
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
    real(RP), intent(in) :: diffCoef_h
    real(RP), intent(in) :: diffCoef_v

    real(RP) :: Fx(elem%Np,PROG_VARS_NUM), Fz(elem%Np,PROG_VARS_NUM)
    real(RP) :: LiftDelFlx(elem%Np,PROG_VARS_NUM)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: dens_(elem%Np), RHOT_(elem%Np), dpres_(elem%Np)
    real(RP) :: pres_(elem%Np), u_(elem%Np), w_(elem%Np)
    real(RP) :: dens_m(elem%Np), theta_m(elem%Np), dpres_m(elem%Np)
    real(RP) :: u_m(elem%Np), w_m(elem%Np)

    integer :: k, p

    integer :: p1, p2, p_
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

      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_DDENS_ID), LiftDelFlx(:,VARS_DDENS_ID))
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_MOMX_ID), LiftDelFlx(:,VARS_MOMX_ID))
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_MOMZ_ID), LiftDelFlx(:,VARS_MOMZ_ID))
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_DRHOT_ID), LiftDelFlx(:,VARS_DRHOT_ID))

      do p2=1, elem%Nfp
      do p1=1, elem%Nfp
        p = p1 + (p2 - 1)*elem%Nfp
        dens_m(:)  = 0.5_RP*(dens_(p)+dens_(:))
        u_m(:)     = 0.5_RP*(u_(p)+u_(:))
        w_m(:)     = 0.5_RP*(w_(p)+w_(:))
        dpres_m(:) = 0.5_RP*(dpres_(p) + dpres_(:))
        theta_m(:)  = 0.5_RP*(RHOT_(p)/dens_(p) + RHOT_(:)/dens_(:))

        ! DENS
        Fx(p,1) = 2.0_RP * sum(elem%Dx1(p,:)*dens_m(:)*u_m(:))
        Fz(p,1) = 2.0_RP * sum(elem%Dx2(p,:)*dens_m(:)*w_m(:))
        
        dQdt(p,k,VARS_DDENS_ID) = &
          - (  lmesh%Escale(p,k,1,1) * Fx(p,1) &
             + lmesh%Escale(p,k,2,2) * Fz(p,1) &
             + LiftDelFlx(p,VARS_DDENS_ID) )        
        
        ! MOMX
        Fx(p,2) = sum( elem%Dx1(p,:) * ( &
            2.0_RP*dens_m(:)*u_m(:)*u_m(:)                 &
          + 2.0_RP*dpres_m(:)                              &
          - diffCoef_h*dens_(:)*GxU_(:,k)                  &
        ))
        Fz(p,2) = sum( elem%Dx2(p,:) * ( &
            2.0_RP*dens_m(:)*u_m(:)*w_m(:)                 &
          - diffCoef_v*dens_(:)*GzU_(:,k)                  &
        ))

        dQdt(p,k,VARS_MOMX_ID) = &
        - (  lmesh%Escale(p,k,1,1) * Fx(p,2) &
           + lmesh%Escale(p,k,2,2) * Fz(p,2) &
           + LiftDelFlx(p,VARS_MOMX_ID))
                
        ! MOMZ
        Fx(p,3) = sum( elem%Dx1(p,:) * ( &
            2.0_RP*dens_m(:)*w_m(:)*u_m(:)                 &
          - diffCoef_h*dens_(:)*GxW_(:,k)                  &
        ))
        Fz(p,3) = sum( elem%Dx2(p,:) * ( &
            2.0_RP*dens_m(:)*w_m(:)*w_m(:)                 &
          + 2.0_RP*dpres_m(:)                              &
          - diffCoef_v*dens_(:)*GzW_(:,k)                  &
        ))  

        dQdt(p,k,VARS_MOMZ_ID) = &
        - (  lmesh%Escale(p,k,1,1) * Fx(p,3) &
           + lmesh%Escale(p,k,2,2) * Fz(p,3) &
           + LiftDelFlx(p,VARS_MOMZ_ID) )    &
        - DDENS_(p,k)*Grav
           !- matmul(IntrpMat_VPOrdM1, DDENS_(:,k)) * Grav

        ! RHOT
        Fx(p,4) = sum( elem%Dx1(p,:) * ( &
             2.0_RP*dens_m(:)*theta_m(:)*u_m(:)         &
           - diffCoef_h*dens_(:)*GxPT_(:,k)             &
        ))
        Fz(p,4) = sum( elem%Dx2(p,:) * ( &
             2.0_RP*dens_m(:)*theta_m(:)*w_m(:)         &
           - diffCoef_v*dens_(:)*GzPT_(:,k)             &
        ))
        dQdt(p,k,VARS_DRHOT_ID) = &
        - (  lmesh%Escale(p,k,1,1) * Fx(p,4) &
           + lmesh%Escale(p,k,2,2) * Fz(p,4) &
           + LiftDelFlx(p,VARS_DRHOT_ID) )        
      end do
      end do
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine dyn_cal_tend2

  subroutine dyn_cal_tend( dQdt, &
    DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, &
    GxU_, GzU_, GxW_, GzW_, GxPT_, GzPT_,             &
    diffCoef_h, diffCoef_v,                           &
    Dx, Dz, Sx, Sz, Lift, lmesh, elem)
    
    use mod_vars, only: &
      DxMOMX, DzMOMZ, LiftDDENS
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dz, Sx, Sz, Lift
    real(RP), intent(out) :: dQdt(elem%Np,lmesh%NeA,PROG_VARS_NUM)
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

      dQdt(:,k,VARS_DDENS_ID) = &
        - (  lmesh%Escale(:,k,1,1) * Fx(:) &
           + lmesh%Escale(:,k,2,2) * Fz(:) &
           + LiftDelFlx(:) )

      DxMOMX%local(1)%val(:,k) = - lmesh%Escale(:,k,1,1) * Fx(:)
      DzMOMZ%local(1)%val(:,k) = - lmesh%Escale(:,k,2,2) * Fz(:)
      LiftDDENS%local(1)%val(:,k) = - LiftDelFlx(:)

      !-- MOMX
      call sparsemat_matmul(Dx, u_(:)*MOMX_(:,k) + dpres_(:) - diffCoef_h*dens_(:)*GxU_(:,k), Fx)
      call sparsemat_matmul(Dz, w_(:)*MOMX_(:,k)             - diffCoef_v*dens_(:)*GzU_(:,k) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_MOMX_ID), LiftDelFlx)

      dQdt(:,k,VARS_MOMX_ID) = &
        - (  lmesh%Escale(:,k,1,1) * Fx(:) &
           + lmesh%Escale(:,k,2,2) * Fz(:) &
           + LiftDelFlx )

      !-- MOMZ
      call sparsemat_matmul(Dx, u_(:)*MOMZ_(:,k)             - diffCoef_h*dens_(:)*GxW_(:,k), Fx)
      call sparsemat_matmul(Dz, w_(:)*MOMZ_(:,k) + dpres_(:) - diffCoef_v*dens_(:)*GzW_(:,k), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_MOMZ_ID), LiftDelFlx)
      
      dQdt(:,k,VARS_MOMZ_ID) = &
        - (  lmesh%Escale(:,k,1,1) * Fx(:) &
           + lmesh%Escale(:,k,2,2) * Fz(:) &
           + LiftDelFlx )                  &
        - matmul(IntrpMat_VPOrdM1, DDENS_(:,k)) * Grav
        !- DDENS_(:,k)*Grav
        

      !-- RHOT
      call sparsemat_matmul(Dx, u_(:)*RHOT_(:) - diffCoef_h*dens_(:)*GxPT_(:,k), Fx)
      call sparsemat_matmul(Dz, w_(:)*RHOT_(:) - diffCoef_v*dens_(:)*GzPT_(:,k), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_DRHOT_ID), LiftDelFlx)
      
      dQdt(:,k,VARS_DRHOT_ID) = &
        - (  lmesh%Escale(:,k,1,1) * Fx(:) &
           + lmesh%Escale(:,k,2,2) * Fz(:) &
           + LiftDelFlx )
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine dyn_cal_tend

  subroutine gen_filter( filter, elem )
    class(elementbase2D), intent(in) :: elem
    real(RP), intent(out) :: Filter(elem%Np,elem%Np)

    real(RP) :: filter1D(elem%Nfp)
    real(RP) :: etac, eta
    integer :: p1, p2, l

    etac = 0.0_RP!(elem%PolyOrder*0.6_RP)/dble(elem%PolyOrder)
    filter1D(:) = 1.0_RP
    do p1=1, elem%Nfp
      eta = dble(p1-1)/dble(elem%PolyOrder)
      if ( eta > etac .and. p1 /= 1) then
        filter1D(p1) = exp( -  0.1_RP*( ((eta - etac)/(1.0_RP - etac))**32 ))
      end if
    end do

    Filter(:,:) = 0.0_RP
    do p2=1, elem%Nfp
    do p1=1, elem%Nfp
      l = p1 + (p2-1)*elem%Nfp
      Filter(l,l) = filter1D(p1) * filter1D(p2)
    end do  
    end do
    Filter(:,:) = matmul(Filter, elem%invV)
    Filter(:,:) = matmul(elem%V, Filter)

  end subroutine gen_filter

  subroutine cal_del_flux_dyn( del_flux, &
    DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,     &
    GxU_, GzU_, GxW_, GzW_, GxPT_, GzPT_,                 &
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
    real(RP) :: mu, diffCoef

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
      mu = (2.0_RP * dble((elem%PolyOrder+1)*(elem%PolyOrder+2)) / 2.0_RP / 1600.0_RP) 
      diffCoef = diffCoef_h*abs(nx(i)) + diffCoef_v*abs(nz(i))
      
      if (diffCoef > 0.0_RP) then
        dDiffFluxU = diffCoef*( &
            (densP*GxU_(iP) - densM*GxU_(iM))*nx(i)                    &
          + (densP*GzU_(iP) - densM*GzU_(iM))*nz(i)                    &
          + mu*(densP + densM)*(MOMX_(iP)/densP - MOMX_(iM)/densM) )
        
        dDiffFluxW = diffCoef*( &
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

  subroutine dyn_cal_grad_diffVars( GxU_, GzU_, GxW_, GzW_, GxPT_, GzPT_,   &
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
  end subroutine dyn_cal_grad_diffVars

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

end module mod_dyn
