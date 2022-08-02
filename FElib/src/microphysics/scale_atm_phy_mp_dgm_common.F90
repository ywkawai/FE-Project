!> module Atmosphere / Physics cloud microphysics / common
!!
!! @par Description
!!      cloud microphysics process
!!      common subroutines
!!  
!! To preserve nonnegativity in precipitation process, 
!! a limiter proposed by Light and Durran (2016, MWR) is used
!!
!! @par Reference
!!  - Light and Durran 2016: 
!!    Preserving Nonnegativity in Discontinuous Galerkin Approximations to Scalar Transport via Truncation and Mass Aware Rescaling (TMAR).
!!    Monthly Weather Review, 144(12), 4771–4786.
!!
!! @author Team SCALE
!!
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_mp_dgm_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    UNDEF => CONST_UNDEF8, &
    GRAV => CONST_GRAV,    &
    PRES00 => CONST_PRE00
  
  use scale_sparsemat
  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
  use scale_localmesh_3d, only: &
    LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  public :: atm_phy_mp_dgm_common_gen_intweight
  public :: atm_phy_mp_dgm_common_gen_vmap
  public :: atm_phy_mp_dgm_precipitation
  public :: atm_phy_mp_dgm_precipitation_momentum

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  private :: atm_phy_mp_dgm_netOutwardFlux
  private :: atm_phy_mp_dgm_precipitation_get_delflux
  private :: atm_phy_mp_dgm_precipitation_momentum_get_delflux

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------


contains

!OCL SERIAL
  subroutine atm_phy_mp_dgm_common_gen_intweight( &
    intWeight,   & ! (out)
    lcmesh       ) ! (in)
    use scale_polynominal, only: Polynominal_GenGaussLobattoPtIntWeight
    implicit none
    class(LocalMesh3D), target :: lcmesh
    real(RP), intent(out) :: IntWeight(lcmesh%refElem3D%Nfaces,lcmesh%refElem3D%NfpTot)

    class(ElementBase3D), pointer :: elem
    real(RP), allocatable :: intWeight_lgl1DPts_h(:)
    real(RP), allocatable :: intWeight_lgl1DPts_v(:)   
    real(RP), allocatable :: intWeight_h(:) 
    real(RP), allocatable :: intWeight_v(:)  
    
    integer :: f
    integer :: i, j, k, l
    integer :: is, ie
    !--------------------------------------------

    elem => lcmesh%refElem3D
    IntWeight(:,:) = 0.0_RP

    allocate( intWeight_lgl1DPts_h(elem%Nnode_h1D) )
    allocate( intWeight_lgl1DPts_v(elem%Nnode_v) )
    allocate( intWeight_h(elem%Nnode_h1D*elem%Nnode_v) )
    allocate( intWeight_v(elem%Nnode_h1D**2) )

    intWeight_lgl1DPts_h(:) = Polynominal_GenGaussLobattoPtIntWeight(elem%PolyOrder_h)
    intWeight_lgl1DPts_v(:) = Polynominal_GenGaussLobattoPtIntWeight(elem%PolyOrder_v)

    do f=1, elem%Nfaces_h
      do k=1, elem%Nnode_v
      do i=1, elem%Nnode_h1D
        l = i + (k-1)*elem%Nnode_h1D
        intWeight_h(l) = intWeight_lgl1DPts_h(i) * intWeight_lgl1DPts_v(k)
      end do
      end do

      is = (f-1)*elem%Nfp_h + 1
      ie = is + elem%Nfp_h - 1
      IntWeight(f,is:ie) = intWeight_h(:)
    end do

    do f=1, elem%Nfaces_v
      do j=1, elem%Nnode_h1D
      do i=1, elem%Nnode_h1D
        l = i + (j-1)*elem%Nnode_h1D
        intWeight_v(l) = intWeight_lgl1DPts_h(i) * intWeight_lgl1DPts_h(j)
      end do
      end do

      is = elem%Nfaces_h*elem%Nfp_h + (f-1)*elem%Nfp_v + 1
      ie = is + elem%Nfp_v - 1
      IntWeight(elem%Nfaces_h+f,is:ie) = intWeight_v(:)
    end do

    return
  end subroutine atm_phy_mp_dgm_common_gen_intweight

!OCL SERIAL
  subroutine atm_phy_mp_dgm_common_gen_vmap( &
    vmapM, vmapP, & ! (out)
    lmesh, elem   ) ! (in)
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(out) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(out) :: vmapP(elem%NfpTot,lmesh%NeZ)    

    integer :: ke_z
    integer :: f
    integer :: vs, ve
    !-------------------------------------------------------

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

    return
  end subroutine atm_phy_mp_dgm_common_gen_vmap

!OCL SERIAL
  subroutine atm_phy_mp_dgm_precipitation( &
    DENS, RHOQ, CPtot, CVtot, RHOE,         & ! (inout)
    FLX_hydro, sflx_rain, sflx_snow, esflx, & ! (inout)
    TEMP, vterm, dt, rnstep,                & ! (in)
    Dz, Lift, nz, vmapM, vmapP, IntWeight,  & ! (in)
    QHA, QLA, QIA, lcmesh, elem             ) ! (in)
    
    use scale_atmos_hydrometeor, only: &
       CV_WATER, &
       CP_WATER, &
       CV_ICE,   &
       CP_ICE
    
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: QHA                   !< hydrometeor (water + ice)
    real(RP), intent(inout) :: DENS (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: RHOQ (elem%Np,lcmesh%NeZ,lcmesh%Ne2D,QHA)
    real(RP), intent(inout) :: CPtot(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: CVtot(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: RHOE (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: FLX_hydro(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(inout) :: sflx_rain(elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(inout) :: sflx_snow(elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(inout) :: esflx    (elem%Nfp_v,lcmesh%Ne2DA)
    real(RP), intent(in) :: TEMP (elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: vterm(elem%Np,lcmesh%NeZ,lcmesh%Ne2D,QHA)
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: rnstep
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    real(RP), intent(in) :: nz(elem%NfpTot,lcmesh%NeZ,lcmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lcmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lcmesh%NeZ)
    real(RP), intent(in) :: IntWeight(elem%Nfaces,elem%NfpTot)
    integer, intent(in) :: QLA, QIA

    real(RP) :: qflx(elem%Np)
    real(RP) :: eflx(elem%Np)
    real(RP) :: DENS0(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOCP(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOCV(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: NDcoefEuler(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: DzRHOQ(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: DzRHOE(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: dDENS(elem%Np)
    real(RP) :: CP(QHA)
    real(RP) :: CV(QHA)

    real(RP) :: fct_coef(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: RHOQ0, RHOQ1, RHOQ_tmp(elem%Np)
    real(RP) :: netOutwardFlux(lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: del_flux(elem%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,2)

    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: RHOQ_save(elem%Np)

    integer :: ke2D
    integer :: ke_z
    integer :: ke
    integer :: iq

    real(RP) :: Q
    real(RP) :: delz
    !-------------------------------------------------------

    do iq = 1, QHA
      if ( iq > QLA + QIA ) then
        CP(iq) = UNDEF 
        CV(iq) = UNDEF
      else if ( iq > QLA ) then ! ice water
        CP(iq) = CP_ICE
        CV(iq) = CV_ICE
      else                      ! liquid water
        CP(iq) = CP_WATER
        CV(iq) = CV_WATER
      end if
    end do

    !$omp parallel do collapse(2)
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ
      DENS0(:,ke_z,ke2D) = DENS(:,ke_z,ke2D)
      RHOCP(:,ke_z,ke2D) = CPtot(:,ke_z,ke2D) * DENS(:,ke_z,ke2D)
      RHOCV(:,ke_z,ke2D) = CVtot(:,ke_z,ke2D) * DENS(:,ke_z,ke2D)
    end do
    end do

    do iq = 1, QHA
      call atm_phy_mp_dgm_precipitation_get_delflux_dq( &
        del_flux(:,:,:,:),                                                      & ! (out)
        DENS0(:,:,:), RHOQ(:,:,:,iq), TEMP(:,:,:), CV(iq), nz(:,:,:), vmapM(:,:), vmapP(:,:), & ! (in)
        lcmesh, elem                                                            ) ! (in)
      
      !$omp parallel do private( &
      !$omp ke2D, ke_z, ke, delz, Fz, LiftDelFlx )
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D
        delz = ( lcmesh%pos_ev(lcmesh%EToV(ke,5),3) - lcmesh%pos_ev(lcmesh%EToV(ke,1),3) ) / dble( elem%Nnode_v )

        call sparsemat_matmul( Dz, RHOQ(:,ke_z,ke2D,iq), Fz )
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,1), LiftDelFlx )
        DzRHOQ(:,ke_z,ke2D) = lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

        call sparsemat_matmul( Dz, RHOQ(:,ke_z,ke2D,iq) * CV(iq) * TEMP(:,ke_z,ke2D), Fz )
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,2), LiftDelFlx )
        DzRHOE(:,ke_z,ke2D) = lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)
        
        NDcoefEuler(:,ke_z,ke2D) = 0.5_RP * delz * abs(vterm(:,ke_z,ke2D,iq))
      end do
      end do

      call atm_phy_mp_dgm_netOutwardFlux( &
        netOutwardFlux(:,:),                                                  & ! (out)
        RHOQ(:,:,:,iq), vterm(:,:,:,iq), DzRHOQ(:,:,:), NDcoefEuler(:,:,:),   & ! (in) 
        lcmesh%J(:,:), lcmesh%Fscale(:,:),                                    & ! (in)
        nz(:,:,:), vmapM(:,:), vmapP(:,:), lcmesh%VMapM(:,:), IntWeight(:,:), & ! (in)
        lcmesh, elem                                                          ) ! (in)
      
      !$omp parallel do collapse(2) private(ke, Q)
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D

        Q = sum( lcmesh%J(:,ke) * elem%IntWeight_lgl(:) * RHOQ(:,ke_z,ke2D,iq) ) / dt      
        fct_coef(:,ke_z,ke2D) = max( 0.0_RP, min( 1.0_RP, Q / ( netOutwardFlux(ke_z,ke2D) + 1.0E-10_RP ) ) )  
      end do ! end loop for ke_z
      end do ! end loop for ke2D

      call atm_phy_mp_dgm_precipitation_get_delflux( &
        del_flux(:,:,:,:),                                           & ! (out)
        DENS0(:,:,:), RHOQ(:,:,:,iq), TEMP(:,:,:), vterm(:,:,:,iq),  & ! (in)
        DzRHOQ(:,:,:), DzRHOE(:,:,:), NDcoefEuler(:,:,:),            & ! (in)
        fct_coef(:,:,:),                                             & ! (in)
        CV(iq), lcmesh%J(:,:), lcmesh%Fscale(:,:), nz(:,:,:),        & ! (in)
        vmapM(:,:), vmapP(:,:), lcmesh%vmapM(:,:), IntWeight(:,:),   & ! (in)
        lcmesh, elem                                                 ) ! (in)

      !$omp parallel do collapse(2) private( &
      !$omp ke2D, ke_z, ke,                            &
      !$omp qflx, eflx, dDENS, RHOQ_tmp, RHOQ0, RHOQ1, &
      !$omp RHOQ_save, &
      !$omp Fz, LiftDelFlx )
      do ke2D = 1, lcmesh%Ne2D
      do ke_z = 1, lcmesh%NeZ
        ke = ke2D + (ke_z-1)*lcmesh%Ne2D

        !--- update falling tracer 

        RHOQ_save(:) = RHOQ(:,ke_z,ke2D,iq)
        qflx(:) = vterm(:,ke_z,ke2D,iq) * RHOQ(:,ke_z,ke2D,iq)   &
                - NDcoefEuler(:,ke_z,ke2D) * DzRHOQ(:,ke_z,ke2D)

        call sparsemat_matmul( Dz, qflx(:), Fz )
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,1), LiftDelFlx )
        
        dDENS(:) =  - dt * ( &
          lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )
        RHOQ_tmp(:) = max( 0.0_RP, RHOQ(:,ke_z,ke2D,iq) + dDENS(:) )
!        RHOQ_tmp(:) = RHOQ(:,ke_z,ke2D,iq) + dDENS(:)

        !
        RHOQ0 = sum( lcmesh%Gsqrt(:,ke) * lcmesh%J(:,ke) * elem%IntWeight_lgl(:) * ( RHOQ(:,ke_z,ke2D,iq) + dDENS(:) ) )
        RHOQ1 = sum( lcmesh%Gsqrt(:,ke) * lcmesh%J(:,ke) * elem%IntWeight_lgl(:) * RHOQ_tmp(:)          )

        dDENS(:) = RHOQ0 / ( RHOQ1 + 1.0E-32_RP ) * RHOQ_tmp(:) &
                 - RHOQ(:,ke_z,ke2D,iq)
        RHOQ(:,ke_z,ke2D,iq) = RHOQ(:,ke_z,ke2D,iq) + dDENS(:)


        ! QTRC(iq; iq>QLA+QLI) is not mass tracer, such as number density
        if ( iq > QLA + QIA ) cycle

        FLX_hydro(:,ke_z,ke2D) = FLX_hydro(:,ke_z,ke2D) &
                               + qflx(:) * rnstep
        if ( ke_z == 1 ) then
          if ( iq > QLA ) then ! ice water
              sflx_snow(:,ke2D) = sflx_snow(:,ke2D)               &
                                + qflx(elem%Hslice(:,1)) * rnstep
          else                 ! liquid water
              sflx_rain(:,ke2D) = sflx_rain(:,ke2D)               &
                                + qflx(elem%Hslice(:,1)) * rnstep
          end if
        end if

        !--- update density

        RHOCP(:,ke_z,ke2D) = RHOCP(:,ke_z,ke2D) + CP(iq) * dDENS(:)
        RHOCV(:,ke_z,ke2D) = RHOCV(:,ke_z,ke2D) + CV(iq) * dDENS(:)
        DENS (:,ke_z,ke2D) = DENS(:,ke_z,ke2D) + dDENS(:)
  
        !--- update internal energy   

        eflx(:) = vterm(:,ke_z,ke2D,iq) * RHOQ_save(:) * TEMP(:,ke_z,ke2D) * CV(iq)  &
                - NDcoefEuler(:,ke_z,ke2D) * DzRHOE(:,ke_z,ke2D)
        
        call sparsemat_matmul( Dz, eflx(:), Fz )
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,2), LiftDelFlx )

        RHOE(:,ke_z,ke2D) = RHOE(:,ke_z,ke2D) - dt * ( &
          + lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) &
          + qflx(:) * Grav                                  )

        if ( ke_z == 1 ) then
          esflx(:,ke2D) = esflx(:,ke2D) &
                        + eflx(elem%Hslice(:,1)) * rnstep
        end if

      end do ! end loop for ke_z
      end do ! end loop for ke2D
    end do ! end loop for iq

    !$omp parallel do collapse(2)
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ  
      CPtot(:,ke_z,ke2D) = RHOCP(:,ke_z,ke2D) / DENS(:,ke_z,ke2D)
      CVtot(:,ke_z,ke2D) = RHOCV(:,ke_z,ke2D) / DENS(:,ke_z,ke2D)
    end do
    end do

    return
  end subroutine atm_phy_mp_dgm_precipitation

  !OCL SERIAL
  subroutine atm_phy_mp_dgm_precipitation_momentum( &
    MOMU_t, MOMV_t, MOMZ_t,                & ! (out)
    DENS, MOMU, MOMV, MOMZ, mflx,         & ! (in)
    Dz, Lift, nz, vmapM, vmapP,            & ! (in)
    lcmesh, elem                           ) ! (in)
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: MOMU_t(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMV_t(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMZ_t(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMU(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMV(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: MOMZ(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: mflx(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    type(SparseMat), intent(in) :: Dz
    type(SparseMat), intent(in) :: Lift
    real(RP), intent(in) :: nz(elem%NfpTot,lcmesh%NeZ,lcmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lcmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lcmesh%NeZ)
    
    integer :: ke2D
    integer :: ke_z
    integer :: ke

    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,3)

    real(RP) :: RDENS(elem%Np)
    !-------------------------------------------------------

    call atm_phy_mp_dgm_precipitation_momentum_get_delflux( &
      del_flux(:,:,:,:),                                  & ! (out)
      DENS(:,:,:), MOMU(:,:,:), MOMV(:,:,:), MOMZ(:,:,:), & ! (in)
      mflx(:,:,:),                                        & ! (in)
      nz(:,:,:), vmapM(:,:), vmapP(:,:),                  & ! (in)
      lcmesh, elem                                        ) ! (in)

    !$omp parallel do collapse(2) private( &
    !$omp ke2D, ke_z, ke, RDENS, Fz, LiftDelFlx )
    do ke2D = 1, lcmesh%Ne2D
    do ke_z = 1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      RDENS(:) = 1.0_RP / DENS(:,ke_z,ke2D)

      call sparsemat_matmul( Dz, mflx(:,ke_z,ke2D) * MOMU(:,ke_z,ke2D) * RDENS(:), Fz )
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,1), LiftDelFlx )
      MOMU_t(:,ke) = - ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )

      call sparsemat_matmul( Dz, mflx(:,ke_z,ke2D) * MOMV(:,ke_z,ke2D) * RDENS(:), Fz )
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,2), LiftDelFlx )
      MOMV_t(:,ke) = - ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )

      call sparsemat_matmul( Dz, mflx(:,ke_z,ke2D) * MOMZ(:,ke_z,ke2D) * RDENS(:), Fz )
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke2D,3), LiftDelFlx )
      MOMZ_t(:,ke) = - ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )
    end do
    end do

    return
  end subroutine atm_phy_mp_dgm_precipitation_momentum

!- private --------------------------------

!OCL SERIAL
  subroutine atm_phy_mp_dgm_netOutwardFlux( &
    net_outward_flux,                     &
    RHOQ_, vterm_,                        &
    DzRHOQ_, NDcoefEuler_,                &
    J, Fscale,                            &
    nz, vmapM, vmapP, vmapM3D, IntWeight, &
    lmesh, elem                           )
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: net_outward_flux(lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: RHOQ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: vterm_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: DzRHOQ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: NDcoefEuler_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: J(elem%Np*lmesh%Ne)
    real(RP), intent(in) :: Fscale(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM3D(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: IntWeight(elem%Nfaces,elem%NfpTot)

    real(RP) :: numflux(elem%NfpTot)
    real(RP) :: outward_flux_tmp(elem%Nfaces)    
    real(RP) :: alpha(elem%NfpTot)
    real(RP) :: velM(elem%NfpTot), velP(elem%NfpTot)
    real(RP) :: RHOQ_M(elem%NfpTot), RHOQ_P(elem%NfpTot)

    integer :: ke
    integer :: ke_z, ke2D
    integer :: p, i
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: iM3D(elem%NfpTot)
    !------------------------------------------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke2D, ke_z, ke, iM3D, iM, iP,      &
    !$omp RHOQ_M, RHOQ_P, velM, velP, alpha, &
    !$omp numflux, outward_flux_tmp          )
    do ke2D=1, lmesh%Ne2D
    do ke_z=1, lmesh%NeZ
      ke = ke2D + (ke_z-1)*lmesh%Ne2D

      iM3D(:) = vmapM3D(:,ke)
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      RHOQ_M(:) = RHOQ_(iM(:),ke2D)
      RHOQ_P(:) = RHOQ_(iP(:),ke2D)
      velM(:) = vterm_(iM(:),ke2D) * nz(:,ke_z,ke2D)
      velP(:) = vterm_(iP(:),ke2D) * nz(:,ke_z,ke2D)
      alpha(:) = nz(:,ke_z,ke2D)**2 * max( abs(velM(:)), abs(velP(:)) )

      where (nz(:,ke_z,ke2D) > 1.0E-10 .and. iP(:) == iM(:) )
        velP(:) = - velM(:)
      end where      

      numflux(:) = 0.5_RP * (  RHOQ_P(:) * velP(:) + RHOQ_M(:) * velM(:)                                                        &
        - ( NDcoefEuler_(iP(:),ke2D) * DzRHOQ_(iP(:),ke2D) + NDcoefEuler_(iM(:),ke2D) * DzRHOQ_(iM(:),ke2D) ) * nz(:,ke_z,ke2D) &
        - alpha(:) * ( RHOQ_P(:) - RHOQ_M(:) )                                                                                  )

      outward_flux_tmp(:) = matmul( IntWeight(:,:), J(iM3D(:)) * Fscale(:,ke) * numflux(:) )
      net_outward_flux(ke_z,ke2D) = sum( max( 0.0_RP, outward_flux_tmp(:) ) )
    end do
    end do

    return
  end subroutine atm_phy_mp_dgm_netOutwardFlux

!OCL SERIAL
  subroutine atm_phy_mp_dgm_precipitation_get_delflux_dq( &
    del_flux,                                             &
    DENS_, RHOQ_,TEMP_, CV, nz, vmapM, vmapP,             &
    lmesh, elem                                           )
    
    implicit none
    
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: del_flux(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D,2)
    real(RP), intent(in) :: DENS_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: RHOQ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: TEMP_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: CV
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)

    integer :: ke
    integer :: ke_z, ke2D
    real(RP) :: RHOQ_P(elem%NfpTot), RHOQ_M(elem%NfpTot)
    integer :: iM(elem%NfpTot), iP(elem%NfpTot)
    !-----------------------------------------

    !$omp parallel do collapse(2) private(       &
    !$omp ke2D, ke_z, ke, iM, iP, RHOQ_M, RHOQ_P )
    do ke2D=1, lmesh%Ne2D
    do ke_z=1, lmesh%NeZ
      ke = ke2D + (ke_z-1)*lmesh%Ne2D
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      RHOQ_M(:) = RHOQ_(iM(:),ke2D) 
      RHOQ_P(:) = RHOQ_(iP(:),ke2D) 
      del_flux(:,ke_z,ke2D,1) = 0.5_RP * ( RHOQ_P(:) - RHOQ_M(:) ) * nz(:,ke_z,ke2D)
      del_flux(:,ke_z,ke2D,2) = 0.5_RP * CV * ( RHOQ_P(:) * TEMP_(iP(:),ke2D) - RHOQ_M(:) * TEMP_(iM(:),ke2D) ) * nz(:,ke_z,ke2D)
    end do
    end do
    
    return
  end subroutine atm_phy_mp_dgm_precipitation_get_delflux_dq


!OCL SERIAL
  subroutine atm_phy_mp_dgm_precipitation_get_delflux( &
    del_flux,                                          &
    DENS_, RHOQ_, TEMP_, vterm_,                       &
    DzRHOQ_, DzRHOE_, NDcoefEuler_,                    &
    fct_coef_, CV,                                     &
    J, Fscale, nz, vmapM, vmapP, vmapM3D, IntWeight,   &
    lmesh, elem                                        )

    implicit none
    
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: del_flux(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D,2)
    real(RP), intent(in) :: DENS_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: RHOQ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: TEMP_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: vterm_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: DzRHOQ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: DzRHOE_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: NDcoefEuler_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: fct_coef_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: CV
    real(RP), intent(in) :: J(elem%Np*lmesh%Ne)
    real(RP), intent(in) :: Fscale(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM3D(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: IntWeight(elem%Nfaces,elem%NfpTot)

    integer :: ke
    integer :: ke_z, ke2D
    integer :: f, p, fp
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: alpha(elem%NfpTot)
    real(RP) :: velM(elem%NfpTot), velP(elem%NfpTot)
    real(RP) :: RHOQ_M(elem%NfpTot), RHOQ_P(elem%NfpTot)
    real(RP) :: TEMP_M(elem%NfpTot), TEMP_P(elem%NfpTot)

    integer :: iM3D(elem%NfpTot)
    real(RP) :: R_M(elem%NfpTot), R_P(elem%NfpTot)
    real(RP) :: NDcoef_M(elem%NfpTot), NDcoef_P(elem%NfpTot)
    real(RP) :: numflux   (elem%NfpTot)
    real(RP) :: numflux_ei(elem%NfpTot)
    real(RP) :: outward_flux_tmp(elem%Nfaces)
    real(RP) :: R
    !-----------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke2D, ke_z, ke, iM3D, iM, iP,                 &
    !$omp R_M, R_P, RHOQ_M, RHOQ_P, TEMP_M, TEMP_P,     &
    !$omp velM, velP, alpha, numflux, numflux_ei,       &
    !$omp NDcoef_M, NDcoef_P,                           &
    !$omp outward_flux_tmp,                             &
    !$omp f, p, fp, R                                   )
    do ke2D=1, lmesh%Ne2D
    do ke_z=1, lmesh%NeZ
      ke = ke2D + (ke_z-1)*lmesh%Ne2D

      iM3D(:) = vmapM3D(:,ke)
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      R_M(:) = fct_coef_(iM(:),ke2D)
      R_P(:) = fct_coef_(iP(:),ke2D) 
      RHOQ_M(:) = RHOQ_(iM(:),ke2D)
      RHOQ_P(:) = RHOQ_(iP(:),ke2D)
      TEMP_M(:) = TEMP_(iM(:),ke2D)
      TEMP_P(:) = TEMP_(iP(:),ke2D)

      velM(:) = vterm_(iM(:),ke2D) * nz(:,ke_z,ke2D)
      velP(:) = vterm_(iP(:),ke2D) * nz(:,ke_z,ke2D)
      alpha(:) = nz(:,ke_z,ke2D)**2 * max( abs(velM(:)), abs(velP(:)) )

      where (nz(:,ke_z,ke2D) > 1.0E-10 .and. iP(:) == iM(:) )
        velP(:) = - velM(:)
      end where      

      NDcoef_M(:) = NDcoefEuler_(iM(:),ke2D)
      NDcoef_P(:) = NDcoefEuler_(iP(:),ke2D)

      numflux(:) = 0.5_RP * (  RHOQ_P(:) * velP(:) + RHOQ_M(:) * velM(:)                              &
        - ( NDcoef_P(:) * DzRHOQ_(iP(:),ke2D) + NDcoef_M(:) * DzRHOQ_(iM(:),ke2D) ) * nz(:,ke_z,ke2D) &
        - alpha(:) * ( RHOQ_P(:) - RHOQ_M(:) )                                                        )

      numflux_ei(:) = 0.5_RP * ( CV * ( RHOQ_P(:) * TEMP_P(:) * velP(:) + RHOQ_M(:) * TEMP_M(:) * velM(:) ) &
        - ( NDcoef_P(:) * DzRHOE_(iP(:),ke2D) + NDcoef_M(:) * DzRHOE_(iM(:),ke2D) ) * nz(:,ke_z,ke2D)       &
        - alpha(:) * CV * ( RHOQ_P(:) * TEMP_P(:) - RHOQ_M(:) * TEMP_M(:) )                                 )


      del_flux(:,ke_z,ke2D,1) = 0.0_RP
      del_flux(:,ke_z,ke2D,2) = 0.0_RP  
      outward_flux_tmp(:) = matmul( IntWeight(:,:), J(iM3D(:)) * Fscale(:,ke) * numflux(:) )
      do f=1, elem%Nfaces_v
      do p=1, elem%Nfp_v
        fp = p + (f-1)*elem%Nfp_v + elem%Nfaces_h * elem%Nfp_h
        R = 0.5_RP * ( R_P(fp) + R_M(fp) - ( R_P(fp) - R_M(fp) ) * sign( 1.0_RP, outward_flux_tmp(elem%Nfaces_h+f) ) )
        del_flux(fp,ke_z,ke2D,1) = numflux   (fp) * R  &
                                 - RHOQ_M(fp) * velM(fp)                                  &
                                 + NDcoef_M(fp) * DzRHOQ_(iM(fp),ke2D) * nz(fp,ke_z,ke2D)
        del_flux(fp,ke_z,ke2D,2) = numflux_ei(fp) * R  &
                                 - RHOQ_M(fp) * velM(fp) * CV * TEMP_M(fp)                &
                                 + NDcoef_M(fp) * DzRHOE_(iM(fp),ke2D) * nz(fp,ke_z,ke2D)
      end do
      end do                     

    end do
    end do

    return
  end subroutine atm_phy_mp_dgm_precipitation_get_delflux

!OCL SERIAL
  subroutine atm_phy_mp_dgm_precipitation_momentum_get_delflux( &
    del_flux,                          & ! (out)
    DENS_, MOMU_, MOMV_, MOMZ_, mflx_, & ! (in)
    nz, vmapM, vmapP, lmesh, elem      ) ! (in)

    implicit none
    
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: del_flux(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D,3)
    real(RP), intent(in) :: DENS_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: MOMU_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: MOMV_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: MOMZ_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: mflx_(elem%Np*lmesh%NeZ,lmesh%Ne2D)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%Ne2D)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)

    integer :: ke_z, ke2D
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: alpha(elem%NfpTot)
    real(RP) :: densM(elem%NfpTot), densP(elem%NfpTot)
    real(RP) :: VelM(elem%NfpTot), VelP(elem%NfpTot)
    !-----------------------------------------

    !$omp parallel do collapse(2) private( &
    !$omp ke2D, ke_z, iM, iP, alpha,       &
    !$omp densM, densP, VelM, VelP         )
    do ke2D=1, lmesh%Ne2D
    do ke_z=1, lmesh%NeZ
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      densM(:) = DENS_(iM(:),ke2D)
      densP(:) = DENS_(iP(:),ke2D)
      VelM(:) = mflx_(iM(:),ke2D) * nz(:,ke_z,ke2D) / densM(:)
      VelP(:) = mflx_(iP(:),ke2D) * nz(:,ke_z,ke2D) / densP(:)
      alpha(:) = nz(:,ke_z,ke2D)**2 * max( abs(VelM(:)), abs(VelP(:)) )

      del_flux(:,ke_z,ke2D,1) = 0.5_RP * ( &
        + ( MOMU_(iP(:),ke2D) * VelP - MOMU_(iM(:),ke2D) * VelM ) &
        - alpha(:) * ( MOMU_(iP(:),ke2D) - MOMU_(iM(:),ke2D) )    )
      
      del_flux(:,ke_z,ke2D,2) = 0.5_RP * ( &
        + ( MOMV_(iP(:),ke2D) * VelP - MOMV_(iM(:),ke2D) * VelM ) &
        - alpha(:) * ( MOMV_(iP(:),ke2D) - MOMV_(iM(:),ke2D) )    )

      del_flux(:,ke_z,ke2D,3) = 0.5_RP * ( &
        + ( MOMZ_(iP(:),ke2D) * VelP - MOMZ_(iM(:),ke2D) * VelM ) &
        - alpha(:) * ( MOMZ_(iP(:),ke2D) - MOMZ_(iM(:),ke2D) )    )
    end do
    end do

    return
  end subroutine atm_phy_mp_dgm_precipitation_momentum_get_delflux

end module scale_atm_phy_mp_dgm_common
