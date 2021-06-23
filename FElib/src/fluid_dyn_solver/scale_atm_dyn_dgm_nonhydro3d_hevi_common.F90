!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics HEVI common
!!
!! @par Description
!!      HEVI DGM scheme for Atmospheric dynamical process
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_hevi_common
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
  use scale_element_modalfilter, only: ModalFilter
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    DENS_VID, MOMX_VID, MOMY_VID, MOMZ_VID, RHOT_VID, &
    PROG_VARS_NUM
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_nonhydro3d_hevi_common_gen_vmap
  public :: atm_dyn_dgm_nonhydro3d_hevi_common_eval_Ax
  public :: atm_dyn_dgm_nonhydro3d_hevi_common_construct_matbnd

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  private :: vi_cal_del_flux_dyn

contains
  !------------------------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_hevi_common_gen_vmap( &
    vmapM, vmapP, & ! (out)
    lmesh, elem   ) ! (in)
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
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
  end subroutine atm_dyn_dgm_nonhydro3d_hevi_common_gen_vmap

!OCL SERIAL  
  subroutine atm_dyn_dgm_nonhydro3d_hevi_common_eval_Ax( &
    Ax, alph,                                   & ! (out)
    PROG_VARS, PROG_VARS0, DENS_hyd, PRES_hyd,  & ! (in)
    Dz, Lift, IntrpMat_VPOrdM1,                 & ! (in)
    GnnM_z, G13_z, G23_z, GsqrtV_z,             & ! (in)
    modalFilterFlag, VModalFilter,              & ! (in)
    impl_fac, dt,                               & ! (in)
    lmesh, elem,                                & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y, cal_tend_flag ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(out) :: alph(elem%NfpTot,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)
    real(RP), intent(in) :: GnnM_z(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: G13_z(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: G23_z(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: GsqrtV_z(elem%Np,lmesh%NeZ)
    logical, intent(in) :: modalFilterFlag
    real(RP), intent(in) :: VModalFilter(elem%Nnode_v,elem%Nnode_v)    
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y
    logical, intent(in) :: cal_tend_flag

    real(RP) :: RGsqrtV(elem%Np)    
    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%NeZ,PROG_VARS_NUM)
    real(RP) :: RHOT_hyd(elem%Np), POT(elem%Np)
    real(RP) :: DPRES(elem%Np)
    integer :: ke_z
    integer :: ke, ke2d
    integer :: v
    integer :: ij
    real(RP) :: gamm, rgamm
    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    call vi_cal_del_flux_dyn( del_flux, alph,              & ! (out)
      PROG_VARS(:,DENS_VID,:), PROG_VARS(:,MOMX_VID,:),    & ! (in)
      PROG_VARS(:,MOMY_VID,:), PROG_VARS(:,MOMZ_VID,:),    & ! (in)
      PROG_VARS(:,RHOT_VID,:),                             & ! (in)
      PROG_VARS0(:,DENS_VID,:), PROG_VARS0(:,MOMX_VID,:),  & ! (in)
      PROG_VARS0(:,MOMY_VID ,:), PROG_VARS0(:,MOMZ_VID,:), & ! (in)
      PROG_VARS0(:,RHOT_VID,:),                            & ! (in)
      DENS_hyd, PRES_hyd,                                  & ! (in)
      GnnM_z, G13_z, G23_z, GsqrtV_z, nz, vmapM, vmapP,    & ! (in)
      lmesh, elem )                                          ! (in)

    !$omp parallel do private( &
    !$omp ke, ke2d, RHOT_hyd, DPRES, POT, Fz, LiftDelFlx, &
    !$omp v, ij, RGsqrtV                                  &                                      
    !$omp )
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
      ke2d = lmesh%EMap3Dto2D(ke)

      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(:,ke_z)/PRES00)**rgamm

      DPRES(:) = PRES_hyd(:,ke_z) * ( (1.0_RP + PROG_VARS(:,RHOT_VID,ke_z)/RHOT_hyd(:))**gamm - 1.0_RP )
      POT(:) = ( RHOT_hyd(:) + PROG_VARS(:,RHOT_VID,ke_z) ) / ( DENS_hyd(:,ke_z) + PROG_VARS(:,DENS_VID,ke_z) )

      RGsqrtV(:)  = 1.0_RP / GsqrtV_z(:,ke_z)

      !- DENS
      call sparsemat_matmul(Dz, PROG_VARS(:,MOMZ_VID,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,DENS_VID), LiftDelFlx)
      Ax(:,DENS_VID,ke_z) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) ) * RGsqrtV(:)

      !- MOMX
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,MOMX_VID), LiftDelFlx)
      Ax(:,MOMX_VID,ke_z) = LiftDelFlx(:) * RGsqrtV(:)

      !-MOMY
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,MOMY_VID), LiftDelFlx)
      Ax(:,MOMY_VID,ke_z) = LiftDelFlx(:) * RGsqrtV(:)

      !-MOMZ
      call sparsemat_matmul(Dz, DPRES(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,MOMZ_VID), LiftDelFlx)
      Ax(:,MOMZ_VID,ke_z) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) ) * RGsqrtV(:) &
                          + Grav * matmul(IntrpMat_VPOrdM1, PROG_VARS(:,DENS_VID,ke_z)) 

      !-RHOT
      call sparsemat_matmul(Dz, POT(:)*PROG_VARS(:,MOMZ_VID,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,RHOT_VID), LiftDelFlx)
      Ax(:,RHOT_VID,ke_z) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) ) * RGsqrtV(:)

      !-- Modal filtering in the vertical direction      
      if ( modalFilterFlag ) then
        do v=1, PROG_VARS_NUM
          do ij=1, elem%Nnode_h1D**2
            Ax(elem%Colmask(:,ij),v,ke_z) = Ax(elem%Colmask(:,ij),v,ke_z)      &
              - matmul(VModalFilter, PROG_VARS(elem%Colmask(:,ij),v,ke_z) ) / dt
          end do
        end do
      end if

      !--
      if ( .not. cal_tend_flag ) then
        Ax(:,:,ke_z) = PROG_VARS(:,:,ke_z) + impl_fac * Ax(:,:,ke_z)
      end if 

    end do    

    return
  end subroutine atm_dyn_dgm_nonhydro3d_hevi_common_eval_Ax


!OCL SERIAL  
  subroutine atm_dyn_dgm_nonhydro3d_hevi_common_construct_matbnd( &
    PmatBnd,                                & ! (out)
    kl, ku, nz_1D,                          & ! (in)
    PROG_VARS0, DENS_hyd, PRES_hyd,         & ! (in)
    G13, G23, GsqrtV, alph,                 & ! (in)
    Dz, Lift, IntrpMat_VPOrdM1,             & ! (in)
    modalFilterFlag, VModalFilter,          & ! (in)
    impl_fac, dt,                           & ! (in)
    lmesh, elem,                            & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y )            ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(in) :: kl, ku, nz_1D
    real(RP), intent(out) :: PmatBnd(2*kl+ku+1,elem%Nnode_v,PROG_VARS_NUM,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in) ::  G13(elem%Np,lmesh%NeZ)
    real(RP), intent(in) ::  G23(elem%Np,lmesh%NeZ)
    real(RP), intent(in) ::  GsqrtV(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: alph(elem%NfpTot,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)
    logical, intent(in) :: modalFilterFlag
    real(RP), intent(in) :: VModalFilter(elem%Nnode_v,elem%Nnode_v)
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y

    real(RP) :: RHOT_hyd(elem%Nnode_v)
    real(RP) :: POT0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: W0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: DENS0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP) :: DPDRHOT0(elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    integer :: ke_z, ke_z2
    integer :: v, ke, p, f1, f2, fp, fp2, FmV
    real(RP) :: gamm, rgamm
    real(RP) :: fac_dz_p(elem%Nnode_v)
    real(RP) :: PmatD(elem%Nnode_v,elem%Nnode_v,PROG_VARS_NUM,PROG_VARS_NUM)
    real(RP) :: PmatL(elem%Nnode_v,elem%Nnode_v,PROG_VARS_NUM,PROG_VARS_NUM)
    real(RP) :: PmatU(elem%Nnode_v,elem%Nnode_v,PROG_VARS_NUM,PROG_VARS_NUM)
    integer :: Colmask(elem%Nnode_v)
    real(RP) :: Id(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: Dd(elem%Nnode_v)
    real(RP) :: tmp1
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
    eval_flag(DENS_VID,MOMZ_VID ) = .true.
    eval_flag(MOMZ_VID,DENS_VID ) = .true.
    eval_flag(MOMZ_VID,RHOT_VID ) = .true.
    eval_flag(RHOT_VID,MOMZ_VID ) = .true.
    eval_flag(RHOT_VID,DENS_VID) = .true.

    Id(:,:) = 0.0_RP
    do p=1, elem%Nnode_v
      Id(p,p) = 1.0_RP
    end do

    !$omp parallel private(RHOT_hyd, Colmask)
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
                  * ( 1.0_RP + PROG_VARS0(Colmask(:),RHOT_VID,ke_z) / RHOT_hyd(:) )**(gamm-1) 

      DENS0(:,ke_z,ij) = DENS_hyd(Colmask(:),ke_z) + PROG_VARS0(Colmask(:),DENS_VID,ke_z)
      POT0(:,ke_z,ij) = ( RHOT_hyd(:) + PROG_VARS0(Colmask(:),RHOT_VID,ke_z) ) / DENS0(:,ke_z,ij)
      W0(:,ke_z,ij) = PROG_VARS0(Colmask(:),MOMZ_VID,ke_z) / DENS0(:,ke_z,ij)
    end do
    end do
    !$omp end parallel

    !$omp parallel do private(ke_z, ke, ColMask, p, fp, fp2, v, f1, f2, ke_z2, fac_dz_p, &
    !$omp fac, tmp1, FmV,                                                                &
    !$omp ij, v1, v2, pv1, pv2, pb1, g_kj, g_kjp1, g_kjm1, bc_flag,                      &
    !$omp Dd                                                                           ) &
    !$omp firstprivate(PmatD, PmatL, PmatU)
    do ij=1, elem%Nnode_h1D**2
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
      Colmask(:) = elem%Colmask(:,ij)

      !-----  
      do p=1, elem%Nnode_v
        fac_dz_p(:) = impl_fac * lmesh%Escale(Colmask(:),ke,3,3) / GsqrtV(Colmask(:),ke_z) &
                    * elem%Dx3(Colmask(:),Colmask(p))
        
        if (modalFilterFlag) then
          Dd(:) = Id(:,p) - VModalFilter(:,p) * impl_fac / dt
        else
          Dd(:) = Id(:,p)
        end if

        ! DDENS
        PmatD(:,p,DENS_VID,DENS_VID) = Dd(:)
        PmatD(:,p,DENS_VID,MOMZ_VID ) = fac_dz_p(:) 

        ! MOMX
        PmatD(:,p,MOMX_VID,MOMX_VID) = Dd(:)

        ! MOMY
        PmatD(:,p,MOMY_VID,MOMY_VID) = Dd(:)

        ! MOMZ
        PmatD(:,p,MOMZ_VID,MOMZ_VID ) = Dd(:) 
        PmatD(:,p,MOMZ_VID,DENS_VID) = impl_fac * Grav * IntrpMat_VPOrdM1(Colmask(:),Colmask(p))
        PmatD(:,p,MOMZ_VID,RHOT_VID) = fac_dz_p(:) * DPDRHOT0(p,ke_z,ij)

        !DRHOT
        PmatD(:,p,RHOT_VID,DENS_VID) = - fac_dz_p(:) * POT0(p,ke_z,ij) * W0(p,ke_z,ij)
        PmatD(:,p,RHOT_VID,MOMZ_VID ) =   fac_dz_p(:) * POT0(p,ke_z,ij)
        PmatD(:,p,RHOT_VID,RHOT_VID) = Dd(:) + fac_dz_p(:) * W0(p,ke_z,ij)
      end do

      do f1=1, 2
        if (f1==1) then
          ke_z2 = max(ke_z-1,1)
          pv1 = 1; pv2 = elem%Nnode_v
          f2 = 2
        else
          ke_z2 = min(ke_z+1,lmesh%NeZ)
          pv1 = elem%Nnode_v; pv2 = 1
          f2 = 1
        end if
        fac  = 0.5_RP * impl_fac / GsqrtV(Colmask(pv1),ke_z)
        if ( (ke_z == 1 .and. f1==1) .or. (ke_z == lmesh%NeZ .and. f1==elem%Nfaces_v) ) then
          bc_flag = .true.
          pv2 = pv1; f2 = f1
        else 
          bc_flag = .false.    
        end if

        FmV  = elem%Fmask_v(ij,f1)
        fp  = elem%Nfp_h * elem%Nfaces_h + (f1-1)*elem%Nfp_v + ij
        fp2 = elem%Nfp_h * elem%Nfaces_h + (f2-1)*elem%Nfp_v + ij

        !--
        tmp1 = fac * elem%Lift(FmV,fp) * lmesh%Fscale(fp,ke) &
             * max( alph(fp,ke_z), alph(fp2,ke_z2) )
        if (bc_flag) then
          PmatD(pv1,pv1,MOMZ_VID,MOMZ_VID) = PmatD(pv1,pv1,MOMZ_VID,MOMZ_VID) + 2.0_RP * tmp1
        else 
          do v=1, PROG_VARS_NUM
            PmatD(pv1,pv1,v,v) = PmatD(pv1,pv1,v,v) + tmp1            
            if (f1 == 1) then
              PmatL(pv1,pv2,v,v) = - tmp1                                
            else
              PmatU(pv1,pv2,v,v) = - tmp1
            end if
          end do
        end if 

        !--
        tmp1 = fac * elem%Lift(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)

        if (bc_flag) then
          PmatD(pv1,pv1,DENS_VID,MOMZ_VID ) = PmatD(pv1,pv1,DENS_VID,MOMZ_VID ) - 2.0_RP * tmp1
          PmatD(pv1,pv1,RHOT_VID,MOMZ_VID ) = PmatD(pv1,pv1,RHOT_VID,MOMZ_VID ) - 2.0_RP * tmp1 * POT0(pv1,ke_z,ij)
          PmatD(pv1,pv1,RHOT_VID,DENS_VID) = PmatD(pv1,pv1,RHOT_VID,DENS_VID) + 2.0_RP * tmp1 * POT0(pv1,ke_z,ij) * W0(pv1,ke_z,ij)
          PmatD(pv1,pv1,RHOT_VID,RHOT_VID) = PmatD(pv1,pv1,RHOT_VID,RHOT_VID) - 2.0_RP * tmp1 * W0(pv1,ke_z,ij)
        else 
          PmatD(pv1,pv1,DENS_VID,MOMZ_VID ) = PmatD(pv1,pv1,DENS_VID,MOMZ_VID ) - tmp1
          PmatD(pv1,pv1,MOMZ_VID ,RHOT_VID) = PmatD(pv1,pv1,MOMZ_VID ,RHOT_VID) - tmp1 * DPDRHOT0(pv1,ke_z,ij)
          PmatD(pv1,pv1,RHOT_VID,MOMZ_VID ) = PmatD(pv1,pv1,RHOT_VID,MOMZ_VID ) - tmp1 * POT0(pv1,ke_z,ij)
          PmatD(pv1,pv1,RHOT_VID,DENS_VID) = PmatD(pv1,pv1,RHOT_VID,DENS_VID) + tmp1 * POT0(pv1,ke_z,ij) * W0(pv1,ke_z,ij)
          PmatD(pv1,pv1,RHOT_VID,RHOT_VID) = PmatD(pv1,pv1,RHOT_VID,RHOT_VID) - tmp1 * W0(pv1,ke_z,ij)

          if (f1 == 1) then
            PmatL(pv1,pv2,DENS_VID,MOMZ_VID ) = + tmp1
            PmatL(pv1,pv2,MOMZ_VID,RHOT_VID ) = + tmp1 * DPDRHOT0(pv2,ke_z2,ij) 
            PmatL(pv1,pv2,RHOT_VID,MOMZ_VID ) = + tmp1 * POT0(pv2,ke_z2,ij)
            PmatL(pv1,pv2,RHOT_VID,DENS_VID) = - tmp1 * POT0(pv2,ke_z2,ij) * W0(pv2,ke_z2,ij)
            PmatL(pv1,pv2,RHOT_VID,RHOT_VID) = PmatL(pv1,pv2,RHOT_VID,RHOT_VID) &
                                                 + tmp1 * W0(pv2,ke_z2,ij)
          else
            PmatU(pv1,pv2,DENS_VID,MOMZ_VID ) = + tmp1
            PmatU(pv1,pv2,MOMZ_VID,RHOT_VID ) = + tmp1 * DPDRHOT0(pv2,ke_z2,ij)     
            PmatU(pv1,pv2,RHOT_VID,MOMZ_VID ) = + tmp1 * POT0(pv2,ke_z2,ij)
            PmatU(pv1,pv2,RHOT_VID,DENS_VID) = - tmp1 * POT0(pv2,ke_z2,ij) * W0(pv2,ke_z2,ij)
            PmatU(pv1,pv2,RHOT_VID,RHOT_VID) = PmatU(pv1,pv2,RHOT_VID,RHOT_VID) &
                                                 + tmp1 * W0(pv2,ke_z2,ij)
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
  end subroutine atm_dyn_dgm_nonhydro3d_hevi_common_construct_matbnd

!-- private ----------------

!OCL SERIAL  
  subroutine vi_cal_del_flux_dyn( del_flux, alph,             & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,                      & ! (in)
    DDENS0_, MOMX0_, MOMY0_, MOMZ0_, DRHOT0_,                 & ! (in)
    DENS_hyd, PRES_hyd,                                       & ! (in)
    Gnn_, G13_, G23_, GsqrtV_, nz, vmapM, vmapP, lmesh, elem  ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%NeZ,PROG_VARS_NUM)
    real(RP), intent(out) :: alph(elem%NfpTot*lmesh%NeZ)
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
    real(RP), intent(in) ::  Gnn_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  G13_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  G23_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  GsqrtV_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%NeZ)
    
    integer :: i, p, ke_z, iP, iM
    real(RP) :: MOMZ_P
    real(RP) :: wt0M, wt0P
    real(RP) :: rhot_hyd_M, rhot_hyd_P
    real(RP) :: dpresM, dpresP, densM, densP,  pottM, pottP
    real(RP) :: pres0M, pres0P, dens0M, dens0P
    real(RP) :: gamm, rgamm
    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry
    
    !$omp parallel do private( p, i, iM, iP,                  &
    !$omp rhot_hyd_M, rhot_hyd_P, densM, densP, pottM, pottP, &
    !$omp dpresM, dpresP, MOMZ_P, wt0M, wt0P,                 &
    !$omp dens0M, dens0P, pres0M, pres0P                      )    
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
   
      dpresM = PRES_hyd(iM) * ( (1.0_RP + DRHOT_(iM)/rhot_hyd_M)**gamm - 1.0_RP ) 
      dpresP = PRES_hyd(iP) * ( (1.0_RP + DRHOT_(iP)/rhot_hyd_P)**gamm - 1.0_RP ) 

      !-
      dens0M = DENS_hyd(iM) + DDENS0_(iM)
      dens0P = DENS_hyd(iP) + DDENS0_(iP)

      pres0M = PRES_hyd(iM) * (1.0_RP + DRHOT0_(iM)/rhot_hyd_M)**gamm
      pres0P = PRES_hyd(iP) * (1.0_RP + DRHOT0_(iP)/rhot_hyd_P)**gamm

      wt0M = ( MOMZ0_(iM) / GsqrtV_(iM) + G13_(iM) * MOMX0_(iM) + G23_(iM) * MOMY0_(iM) ) / dens0M      
      wt0P = ( MOMZ0_(iP) / GsqrtV_(iP) + G13_(iP) * MOMX0_(iP) + G23_(iP) * MOMY0_(iP) ) / dens0P 

      alph(i) = nz(i)**2 * max( abs( wt0M ) + sqrt( Gnn_(iM) * gamm * pres0M / dens0M ), &
                                abs( wt0P ) + sqrt( Gnn_(iP) * gamm * pres0P / dens0P )  )
      
      !----                                
      if ( iM==iP .and. (ke_z == 1 .or. ke_z == lmesh%NeZ) ) then
        MOMZ_P = - GsqrtV_(iM) * ( dens0M * wt0M + G13_(iM) * MOMX0_(iM) + G23_(iM) * MOMY0_(iM) )
        !alph(i) = 0.0_RP
      else
        MOMZ_P = MOMZ0_(iP)
      end if

      del_flux(i,DENS_VID) = 0.5_RP * (                    &
                    + ( MOMZ_P - MOMZ_(iM) ) * nz(i)        &
                    - alph(i) * ( DDENS_(iP) - DDENS_(iM) ) )
      
      del_flux(i,MOMX_VID) = 0.5_RP * (                     &
                    - alph(i) * ( MOMX_(iP) - MOMX_(iM) )  )
      
      del_flux(i,MOMY_VID) = 0.5_RP * (                     &  
                    - alph(i) * ( MOMY_(iP) - MOMY_(iM) )   )               
      
      del_flux(i,MOMZ_VID) = 0.5_RP * (                     &
                    + ( dpresP - dpresM ) * nz(i)           &                    
                    - alph(i) * ( MOMZ_P - MOMZ_(iM) )      )
      
      del_flux(i,RHOT_VID) = 0.5_RP * (                               &
                    + ( pottP * MOMZ_P  - pottM * MOMZ_(iM) ) * nz(i)  &
                    - alph(i) * ( DRHOT_(iP) - DRHOT_(iM) )            )
    end do
    end do

    return
  end subroutine vi_cal_del_flux_dyn

end module scale_atm_dyn_dgm_nonhydro3d_hevi_common
