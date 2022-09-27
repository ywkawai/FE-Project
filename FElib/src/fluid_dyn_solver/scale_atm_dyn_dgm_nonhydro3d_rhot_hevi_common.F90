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
module scale_atm_dyn_dgm_nonhydro3d_rhot_hevi_common
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


  
  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_gen_vmap
  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_eval_Ax_2
  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd_2

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
  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_gen_vmap( &
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
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_gen_vmap

!OCL SERIAL  
  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_eval_Ax_2( &
    DENS_t, MOMX_t, MOMY_t, MOMZ_t, RHOT_t,                  & ! (out)
    alph,                                                    & ! (out)
    PROG_VARS, PROG_VARS0,                                   & ! (in)
    DDENS00, MOMX00, MOMY00, MOMZ00, DRHOT00,                & ! (in)
    DENS_hyd, PRES_hyd,                                      & ! (in)
    Rtot, CPtot_ov_CVtot,                                    & ! (in)
    Dz, Lift, IntrpMat_VPOrdM1,                              & ! (in)
    GnnM, G13, G23, GsqrtV,                                  & ! (in)
    modalFilterFlag, VModalFilter,                           & ! (in)
    impl_fac, dt,                                            & ! (in)
    lmesh, elem,                                             & ! (in)
    nz, vmapM, vmapP,                                        & ! (in)
    b1D_ij, b1D_ij_uv                                        ) ! (out)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: alph(elem%NfpTot,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in)  :: PROG_VARS  (elem%Np,lmesh%NeZ,PROG_VARS_NUM,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in)  :: PROG_VARS0 (elem%Np,lmesh%NeZ,PROG_VARS_NUM,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in)  :: DDENS00(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX00 (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY00 (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ00 (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT00(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in)  :: CPtot_ov_CVtot(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)
    real(RP), intent(in) :: GnnM(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) :: G13 (elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) :: G23 (elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) :: GsqrtV(elem%Np,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    logical, intent(in) :: modalFilterFlag
    real(RP), intent(in) :: VModalFilter(elem%Nnode_v,elem%Nnode_v)    
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    real(RP), intent(out), optional :: b1D_ij(3,elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2,lmesh%NeX*lmesh%NeY)
    real(RP), intent(out), optional :: b1D_ij_uv(elem%Nnode_v,lmesh%NeZ,2,elem%Nnode_h1D**2,lmesh%NeX*lmesh%NeY)

    real(RP) :: RGsqrtV(elem%Np)
    real(RP) :: Fscale(elem%NfpTot), Escale33(elem%Np)
    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%NeZ,lmesh%NeX*lmesh%NeY,PROG_VARS_NUM)
    real(RP) :: MOMZ(elem%Np), DDENS(elem%Np), DPRES(elem%Np), RHOT(elem%Np)
    integer :: ke_xy, ke_z
    integer :: ke, ke2d
    integer :: v
    integer :: ij
    integer :: ColMask(elem%Nnode_v)
    real(RP) :: rdt
    real(RP) :: drho(elem%Np)
    
    integer :: kk, kkk, p, pp
    real(RP) :: vmf_v

    real(RP) :: gamm, rgamm, rP0
    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry
    rdt = 1.0_RP / dt
    rP0 = 1.0_RP / PRES00

    call vi_cal_del_flux_dyn( del_flux, alph,      & ! (out)
      PROG_VARS, PROG_VARS0,                       & ! (in)
      DENS_hyd, PRES_hyd,                          & ! (in)
      Rtot, CPtot_ov_CVtot,                        & ! (in)
      GnnM, G13, G23, GsqrtV, nz, vmapM, vmapP,    & ! (in)
      lmesh, elem )                                  ! (in)

    !$omp parallel private( ke_xy, ke_z, ke, ke2d, ij, v,     &
    !$omp MOMZ, DDENS, DPRES, RHOT, Fz, LiftDelFlx,           &
    !$omp RGsqrtV, ColMask, Fscale, Escale33, drho,           &
    !$omp kk, p, kkk, vmf_v, pp                               )

    !$omp do collapse(2)
    do ke_xy=1, lmesh%NeX*lmesh%NeY
    do ke_z=1, lmesh%NeZ
      ke = Ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
      ke2d = lmesh%EMap3Dto2D(ke)

      DDENS(:) = PROG_VARS(:,ke_z,DENS_VID,ke_xy)
      drho(:) = matmul(IntrpMat_VPOrdM1, DDENS(:))

      MOMZ (:) = PROG_VARS(:,ke_z,MOMZ_VID,ke_xy)

      RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke_z,ke_xy)/PRES00)**rgamm &
              + PROG_VARS(:,ke_z,RHOT_VID,ke_xy)
      DPRES(:) = PRES00 * ( Rtot(:,ke_z,ke_xy) * rP0 * RHOT(:) )**CPtot_ov_CVtot(:,ke_z,ke_xy) &
               - PRES_hyd(:,ke_z,ke_xy)
!      DPRES(:) = PRES_hyd(:,ke_z,ke_xy) * ( ( RHOT(:) / RHOT_hyd(:) )**gamm - 1.0_RP )

      RGsqrtV(:)  = 1.0_RP / GsqrtV(:,ke_z,ke_xy) 
      Fscale(:) = lmesh%Fscale(:,ke)
      Escale33(:) = lmesh%Escale(:,ke,3,3)

      !- DENS
      call sparsemat_matmul(Dz, MOMZ(:), Fz)
      call sparsemat_matmul(Lift, Fscale(:) * del_flux(:,ke_z,ke_xy,DENS_VID), LiftDelFlx)
      DENS_t(:,ke) = - ( Escale33(:) * Fz(:) + LiftDelFlx(:) ) * RGsqrtV(:)

      !- MOMX
      call sparsemat_matmul(Lift, Fscale(:) * del_flux(:,ke_z,ke_xy,MOMX_VID), LiftDelFlx)
      MOMX_t(:,ke) = - LiftDelFlx(:) * RGsqrtV(:)

      !-MOMY
      call sparsemat_matmul(Lift, Fscale(:) * del_flux(:,ke_z,ke_xy,MOMY_VID), LiftDelFlx)
      MOMY_t(:,ke) = - LiftDelFlx(:) * RGsqrtV(:)

      !-MOMZ
!      call sparsemat_matmul(Dz, MOMZ(:)**2 / ( DENS_hyd(:,ke_z,ke_xy) + DDENS(:) ) + DPRES(:), Fz) ! [<- MOMZ x MOMZ / DENS + DPRES ]
      call sparsemat_matmul(Dz, DPRES(:), Fz)
      call sparsemat_matmul(Lift, Fscale(:) * del_flux(:,ke_z,ke_xy,MOMZ_VID), LiftDelFlx)
      MOMZ_t(:,ke)  = - ( Escale33(:) * Fz(:) + LiftDelFlx(:) ) * RGsqrtV(:) &
                      - Grav * drho(:) 

      !-RHOT
      call sparsemat_matmul(Dz, RHOT(:) * MOMZ(:) / ( DENS_hyd(:,ke_z,ke_xy) + DDENS(:) ), Fz)
      call sparsemat_matmul(Lift, Fscale(:) * del_flux(:,ke_z,ke_xy,RHOT_VID), LiftDelFlx)
      RHOT_t(:,ke) = - ( Escale33(:) * Fz(:) + LiftDelFlx(:) ) * RGsqrtV(:)

    end do
    end do
    !$omp end do

    if ( modalFilterFlag ) then
      !$omp do collapse(2)
      do ke_xy=1, lmesh%NeX*lmesh%NeY
      do ke_z=1, lmesh%NeZ
        ke = ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY

          !-- Modal filtering in the vertical direction  
          do kk=1, elem%Nnode_v
          do ij=1, elem%Nnode_h1D**2
            p = elem%Colmask(kk,ij)
            do kkk=1, elem%Nnode_v
              pp = elem%Colmask(kkk,ij)
              vmf_v = rdt * VModalFilter(kk,kkk)
              DENS_t(p,ke) = DENS_t(p,ke) + vmf_v * PROG_VARS(pp,ke_z,DENS_VID,ke_xy)
              MOMX_t(p,ke) = MOMX_t(p,ke) + vmf_v * PROG_VARS(pp,ke_z,MOMX_VID,ke_xy)
              MOMY_t(p,ke) = MOMY_t(p,ke) + vmf_v * PROG_VARS(pp,ke_z,MOMY_VID,ke_xy)
              MOMZ_t(p,ke) = MOMZ_t(p,ke) + vmf_v * PROG_VARS(pp,ke_z,MOMZ_VID,ke_xy)
              RHOT_t(p,ke) = RHOT_t(p,ke) + vmf_v * PROG_VARS(pp,ke_z,RHOT_VID,ke_xy)
            end do
          end do
          end do
      end do
      end do
    !$omp end do 
    end if

    if ( present( b1D_ij ) ) then
      !$omp do collapse(2)
      do ke_xy=1, lmesh%NeX*lmesh%NeY
      do ke_z=1, lmesh%NeZ
        ke = ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
    
        do ij=1, elem%Nnode_h1D**2          
          ColMask(:) = elem%Colmask(:,ij)
          b1D_ij(1,:,ke_z,ij,ke_xy) = impl_fac * DENS_t(ColMask(:),ke)            &
                                    - PROG_VARS  (ColMask(:),ke_z,DENS_VID,ke_xy) &
                                    + DDENS00(ColMask(:),ke)
          b1D_ij(2,:,ke_z,ij,ke_xy) = impl_fac * MOMZ_t(ColMask(:),ke)            &
                                    - PROG_VARS  (ColMask(:),ke_z,MOMZ_VID,ke_xy) &
                                    + MOMZ00(ColMask(:),ke)
          b1D_ij(3,:,ke_z,ij,ke_xy) = impl_fac * RHOT_t(ColMask(:),ke)            &
                                    - PROG_VARS  (ColMask(:),ke_z,RHOT_VID,ke_xy) &
                                    + DRHOT00(ColMask(:),ke)
          b1D_ij_uv(:,ke_z,1,ij,ke_xy) = impl_fac * MOMX_t(ColMask(:),ke)            &
                                       - PROG_VARS  (ColMask(:),ke_z,MOMX_VID,ke_xy) &
                                       + MOMX00(ColMask(:),ke)
          b1D_ij_uv(:,ke_z,2,ij,ke_xy) = impl_fac * MOMY_t(ColMask(:),ke)            &
                                       - PROG_VARS  (ColMask(:),ke_z,MOMY_VID,ke_xy) &
                                       + MOMY00(ColMask(:),ke)
        end do
      end do
      end do
      !$omp end do
    end if
    !$omp end parallel

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_eval_Ax_2

!OCL SERIAL  
  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd_2( &
    PmatBnd, PmatBnd_uv,                    & ! (out)
    kl, ku, nz_1D,                          & ! (in)
    kl_uv, ku_uv, nz_1D_uv,                 & ! (in)
    PROG_VARS0, DENS_hyd, PRES_hyd,         & ! (in)
    G13, G23, GsqrtV, alph,                 & ! (in)
    Rtot, CPtot_ov_CVtot,                   & ! (in)
    Dz, Lift, IntrpMat_VPOrdM1,             & ! (in)
    modalFilterFlag, VModalFilter,          & ! (in)
    impl_fac, dt,                           & ! (in)
    lmesh, elem,                            & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y )            ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(in) :: kl, ku, nz_1D
    real(RP), intent(out) :: PmatBnd(2*kl+ku+1,3,elem%Nnode_v,lmesh%NeZ,elem%Nnode_h1D**2)
    integer, intent(in) :: kl_uv, ku_uv, nz_1D_uv
    real(RP), intent(out) :: PmatBnd_uv(2*kl_uv+ku_uv+1,elem%Nnode_v,1,lmesh%NeZ,elem%Nnode_h1D**2)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,lmesh%NeZ,PROG_VARS_NUM)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in) ::  G13(elem%Np,lmesh%NeZ)
    real(RP), intent(in) ::  G23(elem%Np,lmesh%NeZ)
    real(RP), intent(in) ::  GsqrtV(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: alph(elem%NfpTot,lmesh%NeZ)
    real(RP), intent(in) ::  Rtot(elem%Np,lmesh%NeZ)
    real(RP), intent(in) ::  CPtot_ov_CVtot(elem%Np,lmesh%NeZ)
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
    real(RP) :: PmatD(elem%Nnode_v,elem%Nnode_v,3,3)
    real(RP) :: PmatL(elem%Nnode_v,elem%Nnode_v,3,3)
    real(RP) :: PmatU(elem%Nnode_v,elem%Nnode_v,3,3)
    real(RP) :: PmatD_uv(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: PmatL_uv(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: PmatU_uv(elem%Nnode_v,elem%Nnode_v)

    integer :: Colmask(elem%Nnode_v)
    real(RP) :: Id(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: Dd(elem%Nnode_v)
    real(RP) :: tmp1
    real(RP) :: fac
    integer :: ij, v1, v2, pv1, pv2,  g_kj, g_kjp1, g_kjm1, pb1
    logical :: bc_flag
    logical :: eval_flag(3,3)

    integer, parameter :: DENS_VID_LC = 1
    integer, parameter :: MOMZ_VID_LC = 2
    integer, parameter :: RHOT_VID_LC = 3

    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    eval_flag(:,:) = .false.
    do v=1, 3
      eval_flag(v,v) = .true.
    end do
    eval_flag(DENS_VID_LC,MOMZ_VID_LC) = .true.
    eval_flag(MOMZ_VID_LC,DENS_VID_LC) = .true.
    eval_flag(MOMZ_VID_LC,RHOT_VID_LC) = .true.
    eval_flag(RHOT_VID_LC,MOMZ_VID_LC) = .true.
    eval_flag(RHOT_VID_LC,DENS_VID_LC) = .true.

    Id(:,:) = 0.0_RP
    do p=1, elem%Nnode_v
      Id(p,p) = 1.0_RP
    end do

    !$omp parallel private(RHOT_hyd, Colmask)
    !$omp workshare
    PmatD(:,:,:,:) = 0.0_RP
    PmatL(:,:,:,:) = 0.0_RP
    PmatU(:,:,:,:) = 0.0_RP  
    PmatD_uv(:,:) = 0.0_RP
    PmatL_uv(:,:) = 0.0_RP
    PmatU_uv(:,:) = 0.0_RP  
    !$omp end workshare
    !$omp do
    do ij=1, elem%Nnode_h1D**2
      PmatBnd   (:,:,:,:,ij) = 0.0_RP
      PmatBnd_uv(:,:,:,:,ij) = 0.0_RP
    end do
    !$omp do collapse(2)
    do ij=1, elem%Nnode_h1D**2
    do ke_z=1, lmesh%NeZ
      Colmask(:) = elem%Colmask(:,ij)

      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd(Colmask(:),ke_z)/PRES00)**rgamm
      DENS0(:,ke_z,ij) = DENS_hyd(Colmask(:),ke_z) + PROG_VARS0(Colmask(:),ke_z,DENS_VID)
      POT0(:,ke_z,ij) = ( RHOT_hyd(:) + PROG_VARS0(Colmask(:),ke_z,RHOT_VID) ) / DENS0(:,ke_z,ij)
      W0(:,ke_z,ij) = PROG_VARS0(Colmask(:),ke_z,MOMZ_VID) / DENS0(:,ke_z,ij)

      ! DPDRHOT0(:,ke_z,ij) = gamm * PRES_hyd(Colmask(:),ke_z) / RHOT_hyd(:)                         &
      !             * ( 1.0_RP + PROG_VARS0(Colmask(:),ke_z,RHOT_VID) / RHOT_hyd(:) )**(gamm-1.0_RP) 
      DPDRHOT0(:,ke_z,ij) = CPtot_ov_CVtot(ColMask(:),ke_z) &
                          * PRES00 * ( Rtot(ColMask(:),ke_z) / PRES00 * ( DENS0(:,ke_z,ij) * POT0(:,ke_z,ij) ) )**CPtot_ov_CVtot(ColMask(:),ke_z) &
                          / ( DENS0(:,ke_z,ij) * POT0(:,ke_z,ij) )
    end do
    end do
    !$omp end parallel

    !$omp parallel private(ke_z, ke, ColMask, p, fp, fp2, v, f1, f2, ke_z2, fac_dz_p,    &
    !$omp fac, tmp1, FmV,                                                                &
    !$omp ij, v1, v2, pv1, pv2, pb1, g_kj, g_kjp1, g_kjm1, bc_flag,                      &
    !$omp Dd                                                                           ) &
    !$omp firstprivate(PmatD, PmatL, PmatU, PmatD_uv, PmatL_uv, PmatU_uv)

    !$omp do collapse(2)
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
        PmatD(:,p,DENS_VID_LC,DENS_VID_LC) = Dd(:)
        PmatD(:,p,DENS_VID_LC,MOMZ_VID_LC) = fac_dz_p(:) 

        ! MOMX, MOMY
        PmatD_uv(:,p) = Dd(:)

        ! MOMZ
        PmatD(:,p,MOMZ_VID_LC,MOMZ_VID_LC) = Dd(:) !&
                                           !+ fac_dz_p(:) * 2.0_RP * W0(p,ke_z,ij)  [ <- d_MOMZ ( MOMZ x MOMZ / DENS ) ]
        PmatD(:,p,MOMZ_VID_LC,DENS_VID_LC) = impl_fac * Grav * IntrpMat_VPOrdM1(Colmask(:),Colmask(p)) !&
                                           !- fac_dz_p(:) * W0(p,ke_z,ij)**2        [ <- d_DENS ( MOMZ x MOMZ / DENS ) ]
        PmatD(:,p,MOMZ_VID_LC,RHOT_VID_LC) = fac_dz_p(:) * DPDRHOT0(p,ke_z,ij)

        !DRHOT
        PmatD(:,p,RHOT_VID_LC,DENS_VID_LC) = - fac_dz_p(:) * POT0(p,ke_z,ij) * W0(p,ke_z,ij)
        PmatD(:,p,RHOT_VID_LC,MOMZ_VID_LC) =   fac_dz_p(:) * POT0(p,ke_z,ij)
        PmatD(:,p,RHOT_VID_LC,RHOT_VID_LC) = Dd(:) + fac_dz_p(:) * W0(p,ke_z,ij)
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
          PmatD(pv1,pv1,MOMZ_VID_LC,MOMZ_VID_LC) = PmatD(pv1,pv1,MOMZ_VID_LC,MOMZ_VID_LC) + 2.0_RP * tmp1
        else         
          do v=1, 3
            PmatD(pv1,pv1,v,v) = PmatD(pv1,pv1,v,v) + tmp1            
            if (f1 == 1) then
              PmatL(pv1,pv2,v,v) = - tmp1                                
            else
              PmatU(pv1,pv2,v,v) = - tmp1
            end if
          end do

          PmatD_uv(pv1,pv1) = PmatD_uv(pv1,pv1) + tmp1
          if (f1 == 1) then
            PmatL_uv(pv1,pv2) = - tmp1                                
          else
            PmatU_uv(pv1,pv2) = - tmp1
          end if
        end if 

        !--
        tmp1 = fac * elem%Lift(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)

        if (bc_flag) then
          PmatD(pv1,pv1,DENS_VID_LC,MOMZ_VID_LC) = PmatD(pv1,pv1,DENS_VID_LC,MOMZ_VID_LC) - 2.0_RP * tmp1
          PmatD(pv1,pv1,RHOT_VID_LC,MOMZ_VID_LC) = PmatD(pv1,pv1,RHOT_VID_LC,MOMZ_VID_LC) - 2.0_RP * tmp1 * POT0(pv1,ke_z,ij)
          PmatD(pv1,pv1,RHOT_VID_LC,DENS_VID_LC) = PmatD(pv1,pv1,RHOT_VID_LC,DENS_VID_LC) + 2.0_RP * tmp1 * POT0(pv1,ke_z,ij) * W0(pv1,ke_z,ij)
          PmatD(pv1,pv1,RHOT_VID_LC,RHOT_VID_LC) = PmatD(pv1,pv1,RHOT_VID_LC,RHOT_VID_LC) - 2.0_RP * tmp1 * W0(pv1,ke_z,ij)
        else 
          PmatD(pv1,pv1,DENS_VID_LC,MOMZ_VID_LC) = PmatD(pv1,pv1,DENS_VID_LC,MOMZ_VID_LC) - tmp1

          ! PmatD(pv1,pv1,MOMZ_VID_LC,MOMZ_VID_LC) = PmatD(pv1,pv1,MOMZ_VID_LC,MOMZ_VID_LC) - tmp1 * 2.0_RP * W0(pv1,ke_z,ij) [ <- d_MOMZ ( MOMZ x MOMZ / DENS ) ]
          ! PmatD(pv1,pv1,MOMZ_VID_LC,DENS_VID_LC) = PmatD(pv1,pv1,MOMZ_VID_LC,DENS_VID_LC) + tmp1 * W0(pv1,ke_z,ij)**2       [ <- d_DENS ( MOMZ x MOMZ / DENS ) ]
          PmatD(pv1,pv1,MOMZ_VID_LC,RHOT_VID_LC) = PmatD(pv1,pv1,MOMZ_VID_LC,RHOT_VID_LC) - tmp1 * DPDRHOT0(pv1,ke_z,ij)

          PmatD(pv1,pv1,RHOT_VID_LC,MOMZ_VID_LC) = PmatD(pv1,pv1,RHOT_VID_LC,MOMZ_VID_LC) - tmp1 * POT0(pv1,ke_z,ij)
          PmatD(pv1,pv1,RHOT_VID_LC,DENS_VID_LC) = PmatD(pv1,pv1,RHOT_VID_LC,DENS_VID_LC) + tmp1 * POT0(pv1,ke_z,ij) * W0(pv1,ke_z,ij)
          PmatD(pv1,pv1,RHOT_VID_LC,RHOT_VID_LC) = PmatD(pv1,pv1,RHOT_VID_LC,RHOT_VID_LC) - tmp1 * W0(pv1,ke_z,ij)

          if (f1 == 1) then
            PmatL(pv1,pv2,DENS_VID_LC,MOMZ_VID_LC) = + tmp1

            ! PmatL(pv1,pv2,MOMZ_VID_LC,MOMZ_VID_LC) = PmatL(pv1,pv2,MOMZ_VID_LC,MOMZ_VID_LC) &     [ <- d_MOMZ ( MOMZ x MOMZ / DENS ) ]
            !                                      + tmp1 * 2.0_RP * W0(pv2,ke_z2,ij)
            ! PmatL(pv1,pv2,MOMZ_VID_LC,DENS_VID_LC) = - tmp1 * W0(pv2,ke_z2,ij) * W0(pv2,ke_z2,ij) [ <- d_DENS ( MOMZ x MOMZ / DENS ) ]      
            PmatL(pv1,pv2,MOMZ_VID_LC,RHOT_VID_LC) = + tmp1 * DPDRHOT0(pv2,ke_z2,ij) 

            PmatL(pv1,pv2,RHOT_VID_LC,MOMZ_VID_LC) = + tmp1 * POT0(pv2,ke_z2,ij)
            PmatL(pv1,pv2,RHOT_VID_LC,DENS_VID_LC) = - tmp1 * POT0(pv2,ke_z2,ij) * W0(pv2,ke_z2,ij)
            PmatL(pv1,pv2,RHOT_VID_LC,RHOT_VID_LC) = PmatL(pv1,pv2,RHOT_VID_LC,RHOT_VID_LC) &
                                                   + tmp1 * W0(pv2,ke_z2,ij)
          else
            PmatU(pv1,pv2,DENS_VID_LC,MOMZ_VID_LC) = + tmp1

            ! PmatU(pv1,pv2,MOMZ_VID_LC,MOMZ_VID_LC) = PmatU(pv1,pv2,MOMZ_VID_LC,MOMZ_VID_LC ) &    [ <- d_MOMZ ( MOMZ x MOMZ / DENS ) ]
            !                                     + tmp1 * 2.0_RP * W0(pv2,ke_z2,ij)
            ! PmatU(pv1,pv2,MOMZ_VID_LC,DENS_VID_LC) = - tmp1 * W0(pv2,ke_z2,ij) * W0(pv2,ke_z2,ij) [ <- d_DENS ( MOMZ x MOMZ / DENS ) ]  
            PmatU(pv1,pv2,MOMZ_VID_LC,RHOT_VID_LC) = + tmp1 * DPDRHOT0(pv2,ke_z2,ij)

            PmatU(pv1,pv2,RHOT_VID_LC,MOMZ_VID_LC) = + tmp1 * POT0(pv2,ke_z2,ij)
            PmatU(pv1,pv2,RHOT_VID_LC,DENS_VID_LC) = - tmp1 * POT0(pv2,ke_z2,ij) * W0(pv2,ke_z2,ij)
            PmatU(pv1,pv2,RHOT_VID_LC,RHOT_VID_LC) = PmatU(pv1,pv2,RHOT_VID_LC,RHOT_VID_LC) &
                                                   + tmp1 * W0(pv2,ke_z2,ij)
          end if
        end if
      end do

      do v2=1, 3
      do v1=1, 3
        if ( eval_flag(v1,v2) ) then
          do pv2=1, elem%Nnode_v
            g_kj   = v2 + (pv2-1)*3 + (ke_z-1)*elem%Nnode_v*3
            g_kjm1 = v2 + (pv2-1)*3 + (ke_z-2)*elem%Nnode_v*3
            g_kjp1 = v2 + (pv2-1)*3 + (ke_z  )*elem%Nnode_v*3

            do pv1=1, elem%Nnode_v            
              pb1 = v1 + (pv1-1)*3 + (ke_z-1)*elem%Nnode_v*3

              if (ke_z > 1 .and. pv2 == elem%Nnode_v ) then
                PmatBnd(kl+ku+1+pb1-g_kjm1, v2,pv2,ke_z-1, ij) = PmatL(pv1,pv2,v1,v2)
              end if
              PmatBnd(kl+ku+1+pb1-g_kj, v2,pv2,ke_z, ij) = PmatD(pv1,pv2,v1,v2)
              if (ke_z < lmesh%NeZ .and. pv2 == 1 ) then
                PmatBnd(kl+ku+1+pb1-g_kjp1, v2,pv2,ke_z+1, ij) = PmatU(pv1,pv2,v1,v2)
              end if
            end do
          end do            
        end if
      end do
      end do

      ! uv
      do pv2=1, elem%Nnode_v
        g_kj   = pv2 + (ke_z-1)*elem%Nnode_v
        g_kjm1 = pv2 + (ke_z-2)*elem%Nnode_v
        g_kjp1 = pv2 + (ke_z  )*elem%Nnode_v

        do pv1=1, elem%Nnode_v            
          pb1 = pv1 + (ke_z-1)*elem%Nnode_v
          if (ke_z > 1 .and. pv2 == elem%Nnode_v ) then
            PmatBnd_uv(kl_uv+ku_uv+1+pb1-g_kjm1, pv2,1,ke_z-1, ij) = PmatL_uv(pv1,pv2)
          end if
          PmatBnd_uv(kl_uv+ku_uv+1+pb1-g_kj, pv2,1,ke_z, ij) = PmatD_uv(pv1,pv2)
          if (ke_z < lmesh%NeZ .and. pv2 == 1) then
            PmatBnd_uv(kl_uv+ku_uv+1+pb1-g_kjp1, pv2,1,ke_z+1, ij) = PmatU_uv(pv1,pv2)
          end if
        end do
      end do

    end do  
    end do  
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd_2

!-- private ----------------

!OCL SERIAL  
  subroutine vi_cal_del_flux_dyn( del_flux, alph,             & ! (out)
    PVARS_, PVARS0_,                                          & ! (in)
    DENS_hyd, PRES_hyd, Rtot, CPtot_ov_CVtot,                 & ! (in)
    Gnn_, G13_, G23_, GsqrtV_, nz, vmapM, vmapP, lmesh, elem  ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%NeZ,lmesh%NeX*lmesh%NeY,PROG_VARS_NUM)
    real(RP), intent(out) :: alph(elem%NfpTot*lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) ::  PVARS_ (elem%Np*lmesh%NeZ,PROG_VARS_NUM,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) ::  PVARS0_(elem%Np*lmesh%NeZ,PROG_VARS_NUM,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) ::  Rtot(elem%Np*lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) ::  CPtot_ov_CVtot(elem%Np*lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) ::  Gnn_(elem%Np*lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) ::  G13_(elem%Np*lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) ::  G23_(elem%Np*lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) ::  GsqrtV_(elem%Np*lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%NeZ,lmesh%NeX*lmesh%NeY)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%NeZ)
    
    integer :: i, p, ke_z, iP, iM
    integer :: ij
    integer :: ke
    real(RP) :: MOMZ_P
    real(RP) :: wt0M, wt0P
    real(RP) :: rhot_hyd_M, rhot_hyd_P
    real(RP) :: dpresM, dpresP, densM, densP,  pottM, pottP
    real(RP) :: pres0M, pres0P, dens0M, dens0P

    real(RP) :: gamm, rgamm, PRES0ovRdry, rP0
    !------------------------------------------------------------------------

    gamm = CpDry / CvDry
    rgamm = CvDry / CpDry
    PRES0ovRdry = PRES00 / Rdry
    rP0 = 1.0_RP / PRES00
    
    !$omp parallel private( ij, ke_z, ke, p, i, iM, iP,       &
    !$omp rhot_hyd_M, rhot_hyd_P, densM, densP, pottM, pottP, &
    !$omp dpresM, dpresP, MOMZ_P, wt0M, wt0P,                 &
    !$omp dens0M, dens0P, pres0M, pres0P                      )  

    !$omp do collapse(2)
    do ij=1, lmesh%NeX*lmesh%NeY  
    do ke_z=1, lmesh%NeZ  
      do p=1, elem%NfpTot
        i = p + (ke_z-1)*elem%NfpTot
        iM = vmapM(i); iP = vmapP(i)

        !-
        rhot_hyd_M = PRES0ovRdry * (PRES_hyd(iM,ij)/PRES00)**rgamm
        rhot_hyd_P = PRES0ovRdry * (PRES_hyd(iP,ij)/PRES00)**rgamm

        densM = DENS_hyd(iM,ij) + PVARS_(iM,DENS_VID,ij)
        densP = DENS_hyd(iP,ij) + PVARS_(iP,DENS_VID,ij)
        
        pottM = ( rhot_hyd_M + PVARS_(iM,RHOT_VID,ij) ) / densM
        pottP = ( rhot_hyd_P + PVARS_(iP,RHOT_VID,ij) ) / densP
    
        dpresM = PRES00 * ( Rtot(iM,ij) * rP0 * densM * pottM )**CPtot_ov_CVtot(iM,ij) &
               - PRES_hyd(iM,ij) 
        dpresP = PRES00 * ( Rtot(iP,ij) * rP0 * densP * pottP )**CPtot_ov_CVtot(iP,ij) &
               - PRES_hyd(iP,ij)
        
        !-
        dens0M = DENS_hyd(iM,ij) + PVARS0_(iM,DENS_VID,ij)
        dens0P = DENS_hyd(iP,ij) + PVARS0_(iP,DENS_VID,ij)

        pres0M = PRES00 * ( Rtot(iM,ij) * rP0 * ( rhot_hyd_M + PVARS0_(iM,RHOT_VID,ij) ) )**CPtot_ov_CVtot(iM,ij)        
        pres0P = PRES00 * ( Rtot(iP,ij) * rP0 * ( rhot_hyd_P + PVARS0_(iP,RHOT_VID,ij) ) )**CPtot_ov_CVtot(iP,ij)        
        ! pres0M = PRES_hyd(iM,ij) * (1.0_RP + PVARS0_(iM,RHOT_VID,ij) / rhot_hyd_M)**gamm
        ! pres0P = PRES_hyd(iP,ij) * (1.0_RP + PVARS0_(iP,RHOT_VID,ij) / rhot_hyd_P)**gamm

        wt0M = ( PVARS0_(iM,MOMZ_VID,ij) / GsqrtV_(iM,ij) + G13_(iM,ij) * PVARS0_(iM,MOMX_VID,ij)  &
                                                          + G23_(iM,ij) * PVARS0_(iM,MOMY_VID,ij)  ) / dens0M      
        wt0P = ( PVARS0_(iP,MOMZ_VID,ij) / GsqrtV_(iP,ij) + G13_(iP,ij) * PVARS0_(iP,MOMX_VID,ij)  &
                                                          + G23_(iP,ij) * PVARS0_(iP,MOMY_VID,ij)  ) / dens0P 
        
        alph(i,ij) = nz(i,ij)**2 * max( abs( wt0M ) + sqrt( Gnn_(iM,ij) * gamm * pres0M / dens0M ), &
                                        abs( wt0P ) + sqrt( Gnn_(iP,ij) * gamm * pres0P / dens0P )  )
        
        !----                                
        if ( iM==iP .and. (ke_z == 1 .or. ke_z == lmesh%NeZ) ) then
          ! MOMZ_P = - GsqrtV_(iM,ij) * ( dens0M * wt0M + G13_(iM,ij) * PVARS0_(iM,MOMX_VID,ij) &
          !                                             + G23_(iM,ij) * PVARS0_(iM,MOMY_VID,ij) )
          MOMZ_P = - PVARS_(iM,MOMZ_VID,ij) &
                  - 2.0_RP * GsqrtV_(iM,ij) * (  G13_(iM,ij) * PVARS0_(iM,MOMX_VID,ij) &
                                               + G23_(iM,ij) * PVARS0_(iM,MOMY_VID,ij) )
          ! alph(i,ij) =  0.0_RP
        else
          MOMZ_P = PVARS_(iP,MOMZ_VID,ij)
        end if

        del_flux(i,ij,DENS_VID) = 0.5_RP * ( &
                  + ( MOMZ_P - PVARS_(iM,MOMZ_VID,ij) ) * nz(i,ij)                   &
                  - alph(i,ij) * ( PVARS_(iP,DENS_VID,ij) - PVARS_(iM,DENS_VID,ij) ) )
        
        del_flux(i,ij,MOMX_VID) = 0.5_RP * ( &
                  - alph(i,ij) * ( PVARS_(iP,MOMX_VID,ij) - PVARS_(iM,MOMX_VID,ij) ) )
        
        del_flux(i,ij,MOMY_VID) = 0.5_RP * ( &  
                  - alph(i,ij) * ( PVARS_(iP,MOMY_VID,ij) - PVARS_(iM,MOMY_VID,ij) ) )               
        
        del_flux(i,ij,MOMZ_VID) = 0.5_RP * ( &
  !                + ( MOMZ_P * MOMZ_P / densP - PVARS_(iM,MOMZ_VID,ij) * PVARS_(iM,MOMZ_VID,ij) / densM ) * nz(i,ij) & [<- MOMZ x MOMZ / DENS ]
                  + ( dpresP - dpresM ) * nz(i,ij)                                   &                    
                  - alph(i,ij) * ( MOMZ_P                 - PVARS_(iM,MOMZ_VID,ij) ) )
        
        del_flux(i,ij,RHOT_VID) = 0.5_RP * ( &
                  + ( pottP * MOMZ_P  - pottM * PVARS_(iM,MOMZ_VID,ij) ) * nz(i,ij)  &
                  - alph(i,ij) * ( PVARS_(iP,RHOT_VID,ij) - PVARS_(iM,RHOT_VID,ij) ) )
      end do
    end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine vi_cal_del_flux_dyn

end module scale_atm_dyn_dgm_nonhydro3d_rhot_hevi_common
