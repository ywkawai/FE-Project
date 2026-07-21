!> module FElib / Atmosphere / Physics / boundary layer turbulence
!!
!! @par Description
!!      Boundary layer turbulence process
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_phy_bl_dgm_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    EPS  => CONST_EPS
  use scale_sparsemat  
  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_element_operation_base, only: ElementOperationBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_phy_bl_dgm_common_calc_tendency

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
contains
!OCL SERIAL
  subroutine atm_phy_bl_dgm_common_calc_tendency( &
    RHOU_tp, RHOV_tp, DRHOT_tp,            & ! (out)
    DDENS_, MOMX_, MOMY_, DRHOT_,          & ! (in)
    PT_, DENS_hyd, PRES_hyd, NU, KH,       & ! (in)
    element3D_operation, C_IP, dtsec,      & ! (in)
    lmesh, elem, elem1D, is_bound,         & ! (in)
    use_delta_form                         ) ! (in)
    use scale_atm_dyn_dgm_hevi_common_linalgebra, only: &
      atm_dyn_dgm_hevi_common_linalgebra_get_param
    implicit none
    class(LocalMesh3D), intent(in), target :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(ElementBase1D), intent(in) :: elem1D
    real(RP), intent(out) :: RHOU_tp(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOV_tp(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: DRHOT_tp(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PT_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: NU(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: KH(elem%Np,lmesh%NeA)
    class(ElementOperationBase3D), intent(in) :: element3D_operation
    real(RP), intent(in) :: C_IP
    real(RP), intent(in) :: dtsec
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)
    logical, intent(in) :: use_delta_form
    
    class(LocalMesh2D), pointer :: lmesh2D
    class(ElementBase2D), pointer :: elem2D

    real(RP) :: PROG_VARS (elem%Np,lmesh%NeX*lmesh%NeY,lmesh%NeZ,3)
    real(RP) :: alph_M(elem%NfpTot,lmesh%Ne)
    real(RP) :: alph_H(elem%NfpTot,lmesh%Ne)
    real(RP) :: GsqrtV(elem%Np,lmesh%Ne)

    integer :: vmapM(elem%NfpTot,lmesh%Ne)
    integer :: vmapP(elem%NfpTot,lmesh%Ne)
    integer :: ke_xy, ke_z, ke, ke2d
    integer :: p

    integer :: im, jm

    real(RP) :: DENS(elem%Np,lmesh%Ne)
    real(RP), allocatable :: b1D_ij(:,:,:,:,:)
    real(RP), allocatable :: BndMatL(:,:,:,:,:)
    real(RP), allocatable :: BndMatD(:,:,:,:,:)
    real(RP), allocatable :: G(:,:,:,:,:,:)

    real(RP) :: impl_fac
    !------------------------------------------------------------------------

    lmesh2D => lmesh%lcmesh2D
    elem2D => lmesh2D%refElem2D
    impl_fac = 1.0_RP * dtsec

    call lmesh%GetVmapZ3D( vmapM, vmapP ) ! (out)
    call atm_dyn_dgm_hevi_common_linalgebra_get_param( im, jm, & ! (out)
      elem%Nnode_v )

    allocate( b1D_ij(im*elem%Nnode_v,3,jm,lmesh%Ne2D,lmesh%NeZ) )
    allocate( BndMatD(im*elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D,2) )
    allocate( BndMatL(im*elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D,2) )
    allocate( G(im*elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D,lmesh%NeZ,2) )

    !$omp parallel do collapse(2) private( ke, ke2D )
    do ke_z =1, lmesh%NeZ
    do ke_xy=1, lmesh%NeX * lmesh%NeY
      ke = ke_xy + (ke_z-1)*lmesh%Ne2D
      ke2D = lmesh%EMap3Dto2D(ke)

      DENS(:,ke) = DENS_hyd(:,ke) + DDENS_(:,ke)

      PROG_VARS(:,ke_xy,ke_z,1) = MOMX_ (:,ke)
      PROG_VARS(:,ke_xy,ke_z,2) = MOMY_ (:,ke)
      PROG_VARS(:,ke_xy,ke_z,3) = DENS(:,ke) * PT_(:,ke)

      do p=1, elem%Np
        GsqrtV(p,ke)  = lmesh%Gsqrt(p,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D(p),ke2D)
      end do
    end do
    end do

    call eval_Ax( RHOU_tp, RHOV_tp, DRHOT_tp, alph_M, alph_H,   & !(out)
      PROG_VARS, MOMX_, MOMY_, PT_, NU, KH, DENS, GsqrtV,       & !(in)
      impl_fac, dtsec, lmesh, elem, vmapM, vmapP, is_bound,     & !(in)
      element3D_operation, C_IP, im, jm, b1D_ij, use_delta_form ) !(in)

    call vi_solve( PROG_VARS,                      & ! (inout)
      BndMatL, BndMatD, G, b1D_ij,                 & ! (inout)
      DENS, NU, KH, GsqrtV, C_IP, dtsec, impl_fac, & ! (in)
      im, jm, lmesh, elem, elem1D, use_delta_form  ) ! (in)

    !---
    !$omp parallel do collapse(2) private( ke )
    do ke_z =1, lmesh%NeZ
    do ke_xy=1, lmesh%NeX * lmesh%NeY
      ke = ke_xy + (ke_z-1)*lmesh%Ne2D
      RHOU_tp (:,ke) = ( PROG_VARS(:,ke_xy,ke_z,1) - MOMX_ (:,ke) ) / impl_fac
      RHOV_tp (:,ke) = ( PROG_VARS(:,ke_xy,ke_z,2) - MOMY_ (:,ke) ) / impl_fac
      DRHOT_tp(:,ke) = ( PROG_VARS(:,ke_xy,ke_z,3) - DENS(:,ke) * PT_(:,ke) ) / impl_fac
    end do
    end do

    return
  end subroutine atm_phy_bl_dgm_common_calc_tendency

!- Private subroutines -----------------------------

!OCL SERIAL
  subroutine vi_solve( PROG_VARS,                & ! (inout)
    BndMatL, BndMatD, G, b1D_ij,                 & ! (inout)
    DENS, NU, KH, GsqrtV, C_IP, dtsec, impl_fac, & ! (in)
    im, jm, lmesh, elem, elem1D, use_delta_form  ) ! (in)
    implicit none
    integer, intent(in) :: im, jm
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(ElementBase1D), intent(in) :: elem1D
    real(RP), intent(inout) :: PROG_VARS(elem%Nnode_h1D**2,elem%Nnode_v,lmesh%Ne2D,lmesh%NeZ,3)
    real(RP), intent(inout) :: BndMatL(im*elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D,2)
    real(RP), intent(inout) :: BndMatD(im*elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D,2)
    real(RP), intent(inout) :: G(im*elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D,lmesh%NeZ,2)
    real(RP), intent(inout) :: b1D_ij(im*elem%Nnode_v,3,jm,lmesh%Ne2D,lmesh%NeZ)
    real(RP), intent(in) :: DENS(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: NU(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: KH(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: GsqrtV(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: C_IP
    real(RP), intent(in) :: dtsec
    real(RP), intent(in) :: impl_fac
    logical, intent(in) :: use_delta_form

    integer :: ke_z, ke_xy
    integer :: i, j, ij
    integer :: pv, pv1, pv2, pp, p2

    real(RP) :: tmp(im*elem%Nnode_v,2)
    real(RP) :: tmp_b(im*elem%Nnode_v,3)
    !-------------------------------------------------------

    do ke_z=1, lmesh%NeZ
      call construct_matbnd_sip( &
        BndMatL(:,:,:,:,1), BndMatD(:,:,:,:,1), G(:,:,:,:,ke_z,1), & ! (out)
        DENS, NU, GsqrtV, elem1D%Dx1, elem1D%M, elem1D%invM,       & ! (in)
        impl_fac, dtsec, C_IP, lmesh, elem, im, jm, ke_z           ) ! (in)
    
      call construct_matbnd_sip( &
        BndMatL(:,:,:,:,2), BndMatD(:,:,:,:,2), G(:,:,:,:,ke_z,2), & ! (out)
        DENS, KH, GsqrtV, elem1D%Dx1, elem1D%M, elem1D%invM,       & ! (in)
        impl_fac, dtsec, C_IP, lmesh, elem, im, jm, ke_z           ) ! (in)
    
      if ( ke_z > 1 ) then
        !$omp parallel do collapse(2) private( pp, p2, tmp, tmp_b )
        do ke_xy=1, lmesh%NeX * lmesh%NeY
        do j=1, jm
          !* D_k <- D_k - L_k * G_{k-1} ------------------
          do pv2=1, elem%Nnode_v
            tmp(:,:) = 0.0_RP
            do pv=1, elem%Nnode_v
              do pv1=1, elem%Nnode_v
              do i=1, im
                pp = i + (pv1-1)*im; p2 = i + (pv -1)*im
                tmp(pp,1) = tmp(pp,1) + BndMatL(pp,pv,j,ke_xy,1) * G(p2,pv2,j,ke_xy,ke_z-1,1)
                tmp(pp,2) = tmp(pp,2) + BndMatL(pp,pv,j,ke_xy,2) * G(p2,pv2,j,ke_xy,ke_z-1,2)
              end do
              end do
            end do
            BndMatD(:,pv2,j,ke_xy,1) = BndMatD(:,pv2,j,ke_xy,1) - tmp(:,1)
            BndMatD(:,pv2,j,ke_xy,2) = BndMatD(:,pv2,j,ke_xy,2) - tmp(:,2)
          end do ! loop for pv2

          !* b_k <- b_k - L_k * b_{k-1} ------------------
          tmp_b(:,:) = 0.0_RP
          do pv=1, elem%Nnode_v
          do pv1=1, elem%Nnode_v
          do i=1, im
            pp = i + (pv1-1)*im; p2 = i + (pv -1)*im
            tmp_b(pp,1) = tmp_b(pp,1) + BndMatL(pp,pv,j,ke_xy,1) * b1D_ij(p2,1,j,ke_xy,ke_z-1)
            tmp_b(pp,2) = tmp_b(pp,2) + BndMatL(pp,pv,j,ke_xy,1) * b1D_ij(p2,2,j,ke_xy,ke_z-1)
            tmp_b(pp,3) = tmp_b(pp,3) + BndMatL(pp,pv,j,ke_xy,2) * b1D_ij(p2,3,j,ke_xy,ke_z-1)
          end do
          end do
          end do
          b1D_ij(:,1,j,ke_xy,ke_z) = b1D_ij(:,1,j,ke_xy,ke_z) - tmp_b(:,1)
          b1D_ij(:,2,j,ke_xy,ke_z) = b1D_ij(:,2,j,ke_xy,ke_z) - tmp_b(:,2)
          b1D_ij(:,3,j,ke_xy,ke_z) = b1D_ij(:,3,j,ke_xy,ke_z) - tmp_b(:,3)
        end do
        end do
      end if

      call linkernel_solve_sip( b1D_ij(:,:,:,:,ke_z), G(:,:,:,:,ke_z,1), G(:,:,:,:,ke_z,2),&
        BndMatD(:,:,:,:,1), BndMatD(:,:,:,:,2), elem%Nnode_v, im, jm, lmesh%Ne2D, ke_z == lmesh%NeZ   )

    end do
    do ke_z=lmesh%NeZ-1, 1, -1
      !$omp parallel do collapse(2) private( pp, p2, tmp_b )
      do ke_xy=1, lmesh%NeX * lmesh%NeY
      do j=1, jm
        tmp_b(:,:) = 0.0_RP
        do pv2=1, elem%Nnode_v
        do pv =1, elem%Nnode_v
        do i=1, im
          pp = i + (pv -1)*im; p2 = i + (pv2-1)*im
          tmp_b(pp,1) = tmp_b(pp,1) + G(pp,pv2,j,ke_xy,ke_z,1) * b1D_ij(p2,1,j,ke_xy,ke_z+1)
          tmp_b(pp,2) = tmp_b(pp,2) + G(pp,pv2,j,ke_xy,ke_z,1) * b1D_ij(p2,2,j,ke_xy,ke_z+1)
          tmp_b(pp,3) = tmp_b(pp,3) + G(pp,pv2,j,ke_xy,ke_z,2) * b1D_ij(p2,3,j,ke_xy,ke_z+1)
        end do
        end do
        end do
        b1D_ij(:,1,j,ke_xy,ke_z) = b1D_ij(:,1,j,ke_xy,ke_z) - tmp_b(:,1)
        b1D_ij(:,2,j,ke_xy,ke_z) = b1D_ij(:,2,j,ke_xy,ke_z) - tmp_b(:,2)
        b1D_ij(:,3,j,ke_xy,ke_z) = b1D_ij(:,3,j,ke_xy,ke_z) - tmp_b(:,3)
      end do
      end do      
    end do

    !$omp parallel do collapse(2) private( p2, pp )
    do ke_z=1, lmesh%NeZ
    do ke_xy=1, lmesh%NeX * lmesh%NeY
      do j=1, jm
        if ( use_delta_form ) then 
            do pv=1, elem%Nnode_v
            do i=1, im    
              p2 = i + (pv-1)*im; pp = i + (j-1)*im
              PROG_VARS(pp,pv,ke_xy,ke_z,1) = PROG_VARS(pp,pv,ke_xy,ke_z,1) + b1D_ij(p2,1,j,ke_xy,ke_z)
              PROG_VARS(pp,pv,ke_xy,ke_z,2) = PROG_VARS(pp,pv,ke_xy,ke_z,2) + b1D_ij(p2,2,j,ke_xy,ke_z)
              PROG_VARS(pp,pv,ke_xy,ke_z,3) = PROG_VARS(pp,pv,ke_xy,ke_z,3) + b1D_ij(p2,3,j,ke_xy,ke_z)
            end do
            end do
        else
            do pv=1, elem%Nnode_v
            do i=1, im    
              p2 = i + (pv-1)*im; pp = i + (j-1)*im
              PROG_VARS(pp,pv,ke_xy,ke_z,1) = b1D_ij(p2,1,j,ke_xy,ke_z)
              PROG_VARS(pp,pv,ke_xy,ke_z,2) = b1D_ij(p2,2,j,ke_xy,ke_z)
              PROG_VARS(pp,pv,ke_xy,ke_z,3) = b1D_ij(p2,3,j,ke_xy,ke_z)
            end do
            end do
        end if
      end do
    end do
    end do

    return
  end subroutine vi_solve

!OCL SERIAL
  subroutine linkernel_solve_sip( b, G1, G2, & ! (inout)
    BndMatD1, BndMatD2, nv, im, jm, Ne2D, top_flag  ) ! (in)
    implicit none

    integer, intent(in) :: nv
    integer, intent(in) :: im
    integer, intent(in) :: jm
    integer, intent(in) :: Ne2D
    real(RP), intent(inout) :: b(im*nv,3,jm,Ne2D)
    real(RP), intent(inout) :: G1(im*nv,nv,jm,Ne2D)
    real(RP), intent(inout) :: G2(im*nv,nv,jm,Ne2D)
    real(RP), intent(in) :: BndMatD1(im*nv,nv,jm,Ne2D)
    real(RP), intent(in) :: BndMatD2(im*nv,nv,jm,Ne2D)
    logical, intent(in) :: top_flag

    integer :: ke_xy
    integer :: i, j
    integer :: pv1, pv2
    integer :: iv
    integer :: pp

    integer :: info
    integer :: nrhs1, nrhs2
    integer :: ipiv1(nv), ipiv2(nv)

    real(RP) :: Amat1(nv,nv)
    real(RP) :: RHS1(nv,nv+2)
    real(RP) :: Amat2(nv,nv)
    real(RP) :: RHS2(nv,nv+2)
    !------------------------------------------------------------

    !$omp parallel do collapse(2) &
    !$omp private(ke_xy,j,i,pv1,pv2,iv,pp,info,                 &
    !$omp         nrhs1,nrhs2,ipiv1,ipiv2,Amat1,RHS1,Amat2,RHS2 )
    do ke_xy=1, Ne2D
    do j=1, jm
      do i=1, im

        do pv2=1, nv
        do pv1=1, nv
            pp = i + (pv1-1)*im
            Amat1(pv1,pv2) = BndMatD1(pp,pv2,j,ke_xy)
            Amat2(pv1,pv2) = BndMatD2(pp,pv2,j,ke_xy)
        end do
        end do

        RHS1(:,:) = 0.0_RP
        RHS2(:,:) = 0.0_RP

        if ( top_flag ) then
            nrhs1 = 2
            do iv=1, 2
            do pv1=1, nv
              pp = i + (pv1-1)*im
              RHS1(pv1,iv) = b(pp,iv,j,ke_xy)
            end do
            end do

            nrhs2 = 1
            do pv1=1, nv
              pp = i + (pv1-1)*im
              RHS2(pv1,1) = b(pp,3,j,ke_xy)
            end do
        else
            nrhs1 = nv + 2
            nrhs2 = nv + 1

            do pv2=1, nv
            do pv1=1, nv
              pp = i + (pv1-1)*im
              RHS1(pv1,pv2) = G1(pp,pv2,j,ke_xy)
              RHS2(pv1,pv2) = G2(pp,pv2,j,ke_xy)
            end do
            end do
            do iv=1, 2
            do pv1=1, nv
              pp = i + (pv1-1)*im
              RHS1(pv1,nv+iv) = b(pp,iv,j,ke_xy)
            end do
            end do
            do pv1=1, nv
              pp = i + (pv1-1)*im
              RHS2(pv1,nv+1) = b(pp,3,j,ke_xy)
            end do

        end if

        ! LU factorization:
        call DGETRF( nv, nv, Amat1, nv, ipiv1, info )
        if ( info /= 0 ) then
          LOG_ERROR('linkernel_solve_sip',*)  'NU, DGETRF failed: info=', info, ', i=', i, ', j=', j, ', ke_xy=', ke_xy
          call PRC_abort
        end if
        call DGETRF( nv, nv, Amat2, nv, ipiv2, info )
        if ( info /= 0 ) then
          LOG_ERROR('linkernel_solve_sip',*)  'KH, DGETRF failed: info=', info, ', i=', i, ', j=', j, ', ke_xy=', ke_xy
          call PRC_abort
        end if

        ! Solve all right-hand sides using the same LU factors.
        call DGETRS( 'N', nv, nrhs1, Amat1, nv, ipiv1, RHS1, nv, info )
        if ( info /= 0 ) then
          LOG_ERROR('linkernel_solve_sip',*)  'NU, DGETRS failed: info=', info, ', i=', i, ', j=', j, ', ke_xy=', ke_xy
          call PRC_abort
        end if
        call DGETRS( 'N', nv, nrhs2, Amat2, nv, ipiv2, RHS2, nv, info )
        if ( info /= 0 ) then
          LOG_ERROR('linkernel_solve_sip',*)  'KH, DGETRS failed: info=', info, ', i=', i, ', j=', j, ', ke_xy=', ke_xy
          call PRC_abort
        end if

        if ( top_flag ) then
            ! At the top level RHS(:,1:2) contains D^{-1} b.
            do iv=1, 2
            do pv1=1, nv
              pp = i + (pv1-1)*im
              b(pp,iv,j,ke_xy) = RHS1(pv1,iv)
            end do
            end do
            do pv1=1, nv
              pp = i + (pv1-1)*im
              b(pp,3,j,ke_xy) = RHS2(pv1,1)
            end do

            ! Not used in backward substitution.
            do pv2=1, nv
            do pv1=1, nv
              pp = i + (pv1-1)*im
              G1(pp,pv2,j,ke_xy) = 0.0_RP
              G2(pp,pv2,j,ke_xy) = 0.0_RP
            end do
            end do

        else
            ! RHS(:,1:nv) now contains D^{-1} U.
            do pv2=1, nv
            do pv1=1, nv
              pp = i + (pv1-1)*im
              G1(pp,pv2,j,ke_xy) = RHS1(pv1,pv2)
              G2(pp,pv2,j,ke_xy) = RHS2(pv1,pv2)
            end do
            end do
            ! RHS(:,nv+1:nv+2) contains D^{-1} b.
            do iv=1, 2
            do pv1=1, nv
              pp = i + (pv1-1)*im
              b(pp,iv,j,ke_xy) = RHS1(pv1,nv+iv)
            end do
            end do
            do pv1=1, nv
              pp = i + (pv1-1)*im
              b(pp,3,j,ke_xy) = RHS2(pv1,nv+1)
            end do
        end if

        end do ! i
    end do ! j
    end do ! ke_xy
    !$omp end parallel do

    return
  end subroutine linkernel_solve_sip

!OCL SERIAL
  subroutine construct_matbnd_sip( BndMatL, BndMatD, BndMatU, & ! (out)
    RHO, KDIFF, GsqrtV, Dx1D, M1D,invM1D,                     & ! (in)
    impl_fac, dt, penalty_fac, lmesh, elem, im, jm, ke_z      ) ! (in)
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: im, jm
    real(RP), intent(out) :: BndMatL(im,elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP), intent(out) :: BndMatD(im,elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP), intent(out) :: BndMatU(im,elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP), intent(in) :: RHO(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: KDIFF(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: GsqrtV(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: Dx1D(elem%Nnode_v,elem%Nnode_v)
    real(RP), intent(in) :: M1D(elem%Nnode_v,elem%Nnode_v)
    real(RP), intent(in) :: invM1D(elem%Nnode_v,elem%Nnode_v)
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: penalty_fac    
    integer, intent(in) :: ke_z

    integer :: ke2D, ke, p
    integer :: i, j
    integer :: pv, pv1, pv2
    integer :: f
    integer :: ke_nb, ke_z_nb
    integer :: pvM, pvP

    real(RP) :: lambda

    real(RP) :: mu_loc (elem%Nnode_v)
    real(RP) :: rinv_loc(elem%Nnode_v)
    real(RP) :: rGsqrtV_loc(elem%Nnode_v)

    real(RP) :: mu_nb (elem%Nnode_v)
    real(RP) :: rinv_nb(elem%Nnode_v)

    real(RP) :: Avol (elem%Nnode_v,elem%Nnode_v)
    real(RP) :: AffMM(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: AffMP(elem%Nnode_v,elem%Nnode_v)

    real(RP) :: MinvALoc(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: MinvANb (elem%Nnode_v,elem%Nnode_v)

    real(RP) :: Dz(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: Dz_loc(elem%Nnode_v), Dz_nb(elem%Nnode_v)

    real(RP) :: mu_face_M, mu_face_P
    real(RP) :: sigma
    real(RP) :: Fscale_M

    logical :: boundary_flag

    integer :: Nnode_v
    !---------------------------------------------------------------------------

    call PROF_rapstart('phy_bl_cal_vi_matbnd_sip', 3)

    lambda = impl_fac
    Nnode_v = elem%Nnode_v

    !$omp parallel do collapse(2) private( ke, p, ke_nb, ke_z_nb, pvM, pvP, &
    !$omp mu_loc, rinv_loc, rGsqrtV_loc, mu_nb, rinv_nb,                    &
    !$omp Avol, AffMM, AffMP, MinvALoc, MinvANb, Dz, Dz_loc, Dz_nb,         &
    !$omp mu_face_M, mu_face_P, sigma, Fscale_M, boundary_flag )
    do ke2D = 1, lmesh%Ne2D
    do j = 1, jm
      ke = ke2D + (ke_z-1)*lmesh%Ne2D

      BndMatL(:,:,:,j,ke2D) = 0.0_RP
      BndMatD(:,:,:,j,ke2D) = 0.0_RP
      BndMatU(:,:,:,j,ke2D) = 0.0_RP

      do i=1, im
        do pv=1, Nnode_v
          p = i + (j-1)*im + (pv-1)*im*jm

          mu_loc(pv) = RHO(p,ke) * KDIFF(p,ke)
          rinv_loc(pv) = 1.0_RP / RHO(p,ke)
          rGsqrtV_loc(pv) = 1.0_RP / GsqrtV(p,ke)

          Avol(:,pv) = 0.0_RP
          BndMatD(i,pv,pv,j,ke2D) = 1.0_RP
        end do
        do pv2=1, Nnode_v
        do pv1=1, Nnode_v
          p = i + (j-1)*im + (pv1-1)*im*jm
          Dz(pv1,pv2) = lmesh%Escale(p,ke,3,3) * Dx1D(pv1,pv2)
        end do
        end do

        !- Volume contribution ---

        ! Avol = Dz^T M diag(mu) Dz diag(1/rho)
        do pv2=1, Nnode_v
        do pv1=1, Nnode_v
        do pv=1, Nnode_v
          Avol(pv1,pv2) = Avol(pv1,pv2) &
                        + Dz(pv,pv1) * M1D(pv,pv) * mu_loc(pv) * rGsqrtV_loc(pv) * Dz(pv,pv2) * rinv_loc(pv2)
        end do
        end do
        end do

        ! Minv Avol
        MinvALoc(:,:) = 0.0_RP
        do pv2=1, Nnode_v
          do pv=1, Nnode_v
          do pv1=1, Nnode_v
            MinvALoc(pv1,pv2) = MinvALoc(pv1,pv2) + invM1D(pv1,pv) * Avol(pv,pv2)
          end do
          end do
          MinvALoc(pv1,pv2) = rGsqrtV_loc(pv1) * MinvALoc(pv1,pv2)
        end do
        
        do pv2=1, Nnode_v
        do pv1=1, Nnode_v
          BndMatD(i,pv1,pv2,j,ke2D) = BndMatD(i,pv1,pv2,j,ke2D) + lambda * MinvALoc(pv1,pv2)
        end do
        end do

        !- Bottom and top faces

        do f=1, 2
          if (f == 1) then
            ke_z_nb = max(ke_z-1,1)
            pvM = 1; pvP = Nnode_v
            boundary_flag = (ke_z == 1)
          else
            ke_z_nb = min(ke_z+1,lmesh%NeZ)
            pvM = Nnode_v; pvP = 1
            boundary_flag = (ke_z == lmesh%NeZ)
          end if
          ke_nb = ke2D + (ke_z_nb-1)*lmesh%Ne2D

          ! Homogeneous Neumann boundary condition:
          !   mu dphi/dn = 0
          ! No SIP boundary contribution is added here.

          if (boundary_flag) cycle

          do pv=1, Nnode_v
            p = i + (j-1)*im + (pv-1)*im*jm
            mu_nb(pv) = RHO(p,ke_nb) * KDIFF(p,ke_nb)
            rinv_nb(pv) = 1.0_RP / RHO(p,ke_nb)
          end do

          mu_face_M = mu_loc(pvM)
          mu_face_P = mu_nb (pvP)
          sigma = penalty_fac * real(Nnode_v, kind=RP)**2 * max(mu_face_M, mu_face_P)


          !--

          Dz_loc(:) = Dz(pvM,:)
          p = i + (j-1)*im + (pvP-1)*im*jm
          do pv=1, Nnode_v
            Dz_nb(pv) = lmesh%Escale(p,ke_nb,3,3) * Dx1D(pvP,pv)
          end do

          call construct_sip_face_blocks_lgl( AffMM, AffMP,                & ! (out)
            Dz_loc, Dz_nb, mu_face_M, mu_face_P, rinv_loc, rinv_nb, & ! (in)
            sigma, pvM, pvP, f, elem%Nnode_v                               ) ! (in)

          Fscale_M = lmesh%Fscale(elem%Nfp_h*elem%Nfaces_h+1,ke)
          AffMM(:,:) = Fscale_M * AffMM(:,:)
          AffMP(:,:) = Fscale_M * AffMP(:,:)

          ! M^-1 AffMM, M^-1 AffMP
          MinvALoc(:,:) = 0.0_RP; MinvAnb(:,:) = 0.0_RP
          do pv2=1, Nnode_v
            do pv=1, Nnode_v
            do pv1=1, Nnode_v
                MinvALoc(pv1,pv2) = MinvALoc(pv1,pv2) + invM1D(pv1,pv) * AffMM(pv,pv2)
                MinvANb (pv1,pv2) = MinvANb (pv1,pv2) + invM1D(pv1,pv) * AffMP(pv,pv2)
            end do
            end do
          end do

          do pv2=1, Nnode_v
          do pv1=1, Nnode_v
            BndMatD(i,pv1,pv2,j,ke2D) = BndMatD(i,pv1,pv2,j,ke2D) + lambda * rGsqrtV_loc(pv1) * MinvALoc(pv1,pv2)
            if ( f == 1 ) then
              BndMatL(i,pv1,pv2,j,ke2D) = BndMatL(i,pv1,pv2,j,ke2D) + lambda * rGsqrtV_loc(pv1) * MinvANb(pv1,pv2)
            else
              BndMatU(i,pv1,pv2,j,ke2D) = BndMatU(i,pv1,pv2,j,ke2D) + lambda * rGsqrtV_loc(pv1) * MinvANb(pv1,pv2)
            end if
          end do
          end do

        end do ! end face loop
      end do
    end do    
    end do

    call PROF_rapend('phy_bl_cal_vi_matbnd_sip', 3)
    return
  end subroutine construct_matbnd_sip

!OCL SERIAL
  subroutine construct_sip_face_blocks_lgl( &
    AffMM, AffMP,                           & ! (out)
    dM, dP, muM, muP, rinvM, rinvP,         & ! (in)
    sigma, pvM, pvP, face_id, nv            ) ! (in)
    implicit none
    integer, intent(in) :: nv
    real(RP), intent(out) :: AffMM(nv,nv)
    real(RP), intent(out) :: AffMP(nv,nv)
    real(RP), intent(in) :: dM(nv), dP(nv)
    real(RP), intent(in) :: muM, muP
    real(RP), intent(in) :: rinvM(nv), rinvP(nv)
    real(RP), intent(in) :: sigma
    integer, intent(in) :: pvM, pvP
    integer, intent(in) :: face_id

    integer :: i, j
    real(RP) :: ncom
    !----------------------------------------------------

    if (face_id == 1) then
      ncom = -1.0_RP
    else
      ncom = +1.0_RP
    end if

    AffMM(:,:) = 0.0_RP; AffMP(:,:) = 0.0_RP
    do i=1, nv
    do j=1, nv
      if (i == pvM) then
        AffMM(i,j) = AffMM(i,j)                           &
                   - 0.5_RP * muM * ncom * dM(j) * rinvM(j)
      end if
      if (j == pvM) then
        AffMM(i,j) = AffMM(i,j)                           &
                   - 0.5_RP * muM * ncom * dM(i) * rinvM(j)
      end if
      if (i == pvM .and. j == pvM) then
        AffMM(i,j) = AffMM(i,j) + sigma * rinvM(j)
      end if
      !-  
      if ( i == pvM ) then
        AffMP(i,j) = AffMP(i,j)                           &
                   - 0.5_RP * muP * ncom * dP(j) * rinvP(j)
      end if
      if ( j == pvP ) then
        AffMP(i,j) = AffMP(i,j)                           &
                   + 0.5_RP * muM * ncom * dM(i) * rinvP(j)
      end if
      if ( i == pvM .and. j == pvP ) then
        AffMP(i,j) = AffMP(i,j) - sigma * rinvP(j)
      end if
    end do
    end do
    return
  end subroutine construct_sip_face_blocks_lgl
    
!OCL SERIAL
  subroutine eval_Ax( MOMX_t, MOMY_t, DRHOT_t, alph_M, alph_H, &
    PROG_VARS, MOMX00, MOMY00, PT00, NU, KH, DENS, GsqrtV, impl_fac, dt, &
    lmesh, elem, vmapM, vmapP, is_bound, element3D_operation, C_IP,  &
    im, jm, b, use_delta_form )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: MOMX_t(elem%Np,lmesh%Ne)
    real(RP), intent(out) :: MOMY_t(elem%Np,lmesh%Ne)
    real(RP), intent(out) :: DRHOT_t(elem%Np,lmesh%Ne)
    real(RP), intent(out) :: alph_M(elem%NfpTot,lmesh%Ne)
    real(RP), intent(out) :: alph_H(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: PROG_VARS(elem%Np,lmesh%NeX*lmesh%NeY*lmesh%NeZ,3)
    real(RP), intent(in) :: MOMX00(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY00(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PT00(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: NU(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: KH(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: GsqrtV(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)
    class(ElementOperationBase3D), intent(in) :: element3D_operation
    real(RP), intent(in) :: C_IP
    integer, intent(in) :: im, jm
    real(RP), intent(out), optional :: b(im,elem%Nnode_v,3,jm,lmesh%Ne)
    logical, intent(in), optional :: use_delta_form
    
    real(RP) :: Flux(elem%Np,3), DFlux(elem%Np,2,3)
    real(RP) :: RDENS(elem%Np)
    real(RP) :: RGsqrtV
    real(RP) :: E33

    real(RP) :: DIFF_flux_z_broken(elem%Np,lmesh%NeA,3)
    real(RP) :: DIFF_flux_z(elem%Np,lmesh%NeA,3)
    
    real(RP) :: del_flux(elem%NfpTot,3,lmesh%Ne)

    integer :: ke_xy, ke_z
    integer :: ke, ke2D
    integer :: p, fp
    integer :: iv

    integer :: i, j
    integer :: pv
    logical :: flag_cal_b
    logical :: flag_use_delta_form
    !---------------------------------------------------------

    if ( present(b) .and. present(use_delta_form) ) then
      flag_cal_b = .true.
      flag_use_delta_form = use_delta_form
    else
      flag_cal_b          = .false.
      flag_use_delta_form = .true.
    end if

    if ( flag_use_delta_form ) then
      call cal_grad_del_flux( del_flux,                       & ! (out)
        PROG_VARS, DENS,                                      & ! (in)
        lmesh%normal_fn(:,:,3), lmesh%Fscale, vmapM, vmapP,   & ! (in)
        lmesh, elem, lmesh%lcmesh2D, lmesh%lcmesh2D%refElem2D ) ! (in)

      !$omp parallel do collapse(2) private( ke, ke2D, DFlux, RDENS, RGsqrtV, E33 )
      do ke_z=1, lmesh%NeZ
      do ke_xy=1, lmesh%NeX*lmesh%NeY
        ke = ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
        ke2D = lmesh%EMap3Dto2D(ke)

        RDENS(:) = 1.0_RP / DENS(:,ke)
        do iv=1, 3
          call element3D_operation%Dz( PROG_VARS(:,ke,iv) * RDENS(:), DFlux(:,1,iv) )
          call element3D_operation%Lift( del_flux(:,iv,ke), DFlux(:,2,iv) )
        end do

        do p=1, elem%Np
          RGsqrtV = 1.0_RP / GsqrtV(p,ke)
          E33 = lmesh%Escale(p,ke,3,3)
          DIFF_flux_z_broken(p,ke,1) = DENS(p,ke) * NU(p,ke) * E33 * DFlux(p,1,1) * RGsqrtV
          DIFF_flux_z_broken(p,ke,2) = DENS(p,ke) * NU(p,ke) * E33 * DFlux(p,1,2) * RGsqrtV
          DIFF_flux_z_broken(p,ke,3) = DENS(p,ke) * KH(p,ke) * E33 * DFlux(p,1,3) * RGsqrtV

          DIFF_flux_z(p,ke,1) = DENS(p,ke) * NU(p,ke) * ( E33 * DFlux(p,1,1) + DFlux(p,2,1) ) * RGsqrtV
          DIFF_flux_z(p,ke,2) = DENS(p,ke) * NU(p,ke) * ( E33 * DFlux(p,1,2) + DFlux(p,2,2) ) * RGsqrtV
          DIFF_flux_z(p,ke,3) = DENS(p,ke) * KH(p,ke) * ( E33 * DFlux(p,1,3) + DFlux(p,2,3) ) * RGsqrtV
        end do
      end do
      end do

      !---------------------------------------------------------

      call cal_del_flux( del_flux, alph_M, alph_H,                      & ! (out)
        DIFF_flux_z, DIFF_flux_z_broken, PROG_VARS, DENS, NU, KH, C_IP, & ! (in)
        lmesh%normal_fn(:,:,3), lmesh%Fscale, vmapM, vmapP, is_bound,   & ! (in)
        lmesh, elem, lmesh%lcmesh2D, lmesh%lcmesh2D%refElem2D )           ! (in)

      !$omp parallel do collapse(2) private( ke, ke2D, DFlux, RGsqrtV, E33 )
      do ke_z=1, lmesh%NeZ
      do ke_xy=1, lmesh%NeX*lmesh%NeY
        ke = ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
        ke2D = lmesh%EMap3Dto2D(ke)
        do iv=1, 3
          call element3D_operation%Dz( DIFF_flux_z(:,ke,iv), DFlux(:,1,iv) )
          call element3D_operation%Lift( del_flux(:,iv,ke), DFlux(:,2,iv) )
        end do

        do p=1, elem%Np
          RGsqrtV = 1.0_RP / GsqrtV(p,ke)
          E33 = lmesh%Escale(p,ke,3,3)
          MOMX_t (p,ke) = ( E33 * DFlux(p,1,1) + DFlux(p,2,1) ) * RGsqrtV
          MOMY_t (p,ke) = ( E33 * DFlux(p,1,2) + DFlux(p,2,2) ) * RGsqrtV
          DRHOT_t(p,ke) = ( E33 * DFlux(p,1,3) + DFlux(p,2,3) ) * RGsqrtV
        end do
      end do
      end do
    end if

    if ( flag_cal_b ) then
      if ( use_delta_form ) then
        !$omp parallel do private(p)
        do ke=1, lmesh%Ne
        do pv=1, elem%Nnode_v
        do j=1, jm
        do i=1, im
            p  = i + (j-1)*im + (pv-1)*im*jm
            b(i,pv,1,j,ke) = impl_fac * MOMX_t(p,ke)  &
                        - PROG_VARS  (p,ke,1)      &
                        + MOMX00(p,ke)
            b(i,pv,2,j,ke) = impl_fac * MOMY_t(p,ke)  &
                        - PROG_VARS  (p,ke,2)      &
                        + MOMY00(p,ke)
            b(i,pv,3,j,ke) = impl_fac * DRHOT_t(p,ke) &
                        - PROG_VARS  (p,ke,3)      &
                        + DENS(p,ke) * PT00(p,ke)
        end do
        end do
        end do
        end do
      else
        !$omp parallel do private(p)
        do ke=1, lmesh%Ne
        do pv=1, elem%Nnode_v
        do j=1, jm
        do i=1, im
            p  = i + (j-1)*im + (pv-1)*im*jm
            b(i,pv,1,j,ke) = MOMX00(p,ke)
            b(i,pv,2,j,ke) = MOMY00(p,ke)
            b(i,pv,3,j,ke) = DENS(p,ke) * PT00(p,ke)
        end do
        end do
        end do
        end do
      end if
    end if
    return
  end subroutine eval_Ax

!OCL SERIAL
  subroutine cal_del_flux( del_flux, alph_M, alph_H, &
    DIFF_flux_z, DIFF_flux_z_broken, PROG_VARS, DENS, NU, KH, C_IP,  &
    nz, Fscale, vmapM, vmapP, is_bound, lmesh, elem, lmesh2D, elem2D   )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: del_flux(elem%NfpTot,3,lmesh%Ne)
    real(RP), intent(out) :: alph_M(elem%NfpTot,lmesh%Ne)
    real(RP), intent(out) :: alph_H(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: DIFF_flux_z(elem%Np*lmesh%NeA,3)
    real(RP), intent(in) :: DIFF_flux_z_broken(elem%Np*lmesh%NeA,3)
    real(RP), intent(in) :: PROG_VARS(elem%Np*lmesh%NeX*lmesh%NeY*lmesh%NeZ,3)
    real(RP), intent(in) :: DENS(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: NU(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: KH(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: C_IP
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: Fscale(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)

    integer :: ke, ke_z, ke_xy
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: DIFF_flux_z_P(elem%NfpTot,3)

    real(RP) :: coef(elem%NfpTot)
    real(RP) :: RDENS_M(elem%NfpTot), RDENS_P(elem%NfpTot)
    real(RP) :: numflux(elem%NfpTot,3)
    !------------------------------------------------

    !$omp parallel do collapse(2) &
    !$omp private( ke, iM, iP, alph_M, alph_H, coef, RDENS_M, RDENS_P, numflux )
    do ke_z=1, lmesh%NeZ
    do ke_xy=1, lmesh%NeX*lmesh%NeY
      ke = ke_xy + (ke_z-1)*lmesh%Ne2D
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      
      alph_M(:,ke) = C_IP * (elem%Nnode_v)**2 * max( DENS(iM)*NU(iM), DENS(iP)*NU(iP) )
      alph_H(:,ke) = C_IP * (elem%Nnode_v)**2 * max( DENS(iM)*KH(iM), DENS(iP)*KH(iP) )

      coef(:) = nz(:,ke) * Fscale(:,ke)
      RDENS_M(:) = 1.0_RP / DENS(iM)
      RDENS_P(:) = 1.0_RP / DENS(iP)

      where ( is_bound(:,ke) )
        numflux(:,1) = 0.0_RP
      elsewhere
        numflux(:,1) = 0.5_RP * ( DIFF_flux_z_broken(iP,1) + DIFF_flux_z(iM,1) )
        numflux(:,2) = 0.5_RP * ( DIFF_flux_z_broken(iP,2) + DIFF_flux_z(iM,2) )
        numflux(:,3) = 0.5_RP * ( DIFF_flux_z_broken(iP,3) + DIFF_flux_z(iM,3) )
      end where

      !-
      del_flux(:,1,ke) = coef(:) * ( numflux(:,1) - DIFF_flux_z(iM,1) ) &
                       + alph_M(:,ke) * Fscale(:,ke) * ( PROG_VARS(iP,1) * RDENS_P(:)- PROG_VARS(iM,1) * RDENS_M(:) )

      del_flux(:,2,ke) = coef(:) * ( numflux(:,2) - DIFF_flux_z(iM,2) ) &
                       + alph_M(:,ke) * Fscale(:,ke) * ( PROG_VARS(iP,2) * RDENS_P(:) - PROG_VARS(iM,2) * RDENS_M(:) )
    
      del_flux(:,3,ke) = coef(:) * ( numflux(:,3) - DIFF_flux_z(iM,3) ) &
                       + alph_H(:,ke) * Fscale(:,ke) * ( PROG_VARS(iP,3) * RDENS_P(:) - PROG_VARS(iM,3) * RDENS_M(:) )
    end do
    end do

    return
  end subroutine cal_del_flux

!OCL SERIAL
  subroutine cal_grad_del_flux( del_flux, &
    PROG_VARS, DENS, nz, Fscale,vmapM, vmapP, &
    lmesh, elem, lmesh2D, elem2D              )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: del_flux(elem%NfpTot,3,lmesh%Ne)
    real(RP), intent(in) :: PROG_VARS(elem%Np*lmesh%NeX*lmesh%NeY*lmesh%NeZ,3)
    real(RP), intent(in) :: DENS(elem%Np*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: Fscale(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)

    integer :: ke, ke_z, ke_xy
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)

    real(RP) :: coef(elem%NfpTot)
    real(RP) :: RDENS_M(elem%NfpTot), RDENS_P(elem%NfpTot)
    !------------------------------------------------

    !$omp parallel do collapse(2) &
    !$omp private( ke, iM, iP, coef, RDENS_M, RDENS_P )
    do ke_z=1, lmesh%NeZ
    do ke_xy=1, lmesh%NeX*lmesh%NeY
      ke = ke_xy + (ke_z-1)*lmesh%Ne2D
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)

      coef(:) = 0.5_RP * nz(:,ke) * Fscale(:,ke)
      RDENS_M(:) = 1.0_RP / DENS(iM)
      RDENS_P(:) = 1.0_RP / DENS(iP)

      del_flux(:,1,ke) = coef(:) * ( PROG_VARS(iP,1) * RDENS_P(:) - PROG_VARS(iM,1) * RDENS_M(:) )
      del_flux(:,2,ke) = coef(:) * ( PROG_VARS(iP,2) * RDENS_P(:) - PROG_VARS(iM,2) * RDENS_M(:) )
      del_flux(:,3,ke) = coef(:) * ( PROG_VARS(iP,3) * RDENS_P(:) - PROG_VARS(iM,3) * RDENS_M(:) )
    end do
    end do

    return
  end subroutine cal_grad_del_flux
end module scale_atm_phy_bl_dgm_common