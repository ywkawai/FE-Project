!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Nonhydrostatic model / HEVI / Common
!!
!! @par Description
!!      HEVI DGM scheme for Atmospheric dynamical process
!!      As the thermodynamics equation, a density * potential temperature equation is used.
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_rhot_hevi_common_2
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

  use scale_element_operation_base, only: ElementOperationBase3D

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    DENS_VID => PRGVAR_DDENS_ID, RHOT_VID => PRGVAR_DRHOT_ID, &
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID,                               &
    PRGVAR_NUM
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_Init
  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_Final
  
  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_eval_Ax
  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_solve
  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd

  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_eval_Ax_uv
  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_solve_uv
  public :: atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd_uv

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: VI_use_lapack_flag = .true.
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  private :: vi_cal_del_flux_dyn
  private :: vi_cal_del_flux_dyn_uv

contains
  !------------------------------------------------
  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_Init
    implicit none
    namelist / PARAM_ATMOS_DYN_NONHYDRO3D_RHOT_HEVI_COMMON / &
      VI_use_lapack_flag

    integer :: ierr
    !---------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_NONHYDRO3D_RHOT_HEVI_COMMON,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_nonhydro3d_rhot_hevi_common_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_nonhydro3d_rhot_hevi_common_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_NONHYDRO3D_RHOT_HEVI_COMMON. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN_NONHYDRO3D_RHOT_HEVI_COMMON)

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_Init

  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_Final
    implicit none
    !---------------------------------------------------------
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_Final

!OCL SERIAL  
  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_eval_Ax( &
    DENS_t, MOMZ_t, RHOT_t,                                  & ! (out)
    alph,                                                    & ! (in)
    PROG_VARS, PROG_VARS0,                                   & ! (in)
    DDENS00, MOMX00, MOMY00, MOMZ00, DRHOT00,                & ! (in)
    DENS_hyd, PRES_hyd,                                      & ! (in)
    Rtot, CPtot, CVtot, GsqrtV,                              & ! (in)
    impl_fac, dt,                                            & ! (in)
    lmesh, elem,                                             & ! (in)
    vmapM, vmapP,                                            & ! (in)
    element3D_operation,                                     & ! (in)
    im, jm,                                                  & ! (in)
    b                                                        ) ! (out)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: RHOT_t(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: alph(elem%NfpTot,lmesh%NeX*lmesh%NeY,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS  (elem%Np,lmesh%NeX*lmesh%NeY*lmesh%NeZ,PRGVAR_NUM)
    real(RP), intent(in)  :: PROG_VARS0 (elem%Np,lmesh%NeX*lmesh%NeY*lmesh%NeZ,PRGVAR_NUM)
    real(RP), intent(in)  :: DDENS00(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX00 (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY00 (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ00 (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT00(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: GsqrtV(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)    
    class(ElementOperationBase3D), intent(in) :: element3D_operation
    integer, intent(in) :: im, jm
    real(RP), intent(out), optional :: b(im,3,elem%Nnode_v,jm,lmesh%Ne)

    real(RP) :: Flux(elem%Np,3), DFlux(elem%Np,2,3)
    real(RP) :: del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP) :: DPRES(elem%Np), RHOT(elem%Np)

    real(RP) :: rdens_, u_, v_, w_, pt_
    real(RP) :: RGsqrtV(elem%Np), RGsqrt(elem%Np)
    real(RP) :: Gsqrt_, GsqrtDPRES_, E33

    integer :: ke_xy, ke_z
    integer :: ke, ke2d
    real(RP) :: rdt
    real(RP) :: drho(elem%Np)
    
    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR

    integer :: i, j, p, pp, pv
    !--------------------------------------------------------

    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    rdt = 1.0_RP / dt

    call vi_cal_del_flux_dyn( del_flux, & ! (out)
      alph, PROG_VARS, PROG_VARS0,                                              & ! (in)
      DENS_hyd, PRES_hyd, Rtot, CPtot, CVtot,                                   & ! (in)
      lmesh%Gsqrt, lmesh%GsqrtH, lmesh%gam, GsqrtV, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2), & ! (in)    
      lmesh%normal_fn(:,:,3), vmapM, vmapP, elem%IndexH2Dto3D_bnd,              & ! (in)
      lmesh, elem, lmesh%lcmesh2D, lmesh%lcmesh2D%refElem2D                     ) ! (in)

    !$omp parallel private( ke_xy, ke_z, ke, ke2d, p, &
    !$omp DPRES, RHOT, drho, pt_, Flux, DFlux,        &
    !$omp E33, RGsqrtV                                )

    !$omp do collapse(2)
    do ke_z=1, lmesh%NeZ
    do ke_xy=1, lmesh%NeX*lmesh%NeY
      ke = Ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
      ke2d = lmesh%EMap3Dto2D(ke)

      do p=1, elem%Np
        RGsqrtV(p) = 1.0_RP / GsqrtV(p,ke)
      end do

      !-
      do p=1, elem%Np
        Flux(p,DENS_VID) = PROG_VARS(p,ke,MOMZ_VID) &
                         + GsqrtV(p,ke) * lmesh%GI3(p,ke,1) * PROG_VARS(p,ke,MOMX_VID) &
                         + GsqrtV(p,ke) * lmesh%GI3(p,ke,2) * PROG_VARS(p,ke,MOMY_VID)
      end do

      RHOT(:) = P0ovR * (PRES_hyd(:,ke)/PRES00)**rgamm &
              + PROG_VARS(:,ke,RHOT_VID)
      do p=1, elem%Np
        pt_ = RHOT(p) / ( PROG_VARS(p,ke,DENS_VID) + DENS_hyd(p,ke) )
        Flux(p,RHOT_VID) = pt_ * Flux(p,DENS_VID)
      end do

      do p=1, elem%Np
        Flux(p,MOMZ_VID) = PRES00 * ( Rtot(p,ke) * rP0 * RHOT(p) )**(CPtot(p,ke) / CVtot(p,ke)) &
                         - PRES_hyd(p,ke)
      end do

      !-
      call element3D_operation%Dz( Flux(:,DENS_VID), DFlux(:,1,DENS_VID) )
      call element3D_operation%Lift( del_flux(:,DENS_VID,ke), DFlux(:,2,DENS_VID)  )

      call element3D_operation%Dz( Flux(:,RHOT_VID), DFlux(:,1,RHOT_VID) )
      call element3D_operation%Lift( del_flux(:,RHOT_VID,ke), DFlux(:,2,RHOT_VID)  )

      call element3D_operation%Dz( Flux(:,MOMZ_VID), DFlux(:,1,MOMZ_VID) )
      call element3D_operation%Lift( del_flux(:,MOMZ_VID,ke), DFlux(:,2,MOMZ_VID)  )

      !-
      do p=1, elem%Np
        E33 = lmesh%Escale(p,ke,3,3)

        DENS_t(p,ke) = - ( &
              E33 * DFlux(p,1,DENS_VID)    &
            + DFlux(p,2,DENS_VID) ) * RGsqrtV(p)
          
        RHOT_t(p,ke) = - ( &
              E33 * DFlux(p,1,RHOT_VID)    &
            + DFlux(p,2,RHOT_VID) ) * RGsqrtV(p)
      end do

      call element3D_operation%VFilterPM1( PROG_VARS(:,ke,DENS_VID), &
        drho )
      do p=1, elem%Np
        E33 = lmesh%Escale(p,ke,3,3)

        MOMZ_t(p,ke) = - ( &
              E33 * DFlux(p,1,MOMZ_VID) &
            + DFlux(p,2,MOMZ_VID) ) * RGsqrtV(p) &
            - Grav * drho(p) 
      end do
    end do
    end do
    !$omp end do

    if ( present( b ) ) then
      !$omp do collapse(2)
      do ke=1, lmesh%Ne
      do pv=1, elem%Nnode_v
      do j=1, jm
      do i=1, im
        p  = i + (j-1)*im + (pv-1)*im*jm
        b(i,1,pv,j,ke) = impl_fac * DENS_t(p,ke)         &
                       - PROG_VARS  (p,ke,DENS_VID)      &
                       + DDENS00(p,ke)
        b(i,2,pv,j,ke) = impl_fac * MOMZ_t(p,ke)         &
                       - PROG_VARS  (p,ke,MOMZ_VID)      &
                       + MOMZ00(p,ke)
        b(i,3,pv,j,ke) = impl_fac * RHOT_t(p,ke)         &
                       - PROG_VARS  (p,ke,RHOT_VID)      &
                       + DRHOT00(p,ke)
      end do
      end do
      end do
      end do
      !$omp end do
    end if
    !$omp end parallel

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_eval_Ax

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_solve( PROG_VARS, &
    b, PROG_VARS0, Rtot, CPtot, CVtot, alph, DENS_hyd, PRES_hyd, &
    GsqrtV, IntrpMat_VPOrdM1,                                    &
    impl_fac, dt,                                                &
    vmapM, vmapP, lmesh, elem, im, jm                            )
    use scale_atm_dyn_dgm_hevi_common_linalgebra, only: &
      linkernel_ludecomp => atm_dyn_dgm_hevi_common_linalgebra_ludecomp, &
      linkernel_solve => atm_dyn_dgm_hevi_common_linalgebra_solve
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: im
    integer, intent(in) :: jm
    real(RP), intent(inout) :: PROG_VARS (elem%Nnode_h1D**2,elem%Nnode_v,lmesh%NeX*lmesh%NeY,lmesh%NeZ,PRGVAR_NUM)
    real(RP), intent(inout) :: b(im*3*elem%Nnode_v,jm,lmesh%Ne2D,lmesh%NeZ)
    real(RP), intent(in) :: PROG_VARS0(elem%Np,lmesh%NeX*lmesh%NeY,lmesh%NeZ,PRGVAR_NUM)
    real(RP), intent(in) :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: alph(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: GsqrtV(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)

    real(RP) :: tmp_v3(3)
    real(RP) :: BndMatD(im*3*elem%Nnode_v,3*elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP) :: BndMatL(im*3*elem%Nnode_v,3,jm,lmesh%Ne2D)
    real(RP) :: BndMatU(im*3*elem%Nnode_v,3,jm,lmesh%Ne2D)
    real(RP) :: D(im*3*elem%Nnode_v,3*elem%Nnode_v)
    real(RP) :: G(im*3*elem%Nnode_v,3,jm,lmesh%Ne2D,lmesh%NeZ)
    integer :: ipiv(im*3*elem%Nnode_v)

    integer :: ke_xy, ke_z
    integer :: i, pv, j
    integer :: v
    integer :: pp
    integer :: p1, p2, p3
    !-----------------------------------------------------------------

    ! ke_z=1
    call atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd( &
      BndMatL, BndMatD, BndMatU,                                          & ! (out)
      PROG_VARS0, DENS_hyd, PRES_hyd, Rtot, CPtot, CVtot, alph,           & ! (in)
      lmesh%GI3(:,:,1), lmesh%GI3(:,:,2), GsqrtV, lmesh%normal_fn(:,:,3), & ! (in)
      IntrpMat_VPOrdM1,                                                   & ! (in)
      impl_fac, dt, lmesh, elem, im, jm, vmapM, vmapP,                    & ! (in)
      1 ) ! (in)
    
    !$omp parallel do private( ke_xy,j,v, D, ipiv ) collapse(2)
    do ke_xy=1, lmesh%NeX * lmesh%NeY
    do j=1, jm
      D(:,:) = BndMatD(:,:,j,ke_xy)

      call linkernel_ludecomp( D, ipiv, im, elem%Nnode_v*3 )
      call linkernel_solve( D, b(:,j,ke_xy,1), ipiv, im, elem%Nnode_v*3 )
      if ( lmesh%NeZ > 1) then
        do v=1, 3
          G(:,v,j,ke_xy,1) = BndMatU(:,v,j,ke_xy)
          call linkernel_solve( D, G(:,v,j,ke_xy,1), ipiv, im, elem%Nnode_v*3 )
        end do
      end if
    end do
    end do

    ! ke_z=2..NeZ

    do ke_z=2, lmesh%NeZ
      call atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd( &
        BndMatL, BndMatD, BndMatU,                                          & ! (out)
        PROG_VARS0, DENS_hyd, PRES_hyd, Rtot, CPtot, CVtot, alph,           & ! (in)
        lmesh%GI3(:,:,1), lmesh%GI3(:,:,2), GsqrtV, lmesh%normal_fn(:,:,3), & ! (in)
        IntrpMat_VPOrdM1,                                                   & ! (in)
        impl_fac, dt, lmesh, elem, im, jm, vmapM, vmapP,                    & ! (in)
        ke_z ) ! (in)
    
      !$omp parallel do private( ke_xy,j,pv,i,v, pp,p1,p2,p3,tmp_v3, D, ipiv ) collapse(2)
      do ke_xy=1, lmesh%NeX * lmesh%NeY
      do j=1, jm
        D(:,:) = BndMatD(:,:,j,ke_xy)

        do pv=1, 3*elem%Nnode_v
        do i=1, im
          p1 = i        + (elem%Nnode_v-1)*3*im
          p2 = i +   im + (elem%Nnode_v-1)*3*im
          p3 = i + 2*im + (elem%Nnode_v-1)*3*im
          pp = i + (pv-1)*im

          tmp_v3(:) = BndMatL(pp,:,j,ke_xy)
          D(pp,1) = D(pp,1) - tmp_v3(1) * G(p1,1,j,ke_xy,ke_z-1) &
                            - tmp_v3(2) * G(p2,1,j,ke_xy,ke_z-1) &
                            - tmp_v3(3) * G(p3,1,j,ke_xy,ke_z-1)
          D(pp,2) = D(pp,2) - tmp_v3(1) * G(p1,2,j,ke_xy,ke_z-1) &
                            - tmp_v3(2) * G(p2,2,j,ke_xy,ke_z-1) &
                            - tmp_v3(3) * G(p3,2,j,ke_xy,ke_z-1)
          D(pp,3) = D(pp,3) - tmp_v3(1) * G(p1,3,j,ke_xy,ke_z-1) &
                            - tmp_v3(2) * G(p2,3,j,ke_xy,ke_z-1) &
                            - tmp_v3(3) * G(p3,3,j,ke_xy,ke_z-1)
          b(pp,j,ke_xy,ke_z) = b(pp,j,ke_xy,ke_z) &
            - tmp_v3(1) * b(p1,j,ke_xy,ke_z-1) &
            - tmp_v3(2) * b(p2,j,ke_xy,ke_z-1) &
            - tmp_v3(3) * b(p3,j,ke_xy,ke_z-1)              
        end do
        end do
        call linkernel_ludecomp( D, ipiv, im, elem%Nnode_v*3 )
        call linkernel_solve( D, b(:,j,ke_xy,ke_z), ipiv, im, 3*elem%Nnode_v )

        if ( ke_z < lmesh%NeZ ) then
          do v=1, 3
            G(:,:,j,ke_xy,ke_z) = BndMatU(:,:,j,ke_xy)
            call linkernel_solve( D, G(:,v,j,ke_xy,ke_z), ipiv, im, 3*elem%Nnode_v )
          end do              
        end if
      end do  
      end do
    end do
    do ke_z=lmesh%NeZ-1, 1
      !$omp parallel do private( ke_xy,j,pv,i, pp,p1,p2,p3, tmp_v3 ) collapse(2)          
      do ke_xy=1, lmesh%NeX * lmesh%NeY
      do j=1, jm
        do pv=1, elem%Nnode_v*3
        do i=1, im
          p1 = i        ! + (1-1)*3*im
          p2 = i +   im ! + (1-1)*3*im
          p3 = i + 2*im ! + (1-1)*3*im
          pp = i + (pv-1)*im

          tmp_v3(:) = G(pp,:,j,ke_xy,ke_z)
          b(pp,j,ke_xy,ke_z) = b(pp,j,ke_xy,ke_z) &
            - tmp_v3(1) * b(p1,j,ke_xy,ke_z+1) &
            - tmp_v3(2) * b(p2,j,ke_xy,ke_z+1) &
            - tmp_v3(3) * b(p3,j,ke_xy,ke_z+1)
        end do
        end do            
        do pv=1, elem%Nnode_v
        do i=1, im
          p1 = i        + (pv-1)*3*im
          p2 = i +   im + (pv-1)*3*im
          p3 = i + 2*im + (pv-1)*3*im
          pp = i + (j-1)*im

          PROG_VARS(pp,pv,ke_xy,ke_z,DENS_VID) = PROG_VARS(pp,pv,ke_xy,ke_z,DENS_VID) + b(p1,j,ke_xy,ke_z)
          PROG_VARS(pp,pv,ke_xy,ke_z,MOMZ_VID) = PROG_VARS(pp,pv,ke_xy,ke_z,MOMZ_VID) + b(p2,j,ke_xy,ke_z)
          PROG_VARS(pp,pv,ke_xy,ke_z,RHOT_VID) = PROG_VARS(pp,pv,ke_xy,ke_z,RHOT_VID) + b(p3,j,ke_xy,ke_z)
        end do
        end do
      end do
    end do
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_solve

!OCL SERIAL  
  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_eval_Ax_uv( &
    MOMX_t, MOMY_t,                                          & ! (out)
    alph,                                                    & ! (out)
    PROG_VARS, PROG_VARS0,                                   & ! (in)
    DDENS00, MOMX00, MOMY00, MOMZ00, DRHOT00,                & ! (in)
    DENS_hyd, PRES_hyd,                                      & ! (in)
    Rtot, CPtot, CVtot, GsqrtV,                              & ! (in)
    impl_fac, dt,                                            & ! (in)
    lmesh, elem,                                             & ! (in)
    vmapM, vmapP,                                            & ! (in)
    element3D_operation,                                     & ! (in)
    im, jm,                                                  & ! (in)
    b1D_ij_uv                                                ) ! (out)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: MOMX_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: alph(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in)  :: PROG_VARS  (elem%Np,lmesh%Ne,PRGVAR_NUM)
    real(RP), intent(in)  :: PROG_VARS0 (elem%Np,lmesh%Ne,PRGVAR_NUM)
    real(RP), intent(in)  :: DDENS00(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX00 (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMY00 (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ00 (elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT00(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: GsqrtV(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)   
    class(ElementOperationBase3D), intent(in) :: element3D_operation
    integer, intent(in) :: im, jm
    real(RP), intent(out), optional :: b1D_ij_uv(im*elem%Nnode_v,2,jm,lmesh%Ne)

    real(RP) :: RGsqrtV
    real(RP) :: LiftDelFlx(elem%Np,2)
    real(RP) :: del_flux(elem%NfpTot,2,lmesh%Ne)
    integer :: ke_xy, ke_z
    integer :: ke, ke2D
    integer :: fp

    real(RP) :: rdt

    integer :: i, j, p, pp, pv
    !--------------------------------------------------------

    rdt = 1.0_RP / dt

    call vi_cal_del_flux_dyn_uv( del_flux, alph, & ! (out)
      PROG_VARS, PROG_VARS0,                                                   & ! (in)
      DENS_hyd, PRES_hyd, Rtot, CPtot, CVtot,                                  & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%GsqrtH, lmesh%gam, GsqrtV, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),     & ! (in)    
      lmesh%normal_fn(:,:,3), vmapM, vmapP, elem%IndexH2Dto3D_bnd,             & ! (in)
      lmesh, elem, lmesh%lcmesh2D, lmesh%lcmesh2D%refElem2D                    ) ! (in)

    !$omp parallel private( ke_xy, ke_z, ke, ke2D, fp, pv, i, j, p, pp, &
    !$omp LiftDelFlx, RGsqrtV                         )

    !$omp do collapse(2)
    do ke_z=1, lmesh%NeZ
    do ke_xy=1, lmesh%NeX*lmesh%NeY
      ke = Ke_xy + (ke_z-1)*lmesh%NeX*lmesh%NeY
      ke2D = lmesh%EMap3Dto2D(ke)

      call element3D_operation%Lift( del_flux(:,1,ke), LiftDelFlx(:,1)  )
      call element3D_operation%Lift( del_flux(:,2,ke), LiftDelFlx(:,2)  )

      do fp=1, elem%Np
        RGsqrtV = 1.0_RP / GsqrtV(fp,ke)
        MOMX_t(fp,ke) = - LiftDelFlx(fp,1) * RGsqrtV
        MOMY_t(fp,ke) = - LiftDelFlx(fp,2) * RGsqrtV
      end do
    end do
    end do
    !$omp end do

    if ( present( b1D_ij_uv ) ) then
      !$omp do collapse(2)
      do ke=1, lmesh%Ne
      do pv=1, elem%Nnode_v
      do j=1, jm
      do i=1, im
        pp = i + (pv-1)*im
        p  = i + (j-1)*im + (pv-1)*im*jm
        b1D_ij_uv(pp,1,j,ke) = impl_fac * MOMX_t(p,ke)    &
                             - PROG_VARS  (p,ke,MOMX_VID) &
                             + MOMX00(p,ke)
        b1D_ij_uv(pp,2,j,ke) = impl_fac * MOMY_t(p,ke)    &
                             - PROG_VARS  (p,ke,MOMY_VID) &
                             + MOMY00(p,ke)
      end do
      end do
      end do
      end do
      !$omp end do
    end if
    !$omp end parallel

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_eval_Ax_uv

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_solve_uv( PROG_VARS, &
    b_uv, alph,                                                  &
    GsqrtV,                                                      &
    impl_fac, dt,                                                &
    vmapM, vmapP, lmesh, elem, im, jm                            )
    use scale_atm_dyn_dgm_hevi_common_linalgebra, only: &
      linkernel_ludecomp => atm_dyn_dgm_hevi_common_linalgebra_ludecomp, &
      linkernel_solve => atm_dyn_dgm_hevi_common_linalgebra_solve
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: im
    integer, intent(in) :: jm
    real(RP), intent(inout) :: PROG_VARS(elem%Nnode_h1D**2,elem%Nnode_v,lmesh%NeX*lmesh%NeY,lmesh%NeZ,PRGVAR_NUM)
    real(RP), intent(inout) :: b_uv(im*elem%Nnode_v,2,jm,lmesh%Ne2D,lmesh%NeZ)
    real(RP), intent(in) :: alph(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: GsqrtV(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)

    real(RP) :: tmp
    real(RP) :: BndMatD_uv(im*elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP) :: BndMatL_uv(im*elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP) :: BndMatU_uv(im*elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP) :: D_uv(im*elem%Nnode_v,elem%Nnode_v)
    real(RP) :: G_uv(im*elem%Nnode_v,jm,lmesh%Ne2D,lmesh%NeZ)
    integer :: ipiv_uv(im*elem%Nnode_v)

    integer :: ke_xy, ke_z
    integer :: i, pv, j
    integer :: v
    integer :: p, pp
    integer :: p1, p2, p3
    !-----------------------------------------------------------------

    ! ke_z=1
    call atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd_uv( &
      BndMatL_uv, BndMatD_uv, BndMatU_uv,                            & ! (out)
      alph, GsqrtV, impl_fac, dt, lmesh, elem, im, jm, vmapM, vmapP, & ! (in)
      1 ) ! (in)
        
    !$omp parallel do private( ke_xy,j, D_uv, ipiv_uv ) collapse(2)
    do ke_xy=1, lmesh%NeX * lmesh%NeY
    do j=1, jm
      D_uv(:,:) = BndMatD_uv(:,:,j,ke_xy)
      call linkernel_ludecomp( D_uv, ipiv_uv, im, elem%Nnode_v )
      call linkernel_solve( D_uv, b_uv(:,1,j,ke_xy,1), ipiv_uv, im, elem%Nnode_v )
      call linkernel_solve( D_uv, b_uv(:,2,j,ke_xy,1), ipiv_uv, im, elem%Nnode_v )
      if ( lmesh%NeZ > 1) then
        G_uv(:,j,ke_xy,1) = BndMatU_uv(:,j,ke_xy)
        call linkernel_solve( D_uv, G_uv(:,j,ke_xy,1), ipiv_uv, im, elem%Nnode_v )
      end if
    end do
    end do
    ! ke_z=2..NeZ
    do ke_z=2, lmesh%NeZ
      call atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd_uv( &
        BndMatL_uv, BndMatD_uv, BndMatU_uv,                            & ! (out)
        alph, GsqrtV, impl_fac, dt, lmesh, elem, im, jm, vmapM, vmapP, & ! (in)
        ke_z ) ! (in)
      !$omp parallel do private( ke_xy,j,pv,i, pp,p2,tmp, D_uv, ipiv_uv ) collapse(2)
      do ke_xy=1, lmesh%NeX * lmesh%NeY
      do j=1, jm
        D_uv(:,:) = BndMatD_uv(:,:,j,ke_xy)
        do pv=1, elem%Nnode_v
        do i=1, im
          pp = i + (pv-1)*im
          p2 = i + (elem%Nnode_v-1)*im
          tmp = BndMatL_uv(pp,j,ke_xy)
          D_uv(pp,1)              = D_uv(pp,1)              - tmp * G_uv(p2,j,ke_xy,ke_z-1)
          b_uv(pp,1,j,ke_xy,ke_z) = b_uv(pp,1,j,ke_xy,ke_z) - tmp * b_uv(p2,1,j,ke_xy,ke_z-1)
          b_uv(pp,2,j,ke_xy,ke_z) = b_uv(pp,2,j,ke_xy,ke_z) - tmp * b_uv(p2,2,j,ke_xy,ke_z-1)
        end do
        end do
        call linkernel_ludecomp( D_uv, ipiv_uv, im, elem%Nnode_v )
        call linkernel_solve( D_uv, b_uv(:,1,j,ke_xy,ke_z), ipiv_uv, im, elem%Nnode_v )
        call linkernel_solve( D_uv, b_uv(:,2,j,ke_xy,ke_z), ipiv_uv, im, elem%Nnode_v )

        if ( ke_z < lmesh%NeZ ) then
          G_uv(:,j,ke_xy,ke_z) = BndMatU_uv(:,j,ke_xy)
          call linkernel_solve( D_uv, G_uv(:,j,ke_xy,ke_z), ipiv_uv, im, elem%Nnode_v )
        end if
      end do  
      end do
    end do
    do ke_z=lmesh%NeZ-1, 1
      !$omp parallel do private( ke_xy,j,pv,i, pp,p2,tmp, p ) collapse(2)
      do ke_xy=1, lmesh%NeX * lmesh%NeY
      do j=1, jm
        do pv=1, elem%Nnode_v
        do i=1, im
          pp = i + (pv-1)*im
          p2 = i ! + (1-1)*im        
          tmp = G_uv(pp,j,ke_xy,ke_z)
          b_uv(pp,1,j,ke_xy,ke_z) = b_uv(pp,1,j,ke_xy,ke_z) - tmp * b_uv(p2,1,j,ke_xy,ke_z+1)
          b_uv(pp,2,j,ke_xy,ke_z) = b_uv(pp,2,j,ke_xy,ke_z) - tmp * b_uv(p2,2,j,ke_xy,ke_z+1)

          p  = i + (j-1)*im
          PROG_VARS(p,pv,ke_xy,ke_z,MOMX_VID) = PROG_VARS(p,pv,ke_xy,ke_z,MOMX_VID) + b_uv(pp,1,j,ke_xy,ke_z)
          PROG_VARS(p,pv,ke_xy,ke_z,MOMY_VID) = PROG_VARS(p,pv,ke_xy,ke_z,MOMY_VID) + b_uv(pp,2,j,ke_xy,ke_z)
        end do
        end do
      end do
    end do
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_solve_uv

!--
!OCL SERIAL  
  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd( &
    BndMatL, BndMatD, BndMatU,              & ! (out)
    PROG_VARS0, DENS_hyd, PRES_hyd,         & ! (in)
    Rtot, CPtot, CVtot,                     & ! (in)
    alph, G13, G23, GsqrtV, nz,             & ! (in)
    IntrpMat_VPOrdM1,                       & ! (in)
    impl_fac, dt,                           & ! (in)
    lmesh, elem, im, jm,                    & ! (in)
    vmapM, vmapP, ke_z                      ) ! (in)

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: im, jm
    real(RP), intent(out) :: BndMatL(im,3,elem%Nnode_v,3,jm,lmesh%Ne2D)
    real(RP), intent(out) :: BndMatD(im,3,elem%Nnode_v,3,elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP), intent(out) :: BndMatU(im,3,elem%Nnode_v,3,jm,lmesh%Ne2D)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,lmesh%Ne,PRGVAR_NUM)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: alph(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  G13(elem%Np,lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: GsqrtV(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)    
    integer, intent(in) :: ke_z

    integer :: Colmask(elem%Nnode_v)
    real(RP) :: Id(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: Dd(elem%Nnode_v)
    real(RP) :: fac_dz(im,elem%Nnode_v,elem%Nnode_v,jm)
    real(RP) :: tmp1(im,elem%Nnode_v)

    real(RP) :: RHOT_hyd
    real(RP) :: DENS0(im,elem%Nnode_v,jm)
    real(RP) :: POT0(im,elem%Nnode_v,jm)
    real(RP) :: W0(im,elem%Nnode_v,jm)
    real(RP) :: WT0(im,elem%Nnode_v,jm)
    real(RP) :: DPDRHOT0(im,elem%Nnode_v,jm)
    real(RP) :: CPtot_ov_CVtot

    integer :: p, p2, pv, pv1, pv2
    integer :: fp,fp2
    integer :: ke, ke2D
    integer :: ke_z2, ke_2
    integer :: ij, i, j
    integer :: f1, f2
    integer :: v

    real(RP) :: fac
    logical :: bc_flag

    integer, parameter :: DENS_VID_LC = 1
    integer, parameter :: MOMZ_VID_LC = 2
    integer, parameter :: RHOT_VID_LC = 3

    real(RP) :: gamm, rgamm
    !---------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    Id(:,:) = 0.0_RP
    do p=1, elem%Nnode_v
      Id(p,p) = 1.0_RP
    end do

    do ke2D=1, lmesh%Ne2D
      ke = ke2D + (ke_z-1)*lmesh%Ne2D

      do j=1, jm
      do pv2=1, elem%Nnode_v
      do pv=1, elem%Nnode_v
      do i=1, im
        ij = i + (j-1)*im
        p = ij + (pv-1)*elem%Nnode_h1D**2
        fac_dz(i,pv,pv2,j) = impl_fac * lmesh%Escale(p,ke,3,3) / GsqrtV(p,ke) &
                           * elem%Dx3(p,elem%Colmask(pv2,ij))
      end do
      end do
      end do
      end do
      do j=1, jm
      do pv=1, elem%Nnode_v
      do i=1, im
        p = i + (j-1)*im + (pv-1)*im*jm
        RHOT_hyd = PRES00/Rdry * ( PRES_hyd(p,ke)/PRES00 )**rgamm
        DENS0(i,pv,j) = DENS_hyd(p,ke) + PROG_VARS0(p,ke,DENS_VID)
        POT0(i,pv,j) = ( RHOT_hyd + PROG_VARS0(p,ke,RHOT_VID) ) / DENS0(i,pv,j)
        W0(i,pv,j) = PROG_VARS0(p,ke,MOMZ_VID) / DENS0(i,pv,j)
      end do
      end do
      end do
      do j=1, jm
      do pv=1, elem%Nnode_v
      do i=1, im
        p = i + (j-1)*im + (pv-1)*im*jm
        WT0(i,pv,j) = W0(i,pv,j) + GsqrtV(p,ke) * ( &
              G13(p,ke) * PROG_VARS0(p,ke,MOMX_VID) &
            + G23(p,ke) * PROG_VARS0(p,ke,MOMY_VID) ) / DENS0(i,pv,j)
      end do
      end do
      end do
      do j=1, jm
      do pv=1, elem%Nnode_v
      do i=1, im
        p = i + (j-1)*im + (pv-1)*im*jm
        CPtot_ov_CVtot = CPtot(p,ke) / CVtot(p,ke)
        DPDRHOT0(i,pv,j) = CPtot_ov_CVtot &
                         * PRES00 * ( Rtot(p,ke) / PRES00 * ( DENS0(i,pv,j) * POT0(i,pv,j) ) )**CPtot_ov_CVtot &
                         / ( DENS0(i,pv,j) * POT0(i,pv,j) )
      end do
      end do
      end do

      do j=1, jm
        do pv2=1, elem%Nnode_v
        do pv=1, elem%Nnode_v
        do i=1, im
          ij = i + (j-1)*im
          p  = ij + (pv -1)*elem%Nnode_h1D**2
          p2 = ij + (pv2-1)*elem%Nnode_h1D**2

          ! DDENS
          BndMatD(i,DENS_VID_LC,pv,DENS_VID_LC,pv2,j,ke2D) = Id(pv,pv2)
          BndMatD(i,DENS_VID_LC,pv,MOMZ_VID_LC,pv2,j,ke2D) = fac_dz(i,pv,pv2,j)

          ! MOMZ
          BndMatD(i,MOMZ_VID_LC,pv,MOMZ_VID_LC,pv2,j,ke2D) = Id(pv,pv2) ! &
                                            !+ fac_dz_p(i,pv,pv2,j) * 2.0_RP * W0(i,pv2,j)  [ <- d_MOMZ ( MOMZ x MOMZ / DENS ) ]
          BndMatD(i,MOMZ_VID_LC,pv,DENS_VID_LC,pv2,j,ke2D) = impl_fac * Grav * IntrpMat_VPOrdM1(p,elem%Colmask(pv2,ij)) ! &
                                            !- fac_dz_p(i,pv,pv2,j) * W0(i,pv2,j)**2        [ <- d_DENS ( MOMZ x MOMZ / DENS ) ]
          BndMatD(i,MOMZ_VID_LC,pv,RHOT_VID_LC,pv2,j,ke2D) = fac_dz(i,pv,pv2,j) * DPDRHOT0(i,pv2,j)

          ! DRHOT
          BndMatD(i,RHOT_VID_LC,pv,DENS_VID_LC,pv2,j,ke2D) = - fac_dz(i,pv,pv2,j) * POT0(i,pv2,j) * WT0(i,pv2,j)
          BndMatD(i,RHOT_VID_LC,pv,MOMZ_VID_LC,pv2,j,ke2D) =   fac_dz(i,pv,pv2,j) * POT0(i,pv2,j)
          BndMatD(i,RHOT_VID_LC,pv,RHOT_VID_LC,pv2,j,ke2D) = Id(pv,pv2) + fac_dz(i,pv,pv2,j) * WT0(i,pv2,j)
        end do
        end do
        end do
        do v=1, 3
          BndMatL(:,:,:,v,j,ke2D) = 0.0_RP
          BndMatU(:,:,:,v,j,ke2D) = 0.0_RP
        end do
        do f1=1, 2
          if (f1==1) then
            ke_z2 = max(ke_z-1,1)
            pv1 = 1; f2 = 2
          else
            ke_z2 = min(ke_z+1,lmesh%NeZ)
            pv1 = elem%Nnode_v; f2 = 1
          end if
          ke_2 = ke2D + (ke_z2-1)*lmesh%Ne2D

          if ( (ke_z == 1 .and. f1==1) .or. (ke_z == lmesh%NeZ .and. f1==elem%Nfaces_v) ) then
            bc_flag = .true.
            f2 = f1
          else 
            bc_flag = .false.    
          end if

          do pv=1, elem%Nnode_v
          do i=1, im
            ij = i + (j-1)*im
            fp  = elem%Nfp_h * elem%Nfaces_h + (f1-1)*elem%Nfp_v + ij
            fp2 = elem%Nfp_h * elem%Nfaces_h + (f2-1)*elem%Nfp_v + ij
            p = ij + (pv-1)*im*jm

            !--
            fac  = 0.5_RP * impl_fac / GsqrtV(p,ke)
            tmp1(i,pv) = fac * elem%Lift(p,fp) * lmesh%Fscale(fp,ke) &
                * max( alph(fp,ke), alph(fp2,ke_2) )
          end do
          end do

          if (bc_flag) then
            do pv=1, elem%Nnode_v
            do i=1, im
              BndMatD(i,MOMZ_VID_LC,pv,MOMZ_VID_LC,pv1,j,ke2D) = BndMatD(i,MOMZ_VID_LC,pv,MOMZ_VID_LC,pv1,j,ke2D) + 2.0_RP * tmp1(i,pv)
            end do
            end do
          else
            do v=1, 3
            do pv=1, elem%Nnode_v
            do i=1, im
              BndMatD(i,v,pv,v,pv1,j,ke2D) = BndMatD(i,v,pv,v,pv1,j,ke2D) + tmp1(i,pv)
            end do
            end do
            end do
            if (f1 == 1) then
              do v=1, 3
              do pv=1, elem%Nnode_v
              do i=1, im
                BndMatL(i,v,pv,v,j,ke2D) = - tmp1(i,pv)
              end do
              end do
              end do
            else
              do v=1, 3
              do pv=1, elem%Nnode_v
              do i=1, im
                BndMatU(i,v,pv,v,j,ke2D) = - tmp1(i,pv)
              end do
              end do
              end do
            end if
          end if

          !--
          do pv=1, elem%Nnode_v
          do i=1, im
            ij = i + (j-1)*im
            fp  = elem%Nfp_h * elem%Nfaces_h + (f1-1)*elem%Nfp_v + ij
            p = ij + (pv-1)*im*jm
            !--
            fac  = 0.5_RP * impl_fac / GsqrtV(p,ke)
            tmp1(i,pv) = fac * elem%Lift(p,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke)
          end do
          end do
          if ( bc_flag ) then
            do pv=1, elem%Nnode_v
            do i=1, im
              ij = i + (j-1)*im
              p  = ij + (pv -1)*im*jm
              BndMatD(i,DENS_VID_LC,pv,MOMZ_VID_LC,pv1,j,ke2D) = BndMatD(i,DENS_VID_LC,pv,MOMZ_VID_LC,pv1,j,ke2D) - 2.0_RP * tmp1(i,pv)
              BndMatD(i,RHOT_VID_LC,pv,MOMZ_VID_LC,pv1,j,ke2D) = BndMatD(i,RHOT_VID_LC,pv,MOMZ_VID_LC,pv1,j,ke2D) - 2.0_RP * tmp1(i,pv) * POT0(i,pv1,j)
              BndMatD(i,RHOT_VID_LC,pv,DENS_VID_LC,pv1,j,ke2D) = BndMatD(i,RHOT_VID_LC,pv,DENS_VID_LC,pv1,j,ke2D) + 2.0_RP * tmp1(i,pv) * POT0(i,pv1,j) * WT0(i,pv1,j)
              BndMatD(i,RHOT_VID_LC,pv,RHOT_VID_LC,pv1,j,ke2D) = BndMatD(i,RHOT_VID_LC,pv,RHOT_VID_LC,pv1,j,ke2D) - 2.0_RP * tmp1(i,pv) * WT0(i,pv1,j)
            end do
            end do
          else
            do pv=1, elem%Nnode_v
            do i=1, im
              ij = i + (j-1)*im
              p  = ij + (pv -1)*im*jm
              BndMatD(i,DENS_VID_LC,pv,MOMZ_VID_LC,pv1,j,ke2D) = BndMatD(i,DENS_VID_LC,pv,MOMZ_VID_LC,pv1,j,ke2D) - tmp1(i,pv)

              BndMatD(i,MOMZ_VID_LC,pv,RHOT_VID_LC,pv1,j,ke2D) = BndMatD(i,MOMZ_VID_LC,pv,RHOT_VID_LC,pv1,j,ke2D) - tmp1(i,pv) * DPDRHOT0(i,pv1,j)

              BndMatD(i,RHOT_VID_LC,pv,MOMZ_VID_LC,pv1,j,ke2D) = BndMatD(i,RHOT_VID_LC,pv,MOMZ_VID_LC,pv1,j,ke2D) - tmp1(i,pv) * POT0(i,pv1,j)
              BndMatD(i,RHOT_VID_LC,pv,DENS_VID_LC,pv1,j,ke2D) = BndMatD(i,RHOT_VID_LC,pv,DENS_VID_LC,pv1,j,ke2D) + tmp1(i,pv) * POT0(i,pv1,j) * WT0(i,pv1,j)
              BndMatD(i,RHOT_VID_LC,pv,RHOT_VID_LC,pv1,j,ke2D) = BndMatD(i,RHOT_VID_LC,pv,RHOT_VID_LC,pv1,j,ke2D) - tmp1(i,pv) * WT0(i,pv1,j)
            end do
            end do
            if ( f1 == 1 ) then
              do pv=1, elem%Nnode_v
              do i=1, im
                ij = i + (j-1)*im
                p  = ij + (pv -1)*im*jm
                BndMatL(i,DENS_VID_LC,pv,MOMZ_VID_LC,j,ke2D) = + tmp1(i,pv)
    
                BndMatL(i,MOMZ_VID_LC,pv,RHOT_VID_LC,j,ke2D) = + tmp1(i,pv) * DPDRHOT0(i,pv1,j)
    
                BndMatL(i,RHOT_VID_LC,pv,MOMZ_VID_LC,j,ke2D) = + tmp1(i,pv) * POT0(i,pv1,j)
                BndMatL(i,RHOT_VID_LC,pv,DENS_VID_LC,j,ke2D) = - tmp1(i,pv) * POT0(i,pv1,j) * WT0(i,pv1,j)
                BndMatL(i,RHOT_VID_LC,pv,RHOT_VID_LC,j,ke2D) = BndMatL(i,RHOT_VID_LC,pv,RHOT_VID_LC,j,ke2D) &
                                                             + tmp1(i,pv) * WT0(i,pv1,j)
              end do
              end do            
            else
              do pv=1, elem%Nnode_v
              do i=1, im
                ij = i + (j-1)*im
                p  = ij + (pv -1)*im*jm
                BndMatU(i,DENS_VID_LC,pv,MOMZ_VID_LC,j,ke2D) = + tmp1(i,pv)
    
                BndMatU(i,MOMZ_VID_LC,pv,RHOT_VID_LC,j,ke2D) = + tmp1(i,pv) * DPDRHOT0(i,pv1,j)
    
                BndMatU(i,RHOT_VID_LC,pv,MOMZ_VID_LC,j,ke2D) = + tmp1(i,pv) * POT0(i,pv1,j)
                BndMatU(i,RHOT_VID_LC,pv,DENS_VID_LC,j,ke2D) = - tmp1(i,pv) * POT0(i,pv1,j) * WT0(i,pv1,j)
                BndMatU(i,RHOT_VID_LC,pv,RHOT_VID_LC,j,ke2D) = BndMatU(i,RHOT_VID_LC,pv,RHOT_VID_LC,j,ke2D) &
                                                             + tmp1(i,pv) * WT0(i,pv1,j)
              end do
              end do            
            end if   
          end if

        end do ! end do for f
      end do ! end for j
    end do ! end do for ke2D

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd

!OCL SERIAL  
  subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd_uv( &
    BndMatL_uv, BndMatD_uv, BndMatU_uv,     & ! (out)
    alph, GsqrtV,                           & ! (in)
    impl_fac, dt,                           & ! (in)
    lmesh, elem, im, jm,                    & ! (in)
    vmapM, vmapP, ke_z                      ) ! (in)

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: im, jm
    real(RP), intent(out) :: BndMatL_uv(im,elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP), intent(out) :: BndMatD_uv(im,elem%Nnode_v,elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP), intent(out) :: BndMatU_uv(im,elem%Nnode_v,jm,lmesh%Ne2D)
    real(RP), intent(in) :: alph(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: GsqrtV(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: dt
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)    
    integer, intent(in) :: ke_z

    integer :: Colmask(elem%Nnode_v)
    real(RP) :: Id(elem%Nnode_v,elem%Nnode_v)
    real(RP) :: Dd(elem%Nnode_v)
    real(RP) :: tmp1(im,elem%Nnode_v)

    integer :: p, pv, pv1, pv2
    integer :: fp,fp2
    integer :: ke, ke2D
    integer :: ke_z2, ke_2
    integer :: ij, i, j
    integer :: f1, f2

    real(RP) :: fac
    logical :: bc_flag
    !---------------------------------------------------------------------------

    Id(:,:) = 0.0_RP
    do p=1, elem%Nnode_v
      Id(p,p) = 1.0_RP
    end do

    do ke2D=1, lmesh%Ne2D
    do j=1, jm
      ke = ke2D + (ke_z-1)*lmesh%Ne2D

      do pv2=1, elem%Nnode_v
      do pv=1, elem%Nnode_v
      do i=1, im
        BndMatD_uv(i,pv,pv2,j,ke2D) = Id(pv,pv2)
      end do
      end do
      end do
      BndMatL_uv(:,:,j,ke2D) = 0.0_RP
      BndMatU_uv(:,:,j,ke2D) = 0.0_RP

      do f1=1, 2
        if (f1==1) then
          ke_z2 = max(ke_z-1,1)
          pv1 = 1; f2 = 2
        else
          ke_z2 = min(ke_z+1,lmesh%NeZ)
          pv1 = elem%Nnode_v; f2 = 1
        end if
        ke_2 = ke2D + (ke_z2-1)*lmesh%Ne2D

        if ( (ke_z == 1 .and. f1==1) .or. (ke_z == lmesh%NeZ .and. f1==elem%Nfaces_v) ) then
          bc_flag = .true.
          f2 = f1
        else 
          bc_flag = .false.    
        end if

        do pv=1, elem%Nnode_v
        do i=1, im
          ij = i + (j-1)*im
          fp  = elem%Nfp_h * elem%Nfaces_h + (f1-1)*elem%Nfp_v + ij
          fp2 = elem%Nfp_h * elem%Nfaces_h + (f2-1)*elem%Nfp_v + ij
          p = ij + (pv-1)*im*jm

          !--
          fac  = 0.5_RP * impl_fac / GsqrtV(p,ke)
          tmp1(i,pv) = fac * elem%Lift(p,fp) * lmesh%Fscale(fp,ke) &
              * max( alph(fp,ke), alph(fp2,ke_2) )
        end do
        end do

        if (bc_flag) then
        else
          BndMatD_uv(:,:,pv1,j,ke2D) = BndMatD_uv(:,:,pv1,j,ke2D) + tmp1(:,:)
          if (f1 == 1) then
            BndMatL_uv(:,:,j,ke2D) = - tmp1(:,:)
          else
            BndMatU_uv(:,:,j,ke2D) = - tmp1(:,:)
          end if
        end if
      end do ! end do for f

    end do ! end do for j
    end do ! end do for ke2D

    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_hevi_common_construct_matbnd_uv

!-- private ----------------

!OCL SERIAL  
  subroutine vi_cal_del_flux_dyn_uv( del_flux, alph,          & ! (out)
    PVARS_, PVARS0_,                                          & ! (in)
    DENS_hyd, PRES_hyd, Rtot, CPtot, CVtot,                   & ! (in)
    Gsqrt, G11, G12, G22, GsqrtH, gam, GsqrtV, G13, G23,      & ! (in)
    nz, vmapM, vmapP, iM2Dto3D, lmesh, elem, lmesh2D, elem2D  ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,2,lmesh%Ne)
    real(RP), intent(out) :: alph(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  PVARS_ (elem%Np*lmesh%Ne,PRGVAR_NUM)
    real(RP), intent(in) ::  PVARS0_(elem%Np*lmesh%Ne,PRGVAR_NUM)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Rtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CPtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CVtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G11(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G12(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  G22(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) ::  GsqrtH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: gam(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GsqrtV(elem%Np*lmesh%Ne)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: iM2Dto3D(elem%NfpTot)
    
    integer :: fp, ke, ke2D
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: MOMX0(elem%NfpTot,2)
    real(RP) :: MOMY0(elem%NfpTot,2)
    real(RP) :: DRHOT0(elem%NfpTot,2)
    real(RP) :: rhot_hyd(elem%NfpTot,2)
    real(RP) :: CPtot_ov_CVtot(elem%NfpTot,2)
    real(RP) :: Rtot_(elem%NfpTot,2)
    real(RP) :: wt0(elem%NfpTot,2)
    real(RP) :: pres0(elem%NfpTot,2)
    real(RP) :: rdens0(elem%NfpTot,2)

    real(RP) :: Gsqrt_(elem%NfpTot,2)
    real(RP) :: GsqrtV_(elem%NfpTot,2)
    real(RP) :: RGsqrtV(elem%NfpTot,2)
    real(RP) :: rgam2(elem%NfpTot,2)
    real(RP) :: G13_(elem%NfpTot,2)
    real(RP) :: G23_(elem%NfpTot,2)
    real(RP) :: Gxz_(elem%NfpTot,2)
    real(RP) :: Gyz_(elem%NfpTot,2)
    real(RP) :: G11_, G12_, G22_
    real(RP) :: Gnn_M, Gnn_P

    real(RP) :: tmp1

    real(RP) :: gamm, rgamm, PRES0ovRdry, rP0

    integer, parameter :: IN = 1
    integer, parameter :: EX = 2
    !------------------------------------------------------------------------

    gamm  = CpDry / CvDry
    rgamm = CvDry / CpDry
    PRES0ovRdry = PRES00 / Rdry
    rP0 = 1.0_RP / PRES00
    
    !$omp parallel private( ke, ke2D, fp, iM, iP,                        &
    !$omp MOMX0, MOMY0, DRHOT0, wt0, rdens0, pres0,                      &
    !$omp CPtot_ov_CVtot, Rtot_, rhot_hyd,                               &
    !$omp Gsqrt_, GsqrtV_, RGsqrtV, G11_, G12_, G22_, G13_, G23_, rgam2, &
    !$omp Gxz_, Gyz_, Gnn_M, Gnn_P,                                      &
    !$omp tmp1 )  

    !$omp do
    do ke=lmesh%NeS, lmesh%NeE
      ke2D = lmesh%EMap3Dto2D(ke)
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)

      Gsqrt_(:,IN) = Gsqrt(iM)
      Gsqrt_(:,EX) = Gsqrt(iP)

      rgam2(:,IN) = 1.0_RP / gam(iM)**2
      rgam2(:,EX) = 1.0_RP / gam(iP)**2
      GsqrtV_(:,IN) = GsqrtV(iM)
      GsqrtV_(:,EX) = GsqrtV(iP)
      RGsqrtV(:,:) = 1.0_RP / GsqrtV_(:,:)      

      G13_(:,IN) = G13(iM)
      G13_(:,EX) = G13(iP)
      G23_(:,IN) = G23(iM)
      G23_(:,EX) = G23(iP)

      MOMX0(:,IN) = PVARS0_(iM(:),MOMX_VID)
      MOMX0(:,EX) = PVARS0_(iP(:),MOMX_VID)
      MOMY0(:,IN) = PVARS0_(iM(:),MOMY_VID)
      MOMY0(:,EX) = PVARS0_(iP(:),MOMY_VID)
      DRHOT0(:,IN) = PVARS0_(iM(:),RHOT_VID)
      DRHOT0(:,EX) = PVARS0_(iP(:),RHOT_VID)
      
      Rtot_(:,IN) = Rtot(iM(:))
      Rtot_(:,EX) = Rtot(iP(:))
      CPtot_ov_CVtot(:,IN) = CPtot(iM(:)) / CVtot(iM(:))
      CPtot_ov_CVtot(:,EX) = CPtot(iP(:)) / CVtot(iP(:))

      !-
      rdens0(:,IN) = 1.0_RP / ( DENS_hyd(iM) + PVARS0_(iM,DENS_VID) )
      rdens0(:,EX) = 1.0_RP / ( DENS_hyd(iP) + PVARS0_(iP,DENS_VID) )

      wt0(:,IN) = ( PVARS0_(iM(:),MOMZ_VID) * RGsqrtV(:,IN) + G13_(:,IN) * MOMX0(:,IN)  &
                                                            + G23_(:,IN) * MOMY0(:,IN)  ) * rdens0(:,IN)     
      wt0(:,EX) = ( PVARS0_(iP(:),MOMZ_VID) * RGsqrtV(:,EX) + G13_(:,EX) * MOMX0(:,EX)  &
                                                            + G23_(:,EX) * MOMY0(:,EX)  ) * rdens0(:,EX)
      
      rhot_hyd(:,IN) = PRES0ovRdry * (PRES_hyd(iM(:))/PRES00)**rgamm
      rhot_hyd(:,EX) = PRES0ovRdry * (PRES_hyd(iP(:))/PRES00)**rgamm
      pres0(:,:) = PRES00 * ( Rtot_(:,:) * rP0 * ( rhot_hyd(:,:) + DRHOT0(:,:) ) )**CPtot_ov_CVtot(:,:)

      !-
      do fp=1, elem%NfpTot
        G11_ = G11(iM2Dto3D(fp),ke2D);  G12_ = G12(iM2Dto3D(fp),ke2D);  G22_ = G22(iM2Dto3D(fp),ke2D)

        Gxz_(fp,IN) = rgam2(fp,IN) * ( G11_ * G13_(fp,IN) + G12_ * G23_(fp,IN) )
        Gxz_(fp,EX) = rgam2(fp,EX) * ( G11_ * G13_(fp,EX) + G12_ * G23_(fp,EX) )

        Gyz_(fp,IN) = rgam2(fp,IN) * ( G12_ * G13_(fp,IN) + G22_ * G23_(fp,IN) )
        Gyz_(fp,EX) = rgam2(fp,EX) * ( G12_ * G13_(fp,EX) + G22_ * G23_(fp,EX) )
      end do

      do fp=1, elem%NfpTot
        Gnn_M = ( 1.0_RP * RGsqrtV(fp,IN)**2 + G13_(fp,IN) * Gxz_(fp,IN) + G23_(fp,IN) * Gyz_(fp,IN) ) * abs( nz(fp,ke) )
        Gnn_P = ( 1.0_RP * RGsqrtV(fp,EX)**2 + G13_(fp,EX) * Gxz_(fp,EX) + G23_(fp,EX) * Gyz_(fp,EX) ) * abs( nz(fp,ke) )

        alph(fp,ke) = nz(fp,ke)**2 * max( abs( wt0(fp,IN) ) + sqrt( Gnn_M * gamm * pres0(fp,IN) * rdens0(fp,IN) ), &
                                          abs( wt0(fp,EX) ) + sqrt( Gnn_P * gamm * pres0(fp,EX) * rdens0(fp,EX) )  )
      end do

      !---- 
      do fp=1, elem%NfpTot
        tmp1 = - 0.5_RP * lmesh%Fscale(fp,ke) * alph(fp,ke)
        del_flux(fp,1,ke) = tmp1 * ( PVARS_(iP(fp),MOMX_VID) - PVARS_(iM(fp),MOMX_VID) )
        del_flux(fp,2,ke) = tmp1 * ( PVARS_(iP(fp),MOMY_VID) - PVARS_(iM(fp),MOMY_VID) )                       
      end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine vi_cal_del_flux_dyn_uv

!OCL SERIAL  
  subroutine vi_cal_del_flux_dyn( del_flux, alph,            & ! (out)
    PVARS_, PVARS0_,                                         & ! (in)
    DENS_hyd, PRES_hyd, Rtot, CPtot, CVtot,                  & ! (in)
    Gsqrt, GsqrtH, gam, GsqrtV, G13, G23,                    & ! (in)
    nz, vmapM, vmapP, iM2Dto3D, lmesh, elem, lmesh2D, elem2D ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D 
    real(RP), intent(out) ::  del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP), intent(in) :: alph(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  PVARS_ (elem%Np*lmesh%Ne,PRGVAR_NUM)
    real(RP), intent(in) ::  PVARS0_(elem%Np*lmesh%Ne,PRGVAR_NUM)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Rtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CPtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CVtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GsqrtH(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: gam(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GsqrtV(elem%Np*lmesh%Ne)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: iM2Dto3D(elem%NfpTot)
    
    integer :: fp, ke, ke2D, ke_z, ke_xy
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: DDENS(elem%NfpTot,2)
    real(RP) :: MOMZ(elem%NfpTot,2)
    real(RP) :: MOMW(elem%NfpTot,2)
    real(RP) :: DRHOT(elem%NfpTot,2)
    real(RP) :: rhot_hyd(elem%NfpTot,2)
    real(RP) :: CPtot_ov_CVtot(elem%NfpTot,2)
    real(RP) :: Rtot_(elem%NfpTot,2)
    real(RP) :: pott(elem%NfpTot,2)
    real(RP) :: dpres(elem%NfpTot,2)
    real(RP) :: dens(elem%NfpTot,2)

    real(RP) :: Gsqrt_(elem%NfpTot,2)
    real(RP) :: GsqrtV_(elem%NfpTot,2)
    real(RP) :: rgam2(elem%NfpTot,2)
    real(RP) :: G13_(elem%NfpTot,2)
    real(RP) :: G23_(elem%NfpTot,2)

    real(RP) :: tmp1

    real(RP) :: rgamm, PRES0ovRdry, rP0

    integer, parameter :: IN = 1
    integer, parameter :: EX = 2
    !------------------------------------------------------------------------

    rgamm = CvDry / CpDry
    PRES0ovRdry = PRES00 / Rdry
    rP0 = 1.0_RP / PRES00
    
    !$omp parallel private( ke_z, ke_xy, ke, ke2D, fp, iM, iP, &
    !$omp DDENS, MOMZ, MOMW, DRHOT,                           &
    !$omp CPtot_ov_CVtot, Rtot_, rhot_hyd, pott, dpres, dens, &
    !$omp Gsqrt_, GsqrtV_, G13_, G23_, rgam2                  )  

    !$omp do collapse(2)
    do ke_z=1, lmesh%NeZ
    do ke_xy=1, lmesh%Ne2D
      ke = ke_xy + (ke_z-1)*lmesh%Ne2D
      ke2D = lmesh%EMap3Dto2D(ke)

      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)

      Gsqrt_(:,IN) = Gsqrt(iM)
      Gsqrt_(:,EX) = Gsqrt(iP)

      rgam2(:,IN) = 1.0_RP / gam(iM)**2
      rgam2(:,EX) = 1.0_RP / gam(iP)**2
      GsqrtV_(:,IN) = GsqrtV(iM)
      GsqrtV_(:,EX) = GsqrtV(iP)

      G13_(:,IN) = G13(iM)
      G13_(:,EX) = G13(iP)
      G23_(:,IN) = G23(iM)
      G23_(:,EX) = G23(iP)

      DDENS(:,IN) = PVARS_(iM(:),DENS_VID)
      DDENS(:,EX) = PVARS_(iP(:),DENS_VID)
      MOMZ(:,IN) = PVARS_(iM(:),MOMZ_VID)
      MOMZ(:,EX) = PVARS_(iP(:),MOMZ_VID)
      DRHOT(:,IN) = PVARS_(iM(:),RHOT_VID)
      DRHOT(:,EX) = PVARS_(iP(:),RHOT_VID)

      Rtot_(:,IN) = Rtot(iM(:))
      Rtot_(:,EX) = Rtot(iP(:))
      CPtot_ov_CVtot(:,IN) = CPtot(iM(:)) / CVtot(iM(:))
      CPtot_ov_CVtot(:,EX) = CPtot(iP(:)) / CVtot(iP(:))
      !-

      dens(:,IN) = DENS_hyd(iM(:)) + DDENS(:,IN)
      dens(:,EX) = DENS_hyd(iP(:)) + DDENS(:,EX)

      rhot_hyd(:,IN) = PRES0ovRdry * (PRES_hyd(iM(:))/PRES00)**rgamm
      rhot_hyd(:,EX) = PRES0ovRdry * (PRES_hyd(iP(:))/PRES00)**rgamm      
      pott(:,:) = ( rhot_hyd(:,:) + DRHOT(:,:) ) / dens(:,:)

      dpres(:,:) = PRES00 * ( Rtot_(:,:) * rP0 * dens(:,:) * pott(:,:) )**CPtot_ov_CVtot(:,:)
      dpres(:,IN) = dpres(:,IN) - PRES_hyd(iM(:)) 
      dpres(:,EX) = dpres(:,EX) - PRES_hyd(iP(:)) 
      
      !-
      MOMW(:,IN) = MOMZ(:,IN) &
              + GsqrtV_(:,IN) * G13_(:,IN) * PVARS_(iM,MOMX_VID) &
              + GsqrtV_(:,IN) * G23_(:,IN) * PVARS_(iM,MOMY_VID)
      MOMW(:,EX) = MOMZ(:,EX) &
              + GsqrtV_(:,EX) * G13_(:,EX) * PVARS_(iP,MOMX_VID) &
              + GsqrtV_(:,EX) * G23_(:,EX) * PVARS_(iP,MOMY_VID)
      
      !----                  
                    
      if ( ke_z == 1 .or. ke_z == lmesh%NeZ) then
        do fp=1, elem%NfpTot
          if ( iM(fp) == iP(fp) ) then
            MOMZ(fp,EX) = - PVARS_(fp,MOMZ_VID) &
                          - 2.0_RP * GsqrtV_(fp,IN) * (  G13_(fp,IN) * PVARS_(fp,MOMX_VID) &
                                                       + G23_(fp,IN) * PVARS_(fp,MOMY_VID) )          
            MOMW(fp,EX) = - MOMW(fp,IN)
          end if
        end do
      end if

      do fp=1, elem%NfpTot
        tmp1 = 0.5_RP * lmesh%Fscale(fp,ke)

        del_flux(fp,DENS_VID,ke) = tmp1 * ( &
                  + ( MOMW(fp,EX) - MOMW(fp,IN) ) * nz(fp,ke)    &
                  - alph(fp,ke) * ( DDENS(fp,EX) - DDENS(fp,IN) ) )
        
        
        del_flux(fp,MOMZ_VID,ke) = tmp1 * ( &
                  + ( dpres(fp,EX) - dpres(fp,IN) ) * nz(fp,ke) &                    
                  - alph(fp,ke) * ( MOMZ(fp,EX) - MOMZ(fp,IN) ) )
        
        del_flux(fp,RHOT_VID,ke) = tmp1 * ( &
                  + ( pott(fp,EX) * MOMW(fp,EX)  - pott(fp,IN) * MOMW(fp,IN) ) * nz(fp,ke) &
                  - alph(fp,ke) * ( DRHOT(fp,EX) - DRHOT(fp,IN) ) )
      end do
    end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine vi_cal_del_flux_dyn
end module scale_atm_dyn_dgm_nonhydro3d_rhot_hevi_common_2
