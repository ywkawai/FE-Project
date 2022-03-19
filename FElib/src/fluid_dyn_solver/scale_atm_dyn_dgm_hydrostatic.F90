!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!      Construct hydrostatic state for Atmospheric dynamical process. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_hydrostatic
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc

  use scale_const, only: &
    PI => CONST_PI,       &
    GRAV => CONST_GRAV,   &
    Rdry => CONST_Rdry,   &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00

  use scale_gmres, only: GMRES
  use scale_sparsemat, only: &
    SparseMat,               &
    SparseMat_matmul
  
  use scale_element_base, only: &
    ElementBase, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: hydrostatic_calc_basicstate_constT
  public :: hydrostatic_calc_basicstate_constPT
  public :: hydrostatic_calc_basicstate_constPTLAPS
  public :: hydrostatic_calc_basicstate_constBVFreq
  public :: hydrostaic_build_rho_XYZ

  interface hydrostaic_build_rho_XYZ
    module procedure hydrostaic_build_rho_XYZ_dry
    module procedure hydrostaic_build_rho_XYZ_moist
  end interface
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !

contains
!OCL SERIAL
  subroutine hydrostatic_calc_basicstate_constT( &
    DENS_hyd, PRES_hyd,                          &
    Temp0, PRES_sfc, x, y, z, lcmesh3D, elem     )

    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh3D%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: Temp0
    real(RP), intent(in) :: PRES_sfc

    integer :: ke
    real(RP) :: H0
    !-----------------------------------------------

    H0 = Rdry * Temp0 / Grav

    !$omp parallel do
    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      PRES_hyd(:,ke) = PRES_sfc * exp( - z(:,ke) / H0 )
      DENS_hyd(:,ke) = PRES_hyd(:,ke) / ( Rdry * Temp0 )
    end do

    return
  end subroutine hydrostatic_calc_basicstate_constT

!OCL SERIAL
  subroutine hydrostatic_calc_basicstate_constPT( &
    DENS_hyd, PRES_hyd,                         &
    PotTemp0, PRES_sfc, x, y, z, lcmesh3D, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh3D%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: PotTemp0
    real(RP), intent(in) :: PRES_sfc

    integer :: ke

    real(RP) :: exner(elem%Np)
    real(RP) :: exner_sfc
    real(RP) :: RovCP
    real(RP) :: CPovR
    !-----------------------------------------------

    RovCP = Rdry / CpDry
    CPovR = CpDry / Rdry
    exner_sfc = (PRES_sfc / PRES00)**RovCP

    !$omp parallel do private(exner)
    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      ! Cp * PT0 * d exner / dz = - g
      exner(:) = exner_sfc - Grav / (CpDry * PotTemp0) * z(:,ke)
      PRES_hyd(:,ke) = PRES00 * exner(:)**CPovR
      DENS_hyd(:,ke) =  PRES_hyd(:,ke) / ( Rdry * exner(:) * PotTemp0 )
    end do

    return
  end subroutine hydrostatic_calc_basicstate_constPT

!OCL SERIAL
  subroutine hydrostatic_calc_basicstate_constBVFreq( &
    DENS_hyd, PRES_hyd,                                           &
    BruntVaisalaFreq, PotTemp0, PRES_sfc, x, y, z, lcmesh3D, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh3D%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: BruntVaisalaFreq 
    real(RP), intent(in) :: PotTemp0
    real(RP), intent(in) :: PRES_sfc

    integer :: ke

    real(RP) :: PT(elem%NP)
    real(RP) :: exner(elem%Np)
    real(RP) :: exner_sfc
    real(RP) :: RovCP
    real(RP) :: CPovR
    !-----------------------------------------------

    RovCP = Rdry / CpDry
    CPovR = CpDry / Rdry
    exner_sfc = (PRES_sfc / PRES00)**RovCP

    !$omp parallel do private(PT, exner)
    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      ! d exner / dz = - g / ( Cp * PT0 ) * exp (- N2/g * z)
      ! exner = exner(zs) - g^2 / (Cp * N^2) [ 1/PT (z) - 1/PT(zs) ] 
      PT(:) = PotTemp0 * exp( BruntVaisalaFreq**2 / Grav * z(:,ke) )
      exner(:) = exner_sfc + Grav**2 / ( CpDry * BruntVaisalaFreq**2 ) * ( 1.0_RP / PT(:) - 1.0_RP / PotTemp0 )

      PRES_hyd(:,ke) = PRES00 * exner(:)**CPovR
      DENS_hyd(:,ke) =  PRES_hyd(:,ke) / ( Rdry * exner(:) * PT(:) )
    end do

    return
  end subroutine hydrostatic_calc_basicstate_constBVFreq

!OCL SERIAL
  subroutine hydrostatic_calc_basicstate_constPTLAPS( &
    DENS_hyd, PRES_hyd,                                 &
    PTLAPS, PotTemp0, PRES_sfc, x, y, z, lcmesh3D, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh3D%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh3D%Ne)
    real(RP), intent(in) :: PTLAPS
    real(RP), intent(in) :: PotTemp0
    real(RP), intent(in) :: PRES_sfc

    integer :: ke

    real(RP) :: PT(elem%NP)
    real(RP) :: exner(elem%Np)
    real(RP) :: exner_sfc
    real(RP) :: RovCP
    real(RP) :: CPovR
    !-----------------------------------------------

    RovCP = Rdry / CpDry
    CPovR = CpDry / Rdry
    exner_sfc = (PRES_sfc / PRES00)**RovCP

    !$omp parallel do private(PT, exner)
    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      ! d exner / dz = - g / ( Cp * PT0 ) / (1 + PTLAPS/PT0 * z)
      ! exner = exner(zs) - g / (Cp * PTLAPS ) * log[ 1 + PTLAPS/PT0 * z ] 
      PT(:) = PotTemp0 + PTLAPS * z(:,ke)
      exner(:) = exner_sfc - Grav / ( CpDry * PTLAPS ) * log( 1.0_RP + PTLAPS / PotTemp0 * z(:,ke) )

      PRES_hyd(:,ke) = PRES00 * exner(:)**CPovR
      DENS_hyd(:,ke) =  PRES_hyd(:,ke) / ( Rdry * exner(:) * PT(:) )
    end do

    return
  end subroutine hydrostatic_calc_basicstate_constPTLAPS

!OCL SERIAL
  subroutine hydrostaic_build_rho_XYZ_dry(   &
    DDENS,                                   &
    DENS_hyd, PRES_hyd,                      &
    POT,                                     &
    x, y, z, lcmesh, elem                    )

    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: POT(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)

    real(RP) :: Rtot          (elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: CPtot_ov_CVtot(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    !-----------------------------------------------

    !$omp parallel 
    !$omp workshare
    Rtot(:,:,:,:) = Rdry
    CPtot_ov_CVtot(:,:,:,:) = CPdry / CVdry
    !$omp end workshare
    !$omp end parallel

    call hydrostaic_build_rho_XYZ_moist(   &
      DDENS,                                   &
      DENS_hyd, PRES_hyd,                      &
      POT, Rtot, CPtot_ov_CVtot,               &
      x, y, z, lcmesh, elem                    )
    
    return
  end subroutine hydrostaic_build_rho_XYZ_dry

!OCL SERIAL
  subroutine hydrostaic_build_rho_XYZ_moist( &
    DDENS,                                   &
    DENS_hyd, PRES_hyd,                      &
    POT, Rtot, CPtot_ov_CVtot,               &
    x, y, z, lcmesh, elem                    )

    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: POT(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP), intent(in) :: Rtot(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP), intent(in) :: CPtot_ov_CVtot(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)

    integer :: ke
    integer :: ke_x, ke_y, ke_z
    integer :: itr_lin
    integer :: itr_nlin

    real(RP), parameter :: EPS0 = 1.0E-12_RP
    real(RP), parameter :: EPS  = 1.0E-12_RP

    type(SparseMat) :: Dz, Lift
    type(HexahedralElement) :: elem3D

    type(GMRES) :: gmres_hydro
    real(RP), allocatable :: wj(:)
    real(RP), allocatable :: pinv_v(:)
    integer :: N, m
    integer :: vmapM_z1D(elem%NfpTot,lcmesh%NeZ)
    integer :: vmapP_z1D(elem%NfpTot,lcmesh%NeZ) 
    real(RP) :: VARS  (elem%Np,lcmesh%NeZ)
    real(RP) :: VARS0 (elem%Np,lcmesh%NeZ)
    real(RP) :: VAR_DEL(elem%Np,lcmesh%NeZ)
    real(RP) :: b(elem%Np,lcmesh%NeZ)
    real(RP) :: Ax(elem%Np,lcmesh%NeZ)
    real(RP) :: nz(elem%NfpTot,lcmesh%NeZ)
    real(RP) :: DENS_hyd_z(elem%Np,lcmesh%NeZ)
    real(RP) :: PRES_hyd_z(elem%Np,lcmesh%NeZ)
    real(RP) :: PmatDlu(elem%Np,elem%Np,lcmesh%NeZ)
    integer :: PmatDlu_ipiv(elem%Np,lcmesh%NeZ)
    real(RP) :: PmatL(elem%Np,elem%Np,lcmesh%NeZ)
    real(RP) :: PmatU(elem%Np,elem%Np,lcmesh%NeZ)

    real(RP) :: linf_b

    logical :: is_converged

    real(RP) :: invV_POrdM1(elem%Np,elem%Np)
    real(RP) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)
    integer :: p1, p2, p_
    !-----------------------------------------------

    N = elem%Np * lcmesh%NeZ
    m = min(N / elem%Nnode_h1D**2, 30)
    call gmres_hydro%Init( N, m, EPS, EPS )
    allocate( wj(N), pinv_v(N) )

    call Dz%Init( elem%Dx3, storage_format='ELL' )
    call Lift%Init( elem%Lift, storage_format='ELL' )

    call set_vmapZ1D( vmapM_z1D, vmapP_z1D, & ! (out)
      elem, lcmesh )                          ! (in)

    invV_POrdM1(:,:) = elem%invV
    do p2=1, elem%Nnode_h1D
    do p1=1, elem%Nnode_h1D
      p_ = p1 + (p2-1)*elem%Nnode_h1D + (elem%Nnode_v-1)*elem%Nnode_h1D**2
      invV_POrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem%V, invV_POrdM1)

    call elem3D%Init( elem%PolyOrder_h, elem%PolyOrder_v, .false. )

    !-----------
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX

      do ke_z=1, lcmesh%NeZ
        ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1)*lcmesh%NeX*lcmesh%NeY
        VARS(:,ke_z) = 0.0_RP

        VARS0 (:,ke_z) = VARS(:,ke_z)
        nz(:,ke_z) = lcmesh%normal_fn(:,ke,3)

        DENS_hyd_z(:,ke_z) = DENS_hyd(:,ke)
        PRES_hyd_z(:,ke_z) = PRES_hyd(:,ke)
      end do

      do itr_nlin=1, 10

        do ke_z=1, lcmesh%NeZ
          VAR_DEL(:,ke_z) = 0.0_RP
        end do

        call eval_Ax( Ax(:,:), &
          VARS, VARS0, POT(:,:,ke_x,ke_y), Rtot(:,:,ke_x,ke_y),  &
          CPtot_ov_CVtot(:,:,ke_x,ke_y), DENS_hyd_z, PRES_hyd_z, &
          Dz, Lift, IntrpMat_VPOrdM1, lcmesh, elem,              &
          nz, vmapM_z1D, vmapP_z1D, ke_x, ke_y                   )
        
        do ke_z=1, lcmesh%NeZ
          b(:,ke_z) = - Ax(:,ke_z)
        end do
        if (lcmesh%tileID==1) then
          LOG_PROGRESS(*) ke_x, ke_y, "itr:", itr_nlin, 0, ": VAR", VARS(elem%Colmask(:,1),1)
          LOG_PROGRESS(*) ke_x, ke_y, "itr:", itr_nlin, 0, ": b", b(elem%Colmask(:,1),1)
          if( IO_L ) call flush(IO_FID_LOG)
        end if

        if ( maxval(abs(b(:,:))) < 1.0E-10_RP ) exit
        
        call construct_pmatInv( PmatDlu, PmatDlu_ipiv, PmatL, PmatU,  & ! (out)
          VARS0, POT(:,:,ke_x,ke_y), Rtot(:,:,ke_x,ke_y),             & ! (in)
          CPtot_ov_CVtot(:,:,ke_x,ke_y), DENS_hyd_z, PRES_hyd_z,      & ! (in)
          Dz, Lift, IntrpMat_VPOrdM1, lcmesh, elem,                   & ! (in)
          nz, vmapM_z1D, vmapP_z1D, ke_x, ke_y  )

        do itr_lin=1, 2*int(N/m)
          !
          call GMRES_hydro_core( gmres_hydro, VAR_DEL, wj, is_converged, &
            VARS, b, N, m,                                               &
            PmatDlu,  PmatDlu_ipiv,  PmatL, PmatU, pinv_v,               & ! (in)
            POT(:,:,ke_x,ke_y), Rtot(:,:,ke_x,ke_y),                     & ! (in)
            CPtot_ov_CVtot(:,:,ke_x,ke_y), DENS_hyd_z, PRES_hyd_z,       & ! (in)
            Dz, Lift, IntrpMat_VPOrdM1, lcmesh, elem,                    & ! (in)
            nz, vmapM_z1D, vmapP_z1D, ke_x, ke_y )
          
          if (is_converged) exit            
        end do ! itr_lin
        do ke_z=1, lcmesh%NeZ
          VARS (:,ke_z) = VARS(:,ke_z) + VAR_DEL(:,ke_z)
          VARS0(:,ke_z) = VARS(:,ke_z)
        end do
      end do ! itr_nlin

      do ke_z=1, lcmesh%NeZ
        ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1)*lcmesh%NeX*lcmesh%NeY
        DDENS(:,ke) = VARS(:,ke_z)
      end do
    end do  
    end do
      
    !
    call gmres_hydro%Final()
    call Dz%Final()
    call Lift%Final()

    return
  end subroutine hydrostaic_build_rho_XYZ_moist

!-- private ---------------------------------

!OCL SERIAL
  subroutine set_vmapZ1D( vmapM, vmapP, &
    elem, lmesh )

    implicit none
    type(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(out) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(out) :: vmapP(elem%NfpTot,lmesh%NeZ)    

    integer :: ke_z
    integer :: f
    integer :: vs, ve
    !------------------------------

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

      vs = elem%Nfp_h*elem%Nfaces_h + elem%Nfp_v + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z < lmesh%NeZ) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1) + ke_z*elem%Np
    end do

    return
  end subroutine set_vmapZ1D

!OCL SERIAL
  subroutine GMRES_hydro_core( gmres_hydro, x, wj, is_converged, &
    x0, b, N, m,                                        &
    PmatDlu, PmatDlu_ipiv, PmatL, PmatU, pinv_v,        & ! (in)    
    POT, Rtot, CPtot_ov_CVtot, DENS_hyd, PRES_hyd,      & ! (in)
    Dz, Lift, IntrpMat_VPOrdM1, lmesh, elem,            & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(in) :: N
    integer, intent(in) :: m    

    class(GMRES), intent(inout) :: gmres_hydro
    real(RP), intent(inout) :: x(N)
    real(RP), intent(inout) :: wj(N)
    logical, intent(out) :: is_converged
    real(RP), intent(in) :: x0(N)
    real(RP), intent(in) :: b(N)
    real(RP), intent(in) :: PmatDlu(elem%Np,elem%Np,lmesh%NeZ)
    integer, intent(in) :: PmatDlu_ipiv(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: PmatL(elem%Np,elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: PmatU(elem%Np,elem%Np,lmesh%NeZ)
    real(RP), intent(inout) :: pinv_v(N)
    !---
    real(RP), intent(in) :: POT(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: Rtot(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: CPtot_ov_CVtot(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y  
    
    integer :: j
    !-------------------------------------------------------

    call eval_Ax_lin( wj(:),                   & ! (out)
      x, x0, POT, Rtot, CPtot_ov_CVtot,        & ! (in)
      DENS_hyd, PRES_hyd,                      & ! (in)
      Dz, Lift, IntrpMat_VPOrdM1, lmesh, elem, & ! (in)
      nz, vmapM, vmapP, ke_x, ke_y             ) ! (in)

    call gmres_hydro%Iterate_pre( b, wj, is_converged )
    if (is_converged) return 

    do j=1, min(m, N)
      call matmul_pinv_v( pinv_v, PmatDlu, PmatDlu_ipiv, PmatL, PmatU, gmres_hydro%v(:,j) )

      call eval_Ax_lin( wj(:),                   & ! (out)
        pinv_v, x0, POT, Rtot, CPtot_ov_CVtot,   & ! (in)
        DENS_hyd, PRES_hyd,                      & ! (in)
        Dz, Lift, IntrpMat_VPOrdM1, lmesh, elem, & ! (in)
        nz, vmapM, vmapP, ke_x, ke_y             ) ! (in)
      
      call gmres_hydro%Iterate_step_j( j, wj, is_converged )
      if (is_converged) exit
    end do

    do j=1, N
      wj(j) = 0.0_RP
    end do
    call gmres_hydro%Iterate_post( wj )
    call matmul_pinv_v_plus_x0( x, PmatDlu, PmatDlu_ipiv, PmatL, pmatU, wj)

    return
  contains

!OCL SERIAL
    subroutine matmul_pinv_v( pinv_v_, pDlu_, PmatDlu_ipiv_, pL, pU, v)
      implicit none
      real(RP), intent(out) :: pinv_v_(elem%Np,lmesh%NeZ)      
      real(RP), intent(in) :: pDlu_(elem%Np,elem%Np,lmesh%NeZ)
      integer, intent(in) :: PmatDlu_ipiv_(elem%Np,lmesh%NeZ)
      real(RP), intent(in) :: pL(elem%Np,elem%Np,lmesh%NeZ)
      real(RP), intent(in) :: pU(elem%Np,elem%Np,lmesh%NeZ)
      real(RP), intent(in) :: v(elem%Np,lmesh%NeZ)      

      integer :: k, n
      real(RP) :: tmp(elem%Np)
      integer :: vs, ve
      integer :: info
      !------------------------------------
      
      n = elem%Np

      vs = 1
      ve = vs + elem%Np - 1
      pinv_v_(:,1) = v(vs:ve,1)
      call DGETRS('N', n, 1, pDlu_(:,:,1), n, PmatDlu_ipiv_(:,1), pinv_v_(:,1), n, info)

      do k=2, lmesh%NeZ
        vs = 1; ve = elem%Np
        pinv_v_(:,k) = v(vs:ve,k) &
          - matmul( pL(:,:,k), pinv_v_(:,k-1) )
        
        call DGETRS('N', n, 1, pDlu_(:,:,k), n, PmatDlu_ipiv_(:,k), pinv_v_(:,k), n, info)
      end do

      !
      do k=lmesh%NeZ-1, 1, -1
        vs = 1; ve = elem%Np
        tmp(vs:ve) = matmul( pU(:,:,k), pinv_v_(:,k+1) )
        call DGETRS('N', n, 1, pDlu_(:,:,k), n, PmatDlu_ipiv_(:,k), tmp(:), n, info)

        vs = 1
        ve = vs + elem%Np - 1  
        pinv_v_(:,k) = pinv_v_(:,k) - tmp(vs:ve)
      end do

      return
    end subroutine matmul_pinv_v

!OCL SERIAL
    subroutine matmul_pinv_v_plus_x0( x_, pDlu_, PmatDlu_ipiv_, pL, pU, v)
      implicit none
      real(RP), intent(inout) :: x_(elem%Np,lmesh%NeZ)      
      real(RP), intent(in) :: pDlu_(elem%Np,elem%Np,lmesh%NeZ)
      integer, intent(in) :: PmatDlu_ipiv_(elem%Np,lmesh%NeZ)
      real(RP), intent(in) :: pL(elem%Np,elem%Np,lmesh%NeZ)
      real(RP), intent(in) :: pU(elem%Np,elem%Np,lmesh%NeZ)
      real(RP), intent(in) :: v(elem%Np,lmesh%NeZ)      

      integer :: k
      real(RP) :: tmp(elem%Np,lmesh%NeZ)

      !------------------------------------
      
      call matmul_pinv_v( tmp, pDlu_, PmatDlu_ipiv_, pL, pU, v)
      !$omp parallel do
      do k=1, lmesh%NeZ
        x_(:,k) = x_(:,k) + tmp(:,k)
      end do

      return
    end subroutine matmul_pinv_v_plus_x0    

  end subroutine GMRES_hydro_core

!OCL SERIAL
  subroutine eval_Ax( Ax, &
      DDENS, DENS0, POT, Rtot, CPtot_ov_CVtot,    & ! (in)
      DENS_hyd, PRES_hyd,                         & ! (in)
      Dz, Lift, IntrpMat_VPOrdM1, lmesh, elem,    & ! (in)
      nz, vmapM, vmapP, ke_x, ke_y )
   
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem

    real(RP), intent(out) :: Ax(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: DDENS (elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: DENS0(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: POT(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: Rtot(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: CPtot_ov_CVtot(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift    
    real(RP), intent(in) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)

    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y

    real(RP) :: DPRES(elem%Np), DENS(elem%Np)
    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%NeZ)

    integer :: ke_z
    integer :: ke

    real(RP) :: gamm
    real(RP) :: RdOvP00
    real(RP) :: rP0
    !-------------------------------------------

    gamm = CpDry / CvDry
    RdOvP00 = Rdry / PRES00
    rP0 = 1.0_RP / PRES00

    call cal_del_flux( del_flux,       & ! (out)
      DDENS, POT,Rtot, CPtot_ov_CVtot, & ! (in)
      DENS_hyd, PRES_hyd,              & ! (in)        
      nz, vmapM, vmapP, lmesh, elem    ) ! (in)

    !$omp parallel do private(ke, DPRES, DENS, Fz, LiftDelFlx)
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

      DENS(:) = DENS_hyd(:,ke_z) + DDENS(:,ke_z)
      DPRES(:) = PRES00 * ( Rtot(:,ke_z) * rP0 * DENS(:) * POT(:,ke_z) )**CPtot_ov_CVtot(:,ke_z) !&
               !- PRES_hyd(:,ke_z)
      call sparsemat_matmul(Dz, DPRES, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z), LiftDelFlx)

      Ax(:,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) &
                 + Grav * matmul(IntrpMat_VPOrdM1, DENS(:))!DDENS(:,ke_z))
    end do

    return
  end subroutine eval_Ax

!OCL SERIAL
  subroutine cal_del_flux( del_flux, &
    DDENS_, POT_, Rtot_, CPtot_ov_CVtot_, &
    DENS_hyd, PRES_hyd,                   &
    nz, vmapM, vmapP, lmesh, elem         )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%NeZ)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  POT_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: Rtot_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: CPtot_ov_CVtot_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: DENS_hyd(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: PRES_hyd(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%NeZ)

    integer :: i, iP, iM
    integer :: p, ke_z

    real(RP) :: dpresP, dpresM
    real(RP) :: RtotOvP00M, RtotOvP00P
    real(RP) :: rP0
    !-------------------------------

    rP0 = 1.0_RP / PRES00

    !$omp parallel do private(p, i, iM, iP, dpresM, dpresP, RtotOvP00M, RtotOvP00P)
    do ke_z=1, lmesh%NeZ
    do p=1, elem%NfpTot
      i = p + (ke_z-1)*elem%NfpTot
      iM = vmapM(i); iP = vmapP(i)

      RtotOvP00M = Rtot_(iM) * rP0
      RtotOvP00P = Rtot_(iP) * rP0

      dpresM = PRES00 *  ( RtotOvP00M * (DENS_hyd(iM) + DDENS_(iM)) * POT_(iM) )**CPtot_ov_CVtot_(iM) !- PRES_hyd(iM)
      dpresP = PRES00 *  ( RtotOvP00P * (DENS_hyd(iP) + DDENS_(iP)) * POT_(iP) )**CPtot_ov_CVtot_(iP) !- PRES_hyd(iP)
      if (ke_z==1.and. iM==iP) dpresP = PRES_hyd(iP)!0.0_RP

      del_flux(i) = 0.5_RP * (dpresP - dpresM) * nz(i)
    end do
    end do

    return
  end subroutine cal_del_flux

!OCL SERIAL
  subroutine eval_Ax_lin( Ax, &
    DDENS, DDENS0, POT, Rtot, CPtot_ov_CVtot, & ! (in) 
    DENS_hyd, PRES_hyd,                       & ! (in)
    Dz, Lift, IntrpMat_VPOrdM1, lmesh, elem,  & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y )
 
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem

    real(RP), intent(out) :: Ax(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: DDENS (elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: DDENS0(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: POT(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: Rtot(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: CPtot_ov_CVtot(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y

    real(RP) :: PRES(elem%Np), PRES0(elem%Np)
    real(RP) :: DPRES(elem%Np)
    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%NeZ)

    integer :: ke_z
    integer :: ke

    real(RP) :: gamm
    real(RP) :: rP0
    !-------------------------------------------

    gamm = CpDry/CvDry
    rP0 = 1.0_RP / PRES00

    call cal_del_flux_lin( del_flux,    & ! (out)
      DDENS, DDENS0, POT, Rtot, CPtot_ov_CVtot, & ! (in) 
      DENS_hyd, PRES_hyd,                       & ! (in)        
      nz, vmapM, vmapP, lmesh, elem             ) ! (in)

    !$omp parallel do private(ke, PRES, PRES0, DPRES, Fz, LiftDelFlx)
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

      PRES0(:) = PRES00 * ( Rtot(:,ke_z) * rP0 * (DENS_hyd(:,ke_z) + DDENS0(:,ke_z)) * POT(:,ke_z) )**CPtot_ov_CVtot(:,ke_z)
      DPRES(:) = CPtot_ov_CVtot(:,ke_z) * PRES0(:) / (DENS_hyd(:,ke_z) + DDENS0(:,ke_z)) * DDENS(:,ke_z)

      call sparsemat_matmul(Dz, DPRES(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z), LiftDelFlx)
      
      Ax(:,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) &
                 + Grav * matmul(IntrpMat_VPOrdM1, DDENS(:,ke_z))
    end do

    return
  end subroutine eval_Ax_lin

!OCL SERIAL
  subroutine cal_del_flux_lin( del_flux, &
    DDENS_, DDENS0_, POT_, Rtot_, CPtot_ov_CVtot_, &
    DENS_hyd_, PRES_hyd_,                          &
    nz, vmapM, vmapP, lmesh, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%NeZ)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  DDENS0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  Rtot_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  CPtot_ov_CVtot_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  DENS_hyd_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  PRES_hyd_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  POT_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%NeZ)

    integer :: i, iP, iM
    integer :: p, ke_z

    real(RP) :: dpresP, dpresM
    real(RP) :: pres0P, pres0M
    real(RP) :: gamm
    real(RP) :: rP0
    real(RP) :: RtotOvP00M
    real(RP) :: RtotOvP00P
    !-------------------------------

    gamm = CpDry/CvDry
    rP0 = 1.0_RP / PRES00

    !$omp parallel do private(p, i, iM, iP, &
    !$omp dpresM, dpresP, pres0M, pres0P, RtotOvP00M, RtotOvP00P )
    do ke_z=1, lmesh%NeZ
    do p=1, elem%NfpTot
      i = p + (ke_z-1)*elem%NfpTot
      iM = vmapM(i); iP = vmapP(i)

      RtotOvP00M = Rtot_(iM) * rP0
      RtotOvP00P = Rtot_(iP) * rP0

      pres0M = PRES00 * ( RtotOvP00M * (DENS_hyd_(iM) + DDENS0_(iM)) * POT_(iM) )**CPtot_ov_CVtot_(iM)
      pres0P = PRES00 * ( RtotOvP00P * (DENS_hyd_(iP) + DDENS0_(iP)) * POT_(iP) )**CPtot_ov_CVtot_(iP)

      dpresM = CPtot_ov_CVtot_(iM) * pres0M / (DENS_hyd_(iM) + DDENS0_(iM)) * DDENS_(iM)
      dpresP = CPtot_ov_CVtot_(iP) * pres0P / (DENS_hyd_(iP) + DDENS0_(iP)) * DDENS_(iP)
      if (ke_z==1.and. iM==iP) dpresP = 0.0_RP

      del_flux(i) = 0.5_RP * (dpresP - dpresM) * nz(i)
    end do
    end do

    return
  end subroutine cal_del_flux_lin

!OCL SERIAL
  subroutine construct_pmatInv( PmatDlu, PmatDlu_ipiv, PmatL, PmatU, &
    DDENS0, POT, Rtot, CPtot_ov_CVtot, DENS_hyd, PRES_hyd,           & ! (in)
    Dz, Lift, IntrpMat_VPOrdM1, lmesh, elem,                         & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y )

    use scale_linalgebra, only: linalgebra_LU
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: PmatDlu(elem%Np,elem%Np,lmesh%NeZ)
    integer, intent(out) :: PmatDlu_ipiv(elem%Np,lmesh%NeZ)
    real(RP), intent(out) :: PmatL(elem%Np,elem%Np,lmesh%NeZ)
    real(RP), intent(out) :: PmatU(elem%Np,elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: DDENS0(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: POT(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: CPtot_ov_CVtot(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: IntrpMat_VPOrdM1(elem%Np,elem%Np)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y

    real(RP) :: DENS0(elem%Np,lmesh%NeZ)
    real(RP) :: PRES0(elem%Np,lmesh%NeZ)
    integer :: ke_z, ke_z2
    integer :: ke, p, fp, v
    real(RP) :: gamm, rgamm, rP0
    real(RP) :: dz_p(elem%Np)
    real(RP) :: PmatD(elem%Np,elem%Np)

    integer :: f1, f2, fp_s, fp_e
    integer :: FmV(elem%Nfp_v)
    integer :: FmV2 (elem%Nfp_v)
    real(RP) :: lift_op(elem%Np,elem%NfpTot)
    real(RP) :: lift_(elem%Np,elem%Np)
    real(RP) :: lift_2(elem%Np,elem%Np)
    real(RP) :: tmp(elem%Nfp_v)
    real(RP) :: fac
    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry
    rP0 = 1.0_RP / PRES00

    lift_op(:,:) = elem%Lift

    !$omp parallel do   
    do ke_z=1, lmesh%NeZ
      DENS0(:,ke_z) = DENS_hyd(:,ke_z) + DDENS0(:,ke_z)
      PRES0(:,ke_z) = PRES00 * ( Rtot(:,ke_z) * rP0 * DENS0(:,ke_z) * POT(:,ke_z) )**CPtot_ov_CVtot(:,ke_z)
    end do

!   !$omp parallel do private(ke, p, fp, v, f1, f2, ke_z2, dz_p, &
!   !$omp fac, tmp, lift_, lift_2, &
!   !$omp PmatD, FmV, FmV2, fp_s, fp_e)
    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

      !-----
      PmatD(:,:) = 0.0_RP
      PmatL(:,:,ke_z) = 0.0_RP
      PmatU(:,:,ke_z) = 0.0_RP

      do p=1, elem%Np
        dz_p(:) = lmesh%Escale(p,ke,3,3) * elem%Dx3(p,:) 
        PmatD(p,:) = dz_p(:) * CPtot_ov_CVtot(:,ke_z) * PRES0(:,ke_z) / DENS0(:,ke_z) &
                   + Grav * IntrpMat_VPOrdM1(p,:)
      end do

      do f1=1, 2
        if (f1==1) then
          f2 = 2; ke_z2 = max(ke_z-1, 1)
        else
          f2 = 1; ke_z2 = min(ke_z+1, lmesh%NeZ)
        end if
        fac  = 0.5_RP
        if ( (ke_z == 1 .and. f1==1) .or. (ke_z == lmesh%NeZ .and. f1==elem%Nfaces_v) ) then
          f2 = f1
        end if

        FmV (:) = elem%Fmask_v(:,f1)
        FmV2(:) = elem%Fmask_v(:,f2)

        fp_s = elem%Nfp_h * elem%Nfaces_h + 1 + (f1-1)*elem%Nfp_v
        fp_e = fp_s + elem%Nfp_v - 1
        
        !
        lift_ (:,:) = 0.0_RP
        lift_2(:,:) = 0.0_RP
        do fp=fp_s, fp_e
          p = fp-fp_s+1
          tmp(:) = lift_op(FmV,fp) * lmesh%Fscale(fp,ke) * nz(fp,ke_z)
          lift_ (FmV,FmV (p)) = tmp(:) * CPtot_ov_CVtot(FmV (p),ke_z ) * PRES0(FmV (p),ke_z ) / DENS0(FmV (p),ke_z )
          lift_2(FmV,FmV2(p)) = tmp(:) * CPtot_ov_CVtot(FmV2(p),ke_z2) * PRES0(FmV2(p),ke_z2) / DENS0(FmV2(p),ke_z2)
        end do

        !----        
        if ( (ke_z == 1 .and. f1==1) ) then
          PmatD(:,:) = PmatD(:,:) - fac * lift_(:,:)
        else if ( (ke_z == lmesh%NeZ .and. f1==elem%Nfaces_v) ) then
        else
          PmatD(:,:) = PmatD(:,:) - fac * lift_(:,:)
          if (f1 == 1) then
            PmatL(:,:,ke_z) = fac * lift_2(:,:)
          else
            PmatU(:,:,ke_z) = fac * lift_2(:,:)
          end if  
        end if

      end do
      call get_PmatD_LU( PmatDlu(:,:,ke_z), PmatD(:,:), PmatDlu_ipiv(:,ke_z), elem%Np )
    end do    

    return

  contains 
!OCL SERIAL
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
    
!OCL SERIAL
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
  end subroutine construct_pmatInv

end module scale_atm_dyn_dgm_hydrostatic