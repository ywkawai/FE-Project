!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Global nonhydrostatic model / HEVE
!!
!! @par Description
!!      HEVE DGM scheme for Global Atmospheric Dynamical process. 
!!      The governing equations is a fully compressible nonhydrostatic equations, 
!!      which consist of mass, momentum, and thermodynamics (density * potential temperature conservation) equations. 
!!
!! @author Yuta Kawai, Xuanzhengbo Ren, and Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,    &
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  &
    CVdry => CONST_CVdry,  &
    PRES00 => CONST_PRE00, &
    RPlanet => CONST_RADIUS

  use scale_sparsemat  
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_element_operation_base, only: ElementOperationBase3D
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    atm_dyn_dgm_nonhydro3d_common_Init,                       &
    atm_dyn_dgm_nonhydro3d_common_Final,                      &
    DENS_VID => PRGVAR_DDENS_ID, RHOT_VID => PRGVAR_DRHOT_ID, &
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID,                               &
    PRGVAR_NUM, IntrpMat_VPOrdM1
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_Init
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_Final
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_cal_tend_shallow_atm
  public :: atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_cal_tend_deep_atm

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

contains
!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_Init( mesh )
    implicit none
    class(MeshBase3D), intent(in) :: mesh
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Init( mesh )

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_Init

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_Final()
    implicit none
    !--------------------------------------------
    
    call atm_dyn_dgm_nonhydro3d_common_Final()    
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_Final

  !-------------------------------

  !OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_cal_tend_shallow_atm( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,         & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_,         & ! (in)
    DENS_hyd, PRES_hyd, PRES_hyd_ref, THERM_hyd,         & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,      & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift,   & ! (in)
    lmesh, elem, lmesh2D, elem2D )                         ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_gpu, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc_gpu
    use scale_element_operation_gpu_driver, only: &
      ElementOperationGPUDriver
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    class(ElementOperationBase3D), intent(in) :: element3D_operation
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
    real(RP), intent(in)  :: DPRES_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd_ref(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: THERM_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot (elem%Np,lmesh%NeA)    
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    type(ElementOperationGPUDriver) :: element3D_operation_driver
    real(RP) :: tend_tmp(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne)
    real(RP) :: Flux2D(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne,2)
    real(RP) :: FluxZ_store(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne)
    real(RP) :: del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP) :: drho(elem%Np,lmesh%Ne)

    integer :: IndexH2Dto3D(elem%Np)

    integer :: NeS, NeE, Nnode_h1D, Nnode_v
    !------------------------------------------------------------------------

    IndexH2Dto3D(:) = elem%IndexH2Dto3D(:)
    call element3D_operation_driver%Init( element3D_operation )
    NeS = lmesh%NeS
    NeE = lmesh%NeE
    Nnode_h1D = elem%Nnode_h1D
    Nnode_v = elem%Nnode_v

    !$acc data present(DDENS_,MOMX_,MOMY_,MOMZ_,DRHOT_,DPRES_,    &
    !$acc              DENS_hyd,PRES_hyd,THERM_hyd,               &
    !$acc              CORIOLIS,Rtot,CVtot,CPtot,DPhydDx,DPhydDy, &
    !$acc              DENS_dt,MOMX_dt,MOMY_dt,MOMZ_dt,RHOT_dt, lmesh,elem) &
    !$acc      create(del_flux,Flux2D,FluxZ_store,tend_tmp,drho) copyin(IndexH2Dto3D)

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux,                                                                & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_,                             & ! (in)
      DENS_hyd, PRES_hyd, THERM_hyd, Rtot, CVtot, CPtot,                       & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%GsqrtH, lmesh%gam, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),             & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%vmapM, lmesh%vmapP, elem%IndexH2Dto3D_bnd,                         & ! (in)
      lmesh, elem, lmesh2D, elem2D                                             ) ! (in)
    !$acc wait(1)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)
    call cal_tend_interior_gpu( &
      DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                         & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, del_flux,               & ! (in)
      DENS_hyd, THERM_hyd, DPhydDx, DPhydDy,                               & ! (in)
      lmesh%Gsqrt,lmesh%GsqrtH,lmesh%GI3(:,:,1),lmesh%GI3(:,:,2),          & ! (in)
      lmesh%Escale(:,:,1,1), lmesh%Escale(:,:,2,2), lmesh%Escale(:,:,3,3), & ! (in)
      lmesh%pos_en(:,:,1), lmesh%pos_en(:,:,2),                            & ! (in)
      Flux2D, FluxZ_store, tend_tmp, drho,                                 & ! (in)
      element3D_operation_driver, lmesh%EMap3Dto2D, IndexH2Dto3D,          & ! (in)
      lmesh,elem, lmesh%NeS, lmesh%NeE, lmesh%NeA, lmesh%Ne2DA, elem%Nnode_h1D, elem%Nnode_v ) ! (in)                     ! (in)    
    !$acc wait(1)
    call PROF_rapend('cal_dyn_tend_interior', 3)
    
    !$acc end data
    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_cal_tend_shallow_atm

  subroutine cal_tend_interior_gpu( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                  & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, del_flux,        & ! (in)
    DENS_hyd, THERM_hyd, DPhydDx, DphydDy,                        & ! (in)
    Gsqrt,GsqrtH,G13,G23, E11,E22,E33, alph, beta,                & ! (in)
    Flux2D, FluxZ_store, tend_tmp, drho,                          & ! (in)
    element3D_operation_driver,                                   & ! (in)
    EMap3Dto2D, IndexH2Dto3D, lmesh,elem, NeS, NeE, NeA, Ne2DA, Nnode_h1D, Nnode_v ) ! (in)
    use scale_const, only: &
      OHM => CONST_OHM
    use scale_element_operation_gpu_driver, only: &
      ElementOperationGPUDriver, &
      Div_kplane => ElementOperationGPU_Div_kplane, &
      DivVar5_kplane => ElementOperationGPU_DivVar5_kplane2, &
      DivVar5_z_lift => ElementOperationGPU_DivVar5_z_lift, &
      DivVar5 => ElementOperationGPU_DivVar5, &
      VFilterPM1 => ElementOperationGPU_VFilterPM1
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: NeS, NeE, NeA, Ne2DA
    integer, intent(in) :: Nnode_h1D, Nnode_v
    real(RP), intent(out) :: DENS_dt(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(out) :: MOMX_dt(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(out) :: MOMY_dt(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(out) :: MOMZ_dt(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(out) :: RHOT_dt(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: DDENS_(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: MOMX_(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: MOMY_(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: MOMZ_(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: DRHOT_(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: DPRES_(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP), intent(in)  :: DENS_hyd(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: THERM_hyd(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: DPhydDx(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: DPhydDy(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: Gsqrt(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: GsqrtH(Nnode_h1D**2,lmesh%Ne2D)
    real(RP), intent(in) :: G13(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in) :: G23(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: E11(Nnode_h1D**2,Nnode_v,lmesh%Ne)
    real(RP), intent(in)  :: E22(Nnode_h1D**2,Nnode_v,lmesh%Ne)
    real(RP), intent(in)  :: E33(Nnode_h1D**2,Nnode_v,lmesh%Ne)
    real(RP), intent(in)  :: alph(Nnode_h1D**2,lmesh%Ne2D)
    real(RP), intent(in)  :: beta(Nnode_h1D**2,lmesh%Ne2D)
    real(RP), intent(out) :: Flux2D(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne,2)
    real(RP), intent(out) :: FluxZ_store(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne)
    real(RP), intent(out) :: tend_tmp(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne)
    real(RP), intent(out) :: drho(Nnode_h1D**2,Nnode_v,lmesh%Ne)
    type(ElementOperationGPUDriver), intent(in) :: element3D_operation_driver
    integer, intent(in) :: EMap3Dto2D(lmesh%Ne)
    integer, intent(in) :: IndexH2Dto3D(elem%Nnode_h1D**2,Nnode_v)

    integer :: ke, ke2D
    integer :: ph, pz

    real(RP) :: mflxX, mflxY, mflxZ
    real(RP) :: u_, v_, w_, pt_
    real(RP) :: cor_x, cor_y
    real(RP) :: RDENS_, GsqrtV, RGsqrtV, RGsqrt
    real(RP) :: Gsqrt_, GsqrtDPRES_

    real(RP) :: G11, G12, G22
    real(RP) :: X, Y
    real(RP) :: twoOVdel2

    real(RP) :: s
    logical :: is_panel1to4
    !----------------------------

    s = 1.0_RP 
    is_panel1to4 = .true.
    if ( lmesh%panelID == 5 ) then
      is_panel1to4 = .false.
    else if ( lmesh%panelID == 6 ) then
      is_panel1to4 = .false.
      s = - 1.0_RP
    end if

    !$acc parallel loop gang collapse(2) &
    !$acc   present(DDENS_,MOMX_,MOMY_,MOMZ_,DRHOT_,DPRES_,DENS_hyd,THERM_hyd,Flux2D,FluxZ_store, EMap3Dto2D) async(1)
    do ke = NeS, NeE
    do pz=1, Nnode_v
      ke2d = EMap3Dto2D(ke)
      !$acc loop vector
      do ph=1, Nnode_h1D**2
        Gsqrt_ = Gsqrt(ph,pz,ke)
        GsqrtV  = Gsqrt_ / GsqrtH(ph,ke2d)
        RGsqrtV = 1.0_RP / GsqrtV
        RGsqrt  = 1.0_RP / Gsqrt_
        RDENS_  = 1.0_RP / (DDENS_(ph,pz,ke) + DENS_hyd(ph,pz,ke))

        G11 = lmesh%GIJ(IndexH2Dto3D(ph,pz),ke2d,1,1)
        G12 = lmesh%GIJ(IndexH2Dto3D(ph,pz),ke2d,1,2)
        G22 = lmesh%GIJ(IndexH2Dto3D(ph,pz),ke2d,2,2)

        !-
        mflxX = Gsqrt_ * MOMX_(ph,pz,ke)
        mflxY = Gsqrt_ * MOMY_(ph,pz,ke)
        mflxZ = Gsqrt_ * ( &
            MOMZ_(ph,pz,ke) * RGsqrtV        &
          + G13(ph,pz,ke) * MOMX_(ph,pz,ke) &
          + G23(ph,pz,ke) * MOMY_(ph,pz,ke) )

        Flux2D(ph,DENS_VID,pz,ke,1) = mflxX
        Flux2D(ph,DENS_VID,pz,ke,2) = mflxY
        FluxZ_store(ph,DENS_VID,pz,ke) = mflxZ
        
        !-
        pt_ = ( THERM_hyd(ph,pz,ke) + DRHOT_(ph,pz,ke) ) * RDENS_

        Flux2D(ph,RHOT_VID,pz,ke,1) = mflxX * pt_
        Flux2D(ph,RHOT_VID,pz,ke,2) = mflxY * pt_
        FluxZ_store(ph,RHOT_VID,pz,ke) = mflxZ * pt_

        !-
        GsqrtDPRES_ = Gsqrt_ * DPRES_(ph,pz,ke)
        w_ = MOMZ_(ph,pz,ke) * RDENS_
        Flux2D(ph,MOMZ_VID,pz,ke,1) = mflxX * w_
        Flux2D(ph,MOMZ_VID,pz,ke,2) = mflxY * w_
        FluxZ_store(ph,MOMZ_VID,pz,ke) = mflxZ * w_ + GsqrtDPRES_ * RGsqrtV

        !-        
        u_ = MOMX_(ph,pz,ke) * RDENS_
        Flux2D(ph,MOMX_VID,pz,ke,1) = mflxX * u_ + G11 * GsqrtDPRES_
        Flux2D(ph,MOMX_VID,pz,ke,2) = mflxY * u_ + G12 * GsqrtDPRES_
        FluxZ_store(ph,MOMX_VID,pz,ke) = mflxZ * u_ + GsqrtDPRES_ * ( G11 * G13(ph,pz,ke) + G12 * G23(ph,pz,ke) )

        v_ = MOMY_(ph,pz,ke) * RDENS_
        Flux2D(ph,MOMY_VID,pz,ke,1) = mflxX * v_ + G12 * GsqrtDPRES_
        Flux2D(ph,MOMY_VID,pz,ke,2) = mflxY * v_ + G22 * GsqrtDPRES_
        FluxZ_store(ph,MOMY_VID,pz,ke) = mflxZ * v_ + GsqrtDPRES_ * ( G12 * G13(ph,pz,ke) + G22 * G23(ph,pz,ke) )
      end do
    end do
    end do

    call DivVar5( element3D_operation_driver, Flux2D(:,:,:,:,1), Flux2D(:,:,:,:,2), FluxZ_store, del_flux, &
      E11, E22, E33, Gsqrt, lmesh%Ne, tend_tmp )
    
    call VFilterPM1( element3D_operation_driver, DDENS_, lmesh%NeA, lmesh%Ne, &
      drho )


    !$acc parallel loop gang collapse(2) &
    !$acc   present(DENS_dt,MOMX_dt,MOMY_dt,MOMZ_dt,RHOT_dt, DPhydDx,DPhydDy, tend_tmp, drho, EMap3Dto2D) async(1)
    do ke = NeS, NeE
    do pz=1, Nnode_v
      ke2d = EMap3Dto2D(ke)
      !$acc loop vector
      do ph=1, Nnode_h1D**2
        DENS_dt(ph,pz,ke) = - tend_tmp(ph,DENS_VID,pz,ke)
        RHOT_dt(ph,pz,ke) = - tend_tmp(ph,RHOT_VID,pz,ke)


        MOMZ_dt(ph,pz,ke) = - tend_tmp(ph,MOMZ_VID,pz,ke) &
                          - Grav * drho(ph,pz,ke)
      end do
      !$acc loop vector
      do ph=1, Nnode_h1D**2
        !-
        X = tan(alph(ph,ke2d))
        Y = tan(beta(ph,ke2d))
        twoOVdel2 = 2.0_RP / ( 1.0_RP + X**2 + Y**2 )

        cor_x = s * OHM * twoOVdel2 * ( - X * Y             * MOMX_(ph,pz,ke) + ( 1.0_RP + Y**2 ) * MOMY_(ph,pz,ke) )
        cor_y = s * OHM * twoOVdel2 * ( - ( 1.0_RP + X**2 ) * MOMX_(ph,pz,ke) + X * Y             * MOMY_(ph,pz,ke) )
        if ( is_panel1to4 ) then
          cor_x = s * Y * cor_x
          cor_y = s * Y * cor_y
        end if

        MOMX_dt(ph,pz,ke) = - tend_tmp(ph,MOMX_VID,pz,ke) &
                          - DPhydDx(ph,pz,ke) &
                          + cor_x
        MOMY_dt(ph,pz,ke) = - tend_tmp(ph,MOMY_VID,pz,ke) &
                          - DPhydDy(ph,pz,ke) &
                          - cor_y
      end do
    end do
    end do
    !$acc end parallel
    return
  end subroutine cal_tend_interior_gpu

!OCL SERIAL
  subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_cal_tend_deep_atm( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,         & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_,         & ! (in)
    DENS_hyd, PRES_hyd, PRES_hyd_ref, THERM_hyd,         & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,      & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift,   & ! (in)
    lmesh, elem, lmesh2D, elem2D )                         ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalhvc_asis
    use scale_const, only: &
      OHM => CONST_OHM
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    class(ElementOperationBase3D), intent(in) :: element3D_operation
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
    real(RP), intent(in)  :: DPRES_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd_ref(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: THERM_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in)  :: Rtot (elem%Np,lmesh%NeA)    
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DPRES_hyd(elem%Np), GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: RHOT_(elem%Np)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%np), drho(elem%Np)

    real(RP) :: G11(elem%Np), G12(elem%Np), G22(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np), Rgam2(elem%Np)
    real(RP) :: X2D(elem2D%Np,lmesh2D%Ne), Y2D(elem2D%Np,lmesh2D%Ne)
    real(RP) :: X(elem%Np), Y(elem%Np), twoOVdel2(elem%Np)
    real(RP) :: OM1(elem%Np), OM2(elem%Np), OM3(elem%Np), DEL(elem%Np), R(elem%Np)
    logical :: is_panel1to4
    real(RP) :: s

    integer :: ke, ke2d
    integer :: p
    
    real(RP) :: rgamm    
    real(RP) :: rP0
    real(RP) :: P0ovR
    !------------------------------------------------------------------------

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux, del_flux_hyd,                                                  & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd,         & ! (in)
      Rtot, CVtot, CPtot,                                                      & ! (in)
      lmesh%Gsqrt, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), & ! (in)
      lmesh%GsqrtH, lmesh%gam, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),             & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),  & ! (in)
      lmesh%vmapM, lmesh%vmapP, elem%IndexH2Dto3D_bnd,                         & ! (in)
      lmesh, elem, lmesh2D, elem2D                                             ) ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    P0ovR = PRES00 / Rdry

    s = 2.0_RP * OHM
    is_panel1to4 = .true.
    if ( lmesh%panelID == 5 ) then
      is_panel1to4 = .false.
    else if ( lmesh%panelID == 6 ) then
      is_panel1to4 = .false.
      s = - s
    end if

    !$omp parallel private(                           &
    !$omp RHOT_, rdens_, u_, v_, w_, wt_,             &
    !$omp Fx, Fy, Fz, LiftDelFlx,                     &
    !$omp drho, DPRES_hyd, GradPhyd_x, GradPhyd_y,    &
    !$omp G11, G12, G22, Rgam2, GsqrtV, RGsqrtV,      &
    !$omp X, Y, twoOVdel2,                            &
    !$omp OM1, OM2, OM3, DEL, R, ke, ke2D             )

    !$omp do
    do ke2D = lmesh2D%NeS, lmesh2D%NeE
      X2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,1))
      Y2D(:,ke2d) = tan(lmesh2D%pos_en(:,ke2d,2))
    end do

    !$omp do
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      Rgam2(:) = 1.0_RP / lmesh%gam(:,ke)**2
      G11(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,1) * Rgam2(:)
      G12(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,1,2) * Rgam2(:)
      G22(:) = lmesh%GIJ(elem%IndexH2Dto3D,ke2d,2,2) * Rgam2(:)
      GsqrtV(:)  = lmesh%Gsqrt(:,ke) * Rgam2(:) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2d)
      RGsqrtV(:) = 1.0_RP / GsqrtV(:)

      !--
      RHOT_(:) = P0ovR * ( PRES_hyd(:,ke) * rP0 )**rgamm + DRHOT_(:,ke)
      ! DPRES_(:) = PRES00 * ( Rtot(:,ke) * rP0 * RHOT_(:) )**( CPtot(:,ke) / CVtot(:,ke) ) &
      !           - PRES_hyd(:,ke)

      rdens_(:) = 1.0_RP / ( DDENS_(:,ke) + DENS_hyd(:,ke) )
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)
      wt_(:) = w_(:) * RGsqrtV(:) + lmesh%GI3(:,ke,1) * u_(:) + lmesh%GI3(:,ke,2) * v_(:) 

      X(:) = X2D(elem%IndexH2Dto3D,ke2d)
      Y(:) = Y2D(elem%IndexH2Dto3D,ke2d)
      DEL(:) = sqrt( 1.0_RP + X(:)**2 + Y(:)**2 )
      twoOVdel2(:) = 2.0_RP / ( 1.0_RP + X(:)**2 + Y(:)**2 )

      R(:) = RPlanet * lmesh%gam(:,ke)

      !-  pnl=1~4: OM1:                     0, OM2: s del / (r (1+Y^2)) ,   s OM3 : Y /del
      !-  pnl=5,6: OM1: - s X del/(r (1+X^2)), OM2: -  s Y del/(r(1+Y^2)), OM3 : s/del
      if ( is_panel1to4 ) then
        OM1(:) = 0.0_RP
        OM2(:) = s * DEL(:) / ( R(:) * ( 1.0_RP + Y(:)**2 ) )
        OM3(:) = s * Y(:) / DEL(:)
      else
        OM1(:) = - s * X(:) * DEL(:) / ( R(:) * ( 1.0_RP + X(:)**2 ) )
        OM2(:) = - s * Y(:) * DEL(:) / ( R(:) * ( 1.0_RP + Y(:)**2 ) )
        OM3(:) =   s / DEL(:)
      end if

      drho(:) = matmul(IntrpMat_VPOrdM1, DDENS_(:,ke))

      !-- Gradient hydrostatic pressure
      
      DPRES_hyd(:) = PRES_hyd(:,ke) - PRES_hyd_ref(:,ke)

      call sparsemat_matmul(Dx, GsqrtV(:) * DPRES_hyd(:), Fx)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,1) * DPRES_hyd(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,1), LiftDelFlx)
      GradPhyd_x(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      call sparsemat_matmul(Dy, GsqrtV(:) * DPRES_hyd(:), Fy)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,2) * DPRES_hyd(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,2), LiftDelFlx)
      GradPhyd_y(:) = lmesh%Escale(:,ke,2,2) * Fy(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)
      
      !-- DENS
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * MOMY_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( DDENS_(:,ke) + DENS_hyd(:,ke) ) * wt_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DENS_VID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)     &
          + lmesh%Escale(:,ke,2,2) * Fy(:)     &
          + lmesh%Escale(:,ke,3,3) * Fz(:)     &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      !-- MOMX
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMX_(:,ke) + G11(:) * DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMX_(:,ke) + G12(:) * DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMX_(:,ke)                              &
                                                    + ( lmesh%GI3(:,ke,1) * G11(:) + lmesh%GI3(:,ke,2) * G12(:) ) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                                &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                                &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                                &
            + LiftDelFlx(:)                   ) / lmesh%Gsqrt(:,ke)                         &
          - twoOVdel2(:) * Y(:) *                                                           & !-> metric terms
            ( X(:) * Y(:) * u_(:) - ( 1.0_RP + Y(:)**2 ) * v_(:) ) * MOMX_(:,ke)            & !
          - 2.0_RP * u_(:) * MOMZ_(:,ke) / R(:)                                             & !<-
          - ( G11(:) * GradPhyd_x(:) + G12(:) * GradPhyd_y(:) ) * RGsqrtV(:)                & !-> gradient hydrostatic pressure
          - lmesh%Gsqrt(:,ke) * (  G11(:) * ( OM2(:) * MOMZ_(:,ke) - OM3(:) * MOMY_(:,ke) ) & !-> Coriolis term
                                 - G12(:) * ( OM1(:) * MOMZ_(:,ke) - OM3(:) * MOMX_(:,ke) ) )

      !-- MOMY
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * ( u_ (:) * MOMY_(:,ke) + G12(:) * DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * ( v_ (:) * MOMY_(:,ke) + G22(:) * DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMY_(:,ke)                              &
                                                    + ( lmesh%GI3(:,ke,1) * G12(:) + lmesh%GI3(:,ke,2) * G22(:) ) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = &
            - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                                &
              + lmesh%Escale(:,ke,2,2) * Fy(:)                                                &
              + lmesh%Escale(:,ke,3,3) * Fz(:)                                                &
              + LiftDelFlx(:)                  ) / lmesh%Gsqrt(:,ke)                          &
            - twoOVdel2(:) * X(:) *                                                           & !-> metric terms
              ( - (1.0_RP + X(:)**2) * u_(:) + X(:) * Y(:) * v_(:) ) * MOMY_(:,ke)            & !
            - 2.0_RP * v_(:) * MOMZ_(:,ke) / R(:)                                             & !<-
            - ( G12(:) * GradPhyd_x(:) + G22(:) * GradPhyd_y(:) ) * RGsqrtV(:)                & !-> gradient hydrostatic pressure
            - lmesh%Gsqrt(:,ke) * (  G12(:) * ( OM2(:) * MOMZ_(:,ke) - OM3(:) * MOMY_(:,ke) ) & !-> Coriolis term
                                   - G22(:) * ( OM1(:) * MOMZ_(:,ke) - OM3(:) * MOMX_(:,ke) ) )

      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *   u_ (:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *   v_ (:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMZ_(:,ke) + RGsqrtV(:) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)                                                                    &
            + lmesh%Escale(:,ke,2,2) * Fy(:)                                                                    &
            + lmesh%Escale(:,ke,3,3) * Fz(:)                                                                    &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)                                                               &
          - 0.25_RP * R(:) * twoOVdel2(:)**2 * ( 1.0_RP * X(:)**2 ) * ( 1.0_RP * Y(:)**2 )                      & !-> metric terms
            * ( - ( 1.0_RP + X(:)**2 ) * MOMX_(:,ke) * u_(:)                                                    & !
                + 2.0_RP * X(:) * Y(:) * MOMX_(:,ke) * v_(:)                                                    & !
                - ( 1.0_RP + Y(:)**2 ) * MOMY_(:,ke) * v_(:)  )                                                 & !<-
          + 2.0_RP * DPRES_(:,ke) / R(:)                                                                        & !-> metric term with gradient of pressure deviaition
          - lmesh%Gsqrt(:,ke) * ( OM1(:) * MOMY_(:,ke) - OM2(:) * MOMX_(:,ke) )                                 & !-> Coriolis term
          - Grav * Rgam2(:) * drho(:)                                                                             !-> buoyancy term

      !-- RHOT
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * u_ (:) * RHOT_(:), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * v_ (:) * RHOT_(:), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * wt_(:) * RHOT_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,RHOT_VID), LiftDelFlx)
      
      RHOT_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 
    end do
    !$omp end do
    !$omp end parallel
    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu_cal_tend_deep_atm

end module scale_atm_dyn_dgm_globalnonhydro3d_rhot_heve_gpu
