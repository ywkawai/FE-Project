!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Regional nonhydrostatic model / HEVE
!!
!! @par Description
!!      HEVE DGM scheme for Atmospheric dynamical process which runs on GPU. 
!!      The governing equations is a fully compressible nonhydrostatic equations, 
!!      which consist of mass, momentum, and thermodynamics (density * potential temperature conservation) equations. 
!!
!! @author Yuta Kawai, Xuanzhengbo Ren, and Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_rhot_heve_gpu
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
  use scale_element_operation_base, only: ElementOperationBase3D
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
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
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_gpu_Init
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_gpu_Final
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_gpu
  public :: atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_gpu2

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
  !> Initialize the module for HEVE scheme on GPU
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_gpu_Init( mesh )
    implicit none
    class(MeshBase3D), intent(in) :: mesh
    !--------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_Init( mesh )
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_gpu_Init

  !> Finalize the module for HEVE scheme on GPU
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_gpu_Final()
    implicit none
    !--------------------------------------------
    
    call atm_dyn_dgm_nonhydro3d_common_Final()
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_gpu_Final

  !-------------------------------

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_gpu( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,         & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_,         & ! (in)
    DENS_hyd, PRES_hyd, PRES_hyd_ref, THERM_hyd,         & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,      & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift,   & ! (in)
    lmesh, elem, lmesh2D, elem2D )                         ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_gpu, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_gpu
    use scale_element_operation_gpu_driver, only: &
      ElementOperationGPUDriver, &
      Div_kplane => ElementOperationGPU_Div_kplane, &
      DivVar5_kplane => ElementOperationGPU_DivVar5_kplane, &
      DivVar5_z_lift => ElementOperationGPU_DivVar5_z_lift, &
      VFilterPM1 => ElementOperationGPU_VFilterPM1
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
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    real(RP) :: Flux2D(elem%Nnode_h1D**2,5,2,lmesh%Ne)
    real(RP) :: Flux2D_(elem%Nnode_h1D**2,2,elem%Nnode_v,lmesh%Ne)   
    real(RP) :: mflxX, mflxY, mflxZ
    real(RP) :: FluxZ_store(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne)
    real(RP) :: del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP) :: u_, v_, w_, pt_
    real(RP) :: cor
    real(RP) :: drho(elem%Np,lmesh%Ne)
    real(RP) :: RDENS_, GsqrtH, GsqrtV, RGsqrtV, RGsqrt
    real(RP) :: Gsqrt_, GsqrtDPRES_, E11, E22, E33

    integer :: ke, ke2d
    integer :: p, ph, pz
    integer :: iv
    
    integer :: IndexH2Dto3D(elem%Np)
    integer :: EMap3Dto2D(lmesh%Ne)

    type(ElementOperationGPUDriver) :: element3D_operation_driver
    real(RP) :: tend_tmp(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne)

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
    !$acc      create(del_flux,Flux2D,Flux2D_,FluxZ_store,tend_tmp,drho) copyin(IndexH2Dto3D)

    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux,                                                                   & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, THERM_hyd, & ! (in)
      Rtot, CVtot, CPtot,                                                         & ! (in)
      lmesh%Gsqrt, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),                            & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),     & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                                   & ! (in)
      lmesh, elem, lmesh2D, elem2D )                                                ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)

    call PROF_rapstart('cal_dyn_tend_interior', 3)
    !$acc parallel loop gang collapse(2) present(DDENS_,MOMX_,MOMY_,MOMZ_,DRHOT_,DPRES_,DENS_hyd,THERM_hyd,lmesh%Gsqrt,lmesh%GsqrtH,lmesh%Escale)
    do ke = NeS, NeE
      do pz=1, Nnode_v
        ke2d = lmesh%EMap3Dto2D(ke)

        !$acc loop vector
        do ph=1, Nnode_h1D**2
          p = ph + (pz-1) * Nnode_h1D**2

          Gsqrt_ = lmesh%Gsqrt(p,ke)
          GsqrtV  = Gsqrt_ / lmesh%GsqrtH(IndexH2Dto3D(p),ke2d)
          RGsqrtV = 1.0_RP / GsqrtV
          RGsqrt  = 1.0_RP / Gsqrt_
          
          Flux2D_(ph,1,pz,ke) = Gsqrt_ * MOMX_(p,ke)
          Flux2D_(ph,2,pz,ke) = Gsqrt_ * MOMY_(p,ke)
          FluxZ_store(ph,DENS_VID,pz,ke) = Gsqrt_ * ( &
              MOMZ_(p,ke) * RGsqrtV        &
            + lmesh%GI3(p,ke,1) * MOMX_(p,ke) &
            + lmesh%GI3(p,ke,2) * MOMY_(p,ke) )
        end do
        call Div_kplane( element3D_operation_driver, Flux2D_(:,:,pz,ke), lmesh%Escale(:,ke,1,1), lmesh%Escale(:,ke,2,2), pz, &
          tend_tmp(:,DENS_VID,pz,ke) )
  
        !$acc loop vector
        do ph=1, Nnode_h1D**2
          p = ph + (pz-1) * Nnode_h1D**2

          Gsqrt_ = lmesh%Gsqrt(p,ke)
          RDENS_  = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))
          pt_ = ( THERM_hyd(p,ke) + DRHOT_(p,ke) ) * RDENS_

          Flux2D_(ph,1,pz,ke) = Gsqrt_ * MOMX_(p,ke) * pt_
          Flux2D_(ph,2,pz,ke) = Gsqrt_ * MOMY_(p,ke) * pt_
          FluxZ_store(ph,RHOT_VID,pz,ke)  = FluxZ_store(ph,DENS_VID,pz,ke) * pt_
        end do
        call Div_kplane( element3D_operation_driver, Flux2D_(:,:,pz,ke), lmesh%Escale(:,ke,1,1), lmesh%Escale(:,ke,2,2), pz, &
          tend_tmp(:,RHOT_VID,pz,ke) )
        
        !$acc loop vector
        do ph=1, Nnode_h1D**2
          p = ph + (pz-1) * Nnode_h1D**2

          Gsqrt_ = lmesh%Gsqrt(p,ke)
          GsqrtH  = lmesh%GsqrtH(IndexH2Dto3D(p),ke2d)
          RDENS_  = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))
          w_ = MOMZ_(p,ke) * RDENS_

          Flux2D_(ph,1,pz,ke) = Gsqrt_ * MOMX_(p,ke) * w_
          Flux2D_(ph,2,pz,ke) = Gsqrt_ * MOMY_(p,ke) * w_
          FluxZ_store(ph,MOMZ_VID,pz,ke) = FluxZ_store(ph,DENS_VID,pz,ke) * w_ + GsqrtH * DPRES_(p,ke)
        end do
        call Div_kplane( element3D_operation_driver, Flux2D_(:,:,pz,ke), lmesh%Escale(:,ke,1,1), lmesh%Escale(:,ke,2,2), pz, &
          tend_tmp(:,MOMZ_VID,pz,ke) )  

        !$acc loop vector
        do ph=1, Nnode_h1D**2
          p = ph + (pz-1) * Nnode_h1D**2

          Gsqrt_ = lmesh%Gsqrt(p,ke)
          GsqrtH  = lmesh%GsqrtH(IndexH2Dto3D(p),ke2d)
          RDENS_  = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))

          u_ = MOMX_(p,ke) * RDENS_
          GsqrtDPRES_ = Gsqrt_ * DPRES_(p,ke)

          Flux2D_(ph,1,pz,ke) = Gsqrt_ * MOMX_(p,ke) * u_ + GsqrtDPRES_
          Flux2D_(ph,2,pz,ke) = Gsqrt_ * MOMY_(p,ke) * u_
          FluxZ_store(ph,MOMX_VID,pz,ke) = FluxZ_store(ph,DENS_VID,pz,ke) * u_ + GsqrtDPRES_ * lmesh%GI3(p,ke,1)
        end do
        call Div_kplane( element3D_operation_driver, Flux2D_(:,:,pz,ke), lmesh%Escale(:,ke,1,1), lmesh%Escale(:,ke,2,2), pz, &
          tend_tmp(:,MOMX_VID,pz,ke) )  


        !$acc loop vector
        do ph=1, Nnode_h1D**2
          p = ph + (pz-1) * Nnode_h1D**2

          Gsqrt_ = lmesh%Gsqrt(p,ke)
          GsqrtH  = lmesh%GsqrtH(IndexH2Dto3D(p),ke2d)
          RDENS_  = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))
          
          v_ = MOMY_(p,ke) * RDENS_
          GsqrtDPRES_ = Gsqrt_ * DPRES_(p,ke)

          Flux2D_(ph,1,pz,ke) = Gsqrt_ * MOMX_(p,ke) * v_ + GsqrtDPRES_
          Flux2D_(ph,2,pz,ke) = Gsqrt_ * MOMY_(p,ke) * v_
          FluxZ_store(ph,MOMY_VID,pz,ke) = FluxZ_store(ph,DENS_VID,pz,ke) * v_ + GsqrtDPRES_ * lmesh%GI3(p,ke,2)
        end do
        call Div_kplane( element3D_operation_driver, Flux2D_(:,:,pz,ke), lmesh%Escale(:,ke,1,1), lmesh%Escale(:,ke,2,2), pz, &
          tend_tmp(:,MOMY_VID,pz,ke) )  

        ! !$acc loop vector
        ! do ph=1, Nnode_h1D**2
        !   p = ph + (pz-1) * Nnode_h1D**2

        !   Gsqrt_ = lmesh%Gsqrt(p,ke)
        !   GsqrtV  = Gsqrt_ / lmesh%GsqrtH(IndexH2Dto3D(p),ke2d)
        !   RGsqrtV = 1.0_RP / GsqrtV
        !   RGsqrt  = 1.0_RP / Gsqrt_
        !   RDENS_  = 1.0_RP / (DDENS_(p,ke) + DENS_hyd(p,ke))

        !   !-
        !   mflxX = Gsqrt_ * MOMX_(p,ke)
        !   mflxY = Gsqrt_ * MOMY_(p,ke)
        !   mflxZ = Gsqrt_ * ( &
        !       MOMZ_(p,ke) * RGsqrtV        &
        !     + lmesh%GI3(p,ke,1) * MOMX_(p,ke) &
        !     + lmesh%GI3(p,ke,2) * MOMY_(p,ke) )

        !   Flux2D(ph,DENS_VID,1,ke) = mflxX
        !   Flux2D(ph,DENS_VID,2,ke) = mflxY
        !   FluxZ_store(ph,DENS_VID,pz,ke) = mflxZ
          
        !   !-
        !   pt_ = ( THERM_hyd(p,ke) + DRHOT_(p,ke) ) * RDENS_

        !   Flux2D(ph,RHOT_VID,1,ke) = mflxX * pt_
        !   Flux2D(ph,RHOT_VID,2,ke) = mflxY * pt_
        !   FluxZ_store(ph,RHOT_VID,pz,ke) = mflxZ * pt_

        !   !-
        !   GsqrtDPRES_ = Gsqrt_ * DPRES_(p,ke)
        !   w_ = MOMZ_(p,ke) * RDENS_
        !   Flux2D(ph,MOMZ_VID,1,ke) = mflxX * w_
        !   Flux2D(ph,MOMZ_VID,2,ke) = mflxY * w_
        !   FluxZ_store(ph,MOMZ_VID,pz,ke) = mflxZ * w_ + GsqrtDPRES_ * RGsqrtV

        !   !-
        !   u_ = MOMX_(p,ke) * RDENS_

        !   Flux2D(ph,MOMX_VID,1,ke) = mflxX * u_ + GsqrtDPRES_
        !   Flux2D(ph,MOMX_VID,2,ke) = mflxY * u_ 
        !   FluxZ_store(ph,MOMX_VID,pz,ke) = mflxZ * u_ + GsqrtDPRES_ * lmesh%GI3(p,ke,1)

        !   v_ = MOMY_(p,ke) * RDENS_
        !   Flux2D(ph,MOMY_VID,1,ke) = mflxX * v_
        !   Flux2D(ph,MOMY_VID,2,ke) = mflxY * v_ + GsqrtDPRES_ 
        !   FluxZ_store(ph,MOMY_VID,pz,ke) = mflxZ * v_ + GsqrtDPRES_ * lmesh%GI3(p,ke,2)
        ! end do
        ! call DivVar5_kplane( element3D_operation_driver,              &
        !   Flux2D(:,:,:,ke), lmesh%Escale(:,ke,1,1), lmesh%Escale(:,ke,2,2), pz, &
        !   tend_tmp(:,:,pz,ke) )
      end do
    end do
    !$acc parallel loop gang
    do ke = NeS, NeE
      call DivVar5_z_lift( element3D_operation_driver, &
        FluxZ_store(:,:,:,ke),  del_flux(:,:,ke), lmesh%Escale(:,ke,3,3), lmesh%Gsqrt(:,ke), &
        tend_tmp(:,:,:,ke) )
    end do

    call VFilterPM1( element3D_operation_driver, DDENS_, lmesh%NeA, lmesh%Ne, &
      drho )

    !$acc parallel loop gang collapse(2)
    do ke = NeS, NeE
      !--
      do pz=1, Nnode_v
         ke2d = lmesh%EMap3Dto2D(ke)

        !$acc loop vector
        do ph=1, Nnode_h1D**2
          p = ph + (pz-1) * Nnode_h1D**2

          DENS_dt(p,ke) = - tend_tmp(ph,DENS_VID,pz,ke)
          RHOT_dt(p,ke) = - tend_tmp(ph,RHOT_VID,pz,ke)


          MOMZ_dt(p,ke) = - tend_tmp(ph,MOMZ_VID,pz,ke) &
                          - Grav * drho(p,ke)
          
          cor = CORIOLIS(ph,ke2d)
          MOMX_dt(p,ke) = - tend_tmp(ph,MOMX_VID,pz,ke) &
                          - DPhydDx(p,ke) &
                          + cor * MOMY_(p,ke)
          MOMY_dt(p,ke) = - tend_tmp(ph,MOMY_VID,pz,ke) &
                          - DPhydDy(p,ke) &
                          - cor * MOMX_(p,ke)
        end do
      end do
    end do
    call PROF_rapend('cal_dyn_tend_interior', 3)
    !$acc end data

    call element3D_operation_driver%Final()
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_gpu

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_gpu2( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,         & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_,         & ! (in)
    DENS_hyd, PRES_hyd, PRES_hyd_ref, THERM_hyd,         & ! (in)
    CORIOLIS, Rtot, CVtot, CPtot, DPhydDx, DPhydDy,      & ! (in)
    element3D_operation, Dx, Dy, Dz, Sx, Sy, Sz, Lift,   & ! (in)
    lmesh, elem, lmesh2D, elem2D )                         ! (in)

    use scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_gpu, only: &
      get_ebnd_flux => atm_dyn_dgm_nonhydro3d_rhot_heve_numflux_get_generalvc_gpu
    use scale_element_operation_gpu_driver, only: &
      ElementOperationGPUDriver, &
      Div_kplane => ElementOperationGPU_Div_kplane, &
      DivVar5_kplane => ElementOperationGPU_DivVar5_kplane, &
      DivVar5_z_lift => ElementOperationGPU_DivVar5_z_lift, &
      VFilterPM1 => ElementOperationGPU_VFilterPM1
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
    real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DPhydDy(elem%Np,lmesh%NeA)

    real(RP) :: Flux2D(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne,2)
    real(RP) :: mflxX, mflxY, mflxZ
    real(RP) :: FluxZ_store(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne)
    real(RP) :: del_flux(elem%NfpTot,PRGVAR_NUM,lmesh%Ne)
    real(RP) :: u_, v_, w_, pt_
    real(RP) :: cor
    real(RP) :: drho(elem%Np,lmesh%Ne)
    real(RP) :: RDENS_, GsqrtH, GsqrtV, RGsqrtV, RGsqrt
    real(RP) :: Gsqrt_, GsqrtDPRES_

    integer :: ke, ke2d
    integer :: p, ph, pz
    integer :: iv
    
    integer :: IndexH2Dto3D(elem%Np)
    integer :: EMap3Dto2D(lmesh%Ne)

    type(ElementOperationGPUDriver) :: element3D_operation_driver
    real(RP) :: tend_tmp(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne)

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
      del_flux,                                                                   & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, THERM_hyd, & ! (in)
      Rtot, CVtot, CPtot,                                                         & ! (in)
      lmesh%Gsqrt, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),                            & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3),     & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                                   & ! (in)
      lmesh, elem, lmesh2D, elem2D )
    !!$acc wait(1)                                                      ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)

    call PROF_rapstart('cal_dyn_tend_interior', 3)
    call atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_interior_gpu( &
      DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                  & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, del_flux,        & ! (in)
      DENS_hyd, THERM_hyd, CORIOLIS, DPhydDx, DPhydDy,              & ! (in)
      lmesh%Gsqrt,lmesh%GsqrtH,lmesh%GI3(:,:,1),lmesh%GI3(:,:,2),          &
      lmesh%Escale(:,:,1,1), lmesh%Escale(:,:,2,2), lmesh%Escale(:,:,3,3), & ! (in)
      Flux2D, FluxZ_store, tend_tmp, drho, & ! (in)
      element3D_operation_driver, & ! (in)
      lmesh%EMap3Dto2D, lmesh,elem, lmesh%NeS, lmesh%NeE, lmesh%NeA, lmesh%Ne2DA, elem%Nnode_h1D, elem%Nnode_v )                         ! (in)    
    !$acc wait(1)
    call PROF_rapend('cal_dyn_tend_interior', 3)
    !$acc end data

    call element3D_operation_driver%Final()
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_gpu2

  subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_interior_gpu( &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,                  & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, del_flux,        & ! (in)
    DENS_hyd, THERM_hyd, Coriolis, DPhydDx, DphydDy,              & ! (in)
    Gsqrt,GsqrtH,G13,G23, E11,E22,E33, & ! (in)
    Flux2D, FluxZ_store, tend_tmp, drho, & ! (in)
    element3D_operation_driver, & ! (in)
    EMap3Dto2D, lmesh,elem, NeS, NeE, NeA, Ne2DA, Nnode_h1D, Nnode_v )                         ! (in)
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
    real(RP), intent(in)  :: Coriolis(Nnode_h1D**2,Ne2DA)
    real(RP), intent(in)  :: DPhydDx(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: DPhydDy(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: Gsqrt(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: GsqrtH(Nnode_h1D**2,lmesh%Ne2D)
    real(RP), intent(in) :: G13(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in) :: G23(Nnode_h1D**2,Nnode_v,NeA)
    real(RP), intent(in)  :: E11(Nnode_h1D**2,Nnode_v,lmesh%Ne)
    real(RP), intent(in)  :: E22(Nnode_h1D**2,Nnode_v,lmesh%Ne)
    real(RP), intent(in)  :: E33(Nnode_h1D**2,Nnode_v,lmesh%Ne)
    real(RP), intent(out) :: Flux2D(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne,2)
    real(RP), intent(out) :: FluxZ_store(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne)
    real(RP), intent(out) :: tend_tmp(elem%Nnode_h1D**2,5,elem%Nnode_v,lmesh%Ne)
    real(RP), intent(out) :: drho(Nnode_h1D**2,Nnode_v,lmesh%Ne)
    type(ElementOperationGPUDriver), intent(in) :: element3D_operation_driver
    integer, intent(in) :: EMap3Dto2D(lmesh%Ne)

    integer :: ke, ke2D, ph,pz

    real(RP) :: mflxX, mflxY, mflxZ
    real(RP) :: u_, v_, w_, pt_
    real(RP) :: cor
    real(RP) :: RDENS_, GsqrtV, RGsqrtV, RGsqrt
    real(RP) :: Gsqrt_, GsqrtDPRES_

    !----------------------------
    !!$acc parallel present(DENS_dt,MOMX_dt,MOMY_dt,MOMZ_dt,RHOT_dt,&
    !!$acc DDENS_,MOMX_,MOMY_,MOMZ_,DRHOT_,DPRES_,DENS_hyd,THERM_hyd,&
    !!$acc DPhydDx,DPhydDy,Coriolis, del_flux, tend_tmp,drho, Gsqrt,GsqrtH,G13,G23,E11,E22,E33,Flux2D,FluxZ_store, EMap3Dto2D)

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
        Flux2D(ph,MOMX_VID,pz,ke,1) = mflxX * u_ + GsqrtDPRES_
        Flux2D(ph,MOMX_VID,pz,ke,2) = mflxY * u_ 
        FluxZ_store(ph,MOMX_VID,pz,ke) = mflxZ * u_ + GsqrtDPRES_ * G13(ph,pz,ke)

        v_ = MOMY_(ph,pz,ke) * RDENS_
        Flux2D(ph,MOMY_VID,pz,ke,1) = mflxX * v_
        Flux2D(ph,MOMY_VID,pz,ke,2) = mflxY * v_ + GsqrtDPRES_ 
        FluxZ_store(ph,MOMY_VID,pz,ke) = mflxZ * v_ + GsqrtDPRES_ * G23(ph,pz,ke)
      end do
    end do
    end do

    call DivVar5( element3D_operation_driver, Flux2D(:,:,:,:,1), Flux2D(:,:,:,:,2), FluxZ_store, del_flux, &
      E11, E22, E33, Gsqrt, lmesh%Ne, tend_tmp )
    
    call VFilterPM1( element3D_operation_driver, DDENS_, lmesh%NeA, lmesh%Ne, &
      drho )


    !$acc parallel loop gang collapse(2) &
    !$acc   present(DENS_dt,MOMX_dt,MOMY_dt,MOMZ_dt,RHOT_dt, DPhydDx,DPhydDy,Coriolis, tend_tmp, drho, EMap3Dto2D) async(1)
    do ke = NeS, NeE
    do pz=1, Nnode_v
      ke2d = EMap3Dto2D(ke)
      !$acc loop vector
      do ph=1, Nnode_h1D**2
        DENS_dt(ph,pz,ke) = - tend_tmp(ph,DENS_VID,pz,ke)
        RHOT_dt(ph,pz,ke) = - tend_tmp(ph,RHOT_VID,pz,ke)


        MOMZ_dt(ph,pz,ke) = - tend_tmp(ph,MOMZ_VID,pz,ke) &
                        - Grav * drho(ph,pz,ke)
        
        cor = CORIOLIS(ph,ke2d)
        MOMX_dt(ph,pz,ke) = - tend_tmp(ph,MOMX_VID,pz,ke) &
                        - DPhydDx(ph,pz,ke) &
                        + cor * MOMY_(ph,pz,ke)
        MOMY_dt(ph,pz,ke) = - tend_tmp(ph,MOMY_VID,pz,ke) &
                        - DPhydDy(ph,pz,ke) &
                        - cor * MOMX_(ph,pz,ke)
      end do
    end do
    end do
    !$acc end parallel
    return
  end subroutine atm_dyn_dgm_nonhydro3d_rhot_heve_cal_tend_interior_gpu

end module scale_atm_dyn_dgm_nonhydro3d_rhot_heve_gpu
