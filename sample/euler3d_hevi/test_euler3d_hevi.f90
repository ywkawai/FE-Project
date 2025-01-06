#include "scaleFElib.h"
program test_euler3d_hevi
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof
  use scale_const, only: &
    PI => CONST_PI, &
    GRAV => CONST_GRAV,  &
    Rdry => CONST_Rdry,  &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00
  
  use scale_sparsemat  
  use scale_element_base
  use scale_element_hexahedral
  use scale_localmesh_2d
  use scale_localmesh_3d
  use scale_mesh_cubedom3d

  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_put,   &
    FILE_HISTORY_meshfield_write

  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate

  use scale_time_manager, only: &
    TIME_manager_advance,                              &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,        &
    TIME_DTSEC, TIME_NSTEP 
  use scale_timeint_rk, only: &
    timeint_rk
  
  use mod_fieldutil, only: &
    get_upwind_pos1d => fieldutil_get_upwind_pos1d,         &
    get_profile3d_tracer => fieldutil_get_profile3d_tracer, &
    get_profile3d_flow => fieldutil_get_profile3d_flow

  !-----------------------------------------------------------------------------
  implicit none

  integer :: NprcX, NprcY
  integer :: NeX, NeY, NeGZ

  real(RP), parameter :: dom_xmin =  0.0_RP
  real(RP), parameter :: dom_xmax = +10.0E3_RP
  real(RP), parameter :: dom_ymin =  0.0_RP
  real(RP), parameter :: dom_ymax = +10.0E3_RP
  real(RP), parameter :: dom_zmin =  0.0_RP
  real(RP), parameter :: dom_zmax = +10.0E3_RP
  
  type(HexahedralElement) :: refElem, refElem_l
  integer :: PolyOrder_h, PolyOrder_v
  logical :: LumpedMassMatFlag
  type(sparsemat) :: Dx, Dy, Dz, Lift
  type(sparsemat) :: Dx_l, Lift_l  
  real(RP), allocatable :: IntrpMat_VPOrdM1(:,:)


  integer, parameter :: PROG_VARS_NUM = 5
  integer, parameter :: VARID_DDENS = 1
  integer, parameter :: VARID_MOMX  = 2
  integer, parameter :: VARID_MOMY  = 3
  integer, parameter :: VARID_MOMZ  = 4
  integer, parameter :: VARID_DRHOT = 5

  integer, parameter :: AUX_VARS_NUM = 2
  integer, parameter :: AUXVARID_PRES_hyd  = 1
  integer, parameter :: AUXVARID_DENS_hyd  = 2
  
  type(MeshCubeDom3D), target :: mesh
  type(MeshField3D), target :: DDENS, MOMX, MOMY, MOMZ, DRHOT
  type(MeshField3D), target :: PRES_hyd, DENS_hyd
  type(MeshFieldCommCubeDom3D) :: prgvar_comm
  type(MeshFieldContainer), save :: prgvar_list(PROG_VARS_NUM)  
  type(MeshFieldCommCubeDom3D) :: auxvar_comm
  type(MeshFieldContainer), save :: auxvar_list(AUX_VARS_NUM)  

  integer, save :: PRGVAR_HST_ID(PROG_VARS_NUM)
  integer, save :: AUXVAR_HST_ID(AUX_VARS_NUM)

  integer :: n
  type(LocalMesh3D), pointer :: lcmesh
  
  character(len=H_SHORT) :: TINTEG_SCHEME_TYPE
  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: nowstep
  integer :: rkstage
  integer :: tintbuf_ind

  real(RP) :: tsec_

  real(RP) :: impl_fac
  !-------------------------------------------------------

  call init()
  call set_initcond()

  !---------------------

  do nowstep=1, TIME_NSTEP
    do rkstage=1, tinteg_lc(1)%nstage
      tsec_ =  TIME_NOWDATE(6) + TIME_NOWSUBSEC
      
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)
        impl_fac = tinteg_lc(n)%Get_implicit_diagfac(rkstage)

        call PROF_rapstart( 'cal_dyn_tend_vi', 1)
        call cal_dyn_tend_vi( &
           tinteg_lc(n)%tend_buf2D_im(:,:,VARID_DDENS,tintbuf_ind),         &
           tinteg_lc(n)%tend_buf2D_im(:,:,VARID_MOMX,tintbuf_ind),          &
           tinteg_lc(n)%tend_buf2D_im(:,:,VARID_MOMY,tintbuf_ind),          &
           tinteg_lc(n)%tend_buf2D_im(:,:,VARID_MOMZ,tintbuf_ind),          &
           tinteg_lc(n)%tend_buf2D_im(:,:,VARID_DRHOT,tintbuf_ind),         &
           DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val,        &
           MOMZ%local(n)%val,  DRHOT%local(n)%val,                          &
           DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,                    &
           impl_fac, lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D ) 
        call PROF_rapend( 'cal_dyn_tend_vi', 1)
        
        call PROF_rapstart( 'update_var_vi', 1)
        call tinteg_lc(n)%StoreImplicit2D( rkstage, DDENS%local(n)%val, VARID_DDENS,    &
                                           1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call tinteg_lc(n)%StoreImplicit2D( rkstage, MOMX%local(n)%val, VARID_MOMX,      &
                                           1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call tinteg_lc(n)%StoreImplicit2D( rkstage, MOMY%local(n)%val, VARID_MOMY,      &
                                           1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call tinteg_lc(n)%StoreImplicit2D( rkstage, MOMZ%local(n)%val, VARID_MOMZ,      &
                                           1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call tinteg_lc(n)%StoreImplicit2D( rkstage, DRHOT%local(n)%val, VARID_DRHOT,    &
                                           1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )

        call PROF_rapend('update_var_vi', 1)      
      end do

      !* Exchange halo data
      call PROF_rapstart( 'exchange_halo', 1)
      call prgvar_comm%Put(prgvar_list, 1)
      call prgvar_comm%Exchange()
      call prgvar_comm%Get(prgvar_list, 1)
      call PROF_rapend( 'exchange_halo', 1)

      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'cal_dyn_tend', 1)
        call cal_dyn_tend_he( &
          tinteg_lc(n)%tend_buf2D_ex(:,:,VARID_DDENS,tintbuf_ind),         &
          tinteg_lc(n)%tend_buf2D_ex(:,:,VARID_MOMX,tintbuf_ind),          &
          tinteg_lc(n)%tend_buf2D_ex(:,:,VARID_MOMY,tintbuf_ind),          &
          tinteg_lc(n)%tend_buf2D_ex(:,:,VARID_MOMZ,tintbuf_ind),          &
          tinteg_lc(n)%tend_buf2D_ex(:,:,VARID_DRHOT,tintbuf_ind),         &
          DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val,        &
          MOMZ%local(n)%val,  DRHOT%local(n)%val,                          &
          DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,                    &
          lcmesh, lcmesh%refElem3D ) 
        call PROF_rapend( 'cal_dyn_tend', 1)

        call PROF_rapstart( 'update_var', 1)
        call tinteg_lc(n)%Advance( rkstage, DDENS%local(n)%val, VARID_DDENS,    &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call tinteg_lc(n)%Advance( rkstage, MOMX%local(n)%val, VARID_MOMX,      &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call tinteg_lc(n)%Advance( rkstage, MOMY%local(n)%val, VARID_MOMY,      &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call tinteg_lc(n)%Advance( rkstage, MOMZ%local(n)%val, VARID_MOMZ,      &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call tinteg_lc(n)%Advance( rkstage, DRHOT%local(n)%val, VARID_DRHOT,    &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call PROF_rapend('update_var', 1)      
      end do
    end do
    
    !* Advance time
    call TIME_manager_advance()

    tsec_ = TIME_NOWDATE(6) + TIME_NOWSUBSEC
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )

    !* Output
    call FILE_HISTORY_meshfield_put(PRGVAR_HST_ID(VARID_DDENS), DDENS)
    call FILE_HISTORY_meshfield_put(PRGVAR_HST_ID(VARID_MOMX), MOMX)
    call FILE_HISTORY_meshfield_put(PRGVAR_HST_ID(VARID_MOMY), MOMY)
    call FILE_HISTORY_meshfield_put(PRGVAR_HST_ID(VARID_MOMZ), MOMZ)
    call FILE_HISTORY_meshfield_put(PRGVAR_HST_ID(VARID_DRHOT), DRHOT)
    call FILE_HISTORY_meshfield_write()
  end do

  call final()

contains
  subroutine cal_dyn_tend_he( &
    DDENS_t, MOMX_t, MOMY_t, MOMZ_t, DRHOT_t, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,      &
    DENS_hyd_, PRES_hyd_, lmesh, elem)
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: DDENS_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_t(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: DRHOT_t(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: dens_(elem%Np), RHOT_hyd_(elem%Np), RHOT_(elem%Np), dpres_(elem%Np)
    real(RP) :: pres_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np)

    integer :: ke
    real(RP) :: gamm, rgamm

    !------------------------------------------------------------------------

    gamm = CPdry / CVdry
    rgamm = CVdry / CPdry

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 2)
    call cal_del_flux_dyn_he( del_flux,                                       & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,                                    & ! (in)
      DENS_hyd_, PRES_hyd_,                                                   & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem )                                                             ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 2)

    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 2)
    do ke = lmesh%NeS, lmesh%NeE

      RHOT_hyd_(:) = PRES00/Rdry * (PRES_hyd_(:,ke)/PRES00)**rgamm
      RHOT_(:) = RHOT_hyd_(:) + DRHOT_(:,ke)
      pres_(:) = PRES_hyd_(:,ke) * (1.0_RP + DRHOT_(:,ke)/RHOT_hyd_(:))**gamm
      dpres_(:) = pres_(:) - PRES_hyd_(:,ke)
      dens_(:) = DDENS_(:,ke) + DENS_hyd_(:,ke)

      u_(:) = MOMX_(:,ke)/dens_(:)
      v_(:) = MOMY_(:,ke)/dens_(:)
      w_(:) = MOMZ_(:,ke)/dens_(:)

      !-
      call sparsemat_matmul(Dx, MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, MOMY_(:,ke), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARID_DDENS), LiftDelFlx)

      DDENS_t(:,ke) = - 0d0*( lmesh%Escale(:,ke,1,1) * Fx(:) &
                        + lmesh%Escale(:,ke,2,2) * Fy(:) &
                        + LiftDelFlx )
      if(ke==1) then
        write(*,*) "Ex:", - LiftDelFlx(:)
        write(*,*) "Ex DDENS:", DDENS_(:,ke)
      end if

      !-
      call sparsemat_matmul(Dx, u_(:)*MOMX_(:,ke) + pres_(:), Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMX_(:,ke)           , Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMX_(:,ke)           , Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARID_MOMX), LiftDelFlx)
                  
      MOMX_t(:,ke) = - ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                       + lmesh%Escale(:,ke,2,2) * Fy(:) &
                       + lmesh%Escale(:,ke,3,3) * Fz(:) &
                       + LiftDelFlx )
      !-
      call sparsemat_matmul(Dx, u_(:)*MOMY_(:,ke)           , Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMY_(:,ke) + pres_(:), Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMY_(:,ke)           , Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARID_MOMY), LiftDelFlx)
                  
      MOMY_t(:,ke) = - ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                       + lmesh%Escale(:,ke,2,2) * Fy(:) &
                       + lmesh%Escale(:,ke,3,3) * Fz(:) &
                       + LiftDelFlx )
      
      !-
      call sparsemat_matmul(Dx, u_(:)*MOMZ_(:,ke)           , Fx)
      call sparsemat_matmul(Dy, v_(:)*MOMZ_(:,ke)           , Fy)
      call sparsemat_matmul(Dz, w_(:)*MOMZ_(:,ke)           , Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARID_MOMZ), LiftDelFlx)
                  
      MOMZ_t(:,ke) = - ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                       + lmesh%Escale(:,ke,2,2) * Fy(:) &
                       + lmesh%Escale(:,ke,3,3) * Fz(:) &
                       + LiftDelFlx ) 
                                        
      !-
      call sparsemat_matmul(Dx, u_(:)*RHOT_(:), Fx)
      call sparsemat_matmul(Dy, v_(:)*RHOT_(:), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke,VARID_DRHOT), LiftDelFlx)
                  
      DRHOT_t(:,ke) = - ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                        + lmesh%Escale(:,ke,2,2) * Fy(:) &
                        + LiftDelFlx )   
                          
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine cal_dyn_tend_he

  subroutine cal_del_flux_dyn_he( del_flux, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd_, PRES_hyd_, &
    nx, ny, nz, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne,PROG_VARS_NUM)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, MOMZ_P, alpha, swV
    real(RP) :: presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P
    real(RP) :: DENS_M, DENS_P

    real(RP) :: gamm, rgamm
    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)
      rhot_hyd_M = PRES00/Rdry * (PRES_hyd_(iM)/PRES00)**rgamm
      rhot_hyd_P = PRES00/Rdry * (PRES_hyd_(iP)/PRES00)**rgamm
      
      rhotM = rhot_hyd_M + DRHOT_(iM)
      presM = PRES_hyd_(iM) * (1.0_RP + DRHOT_(iM)/rhot_hyd_M)**gamm
      dpresM = presM - PRES_hyd_(iM)*abs(nz(i))

      rhotP = rhot_hyd_P + DRHOT_(iP) 
      presP = PRES_hyd_(iP) * (1.0_RP + DRHOT_(iP)/rhot_hyd_P)**gamm
      dpresP = presP - PRES_hyd_(iP)*abs(nz(i))

      densM = DDENS_(iM) + DENS_hyd_(iM)
      densP = DDENS_(iP) + DENS_hyd_(iP)

      swV = 1.0_RP - nz(i)**2
      VelM = (MOMX_(iM)*nx(i) + MOMY_(iM)*ny(i) + MOMZ_(iM)*nz(i))/densM
      VelP = (MOMX_(iP)*nx(i) + MOMY_(iP)*ny(i) + MOMZ_(iP)*nz(i))/densP
      MOMZ_P = MOMZ_(iP)

      alpha = swV*max( sqrt(gamm*presM/densM) + abs(VelM), &
                       sqrt(gamm*presP/densP) + abs(VelP)  )
    
      del_flux(i,VARID_DDENS) = 0.5_RP*(                   &
                       (MOMX_(iP) - MOMX_(iM))*nx(i)       &
                     + (MOMY_(iP) - MOMY_(iM))*ny(i)       &
                     - alpha * (DDENS_(iP) - DDENS_(iM))   )
       
      del_flux(i,VARID_MOMX) = 0.5_RP*(                    &
                       ( MOMX_(iP)*VelP - MOMX_(iM)*VelM ) &
                     + (dpresP - dpresM )*nx(i)            &
                     - alpha * (MOMX_(iP) - MOMX_(iM))     )
 
      del_flux(i,VARID_MOMY) = 0.5_RP*(                     &
                       ( MOMY_(iP)*VelP - MOMY_(iM)*VelM )  &
                     + ( dpresP - dpresM )*ny(i)            &
                     - alpha * (MOMY_(iP) - MOMY_(iM))      )               
       
      del_flux(i,VARID_MOMZ) = 0.5_RP*(                    &
                       ( MOMZ_P*VelP - MOMZ_(iM)*VelM )    &
                     - alpha * (MOMZ_(iP) - MOMZ_(iM))     )
       
      del_flux(i,VARID_DRHOT) = 0.5_RP*(                    &
                       swV*( rhotP*VelP - rhotM*VelM )     &
                     - alpha *(DRHOT_(iP) - DRHOT_(iM))    )
    end do

    return
  end subroutine cal_del_flux_dyn_he

  subroutine cal_dyn_tend_vi(  &
    DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,             & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd_, PRES_hyd_, & ! (in)
    impl_fac_, lmesh, elem, lmesh2D, elem2D )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D
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
    real(RP), intent(in)  :: DENS_hyd_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: PRES_hyd_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: impl_fac_

    real(RP) :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: PROG_DEL(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: b(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: tend(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: DENS_hyd_z(elem%Np,lmesh%NeZ)
    real(RP) :: PRES_hyd_z(elem%Np,lmesh%NeZ)
    real(RP) :: nz(elem%NfpTot,lmesh%NeZ)
    integer :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer :: ke_x, ke_y, ke_z, ke, p, v
    integer :: itr, m, N
    integer :: f, vs, ve
    logical :: is_converged
    !------------------------------------------------------------------------
    
    N = elem%Np * PROG_VARS_NUM * lmesh%NeZ
    m = N

    do ke_z=1, lmesh%NeZ
      do f=1, elem%Nfaces_h
        vs = 1 + (f-1)*elem%Nfp_h
        ve = vs + elem%Nfp_h - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_h(:,f) + (ke_z-1)*elem%Np
      end do
      do f=1, elem%Nfaces_v
        vs = elem%Nfaces_h*elem%Nfp_h + (f-1)*elem%Nfp_v + 1
        ve = vs + elem%Nfp_v - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_v(:,f) + (ke_z-1)*elem%Np
      end do
      vmapP(:,ke_z) = vmapM(:,ke_z)
    end do

    do ke_z=1, lmesh%NeZ
      vs = elem%Nfp_h*elem%Nfaces_h + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z > 1) then
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,2) + (ke_z-2)*elem%Np
      else
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,2) + (lmesh%NeZ-1)*elem%Np
      end if

      vs = elem%Nfp_h*elem%Nfaces_h + elem%Nfp_v + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z < lmesh%NeZ) then
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1) + ke_z*elem%Np
      else
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1)
      end if
    end do
    
    do ke_y=1, lmesh%NeY
    do ke_x=1, lmesh%NeX

      do ke_z=1, lmesh%NeZ
        ke = ke_x + (ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
        PROG_VARS(:,VARID_DDENS,ke_z) = DDENS_(:,ke)
        PROG_VARS(:,VARID_MOMX,ke_z) = MOMX_(:,ke)
        PROG_VARS(:,VARID_MOMY,ke_z) = MOMY_(:,ke)
        PROG_VARS(:,VARID_MOMZ,ke_z) = MOMZ_(:,ke)
        PROG_VARS(:,VARID_DRHOT,ke_z) = DRHOT_(:,ke)
        PROG_VARS0(:,:,ke_z) = PROG_VARS(:,:,ke_z)
        PROG_DEL(:,:,ke_z) = 0.0_RP

        DENS_hyd_z(:,ke_z) = DENS_hyd_(:,ke)
        PRES_hyd_z(:,ke_z) = PRES_hyd_(:,ke)
        nz(:,ke_z) = lmesh%normal_fn(:,ke,3)
      end do
      
      if ( impl_fac_ > 0.0_RP ) then
        call vi_eval_Ax( Ax(:,:,:),                       & ! (out)
          PROG_VARS, PROG_VARS0, DENS_hyd_z, PRES_hyd_z,  & ! (in)
          Dz, Lift, impl_fac_, lmesh, elem,               & ! (in)
          nz, vmapM, vmapP, ke_x, ke_y, .false. )

        do ke_z=1, lmesh%NeZ
          b(:,:,ke_z) = - Ax(:,:,ke_z) + PROG_VARS0(:,:,ke_z)
        end do

        is_converged = .false.
        do itr=1, 2*int(N/m)

          call vi_GMRES_core( PROG_DEL(:,:,:),             & ! (inout)
            is_converged,                                  & ! (out)
            PROG_VARS(:,:,:), b(:,:,:), N, m,              & ! (in)
            DENS_hyd_z, PRES_hyd_z,                        & ! (in)
            Dz, Lift, impl_fac_, lmesh, elem,               & ! (in)
            nz, vmapM, vmapP, ke_x, ke_y ) 
          
          if (is_converged) exit
        end do ! itr

        do ke_z=1, lmesh%NeZ
          PROG_VARS(:,:,ke_z) = PROG_VARS(:,:,ke_z) + PROG_DEL(:,:,ke_z)
        end do        
      end if

      call vi_eval_Ax( tend(:,:,:),                     & ! (out)
        PROG_VARS, PROG_VARS0, DENS_hyd_z, PRES_hyd_z,  & ! (in)
        Dz, Lift, impl_fac_, lmesh, elem,                & ! (in)
        nz, vmapM, vmapP, ke_x, ke_y, .true. )
      
      do ke_z=1, lmesh%NeZ
        ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
        DENS_dt(:,ke) = - tend(:,VARID_DDENS,ke_z)
        MOMX_dt(:,ke) = - tend(:,VARID_MOMX,ke_z)
        MOMY_dt(:,ke) = - tend(:,VARID_MOMY,ke_z)
        MOMZ_dt(:,ke) = - tend(:,VARID_MOMZ,ke_z)
        RHOT_dt(:,ke) = - tend(:,VARID_DRHOT,ke_z)
      end do

    end do
    end do
    if (impl_fac == 0.0_RP) write(*,*) "------------------------"
    ! write(*,*) "DDENS1", DDENS_(:,1)
    ! write(*,*) "DENS_dt1", DENS_dt(:,1)
    write(*,*) "DRHOT1", DRHOT_(:,10)
    write(*,*) "DRHOT_dt", RHOT_dt(:,10)
!     write(*,*) "DDENS2", DDENS_(:,2)
!     write(*,*) "DENS_dt2", DENS_dt(:,2)
!     write(*,*) "DDENS3", DDENS_(:,10)
!     write(*,*) "DENS_dt3", DENS_dt(:,10)
!     stop

    return
  end subroutine cal_dyn_tend_vi

  !------------------------------------------------

  subroutine vi_GMRES_core( x, is_converged,  & ! (inout)
    x0, b, N, m,                              & ! (in)
    DENS_hyd, PRES_hyd,                       & ! (in)
    Dz, Lift, impl_fac, lmesh, elem,          & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y  )
    
    implicit none
    integer, intent(in) :: N
    integer, intent(in) :: m    
    real(RP), intent(inout) :: x(N)
    logical, intent(out) :: is_converged
    real(RP), intent(in) :: x0(N)
    real(RP), intent(in) :: b(N)
    !---
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y    

    real(RP) :: r0(N)
    real(RP) :: beta
    real(RP) :: v(N,m+1)
    real(RP) :: hj(m+1)
    real(RP) :: g(m+1)
    real(RP) :: wj(N)
    integer :: i, j
    integer :: m_out
    real(RP) :: r(m+1,m)
    real(RP) :: co(m), si(m)
    real(RP) :: tmp1, tmp2
    real(RP) :: y(m)
    real(RP), parameter :: EPS0 = 1.0E-16_RP
    real(RP), parameter :: EPS = 1.0E-16_RP

    !--------------------------------------

    call vi_eval_Ax_lin( wj(:),             & ! (out)
      x, x0, DENS_hyd, PRES_hyd,            & ! (in)
      Dz, Lift, impl_fac, lmesh, elem,      & ! (in)
      nz, vmapM, vmapP, ke_x, ke_y, .false. ) ! (in)

    r0(:) = b(:) - wj(:)
    beta = sqrt(sum(r0(:)**2))
    if (beta < EPS0*N) then
      is_converged = .true.
      return
    end if

    v(:,1) = r0(:)/beta

    g(1) = beta

    m_out = min(m, N)
    is_converged = .false.
    do j=1, min(m, N)
      call vi_eval_Ax_lin( wj(:),            & ! (out)
       v(:,j), x0, DENS_hyd, PRES_hyd,       & ! (in)
       Dz, Lift, impl_fac, lmesh, elem,      & ! (in)
       nz, vmapM, vmapP, ke_x, ke_y, .false. ) ! (in)
      
      do i=1, j
        hj(i) = sum( wj(:)*v(:,i) )
        wj(:) = wj(:) - hj(i)*v(:,i)
      end do
      hj(j+1) = sqrt(sum(wj(:)**2))

      if (abs(hj(j+1)) < EPS0) then
        m_out = j
        if (ke_x==1 .and. ke_y==3) write(*,*) m_out, "small hj=", abs(hj(j+1))
        is_converged = .true.
        exit
      else
        v(:,j+1) = wj(:)/hj(j+1)
      end if

      r(1,j) = hj(1)
      do i=1, j-1
        tmp1 =  co(i)*r(i,j) + si(i)*hj(i+1)
        tmp2 = -si(i)*r(i,j) + co(i)*hj(i+1)
        r(i,j) = tmp1
        r(i+1,j) = tmp2
      end do

      tmp1 = 1.0_RP / sqrt(r(j,j)**2 + hj(j+1)**2)
      co(j) = tmp1 * r(j,j)
      si(j) = tmp1 * hj(j+1)

      g(j+1) = - si(j)*g(j)
      g(j)   =   co(j)*g(j)

      r(j,j) = co(j)*r(j,j) + si(j)*hj(j+1)
      r(j+1,j) = 0.0_RP
      if ( abs(si(j)*g(j)) < EPS ) then
        m_out = j
        if (ke_x==1 .and. ke_y==3) write(*,*) m_out, "converge: RES=",  abs(si(j)*g(j)), si(j), g(j), EPS
        is_converged = .true.
        exit
      end if
    end do

    do j=m_out,1,-1
      y(j) = g(j)
      do i=j+1,m_out
        y(j) = y(j) - r(j,i)*y(i)
      end do
      y(j) = y(j)/r(j,j)
    end do

    x(:) = x(:) + matmul(v(:,1:m_out),y(1:m_out))
    return
  end subroutine vi_GMRES_core

  !---
  subroutine vi_eval_Ax( Ax, &
    PROG_VARS, PROG_VARS0, DENS_hyd_, PRES_hyd_,      & ! (in)
    Dz, Lift, impl_fac, lmesh, elem,                  & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y, cal_tend_flag )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: DENS_hyd_(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd_(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y
    logical, intent(in) :: cal_tend_flag

    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%NeZ,PROG_VARS_NUM)
    real(RP) :: RHOT_hyd(elem%Np), POT(elem%Np)
    real(RP) :: DPRES(elem%Np)
    integer :: ke_z
    integer :: ke
    real(RP) :: gamm, rgamm

    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    call vi_cal_del_flux_dyn( del_flux,                    & ! (out)
      PROG_VARS(:,VARID_DDENS,:), PROG_VARS(:,VARID_MOMX,:),   & ! (in)
      PROG_VARS(:,VARID_MOMY ,:), PROG_VARS(:,VARID_MOMZ,:),   & ! (in)
      PROG_VARS(:,VARID_DRHOT,:),                            & ! (in)
      PROG_VARS0(:,VARID_DDENS,:), PROG_VARS0(:,VARID_MOMX,:), & ! (in)
      PROG_VARS0(:,VARID_MOMY ,:), PROG_VARS0(:,VARID_MOMZ,:), & ! (in)
      PROG_VARS0(:,VARID_DRHOT,:),                           & ! (in)
      DENS_hyd_, PRES_hyd_, nz, vmapM, vmapP,              & ! (in)
      lmesh, elem )                                          ! (in)

    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd_(:,ke_z)/PRES00)**rgamm

      DPRES(:) = PRES_hyd_(:,ke_z) * ((1.0_RP + PROG_VARS(:,VARID_DRHOT,ke_z)/RHOT_hyd(:))**gamm - 1.0_RP)
      POT(:) = (RHOT_hyd(:) + PROG_VARS(:,VARID_DRHOT,ke_z))/(DENS_hyd_(:,ke_z) + PROG_VARS(:,VARID_DDENS,ke_z))

      !- DENS
      call sparsemat_matmul(Dz, PROG_VARS(:,VARID_MOMZ,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,VARID_DDENS), LiftDelFlx)
      Ax(:,VARID_DDENS,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      !- MOMX
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,VARID_MOMX), LiftDelFlx)
      Ax(:,VARID_MOMX,ke_z) = LiftDelFlx(:)

      !-MOMY
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,VARID_MOMY), LiftDelFlx)
      Ax(:,VARID_MOMY,ke_z) = LiftDelFlx(:)

      !-MOMZ
      call sparsemat_matmul(Dz, DPRES(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,VARID_MOMZ), LiftDelFlx)
      Ax(:,VARID_MOMZ,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)             &
                          + Grav * matmul(IntrpMat_VPOrdM1, PROG_VARS(:,VARID_DDENS,ke_z))

      !-RHOT
      call sparsemat_matmul(Dz, POT(:)*PROG_VARS(:,VARID_MOMZ,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,VARID_DRHOT), LiftDelFlx)
      Ax(:,VARID_DRHOT,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      !--
      if ( .not. cal_tend_flag ) then
        Ax(:,:,ke_z) =  PROG_VARS(:,:,ke_z) + impl_fac * Ax(:,:,ke_z)
      end if 

    end do    

    return
  end subroutine vi_eval_Ax

  subroutine vi_cal_del_flux_dyn( del_flux, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,                &
    DDENS0_, MOMX0_, MOMY0_, MOMZ0_, DRHOT0_,           &
    DENS_hyd_, PRES_hyd_, nz, vmapM, vmapP, lmesh, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%NeZ,PROG_VARS_NUM)
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
    real(RP), intent(in) ::  DENS_hyd_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  PRES_hyd_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%NeZ)
    
    integer :: i, p, ke_z, iP, iM
    real(RP) :: alpha0, swV
    real(RP) ::  MOMZ_P
    real(RP) :: rhot_hyd_M, rhot_hyd_P
    real(RP) :: dpresM, dpresP, densM, densP,  pottM, pottP
    real(RP) :: pres0M, pres0P, dens0M, dens0P
    real(RP) :: gamm, rgamm
    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry
    
    do ke_z=1, lmesh%NeZ
    do p=1, elem%NfpTot
      i = p + (ke_z-1)*elem%NfpTot
      iM = vmapM(i); iP = vmapP(i)

      !-
      densM = DENS_hyd_(iM) + DDENS_(iM)
      densP = DENS_hyd_(iP) + DDENS_(iP)
      
      rhot_hyd_M = PRES00/Rdry * (PRES_hyd_(iM)/PRES00)**rgamm
      rhot_hyd_P = PRES00/Rdry * (PRES_hyd_(iP)/PRES00)**rgamm

      pottM = (rhot_hyd_M + DRHOT_(iM)) / densM
      pottP = (rhot_hyd_P + DRHOT_(iP)) / densP

      dpresM = PRES_hyd_(iM) * ((1.0_RP + DRHOT_(iM)/rhot_hyd_M)**gamm - 1.0_RP) 
      dpresP = PRES_hyd_(iP) * ((1.0_RP + DRHOT_(iP)/rhot_hyd_P)**gamm - 1.0_RP) 

      !-
      dens0M = DENS_hyd_(iM) + DDENS0_(iM)
      dens0P = DENS_hyd_(iP) + DDENS0_(iP)

      pres0M = PRES_hyd_(iM) * (1.0_RP + DRHOT0_(iM)/rhot_hyd_M)**gamm
      pres0P = PRES_hyd_(iP) * (1.0_RP + DRHOT0_(iP)/rhot_hyd_P)**gamm

      swV = nz(i)**2      
      alpha0 = swV * max( abs(MOMZ0_(iM)/dens0M) + sqrt(gamm * pres0M/dens0M), &
                          abs(MOMZ0_(iP)/dens0P) + sqrt(gamm * pres0P/dens0P)  )
      
      ! if (iM==iP .and. (ke_z == 1 .or. ke_z == lmesh%NeZ)) then
      !   MOMZ_P = - MOMZ_(iM)
      ! else
        MOMZ_P = MOMZ_(iP)
      !end if

      del_flux(i,VARID_DDENS) = 0.5_RP * (                   &
                    + ( MOMZ_P - MOMZ_(iM) ) * nz(i)       &
                    - alpha0 * ( DDENS_(iP) - DDENS_(iM) ) )
      !if (swV - 1.0_RP < 1.0E-14_RP) then
      !  write(*,*) p, ke_z, del_flux(i,VARID_DDENS), DDENS_(iP), DDENS_(iM), MOMZ_P, MOMZ_(iM)
      !end if
      del_flux(i,VARID_MOMX) = 0.5_RP * (                   &
                    - alpha0 * (  MOMX_(iP) - MOMX_(iM) ) )
      
      del_flux(i,VARID_MOMY) = 0.5_RP * (                   &  
                    - alpha0 * ( MOMY_(iP) - MOMY_(iM) )  )               
      
      del_flux(i,VARID_MOMZ) = 0.5_RP * (                   &
                    + ( dpresP - dpresM ) * nz(i)         &                    
                    - alpha0 * ( MOMZ_P - MOMZ_(iM) )     )
      
      del_flux(i,VARID_DRHOT) = 0.5_RP * (                               &
                    + ( pottP * MOMZ_P  - pottM * MOMZ_(iM) ) * nz(i)  &
                    - alpha0 * ( DRHOT_(iP) - DRHOT_(iM) )             )
    end do
    end do

    return
  end subroutine vi_cal_del_flux_dyn
!-
  subroutine vi_eval_Ax_lin( Ax, &
    PROG_VARS, PROG_VARS0, DENS_hyd_, PRES_hyd_,        & ! (in)
    Dz, Lift, impl_fac, lmesh, elem,                  & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y, cal_tend_flag )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: DENS_hyd_(elem%Np,lmesh%NeZ)
    real(RP), intent(in)  :: PRES_hyd_(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y
    logical, intent(in) :: cal_tend_flag

    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%NeZ,PROG_VARS_NUM)
    real(RP) :: RHOT_hyd(elem%Np)
    real(RP) :: POT0(elem%Np), DENS0(elem%Np)
    real(RP) :: DPRES(elem%Np)
    integer :: ke_z
    integer :: ke
    real(RP) :: gamm, rgamm

    !--------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    call vi_cal_del_flux_dyn_lin( del_flux,                & ! (out)
      PROG_VARS(:,VARID_DDENS,:), PROG_VARS(:,VARID_MOMX,:),   & ! (in)
      PROG_VARS(:,VARID_MOMY ,:), PROG_VARS(:,VARID_MOMZ,:),   & ! (in)
      PROG_VARS(:,VARID_DRHOT,:),                            & ! (in)
      PROG_VARS0(:,VARID_DDENS,:), PROG_VARS0(:,VARID_MOMX,:), & ! (in)
      PROG_VARS0(:,VARID_MOMY ,:), PROG_VARS0(:,VARID_MOMZ,:), & ! (in)
      PROG_VARS0(:,VARID_DRHOT,:),                           & ! (in)
      DENS_hyd_, PRES_hyd_, nz, vmapM, vmapP,                & ! (in)
      lmesh, elem )                                          ! (in)

    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd_(:,ke_z)/PRES00)**rgamm

      DPRES(:) = gamm * PRES_hyd_(:,ke_z) / RHOT_hyd(:)                               &
                 * ( 1.0_RP + PROG_VARS0(:,VARID_DRHOT,ke_z) / RHOT_hyd(:) )**(gamm-1) &
                 * PROG_VARS(:,VARID_DRHOT,ke_z)

      DENS0(:) = DENS_hyd_(:,ke_z) + PROG_VARS0(:,VARID_DDENS,ke_z)
      POT0(:) = ( RHOT_hyd(:) + PROG_VARS0(:,VARID_DRHOT,ke_z) ) / DENS0(:)

      !- DENS
      call sparsemat_matmul(Dz, PROG_VARS(:,VARID_MOMZ,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,VARID_DDENS), LiftDelFlx)
      Ax(:,VARID_DDENS,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)  

      !- MOMX
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,VARID_MOMX), LiftDelFlx)
      Ax(:,VARID_MOMX,ke_z) = LiftDelFlx(:)

      !-MOMY
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,VARID_MOMY), LiftDelFlx)
      Ax(:,VARID_MOMY,ke_z) = LiftDelFlx(:)

      !-MOMZ
      call sparsemat_matmul(Dz, DPRES(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,VARID_MOMZ), LiftDelFlx)
      Ax(:,VARID_MOMZ,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)             &
                          + Grav * matmul(IntrpMat_VPOrdM1, PROG_VARS(:,VARID_DDENS,ke_z))

      !-RHOT
      call sparsemat_matmul(Dz,   POT0(:) * PROG_VARS(:,VARID_MOMZ,ke_z)                                            &
                                + PROG_VARS0(:,VARID_MOMZ,ke_z) / DENS0(:) * PROG_VARS(:,VARID_DRHOT,ke_z)            &
                                - POT0(:) * PROG_VARS0(:,VARID_MOMZ,ke_z) / DENS0(:) * PROG_VARS(:,VARID_DDENS,ke_z), & 
                           Fz )
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,VARID_DRHOT), LiftDelFlx)
      Ax(:,VARID_DRHOT,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      !--
      if ( .not. cal_tend_flag ) then
        Ax(:,:,ke_z) =  PROG_VARS(:,:,ke_z) + impl_fac * Ax(:,:,ke_z)
      end if 

    end do    

    return
  end subroutine vi_eval_Ax_lin

  subroutine vi_cal_del_flux_dyn_lin( del_flux, &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,         &
    DDENS0_, MOMX0_, MOMY0_, MOMZ0_, DRHOT0_,    &
    DENS_hyd_, PRES_hyd_, nz, vmapM, vmapP, lmesh, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%NeZ,PROG_VARS_NUM)
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
    real(RP), intent(in) ::  DENS_hyd_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  PRES_hyd_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%NeZ)
    
    integer :: i, p, ke_z, iP, iM
    real(RP) :: alpha0, swV
    real(RP) :: rhot_hyd_M, rhot_hyd_P
    real(RP) :: dpresM, dpresP, densM, desnP, MOMZ_P
    real(RP) :: pres0M, pres0P, dens0M, dens0P, pott0M, pott0P, MOMZ0_P

    real(RP) :: gamm, rgamm
    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry
    
    do ke_z=1, lmesh%NeZ
    do p=1, elem%NfpTot
      i = p + (ke_z-1)*elem%NfpTot
      iM = vmapM(i); iP = vmapP(i)
      
      rhot_hyd_M = PRES00/Rdry * (PRES_hyd_(iM)/PRES00)**rgamm
      rhot_hyd_P = PRES00/Rdry * (PRES_hyd_(iP)/PRES00)**rgamm

      dpresM = gamm * PRES_hyd_(iM) / rhot_hyd_M * ( 1.0_RP + DRHOT0_(iM) / rhot_hyd_M )**(gamm-1) * DRHOT_(iM)
      dpresP = gamm * PRES_hyd_(iP) / rhot_hyd_P * ( 1.0_RP + DRHOT0_(iP) / rhot_hyd_P )**(gamm-1) * DRHOT_(iP)

      !-
      dens0M = DENS_hyd_(iM) + DDENS0_(iM)
      dens0P = DENS_hyd_(iP) + DDENS0_(iP)

      pott0M = ( rhot_hyd_M + DRHOT0_(iM) ) / dens0M
      pott0P = ( rhot_hyd_P + DRHOT0_(iP) ) / dens0P

      pres0M = PRES_hyd_(iM) * ( 1.0_RP + DRHOT0_(iM) / rhot_hyd_M )**gamm
      pres0P = PRES_hyd_(iP) * ( 1.0_RP + DRHOT0_(iP) / rhot_hyd_P )**gamm

      swV = nz(i)**2      
      alpha0 = swV * max( abs(MOMZ0_(iM)/dens0M) + sqrt(gamm*pres0M/dens0M), &
                          abs(MOMZ0_(iP)/dens0P) + sqrt(gamm*pres0P/dens0P)  )

      ! if (iM==iP .and. (ke_z == 1 .or. ke_z == lmesh%NeZ)) then
      !   MOMZ_P = - MOMZ_(iM)
      !   MOMZ0_P = - MOMZ0_(iM)
      ! else
        MOMZ_P = MOMZ_(iP)
        MOMZ0_P = MOMZ0_(iP)
      ! end if

      del_flux(i,VARID_DDENS) = 0.5_RP * (                    &
                    + ( MOMZ_P - MOMZ_(iM) ) * nz(i)        &
                    - alpha0 * ( DDENS_(iP) - DDENS_(iM) )  )
      
      del_flux(i,VARID_MOMX) = 0.5_RP * (                     &
                    - alpha0 * ( MOMX_(iP) - MOMX_(iM) )    )
      
      del_flux(i,VARID_MOMY) = 0.5_RP * (                     &  
                    - alpha0 * ( MOMY_(iP) - MOMY_(iM) )    )               
      
      del_flux(i,VARID_MOMZ) = 0.5_RP * (                     &
                    + ( dpresP - dpresM ) * nz(i)           &                    
                    - alpha0 * ( MOMZ_P - MOMZ_(iM) )       )
      
      del_flux(i,VARID_DRHOT) = 0.5_RP * (                             &
                       (  pott0P * MOMZ_P                            &
                        - pott0M * MOMZ_(iM)                         &
                        + MOMZ0_P    / dens0P * DRHOT_(iP)           &
                        - MOMZ0_(iM) / dens0M * DRHOT_(iM)           &
                        - pott0P * MOMZ0_P    / dens0P * DDENS_(iP)  &
                        + pott0M * MOMZ0_(iM) / dens0M * DDENS_(iM)  &
                       ) * nz(i)                                     &
                       - alpha0 * ( DRHOT_(iP) - DRHOT_(iM) )        )
    end do
    end do

    return
  end subroutine vi_cal_del_flux_dyn_lin


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_initcond()
    use scale_linalgebra, only: linalgebra_inv
    use scale_polynominal, only: &
      Polynominal_GenLagrangePoly, Polynominal_GenGaussLobattoPt, Polynominal_GenGaussLegendrePt
    implicit none


    integer :: ke, n

    !------------------------------------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do ke=lcmesh%NeS, lcmesh%NeE
        call set_initcond_lc( DENS_hyd%local(n)%val(:,ke), PRES_hyd%local(n)%val(:,ke),                         &
          DDENS%local(n)%val(:,ke), MOMX%local(n)%val(:,ke), MOMY%local(n)%val(:,ke), MOMZ%local(n)%val(:,ke), DRHOT%local(n)%val(:,ke), &
          lcmesh%pos_en(:,ke,1), lcmesh%pos_en(:,ke,2), lcmesh%pos_en(:,ke,3), lcmesh%refElem3D )
      end do
    end do
    write(*,*) "Allowable Dt=", 10.0E3_RP / sqrt(1.0E5_RP/1.0_RP * (CpDry/CvDry)) / 3.0_RP / sqrt(3.0_RP)

    !--
    call auxvar_comm%Put(auxvar_list, 1)
    call auxvar_comm%Exchange()
    call auxvar_comm%Get(auxvar_list, 1)

    !--
    call FILE_HISTORY_meshfield_put(PRGVAR_HST_ID(VARID_DDENS), DDENS)
    call FILE_HISTORY_meshfield_put(PRGVAR_HST_ID(VARID_MOMX), MOMX)
    call FILE_HISTORY_meshfield_put(PRGVAR_HST_ID(VARID_MOMY), MOMY)
    call FILE_HISTORY_meshfield_put(PRGVAR_HST_ID(VARID_MOMZ), MOMZ)
    call FILE_HISTORY_meshfield_put(PRGVAR_HST_ID(VARID_DRHOT), DRHOT)
    call FILE_HISTORY_meshfield_put(AUXVAR_HST_ID(AUXVARID_DENS_hyd), DENS_hyd)
    call FILE_HISTORY_meshfield_put(AUXVAR_HST_ID(AUXVARID_PRES_hyd), PRES_hyd)

    call FILE_HISTORY_meshfield_write()   
  
    !--
    LOG_PROGRESS('(A,F13.5,A)') "t=", real(0.0_RP), "[s]"

    return
  end subroutine set_initcond

  subroutine set_initcond_lc( &
      DENS_hyd_lc, PRES_hyd_lc, &
      DDENS_lc, MOMX_lc, MOMY_lc, MOMZ_lc, DRHOT_lc, x, y, z, elem  )
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd_lc(elem%Np)
    real(RP), intent(out) :: PRES_hyd_lc(elem%Np)
    real(RP), intent(out) :: DDENS_lc(elem%Np)
    real(RP), intent(out) :: MOMX_lc(elem%Np)
    real(RP), intent(out) :: MOMY_lc(elem%Np)
    real(RP), intent(out) :: MOMZ_lc(elem%Np)
    real(RP), intent(out) :: DRHOT_lc(elem%Np)
    real(RP), intent(in) :: x(elem%Np)
    real(RP), intent(in) :: y(elem%Np)
    real(RP), intent(in) :: z(elem%Np)

    real(RP) :: r(elem%Np)
    !------------------------------------------

    DENS_hyd_lc(:) = 1.0_RP
    PRES_hyd_lc(:) = 1.0E5_RP

    ! Advection test
    ! DDENS_lc(:) = 1.E-15_RP * sin(2.0_RP * PI * lcmesh%pos_en(:,ke,3)/(dom_zmax - dom_zmin) )
    ! DRHOT_lc(:) = 0.0_RP
    ! MOMX_lc(:) = 1.0_RP * (DENS_hyd_lc(:) + DDENS_lc(:))
    ! MOMY_lc(:) = 0.0_RP
    ! MOMZ_lc(:) = 0.0_RP

    ! Sound wave test
    DDENS_lc(:) = 0.0_RP

    r(:) = (z(:) - 0.5_RP*(dom_zmax + dom_zmin)) / ((dom_zmax - dom_zmin)*0.1_RP)
    where( abs(r) <= 1.0_RP) 
      DRHOT_lc(:) = 1.E-12_RP * cos(0.5_RP*PI*r(:))
    elsewhere
      DRHOT_lc(:) = 0.0_RP
    end where

    MOMX_lc(:) = 0.0_RP * (DENS_hyd_lc(:) + DDENS_lc(:))
    MOMY_lc(:) = 0.0_RP
    MOMZ_lc(:) = 0.0_RP

    return
  end subroutine set_initcond_lc

  subroutine init()

    use scale_const, only: CONST_setup
    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only: TIME_manager_Init 
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_setup  
    use scale_file_history, only: FILE_HISTORY_reg 
        
    implicit none

    integer, parameter :: NLocalMeshPerPrc = 1
    namelist /PARAM_TEST/ &
      NprcX, NeX, NprcY, NeY, NeGZ,   & 
      PolyOrder_h, PolyOrder_v,       &
      TINTEG_SCHEME_TYPE, &
      LumpedMassMatFlag
    
    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr
    integer :: p1, p2, p_
    type(ElementBase3D), pointer :: elem
    real(RP), allocatable :: invV_VPOrdM1(:,:)
    real(RP), allocatable :: lift_tmp(:,:)

    !------------------------------------------------------------------------

    call PRC_MPIstart( comm )
     
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    call IO_setup( "test_euler3d", "test.conf" )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
  

    !--- read namelist

    NeX = 2; NeY = 2; NeGZ = 2
    PolyOrder_h = 1; PolyOrder_v = 1
    TINTEG_SCHEME_TYPE = 'IMEX_ARK232'
    LumpedMassMatFlag = .false.

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TEST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_TEST. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_TEST)

    ! setup profiler
    call PROF_setup
    call PROF_rapstart( "total", 0 )
    call PROF_rapstart( "init", 1 )

    ! setup constants
    call CONST_setup

    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init

    !------   
    call refElem%Init(PolyOrder_h, PolyOrder_v, LumpedMassMatFlag)
    call refElem_l%Init(PolyOrder_h, PolyOrder_v, .true.)    
    call Dx%Init(refElem%Dx1)
    call Dx_l%Init(refElem_l%Dx1)    
    call Dy%Init(refElem%Dx2)
    call Dz%Init(refElem%Dx3)
    allocate( lift_tmp(refElem%Np,refElem%NfpTot) )
    call Lift%Init(refElem%Lift)
    ! lift_tmp(:,:) = matmul(refElem%M, refElem%Lift)
    ! call Lift%Init(matmul(refElem_l%invM, lift_tmp))
    call Lift_l%Init(refElem_l%Lift)

    call mesh%Init( &
      NeX*NprcX, NeY*NprcY, NeGZ,                                 &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
      .true., .true., .true.,                                     &
      refElem, NLocalMeshPerPrc, NprcX, NprcY )
    
    call mesh%Generate()
    
    ! setup for time integrator
    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call tinteg_lc(n)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, PROG_VARS_NUM,  &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)          )
    end do
    
    !---
    call DDENS%Init( "DDENS", "kg.m-3", mesh )
    call MOMX%Init( "MOMX", "kg.m-2.s-1", mesh )
    call MOMY%Init( "MOMY", "kg.m-2.s-1", mesh )
    call MOMZ%Init( "MOMZ", "kg.m-2.s-1", mesh )
    call DRHOT%Init( "DRHOT", "kg.m-3.K", mesh )
   
    call prgvar_comm%Init(PROG_VARS_NUM, 0, 0, mesh)
    prgvar_list(VARID_DDENS)%field3d => DDENS
    prgvar_list(VARID_MOMX)%field3d => MOMX
    prgvar_list(VARID_MOMY)%field3d => MOMY
    prgvar_list(VARID_MOMZ)%field3d => MOMZ
    prgvar_list(VARID_DRHOT)%field3d => DRHOT

    !--
    call PRES_hyd%Init( "PRES_hyd", "Pa", mesh )
    call DENS_hyd%Init( "DENS_hyd", "kg.m-3", mesh )
    
    call auxvar_comm%Init(AUX_VARS_NUM, 0, 0, mesh)
    auxvar_list(AUXVARID_DENS_hyd)%field3d => DENS_hyd
    auxvar_list(AUXVARID_PRES_hyd)%field3d => PRES_hyd

    !--
    call FILE_HISTORY_meshfield_setup( mesh3d_=mesh )

    call FILE_HISTORY_reg( DDENS%varname, "DDENS", DDENS%unit, PRGVAR_HST_ID(VARID_DDENS), dim_type='XYZ')
    call FILE_HISTORY_reg( MOMX%varname, "MOMX", MOMX%unit, PRGVAR_HST_ID(VARID_MOMX), dim_type='XYZ')
    call FILE_HISTORY_reg( MOMY%varname, "MOMY", MOMY%unit, PRGVAR_HST_ID(VARID_MOMY), dim_type='XYZ')
    call FILE_HISTORY_reg( MOMZ%varname, "MOMZ", MOMZ%unit, PRGVAR_HST_ID(VARID_MOMZ), dim_type='XYZ')
    call FILE_HISTORY_reg( DRHOT%varname, "DRHOT", DRHOT%unit, PRGVAR_HST_ID(VARID_DRHOT), dim_type='XYZ')

    call FILE_HISTORY_reg( DENS_hyd%varname, "DENS_hyd", DENS_hyd%unit, AUXVAR_HST_ID(AUXVARID_DENS_hyd), dim_type='XYZ')
    call FILE_HISTORY_reg( PRES_hyd%varname, "PRES_hyd", PRES_hyd%unit, AUXVAR_HST_ID(AUXVARID_PRES_hyd), dim_type='XYZ')

    !---
    elem => mesh%refElem3D
    allocate( invV_VPOrdM1(elem%Np,elem%Np) )
    allocate( IntrpMat_VPOrdM1(elem%Np,elem%Np) )
    
    InvV_VPOrdM1(:,:) = elem%invV
    do p2=1, elem%Nnode_h1D
    do p1=1, elem%Nnode_h1D
      p_ = p1 + (p2-1)*elem%Nnode_h1D + (elem%Nnode_v-1)*elem%Nnode_h1D**2
      InvV_VPOrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem%V, invV_VPOrdM1)

    !---

    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()

    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    use scale_time_manager, only: TIME_manager_Final    
    implicit none
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )

    call FILE_HISTORY_meshfield_finalize()

    do n=1, mesh%LOCAL_MESH_NUM
      call tinteg_lc(n)%Final()
    end do

    call DDENS%Final()
    call MOMX%Final()
    call MOMY%Final()
    call MOMZ%Final()
    call DRHOT%Final()

    call DENS_hyd%Final()
    call PRES_hyd%Final()

    call prgvar_comm%Final()
    call auxvar_comm%Final()

    call mesh%Final()
    
    call Dx%Final(); call Dy%Final(); call Dz%Final()
    call Lift%Final()
    call refElem%Final()
    
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()

    return
  end subroutine final
end program test_euler3d_hevi
