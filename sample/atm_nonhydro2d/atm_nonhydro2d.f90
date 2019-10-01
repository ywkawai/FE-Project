#include "scaleFElib.h"
program atm_nonhydro2d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_prc
  use scale_io
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

  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate

  use scale_time_manager, only: &
    TIME_manager_advance,                              &
    TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP,            &
    TIME_NOWDAYSEC, TIME_DTSEC, TIME_NSTEP
  use scale_time, only: TIME_time2label

  use mod_timeint_rk, only: &
    timeint_rk

  use mod_vars, only: &
    vars_Init, vars_Final, vars_Output,                       &
    DDENS, MOMX, MOMZ, DRHOT,                                 &
    VARS_DDENS_ID, VARS_MOMX_ID, VARS_MOMZ_ID, VARS_DRHOT_ID, &
    PROG_VARS_NUM,                                            &
    PROG_VARS_comm, PROG_VARS_list,                           &
    GxU, GzU, GxW, GzW, GxTHETA, GzTHETA,                     &
    VARS_GxU_ID, VARS_GzU_ID, VARS_GxW_ID, VARS_GzW_ID,       &
    VARS_GxTHETA_ID, VARS_GzTHETA_ID,                         &
    AUX_DIFFVARS_NUM,                                         &
    AUX_DIFFVARS_comm, AUX_DIFFVARS_list,                     &
    DPRES, DENS_hydro, PRES_hydro      

  use mod_exp, only: &
    exp_Init, exp_Final, exp_SetInitCond

  !-----------------------------------------------------------------------------
  implicit none

  integer :: NeGX, NeGZ
  integer, parameter :: NLocalMeshPerPrc = 1

  ! The type of initial q (sin, gaussian-hill, cosine-bell, top-hat)
  character(len=H_SHORT) :: InitShapeName
  real(RP) :: InitShapeParams(4)
  ! The type of specified velocify field (constant, rigid-body-rot)
  character(len=H_SHORT) :: VelTypeName 
  real(RP) :: VelTypeParams(4)

  real(RP) :: dom_xmin, dom_xmax
  real(RP) :: dom_zmin, dom_zmax
  logical :: isPeriodicX, isPeriodicZ

  real(RP), parameter :: diffCoef = 75.0_RP*0.0_RP


  type(QuadrilateralElement) :: refElem
  integer :: PolyOrder
  logical, parameter :: LumpedMassMatFlag = .false.
  logical :: InitCond_GalerkinProjFlag 
  integer, parameter :: PolyOrderErrorCheck = 6
  type(sparsemat) :: Dx, Sx, Dz, Sz, Lift
  
  type(MeshRectDom2D), target :: mesh

  integer :: n, k, p
  type(LocalMesh2D), pointer :: lcmesh
  
  character(len=H_SHORT) :: TINTEG_SCHEME_TYPE
  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: nowstep
  integer :: rkstage
  integer :: tintbuf_ind
  character(len=H_MID) :: timelabel

  real(RP), allocatable :: IntrpMat(:,:)
  real(RP) :: intw_intrp(PolyOrderErrorCheck**2)
  real(RP) :: x_intrp(PolyOrderErrorCheck**2)
  real(RP) :: y_intrp(PolyOrderErrorCheck**2)

  integer :: nstep_eval_error
  !-------------------------------------------------------

  call init()
  call set_initcond()
  
  do nowstep=1, TIME_NSTEP
    do rkstage=1, tinteg_lc(1)%nstage
      
      !* Exchange halo data
      call PROF_rapstart( 'exchange_PROGVARS', 1)
      call PROG_VARS_comm%Put(PROG_VARS_list, 1)
      call PROG_VARS_comm%Exchange()
      call PROG_VARS_comm%Get(PROG_VARS_list, 1)
      call PROF_rapend( 'exchange_PROGVARS', 1)
      
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        call PROF_rapstart( 'cal_grad_diffVars', 1)
        call cal_grad_diffVars( &
          GxU%local(n)%val, GzU%local(n)%val, GxW%local(n)%val, GzW%local(n)%val,        &
          GxTheta%local(n)%val, GzTheta%local(n)%val,                                    &
          DDENS%local(n)%val, MOMX%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,  &
          DENS_hydro%local(n)%val, PRES_hydro%local(n)%val,                              &
          lcmesh, lcmesh%refElem2D ) 
        call PROF_rapend( 'cal_grad_diffVars', 1)
      end do

      !* Exchange halo data
      call PROF_rapstart( 'exchange_DIFFVARS', 1)
      call AUX_DIFFVARS_comm%Put(AUX_DIFFVARS_list, 1)
      call AUX_DIFFVARS_comm%Exchange()
      call AUX_DIFFVARS_comm%Get(AUX_DIFFVARS_list, 1)
      call PROF_rapend( 'exchange_DIFFVARS', 1)

      !* Update prognostic variables
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'cal_dyn_tend', 1)
        call cal_dyn_tend( &
           tinteg_lc(n)%tend_buf2D(:,:,:,tintbuf_ind),                                    &
           DDENS%local(n)%val, MOMX%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,  &
           DENS_hydro%local(n)%val, PRES_hydro%local(n)%val,                              &
           GxU%local(n)%val, GzU%local(n)%val, GxW%local(n)%val, GzW%local(n)%val,        &
           GxTheta%local(n)%val, GzTheta%local(n)%val,                                    & 
           lcmesh, lcmesh%refElem2D ) 
        call PROF_rapend( 'cal_dyn_tend', 1)
        
        call PROF_rapstart( 'update_var', 1)
        call tinteg_lc(n)%Advance( rkstage, DDENS%local(n)%val, VARS_DDENS_ID,     &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        call tinteg_lc(n)%Advance( rkstage, MOMX%local(n)%val, VARS_MOMX_ID,       &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )  

        call tinteg_lc(n)%Advance( rkstage, MOMZ%local(n)%val, VARS_MOMZ_ID,       &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        call tinteg_lc(n)%Advance( rkstage, DRHOT%local(n)%val, VARS_DRHOT_ID,     &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )                                                                     
        call PROF_rapend('update_var', 1)      
      end do
    end do

    !* Advance time
    call TIME_manager_advance()

    if (mod(nowstep,1000) == 0) then 
      call TIME_time2label( TIME_NOWDATE, TIME_NOWMS, timelabel )
      LOG_PROGRESS('(A,A)') "time=", trim(timelabel)
    !  call evaluate_error(tsec_)
    end if
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP )
    
    !* Output
    call vars_output( mesh, TIME_NOWDAYSEC )

    if( IO_L ) call flush(IO_FID_LOG)
  end do

  call final()

contains
  subroutine cal_dyn_tend( dQdt, &
    DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, &
    GxU_, GzU_, GxW_, GzW_, GxTHETA_, GzTHETA_,       &
    lmesh, elem)
    
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem
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
    real(RP), intent(in)  :: GxTHETA_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzTHETA_(elem%Np,lmesh%NeA)  

    real(RP) :: Fx(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PROG_VARS_NUM)
    real(RP) :: dens_(elem%Np), RHOT_(elem%Np), dpres_(elem%Np)
    real(RP) :: pres_(elem%Np), u_(elem%Np), w_(elem%Np)

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
      GxU_, GzU_, GxW_, GzW_, GxTHETA_, GzTHETA_,                 & ! (in)
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

      !-- MOMX
      call sparsemat_matmul(Dx, u_(:)*MOMX_(:,k) + dpres_(:) - diffCoef*dens_(:)*GxU_(:,k), Fx)
      call sparsemat_matmul(Dz, w_(:)*MOMX_(:,k)             - diffCoef*dens_(:)*GzU_(:,k) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_MOMX_ID), LiftDelFlx)

      dQdt(:,k,VARS_MOMX_ID) = &
        - (  lmesh%Escale(:,k,1,1) * Fx(:) &
           + lmesh%Escale(:,k,2,2) * Fz(:) &
           + LiftDelFlx )

      !-- MOMZ
      call sparsemat_matmul(Dx, u_(:)*MOMZ_(:,k)             - diffCoef*dens_(:)*GxW_(:,k), Fx)
      call sparsemat_matmul(Dz, w_(:)*MOMZ_(:,k) + dpres_(:) - diffCoef*dens_(:)*GzW_(:,k), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_MOMZ_ID), LiftDelFlx)
      
      dQdt(:,k,VARS_MOMZ_ID) = &
        - (  lmesh%Escale(:,k,1,1) * Fx(:) &
           + lmesh%Escale(:,k,2,2) * Fz(:) &
           + LiftDelFlx )                  &
        - DDENS_(:,k)*Grav
        !- matmul(IntrpMat_VPOrdM1, DDENS_(:,k)) * Grav

      !-- RHOT
      call sparsemat_matmul(Dx, u_(:)*RHOT_(:) - diffCoef*dens_(:)*GxTHETA_(:,k), Fx)
      call sparsemat_matmul(Dz, w_(:)*RHOT_(:) - diffCoef*dens_(:)*GzTHETA_(:,k), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_DRHOT_ID), LiftDelFlx)
      
      dQdt(:,k,VARS_DRHOT_ID) = &
        - (  lmesh%Escale(:,k,1,1) * Fx(:) &
           + lmesh%Escale(:,k,2,2) * Fz(:) &
           + LiftDelFlx )

    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine cal_dyn_tend

  subroutine cal_del_flux_dyn( del_flux, &
      DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, &
      GxU_, GzU_, GxW_, GzW_, GxTHETA_, GzTHETA_,       &
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
    real(RP), intent(in) ::  GxTHETA_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzTHETA_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    real(RP) :: uM, uP, wM, wP, presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd
    real(RP) :: MOMZ_P, MOMX_P, GxU_P, GzU_P, GxW_P, GzW_P, GxTHETA_P, GzTHETA_P
    real(RP) :: gamm, rgamm
    
    integer :: is_topBC, ie_topBC
    integer :: is_btmBC, ie_btmBC
    integer :: is_leftBC, ie_leftBC
    integer :: is_rightBC, ie_rightBC

    real(RP) :: mu

    !------------------------------------------------------------------------

    gamm = CpDry/CvDry
    rgamm = CvDry/CpDry

    is_topBC = elem%Np*lmesh%Ne + elem%Nfp*(lmesh%NeX + lmesh%NeY) + 1
    ie_topBC = is_topBC + elem%Nfp*lmesh%NeX - 1
    is_btmBC = elem%Np*lmesh%Ne + 1
    ie_btmBC = is_btmBC + elem%Nfp*lmesh%NeX - 1

    is_leftBC = elem%Np*lmesh%Ne + elem%Nfp*(2*lmesh%NeX + lmesh%NeY) + 1
    ie_leftBC = is_leftBC + elem%Nfp*lmesh%NeY - 1
    is_rightBC = elem%Np*lmesh%Ne + elem%Nfp*lmesh%NeX + 1
    ie_rightBC = is_rightBC + elem%Nfp*lmesh%NeY - 1

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


      MOMZ_P = MOMZ_(iP)
      MOMX_P = MOMX_(iP)

      GxU_P = GxU_(iP); GzU_P = GzU_(iP)
      GxW_P = GxW_(iP); GzW_P = GzW_(iP)
      GxTHETA_P = GxTHETA_(iP); GzTHETA_P = GzTHETA_(iP)
      if ( (is_topBC <= iP .and. iP <= ie_topBC) .or. (is_btmBC <= iP .and. iP <= ie_btmBC) ) then
        MOMZ_P = - MOMZ_(iM)
        GzU_P = - GzU_(iM)
        GzW_P = - GzW_(iM)
        GzTHETA_P = - GzTHETA_(iM)
      end if
      if ( (is_leftBC <= iP .and. iP <= ie_leftBC) .or. (is_rightBC <= iP .and. iP <= ie_rightBC) ) then
        !MOMX_P = - MOMX_(iM)
        GxU_P = - GxU_(iM)
        GxW_P = - GxW_(iM)
        GxTHETA_P = - GxTHETA_(iM)
      end if

      VelM = (MOMX_(iM)*nx(i) + MOMZ_(iM)*nz(i))/densM
      VelP = (MOMX_P   *nx(i) + MOMZ_P   *nz(i))/densP

      alpha = max( sqrt(gamm*presM/densM) + abs(VelM), sqrt(gamm*presP/densP) + abs(VelP)  )
      mu = 2.0_RP * (PolyOrder+1)*(PolyOrder+2) / 2.0_RP / 400.0_RP
      
      del_flux(i,VARS_DDENS_ID) = 0.5_RP*(                    &
                    ( densP*VelP - densM*VelM )               &
                    - alpha*(DDENS_(iP) - DDENS_(iM))   )
      
      del_flux(i,VARS_MOMX_ID) = 0.5_RP*(                                &
                    ( MOMX_(iP)*VelP - MOMX_(iM)*VelM )                  &
                    + ( dpresP - dpresM )*nx(i)                          &
                    - diffCoef*( densP*GxU_P - densM*GxU_(iM) )*nx(i)    &
                    - diffCoef*( densP*GzU_P - densM*GzU_(iM) )*nz(i)    &
                    - alpha*(MOMX_P - MOMX_(iM))                         &
                    - mu*diffCoef*(densP + densM)*(MOMX_(iP)/densP - MOMX_(iM)/densM) )
      
      del_flux(i,VARS_MOMZ_ID) = 0.5_RP*(                                      &
                    ( MOMZ_(iP)*VelP - MOMZ_(iM)*VelM)                         &
                    + ( dpresP - dpresM )*nz(i)                                &                    
                    - diffCoef*( densP*GxW_P - densM*GxW_(iM) )*nx(i)          &
                    - diffCoef*( densP*GzW_P - densM*GzW_(iM) )*nz(i)          &
                    - alpha*(MOMZ_P - MOMZ_(iM))                               &
                    - mu*diffCoef*(densP + densM)*(MOMZ_(iP)/densP - MOMZ_(iM)/densM) )
      
      del_flux(i,VARS_DRHOT_ID) = 0.5_RP*(                                     &
                    ( rhotP*VelP - rhotM*VelM )                                &
                  - diffCoef*( densP*GxTHETA_P - densM*GxTHETA_(iM) )*nx(i)    &
                  - diffCoef*( densP*GzTHETA_P - densM*GzTHETA_(iM) )*nz(i)    &
                  - alpha*(DRHOT_(iP) - DRHOT_(iM))                            &
                  - mu*diffCoef*(densP + densM)*(rhotP/densP - rhotM/densM)    )

    end do

    return
  end subroutine cal_del_flux_dyn

  subroutine cal_grad_diffVars( GxU_, GzU_, GxW_, GzW_, GxTHETA_, GzTHETA_, &
    DDENS_, MOMX_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, lmesh, elem )
    
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem
    real(RP), intent(out)  :: GxU_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzU_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GxW_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzW_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GxTheta_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzTheta_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%Ne)
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
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_GxTHETA_ID), LiftDelFlx)
      GxTHETA_(:,k) = lmesh%Escale(:,k,1,1)*Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, DTHETA_, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k,VARS_GzTHETA_ID), LiftDelFlx)
      GzTHETA_(:,k) = lmesh%Escale(:,k,2,2)*Fz(:) + LiftDelFlx(:)

    end do
  end subroutine cal_grad_diffVars

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
    real(RP) :: delU, delW, delTHETA
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
      delTHETA = 0.5_RP*(rhotP/densP - rhotM/densM)

      del_flux(i,VARS_GxU_ID) = delU * nx(i)
      del_flux(i,VARS_GzU_ID) = delU * nz(i)
      del_flux(i,VARS_GxW_ID) = delW * nx(i)
      del_flux(i,VARS_GzW_ID) = delW * nz(i)
      del_flux(i,VARS_GxTHETA_ID) = delTHETA * nx(i)
      del_flux(i,VARS_GzTHETA_ID) = delTHETA * nz(i)
    end do

    return
  end subroutine cal_del_gradDiffVar

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_initcond()
    implicit none
    !------------------------------------------------------------------------

    call exp_SetInitCond( &
      DENS_hydro, PRES_hydro,  DDENS, MOMX, MOMZ, DRHOT, & ! (out)
      mesh )                                               ! (in)
    
    call vars_output( mesh, TIME_NOWDAYSEC )
    LOG_PROGRESS('(A,F13.5,A)') "time=", real(0.0_RP), "[s]"

    return
  end subroutine set_initcond

  subroutine init()
    use scale_const, only: CONST_setup
    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only: TIME_manager_Init 
        
    implicit none

    character(len=H_SHORT) :: exp_name

    namelist /PARAM_ATM_NONHYDRO2D/ &
      dom_xmin, dom_xmax,             &
      dom_zmin, dom_zmax,             &
      isPeriodicX, isPeriodicZ,       &
      NeGX, NeGZ, PolyOrder,          &
      TINTEG_SCHEME_TYPE,             &
      InitCond_GalerkinProjFlag,      &
      exp_name
    
    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr

    character(len=H_MID) :: conf_name
    !------------------------------------------------------------------------

    call PRC_MPIstart( comm )
     
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    call get_command_argument(1, conf_name)
    call IO_setup( "atm_nonhydro2d", trim(conf_name) )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
 

    !--- read namelist

    dom_xmin = 0.0_RP; dom_xmax = 100.0E3_RP
    dom_zmin = 0.0_RP; dom_zmax =  10.0E3_RP
    isPeriodicX = .true.; isPeriodicZ = .false.

    NeGX = 2; NeGZ = 2; PolyOrder = 1 
    TINTEG_SCHEME_TYPE = 'RK_TVD_3'
    InitCond_GalerkinProjFlag = .false.

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATM_NONHYDRO2D,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_ATM_NONHYDRO2D. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATM_NONHYDRO2D)
    
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
    
    call refElem%Init(PolyOrder, LumpedMassMatFlag)
    call Dx%Init(refElem%Dx1)
    call Sx%Init(refElem%Sx1)
    call Dz%Init(refElem%Dx2)
    call Sz%Init(refElem%Sx2)
    call Lift%Init(refElem%Lift)

    call mesh%Init( &
      NeGX, NeGZ,                             &
      dom_xmin, dom_xmax, dom_zmin, dom_zmax, &
      isPeriodicX, isPeriodicZ,               &
      refElem, NLocalMeshPerPrc )

    call mesh%Generate()
    
    ! setup for time integrator
    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call tinteg_lc(n)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, PROG_VARS_NUM,  &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)  )
    end do
    
    !---
    call vars_Init( mesh )
    call exp_Init( exp_name )

    !---
    allocate( IntrpMat(PolyOrderErrorCheck**2,(PolyOrder+1)**2) )
    IntrpMat(:,:) = refElem%GenIntGaussLegendreIntrpMat( PolyOrderErrorCheck,          & ! (in)
                                                         intw_intrp, x_intrp, y_intrp )  ! (out)

    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()

    use scale_time_manager, only: TIME_manager_Final    
    implicit none
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )

    do n=1, mesh%LOCAL_MESH_NUM
      call tinteg_lc(n)%Final()
    end do

    call exp_Final()
    call vars_Final()

    call mesh%Final()
    
    call Dx%Final()
    call Sx%Final()
    call Dz%Final()
    call Sz%Final()
    call Lift%Final()
    call refElem%Final()
    
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()

    return
  end subroutine final
end program atm_nonhydro2d
