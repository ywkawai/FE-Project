#include "scaleFElib.h"
program test_advect3dGlobal
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof
  use scale_const, only: &
    RPlanet => CONST_radius
    
  use scale_sparsemat  
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d
  use scale_mesh_cubedspheredom3d
  use scale_cubedsphere_coord_cnv, only: &
    CubedSphereCnv_CS2LonLatCoord, &
    CubedSphereCnv_LonLat2CSVec
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_cubedspheredom3d, only: MeshFieldCommCubedSphereDom3D

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
    get_upwind_pos1d => fieldutil_get_upwind_pos1d,               &
    get_profile3d_tracer => fieldutil_get_profile3dGlobal_tracer, &
    get_profile3d_flow => fieldutil_get_profile3dGlobal_flow


  !-----------------------------------------------------------------------------
  implicit none

  integer :: NeGX
  integer :: NeGY
  integer :: NeGZ
  integer :: NLocalMeshPerPrc

  real(RP) :: dom_zmin
  real(RP) :: dom_zmax

  ! The type of initial q (gaussian-hill, cosine-bell, top-hat)
  character(len=H_SHORT) :: InitShapeName
  real(RP) :: InitShapeParams(4)
  ! The type of specified velocify field (rigid-body-rot)
  character(len=H_SHORT) :: VelTypeName 
  real(RP) :: VelTypeParams(4)
  
  type(HexahedralElement) :: refElem
  integer :: PolyOrder_h
  integer :: PolyOrder_v
  logical, parameter :: LumpedMassMatFlag = .false.
  logical :: InitCond_GalerkinProjFlag 
  integer, parameter :: PolyOrderErrorCheck = 6
  type(sparsemat) :: Dx, Sx, Dy, Sy, Dz, Sz, Lift
  
  type(MeshCubedSphereDom3D), target :: mesh
  type(MeshField3D), target :: q, qexact  
  type(MeshField3D), target :: U, V, W
  type(MeshField3D), target :: Vellon, Vellat
  type(MeshFieldCommCubedSphereDom3D) :: prgvars_comm
  type(MeshFieldContainer) :: prgvars_comm_vars(1)  
  type(MeshFieldCommCubedSphereDom3D) :: auxvars_comm
  type(MeshFieldContainer) :: auxvars_comm_vars(3)  

  integer :: HST_ID(7)

  integer :: n, ke, p
  type(LocalMesh3D), pointer :: lcmesh
  
  character(len=H_SHORT) :: TINTEG_SCHEME_TYPE
  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: nowstep
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1
  real(RP) :: tsec_

  real(RP), allocatable :: IntrpMat(:,:)
  real(RP) :: intw_intrp(PolyOrderErrorCheck**3)
  real(RP) :: x_intrp(PolyOrderErrorCheck**3)
  real(RP) :: y_intrp(PolyOrderErrorCheck**3)
  real(RP) :: z_intrp(PolyOrderErrorCheck**3)

  integer :: nstep_eval_error
  !-------------------------------------------------------

  call init()
  call set_initcond()
  
  prgvars_comm_vars(1)%field3d => q
  auxvars_comm_vars(1)%field3d => W
  auxvars_comm_vars(2)%field3d => U
  auxvars_comm_vars(3)%field3d => V
  
  do nowstep=1, TIME_NSTEP
    do rkstage=1, tinteg_lc(1)%nstage
      tsec_ =  ( dble(nowstep-1) + tinteg_lc(1)%coef_c_ex(rkstage) ) * TIME_DTSEC
      
      call PROF_rapstart( 'set_velocity', 1)
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        call set_velocity_lc( U%local(n)%val, V%local(n)%val, W%local(n)%val, &
         Vellon%local(n)%val, Vellat%local(n)%val, tsec_,                     &
         lcmesh, lcmesh%refElem3D )
      end do
      call PROF_rapend( 'set_velocity', 1)  

      !* Exchange halo data
      call PROF_rapstart( 'exchange_halo', 1)
      call prgvars_comm%Put(prgvars_comm_vars, 1)
      call prgvars_comm%Exchange()
      call prgvars_comm%Get(prgvars_comm_vars, 1)

      call auxvars_comm%Put(auxvars_comm_vars, 1)
      call auxvars_comm%Exchange()
      call auxvars_comm%Get(auxvars_comm_vars, 1)    
      call PROF_rapend( 'exchange_halo', 1)


      !* Update prognostic variables
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'cal_dyn_tend', 1)
        call cal_dyn_tend( &
           tinteg_lc(n)%tend_buf2D_ex(:,:,RKVAR_Q,tintbuf_ind),               &
           q%local(n)%val, U%local(n)%val, V%local(n)%val,  W%local(n)%val,   &
           lcmesh, lcmesh%refElem3D ) 
        call PROF_rapend( 'cal_dyn_tend', 1)

        call PROF_rapstart( 'update_var', 1)
        call tinteg_lc(n)%Advance( rkstage, q%local(n)%val, RKVAR_Q,            &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call PROF_rapend('update_var', 1)      
      end do
    end do
    
    !* Advance time
    call TIME_manager_advance()

    tsec_ =  dble(nowstep) * TIME_DTSEC
    if (mod(nowstep,nstep_eval_error) == 0) then 
      LOG_PROGRESS('(A,F13.5,A)') "t=", real(tsec_), "[s]"
      call evaluate_error(tsec_)
    end if
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )

    !* Output
    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
    call FILE_HISTORY_meshfield_put(HST_ID(3), U)
    call FILE_HISTORY_meshfield_put(HST_ID(4), V)
    call FILE_HISTORY_meshfield_put(HST_ID(5), W)
    call FILE_HISTORY_meshfield_put(HST_ID(6), Vellon)
    call FILE_HISTORY_meshfield_put(HST_ID(7), Vellat)
    call FILE_HISTORY_meshfield_write()

    if( IO_L ) call flush(IO_FID_LOG)
  end do

  call final()

contains
  subroutine cal_dyn_tend( dqdt, q_, U_, V_, W_, lmesh, elem)
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: U_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: V_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: W_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np)
    real(RP) :: LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)

    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 2)
    call cal_del_flux_dyn( del_flux,                                                       & ! (out)
      q_, U_, V_, W_,                                                                      & ! (in)
      lmesh%Gsqrt, lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                                            & ! (in)
      lmesh, elem )                                                                          ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 2)

    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 2)
    !$omp parallel do private(Fx, Fy, Fz, LiftDelFlx)
    do ke = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * q_(:,ke) * U_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * q_(:,ke) * V_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * q_(:,ke) * W_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx)

      dqdt(:,ke) = - ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                     + lmesh%Escale(:,ke,2,2) * Fy(:) &
                     + lmesh%Escale(:,ke,3,3) * Fz(:) &
                     + LiftDelFlx ) / lmesh%Gsqrt(:,ke)
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine cal_dyn_tend

  subroutine cal_del_flux_dyn( del_flux, q_, U_, V_, W_, Gsqrt_, nx, ny, nz, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  U_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  V_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  W_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  Gsqrt_(elem%Np*lmesh%Ne)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    !------------------------------------------------------------------------

    !$omp parallel do private(iM, iP, VelM, VelP, alpha)
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      VelM = U_(iM) * nx(i) + V_(iM) * ny(i) + W_(iM) * nz(i)
      VelP = U_(iP) * nx(i) + V_(iP) * ny(i) + W_(iP) * nz(i)
      
      alpha = 0.5_RP * abs( VelM + VelP )
      del_flux(i) = 0.5_RP * (              &
          ( q_(iP) * VelP - q_(iM) * VelM ) &
        - alpha * ( q_(iP) - q_(iM) )       &
        ) * Gsqrt_(iM)
    end do
    
    return
  end subroutine cal_del_flux_dyn

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine evaluate_error(tsec)

    implicit none

    real(DP), intent(in) :: tsec

    real(RP) :: q_intrp(PolyOrderErrorCheck**3)
    real(RP) :: qexact_intrp(PolyOrderErrorCheck**3)
    real(RP) :: x_uwind(refElem%Np), y_vwind(refElem%Np), z_vwind(refElem%Np)
    real(RP) :: x_uwind_intrp(PolyOrderErrorCheck**3)
    real(RP) :: y_vwind_intrp(PolyOrderErrorCheck**3)
    real(RP) :: z_vwind_intrp(PolyOrderErrorCheck**3)
    real(RP) :: pos_intrp(PolyOrderErrorCheck**3,3)
    real(RP) vx(6), vy(6), vz(6)

    real(RP) :: l2error
    real(RP) :: linferror
    real(RP) :: ADV_VELX, ADV_VELY, ADV_VELZ
    !------------------------------------------------------------------------

    l2error   = 0.0_RP   
    linferror = 0.0_RP
    ADV_VELX = VelTypeParams(1); ADV_VELY = VelTypeParams(2)

    ! do n=1, mesh%LOCAL_MESH_NUM
    !   lcmesh => mesh%lcmesh_list(n)
    !   do k=lcmesh%NeS, lcmesh%NeE

    !     x_uwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,k,1), ADV_VELX, tsec, dom_xmin, dom_xmax)
    !     y_vwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,k,2), ADV_VELY, tsec, dom_ymin, dom_ymax)

    !     vx(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),1)
    !     vy(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),2)
    !     pos_intrp(:,1) = vx(1) + 0.5_RP*(x_intrp(:) + 1.0_RP)*(vx(2) - vx(1))
    !     pos_intrp(:,2) = vy(1) + 0.5_RP*(y_intrp(:) + 1.0_RP)*(vy(3) - vy(1))
    !     x_uwind_intrp(:) = get_upwind_pos1d(pos_intrp(:,1), ADV_VELX, tsec, dom_xmin, dom_xmax)
    !     y_vwind_intrp(:) = get_upwind_pos1d(pos_intrp(:,2), ADV_VELY, tsec, dom_ymin, dom_ymax)

    !     call get_profile2d_tracer( qexact%local(n)%val(:,k),             & ! (out)
    !       InitShapeName, x_uwind, y_vwind, InitShapeParams, refElem%Np )   ! (in)

    !     call get_profile2d_tracer( qexact_intrp(:),                                              & ! (out) 
    !       InitShapeName, x_uwind_intrp, y_vwind_intrp, InitShapeParams, PolyOrderErrorCheck**2 )   ! (in)

    !     q_intrp(:) = matmul(IntrpMat, q%local(n)%val(:,k))

    !     l2error = l2error &
    !         + sum(   lcmesh%J(1,k) * intw_intrp(:) * ( q_intrp(:) - qexact_intrp(:) )**2 )
        
    !     linferror = max(linferror, maxval(abs(q%local(n)%val(:,k) - qexact%local(n)%val(:,k))))
    !   end do
    ! end do

    ! LOG_INFO("evaluate_error_l2",*) sqrt(l2error)/( (dom_xmax - dom_xmin) * (dom_ymax - dom_ymin) )
    ! LOG_INFO("evaluate_error_linf",*) linferror

    return
  end subroutine evaluate_error

  subroutine set_velocity_lc( U_, V_, W_, Vellon_, Vellat_, &
      tsec, lmesh, elem )
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: U_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: V_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: W_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: Vellon_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: Vellat_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: tsec

    integer :: ke_
    real(RP) :: svec(elem%Np,lmesh%Ne,2)
    real(RP) :: lon3D(elem%Np), lat3D(elem%Np)

    !----------------------------------------

    VelTypeParams(4) = tsec

    !$omp parallel do private(lon3D, lat3D)
    do ke_=lmesh%NeS, lmesh%NeE
      lon3D(:) = lmesh%lon2D(elem%IndexH2Dto3D(:),lmesh%EMap3Dto2D(ke_))
      lat3D(:) = lmesh%lat2D(elem%IndexH2Dto3D(:),lmesh%EMap3Dto2D(ke_))

      call get_profile3d_flow( svec(:,ke_,1), svec(:,ke_,2), W_(:,ke_), & ! (out)
        VelTypeName, lon3D(:), lat3D(:), lmesh%pos_en(:,ke_,3),         & ! (in)
        VelTypeParams, elem%Np                                          ) ! (in)
      
      Vellon_(:,ke_) = svec(:,ke_,1)
      Vellat_(:,ke_) = svec(:,ke_,2)

      svec(:,ke_,1) = svec(:,ke_,1) / cos(lat3D(:))
    end do

    call CubedSphereCnv_LonLat2CSVec( &
      lmesh%panelID, lmesh%pos_en(:,:,1), lmesh%pos_en(:,:,2),    &
      lmesh%Ne * elem%Np, RPlanet, svec(:,:,1), svec(:,:,2),      &
      U_(:,lmesh%NeS:lmesh%NeE), V_(:,lmesh%NeS:lmesh%NeE)        )

    return
  end subroutine set_velocity_lc

  subroutine set_initcond()
    use scale_linalgebra, only: linalgebra_inv
    use scale_polynominal, only: &
      Polynominal_GenLagrangePoly, Polynominal_GenGaussLobattoPt, Polynominal_GenGaussLegendrePt
    implicit none

    real(RP) :: q_intrp(PolyOrderErrorCheck**3)
    real(RP) :: lgl1D_h(refElem%PolyOrder_h+1), lgl1D_v(refElem%PolyOrder_v+1)
    real(RP) :: r_int1D_i(PolyOrderErrorCheck)
    real(RP) :: lagrange_intrp1D(PolyOrderErrorCheck,refElem%PolyOrder_h+1)    
    real(RP) :: lagrange_intrp(PolyOrderErrorCheck**2,refElem%Np)
    real(RP) :: pos_intrp(PolyOrderErrorCheck**3,3)
    real(RP) vx(6), vy(6), vz(6)
    integer :: p1, p2, p3, p1_, p2_, p3_
    integer :: n_, l_
    real(RP) int_gphi(refElem%Np)

    real(RP) :: lon3D(refElem%Np), lat3D(refElem%Np)
    integer :: ke_
    type(ElementBase3D), pointer :: elem
    !------------------------------------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      elem => lcmesh%refElem3D

      !$omp parallel do private(ke_, lon3D, lat3D)
      do ke_=lcmesh%NeS, lcmesh%NeE
        lon3D(:) = lcmesh%lon2D(elem%IndexH2Dto3D(:),lcmesh%EMap3Dto2D(ke_))
        lat3D(:) = lcmesh%lat2D(elem%IndexH2Dto3D(:),lcmesh%EMap3Dto2D(ke_))
  
        call get_profile3d_tracer( qexact%local(n)%val(:,ke_),        & ! (out)
          InitShapeName, lon3D(:), lat3D(:), lcmesh%pos_en(:,ke_,3),  & ! (in)
          RPlanet, InitShapeParams, elem%Np                           ) ! (in)
        
        q%local(n)%val(:,ke_) = qexact%local(n)%val(:,ke_)
      end do

      call set_velocity_lc( U%local(n)%val, V%local(n)%val, W%local(n)%val, &
        Vellon%local(n)%val, Vellat%local(n)%val,                           &
        0.0_RP, lcmesh, elem )
    end do

    ! if (InitCond_GalerkinProjFlag) then

    !   lgl1D(:) = Polynominal_GenGaussLobattoPt(refElem%PolyOrder)
    !   r_int1D_i(:) = Polynominal_GenGaussLegendrePt( PolyOrderErrorCheck )
    !   lagrange_intrp1D(:,:) = Polynominal_GenLagrangePoly(refElem%PolyOrder, lgl1D, r_int1D_i)
    !   do p2_=1, PolyOrderErrorCheck
    !   do p1_=1, PolyOrderErrorCheck
    !     n_= p1_ + (p2_-1)*PolyOrderErrorCheck
    !     do p2=1, refElem%Nfp
    !     do p1=1, refElem%Nfp
    !       l_ = p1 + (p2-1)*refElem%Nfp
    !       lagrange_intrp(n_,l_) = lagrange_intrp1D(p1_,p1)*lagrange_intrp1D(p2_,p2)
    !     end do
    !     end do
    !   end do
    !   end do

    !   do n=1, mesh%LOCAL_MESH_NUM
    !     lcmesh => mesh%lcmesh_list(n)
    !     do k=lcmesh%NeS, lcmesh%NeE      
    !       vx(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),1)
    !       vy(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),2)                                                  
    !       pos_intrp(:,1) = vx(1) + 0.5_RP*(x_intrp(:) + 1.0_RP)*(vx(2) - vx(1))
    !       pos_intrp(:,2) = vy(1) + 0.5_RP*(y_intrp(:) + 1.0_RP)*(vy(3) - vy(1))       
  
    !       call get_profile2d_tracer( q_intrp(:),                                                     & ! (out)
    !         InitShapeName, pos_intrp(:,1), pos_intrp(:,2), InitShapeParams, PolyOrderErrorCheck**2 )   ! (in)        
  
    !       do l_=1, refElem%Np
    !         int_gphi(l_) = sum(intw_intrp(:)*lagrange_intrp(:,l_)*q_intrp(:)) 
    !       end do
    !       q%local(n)%val(:,k) = matmul(refElem%invM, int_gphi)
    !     end do
    !   end do
    ! end if

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
    call FILE_HISTORY_meshfield_put(HST_ID(3), U)
    call FILE_HISTORY_meshfield_put(HST_ID(4), V)
    call FILE_HISTORY_meshfield_put(HST_ID(5), W)
    call FILE_HISTORY_meshfield_put(HST_ID(6), Vellon)
    call FILE_HISTORY_meshfield_put(HST_ID(7), Vellat)
    call FILE_HISTORY_meshfield_write()   
  
    LOG_PROGRESS('(A,F13.5,A)') "t=", real(0.0_RP), "[s]"
    call evaluate_error(0.0_RP)

    return
  end subroutine set_initcond

  subroutine init()

    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only: TIME_manager_Init 
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_setup  
    use scale_file_history, only: FILE_HISTORY_reg 
    use scale_const, only: &
      CONST_setup
        
    implicit none

    namelist /PARAM_TEST/ &
      NeGX, NeGY, NeGZ, NLocalMeshPerPrc, &
      PolyOrder_h, PolyOrder_v,           &
      dom_zmin, dom_zmax,                 &
      TINTEG_SCHEME_TYPE,                 &
      InitShapeName, InitShapeParams,     &
      InitCond_GalerkinProjFlag,          &      
      VelTypeName, VelTypeParams,         &
      nstep_eval_error
    
    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr
    !------------------------------------------------------------------------

    call PRC_MPIstart( comm )
     
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    call IO_setup( "test_advect3d", "test.conf" )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
      
    !--- read namelist

    NeGX = 2; NeGY = 2; NeGZ = 2
    PolyOrder_h = 1
    PolyOrder_v = 1 
    TINTEG_SCHEME_TYPE = 'ERK_SSP_3s3o'
    InitShapeName      = 'sin'
    InitShapeParams(:) = (/ 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP /)
    VelTypeName        = 'const'
    InitCond_GalerkinProjFlag = .false.
    VelTypeParams(:)   = (/ 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP /)
    nstep_eval_error = 5

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

    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init

    ! setup constat
    call CONST_setup

    !------   
    
    call refElem%Init(PolyOrder_h, PolyOrder_v, LumpedMassMatFlag)
    call Dx%Init(refElem%Dx1)
    call Sx%Init(refElem%Sx1)
    call Dy%Init(refElem%Dx2)
    call Sy%Init(refElem%Sx2)
    call Dz%Init(refElem%Dx3)
    call Sz%Init(refElem%Sx3)
    call Lift%Init(refElem%Lift)

    call mesh%Init( &
      NeGX, NeGY, NeGz, RPlanet, dom_zmin, dom_zmax, &
      refElem, NLocalMeshPerPrc                      )
    
    call mesh%Generate()
    
    ! setup for time integrator
    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call tinteg_lc(n)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,      &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)  )
    end do

    !---
    call q%Init( "q", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )
    call U%Init( "U", "s-1", mesh )
    call V%Init( "V", "s-1", mesh )
    call W%Init( "W", "s-1", mesh )

    call prgvars_comm%Init(1, 0, mesh)
    call auxvars_comm%Init(1, 1, mesh)
    
    call Vellon%Init( "Vellon", "m/s", mesh )
    call Vellat%Init( "Vellat", "m/s", mesh )

    call FILE_HISTORY_meshfield_setup( meshCubedSphere3d_=mesh )
    
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='XYZ')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='XYZ')
    call FILE_HISTORY_reg( U%varname, "U", U%unit, HST_ID(3), dim_type='XYZ')
    call FILE_HISTORY_reg( V%varname, "V", V%unit, HST_ID(4), dim_type='XYZ')
    call FILE_HISTORY_reg( W%varname, "W", W%unit, HST_ID(5), dim_type='XYZ')
    call FILE_HISTORY_reg( Vellon%varname, "Vellon", Vellon%unit, HST_ID(6), dim_type='XYZ')
    call FILE_HISTORY_reg( Vellat%varname, "Vellat", Vellat%unit, HST_ID(7), dim_type='XYZ')

    allocate( IntrpMat(PolyOrderErrorCheck**3,(PolyOrder_h+1)**2*(PolyOrder_v+1)) )
    IntrpMat(:,:) = refElem%GenIntGaussLegendreIntrpMat( PolyOrderErrorCheck,                   & ! (in)
                                                         intw_intrp, x_intrp, y_intrp, z_intrp )  ! (out)

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

    call q%Final()
    call qexact%Final()
    call U%Final()
    call V%Final()
    call W%Final()

    call prgvars_comm%Final()
    call auxvars_comm%Final()
    call mesh%Final()
    
    call Dx%Final()
    call Sx%Final()
    call Dy%Final()
    call Sy%Final()
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
end program test_advect3dGlobal
