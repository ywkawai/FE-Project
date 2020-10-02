#include "scaleFElib.h"
program test_advect3d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof

  use scale_sparsemat  
  use scale_element_base
  use scale_element_hexahedral
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
  integer, parameter :: NLocalMeshPerPrc = 1

  ! The type of initial q (sin, gaussian-hill, cosine-bell, top-hat)
  character(len=H_SHORT) :: InitShapeName
  real(RP) :: InitShapeParams(6)
  ! The type of specified velocify field (constant)
  character(len=H_SHORT) :: VelTypeName 
  real(RP) :: VelTypeParams(6)

  real(RP), parameter :: dom_xmin =  0.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: dom_ymin =  0.0_RP
  real(RP), parameter :: dom_ymax = +1.0_RP
  real(RP), parameter :: dom_zmin =  0.0_RP
  real(RP), parameter :: dom_zmax = +1.0_RP
  
  type(HexahedralElement) :: refElem
  integer :: PolyOrder_h, PolyOrder_v
  logical, parameter :: LumpedMassMatFlag = .false.
  logical :: InitCond_GalerkinProjFlag 
  integer, parameter :: PolyOrderErrorCheck = 6
  type(sparsemat) :: Dx, Sx, Dy, Sy, Dz, Sz, Lift
  
  type(MeshCubeDom3D), target :: mesh
  type(MeshField3D), target :: q, qexact  
  type(MeshField3D), target :: u, v, w
  integer, parameter :: PROG_VARS_NUM = 1
  type(MeshFieldCommCubeDom3D) :: fields_comm
  type(MeshFieldContainer) :: field_list(4)  
  integer :: HST_ID(2)

  integer :: n
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
  real(RP) :: impl_fac
  !-------------------------------------------------------

  call init()
  call set_initcond()

  field_list(1)%field3d => q
  field_list(2)%field3d => u
  field_list(3)%field3d => v
  field_list(4)%field3d => w

  do nowstep=1, TIME_NSTEP
    do rkstage=1, tinteg_lc(1)%nstage
      tsec_ =  TIME_NOWDATE(6) + TIME_NOWSUBSEC
      
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)
        impl_fac = tinteg_lc(n)%Get_implicit_diagfac(rkstage)

        call PROF_rapstart( 'cal_dyn_tend_vi', 1)
        call cal_dyn_tend_vi( &
           tinteg_lc(n)%tend_buf2D_im(:,:,RKVAR_Q,tintbuf_ind),              &
           q%local(n)%val, u%local(n)%val, v%local(n)%val, w%local(n)%val,   &
           impl_fac, lcmesh, lcmesh%refElem3D ) 
        call PROF_rapend( 'cal_dyn_tend_vi', 1)
        call PROF_rapstart( 'update_var_vi', 1)
        call tinteg_lc(n)%StoreImplicit2D( rkstage, q%local(n)%val, RKVAR_Q,    &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call PROF_rapend('update_var_vi', 1)      
      end do

      !* Exchange halo data
      call PROF_rapstart( 'exchange_halo', 1)
      call fields_comm%Put(field_list, 1)
      call fields_comm%Exchange()
      call fields_comm%Get(field_list, 1)
      call PROF_rapend( 'exchange_halo', 1)

      call PROF_rapstart( 'set_velocity', 1)
      call set_velocity( u, v, w, tsec_ )
      call PROF_rapend( 'set_velocity', 1)  

      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'cal_dyn_tend', 1)
        call cal_dyn_tend_he( &
           tinteg_lc(n)%tend_buf2D_ex(:,:,RKVAR_Q,tintbuf_ind),              &
           q%local(n)%val, u%local(n)%val, v%local(n)%val, w%local(n)%val,   &
           lcmesh, lcmesh%refElem3D ) 
        call PROF_rapend( 'cal_dyn_tend', 1)

        call PROF_rapstart( 'update_var', 1)
        call tinteg_lc(n)%Advance( rkstage, q%local(n)%val, RKVAR_Q,              &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call PROF_rapend('update_var', 1)      
      end do
    end do
    
    !* Advance time
    call TIME_manager_advance()

    tsec_ = TIME_NOWDATE(6) + TIME_NOWSUBSEC
    if (mod(nowstep,nstep_eval_error) == 0) then 
      LOG_PROGRESS('(A,F13.5,A)') "t=", real(tsec_), "[s]"
      call evaluate_error(tsec_)
    end if
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )

    !* Output
    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
    call FILE_HISTORY_meshfield_write()
  end do

  call final()

contains
  subroutine cal_dyn_tend_he( dqdt, q_, u_, v_, w_, lmesh, elem)
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: v_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: w_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)

    integer :: ke
    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 2)
    call cal_del_flux_dyn_he( del_flux,                                       & ! (out)
      q_, u_, v_, w_,                                                         & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem )                                                             ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 2)

    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 2)
    do ke = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,ke)*u_(:,ke), Fx)
      call sparsemat_matmul(Dy, q_(:,ke)*v_(:,ke), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke), LiftDelFlx)

      dqdt(:,ke) = - ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                     + lmesh%Escale(:,ke,2,2) * Fy(:) &
                     + LiftDelFlx )
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine cal_dyn_tend_he

  subroutine cal_del_flux_dyn_he( del_flux, q_, u_, v_, w_, nx, ny, nz, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  v_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  w_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    !------------------------------------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      VelM = u_(iM)*nx(i) + v_(iM)*ny(i)
      VelP = u_(iP)*nx(i) + v_(iP)*ny(i)

      alpha = 0.5_RP*abs(VelM + VelP) * (1.0_RP - nz(i)**2)
      del_flux(i) = 0.5_RP*(               &
          ( q_(iP)*VelP - q_(iM)*VelM )    &
        - alpha*(q_(iP) - q_(iM))        )
    end do

    return
  end subroutine cal_del_flux_dyn_he

  subroutine cal_dyn_tend_vi( tend_q, q_, u_, v_, w_, impl_fac, lmesh, elem)
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: tend_q(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: v_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: w_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: impl_fac

    !------------------------------------------------------------------------
    real(RP) :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: PROG_DEL(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: b(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: tend(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP) :: nz(elem%NfpTot,lmesh%NeZ)
    real(RP) :: w_l(elem%Np,lmesh%NeZ)
    integer :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer :: ke_x, ke_y, ke_z, ke, p
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
        vs = elem%Nfp_h*elem%Nfaces_h + 1 + (f-1)*elem%Nfp_v
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

      vs = elem%Nfp_h*elem%Nfaces_h + 1 + elem%Nfp_v
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
        PROG_VARS(:,1,ke_z) = q_(:,ke)
        PROG_VARS0(:,:,ke_z) = PROG_VARS(:,:,ke_z)
        PROG_DEL(:,:,ke_z) = 0.0_RP

        w_l(:,ke_z) = w_(:,ke)
        nz(:,ke_z) = lmesh%normal_fn(:,ke,3)
      end do
      
      if ( impl_fac > 0.0_RP ) then
        call vi_eval_Ax( Ax(:,:,:),                       & ! (out)
          PROG_VARS, PROG_VARS0, w_l,                     & ! (in)
          Dz, Lift, impl_fac, lmesh, elem,                & ! (in)
          nz, vmapM, vmapP, ke_x, ke_y, .false. )

        do ke_z=1, lmesh%NeZ
          b(:,:,ke_z) = - Ax(:,:,ke_z) + PROG_VARS0(:,:,ke_z)
        end do

        is_converged = .false.
        do itr=1, 2*int(N/m)

          call vi_GMRES_core( PROG_DEL(:,:,:),             & ! (inout)
            is_converged,                                  & ! (out)
            PROG_VARS(:,:,:), b(:,:,:), N, m,              & ! (in)
            w_l,                                           & ! (in)
            Dz, Lift, impl_fac, lmesh, elem,               & ! (in)
            nz, vmapM, vmapP, ke_x, ke_y ) 

          if (is_converged) exit
        end do ! itr

        do ke_z=1, lmesh%NeZ
          PROG_VARS(:,:,ke_z) = PROG_VARS(:,:,ke_z) + PROG_DEL(:,:,ke_z)
        end do        
      end if

      call vi_eval_Ax( tend(:,:,:),                     & ! (out)
        PROG_VARS, PROG_VARS0, w_l,                     & ! (in)
        Dz, Lift, impl_fac, lmesh, elem,                & ! (in)
        nz, vmapM, vmapP, ke_x, ke_y, .true. )
      
      !tend(:,:,:) = 0.0_RP
      do ke_z=1, lmesh%NeZ
        ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
        tend_q(:,ke) = - tend(:,1,ke_z)
      end do

    end do
    end do

    return
  end subroutine cal_dyn_tend_vi

  subroutine cal_del_flux_dyn_vi( del_flux, q_, u_, v_, w_, nx, ny, nz, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  v_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  w_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    !------------------------------------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      VelM = u_(iM)*nx(i) + v_(iM)*ny(i)
      VelP = u_(iP)*nx(i) + v_(iP)*ny(i)

      alpha = 0.5_RP*abs(VelM + VelP) * (1.0_RP - nz(i)**2)
      del_flux(i) = 0.5_RP*(               &
          ( q_(iP)*VelP - q_(iM)*VelM )    &
        - alpha*(q_(iP) - q_(iM))        )
    end do

    return
  end subroutine cal_del_flux_dyn_vi

  subroutine vi_GMRES_core( x, is_converged,  & ! (inout)
    x0, b, N, m,                              & ! (in)
    w_,                                       & ! (in)
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
    real(RP), intent(in)  :: w_(elem%Np,lmesh%NeZ)
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
      x, x0, w_,                            & ! (in)
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
       v(:,j), x0, w_,                       & ! (in)
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


  subroutine vi_eval_Ax( Ax, &
    PROG_VARS, PROG_VARS0, w_,                        & ! (in)
    Dz, Lift, impl_fac, lmesh, elem,                  & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y, cal_tend_flag )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: w_(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y
    logical, intent(in) :: cal_tend_flag

    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%NeZ,PROG_VARS_NUM)
    integer :: ke_z
    integer :: ke

    !--------------------------------------------------------

    call vi_cal_del_flux_dyn( del_flux,                    & ! (out)
      PROG_VARS(:,1,:), PROG_VARS0(:,1,:),                 & ! (in)
      w_, nz, vmapM, vmapP,                                & ! (in)
      lmesh, elem )                                          ! (in)

    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

      !- Q
      call sparsemat_matmul(Dz, w_(:,ke_z) * PROG_VARS(:,1,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke_z,1), LiftDelFlx)
      Ax(:,1,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      !--
      if ( .not. cal_tend_flag ) then
        Ax(:,:,ke_z) =  PROG_VARS(:,:,ke_z) + impl_fac * Ax(:,:,ke_z)
      end if 
    end do    

    return
  end subroutine vi_eval_Ax

  subroutine vi_cal_del_flux_dyn( del_flux, &
    q_, q0_,                                &
    w_, nz, vmapM, vmapP, lmesh, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%NeZ,PROG_VARS_NUM)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  q0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  w_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%NeZ)
    
    integer :: i, p, ke_z, iP, iM
    real(RP) :: alpha0, swV
    !------------------------------------------------------------------------

    
    do ke_z=1, lmesh%NeZ
    do p=1, elem%NfpTot
      i = p + (ke_z-1)*elem%NfpTot
      iM = vmapM(i); iP = vmapP(i)

      !-
      swV = nz(i)**2      
      alpha0 = swV * 0.5_RP * ( w_(iP) + w_(iM) )
      
      del_flux(i,1) = 0.5_RP * (                   &
                    + ( w_(iP)*q_(iP) - w_(iM)*q_(iM) ) * nz(i)  &
                    - alpha0 * ( q_(iP) - q_(iM) )               )
    end do
    end do

    return
  end subroutine vi_cal_del_flux_dyn
!-
  subroutine vi_eval_Ax_lin( Ax, &
    PROG_VARS, PROG_VARS0, w_,                        & ! (in)
    Dz, Lift, impl_fac, lmesh, elem,                  & ! (in)
    nz, vmapM, vmapP, ke_x, ke_y, cal_tend_flag )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: Ax(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: PROG_VARS0(elem%Np,PROG_VARS_NUM,lmesh%NeZ)
    real(RP), intent(in)  :: w_(elem%Np,lmesh%NeZ)
    class(SparseMat), intent(in) :: Dz, Lift
    real(RP), intent(in) :: impl_fac
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(in) :: ke_x, ke_y
    logical, intent(in) :: cal_tend_flag

    real(RP) :: Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%NeZ,PROG_VARS_NUM)
    integer :: ke_z
    integer :: ke

    !--------------------------------------------------------


    call vi_cal_del_flux_dyn_lin( del_flux,                & ! (out)
      PROG_VARS(:,1,:), PROG_VARS0(:,1,:),                 & ! (in)
      w_, nz, vmapM, vmapP,                                & ! (in)
      lmesh, elem )                                          ! (in)

    do ke_z=1, lmesh%NeZ
      ke = Ke_x + (Ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY

      !- DENS
      call sparsemat_matmul(Dz, PROG_VARS(:,1,ke_z), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke_z,1), LiftDelFlx)
      Ax(:,1,ke_z) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)

      !--
      if ( .not. cal_tend_flag ) then
        Ax(:,:,ke_z) =  PROG_VARS(:,:,ke_z) + impl_fac * Ax(:,:,ke_z)
      end if 
    end do    

    return
  end subroutine vi_eval_Ax_lin

  subroutine vi_cal_del_flux_dyn_lin( del_flux, &
    q_, q0_, w_, nz, vmapM, vmapP, lmesh, elem )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%NeZ,PROG_VARS_NUM)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  q0_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) ::  w_(elem%Np*lmesh%NeZ)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%NeZ)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%NeZ)
    
    integer :: i, p, ke_z, iP, iM
    real(RP) :: alpha0, swV
    !------------------------------------------------------------------------
    
    do ke_z=1, lmesh%NeZ
    do p=1, elem%NfpTot
      i = p + (ke_z-1)*elem%NfpTot
      iM = vmapM(i); iP = vmapP(i)
      
      !-
      swV = nz(i)**2      
      alpha0 = swV * 0.5_RP * ( w_(iP) + w_(iM) )
      
      del_flux(i,1) = 0.5_RP * (                   &
                    + ( w_(iP)*q_(iP) - w_(iM)*q_(iM) ) * nz(i)  &
                    - alpha0 * ( q_(iP) - q_(iM) )               )
    end do
    end do

    return
  end subroutine vi_cal_del_flux_dyn_lin

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine evaluate_error(tsec)

    implicit none

    real(DP), intent(in) :: tsec

    real(RP) :: q_intrp(PolyOrderErrorCheck**3)
    real(RP) :: qexact_intrp(PolyOrderErrorCheck**3)
    real(RP) :: x_uwind(refElem%Np), y_vwind(refElem%Np), z_wwind(refElem%Np)
    real(RP) :: x_uwind_intrp(PolyOrderErrorCheck**3), y_vwind_intrp(PolyOrderErrorCheck**3), z_wwind_intrp(PolyOrderErrorCheck**3)
    real(RP) :: pos_intrp(PolyOrderErrorCheck**3,3)
    real(RP) vx(8), vy(8), vz(8)

    real(RP) :: l2error
    real(RP) :: linferror
    real(RP) :: ADV_VELX, ADV_VELY, ADV_VELZ
    integer :: ke
    !------------------------------------------------------------------------

    l2error   = 0.0_RP   
    linferror = 0.0_RP
    ADV_VELX = VelTypeParams(1); ADV_VELY = VelTypeParams(2); ADV_VELZ = VelTypeParams(3)

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do ke=lcmesh%NeS, lcmesh%NeE

        x_uwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,ke,1), ADV_VELX, tsec, dom_xmin, dom_xmax)
        y_vwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,ke,2), ADV_VELY, tsec, dom_ymin, dom_ymax)
        z_wwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,ke,3), ADV_VELZ, tsec, dom_zmin, dom_zmax)

        vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
        vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
        vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)

        pos_intrp(:,1) = vx(1) + 0.5_RP*(x_intrp(:) + 1.0_RP)*(vx(2) - vx(1))
        pos_intrp(:,2) = vy(1) + 0.5_RP*(y_intrp(:) + 1.0_RP)*(vy(3) - vy(1))
        pos_intrp(:,3) = vz(1) + 0.5_RP*(z_intrp(:) + 1.0_RP)*(vz(5) - vz(1))

        x_uwind_intrp(:) = get_upwind_pos1d(pos_intrp(:,1), ADV_VELX, tsec, dom_xmin, dom_xmax)
        y_vwind_intrp(:) = get_upwind_pos1d(pos_intrp(:,2), ADV_VELY, tsec, dom_ymin, dom_ymax)
        z_wwind_intrp(:) = get_upwind_pos1d(pos_intrp(:,3), ADV_VELZ, tsec, dom_zmin, dom_zmax)

        call get_profile3d_tracer( qexact%local(n)%val(:,ke),                      & ! (out)
          InitShapeName, x_uwind, y_vwind, z_wwind, InitShapeParams, refElem%Np )    ! (in)

        call get_profile3d_tracer( qexact_intrp(:),                                                             & ! (out) 
          InitShapeName, x_uwind_intrp, y_vwind_intrp, z_wwind_intrp, InitShapeParams, PolyOrderErrorCheck**3 )   ! (in)

        q_intrp(:) = matmul(IntrpMat, q%local(n)%val(:,ke))

        l2error = l2error &
            + sum(   lcmesh%J(1,ke) * intw_intrp(:) * ( q_intrp(:) - qexact_intrp(:) )**2 )
        
        linferror = max(linferror, maxval(abs(q%local(n)%val(:,ke) - qexact%local(n)%val(:,ke))))
      end do
    end do

    LOG_INFO("evaluate_error_l2",*) sqrt(l2error)/( (dom_xmax - dom_xmin) * (dom_ymax - dom_ymin) * (dom_zmax - dom_zmin) )
    LOG_INFO("evaluate_error_linf",*) linferror
    return
  end subroutine evaluate_error

  subroutine set_velocity( u_, v_, w_, tsec )
    type(MeshField3D), intent(inout) :: u_
    type(MeshField3D), intent(inout) :: v_ 
    type(MeshField3D), intent(inout) :: w_ 
    real(RP), intent(in) :: tsec
    
    integer :: n, ke
    !----------------------------------------
    
    VelTypeParams(5) = tsec
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do ke=lcmesh%NeS, lcmesh%NeE
        call get_profile3d_flow( u%local(n)%val(:,ke), v%local(n)%val(:,ke), w%local(n)%val(:,ke),  & ! (out)
          VelTypeName, lcmesh%pos_en(:,ke,1), lcmesh%pos_en(:,ke,2), lcmesh%pos_en(:,ke,3),         & ! (in)
          VelTypeParams, refElem%Np )                                                                 ! (in)
      end do
    end do
    
    return
  end subroutine set_velocity

  subroutine set_initcond()
    use scale_linalgebra, only: linalgebra_inv
    use scale_polynominal, only: &
      Polynominal_GenLagrangePoly, Polynominal_GenGaussLobattoPt, Polynominal_GenGaussLegendrePt
    implicit none

    real(RP) :: q_intrp(PolyOrderErrorCheck**3)
    real(RP) :: lgl1D_h(refElem%PolyOrder_h+1), lgl1D_v(refElem%PolyOrder_v+1)
    real(RP) :: r_int1D_i(PolyOrderErrorCheck)
    real(RP) :: lagrange_intrp1D_h(PolyOrderErrorCheck,refElem%PolyOrder_h+1)    
    real(RP) :: lagrange_intrp1D_v(PolyOrderErrorCheck,refElem%PolyOrder_v+1)    
    real(RP) :: lagrange_intrp(PolyOrderErrorCheck**3,refElem%Np)
    real(RP) :: pos_intrp(PolyOrderErrorCheck**3,3)
    real(RP) vx(8), vy(8), vz(8)
    integer :: p1, p2, p3, p1_, p2_, p3_
    integer :: ke, n_, l_
    real(RP) int_gphi(refElem%Np)
    !------------------------------------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do ke=lcmesh%NeS, lcmesh%NeE
        call get_profile3d_tracer( qexact%local(n)%val(:,ke),                                        & ! (out)
          InitShapeName, lcmesh%pos_en(:,ke,1), lcmesh%pos_en(:,ke,2), lcmesh%pos_en(:,ke,3),        & ! (in)
          InitShapeParams, refElem%Np )                                                                ! (in)
        
        q%local(n)%val(:,ke) = qexact%local(n)%val(:,ke)
      end do
    end do
    call set_velocity( u, v, w, 0.0_RP )

    if (InitCond_GalerkinProjFlag) then

      lgl1D_h(:) = Polynominal_GenGaussLobattoPt(refElem%PolyOrder_h)
      lgl1D_v(:) = Polynominal_GenGaussLobattoPt(refElem%PolyOrder_v)
      r_int1D_i(:) = Polynominal_GenGaussLegendrePt( PolyOrderErrorCheck )
      lagrange_intrp1D_h(:,:) = Polynominal_GenLagrangePoly(refElem%PolyOrder_h, lgl1D_h, r_int1D_i)
      lagrange_intrp1D_v(:,:) = Polynominal_GenLagrangePoly(refElem%PolyOrder_v, lgl1D_v, r_int1D_i)

      do p3_=1, PolyOrderErrorCheck
      do p2_=1, PolyOrderErrorCheck
      do p1_=1, PolyOrderErrorCheck
        n_= p1_ + (p2_-1)*PolyOrderErrorCheck
        do p3=1, refElem%Nnode_v
        do p2=1, refElem%Nnode_h1D
        do p1=1, refElem%Nnode_h1D
          l_ = p1 + (p2-1)*refElem%Nnode_h1D + (p3-1)*refElem%Nnode_h1D**2
          lagrange_intrp(n_,l_) = lagrange_intrp1D_h(p1_,p1)*lagrange_intrp1D_h(p2_,p2)*lagrange_intrp1D_v(p3_,p3)
        end do
        end do
        end do
      end do
      end do
      end do
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        do ke=lcmesh%NeS, lcmesh%NeE      
          vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
          vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)                                       
          vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)
          pos_intrp(:,1) = vx(1) + 0.5_RP*(x_intrp(:) + 1.0_RP)*(vx(2) - vx(1))
          pos_intrp(:,2) = vy(1) + 0.5_RP*(y_intrp(:) + 1.0_RP)*(vy(3) - vy(1))       
          pos_intrp(:,3) = vz(1) + 0.5_RP*(z_intrp(:) + 1.0_RP)*(vz(5) - vz(1))       
  
          call get_profile3d_tracer( q_intrp(:),                            & ! (out)
            InitShapeName, pos_intrp(:,1), pos_intrp(:,2), pos_intrp(:,3), &  ! (in)
            InitShapeParams, PolyOrderErrorCheck**3 )                         ! (in)        
  
          do l_=1, refElem%Np
            int_gphi(l_) = sum(intw_intrp(:)*lagrange_intrp(:,l_)*q_intrp(:)) 
          end do
          q%local(n)%val(:,ke) = matmul(refElem%invM, int_gphi)
        end do
      end do
    end if

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
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
        
    implicit none

    namelist /PARAM_TEST/ &
      NprcX, NeX, NprcY, NeY, NeGZ,   & 
      PolyOrder_h, PolyOrder_v,       &
      TINTEG_SCHEME_TYPE,             &
      InitShapeName, InitShapeParams, &
      InitCond_GalerkinProjFlag,      &      
      VelTypeName, VelTypeParams,     &
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

    NeX = 2; NeY = 2; NeGZ = 2
    PolyOrder_h = 1; PolyOrder_v = 1
    TINTEG_SCHEME_TYPE = 'RK_TVD_3'
    InitShapeName      = 'sin'
    InitShapeParams(:) = (/ 1.0_RP, 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP, 0.0_RP /)
    VelTypeName        = 'const'
    InitCond_GalerkinProjFlag = .false.
    VelTypeParams(:)   = (/ 1.0_RP, 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP, 0.0_RP /)
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
      NeX*NprcX, NeY*NprcY, NeGZ,                                 &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
      .true., .true., .true.,                                     &
      refElem, NLocalMeshPerPrc, NprcX, NprcY )
    
    call mesh%Generate()
    
    ! setup for time integrator
    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call tinteg_lc(n)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,        &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)  )
    end do
    
    !---
    call q%Init( "q", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )
    call u%Init( "u", "m/s", mesh )
    call v%Init( "v", "m/s", mesh )
    call w%Init( "w", "m/s", mesh )
    call fields_comm%Init(4, 0, mesh)

    call FILE_HISTORY_meshfield_setup( mesh3d_=mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='XYZ')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='XYZ')
    
    !---
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
    call u%Final()
    call v%Final()
    call w%Final()

    call fields_comm%Final()
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
end program test_advect3d
