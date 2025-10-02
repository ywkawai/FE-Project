!-------------------------------------------------------------------------------
!> Program A sample program: 3-dimensional linear advection test
!! 
!! 
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
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

  ! The type of initial q (sin, gaussian-hill, cosine-bell, top-hat)
  character(len=H_SHORT) :: InitShapeName
  real(RP), save :: InitShapeParams(6)
  ! The type of specified velocify field (constant)
  character(len=H_SHORT) :: VelTypeName 
  real(RP), save :: VelTypeParams(6)

  real(RP), parameter :: dom_xmin =  0.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: dom_ymin =  0.0_RP
  real(RP), parameter :: dom_ymax = +1.0_RP
  real(RP), parameter :: dom_zmin =  0.0_RP
  real(RP), parameter :: dom_zmax = +1.0_RP
  
  type(HexahedralElement) :: refElem
  integer :: PolyOrder_h, PolyOrder_v
  logical :: InitCond_GalerkinProjFlag 
  type(sparsemat) :: Dx, Dy, Dz, Lift
  
  type(MeshCubeDom3D), target :: mesh
  type(MeshField3D), target :: q, qexact  
  type(MeshField3D), target :: u, v, w
  type(MeshFieldCommCubeDom3D) :: fields_comm
  type(MeshFieldContainer), save :: field_list(4)  
  integer, save :: HST_ID(2)

  integer :: n
  type(LocalMesh3D), pointer :: lcmesh
  
  character(len=H_SHORT) :: TINTEG_SCHEME_TYPE
  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: nowstep
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1
  real(RP) :: tsec_

  integer :: PolyOrderErrorCheck
  real(RP), allocatable :: IntrpMat(:,:)
  real(RP), allocatable :: intw_intrp(:)
  real(RP), allocatable :: x_intrp(:)
  real(RP), allocatable :: y_intrp(:)
  real(RP), allocatable :: z_intrp(:)
  integer :: nstep_eval_error
  !-------------------------------------------------------

  call init()
  call set_initcond()

  do nowstep=1, TIME_NSTEP
    do rkstage=1, tinteg_lc(1)%nstage
      tsec_ =  TIME_NOWDATE(6) + TIME_NOWSUBSEC
      
      !* Exchange halo data
      call PROF_rapstart( 'exchange_halo', 1)
      call fields_comm%Put(field_list, 1)
      call fields_comm%Exchange()
      call fields_comm%Get(field_list, 1)
      call PROF_rapend( 'exchange_halo', 1)

      call PROF_rapstart( 'set_velocity', 1)
      call set_velocity( u, v, w, tsec_ )
      call PROF_rapend( 'set_velocity', 1)  

      !* Update prognostic variables
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'cal_dyn_tend', 1)
        call cal_dyn_tend( &
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
  subroutine cal_dyn_tend( dqdt, q_, u_, v_, w_, lmesh, elem)
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
    call cal_bnd_flux_dyn( del_flux,                                          & ! (out)
      q_, u_, v_, w_,                                                         & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem )                                                             ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 2)

    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 2)
    !$omp parallel do private(Fx, Fy, Fz, LiftDelFlx)
    do ke = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,ke) * u_(:,ke), Fx)
      call sparsemat_matmul(Dy, q_(:,ke) * v_(:,ke), Fy)
      call sparsemat_matmul(Dz, q_(:,ke) * w_(:,ke), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke), LiftDelFlx)

      dqdt(:,ke) = - ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                     + lmesh%Escale(:,ke,2,2) * Fy(:) &
                     + lmesh%Escale(:,ke,3,3) * Fz(:) &
                     + LiftDelFlx )
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine cal_dyn_tend

  subroutine cal_bnd_flux_dyn( del_flux, q_, u_, v_, w_, nx, ny, nz, vmapM, vmapP, lmesh, elem )
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

    !$omp parallel do private(i, iM, iP, VelM, VelP, alpha)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      VelM = u_(iM) * nx(i) + v_(iM) * ny(i) + w_(iM) * nz(i)
      VelP = u_(iP) * nx(i) + v_(iP) * ny(i) + w_(iP) * nz(i)

      alpha = 0.5_RP * abs( VelM + VelP )
      del_flux(i) = 0.5_RP * (               &
          ( q_(iP) * VelP - q_(iM) * VelM )  &
        - alpha * ( q_(iP) - q_(iM) )        )
    end do

    return
  end subroutine cal_bnd_flux_dyn

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
      
      !$omp parallel do private( &
      !$omp ke, x_uwind, y_vwind, z_wwind, vx, vy, vz, pos_intrp, x_uwind_intrp, y_vwind_intrp, z_wwind_intrp, &
      !$omp qexact_intrp, q_intrp ) reduction(+: l2error) reduction(max: linferror)
      do ke=lcmesh%NeS, lcmesh%NeE

        x_uwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,ke,1), ADV_VELX, tsec, dom_xmin, dom_xmax)
        y_vwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,ke,2), ADV_VELY, tsec, dom_ymin, dom_ymax)
        z_wwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,ke,3), ADV_VELZ, tsec, dom_zmin, dom_zmax)

        vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
        vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
        vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)

        pos_intrp(:,1) = vx(1) + 0.5_RP * ( x_intrp(:) + 1.0_RP ) * ( vx(2) - vx(1) )
        pos_intrp(:,2) = vy(1) + 0.5_RP * ( y_intrp(:) + 1.0_RP ) * ( vy(3) - vy(1) )
        pos_intrp(:,3) = vz(1) + 0.5_RP * ( z_intrp(:) + 1.0_RP ) * ( vz(5) - vz(1) )

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

    logical, parameter :: LumpedMassMatFlag = .false.
    integer, parameter :: NLocalMeshPerPrc = 1
    namelist /PARAM_TEST/ &
      NprcX, NeX, NprcY, NeY, NeGZ,   & 
      PolyOrder_h, PolyOrder_v,       &
      TINTEG_SCHEME_TYPE,             &
      InitShapeName, InitShapeParams, &
      InitCond_GalerkinProjFlag,      &      
      VelTypeName, VelTypeParams,     &
      PolyOrderErrorCheck,            &
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
    TINTEG_SCHEME_TYPE = 'ERK_SSP_3s3o'
    InitShapeName      = 'sin'
    InitShapeParams(:) = (/ 1.0_RP, 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP, 0.0_RP /)
    VelTypeName        = 'const'
    InitCond_GalerkinProjFlag = .false.
    VelTypeParams(:)   = (/ 1.0_RP, 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP, 0.0_RP /)
    PolyOrderErrorCheck = 6
    nstep_eval_error    = 5

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
    call Dx%Init(refElem%Dx1, storage_format='ELL')
    call Dy%Init(refElem%Dx2, storage_format='ELL')
    call Dz%Init(refElem%Dx3, storage_format='ELL')
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
    
    call fields_comm%Init(4, 0, 0, mesh)
    field_list(1)%field3d => q
    field_list(2)%field3d => u
    field_list(3)%field3d => v
    field_list(4)%field3d => w
  
    call FILE_HISTORY_meshfield_setup( mesh3d_=mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='XYZ')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='XYZ')
    
    !---
    allocate( IntrpMat(PolyOrderErrorCheck**3,(PolyOrder_h+1)**2*(PolyOrder_v+1)) )
    allocate( intw_intrp(PolyOrderErrorCheck**3), x_intrp(PolyOrderErrorCheck**3), y_intrp(PolyOrderErrorCheck**3), z_intrp(PolyOrderErrorCheck**3) )
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
end program test_advect3d
