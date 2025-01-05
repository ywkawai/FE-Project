#include "scalelib.h"
program test_advdiff1d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_const, only: &
    PI  => CONST_PI
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof

  use scale_sparsemat  
  use scale_element_base
  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d

  use scale_localmeshfield_base, only: LocalMeshField1D
  use scale_meshfield_base, only: MeshField1D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_1d, only: MeshFieldComm1D

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_put,           &
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
    get_upwind_pos1d => fieldutil_get_upwind_pos1d,        &
    get_profile1d_tracer => fieldutil_get_profile1d_tracer 
  
  !-----------------------------------------------------------------------------
  implicit none

  ! The type of initial q (sin, gaussian-hill, cosine-bell, top-hat)
  character(len=H_SHORT) :: InitShapeName
  real(RP), save :: InitShapeParams(2)
  ! The type of specified velocify field (constant)
  real(RP) :: ADV_VEL
  ! The coefficient of diffusion
  real(RP) :: DIFF_COEF

  real(RP), parameter :: dom_xmin =  0.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP

  type(LineElement)  :: refElem
  type(sparsemat) :: Dx, Sx, Lift

  type(MeshLineDom1D), target :: mesh
  type(MeshField1D), target :: q, qexact 
  type(MeshField1D), target :: dqdx
  type(MeshField1D), target :: u
  type(MeshFieldComm1D) :: fields_comm
  type(MeshFieldComm1D) :: auxvars_comm
  type(MeshFieldContainer), save :: field_list(2)
  type(MeshFieldContainer), save :: auxvars_list(1)

  integer, save :: HST_ID(2)

  integer :: n, k, p
  type(LocalMesh1D), pointer :: lcmesh
  
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
  integer :: nstep_eval_error
  !-------------------------------------------------------

  call init()
  call set_initcond()

  do nowstep=1, TIME_NSTEP
    do rkstage=1, tinteg_lc(1)%nstage
      !* Exchange halo data
      call PROF_rapstart( 'exchange_halo', 1)
      
      call fields_comm%Put(field_list, 1)
      call fields_comm%Exchange()
      call fields_comm%Get(field_list, 1)

      call PROF_rapend( 'exchange_halo', 1)

      !* Update auxiliary variables
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)

        call PROF_rapstart( 'cal_prg_tend', 1)
        call cal_aux( dqdx%local(n)%val, &
          q%local(n)%val,                &
          lcmesh, lcmesh%refElem1D       )
        call PROF_rapend( 'cal_prg_tend', 1) 
      end do

      call PROF_rapstart( 'exchange_halo_aux', 1)

      call auxvars_comm%Put(auxvars_list, 1)
      call auxvars_comm%Exchange()
      call auxvars_comm%Get(auxvars_list, 1)

      call PROF_rapend( 'exchange_halo_aux', 1)

      !* Update prognostic variables
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)      
        
        call PROF_rapstart( 'cal_prg_tend', 1)
        call cal_prg_tend( &
          tinteg_lc(n)%tend_buf2D_ex(:,:,RKVAR_Q,tintbuf_ind), &
          q%local(n)%val, dqdx%local(n)%val, u%local(n)%val,   &
          lcmesh, lcmesh%refElem1D )
        call PROF_rapend( 'cal_prg_tend', 1) 

        call PROF_rapstart( 'update_var', 1)
        call tinteg_lc(n)%Advance( rkstage, q%local(n)%val, RKVAR_Q, &
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
  subroutine cal_prg_tend( dqdt, q_, dqdx_, u_, lmesh, elem)
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(elementbase1D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dqdx_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)

    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 2)
    call cal_del_flux_prg( del_flux,           & ! (out)
      q_, dqdx_, u_, lmesh%normal_fn(:,:,1),   & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem )    ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 2)

    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 2)
    do k = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,k)*u_(:,k) - DIFF_COEF * dqdx_(:,k), Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k), LiftDelFlx)

      dqdt(:,k) = - (  lmesh%Escale(:,k,1,1) * Fx(:) &
                    + LiftDelFlx )
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine cal_prg_tend

  subroutine cal_del_flux_prg( del_flux, q_, dqdx_, u_, nx, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(elementbase1D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  dqdx_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: alpha
    !------------------------------------------------------------------------

    !$omp parallel do private(iM, iP, alpha)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      alpha = 0.5_RP*abs(u_(iP) + u_(iM))
      del_flux(i) = 0.5_RP * (                      &
          (q_(iP)*u_(iP) - q_(iM)*u_(iM)) * nx(i)   &
        - alpha * (q_(iP) - q_(iM))                 & 
        - DIFF_COEF * (1.0_RP - nx(i))              &
           * (dqdx_(iP) - dqdx_(iM)) * nx(i)        )
    end do

    return
  end subroutine cal_del_flux_prg

  subroutine cal_aux( dqdx_, q_, lmesh, elem)
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(elementbase1D), intent(in) :: elem
    real(RP), intent(out)  :: dqdx_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)

    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 2)
    call cal_del_flux_aux( del_flux,           & ! (out)
      q_, lmesh%normal_fn(:,:,1),              & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem )    ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 2)

    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 2)
    do k = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,k), Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k), LiftDelFlx)

      dqdx_(:,k) = lmesh%Escale(:,k,1,1) * Fx(:) &
                 + LiftDelFlx(:)
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine cal_aux

  subroutine cal_del_flux_aux( del_flux, q_, nx, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(elementbase1D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: delVar
    !------------------------------------------------------------------------

    !$omp parallel do private(iM, iP)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)
      delVar= 0.5_RP * (q_(iP) - q_(iM))
      del_flux(i) = (1.0_RP + nx(i)) * delVar * nx(i)
    end do

    return
  end subroutine cal_del_flux_aux  

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine evaluate_error(tsec)
    use scale_polynominal, only: Polynominal_genLegendrePoly
    implicit none
    real(DP), intent(in) :: tsec

    real(RP) :: q_intrp(PolyOrderErrorCheck)
    real(RP) :: qexact_intrp(PolyOrderErrorCheck)
    real(RP) :: x_uwind(refElem%Np)
    real(RP) :: x_uwind_intrp(PolyOrderErrorCheck)
    real(RP) :: pos_intrp(PolyOrderErrorCheck)
    real(RP) vx(2)

    real(RP) :: l2error   
    real(RP) :: linferror   
    !------------------------------------------------------------------------

    l2error   = 0.0_RP   
    linferror = 0.0_RP      
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE

        x_uwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,k,1), ADV_VEL, tsec, dom_xmin, dom_xmax)

        vx(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),1)
        pos_intrp(:) = vx(1) + 0.5_RP*(x_intrp(:) + 1.0_RP)*(vx(2) - vx(1))
        x_uwind_intrp(:) = get_upwind_pos1d(pos_intrp(:), ADV_VEL, tsec, dom_xmin, dom_xmax)

        call get_profile1d_tracer( qexact%local(n)%val(:,k),    & ! (out)
          InitShapeName, x_uwind, InitShapeParams, refElem%Np )   ! (in)

        call get_profile1d_tracer( qexact_intrp(:),                              & ! (out)
          InitShapeName, x_uwind_intrp, InitShapeParams(:), PolyOrderErrorCheck)   ! (in)
        
        q_intrp(:) = matmul(IntrpMat, q%local(n)%val(:,k)) 
        qexact_intrp(:) = qexact_intrp(:) * exp(- DIFF_COEF * (2.0_RP*PI)**2 * tsec) 

        l2error = l2error &
          + sum(   lcmesh%J(1,k) * intw_intrp(:) * ( q_intrp(:) - qexact_intrp(:) )**2 ) 
        
        linferror = max(linferror, maxval(abs(q%local(n)%val(:,k) - qexact%local(n)%val(:,k))))
      end do
    end do

    LOG_INFO("evaluate_error_l2",*) sqrt(l2error)/(dom_xmax - dom_xmin) 
    LOG_INFO("evaluate_error_linf",*) linferror

  end subroutine evaluate_error

  subroutine set_initcond()
    implicit none

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE
        call get_profile1d_tracer( q%local(n)%val(:,k),                      & ! (out)
          InitShapeName, lcmesh%pos_en(:,k,1), InitShapeParams, refElem%Np )   ! (in)   
             
        qexact%local(n)%val(:,k) = q%local(n)%val(:,k)
        u%local(n)%val(:,k) = ADV_VEL
      end do
    end do

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
  
    integer :: NeGX = 8
    integer, parameter :: NLocalMeshPerPrc = 1
    integer            :: PolyOrder
    logical, parameter :: DumpedMassMatFlag = .false.
    character(len=H_SHORT) :: TINTEG_SCHEME_TYPE

    namelist /PARAM_TEST/ &
      NeGX, PolyOrder,                &
      TINTEG_SCHEME_TYPE,             &
      InitShapeName, InitShapeParams, &
      ADV_VEL,                        &
      DIFF_COEF,                      &
      PolyOrderErrorCheck,            &
      nstep_eval_error
        
    character(len=H_LONG) :: cnf_fname  ! config file for launcher

    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr
    !----------------------------------------------
    
    call PRC_MPIstart( comm )
     
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    cnf_fname = IO_ARG_getfname( ismaster )
    call IO_setup( "test_advdiff1d", cnf_fname )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   

    !--- read namelist

    NeGX = 2; PolyOrder = 1 
    InitShapeName    = 'sin'; 
    InitShapeParams  = (/ 1.0_RP, 0.0_RP /)
    ADV_VEL          = 0.0_RP
    DIFF_COEF        = 0.0_RP
    PolyOrderErrorCheck = 6
    nstep_eval_error    = 5
    TINTEG_SCHEME_TYPE = 'ERK_SSP_3s3o'
    
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
    call refElem%Init(PolyOrder, DumpedMassMatFlag)
    call Dx%Init(refElem%Dx1)
    call Sx%Init(refElem%Sx1)
    call Lift%Init(refElem%Lift)

    call mesh%Init( &
      NeGX,                              &
      dom_xmin, dom_xmax,                &
      refElem, NLocalMeshPerPrc )
    
    call mesh%Generate()

    !---
    call q%Init( "q", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )
    call dqdx%Init( "q", "1", mesh )
    call u%Init( "u", "m/s", mesh )

    call fields_comm%Init(2, 0, mesh)
    call auxvars_comm%Init(1, 0, mesh)
    field_list(1)%field1d => q
    field_list(2)%field1d => u
    auxvars_list(1)%field1d => dqdx
      
    call FILE_HISTORY_meshfield_setup( mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='X')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='X')

    ! setup for time integrator
    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call tinteg_lc(n)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,  &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)  )
    end do

    !------------------

    ! calculate the information about GL points, weights, and matrix for interpolation

    allocate(IntrpMat(PolyOrderErrorCheck,PolyOrder+1))
    allocate(intw_intrp(PolyOrderErrorCheck), x_intrp(PolyOrderErrorCheck))

    IntrpMat(:,:) = refElem%GenIntGaussLegendreIntrpMat( PolyOrderErrorCheck,  & ! (in)
                                                         intw_intrp, x_intrp  )  ! (out)

    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    use scale_time_manager, only: TIME_manager_Final    
    implicit none

    call PROF_rapstart( "final", 1 )
    call FILE_HISTORY_meshfield_finalize()

    call q%Final()
    call qexact%Final()
    call dqdx%Final()
    call u%Final()

    call fields_comm%Final()
    call auxvars_comm%Final()
    call mesh%Final()
    
    call Dx%Final()
    call Sx%Final()
    call Lift%Final()
    call refElem%Final()
    
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()

    return
  end subroutine final

end program test_advdiff1d
