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

  use scale_sparsemat, only: SparseMat
  use scale_element_base, only: ElementBase1D
  use scale_element_line, only: LineElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_mesh_linedom1d, only: MeshLineDom1D

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
    TIME_manager_checkstate, TIME_manager_advance,     &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,        &
    TIME_DTSEC, TIME_NSTEP, TIME_DOresume, TIME_DOend

  use scale_timeint_rk, only: timeint_rk

  use mod_advdiff1d_numerror, only: advdiff1d_numerror_eval
  !-----------------------------------------------------------------------------
  implicit none

  character(len=H_SHORT) :: InitShapeName   !< The type of initial profile (sin, gaussian-hill, cosine-bell, top-hat)
  real(RP), save :: InitShapeParams(2)
  integer :: InitGPMatPolyOrder
  real(RP) :: ADV_VEL                       !< The constant speed of advection
  real(RP) :: DIFF_COEF                     !< The coefficient of diffusion

  type(LineElement)  :: refElem
  type(sparsemat) :: Dx, Lift

  type(MeshLineDom1D), target :: mesh
  type(MeshField1D), target :: q, qexact 
  type(MeshField1D), target :: dqdx
  type(MeshField1D), target :: u
  type(MeshFieldComm1D) :: fields_comm
  type(MeshFieldComm1D) :: auxvars_comm
  type(MeshFieldContainer), save :: field_list(2)
  type(MeshFieldContainer), save :: auxvars_list(1)

  integer, save :: HST_ID(2)

  integer :: domid, k, p
  type(LocalMesh1D), pointer :: lcmesh
  
  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: nowstep
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1
  real(RP) :: tsec_
  !-------------------------------------------------------

  call init()

  do 
    !* Report current time
    call TIME_manager_checkstate

    if (TIME_DOresume) call set_initcond()

    !* Advance time
    call TIME_manager_advance()
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )

    do rkstage=1, tinteg_lc(1)%nstage

      !* Exchange halo data (prognostic variables)

      call PROF_rapstart( 'exchange_halo', 1)
      call fields_comm%Put(field_list, 1)
      call fields_comm%Exchange()
      call fields_comm%Get(field_list, 1)
      call PROF_rapend( 'exchange_halo', 1)

      !* Update auxiliary variables

      do domid=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(domid)

        call PROF_rapstart( 'cal_prg_tend', 1)
        call cal_aux( dqdx%local(domid)%val, &
          q%local(domid)%val,                &
          lcmesh, lcmesh%refElem1D       )
        call PROF_rapend( 'cal_prg_tend', 1) 
      end do

      !* Exchange halo data (auxiliary variables)

      call PROF_rapstart( 'exchange_halo_aux', 1)
      call auxvars_comm%Put(auxvars_list, 1)
      call auxvars_comm%Exchange()
      call auxvars_comm%Get(auxvars_list, 1)
      call PROF_rapend( 'exchange_halo_aux', 1)

      !* Update prognostic variables

      do domid=1, mesh%LOCAL_MESH_NUM

        lcmesh => mesh%lcmesh_list(domid)
        tintbuf_ind = tinteg_lc(domid)%tend_buf_indmap(rkstage)      
        
        call PROF_rapstart( 'cal_prg_tend', 1)
        call cal_prg_tend( &
          tinteg_lc(domid)%tend_buf2D_ex(:,:,RKVAR_Q,tintbuf_ind), &
          q%local(domid)%val, dqdx%local(domid)%val, u%local(domid)%val,   &
          lcmesh, lcmesh%refElem1D )
        call PROF_rapend( 'cal_prg_tend', 1) 

        call PROF_rapstart( 'update_var', 1)
        call tinteg_lc(domid)%Advance( rkstage, q%local(domid)%val, RKVAR_Q,    &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call PROF_rapend('update_var', 1)
      end do
    end do
    
    tsec_ = TIME_DTSEC * real(TIME_NOWSTEP-1, kind=RP)
    call advdiff1d_numerror_eval( qexact, & ! (out)
      q, TIME_NOWSTEP, tsec_, ADV_VEL, DIFF_COEF, InitShapeName, InitShapeParams, & ! (in)
      mesh, mesh%refElem1D                                                        ) ! (in)

    !* Output history file

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
    call FILE_HISTORY_meshfield_write()

    if (TIME_DOend) exit
  end do

  call final()

contains
  subroutine cal_prg_tend( dqdt, q_, dqdx_, u_, lmesh, elem)
    use scale_sparsemat, only: sparsemat_matmul
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dqdx_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)

    integer :: ke
    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 2)
    call cal_del_flux_prg( del_flux,           & ! (out)
      q_, dqdx_, u_, lmesh%normal_fn(:,:,1),   & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem )    ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 2)

    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 2)
    !$omp parallel do private(Fx, LiftDelFlx)
    do ke = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,ke) * u_(:,ke) - DIFF_COEF * dqdx_(:,ke), Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx)

      dqdt(:,ke) = - (  lmesh%Escale(:,ke,1,1) * Fx(:) &
                     + LiftDelFlx(:) )
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine cal_prg_tend

  subroutine cal_del_flux_prg( del_flux, q_, dqdx_, u_, nx, vmapM, vmapP, lmesh, elem )    
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem  
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

      alpha = 0.5_RP * abs( u_(iP) + u_(iM) )
      del_flux(i) = 0.5_RP * (                          &
          ( q_(iP) * u_(iP) - q_(iM) * u_(iM) ) * nx(i) &
        - alpha * ( q_(iP) - q_(iM) )                   & 
        - DIFF_COEF * ( 1.0_RP - nx(i) )                &
           * ( dqdx_(iP) - dqdx_(iM) ) * nx(i)          )
    end do

    return
  end subroutine cal_del_flux_prg

  subroutine cal_aux( dqdx_, q_, lmesh, elem)
    use scale_sparsemat, only: sparsemat_matmul
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dqdx_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)

    integer :: ke
    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_dyn_tend_bndflux', 2)
    call cal_bnd_flux_aux( del_flux,           & ! (out)
      q_, lmesh%normal_fn(:,:,1),              & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem )    ! (in)
    call PROF_rapend( 'cal_dyn_tend_bndflux', 2)

    !-----
    call PROF_rapstart( 'cal_dyn_tend_interior', 2)
    !$omp parallel do private(Fx, LiftDelFlx)
    do ke = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,ke), Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke), LiftDelFlx)

      dqdx_(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) &
                 + LiftDelFlx(:)
    end do
    call PROF_rapend( 'cal_dyn_tend_interior', 2)

    return
  end subroutine cal_aux

  subroutine cal_bnd_flux_aux( del_flux, q_, nx, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: delVar
    !------------------------------------------------------------------------

    !$omp parallel do private(i, iM, iP, delVar)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)
      delVar = 0.5_RP * ( q_(iP) - q_(iM) )
      del_flux(i) = ( 1.0_RP + nx(i) ) * delVar * nx(i)
    end do

    return
  end subroutine cal_bnd_flux_aux  

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_initcond()
    use mod_fieldutil, only: fieldutil_get_profile1d_tracer 
    implicit none

    class(LocalMeshBase), pointer :: lmesh
    integer :: idom
    integer :: ke

    type(LineElement) :: intrpElem
    real(RP) :: InvV_intrp(refElem%PolyOrder+1,InitGPMatPolyOrder+1)    
    real(RP) :: GPMat(refElem%PolyOrder+1,InitGPMatPolyOrder+1)

    real(RP) :: pos_intrp(InitGPMatPolyOrder+1)
    real(RP) :: vx(2)

    real(RP) :: q_intrp(InitGPMatPolyOrder+1)
    !--------------------------------

    call intrpElem%Init( InitGPMatPolyOrder, .false. )
    InvV_intrp(:,:) = intrpElem%invV(1:refElem%PolyOrder+1,:)
    GPMat(:,:) = matmul( refElem%V, InvV_intrp )

    do idom=1, mesh%LOCAL_MESH_NUM
      call mesh%GetLocalMesh( idom, lmesh )
      do ke=lmesh%NeS, lmesh%NeE
        vx(:) = lmesh%pos_ev(lmesh%EToV(ke,:),1)
        pos_intrp(:) = vx(1) + 0.5_RP*( intrpElem%x1(:) + 1.0_RP ) * ( vx(2) - vx(1) )

        call fieldutil_get_profile1d_tracer( q_intrp(:),             & ! (out)
          InitShapeName, pos_intrp(:), InitShapeParams, intrpElem%Np ) ! (in)   
             
        q%local(idom)%val(:,ke) = matmul( GPMat, q_intrp )
        u%local(idom)%val(:,ke) = ADV_VEL
      end do
    end do

    call advdiff1d_numerror_eval( qexact, & ! (out)
      q, 1, 0.0_RP, ADV_VEL, DIFF_COEF, InitShapeName, InitShapeParams, & ! (in)
      mesh, mesh%refElem1D                                              ) ! (in)
  
    call FILE_HISTORY_meshfield_put( HST_ID(1), q )
    call FILE_HISTORY_meshfield_put( HST_ID(2), qexact )
    call FILE_HISTORY_meshfield_write()   

    call intrpElem%Final()

    return
  end subroutine set_initcond

  !> Initialization
  subroutine init()
    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only:           &
      TIME_manager_Init,                    &
      TIME_manager_report_timeintervals
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_setup  
    use scale_file_history, only: FILE_HISTORY_reg       
    use mod_advdiff1d_numerror, only: advdiff1d_numerror_Init     
    implicit none
  
    real(RP), parameter :: dom_xmin =  0.0_RP
    real(RP), parameter :: dom_xmax = +1.0_RP
  
    integer :: NeGX                         = 8
    integer, parameter :: NLocalMeshPerPrc  = 1
    integer            :: PolyOrder         = 3
    logical, parameter :: LumpedMassMatFlag = .false.
    character(len=H_SHORT) :: TINTEG_SCHEME_TYPE

    namelist /PARAM_TEST/ &
      NeGX, PolyOrder,                &
      TINTEG_SCHEME_TYPE,             &
      InitShapeName, InitShapeParams, &
      InitGPMatPolyOrder,             &      
      ADV_VEL,                        &
      DIFF_COEF
    
    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr

    class(LocalMeshBase), pointer :: lmesh    
    integer :: idom

    character(len=H_LONG) :: cnf_fname  ! config file for launcher
    !----------------------------------------------
    
    !-- setup MPI

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
    InitShapeName      = 'sin'; 
    InitShapeParams    = (/ 1.0_RP, 0.0_RP /)
    InitGPMatPolyOrder = 7
    ADV_VEL            = 0.0_RP
    DIFF_COEF          = 0.05_RP
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

    !-- setup reference element and spatial operators

    call refElem%Init(PolyOrder, LumpedMassMatFlag)
    call Dx%Init(refElem%Dx1)
    call Lift%Init(refElem%Lift)

    !-- setup mesh

    call mesh%Init( NeGX, dom_xmin, dom_xmax, refElem, NLocalMeshPerPrc )
    call mesh%Generate()

    !-- seup fields

    call q%Init( "q", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )
    call dqdx%Init( "q", "1", mesh )
    call u%Init( "u", "m/s", mesh )

    !-- seup data communicators

    call fields_comm%Init(2, 0, mesh)
    call auxvars_comm%Init(1, 0, mesh)
    field_list(1)%field1d => q
    field_list(2)%field1d => u
    auxvars_list(1)%field1d => dqdx
      
    !-- setup history files    

    call FILE_HISTORY_meshfield_setup( mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='X')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='X')

    !-- setup for time integrator

    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do domid=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(domid)
      call tinteg_lc(domid)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,  &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)  )
    end do

    !-- setup a module for evaluating numerical errors 
    call advdiff1d_numerror_Init( refElem )

    !-- report information of time intervals
    call TIME_manager_report_timeintervals

    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  !> Finalization
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
    
    call Dx%Final(); call Lift%Final()
    call refElem%Final()
    
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()

    return
  end subroutine final

end program test_advdiff1d
