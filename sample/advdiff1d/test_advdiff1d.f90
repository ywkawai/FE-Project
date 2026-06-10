!-------------------------------------------------------------------------------
!> Program A sample program: 1-dimensional linear advection and diffusion test
!! 
!! 
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
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

  use mod_advdiff1d_kernel, only: &
    advdiff1d_kernel_cal_aux, advdiff1d_kernel_cal_tend
  use mod_advdiff1d_numerror, only: AdvDiff1DNumErrorAnalysis
  !-----------------------------------------------------------------------------
  implicit none

  character(len=H_SHORT) :: InitShapeName   !< The type of initial profile (sin, gaussian-hill, cosine-bell, top-hat)
  real(RP), save :: InitShapeParams(2)
  integer :: InitGPMatPolyOrder
  real(RP) :: ADV_VEL                       !< The constant speed of advection
  real(RP) :: DIFF_COEF                     !< The coefficient of diffusion
  logical :: Do_NumErrorAnalysis            !< Flag wheter analysis of numerical error is performed

  type(LineElement)  :: refElem
  type(sparsemat) :: Dx, Lift

  type(MeshLineDom1D), target :: mesh
  type(LocalMesh1D), pointer :: lcmesh
  integer :: domid

  type(MeshField1D), target :: q, qexact 
  type(MeshField1D), target :: dqdx
  type(MeshField1D), target :: u
  type(MeshFieldComm1D) :: fields_comm
  type(MeshFieldComm1D) :: auxvars_comm
  type(MeshFieldContainer), save :: field_list(2)
  type(MeshFieldContainer), save :: auxvars_list(1)

  integer, save :: HST_ID(3)

  real(RP) :: tsec_
  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1

  type(AdvDiff1DNumErrorAnalysis) :: numerror_analysis
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
        call advdiff1d_kernel_cal_aux( dqdx%local(domid)%val,    & ! (out)
          q%local(domid)%val, Dx, Lift, lcmesh, lcmesh%refElem1D ) ! (in)
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
        call advdiff1d_kernel_cal_tend( &
          tinteg_lc(domid)%tend_buf2D_ex(:,:,RKVAR_Q,tintbuf_ind),                  & ! (out)
          q%local(domid)%val, dqdx%local(domid)%val, u%local(domid)%val, DIFF_COEF, & ! (in)
          Dx, Lift, lcmesh, lcmesh%refElem1D )                                        ! (in)
        call PROF_rapend( 'cal_prg_tend', 1) 

        call PROF_rapstart( 'update_var', 1)
        call tinteg_lc(domid)%Advance( rkstage, q%local(domid)%val, RKVAR_Q,    &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call PROF_rapend('update_var', 1)
      end do
    end do
    !$acc update host( q%local(1)%val )
    write(*,*) "q=", q%local(1)%val(1:refElem%Np,lcmesh%NeS:lcmesh%NeE)

    tsec_ = TIME_DTSEC * real(TIME_NOWSTEP-1, kind=RP)
    if ( Do_NumErrorAnalysis ) then
      call numerror_analysis%Eval( qexact, & ! (out)
        q, TIME_NOWSTEP, tsec_             ) ! (in)
    end if

    !* Output history file

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), dqdx)
    call FILE_HISTORY_meshfield_put(HST_ID(3), qexact)
    call FILE_HISTORY_meshfield_write()

    if (TIME_DOend) exit
  end do

  call final()

contains
  !> Set inital data
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
      !$acc update device(q%local(idom)%val, u%local(idom)%val)
    end do

    if ( Do_NumErrorAnalysis ) then
      call numerror_analysis%Eval( qexact, & ! (out)
        q, 1, 0.0_RP                       ) ! (in)
    end if

    call FILE_HISTORY_meshfield_put( HST_ID(1), q )
    call FILE_HISTORY_meshfield_put( HST_ID(2), dqdx )
    call FILE_HISTORY_meshfield_put( HST_ID(3), qexact )
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
      DIFF_COEF,                      &
      Do_NumErrorAnalysis
    
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
    InitShapeName       = 'sin'; 
    InitShapeParams     = (/ 1.0_RP, 0.0_RP /)
    InitGPMatPolyOrder  = 7
    ADV_VEL             = 0.0_RP
    DIFF_COEF           = 0.05_RP
    TINTEG_SCHEME_TYPE  = 'ERK_SSP_3s3o'
    Do_NumErrorAnalysis = .false.
    
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

    !-- setup fields

    call q%Init( "q", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )
    call dqdx%Init( "dqdx", "1", mesh )
    call u%Init( "u", "m/s", mesh )

    !-- setup data communicators

    field_list(1)%field1d => q
    field_list(2)%field1d => u
    call fields_comm%Init( size(field_list), 0, mesh )

    auxvars_list(1)%field1d => dqdx
    call auxvars_comm%Init( size(auxvars_list), 0, mesh )
      
    !-- setup history files    

    call FILE_HISTORY_meshfield_setup( mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='X')
    call FILE_HISTORY_reg( dqdx%varname, "dqdx", dqdx%unit, HST_ID(2), dim_type='X')
    call FILE_HISTORY_reg( qexact%varname, "qexact", qexact%unit, HST_ID(3), dim_type='X')

    !-- setup for time integrator

    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do domid=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(domid)
      call tinteg_lc(domid)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,  &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)  )
    end do

    !-- setup a module for evaluating numerical errors 
    if ( Do_NumErrorAnalysis ) &
      call numerror_analysis%Init( ADV_VEL, DIFF_COEF, InitShapeName, InitShapeParams, mesh, refElem )

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
    if ( Do_NumErrorAnalysis ) &
      call numerror_analysis%Final()

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
