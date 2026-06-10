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

  use scale_sparsemat, only: SparseMat
  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D

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
    TIME_manager_checkstate, TIME_manager_advance,     &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,        &
    TIME_DTSEC, TIME_NSTEP, TIME_DOresume, TIME_DOend

  use scale_timeint_rk, only: &
    timeint_rk
  
  use mod_advect3d_kernel, only: advect3d_kernel_cal_tend
  use mod_advect3d_numerror, only: Advect3DNumErrorAnalysis  
  !-----------------------------------------------------------------------------
  implicit none

  character(len=H_SHORT) :: InitShapeName   !< The type of initial profile (sin, gaussian-hill, cosine-bell, top-hat)
  real(RP), save :: InitShapeParams(6)
  integer :: InitGPMatPolyOrder
  character(len=H_SHORT) :: VelTypeName     !< The type of specified velocify field (constant, rigid-body-rot)
  real(RP), save :: VelTypeParams(6)
  logical :: Do_NumErrorAnalysis            !< Flag wheter analysis of numerical error is performed
  
  type(HexahedralElement) :: refElem
  type(sparsemat) :: Dx, Dy, Dz, Lift
  
  type(MeshCubeDom3D), target :: mesh
  type(MeshField3D), target :: q, qexact  
  type(MeshField3D), target :: u, v, w
  type(MeshFieldCommCubeDom3D) :: fields_comm
  type(MeshFieldContainer), save :: field_list(4)  
  integer, save :: HST_ID(2)

  integer :: domid
  type(LocalMesh3D), pointer :: lcmesh

  integer :: nowstep
  real(RP) :: tsec_
  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1

  type(Advect3DNumErrorAnalysis) :: numerror_analysis
  !-------------------------------------------------------

  call init()

  tsec_ = 0.0_RP  
  do
    
    !* Report current time
    call TIME_manager_checkstate

    if (TIME_DOresume) call set_initcond()

    !* Advance time
    call TIME_manager_advance()
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )
      
    do rkstage=1, tinteg_lc(1)%nstage

      !* Set velocity field

      call PROF_rapstart( 'set_velocity', 1)
      call set_velocity( u, v, w, tsec_ )
      call PROF_rapend( 'set_velocity', 1)  

      !* Exchange halo data

      call PROF_rapstart( 'exchange_halo', 1)
      call fields_comm%Put(field_list, 1)
      call fields_comm%Exchange()
      call fields_comm%Get(field_list, 1)
      call PROF_rapend( 'exchange_halo', 1)

      !* Update prognostic variables

      do domid=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(domid)
        tintbuf_ind = tinteg_lc(domid)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'cal_dyn_tend', 1)
        call advect3d_kernel_cal_tend( &
           tinteg_lc(domid)%tend_buf2D_ex(:,:,RKVAR_Q,tintbuf_ind),                        & ! (out)
           q%local(domid)%val, u%local(domid)%val, v%local(domid)%val, w%local(domid)%val, & ! (in)
           Dx, Dy, Dz, Lift, lcmesh, lcmesh%refElem3D )                                      ! (in)
        call PROF_rapend( 'cal_dyn_tend', 1)

        call PROF_rapstart( 'update_var', 1)
        call tinteg_lc(domid)%Advance( rkstage, q%local(domid)%val, RKVAR_Q,              &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        call PROF_rapend('update_var', 1)
      end do
    end do

    tsec_ = TIME_DTSEC * real(TIME_NOWSTEP-1, kind=RP)
    if ( Do_NumErrorAnalysis ) then
      call numerror_analysis%Eval( qexact, & ! (inout)
        q, TIME_NOWSTEP, tsec_             ) ! (in)
    end if
    
    !* Output history file

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
    call FILE_HISTORY_meshfield_write()

    if( IO_L ) call flush(IO_FID_LOG)

    if (TIME_DOend) exit
  end do

  call final()

contains
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_velocity( u_, v_, w_, tsec )
    use mod_fieldutil, only: &
      get_profile3d_flow => fieldutil_get_profile3d_flow
    implicit none
    type(MeshField3D), intent(inout) :: u_
    type(MeshField3D), intent(inout) :: v_ 
    type(MeshField3D), intent(inout) :: w_ 
    real(RP), intent(in) :: tsec
    
    class(LocalMesh3D), pointer :: lmesh
    integer :: idom
    integer :: ke
    !----------------------------------------
    
    VelTypeParams(5) = tsec
    do idom=1, mesh%LOCAL_MESH_NUM
      lmesh => mesh%lcmesh_list(idom)
      call get_profile3d_flow( &
        u%local(idom)%val(:,lmesh%NeS:lmesh%NeE), v%local(idom)%val(:,lmesh%NeS:lmesh%NeE), w%local(idom)%val(:,lmesh%NeS:lmesh%NeE), & ! (out)
        VelTypeName, lmesh%pos_en(:,:,1), lmesh%pos_en(:,:,2), lmesh%pos_en(:,:,3),                                                   & ! (in)
        VelTypeParams, 1,refElem%Np,refElem%Np, lmesh%NeS,lmesh%NeE,lmesh%NeA, 1,1,1 )                                                  ! (in)
    end do    
    return
  end subroutine set_velocity

  subroutine set_initcond()
    use mod_fieldutil, only: &
      get_profile3d_tracer => fieldutil_get_profile3d_tracer
    use scale_linalgebra, only: linalgebra_inv
    implicit none

    class(LocalMesh3D), pointer :: lmesh
    integer :: idom
    integer :: ke
    integer :: p, p_, p1, p2, p3, p1_, p2_, p3_

    type(HexahedralElement) :: intrpElem
    real(RP) :: InvV_intrp(refElem%Np,(InitGPMatPolyOrder+1)**3)
    real(RP) :: GPMat(refElem%Np,(InitGPMatPolyOrder+1)**3)

    real(RP) :: pos_intrp((InitGPMatPolyOrder+1)**3,3)
    real(RP) vx(8), vy(8), vz(8)

    real(RP) :: q_intrp((InitGPMatPolyOrder+1)**3)
    !------------------------------------------------------------------------

    call intrpElem%Init( InitGPMatPolyOrder, InitGPMatPolyOrder, .false. )
    InvV_intrp(:,:) = 0.0_RP
    do p3=1, refElem%Nnode_v
    do p2=1, refElem%Nnode_h1D
    do p1=1, refElem%Nnode_h1D
      p = p1 + (p2-1)*refElem%Nnode_h1D + (p3-1)*refElem%Nnode_h1D**2
      p_ = p1 + (p2-1)*intrpElem%Nnode_h1D + (p3-1)*intrpElem%Nnode_h1D**2
      InvV_intrp(p,:) = intrpElem%invV(p_,:)
    end do
    end do
    end do
    GPMat(:,:) = matmul( refElem%V, InvV_intrp )

    do idom=1, mesh%LOCAL_MESH_NUM
      lmesh => mesh%lcmesh_list(idom)
      !$omp parallel do private(ke, vx, vy, vz, pos_intrp, q_intrp)
      do ke=lmesh%NeS, lmesh%NeE
        vx(:) = lmesh%pos_ev(lmesh%EToV(ke,:),1)
        vy(:) = lmesh%pos_ev(lmesh%EToV(ke,:),2)
        vz(:) = lmesh%pos_ev(lmesh%EToV(ke,:),3)
        pos_intrp(:,1) = vx(1) + 0.5_RP*( intrpElem%x1(:) + 1.0_RP ) * ( vx(2) - vx(1) )
        pos_intrp(:,2) = vy(1) + 0.5_RP*( intrpElem%x2(:) + 1.0_RP ) * ( vy(3) - vy(1) )
        pos_intrp(:,3) = vz(1) + 0.5_RP*( intrpElem%x3(:) + 1.0_RP ) * ( vz(5) - vz(1) )

        call get_profile3d_tracer( q_intrp,                               & ! (out)
          InitShapeName, pos_intrp(:,1), pos_intrp(:,2), pos_intrp(:,3),  & ! (in)
          InitShapeParams, intrpElem%Np )                                   ! (in)
        
        q%local(idom)%val(:,ke) = matmul( GPMat, q_intrp )
      end do
      !$acc update device( q%local(idom)%val )      
    end do
    call set_velocity( u, v, w, 0.0_RP )

    if ( Do_NumErrorAnalysis ) then
      call numerror_analysis%Eval( qexact, & ! (inout)
        q, 1, 0.0_RP             ) ! (in)
    end if

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
    call FILE_HISTORY_meshfield_write()   
  
    return
  end subroutine set_initcond

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
    real(RP), parameter :: dom_ymin =  0.0_RP
    real(RP), parameter :: dom_ymax = +1.0_RP
    real(RP), parameter :: dom_zmin =  0.0_RP
    real(RP), parameter :: dom_zmax = +1.0_RP

    integer :: NprcX
    integer :: NprcY
    integer :: NeX
    integer :: NeY
    integer :: NeGZ    
    integer :: PolyOrder_h
    integer :: PolyOrder_v
    integer, parameter :: NLocalMeshPerPrc = 1
    logical, parameter :: LumpedMassMatFlag = .false.
    character(len=H_SHORT) :: TINTEG_SCHEME_TYPE

    namelist /PARAM_TEST/ &
      NprcX, NeX, NprcY, NeY, NeGZ,   & 
      PolyOrder_h, PolyOrder_v,       &
      TINTEG_SCHEME_TYPE,             &
      InitShapeName, InitShapeParams, &
      InitGPMatPolyOrder,             &
      VelTypeName, VelTypeParams,     &
      Do_NumErrorAnalysis
    
    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr

    character(len=H_LONG) :: cnf_fname  ! config file for launcher
    !------------------------------------------------------------------------

    !-- setup MPI

    call PRC_MPIstart( comm )
     
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    !-- setup scale_io
    cnf_fname = IO_ARG_getfname( ismaster )
    call IO_setup( "test_advect3d", cnf_fname )
    
    !-- setup log
    call IO_LOG_setup( myrank, ismaster )   
  
    !--- read namelist

    NeX = 2; NeY = 2; NeGZ = 2
    PolyOrder_h = 1; PolyOrder_v = 1
    TINTEG_SCHEME_TYPE = 'ERK_SSP_3s3o'
    InitShapeName      = 'sin'
    InitShapeParams(:) = (/ 1.0_RP, 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP, 0.0_RP /)
    VelTypeName        = 'const'
    InitGPMatPolyOrder = max(PolyOrder_h, PolyOrder_v)
    VelTypeParams(:)   = (/ 1.0_RP, 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP, 0.0_RP /)
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

    !-- setup profiler

    call PROF_setup
    call PROF_rapstart( "total", 0 )
    call PROF_rapstart( "init", 1 )

    !-- setup calendar & initial time

    call CALENDAR_setup
    call TIME_manager_Init

    !-- setup reference element and spatial operators
    
    call refElem%Init(PolyOrder_h, PolyOrder_v, LumpedMassMatFlag)
    call Dx%Init(refElem%Dx1, storage_format='ELL')
    call Dy%Init(refElem%Dx2, storage_format='ELL')
    call Dz%Init(refElem%Dx3, storage_format='ELL')
    call Lift%Init(refElem%Lift, storage_format='ELL')

    !-- setup mesh

    call mesh%Init( &
      NeX*NprcX, NeY*NprcY, NeGZ,                                 &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
      .true., .true., .true.,                                     &
      refElem, NLocalMeshPerPrc, NprcX, NprcY )

    call mesh%Generate()
    
    !-- setup fields

    call q%Init( "q", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )
    call u%Init( "u", "m/s", mesh )
    call v%Init( "v", "m/s", mesh )
    call w%Init( "w", "m/s", mesh )
    
    !-- setup data communicators

    field_list(1)%field3d => q
    field_list(2)%field3d => u
    field_list(3)%field3d => v
    field_list(4)%field3d => w
    call fields_comm%Init( size(field_list), 0, 0, mesh)
    
    !-- setup history files    

    call FILE_HISTORY_meshfield_setup( mesh3d_=mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='XYZ')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='XYZ')

    !-- setup for time integrator

    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do domid=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(domid)
      call tinteg_lc(domid)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,        &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)  )
    end do

    !-- setup a module for evaluating numerical errors 
    if ( Do_NumErrorAnalysis ) then
      call numerror_analysis%Init( &
         VelTypeName, VelTypeParams, InitShapeName, InitShapeParams, & ! (in)
         mesh, refElem ) ! (in)
    end if

    !-- report information of time intervals
    call TIME_manager_report_timeintervals

    
    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    use scale_time_manager, only: TIME_manager_Final    
    implicit none
    integer :: idom
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )
    if ( Do_NumErrorAnalysis ) &
      call numerror_analysis%Final()
    
    call FILE_HISTORY_meshfield_finalize()

    do domid=1, mesh%LOCAL_MESH_NUM
      call tinteg_lc(domid)%Final()
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
