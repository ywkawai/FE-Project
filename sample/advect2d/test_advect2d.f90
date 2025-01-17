!-------------------------------------------------------------------------------
!> Program A sample program: 1-dimensional linear advection test
!! 
!! 
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
program test_advect2d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof

  use scale_sparsemat, only: SparseMat
  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D

  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_put,   &
    FILE_HISTORY_meshfield_write
  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate

    use scale_time_manager, only: &
    TIME_manager_checkstate, TIME_manager_advance,     &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,        &
    TIME_DTSEC, TIME_NSTEP, TIME_DOresume, TIME_DOend
  
  use scale_timeint_rk, only: timeint_rk  

  use mod_advect2d_numerror, only: advect2d_numerror_eval
  
  use mod_fieldutil, only: &
    get_upwind_pos1d => fieldutil_get_upwind_pos1d,         &
    get_profile2d_tracer => fieldutil_get_profile2d_tracer, &
    get_profile2d_flow => fieldutil_get_profile2d_flow

  !-----------------------------------------------------------------------------
  implicit none

  character(len=H_SHORT) :: InitShapeName   !< The type of initial profile (sin, gaussian-hill, cosine-bell, top-hat)
  real(RP), save :: InitShapeParams(4)
  integer :: InitGPMatPolyOrder
  character(len=H_SHORT) :: VelTypeName     !< The type of specified velocify field (constant, rigid-body-rot)
  real(RP), save :: VelTypeParams(4)
  logical :: Do_NumErrorAnalysis            !< Flag wheter analysis of numerical error is performed
  
  type(QuadrilateralElement) :: refElem
  type(sparsemat) :: Dx, Dy, Lift
  
  type(MeshRectDom2D), target :: mesh
  type(MeshField2D), target :: q, qexact  
  type(MeshField2D), target :: u, v
  type(MeshFieldCommRectDom2D) :: fields_comm
  type(MeshFieldContainer), save :: field_list(3)  
  integer, save :: HST_ID(2)

  integer :: domid
  type(LocalMesh2D), pointer :: lcmesh
  
  integer :: nowstep
  real(RP) :: tsec_
  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1
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

      !* Exchange halo data

      call PROF_rapstart( 'exchange_halo', 1)
      call fields_comm%Put(field_list, 1)
      call fields_comm%Exchange()
      call fields_comm%Get(field_list, 1)
      call PROF_rapend( 'exchange_halo', 1)

      !* Set velocity field

      call PROF_rapstart( 'set_velocity', 1)
      call set_velocity( u, v, tsec_ )
      call PROF_rapend( 'set_velocity', 1)  

      !* Update prognostic variables
      
      do domid=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(domid)
        tintbuf_ind = tinteg_lc(domid)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'cal_tend', 1)
        call cal_tend( &
           tinteg_lc(domid)%tend_buf2D_ex(:,:,RKVAR_Q,tintbuf_ind),    & ! (out)
           q%local(domid)%val, u%local(domid)%val, v%local(domid)%val, & ! (in)
           lcmesh, lcmesh%refElem2D )                                    ! (in)
        call PROF_rapend( 'cal_tend', 1)

        call PROF_rapstart( 'update_var', 1)
        call tinteg_lc(domid)%Advance( rkstage, q%local(domid)%val, RKVAR_Q,    & ! (out)
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE ) ! (in)
        call PROF_rapend('update_var', 1)      
      end do
    end do
    
    tsec_ = TIME_DTSEC * real(TIME_NOWSTEP-1, kind=RP)
    if ( Do_NumErrorAnalysis ) then
      call advect2d_numerror_eval( qexact, & ! (out)
        q, TIME_NOWSTEP, tsec_, VelTypeName, VelTypeParams, InitShapeName, InitShapeParams, & ! (in)
        mesh, mesh%refElem2D                                                                ) ! (in)
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
  !> Calculate the tendency
  !! dqdt = - Dx ( uq ) - Dy ( vq ) + L ( < vec q>_numflx - vec q ).n
  !!
  subroutine cal_tend( dqdt, & ! (out)
    q_, u_, v_, lmesh, elem  ) ! (in)

    use scale_sparsemat, only: sparsemat_matmul
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: v_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), LiftBndFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)

    integer :: ke
    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_tend_bndflux', 2)
    call cal_elembnd_flux( del_flux,                              & ! (out)
      q_, u_, v_, lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                   & ! (in)
      lmesh, elem )                                                 ! (in)
    call PROF_rapend( 'cal_tend_bndflux', 2)

    !-----
    call PROF_rapstart( 'cal_tend_interior', 2)
    !$omp parallel do private(ke, Fx, Fy, LiftBndFlx)
    do ke = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,ke)*u_(:,ke), Fx)
      call sparsemat_matmul(Dy, q_(:,ke)*v_(:,ke), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke), LiftBndFlx)

      dqdt(:,ke) = - ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                     + lmesh%Escale(:,ke,2,2) * Fy(:) &
                     + LiftBndFlx(:) )
    end do
    call PROF_rapend( 'cal_tend_interior', 2)

    return
  end subroutine cal_tend

  subroutine cal_elembnd_flux( ebnd_flux, q_, u_, v_, nx, ny, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem  
    real(RP), intent(out) ::  ebnd_flux(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  v_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
     
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: VelP(elem%NfpTot), VelM(elem%NfpTot)
    real(RP) :: alpha(elem%NfpTot)

    integer :: ke
    !------------------------------------------------------------------------

    !$omp parallel do private(ke, iM, iP, VelM, VelP, alpha)
    do ke=1, lmesh%Ne
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)

      VelM(:) = u_(iM(:)) * nx(:,ke) + v_(iM(:)) * ny(:,ke)
      VelP(:) = u_(iP(:)) * nx(:,ke) + v_(iP(:)) * ny(:,ke)

      alpha(:) = 0.5_RP * abs( VelM(:) + VelP(:) )
      ebnd_flux(:,ke) = 0.5_RP * ( &
          ( q_(iP(:)) * VelP(:) - q_(iM(:)) * VelM(:) ) &
         - alpha(:) * ( q_(iP(:)) - q_(iM(:)) )         )
    end do

    return
  end subroutine cal_elembnd_flux

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_velocity( u_, v_, tsec )
    implicit none
    type(MeshField2D), intent(inout) :: u_
    type(MeshField2D), intent(inout) :: v_ 
    real(RP), intent(in) :: tsec

    integer :: idom, ke
    !----------------------------------------

    VelTypeParams(4) = tsec

    do idom=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(idom)
      !$omp parallel do private(ke)
      do ke=lcmesh%NeS, lcmesh%NeE
        call get_profile2d_flow( u%local(idom)%val(:,ke), v%local(idom)%val(:,ke),                   & ! (out)
          VelTypeName, lcmesh%pos_en(:,ke,1), lcmesh%pos_en(:,ke,2), VelTypeParams, refElem%Np )       ! (in)
      end do
    end do
    return
  end subroutine set_velocity

  subroutine set_initcond()
    use mod_fieldutil, only: fieldutil_get_profile2d_tracer 
    implicit none

    class(LocalMeshBase), pointer :: lmesh
    integer :: idom
    integer :: ke
    integer :: p, p_, p1, p2

    type(QuadrilateralElement) :: intrpElem
    real(RP) :: InvV_intrp(refElem%Np,(InitGPMatPolyOrder+1)**2)
    real(RP) :: GPMat(refElem%Np,(InitGPMatPolyOrder+1)**2)

    real(RP) :: pos_intrp((InitGPMatPolyOrder+1)**2,2)
    real(RP) :: vx(4), vy(4)

    real(RP) :: q_intrp((InitGPMatPolyOrder+1)**2)
    !------------------------------------------------------------------------

    call intrpElem%Init( InitGPMatPolyOrder, .false. )
    InvV_intrp(:,:) = 0.0_RP
    do p2=1, refElem%Nfp
    do p1=1, refElem%Nfp
      p = p1 + (p2-1)*refElem%Nfp
      p_ = p1 + (p2-1)*intrpElem%Nfp
      InvV_intrp(p,:) = intrpElem%invV(p_,:)
    end do
    end do
    GPMat(:,:) = matmul( refElem%V, InvV_intrp )

    do idom=1, mesh%LOCAL_MESH_NUM
      call mesh%GetLocalMesh( idom, lmesh )
      !$omp parallel do private(vx, vy, pos_intrp, q_intrp)
      do ke=lmesh%NeS, lmesh%NeE
        vx(:) = lmesh%pos_ev(lmesh%EToV(ke,:),1)
        vy(:) = lmesh%pos_ev(lmesh%EToV(ke,:),2)
        pos_intrp(:,1) = vx(1) + 0.5_RP*( intrpElem%x1(:) + 1.0_RP ) * ( vx(2) - vx(1) )
        pos_intrp(:,2) = vy(1) + 0.5_RP*( intrpElem%x2(:) + 1.0_RP ) * ( vy(3) - vy(1) )

        call fieldutil_get_profile2d_tracer( q_intrp(:),                               & ! (out)
          InitShapeName, pos_intrp(:,1), pos_intrp(:,2), InitShapeParams, intrpElem%Np ) ! (in)   
             
        q%local(idom)%val(:,ke) = matmul( GPMat, q_intrp )
      end do
    end do
    call set_velocity( u, v, 0.0_RP )

    if ( Do_NumErrorAnalysis ) then
      call advect2d_numerror_eval( qexact, & ! (out)
        q, 1, 0.0_RP, VelTypeName, VelTypeParams, InitShapeName, InitShapeParams, & ! (in)
        mesh, mesh%refElem2D                                                      ) ! (in)
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
    use mod_advect2d_numerror, only: advect2d_numerror_Init             
    implicit none

    real(RP), parameter :: dom_xmin =  0.0_RP
    real(RP), parameter :: dom_xmax = +1.0_RP
    real(RP), parameter :: dom_ymin =  0.0_RP
    real(RP), parameter :: dom_ymax = +1.0_RP

    integer :: NeGX                       
    integer :: NeGY                        
    integer :: PolyOrder                   
    integer, parameter :: NLocalMeshPerPrc = 1
    logical :: InitCond_GalerkinProjFlag         
    logical, parameter :: LumpedMassMatFlag = .false.
    character(len=H_SHORT) :: TINTEG_SCHEME_TYPE

    namelist /PARAM_TEST/ &
      NeGX, NeGY, PolyOrder,          &
      TINTEG_SCHEME_TYPE,             &
      InitShapeName, InitShapeParams, &
      InitCond_GalerkinProjFlag,      &
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
    call IO_setup( "test_advect2d", cnf_fname )
    
    !-- setup log
    call IO_LOG_setup( myrank, ismaster )   
  
    !--- read namelist

    NeGX = 2; NeGY = 2; PolyOrder = 1 
    TINTEG_SCHEME_TYPE = 'ERK_SSP_3s3o'
    InitShapeName      = 'sin'
    InitShapeParams(:) = (/ 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP /)
    VelTypeName        = 'const'
    InitCond_GalerkinProjFlag = .false.
    InitGPMatPolyOrder = 7
    VelTypeParams(:)   = (/ 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP /)
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
    
    call refElem%Init(PolyOrder, LumpedMassMatFlag)
    call Dx%Init(refElem%Dx1, storage_format='ELL')
    call Dy%Init(refElem%Dx2, storage_format='ELL')
    call Lift%Init(refElem%Lift, storage_format='ELL')

    !-- setup mesh

    call mesh%Init( &
      NeGX, NeGY,                             &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
      .true., .true.,                         &
      refElem, NLocalMeshPerPrc, 1, 1 )
    
    call mesh%Generate()
    
    !-- seup fields

    call q%Init( "q", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )
    call u%Init( "u", "m/s", mesh )
    call v%Init( "v", "m/s", mesh )
    
    !-- seup data communicators

    call fields_comm%Init(3, 0, 0, mesh)
    field_list(1)%field2d => q
    field_list(2)%field2d => u
    field_list(3)%field2d => v
  
    !-- setup history files    

    call FILE_HISTORY_meshfield_setup( mesh2d_=mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='XY')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='XY')
    
    !-- setup for time integrator

    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do domid=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(domid)
      call tinteg_lc(domid)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,        &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)  )
    end do

    !-- setup a module for evaluating numerical errors 
    if ( Do_NumErrorAnalysis ) &
      call advect2d_numerror_Init( refElem )

    !-- report information of time intervals
    call TIME_manager_report_timeintervals

    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    use scale_time_manager, only: TIME_manager_Final  
    use mod_advect2d_numerror, only: advect2d_numerror_Final   
    implicit none
    integer :: idom
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )
    call advect2d_numerror_Final()

    call FILE_HISTORY_meshfield_finalize()

    do idom=1, mesh%LOCAL_MESH_NUM
      call tinteg_lc(idom)%Final()
    end do

    call q%Final()
    call qexact%Final()
    call u%Final()
    call V%Final()

    call fields_comm%Final()
    call mesh%Final()
    
    call Dx%Final(); call Dy%Final(); call Lift%Final()
    call refElem%Final()
    
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()

    return
  end subroutine final
end program test_advect2d
