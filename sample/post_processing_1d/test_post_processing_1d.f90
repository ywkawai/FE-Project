!-------------------------------------------------------------------------------
!> A sample program: post-processing for 1-dimensional DG data
!! 
!! 
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
program test_post_processing_1d
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

  use mod_advect1d_numerror, only: advect1d_numerror_eval

  use scale_element_SIACfilter, only: SIAC_filter
  use scale_meshfield_spectral_transform, only: &
    MeshField_SpetralTransform1D, &
    ST_EVALTYPE_SAMPLE_UNIFORM_PTS, &
    ST_EVALTYPE_L2PROJECTION_1

  use netcdf
  implicit none
  !-----------------------------------------------------------------------------

  character(len=H_SHORT) :: InitShapeName   !< The type of initial profile (sin, gaussian-hill, cosine-bell, top-hat)
  real(RP), save :: InitShapeParams(2)
  integer :: InitGPMatPolyOrder
  real(RP) :: ADV_VEL                       !< The constant speed of advection
  logical :: Do_NumErrorAnalysis            !< Flag wheter analysis of numerical error is performed

  type(LineElement)  :: refElem
  type(sparsemat) :: Dx, Lift

  type(MeshLineDom1D), target :: mesh
  type(LocalMesh1D), pointer :: lcmesh
  integer :: domid

  type(MeshField1D), target :: q, u, qexact  
  type(MeshFieldComm1D) :: fields_comm
  type(MeshFieldContainer), save :: field_list(2)  
  integer, save :: HST_ID(2)
  
  real(RP) :: tsec_
  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1

  type(SIAC_filter) :: filter
  integer :: NsamplePtTot
  integer :: NsamplePtPerElem  
  real(RP), allocatable :: sampling_x(:)
  real(RP), allocatable :: IntrpMatSampling(:,:)

  type(MeshField_SpetralTransform1D) :: STtool
  integer :: ST_ks, ST_ke

  integer :: ncido
  integer :: sampling_x_dim_id
  integer :: sampling_x_id
  integer :: time_dim_id
  integer :: time_id
  integer :: q_ori_id
  integer :: q_exact_id
  integer :: filtered_q_id

  integer :: st_k_dim_id
  integer :: st_k_id
  integer :: st_q_r_id
  integer :: st_q_i_id
  !-----------------------------------------------------------------------------

  call init()
  !-
  call set_profile()
  call apply_SIAC_filter(1, 0.0_RP)
  call apply_spectral_transform(1, 0.0_RP)

  do
    !* Report current time
    call TIME_manager_checkstate

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

      !* Update prognostic variables

      do domid=1, mesh%LOCAL_MESH_NUM

        lcmesh => mesh%lcmesh_list(domid)
        tintbuf_ind = tinteg_lc(domid)%tend_buf_indmap(rkstage)      
        
        call PROF_rapstart( 'cal_tend', 1)
        call cal_tend( &
          tinteg_lc(domid)%tend_buf2D_ex(:,:,RKVAR_Q,tintbuf_ind), & ! (out)
          q%local(domid)%val, u%local(domid)%val,                  & ! (in)
          lcmesh, lcmesh%refElem1D )                                 ! (in)
        call PROF_rapend( 'cal_tend', 1) 

        call PROF_rapstart( 'update_var', 1)
        call tinteg_lc(domid)%Advance( rkstage, q%local(domid)%val, & ! (out) 
          RKVAR_Q, 1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE     ) ! (in)
        call PROF_rapend('update_var', 1)
      end do
    end do

    tsec_ = TIME_DTSEC * real(TIME_NOWSTEP-1, kind=RP)
    if ( Do_NumErrorAnalysis ) then
      call advect1d_numerror_eval( qexact, & ! (out)
        q, TIME_NOWSTEP, tsec_, ADV_VEL, InitShapeName, InitShapeParams, & ! (in)
        mesh, mesh%refElem1D                                             ) ! (in)
    end if

    !* Output history file

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
    call FILE_HISTORY_meshfield_write()

    if (TIME_DOend) exit
  end do
  LOG_INFO("tsec=",*) tsec_

  call apply_SIAC_filter(2, 1.0_RP)
  call apply_spectral_transform(2, 1.0_RP)

  !-
  call final()

contains

  !> Calculate the tendency
  !! dqdt = - Dx ( uq ) + L ( <u q>_numflx - uq )
  !!
  subroutine cal_tend( dqdt, & ! (out)
    q_, u_, lmesh, elem      ) ! (in)

    use scale_sparsemat, only: sparsemat_matmul
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), LiftBndFlux(elem%Np)
    real(RP) :: ebnd_flux(elem%NfpTot,lmesh%Ne)

    integer :: ke
    !------------------------------------------------------------------------

    call PROF_rapstart( 'cal_tend_bndflux', 2)
    call cal_elembnd_flux( ebnd_flux,        & ! (out)
      q_, u_, lmesh%normal_fn(:,:,1),        & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem  ) ! (in)
    call PROF_rapend( 'cal_tend_bndflux', 2)

    call PROF_rapstart( 'cal_tend_interior', 2)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul( Dx, q_(:,ke) * u_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * ebnd_flux(:,ke), LiftBndFlux )

      dqdt(:,ke) = - (  lmesh%Escale(:,ke,1,1) * Fx(:) &
                      + LiftBndFlux )
    end do
    call PROF_rapend( 'cal_tend_interior', 2)

    return
  end subroutine cal_tend

  !> Calculate the contribution at element boundaries: 
  !! 0.5 * [ ( [qu]^+ [qu]^- ) - ( [qu]^+ [qu]^- ) ] - [qu]^-
  subroutine cal_elembnd_flux( ebnd_flux,   & ! (out)
      q_, u_, nx, vmapM, vmapP, lmesh, elem ) ! (in)
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem  
    real(RP), intent(out) ::  ebnd_flux(elem%NfpTot,lmesh%Ne) !< Flux at element boundaries
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)          !< Normal vector at element boundaries
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (outside own element) into that of all nodes
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (inside) own element) into that of all nodes
     
    integer :: ke
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: uP(elem%NfpTot), uM(elem%NfpTot)
    real(RP) :: qP(elem%NfpTot), qM(elem%NfpTot)
    real(RP) :: alpha(elem%NfpTot)
    !------------------------------------------------------------------------

    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      uM(:) = u_(iM(:)); uP(:) = u_(iP(:))
      qM(:) = q_(iM(:)); qP(:) = q_(iP(:))

      alpha = 0.5_RP * abs( uP(:) + uM(:) )
      ebnd_flux(:,ke) = 0.5_RP * ( &  
          ( qP(:) * uP(:) - qM(:) * uM(:) ) * nx(:,ke) &
           - alpha(:) * ( qP(:) - qM(:) )              )  
    end do

    return
  end subroutine cal_elembnd_flux

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!OCL SERIAL
  subroutine apply_SIAC_filter(istep, output_dt)
    use mod_fieldutil, only: fieldutil_get_profile1d_tracer
    implicit none
    integer, intent(in) :: istep
    real(RP), intent(in) :: output_dt

    real(RP), allocatable :: q_ori(:,:,:)
    real(RP), allocatable :: q_out_ori(:,:)
    real(RP), allocatable :: q_out_exact(:,:)
    real(RP), allocatable :: filtered_q(:,:)

    integer :: NmeshHalo
    integer :: ke
    integer :: ldomID
    integer :: i

    integer :: meshID
    !---------------------------------------

    ldomID = 1
    lcmesh => mesh%lcmesh_list(ldomID)

    NmeshHalo = ceiling( (3 * refElem%PolyOrder + 1) / real(lcmesh%Ne, kind=RP) )

    allocate( q_ori(refElem%Np,lcmesh%Ne,1+2*NmeshHalo) )
    allocate( q_out_exact(filter%Npts_per_elem,lcmesh%Ne) )
    allocate( q_out_ori(filter%Npts_per_elem,lcmesh%Ne) )
    allocate( filtered_q(filter%Npts_per_elem,lcmesh%Ne) )

    call fieldutil_get_profile1d_tracer( q_out_exact,          & ! (out)
      InitShapeName, sampling_x, InitShapeParams, NsamplePtTot ) ! (in)   

    meshID = NmeshHalo + 1
    do ke=1, lcmesh%Ne
      q_ori(:,ke,meshID) = q%local(ldomID)%val(:,ke)
      q_out_ori(:,ke) = matmul( IntrpMatSampling, q_ori(:,ke,meshID) )
    end do
    do i=1, NmeshHalo
      do ke=1, lcmesh%Ne
        q_ori(:,ke,meshID-i) = q_ori(:,ke,meshID)
        q_ori(:,ke,meshID+i) = q_ori(:,ke,meshID)
      end do
    end do

    call filter%Apply1D( filtered_q, &
      q_ori, refElem%Np, lcmesh%Ne, size(q_ori,3), NmeshHalo )
    
    !-
    call output_filtered_q( istep, output_dt, filtered_q, q_out_ori, q_out_exact )
    return
  end subroutine apply_SIAC_filter
!OCL SERIAL
  subroutine output_filtered_q( istep, output_dt, filtered_q, q_ori, q_exact )
    implicit none
    integer, intent(in) :: istep
    real(RP), intent(in) :: output_dt
    real(RP), intent(in) :: filtered_q(NsamplePtTot)
    real(RP), intent(in) :: q_ori(NsamplePtTot)
    real(RP), intent(in) :: q_exact(NsamplePtTot)
    !------------------------------------------

    call nc_check( nf90_put_var( ncido, time_id, (/ (istep-1) * output_dt /), &
      start=[istep], count=[1]) )      

    call nc_check( nf90_put_var( ncido, q_ori_id, q_ori, &
      start=[1,istep], count=[NsamplePtTot,1]) )    
    call nc_check( nf90_put_var( ncido, q_exact_id, q_exact, &
      start=[1,istep], count=[NsamplePtTot,1]) )    
    call nc_check( nf90_put_var( ncido, filtered_q_id, filtered_q, &
      start=[1,istep], count=[NsamplePtTot,1]) )
    return
  end subroutine output_filtered_q

!OCL SERIAL
  subroutine apply_spectral_transform(istep, output_dt)
    use scale_meshfield_base, only: MeshField1DList
    implicit none
    integer, intent(in) :: istep
    real(RP), intent(in) :: output_dt

    type(MeshField1DList) :: target_var_list(1,1)
    !-----------------------
    target_var_list(1,1)%ptr => q
    
    call STtool%Transform( target_var_list, 1 )
    call output_spectral_data( istep, output_dt, STtool%spectral_coef, 1 )
    return
  end subroutine apply_spectral_transform
  subroutine output_spectral_data( istep, output_dt, s_q, var_num )
    implicit none
    integer, intent(in) :: istep
    real(RP), intent(in) :: output_dt
    integer, intent(in) :: var_num
    real(RP), intent(in) :: s_q(STtool%ks:STtool%ke,2,var_num)
    !------------------------------------------

    call nc_check( nf90_put_var( ncido, st_q_r_id, s_q(:,1,1), &
      start=[1,istep], count=[STtool%kall,1]) )    
    call nc_check( nf90_put_var( ncido, st_q_i_id, s_q(:,2,1), &
      start=[1,istep], count=[STtool%kall,1]) )    
    return
  end subroutine output_spectral_data

  !> Set inital data
  subroutine set_profile()
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
    !------------------------------------------------------------------------

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

    if ( Do_NumErrorAnalysis ) then
      call advect1d_numerror_eval( qexact, & ! (out)
        q, 1, 0.0_RP, ADV_VEL, InitShapeName, InitShapeParams, & ! (in)
        mesh, mesh%refElem1D                                   ) ! (in)
    end if

    call FILE_HISTORY_meshfield_put( HST_ID(1), q )
    call FILE_HISTORY_meshfield_put( HST_ID(2), qexact )
    call FILE_HISTORY_meshfield_write()   
  
    call intrpElem%Final()
    return
  end subroutine set_profile

  !> Initialization
  subroutine init()
    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only:           &
      TIME_manager_Init,                    &
      TIME_manager_report_timeintervals
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_setup  
    use scale_file_history, only: FILE_HISTORY_reg  
    use mod_advect1d_numerror, only: advect1d_numerror_Init     
    use scale_polynominal, only: Polynominal_GenLagrangePoly, Polynominal_GenLegendrePoly
    implicit none
    real(RP), parameter :: dom_xmin =   0.0_RP
    real(RP), parameter :: dom_xmax = + 1.0_RP
  
    integer            :: NeGX                = 2
    integer            :: PolyOrder           = 1
    logical, parameter :: LumpedMassMatFlag   = .false.
    character(len=H_SHORT) :: TINTEG_SCHEME_TYPE

    integer :: ST_GPMatPolyOrder
    character(len=H_MID) :: ST_TYPE ! L2_Projection, Uniform_Sampling
    integer :: ST_typeID
    
    namelist /PARAM_TEST/ &
      NeGX, PolyOrder,                &
      TINTEG_SCHEME_TYPE,             &
      InitShapeName, InitShapeParams, &
      InitGPMatPolyOrder,             &
      ADV_VEL,                        &
      Do_NumErrorAnalysis,            &
      NsamplePtPerElem,               &
      ST_GPMatPolyOrder,              &
      ST_TYPE

    
    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr
    
    class(LocalMeshBase), pointer :: lmesh    
    integer :: idom

    character(len=H_MID) :: conf_name

    type(LineElement) :: tmp_elem
    real(RP), allocatable :: sampling_x_per_elem(:)

    integer :: p, ke
    real(RP) :: delx

    real(RP), allocatable :: P_1D_ori(:,:)
    !----------------------------------------------

    !-- setup MPI
    call PRC_MPIstart( comm )

    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    !-- setup scale_io
    conf_name = IO_ARG_getfname( ismaster )
    call IO_setup( "test_advect1d", conf_name )
    
    !-- setup log
    call IO_LOG_setup( myrank, ismaster )   

    !-- read namelist

    InitShapeName      = 'sin'; 
    InitShapeParams    = (/ 1.0_RP, 0.0_RP /)
    InitGPMatPolyOrder = 7
    ADV_VEL            = 1.0_RP
    TINTEG_SCHEME_TYPE = 'ERK_SSP_3s3o'
    Do_NumErrorAnalysis = .false.
    NsamplePtPerElem   = 3
    ST_TYPE            = "L2_Projection"
    
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
    call Dx%Init(refElem%Dx1)
    call Lift%Init(refElem%Lift)

    !-- setup mesh

    call mesh%Init( NeGX, dom_xmin, dom_xmax, refElem, 1 )
    call mesh%Generate()

    !-- seup fields

    call q%Init( "q", "1", mesh )
    call u%Init( "u", "m/s", mesh )
    call qexact%Init( "qexact", "1", mesh )

    !-- seup data communicators

    field_list(1)%field1d => q
    field_list(2)%field1d => u  
    call fields_comm%Init( size(field_list), 0, mesh )
    
    !-- setup history files    

    call FILE_HISTORY_meshfield_setup( mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='X')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='X')

    !-- setup for time integrator

    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do idom=1, mesh%LOCAL_MESH_NUM
      call mesh%GetLocalMesh( idom, lmesh )
      call tinteg_lc(idom)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,  &
                               2, (/ lmesh%refElem%Np, lmesh%NeA /)  )
    end do

    !-- setup a module for evaluating numerical errors 
    if ( Do_NumErrorAnalysis ) &
      call advect1d_numerror_Init( refElem )

    !-- report information of time intervals
    call TIME_manager_report_timeintervals
    
    !----
    call mesh%GetLocalMesh( 1, lmesh )

    allocate( sampling_x_per_elem(NsamplePtPerElem), sampling_x(NsamplePtPerElem*lmesh%Ne) )
    NsamplePtTot = size(sampling_x)

    delx = 2.0_RP / real(NsamplePtPerElem, kind=RP)
    do p=1, NsamplePtPerElem
      sampling_x_per_elem(p) = - 1.0_RP + ( p - 0.5_RP )*delx
    end do

    allocate( IntrpMatSampling(NsamplePtPerElem,refElem%Np) )
    ! IntrpMatSampling(:,:) = Polynominal_GenLagrangePoly( refElem%PolyOrder, refElem%x1, sampling_x_per_elem(:) )
    allocate( P_1D_ori(NsamplePtPerElem,refElem%Np) )
    P_1D_ori(:,:) = Polynominal_GenLegendrePoly( refElem%PolyOrder, sampling_x_per_elem(:) )
    do p=1, refElem%Np
      P_1D_ori(:,p) = P_1D_ori(:,p) * sqrt(real(p-1,kind=RP) + 0.5_RP)
    end do
    IntrpMatSampling(:,:) = matmul( P_1D_ori(:,:), refElem%invV )

    delx = ( dom_xmax - dom_xmin ) / real(lmesh%Ne, kind=RP)
    do ke=1, lmesh%Ne
    do p=1, NsamplePtPerElem
      sampling_x(p+(ke-1)*NsamplePtPerElem) = dom_xmin + ( ke-1 ) * delx + 0.5_RP * delx * ( 1.0_RP + sampling_x_per_elem(p) )
    end do
    end do

    call tmp_elem%Init( refElem%PolyOrder, .false. )
    call filter%Init( 2*tmp_elem%PolyOrder, tmp_elem%PolyOrder+1, sampling_x_per_elem, tmp_elem )
    call tmp_elem%Final()
  
    !-
    ST_ke = refElem%Np * NeGX / 2
    ST_ks = - ST_ke
    select case(trim(ST_TYPE))
    case( 'L2_Projection')
      ST_typeID = ST_EVALTYPE_L2PROJECTION_1
      call STtool%Init( ST_typeID, ST_ks, ST_ke, mesh, 1, GLQuadOrd=ST_GPMatPolyOrder )
    case( 'Uniform_Sampling' )
      ST_typeID = ST_EVALTYPE_SAMPLE_UNIFORM_PTS
      call STtool%Init( ST_typeID, ST_ks, ST_ke, mesh, 1, NsamplePtPerElem=refElem%Np )
    end select

    !--
    call nc_check( nf90_create("filtered_q.nc", NF90_clobber, ncido) )
    call nc_check( nf90_def_dim(ncido, 'x', NsamplePtTot, sampling_x_dim_id) )
    call nc_check( nf90_def_dim(ncido, 'k', size(STtool%k), st_k_dim_id) )
    call nc_check( nf90_def_dim( ncido, 'time', NF90_UNLIMITED, time_dim_id ) )

    call nc_check( nf90_def_var(ncido, 'x', NF90_DOUBLE, sampling_x_dim_id, sampling_x_id) )
    call nc_check( nf90_def_var(ncido, 'k', NF90_DOUBLE, st_k_dim_id, st_k_id) )
    call nc_check( nf90_def_var(ncido, 'time', NF90_DOUBLE, time_dim_id, time_id) )

    call nc_check( nf90_def_var(ncido, 'q', NF90_DOUBLE, (/ sampling_x_dim_id, time_dim_id /), q_ori_id ) )
    call nc_check( nf90_def_var(ncido, 'qexact', NF90_DOUBLE, (/ sampling_x_dim_id, time_dim_id /), q_exact_id ) )
    call nc_check( nf90_def_var(ncido, 'filtered_q', NF90_DOUBLE, (/ sampling_x_dim_id, time_dim_id /), filtered_q_id ) )

    call nc_check( nf90_def_var(ncido, 'st_q_r', NF90_DOUBLE, (/ st_k_dim_id, time_dim_id /), st_q_r_id ) )
    call nc_check( nf90_def_var(ncido, 'st_q_i', NF90_DOUBLE, (/ st_k_dim_id, time_dim_id /), st_q_i_id ) )

    call nc_check( nf90_enddef(ncido) )

    call nc_check( nf90_put_var(ncido, sampling_x_id, sampling_x ) )    
    call nc_check( nf90_put_var(ncido, st_k_id, STtool%k ) )    
    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  !> Finalization
  subroutine final()
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    use scale_time_manager, only: TIME_manager_Final 
    use mod_advect1d_numerror, only: advect1d_numerror_Final   
    implicit none
    integer :: idom
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )
    if ( Do_NumErrorAnalysis ) &
      call advect1d_numerror_Final()

    call FILE_HISTORY_meshfield_finalize()

    call nc_check( nf90_close(ncido) )

    do idom=1, mesh%LOCAL_MESH_NUM
      call tinteg_lc(idom)%Final()
    end do

    call q%Final()
    call qexact%Final()
    call u%Final()

    call fields_comm%Final()
    call mesh%Final()
    
    call filter%Final()
    call STtool%Final()

    call Dx%Final(); call Lift%Final()
    call refElem%Final()
    
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()
    return
  end subroutine final  

!-----------------------------

!OCL SERIAL
  subroutine nc_check( status )
    implicit none
    integer, intent (in) :: status
    if(status /= nf90_noerr) then 
      LOG_INFO('NetCDF_check',*) trim(nf90_strerror(status))
      call flush(IO_FID_CONF)
      call PRC_abort
    end if

    return
  end subroutine nc_check  
end program test_post_processing_1d
