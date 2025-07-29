!-------------------------------------------------------------------------------
!> A sample program: post-processing for 2-dimensional DG data
!! 
!! 
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
program test_post_processing_2d
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
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_linedom1d, only: MeshLineDom1D
  use scale_mesh_rectdom2d, only: MeshRectDom2D

  use scale_localmeshfield_base, only: LocalMeshField1D
  use scale_meshfield_base, only: MeshField1D, MeshField2D

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

  use scale_element_SIACfilter, only: SIAC_filter
  use scale_meshfield_spectral_transform, only: &
    MeshField_SpetralTransform2D,   &
    ST_EVALTYPE_SAMPLE_UNIFORM_PTS, &
    ST_EVALTYPE_L2PROJECTION_1

  use netcdf
  implicit none
  !-----------------------------------------------------------------------------

  character(len=H_SHORT) :: InitShapeName   !< The type of initial profile (sin, gaussian-hill, cosine-bell, top-hat)
  real(RP), save :: InitShapeParams(4)
  integer :: InitGPMatPolyOrder
  real(RP) :: ADV_VEL                       !< The constant speed of advection
  logical :: Do_NumErrorAnalysis            !< Flag wheter analysis of numerical error is performed

  type(QuadrilateralElement)  :: refElem
  type(SparseMat) :: Dx, Dy, Lift

  type(MeshRectDom2D), target :: mesh
  type(LocalMesh2D), pointer :: lcmesh
  integer :: domid

  type(MeshField2D), target :: q, qexact  
  type(MeshField2D), target :: u, v
  integer, save :: HST_ID(2)
  
  real(RP) :: tsec_
  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_Q = 1

  type(SIAC_filter) :: filter
  integer :: NsamplePtTot_x
  integer :: NsamplePtTot_y
  integer :: NsamplePtPerElem1D  
  real(RP), allocatable :: sampling_x(:)
  real(RP), allocatable :: sampling_y(:)
  real(RP), allocatable :: IntrpMatSampling(:,:)

  type(MeshField_SpetralTransform2D) :: STtool
  integer :: ST_ks, ST_ke
  integer :: ST_ls, ST_le

  integer :: ncido
  integer :: time_id
  integer :: q_ori_id
  integer :: q_exact_id
  integer :: filtered_q_id

  integer :: st_q_r_id
  integer :: st_q_i_id

  integer :: istep
  !-----------------------------------------------------------------------------

  call init()
  !-
  call set_profile()

  istep = 1
  call nc_check( nf90_put_var( ncido, time_id, (/ istep * 0.0_RP /), &
    start=[istep], count=[1]) ) 
       
  ! call apply_SIAC_filter()
  call apply_spectral_transform()

  !-
  call final()

contains
!OCL SERIAL
  subroutine apply_SIAC_filter()
    use mod_fieldutil, only: fieldutil_get_profile1d_tracer
    implicit none

    real(RP), allocatable :: q_ori(:,:,:)
    real(RP), allocatable :: q_out_ori(:,:)
    real(RP), allocatable :: q_out_exact(:,:)
    real(RP), allocatable :: filtered_q(:,:)

    integer, parameter :: NmeshHalo = 1
    integer :: ke
    integer :: ldomID

    integer :: meshID
    !---------------------------------------

    ldomID = 1
    lcmesh => mesh%lcmesh_list(ldomID)

    allocate( q_ori(refElem%Np,lcmesh%Ne,1+2*NmeshHalo) )
    allocate( q_out_exact(filter%Npts_per_elem,lcmesh%Ne) )
    allocate( q_out_ori(filter%Npts_per_elem,lcmesh%Ne) )
    allocate( filtered_q(filter%Npts_per_elem,lcmesh%Ne) )

    call fieldutil_get_profile1d_tracer( q_out_exact,          & ! (out)
      InitShapeName, sampling_x, InitShapeParams, NsamplePtTot_x ) ! (in)   

    meshID = NmeshHalo + 1
    do ke=1, lcmesh%Ne
      q_ori(:,ke,meshID) = q%local(ldomID)%val(:,ke)
      q_out_ori(:,ke) = matmul( IntrpMatSampling, q_ori(:,ke,meshID) )
    end do
    do ke=1, lcmesh%Ne
      q_ori(:,ke,meshID-1) = q_ori(:,ke,meshID)
      q_ori(:,ke,meshID+1) = q_ori(:,ke,meshID)
    end do
    
    call filter%Apply1D( filtered_q, &
      q_ori, refElem%Np, lcmesh%Ne, size(q_ori,3), NmeshHalo )
  
    !-
    call output_filtered_q( filtered_q, q_out_ori, q_out_exact )
    return
  end subroutine apply_SIAC_filter
  subroutine output_filtered_q( filtered_q, q_ori, q_exact )
    implicit none
    real(RP), intent(in) :: filtered_q(NsamplePtTot_x,NsamplePtTot_y)
    real(RP), intent(in) :: q_ori(NsamplePtTot_x,NsamplePtTot_y)
    real(RP), intent(in) :: q_exact(NsamplePtTot_x,NsamplePtTot_y)
    !------------------------------------------

    call nc_check( nf90_put_var( ncido, q_ori_id, q_ori, &
      start=[1,1,istep], count=[NsamplePtTot_x,NsamplePtTot_y,1]) )    
    call nc_check( nf90_put_var( ncido, q_exact_id, q_exact, &
      start=[1,1,istep], count=[NsamplePtTot_x,NsamplePtTot_y,1]) )    
    call nc_check( nf90_put_var( ncido, filtered_q_id, filtered_q, &
      start=[1,1,istep], count=[NsamplePtTot_x,NsamplePtTot_y,1]) )
    return
  end subroutine output_filtered_q

!OCL SERIAL
  subroutine apply_spectral_transform()
    use scale_meshfield_base, only: MeshField2DList
    implicit none

    type(MeshField2DList) :: target_var_list(1,1)
    !-----------------------
    target_var_list(1,1)%ptr => q
    
    call STtool%Transform( target_var_list, 1, 1 )
    call output_spectral_data( STtool%spectral_coef, 1 )
    return
  end subroutine apply_spectral_transform
  subroutine output_spectral_data( s_q, var_num )
    implicit none
    integer, intent(in) :: var_num
    real(RP), intent(in) :: s_q(STtool%ks:STtool%ke,STtool%ls:STtool%le,2,var_num)
    !------------------------------------------
    call nc_check( nf90_put_var( ncido, st_q_r_id, s_q(:,:,1,1), &
      start=[1,1,istep], count=[STtool%kall,STtool%lall,1]) )    
    call nc_check( nf90_put_var( ncido, st_q_i_id, s_q(:,:,2,1), &
      start=[1,1,istep], count=[STtool%kall,STtool%lall,1]) )    
    return
  end subroutine output_spectral_data

  !> Set inital data
  subroutine set_profile()
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
        u%local(idom)%val(:,ke) = 0.0_RP
        v%local(idom)%val(:,ke) = 0.0_RP
      end do
    end do
  
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
!    use mod_advect1d_numerror, only: advect1d_numerror_Init     
    use scale_polynominal, only: Polynominal_GenLagrangePoly
    implicit none
    real(RP), parameter :: dom_xmin =   0.0_RP
    real(RP), parameter :: dom_xmax = + 1.0_RP
    real(RP), parameter :: dom_ymin =   0.0_RP
    real(RP), parameter :: dom_ymax = + 1.0_RP

    integer            :: NprcX               = 2
    integer            :: NprcY               = 2
    integer            :: NeGX                = 2
    integer            :: NeGY                = 2
    integer            :: PolyOrder           = 1
    logical, parameter :: LumpedMassMatFlag   = .false.
    character(len=H_SHORT) :: TINTEG_SCHEME_TYPE

    integer :: ST_GPMatPolyOrder
    character(len=H_MID) :: ST_TYPE ! L2_Projection, Uniform_Sampling
    integer :: ST_typeID

    namelist /PARAM_TEST/ &
      NprcX, NprcY,                   &
      NeGX, NeGY, PolyOrder,          &
      TINTEG_SCHEME_TYPE,             &
      InitShapeName, InitShapeParams, &
      InitGPMatPolyOrder,             &
      ADV_VEL,                        &
      Do_NumErrorAnalysis,            &
      NsamplePtPerElem1D,             &
      ST_GPMatPolyOrder,              &
      ST_TYPE

    integer :: comm, myrank, nprocs
    logical :: ismaster
    integer :: ierr
    
    class(LocalMesh2D), pointer :: lmesh    
    integer :: idom

    character(len=H_MID) :: conf_name

    type(LineElement) :: tmp_elem1D
    real(RP), allocatable :: sampling_1d_per_elem(:)

    integer :: p, ke
    real(RP) :: delx, dely

    integer :: sampling_x_dim_id
    integer :: sampling_x_id
    integer :: sampling_y_dim_id
    integer :: sampling_y_id
    integer :: time_dim_id
    integer :: st_k_dim_id
    integer :: st_k_id
    integer :: st_l_dim_id
    integer :: st_l_id
    character(len=6) :: pe_str
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
    InitShapeParams(:) = (/ 1.0_RP, 1.0_RP, 0.0_RP, 0.0_RP /)
    InitGPMatPolyOrder = 7
    ADV_VEL            = 1.0_RP
    TINTEG_SCHEME_TYPE = 'ERK_SSP_3s3o'
    Do_NumErrorAnalysis = .false.
    NsamplePtPerElem1D = 3
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
    call Dy%Init(refElem%Dx2)
    call Lift%Init(refElem%Lift)

    !-- setup mesh

    call mesh%Init( &
      NeGX, NeGY,                             &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
      .true., .true.,                         &
      refElem, 1, NprcX, NprcY )
    call mesh%Generate()

    !-- seup fields

    call q%Init( "q", "1", mesh )
    call u%Init( "u", "m/s", mesh )
    call v%Init( "v", "m/s", mesh )
    call qexact%Init( "qexact", "1", mesh )
    
    !-- setup history files    

    call FILE_HISTORY_meshfield_setup( mesh2D_=mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='XY')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='XY')

    !-- setup for time integrator

    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do idom=1, mesh%LOCAL_MESH_NUM
      lmesh => mesh%lcmesh_list(idom)
      call tinteg_lc(idom)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,  &
                               2, (/ lmesh%refElem%Np, lmesh%NeA /)  )
    end do

    !-- setup a module for evaluating numerical errors 
    ! if ( Do_NumErrorAnalysis ) &
    !   call advect1d_numerror_Init( refElem )

    !-- report information of time intervals
    call TIME_manager_report_timeintervals
    
    !----
    lmesh => mesh%lcmesh_list(1)

    call tmp_elem1D%Init( refElem%PolyOrder, .false. )

    allocate( sampling_1d_per_elem(NsamplePtPerElem1D) )
    allocate( sampling_x(NsamplePtPerElem1D*lmesh%NeX) )
    allocate( sampling_y(NsamplePtPerElem1D*lmesh%NeY) )

    NsamplePtTot_x = size(sampling_x)
    NsamplePtTot_y = size(sampling_y)

    delx = 2.0_RP / real(NsamplePtPerElem1D, kind=RP)
    do p=1, NsamplePtPerElem1D
      sampling_1d_per_elem(p) = - 1.0_RP + ( p - 0.5_RP )*delx
    end do

    allocate( IntrpMatSampling(NsamplePtPerElem1D,refElem%Nfp) )
    IntrpMatSampling(:,:) = Polynominal_GenLagrangePoly( refElem%PolyOrder, refElem%x1, sampling_1d_per_elem(:) )

    delx = ( dom_xmax - dom_xmin ) / real(lmesh%NeX, kind=RP)
    do ke=1, lmesh%NeX
    do p=1, NsamplePtPerElem1D
      sampling_x(p+(ke-1)*NsamplePtPerElem1D) = dom_xmin + ( ke-1 ) * delx + 0.5_RP * delx * ( 1.0_RP + sampling_1d_per_elem(p) )
    end do
    end do

    dely = ( dom_ymax - dom_ymin ) / real(lmesh%NeY, kind=RP)
    do ke=1, lmesh%NeY
    do p=1, NsamplePtPerElem1D
      sampling_y(p+(ke-1)*NsamplePtPerElem1D) = dom_ymin + ( ke-1 ) * dely + 0.5_RP * dely * ( 1.0_RP + sampling_1d_per_elem(p) )
    end do
    end do

    call filter%Init( 2*tmp_elem1D%PolyOrder, tmp_elem1D%PolyOrder+1, sampling_1d_per_elem, tmp_elem1D )
    call tmp_elem1D%Final()

    !-
    ST_ke = refElem%Nfp * NeGX / 2
    ST_ks = - ST_ke
    ST_le = refElem%Nfp * NeGY / 2
    ST_ls = - ST_le
    select case(trim(ST_TYPE))
    case( 'L2_Projection')
      ST_typeID = ST_EVALTYPE_L2PROJECTION_1
      call STtool%Init( ST_typeID, ST_ks, ST_ke, ST_ls, ST_le, mesh, 1, GLQuadOrd=ST_GPMatPolyOrder )
    case( 'Uniform_Sampling' )
      ST_typeID = ST_EVALTYPE_SAMPLE_UNIFORM_PTS
      call STtool%Init( ST_typeID, ST_ks, ST_ke, ST_ls, ST_le, mesh, 1, NsamplePtPerElem1D=refElem%Nfp )
    end select

    !--
    write(pe_str,'(I6.6)') PRC_myrank
    call nc_check( nf90_create("filtered_q.pe"//pe_str//".nc", NF90_clobber, ncido) )
    call nc_check( nf90_def_dim(ncido, 'x', NsamplePtTot_x, sampling_x_dim_id) )
    call nc_check( nf90_def_dim(ncido, 'y', NsamplePtTot_y, sampling_y_dim_id) )
    call nc_check( nf90_def_dim(ncido, 'k', size(STtool%k), st_k_dim_id) )
    call nc_check( nf90_def_dim(ncido, 'l', size(STtool%l), st_l_dim_id) )
    call nc_check( nf90_def_dim( ncido, 'time', NF90_UNLIMITED, time_dim_id ) )

    call nc_check( nf90_def_var(ncido, 'x', NF90_DOUBLE, sampling_x_dim_id, sampling_x_id) )
    call nc_check( nf90_def_var(ncido, 'y', NF90_DOUBLE, sampling_y_dim_id, sampling_y_id) )
    call nc_check( nf90_def_var(ncido, 'k', NF90_DOUBLE, st_k_dim_id, st_k_id) )
    call nc_check( nf90_def_var(ncido, 'l', NF90_DOUBLE, st_l_dim_id, st_l_id) )
    call nc_check( nf90_def_var(ncido, 'time', NF90_DOUBLE, time_dim_id, time_id) )

    call nc_check( nf90_def_var(ncido, 'q', NF90_DOUBLE, (/ sampling_x_dim_id, sampling_y_dim_id, time_dim_id /), q_ori_id ) )
    call nc_check( nf90_def_var(ncido, 'qexact', NF90_DOUBLE, (/ sampling_x_dim_id, sampling_y_dim_id, time_dim_id /), q_exact_id ) )
    call nc_check( nf90_def_var(ncido, 'filtered_q', NF90_DOUBLE, (/ sampling_x_dim_id, sampling_y_dim_id, time_dim_id /), filtered_q_id ) )

    call nc_check( nf90_def_var(ncido, 'st_q_r', NF90_DOUBLE, (/ st_k_dim_id, st_l_dim_id, time_dim_id /), st_q_r_id ) )
    call nc_check( nf90_def_var(ncido, 'st_q_i', NF90_DOUBLE, (/ st_k_dim_id, st_l_dim_id, time_dim_id /), st_q_i_id ) )

    call nc_check( nf90_enddef(ncido) )

    call nc_check( nf90_put_var(ncido, sampling_x_id, sampling_x ) )    
    call nc_check( nf90_put_var(ncido, sampling_y_id, sampling_y ) )
    call nc_check( nf90_put_var(ncido, st_k_id, STtool%k ) )    
    call nc_check( nf90_put_var(ncido, st_l_id, STtool%l ) )
    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  !> Finalization
  subroutine final()
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    use scale_time_manager, only: TIME_manager_Final 
!    use mod_advect1d_numerror, only: advect1d_numerror_Final   
    implicit none
    integer :: idom
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )
    ! if ( Do_NumErrorAnalysis ) &
    !   call advect1d_numerror_Final()

    call FILE_HISTORY_meshfield_finalize()

    call nc_check( nf90_close(ncido) )

    do idom=1, mesh%LOCAL_MESH_NUM
      call tinteg_lc(idom)%Final()
    end do

    call q%Final()
    call qexact%Final()
    call u%Final()
    call v%Final()

    call mesh%Final()
    
    call filter%Final()
    call STtool%Final()

    call Dx%Final(); call Dy%Final(); call Lift%Final()
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
    integer, intent (in) :: status
    implicit none
    if(status /= nf90_noerr) then 
      LOG_INFO('NetCDF_check',*) trim(nf90_strerror(status))
      call flush(IO_FID_CONF)
      call PRC_abort
    end if

    return
  end subroutine nc_check  
end program test_post_processing_2d
