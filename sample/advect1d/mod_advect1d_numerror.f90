#include "scalelib.h"
module mod_advect1d_numerror
  use scale_precision
  use scale_prc
  use scale_io

  use scale_element_base, only: ElementBase1D
  use scale_element_line, only: LineElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_mesh_linedom1d, only: MeshLineDom1D
  use scale_meshfield_base, only: MeshField1D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  public :: advect1d_numerror_Init
  public :: advect1d_numerror_do_step
  public :: advect1d_numerror_eval
  public :: advect1d_numerror_Final

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  integer :: polyOrderErrorCheck
  integer :: LOG_STEP_INTERVAL
  integer :: LOG_RSTEP_INTERVAL

  real(RP), allocatable :: IntrpMat(:,:)
  real(RP), allocatable :: intw_intrp(:)
  real(RP), allocatable :: x_intrp(:)

  integer               :: NUMERROR_LOG_FID = -1   !< fileID for monitor output file

contains
  !> Initialization
  subroutine advect1d_numerror_Init( elem )
    use scale_prc, only: &
       PRC_abort    
    implicit none
    class(LineElement), intent(in) :: elem

    integer :: polyorder
    integer :: ierr

    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    character(len=H_MID) :: fname

    namelist / PARAM_ADVECT1D_NUMERROR / &
      polyOrderErrorCheck, &
      LOG_OUT_BASENAME,    &
      LOG_STEP_INTERVAL    
    !----------------------------------------------

    LOG_NEWLINE
    LOG_INFO("advect1d_numerror_init",*) 'Setup'

    polyOrderErrorCheck = 6
    LOG_STEP_INTERVAL   = 5

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ADVECT1D_NUMERROR,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO('advect1d_numerror_init',*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR('advect1d_numerror_init',*) 'Not appropriate names in namelist PARAM_ADVECT1D_NUMERROR. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ADVECT1D_NUMERROR)


    !--
    allocate( IntrpMat(polyOrderErrorCheck, elem%PolyOrder+1) )
    allocate( intw_intrp(polyOrderErrorCheck), x_intrp(polyOrderErrorCheck) )

    IntrpMat(:,:) = elem%GenIntGaussLegendreIntrpMat( &
      PolyOrderErrorCheck,  & ! (in)
      intw_intrp, x_intrp   ) ! (out)

    !--- Open logfile
    
    NUMERROR_LOG_FID = IO_get_available_fid()
    fname = trim(LOG_OUT_BASENAME)//".peall"
    open( unit   = NUMERROR_LOG_FID,                 &
          file   = fname,                            &
          form   = 'formatted',                      &
          iostat = ierr                              )
    if ( ierr /= 0 ) then
     LOG_ERROR('advect1d_numerror_init',*) 'File open error! :', trim(fname)
     call PRC_abort
    endif      

    write(NUMERROR_LOG_FID,'(A)',advance='no') '                   '
    write(NUMERROR_LOG_FID,'(A16)',advance='no') 'L1_error'
    write(NUMERROR_LOG_FID,'(A16)',advance='no') 'L2_error'
    write(NUMERROR_LOG_FID,'(A16)',advance='no') 'Linf_error'
    write(NUMERROR_LOG_FID,*)

    LOG_RSTEP_INTERVAL = LOG_STEP_INTERVAL

    return
  end subroutine advect1d_numerror_Init

  !> Inquire when we should evaluate numerical error
  function advect1d_numerror_do_step( step ) result(do_flag)
    integer, intent(in) :: step
    logical :: do_flag
    !------------------------------------------------------------------------

    LOG_RSTEP_INTERVAL = LOG_RSTEP_INTERVAL - 1
    if ( LOG_RSTEP_INTERVAL == 0 ) then
      do_flag = .true.
      LOG_RSTEP_INTERVAL = LOG_STEP_INTERVAL
    else
      do_flag = .false. 
    end if

    return
  end function advect1d_numerror_do_step

  !> Evaluate numerical errors
  subroutine advect1d_numerror_eval( qexact,            & ! (inout)
      q, tsec, ADV_VEL, InitShapeName, InitShapeParams, & ! (in)
      mesh, elem                                        ) ! (in)

    use scale_polynominal, only: Polynominal_genLegendrePoly
    use mod_fieldutil, only: &
      get_upwind_pos1d => fieldutil_get_upwind_pos1d,        &
      get_profile1d_tracer => fieldutil_get_profile1d_tracer 

    implicit none
    class(MeshField1D), intent(inout) :: qexact
    class(MeshField1D), intent(in) :: q
    real(RP), intent(in) :: tsec
    real(RP), intent(in) :: ADV_VEL
    character(*), intent(in) :: InitShapeName
    real(RP), intent(in) :: InitShapeParams(2)
    class(MeshLineDom1D), intent(in) :: mesh
    class(ElementBase1D), intent(in) :: elem

    class(LocalMeshBase), pointer :: lcmesh

    real(RP) :: q_intrp(polyOrderErrorCheck)
    real(RP) :: qexact_intrp(polyOrderErrorCheck)
    real(RP) :: x_uwind(elem%Np)
    real(RP) :: x_uwind_intrp(polyOrderErrorCheck)
    real(RP) :: pos_intrp(polyOrderErrorCheck)
    real(RP) :: vx(2)

    real(RP) :: l1error
    real(RP) :: l2error   
    real(RP) :: linferror
    real(RP) :: dom_xmin, dom_xmax

    integer :: idom, ke    
    !------------------------------------------------------------------------

    l1error   = 0.0_RP
    l2error   = 0.0_RP   
    linferror = 0.0_RP    
    dom_xmin = mesh%xmin_gl
    dom_xmax = mesh%xmax_gl

    do idom=1, mesh%LOCAL_MESH_NUM
      call mesh%GetLocalMesh( idom, lcmesh )
      do ke=lcmesh%NeS, lcmesh%NeE

        x_uwind(:) = get_upwind_pos1d(lcmesh%pos_en(:,ke,1), ADV_VEL, tsec, dom_xmin, dom_xmax)

        vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
        pos_intrp(:) = vx(1) + 0.5_RP*( x_intrp(:) + 1.0_RP ) * ( vx(2) - vx(1) )
        x_uwind_intrp(:) = get_upwind_pos1d(pos_intrp(:), ADV_VEL, tsec, dom_xmin, dom_xmax)

        call get_profile1d_tracer( qexact%local(idom)%val(:,ke),  & ! (out)
          InitShapeName, x_uwind, InitShapeParams, elem%Np        ) ! (in)

        call get_profile1d_tracer( qexact_intrp(:),                              & ! (out)
          InitShapeName, x_uwind_intrp, InitShapeParams(:), PolyOrderErrorCheck)   ! (in)
        
        q_intrp(:) = matmul(IntrpMat, q%local(idom)%val(:,ke))

        l1error = l1error &
          + sum( lcmesh%J(1,ke) * intw_intrp(:) * abs( q_intrp(:) - qexact_intrp(:) ) ) 

        l2error = l2error &
          + sum( lcmesh%J(1,ke) * intw_intrp(:) * ( q_intrp(:) - qexact_intrp(:) )**2 ) 
        
        linferror = max(linferror, maxval(abs(q%local(idom)%val(:,ke) - qexact%local(idom)%val(:,ke))))
      end do
    end do

    ! LOG_INFO("evaluate_error_l2",*) sqrt(l2error)/(dom_xmax - dom_xmin) 
    ! LOG_INFO("evaluate_error_linf",*) linferror

    write(NUMERROR_LOG_FID,'(A,ES15.8)',advance='no') 'tsec=', tsec
    write(NUMERROR_LOG_FID,'(A,ES15.8)',advance='no') ' ', l1error / (dom_xmax - dom_xmin) 
    write(NUMERROR_LOG_FID,'(A,ES15.8)',advance='no') ' ', sqrt(l2error) / (dom_xmax - dom_xmin) 
    write(NUMERROR_LOG_FID,'(A,ES15.8)',advance='no') ' ', linferror
    write(NUMERROR_LOG_FID,*)

    return
  end subroutine advect1d_numerror_eval

  !> Finalization
  subroutine advect1d_numerror_Final()
    implicit none
    !---------------------------------
    
    close( NUMERROR_LOG_FID )
    deallocate( IntrpMat, intw_intrp, x_intrp )

    return
  end subroutine advect1d_numerror_Final

end module mod_advect1d_numerror