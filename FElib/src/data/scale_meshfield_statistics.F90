#include "scaleFElib.h"
module scale_meshfield_statistics

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi ! TODO: to replace functions in scale_comm module
  use scale_precision
  use scale_io
  use scale_prc, only: &
     PRC_LOCAL_COMM_WORLD
  use scale_prof

  use scale_element_base, only: ElementBase
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField1D, MeshField2D, MeshField3D
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  public :: MeshField_statistics_setup

  interface MeshField_statistics_detail
    module procedure :: MeshField_statistics_detail_1D
    module procedure :: MeshField_statistics_detail_2D
    module procedure :: MeshField_statistics_detail_3D
  end interface
  public :: MeshField_statistics_detail

  interface MeshField_statistics_total
    module procedure :: MeshField_statistics_total_1D
    module procedure :: MeshField_statistics_total_2D
    module procedure :: MeshField_statistics_total_3D
  end interface
  public :: MeshField_statistics_total

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !
  type :: MeshField_statistics_base
    logical :: use_globalcomm
    integer :: comm_datatype
  end type MeshField_statistics_base

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  type(MeshField_statistics_base), private :: base

contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MeshField_statistics_setup()
    use scale_prc, only: &
      PRC_abort
    implicit none

    logical :: use_globalcomm = .false.    
    namelist / PARAM_MESHFIELD_statistics / &
      use_globalcomm

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MeshField_statistics_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MESHFIELD_STATISTICS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("MeshField_statistics_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("MeshField_statistics_setup",*) 'Not appropriate names in namelist PARAM_MeshField_statistics. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_MESHFIELD_STATISTICS)

    base%use_globalcomm = use_globalcomm
    
    if ( base%use_globalcomm ) then
      LOG_INFO_CONT(*) '=> The total is calculated for the global domain.'
    else
      LOG_INFO_CONT(*) '=> The total is calculated for the local domain.'
    endif

    if ( RP == kind(0.D0) ) then
      base%comm_datatype = MPI_DOUBLE_PRECISION
   elseif( RP == kind(0.0) ) then
      base%comm_datatype = MPI_REAL
   else
      LOG_ERROR("MeshField_statistics_setup",*) 'precision is not supportd'
      call PRC_abort
   endif

    return
  end subroutine MeshField_statistics_setup

  !-----------------------------------------------------------------------------
  !> Calculate domain sum and area-weighted mean for 1D field
  subroutine MeshField_statistics_total_1D( field, & ! (in)
    log_suppress, global,                          & ! (in)
    mean, sum                                      ) ! (out)

    use scale_prc, only: &
        PRC_myrank, &
        PRC_abort
    use scale_const, only: &
        EPS   => CONST_EPS, &
        UNDEF => CONST_UNDEF
    implicit none

    class(MeshField1D), intent(in) :: field
    logical,  intent(in),  optional :: log_suppress !< suppress log output
    logical,  intent(in),  optional :: global       !< global or local sum
    real(RP), intent(out), optional :: mean !< area-weighted mean
    real(DP), intent(out), optional :: sum  !< domain sum

    real(DP) :: statval
    real(DP) :: total    !< total volume
    !---------------------------------------------------------------------------

    call calculate_statval( field, field%mesh%lcmesh_list(:), &  ! (in)
                            statval, total                    )  ! (out)

    call statistics_total_core(       &
      field%varname, statval, total,  & ! (in)
      log_suppress, global,           & ! (in)
      mean, sum     )                   ! (out)
    
    return
  end subroutine MeshField_statistics_total_1D

  !-----------------------------------------------------------------------------
  !> Calculate domain sum and area-weighted mean for 2D field
  subroutine MeshField_statistics_total_2D( field, & ! (in)
    log_suppress, global,                          & ! (in)
    mean, sum                                      ) ! (out)

    use scale_prc, only: &
        PRC_myrank, &
        PRC_abort
    use scale_const, only: &
        EPS   => CONST_EPS, &
        UNDEF => CONST_UNDEF
    implicit none

    class(MeshField2D), intent(in) :: field
    logical,  intent(in),  optional :: log_suppress !< suppress log output
    logical,  intent(in),  optional :: global       !< global or local sum
    real(RP), intent(out), optional :: mean !< area-weighted mean
    real(DP), intent(out), optional :: sum  !< domain sum

    real(DP) :: statval
    real(DP) :: total    !< total volume
    !---------------------------------------------------------------------------

    call calculate_statval( field, field%mesh%lcmesh_list(:), &  ! (in)
                            statval, total                    )  ! (out)

    call statistics_total_core(       &
      field%varname, statval, total,  & ! (in)
      log_suppress, global,           & ! (in)
      mean, sum     )                   ! (out)
    
    return
  end subroutine MeshField_statistics_total_2D

  !-----------------------------------------------------------------------------
  !> Calculate domain sum and area-weighted mean for 3D field
  subroutine MeshField_statistics_total_3D( field, & ! (in)
    log_suppress, global,                          & ! (in)
    mean, sum                                      ) ! (out)

    use scale_prc, only: &
        PRC_myrank, &
        PRC_abort
    use scale_const, only: &
        EPS   => CONST_EPS, &
        UNDEF => CONST_UNDEF
    implicit none

    class(MeshField3D), intent(in) :: field
    logical,  intent(in),  optional :: log_suppress !< suppress log output
    logical,  intent(in),  optional :: global       !< global or local sum
    real(RP), intent(out), optional :: mean !< area-weighted mean
    real(DP), intent(out), optional :: sum  !< domain sum

    real(DP) :: statval
    real(DP) :: total    !< total volume
    !---------------------------------------------------------------------------

    call calculate_statval( field, field%mesh%lcmesh_list(:), &  ! (in)
                            statval, total                    )  ! (out)

    call statistics_total_core(       &
      field%varname, statval, total,  & ! (in)
      log_suppress, global,           & ! (in)
      mean, sum     )                   ! (out)
    
    return
  end subroutine MeshField_statistics_total_3D

  !-----------------------------------------------------------------------------
  !> Search global maximum & minimum value for 1D field
  subroutine MeshField_statistics_detail_1D( field_list, local )
    implicit none

    class(MeshField1D), intent(in) :: field_list(:)
    logical,          intent(in), optional :: local  !< calculate in local node

    real(RP) :: statval_l(  size(field_list),2)
    integer  :: statidx_l(3,size(field_list),2)
    integer :: VA
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MeshField_STATISTICS_detail_1D",*) 'Variable Statistics '

    VA = size(field_list)
    call search_maxmin_local( &
      VA, field_list, field_list(1)%mesh%lcmesh_list, & ! (in)
      statval_l, statidx_l                            ) ! (out)

    call statistics_detail_core( &
      VA, field_list, statval_l, statidx_l, local ) ! (in)
    
    LOG_NEWLINE

    return
  end subroutine MeshField_statistics_detail_1D

  !-----------------------------------------------------------------------------
  !> Search global maximum & minimum value for 2D field
  subroutine MeshField_statistics_detail_2D( field_list, local )
    implicit none

    class(MeshField2D), intent(in) :: field_list(:)
    logical,          intent(in), optional :: local  !< calculate in local node

    real(RP) :: statval_l(  size(field_list),2)
    integer  :: statidx_l(3,size(field_list),2)
    integer :: VA
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MeshField_STATISTICS_detail_2D",*) 'Variable Statistics '

    VA = size(field_list)
    call search_maxmin_local( &
      VA, field_list, field_list(1)%mesh%lcmesh_list, & ! (in)
      statval_l, statidx_l                            ) ! (out)

    call statistics_detail_core( &
      VA, field_list, statval_l, statidx_l, local ) ! (in)
    
    LOG_NEWLINE

    return
  end subroutine MeshField_statistics_detail_2D  
  !-----------------------------------------------------------------------------
  !> Search global maximum & minimum value for 3D field
  subroutine MeshField_statistics_detail_3D( field_list, local )
    implicit none

    class(MeshField3D), intent(in) :: field_list(:)
    logical,          intent(in), optional :: local  !< calculate in local node

    real(RP) :: statval_l(  size(field_list),2)
    integer  :: statidx_l(3,size(field_list),2)
    integer :: VA
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MeshField_STATISTICS_detail_3D",*) 'Variable Statistics '

    VA = size(field_list)
    call search_maxmin_local( &
      VA, field_list, field_list(1)%mesh%lcmesh_list, & ! (in)
      statval_l, statidx_l                            ) ! (out)

    call statistics_detail_core( &
      VA, field_list, statval_l, statidx_l, local ) ! (in)
    
    LOG_NEWLINE

    return
  end subroutine MeshField_statistics_detail_3D

!-- private --------------------------------------------------------------------

  subroutine search_maxmin_local( VA, field_list, lcmesh_list, &
    statval_l, statidx_l )

    implicit none
    integer, intent(in) :: VA
    class(MeshFieldBase), intent(in) :: field_list(VA)
    class(LocalMeshBase), intent(in), target :: lcmesh_list(:) 
    real(RP), intent(out) :: statval_l (  VA,2)
    integer,  intent(out) :: statidx_l (3,VA,2)

    type(LocalMeshBase), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: lcfield
    class(ElementBase), pointer :: refElem

    integer :: v
    integer :: n
    integer :: ke, p
    !---------------------------------------------------------------------------

    do n=1, size(lcmesh_list)
      lcmesh => lcmesh_list(n)
      refElem => lcmesh%refElem
      do v=1, VA
        call field_list(v)%GetLocalMeshField(n, lcfield)
        if (n==1) then
          statval_l(  v,:) = lcfield%val(1,lcmesh%NeS)
          statidx_l(1,v,:) = 1
          statidx_l(2,v,:) = lcmesh%NeS
          statidx_l(3,v,:) = n
        end if
        do ke=lcmesh%NeS, lcmesh%NeE
        do p =1, refElem%Np
          if ( lcfield%val(p,ke) > statval_l(v,1) ) then
            statval_l(  v,1) = lcfield%val(p,ke)
            statidx_l(:,v,1) = (/ p, ke, n /)
          end if
          if ( lcfield%val(p,ke) < statval_l(v,2) ) then
            statval_l(  v,2) = lcfield%val(p,ke)
            statidx_l(:,v,2) = (/ p, ke, n /)
          end if
        end do
        end do
      end do
    end do

    return
  end subroutine search_maxmin_local

  subroutine statistics_detail_core(      &
    VA, field_list, statval_l, statidx_l, &
    local )

    use scale_prc, only: &
       PRC_nprocs
    implicit none

    integer, intent(in) :: VA
    class(MeshFieldBase), intent(in) :: field_list(VA)
    real(DP), intent(in) :: statval_l(  VA,2)
    integer,  intent(in) :: statidx_l(3,VA,2)
    logical,  intent(in), optional :: local  !< calculate in local node

    real(RP) :: statval   (  VA,2,0:PRC_nprocs-1)
    integer  :: statidx   (3,VA,2,0:PRC_nprocs-1)
    real(RP) :: allstatval(VA,2)
    integer  :: allstatidx(VA,2)
    logical :: do_globalcomm

    integer :: v
    integer :: p
    integer :: ierr
    !---------------------------------------------------------------------------

    do_globalcomm = base%use_globalcomm
    if ( present(local) ) do_globalcomm = ( .not. local )

    if ( do_globalcomm ) then
      call PROF_rapstart('COMM_Bcast', 2)

      call MPI_AllGather( statval_l(:,:),       &
                          VA*2,                 &
                          base%comm_datatype,   &
                          statval(:,:,:),       &
                          VA*2,                 &
                          base%comm_datatype,   &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr )

      call MPI_AllGather( statidx_l(:,:,:),     &
                          3*VA*2,               &
                          MPI_INTEGER,          &
                          statidx(:,:,:,:),     &
                          3*VA*2,               &
                          MPI_INTEGER,          &
                          PRC_LOCAL_COMM_WORLD, &
                          ierr )

      call PROF_rapend  ('COMM_Bcast', 2)

      do v = 1, VA
        allstatval(v,1) = statval(v,1,0)
        allstatval(v,2) = statval(v,2,0)
        allstatidx(v,:) = 0
        do p = 1, PRC_nprocs-1
            if ( statval(v,1,p) > allstatval(v,1) ) then
              allstatval(v,1) = statval(v,1,p)
              allstatidx(v,1) = p
            end if
            if ( statval(v,2,p) < allstatval(v,2) ) then
              allstatval(v,2) = statval(v,2,p)
              allstatidx(v,2) = p
            end if
        end do
        LOG_INFO_CONT(*) '[', trim(field_list(v)%varname), ']'
        LOG_INFO_CONT('(1x,A,ES17.10,A,4(I5,A))') '  MAX =', &
                                                      allstatval(v,1), ' (rank=', &
                                                      allstatidx(v,1), '; ', &
                                        statidx(1,v,1,allstatidx(v,1)),',', &
                                        statidx(2,v,1,allstatidx(v,1)),',', &
                                        statidx(3,v,1,allstatidx(v,1)),')'
        LOG_INFO_CONT('(1x,A,ES17.10,A,4(I5,A))') '  MIN =', &
                                                      allstatval(v,2), ' (rank=', &
                                                      allstatidx(v,2), '; ', &
                                        statidx(1,v,2,allstatidx(v,2)),',', &
                                        statidx(2,v,2,allstatidx(v,2)),',', &
                                        statidx(3,v,2,allstatidx(v,2)),')'
      enddo
    else
        ! statistics on each node
        do v = 1, VA
          LOG_INFO_CONT(*) '[', trim(field_list(v)%varname), ']'
          LOG_INFO_CONT('(1x,A,ES17.10,A,3(I5,A))') 'MAX = ', &
                                                statval_l(  v,1),' (', &
                                                statidx_l(1,v,1),',', &
                                                statidx_l(2,v,1),',', &
                                                statidx_l(3,v,1),')'
          LOG_INFO_CONT('(1x,A,ES17.10,A,3(I5,A))') 'MIN = ', &
                                                statval_l(  v,2),' (', &
                                                statidx_l(1,v,2),',', &
                                                statidx_l(2,v,2),',', &
                                                statidx_l(3,v,2),')'
        enddo
    endif

    return
  end subroutine statistics_detail_core

  !-----

  subroutine calculate_statval( field, lcmesh_list, &                   
    statval, total )
    
    implicit none

    class(MeshFieldBase), intent(in) :: field
    class(LocalMeshBase), intent(in), target :: lcmesh_list(:) 
    real(DP), intent(out) :: statval
    real(DP), intent(out) :: total    !< total volume

    type(LocalMeshBase), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: lcfield
    class(ElementBase), pointer :: refElem

    real(DP) :: statval_lc
    real(DP) :: total_lc
    real(DP) :: weight

    integer :: n
    integer :: ke, p
    !---------------------------------------------------------------------------

    statval = 0.0_DP
    total   = 0.0_DP
    do n=1, size(lcmesh_list)
      lcmesh => lcmesh_list(n)
      refElem => lcmesh%refElem
      call field%GetLocalMeshField(n, lcfield)

      total_lc   = 0.0_DP
      statval_lc = 0.0_DP
      !$omp parallel do private(p, weight) reduction(+:total_lc, statval_lc)
      do ke=lcmesh%NeS, lcmesh%NeE
      do p=1, refElem%Np
        weight = lcmesh%J(p,ke) * refElem%IntWeight_lgl(p)
        total_lc = total_lc + weight
        statval_lc = statval_lc + weight * lcfield%val(p,ke)
      end do
      end do
      total = total + total_lc
      statval = statval + statval_lc
    end do

    return
  end subroutine calculate_statval
  
  subroutine statistics_total_core( &
    varname, statval, total,        &
    log_suppress, &
    global,       &
    mean, sum     )
    use scale_prc, only: &
        PRC_myrank, &
        PRC_abort
    use scale_const, only: &
        EPS   => CONST_EPS, &
        UNDEF => CONST_UNDEF

    implicit none

    character(len=*), intent(in) :: varname
    real(DP), intent(in) :: statval
    real(DP), intent(in) :: total    !< total volume
    logical,  intent(in),  optional :: log_suppress !< suppress log output
    logical,  intent(in),  optional :: global       !< global or local sum
    real(RP), intent(out), optional :: mean !< area-weighted mean
    real(DP), intent(out), optional :: sum  !< domain sum

    real(DP) :: sendbuf(2), recvbuf(2)
    real(DP) :: sum_, mean_

    logical :: suppress_, global_
    integer :: ierr
    !---------------------------------------------------------------------------

    if ( .NOT. ( statval > -1.0_DP .OR. statval < 1.0_DP ) ) then ! must be NaN
        LOG_ERROR("MeshField_STATISTICS_total",*) 'NaN is detected for ', trim(varname), ' in rank ', PRC_myrank
        call PRC_abort
    endif

    if ( present(log_suppress) ) then
        suppress_ = log_suppress
    else
        suppress_ = .false.
    end if

    if ( present(global) ) then
        global_ = global
    else
        global_ = base%use_globalcomm
    end if

    if ( global_ ) then
        call PROF_rapstart('COMM_Allreduce', 2)
        sendbuf(1) = statval
        sendbuf(2) = total
        ! All reduce
        call MPI_Allreduce( sendbuf(:), recvbuf(:), &
                            2,                      &
                            MPI_DOUBLE_PRECISION,   &
                            MPI_SUM,                &
                            PRC_LOCAL_COMM_WORLD,   &
                            ierr                    )
        call PROF_rapend  ('COMM_Allreduce', 2)

        if ( recvbuf(2) < EPS ) then
          sum_  = UNDEF
          mean_ = UNDEF
        else
          sum_  = recvbuf(1)
          mean_ = recvbuf(1) / recvbuf(2)
        end if
        ! statistics over the all node
        if ( .not. suppress_ ) then ! if varname is empty, suppress output
          LOG_INFO("MeshField_STATISTICS_total_3D",'(1x,A,A24,A,ES24.17)') &
                      '[', trim(varname), '] MEAN(global) = ', mean_
        endif
    else
        if ( total < EPS ) then
          sum_  = UNDEF
          mean_ = UNDEF
        else
          sum_  = statval
          mean_ = statval / total
        end if

        ! statistics on each node
        if ( .not. suppress_ ) then ! if varname is empty, suppress output
          LOG_INFO("MeshField_STATISTICS_total_3D",'(1x,A,A24,A,ES24.17)') &
                      '[', trim(varname), '] MEAN(local)  = ', mean_
        endif
    endif

    if ( present(mean) ) mean = mean_
    if ( present(sum ) ) sum  = sum_

    return
  end subroutine statistics_total_core

end module scale_meshfield_statistics
