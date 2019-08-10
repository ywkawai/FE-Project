!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_file_history_meshfield
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io

  use scale_file_history, only: &
    FILE_HISTORY_Setup, &
    FILE_HISTORY_Set_NowDate, &
    FILE_HISTORY_truncate_1D, &
    FILE_HISTORY_truncate_2D, &
    FILE_HISTORY_truncate_3D, &
    FILE_HISTORY_put,         &
    FILE_HISTORY_write,       &
    FILE_HISTORY_Set_Dim,     &
    FILE_HISTORY_Set_Axis,    &
    FILE_HISTORY_finalize

  use scale_element_base, only: elementbase1D, elementbase2D
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_localmesh_2d, only: LocalMesh2D

  use scale_meshfield_base, only: MeshField1D, MeshField2D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: FILE_HISTORY_meshfield_setup
  public :: FILE_HISTORY_meshfield_put
  public :: FILE_HISTORY_meshfield_write  
  public :: FILE_HISTORY_meshfield_finalize

  interface FILE_HISTORY_meshfield_put
    module procedure FILE_HISTORY_meshfield_put1D
    module procedure FILE_HISTORY_meshfield_put2D
  end interface
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  integer,          parameter :: nzs = 1
  character(len=8), parameter :: zs(nzs) = (/  "model   " /)

  integer  :: FILE_HISTORY_MESHFIELD_STARTDATE(6) !< start time [YYYY MM DD HH MM SS]
  real(DP) :: FILE_HISTORY_MESHFIELD_STARTMS      !< subsecond part of start time [millisec]

  class(MeshBase1D), pointer :: mesh1D
  
  class(MeshRectDom2D), pointer :: mesh2D
  integer :: mesh2D_icount, mesh2D_jcount

contains

!----------------

subroutine FILE_HISTORY_meshfield_setup( mesh1D_, mesh2D_ )

  use scale_file_h, only: &
    FILE_HSHORT
  use scale_prc, only: &
    PRC_masterrank, &
    PRC_myrank,     &
    PRC_abort
  use scale_time, only: &
    TIME_NOWDATE,       &
    TIME_NOWMS,         &
    TIME_STARTDAYSEC,   &
    TIME_DTSEC,         &
    TIME_NOWSTEP
  use scale_calendar, only: &
    CALENDAR_get_name
  implicit none

  class(Meshbase1d), intent(in), target, optional :: mesh1D_
  class(MeshRectDom2d), intent(in), target, optional :: mesh2D_

  character(len=H_MID) :: FILE_HISTORY_MESHFILED_H_TITLE = 'SCALE-FEM FILE_HISTORY_MESHFIELD' !< title of the output file
  character(len=H_MID) :: FILE_HISTORY_MESHFIELD_T_SINCE

  
  character(len=FILE_HSHORT) :: calendar
  real(DP) :: start_daysec
  integer  :: ierr
  integer  :: k

  !---------------------------------------------------------------------------

  FILE_HISTORY_MESHFIELD_STARTDATE(:) = TIME_NOWDATE
  FILE_HISTORY_MESHFIELD_STARTMS      = TIME_NOWMS

  start_daysec = TIME_STARTDAYSEC
  if ( TIME_NOWDATE(1) > 0 ) then
     write(FILE_HISTORY_MESHFIELD_T_SINCE,'(I4.4,5(A1,I2.2))') TIME_NOWDATE(1), &
                                                        '-', TIME_NOWDATE(2), &
                                                        '-', TIME_NOWDATE(3), &
                                                        ' ', TIME_NOWDATE(4), &
                                                        ':', TIME_NOWDATE(5), &
                                                        ':', TIME_NOWDATE(6)
     start_daysec = TIME_NOWMS
  else
     FILE_HISTORY_MESHFIELD_T_SINCE = ''
  endif

  ! get calendar name
  call CALENDAR_get_name( calendar )

  call FILE_HISTORY_Setup( FILE_HISTORY_MESHFILED_H_TITLE,        & ! [IN]
    H_SOURCE, H_INSTITUTE,                                        & ! [IN]
    start_daysec, TIME_DTSEC,                                     & ! [IN]
    time_since = FILE_HISTORY_MESHFIELD_T_SINCE,                  & ! [IN]
    calendar = calendar,                                          & ! [IN]
    default_zcoord = 'model',                                     & ! [IN]
    myrank = PRC_myrank                        )                    ! [IN]

  call FILE_HISTORY_Set_NowDate( TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP )

  if ( present(mesh1D_) ) then
    mesh1D => mesh1D_
    call set_dims1D()
    call set_axis1D()
  else if ( present(mesh2D_) ) then
    mesh2D => mesh2D_
    call set_dims2D()
    call set_axis2D()
  else
    LOG_ERROR("FILE_HISTORY_meshfield_setup",*)   "Any mesh (mesh1d, mesh2d) are not specified."
    call PRC_abort
  end if

  return
end subroutine FILE_HISTORY_meshfield_setup

subroutine FILE_HISTORY_meshfield_write
  call FILE_HISTORY_write
end subroutine FILE_HISTORY_meshfield_write

subroutine FILE_HISTORY_meshfield_finalize
  call FILE_HISTORY_finalize
end subroutine FILE_HISTORY_meshfield_finalize

subroutine FILE_HISTORY_meshfield_put1D(hstid, field1d)
  implicit none
  integer, intent(in) :: hstid
  class(MeshField1D), intent(in) :: field1d

  integer :: n, k, p
  class(Meshbase1d), pointer :: mesh
  type(LocalMesh1D), pointer :: lcmesh
  integer :: bufsize, ptr
  real(RP), allocatable :: buf(:)

  !--------------------
    
  bufsize = 0
  mesh => field1d%mesh

  do n=1, mesh%LOCAL_MESH_NUM
    lcmesh => mesh%lcmesh_list(n)
    bufsize = bufsize + lcmesh%refElem%Np * lcmesh%Ne 
  end do

  allocate( buf(bufsize) )
  ptr = 0
  do n=1, mesh%LOCAL_MESH_NUM
    lcmesh => mesh%lcmesh_list(n)
    do k=lcmesh%NeS, lcmesh%NeE
    do p=1, lcmesh%refElem%Np
      ptr = ptr + 1
      buf(ptr) =  field1d%local(n)%val(p,k)
    end do
    end do
  end do
  call FILE_HISTORY_put(hstid, buf)

  return
end subroutine FILE_HISTORY_meshfield_put1D

subroutine FILE_HISTORY_meshfield_put2D(hstid, field2d)
  implicit none
  integer, intent(in) :: hstid
  class(MeshField2D), intent(in) :: field2d

  integer :: n, k1, p
  integer :: i0, j0, i1, j1, i2, j2, i, j
  class(Meshbase2d), pointer :: mesh
  type(LocalMesh2D), pointer :: lcmesh
  type(elementbase2D), pointer :: refElem
  integer :: i0_s, j0_s

  real(RP), allocatable :: buf(:,:)

  !--------------------

  mesh => field2d%mesh
  allocate( buf(mesh2D_icount, mesh2D_jcount) )

  i0_s = 0; j0_s = 0
  
  do j0=1, size(mesh2D%rcdomIJ2LCMeshID,2)
  do i0=1, size(mesh2D%rcdomIJ2LCMeshID,1)
    n =  mesh2D%rcdomIJ2LCMeshID(i0,j0)

    lcmesh => mesh%lcmesh_list(n)
    refElem => lcmesh%refElem2D
      
    do j1=1, lcmesh%NeY
    do i1=1, lcmesh%NeX
       k1 = i1 + (j1-1)*lcmesh%NeX
       do j2=1, refElem%Nfp
       do i2=1, refElem%Nfp
         i = i0_s + i2 + (i1-1)*refElem%Nfp 
         j = j0_s + j2 + (j1-1)*refElem%Nfp
         buf(i,j) = field2d%local(n)%val(i2+(j2-1)*refElem%Nfp,k1)
       end do
       end do
    end do
    end do

    i0_s = i0_s + lcmesh%NeX * refElem%Nfp
    j0_s = j0_s + lcmesh%NeY * refElem%Nfp
  end do
  end do

  call FILE_HISTORY_put(hstid, buf)

  return
end subroutine FILE_HISTORY_meshfield_put2D

!----------------

subroutine set_dims1D()
  implicit none

  character(len=H_SHORT) :: dims(3,3)
  integer :: start(3,3), count(3,3)  

  start(1,1) = 1
  dims(1,1)  = "x"
  count(1,1) = mesh1D%NeG * mesh1D%refElem1D%Np
  call FILE_HISTORY_Set_Dim( "X", 1, 1, dims(:,:), zs(:), start(:,:), count(:,:))

  return
end subroutine set_dims1D

subroutine set_axis1D()
  implicit none

  real(DP) :: x(mesh1D%refElem1D%Np * mesh1D%NeG)
  integer :: n
  integer :: i
  integer :: i2
  type(ElementBase1D), pointer :: refElem
  type(LocalMesh1D), pointer :: lcmesh

  !--------------------
  
  do n=1,mesh1D%LOCAL_MESH_NUM
    lcmesh => mesh1D%lcmesh_list(n)
    refElem => lcmesh%refElem1D

    do i=1,mesh1D%NeG
    do i2=1, refElem%Np
      x(i2 + (i-1)*refElem%Np + (n-1)*refElem%Np*lcmesh%Ne) = mesh1D%lcmesh_list(n)%pos_en(i2,i,1)
    end do
    end do
  end do
  call FILE_HISTORY_Set_Axis( 'x', 'X-coordinate', '1', 'x', x)
  
  return
end subroutine set_axis1D

!----------------

subroutine set_dims2D()
  implicit none

  character(len=H_SHORT) :: dims(3,3)
  integer :: start(3,3), count(3,3)  

  type(ElementBase2D), pointer :: refElem
  type(LocalMesh2D), pointer :: lcmesh
  integer :: i, j, n
  integer :: icount, jcount

  !-------------------------------------
  
  icount = 0
  do i=1, size(mesh2D%rcdomIJ2LCMeshID,1)
    n = mesh2D%rcdomIJ2LCMeshID(i,1)
    lcmesh => mesh2D%lcmesh_list(n)
    icount = icount + lcmesh%NeX * lcmesh%refElem2D%Nfp
  end do

  jcount = 0
  do j=1, size(mesh2D%rcdomIJ2LCMeshID,2)
    n = mesh2D%rcdomIJ2LCMeshID(1,j)
    lcmesh => mesh2D%lcmesh_list(n)
    jcount = jcount + lcmesh%NeY * lcmesh%refElem2D%Nfp
  end do
  mesh2D_icount = icount
  mesh2D_jcount = jcount

  !--

  start(1,1) = 1
  dims(1,1)  = "x"
  count(1,1) = icount
  call FILE_HISTORY_Set_Dim( "X", 1, 1, dims(:,:), zs(:), start(:,:), count(:,:))

  start(1,1) = 1
  dims(1,1)  = "y"
  count(1,1) = jcount
  call FILE_HISTORY_Set_Dim( "Y", 1, 1, dims(:,:), zs(:), start(:,:), count(:,:))

  start(1,1) = 1
  start(2,1) = 1
  dims(1,1)  = "x"
  dims(2,1)  = "y"
  count(1,1) = icount
  count(2,1) = jcount
  call FILE_HISTORY_Set_Dim( "XY", 2, 1, dims(:,:), zs(:), start(:,:), count(:,:))

  return
end subroutine set_dims2D

subroutine set_axis2D()
  implicit none

  real(RP), allocatable :: x(:)
  real(RP), allocatable :: y(:)
  integer :: n
  integer :: k
  integer :: i, j 
  integer :: i2, j2
  type(ElementBase2D), pointer :: refElem
  type(LocalMesh2D), pointer :: lcmesh

  integer :: is, js, ie, je, igs, jgs

  !--------------------
  
  allocate( x(mesh2d_icount), y(mesh2d_jcount) )
  
  !--

  igs = 0; jgs = 0
  do n=1 ,mesh2D%LOCAL_MESH_NUM
    lcmesh => mesh2D%lcmesh_list(n)
    refElem => lcmesh%refElem2D

    do j=1, lcmesh%NeY
    do i=1, lcmesh%NeX
      k = i + (j-1) * lcmesh%NeX
      if ( j==1 ) then
        is = igs + 1 + (i-1)*refElem%Nfp
        ie = is + refElem%Nfp - 1
        x(is:ie) = lcmesh%pos_en(refElem%Fmask(:,1),k,1)
      end if
      if ( i==1 ) then
        js = jgs + 1 + (j-1)*refElem%Nfp
        je = js + refElem%Nfp - 1
        y(js:je) = lcmesh%pos_en(refElem%Fmask(:,4),k,2)
        write(*,*) js, je, ":", refElem%Fmask(:,4)
      end if
    end do
    end do

    igs = ie; jgs = je
  end do

  call FILE_HISTORY_Set_Axis( 'x', 'X-coordinate', '1', 'x', x)
  call FILE_HISTORY_Set_Axis( 'y', 'Y-coordinate', '1', 'y', y)

end subroutine set_axis2D

!----------------

end module scale_file_history_meshfield
