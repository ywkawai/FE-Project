#include "scalelib.h"
program test_field2d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale_file_history
  use scale

  use scale_element_quadrilateral
  use scale_localmesh_2d
  use scale_mesh_rectdom2d

  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  use scale_time_manager, only: &
    TIME_manager_Init , TIME_manager_Final,            &
    TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP,            &
    TIME_DTSEC, TIME_NSTEP 
    
  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_setup,   &
    FILE_HISTORY_meshfield_put,     &
    FILE_HISTORY_meshfield_write,   &
    FILE_HISTORY_meshfield_finalize

  implicit none

  integer, parameter :: NeGX = 3
  integer, parameter :: NeGY = 3
  integer, parameter :: NLocalMeshPerPrc = 1

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: dom_ymin = -1.0_RP
  real(RP), parameter :: dom_ymax = +1.0_RP

  type(QuadrilateralElement) :: refElem
  integer, parameter :: PolyOrder = 2
  
  type(MeshRectDom2D), target :: mesh

  type(MeshField2D), target :: q
  character(len=H_short), parameter :: q_VARNAME = "q"
  character(len=H_short), parameter :: q_DESC    = "q"
  character(len=H_short), parameter :: q_UNITS   = "K"
  integer :: HST_ID(1)

  type(MeshFieldCommRectDom2D) :: fields_comm
  type(MeshFieldContainer) :: field_list(1)  

  integer :: n
  type(LocalMesh2D), pointer :: lcmesh

  !---------------
  
  call init()

  !----
  write(*,*) "* Check data communication.."
  call perform_comm()

  do n=1, mesh%LOCAL_MESH_NUM
    write(*,*) "Check interior & halo data.."
    write(*,*) "tileID=", lcmesh%tileID
    call check_interior_data(n, mesh%lcmesh_list(n), q%local(n)%val, q%local(n)%val)
    call check_halo_data(n, mesh%lcmesh_list(n), q%local(n)%val, q%local(n)%val)
  end do
  !----

  call final()

contains

  elemental function get_field_val(tileID, k, p) result(val)
    integer, intent(in) :: tileID, k, p
    real(RP) :: val
    !----------------------------
    val = tileID*1000000 + k*1000 + p
  end function get_field_val

  subroutine init()
    use scale_calendar, only: CALENDAR_setup

    integer :: comm, myrank, nprocs
    logical :: ismaster
    
    integer :: k, p
    !---------

    call PRC_MPIstart( comm )
    
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    call IO_setup( "test", allow_noconf = .true. )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
    
    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init

    !------
    call refElem%Init(PolyOrder, .true.)

    call mesh%Init( &
      NeGX, NeGY,                             &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
      .true., .true.,                         &
      refElem, NLocalMeshPerPrc )
    
    call mesh%Generate()

    !---
    call q%Init( q_VARNAME, q_UNITS, mesh )
    call fields_comm%Init(1, 0, mesh)
    
    call FILE_HISTORY_meshfield_setup( mesh2d_=mesh )
    call FILE_HISTORY_reg( q_VARNAME, q_DESC, q_UNITS, HST_ID(1), dim_type='XY')

    !---

    !---
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE
      do p=1, refElem%Np
        q%local(n)%val(p,k) = get_field_val(lcmesh%tileID, k, p)
      end do
      end do
    end do

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_write()    

    !---


  end subroutine init

  subroutine final()

    call FILE_HISTORY_meshfield_finalize()

    call q%Final()
    call fields_comm%Final()
    call mesh%Final()
    call refElem%Final()
    
    call TIME_manager_Final()
    call PRC_MPIfinish()

  end subroutine final

  subroutine perform_comm()
    type(MeshFieldContainer) :: field_list(1)
    !---------------------------

    field_list(1)%field2d => q
    write(*,*) " - Put.."
    call fields_comm%Put(field_list, 1)
    write(*,*) " - Exchange.."
    call fields_comm%Exchange()
    write(*,*) " - Get.."
    call fields_comm%Get(field_list, 1)

  end subroutine perform_comm
    
  subroutine check_interior_data(n, lcmesh_, lcfield, lcfield_fl)
    integer, intent(in) :: n
    type(LocalMesh2D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: k, f, p
    integer :: iM, iP, fnid
    integer :: ans(refElem%Np)
    !------------------------------

    do k=lcmesh_%NeS, lcmesh_%NeE
      do p=1, refElem%Np
        ans(p) = get_field_val(lcmesh_%tileID,k,p)
      end do

      call assert(k, int(lcfield(:,k)), ans(:), 'check_interior_data', 'q', refElem%Np)
      do f=1, refElem%Nfaces
      do p=1, refElem%Nfp
        fnid = p+(f-1)*refElem%Nfp
        iM = lcmesh_%VMapM(fnid,k)
        iP = lcmesh_%VMapP(fnid,k)
        write(*,'(a,3(i3),a,i10,a,i10,a,3i4)') "k , f, p =", k,f,p,           &
          " : q- = ", int(lcfield_fl(iM)), ",  q+ =", int(lcfield_fl(iP)),  &
          ", normal=", int(lcmesh_%normal_fn(fnid,k,:))
      end do
      end do

    end do

  end subroutine check_interior_data

  subroutine check_halo_data(n, lcmesh_, lcfield, lcfield_fl)
    implicit none

    integer, intent(in) :: n
    type(LocalMesh2D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ans_bx(refElem%Nfp,lcmesh_%NeX)
    integer :: ans_by(refElem%Nfp,lcmesh_%NeY)    
    integer :: haloInd_s, haloInd_e
    integer :: i, j, k, f, fp
    integer :: tileID

    !------------------------------

    tileID = lcmesh_%tileID

    write(*,*) "Check values in halo.."
    write(*,*) "tileID_globalMap=", mesh%tileID_globalMap(:,tileID)
    write(*,*) "tileFaceID_globalMap=", mesh%tileFaceID_globalMap(:,tileID)


    haloInd_s = refElem%Np*lcmesh_%Ne + 1    
    haloInd_e = haloInd_s + refElem%Nfp*lcmesh_%NeX - 1
    do i=1, lcmesh_%NeX
      k = i + (lcmesh_%NeY-1)*lcmesh_%NeX
      do fp=1, refElem%Nfp
        ans_bx(fp,i) = get_field_val(mesh%tileID_globalMap(1,tileID), k, refElem%Fmask(fp,3))
      end do
    end do
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans_bx(:,:),     &
                 'check_halo_data', 'q_halo_south', refElem%Nfp*lcmesh_%NeX)

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + refElem%Nfp*lcmesh_%NeY - 1
    do j=1, lcmesh_%NeY
      k = 1 + (j-1)*lcmesh_%NeX
      do fp=1, refElem%Nfp
        ans_by(fp,j) = get_field_val(mesh%tileID_globalMap(2,tileID), k, refElem%Fmask(fp,4))
      end do
    end do
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans_by(:,:),     &
      'check_halo_data', 'q_halo_west', refElem%Nfp*lcmesh_%NeY)

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + refElem%Nfp*lcmesh_%NeX - 1    
    do i=1, lcmesh_%NeX
      k = i
      do fp=1, refElem%Nfp
        ans_bx(fp,i) = get_field_val(mesh%tileID_globalMap(3,tileID), k, refElem%Fmask(fp,1))
      end do
    end do
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans_bx(:,:),     &
                 'check_halo_data', 'q_halo_north', refElem%Nfp*lcmesh_%NeX)

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + refElem%Nfp*lcmesh_%NeY - 1    
    do j=1, lcmesh_%NeY
      k = lcmesh_%NeX + (j-1)*lcmesh_%NeX
      do fp=1, refElem%Nfp
        ans_by(fp,j) = get_field_val(mesh%tileID_globalMap(4,tileID), k, refElem%Fmask(fp,2))
      end do
    end do    
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans_by(:,:),     &
      'check_halo_data', 'q_halo_east', refElem%Nfp*lcmesh_%NeY)    

  end subroutine check_halo_data

  subroutine assert(k, vals, ans, assert_name, var_name, val_size)
    integer, intent(in) :: val_size
    integer, intent(in) :: k
    integer, intent(in) :: vals(val_size)
    integer, intent(in) :: ans(val_size)
    character(*), intent(in) :: assert_name
    character(*), intent(in) :: var_name

    real(RP), parameter :: EPS = 1.0E-15_RP
    !--------------------------------------

    write(*,*) trim(var_name), "=", vals(:)
    if ( sum((vals(:) - ans(:))**2) > EPS ) then
      LOG_ERROR(assert_name,*) 'The value of '//trim(var_name)//' is unexcepted!', &
        ' k=', k, ": val=", vals(:), " ans=", ans(:)
      call PRC_abort
    end if    
  end subroutine assert

end program test_field2d
