#include "scalelib.h"
program test_field1d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale_file_history
  use scale

  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d

  use scale_localmeshfield_base, only: LocalMeshField1D
  use scale_meshfield_base, only: MeshField1D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_1d, only: MeshFieldComm1D

  use scale_time_manager, only: &
    TIME_manager_Init , TIME_manager_Final,            &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,        &
    TIME_DTSEC, TIME_NSTEP 
    
  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_setup,   &
    FILE_HISTORY_meshfield_put,     &
    FILE_HISTORY_meshfield_write,   &
    FILE_HISTORY_meshfield_finalize

  implicit none

  integer, parameter :: NeGX = 4
  integer, parameter :: NLocalMeshPerPrc = 1

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP

  type(LineElement) :: refElem
  integer, parameter :: PolyOrder = 2
  
  type(MeshLineDom1D), target :: mesh

  type(MeshField1D), target :: q
  character(len=H_short), parameter :: q_VARNAME = "q"
  character(len=H_short), parameter :: q_DESC    = "q"
  character(len=H_short), parameter :: q_UNITS   = "K"
  integer, save :: HST_ID(1)

  type(MeshFieldComm1D) :: fields_comm
  type(MeshFieldContainer), save :: field_list(1)  

  type(LocalMesh1D), pointer :: lcmesh
  integer :: n
  
  !------------------------------------------------------------------

  call init()

  !----
  write(*,*) "* Check data communication.."
  call perform_comm()

  do n=1, mesh%LOCAL_MESH_NUM
    lcmesh => mesh%lcmesh_list(n)

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
    integer :: n, k, p

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
      NeGX,                              &
      dom_xmin, dom_xmax,                &
      refElem, NLocalMeshPerPrc )
    
    call mesh%Generate()

    !---
    call q%Init( q_VARNAME, q_UNITS, mesh )
    call fields_comm%Init(1, 0, mesh)
    
    call FILE_HISTORY_meshfield_setup( mesh )
    call FILE_HISTORY_reg( q_VARNAME, q_DESC, q_UNITS, HST_ID(1), dim_type='X')

    !---
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE
      do p=1, refElem%Np
        q%local(n)%val(p,k) = get_field_val(lcmesh%tileID,k,p)
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

    field_list(1)%field1d => q
    write(*,*) " - Put.."
    call fields_comm%Put(field_list, 1)
    write(*,*) " - Exchange.."
    call fields_comm%Exchange()
    write(*,*) " - Get.."
    call fields_comm%Get(field_list, 1)

  end subroutine perform_comm

  subroutine check_interior_data(n, lcmesh_, lcfield, lcfield_fl)
    integer, intent(in) :: n
    type(LocalMesh1D), intent(in) :: lcmesh_
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
        write(*,'(a,3(i3),a,i10,a,i10,a,i10)') "k , f, p =", k,f,p,           &
          " : q- = ", int(lcfield_fl(iM)), ",  q+ =", int(lcfield_fl(iP)),  &
          ", nx=", int(lcmesh_%normal_fn(fnid,k,1))
      end do
      end do

    end do

  end subroutine check_interior_data

  subroutine check_halo_data(n, lcmesh_, lcfield, lcfield_fl)
    implicit none

    integer, intent(in) :: n
    type(LocalMesh1D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ans(2)
    integer :: haloInd_s, haloInd_e

    !------------------------------

    write(*,*) "Check values in halo.."
    write(*,*) "tileID_globalMap=", mesh%tileID_globalMap(:,n)
    write(*,*) "tileFaceID_globalMap=", mesh%tileFaceID_globalMap(:,n)

    haloInd_s = refElem%Np*lcmesh_%Ne + 1
    haloInd_e = haloInd_s + 2 - 1

    ans(1) = get_field_val(mesh%tileID_globalMap(1,n), lcmesh_%Ne, refElem%Np) 
    ans(2) = get_field_val(mesh%tileID_globalMap(2,n), 1, 1) 
    
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans(:), &
                 'check_halo_data', 'q_halo', 2)

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

end program test_field1d
