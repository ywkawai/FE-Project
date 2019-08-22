program test_field2d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale_file_history
  use scale

  use scale_element_quadrial
  use scale_localmesh_2d
  use scale_mesh_rectdom2d

  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

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

  type(QuadrialElement) :: refElem
  integer, parameter :: PolyOrder = 2
  
  type(MeshRectDom2D), target :: mesh

  type(MeshField2D), target :: q
  character(len=H_short), parameter :: q_VARNAME = "q"
  character(len=H_short), parameter :: q_DESC    = "q"
  character(len=H_short), parameter :: q_UNITS   = "K"
  integer :: HST_ID(1)

  type(MeshFieldCommRectDom2D) :: fields_comm
  type(MeshFieldContainer) :: field_list(1)  

  integer :: n, k, p
  type(LocalMesh2D), pointer :: lcmesh

  !---------------
  
  call init()
  !----
  call check_field_info()
  call check_comm()
  !----
  call final()

contains
  subroutine init()
    use scale_calendar, only: &
      CALENDAR_setup,        &
      CALENDAR_date2daysec,  &
      CALENDAR_combine_daysec
    use scale_time, only: &
      TIME_STARTDAYSEC, &
      TIME_NOWDATE,     &
      TIME_NOWMS,       &
      TIME_NOWSTEP,     &
      TIME_DTSEC,       &
      TIME_OFFSET_YEAR, &
      TIME_NOWSTEP

    integer :: TIME_STARTDATE(6) = (/ 0, 1, 1, 0, 0, 0 /)
    real(DP) :: TIME_STARTMS      = 0.0_DP !< [millisec]
    integer :: TIME_STARTDAY
    real(DP) :: TIME_STARTSEC
    
    integer :: comm, myrank, nprocs
    logical :: ismaster

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

    TIME_OFFSET_YEAR = TIME_STARTDATE(1)
    
    call CALENDAR_date2daysec( TIME_STARTDAY,     & ! [OUT]
                               TIME_STARTSEC,     & ! [OUT]
                               TIME_STARTDATE(:), & ! [IN]
                               TIME_STARTMS,      & ! [IN]
                               TIME_OFFSET_YEAR   ) ! [IN]

    TIME_STARTDAYSEC  = CALENDAR_combine_daysec( TIME_STARTDAY, TIME_STARTSEC )
    TIME_NOWDATE(:) = TIME_STARTDATE(:)
    TIME_NOWMS      = TIME_STARTMS
    TIME_DTSEC      = 1.0_DP
    TIME_NOWSTEP    = 1

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
        q%local(n)%val(p,k) = n*1000000 + k*1000 + p
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
    
    call PRC_MPIfinish()

  end subroutine final

  subroutine check_field_info()
    
    write(*,*) "Check values of q.."
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE
        write(*,*) "n,k=", n, k, ":", int(q%local(n)%val(:,k))
      end do
    end do

  end subroutine check_field_info

  subroutine check_comm()
    type(MeshFieldContainer) :: field_list(1)
    integer :: f
    !---------------------------

    write(*,*) "* Check data communication.."

    field_list(1)%field2d => q
    write(*,*) " - Put.."
    call fields_comm%Put(field_list, 1)
    write(*,*) " - Exchange.."
    call fields_comm%Exchange()
    write(*,*) " - Get.."
    call fields_comm%Get(field_list, 1)

    write(*,*) " - Check halo.."
    do n=1, mesh%LOCAL_MESH_NUM
      write(*,*) "n=", n
      lcmesh => mesh%lcmesh_list(n)
      call check_halo_data( q%local(n)%val, lcmesh)
    end do

  end subroutine check_comm

  subroutine check_halo_data( q_lc, lcmesh_ )
    type(LocalMesh2D), intent(in) :: lcmesh_
    real(RP), intent(in) :: q_lc(lcmesh%refElem2D%Np*lcmesh%NeA)

    integer :: f, iM, iP, fnid
    real(RP) :: qM, qP

    do k=lcmesh%NeS, lcmesh%NeE
      do f=1, refElem%Nfaces
      do p=1, refElem%Nfp
        fnid = p + (f-1) * refElem%Nfp
        iM = lcmesh_%VMapM(fnid,k)
        iP = lcmesh_%VMapP(fnid,k)
        write(*,*) iM, iP
        write(*,'(a,3(i3),a,i10,a,i10,2(a,i10))') "k , f, p =", k,f,p, &
          " : q- = ", int(q_lc(iM)), ",  q+ =", int(q_lc(iP)),        &
          ", nx=", int(lcmesh%normal_fn(fnid,k,1)),                   &           
          ", nx=", int(lcmesh%normal_fn(fnid,k,2))

      end do
      end do
      write(*,*) "-"
    end do
  end subroutine check_halo_data

end program test_field2d
