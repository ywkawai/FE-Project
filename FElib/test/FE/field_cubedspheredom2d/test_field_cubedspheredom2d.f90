#include "scalelib.h"
program test_field_cubedspheredom2d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale_file_history
  use scale
  use scale_const, only: &
    RPlanet => CONST_Radius

  use scale_element_quadrilateral
  use scale_localmesh_2d
  use scale_mesh_cubedspheredom2d

  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_cubedspheredom2d, only: MeshFieldCommCubedSphereDom2D

  use scale_time_manager, only: &
    TIME_manager_Init , TIME_manager_Final,            &
    TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP,            &
    TIME_DTSEC, TIME_NSTEP 
    
  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_setup,   &
    FILE_HISTORY_meshfield_put,     &
    FILE_HISTORY_meshfield_write,   &
    FILE_HISTORY_meshfield_finalize

  implicit none

  integer, parameter :: PolyOrder = 2
  integer, parameter :: NeGX = 3
  integer, parameter :: NeGY = 3
  integer, parameter :: NLocalMeshPerPrc = 6

  type(QuadrilateralElement) :: refElem  
  type(MeshCubedSphereDom2D), target :: mesh

  type(MeshField2D), target :: lon
  type(MeshField2D), target :: lat
  type(MeshField2D), target :: q
  character(len=H_short), parameter :: q_VARNAME = "q"
  character(len=H_short), parameter :: q_DESC    = "q"
  character(len=H_short), parameter :: q_UNITS   = "K"
  integer :: HST_ID(3)

  type(MeshFieldCommCubedSphereDom2D) :: fields_comm
  type(MeshFieldContainer) :: field_list(1)  

  integer :: n
  type(LocalMesh2D), pointer :: lcmesh

  !---------------
  
  call init()

  !----
  write(*,*) "* Check data communication.."
  call perform_comm()

  do n=1, mesh%LOCAL_MESH_NUM
    lcmesh => mesh%lcmesh_list(n)
    
    write(*,*) "Check interior & halo data.."
    write(*,*) "tileID=", lcmesh%tileID
    ! call check_interior_data(n, mesh%lcmesh_list(n), q%local(n)%val, q%local(n)%val)
    ! call check_halo_data(n, mesh%lcmesh_list(n), q%local(n)%val, q%local(n)%val)
  end do
  !----

  call final()

contains
!OCL SERIAL
  elemental function get_field_val(tileID, k, p) result(val)
    implicit none
    integer, intent(in) :: tileID, k, p
    real(RP) :: val
    !----------------------------
    val = tileID*1000000 + k*1000 + p

    return
  end function get_field_val

!OCL SERIAL
  subroutine init()
    use scale_calendar, only: CALENDAR_setup
    use scale_const, only: CONST_setup
    implicit none

    integer :: comm, myrank, nprocs
    logical :: ismaster
    
    integer :: domid
    integer :: ke, p
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

    ! setup constant
    call CONST_setup

    !------
    call refElem%Init(PolyOrder, .true.)

    call mesh%Init( &
      NeGX, NeGY,                                        &
      RPlanet, refElem, NLocalMeshPerPrc, nprocs, myrank )
    
    call mesh%Generate()

    !---
    call q%Init( q_VARNAME, q_UNITS, mesh )
    call lon%Init( "lon", "rad", mesh )
    call lat%Init( "lat", "rad", mesh )
    call fields_comm%Init(1, 0, 0, mesh)
    
    call FILE_HISTORY_meshfield_setup( meshCubedSphere2D_ = mesh )
    call FILE_HISTORY_reg( q_VARNAME, q_DESC, q_UNITS, HST_ID(1), dim_type='XY')
    call FILE_HISTORY_reg( lon%varname, "longitude", lon%unit, HST_ID(2), dim_type='XY')
    call FILE_HISTORY_reg( lat%varname,  "latitude", lat%unit, HST_ID(3), dim_type='XY')

    !---
    do domid=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(domid)

      !$omp parallel do private(p)
      do ke=lcmesh%NeS, lcmesh%NeE
        do p=1, refElem%Np
          q%local(domid)%val(p,ke) = get_field_val(lcmesh%tileID, ke, p)
        end do
        lon%local(domid)%val(:,ke) = lcmesh%lon(:,ke)
        lat%local(domid)%val(:,ke) = lcmesh%lat(:,ke)
      end do
    end do

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), lon)
    call FILE_HISTORY_meshfield_put(HST_ID(3), lat)
    call FILE_HISTORY_meshfield_write()    

    !---

    return
  end subroutine init

!OCL SERIAL
  subroutine final()
    implicit none

    call FILE_HISTORY_meshfield_finalize()

    call q%Final()
    call lon%Final()
    call lat%Final()
    
    call fields_comm%Final()
    call mesh%Final()
    call refElem%Final()
    
    call TIME_manager_Final()
    call PRC_MPIfinish()

    return
  end subroutine final

!OCL SERIAL
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
    
!OCL SERIAL
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

!OCL SERIAL
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
      'check_halo_data', 'q_halo_east', refElem%Nfp*lcmesh_%NeY)

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
      'check_halo_data', 'q_halo_west', refElem%Nfp*lcmesh_%NeY)    

  end subroutine check_halo_data

!OCL SERIAL
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

end program test_field_cubedspheredom2d
