#include "scalelib.h"
program test_field3d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale_file_history
  use scale

  use scale_element_hexahedral
  use scale_localmesh_3d
  use scale_mesh_cubedom3d

  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D

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

  integer, parameter :: NeGX = 3
  integer, parameter :: NeGY = 3
  integer, parameter :: NeGZ = 3
  integer, parameter :: NLocalMeshPerPrc = 1

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: dom_ymin = -1.0_RP
  real(RP), parameter :: dom_ymax = +1.0_RP
  real(RP), parameter :: dom_zmin = -1.0_RP
  real(RP), parameter :: dom_zmax = +1.0_RP

  type(HexahedralElement) :: refElem
  integer, parameter :: PolyOrder_h = 2
  integer, parameter :: PolyOrder_v = 2
  
  type(MeshCubeDom3D), target :: mesh

  type(MeshField3D), target :: q
  character(len=H_short), parameter :: q_VARNAME = "q"
  character(len=H_short), parameter :: q_DESC    = "q"
  character(len=H_short), parameter :: q_UNITS   = "K"
  integer :: HST_ID(1)

  type(MeshFieldCommCubeDom3D) :: fields_comm

  integer :: n
  type(LocalMesh3D), pointer :: lcmesh

  !---------------
  
  call init()

  !----
  ! write(*,*) "* Check data communication.."
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
    val = dble(tileID*1000000 + k*1000 + p)
    return
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
    call refElem%Init(PolyOrder_h, PolyOrder_v, .true.)

    call mesh%Init( &
      NeGX, NeGY, NeGZ,                                           &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
      .true., .true., .true.,                                     &
      refElem, NLocalMeshPerPrc, nprocs, 1 )
    
    call mesh%Generate()

    !---
    call q%Init( q_VARNAME, q_UNITS, mesh )
    call fields_comm%Init(1, 0, mesh)
    
    call FILE_HISTORY_meshfield_setup( mesh3d_=mesh )
    call FILE_HISTORY_reg( q_VARNAME, q_DESC, q_UNITS, HST_ID(1), dim_type='XYZ')

    !---

    !---
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      q%local(n)%val(:,:) = -1.0E2_RP
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

    field_list(1)%field3d => q
    write(*,*) " - Put.."
    call fields_comm%Put(field_list, 1)
    write(*,*) " - Exchange.."
    call fields_comm%Exchange()
    write(*,*) " - Get.."
    call fields_comm%Get(field_list, 1)

  end subroutine perform_comm
    
  subroutine check_interior_data(n, lcmesh_, lcfield, lcfield_fl)
    integer, intent(in) :: n
    type(LocalMesh3D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ke, f, p
    integer :: iM, iP, fnid
    integer :: ans(refElem%Np)
    !------------------------------

    do ke=lcmesh_%NeS, lcmesh_%NeE
      do p=1, refElem%Np
        ans(p) = get_field_val(lcmesh_%tileID,ke,p)
      end do

      call assert(ke, int(lcfield(:,ke)), ans(:), 'check_interior_data', 'q', refElem%Np)
      do f=1, refElem%Nfaces_h
      do p=1, refElem%Nfp_h
        fnid = p + (f-1)*refElem%Nfp_h
        iM = lcmesh_%VMapM(fnid,ke)
        iP = lcmesh_%VMapP(fnid,ke)
        write(*,'(a,3(i3),a,i10,a,i10,a,3i4)') "ke, f, p =", ke,f,p,           &
          " : q- = ", int(lcfield_fl(iM)), ",  q+ =", int(lcfield_fl(iP)),    &
          ", normal=", int(lcmesh_%normal_fn(fnid,ke,:))
      end do
      end do
      do f=1, refElem%Nfaces_v
      do p=1, refElem%Nfp_v
        fnid = p + refElem%Nfaces_h*refElem%Nfp_h + (f-1)*refElem%Nfp_v
        iM = lcmesh_%VMapM(fnid,ke)
        iP = lcmesh_%VMapP(fnid,ke)
        write(*,'(a,3(i3),a,i10,a,i10,a,3i4)') "ke, f, p =", ke,f+refElem%Nfaces_h,p,  &
          " : q- = ", int(lcfield_fl(iM)), ",  q+ =", int(lcfield_fl(iM)),            &
          ", normal=", int(lcmesh_%normal_fn(fnid,ke,:))
      end do
      end do
    end do

  end subroutine check_interior_data

  subroutine check_halo_data(n, lcmesh_, lcfield, lcfield_fl)
    implicit none

    integer, intent(in) :: n
    type(LocalMesh3D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ans_bx(refElem%Nfp_h,lcmesh_%NeX,lcmesh_%NeZ)
    integer :: ans_by(refElem%Nfp_h,lcmesh_%NeY,lcmesh_%NeZ)    
    integer :: ans_bz(refElem%Nfp_v,lcmesh_%NeX,lcmesh_%NeY) 
    integer :: haloInd_s, haloInd_e
    integer :: ke, f, fp
    integer :: i, j, k
    integer :: tileID

    integer :: Ne, NeX, NeY, NeZ
    integer :: Np, Nfp_h, Nfp_v

    !------------------------------

    tileID = lcmesh_%tileID
    Ne  = lcmesh_%Ne
    NeX = lcmesh_%NeX
    NeY = lcmesh_%NeY
    NeZ = lcmesh_%NeZ
    Np = refElem%Np
    Nfp_h = refElem%Nfp_h
    Nfp_v = refElem%Nfp_v

    write(*,*) "Check values in halo.."
    write(*,*) "tileID_globalMap=", mesh%tileID_globalMap(:,tileID)
    write(*,*) "tileFaceID_globalMap=", mesh%tileFaceID_globalMap(:,tileID)


    haloInd_s = Np*Ne + 1    
    haloInd_e = haloInd_s + Nfp_h*NeX*NeZ - 1
    do k=1, NeZ
    do i=1, NeX
      ke = i + (NeY-1)*NeX  + (k-1)*NeX*NeY
      do fp=1, Nfp_h
        ans_bx(fp,i,k) = get_field_val(mesh%tileID_globalMap(1,tileID), ke, refElem%Fmask_h(fp,3))
      end do
    end do
    end do
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans_bx(:,:,:),     &
                 'check_halo_data', 'q_halo_south', Nfp_h*NeX*NeZ )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + Nfp_h*NeY*NeZ - 1
    do k=1, NeZ
    do j=1, NeY
      ke = 1 + (j-1)*NeX + (k-1)*NeX*NeY
      do fp=1, Nfp_h
        ans_by(fp,j,k) = get_field_val(mesh%tileID_globalMap(2,tileID), ke, refElem%Fmask_h(fp,4))
      end do
    end do
    end do
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans_by(:,:,:),     &
      'check_halo_data', 'q_halo_east', Nfp_h*NeY*NeZ )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + Nfp_h*NeX*NeZ - 1    
    do k=1, NeZ
    do i=1, NeX
      ke = i + (k-1)*NeX*NeY
      do fp=1, Nfp_h
        ans_bx(fp,i,k) = get_field_val(mesh%tileID_globalMap(3,tileID), ke, refElem%Fmask_h(fp,1))
      end do
    end do
    end do
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans_bx(:,:,:),     &
                 'check_halo_data', 'q_halo_north', Nfp_h*NeX*NeZ )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + Nfp_h*NeY*NeZ - 1    
    do k=1, NeZ
    do j=1, NeY
      ke = NeX + (j-1)*NeX + (k-1)*NeX*NeY
      do fp=1, Nfp_h
        ans_by(fp,j,k) = get_field_val(mesh%tileID_globalMap(4,tileID), ke, refElem%Fmask_h(fp,2))
      end do
    end do
    end do
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans_by(:,:,:),     &
      'check_halo_data', 'q_halo_west', Nfp_h*NeY*NeZ )    

    haloInd_s = haloInd_e + 1   
    haloInd_e = haloInd_s + Nfp_v*NeX*NeY - 1
    do j=1, NeY
    do i=1, NeX
      ke = i + (j-1)*NeX + (NeZ-1)*NeX*NeY
      do fp=1, Nfp_v
        ans_bz(fp,i,j) = get_field_val(mesh%tileID_globalMap(5,tileID), ke, refElem%Fmask_v(fp,2))
      end do
    end do
    end do
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans_bz(:,:,:),     &
                  'check_halo_data', 'q_halo_bottom', Nfp_v*NeX*NeY )

    haloInd_s = haloInd_e + 1   
    haloInd_e = haloInd_s + Nfp_v*NeX*NeY - 1
    do j=1, NeY
    do i=1, NeX
      ke = i + (j-1)*NeX
      do fp=1, Nfp_v
        ans_bz(fp,i,j) = get_field_val(mesh%tileID_globalMap(6,tileID), ke, refElem%Fmask_v(fp,1))
      end do
    end do
    end do
    call assert( -1, int(lcfield_fl(haloInd_s:haloInd_e)), ans_bz(:,:,:),     &
                  'check_halo_data', 'q_halo_top', Nfp_v*NeX*NeY ) 

    return
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
    if ( sum(dble(vals(:) - ans(:))**2) > EPS ) then
      LOG_ERROR(assert_name,*) 'The value of '//trim(var_name)//' is unexcepted!', &
        ' k=', k, ": val=", vals(:), " ans=", ans(:)
      call PRC_abort
    end if    
  end subroutine assert

end program test_field3d
