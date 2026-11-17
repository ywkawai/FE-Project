#include "scalelib.h"
program test_field3d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale_file_history
  use scale_const, only: &
    PI => CONST_PI
  use scale

  use scale_sparsemat
  use scale_element_base
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

  use mod_test_meshfield_util, only: &
    TEST_field_val, TEST_perform_comm, TEST_check_face_node
  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: UNITTEST_assert_equal, UNITTEST_assert_less
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
  type(SparseMat) :: Dx, Dy, Dz, Lift
  
  type(MeshCubeDom3D), target :: mesh

  type(MeshField3D), target :: q
  character(len=H_short), parameter :: q_VARNAME = "q"
  character(len=H_short), parameter :: q_DESC    = "q"
  character(len=H_short), parameter :: q_UNITS   = "K"
  type(MeshField3D), target :: u
  type(MeshField3D), target :: dudx
  type(MeshField3D), target :: dudy
  type(MeshField3D), target :: dudz
  integer, save :: HST_ID(5)

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
    !$acc update host(q%local(n)%val)
    call check_interior_data(n, mesh%lcmesh_list(n), q%local(n)%val, q%local(n)%val)
    call check_halo_data(n, mesh%lcmesh_list(n), q%local(n)%val, q%local(n)%val)

    write(*,*) "Check gradient operator.."
    write(*,*) "tileID=", lcmesh%tileID
    call calc_grad( dudx%local(n)%val, dudy%local(n)%val, dudz%local(n)%val, &
      u%local(n)%val, lcmesh, lcmesh%refElem3D )
    !$acc update host(dudx%local(n)%val, dudy%local(n)%val, dudz%local(n)%val)
    call check_grad( dudx%local(n)%val, dudy%local(n)%val, dudz%local(n)%val, &
      lcmesh, refElem )

  end do
  !----
  
  call final()

contains

  subroutine get_field_val(val, lmesh, elem)
    class(LocalMesh3D), intent(in) :: lmesh
    class(HexahedralElement), intent(in) :: elem
    real(RP), intent(out) :: val(elem%Np,lmesh%NeA)

    integer :: tileID, ke, p
    !----------------------------

    tileID = lmesh%tileID
    !$acc parallel loop collapse(2) present(lmesh,val)
    do ke=lmesh%NeS, lmesh%NeE
    do p=1, elem%Np
      val(p,ke) = TEST_field_val(tileID, ke, p)
    end do
    end do
  end subroutine get_field_val

  subroutine get_field_u(val, lmesh, elem)
    class(LocalMesh3D), intent(in) :: lmesh
    class(HexahedralElement), intent(in) :: elem
    real(RP), intent(out) :: val(elem%Np,lmesh%NeA)

    integer :: ke, p
    !----------------------------

    !$acc parallel loop collapse(2) present(lmesh,val)
    do ke=lmesh%NeS, lmesh%NeE
    do p=1, elem%Np
      val(p,ke) = sin( PI * lmesh%pos_en(p,ke,1) ) &
                * sin( PI * lmesh%pos_en(p,ke,2) ) &
                * sin( PI * lmesh%pos_en(p,ke,3) ) 
    end do
    end do
  end subroutine get_field_u


  subroutine init()
    implicit none

    integer :: comm, myrank, nprocs
    logical :: ismaster
    !---------

    call UNITTEST_init( "test_field_cubedom3d_hexahedral",          &
      comm=comm, myrank=myrank, nprocs=nprocs, ismaster=ismaster, use_time=.true. )

    !------
    call refElem%Init(PolyOrder_h, PolyOrder_v, .true.)
    call Dx%Init( refElem%Dx1 )
    call Dy%Init( refElem%Dx2 )
    call Dz%Init( refElem%Dx3 )
    call Lift%Init( refElem%Lift )

    call mesh%Init( &
      NeGX, NeGY, NeGZ,                                           &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
      .true., .true., .true.,                                     &
      refElem, NLocalMeshPerPrc, nprocs, 1 )
    
    call mesh%Generate()

    !---
    call q%Init( q_VARNAME, q_UNITS, mesh )
    call u   %Init( "u"   , "1", mesh )
    call dudx%Init( "dudx", "1", mesh )
    call dudy%Init( "dudy", "1", mesh )
    call dudz%Init( "dudz", "1", mesh )
    call fields_comm%Init(2, 0, 0, mesh)

    call FILE_HISTORY_meshfield_setup( mesh3d_=mesh )
    call FILE_HISTORY_reg( q_VARNAME, q_DESC, q_UNITS, HST_ID(1), dim_type='XYZ')
    call FILE_HISTORY_reg( "u"   ,    "u", "1", HST_ID(2), dim_type='XYZ')
    call FILE_HISTORY_reg( "dudx", "dudx", "1", HST_ID(3), dim_type='XYZ')
    call FILE_HISTORY_reg( "dudy", "dudy", "1", HST_ID(4), dim_type='XYZ')
    call FILE_HISTORY_reg( "dudz", "dudz", "1", HST_ID(5), dim_type='XYZ')

    !---
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)

      call get_field_val( q%local(n)%val, &
        mesh%lcmesh_list(n), refElem )
      call get_field_u( u%local(n)%val, &
        mesh%lcmesh_list(n), refElem )
    end do
    !---

    return
  end subroutine init

  subroutine final()
    implicit none
    !------------------------------------------------------

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), u)
    call FILE_HISTORY_meshfield_put(HST_ID(3), dudx)
    call FILE_HISTORY_meshfield_put(HST_ID(4), dudy)
    call FILE_HISTORY_meshfield_put(HST_ID(5), dudz)

    call FILE_HISTORY_meshfield_write()   
    call FILE_HISTORY_meshfield_finalize()

    call q%Final()
    call u%Final()
    call dudx%Final()
    call dudy%Final()
    call dudz%Final()
    call fields_comm%Final()

    call mesh%Final()
    call Dx%Final(); call Dy%Final(); call Dz%Final()
    call Lift%Final()    
    call refElem%Final()
    
    call UNITTEST_final()

    return
  end subroutine final

  subroutine perform_comm()
    implicit none
    type(MeshFieldContainer) :: field_list(2)
    !---------------------------

    field_list(1)%field3d => q
    field_list(2)%field3d => u
    call TEST_perform_comm( fields_comm, field_list )

    return
  end subroutine perform_comm
    
  subroutine check_interior_data(n, lcmesh_, lcfield, lcfield_fl)
    implicit none

    integer, intent(in) :: n
    type(LocalMesh3D), intent(in) :: lcmesh_
    real(RP), intent(in) :: lcfield(refElem%Np,lcmesh_%NeA)
    real(RP), intent(in) :: lcfield_fl(refElem%Np*lcmesh_%NeA)
    
    integer :: ke, f, p
    integer :: fnid
    real(RP) :: ans(refElem%Np,lcmesh_%NeA)
    !------------------------------

    !$acc data create(ans)
    call get_field_val( ans, &
      lcmesh_, refElem )
    !$acc update host(ans)
    
    do ke=lcmesh_%NeS, lcmesh_%NeE
      call UNITTEST_assert_equal( int(lcfield(:,ke)), int(ans(:,ke)), 'check_interior_data: q', idx=ke )
      do f=1, refElem%Nfaces_h
      do p=1, refElem%Nfp_h
        fnid = p + (f-1)*refElem%Nfp_h
        call TEST_check_face_node( lcmesh_, lcfield_fl, ke, f, fnid )
      end do
      end do
      do f=1, refElem%Nfaces_v
      do p=1, refElem%Nfp_v
        fnid = p + refElem%Nfaces_h*refElem%Nfp_h + (f-1)*refElem%Nfp_v
        call TEST_check_face_node( lcmesh_, lcfield_fl, ke, f+refElem%Nfaces_h, fnid )
      end do
      end do
    end do
    !$acc end data
    return
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
        ans_bx(fp,i,k) = TEST_field_val(mesh%tileID_globalMap(1,tileID), ke, refElem%Fmask_h(fp,3))
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_bx(:,:,:), &
      'check_halo_data: q_halo_south' )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + Nfp_h*NeY*NeZ - 1
    do k=1, NeZ
    do j=1, NeY
      ke = 1 + (j-1)*NeX + (k-1)*NeX*NeY
      do fp=1, Nfp_h
        ans_by(fp,j,k) = TEST_field_val(mesh%tileID_globalMap(2,tileID), ke, refElem%Fmask_h(fp,4))
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_by(:,:,:), &
      'check_halo_data: q_halo_east' )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + Nfp_h*NeX*NeZ - 1    
    do k=1, NeZ
    do i=1, NeX
      ke = i + (k-1)*NeX*NeY
      do fp=1, Nfp_h
        ans_bx(fp,i,k) = TEST_field_val(mesh%tileID_globalMap(3,tileID), ke, refElem%Fmask_h(fp,1))
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_bx(:,:,:), &
      'check_halo_data: q_halo_north' )

    haloInd_s = haloInd_e + 1    
    haloInd_e = haloInd_s + Nfp_h*NeY*NeZ - 1    
    do k=1, NeZ
    do j=1, NeY
      ke = NeX + (j-1)*NeX + (k-1)*NeX*NeY
      do fp=1, Nfp_h
        ans_by(fp,j,k) = TEST_field_val(mesh%tileID_globalMap(4,tileID), ke, refElem%Fmask_h(fp,2))
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_by(:,:,:), &
      'check_halo_data: q_halo_west' )    

    haloInd_s = haloInd_e + 1   
    haloInd_e = haloInd_s + Nfp_v*NeX*NeY - 1
    do j=1, NeY
    do i=1, NeX
      ke = i + (j-1)*NeX + (NeZ-1)*NeX*NeY
      do fp=1, Nfp_v
        ans_bz(fp,i,j) = TEST_field_val(mesh%tileID_globalMap(5,tileID), ke, refElem%Fmask_v(fp,2))
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_bz(:,:,:), &
      'check_halo_data: q_halo_bottom' )

    haloInd_s = haloInd_e + 1   
    haloInd_e = haloInd_s + Nfp_v*NeX*NeY - 1
    do j=1, NeY
    do i=1, NeX
      ke = i + (j-1)*NeX
      do fp=1, Nfp_v
        ans_bz(fp,i,j) = TEST_field_val(mesh%tileID_globalMap(6,tileID), ke, refElem%Fmask_v(fp,1))
      end do
    end do
    end do
    call UNITTEST_assert_equal( int(lcfield_fl(haloInd_s:haloInd_e)), ans_bz(:,:,:), &
      'check_halo_data: q_halo_top' ) 

    return
  end subroutine check_halo_data

  !> Compare the DG gradient of u = sin(pi x) sin(pi y) sin(pi z) with the analytical solution
  subroutine check_grad( dudx_, dudy_, dudz_, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(HexahedralElement), intent(in) :: elem
    real(RP), intent(in) :: dudx_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: dudy_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: dudz_(elem%Np,lmesh%NeA)

    integer :: ke, d
    real(RP) :: x(elem%Np), y(elem%Np), z(elem%Np), w(elem%Np)
    real(RP) :: grad_ans(elem%Np,3)
    real(RP) :: l2err(3), l2norm(3), rel_err(3)
    character(len=4), parameter :: lbl(3) = (/ 'dudx', 'dudy', 'dudz' /)

    ! Relative L2 error tolerance for PolyOrder=2, NeGX=NeGY=NeGZ=3
    ! (measured value: 2.413E-01 for each component)
    real(RP), parameter :: GRAD_ERR_TOL = 0.25_RP
    !------------------------------------------------

    l2err(:) = 0.0_RP; l2norm(:) = 0.0_RP
    do ke=lmesh%NeS, lmesh%NeE
      x(:) = lmesh%pos_en(:,ke,1); y(:) = lmesh%pos_en(:,ke,2); z(:) = lmesh%pos_en(:,ke,3)
      w(:) = lmesh%J(:,ke) * elem%IntWeight_lgl(:)
      grad_ans(:,1) = PI * cos( PI * x(:) ) * sin( PI * y(:) ) * sin( PI * z(:) )
      grad_ans(:,2) = PI * sin( PI * x(:) ) * cos( PI * y(:) ) * sin( PI * z(:) )
      grad_ans(:,3) = PI * sin( PI * x(:) ) * sin( PI * y(:) ) * cos( PI * z(:) )

      l2err(1) = l2err(1) + sum( w(:) * ( dudx_(:,ke) - grad_ans(:,1) )**2 )
      l2err(2) = l2err(2) + sum( w(:) * ( dudy_(:,ke) - grad_ans(:,2) )**2 )
      l2err(3) = l2err(3) + sum( w(:) * ( dudz_(:,ke) - grad_ans(:,3) )**2 )
      do d=1, 3
        l2norm(d) = l2norm(d) + sum( w(:) * grad_ans(:,d)**2 )
      end do
    end do
    rel_err(:) = sqrt( l2err(:) / l2norm(:) )

    do d=1, 3
      write(*,'(a,a,a,es12.5)') " ", lbl(d), ": relative L2 error=", rel_err(d)
      call UNITTEST_assert_less( rel_err(d), GRAD_ERR_TOL, 'check_grad: relative L2 error of '//lbl(d) )
    end do

    return
  end subroutine check_grad

  subroutine calc_grad( dudx_, dudy_, dudz_, &
      u_, lmesh, elem                        )
    
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: dudx_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: dudy_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: dudz_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: u_(elem%Np,lmesh%NeA)

    integer :: ke
    real(RP) :: del_flux(elem%NfpTot,lcmesh%Ne,3)
    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlux(elem%Np)
    !------------------------------------------------

    !$acc data create(del_flux)

    call calc_grad_delflx( del_flux,                                              &
      u_, lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), &
      lmesh%Fscale, lmesh%VMapM, lmesh%VMapP, lmesh, elem )
    
    !$omp parallel do private(Fx, Fy, Fz, LiftDelFlux)
    !$acc parallel loop private(Fx, Fy, Fz, LiftDelFlux) present(lmesh,u_,del_flux,dudx_,dudy_,dudz_)
    do ke=1, lmesh%Ne
      call sparsemat_matmul(Dx, u_(:,ke), Fx)
      call sparsemat_matmul(Lift, del_flux(:,ke,1), LiftDelFlux )
      dudx_(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlux(:)

      call sparsemat_matmul(Dy, u_(:,ke), Fy)
      call sparsemat_matmul(Lift, del_flux(:,ke,2), LiftDelFlux )
      dudy_(:,ke) = lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlux(:)

      call sparsemat_matmul(Dz, u_(:,ke), Fz)
      call sparsemat_matmul(Lift, del_flux(:,ke,3), LiftDelFlux )
      dudz_(:,ke) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlux(:)
    end do

    !$acc end data
    return
  end subroutine calc_grad

  subroutine calc_grad_delflx( del_flux,                &
      u_, nx, ny, nz, Fscale, VMapM, VMapP, lmesh, elem )
          
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: del_flux(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(in) :: u_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: Fscale(elem%NfpTot*lmesh%Ne) 
    integer, intent(in) :: VMapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: VMapP(elem%NfpTot*lmesh%Ne)
  
    integer :: i
    integer :: iM, iP
    real(RP) :: del
    !-----------------------------------------

    !$omp parallel do private(iM, iP, del)
    !$acc parallel loop private(iM, iP, del) present(u_,nx,ny,nz,Fscale,VMapM,VMapP,del_flux)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = VMapM(i); iP = VMapP(i)
      del = 0.5_RP * (u_(iP) - u_(iM)) * Fscale(i)

      del_flux(i,1) = del * nx(i)
      del_flux(i,2) = del * ny(i)
      del_flux(i,3) = del * nz(i)
    end do

    return
  end subroutine calc_grad_delflx


end program test_field3d
