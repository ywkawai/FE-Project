!-------------------------------------------------------------------------------
!> module FElib / Mesh / utility for 3D mesh
!!
!! @par Description
!!          A module useful for generating 3D mesh 
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_meshutil_3d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !  
  public :: MeshUtil3D_genCubeDomain
  public :: MeshUtil3D_genConnectivity
  public :: MeshUtil3D_buildInteriorMap
  public :: MeshUtil3D_genPatchBoundaryMap
  public :: MeshUtil3D_buildGlobalMap

contains

!OCL SERIAL
  subroutine MeshUtil3D_genCubeDomain( pos_v, EToV,        &
    Ke_x, xmin, xmax, Ke_y, ymin, ymax, Ke_z, zmin, zmax,  &
    Fz )

    implicit none

    integer, intent(in) :: Ke_x
    real(RP), intent(in) :: xmin, xmax
    integer, intent(in) :: Ke_y
    real(RP), intent(in) :: ymin, ymax
    integer, intent(in) :: Ke_z
    real(RP), intent(in) :: zmin, zmax
    real(RP), intent(in), optional :: Fz(Ke_z+1)
    
    real(RP), intent(out) :: pos_v((Ke_x+1)*(Ke_y+1)*(Ke_z+1),3)
    integer, intent(out) :: EToV(Ke_x*Ke_y*Ke_z,8)

    integer :: i, j, k
    integer :: n
    integer :: k_
    integer :: NvX, NvY, NvZ
    !-----------------------------------------------------------------------------

    NvX = Ke_X + 1
    NvY = Ke_Y + 1
    NvZ = Ke_Z + 1

    !$omp parallel private(i,j,k,n)

    !$omp do collapse(2) 
    do k=1, NvZ
    do j=1, NvY
    do i=1, NvX
      n = i + (j-1)*NvX + (k-1)*NvX*NvY
      pos_v(n,1) = (xmax - xmin)*dble(i - 1)/dble(NvX - 1) + xmin
      pos_v(n,2) = (ymax - ymin)*dble(j - 1)/dble(NvY - 1) + ymin
      if ( present(Fz) ) then
        pos_v(n,3) = Fz(k)
      else
        pos_v(n,3) = (zmax - zmin)*dble(k - 1)/dble(NvZ - 1) + zmin
      end if
    end do
    end do
    end do

    !$omp do collapse(2)
    do k=1, Ke_z
    do j=1, Ke_y
    do i=1, Ke_x
      n = i + (j-1)*Ke_x + (k-1)*Ke_x*Ke_y
      EToV(n,1:4) = (k-1)*NvX*NvY + (j-1)*NvX + &
            & (/ i, i+1,              &
            &    i + NvX, i+1 + NvX /)
      EToV(n,5:8) = k    *NvX*NvY + (j-1)*NvX + &
            & (/ i, i+1,              &
            &    i + NvX, i+1 + NvX /)              
    end do
    end do
    end do

    !$omp end parallel
    
    !---
  !!$
  !!$    write(*,*) "-- vx, vy --"
  !!$    do j=1, NvY
  !!$    do i=1, NvX
  !!$       n = i + (j-1)*NvY
  !!$       write(*,*) i, j, " :", mesh%vx(n), mesh%vy(n)
  !!$    end do
  !!$    end do
  !!$
  !!$    write(*,*) "-- EToV --"
  !!$    do j=1, Ke_y
  !!$    do i=1, Ke_x
  !!$       n = i + (j-1)*Ke_y
  !!$       write(*,*) i, j, " :", mesh%EToV(n,:)
  !!$    end do
  !!$    end do

    return
  end subroutine MeshUtil3D_genCubeDomain

!OCL SERIAL
  subroutine MeshUtil3D_genConnectivity( EToE, EToF, &
    EToV, Ne, Nfaces )
    
    use scale_quicksort, only: QUICKSORT_exec_with_idx
    implicit none

    integer, intent(in) :: Ne
    integer, intent(in) :: Nfaces
    integer, intent(out) :: EToE(Ne, Nfaces)
    integer, intent(out) :: EToF(Ne, Nfaces)
    integer, intent(in) :: EToV(Ne,8)

    integer :: nodes(Ne*Nfaces,4)
    integer(kind=8) :: face_ids(Ne*Nfaces)
    integer :: ke
    integer :: f
    integer :: n
    integer :: n1, n2
    integer :: Nnodes
    integer :: tmp
    integer :: Nnodes_row
    integer :: matchL(2,3), matchR(2,3)

    real(RP) :: EToE_1d(Ne*Nfaces)
    real(RP) :: EToF_1d(Ne*Nfaces)

    integer :: spNodeToNode(3,Ne*Nfaces)
    integer :: spNodeToNodeRowTmp(3,Ne*Nfaces)
    integer :: sort_indx(Ne*Nfaces)
    integer(kind=8) :: sorted_faceid(Ne*Nfaces)

    real(RP) :: vtmp(4)
    
    !-----------------------------------------------------------------------------

    Nnodes = maxval( EToV )
    Nnodes_row = size(nodes,1)
    
    !---------
    do ke=1, Ne
       nodes(ke     ,:) = EToV(ke,(/ 1, 2, 6, 5 /))
       nodes(ke+Ne  ,:) = EToV(ke,(/ 2, 4, 8, 6 /))
       nodes(ke+2*Ne,:) = EToV(ke,(/ 4, 3, 7, 8 /))
       nodes(ke+3*Ne,:) = EToV(ke,(/ 3, 1, 5, 7 /))
       nodes(ke+4*Ne,:) = EToV(ke,(/ 1, 3, 4, 2 /))       
       nodes(ke+5*Ne,:) = EToV(ke,(/ 5, 6, 8, 7 /))       
    end do
   
    ! Sort
    do n=1, Nnodes_row
      vtmp(:) = nodes(n,:)
      call bubbleSort( vtmp )
      nodes(n,:) = vtmp(:)
      ! write(*,*) n, ":", nodes(n,:)
    end do
    nodes = nodes - 1
    !---------

    do ke=1, Ne
       EToE(ke,:) = ke
       EToF(ke,:) = (/ 1, 2, 3, 4, 5, 6 /)
    end do

    face_ids(:) =   nodes(:,1)*Nnodes**3 + nodes(:,2)*Nnodes**2 &
                + nodes(:,3)*Nnodes + nodes(:,4) + 1

    do f=1, Nfaces
    do ke=1, Ne
      n = ke + (f-1)*Ne
      spNodeToNode(:,n) = (/ n, EToE(ke,f), EToF(ke,f) /)
      sorted_faceid(n)  = face_ids(n)
      sort_indx(n)      = n
       ! write(*,*) "face_id, n, EToE, EToF:", spNodeToNode(:,n)
    end do
    end do
    
    !- sort row
    call QUICKSORT_exec_with_idx( Ne*Nfaces, sorted_faceid, sort_indx )
    spNodeToNodeRowTmp(:,:) = spNodeToNode(:,:)

    do n=1, Nnodes_row
      spNodeToNode(:,n) = spNodeToNodeRowTmp(:,sort_indx(n))
      !  write(*,'(a,4i10)') "(sorted) face_id, n, EToE, EToF:", spNodeToNode(:,n)
    end do

    EToE_1d(:) = -1
    EToF_1d(:) = -1
    do n=1, Nnodes_row-1
       if ( sorted_faceid(n) - sorted_faceid(n+1) == 0 ) then
          matchL(:,:) = transpose( spNodeToNode(:,(/ n, n+1 /)) )
          matchR(:,:) = transpose( spNodeToNode(:,(/ n+1, n /)) )

          EToE_1d(matchL(:,1)) = matchR(:,2)
          EToF_1d(matchL(:,1)) = matchR(:,3)
       end if
    end do

    do f=1, Nfaces
    do ke=1, Ne
       n = ke + (f-1)*Ne
       if ( EToE_1d(n) /= -1 ) then
          EToE(ke,f) = EToE_1d(n)
          EToF(ke,f) = EToF_1d(n)
       end if
    end do
    end do

    !-------------------

!    write(*,*) "EToE------"
!    do k=1, mesh%Ne
!      write(*,*) "k=", k, ":", mesh%EToE(k,:)
!    end do
!    write(*,*) "EToF------"
!    do k=1, mesh%Ne
!       write(*,*) "k=", k, ":", mesh%EToF(k,:)
!    end do

    return
  end subroutine MeshUtil3D_genConnectivity

!OCL SERIAL
  subroutine bubbleSort( array )
    implicit none
    real(RP), intent(inout) :: array(:)
    
    integer :: i, j, N
    real(RP) :: t
    !--------------------

    N = size(array)
    do i=1, N-1
      do j=i+1,N
        t = array(i)      
        if(t > array(j)) then
          array(i) = array(j)
          array(j) = t
        end if
      end do    
    end do

    return
  end subroutine bubbleSort

!OCL SERIAL
  subroutine MeshUtil3D_BuildInteriorMap( VMapM, VMapP, MapM, MapP, &
    pos_en, pos_ev, EtoE, EtoF, EtoV, Fmask_h, Fmask_v,             &
    Ne, Nv, Np, Nfp_h, Nfp_v, NfpTot, Nfaces_h, Nfaces_v, Nfaces)

    implicit none

    integer, intent(in) :: Ne
    integer, intent(in) :: Nv
    integer, intent(in) :: Np
    integer, intent(in) :: Nfp_h
    integer, intent(in) :: Nfp_v
    integer, intent(in) :: NfpTot
    integer, intent(in) :: Nfaces_h
    integer, intent(in) :: Nfaces_v
    integer, intent(in) :: Nfaces

    integer, intent(out) :: VMapM(NfpTot,Ne)
    integer, intent(out) :: VMapP(NfpTot,Ne)
    integer, intent(out) :: MapM(NfpTot,Ne)
    integer, intent(out) :: MapP(NfpTot,Ne)

    real(RP), intent(in) :: pos_en(Np,Ne,3)
    real(RP), intent(in) :: pos_ev(Nv,3)
    integer, intent(in) :: EToE(Ne,Nfaces)
    integer, intent(in) :: EToF(Ne,Nfaces)
    integer, intent(in) :: EToV(Ne,Nv)
    integer, intent(in) :: Fmask_h(Nfp_h,Nfaces_h)
    integer, intent(in) :: Fmask_v(Nfp_v,Nfaces_v)

    integer :: ke, ke1, ke2
    integer :: f, f1, f2
    integer :: p
    integer :: n
    integer :: i
    integer :: idP, idM
    integer :: v1, v2
    integer :: nodeids(Np,Ne)
    real(RP) :: r_h(Nfp_h,Nfp_h,3,2)
    real(RP) :: r_v(Nfp_v,Nfp_v,3,2)
    real(RP) :: dist_h(Nfp_h,Nfp_h)
    real(RP) :: dist_v(Nfp_v,Nfp_v)
    real(RP) :: x(Np*Ne), y(Np*Ne), z(Np*Ne)

    integer :: VMapM_h(Nfp_h,Nfaces_h,Ne)
    integer :: VMapP_h(Nfp_h,Nfaces_h,Ne)
    integer :: VMapM_v(Nfp_v,Nfaces_v,Ne)
    integer :: VMapP_v(Nfp_v,Nfaces_v,Ne)
    integer :: MapM_h(Nfp_h,Nfaces_h,Ne)
    integer :: MapP_h(Nfp_h,Nfaces_h,Ne)
    integer :: MapM_v(Nfp_v,Nfaces_v,Ne)
    integer :: MapP_v(Nfp_v,Nfaces_v,Ne)

    integer :: mindist_indx(1)
    !-----------------------------------------------------------------------------

    !$omp parallel private(ke,f,p,n)

    !$omp do
    do ke=1, Ne
    do p=1, Np
      n = p + (ke-1)*Np
      nodeids(p,ke) = n
      x(n) = pos_en(p,ke,1)
      y(n) = pos_en(p,ke,2)
      z(n) = pos_en(p,ke,3)
    end do
    end do
    !$omp end do

    !$omp do
    do ke=1, Ne
      do f=1, Nfaces_h
      do p=1, Nfp_h
        n = p + (f-1)*Nfp_h + (ke-1)*NfpTot
        MapM_h(p,f,ke) = n
        MapP_h(p,f,ke) = n
        VMapM_h(p,f,ke) = nodeids(Fmask_h(p,f),ke)
      end do
      end do
      do f=1, Nfaces_v
      do p=1, Nfp_v
        n = p + Nfaces_h*Nfp_h + (f-1)*Nfp_v + (ke-1)*NfpTot
        MapM_v(p,f,ke) = n
        MapP_v(p,f,ke) = n
        VMapM_v(p,f,ke) = nodeids(Fmask_v(p,f),ke)
      end do
      end do    
    end do
    !$omp end do

    !$omp workshare
    VMapP_h(:,:,:) = -1
    VMapP_v(:,:,:) = -1
    !$omp end workshare
    !$omp end parallel

    !$omp parallel private( &
    !$omp ke1, f1, ke2, f2, v1, v2,                         &
    !$omp r_h, r_v, dist_h, dist_v, mindist_indx, idP, idM  )

    !$omp do
    do ke1=1, Ne
    do f1=1, Nfaces_h
      ke2 = EToE(ke1,f1); f2 = EToF(ke1,f1)

      v1 = EToV(ke1,f1); v2 = EToV(ke1,1+mod(f1,Nfaces_h))

      r_h(:,:,1,1) = spread( x(VMapM_h(:,f1,ke1)), 2, Nfp_h )
      r_h(:,:,1,2) = spread( x(VMapM_h(:,f2,ke2)), 1, Nfp_h )
      r_h(:,:,2,1) = spread( y(VMapM_h(:,f1,ke1)), 2, Nfp_h )
      r_h(:,:,2,2) = spread( y(VMapM_h(:,f2,ke2)), 1, Nfp_h )
      r_h(:,:,3,1) = spread( z(VMapM_h(:,f1,ke1)), 2, Nfp_h )
      r_h(:,:,3,2) = spread( z(VMapM_h(:,f2,ke2)), 1, Nfp_h )
      
      dist_h(:,:) =   (r_h(:,:,1,1) - r_h(:,:,1,2))**2  &
                    + (r_h(:,:,2,1) - r_h(:,:,2,2))**2  &
                    + (r_h(:,:,3,1) - r_h(:,:,3,2))**2
      do idM=1, Nfp_h
        mindist_indx(:) = minloc(dist_h(idM,:))
        idP = mindist_indx(1)
        VMapP_h(idM,f1,ke1) = VMapM_h(idP,f2,ke2)
        MapP_h(idM,f1,ke1) = idP + (f2-1)*Nfp_h + (ke2-1)*NfpTot
      end do
    end do
    end do
    !omp end do

    !$omp do
    do ke1=1, Ne
    do f1=1, Nfaces_v
      ke2 = EToE(ke1,Nfaces_h+f1); f2 = EToF(ke1,Nfaces_h+f1) - Nfaces_h
      
      v1 = EToV(ke1,1); v2 = EToV(ke1,Nfaces_h+1)

      r_v(:,:,1,1) = spread( x(VMapM_v(:,f1,ke1)), 2, Nfp_v )
      r_v(:,:,1,2) = spread( x(VMapM_v(:,f2,ke2)), 1, Nfp_v )
      r_v(:,:,2,1) = spread( y(VMapM_v(:,f1,ke1)), 2, Nfp_v )
      r_v(:,:,2,2) = spread( y(VMapM_v(:,f2,ke2)), 1, Nfp_v )
      r_v(:,:,3,1) = spread( z(VMapM_v(:,f1,ke1)), 2, Nfp_v )
      r_v(:,:,3,2) = spread( z(VMapM_v(:,f2,ke2)), 1, Nfp_v )
      
      dist_v(:,:) =   (r_v(:,:,1,1) - r_v(:,:,1,2))**2  &
                    + (r_v(:,:,2,1) - r_v(:,:,2,2))**2  &
                    + (r_v(:,:,3,1) - r_v(:,:,3,2))**2
      do idM=1, Nfp_v
        mindist_indx(:) = minloc(dist_v(idM,:))
        idP = mindist_indx(1)
        VMapP_v(idM,f1,ke1) = VMapM_v(idP,f2,ke2)
        MapP_v(idM,f1,ke1) = idP + Nfaces_h*Nfp_h + (f2-1)*Nfp_v + (ke2-1)*NfpTot
      end do
    end do
    end do
    !omp end do
    !$omp end parallel

    !$omp parallel do private(ke,f,n,i)
    do ke=1, Ne
      do f=1, Nfaces_h
      do n=1, Nfp_h
        i = n + (f-1)*Nfp_h
        VMapM(i,ke) = VMapM_h(n,f,ke)
        MapM(i,ke) = MapM_h(n,f,ke)
        VMapP(i,ke) = VMapP_h(n,f,ke)
        MapP(i,ke) = MapP_h(n,f,ke)
      end do
      end do
      do f=1, Nfaces_v
      do n=1, Nfp_v
        i = n + Nfaces_h*Nfp_h + (f-1)*Nfp_v
        VMapM(i,ke) = VMapM_v(n,f,ke)
        MapM(i,ke) = MapM_v(n,f,ke)
        VMapP(i,ke) = VMapP_v(n,f,ke)
        MapP(i,ke) = MapP_v(n,f,ke)
      end do
      end do
    end do

    !-----
  !    mapB_counter = 0
  !    do k=1,mesh%Ne
  !    do f=1,elem%Nfaces
  !    do p=1,elem%Nfp
  !      n = p + (f-1)*elem%Nfp + (k-1)*elem%Nfp*elem%Nfaces
  !      if (mesh%VMapM(p,f,k) == mesh%VMapP(p,f,k)) then
  !        mapB_counter = mapB_counter + 1
  !        mapB_tmp(mapB_counter) = n
  !
  !        vmapB_tmp(mapB_counter) = mesh%VMapM(p,f,k)
  !        mesh%VMapP(p,f,k) = elem%Np*mesh%NeE + mapB_counter
  !      end if
  !    end do
  !    end do
  !    end do
  !
  !    allocate( mesh%mapB(mapB_counter) )
  !    allocate( mesh%vmapB(mapB_counter) )
  !    mesh%mapB(:) = mapB_tmp(1:mapB_counter)
  !    mesh%vmapB(:) = vmapB_tmp(1:mapB_counter)

    !-------
  !    write(*,*) "Build MapInfo: "
  !    do k=mesh%NeS,mesh%NeE
  !       write(*,*) "k=", k, "---"
  !       write(*,*) " VMapM:", mesh%VMapM(:,:,k)
  !       write(*,*) " VMapP:", mesh%VMapP(:,:,k)
  !    end do

  !    write(*,*) "mapB:", mesh%mapB(:)
  !    write(*,*) "vmapB:", mesh%vmapB(:)

    return
  end subroutine MeshUtil3D_BuildInteriorMap

!OCL SERIAL
  subroutine MeshUtil3D_genPatchBoundaryMap(  VMapB, MapB, VMapP,                  &
    pos_en, xmin, xmax, ymin, ymax, zmin, zmax,                                    &
    Fmask_h, Fmask_v, Ne, Nv, Np, Nfp_h, Nfp_v, NfpTot, Nfaces_h, Nfaces_v, Nfaces )

    implicit none

    integer, intent(in) :: Ne
    integer, intent(in) :: Nv
    integer, intent(in) :: Np
    integer, intent(in) :: Nfp_h
    integer, intent(in) :: Nfp_v
    integer, intent(in) :: NfpTot
    integer, intent(in) :: Nfaces_h
    integer, intent(in) :: Nfaces_v
    integer, intent(in) :: Nfaces

    integer, intent(inout), allocatable :: VMapB(:)
    integer, intent(inout), allocatable :: MapB(:)
    integer, intent(inout) :: VMapP(NfpTot,Ne)

    real(RP), intent(in) :: pos_en(Np,Ne,3)
    real(RP), intent(in) :: xmin, xmax
    real(RP), intent(in) :: ymin, ymax
    real(RP), intent(in) :: zmin, zmax
    integer, intent(in) :: Fmask_h(Nfp_h,Nfaces_h)
    integer, intent(in) :: Fmask_v(Nfp_v,Nfaces_v)

    integer :: ke, n
    integer :: b
    integer :: f
    integer :: i, j
    real(RP) :: x, y, z

    real(RP), parameter :: NODE_TOL = 1.0E-12_RP

    integer :: elemIDs_h(Nfp_h*Ne,Nfaces_h)
    real(RP) :: ordInfo_h(Nfp_h*Ne,Nfaces_h)
    integer :: faceIDs_h(Nfp_h*Ne,Nfaces_h)
    integer :: counterB_h(Nfaces_h)

    integer :: elemIDs_v(Nfp_v*Ne,Nfaces_v)
    real(RP) :: ordInfo_v(Nfp_v*Ne,Nfaces_v)
    integer :: faceIDs_v(Nfp_v*Ne,Nfaces_v)
    integer :: counterB_v(Nfaces_v)

    integer :: mapB_counter
    real(RP) :: rdomx, rdomy, rdomz
    !-----------------------------------------------------------------------------

    counterB_h(:) = 0
    counterB_v(:) = 0

    rdomx = 1.0_RP/(xmax - xmin)
    rdomy = 1.0_RP/(ymax - ymin)
    rdomz = 1.0_RP/abs(zmax - zmin)
    
    do ke=1, Ne
      do f=1, Nfaces_h
        x = sum(pos_en(Fmask_h(:,f),ke,1)) / dble(Nfp_h)
        y = sum(pos_en(Fmask_h(:,f),ke,2)) / dble(Nfp_h)

        call eval_domain_boundary( &
          elemIDs_h, ordInfo_h, faceIDs_h, counterB_h, & ! (inout)
          1, y, ymin, x, ke, f, rdomy                  ) ! (in)
        call eval_domain_boundary( &
          elemIDs_h, ordInfo_h, faceIDs_h, counterB_h, & ! (inout)
          2, x, xmax, y, ke, f, rdomx                  ) ! (in)    
        call eval_domain_boundary( &
          elemIDs_h, ordInfo_h, faceIDs_h, counterB_h, & ! (inout)
          3, y, ymax, x, ke, f, rdomy                  ) ! (in)
        call eval_domain_boundary( &
          elemIDs_h, ordInfo_h, faceIDs_h, counterB_h, & ! (inout)
          4, x, xmin, y, ke, f, rdomx                  ) ! (in)
      end do
      do f=1, Nfaces_v
        x = sum(pos_en(Fmask_v(:,f),ke,1)) / dble(Nfp_v)
        z = sum(pos_en(Fmask_v(:,f),ke,3)) / dble(Nfp_v)

        call eval_domain_boundary( &
          elemIDs_v, ordInfo_v, faceIDs_v, counterB_v, & ! (inout)
          1, z, zmin, x, ke, f, rdomz                  ) ! (in)
        call eval_domain_boundary( &
          elemIDs_v, ordInfo_v, faceIDs_v, counterB_v, & ! (inout)
          2, z, zmax, x, ke, f, rdomz                  ) ! (in)  
      end do    
    end do


    allocate( mapB(sum(counterB_h(:)*Nfp_h)+sum(counterB_v(:)*Nfp_v)) )
    allocate( vmapB(size(mapB)) )

    mapB_counter = 1
    do b = 1, Nfaces_h
      ! write(*,*) "LocalMesh boundary ID:", b
      ! write(*,*) counterB_h(b)
      ! write(*,*) ordInfo_h(1:counterB_h(1),b)
      ! write(*,*) elemIds_h(1:counterB_h(1),b)
      ! write(*,*) faceIds_h(1:counterB_h(1),b)

      do i=1, counterB_h(b)
        ke = elemIDs_h(i,b); f = faceIDs_h(i,b)
        do j=1, Nfp_h
          n = j + (f-1)*Nfp_h
          VMapP(n,ke) = Np*Ne + mapB_counter
          VmapB(mapB_counter) = Fmask_h(j,f) + (ke-1)*Np
          mapB_counter = mapB_counter + 1
        end do
      end do
    end do

    do b = 1, Nfaces_v
      do i=1, counterB_v(b)
        ke = elemIDs_v(i,b); f = faceIDs_v(i,b)
        do j=1, Nfp_v
          n = j + Nfp_h*Nfaces_h + Nfp_v*(f-1)
          VMapP(n,ke) = Np*Ne + mapB_counter
          VmapB(mapB_counter) = Fmask_v(j,f) + (ke-1)*Np
          mapB_counter = mapB_counter + 1
        end do
      end do
    end do

    ! write(*,*) "VMapP:-----"
    ! do b=1, 4
    !   do i=1, counterB(b)
    !     k = elemIds(i,b); f = faceIDs(i,b)
    !     write(*,*) "bid=", b, ":", mesh%VmapP(:,f,k)
    !    end do
    ! end do
    ! write(*,*) "-----"
    ! write(*,*) "VMapB:", mesh%VmapB(:)
    ! write(*,*) "NeA=", mesh%NeA

    return

  contains
!OCL SERIAL
    subroutine eval_domain_boundary( &
        elemIDs, ordInfo, faceIDs, counterB,              &
        domb_id, r, rbc, ord_info, k_, f_, normalized_fac )
      implicit none

      integer, intent(inout) :: elemIDs(:,:)
      real(RP), intent(inout) :: ordInfo(:,:)
      integer, intent(inout) :: counterB(:)
      integer, intent(inout) :: faceIDs(:,:)
      integer, intent(in) :: domb_id
      real(RP), intent(in) :: r
      real(RP), intent(in) :: rbc
      real(RP), intent(in) :: ord_info
      integer, intent(in) :: k_, f_
      real(RP), intent(in) :: normalized_fac
      !-------------------------------------------------------------

      if ( abs(r - rbc)*normalized_fac < NODE_TOL ) then
        counterB(domB_ID) = counterB(domB_ID) + 1
        ordInfo(counterB(domB_ID),domB_ID) = ord_info
        elemIDs(counterB(domB_ID),domB_ID) = k_
        faceIDs(counterB(domB_ID),domB_ID) = f_
      end if

      return
    end subroutine eval_domain_boundary
  end subroutine MeshUtil3D_genPatchBoundaryMap

!OCL SERIAL
  subroutine MeshUtil3D_buildGlobalMap( &
    panelID_table, pi_table, pj_table, pk_table,    &
    tileID_map, tileFaceID_map, tilePanelID_map,    &
    Ntile, NtileFace, NtileVertex,                  &
    isPeriodicX, isPeriodicY, isPeriodicZ,          &
    Ne_x, Ne_y, Ne_z )
    
    ! use scale_prc, only: PRC_isMaster
    implicit none

    integer, intent(in) :: Ntile
    integer, intent(in) :: NtileFace
    integer, intent(in) :: NtileVertex
    integer, intent(out) :: panelID_table(Ntile)
    integer, intent(out) :: pi_table(Ntile)
    integer, intent(out) :: pj_table(Ntile)
    integer, intent(out) :: pk_table(Ntile)
    integer, intent(out) :: tileID_map(NtileFace,Ntile)
    integer, intent(out) :: tileFaceID_map(NtileFace,Ntile)
    integer, intent(out) :: tilePanelID_map(NtileFace,Ntile)
    logical, intent(in) :: isPeriodicX
    logical, intent(in) :: isPeriodicY
    logical, intent(in) :: isPeriodicZ
    integer, intent(in) :: Ne_x
    integer, intent(in) :: Ne_y
    integer, intent(in) :: Ne_z

    integer :: NtilePerPanel
    integer :: Nv_x, Nv_y, Nv_z
    integer, allocatable :: nodesID_3d(:,:,:)
    integer, allocatable :: EToV(:,:)
    integer, allocatable :: EToE(:,:)
    integer, allocatable :: EToF(:,:)
    integer :: i, j, k, f
    integer :: tileID, tileID_R
    integer :: counter

    !-----------------------------------------------------------------------------

    NtilePerPanel = Ntile/1
    
    Nv_x = Ne_x + 1
    Nv_y = Ne_y + 1
    Nv_z = Ne_z + 1

    allocate( nodesID_3d(Nv_x,Nv_y,Nv_z) )
    allocate( EToV(Ntile,NtileVertex), EToE(Ntile,NtileFace), EToF(Ntile,NtileFace) )

    counter = 0
    do k = 1, Nv_z
    do j = 1, Nv_y
    do i = 1, Nv_x
      counter = counter + 1
      nodesID_3d(i,j,k) = counter
    end do
    end do
    end do


    !----

    tileID = 0
    do k = 1, Ne_z
    do j = 1, Ne_y
    do i = 1, Ne_x
      tileID = tileID + 1
      panelID_table(tileID) = 1
      pi_table(tileID) = i; pj_table(tileID) = j; pk_table(tileID) = k
      EToV(tileID,:) = (/ nodesID_3d(i,j  ,k  ), nodesID_3d(i+1,j  ,k  ),   &
                          nodesID_3d(i,j+1,k  ), nodesID_3d(i+1,j+1,k  ),   &
                          nodesID_3d(i,j  ,k+1), nodesID_3d(i+1,j  ,k+1),   &
                          nodesID_3d(i,j+1,k+1), nodesID_3d(i+1,j+1,k+1) /)
    end do
    end do
    end do

    call MeshUtil3D_genConnectivity( EToE, EToF, &
      EToV, Ntile, NtileFace )
    tileID_map(:,:) = transpose(EToE)
    tileFaceID_map(:,:) = transpose(EToF)

    do tileID=1, Ntile
    do f=1, NtileFace
      tileID_R = tileID_map(f,tileID)
      tilePanelID_map(f,tileID) = panelID_table(tileID_R)
    end do
    end do

    if (isPeriodicX) then
      do tileID=1, Ntile
        if (pi_table(tileID) == 1 .and. tileFaceID_map(4,tileID) == 4) then
          tileID_map(4,tileID) = Ne_x + (pj_table(tileID) - 1)*Ne_x + (pk_table(tileID) - 1)*Ne_x*Ne_y
          tileFaceID_map(4,tileID) = 2
        end if
        if (pi_table(tileID) == Ne_x .and. tileFaceID_map(2,tileID) == 2) then
          tileID_map(2,tileID) = 1 + (pj_table(tileID) - 1)*Ne_x + (pk_table(tileID) - 1)*Ne_x*Ne_y
          tileFaceID_map(2,tileID) = 4
        end if
      end do
    end if

    if (isPeriodicY) then
      do tileID=1, Ntile
        if (pj_table(tileID) == 1 .and. tileFaceID_map(1,tileID) == 1) then
          tileID_map(1,tileID) = pi_table(tileID) + (Ne_y - 1)*Ne_x + (pk_table(tileID) - 1)*Ne_x*Ne_y
          tileFaceID_map(1,tileID) = 3
        end if
        if (pj_table(tileID) == Ne_y .and. tileFaceID_map(3,tileID) == 3) then
          tileID_map(3,tileID) = pi_table(tileID) + (pk_table(tileID) - 1)*Ne_x*Ne_y
          tileFaceID_map(3,tileID) = 1
        end if
      end do
    end if

    if (isPeriodicZ) then
      do tileID=1, Ntile
        if (pk_table(tileID) == 1 .and. tileFaceID_map(5,tileID) == 5) then
          tileID_map(5,tileID) = pi_table(tileID) + (pj_table(tileID) - 1)*Ne_x  + (Ne_z - 1)*Ne_x*Ne_y
          tileFaceID_map(5,tileID) = 6
        end if
        if (pk_table(tileID) == Ne_z .and. tileFaceID_map(6,tileID) == 6) then
          tileID_map(6,tileID) = pi_table(tileID) + (pj_table(tileID) - 1)*Ne_x
          tileFaceID_map(6,tileID) = 5
        end if
      end do
    end if

    return
    
    !--
    ! if (PRC_isMaster) then
    !   write(*,*) "TotTile", Ntile
    !   do tileID = 1, Ntile
    !     write(*,*) "tileID:", tileID, ", EtoV:", EtoV(tileID,:)
    !   end do
    !   write(*,*) "-----------"
    !   do tileID = 1, Ntile
    !     write(*,*) "tileID:", tileID, ", EtoE:", EtoE(tileID,:)
    !   end do
    !   write(*,*) "-----------"
    !   do tileID = 1, Ntile
    !     write(*,*) "tileID:", tileID, ", EtoF:", EtoF(tileID,:)
    !   end do
    ! end if
    
  end subroutine MeshUtil3D_buildGlobalMap

end module scale_meshutil_3d
