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

  subroutine MeshUtil3D_genCubeDomain( pos_v, EToV,        &
    Ke_x, xmin, xmax, Ke_y, ymin, ymax, Ke_z, zmin, zmax )

    implicit none

    integer, intent(in) :: Ke_x
    real(RP), intent(in) :: xmin, xmax
    integer, intent(in) :: Ke_y
    real(RP), intent(in) :: ymin, ymax
    integer, intent(in) :: Ke_z
    real(RP), intent(in) :: zmin, zmax
    
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

    do k=1, NvZ
    do j=1, NvY
    do i=1, NvX
      n = i + (j-1)*NvX + (k-1)*NvX*NvY
      pos_v(n,1) = (xmax - xmin)*dble(i - 1)/dble(NvX - 1) + xmin
      pos_v(n,2) = (ymax - ymin)*dble(j - 1)/dble(NvY - 1) + ymin
      pos_v(n,3) = (zmax - zmin)*dble(k - 1)/dble(NvZ - 1) + zmin
    end do
    end do
    end do

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
    integer :: face_ids(Ne*Nfaces)
    integer :: k
    integer :: f
    integer :: n
    integer :: n1, n2
    integer :: Nnodes
    integer :: tmp
    integer :: Nnodes_row
    integer :: matchL(2,4), matchR(2,4)

    real(RP) :: EToE_1d(Ne*Nfaces)
    real(RP) :: EToF_1d(Ne*Nfaces)

    integer :: spNodeToNode(4,Ne*Nfaces)
    integer :: spNodeToNodeRowTmp(4,Ne*Nfaces)
    integer :: sort_indx(Ne*Nfaces)
    integer :: sort_val(Ne*Nfaces)

    real(RP) :: vtmp(4)

    !-----------------------------------------------------------------------------

    Nnodes = maxval( EToV )
    Nnodes_row = size(nodes,1)

    !---------
    do n=1, Ne
       nodes(n     ,:) = EToV(n,(/ 1, 2, 6, 5 /))
       nodes(n+Ne  ,:) = EToV(n,(/ 2, 4, 8, 6 /))
       nodes(n+2*Ne,:) = EToV(n,(/ 4, 3, 7, 8 /))
       nodes(n+3*Ne,:) = EToV(n,(/ 3, 1, 5, 7 /))
       nodes(n+4*Ne,:) = EToV(n,(/ 1, 3, 4, 2 /))       
       nodes(n+5*Ne,:) = EToV(n,(/ 5, 6, 8, 7 /))       
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

    do n=1, Ne
       EToE(n,:) = n
       EToF(n,:) = (/ 1, 2, 3, 4, 5, 6 /)
    end do

    face_ids(:) =   nodes(:,1)*Nnodes**3 + nodes(:,2)*Nnodes**2 &
                  + nodes(:,3)*Nnodes + nodes(:,4) + 1

    do f=1, Nfaces
    do k=1, Ne
       n = k + (f-1)*Ne
       spNodeToNode(:,n) = (/ face_ids(n), n, EToE(k,f), EToF(k,f) /)
       sort_val(n) = face_ids(n)
       sort_indx(n) = n
       ! write(*,*) "face_id, n, EToE, EToF:", spNodeToNode(:,n)
    end do
    end do

    
    !- sort row
    call QUICKSORT_exec_with_idx( Ne*Nfaces, sort_val, sort_indx )
    spNodeToNodeRowTmp(:,:) = spNodeToNode(:,:)
    do n=1, Nnodes_row
      spNodeToNode(:,n) = spNodeToNodeRowTmp(:,sort_indx(n))
      ! write(*,'(a,4i10)') "(sorted) face_id, n, EToE, EToF:", spNodeToNode(:,n)
    end do

    EToE_1d(:) = -1
    EToF_1d(:) = -1
    do n=1, Nnodes_row-1
       if ( spNodeToNode(1,n) - spNodeToNode(1,n+1) == 0 ) then
          matchL(:,:) = transpose( spNodeToNode(:,(/ n, n+1 /)) )
          matchR(:,:) = transpose( spNodeToNode(:,(/ n+1, n /)) )

          EToE_1d(matchL(:,2)) = matchR(:,3)
          EToF_1d(matchL(:,2)) = matchR(:,4)
       end if
    end do

    do f=1, Nfaces
    do k=1, Ne
       n = k + (f-1)*Ne
       if ( EToE_1d(n) /= -1 ) then
          EToE(k,f) = EToE_1d(n)
          EToF(k,f) = EToF_1d(n)
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

    integer :: k, k1, k2
    integer :: f, f1, f2
    integer :: p
    integer :: n
    integer :: i
    integer :: idP, idM
    integer :: v1, v2
    integer :: nodeids(Np,Ne)
    real(RP) :: refd2
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
    !-----------------------------------------------------------------------------

    do k=1, Ne
    do p=1, Np
      n = p + (k-1)*Np
      nodeids(p,k) = n
      x(n) = pos_en(p,k,1)
      y(n) = pos_en(p,k,2)
      z(n) = pos_en(p,k,3)
    end do
    end do

    do k=1, Ne
      do f=1, Nfaces_h
      do p=1, Nfp_h
        n = p + (f-1)*Nfp_h + (k-1)*NfpTot
        MapM_h(p,f,k) = n
        MapP_h(p,f,k) = n
        VMapM_h(p,f,k) = nodeids(Fmask_h(p,f),k)
      end do
      end do
      do f=1, Nfaces_v
      do p=1, Nfp_v
        n = p + Nfaces_h*Nfp_h + (f-1)*Nfp_v + (k-1)*NfpTot
        MapM_v(p,f,k) = n
        MapP_v(p,f,k) = n
        VMapM_v(p,f,k) = nodeids(Fmask_v(p,f),k)
      end do
      end do    
    end do

    VMapP_h(:,:,:) = -1
    do k1=1, Ne
    do f1=1, Nfaces_h
      k2 = EToE(k1,f1); f2 = EToF(k1,f1)

      v1 = EToV(k1,f1); v2 = EToV(k1,1+mod(f1,Nfaces_h))

      refd2 =    (pos_ev(v1,1) - pos_ev(v2,1))**2   &
               + (pos_ev(v1,2) - pos_ev(v2,2))**2   &
               + (pos_ev(v1,3) - pos_ev(v2,3))**2

      r_h(:,:,1,1) = spread( x(VMapM_h(:,f1,k1)), 2, Nfp_h )
      r_h(:,:,1,2) = spread( x(VMapM_h(:,f2,k2)), 1, Nfp_h )
      r_h(:,:,2,1) = spread( y(VMapM_h(:,f1,k1)), 2, Nfp_h )
      r_h(:,:,2,2) = spread( y(VMapM_h(:,f2,k2)), 1, Nfp_h )
      r_h(:,:,3,1) = spread( z(VMapM_h(:,f1,k1)), 2, Nfp_h )
      r_h(:,:,3,2) = spread( z(VMapM_h(:,f2,k2)), 1, Nfp_h )
      
      dist_h(:,:) =  (r_h(:,:,1,1) - r_h(:,:,1,2))**2  &
                    + (r_h(:,:,2,1) - r_h(:,:,2,2))**2  &
                    + (r_h(:,:,3,1) - r_h(:,:,3,2))**2
      do idP=1, Nfp_h
      do idM=1, Nfp_h
          if (dist_h(idM,idP)/refd2 < 1d-14) then
            VMapP_h(idM,f1,k1) = VMapM_h(idP,f2,k2)
            MapP_h(idM,f1,k1) = idP + (f2-1)*Nfp_h + (k2-1)*NfpTot
          end if
      end do
      end do
    end do
    end do

    VMapP_v(:,:,:) = -1
    do k1=1, Ne
    do f1=1, Nfaces_v
      k2 = EToE(k1,Nfaces_h+f1); f2 = EToF(k1,Nfaces_h+f1) - Nfaces_h
      
      v1 = EToV(k1,1); v2 = EToV(k1,Nfaces_h+1)
      refd2 =  sum( (pos_ev(v1,:) - pos_ev(v2,:))**2 )

      r_v(:,:,1,1) = spread( x(VMapM_v(:,f1,k1)), 2, Nfp_v )
      r_v(:,:,1,2) = spread( x(VMapM_v(:,f2,k2)), 1, Nfp_v )
      r_v(:,:,2,1) = spread( y(VMapM_v(:,f1,k1)), 2, Nfp_v )
      r_v(:,:,2,2) = spread( y(VMapM_v(:,f2,k2)), 1, Nfp_v )
      r_v(:,:,3,1) = spread( z(VMapM_v(:,f1,k1)), 2, Nfp_v )
      r_v(:,:,3,2) = spread( z(VMapM_v(:,f2,k2)), 1, Nfp_v )
      
      dist_v(:,:) =   (r_v(:,:,1,1) - r_v(:,:,1,2))**2  &
                    + (r_v(:,:,2,1) - r_v(:,:,2,2))**2  &
                    + (r_v(:,:,3,1) - r_v(:,:,3,2))**2
      do idP=1, Nfp_v
      do idM=1, Nfp_v
          if (dist_v(idM,idP)/refd2 < 1d-14) then
            VMapP_v(idM,f1,k1) = VMapM_v(idP,f2,k2)
            MapP_v(idM,f1,k1) = idP + Nfaces_h*Nfp_h + (f2-1)*Nfp_v + (k2-1)*NfpTot
          end if
      end do
      end do
    end do
    end do

    do k=1, Ne
      do f=1, Nfaces_h
      do n=1, Nfp_h
        i = n + (f-1)*Nfp_h
        VMapM(i,k) = VMapM_h(n,f,k)
        MapM(i,k) = MapM_h(n,f,k)
        VMapP(i,k) = VMapP_h(n,f,k)
        MapP(i,k) = MapP_h(n,f,k)
      end do
      end do
      do f=1, Nfaces_v
      do n=1, Nfp_v
        i = n + Nfaces_h*Nfp_h + (f-1)*Nfp_v
        VMapM(i,k) = VMapM_v(n,f,k)
        MapM(i,k) = MapM_v(n,f,k)
        VMapP(i,k) = VMapP_v(n,f,k)
        MapP(i,k) = MapP_v(n,f,k)
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

    integer :: k, n
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
    rdomz = 1.0_RP/(zmax - zmin)

    do k=1, Ne
      do f=1, Nfaces_h
        x = sum(pos_en(Fmask_h(:,f),k,1)/dble(Nfp_h))
        y = sum(pos_en(Fmask_h(:,f),k,2)/dble(Nfp_h))

        call eval_domain_boundary( &
          elemIDs_h, ordInfo_h, faceIDs_h, counterB_h, & ! (inout)
          1, y, ymin, x, k, f, rdomy                   ) ! (in)
        call eval_domain_boundary( &
          elemIDs_h, ordInfo_h, faceIDs_h, counterB_h, & ! (inout)
          2, x, xmax, y, k, f, rdomx                   ) ! (in)    
        call eval_domain_boundary( &
          elemIDs_h, ordInfo_h, faceIDs_h, counterB_h, & ! (inout)
          3, y, ymax, x, k, f, rdomy                   ) ! (in)
        call eval_domain_boundary( &
          elemIDs_h, ordInfo_h, faceIDs_h, counterB_h, & ! (inout)
          4, x, xmin, y, k, f, rdomx                   ) ! (in)
      end do
      do f=1, Nfaces_v
        x = sum(pos_en(Fmask_v(:,f),k,1)/dble(Nfp_v))
        z = sum(pos_en(Fmask_v(:,f),k,3)/dble(Nfp_v))

        call eval_domain_boundary( &
          elemIDs_v, ordInfo_v, faceIDs_v, counterB_v, & ! (inout)
          1, z, zmin, x, k, f, rdomz                   ) ! (in)
        call eval_domain_boundary( &
          elemIDs_v, ordInfo_v, faceIDs_v, counterB_v, & ! (inout)
          2, z, zmax, x, k, f, rdomz                   ) ! (in)  
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
        k = elemIDs_h(i,b); f = faceIDs_h(i,b)
        do j=1, Nfp_h
          n = j + (f-1)*Nfp_h
          VMapP(n,k) = Np*Ne + mapB_counter
          VmapB(mapB_counter) = Fmask_h(j,f) + (k-1)*Np
          mapB_counter = mapB_counter + 1
        end do
      end do
    end do

    do b = 1, Nfaces_v
      do i=1, counterB_v(b)
        k = elemIDs_v(i,b); f = faceIDs_v(i,b)
        do j=1, Nfp_v
          n = j + Nfp_h*Nfaces_h + Nfp_v*(f-1)
          VMapP(n,k) = Np*Ne + mapB_counter
          VmapB(mapB_counter) = Fmask_v(j,f) + (k-1)*Np
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

  subroutine MeshUtil3D_buildGlobalMap( &
    panelID_table, pi_table, pj_table, pk_table,    &
    tileID_map, tileFaceID_map, tilePanelID_map,    &
    Ntile, NtileFace,                               &
    isPeriodicX, isPeriodicY, isPeriodicZ )
    
    use scale_prc, only: PRC_isMaster
    implicit none

    integer, intent(in) :: Ntile
    integer, intent(in) :: NtileFace
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

    integer :: NtilePerPanel
    integer :: Ne_h, Ne_v, Nv_h, Nv_v
    integer, allocatable :: nodesID_3d(:,:,:)
    integer, allocatable :: EToV(:,:)
    integer, allocatable :: EToE(:,:)
    integer, allocatable :: EToF(:,:)
    integer :: i, j, k, f
    integer :: panelID
    integer :: tileID, tileID_R
    integer :: counter

    !-----------------------------------------------------------------------------

    NtilePerPanel = Ntile/1
    
    Ne_h = sqrt(dble(NtilePerPanel))
    Nv_h = Ne_h + 1
    
    Ne_v = 1
    Nv_v = Ne_v + 1

    allocate( nodesID_3d(Nv_h,Nv_h,Nv_v) )
    allocate( EToV(Ntile,NtileFace), EToE(Ntile,NtileFace), EToF(Ntile,NtileFace) )

    counter = 0
    do k = 1, Nv_v
    do j = 1, Nv_h
    do i = 1, Nv_h
      counter = counter + 1
      nodesID_3d(i,j,k) = counter
    end do
    end do
    end do


    !----

    tileID = 0
    do k = 1, Ne_v
    do j = 1, Ne_h
    do i = 1, Ne_h
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
      & EToV, Ntile, NtileFace )
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
          tileID_map(4,tileID) = Ne_h + (pj_table(tileID) - 1)*Ne_h + (pk_table(tileID) - 1)*Ne_h**2
          tileFaceID_map(4,tileID) = 2
        end if
        if (pi_table(tileID) == Ne_h .and. tileFaceID_map(2,tileID) == 2) then
          tileID_map(2,tileID) = 1 + (pj_table(tileID) - 1)*Ne_h + (pk_table(tileID) - 1)*Ne_h**2
          tileFaceID_map(2,tileID) = 4
        end if
      end do
    end if

    if (isPeriodicY) then
      do tileID=1, Ntile
        if (pj_table(tileID) == 1 .and. tileFaceID_map(1,tileID) == 1) then
          tileID_map(1,tileID) = pi_table(tileID) + (Ne_h - 1)*Ne_h + (pk_table(tileID) - 1)*Ne_h**2
          tileFaceID_map(1,tileID) = 3
        end if
        if (pj_table(tileID) == Ne_h .and. tileFaceID_map(3,tileID) == 3) then
          tileID_map(3,tileID) = pi_table(tileID) + (pk_table(tileID) - 1)*Ne_h**2
          tileFaceID_map(3,tileID) = 1
        end if
      end do
    end if

    if (isPeriodicZ) then
      do tileID=1, Ntile
        if (pk_table(tileID) == 1 .and. tileFaceID_map(5,tileID) == 5) then
          tileID_map(5,tileID) = pi_table(tileID) + (pj_table(tileID) - 1)*Ne_h  + (Ne_v - 1)*Ne_h**2
          tileFaceID_map(5,tileID) = 6
        end if
        if (pk_table(tileID) == Ne_v .and. tileFaceID_map(6,tileID) == 6) then
          tileID_map(6,tileID) = pi_table(tileID) + (pj_table(tileID) - 1)*Ne_h
          tileFaceID_map(6,tileID) = 5
        end if
      end do
    end if

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

    return
  end subroutine MeshUtil3D_buildGlobalMap

end module scale_meshutil_3d
