!-------------------------------------------------------------------------------
!> module FElib / Mesh / utility for 2D mesh
!!
!! @par Description
!!          A module useful for generating 2D mesh 
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_meshutil_2d
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
  public :: MeshUtil2D_genRectDomain
  public :: MeshUtil2D_genConnectivity
  public :: MeshUtil2D_buildInteriorMap
  public :: MeshUtil2D_genPatchBoundaryMap
  public :: MeshUtil2D_buildGlobalMap

contains
!OCL SERIAL
  subroutine MeshUtil2D_genRectDomain( pos_v, EToV, &
    Ke_x, xmin, xmax, Ke_y, ymin, ymax )

    implicit none

    integer, intent(in) :: Ke_x
    real(RP), intent(in) :: xmin, xmax
    integer, intent(in) :: Ke_y
    real(RP), intent(in) :: ymin, ymax

    real(RP), intent(out) :: pos_v((Ke_x+1)*(Ke_y+1),2)
    integer, intent(out) :: EToV(Ke_x*Ke_y,4)

    integer :: i, j
    integer :: n
    integer :: NvX, NvY
    !-----------------------------------------------------------------------------

    NvX = Ke_X + 1
    NvY = Ke_Y + 1

    do j=1, NvY
    do i=1, NvX
        n = i + (j-1)*NvX
        pos_v(n,1) = (xmax - xmin)*dble(i - 1)/dble(NvX - 1) + xmin
        pos_v(n,2) = (ymax - ymin)*dble(j - 1)/dble(NvY - 1) + ymin
    end do
    end do

    do j=1, Ke_y
    do i=1, Ke_x
        n = i + (j-1)*Ke_x
        EToV(n,1:4) = (j-1)*NvX + &
            & (/ i, i+1,              &
            &    i + NvX, i+1 + NvX /)
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
  end subroutine MeshUtil2D_genRectDomain

!OCL SERIAL
  subroutine MeshUtil2D_genConnectivity( EToE, EToF, &
    EToV, Ne, Nfaces )
    
    use scale_quicksort, only: QUICKSORT_exec_with_idx
    implicit none

    integer, intent(in) :: Ne
    integer, intent(in) :: Nfaces
    integer, intent(out) :: EToE(Ne, Nfaces)
    integer, intent(out) :: EToF(Ne, Nfaces)
    integer, intent(in) :: EToV(Ne,4)

    integer :: nodes(Ne*Nfaces,2)
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
    !-----------------------------------------------------------------------------

    Nnodes = maxval( EToV )
    Nnodes_row = size(nodes,1)

    !---------
    do ke=1, Ne
       nodes(ke     ,:) = EToV(ke,(/ 1, 2 /))
       nodes(ke+Ne  ,:) = EToV(ke,(/ 2, 4 /))
       nodes(ke+2*Ne,:) = EToV(ke,(/ 4, 3 /))
       nodes(ke+3*Ne,:) = EToV(ke,(/ 3, 1 /))
    end do
   
    ! Sort
    do n=1, Nnodes_row
       if (nodes(n,1) > nodes(n,2) ) then
          tmp = nodes(n,1);
          nodes(n,1) = nodes(n,2)
          nodes(n,2) = tmp
       end if
       ! write(*,*) n, ":", nodes(n,:)
    end do
    nodes = nodes - 1
    !---------

    do ke=1, Ne
       EToE(ke,:) = ke
       EToF(ke,:) = (/ 1, 2, 3, 4 /)
    end do

    face_ids(:) = nodes(:,1)*Nnodes + nodes(:,2) + 1

    do f=1, Nfaces
    do ke=1, Ne
       n = ke + (f-1)*Ne
       spNodeToNode(:,n) = (/ n, EToE(ke,f), EToF(ke,f) /)
       sorted_faceid(n) = face_ids(n)
       sort_indx(n) = n
       ! write(*,*) "face_id, n, EToE, EToF:", spNodeToNode(:,n)
    end do
    end do

    !- sort row
    call QUICKSORT_exec_with_idx( Ne*Nfaces, sorted_faceid, sort_indx )
    spNodeToNodeRowTmp(:,:) = spNodeToNode(:,:)
    do n=1, Nnodes_row
      spNodeToNode(:,n) = spNodeToNodeRowTmp(:,sort_indx(n))
      ! write(*,'(a,4i10)') "(sorted) face_id, n, EToE, EToF:", spNodeToNode(:,n)
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
  end subroutine MeshUtil2D_genConnectivity

!OCL SERIAL
  subroutine MeshUtil2D_BuildInteriorMap( VMapM, VMapP, MapM, MapP, &
    pos_en, pos_ev, EtoE, EtoF, EtoV, Fmask, Ne, Np, Nfp, Nfaces, Nv)

    implicit none

    integer, intent(in) :: Ne
    integer, intent(in) :: Np
    integer, intent(in) :: Nfp
    integer, intent(in) :: Nfaces
    integer, intent(in) :: Nv

    integer, intent(out) :: VMapM(Nfp,Nfaces,Ne)
    integer, intent(out) :: VMapP(Nfp,Nfaces,Ne)
    integer, intent(out) :: MapM(Nfp,Nfaces,Ne)
    integer, intent(out) :: MapP(Nfp,Nfaces,Ne)

    real(RP), intent(in) :: pos_en(Np,Ne,2)
    real(RP), intent(in) :: pos_ev(Nv,2)
    integer, intent(in) :: EToE(Ne, Nfaces)
    integer, intent(in) :: EToF(Ne, Nfaces)
    integer, intent(in) :: EToV(Ne, 4)
    integer, intent(in) :: Fmask(Nfp,4)


    integer :: ke, ke1, ke2
    integer :: f, f1, f2
    integer :: p
    integer :: n
    integer :: idP, idM
    integer :: v1, v2
    integer :: nodeids(Np,Ne)
    real(RP) :: x1(Nfp,Nfp), x2(Nfp,Nfp)
    real(RP) :: y1(Nfp,Nfp), y2(Nfp,Nfp)
    real(RP) :: dist(Nfp,Nfp)
    real(RP) :: x(Np*Ne), y(Np*Ne)

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
    end do
    end do
    !$omp end do
    !$omp do collapse(2)
    do ke=1, Ne
    do f=1, Nfaces
    do p=1, Nfp
      n = p + (f-1)*Nfp + (ke-1)*Nfp*Nfaces
      MapM(p,f,ke) = n
      MapP(p,f,ke) = n
      VMapM(p,f,ke) = nodeids(Fmask(p,f),ke)
    end do
    end do
    end do

    !$omp workshare
    VMapP(:,:,:) = -1
    !$omp end workshare
    !$omp end parallel

    !$omp parallel private( &
    !$omp ke1, f1, ke2, f2, v1, v2,                     &
    !$omp x1, x2, y1, y2, dist, mindist_indx, idP, idM  )
    
    !$omp do
    do ke1=1, Ne
    do f1=1, Nfaces
      ke2 = EToE(ke1,f1); f2 = EToF(ke1,f1)

      v1 = EToV(ke1,f1); v2 = EToV(ke1,1+mod(f1,Nfaces))

      x1(:,:) = spread( x(VMapM(:,f1,ke1)), 2, Nfp )
      x2(:,:) = spread( x(VMapM(:,f2,ke2)), 1, Nfp )
      y1(:,:) = spread( y(VMapM(:,f1,ke1)), 2, Nfp )
      y2(:,:) = spread( y(VMapM(:,f2,ke2)), 1, Nfp )

      dist(:,:) = (x1(:,:) - x2(:,:))**2 &
                + (y1(:,:) - y2(:,:))**2
      do idM=1, Nfp
        mindist_indx(:) = minloc(dist(idM,:))
        idP = mindist_indx(1)
        VMapP(idM,f1,ke1) = VMapM(idP,f2,ke2)
        MapP(idM,f1,ke1) = idP + (f2-1)*Nfp + (ke2-1)*Nfp*Nfaces
      end do
    end do
    end do
    !$omp end do
    !$omp end parallel

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

  end subroutine MeshUtil2D_BuildInteriorMap

!OCL SERIAL
  subroutine MeshUtil2D_genPatchBoundaryMap(  VMapB, MapB, VMapP, &
    pos_en, xmin, xmax, ymin, ymax, Fmask, Ne, Np, Nfp, Nfaces, Nv)

    implicit none

    integer, intent(in) :: Ne
    integer, intent(in) :: Np
    integer, intent(in) :: Nfp
    integer, intent(in) :: Nfaces
    integer, intent(in) :: Nv

    integer, intent(inout), allocatable :: VMapB(:)
    integer, intent(inout), allocatable :: MapB(:)
    integer, intent(inout) :: VMapP(Nfp,Nfaces,Ne)

    real(RP), intent(in) :: pos_en(Np,Ne,2)
    real(RP), intent(in) :: xmin, xmax
    real(RP), intent(in) :: ymin, ymax
    integer, intent(in) :: Fmask(Nfp,4)


    integer :: ke
    integer :: b
    integer :: f
    integer :: i, j
    real(RP) :: x, y

    real(RP), parameter :: NODE_TOL = 1.0E-12_RP

    real(RP) :: ordInfo(Nfp*Ne,4)
    integer :: elemIDs(Nfp*Ne,4)
    integer :: faceIDs(Nfp*Ne,4)
    integer :: counterB(4)
    integer :: mapB_counter
    real(RP) :: rdomx, rdomy
    !-----------------------------------------------------------------------------


    counterB(:) = 0
    rdomx = 1.0_RP/(xmax - xmin)
    rdomy = 1.0_RP/(ymax - ymin)

    do ke=1, Ne
    do f=1, Nfaces
      x = sum(pos_en(Fmask(:,f),ke,1)/dble(Nfp))
      y = sum(pos_en(Fmask(:,f),ke,2)/dble(Nfp))

      call eval_domain_boundary(1, y, ymin, x, ke, f, rdomy)
      call eval_domain_boundary(2, x, xmax, y, ke, f, rdomx)
      call eval_domain_boundary(3, y, ymax, x, ke, f, rdomy)
      call eval_domain_boundary(4, x, xmin, y, ke, f, rdomx)
    end do
    end do

    allocate( mapB(sum(counterB*Nfp))    )
    allocate( vmapB(size(mapB)) )

    mapB_counter = 1
    do b = 1, 4
    !  write(*,*) "LocalMesh boundary ID:", b
    !  write(*,*) counterB(b)
    !  write(*,*) ordInfo(1:counterB(b),b)
    !  write(*,*) elemIds(1:counterB(b),b)
    !  write(*,*) faceIds(1:counterB(b),b)

      do i=1, counterB(b)
        ke = elemIds(i,b); f = faceIDs(i,b)
        do j=1, Nfp
          VMapP(j,f,ke) = Np*Ne + mapB_counter
          VmapB(mapB_counter) = Fmask(j,f) + (ke-1)*Np
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
     subroutine eval_domain_boundary(domb_id, r, rbc, ord_info, ke_, f_, normalized_fac)
        integer, intent(in) :: domb_id
        real(RP), intent(in) :: r
        real(RP), intent(in) :: rbc
        real(RP), intent(in) :: ord_info
        integer, intent(in) :: ke_, f_
        real(RP), intent(in) :: normalized_fac
        !-------------------------------------------------------------

        if ( abs(r - rbc)*normalized_fac < NODE_TOL ) then
          counterB(domB_ID) = counterB(domB_ID) + 1
          ordInfo(counterB(domB_ID),domB_ID) = ord_info
          elemIds(counterB(domB_ID),domB_ID) = ke_
          faceIds(counterB(domB_ID),domB_ID) = f_
        end if

        return
     end subroutine eval_domain_boundary
     
  end subroutine MeshUtil2D_genPatchBoundaryMap

!OCL SERIAL
  subroutine MeshUtil2D_buildGlobalMap( &
    panelID_table, pi_table, pj_table,           &
    tileID_map, tileFaceID_map, tilePanelID_map, &
    Ntile, isPeriodicX, isPeriodicY,             &
    Ne_x, Ne_y )

    implicit none

    integer, intent(in) :: Ntile
    integer, intent(out) :: panelID_table(Ntile)
    integer, intent(out) :: pi_table(Ntile)
    integer, intent(out) :: pj_table(Ntile)
    integer, intent(out) :: tileID_map(4,Ntile)
    integer, intent(out) :: tileFaceID_map(4,Ntile)
    integer, intent(out) :: tilePanelID_map(4,Ntile)
    logical, intent(in) :: isPeriodicX
    logical, intent(in) :: isPeriodicY
    integer, intent(in) :: Ne_x
    integer, intent(in) :: Ne_y

    integer :: NtilePerPanel
    integer :: NvX, NvY
    integer, allocatable :: nodesID_2d(:,:)
    integer, allocatable :: EToV(:,:)
    integer, allocatable :: EToE(:,:)
    integer, allocatable :: EToF(:,:)
    integer :: i, j, f
    integer :: panelID
    integer :: tileID, tileID_R
    integer :: counter

    !-----------------------------------------------------------------------------

    NtilePerPanel = Ntile/1
    NvX = Ne_x + 1
    NvY = Ne_y + 1
    allocate( nodesID_2d(NvX, NvY) )
    allocate( EToV(Ntile,4), EToE(Ntile,4), EToF(Ntile,4) )

    counter = 0
    do j = 1, NvY
    do i = 1, NvX
      counter = counter + 1
      nodesID_2d(i,j) = counter
    end do
    end do

    !----

    tileID = 0
    do j = 1, Ne_y
    do i = 1, Ne_x
      tileID = tileID + 1
      panelID_table(tileID) = 1
      pi_table(tileID) = i; pj_table(tileID) = j
      EToV(tileID,:) = (/ nodesID_2d(i,j  ), nodesID_2d(i+1,j  ),    &
        &                 nodesID_2d(i,j+1), nodesID_2d(i+1,j+1)/)
    end do
    end do

    call MeshUtil2D_genConnectivity( EToE, EToF, &
      & EToV, Ntile, 4 )
    tileID_map(:,:) = transpose(EToE)
    tileFaceID_map(:,:) = transpose(EToF)

    do tileID=1, Ntile
    do f=1, 4
      tileID_R = tileID_map(f,tileID)
      tilePanelID_map(f,tileID) = panelID_table(tileID_R)
    end do
    end do

    if (isPeriodicX) then
      do tileID=1, Ntile
        if (pi_table(tileID) == 1 .and. tileFaceID_map(4,tileID) == 4) then
          tileID_map(4,tileID) = Ne_x + (pj_table(tileID) - 1)*Ne_x
          tileFaceID_map(4,tileID) = 2
        end if
        if (pi_table(tileID) == Ne_x .and. tileFaceID_map(2,tileID) == 2) then
          tileID_map(2,tileID) = 1 + (pj_table(tileID) - 1)*Ne_x
          tileFaceID_map(2,tileID) = 4
        end if
      end do
    end if

    if (isPeriodicY) then
      do tileID=1, Ntile
        if (pj_table(tileID) == 1 .and. tileFaceID_map(1,tileID) == 1) then
          tileID_map(1,tileID) = pi_table(tileID) + (Ne_y - 1)*Ne_x
          tileFaceID_map(1,tileID) = 3
        end if
        if (pj_table(tileID) == Ne_y .and. tileFaceID_map(3,tileID) == 3) then
          tileID_map(3,tileID) = pi_table(tileID)
          tileFaceID_map(3,tileID) = 1
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
  end subroutine MeshUtil2D_buildGlobalMap

end module scale_meshutil_2d
