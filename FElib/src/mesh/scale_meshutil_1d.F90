#include "scaleFElib.h"
module scale_meshutil_1d
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
  public :: MeshUtil1D_genLineDomain
  public :: MeshUtil1D_genConnectivity
  public :: MeshUtil1D_buildInteriorMap
  public :: MeshUtil1D_genPatchBoundaryMap
  public :: MeshUtil1D_genPeriodicBoundaryMap
  public :: MeshUtil1D_buildGlobalMap

contains

  subroutine MeshUtil1D_genLineDomain( pos_v, EToV, &
      Ke_x, xmin, xmax )

    implicit none

    integer, intent(in) :: Ke_x
    real(RP), intent(in) :: xmin, xmax

    real(RP), intent(out) :: pos_v(Ke_x+1,1)
    integer, intent(out) :: EToV(Ke_x,2)

    integer :: i
    integer :: n
    integer :: NvX

    !-----------------------------------------------------------------------------

    NvX = Ke_X + 1

    do i=1, NvX
      pos_v(i,1) = (xmax - xmin)*(i - 1)/dble(NvX - 1) + xmin
    end do

    do i=1, Ke_x
      EToV(i,1:2) = (/ i, i+1 /)
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

  end subroutine MeshUtil1D_genLineDomain


  subroutine MeshUtil1D_genConnectivity( EToE, EToF, &
    EToV, Ne, Nfaces )
    
    implicit none

    integer, intent(in) :: Ne
    integer, intent(in) :: Nfaces
    integer, intent(out) :: EToE(Ne, Nfaces)
    integer, intent(out) :: EToF(Ne, Nfaces)
    integer, intent(in) :: EToV(Ne,2)

    integer :: nodes(Ne*Nfaces,1)
    integer :: face_ids(Ne*Nfaces)
    integer :: spNodeToNode(Ne*Nfaces,4)
    integer :: k
    integer :: f
    integer :: n
    integer :: n1, n2
    integer :: Nnodes
    integer :: tmp
    integer :: spNodeToNodeRowTmp(4)
    integer :: Nnodes_row
    integer :: matchL(2,4), matchR(2,4)

    real(RP) :: EToE_1d(Ne*Nfaces)
    real(RP) :: EToF_1d(Ne*Nfaces)

    !-----------------------------------------------------------------------------

    Nnodes = maxval( EToV )
    Nnodes_row = size(nodes,1)

    !---------
    do n=1, Ne
       nodes(n     ,:) = EToV(n,(/ 1 /))
       nodes(n+Ne  ,:) = EToV(n,(/ 2 /))
    end do
    nodes = nodes - 1

    !---------

    do n=1, Ne
       EToE(n,:) = n
       EToF(n,:) = (/ 1, 2 /)
    end do

    face_ids(:) = nodes(:,1) + 1

    do f=1, Nfaces
    do k=1, Ne
       n = k + (f-1)*Ne
       spNodeToNode(n,:) = (/ face_ids(n), n, EToE(k,f), EToF(k,f) /)
       ! write(*,*) "face_id, n, EToE, EToF:", spNodeToNode(n,:)
    end do
    end do

    ! Sort row
    do n1=1, Nnodes_row-1
    do n2=n1+1, Nnodes_row
       if (spNodeToNode(n1,1) > spNodeToNode(n2,1)) then
          spNodeToNodeRowTmp(:) = spNodeToNode(n1,:)
          spNodeToNode(n1,:) = spNodeToNode(n2,:)
          spNodeToNode(n2,:) = spNodeToNodeRowTmp
       end if
    end do
    end do

    ! do n=1, Nnodes_row
    !    write(*,*) "(sorted) face_id, n, EToE, EToF:", spNodeToNode(n,:)
    ! end do

    !
    EToE_1d(:) = -1
    EToF_1d(:) = -1
    do n=1, Nnodes_row-1
       if ( spNodeToNode(n,1) - spNodeToNode(n+1,1) == 0 ) then
          matchL(:,:) = spNodeToNode((/ n, n+1 /), :)
          matchR(:,:) = spNodeToNode((/ n+1, n /), :)

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

end subroutine MeshUtil1D_genConnectivity

subroutine MeshUtil1D_BuildInteriorMap( VMapM, VMapP, MapM, MapP, &
  & pos_en, pos_ev, EtoE, EtoF, EtoV, Fmask, Ne, Np, Nfp, Nfaces, Nv)


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

  real(RP), intent(in) :: pos_en(Np,Ne,1)
  real(RP), intent(in) :: pos_ev(Nv,1)
  integer, intent(in) :: EToE(Ne, Nfaces)
  integer, intent(in) :: EToF(Ne, Nfaces)
  integer, intent(in) :: EToV(Ne,2)
  integer, intent(in) :: Fmask(Nfp,2)


  integer :: k, k1, k2
  integer :: f, f1, f2
  integer :: p
  integer :: n
  integer :: idP, idM
  integer :: v1, v2
  integer :: nodeids(Np,Ne)
  real(RP) :: refd2
  real(RP) :: x1(Nfp,Nfp), x2(Nfp,Nfp)
  real(RP) :: dist(Nfp,Nfp)
  real(RP) :: x(Np*Ne)

  !-----------------------------------------------------------------------------

  do k=1, Ne
  do p=1, Np
     n = p + (k-1)*Np
     nodeids(p,k) = n
     x(n) = pos_en(p,k,1)
  end do
  end do

  do k=1, Ne
  do f=1, Nfaces
  do p=1, Nfp
     n = p + (f-1)*Nfp + (k-1)*Nfp*Nfaces
     MapM(p,f,k) = n
     MapP(p,f,k) = n
     VMapM(p,f,k) = nodeids(Fmask(p,f),k)
  end do
  end do
  end do

  VMapP = -1
  do k1=1, Ne
  do f1=1, Nfaces
     k2 = EToE(k1,f1); f2 = EToF(k1,f1)

     v1 = EToV(k1,f1); v2 = EToV(k1,1+mod(f1,Nfaces))
     refd2 =    (pos_ev(v1,1) - pos_ev(v2,1))**2

     x1(:,:) = spread( x(VMapM(:,f1,k1)), 2, Nfp )
     x2(:,:) = spread( x(VMapM(:,f2,k2)), 1, Nfp )
  
     dist(:,:) = (x1 - x2)**2
     do idP=1, Nfp
     do idM=1, Nfp
        if (dist(idM,idP)/refd2 < 1.0E-14_RP) then
           VMapP(idM,f1,k1) = VMapM(idP,f2,k2)
           MapP(idM,f1,k1) = idP + (f2-1)*Nfp + (k2-1)*Nfp*Nfaces
        end if
     end do
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

end subroutine MeshUtil1D_BuildInteriorMap

subroutine MeshUtil1D_genPatchBoundaryMap(  VMapB, MapB, VMapP, &
  & pos_en, xmin, xmax, Fmask, Ne, Np, Nfp, Nfaces, Nv)

  implicit none

  integer, intent(in) :: Ne
  integer, intent(in) :: Np
  integer, intent(in) :: Nfp
  integer, intent(in) :: Nfaces
  integer, intent(in) :: Nv

  integer, intent(inout), allocatable :: VMapB(:)
  integer, intent(inout), allocatable :: MapB(:)
  integer, intent(inout) :: VMapP(Nfp,Nfaces,Ne)

  real(RP), intent(in) :: pos_en(Np,Ne,1)
  real(RP), intent(in) :: xmin, xmax
  integer, intent(in) :: Fmask(Nfp,2)


  integer :: k
  integer :: b
  integer :: f
  integer :: i, j
  real(RP) :: x

   real(RP), parameter :: NODE_TOL = 1.0E-12_RP

   integer :: face

   integer :: elemIDs(Nfp*Ne,2)
   integer :: faceIDs(Nfp*Ne,2)
   integer :: counterB(2)
   integer :: mapB_counter

   !-----------------------------------------------------------------------------

   counterB(:) = 0d0

   do k=1, Ne
   do f=1, Nfaces
     x = sum(pos_en(Fmask(:,f),k,1)/dble(Nfp))

     call eval_domain_boundary(1, x, xmin, k, f)
     call eval_domain_boundary(2, x, xmax, k, f)
   end do
   end do


   allocate( mapB(sum(counterB*Nfp))    )
   allocate( vmapB(size(mapB)) )

   mapB_counter = 1
   do b = 1, 2
     ! write(*,*) "LocalMesh boundary ID:", b
     ! write(*,*) counterB(b)
     ! write(*,*) ordInfo(1:counterB(1),b)
     ! write(*,*) elemIds(1:counterB(1),b)
     ! write(*,*) faceIds(1:counterB(1),b)

     do i=1, counterB(b)
       k = elemIds(i,b); f = faceIDs(i,b)
       do j=1, Nfp
         VMapP(j,f,k) = Np*Ne + mapB_counter
         VmapB(mapB_counter) = Fmask(j,f) + (k-1)*Np
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

   contains
     subroutine eval_domain_boundary(domb_id, r, rbc, k, f)
       integer, intent(in) :: domb_id
       real(RP), intent(in) :: r
       real(RP), intent(in) :: rbc
       integer, intent(in) :: k, f

       if ( abs(r - rbc) < NODE_TOL ) then
         counterB(domB_ID) = counterB(domB_ID) + 1
         elemIds(counterB(domB_ID),domB_ID) = k
         faceIds(counterB(domB_ID),domB_ID) = f
       end if
     end subroutine eval_domain_boundary
 end subroutine MeshUtil1D_genPatchBoundaryMap

 subroutine MeshUtil1D_genPeriodicBoundaryMap( EToE, EToF, VMapP, &
    & xperiod, x, fx, Fmask, Ne, Np, Nfp, Nfaces, Nv)

   implicit none
   integer, intent(in) :: Ne
   integer, intent(in) :: Np
   integer, intent(in) :: Nfp
   integer, intent(in) :: Nfaces
   integer, intent(in) :: Nv   
   integer, intent(inout) :: EToE(Ne, Nfaces)
   integer, intent(inout) :: EToF(Ne, Nfaces)  
   integer, intent(inout) :: VMapP(Nfp,Nfaces,Ne)    
   real(RP), intent(in) :: xperiod
   real(RP), intent(in) :: x(Np,Ne)
   real(RP), intent(in) :: fx(Nfp*Ne)
   integer, intent(in) :: Fmask(Nfp,4)

   integer :: k, k1, k2
   integer :: f, f1, f2
   integer :: p
   integer :: n

   real(RP) :: cx1, cx2
   real(RP) :: dx_pbc
   integer :: idP, idM
   integer :: vidL(Nfp), vidR(Nfp)
   integer :: fidL(Nfp), fidR(Nfp)

   real(RP) :: refd2
   real(RP) :: x1(Nfp,Nfp), x2(Nfp,Nfp)
   real(RP) :: dist(Nfp,Nfp)

   real(RP), parameter :: NODE_TOL = 1.0E-12_RP

    !-----------------------------------------------------------------------------


    do k1=1, Ne
    do f1=1, Nfaces
     cx1 = sum(x(Fmask(:,f1),k1)/dble(Nfp))
     
     k2 = EToE(k1,f1); f2 = EToF(k1,f1)
     if (k2==k1) then
       do k2=1, Ne
         if (k1 /= k2) then
           do f2=1, Nfaces
             if (EToE(k2,f2) == k2) then
               cx2 = sum(x(Fmask(:,f2),k2)/dble(Nfp))
               
               dx_pbc = sqrt( (abs(cx1-cx2)-xperiod)**2 )
               
               if (dx_pbc < NODE_TOL) then

                 EtoE(k1,f1) = k2; EToE(k2,f2) = k1
                 EtoF(k1,f1) = f2; EToF(k2,f2) = f1

                 !mesh%BCType(f1,k1) = BCTYPE_PERIODIC
                 !mesh%BCType(f2,k2) = BCTYPE_PERIODIC

                 do p=1, Nfp
                   vidL(p) = Fmask(p,f1) + (k1-1)*Np
                   vidR(p) = Fmask(p,f2) + (k2-1)*Np
                   fidL(p) = p + (f1-1)*Nfp + (k1-1)*Nfp*Nfaces
                   fidR(p) = p + (f2-1)*Nfp + (k2-1)*Nfp*Nfaces
                 end do

                 x1(:,:) = spread( fx(fidL(:)), 2, Nfp )
                 x2(:,:) = spread( fx(fidR(:)), 1, Nfp )
                 
                 dist(:,:) = (abs(x1-x2) -xperiod)**2

                 do idP=1,Nfp
                 do idM=1,Nfp
                     if (dist(idM,idP) < NODE_TOL) then
                       VMapP(idM,f1,k1) = vidR(idP)
                       VMapP(idP,f2,k2) = vidL(idM)
                     end if
                 end do
                 end do

               end if
             end if
           end do
         end if
       end do
     end if

   end do
   end do

   !-------

!    write(*,*) "Build MapInfo (Periodic BC): "
!    do k=1,mesh%Ne
!       write(*,*) "k=", k, "---"
!       write(*,*) " VMapM:", mesh%VMapM(:,:,k)
!       write(*,*) " VMapP:", mesh%VMapP(:,:,k)
!       write(*,*) "BCType:", mesh%BCType(:,k)
!    end do

  end subroutine MeshUtil1D_genPeriodicBoundaryMap

  subroutine MeshUtil1D_buildGlobalMap( &
    & panelID_table, pi_table,    &
    & tileID_map, tileFaceID_map, tilePanelID_map, &
    & Ntile )
    
    use scale_prc, only: PRC_isMaster

    integer, intent(in) :: Ntile
    integer, intent(out) :: panelID_table(Ntile)
    integer, intent(out) :: pi_table(Ntile)
    integer, intent(out) :: tileID_map(2,Ntile)
    integer, intent(out) :: tileFaceID_map(2,Ntile)
    integer, intent(out) :: tilePanelID_map(2,Ntile)

    integer :: NtilePerPanel
    integer :: Ne, Nv
    integer, allocatable :: nodesID_1d(:)
    integer, allocatable :: EToV(:,:)
    integer, allocatable :: EToE(:,:)
    integer, allocatable :: EToF(:,:)
    integer :: i, f
    integer :: panelID
    integer :: tileID, tileID_R
    integer :: counter
    real(DP) :: del

    !-----------------------------------------------------------------------------

    NtilePerPanel = Ntile
    Ne = NtilePerPanel
    Nv = Ne + 1
    allocate( nodesID_1d(Nv) )
    allocate( EToV(Ntile,2), EToE(Ntile,2), EToF(Ntile,2) )

    nodesID_1d(:) = -1
    counter = 0
    do i = 2, Nv-1
      counter = counter + 1
      nodesID_1d(i) = counter
    end do

    !----

    tileID = 0

    do i = 1, Ne
      tileID = tileID + 1
      panelID_table(tileID) = 1
      pi_table(tileID) = i
      EToV(tileID,:) = (/ nodesID_1d(i), nodesID_1d(i+1) /)
    end do

    call MeshUtil1D_genConnectivity( EToE, EToF, &
      & EToV, Ntile, 2 )
    tileID_map(:,:) = transpose(EToE)
    tileFaceID_map(:,:) = transpose(EToF)

    do tileID=1, Ntile
    do f=1, 2
      tileID_R = tileID_map(f,tileID)
      tilePanelID_map(f,tileID) = panelID_table(tileID_R)
    end do
    end do

    !--

    ! if (PRC_isMaster) then
    !
    !   write(*,*) "TotTile", tileID, Ntile
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
    !
    ! end if
  end subroutine MeshUtil1D_buildGlobalMap

end module scale_meshutil_1d
