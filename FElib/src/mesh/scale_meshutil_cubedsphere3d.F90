#include "scaleFElib.h"
module scale_meshutil_cubedsphere3d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_const, only: &
    PI => CONST_PI,      &
    EPS => CONST_EPS
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prc

  use scale_meshutil_3d, only: &
    MeshUtilCubedSphere3D_genCubeDomain    => MeshUtil3D_genCubeDomain,    &
    MeshUtilCubedSphere3D_genConnectivity  => MeshUtil3D_genConnectivity,  &
    MeshUtilCubedSphere3D_BuildInteriorMap => MeshUtil3D_BuildInteriorMap, &
    MeshUtilCubedSphere3D_genPatchBoundaryMap => MeshUtil3D_genPatchBoundaryMap
  !-----------------------------------------------------------------------------
  implicit none
  private
  
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MeshUtilCubedSphere3D_genCubeDomain
  public :: MeshUtilCubedSphere3D_genConnectivity  
  public :: MeshUtilCubedSphere3D_buildInteriorMap
  public :: MeshUtilCubedSphere3D_buildGlobalMap
  public :: MeshUtilCubedSphere3D_genPatchBoundaryMap
  
contains
  subroutine MeshUtilCubedSphere3D_buildGlobalMap( &
    panelID_table, pi_table, pj_table, pk_table,    &
    tileID_map, tileFaceID_map, tilePanelID_map,    &
    Ntile,                                          &
    NeZ  )
    
    ! use scale_prc, only: PRC_isMaster
    use scale_meshutil_3d, only: &
      MeshUtil3D_genConnectivity
    use scale_meshutil_cubedsphere2d, only: &
      MeshUtilCubedSphere2D_getPanelConnectivity
    implicit none

    integer, intent(in) :: Ntile
    integer, intent(out) :: panelID_table(Ntile)
    integer, intent(out) :: pi_table(Ntile)
    integer, intent(out) :: pj_table(Ntile)
    integer, intent(out) :: pk_table(Ntile)
    integer, intent(out) :: tileID_map(6,Ntile)
    integer, intent(out) :: tileFaceID_map(6,Ntile)
    integer, intent(out) :: tilePanelID_map(6,Ntile)
    integer, intent(in) :: NeZ

    integer :: NtilePerPanel
    integer :: NeX, NeY, NvX, NvY, NvZ
    integer, allocatable :: nodesID_3d(:,:,:,:)
    integer, allocatable :: EToV(:,:)
    integer, allocatable :: EToE(:,:)
    integer, allocatable :: EToF(:,:)
    integer :: i, j, k, f
    integer :: panelID
    integer :: tileID, tileID_R
    integer :: counter

    integer :: panel_connectivity(4,6)
    integer :: face_connectivity (4,6)
    
    integer :: pi_, pj_
    !-----------------------------------------------------------------------------

    NtilePerPanel = Ntile / 6
    NeY = int( sqrt(dble(NtilePerPanel)) )
    NeX = NtilePerPanel/NeY
    NvX = NeX + 1
    NvY = NeY + 1
    NvZ = NeZ + 1
    allocate( nodesID_3d(NvX,NvY,NvZ,6) )
    allocate( EToV(Ntile,8), EToE(Ntile,8), EToF(Ntile,8) )

    counter = 0
    do panelID = 1, 6
    do k = 1, NvZ
    do j = 1, NvY
    do i = 1, NvX
      counter = counter + 1
      nodesID_3d(i,j,k,panelID) = counter
    end do
    end do
    end do
    end do

    !----

    tileID = 0
    do panelID = 1, 6
    do k = 1, NeZ
    do j = 1, NeY
    do i = 1, NeX
      tileID = tileID + 1
      panelID_table(tileID) = panelID
      pi_table(tileID) = i; pj_table(tileID) = j; pk_table(tileID) = k
      EToV(tileID,:) = (/ nodesID_3d(i,j  ,k  ,panelID), nodesID_3d(i+1,j  ,k  ,panelID), &
                          nodesID_3d(i,j+1,k  ,panelID), nodesID_3d(i+1,j+1,k  ,panelID), &
                          nodesID_3d(i,j  ,k+1,panelID), nodesID_3d(i+1,j  ,k+1,panelID), &
                          nodesID_3d(i,j+1,k+1,panelID), nodesID_3d(i+1,j+1,k+1,panelID)  /)
    end do
    end do
    end do
    end do

    call MeshUtil3D_genConnectivity( EToE, EToF, &
      EToV, Ntile, 6 )
    tileID_map(:,:) = transpose(EToE)
    tileFaceID_map(:,:) = transpose(EToF)

    call MeshUtilCubedSphere2D_getPanelConnectivity( &
      panel_connectivity, face_connectivity )
    
    do tileID=1, Ntile
    do f=1, 6
      tileID_R = tileID_map(f,tileID)
      tilePanelID_map(f,tileID) = panelID_table(tileID_R)
    end do
    end do

    !-
    do tileID=1, Ntile
      panelID = panelID_table(tileID)

      do f=1, 4
        if ( tileFaceID_map(f,tileID) /= f ) cycle ! Does the face correspond the boundary of panel of cubed sphere?

        pi_ = pi_table(tileID)          
        pj_ = pj_table(tileID)  
        
        select case( panelID )
        case ( 1, 2, 3, 4 )
          if ( pi_table(tileID) == 1   ) pi_ = NeX
          if ( pi_table(tileID) == NeX ) pi_ = 1
          if ( pj_table(tileID) == 1   ) pj_ = NeY
          if ( pj_table(tileID) == NeY ) pj_ = 1          
        case ( 5 )
          pj_ = NeY
          if ( pi_table(tileID) == 1   ) pi_ = NeY - pj_table(tileID) + 1  ! West
          if ( pi_table(tileID) == NeX ) pi_ = pj_table(tileID)            ! East
          if ( pj_table(tileID) == 1   ) pi_ = pi_table(tileID)            ! South
          if ( pj_table(tileID) == NeY ) pi_ = NeX - pi_table(tileID) + 1  ! North
        case ( 6 )
          pj_ = 1
          if ( pi_table(tileID) == 1   ) pi_ = pj_table(tileID)            ! West
          if ( pi_table(tileID) == NeX ) pi_ = NeY - pj_table(tileID) + 1  ! East
          if ( pj_table(tileID) == 1   ) pi_ = NeX - pi_table(tileID) + 1  ! South
          if ( pj_table(tileID) == NeY ) pi_ = pi_table(tileID)            ! North
        end select

        tilePanelID_map(f,tileID) = panel_connectivity(f,panelID)
        tileID_map(f,tileID) = pi_ + (pj_ - 1) * NeX + (tilePanelID_map(f,tileID) - 1) * NeX * NeY
        tileFaceID_map(f,tileID) = face_connectivity(f,panelID)
      end do ! loop for f
    end do ! loop for tile

    return
  end subroutine MeshUtilCubedSphere3D_buildGlobalMap

end module scale_meshutil_cubedsphere3d
