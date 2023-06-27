#include "scaleFElib.h"
module scale_meshutil_cubedsphere2d
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

  use scale_meshutil_2d, only: &
    MeshUtilCubedSphere2D_genRectDomain    => MeshUtil2D_genRectDomain,    &
    MeshUtilCubedSphere2D_genConnectivity  => MeshUtil2D_genConnectivity,  &
    MeshUtilCubedSphere2D_BuildInteriorMap => MeshUtil2D_BuildInteriorMap, &
    MeshUtilCubedSphere2D_genPatchBoundaryMap => MeshUtil2D_genPatchBoundaryMap
  !-----------------------------------------------------------------------------
  implicit none
  private
  
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MeshUtilCubedSphere2D_genRectDomain
  public :: MeshUtilCubedSphere2D_genConnectivity  
  public :: MeshUtilCubedSphere2D_buildInteriorMap
  public :: MeshUtilCubedSphere2D_buildGlobalMap
  public :: MeshUtilCubedSphere2D_genPatchBoundaryMap
  public :: MeshUtilCubedSphere2D_modifyConnectivity
  public :: MeshUtilCubedSphere2D_GetPanelConnectivity
  public :: MeshUtilCubedSphere2D_getPanelID
  
contains
!OCL SERIAL
  subroutine MeshUtilCubedSphere2D_buildGlobalMap( &
    panelID_table, pi_table, pj_table,              &
    tileID_map, tileFaceID_map, tilePanelID_map,    &
    Ntile  )
    
    use scale_meshutil_2d, only: &
      MeshUtil2D_genConnectivity
    implicit none

    integer, intent(in) :: Ntile
    integer, intent(out) :: panelID_table(Ntile)
    integer, intent(out) :: pi_table(Ntile)
    integer, intent(out) :: pj_table(Ntile)
    integer, intent(out) :: tileID_map(4,Ntile)
    integer, intent(out) :: tileFaceID_map(4,Ntile)
    integer, intent(out) :: tilePanelID_map(4,Ntile)

    integer :: NtilePerPanel
    integer :: NeX, NeY, NvX, NvY
    integer, allocatable :: nodesID_2d(:,:,:)
    integer, allocatable :: EToV(:,:)
    integer, allocatable :: EToE(:,:)
    integer, allocatable :: EToF(:,:)
    integer :: i, j, k, f
    integer :: panelID
    integer :: tileID, tileID_R
    integer :: counter
    !-----------------------------------------------------------------------------

    NtilePerPanel = Ntile / 6
    NeY = int( sqrt(dble(NtilePerPanel)) )
    NeX = NtilePerPanel/NeY
    NvX = NeX + 1
    NvY = NeY + 1
    allocate( nodesID_2d(NvX,NvY,6) )
    allocate( EToV(Ntile,4), EToE(Ntile,4), EToF(Ntile,4) )

    counter = 0
    do panelID = 1, 6
    do j = 1, NvY
    do i = 1, NvX
      counter = counter + 1
      nodesID_2d(i,j,panelID) = counter
    end do
    end do
    end do

    !----

    tileID = 0
    do panelID = 1, 6
    do j = 1, NeY
    do i = 1, NeX
      tileID = tileID + 1
      panelID_table(tileID) = panelID
      pi_table(tileID) = i; pj_table(tileID) = j;
      EToV(tileID,:) = (/ nodesID_2d(i,j  ,panelID), nodesID_2d(i+1,j  ,panelID),   &
                          nodesID_2d(i,j+1,panelID), nodesID_2d(i+1,j+1,panelID)   /)
    end do
    end do
    end do

    call MeshUtil2D_genConnectivity( EToE, EToF, &
      EToV, Ntile, 4 )
    tileID_map(:,:) = transpose(EToE)
    tileFaceID_map(:,:) = transpose(EToF)
    
    do tileID=1, Ntile
    do f=1, 4
      tileID_R = tileID_map(f,tileID)
      tilePanelID_map(f,tileID) = panelID_table(tileID_R)
    end do
    end do

    call MeshUtilCubedSphere2D_modifyConnectivity( &
      tilePanelID_map, tileID_map, tileFaceID_map,          & ! (inout)
      panelID_table, pi_table, pj_table, NeX, NeY, Ntile, 4 ) ! (in)

    return
  end subroutine MeshUtilCubedSphere2D_buildGlobalMap

  !----
!OCL SERIAL
  subroutine MeshUtilCubedSphere2D_modifyConnectivity( tilePanelID_map, tileID_map, tileFaceID_map, &
    panelID_table, pi_table, pj_table, NeX, NeY, Ntile, Nface )

    integer, intent(in) :: Ntile
    integer, intent(in) :: Nface
    integer, intent(out) :: tileID_map(Nface,Ntile)
    integer, intent(out) :: tileFaceID_map(Nface,Ntile)
    integer, intent(out) :: tilePanelID_map(Nface,Ntile)
    integer, intent(in) :: panelID_table(Ntile)
    integer, intent(in) :: pi_table(Ntile)
    integer, intent(in) :: pj_table(Ntile)
    integer, intent(in) :: NeX, NeY

    integer :: panel_connectivity(4,6)
    integer :: face_connectivity (4,6)

    integer :: tileID
    integer :: panelID
    integer :: f
    integer :: pi_, pj_    
    !-----------------------------------------------------------------------------

    call MeshUtilCubedSphere2D_getPanelConnectivity( &
      panel_connectivity, face_connectivity )

    
    do tileID=1, Ntile
      panelID = panelID_table(tileID)

      do f=1, 4
        if ( tileFaceID_map(f,tileID) /= f ) cycle ! Does the face correspond the boundary of panel of cubed sphere?

        pi_ = pi_table(tileID)          
        pj_ = pj_table(tileID)  
        
        select case( panelID )
        case ( 1 )
          if ( mod(f,2) == 0 ) then ! West / East
            if ( pi_table(tileID) == 1   ) pi_ = NeX
            if ( pi_table(tileID) == NeX ) pi_ = 1
          else                      ! North / South
            if ( pj_table(tileID) == 1   ) pj_ = NeY
            if ( pj_table(tileID) == NeY ) pj_ = 1  
          end if
        case ( 2 )
          if ( mod(f,2) == 0 ) then ! West / East
            if ( pi_table(tileID) == 1   ) pi_ = NeX
            if ( pi_table(tileID) == NeX ) pi_ = 1
          else                      ! North / South
            if ( pj_table(tileID) == 1   ) then
              pi_ = NeX; pj_ = NeY - pi_table(tileID) + 1
            end if
            if ( pj_table(tileID) == NeY ) then
              pi_ = NeX ; pj_ = pi_table(tileID)
            end if
          end if
        case ( 3 )
          if ( mod(f,2) == 0 ) then ! West / East
            if ( pi_table(tileID) == 1   ) pi_ = NeX
            if ( pi_table(tileID) == NeX ) pi_ = 1
          else                      ! North / South
            if ( pj_table(tileID) == 1   ) then
              pi_ = NeX - pi_table(tileID) + 1; pj_ = 1
            end if
            if ( pj_table(tileID) == NeY ) then
              pi_ = NeX - pi_table(tileID) + 1; pj_ = NeY
            end if
          end if
        case ( 4 )
          if ( mod(f,2) == 0 ) then ! West / East
            if ( pi_table(tileID) == 1   ) pi_ = NeX
            if ( pi_table(tileID) == NeX ) pi_ = 1
          else                      ! North / South 
            if ( pj_table(tileID) == 1   ) then
              pi_ = 1; pj_ = pi_table(tileID)
            end if
            if ( pj_table(tileID) == NeY ) then
              pi_ = 1; pj_ = NeY - pi_table(tileID) + 1; 
            end if
          end if          
        case ( 5 )
          pj_ = NeY
          if ( mod(f,2) == 0 ) then 
            if ( pi_table(tileID) == 1   ) pi_ = NeY - pj_table(tileID) + 1  ! West            
            if ( pi_table(tileID) == NeX ) pi_ = pj_table(tileID)            ! East
          else
            if ( pj_table(tileID) == 1   ) pi_ = pi_table(tileID)            ! South            
            if ( pj_table(tileID) == NeY ) pi_ = NeX - pi_table(tileID) + 1  ! North
          end if
        case ( 6 )
          pj_ = 1
          if ( mod(f,2) == 0 ) then 
            if ( pi_table(tileID) == 1   ) pi_ = pj_table(tileID)            ! West            
            if ( pi_table(tileID) == NeX ) pi_ = NeY - pj_table(tileID) + 1  ! East
          else
            if ( pj_table(tileID) == 1   ) pi_ = NeX - pi_table(tileID) + 1  ! South            
            if ( pj_table(tileID) == NeY ) pi_ = pi_table(tileID)            ! North
          end if
        end select

        tilePanelID_map(f,tileID) = panel_connectivity(f,panelID)
        tileID_map(f,tileID) = pi_ + (pj_ - 1) * NeX + (tilePanelID_map(f,tileID) - 1) * NeX * NeY
        tileFaceID_map(f,tileID) = face_connectivity(f,panelID)
      end do ! loop for f
    end do ! loop for tile

    ! if (PRC_myrank==0) then
    !   write(*,*) " MeshUtilCubedSphere2D_modifyConnectivity"
    !   do tileID=1, Ntile
    !     write(*,'(a,i3,a,2i3)') "tileID=", tileID, ":", pi_table(tileID), pj_table(tileID)
    !     write(*,'(a,6i4)') "map_tile:", tileID_map(:,tileID)
    !     write(*,'(a,6i4)') "map_panel:", tilePanelID_map(:,tileID)
    !     write(*,'(a,6i4)') "mac_face:", tileFaceID_map(:,tileID)
    !   end do
    ! end if

    return
  end subroutine MeshUtilCubedSphere2D_modifyConnectivity

!OCL SERIAL
  subroutine MeshUtilCubedSphere2D_getPanelConnectivity( panel_connectivity, face_connectivity )

    implicit none

    integer, intent(out) :: panel_connectivity(4,6)
    integer, intent(out) :: face_connectivity(4,6)

    integer :: n
    integer :: zonal_panelID_list(0:5)

    zonal_panelID_list(:) = (/ 4, 1, 2, 3, 4, 1/)
    do n=1, 4
      panel_connectivity(1,n) = 6
      if (zonal_panelID_list(4-n) > 2) then
        face_connectivity(1,n)  = +zonal_panelID_list(4-n)
      else
        face_connectivity(1,n) =  -zonal_panelID_list(4-n)
      end if

      panel_connectivity(2,n) = zonal_panelID_list(n+1)
      face_connectivity(2,n)  = 4

      panel_connectivity(3,n) = 5
      face_connectivity(3,n)  = n
      if (n > 2) then
        face_connectivity(3,n)  = -n
      else
        face_connectivity(3,n) = n
      end if

      panel_connectivity(4,n) = zonal_panelID_list(n-1)
      face_connectivity(4,n)  = 2
    end do

    panel_connectivity(:,5) = (/ 1, 2, 3, 4 /)
    face_connectivity(:,5) = (/ 3, 3, -3, -3 /)
    panel_connectivity(:,6) = (/ 3, 2, 1, 4 /)
    face_connectivity(:,6) = (/ -1, -1, 1, 1 /)

  end subroutine MeshUtilCubedSphere2D_getPanelConnectivity

!OCL SERIAL
  subroutine MeshUtilCubedSphere2D_getPanelID( panelID, lon, lat, Np )
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSPos
    implicit none

    integer, intent(in) :: Np
    integer, intent(out) :: panelID(Np)
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)

    integer :: p
    integer :: pnl
    real(RP) :: alph(Np), beta(Np)
    real(RP) :: lon_(Np), lat_(Np)

    real(RP), parameter :: EPS = 1.0E-64_RP
    !------------------------------------------

    !$omp parallel do
    do p=1, Np
      if ( abs(cos(lon(p))) < EPS ) then
        lon_(p) = lon(p) + EPS
      else
        lon_(p) = lon(p)
      end if
      if ( abs( lat(p) - 0.5_RP * PI ) < EPS ) then
        lat_(p) = lat(p) - sign(EPS, lat(p))
      else
        lat_(p) = lat(p)
      end if
    end do

    panelID(:) = -1
    do pnl=1, 6
      call CubedSphereCoordCnv_LonLat2CSPos( pnl, lon_, lat_, Np, &
        alph, beta )
    
      select case(pnl)
      case (5)
        where ( lat(:) > 0.0_RP .and. abs(alph(:)) <= 0.25_RP * PI .and. abs(beta(:)) <= 0.25_RP * PI  )
          panelID(:) = pnl
        end where
      case (6)
        where ( lat(:) < 0.0_RP .and. abs(alph(:)) <= 0.25_RP * PI .and. abs(beta(:)) <= 0.25_RP * PI  )
          panelID(:) = pnl
        end where
      case default 
        where ( abs(alph(:)) <= 0.25_RP * PI .and. abs(beta(:)) <= 0.25_RP * PI  )
          panelID(:) = pnl
        end where
      end select     
    end do

    do p=1, Np
      if (panelID(p) <  0) then
        LOG_ERROR("MeshUtilCubedSphere2D_getPanelID",*) 'Fail to search a panel ID of cubed sphere grid!'
        write(*,*) "p=", p, ": (lon,lat)=", lon(p), lat(p), "(alpha,beta)=", alph(p), beta(p) 
        call PRC_abort
      end if 
    end do

    return
  end subroutine MeshUtilCubedSphere2D_getPanelID

end module scale_meshutil_cubedsphere2d
