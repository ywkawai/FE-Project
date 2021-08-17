!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_regrid_nodemap
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort
  use scale_const, only: &
    PI => CONST_PI

  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D  

  use mod_regrid_mesh_base, only: &
    regrid_mesh_base,                    &
    MESH_INID => REGRID_MESH_BASE_IN_ID, &
    REGRID_MESHTYPE_CUBEDSPHERE2D_ID,    &
    REGRID_MESHTYPE_CUBEDSPHERE3D_ID,    &
    REGRID_MESHTYPE_LONLAT2D_ID,         &
    REGRID_MESHTYPE_LONLAT3D_ID

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  type, public :: regrid_nodemap
    integer, allocatable :: local_domID(:,:)
    integer, allocatable :: lcprc(:,:)
    integer, allocatable :: inPanelID(:,:)
    integer, allocatable :: elem_i(:,:)
    integer, allocatable :: elem_j(:,:)
    integer, allocatable :: elem_k(:,:)
    real(RP), allocatable :: elem_x(:,:)
    real(RP), allocatable :: elem_y(:,:)
    real(RP), allocatable :: elem_z(:,:)

    type(regrid_mesh_base), allocatable :: in_mesh_list(:)
    integer, allocatable :: in_tileID_list(:)
  contains
    procedure :: Init => Nodemap_Init
    procedure :: Final => Nodemap_Final
  end type

contains
!OCL SERIAL
  subroutine NodeMap_Init( this, &
      inmesh_type, out_mesh      )
    implicit none
    class(regrid_nodemap), intent(inout) :: this
    character(len=*), intent(in) :: inmesh_type
    class(regrid_mesh_base), intent(in), target :: out_mesh

    type(regrid_mesh_base) :: inmesh_dummy

    integer :: n
    class(LocalMesh2D), pointer :: lcmesh2D
    !------------------------------------

    call inmesh_dummy%Init( MESH_INID, inmesh_type )

    if ( associated(out_mesh%ptr_mesh2D ) ) then
      do n=1, out_mesh%ptr_mesh2D%LOCAL_MESH_NUM
        lcmesh2D => out_mesh%ptr_mesh2D%lcmesh_list(n)
        call NodeMap_construct_nodemap_2D( this, &
          inmesh_dummy%mesh_type_id, out_mesh%mesh_type_id,                  &
          lcmesh2D%refElem2D%Nfp, lcmesh2D%NeX * lcmesh2D%NeY, lcmesh2D,     &
          inmesh_dummy%Nprc, inmesh_dummy%NprcX, inmesh_dummy%NprcY,         &
          inmesh_dummy%NeX, inmesh_dummy%NeY, inmesh_dummy%NLocalMeshPerPRC  )
      end do

    end if
    if ( associated(out_mesh%ptr_mesh3D ) ) then
    end if

    call inmesh_dummy%Final()

    return
  end subroutine NodeMap_Init

!OCL SERIAL  
  subroutine NodeMap_Final( this )
    implicit none
    class(regrid_nodemap), intent(inout) :: this
    !------------------------------------

    return
  end subroutine NodeMap_Final

!-- private -----------------

!OCL SERIAL  
  subroutine NodeMap_construct_nodemap_2D( this, &
    in_meshtype_id, out_meshtype_id,             &
    Np1D, Ne2D, lcmesh,                          &
    in_Nprc, in_NprcX, in_NprcY, in_NeX, in_NeY, &   
    in_NLocalMeshPerPrc                 )

    use scale_polygon, only: &
      polygon_inpoly    
        
    implicit none

    class(regrid_nodemap), intent(inout) :: this
    integer, intent(in) :: Np1D
    integer, intent(in) :: Ne2D
    class(LocalMesh2D), intent(in) :: lcmesh
    integer, intent(in) :: in_meshtype_id
    integer, intent(in) :: out_meshtype_id
    integer, intent(in) :: in_Nprc
    integer, intent(in) :: in_NprcX
    integer, intent(in) :: in_NprcY
    integer, intent(in) :: in_NeX
    integer, intent(in) :: in_NeY
    integer, intent(in) :: in_NLocalMeshPerPrc

    real(RP) :: tile_x(4,in_NLocalMeshPerPrc,in_Nprc)
    real(RP) :: tile_y(4,in_NLocalMeshPerPrc,in_Nprc)
    integer  :: tileID_table(in_NLocalMeshPerPrc,in_Nprc)
    integer  :: panelID_table(in_NLocalMeshPerPrc,in_Nprc)
    
    integer :: lc_domID, prcID, in_prc, in_n
    integer :: i, j
    integer :: p, ke
    integer :: ke_h
    integer :: p_h, p_h_x, p_h_y
    integer :: ke_z, ke_z2, p_z, p_z2
    logical :: is_inside_tile(Np1D**2)
    logical :: is_inside_elem
    real(RP) :: in_elem_x(4)
    real(RP) :: in_elem_y(4)
    real(RP) :: delx, dely
    logical :: target_tile_flag(in_NLocalMeshPerPrc,in_Nprc)
    logical :: target_prc_flag(in_Nprc)
    integer :: in_tileID_tmp(in_NLocalMeshPerPrc*in_Nprc)
    integer :: in_prc2lcprc_tmp(in_Nprc)
    integer :: in_lcprc2prc_tmp(in_Nprc)
    integer :: in_tile_num
    integer :: in_prc_num    
    type(LocalMesh2D), pointer :: in_lcmesh

    real(RP), allocatable :: out_x(:,:)
    real(RP), allocatable :: out_y(:,:)

    integer :: out_panel
    real(RP) :: out_x0, del_x
    real(RP) :: out_y0, del_y

    integer :: Np2D
    !-----------------------------------------------------------------

    Np2D = Np1D**2
    allocate( this%local_domID(Np2D, Ne2D) )
    allocate( this%lcprc(Np2D, Ne2D) )
    allocate( this%inPanelID(Np2D, Ne2D) )
    allocate( this%elem_i(Np2D, Ne2D) )
    allocate( this%elem_j(Np2D, Ne2D) )
    allocate( this%elem_x(Np2D, Ne2D) )
    allocate( this%elem_y(Np2D, Ne2D) )  
    
    allocate( out_x(Np2D,Ne2D), out_y(Np2D,Ne2D) )

    target_tile_flag(:,:) = .false.
    target_prc_flag(:)    = .false.
    in_lcprc2prc_tmp(:)   = -1
    in_prc2lcprc_tmp(:)   = -1

    call get_panelID( this%inPanelID,  & ! (out)
      in_meshtype_id, out_meshtype_id, & ! (in)
      Np1D, Ne2D, lcmesh               ) ! (in)

    call get_out_xy( out_x, out_y,       & ! (out)
      in_meshtype_id, out_meshtype_id,   & ! (in)
      this%inPanelID,                    & ! (in)
      Np1D, Ne2D, lcmesh )                 ! (in)

    in_tile_num = 0
    in_prc_num  = 0
    do ke_h=1, Ne2D

      do p_h_y=1, Np1D
      do p_h_x=1, Np1D
        p_h = p_h_x + (p_h_y-1)*Np1D

        this%local_domID(p_h,ke_h) = -1
        this%lcprc(p_h,ke_h)       = -1

        out_panel = this%inPanelID(p_h,ke_h)

        loop_prc: do in_prc=1, in_Nprc
        do in_n=1, in_NLocalMeshPerPrc
          
          if ( out_panel /= panelID_table(in_n,in_prc) ) cycle

          is_inside_tile(p_h) = polygon_inpoly( out_x(p_h,ke_h), out_y(p_h,ke_h),               &
                                                4, tile_x(:,in_n,in_prc), tile_y(:,in_n,in_prc) )

          if ( is_inside_tile(p_h) ) then
            this%local_domID(p_h,ke_h) = in_n

            if ( .not. target_tile_flag(in_n,in_prc) ) then
              target_tile_flag(in_n,in_prc) = .true.
              in_tile_num = in_tile_num + 1
              in_tileID_tmp(in_tile_num) = tileID_table(in_n,in_prc)
            end if
            if ( .not. target_prc_flag(in_prc) ) then
              target_prc_flag(in_prc) = .true.
              in_prc_num = in_prc_num + 1
              in_prc2lcprc_tmp(in_prc) = in_prc_num
              in_lcprc2prc_tmp(in_prc_num) = in_prc
            end if
            this%lcprc(p_h,ke_h) = in_prc2lcprc_tmp(in_prc)

            exit loop_prc
          end if

        end do
        end do loop_prc

      end do
      end do

    end do


    allocate( this%in_tileID_list(in_tile_num) )
    this%in_tileID_list(:) = in_tileID_tmp(1:in_tile_num)
  
    do ke_h=1, Ne2D
      do p_h_y=1, Np1D
      do p_h_x=1, Np1D

        p_h = p_h_x + (p_h_y-1)*Np1D

        lc_domID = this%local_domID(p_h,ke_h)
        prcID = in_lcprc2prc_tmp( this%lcprc(p_h,ke_h) )
        out_panel = this%inPanelID(p_h,ke_h)

        this%elem_i(p_h,ke_h) = -1
        this%elem_j(p_h,ke_h) = -1
        if ( lc_domID > 0 .and. prcID > 0) then
          delx = ( tile_x(2,lc_domID,prcID) - tile_x(1,lc_domID,prcID) ) / dble(in_NeX)
          dely = ( tile_y(4,lc_domID,prcID) - tile_y(1,lc_domID,prcID) ) / dble(in_NeY) 
  
          loop_ne: do j=1, in_NeY
          do i=1, in_NeX        
            in_elem_x(:) = tile_x(1,lc_domID,prcID) + delx * dble( (/ i-1, i, i, i-1 /) )
            in_elem_y(:) = tile_y(1,lc_domID,prcID) + dely * dble( (/ j-1, j-1, j, j /) )

            if (i==in_NeX) then
              in_elem_x(2:3) = in_elem_x(2:3) + 1.0E-12_RP * delx       
            end if
            if (j==in_NeY) then
              in_elem_y(3:4) = in_elem_y(3:4) + 1.0E-12_RP * dely
            end if
                  
            is_inside_elem  = polygon_inpoly( out_x(p_h,ke_h), out_y(p_h,ke_h),  &
                                              4, in_elem_x(:), in_elem_y(:)      )
                        
            if (is_inside_elem) then
              this%elem_i(p_h,ke_h) = i
              this%elem_j(p_h,ke_h) = j
              this%elem_x(p_h,ke_h) = out_x(p_h,ke_h)
              this%elem_y(p_h,ke_h) = out_y(p_h,ke_h)

              exit loop_ne
            end if
          end do  
          end do loop_ne
        end if
      end do
      end do
    end do

    !-- prepair mesh for input data --------------------------------


    return
  end subroutine NodeMap_construct_nodemap_2D

!OCL SERIAL
  subroutine get_inmesh_mapinfo( &
    in_meshtype_id, in_Nprc, in_NLocalMeshPerPrc )
    implicit none

    integer, intent(in) :: in_meshtype_id
    integer, intent(in) :: in_Nprc
    integer, intent(in) :: in_NLocalMeshPerPrc
    !-----------------------------------------------------
    
    return
  end subroutine get_inmesh_mapinfo

!OCL SERIAL  
  subroutine get_panelID( inPanelID, &
    in_meshtype_id, out_meshtype_id, &
    Np1D, Ne2D, lcmesh ) 

    use scale_meshutil_cubedsphere2d, only: &
      MeshUtilCubedSphere2D_getPanelID     
    implicit none

    integer, intent(in) :: Np1D
    integer, intent(in) :: Ne2D
    class(LocalMesh2D), intent(in) :: lcmesh
    integer, intent(out) :: inPanelID(Np1D**2,Ne2D)
    integer, intent(in) :: in_meshtype_id
    integer, intent(in) :: out_meshtype_id
    !----------------------------------------------

    if (     (       in_meshtype_id == REGRID_MESHTYPE_CUBEDSPHERE2D_ID    &
              .and. out_meshtype_id == REGRID_MESHTYPE_LONLAT2D_ID      )  &
        .or. (       in_meshtype_id == REGRID_MESHTYPE_CUBEDSPHERE3D_ID    &
              .and. out_meshtype_id == REGRID_MESHTYPE_LONLAT3D_ID      )  ) then
      

      call MeshUtilCubedSphere2D_getPanelID ( &
        inPanelID(:,:),                             & ! (out)
        lcmesh%pos_en(:,:,1) * PI / 180.0_RP,       & ! (in)
        lcmesh%pos_en(:,:,2) * PI / 180.0_RP,       & ! (in)
        Np1D**2 * Ne2D                              ) ! (in)     
    else
      inPanelID(:,:) = 1
    end if

    return
  end subroutine get_panelID

!OCL SERIAL  
  subroutine get_out_xy( out_x, out_y, &
    in_meshtype_id, out_meshtype_id,   &
    inPanelID,                         &
    Np1D, Ne2D, lcmesh ) 

    use scale_cubedsphere_cnv, only: &
      CubedSphereCnv_LonLat2CSPos   
    implicit none

    integer, intent(in) :: Np1D
    integer, intent(in) :: Ne2D
    class(LocalMesh2D), intent(in) :: lcmesh
    real(RP), intent(out) :: out_x(Np1D**2,Ne2D)
    real(RP), intent(out) :: out_y(Np1D**2,Ne2D)
    integer, intent(in) :: in_meshtype_id
    integer, intent(in) :: out_meshtype_id
    integer, intent(in) :: inPanelID(Np1D**2,Ne2D)

    integer :: p_h, ke_h

    real(RP) :: out_lon(1)
    real(RP) :: out_lat(1)

    !----------------------------------------------

    if (     (       in_meshtype_id == REGRID_MESHTYPE_CUBEDSPHERE2D_ID    &
              .and. out_meshtype_id == REGRID_MESHTYPE_LONLAT2D_ID      )  &
        .or. (       in_meshtype_id == REGRID_MESHTYPE_CUBEDSPHERE3D_ID    &
              .and. out_meshtype_id == REGRID_MESHTYPE_LONLAT3D_ID      )  ) then
      
      do ke_h=1, Ne2D
      do p_h=1, Np1D**2
        out_lon(1) = lcmesh%pos_en(p_h,ke_h,1) * PI / 180.0_RP
        out_lat(1) = lcmesh%pos_en(p_h,ke_h,2) * PI / 180.0_RP

        call CubedSphereCnv_LonLat2CSPos( &
          inPanelID(p_h,ke_h),                       & ! (in)
          out_lon(1), out_lat(1), 1,                 & ! (in)
          out_x(p_h,ke_h), out_y(p_h,ke_h)           ) ! (out)
      end do
      end do
    else
      !$omp parallel do private(p_h)
      do ke_h=1, Ne2D
      do p_h=1, Np1D**2
        out_x(p_h,ke_h) = lcmesh%pos_en(p_h,ke_h,1)
        out_y(p_h,ke_h) = lcmesh%pos_en(p_h,ke_h,2)
      end do
      end do
    end if

    return
  end subroutine get_out_xy

end module mod_regrid_nodemap