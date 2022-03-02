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
    class(LocalMesh3D), pointer :: lcmesh3D

    real(RP), allocatable :: in_tiles_x(:,:,:) 
    real(RP), allocatable :: in_tiles_y(:,:,:) 
    integer, allocatable  :: tileID_table(:,:) 
    integer, allocatable  :: panelID_table(:,:) 
    integer, allocatable :: pi_table(:)
    integer, allocatable :: pj_table(:)
    
    integer :: in_p, in_n
    integer :: i, j
    integer :: tileID    

    real(RP) :: delx, dely
    !------------------------------------

    call inmesh_dummy%Init( MESH_INID, inmesh_type )

    allocate( in_tiles_x(4,inmesh_dummy%NLocalMeshPerPrc,inmesh_dummy%Nprc) )
    allocate( in_tiles_y(4,inmesh_dummy%NLocalMeshPerPrc,inmesh_dummy%Nprc) )
    allocate( tileID_table(inmesh_dummy%NLocalMeshPerPrc,inmesh_dummy%Nprc) )
    allocate( panelID_table(inmesh_dummy%NLocalMeshPerPrc,inmesh_dummy%Nprc) )
    allocate( pi_table(inmesh_dummy%NLocalMeshPerPrc*inmesh_dummy%Nprc) )
    allocate( pj_table(inmesh_dummy%NLocalMeshPerPrc*inmesh_dummy%Nprc) )

    !---

    call inmesh_dummy%Get_inmesh_hmapinfo( &
      tileID_table, panelID_table, pi_table, pj_table ) ! (out)

    delx = ( inmesh_dummy%dom_xmax - inmesh_dummy%dom_xmin ) / dble(inmesh_dummy%NprcX)
    dely = ( inmesh_dummy%dom_ymax - inmesh_dummy%dom_ymin ) / dble(inmesh_dummy%NprcY)

    do in_p=1, inmesh_dummy%Nprc
    do in_n=1, inmesh_dummy%NLocalMeshPerPRC
      tileID = tileID_table(in_n,in_p)
      i = pi_table(tileID); j = pj_table(tileID)

      in_tiles_x(:,in_n,in_p) = inmesh_dummy%dom_xmin + delx * dble( (/ i-1, i, i, i-1 /) )
      in_tiles_y(:,in_n,in_p) = inmesh_dummy%dom_ymin + dely * dble( (/ j-1, j-1, j, j /) )

      if (i==inmesh_dummy%NprcX) then
        in_tiles_x(2:3,in_n,in_p) = in_tiles_x(2:3,in_n,in_p) + 1.0E-12_RP * delx       
      end if
      if (j==inmesh_dummy%NprcY) then
        in_tiles_y(3:4,in_n,in_p) = in_tiles_y(3:4,in_n,in_p) + 1.0E-12_RP * dely
      end if
    end do
    end do

    if ( associated(out_mesh%ptr_mesh2D ) ) then

      do n=1, out_mesh%ptr_mesh2D%LOCAL_MESH_NUM

        lcmesh2D => out_mesh%ptr_mesh2D%lcmesh_list(n)

        call NodeMap_construct_nodemap_2D( this, &
          inmesh_dummy%mesh_type_id, out_mesh%mesh_type_id,                  &
          inmesh_type,                                                       &
          lcmesh2D%refElem2D%Nfp, lcmesh2D%NeX * lcmesh2D%NeY, lcmesh2D,     &
          in_tiles_x, in_tiles_y, tileID_table, panelID_table,               &
          inmesh_dummy%Nprc, inmesh_dummy%NprcX, inmesh_dummy%NprcY,         &
          inmesh_dummy%NeX, inmesh_dummy%NeY, inmesh_dummy%NLocalMeshPerPRC, &
          .true.  )
        
      end do

    end if

    if ( associated(out_mesh%ptr_mesh3D ) ) then

      do n=1, out_mesh%ptr_mesh3D%LOCAL_MESH_NUM
            
        lcmesh3D => out_mesh%ptr_mesh3D%lcmesh_list(n)
        lcmesh2D => lcmesh3D%lcmesh2D

        call NodeMap_construct_nodemap_3D( this, &
          inmesh_dummy%mesh_type_id, out_mesh%mesh_type_id,                    &
          inmesh_type,                                                         &
          lcmesh3D%refElem3D%Nnode_h1D, lcmesh3D%NeX * lcmesh3D%NeY, lcmesh2D, &
          lcmesh3D%refElem3D%Np, lcmesh3D%Ne, lcmesh3D,                        &
          in_tiles_x, in_tiles_y, tileID_table, panelID_table,                 &
          inmesh_dummy%Nprc, inmesh_dummy%NprcX, inmesh_dummy%NprcY,           &
          inmesh_dummy%NeX, inmesh_dummy%NeY, inmesh_dummy%NLocalMeshPerPRC    )
        
      end do

    end if

    call inmesh_dummy%Final()

    return
  end subroutine NodeMap_Init

!OCL SERIAL  
  subroutine NodeMap_Final( this )
    implicit none
    class(regrid_nodemap), intent(inout) :: this

    integer :: n
    !------------------------------------
    
    deallocate( this%local_domID )
    deallocate( this%lcprc )
    deallocate( this%inPanelID )
    deallocate( this%elem_i, this%elem_j )
    deallocate( this%elem_x, this%elem_y )
    if ( allocated(this%elem_k) ) deallocate( this%elem_k, this%elem_z )

    do n=1, size(this%in_mesh_list)
      call this%in_mesh_list(n)%Final()
    end do
    deallocate( this%in_mesh_list )

    deallocate( this%in_tileID_list )

    return
  end subroutine NodeMap_Final

!-- private -----------------

!OCL SERIAL  
  subroutine NodeMap_construct_nodemap_2D( this, &
    in_meshtype_id, out_meshtype_id,             &
    in_meshtype_name,                            &
    Np1D, Ne2D, lcmesh,                          &
    tile_x, tile_y, tileID_table, panelID_table, &
    in_Nprc, in_NprcX, in_NprcY, in_NeX, in_NeY, &   
    in_NLocalMeshPerPrc,                         &
    do_mesh_generation,                          &
    in_lcprc2prc, in_prcnum_out                  )

    use scale_polygon, only: &
      polygon_inpoly    
        
    implicit none

    class(regrid_nodemap), intent(inout) :: this
    integer, intent(in) :: Np1D
    integer, intent(in) :: Ne2D
    class(LocalMesh2D), intent(in) :: lcmesh
    integer, intent(in) :: in_meshtype_id
    integer, intent(in) :: out_meshtype_id
    character(len=*), intent(in) :: in_meshtype_name
    integer, intent(in) :: in_Nprc
    integer, intent(in) :: in_NprcX
    integer, intent(in) :: in_NprcY
    integer, intent(in) :: in_NeX
    integer, intent(in) :: in_NeY
    integer, intent(in) :: in_NLocalMeshPerPrc
    real(RP), intent(in) :: tile_x(4,in_NLocalMeshPerPrc,in_Nprc)
    real(RP), intent(in) :: tile_y(4,in_NLocalMeshPerPrc,in_Nprc)
    integer, intent(in) :: tileID_table(in_NLocalMeshPerPrc,in_Nprc)
    integer, intent(in) :: panelID_table(in_NLocalMeshPerPrc,in_Nprc)
    logical, intent(in) :: do_mesh_generation
    integer, intent(out), optional :: in_lcprc2prc(in_Nprc)
    integer, intent(out), optional :: in_prcnum_out

    integer :: lc_domID, prcID, in_prc, in_n
    integer :: i, j
    integer :: ke_h
    integer :: p_h, p_h_x, p_h_y
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

    real(RP), allocatable :: out_x(:,:)
    real(RP), allocatable :: out_y(:,:)

    integer :: out_panel
    integer :: Np2D
    !-----------------------------------------------------------------

    LOG_INFO("regrid_nodemap",*) "NodeMap_construct_nodemap_2D_prep"    
    if( IO_L ) call flush(IO_FID_LOG)

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
   
    LOG_INFO("regrid_nodemap",*) "NodeMap_construct_nodemap_2D_mapping"    
    if( IO_L ) call flush(IO_FID_LOG)

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

    if ( present(in_lcprc2prc) ) in_lcprc2prc(:) = in_lcprc2prc_tmp(:)
    if ( present(in_prcnum_out) ) in_prcnum_out  = in_prc_num

    !-- prepair mesh for input data --------------------------------
    LOG_INFO("regrid_nodemap",*) "NodeMap_construct_nodemap_2D_mesh_gen"    
    if( IO_L ) call flush(IO_FID_LOG)

    if ( do_mesh_generation ) then

      allocate( this%in_mesh_list(in_prc_num)  )
      do i=1, in_prc_num
        call this%in_mesh_list(i)%Init( MESH_INID, in_meshtype_name )
        call this%in_mesh_list(i)%Generate( myrank=in_lcprc2prc_tmp(i)-1 )
      end do

    end if

    return
  end subroutine NodeMap_construct_nodemap_2D


!OCL SERIAL  
  subroutine NodeMap_construct_nodemap_3D( this, &
    in_meshtype_id, out_meshtype_id,             &
    in_meshtype_name,                            &
    Np1D, Ne2D, lcmesh2D,                        &
    Np3D, Ne3D, lcmesh3D,                        &
    tile_x, tile_y, tileID_table, panelID_table, &
    in_Nprc, in_NprcX, in_NprcY, in_NeX, in_NeY, &   
    in_NLocalMeshPerPrc                          )

    use scale_polygon, only: &
      polygon_inpoly    
        
    implicit none

    class(regrid_nodemap), intent(inout) :: this
    integer, intent(in) :: Np1D
    integer, intent(in) :: Ne2D
    class(LocalMesh2D), intent(in) :: lcmesh2D    
    integer, intent(in) :: Np3D
    integer, intent(in) :: Ne3D
    class(LocalMesh3D), intent(in) :: lcmesh3D    
    integer, intent(in) :: in_meshtype_id
    character(len=*), intent(in) :: in_meshtype_name
    integer, intent(in) :: out_meshtype_id
    integer, intent(in) :: in_Nprc
    integer, intent(in) :: in_NprcX
    integer, intent(in) :: in_NprcY
    integer, intent(in) :: in_NeX
    integer, intent(in) :: in_NeY
    integer, intent(in) :: in_NLocalMeshPerPrc
    real(RP), intent(in) :: tile_x(4,in_NLocalMeshPerPrc,in_Nprc)
    real(RP), intent(in) :: tile_y(4,in_NLocalMeshPerPrc,in_Nprc)
    integer, intent(in) :: tileID_table(in_NLocalMeshPerPrc,in_Nprc)
    integer, intent(in) :: panelID_table(in_NLocalMeshPerPrc,in_Nprc)
 
    integer :: in_lcprc2prc(in_Nprc)
    integer :: in_prc_num
    integer :: i
    !-----------------------------------------------------------------

    allocate( this%elem_k(Np3D, Ne3D) )
    allocate( this%elem_z(Np3D, Ne3D) )

    LOG_INFO("regrid_nodemap",*) "NodeMap_construct_nodemap_2D"
    if( IO_L ) call flush(IO_FID_LOG)

    call NodeMap_construct_nodemap_2D( this, &
      in_meshtype_id, out_meshtype_id,             &
      in_meshtype_name,                            &
      Np1D, Ne2D, lcmesh2D,                        &
      tile_x, tile_y, tileID_table, panelID_table, &
      in_Nprc, in_NprcX, in_NprcY, in_NeX, in_NeY, &   
      in_NLocalMeshPerPrc, .false.,                &
      in_lcprc2prc=in_lcprc2prc,                   &
      in_prcnum_out=in_prc_num                     )    

    !-- prepair mesh for input data --------------------------------

    LOG_INFO("regrid_nodemap",*) "NodeMap_construct_nodemap_3D"
    LOG_INFO("regrid_nodemap",*) "in_prc_num=", in_prc_num
    if( IO_L ) call flush(IO_FID_LOG)

    allocate( this%in_mesh_list(in_prc_num) )
    do i=1, in_prc_num
      LOG_INFO("regrid_nodemap",*) "i=", i
      if( IO_L ) call flush(IO_FID_LOG)
  
      LOG_INFO("regrid_nodemap",*) "Init in_mesh .."
      call this%in_mesh_list(i)%Init( MESH_INID, in_meshtype_name )
      LOG_INFO("regrid_nodemap",*) "Generate in_mesh .."
      call this%in_mesh_list(i)%Generate( myrank=in_lcprc2prc(i)-1 )
    end do

    return
  end subroutine NodeMap_construct_nodemap_3D  

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
        
    else if (     (       in_meshtype_id == REGRID_MESHTYPE_CUBEDSPHERE2D_ID    &
                   .and. out_meshtype_id == REGRID_MESHTYPE_CUBEDSPHERE2D_ID )  &
             .or. (       in_meshtype_id == REGRID_MESHTYPE_CUBEDSPHERE3D_ID    &
                   .and. out_meshtype_id == REGRID_MESHTYPE_CUBEDSPHERE3D_ID )  ) then
        
      ! We assume that panelID of input local mesh is equal to that of the corresponding local mesh. 
      inPanelID(:,:) = lcmesh%panelID

      ! Memo: 
      ! If we allow the panelID of input local mesh to be different from that of the corresponding local mesh, 
      ! the procedure should be the following. In addition, vector quantities in cubed sphere mesh need to be 
      ! temporally converted to that in a coordinate system that is commonly used between input and output local meshes.
      !-----
      ! call MeshUtilCubedSphere2D_getPanelID ( &
      !   inPanelID(:,:),                       & ! (out)
      !   lcmesh%lon(:,:), lcmesh%lat(:,:),     & ! (in)
      !   Np1D**2 * Ne2D                        ) ! (in)
      !------

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