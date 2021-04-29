!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_interp_mesh
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  public :: interp_mesh_Init
  public :: interp_mesh_Final

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, parameter :: ATMOS_MESH_NLocalMeshPerPrc = 1
  type(MeshCubeDom3D), public, target :: out_mesh

  type, public :: NodeMappingInfo
    type(LocalMesh3D), pointer :: lcmesh
    integer, allocatable :: prc_x(:,:)
    integer, allocatable :: prc_y(:,:)
    integer, allocatable :: elem_i(:,:)
    integer, allocatable :: elem_j(:,:)
    integer, allocatable :: elem_k(:,:)

    type(MeshCubeDom3D), allocatable :: in_mesh_list(:)
    integer, allocatable :: in_tileID_list(:)
    integer, allocatable :: prcXY2inListID(:,:)
    integer :: cubedom_icount
    integer :: cubedom_jcount
    integer :: cubedom_kcount

  contains
    procedure :: Init => NodeMappingInfo_Init
    procedure :: Final => NodeMappingInfo_Final
  end type
  type(NodeMappingInfo), public :: nodeMap_list(ATMOS_MESH_NLocalMeshPerPrc)

  integer, public :: in_NprcX         = 1        ! x length of 2D processor topology (input)
  integer, public :: in_NprcY         = 1        ! y length of 2D processor topology (input)  
  integer, public :: in_NeX           = 1
  integer, public :: in_NeY           = 1
  integer, public :: in_NeZ           = 1
  type(HexahedralElement), public :: in_elem3D

  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private :: out_NprcX        = 1        ! x length of 2D processor topology (output)
  integer, private :: out_NprcY        = 1        ! y length of 2D processor topology (output)
  type(HexahedralElement), private :: out_elem3D

contains
  subroutine interp_mesh_Init()
    implicit none

    integer :: out_NeX          = 1
    integer :: out_NeY          = 1
    integer :: out_NeZ          = 1
    integer :: in_PolyOrder_h   = 1
    integer :: in_PolyOrder_v   = 1  
    integer :: out_PolyOrder_h  = 1
    integer :: out_PolyOrder_v  = 1
    real(RP) :: dom_xmin        = 0.0_RP
    real(RP) :: dom_xmax        = 0.0_RP
    logical  :: isPeriodicX     = .false.
    real(RP) :: dom_ymin        = 0.0_RP
    real(RP) :: dom_ymax        = 0.0_RP
    logical  :: isPeriodicY     = .false.
    real(RP) :: dom_zmin        = 0.0_RP
    real(RP) :: dom_zmax        = 0.0_RP
    logical  :: isPeriodicZ     = .false.
 
    integer, parameter :: FZ_nmax = 1000
    real(RP) :: in_FZ(FZ_nmax)
    real(RP) :: out_FZ(FZ_nmax)

    namelist / PARAM_INTERP_MESH / &
      in_NprcX,        &
      in_NeX,          &
      in_NprcY,        &
      in_NeY,          &
      in_NeZ,          &
      in_PolyOrder_h,  &
      in_PolyOrder_v,  &    
      out_NprcX,       &
      out_NeX,         &
      out_NprcY,       &
      out_NeY,         &
      out_NeZ,         &
      out_PolyOrder_h, &
      out_PolyOrder_v, & 
      dom_xmin, dom_xmax, &
      isPeriodicX,        &
      dom_ymin, dom_ymax, &
      isPeriodicY,        &
      dom_zmin, dom_zmax, &
      isPeriodicZ,        &
      in_Fz, out_Fz

    integer :: ierr
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("interp_mesh",*) 'Setup'
  
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_INTERP_MESH,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
        LOG_INFO("interp_mesh",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("interp_mesh",*) 'Not appropriate names in namelist PARAM_INTERP_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_INTERP_MESH)

    
    !-
    call in_elem3D%Init( in_PolyOrder_h, in_PolyOrder_v, .true. )
    call out_elem3D%Init( out_PolyOrder_h, out_PolyOrder_v, .true. )

    call out_mesh%Init( out_NprcX * out_NeX, out_NprcY * out_NeY, out_NeZ,             &
       dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,                     &
       isPeriodicX, isPeriodicY, isPeriodicZ, out_elem3D, ATMOS_MESH_NLocalMeshPerPrc, &
       out_NprcX, out_NprcY,                                                           &
       FZ=out_FZ(1:out_NeZ+1) )
    
    call out_mesh%Generate()

    !-
    call construct_map( in_FZ(1:in_NeZ+1) )


    return
  end subroutine interp_mesh_Init

  subroutine interp_mesh_Final()
    implicit none

    integer :: n
    !-------------------------------------------

    call in_elem3D%Final()

    call out_mesh%Final()
    call out_elem3D%Final()

    do n=1, ATMOS_MESH_NLocalMeshPerPrc
      call nodeMap_list(n)%Final()
    end do
    
    return
  end subroutine interp_mesh_Final

  !- private -------------------------------------

  subroutine construct_map( in_Fz )
    implicit none

    real(RP), intent(in) :: in_Fz(1:in_NeZ+1)
    integer :: n
    integer :: i, j
    real(RP) :: in_tiles_x(4,in_NprcX,in_NprcY)
    real(RP) :: in_tiles_y(4,in_NprcX,in_NprcY)
    real(RP) :: delx, dely
    type(LocalMesh3D), pointer :: lcmesh
    !-------------------------------------------

    delx = ( out_mesh%xmax_gl - out_mesh%xmin_gl ) / dble(in_NprcX)
    dely = ( out_mesh%ymax_gl - out_mesh%ymin_gl ) / dble(in_NprcY)
    do j=1, in_NprcY
    do i=1, in_NprcX
      in_tiles_x(:,i,j) = out_mesh%xmin_gl + delx * dble( (/ i-1, i, i, i-1 /) )
      in_tiles_y(:,i,j) = out_mesh%ymin_gl + dely * dble( (/ j-1, j-1, j, j /) )
    end do
    end do
    in_tiles_x(2:3,in_NprcX,:) = in_tiles_x(2:3,in_NprcX,:) + 1.0E-12_RP * delx
    in_tiles_y(3:4,:,in_NprcY) = in_tiles_y(3:4,:,in_NprcY) + 1.0E-12_RP * dely

    do n=1, out_mesh%LOCAL_MESH_NUM
      lcmesh => out_mesh%lcmesh_list(n)
      call nodeMap_list(n)%Init( lcmesh, lcmesh%refElem3D, in_tiles_x, in_tiles_y, in_Fz )
    end do

    return
  end subroutine construct_map

  subroutine NodeMappingInfo_Init( this, lcmesh, elem, tile_x, tile_y, in_Fz )
    use scale_prc
    implicit none

    class(NodeMappingInfo), intent(inout), target :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    type(ElementBase3D), intent(in) :: elem
    real(RP), intent(in) :: tile_x(4,in_NprcX,in_NprcY)
    real(RP), intent(in) :: tile_y(4,in_NprcX,in_NprcY)
    real(RP), intent(in) :: in_Fz(1:in_NeZ+1)

    integer :: prc_i, prc_j
    integer :: i, j
    integer :: p, ke
    integer :: ke_h
    integer :: p_h, p_h_x, p_h_y
    integer :: ke_z, ke_z2, p_z, p_z2
    logical :: is_inside_tile(elem%Nnode_h1D**2)
    logical :: is_inside_elem
    real(RP) :: in_elem_x(4)
    real(RP) :: in_elem_y(4)
    real(RP) :: delx, dely
    logical :: target_tile_flag(in_NprcX,in_NprcY)
    integer :: in_tileID_tmp(in_NprcX*in_NprcY)
    integer :: in_tile_num
    integer :: in_rank
    type(LocalMesh3D), pointer :: in_lcmesh
    real(RP) :: in_Z0, in_Z1
    integer :: in_ke3D

    real(RP) :: out_x, out_y
    real(RP) :: out_x0, del_x
    real(RP) :: out_y0, del_y
    !-------------------------------------------

    allocate( this%prc_x(elem%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%prc_y(elem%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%prcXY2inListID(in_NprcX,in_NprcY) )
    allocate( this%elem_i(elem%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_j(elem%Nnode_h1D**2, lcmesh%NeX*lcmesh%NeY) )
    allocate( this%elem_k(elem%Np, lcmesh%NeX*lcmesh%NeY*lcmesh%NeZ) )

    target_tile_flag(:,:) = .false.
    this%prcXY2inListID(:,:) = -1

    in_tile_num = 0
    do ke_h=1, lcmesh%NeX * lcmesh%NeY
      do p_h_y=1, elem%Nnode_h1D
      do p_h_x=1, elem%Nnode_h1D
        p_h = p_h_x + (p_h_y-1)*elem%Nnode_h1D
        this%prc_x(p_h,ke_h) = -1
        this%prc_y(p_h,ke_h) = -1
        
        out_x = lcmesh%pos_en(elem%Hslice(p_h,1),ke_h,1)
        out_y = lcmesh%pos_en(elem%Hslice(p_h,1),ke_h,2)

        loop_prc: do prc_j=1, in_NprcY
        do prc_i=1, in_NprcX
          is_inside_tile(p_h) = inpoly( out_x, out_y,                                  &
                                        4, tile_x(:,prc_i,prc_j), tile_y(:,prc_i,prc_j))
          if ( is_inside_tile(p_h) ) then
            this%prc_x(p_h,ke_h) =  prc_i
            this%prc_y(p_h,ke_h) =  prc_j
            if ( .not. target_tile_flag(prc_i,prc_j) ) then
              target_tile_flag(prc_i,prc_j) = .true.
              in_tile_num = in_tile_num + 1
              in_tileID_tmp(in_tile_num) = prc_i + in_NprcX * (prc_j - 1)
              this%prcXY2inListID(prc_i,prc_j) = in_tile_num
            end if
            exit loop_prc
          end if

        end do
        end do loop_prc
      end do
      end do

    end do
    allocate( this%in_tileID_list(in_tile_num) )
    this%in_tileID_list(:) = in_tileID_tmp(1:in_tile_num)

    do ke_h=1, lcmesh%NeX * lcmesh%NeY
      do p_h_y=1, elem%Nnode_h1D
      do p_h_x=1, elem%Nnode_h1D

        p_h = p_h_x + (p_h_y-1)*elem%Nnode_h1D
        out_x = lcmesh%pos_en(elem%Hslice(p_h,1),ke_h,1)
        out_y = lcmesh%pos_en(elem%Hslice(p_h,1),ke_h,2)

        prc_i = this%prc_x(p_h,ke_h)
        prc_j = this%prc_y(p_h,ke_h)
        this%elem_i(p_h,ke_h) = -1
        this%elem_j(p_h,ke_h) = -1
        if ( prc_i > 0 .and. prc_j > 0) then
          delx = ( tile_x(2,prc_i,prc_j) - tile_x(1,prc_i,prc_j) ) / dble(in_NeX)
          dely = ( tile_y(4,prc_i,prc_j) - tile_y(1,prc_i,prc_j) ) / dble(in_NeY)
          loop_ne: do j=1, in_NeY
          do i=1, in_NeX        
            in_elem_x(:) = tile_x(1,prc_i,prc_j) + delx * dble( (/ i-1, i, i, i-1 /) )
            in_elem_y(:) = tile_y(1,prc_i,prc_j) + dely * dble( (/ j-1, j-1, j, j /) )
            is_inside_elem  = inpoly( out_x, out_y,                 &
                                      4, in_elem_x(:), in_elem_y(:) )
            if (is_inside_elem) then
              this%elem_i(p_h,ke_h) = i
              this%elem_j(p_h,ke_h) = j
              exit loop_ne
            end if
          end do  
          end do loop_ne
        end if
      end do
      end do
    end do

    !-- prepair mesh for input data --------------------------------

    allocate( this%in_mesh_list(in_tile_num) )
    do i=1, in_tile_num
      in_rank = this%in_tileID_list(i) - 1
      call this%in_mesh_list(i)%Init( in_NprcX*in_NeX, in_NprcY*in_NeY, in_NeZ,  &
        out_mesh%xmin_gl, out_mesh%xmax_gl, out_mesh%ymin_gl, out_mesh%ymax_gl,  &
        out_mesh%zmin_gl, out_mesh%zmax_gl,                                      &
        out_mesh%isPeriodicX, out_mesh%isPeriodicY, out_mesh%isPeriodicZ,        &
        in_elem3D, 1, in_NprcX, in_NprcY,                                        &
        nproc=in_NprcX*in_NprcY, myrank=in_rank,                                 &
        FZ=in_FZ(:)                                                              )
      
      call this%in_mesh_list(i)%Generate()
    end do

    !--
    this%elem_k(:,:) = -1

    do ke_h=1, lcmesh%NeX * lcmesh%NeY
    do p_h=1, elem%Nnode_h1D**2
      prc_i = this%prc_x(p_h,ke_h)
      prc_j = this%prc_y(p_h,ke_h)

      if ( prc_i > 0 .and. prc_j > 0 .and.                          &
          this%elem_i(p_h,ke_h) > 0 .and. this%elem_j(p_h,ke_h) > 0 ) then
        
        in_lcmesh => this%in_mesh_list( this%prcXY2inListID(prc_i,prc_j) )%lcmesh_list(1)
        do ke_z=1, lcmesh%NeZ
        do p_z=1, out_elem3D%Nnode_v
          ke = ke_h + (ke_z-1)*lcmesh%NeX*lcmesh%NeY
          p = p_h + (p_z-1)*elem%Nnode_h1D**2
          do ke_z2=1, in_NeZ
            in_ke3D = this%elem_i(p_h,ke_h) + (this%elem_j(p_h,ke_h)-1)*in_NeX &
                    + (ke_z2-1)*in_NeX*in_NeY
            in_Z0 = in_lcmesh%pos_ev(in_lcmesh%EToV(in_ke3D,1),3)
            in_Z1 = in_lcmesh%pos_ev(in_lcmesh%EToV(in_ke3D,5),3)
            if ( in_Z0 <= lcmesh%pos_en(p,ke,3) .and. lcmesh%pos_en(p,ke,3) <= in_Z1 ) then
              this%elem_k(p,ke) = ke_z2
              exit
            end if
          end do
        end do
        end do
      end if

    end do
    end do


    return
  end subroutine NodeMappingInfo_Init

  subroutine NodeMappingInfo_Final( this )
    implicit none
    class(NodeMappingInfo), intent(inout) :: this

    integer :: i
    !-------------------------------------------

    if ( allocated(this%prc_x) ) then
      deallocate( this%prc_x, this%prc_y )
      deallocate( this%elem_i, this%elem_j )
      deallocate( this%in_tileID_list )

      do i=1, size(this%in_mesh_list)
        call this%in_mesh_list(i)%Final()
      end do   
      deallocate( this%in_mesh_list )
    end if

    return
  end subroutine NodeMappingInfo_Final

  !> Check whether the point is located inside a polyngon
  function inpoly( pt_x, pt_y, num_node, v_x, v_y ) result(ret)
    implicit none
    real(RP), intent(in) :: pt_x
    real(RP), intent(in) :: pt_y
    integer, intent(in) :: num_node
    real(RP), intent(in) :: v_x(num_node)
    real(RP), intent(in) :: v_y(num_node)
    logical :: ret

    integer :: wn
    integer :: i, ii
    !------------------------------------------

    wn = 0
    do i=1, num_node
      ii = mod(i, num_node) + 1
      if ( v_y(i) <= pt_y .and. pt_y < v_y(ii)) then
        if( pt_x < v_x(i) + (pt_y - v_y(i)) * (v_x(ii) - v_x(i))/(v_y(ii) - v_y(i)) ) then
          wn = wn + 1
        end if
      else if ( v_y(i) > pt_y .and. v_y(ii) <= pt_y ) then
        if( pt_x < v_x(i) + (pt_y - v_y(i)) * (v_x(ii) - v_x(i))/(v_y(ii) - v_y(i)) ) then
          wn = wn - 1
        end if
      end if
    end do

    if (wn == 0) then
      ret = .false.
    else
      ret = .true.
    end if

    return
  end function inpoly


end module mod_interp_mesh