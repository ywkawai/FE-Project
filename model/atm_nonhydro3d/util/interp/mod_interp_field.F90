!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_interp_field
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_file_h
  use scale_prof
  use scale_prc, only: &
    PRC_myrank, PRC_abort

  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase3D
  use scale_meshfield_base, only: MeshField3D
  use scale_file_base_meshfield, only: &
    FILE_base_meshfield
  
  use mod_interp_mesh, only: &
    NodeMappingInfo

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  public :: interp_field_Init
  public :: interp_field_Final
  public :: interp_field_Interpolate

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type(MeshField3D), public :: out_var3D
  type, public :: OutVarInfo
    character(FILE_HSHORT) :: varname
    character(FILE_HSHORT) :: units
    integer :: num_step
    real(DP) :: dt
    real(DP) :: start_sec
  end type OutVarInfo
  integer, public :: out_var3D_num
  type(OutVarInfo), public, allocatable :: out_vinfo(:)
  
  character(len=H_LONG), public   :: in_basename      = ''       ! Basename of the input  file 

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: ITEM_MAX_NUM = 128

  type :: in_local_val
    real(RP), allocatable :: val(:,:)
    real(RP), allocatable :: spectral_coef(:,:,:,:)
  end type in_local_val

  type :: in_local_file_list
    type(FILE_base_meshfield), allocatable :: in_files(:)
  end type in_local_file_list
  type(in_local_file_list), private, allocatable, target :: in_file_list(:)
  logical, private :: is_cached_in_files

contains
  subroutine interp_field_Init( out_mesh, nodeMap_list )
    use scale_file_h
    implicit none
    class(MeshCubeDom3D), intent(in) :: out_mesh
    type(NodeMappingInfo), intent(in) :: nodeMap_list(:)

    integer :: nn
    character(len=H_SHORT)  :: vars(ITEM_MAX_NUM) = ''       ! name of variables

    namelist /PARAM_INTERP_FIELD/ &
      in_basename,     &
      vars
    integer :: ierr

    type(FILE_base_meshfield) :: in_file
    real(DP) :: time_endsec
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("interp_field",*) 'Setup'
  
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_INTERP_FIELD,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
        LOG_INFO("interp_field",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("interp_field",*) 'Not appropriate names in namelist PARAM_INTERP_FIELD. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_INTERP_FIELD)

    out_var3D_num = 0
    do nn= 1, ITEM_MAX_NUM
      if ( vars(nn) == '' ) then
        exit
      else
        out_var3D_num = out_var3D_num + 1
      end if
    end do
    allocate( out_vinfo(out_var3D_num) )

    !--

    call in_file%Init( out_var3D_num, mesh3D=out_mesh )
    call in_file%Open( in_basename, myrank=0 )
    do nn = 1, out_var3D_num
      out_vinfo(nn)%varname = vars(nn)
      call in_file%Get_dataInfo( vars(nn), istep=1, & ! (in)
        units=out_vinfo(nn)%units,                  & ! (out)
        time_start=out_vinfo(nn)%start_sec,         & ! (out)
        time_end=time_endsec                        ) ! (out)

      call in_file%Get_VarStepSize( vars(nn), & ! (in)
        out_vinfo(nn)%num_step                ) ! (out)
      
        out_vinfo(nn)%dt = time_endsec - out_vinfo(nn)%start_sec
    end do
    call in_file%Close()
    call in_file%Final()

    !--
    is_cached_in_files = .false.

    !
    call out_var3D%Init( "tmp", "1", out_mesh )

    return
  end subroutine interp_field_Init

  subroutine interp_field_Final()
    implicit none

    integer :: vid
    integer :: domID
    integer :: in_ii
    type(FILE_base_meshfield), pointer :: in_file_ptr
    !-------------------------------------------

    call out_var3D%Final()

    if (is_cached_in_files) then
      do domID=1, size(in_file_list)
        do in_ii=1, size(in_file_list(domID)%in_files)
          in_file_ptr => in_file_list(domID)%in_files(in_ii)
          call in_file_ptr%Close()
          call in_file_ptr%Final()
        end do
        deallocate( in_file_list(domID)%in_files )
      end do
      deallocate( in_file_list )
    end if

    return
  end subroutine interp_field_Final

  subroutine interp_field_Interpolate( istep, varname, out_mesh, out_field, nodeMap_list )
    use scale_mesh_cubedom3d, only: MeshCubeDom3D
    implicit none

    integer, intent(in) :: istep
    character(*), intent(in) :: varname
    class(MeshCubeDom3D), intent(in), target :: out_mesh
    class(MeshField3D), intent(inout) :: out_field
    type(NodeMappingInfo), intent(in) :: nodeMap_list(:)

    class(LocalMesh3D), pointer :: lcmesh
    integer :: n

    integer :: domID
    integer :: in_mesh_num
    integer :: in_ii
    integer :: in_rank
    type(FILE_base_meshfield), pointer :: in_file_ptr
    !-------------------------------------------

    call PROF_rapstart('INTERP_field_interpolate', 0)

    if (.not. is_cached_in_files) then
      allocate( in_file_list(out_mesh%LOCAL_MESH_NUM)  )
      do domID=1, out_mesh%LOCAL_MESH_NUM
        in_mesh_num = size(nodeMap_list(domID)%in_mesh_list)
        allocate( in_file_list(domID)%in_files(in_mesh_num) )
        do in_ii=1, in_mesh_num
          in_file_ptr => in_file_list(domID)%in_files(in_ii)
          in_rank = nodeMap_list(domID)%in_tileID_list(in_ii) - 1

          LOG_INFO("interp_field",'(a,i4,a,a,i6)') 'domID=', domID, 'Open in_file:', trim(in_basename), in_rank        
          call in_file_ptr%Init( out_var3D_num,            &
            mesh3D=nodeMap_list(domID)%in_mesh_list(in_ii) )
          call in_file_ptr%Open( in_basename, in_rank )
        end do
      end do
      is_cached_in_files = .true.
    end if

    do n=1, out_mesh%LOCAL_MESH_NUM
      lcmesh => out_mesh%lcmesh_list(n)
      call interpolate_local( out_field%local(n)%val(:,:),            &
        n, istep, varname, lcmesh, lcmesh%refElem3D, nodeMap_list(n), &
        out_mesh   )
    end do

    call PROF_rapend('INTERP_field_interpolate', 0)

    return
  end subroutine interp_field_interpolate

  subroutine interpolate_local( out_val, &
      out_domID, istep, varname, out_lcmesh, out_elem, mappingInfo, out_mesh3D )

    use scale_const, only: &
      UNDEF8 => CONST_UNDEF8
    use scale_file_base_meshfield, only: &
      FILE_base_meshfield
    use scale_file_common_meshfield, only: &
      MF3D_XYZT => FILE_COMMON_MESHFILED3D_DIMTYPEID_XYZT
    use scale_polynominal, only: &
      polynominal_genLegendrePoly
    use mod_interp_mesh, only: &
      in_NprcX, in_NprcY, in_NeX, in_NeY, in_NeZ,  &
      in_elem3D,                                   &
      LOCAL_MESH_NUM => ATMOS_MESH_NLocalMeshPerPrc

    implicit none
    class(LocalMesh3D), intent(in) :: out_lcmesh
    class(ElementBase3D), intent(in) :: out_elem
    real(DP), intent(out) :: out_val(out_elem%Np,out_lcmesh%NeA)
    integer, intent(in) :: out_domID
    integer, intent(in) :: istep
    character(*), intent(in) :: varname
    type(NodeMappingInfo), intent(in), target :: mappingInfo
    class(MeshCubeDom3D), intent(in) :: out_mesh3D

    integer :: ke3D, ke_h
    integer :: p, p_h, pX, pY, pZ
    integer :: p1, p2, p3, l
    integer :: elem_i, elem_j
    integer :: in_ke3D, in_ex, in_ey, in_ez, in_p3D, in_px, in_py, in_listID
    integer :: i
    integer :: in_tile_num
    integer :: in_rank
    type(LocalMesh3D), pointer :: in_lcmesh
    integer :: n
    integer :: ii, jj, kk, i0_s, j0_s, k0_s

    type(in_local_val) :: in_val_list(size(mappingInfo%in_tileID_list))
    type(MeshCubeDom3D), pointer :: in_mesh
    real(RP) :: in_spectral_coef(in_elem3D%Np,in_NeX,in_NeY,in_NeZ)
    real(RP) :: P1D_ori_x(1,in_elem3D%Nnode_h1D)
    real(RP) :: P1D_ori_y(1,in_elem3D%Nnode_h1D)
    real(RP) :: P1D_ori_z(1,in_elem3D%Nnode_v)
    real(RP) :: ox, oy, oz
    real(RP) :: vx(in_elem3D%Nv), vy(in_elem3D%Nv), vz(in_elem3D%Nv) 
    integer :: node_ids(in_elem3D%Nv)

    !---------------------------------------------

    in_tile_num = size(mappingInfo%in_tileID_list)
    do i=1, in_tile_num
      in_rank = mappingInfo%in_tileID_list(i) - 1
      in_mesh => mappingInfo%in_mesh_list(i)

      n = 1
      in_lcmesh => in_mesh%lcmesh_list(n)
      allocate( in_val_list(i)%val(in_lcmesh%refElem3D%Np,in_lcmesh%NeA) )
      allocate( in_val_list(i)%spectral_coef(in_lcmesh%refElem3D%Np,in_NeX,in_NeY,in_NeZ) )

      i0_s = 0; j0_s = 0; k0_s = 0
      do kk = 1, size(in_mesh%rcdomIJK2LCMeshID,3)
      do jj = 1, size(in_mesh%rcdomIJK2LCMeshID,2)
      do ii = 1, size(in_mesh%rcdomIJK2LCMeshID,1)
        if ( in_mesh%rcdomIJK2LCMeshID(ii,jj,kk) == n ) then
          call in_file_list(out_domID)%in_files(i)%Read_Var( &
            MF3D_XYZT, varname, in_lcmesh, &
            i0_s, j0_s, k0_s, in_val_list(i)%val(:,:), step=istep )
        end if
        i0_s = i0_s + in_lcmesh%NeX * in_lcmesh%refElem3D%Nnode_h1D
        j0_s = j0_s + in_lcmesh%NeY * in_lcmesh%refElem3D%Nnode_h1D
        k0_s = k0_s + in_lcmesh%NeZ * in_lcmesh%refElem3D%Nnode_v
      end do
      end do
      end do

      !$omp parallel do collapse(2) private(in_ke3D,ii) 
      do kk=1, in_NeZ
      do jj=1, in_NeY
      do ii=1, in_NeX
        in_ke3D = ii + (jj-1)*in_NeX + (kk-1)*in_NeX*in_NeY
        in_val_list(i)%spectral_coef(:,ii,jj,kk) = matmul( in_elem3D%invV, in_val_list(i)%val(:,in_ke3D) )
      end do
      end do
      end do    
    end do


    !$omp parallel do collapse(3) private( ke_h, ke3D, &
    !$omp pX, pY, pZ, p_h, p, in_px, in_py,            &
    !$omp in_ex, in_ey, in_ez, in_listID, in_ke3D,     &
    !$omp in_lcmesh, node_ids, vx, vy, vz, ox, oy, oz, &
    !$omp P1D_ori_x, P1D_ori_y, P1D_ori_z,             &
    !$omp p1, p2, p3, l                                )    
    do kk=1, out_lcmesh%NeZ
    do jj=1, out_lcmesh%NeY
    do ii=1, out_lcmesh%NeX
      ke_h = ii + (jj-1)*out_lcmesh%NeX
      ke3D = ke_h + (kk-1)*out_lcmesh%NeX*out_lcmesh%NeY

      do pZ=1, out_elem%Nnode_v
      do pY=1, out_elem%Nnode_h1D
      do pX=1, out_elem%Nnode_h1D
        p_h = pX + (pY-1)*out_elem%Nnode_h1D
        p = p_h + (pZ-1)*out_elem%Nnode_h1D**2
        out_val(p,ke3D) = UNDEF8

        in_px = mappingInfo%prc_x(p_h,ke_h)
        in_py = mappingInfo%prc_y(p_h,ke_h)
        if (in_px > 0 .and. in_py > 0) then
          in_ex = mappingInfo%elem_i(p_h,ke_h)
          in_ey = mappingInfo%elem_j(p_h,ke_h)
          in_ez = mappingInfo%elem_k(p,ke3D)
          in_listID = mappingInfo%prcXY2inListID(in_px,in_py)
          
          if ( in_ex > 0 .and. in_ey > 0 .and. in_ez > 0 ) then
            in_ke3D = in_ex + (in_ey-1)*in_NeX + (in_ez-1)*in_NeX*in_NeY
            in_lcmesh => mappingInfo%in_mesh_list(in_listID)%lcmesh_list(1)

            node_ids(:) = in_lcmesh%EToV(in_ke3D,:)
            vx(:) = in_lcmesh%pos_ev(node_ids(:),1)
            vy(:) = in_lcmesh%pos_ev(node_ids(:),2)
            vz(:) = in_lcmesh%pos_ev(node_ids(:),3)

            ox = - 1.0_RP + 2.0_RP * (out_lcmesh%pos_en(p,ke3D,1) - vx(1)) / (vx(2) - vx(1))
            oy = - 1.0_RP + 2.0_RP * (out_lcmesh%pos_en(p,ke3D,2) - vy(1)) / (vy(3) - vy(1))
            oz = - 1.0_RP + 2.0_RP * (out_lcmesh%pos_en(p,ke3D,3) - vz(1)) / (vz(5) - vz(1))

            P1D_ori_x(:,:) = polynominal_genLegendrePoly( in_elem3D%PolyOrder_h, (/ ox /) )
            P1D_ori_y(:,:) = polynominal_genLegendrePoly( in_elem3D%PolyOrder_h, (/ oy /) )
            P1D_ori_z(:,:) = polynominal_genLegendrePoly( in_elem3D%PolyOrder_v, (/ oz /) )

            out_val(p,ke3D) = 0.0_RP
            do p3=1, in_elem3D%Nnode_v
            do p2=1, in_elem3D%Nnode_h1D
            do p1=1, in_elem3D%Nnode_h1D
              l = p1 + (p2-1)*in_elem3D%Nnode_h1D + (p3-1)*in_elem3D%Nnode_h1D**2
              out_val(p,ke3D) = out_val(p,ke3D) + &
                  ( P1D_ori_x(1,p1) * P1D_ori_y(1,p2) * P1D_ori_z(1,p3) )                 &
                * sqrt((dble(p1-1) + 0.5_RP)*(dble(p2-1) + 0.5_RP)*(dble(p3-1) + 0.5_RP)) &
                * in_val_list(in_listID)%spectral_coef(l,in_ex,in_ey,in_ez)
            end do
            end do
            end do
          end if
        end if

      end do
      end do
      end do
    end do
    end do
    end do

    do i=1, in_tile_num
     deallocate( in_val_list(i)%val )
     deallocate( in_val_list(i)%spectral_coef )
    end do

    return
  end subroutine interpolate_local

end module mod_interp_field