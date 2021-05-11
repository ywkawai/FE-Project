!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_cs2lonlat_interp_field
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

  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_element_base, only: ElementBase2D
  use scale_meshfield_base, only: MeshField2D
  use scale_file_base_meshfield, only: &
    FILE_base_meshfield
  
  use mod_cs2lonlat_interp_mesh, only: &
    ATMOS_MESH_NLocalMeshPerPrc,       &
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
  type(MeshField2D), public :: out_var2D
  type, public :: OutVarInfo
    character(FILE_HSHORT) :: varname
    character(FILE_HSHORT) :: units
    integer :: num_step
    real(DP) :: dt
    real(DP) :: start_sec
    integer :: out_tintrv
  end type OutVarInfo
  integer, public :: out_var2D_num
  type(OutVarInfo), public, allocatable, target :: out_vinfo(:)
  
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
    real(RP), allocatable :: spectral_coef(:,:,:,:)
  end type in_local_val

  type :: in_local_file_list
    type(FILE_base_meshfield), allocatable :: in_files(:)
  end type in_local_file_list
  type(in_local_file_list), private, allocatable, target :: in_file_list(:)
  logical, private :: is_cached_in_files

contains
!OCL SERIAL
  subroutine interp_field_Init( out_mesh, nodeMap_list )
    use scale_file_h
    implicit none
    class(MeshRectDom2D), intent(in) :: out_mesh
    type(NodeMappingInfo), intent(in) :: nodeMap_list(:)

    integer :: nn
    character(len=H_SHORT)  :: vars(ITEM_MAX_NUM) = ''       ! name of variables
    integer :: out_tinterval(ITEM_MAX_NUM)

    namelist /PARAM_INTERP_FIELD/ &
      in_basename,     &
      vars,            &
      out_tinterval
    integer :: ierr

    type(FILE_base_meshfield) :: in_file
    real(DP) :: time_endsec
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("interp_field",*) 'Setup'
  
    out_tinterval(:) = 1

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

    out_var2D_num = 0
    do nn= 1, ITEM_MAX_NUM
      if ( vars(nn) == '' ) then
        exit
      else
        out_var2D_num = out_var2D_num + 1
      end if
    end do
    allocate( out_vinfo(out_var2D_num) )

    !--

    call in_file%Init( out_var2D_num, mesh2D=out_mesh )
    call in_file%Open( in_basename, myrank=0 )
    do nn = 1, out_var2D_num
      out_vinfo(nn)%varname = vars(nn)
      call in_file%Get_dataInfo( vars(nn), istep=1, & ! (in)
        units=out_vinfo(nn)%units,                  & ! (out)
        time_start=out_vinfo(nn)%start_sec,         & ! (out)
        time_end=time_endsec                        ) ! (out)

      call in_file%Get_VarStepSize( vars(nn), & ! (in)
        out_vinfo(nn)%num_step                ) ! (out)
      
      out_vinfo(nn)%dt = time_endsec - out_vinfo(nn)%start_sec
      out_vinfo(nn)%out_tintrv = out_tinterval(nn)
    end do
    call in_file%Close()
    call in_file%Final()

    !--
    is_cached_in_files = .false.

    !
    call out_var2D%Init( "tmp", "1", out_mesh )

    return
  end subroutine interp_field_Init

!OCL SERIAL
  subroutine interp_field_Final()
    implicit none

    integer :: vid
    integer :: domID
    integer :: in_ii
    type(FILE_base_meshfield), pointer :: in_file_ptr
    !-------------------------------------------

    call out_var2D%Final()

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

!OCL SERIAL
  subroutine interp_field_Interpolate( istep, varname, out_mesh, out_field, nodeMap_list )
    use scale_mesh_rectdom2d, only: MeshRectDom2D
    implicit none

    integer, intent(in) :: istep
    character(*), intent(in) :: varname
    class(MeshRectDom2D), intent(in), target :: out_mesh
    class(MeshField2D), intent(inout) :: out_field
    type(NodeMappingInfo), intent(in) :: nodeMap_list(:)

    class(LocalMesh2D), pointer :: lcmesh
    integer :: n

    integer :: domID
    integer :: in_mesh_num
    integer :: in_ii
    integer :: in_tileID
    integer :: in_rank
    type(FILE_base_meshfield), pointer :: in_file_ptr
    !-------------------------------------------

    call PROF_rapstart('INTERP_field_interpolate', 0)

    if (.not. is_cached_in_files) then
      allocate( in_file_list(out_mesh%LOCAL_MESH_NUM)  )
      do domID=1, out_mesh%LOCAL_MESH_NUM
        in_mesh_num = size(nodeMap_list(domID)%in_mesh2D_list)
        allocate( in_file_list(domID)%in_files(in_mesh_num) )
        do in_ii=1, in_mesh_num
          in_file_ptr => in_file_list(domID)%in_files(in_ii)
          in_tileID = nodeMap_list(domID)%in_tileID_list(in_ii)
          in_rank = nodeMap_list(domID)%in_mesh2D_list(in_ii)%PRCrank_globalMap(in_tileID)

          LOG_INFO("interp_field",'(a,i4,a,a,i6)') 'domID=', domID, 'Open in_file:', trim(in_basename), in_rank        
          call in_file_ptr%Init( out_var2D_num,                         &
            meshCubedSphere2D=nodeMap_list(domID)%in_mesh2D_list(in_ii) )
          call in_file_ptr%Open( in_basename, in_rank )
        end do
      end do
      is_cached_in_files = .true.
    end if

    do n=1, out_mesh%LOCAL_MESH_NUM
      lcmesh => out_mesh%lcmesh_list(n)
      call interpolate_local( out_field%local(n)%val(:,:),            &
        n, istep, varname, lcmesh, lcmesh%refElem2D, nodeMap_list(n), &
        out_mesh   )
    end do

    call PROF_rapend('INTERP_field_interpolate', 0)

    return
  end subroutine interp_field_interpolate

!OCL SERIAL
  subroutine interpolate_local( out_val, &
      out_domID, istep, varname, out_lcmesh, out_elem2D, mappingInfo, out_mesh2D )

    use scale_const, only: &
      UNDEF8 => CONST_UNDEF8, &
      PI => CONST_PI
    use scale_file_base_meshfield, only: &
      FILE_base_meshfield
    use scale_mesh_base2d, only: &
      MF2D_XYT => MeshBase2D_DIMTYPEID_XYT
    use scale_polynominal, only: &
      polynominal_genLegendrePoly
    use mod_cs2lonlat_interp_mesh, only: &
      nodeMap_list, &
      in_Nprc, in_elem2D, in_NLocalMeshPerPrc

    implicit none
    class(LocalMesh2D), intent(in) :: out_lcmesh
    class(ElementBase2D), intent(in) :: out_elem2D
    real(DP), intent(out) :: out_val(out_elem2D%Np,out_lcmesh%NeA)
    integer, intent(in) :: out_domID
    integer, intent(in) :: istep
    character(*), intent(in) :: varname
    type(NodeMappingInfo), intent(in), target :: mappingInfo
    class(MeshRectDom2D), intent(in) :: out_mesh2D

    integer :: ke2D, ke_h
    integer :: p, p_h, pX, pY
    integer :: p1, p2, l
    integer :: elem_i, elem_j
    integer :: in_ke2D, in_ex, in_ey, in_domID, in_prc
    integer :: in_rank
    type(LocalMesh2D), pointer :: in_lcmesh
    integer :: n
    integer :: ii, jj, pp, i0_s, j0_s, p0_s

    type(in_local_val), allocatable :: in_val_list(:)
    type(MeshCubedSphereDom2D), pointer :: in_csmesh
    real(RP) :: P1D_ori_x(1,in_elem2D%Nfp)
    real(RP) :: P1D_ori_y(1,in_elem2D%Nfp)
    real(RP) :: ox, oy
    real(RP) :: vx(in_elem2D%Nv), vy(in_elem2D%Nv)
    integer :: node_ids(in_elem2D%Nv)
    type(MeshField2D) :: tmp_field2D

    !---------------------------------------------

    allocate( in_val_list(in_Nprc) )

    do in_prc=1, in_Nprc
      in_rank = in_prc - 1
      in_csmesh => mappingInfo%in_mesh2D_list(in_prc)

      n = 1
      in_lcmesh => in_csmesh%lcmesh_list(n)
      allocate( in_val_list(in_prc)%spectral_coef(in_lcmesh%refElem2D%Np,in_lcmesh%NeX,in_lcmesh%NeY,in_csmesh%LOCAL_MESH_NUM) )

      call tmp_field2D%Init( varname, "", in_csmesh )
      call in_file_list(out_domID)%in_files(in_prc)%Read_Var( &
        MF2D_XYT, varname, tmp_field2D, step=istep            )

      !$omp parallel do collapse(2) private(in_ke2D,ii,jj) 
      do in_domID=1, in_csmesh%LOCAL_MESH_NUM
      do jj=1, in_lcmesh%NeY
      do ii=1, in_lcmesh%NeX
        in_ke2D = ii + (jj-1)*in_lcmesh%NeX
        in_val_list(in_prc)%spectral_coef(:,ii,jj,in_domID) = &
          matmul( in_elem2D%invV, tmp_field2D%local(in_domID)%val(:,in_ke2D) )
      end do
      end do
      end do

      call tmp_field2D%Final()
    end do


   !$omp parallel do collapse(2) private( ke_h, ke2D, &
   !$omp pX, pY, p_h, p, in_domID, in_prc,            &
   !$omp in_ex, in_ey, in_ke2D,                       &
   !$omp in_lcmesh, node_ids, vx, vy, ox, oy,         &
   !$omp P1D_ori_x, P1D_ori_y,                        &
   !$omp p1, p2, l                                )    
    do jj=1, out_lcmesh%NeY
    do ii=1, out_lcmesh%NeX
      ke_h = ii + (jj-1)*out_lcmesh%NeX
      ke2D = ke_h

      do pY=1, out_elem2D%Nfp
      do pX=1, out_elem2D%Nfp
        p_h = pX + (pY-1)*out_elem2D%Nfp
        p = p_h !+ (pZ-1)*out_elem2D%Nfp**2
        out_val(p,ke2D) = UNDEF8

        in_domID = mappingInfo%local_domID(p_h,ke_h)
        in_prc = mappingInfo%lcprc(p_h,ke_h)

        if (in_domID > 0 .and. in_prc > 0) then
          in_ex = mappingInfo%elem_i(p_h,ke_h)
          in_ey = mappingInfo%elem_j(p_h,ke_h)  
          
          if ( in_ex > 0 .and. in_ey  > 0 ) then
            in_lcmesh => mappingInfo%in_mesh2D_list(in_prc)%lcmesh_list(in_domID)
            in_ke2D = in_ex + (in_ey - 1)*in_lcmesh%NeX

            node_ids(:) = in_lcmesh%EToV(in_ke2D,:)
            vx(:) = in_lcmesh%pos_ev(node_ids(:),1)
            vy(:) = in_lcmesh%pos_ev(node_ids(:),2)

            ox = - 1.0_RP + 2.0_RP * (mappingInfo%elem_x(p,ke2D) - vx(1)) / (vx(2) - vx(1))
            oy = - 1.0_RP + 2.0_RP * (mappingInfo%elem_y(p,ke2D) - vy(1)) / (vy(3) - vy(1))

            if ( abs(ox) > 1.001_RP .or. abs(oy) > 1.001_RP ) then
              write(*,*) ii,jj,px,py,":",ox,oy
              stop
            end if
            P1D_ori_x(:,:) = polynominal_genLegendrePoly( in_elem2D%PolyOrder, (/ ox /) )
            P1D_ori_y(:,:) = polynominal_genLegendrePoly( in_elem2D%PolyOrder, (/ oy /) )

            out_val(p,ke2D) = 0.0_RP
            do p2=1, in_elem2D%Nfp
            do p1=1, in_elem2D%Nfp
              l = p1 + (p2-1)*in_elem2D%Nfp
              out_val(p,ke2D) = out_val(p,ke2D) + &
                  ( P1D_ori_x(1,p1) * P1D_ori_y(1,p2)  )                     &
                * sqrt((dble(p1-1) + 0.5_RP)*(dble(p2-1) + 0.5_RP))          &
                * in_val_list(in_prc)%spectral_coef(l,in_ex,in_ey,in_domID)
            end do
            end do
          end if
        end if

      end do
      end do
    end do
    end do

    do in_prc=1, in_Nprc
     deallocate( in_val_list(in_prc)%spectral_coef )
    end do

    return
  end subroutine interpolate_local

end module mod_cs2lonlat_interp_field