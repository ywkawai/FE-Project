!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_regrid_interp_field
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
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_meshfield_base, only: MeshField2D, MeshField3D
  use scale_file_base_meshfield, only: &
    FILE_base_meshfield
  
  use mod_regrid_mesh_base, only: &
    regrid_mesh_base
  use mod_regrid_nodemap, only: &
    regrid_nodemap
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  public :: regrid_interp_field_Init
  public :: regrid_interp_field_Final
  interface regrid_interp_field_Interpolate
    module procedure :: regrid_field_Interpolate_2D
    module procedure :: regrid_field_Interpolate_3D
  end interface
  public :: regrid_interp_field_Interpolate

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type(MeshField2D), public :: out_var2D
  type(MeshField3D), public :: out_var3D
  type, public :: OutVarInfo
    character(FILE_HSHORT) :: varname
    character(FILE_HSHORT) :: units
    integer :: num_step
    real(DP) :: dt
    real(DP) :: start_sec
    integer :: out_tintrv
  end type OutVarInfo
  integer, public :: out_var_num
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
    real(RP), allocatable :: spectral_coef2D(:,:,:,:)
    real(RP), allocatable :: spectral_coef3D(:,:,:,:,:)
  end type in_local_val

  type :: in_local_file_list
    type(FILE_base_meshfield), allocatable :: in_files(:)
  end type in_local_file_list
  type(in_local_file_list), private, allocatable, target :: in_file_list(:)
  logical, private :: is_cached_in_files

contains
!OCL SERIAL
  subroutine regrid_interp_field_Init( out_mesh )
    use scale_const, only: &
      EPS => CONST_EPS
    use scale_file_h
    implicit none
    class(regrid_mesh_base), intent(in), target :: out_mesh

    integer :: nn
    character(len=H_SHORT)  :: vars(ITEM_MAX_NUM) = ''       ! name of variables
    integer :: out_tinterval(ITEM_MAX_NUM)

    namelist /PARAM_REGRID_INTERP_FIELD/ &
      in_basename,     &
      vars,            &
      out_tinterval
    integer :: ierr

    type(FILE_base_meshfield) :: in_file
    real(DP) :: time_endsec

    class(MeshBase2D), pointer :: ptr_mesh2D
    class(MeshBase3D), pointer :: ptr_mesh3D
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("regrid_interp_field",*) 'Setup'
  
    out_tinterval(:) = 1

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_REGRID_INTERP_FIELD,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_interp_field",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_interp_field",*) 'Not appropriate names in namelist PARAM_REGRID_INTERP_FIELD. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_REGRID_INTERP_FIELD)

    out_var_num = 0
    do nn= 1, ITEM_MAX_NUM
      if ( vars(nn) == '' ) then
        exit
      else
        out_var_num = out_var_num + 1
      end if
    end do
    allocate( out_vinfo(out_var_num) )

    !--

    if ( associated( out_mesh%ptr_mesh2D ) ) then

      select type( ptr_mesh2D => out_mesh%ptr_mesh2D  )
      class is ( MeshRectDom2D )
        call in_file%Init( out_var_num, mesh2D=ptr_mesh2D )
      class is ( MeshCubedSphereDom2D )
        call in_file%Init( out_var_num, meshCubedSphere2D=ptr_mesh2D )
      end select

      call out_var2D%Init( "tmp", "1", out_mesh%ptr_mesh2D )

    else if ( associated( out_mesh%ptr_mesh3D ) ) then

      select type( ptr_mesh3D => out_mesh%ptr_mesh3D  )
      class is ( MeshCubeDom3D )
        call in_file%Init( out_var_num, mesh3D=ptr_mesh3D )
      class is ( MeshCubedSphereDom3D )
        call in_file%Init( out_var_num, meshCubedSphere3D=ptr_mesh3D )
      end select

      call out_var3D%Init( "tmp", "1", out_mesh%ptr_mesh3D )
    end if

    call in_file%Open( in_basename, myrank=0 )
    do nn = 1, out_var_num
      out_vinfo(nn)%varname = vars(nn)
      call in_file%Get_dataInfo( vars(nn), istep=1, & ! (in)
        units=out_vinfo(nn)%units,                  & ! (out)
        time_start=out_vinfo(nn)%start_sec,         & ! (out)
        time_end=time_endsec                        ) ! (out)

      call in_file%Get_VarStepSize( vars(nn), & ! (in)
        out_vinfo(nn)%num_step                ) ! (out)
      
      out_vinfo(nn)%dt = time_endsec - out_vinfo(nn)%start_sec
      out_vinfo(nn)%out_tintrv = out_tinterval(nn)

      if (       abs(out_vinfo(nn)%dt) < EPS &
           .and. out_vinfo(nn)%num_step == 0 ) then
        out_vinfo(nn)%num_step   = 1
        out_vinfo(nn)%out_tintrv = 1
      end if

      LOG_INFO("regrid_interp_field_Init", '(3a,i4,a,i4)') &
        " Regist: name=", trim(vars(nn)), ", out_nstep=", out_vinfo(nn)%num_step, &
        ", out_tinterval=", out_vinfo(nn)%out_tintrv
    end do
    
    call in_file%Close()
    call in_file%Final()

    !--
    is_cached_in_files = .false.

    return
  end subroutine regrid_interp_field_Init

!OCL SERIAL
  subroutine regrid_interp_field_Final( out_mesh )
    implicit none

    class(regrid_mesh_base), intent(in), target :: out_mesh

    integer :: domID
    integer :: in_ii
    type(FILE_base_meshfield), pointer :: in_file_ptr
    !-------------------------------------------

    if ( associated( out_mesh%ptr_mesh2D ) ) then
      call out_var2D%Final()
    end if
    if ( associated( out_mesh%ptr_mesh3D ) ) then
      call out_var3D%Final()
    end if

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
  end subroutine regrid_interp_field_Final

!OCL SERIAL
  subroutine regrid_field_Interpolate_2D( istep, varname, out_mesh, out_field, nodeMap_list )
    use scale_mesh_rectdom2d, only: MeshRectDom2D
    implicit none

    integer, intent(in) :: istep
    character(*), intent(in) :: varname
    class(regrid_mesh_base), intent(in), target :: out_mesh
    class(MeshField2D), intent(inout) :: out_field
    type(regrid_nodemap), intent(in) :: nodeMap_list(:)

    class(LocalMesh2D), pointer :: lcmesh
    integer :: n

    integer :: domID
    integer :: in_mesh_num
    integer :: in_ii
    integer :: in_rank
    type(FILE_base_meshfield), pointer :: in_file_ptr

    class(MeshBase2D), pointer :: ptr_outmesh2D
    class(MeshBase2D), pointer :: ptr_inmesh2D
    !-------------------------------------------

    call PROF_rapstart('regrid_field_interpolate_2D', 0)

    ptr_outmesh2D => out_mesh%ptr_mesh2D

    if (.not. is_cached_in_files) then
      allocate( in_file_list(ptr_outmesh2D%LOCAL_MESH_NUM)  )
      do domID=1, ptr_outmesh2D%LOCAL_MESH_NUM
        in_mesh_num = size(nodeMap_list(domID)%in_mesh_list)
        allocate( in_file_list(domID)%in_files(in_mesh_num) )
        do in_ii=1, in_mesh_num
          in_file_ptr => in_file_list(domID)%in_files(in_ii)
          ptr_inmesh2D => nodeMap_list(domID)%in_mesh_list(in_ii)%ptr_mesh2D
          in_rank = ptr_inmesh2D%lcmesh_list(1)%PRC_myrank
          LOG_INFO("regrid_field",'(a,i4,a,a,i6)') 'domID=', domID, ', Open in_file:', trim(in_basename), in_rank        

          select type( ptr_inmesh2D )
          class is (MeshRectDom2D)
            call in_file_ptr%Init( out_var_num, mesh2D=ptr_inmesh2D )
          class is ( MeshCubedSphereDom2D )
            call in_file_ptr%Init( out_var_num, meshCubedSphere2D=ptr_inmesh2D )
          end select
          call in_file_ptr%Open( in_basename, in_rank )
        end do
      end do
      is_cached_in_files = .true.
    end if

    do n=1, ptr_outmesh2D%LOCAL_MESH_NUM
      lcmesh => ptr_outmesh2D%lcmesh_list(n)
      call interpolate_local_2D( out_field%local(n)%val(:,:),                                 &
        n, istep, varname, lcmesh, lcmesh%refElem2D, nodeMap_list(n)%in_mesh_list(1)%elem2D,  &
        nodeMap_list(n)  )
    end do

    call PROF_rapend('regrid_field_interpolate_2D', 0)

    return
  end subroutine regrid_field_Interpolate_2D

!OCL SERIAL
  subroutine regrid_field_Interpolate_3D( istep, varname, out_mesh, out_field, nodeMap_list )
    use scale_mesh_cubedom3d, only: MeshCubeDom3D
    implicit none

    integer, intent(in) :: istep
    character(*), intent(in) :: varname
    class(regrid_mesh_base), intent(in), target :: out_mesh
    class(MeshField3D), intent(inout) :: out_field
    type(regrid_nodemap), intent(in), target :: nodeMap_list(:)

    class(LocalMesh3D), pointer :: lcmesh
    integer :: n

    integer :: domID
    integer :: in_mesh_num
    integer :: in_ii
    integer :: in_rank
    type(FILE_base_meshfield), pointer :: in_file_ptr

    class(MeshBase3D), pointer :: ptr_outmesh3D
    class(MeshBase3D), pointer :: ptr_inmesh3D
    !-------------------------------------------

    call PROF_rapstart('regrid_field_interpolate_3D', 0)

    ptr_outmesh3D => out_mesh%ptr_mesh3D

    if (.not. is_cached_in_files) then
      allocate( in_file_list(ptr_outmesh3D%LOCAL_MESH_NUM)  )
      do domID=1, ptr_outmesh3D%LOCAL_MESH_NUM
        in_mesh_num = size(nodeMap_list(domID)%in_mesh_list)
        allocate( in_file_list(domID)%in_files(in_mesh_num) )
        do in_ii=1, in_mesh_num
          in_file_ptr => in_file_list(domID)%in_files(in_ii)
          ptr_inmesh3D => nodeMap_list(domID)%in_mesh_list(in_ii)%ptr_mesh3D
          in_rank = ptr_inmesh3D%lcmesh_list(1)%PRC_myrank

          LOG_INFO("regrid_interp_field",'(a,i4,a,a,i6)') 'domID=', domID, ', Open in_file:', trim(in_basename), in_rank

          select type( ptr_inmesh3D )
          class is (MeshCubeDom3D)
            call in_file_ptr%Init( out_var_num, mesh3D=ptr_inmesh3D )
          class is ( MeshCubedSphereDom3D )
            call in_file_ptr%Init( out_var_num, meshCubedSphere3D=ptr_inmesh3D )
          end select
          call in_file_ptr%Open( in_basename, in_rank )
        end do
      end do
      is_cached_in_files = .true.
    end if

    do n=1, ptr_outmesh3D%LOCAL_MESH_NUM
      lcmesh => ptr_outmesh3D%lcmesh_list(n)
      call interpolate_local_3D( out_field%local(n)%val(:,:),                                &
        n, istep, varname, lcmesh, lcmesh%refElem3D, nodeMap_list(n)%in_mesh_list(1)%elem3D, &
        nodeMap_list(n)  )
    end do

    call PROF_rapend('regrid_field_interpolate_3D', 0)

    return
  end subroutine regrid_field_Interpolate_3D

!------------------------

!OCL SERIAL
  subroutine interpolate_local_2D( out_val, &
      out_domID, istep, varname,            &
      out_lcmesh, out_elem2D, in_elem2D,    &
      mappingInfo )

    use scale_const, only: &
      UNDEF8 => CONST_UNDEF8, &
      PI => CONST_PI
    use scale_file_base_meshfield, only: &
      FILE_base_meshfield
    use scale_mesh_base2d, only: &
      MF2D_XYT => MeshBase2D_DIMTYPEID_XYT
    use scale_polynominal, only: &
      Polynominal_GenLegendrePoly_sub

    implicit none
    class(LocalMesh2D), intent(in) :: out_lcmesh
    class(ElementBase2D), intent(in) :: out_elem2D
    class(ElementBase2D), intent(in) :: in_elem2D
    real(DP), intent(out) :: out_val(out_elem2D%Np,out_lcmesh%NeA)
    integer, intent(in) :: out_domID
    integer, intent(in) :: istep
    character(*), intent(in) :: varname
    type(regrid_nodemap), intent(in), target :: mappingInfo

    integer :: nprc_local
    integer :: ke2D, ke_h
    integer :: p, p_h, pX, pY
    integer :: p1, p2, l
    integer :: in_ke2D, in_ex, in_ey, in_domID, in_prc
    integer :: in_rank
    type(LocalMesh2D), pointer :: in_lcmesh
    integer :: n
    integer :: ii, jj

    type(in_local_val), allocatable :: in_val_list(:)
    class(MeshBase2D), pointer :: in_mesh
    real(RP) :: P1D_ori_x(1,in_elem2D%Nfp)
    real(RP) :: P1D_ori_y(1,in_elem2D%Nfp)
    real(RP) :: ox(1), oy(1)
    real(RP) :: vx(in_elem2D%Nv), vy(in_elem2D%Nv)
    integer :: node_ids(in_elem2D%Nv)
    type(MeshField2D) :: tmp_field2D

    !---------------------------------------------

    nprc_local = size(mappingInfo%in_mesh_list)

    allocate( in_val_list(nprc_local) )

    do in_prc=1, nprc_local
      in_rank = in_prc - 1
      in_mesh => mappingInfo%in_mesh_list(in_prc)%ptr_mesh2D

      n = 1
      in_lcmesh => in_mesh%lcmesh_list(n)
      allocate( in_val_list(in_prc)%spectral_coef2D(in_lcmesh%refElem2D%Np,in_lcmesh%NeX,in_lcmesh%NeY,in_mesh%LOCAL_MESH_NUM) )

      call tmp_field2D%Init( varname, "", in_mesh )
      call in_file_list(out_domID)%in_files(in_prc)%Read_Var( &
        MF2D_XYT, varname, tmp_field2D, step=istep            )

      !$omp parallel do collapse(2) private(in_ke2D,ii,jj) 
      do in_domID=1, in_mesh%LOCAL_MESH_NUM
      do jj=1, in_lcmesh%NeY
      do ii=1, in_lcmesh%NeX
        in_ke2D = ii + (jj-1)*in_lcmesh%NeX
        in_val_list(in_prc)%spectral_coef2D(:,ii,jj,in_domID) = &
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
   !$omp p1, p2, l                                    )    
    do jj=1, out_lcmesh%NeY
    do ii=1, out_lcmesh%NeX
      ke_h = ii + (jj-1)*out_lcmesh%NeX
      ke2D = ke_h

      do pY=1, out_elem2D%Nfp
      do pX=1, out_elem2D%Nfp
        p_h = pX + (pY-1)*out_elem2D%Nfp
        p = p_h
        out_val(p,ke2D) = UNDEF8

        in_domID = mappingInfo%local_domID(p_h,ke_h)
        in_prc = mappingInfo%lcprc(p_h,ke_h)

        if (in_domID > 0 .and. in_prc > 0) then
          in_ex = mappingInfo%elem_i(p_h,ke_h)
          in_ey = mappingInfo%elem_j(p_h,ke_h)  
          
          if ( in_ex > 0 .and. in_ey  > 0 ) then
            in_lcmesh => mappingInfo%in_mesh_list(in_prc)%ptr_mesh2D%lcmesh_list(in_domID)
            in_ke2D = in_ex + (in_ey - 1)*in_lcmesh%NeX

            node_ids(:) = in_lcmesh%EToV(in_ke2D,:)
            vx(:) = in_lcmesh%pos_ev(node_ids(:),1)
            vy(:) = in_lcmesh%pos_ev(node_ids(:),2)

            ox(1) = - 1.0_RP + 2.0_RP * (mappingInfo%elem_x(p,ke2D) - vx(1)) / (vx(2) - vx(1))
            oy(1) = - 1.0_RP + 2.0_RP * (mappingInfo%elem_y(p,ke2D) - vy(1)) / (vy(3) - vy(1))

            call Polynominal_GenLegendrePoly_sub( in_elem2D%PolyOrder, ox, P1D_ori_x(:,:) )
            call Polynominal_GenLegendrePoly_sub( in_elem2D%PolyOrder, oy, P1D_ori_y(:,:) )

            out_val(p,ke2D) = 0.0_RP
            do p2=1, in_elem2D%Nfp
            do p1=1, in_elem2D%Nfp
              l = p1 + (p2-1)*in_elem2D%Nfp
              out_val(p,ke2D) = out_val(p,ke2D) + &
                  ( P1D_ori_x(1,p1) * P1D_ori_y(1,p2)  )                     &
                * sqrt((dble(p1-1) + 0.5_RP)*(dble(p2-1) + 0.5_RP))          &
                * in_val_list(in_prc)%spectral_coef2D(l,in_ex,in_ey,in_domID)
            end do
            end do
          end if
        end if

      end do
      end do
    end do
    end do

    do in_prc=1, nprc_local
     deallocate( in_val_list(in_prc)%spectral_coef2D )
    end do

    return
  end subroutine interpolate_local_2D


!OCL SERIAL
  subroutine interpolate_local_3D( out_val, &
    out_domID, istep, varname,              &
    out_lcmesh, out_elem3D, in_elem3D,      &
    mappingInfo )

    use scale_const, only: &
      UNDEF8 => CONST_UNDEF8, &
      PI => CONST_PI
    use scale_file_base_meshfield, only: &
      FILE_base_meshfield
    use scale_mesh_base3d, only: &
      MF3D_XYZT => MeshBase3D_DIMTYPEID_XYZT
    use scale_polynominal, only: &
      Polynominal_GenLegendrePoly_sub
    implicit none
    
    class(LocalMesh3D), intent(in) :: out_lcmesh
    class(ElementBase3D), intent(in) :: out_elem3D
    class(ElementBase3D), intent(in) :: in_elem3D
    real(DP), intent(out) :: out_val(out_elem3D%Np,out_lcmesh%NeA)
    integer, intent(in) :: out_domID
    integer, intent(in) :: istep
    character(*), intent(in) :: varname
    type(regrid_nodemap), intent(in), target :: mappingInfo

    integer :: nprc_local
    integer :: ke3D, ke_h
    integer :: p, p_h, pX, pY, pZ
    integer :: p1, p2, p3, l
    integer :: in_ke3D, in_ex, in_ey, in_ez, in_domID, in_prc
    integer :: in_rank
    type(LocalMesh3D), pointer :: in_lcmesh
    integer :: n
    integer :: ii, jj, kk

    type(in_local_val), allocatable :: in_val_list(:)
    class(MeshBase3D), pointer :: in_mesh
    real(RP) :: P1D_ori_x(1,in_elem3D%Nnode_h1D)
    real(RP) :: P1D_ori_y(1,in_elem3D%Nnode_h1D)
    real(RP) :: P1D_ori_z(1,in_elem3D%Nnode_v)
    real(RP) :: ox(1), oy(1), oz(1)
    real(RP) :: vx(in_elem3D%Nv), vy(in_elem3D%Nv), vz(in_elem3D%Nv)
    integer :: node_ids(in_elem3D%Nv)
    type(MeshField3D) :: tmp_field3D
    !---------------------------------------------

    nprc_local = size(mappingInfo%in_mesh_list)

    allocate( in_val_list(nprc_local) )

    do in_prc=1, nprc_local
      in_rank = in_prc - 1
      in_mesh => mappingInfo%in_mesh_list(in_prc)%ptr_mesh3D

      n = 1
      in_lcmesh => in_mesh%lcmesh_list(n)
      allocate( in_val_list(in_prc)%spectral_coef3D(in_lcmesh%refElem3D%Np,in_lcmesh%NeX,in_lcmesh%NeY,in_lcmesh%NeZ,in_mesh%LOCAL_MESH_NUM) )

      call tmp_field3D%Init( varname, "", in_mesh )
      call in_file_list(out_domID)%in_files(in_prc)%Read_Var( &
        MF3D_XYZT, varname, tmp_field3D, step=istep           )

      !$omp parallel do collapse(2) private(in_ke3D,ii,jj,kk) 
      do in_domID=1, in_mesh%LOCAL_MESH_NUM
      do kk=1, in_lcmesh%NeZ
      do jj=1, in_lcmesh%NeY
      do ii=1, in_lcmesh%NeX
        in_ke3D = ii + (jj-1)*in_lcmesh%NeX + (kk-1)*in_lcmesh%NeX*in_lcmesh%NeY
        in_val_list(in_prc)%spectral_coef3D(:,ii,jj,kk,in_domID) = &
          matmul( in_elem3D%invV, tmp_field3D%local(in_domID)%val(:,in_ke3D) )
      end do
      end do
      end do
      end do

      call tmp_field3D%Final()
    end do
    
    !$omp parallel do collapse(3) private( ke_h, ke3D, &
    !$omp pX, pY, pZ, p_h, p, in_domID, in_prc,        &
    !$omp in_ex, in_ey, in_ez, in_ke3D,                &
    !$omp in_lcmesh, node_ids, vx, vy, vz, ox, oy, oz, &
    !$omp P1D_ori_x, P1D_ori_y, P1D_ori_z,             &
    !$omp p1, p2, p3, l                                )    
    do kk=1, out_lcmesh%NeZ
    do jj=1, out_lcmesh%NeY
    do ii=1, out_lcmesh%NeX
      ke_h = ii + (jj-1)*out_lcmesh%NeX
      ke3D = ke_h + (kk-1)*out_lcmesh%NeX*out_lcmesh%NeY

      do pZ=1, out_elem3D%Nnode_v
      do pY=1, out_elem3D%Nnode_h1D
      do pX=1, out_elem3D%Nnode_h1D
        p_h = pX + (pY-1)*out_elem3D%Nnode_h1D
        p = p_h + (pZ-1)*out_elem3D%Nnode_h1D**2
        out_val(p,ke3D) = UNDEF8

        in_domID = mappingInfo%local_domID(p_h,ke_h)
        in_prc = mappingInfo%lcprc(p_h,ke_h)

        if (in_domID > 0 .and. in_prc > 0) then
          in_ex = mappingInfo%elem_i(p_h,ke_h)
          in_ey = mappingInfo%elem_j(p_h,ke_h)  
          in_ez = mappingInfo%elem_k(p  ,ke3D)  
          
          if ( in_ex > 0 .and. in_ey  > 0 .and. in_ez  > 0 ) then
            in_lcmesh => mappingInfo%in_mesh_list(in_prc)%ptr_mesh3D%lcmesh_list(in_domID)
            in_ke3D = in_ex + (in_ey - 1)*in_lcmesh%NeX      &
                    + (in_ez - 1)*in_lcmesh%NeX*in_lcmesh%NeY

            node_ids(:) = in_lcmesh%EToV(in_ke3D,:)
            vx(:) = in_lcmesh%pos_ev(node_ids(:),1)
            vy(:) = in_lcmesh%pos_ev(node_ids(:),2)
            vz(:) = in_lcmesh%pos_ev(node_ids(:),3)

            ox(1) = - 1.0_RP + 2.0_RP * (mappingInfo%elem_x(p_h,ke_h) - vx(1)) / (vx(2) - vx(1))
            oy(1) = - 1.0_RP + 2.0_RP * (mappingInfo%elem_y(p_h,ke_h) - vy(1)) / (vy(3) - vy(1))
            oz(1) = - 1.0_RP + 2.0_RP * (mappingInfo%elem_z(p  ,ke3D) - vz(1)) / (vz(5) - vz(1))

            call Polynominal_GenLegendrePoly_sub( in_elem3D%PolyOrder_h, ox, P1D_ori_x(:,:) )
            call Polynominal_GenLegendrePoly_sub( in_elem3D%PolyOrder_h, oy, P1D_ori_y(:,:) )
            call Polynominal_GenLegendrePoly_sub( in_elem3D%PolyOrder_v, oz, P1D_ori_z(:,:) )

            out_val(p,ke3D) = 0.0_RP
            do p3=1, in_elem3D%Nnode_v
            do p2=1, in_elem3D%Nnode_h1D
            do p1=1, in_elem3D%Nnode_h1D
              l = p1 + (p2-1)*in_elem3D%Nnode_h1D + (p3-1)*in_elem3D%Nnode_h1D**2
              out_val(p,ke3D) = out_val(p,ke3D) + &
                  ( P1D_ori_x(1,p1) * P1D_ori_y(1,p2) * P1D_ori_z(1,p3)  )                  &
                * sqrt( (dble(p1-1) + 0.5_RP)*(dble(p2-1) + 0.5_RP)*(dble(p3-1) + 0.5_RP) ) &
                * in_val_list(in_prc)%spectral_coef3D(l,in_ex,in_ey,in_ez,in_domID)
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

    do in_prc=1, nprc_local
      deallocate( in_val_list(in_prc)%spectral_coef3D )
    end do

    return
  end subroutine interpolate_local_3D

end module mod_regrid_interp_field