!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_file_common_meshfield
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: elementbase1D, elementbase2D, elementbase3D
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D

  use scale_meshfield_base, only: MeshField1D, MeshField2D, MeshField3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !

  public :: File_common_meshfield_get_dims1D
  public :: File_common_meshfield_get_dims2D
  public :: File_common_meshfield_get_dims3D
  interface File_common_meshfield_get_dims
    module procedure File_common_meshfield_get_dims1D
    module procedure File_common_meshfield_get_dims2D
    module procedure File_common_meshfield_get_dims3D
  end interface
  public :: FILE_common_meshfield_get_dims


  public :: File_common_meshfield_get_axis1D
  public :: File_common_meshfield_get_axis2D
  public :: File_common_meshfield_get_axis3D
  interface File_common_meshfield_get_axis
    module procedure File_common_meshfield_get_axis1D
    module procedure File_common_meshfield_get_axis2D
    module procedure File_common_meshfield_get_axis3D
  end interface
  public :: FILE_common_meshfield_get_axis

  public :: File_common_meshfield_put_field1D_cartesbuf
  public :: File_common_meshfield_put_field2D_cartesbuf
  public :: File_common_meshfield_put_field3D_cartesbuf

  public :: File_common_meshfield_set_cartesbuf_field3D
  public :: File_common_meshfield_set_cartesbuf_field3D_local

  type, public :: FILE_common_meshfield_diminfo
    character(len=H_SHORT) :: type
    integer :: ndim
    character(len=H_SHORT) :: dims(3)
    character(len=H_SHORT) :: name
    character(len=H_MID) :: desc
    character(len=H_SHORT) :: unit
    integer :: count(3)
    integer :: size
  end type FILE_common_meshfield_diminfo

  public :: File_common_meshfield_get_dtype
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------

  integer, public :: FILE_COMMON_MESHFILED1D_DIMTYPE_NUM   = 2
  integer, public :: FILE_COMMON_MESHFILED1D_DIMTYPEID_X   = 1
  integer, public :: FILE_COMMON_MESHFILED1D_DIMTYPEID_XT  = 2

  integer, public :: FILE_COMMON_MESHFILED2D_DIMTYPE_NUM   = 4
  integer, public :: FILE_COMMON_MESHFILED2D_DIMTYPEID_X   = 1
  integer, public :: FILE_COMMON_MESHFILED2D_DIMTYPEID_Y   = 2
  integer, public :: FILE_COMMON_MESHFILED2D_DIMTYPEID_XY  = 3
  integer, public :: FILE_COMMON_MESHFILED2D_DIMTYPEID_XYT = 4

  integer, public :: FILE_COMMON_MESHFILED3D_DIMTYPE_NUM    = 5
  integer, public :: FILE_COMMON_MESHFILED3D_DIMTYPEID_X    = 1
  integer, public :: FILE_COMMON_MESHFILED3D_DIMTYPEID_Y    = 2
  integer, public :: FILE_COMMON_MESHFILED3D_DIMTYPEID_Z    = 3
  integer, public :: FILE_COMMON_MESHFILED3D_DIMTYPEID_XYZ  = 4
  integer, public :: FILE_COMMON_MESHFILED3D_DIMTYPEID_XYZT = 5

  !--------------------
  !
  !++ Private procedures
  !
  !-------------------


contains

  !- 1D ---------------

  subroutine File_common_meshfield_get_dims1D( mesh1D, dimsinfo )
    implicit none

    class(MeshBase1D), target, intent(in) :: mesh1D
    type(FILE_common_meshfield_diminfo), intent(out) :: dimsinfo(FILE_COMMON_MESHFILED1D_DIMTYPE_NUM)

    integer :: i_size
    !-------------------------------------------------

    i_size = mesh1D%NeG * mesh1D%refElem1D%Np
    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED1D_DIMTYPEID_X),  &
      "x", "X-coordinate", "X", 1, (/ "x" /), (/ i_size /)              )
    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED1D_DIMTYPEID_XT), &
      "xt", "X-coordinate", "XT", 1, (/ "x" /), (/ i_size  /)           )

    return
  end subroutine File_common_meshfield_get_dims1D

  subroutine File_common_meshfield_get_axis1D( mesh1D, dimsinfo, x )
    implicit none

    class(MeshBase1D), target, intent(in) :: mesh1D
    type(FILE_common_meshfield_diminfo), intent(in) :: dimsinfo(FILE_COMMON_MESHFILED1D_DIMTYPE_NUM)
    real(DP), intent(out) :: x(dimsinfo(FILE_COMMON_MESHFILED1D_DIMTYPEID_X)%size)

    integer :: n
    integer :: i
    integer :: i2
    type(ElementBase1D), pointer :: refElem
    type(LocalMesh1D), pointer :: lcmesh
    !-------------------------------------------------
    
    do n=1,mesh1D%LOCAL_MESH_NUM
      lcmesh => mesh1D%lcmesh_list(n)
      refElem => lcmesh%refElem1D

      do i=1,mesh1D%NeG
      do i2=1, refElem%Np
        x(i2 + (i-1)*refElem%Np + (n-1)*refElem%Np*lcmesh%Ne) = mesh1D%lcmesh_list(n)%pos_en(i2,i,1)
      end do
      end do
    end do
    
    return
  end subroutine File_common_meshfield_get_axis1D

  subroutine File_common_meshfield_put_field1D_cartesbuf( mesh1D, field1D, &
    buf )
    implicit none
    class(MeshBase1D), target, intent(in) :: mesh1D
    class(MeshField1D), intent(in) :: field1d
    real(RP), intent(inout) :: buf(:)

    integer :: n, k, p
    type(LocalMesh1D), pointer :: lcmesh
    integer :: bufsize, ptr
    !------------------------------------------------

    ptr = 0
    do n=1, mesh1D%LOCAL_MESH_NUM
      lcmesh => mesh1D%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE
      do p=1, lcmesh%refElem%Np
        ptr = ptr + 1
        buf(ptr) =  field1d%local(n)%val(p,k)
      end do
      end do
    end do

    return
  end subroutine File_common_meshfield_put_field1D_cartesbuf

  !- 2D ---------------

  subroutine File_common_meshfield_get_dims2D( mesh2D, dimsinfo )
    implicit none

    class(MeshRectDom2D), target, intent(in) :: mesh2D
    type(FILE_common_meshfield_diminfo), intent(out) :: dimsinfo(FILE_COMMON_MESHFILED2D_DIMTYPE_NUM)

    type(ElementBase2D), pointer :: refElem
    type(LocalMesh2D), pointer :: lcmesh
    integer :: i, j, n

    integer :: i_size, j_size
    !-------------------------------------------------
    
    i_size = 0
    do i=1, size(mesh2D%rcdomIJ2LCMeshID,1)
      n = mesh2D%rcdomIJ2LCMeshID(i,1)
      lcmesh => mesh2D%lcmesh_list(n)
      i_size =i_size + lcmesh%NeX * lcmesh%refElem2D%Nfp
    end do

    j_size = 0
    do j=1, size(mesh2D%rcdomIJ2LCMeshID,2)
      n = mesh2D%rcdomIJ2LCMeshID(1,j)
      lcmesh => mesh2D%lcmesh_list(n)
      j_size = j_size + lcmesh%NeY * lcmesh%refElem2D%Nfp
    end do
  
    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED2D_DIMTYPEID_X),  &
      "x", "X-coordinate", "X", 1, (/ "x" /), (/ i_size /)              )
    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED2D_DIMTYPEID_Y),  &
      "y", "Y-coordinate", "Y", 1, (/ "y" /), (/ j_size /)              )

    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED2D_DIMTYPEID_XY),    &
      "xy", "XY-coordinate", "XY", 2, (/ "x", "y" /), (/ i_size, j_size /) )
    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED2D_DIMTYPEID_XYT),  &
      "xyt", "XY-coordinate", "XYT", 2, (/ "x", "y" /), (/ i_size, j_size /) )

    return
  end subroutine File_common_meshfield_get_dims2D

  subroutine File_common_meshfield_get_axis2D( mesh2D, dimsinfo, x, y  )
    implicit none

    class(MeshRectDom2D), target, intent(in) :: mesh2D  
    type(FILE_common_meshfield_diminfo), intent(in) :: dimsinfo(FILE_COMMON_MESHFILED2D_DIMTYPE_NUM)
    real(DP), intent(out) :: x(dimsinfo(FILE_COMMON_MESHFILED2D_DIMTYPEID_X)%size)
    real(DP), intent(out) :: y(dimsinfo(FILE_COMMON_MESHFILED2D_DIMTYPEID_Y)%size)

    integer :: n
    integer :: k
    integer :: i, j 
    integer :: i2, j2
    type(ElementBase2D), pointer :: refElem
    type(LocalMesh2D), pointer :: lcmesh

    integer :: is, js, ie, je, igs, jgs
    !-------------------------------------------------
    

    igs = 0; jgs = 0
    do n=1 ,mesh2D%LOCAL_MESH_NUM
      lcmesh => mesh2D%lcmesh_list(n)
      refElem => lcmesh%refElem2D

      do j=1, lcmesh%NeY
      do i=1, lcmesh%NeX
        k = i + (j-1) * lcmesh%NeX
        if ( j==1 ) then
          is = igs + 1 + (i-1)*refElem%Nfp
          ie = is + refElem%Nfp - 1
          x(is:ie) = lcmesh%pos_en(refElem%Fmask(:,1),k,1)
        end if
        if ( i==1 ) then
          js = jgs + 1 + (j-1)*refElem%Nfp
          je = js + refElem%Nfp - 1
          y(js:je) = lcmesh%pos_en(refElem%Fmask(:,4),k,2)
        end if
      end do
      end do

      igs = ie; jgs = je
    end do

    return
  end subroutine File_common_meshfield_get_axis2D

  subroutine File_common_meshfield_put_field2D_cartesbuf( mesh2D, field2D, &
    buf )
    implicit none
    class(MeshRectDom2D), target, intent(in) :: mesh2D
    class(MeshField2D), intent(in) :: field2d
    real(RP), intent(inout) :: buf(:,:)

    integer :: n, k1, p
    integer :: i0, j0, i1, j1, i2, j2, i, j
    type(LocalMesh2D), pointer :: lcmesh
    type(elementbase2D), pointer :: refElem
    integer :: i0_s, j0_s
    !------------------------------------------------

    i0_s = 0; j0_s = 0
    
    do j0=1, size(mesh2D%rcdomIJ2LCMeshID,2)
    do i0=1, size(mesh2D%rcdomIJ2LCMeshID,1)
      n =  mesh2D%rcdomIJ2LCMeshID(i0,j0)

      lcmesh => mesh2D%lcmesh_list(n)
      refElem => lcmesh%refElem2D
        
      do j1=1, lcmesh%NeY
      do i1=1, lcmesh%NeX
        k1 = i1 + (j1-1)*lcmesh%NeX
        do j2=1, refElem%Nfp
        do i2=1, refElem%Nfp
          i = i0_s + i2 + (i1-1)*refElem%Nfp 
          j = j0_s + j2 + (j1-1)*refElem%Nfp
          buf(i,j) = field2d%local(n)%val(i2+(j2-1)*refElem%Nfp,k1)
        end do
        end do
      end do
      end do

      i0_s = i0_s + lcmesh%NeX * refElem%Nfp
      j0_s = j0_s + lcmesh%NeY * refElem%Nfp
    end do
    end do

    return
  end subroutine File_common_meshfield_put_field2D_cartesbuf

  !- 3D ------------

  subroutine File_common_meshfield_get_dims3D( mesh3D, dimsinfo )
    implicit none

    class(MeshCubeDom3D), target, intent(in) :: mesh3D
    type(FILE_common_meshfield_diminfo), intent(out) :: dimsinfo(FILE_COMMON_MESHFILED3D_DIMTYPE_NUM)

    type(LocalMesh3D), pointer :: lcmesh
    integer :: i, j, k, n
    integer :: icount, jcount, kcount

    integer :: i_size, j_size, k_size
    !-------------------------------------------------
    
    i_size = 0
    do i=1, size(mesh3D%rcdomIJK2LCMeshID,1)
      n = mesh3D%rcdomIJK2LCMeshID(i,1,1)
      lcmesh => mesh3D%lcmesh_list(n)
      i_size = i_size + lcmesh%NeX * lcmesh%refElem3D%Nnode_h1D
    end do

    j_size = 0
    do j=1, size(mesh3D%rcdomIJK2LCMeshID,2)
      n = mesh3D%rcdomIJK2LCMeshID(1,j,1)
      lcmesh => mesh3D%lcmesh_list(n)
      j_size = j_size + lcmesh%NeY * lcmesh%refElem3D%Nnode_h1D
    end do

    k_size = 0
    do k=1, size(mesh3D%rcdomIJK2LCMeshID,3)
      n = mesh3D%rcdomIJK2LCMeshID(1,1,k)
      lcmesh => mesh3D%lcmesh_list(n)
      k_size = k_size + lcmesh%NeZ * lcmesh%refElem3D%Nnode_v
    end do  

    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED3D_DIMTYPEID_X),  &
      "x", "X-coordinate", "X", 1, (/ "x" /), (/ i_size /)              )
    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED3D_DIMTYPEID_Y),  &
      "y", "Y-coordinate", "Y", 1, (/ "y" /), (/ j_size /)              )
    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED3D_DIMTYPEID_Z),  &
      "z", "Z-coordinate", "Z", 1, (/ "z" /), (/ k_size /)              )

    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED3D_DIMTYPEID_XYZ),                   &
      "xyz", "XYZ-coordinate", "XYZ", 3, (/ "x", "y", "z" /), (/ i_size, j_size, k_size /) )
    call set_dimension( dimsinfo(FILE_COMMON_MESHFILED3D_DIMTYPEID_XYZT),                &
      "xyzt", "XYZ-coordinate", "XYZT", 3, (/ "x", "y", "z" /), (/ i_size, j_size, k_size /) )

    return
  end subroutine File_common_meshfield_get_dims3D

  subroutine File_common_meshfield_get_axis3D( mesh3D, dimsinfo, x, y, z )
    implicit none

    class(MeshCubeDom3D), target, intent(in) :: mesh3D  
    type(FILE_common_meshfield_diminfo), intent(in) :: dimsinfo(FILE_COMMON_MESHFILED3D_DIMTYPE_NUM)
    real(DP), intent(out) :: x(dimsinfo(FILE_COMMON_MESHFILED3D_DIMTYPEID_X)%size)
    real(DP), intent(out) :: y(dimsinfo(FILE_COMMON_MESHFILED3D_DIMTYPEID_Y)%size)
    real(DP), intent(out) :: z(dimsinfo(FILE_COMMON_MESHFILED3D_DIMTYPEID_Z)%size)

    integer :: n, kelem
    integer :: i, j, k
    integer :: i2, j2, k2
    type(ElementBase3D), pointer :: refElem
    type(LocalMesh3D), pointer :: lcmesh

    integer :: is, js, ks, ie, je, ke, igs, jgs, kgs
    integer :: Nnode_h1D

    !--------------------
    
    igs = 0; jgs = 0; kgs = 0

    do n=1 ,mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      refElem => lcmesh%refElem3D
      Nnode_h1D = refElem%Nnode_h1D

      do k=1, lcmesh%NeZ
      do j=1, lcmesh%NeY
      do i=1, lcmesh%NeX
        kelem = i + (j-1)*lcmesh%NeX + (k-1)*lcmesh%NeX*lcmesh%NeY
        if ( j==1 .and. k==1) then
          is = igs + 1 + (i-1)*Nnode_h1D
          ie = is + Nnode_h1D - 1
          x(is:ie) = lcmesh%pos_en(refElem%Fmask_h(1:Nnode_h1D,1),kelem,1)
        end if
        if ( i==1 .and. k==1) then
          js = jgs + 1 + (j-1)*Nnode_h1D
          je = js + Nnode_h1D - 1
          y(js:je) = lcmesh%pos_en(refElem%Fmask_h(1:Nnode_h1D,4),kelem,2)
        end if
        if ( i==1 .and. j==1) then
          ks = kgs + 1 + (k-1)*refElem%Nnode_v
          ke = ks + refElem%Nnode_v - 1
          z(ks:ke) = lcmesh%pos_en(refElem%Colmask(:,1),kelem,3)
        end if
      end do
      end do
      end do

      igs = ie; jgs = je; kgs = ke
    end do

    return
  end subroutine File_common_meshfield_get_axis3D

  subroutine File_common_meshfield_put_field3D_cartesbuf( mesh3D, field3D, &
    buf )
    implicit none
    class(MeshCubeDom3D), target, intent(in) :: mesh3D
    class(MeshField3D), intent(in) :: field3d
    real(RP), intent(inout) :: buf(:,:,:)

    integer :: n, kelem1, p
    integer :: i0, j0, k0, i1, j1, k1, i2, j2, k2, i, j, k
    type(LocalMesh3D), pointer :: lcmesh
    type(elementbase3D), pointer :: refElem
    integer :: i0_s, j0_s, k0_s, indx
    !----------------------------------------------------

    i0_s = 0; j0_s = 0; k0_s = 0

    do k0=1, size(mesh3D%rcdomIJK2LCMeshID,3)  
    do j0=1, size(mesh3D%rcdomIJK2LCMeshID,2)
    do i0=1, size(mesh3D%rcdomIJK2LCMeshID,1)
      n =  mesh3D%rcdomIJK2LCMeshID(i0,j0,k0)

      lcmesh => mesh3D%lcmesh_list(n)
      refElem => lcmesh%refElem3D

      do k1=1, lcmesh%NeZ
      do j1=1, lcmesh%NeY
      do i1=1, lcmesh%NeX
        kelem1 = i1 + (j1-1)*lcmesh%NeX + (k1-1)*lcmesh%NeX*lcmesh%NeY
        do k2=1, refElem%Nnode_v
        do j2=1, refElem%Nnode_h1D
        do i2=1, refElem%Nnode_h1D
          i = i0_s + i2 + (i1-1)*refElem%Nnode_h1D
          j = j0_s + j2 + (j1-1)*refElem%Nnode_h1D
          k = k0_s + k2 + (k1-1)*refElem%Nnode_v
          indx = i2 + (j2-1)*refElem%Nnode_h1D + (k2-1)*refElem%Nnode_h1D**2
          buf(i,j,k) = field3d%local(n)%val(indx,kelem1)
        end do
        end do
        end do
      end do
      end do
      end do

      i0_s = i0_s + lcmesh%NeX * refElem%Nnode_h1D
      j0_s = j0_s + lcmesh%NeY * refElem%Nnode_h1D
      k0_s = k0_s + lcmesh%NeZ * refElem%Nnode_v
    end do
    end do
    end do

    return
  end subroutine File_common_meshfield_put_field3D_cartesbuf

  subroutine File_common_meshfield_set_cartesbuf_field3D( mesh3D, buf, &
    field3D )
    implicit none
    class(MeshCubeDom3D), target, intent(in) :: mesh3D
    real(RP), intent(in) :: buf(:,:,:)
    class(MeshField3D), intent(inout) :: field3d

    integer :: n
    integer :: i0, j0, k0
    type(LocalMesh3D), pointer :: lcmesh
    type(elementbase3D), pointer :: refElem
    integer :: i0_s, j0_s, k0_s
    !----------------------------------------------------

    i0_s = 0; j0_s = 0; k0_s = 0

    do k0=1, size(mesh3D%rcdomIJK2LCMeshID,3)  
    do j0=1, size(mesh3D%rcdomIJK2LCMeshID,2)
    do i0=1, size(mesh3D%rcdomIJK2LCMeshID,1)
      n = mesh3D%rcdomIJK2LCMeshID(i0,j0,k0)
      lcmesh => mesh3D%lcmesh_list(n)
      refElem => lcmesh%refElem3D

      call File_common_meshfield_set_cartesbuf_field3D_local(  &
        lcmesh, buf(:,:,:), i0_s, j0_s, k0_s,                  &
        field3d%local(n)%val(:,:)                              )

      i0_s = i0_s + lcmesh%NeX * refElem%Nnode_h1D
      j0_s = j0_s + lcmesh%NeY * refElem%Nnode_h1D
      k0_s = k0_s + lcmesh%NeZ * refElem%Nnode_v
    end do
    end do
    end do

    return
  end subroutine File_common_meshfield_set_cartesbuf_field3D

  subroutine File_common_meshfield_set_cartesbuf_field3D_local( &
    lcmesh, buf, i0_s, j0_s, k0_s,                              &
    val )
    implicit none
    type(LocalMesh3D), intent(in) :: lcmesh
    real(RP), intent(in) :: buf(:,:,:)
    integer, intent(in) :: i0_s, j0_s, k0_s
    real(RP), intent(inout) :: val(lcmesh%refElem3D%Np,lcmesh%NeA)

    integer :: n, kelem1, p
    integer :: i1, j1, k1, i2, j2, k2, i, j, k
    type(elementbase3D), pointer :: refElem
    integer :: indx
    !----------------------------------------------------

    refElem => lcmesh%refElem3D

    do k1=1, lcmesh%NeZ
    do j1=1, lcmesh%NeY
    do i1=1, lcmesh%NeX
      kelem1 = i1 + (j1-1)*lcmesh%NeX + (k1-1)*lcmesh%NeX*lcmesh%NeY
      do k2=1, refElem%Nnode_v
      do j2=1, refElem%Nnode_h1D
      do i2=1, refElem%Nnode_h1D
        i = i0_s + i2 + (i1-1)*refElem%Nnode_h1D
        j = j0_s + j2 + (j1-1)*refElem%Nnode_h1D
        k = k0_s + k2 + (k1-1)*refElem%Nnode_v
        indx = i2 + (j2-1)*refElem%Nnode_h1D + (k2-1)*refElem%Nnode_h1D**2
        val(indx,kelem1) = buf(i,j,k)
      end do
      end do
      end do
    end do
    end do
    end do

    return
  end subroutine File_common_meshfield_set_cartesbuf_field3D_local

  function File_common_meshfield_get_dtype( datatype ) result( dtype )

    use scale_file_h, only: &
      FILE_REAL8, FILE_REAL4
    implicit none

    character(*), intent(in) :: datatype
    integer :: dtype
    !--------------------------

    ! dtype is used to define the data type of axis variables in file
    if    ( datatype == 'REAL8' ) then
      dtype = FILE_REAL8
    elseif( datatype == 'REAL4' ) then
        dtype = FILE_REAL4
    else
      if    ( RP == 8 ) then
        dtype = FILE_REAL8
      elseif( RP == 4 ) then
        dtype = FILE_REAL4
      else
        LOG_ERROR("file_restart_meshfield_get_dtype",*) 'unsupported data type. Check!', trim(datatype)
        call PRC_abort
      endif
    endif

    return
  end function File_common_meshfield_get_dtype

  !- private -----------------------------------------------------------------------
  subroutine set_dimension( dim, name, desc, dim_type, ndims, dims, count )
    implicit none

    type(FILE_common_meshfield_diminfo), intent(out) :: dim
    character(*), intent(in) :: name
    character(*), intent(in) :: desc
    character(*), intent(in) :: dim_type
    integer, intent(in) :: ndims
    character(len=*), intent(in) :: dims(ndims)
    integer, intent(in) :: count(ndims)

    integer :: d
    !----------------------------------------------------

    dim%name = name
    dim%unit = "m"
    dim%desc = desc
    dim%type = dim_type
    dim%ndim = ndims
    dim%size = 1
    do d=1, ndims
      dim%dims(d) = dims(d)
      dim%count(d) = count(d)
      dim%size = dim%size * count(d)
    end do

    return
  end subroutine set_dimension

end module scale_file_common_meshfield
