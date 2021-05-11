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
  use scale_mesh_base, only: MeshDimInfo
  use scale_mesh_base1d, only: MeshBase1D, &
    MeshBase1D_DIMTYPEID_X, MeshBase1D_DIMTYPEID_XT, &
    MeshBase1D_DIMTYPE_NUM 
  use scale_mesh_base2d, only: MeshBase2D, &
    MeshBase2D_DIMTYPEID_X, MeshBase2D_DIMTYPEID_Y,    &
    MeshBase2D_DIMTYPEID_XY, MeshBase2D_DIMTYPEID_XYT, &
    MeshBase2D_DIMTYPE_NUM 
  use scale_mesh_base3d, only: MeshBase3D, &
    MeshBase3D_DIMTYPEID_X, MeshBase3D_DIMTYPEID_Y, MeshBase3D_DIMTYPEID_Z,       &
    MeshBase3D_DIMTYPEID_ZT, MeshBase3D_DIMTYPEID_XYZ, MeshBase3D_DIMTYPEID_XYZT, &
    MeshBase3D_DIMTYPE_NUM

  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
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
    module procedure File_common_meshfield_get_dims2D_cubedsphere
    module procedure File_common_meshfield_get_dims3D
  end interface
  public :: FILE_common_meshfield_get_dims


  public :: File_common_meshfield_get_axis1D
  public :: File_common_meshfield_get_axis2D
  public :: File_common_meshfield_get_axis3D
  interface File_common_meshfield_get_axis
    module procedure File_common_meshfield_get_axis1D
    module procedure File_common_meshfield_get_axis2D
    module procedure File_common_meshfield_get_axis2D_cubedsphere
    module procedure File_common_meshfield_get_axis3D
  end interface
  public :: FILE_common_meshfield_get_axis

  public :: File_common_meshfield_put_field1D_cartesbuf
  public :: File_common_meshfield_put_field2D_cartesbuf
  public :: File_common_meshfield_put_field2D_cubedsphere_cartesbuf  
  public :: File_common_meshfield_put_field3D_cartesbuf

  public :: File_common_meshfield_set_cartesbuf_field1D
  public :: File_common_meshfield_set_cartesbuf_field2D
  public :: File_common_meshfield_set_cartesbuf_field2D_cubedsphere
  public :: File_common_meshfield_set_cartesbuf_field3D

  public :: File_common_meshfield_set_cartesbuf_field1D_local
  public :: File_common_meshfield_set_cartesbuf_field2D_local
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
    type(FILE_common_meshfield_diminfo), intent(out) :: dimsinfo(MeshBase1D_DIMTYPE_NUM)

    integer :: i_size
    type(MeshDimInfo), pointer :: diminfo
    !-------------------------------------------------

    i_size = mesh1D%NeG * mesh1D%refElem1D%Np

    dimInfo => mesh1D%dimInfo(MeshBase1D_DIMTYPEID_X)
    call set_dimension( dimsinfo(MeshBase1D_DIMTYPEID_X),  &
      dimInfo, "X", 1, (/ dimInfo%name /), (/ i_size /)    )

    dimInfo => mesh1D%dimInfo(MeshBase1D_DIMTYPEID_XT)
    call set_dimension( dimsinfo(MeshBase1D_DIMTYPEID_XT), &
      dimInfo, "XT", 1, (/ dimInfo%name /), (/ i_size  /)  )

    return
  end subroutine File_common_meshfield_get_dims1D

  subroutine File_common_meshfield_get_axis1D( mesh1D, dimsinfo, x, &
    force_uniform_grid )
    implicit none

    class(MeshBase1D), target, intent(in) :: mesh1D
    type(FILE_common_meshfield_diminfo), intent(in) :: dimsinfo(MeshBase1D_DIMTYPE_NUM)
    real(DP), intent(out) :: x(dimsinfo(MeshBase1D_DIMTYPEID_X)%size)
    logical, intent(in), optional :: force_uniform_grid

    integer :: n
    integer :: i
    integer :: is, ie
    type(ElementBase1D), pointer :: refElem
    type(LocalMesh1D), pointer :: lcmesh

    logical :: uniform_grid = .false.
    real(RP), allocatable :: x_local(:)
    !-------------------------------------------------
    
    if ( present(force_uniform_grid) ) uniform_grid = force_uniform_grid

    do n=1,mesh1D%LOCAL_MESH_NUM
      lcmesh => mesh1D%lcmesh_list(n)
      refElem => lcmesh%refElem1D

      allocate( x_local(refElem%Np) )

      do i=1,mesh1D%NeG
        x_local(:) = lcmesh%pos_en(:,i,1)
        if ( uniform_grid ) call get_uniform_grid1D( x_local, refElem%Nfp )

        is = 1 + (i-1)*refElem%Np + (n-1)*refElem%Np*lcmesh%Ne
        ie = is + refElem%Np
        x(is:ie) = x_local(:)
      end do

      deallocate( x_local )
    end do
    
    return
  end subroutine File_common_meshfield_get_axis1D

  subroutine File_common_meshfield_put_field1D_cartesbuf( mesh1D, field1D, &
    buf, force_uniform_grid )
    use scale_polynominal, only: &
      polynominal_genLegendrePoly
    implicit none
    class(MeshBase1D), target, intent(in) :: mesh1D
    class(MeshField1D), intent(in) :: field1d
    real(RP), intent(inout) :: buf(:)
    logical, intent(in), optional :: force_uniform_grid

    integer :: n, kelem1, p
    integer :: i0, i1, i2, i
    type(LocalMesh1D), pointer :: lcmesh
    type(elementbase1D), pointer :: refElem
    integer :: i0_s

    logical :: uniform_grid = .false.
    integer :: Np
    real(RP), allocatable :: x_local(:)
    real(RP) :: x_local0, delx
    real(RP) :: ox
    real(RP), allocatable :: spectral_coef(:)
    real(RP), allocatable :: P1D_ori_x(:,:)
    !------------------------------------------------

    if ( present(force_uniform_grid) ) uniform_grid = force_uniform_grid

    i0_s = 0
    do n=1, mesh1D%LOCAL_MESH_NUM
      lcmesh => mesh1D%lcmesh_list(n)
      refElem => lcmesh%refElem1D
      Np = refElem%Np
      
      if ( uniform_grid ) then
        allocate( x_local(Np) )
        allocate( spectral_coef(Np) )
        allocate( P1D_ori_x(1,Np) )
      end if

      do kelem1=lcmesh%NeS, lcmesh%NeE
        if ( uniform_grid ) then
          x_local(:) = lcmesh%pos_en(:,kelem1,1)
          x_local0 = x_local(1); delx = x_local(Np) - x_local0
          call get_uniform_grid1D( x_local, Np )
  
          spectral_coef(:) = matmul(refElem%invV(:,:), field1d%local(n)%val(:,kelem1))
          do i2=1, Np
            ox = - 1.0_RP + 2.0_RP * (x_local(i2) - x_local0) / delx  
            P1D_ori_x(:,:) = polynominal_genLegendrePoly( refElem%PolyOrder, (/ ox /) )
  
            i = i0_s + i2 + (kelem1-1)*Np 
            buf(i) = 0.0_RP 
            do p=1, Np
              buf(i) = buf(i) + &
                  P1D_ori_x(1,p) * sqrt( dble(p-1) + 0.5_RP ) * spectral_coef(p)
            end do
          end do
        else
          do i2=1, Np
            i = i0_s + i2 + (kelem1-1)*Np 
            buf(i) =  field1d%local(n)%val(i2,kelem1)
          end do
        end if
      end do

      i0_s = i0_s + lcmesh%Ne * refElem%Np
      if ( uniform_grid ) then
        deallocate( x_local )
        deallocate( spectral_coef )
        deallocate( P1D_ori_x )
      end if
    end do

    return
  end subroutine File_common_meshfield_put_field1D_cartesbuf

  subroutine File_common_meshfield_set_cartesbuf_field1D( mesh1D, buf, &
    field1D )
    implicit none
    class(MeshBase1D), target, intent(in) :: mesh1D
    real(RP), intent(in) :: buf(:)
    class(MeshField1D), intent(inout) :: field1d

    integer :: n
    integer :: i0
    type(LocalMesh1D), pointer :: lcmesh
    type(elementbase1D), pointer :: refElem
    integer :: i0_s
    !----------------------------------------------------

    i0_s = 0

    do i0=1, mesh1D%LOCAL_MESH_NUM
      n = i0
      lcmesh => mesh1D%lcmesh_list(n)
      refElem => lcmesh%refElem1D

      call File_common_meshfield_set_cartesbuf_field1D_local(  &
        lcmesh, buf(:), i0_s,                                  &
        field1d%local(n)%val(:,:)                              )

      i0_s = i0_s + lcmesh%Ne * refElem%Np
    end do

    return
  end subroutine File_common_meshfield_set_cartesbuf_field1D

  subroutine File_common_meshfield_set_cartesbuf_field1D_local( &
    lcmesh, buf, i0_s,                                          &
    val )
    implicit none
    type(LocalMesh1D), intent(in) :: lcmesh
    real(RP), intent(in) :: buf(:)
    integer, intent(in) :: i0_s
    real(RP), intent(inout) :: val(lcmesh%refElem1D%Np,lcmesh%NeA)

    integer :: n, kelem1, p
    integer :: i, i1, i2
    type(elementbase1D), pointer :: refElem
    integer :: indx
    !----------------------------------------------------

    refElem => lcmesh%refElem1D

    do i1=1, lcmesh%Ne
      kelem1 = i1
      do i2=1, refElem%Np
        i = i0_s + i2 + (i1-1)*refElem%Np
        indx = i2
        val(indx,kelem1) = buf(i)
      end do
    end do

    return
  end subroutine File_common_meshfield_set_cartesbuf_field1D_local

  !- 2D ---------------

  subroutine File_common_meshfield_get_dims2D( mesh2D, dimsinfo )
    implicit none

    class(MeshRectDom2D), target, intent(in) :: mesh2D
    type(FILE_common_meshfield_diminfo), intent(out) :: dimsinfo(MeshBase2D_DIMTYPE_NUM)

    type(ElementBase2D), pointer :: refElem
    type(LocalMesh2D), pointer :: lcmesh
    integer :: i, j, n

    integer :: i_size, j_size
    type(MeshDimInfo), pointer :: diminfo
    type(MeshDimInfo), pointer :: diminfo_x
    type(MeshDimInfo), pointer :: diminfo_y
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
  
    diminfo_x => mesh2D%dimInfo(MeshBase2D_DIMTYPEID_X)
    call set_dimension( dimsinfo(MeshBase2D_DIMTYPEID_X),   &
      dimInfo_x, "X", 1, (/ diminfo_x%name /), (/ i_size /) )

    diminfo_y => mesh2D%dimInfo(MeshBase2D_DIMTYPEID_Y)
    call set_dimension( dimsinfo(MeshBase2D_DIMTYPEID_Y),   &
      diminfo_y, "Y", 1, (/ diminfo_y%name /), (/ j_size /) )

    diminfo => mesh2D%dimInfo(MeshBase2D_DIMTYPEID_XY)
    call set_dimension( dimsinfo(MeshBase2D_DIMTYPEID_XY),    &
      diminfo, "XY", 2, (/ diminfo_x%name, diminfo_y%name /), &
      (/ i_size, j_size /) )

    diminfo => mesh2D%dimInfo(MeshBase2D_DIMTYPEID_XYT)
    call set_dimension( dimsinfo(MeshBase2D_DIMTYPEID_XYT),    &
      diminfo, "XYT", 2, (/ diminfo_x%name, diminfo_y%name /), &
      (/ i_size, j_size /) )
  
    return
  end subroutine File_common_meshfield_get_dims2D

  subroutine File_common_meshfield_get_dims2D_cubedsphere( mesh2D, dimsinfo )
    implicit none

    class(MeshCubedSphereDom2D), target, intent(in) :: mesh2D
    type(FILE_common_meshfield_diminfo), intent(out) :: dimsinfo(MeshBase2D_DIMTYPE_NUM)

    type(ElementBase2D), pointer :: refElem
    type(LocalMesh2D), pointer :: lcmesh
    integer :: i, j, n

    integer :: i_size, j_size

    type(MeshDimInfo), pointer :: diminfo
    type(MeshDimInfo), pointer :: diminfo_x
    type(MeshDimInfo), pointer :: diminfo_y
    !-------------------------------------------------
    
    i_size = 0
    do i=1, size(mesh2D%rcdomIJP2LCMeshID,1)
      n = mesh2D%rcdomIJP2LCMeshID(i,1,1)
      lcmesh => mesh2D%lcmesh_list(n)
      i_size =i_size + lcmesh%NeX * lcmesh%refElem2D%Nfp
    end do

    j_size = 0
    do j=1, size(mesh2D%rcdomIJP2LCMeshID,2)
      n = mesh2D%rcdomIJP2LCMeshID(1,j,1)
      lcmesh => mesh2D%lcmesh_list(n)
      j_size = j_size + lcmesh%NeY * lcmesh%refElem2D%Nfp
    end do

    j_size = j_size * size(mesh2D%rcdomIJP2LCMeshID,3)
  
    diminfo_x => mesh2D%dimInfo(MeshBase2D_DIMTYPEID_X)
    call set_dimension( dimsinfo(MeshBase2D_DIMTYPEID_X),   &
      dimInfo_x, "X", 1, (/ diminfo_x%name /), (/ i_size /) )

    diminfo_y => mesh2D%dimInfo(MeshBase2D_DIMTYPEID_Y)
    call set_dimension( dimsinfo(MeshBase2D_DIMTYPEID_Y),   &
      diminfo_y, "Y", 1, (/ diminfo_y%name /), (/ j_size /) )

    diminfo => mesh2D%dimInfo(MeshBase2D_DIMTYPEID_XY)
    call set_dimension( dimsinfo(MeshBase2D_DIMTYPEID_XY),    &
      diminfo, "XY", 2, (/ diminfo_x%name, diminfo_y%name /), &
      (/ i_size, j_size /) )

    diminfo => mesh2D%dimInfo(MeshBase2D_DIMTYPEID_XYT)
    call set_dimension( dimsinfo(MeshBase2D_DIMTYPEID_XYT),    &
      diminfo, "XYT", 2, (/ diminfo_x%name, diminfo_y%name /), &
      (/ i_size, j_size /) )

    return
  end subroutine File_common_meshfield_get_dims2D_cubedsphere

  subroutine File_common_meshfield_get_axis2D( mesh2D, dimsinfo, x, y, &
    force_uniform_grid  )
    implicit none

    class(MeshRectDom2D), target, intent(in) :: mesh2D  
    type(FILE_common_meshfield_diminfo), intent(in) :: dimsinfo(MeshBase2D_DIMTYPE_NUM)
    real(DP), intent(out) :: x(dimsinfo(MeshBase2D_DIMTYPEID_X)%size)
    real(DP), intent(out) :: y(dimsinfo(MeshBase2D_DIMTYPEID_Y)%size)
    logical, intent(in), optional :: force_uniform_grid

    integer :: n
    integer :: ni, nj
    integer :: k
    integer :: i, j 
    integer :: i2, j2
    type(ElementBase2D), pointer :: refElem
    type(LocalMesh2D), pointer :: lcmesh

    integer :: is, js, ie, je, igs, jgs

    logical :: uniform_grid = .false.
    real(RP), allocatable :: x_local(:)
    real(RP), allocatable :: y_local(:)
    !-------------------------------------------------
    
    if ( present(force_uniform_grid) ) uniform_grid = force_uniform_grid

    igs = 0; jgs = 0
    do nj=1, size(mesh2D%rcdomIJ2LCMeshID,2) 
    do ni=1, size(mesh2D%rcdomIJ2LCMeshID,1)
      n = mesh2D%rcdomIJ2LCMeshID(ni,nj)
      lcmesh => mesh2D%lcmesh_list(n)
      refElem => lcmesh%refElem2D

      allocate( x_local(refElem%Nfp), y_local(refElem%Nfp) )

      do j=1, lcmesh%NeY
      do i=1, lcmesh%NeX
        k = i + (j-1) * lcmesh%NeX
        if ( j==1 .and. nj == 1 ) then
          x_local(:) = lcmesh%pos_en(refElem%Fmask(:,1),k,1)
          if ( uniform_grid ) call get_uniform_grid1D( x_local, refElem%Nfp )

          is = igs + 1 + (i-1)*refElem%Nfp
          ie = is + refElem%Nfp - 1
          x(is:ie) = x_local(:)
        end if
        if ( i==1 .and. ni == 1 ) then
          y_local(:) = lcmesh%pos_en(refElem%Fmask(:,4),k,2)
          if ( uniform_grid ) call get_uniform_grid1D( y_local, refElem%Nfp )

          js = jgs + 1 + (j-1)*refElem%Nfp
          je = js + refElem%Nfp - 1
          y(js:je) = y_local(:)
        end if
      end do
      end do

      igs = ie; jgs = je
      deallocate( x_local, y_local )
    end do
    end do

    return
  end subroutine File_common_meshfield_get_axis2D

  subroutine File_common_meshfield_get_axis2D_cubedsphere( mesh2D, dimsinfo, x, y  )

    use scale_const, only: &
      PI => CONST_PI
    use scale_prc
    implicit none

    class(MeshCubedSphereDom2D), target, intent(in) :: mesh2D  
    type(FILE_common_meshfield_diminfo), intent(in) :: dimsinfo(MeshBase2D_DIMTYPE_NUM)
    real(DP), intent(out) :: x(dimsinfo(MeshBase2D_DIMTYPEID_X)%size)
    real(DP), intent(out) :: y(dimsinfo(MeshBase2D_DIMTYPEID_Y)%size)

    integer :: ni, nj, np, n
    integer :: k
    integer :: i, j, p
    integer :: i2, j2, p2
    type(ElementBase2D), pointer :: refElem
    type(LocalMesh2D), pointer :: lcmesh

    integer :: is, js, ie, je, igs, jgs

    logical :: uniform_grid = .false.
    real(RP), allocatable :: x_local(:)
    real(RP), allocatable :: y_local(:)
    !-------------------------------------------------
    
    igs = 0; jgs = 0
    
    do np=1, size(mesh2D%rcdomIJP2LCMeshID,3) 
    do nj=1, size(mesh2D%rcdomIJP2LCMeshID,2) 
    do ni=1, size(mesh2D%rcdomIJP2LCMeshID,1)
      n = mesh2D%rcdomIJP2LCMeshID(ni,nj,np)
      lcmesh => mesh2D%lcmesh_list(n)
      refElem => lcmesh%refElem2D

      allocate( x_local(refElem%Nfp), y_local(refElem%Nfp) )

      do j=1, lcmesh%NeY
      do i=1, lcmesh%NeX
        k = i + (j-1) * lcmesh%NeX
        if ( j==1 .and. nj == 1 .and. np == 1) then
          x_local(:) = lcmesh%pos_en(refElem%Fmask(:,1),k,1)

          is = igs + 1 + (i-1)*refElem%Nfp
          ie = is + refElem%Nfp - 1
          x(is:ie) = x_local(:)
        end if
        if ( i==1 .and. ni == 1 ) then
          y_local(:) = lcmesh%pos_en(refElem%Fmask(:,4),k,2) &
                     + ( lcmesh%panelID - 1.0_RP ) * 0.5_RP * PI

          js = jgs + 1 + (j-1)*refElem%Nfp
          je = js + refElem%Nfp - 1
          y(js:je) = y_local(:)
        end if
      end do
      end do

      igs = ie; jgs = je
      deallocate( x_local, y_local )
    end do
    end do
    end do

    return
  end subroutine File_common_meshfield_get_axis2D_cubedsphere

  subroutine File_common_meshfield_put_field2D_cartesbuf( mesh2D, field2D, &
    buf, force_uniform_grid )
    use scale_polynominal, only: &
      polynominal_genLegendrePoly
    implicit none
    class(MeshRectDom2D), target, intent(in) :: mesh2D
    class(MeshField2D), intent(in) :: field2d
    real(RP), intent(inout) :: buf(:,:)
    logical, intent(in), optional :: force_uniform_grid

    integer :: n, kelem1, p
    integer :: i0, j0, i1, j1, i2, j2, i, j
    type(LocalMesh2D), pointer :: lcmesh
    type(elementbase2D), pointer :: refElem
    integer :: i0_s, j0_s

    logical :: uniform_grid = .false.
    integer :: Nfp
    real(RP), allocatable :: x_local(:)
    real(RP) :: x_local0, delx
    real(RP), allocatable :: y_local(:)
    real(RP) :: y_local0, dely
    real(RP) :: ox, oy
    real(RP), allocatable :: spectral_coef(:)
    real(RP), allocatable :: P1D_ori_x(:,:)
    real(RP), allocatable :: P1D_ori_y(:,:)
    integer :: l, p1, p2
    !------------------------------------------------

    if ( present(force_uniform_grid) ) uniform_grid = force_uniform_grid

    i0_s = 0; j0_s = 0
    
    do j0=1, size(mesh2D%rcdomIJ2LCMeshID,2)
    do i0=1, size(mesh2D%rcdomIJ2LCMeshID,1)
      n =  mesh2D%rcdomIJ2LCMeshID(i0,j0)

      lcmesh => mesh2D%lcmesh_list(n)
      refElem => lcmesh%refElem2D
      Nfp = refElem%Nfp
      
      if ( uniform_grid ) then
        allocate( x_local(Nfp), y_local(Nfp) )
        allocate( spectral_coef(refElem%Np) )
        allocate( P1D_ori_x(1,Nfp), P1D_ori_y(1,Nfp) )
      end if
  
      do j1=1, lcmesh%NeY
      do i1=1, lcmesh%NeX
        kelem1 = i1 + (j1-1)*lcmesh%NeX

        if ( uniform_grid ) then
          x_local(:) = lcmesh%pos_en(refElem%Fmask(1:Nfp,1),kelem1,1)
          x_local0 = x_local(1); delx = x_local(Nfp) - x_local0
          y_local(:) = lcmesh%pos_en(refElem%Fmask(1:Nfp,4),kelem1,2)
          y_local0 = y_local(1); dely = y_local(Nfp) - y_local0
          call get_uniform_grid1D( x_local, Nfp )
          call get_uniform_grid1D( y_local, Nfp )
  
          spectral_coef(:) = matmul(refElem%invV(:,:), field2d%local(n)%val(:,kelem1))
          do j2=1, Nfp
          do i2=1, Nfp
            ox = - 1.0_RP + 2.0_RP * (x_local(i2) - x_local0) / delx
            oy = - 1.0_RP + 2.0_RP * (y_local(j2) - y_local0) / dely
  
            P1D_ori_x(:,:) = polynominal_genLegendrePoly( refElem%PolyOrder, (/ ox /) )
            P1D_ori_y(:,:) = polynominal_genLegendrePoly( refElem%PolyOrder, (/ oy /) )
  
            i = i0_s + i2 + (i1-1)*Nfp
            j = j0_s + j2 + (j1-1)*Nfp
            buf(i,j) = 0.0_RP 
            do p2=1, Nfp
            do p1=1, Nfp
              l = p1 + (p2-1)*Nfp
              buf(i,j) = buf(i,j) + &
                  ( P1D_ori_x(1,p1) * P1D_ori_y(1,p2) )             &
                * sqrt((dble(p1-1) + 0.5_RP)*(dble(p2-1) + 0.5_RP)) &
                * spectral_coef(l)
            end do
            end do
          end do
          end do
        
        else

          do j2=1, Nfp
          do i2=1, Nfp
            i = i0_s + i2 + (i1-1)*Nfp 
            j = j0_s + j2 + (j1-1)*Nfp
            buf(i,j) = field2d%local(n)%val(i2+(j2-1)*Nfp,kelem1)
          end do
          end do
        
        end if
      end do
      end do
      
      i0_s = i0_s + lcmesh%NeX * refElem%Nfp
      j0_s = j0_s + lcmesh%NeY * refElem%Nfp
      if ( uniform_grid ) then
        deallocate( x_local, y_local  )
        deallocate( spectral_coef )
        deallocate( P1D_ori_x, P1D_ori_y )
      end if
    end do
    end do

    return
  end subroutine File_common_meshfield_put_field2D_cartesbuf

  subroutine File_common_meshfield_put_field2D_cubedsphere_cartesbuf( mesh2D, field2D, &
    buf )
    use scale_prc, only: PRC_abort
    use scale_polynominal, only: &
      polynominal_genLegendrePoly
    implicit none
    class(MeshCubedSphereDom2D), target, intent(in) :: mesh2D
    class(MeshField2D), intent(in) :: field2d
    real(RP), intent(inout) :: buf(:,:)

    integer :: n, kelem1, p
    integer :: i0, j0, p0, i1, j1, i2, j2, i, j
    type(LocalMesh2D), pointer :: lcmesh
    type(elementbase2D), pointer :: refElem
    integer :: i0_s, j0_s
    integer :: Nfp
    !------------------------------------------------

    i0_s = 0; j0_s = 0

    do p0=1, size(mesh2D%rcdomIJP2LCMeshID,3)    
      do j0=1, size(mesh2D%rcdomIJP2LCMeshID,2)
      do i0=1, size(mesh2D%rcdomIJP2LCMeshID,1)
        n = mesh2D%rcdomIJP2LCMeshID(i0,j0,p0)

        lcmesh => mesh2D%lcmesh_list(n)
        refElem => lcmesh%refElem2D
        Nfp = refElem%Nfp
            
        do j1=1, lcmesh%NeY
        do i1=1, lcmesh%NeX
          kelem1 = i1 + (j1-1)*lcmesh%NeX

          do j2=1, Nfp
          do i2=1, Nfp
            i = i0_s + i2 + (i1-1)*Nfp
            j = j0_s + j2 + (j1-1)*Nfp
            buf(i,j) = field2d%local(n)%val(i2+(j2-1)*Nfp,kelem1)
          end do
          end do
          
        end do
        end do
        
        i0_s = i0_s + lcmesh%NeX * refElem%Nfp
        j0_s = j0_s + lcmesh%NeY * refElem%Nfp
      end do
      end do
      i0_s = 0
    end do
    
    return
  end subroutine File_common_meshfield_put_field2D_cubedsphere_cartesbuf


  subroutine File_common_meshfield_set_cartesbuf_field2D( mesh2D, buf, &
    field2D )
    implicit none
    class(MeshRectDom2D), target, intent(in) :: mesh2D
    real(RP), intent(in) :: buf(:,:)
    class(MeshField2D), intent(inout) :: field2d

    integer :: n
    integer :: i0, j0
    type(LocalMesh2D), pointer :: lcmesh
    type(elementbase2D), pointer :: refElem
    integer :: i0_s, j0_s
    !----------------------------------------------------

    i0_s = 0; j0_s = 0

    do j0=1, size(mesh2D%rcdomIJ2LCMeshID,2)
    do i0=1, size(mesh2D%rcdomIJ2LCMeshID,1)
      n = mesh2D%rcdomIJ2LCMeshID(i0,j0)
      lcmesh => mesh2D%lcmesh_list(n)
      refElem => lcmesh%refElem2D

      call File_common_meshfield_set_cartesbuf_field2D_local(  &
        lcmesh, buf(:,:), i0_s, j0_s,                          &
        field2d%local(n)%val(:,:)                              )

      i0_s = i0_s + lcmesh%NeX * refElem%Nfp
      j0_s = j0_s + lcmesh%NeY * refElem%Nfp
    end do
    end do

    return
  end subroutine File_common_meshfield_set_cartesbuf_field2D

  subroutine File_common_meshfield_set_cartesbuf_field2D_local( &
    lcmesh, buf, i0_s, j0_s,                                    &
    val )
    implicit none
    type(LocalMesh2D), intent(in) :: lcmesh
    real(RP), intent(in) :: buf(:,:)
    integer, intent(in) :: i0_s, j0_s
    real(RP), intent(inout) :: val(lcmesh%refElem2D%Np,lcmesh%NeA)

    integer :: n, kelem1, p
    integer :: i1, j1, i2, j2, i, j
    type(elementbase2D), pointer :: refElem
    integer :: indx
    !----------------------------------------------------

    refElem => lcmesh%refElem2D

    do j1=1, lcmesh%NeY
    do i1=1, lcmesh%NeX
      kelem1 = i1 + (j1-1)*lcmesh%NeX
      do j2=1, refElem%Nfp
      do i2=1, refElem%Nfp
        i = i0_s + i2 + (i1-1)*refElem%Nfp
        j = j0_s + j2 + (j1-1)*refElem%Nfp
        indx = i2 + (j2-1)*refElem%Nfp
        val(indx,kelem1) = buf(i,j)
      end do
      end do
    end do
    end do

    return
  end subroutine File_common_meshfield_set_cartesbuf_field2D_local

  subroutine File_common_meshfield_set_cartesbuf_field2D_cubedsphere( mesh2D, buf, &
    field2D )
    implicit none
    class(MeshCubedSphereDom2D), target, intent(in) :: mesh2D
    real(RP), intent(in) :: buf(:,:)
    class(MeshField2D), intent(inout) :: field2d

    integer :: n
    integer :: i0, j0, p0
    type(LocalMesh2D), pointer :: lcmesh
    type(elementbase2D), pointer :: refElem
    integer :: i0_s, j0_s, p0_s
    !----------------------------------------------------

    i0_s = 0; j0_s = 0; p0_s = 0

    do p0=1, size(mesh2D%rcdomIJP2LCMeshID,3)
      do j0=1, size(mesh2D%rcdomIJP2LCMeshID,2)
      do i0=1, size(mesh2D%rcdomIJP2LCMeshID,1)
        n = mesh2D%rcdomIJP2LCMeshID(i0,j0,p0)
        lcmesh => mesh2D%lcmesh_list(n)
        refElem => lcmesh%refElem2D

        call File_common_meshfield_set_cartesbuf_field2D_local(  &
          lcmesh, buf(:,:), i0_s, j0_s,                          &
          field2d%local(n)%val(:,:)                              )

        i0_s = i0_s + lcmesh%NeX * refElem%Nfp
        j0_s = j0_s + lcmesh%NeY * refElem%Nfp
      end do
      end do
      i0_s = 0
    end do

    return
  end subroutine File_common_meshfield_set_cartesbuf_field2D_cubedsphere

  !- 3D ------------

  subroutine File_common_meshfield_get_dims3D( mesh3D, dimsinfo )
    implicit none

    class(MeshCubeDom3D), target, intent(in) :: mesh3D
    type(FILE_common_meshfield_diminfo), intent(out) :: dimsinfo(MeshBase3D_DIMTYPE_NUM)

    type(LocalMesh3D), pointer :: lcmesh
    integer :: i, j, k, n
    integer :: icount, jcount, kcount

    integer :: i_size, j_size, k_size

    type(MeshDimInfo), pointer :: dimInfo
    type(MeshDimInfo), pointer :: dimInfo_x
    type(MeshDimInfo), pointer :: dimInfo_y
    type(MeshDimInfo), pointer :: dimInfo_z
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

    diminfo_x => mesh3D%dimInfo(MeshBase3D_DIMTYPEID_X)
    call set_dimension( dimsinfo(MeshBase3D_DIMTYPEID_X),   &
      dimInfo_x, "X", 1, (/ diminfo_x%name /), (/ i_size /) )

    diminfo_y => mesh3D%dimInfo(MeshBase3D_DIMTYPEID_Y)
    call set_dimension( dimsinfo(MeshBase3D_DIMTYPEID_Y),   &
      diminfo_y, "Y", 1, (/ diminfo_y%name /), (/ j_size /) )

    diminfo_z => mesh3D%dimInfo(MeshBase3D_DIMTYPEID_Z)
    call set_dimension( dimsinfo(MeshBase3D_DIMTYPEID_Z),   &
      diminfo_z, "Z", 1, (/ diminfo_z%name /), (/ k_size /) )
  
    diminfo => mesh3D%dimInfo(MeshBase3D_DIMTYPEID_ZT)
    call set_dimension( dimsinfo(MeshBase3D_DIMTYPEID_ZT), &
      diminfo, "ZT", 1, (/ diminfo_z%name /), (/ k_size /) )

    diminfo => mesh3D%dimInfo(MeshBase3D_DIMTYPEID_XYZ)
    call set_dimension( dimsinfo(MeshBase3D_DIMTYPEID_XYZ), &
      diminfo, "XYZ", 3, (/ diminfo_x%name, diminfo_y%name, diminfo_z%name /), &
      (/ i_size, j_size, k_size /) )
  
    diminfo => mesh3D%dimInfo(MeshBase3D_DIMTYPEID_XYZT)
    call set_dimension( dimsinfo(MeshBase3D_DIMTYPEID_XYZT), &
      diminfo, "XYZT", 3, (/ diminfo_x%name, diminfo_y%name, diminfo_z%name /), &
      (/ i_size, j_size, k_size /) )

    return
  end subroutine File_common_meshfield_get_dims3D

  subroutine File_common_meshfield_get_axis3D( mesh3D, dimsinfo, x, y, z, &
    force_uniform_grid )
    implicit none

    class(MeshCubeDom3D), target, intent(in) :: mesh3D  
    type(FILE_common_meshfield_diminfo), intent(in) :: dimsinfo(MeshBase3D_DIMTYPE_NUM)
    real(DP), intent(out) :: x(dimsinfo(MeshBase3D_DIMTYPEID_X)%size)
    real(DP), intent(out) :: y(dimsinfo(MeshBase3D_DIMTYPEID_Y)%size)
    real(DP), intent(out) :: z(dimsinfo(MeshBase3D_DIMTYPEID_Z)%size)
    logical, intent(in), optional :: force_uniform_grid

    integer :: n, kelem
    integer :: i, j, k
    integer :: i2, j2, k2
    type(ElementBase3D), pointer :: refElem
    type(LocalMesh3D), pointer :: lcmesh

    integer :: is, js, ks, ie, je, ke, igs, jgs, kgs
    integer :: Nnode_h1D, Nnode_v

    logical :: uniform_grid = .false.
    real(RP), allocatable :: x_local(:)
    real(RP), allocatable :: y_local(:)
    real(RP), allocatable :: z_local(:)
    !------------------------------------------------------------------------------------------
    
    if ( present(force_uniform_grid) ) uniform_grid = force_uniform_grid

    igs = 0; jgs = 0; kgs = 0

    do n=1 ,mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      refElem => lcmesh%refElem3D
      Nnode_h1D = refElem%Nnode_h1D
      Nnode_v   = refElem%Nnode_v

      allocate( x_local(Nnode_h1D), y_local(Nnode_h1D) )
      allocate( z_local(Nnode_v) )

      do k=1, lcmesh%NeZ
      do j=1, lcmesh%NeY
      do i=1, lcmesh%NeX
        kelem = i + (j-1)*lcmesh%NeX + (k-1)*lcmesh%NeX*lcmesh%NeY
        if ( j==1 .and. k==1) then
          x_local(:) = lcmesh%pos_en(refElem%Fmask_h(1:Nnode_h1D,1),kelem,1)
          if ( uniform_grid ) call get_uniform_grid1D( x_local, Nnode_h1D )

          is = igs + 1 + (i-1)*Nnode_h1D
          ie = is + Nnode_h1D - 1
          x(is:ie) = x_local(:)
        end if
        if ( i==1 .and. k==1) then
          y_local(:) = lcmesh%pos_en(refElem%Fmask_h(1:Nnode_h1D,4),kelem,2)
          if ( uniform_grid ) call get_uniform_grid1D( y_local, Nnode_h1D )

          js = jgs + 1 + (j-1)*Nnode_h1D
          je = js + Nnode_h1D - 1
          y(js:je) = y_local(:)
        end if
        if ( i==1 .and. j==1) then
          z_local(:) = lcmesh%pos_en(refElem%Colmask(:,1),kelem,3)
          if ( uniform_grid ) call get_uniform_grid1D( z_local, Nnode_v )

          ks = kgs + 1 + (k-1)*Nnode_v
          ke = ks + Nnode_v - 1
          z(ks:ke) = z_local(:)
        end if
      end do
      end do
      end do

      igs = ie; jgs = je; kgs = ke
      deallocate( x_local, y_local )
      deallocate( z_local )
    end do

    return
  end subroutine File_common_meshfield_get_axis3D

  subroutine File_common_meshfield_put_field3D_cartesbuf( mesh3D, field3D, &
    buf, force_uniform_grid )
    use scale_polynominal, only: &
      polynominal_genLegendrePoly    
    implicit none
    class(MeshCubeDom3D), target, intent(in) :: mesh3D
    class(MeshField3D), intent(in) :: field3d
    real(RP), intent(inout) :: buf(:,:,:)
    logical, intent(in), optional :: force_uniform_grid

    integer :: n, kelem1, p
    integer :: i0, j0, k0, i1, j1, k1, i2, j2, k2, i, j, k
    type(LocalMesh3D), pointer :: lcmesh
    type(elementbase3D), pointer :: refElem
    integer :: i0_s, j0_s, k0_s, indx

    logical :: uniform_grid = .false.
    integer :: Nnode_h1D, Nnode_v
    real(RP), allocatable :: x_local(:)
    real(RP) :: x_local0, delx
    real(RP), allocatable :: y_local(:)
    real(RP) :: y_local0, dely
    real(RP), allocatable :: z_local(:)
    real(RP) :: z_local0, delz
    real(RP) :: ox, oy, oz
    real(RP), allocatable :: spectral_coef(:)
    real(RP), allocatable :: P1D_ori_x(:,:)
    real(RP), allocatable :: P1D_ori_y(:,:)
    real(RP), allocatable :: P1D_ori_z(:,:)
    integer :: l, p1, p2, p3
    !----------------------------------------------------

    if ( present(force_uniform_grid) ) uniform_grid = force_uniform_grid

    i0_s = 0; j0_s = 0; k0_s = 0

    do k0=1, size(mesh3D%rcdomIJK2LCMeshID,3)  
    do j0=1, size(mesh3D%rcdomIJK2LCMeshID,2)
    do i0=1, size(mesh3D%rcdomIJK2LCMeshID,1)
      n =  mesh3D%rcdomIJK2LCMeshID(i0,j0,k0)

      lcmesh => mesh3D%lcmesh_list(n)
      refElem => lcmesh%refElem3D
      Nnode_h1D = refElem%Nnode_h1D
      Nnode_v   = refElem%Nnode_v

      if ( uniform_grid ) then
        allocate( x_local(Nnode_h1D), y_local(Nnode_h1D) )
        allocate( z_local(Nnode_v) ) 
        allocate( spectral_coef(refElem%Np) )
        allocate( P1D_ori_x(1,Nnode_h1D), P1D_ori_y(1,Nnode_h1D) )
        allocate( P1D_ori_z(1,Nnode_v) )     
      end if

    !$omp parallel do collapse(2) private( kelem1, &
    !$omp i, i1, i2, j, j2, k, k2, indx,                               &
    !$omp x_local, x_local0, y_local, y_local0, z_local, z_local0,     &
    !$omp delx, dely, delz, ox, oy, oz,                                &
    !$omp spectral_coef, P1D_ori_x, P1D_ori_y, P1D_ori_z,              &
    !$omp p1, p2, p3, l                                                )        
      do k1=1, lcmesh%NeZ
      do j1=1, lcmesh%NeY
      do i1=1, lcmesh%NeX
        kelem1 = i1 + (j1-1)*lcmesh%NeX + (k1-1)*lcmesh%NeX*lcmesh%NeY

        if ( uniform_grid ) then
          x_local(:) = lcmesh%pos_en(refElem%Fmask_h(1:Nnode_h1D,1),kelem1,1)
          x_local0 = x_local(1); delx = x_local(Nnode_h1D) - x_local0
          y_local(:) = lcmesh%pos_en(refElem%Fmask_h(1:Nnode_h1D,4),kelem1,2)
          y_local0 = y_local(1); dely = y_local(Nnode_h1D) - y_local0
          z_local(:) = lcmesh%pos_en(refElem%Colmask(:,1),kelem1,3)
          z_local0 = z_local(1); delz = z_local(Nnode_v  ) - z_local0
          call get_uniform_grid1D( x_local, Nnode_h1D )
          call get_uniform_grid1D( y_local, Nnode_h1D )
          call get_uniform_grid1D( z_local, Nnode_v   )

          spectral_coef(:) = matmul(refElem%invV(:,:), field3d%local(n)%val(:,kelem1))
          do k2=1, Nnode_v
          do j2=1, Nnode_h1D
          do i2=1, Nnode_h1D
            ox = - 1.0_RP + 2.0_RP * (x_local(i2) - x_local0) / delx
            oy = - 1.0_RP + 2.0_RP * (y_local(j2) - y_local0) / dely
            oz = - 1.0_RP + 2.0_RP * (z_local(k2) - z_local0) / delz
  
            P1D_ori_x(:,:) = polynominal_genLegendrePoly( refElem%PolyOrder_h, (/ ox /) )
            P1D_ori_y(:,:) = polynominal_genLegendrePoly( refElem%PolyOrder_h, (/ oy /) )
            P1D_ori_z(:,:) = polynominal_genLegendrePoly( refElem%PolyOrder_v, (/ oz /) )
  
            i = i0_s + i2 + (i1-1)*Nnode_h1D
            j = j0_s + j2 + (j1-1)*Nnode_h1D
            k = k0_s + k2 + (k1-1)*Nnode_v
            buf(i,j,k) = 0.0_RP 
            do p3=1, Nnode_v
            do p2=1, Nnode_h1D
            do p1=1, Nnode_h1D
              l = p1 + (p2-1)*Nnode_h1D + (p3-1)*Nnode_h1D**2
              buf(i,j,k) = buf(i,j,k) + &
                  ( P1D_ori_x(1,p1) * P1D_ori_y(1,p2) * P1D_ori_z(1,p3) )                 &
                * sqrt((dble(p1-1) + 0.5_RP)*(dble(p2-1) + 0.5_RP)*(dble(p3-1) + 0.5_RP)) &
                * spectral_coef(l)
            end do
            end do
            end do            
          end do
          end do
          end do
        else
          do k2=1, Nnode_v
          do j2=1, Nnode_h1D
          do i2=1, Nnode_h1D
            i = i0_s + i2 + (i1-1)*Nnode_h1D
            j = j0_s + j2 + (j1-1)*Nnode_h1D
            k = k0_s + k2 + (k1-1)*Nnode_v
            indx = i2 + (j2-1)*Nnode_h1D + (k2-1)*Nnode_h1D**2
            buf(i,j,k) = field3d%local(n)%val(indx,kelem1)
          end do
          end do
          end do
        end if
      end do
      end do
      end do

      i0_s = i0_s + lcmesh%NeX * refElem%Nnode_h1D
      j0_s = j0_s + lcmesh%NeY * refElem%Nnode_h1D
      k0_s = k0_s + lcmesh%NeZ * refElem%Nnode_v
      if ( uniform_grid ) then
        deallocate( x_local, y_local )
        deallocate( z_local )
        deallocate( spectral_coef )
      end if
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

  subroutine get_uniform_grid1D( pos1D, Np )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(inout) :: pos1D(Np)

    real(RP) :: del
    integer :: i
    !-----------------------------------------------

    del = ( pos1D(Np) - pos1D(1) ) / dble(Np)
    pos1D(1) = pos1D(1) + 0.5_RP * del
    do i=2, Np
      pos1D(i) = pos1D(i-1) + del
    end do

    return
  end subroutine get_uniform_grid1D

  subroutine set_dimension( dim, diminfo, dim_type, ndims, dims, count )
    implicit none

    type(FILE_common_meshfield_diminfo), intent(out) :: dim
    type(MeshDimInfo), intent(in) :: diminfo
    character(*), intent(in) :: dim_type
    integer, intent(in) :: ndims
    character(len=*), intent(in) :: dims(ndims)
    integer, intent(in) :: count(ndims)

    integer :: d
    !----------------------------------------------------

    dim%name = diminfo%name
    dim%unit = diminfo%unit
    dim%desc = diminfo%desc
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
