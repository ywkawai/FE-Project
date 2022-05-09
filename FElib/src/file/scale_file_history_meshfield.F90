!-------------------------------------------------------------------------------
!> module file_history_meshfield
!!
!! @par Description
!!           managing file history 
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_file_history_meshfield
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io

  use scale_file_history, only: &
    FILE_HISTORY_Setup, &
    FILE_HISTORY_Set_NowDate, &
    FILE_HISTORY_truncate_1D, &
    FILE_HISTORY_truncate_2D, &
    FILE_HISTORY_truncate_3D, &
    FILE_HISTORY_put,         &
    FILE_HISTORY_write,       &
    FILE_HISTORY_Set_Dim,     &
    FILE_HISTORY_Set_Axis,    &
    FILE_HISTORY_finalize

  use scale_element_base, only: elementbase1D, elementbase2D, elementbase3D
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D

  use scale_meshfield_base, only: &
    MeshFieldBase,                        &
    MeshField1D, MeshField2D, MeshField3D    

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: FILE_HISTORY_meshfield_setup
  public :: FILE_HISTORY_meshfield_put
  public :: FILE_HISTORY_meshfield_in
  public :: FILE_HISTORY_meshfield_write  
  public :: FILE_HISTORY_meshfield_finalize

  interface FILE_HISTORY_meshfield_put
    module procedure FILE_HISTORY_meshfield_put1D
    module procedure FILE_HISTORY_meshfield_put2D
    module procedure FILE_HISTORY_meshfield_put3D
  end interface
  
  interface FILE_HISTORY_meshfield_in
    module procedure FILE_HISTORY_meshfield_in1D
    module procedure FILE_HISTORY_meshfield_in2D
    module procedure FILE_HISTORY_meshfield_in3D
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !-------------------
  integer,          parameter :: nzs = 1
  character(len=8), parameter :: zs(nzs) = (/  "model   " /)

  integer  :: FILE_HISTORY_MESHFIELD_STARTDATE(6) !< start time [YYYY MM DD HH MM SS]
  real(DP) :: FILE_HISTORY_MESHFIELD_STARTSUBSEC      !< subsecond part of start time [millisec]

  class(MeshBase1D), pointer :: mesh1D
  class(MeshRectDom2D), pointer :: mesh2D
  class(MeshCubeDom3D), pointer :: mesh3D
  class(MeshCubedSphereDom2D), pointer :: meshCubedSphere2D
  class(MeshCubedSphereDom3D), pointer :: meshCubedSphere3D

  integer :: dims1D_size(1)
  integer :: dims2D_size(2)
  integer :: dims3D_size(3,nzs)

  logical, private :: FILE_HISTORY_FILEMESHFILED_DISABLE = .true.

contains

!----------------

  subroutine FILE_HISTORY_meshfield_setup( &
    mesh1D_, mesh2D_, mesh3D_,             &
    meshcubedsphere2D_, meshcubedsphere3D_ )

    use scale_file_h, only: &
      FILE_HSHORT
    use scale_prc, only: &
      PRC_masterrank, &
      PRC_myrank,     &
      PRC_abort
    use scale_time, only: &
      TIME_NOWDATE,       &
      TIME_NOWSUBSEC,     &
      TIME_STARTDAYSEC,   &
      TIME_DTSEC,         &
      TIME_NOWSTEP
    use scale_calendar, only: &
      CALENDAR_get_name
    implicit none

    class(Meshbase1d), intent(in), target, optional :: mesh1D_
    class(MeshRectDom2d), intent(in), target, optional :: mesh2D_
    class(MeshCubeDom3D), intent(in), target, optional :: mesh3D_
    class(MeshCubedSphereDom2D), intent(in), target, optional :: meshCubedsphere2D_
    class(MeshCubedSphereDom3D), intent(in), target, optional :: meshCubedsphere3D_

    character(len=H_MID) :: FILE_HISTORY_MESHFILED_H_TITLE = 'SCALE-FEM FILE_HISTORY_MESHFIELD' !< title of the output file
    character(len=H_MID) :: FILE_HISTORY_MESHFIELD_T_SINCE

    
    character(len=FILE_HSHORT) :: calendar
    real(DP) :: start_daysec
    !---------------------------------------------------------------------------

    FILE_HISTORY_MESHFIELD_STARTDATE(:) = TIME_NOWDATE
    FILE_HISTORY_MESHFIELD_STARTSUBSEC  = TIME_NOWSUBSEC

    start_daysec = TIME_STARTDAYSEC
    if ( TIME_NOWDATE(1) > 0 ) then
      write(FILE_HISTORY_MESHFIELD_T_SINCE,'(I4.4,5(A1,I2.2))') TIME_NOWDATE(1), &
                                                          '-', TIME_NOWDATE(2), &
                                                          '-', TIME_NOWDATE(3), &
                                                          ' ', TIME_NOWDATE(4), &
                                                          ':', TIME_NOWDATE(5), &
                                                          ':', TIME_NOWDATE(6)
      start_daysec = TIME_NOWSUBSEC
    else
      FILE_HISTORY_MESHFIELD_T_SINCE = ''
    endif

    ! get calendar name
    call CALENDAR_get_name( calendar )

    call FILE_HISTORY_Setup( FILE_HISTORY_MESHFILED_H_TITLE,        & ! [IN]
      H_SOURCE, H_INSTITUTE,                                        & ! [IN]
      start_daysec, TIME_DTSEC,                                     & ! [IN]
      time_since = FILE_HISTORY_MESHFIELD_T_SINCE,                  & ! [IN]
      calendar = calendar,                                          & ! [IN]
      default_zcoord = 'model',                                     & ! [IN]
      myrank = PRC_myrank                        )                    ! [IN]

    call FILE_HISTORY_Set_NowDate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )

    !- Set a pointer to the  variable of mesh 
    
    nullify( mesh1D, mesh2D, mesh3D )
    nullify( meshCubedSphere2D, meshCubedsphere3D )

    if ( present(mesh1D_) ) then
      mesh1D => mesh1D_
      call set_dim_axis1D()
    else if ( present(mesh2D_) ) then
      mesh2D => mesh2D_
      call set_dim_axis2D()
    else if ( present(mesh3D_) ) then
      mesh3D => mesh3D_
      mesh2D => mesh3D%mesh2D
      call set_dim_axis3D()
    else if ( present(meshCubedsphere2D_) ) then
      meshCubedSphere2D => meshCubedsphere2D_
      call set_dim_axis2D_cubedsphere()
    else if ( present(meshCubedsphere3D_) ) then
      meshCubedSphere3D => meshCubedsphere3D_
      meshCubedSphere2D => meshCubedSphere3D%mesh2D
      call set_dim_axis3D_cubedsphere()
    else
      LOG_ERROR("FILE_HISTORY_meshfield_setup",*)   "Any mesh (mesh1d/2d/3d) are not specified."
      call PRC_abort
    end if

    FILE_HISTORY_FILEMESHFILED_DISABLE = .false. 

    return
  end subroutine FILE_HISTORY_meshfield_setup

  subroutine FILE_HISTORY_meshfield_write()
    implicit none
    !-------------------------------------------------

    call FILE_HISTORY_write
    return
  end subroutine FILE_HISTORY_meshfield_write

  subroutine FILE_HISTORY_meshfield_finalize()
    implicit none
    !-------------------------------------------------

    call FILE_HISTORY_finalize()
    return
  end subroutine FILE_HISTORY_meshfield_finalize

  !-- 1D

  subroutine FILE_HISTORY_meshfield_put1D(hstid, field1d)
    use scale_file_common_meshfield, only: &
      File_common_meshfield_put_field1D_cartesbuf
    implicit none
    integer, intent(in) :: hstid
    class(MeshField1D), intent(in) :: field1d

    real(RP), allocatable :: buf(:)
    !-------------------------------------------------
      
    allocate( buf(dims1D_size(1)) )

    call File_common_meshfield_put_field1D_cartesbuf( mesh1D, field1d, buf(:) )
    call FILE_HISTORY_put(hstid, buf)

    return
  end subroutine FILE_HISTORY_meshfield_put1D

  subroutine FILE_HISTORY_meshfield_in1D( field1d, desc, standard_name )
    implicit none
    class(MeshField1D), intent(in) :: field1d
    character(len=*), intent(in) :: desc                     !< description of the item
    character(len=*), intent(in), optional :: standard_name

    integer :: hstid
    logical :: do_put
    !-------------------------------------------------

    call history_in_regvar( hstid, do_put,    & ! (out)
      field1d, desc, 1, standard_name, 'XYZ'  ) ! (in)

    if ( do_put ) call FILE_HISTORY_meshfield_put( hstid, field1d )

    return
  end subroutine FILE_HISTORY_meshfield_in1D
  
  !-- 2D

  subroutine FILE_HISTORY_meshfield_put2D(hstid, field2d)
    use scale_file_common_meshfield, only: &
      File_common_meshfield_put_field2D_cartesbuf,            &
      File_common_meshfield_put_field2D_cubedsphere_cartesbuf

    implicit none
    integer, intent(in) :: hstid
    class(MeshField2D), intent(in) :: field2d

    real(RP), allocatable :: buf(:,:)
    !-------------------------------------------------

    allocate( buf(dims2D_size(1),dims2D_size(2)) )
    
    if ( associated(mesh2D) ) then
      call File_common_meshfield_put_field2D_cartesbuf( mesh2D, field2d, buf(:,:) )
    else if ( associated(meshCubedSphere2D) ) then
      call File_common_meshfield_put_field2D_cubedsphere_cartesbuf( &
        meshCubedSphere2D, field2d, buf(:,:) )
    end if
    call FILE_HISTORY_put(hstid, buf)

    return
  end subroutine FILE_HISTORY_meshfield_put2D

  subroutine FILE_HISTORY_meshfield_in2D( field2d, desc, standard_name )
    implicit none
    class(MeshField2D), intent(in) :: field2d
    character(len=*), intent(in) :: desc                     !< description of the item
    character(len=*), intent(in), optional :: standard_name

    integer :: hstid
    logical :: do_put
    !-------------------------------------------------

    call history_in_regvar( hstid, do_put,  & ! (out)
      field2d, desc, 2, standard_name, 'XY' ) ! (in)

    if ( do_put ) call FILE_HISTORY_meshfield_put( hstid, field2d )

    return
  end subroutine FILE_HISTORY_meshfield_in2D

  !-- 3D

  subroutine FILE_HISTORY_meshfield_put3D(hstid, field3d)
    use scale_file_common_meshfield, only: &
      File_common_meshfield_put_field3D_cartesbuf,            &
      File_common_meshfield_put_field3D_cubedsphere_cartesbuf

    implicit none
    integer, intent(in) :: hstid
    class(MeshField3D), intent(in) :: field3d

    real(RP), allocatable :: buf(:,:,:)

    !-------------------------------------------------

    allocate( buf(dims3D_size(1,1),dims3D_size(2,1),dims3D_size(3,1)) )
    if ( associated(mesh3D) ) then
      call File_common_meshfield_put_field3D_cartesbuf( mesh3D, field3d, buf(:,:,:) )
    else if ( associated(meshCubedSphere3D) ) then
      call File_common_meshfield_put_field3D_cubedsphere_cartesbuf( &
        meshCubedSphere3D, field3d, buf(:,:,:) )
    end if
    call FILE_HISTORY_put(hstid, buf)

    return
  end subroutine FILE_HISTORY_meshfield_put3D
  
  subroutine FILE_HISTORY_meshfield_in3D( field3d, desc, standard_name )
    implicit none
    class(MeshField3D), intent(in) :: field3d
    character(len=*), intent(in) :: desc                     !< description of the item
    character(len=*), intent(in), optional :: standard_name

    integer :: hstid
    logical :: do_put
    !-------------------------------------------------

    call history_in_regvar( hstid, do_put,   & ! (out)
      field3d, desc, 3, standard_name, 'XYZ' ) ! (in)

    if ( do_put ) call FILE_HISTORY_meshfield_put( hstid, field3d )

    return
  end subroutine FILE_HISTORY_meshfield_in3D

!----------------

  subroutine history_in_regvar( hstid, do_put,  &
    field, desc, ndim, standard_name, dim_type  )

    use scale_file_history, only: &
      FILE_HISTORY_reg,           &
      FILE_HISTORY_query
    
    implicit none

    integer, intent(out) :: hstid
    logical, intent(out) :: do_put
    class(MeshFieldBase), intent(in) :: field  
    character(len=*), intent(in) :: desc       !< description of the item
    integer, intent(in)          :: ndim
    character(len=*), intent(in), optional :: standard_name
    character(len=*), intent(in), optional :: dim_type

    logical, parameter     :: fill_halo = .false.
    !------------------------------------------------------

    hstid = -1
    do_put = .false.

    if ( FILE_HISTORY_FILEMESHFILED_DISABLE ) return

    ! Check whether the item has been already registered
    call FILE_HISTORY_reg( field%varname, desc, field%unit, & ! [IN]
                           hstid,                           & ! [OUT]
                           standard_name=standard_name,     & ! [IN]
                           ndims=ndim,                      & ! [IN]
                           dim_type=dim_type,               & ! [IN]
                           fill_halo=fill_halo              ) ! [IN]
    
    if ( hstid < 0 ) return

    ! Check whether it is time to input the item
    call FILE_HISTORY_query( hstid, do_put ) ! [IN], [OUT]
    
    return
  end subroutine history_in_regvar

  subroutine set_dim_axis1D()
    use scale_file_common_meshfield, only: &
      FILE_common_meshfield_diminfo,    &
      File_common_meshfield_get_dims1D, &
      File_common_meshfield_get_axis1D
    use scale_mesh_base1d, only: &
      DIMTYPE_NUM => MeshBase1D_DIMTYPE_NUM, &
      DIMTYPE_X   => MeshBase1D_DIMTYPEID_X
  
    implicit none

    type(FILE_common_meshfield_diminfo) :: dimsinfo(DIMTYPE_NUM)
    real(RP), allocatable :: x(:)
    integer :: start(1,1), count(1,1)
    character(len=H_SHORT) :: dims(1,1)
    integer :: n, ndim
    !-------------------------------------------------
    
    call File_common_meshfield_get_dims1D( mesh1D, & ! (in)
      dimsinfo(:) )                                  ! (out)
    
    dims1D_size(1) = dimsinfo(DIMTYPE_X)%size
    allocate( x(dims1D_size(1)) )
    call File_common_meshfield_get_axis1D( mesh1D, dimsinfo, & ! (in)
      x )                                                      ! (out)

    start(:,:) = 1
    do n=1, DIMTYPE_NUM
      ndim = dimsinfo(n)%ndim
      dims(1:ndim,1)  = dimsinfo(n)%dims(1:ndim)
      count(1:ndim,1) = dimsinfo(n)%count(1:ndim)
      call FILE_HISTORY_Set_Dim ( dimsinfo(n)%type, ndim, 1, dims(1:ndim,:), zs(:), start(1:ndim,:), count(1:ndim,:))
    end do
    
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_X)%name, dimsinfo(DIMTYPE_X)%desc, &
      dimsinfo(DIMTYPE_X)%unit, dimsinfo(DIMTYPE_X)%name, x(:))

    return
  end subroutine set_dim_axis1D

  subroutine set_dim_axis2D()
    use scale_file_common_meshfield, only: &
      FILE_common_meshfield_diminfo,       &
      File_common_meshfield_get_dims2D,    &
      File_common_meshfield_get_axis2D
    use scale_mesh_base2d, only: &
      DIMTYPE_NUM => MeshBase2D_DIMTYPE_NUM, &
      DIMTYPE_X   => MeshBase2D_DIMTYPEID_X, &
      DIMTYPE_Y   => MeshBase2D_DIMTYPEID_Y
   
    implicit none

    type(FILE_common_meshfield_diminfo) :: dimsinfo(DIMTYPE_NUM)
    real(RP), allocatable :: x(:), y(:)
    integer :: start(2,1), count(2,1)
    character(len=H_SHORT) :: dims(2,1)
    integer :: n, ndim
    !-------------------------------------------------
    
    call File_common_meshfield_get_dims2D( mesh2D, & ! (in)
      dimsinfo(:) )                                  ! (out)
    
    dims2D_size(1) = dimsinfo(DIMTYPE_X)%size
    dims2D_size(2) = dimsinfo(DIMTYPE_Y)%size
    allocate( x(dims2D_size(1)), y(dims2D_size(2)) )
    call File_common_meshfield_get_axis2D( mesh2D, dimsinfo, & ! (in)
      x, y )                                                   ! (out)

    start(:,:) = 1

    do n=1, DIMTYPE_NUM
      ndim = dimsinfo(n)%ndim
      dims(1:ndim,1)  = dimsinfo(n)%dims(1:ndim)
      count(1:ndim,1) = dimsinfo(n)%count(1:ndim)
      call FILE_HISTORY_Set_Dim ( dimsinfo(n)%type, ndim, 1, dims(1:ndim,:), zs(:), start(1:ndim,:), count(1:ndim,:))
    end do
    
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_X)%name, dimsinfo(DIMTYPE_X)%desc, dimsinfo(DIMTYPE_X)%unit, dimsinfo(DIMTYPE_X)%name, x(:))
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_Y)%name, dimsinfo(DIMTYPE_Y)%desc, dimsinfo(DIMTYPE_Y)%unit, dimsinfo(DIMTYPE_Y)%name, y(:))
    
    return
  end subroutine set_dim_axis2D

  subroutine set_dim_axis3D()
    use scale_file_common_meshfield, only: &
      FILE_common_meshfield_diminfo,       &
      File_common_meshfield_get_dims3D,    &
      File_common_meshfield_get_axis3D
    use scale_mesh_base3d, only: &
      DIMTYPE_NUM => MeshBase3D_DIMTYPE_NUM, &
      DIMTYPE_X   => MeshBase3D_DIMTYPEID_X, &
      DIMTYPE_Y   => MeshBase3D_DIMTYPEID_Y, &
      DIMTYPE_Z   => MeshBase3D_DIMTYPEID_Z

    implicit none

    type(FILE_common_meshfield_diminfo) :: dimsinfo(DIMTYPE_NUM)
    real(RP), allocatable :: x(:), y(:), z(:)
    integer :: start(3,1), count(3,1)
    character(len=H_SHORT) :: dims(3,1)
    integer :: n, ndim
    !-------------------------------------------------
    
    call File_common_meshfield_get_dims3D( mesh3D, & ! (in)
      dimsinfo(:) )                                  ! (out)
  
    dims2D_size(1) = dimsinfo(DIMTYPE_X)%size
    dims2D_size(2) = dimsinfo(DIMTYPE_Y)%size  
    dims3D_size(1,1) = dimsinfo(DIMTYPE_X)%size
    dims3D_size(2,1) = dimsinfo(DIMTYPE_Y)%size
    dims3D_size(3,1) = dimsinfo(DIMTYPE_Z)%size
    allocate( x(dims3D_size(1,1)), y(dims3D_size(2,1)), z(dims3D_size(3,1)) )
    call File_common_meshfield_get_axis3D( mesh3D, dimsinfo, & ! (in)
      x, y, z )                                                ! (out)

    start(:,:) = 1
    do n=1, DIMTYPE_NUM
      ndim = dimsinfo(n)%ndim
      dims(1:ndim,1)  = dimsinfo(n)%dims(1:ndim)
      count(1:ndim,1) = dimsinfo(n)%count(1:ndim)
      call FILE_HISTORY_Set_Dim ( dimsinfo(n)%type, ndim, 1, dims(1:ndim,:), zs(:), start(1:ndim,:), count(1:ndim,:))
    end do
    
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_X)%name, dimsinfo(DIMTYPE_X)%desc, dimsinfo(DIMTYPE_X)%unit, dimsinfo(DIMTYPE_X)%name, x(:) )
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_Y)%name, dimsinfo(DIMTYPE_Y)%desc, dimsinfo(DIMTYPE_Y)%unit, dimsinfo(DIMTYPE_Y)%name, y(:) )
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_Z)%name, dimsinfo(DIMTYPE_Z)%desc, dimsinfo(DIMTYPE_Z)%unit, dimsinfo(DIMTYPE_Z)%name, z(:), &
      down=dimsinfo(DIMTYPE_Z)%positive_down(1) )

    return
  end subroutine set_dim_axis3D

  subroutine set_dim_axis2D_cubedsphere()
    use scale_file_common_meshfield, only: &
      FILE_common_meshfield_diminfo,  &
      File_common_meshfield_get_dims, &
      File_common_meshfield_get_axis
    use scale_mesh_base2d, only: &
      DIMTYPE_NUM => MeshBase2D_DIMTYPE_NUM, &
      DIMTYPE_X  => MeshBase2D_DIMTYPEID_X,  &
      DIMTYPE_Y  => MeshBase2D_DIMTYPEID_Y

    implicit none

    type(FILE_common_meshfield_diminfo) :: dimsinfo(DIMTYPE_NUM)
    real(RP), allocatable :: x(:), y(:)
    integer :: start(2,1), count(2,1)
    character(len=H_SHORT) :: dims(2,1)
    integer :: n, ndim
    !-------------------------------------------------
    
    call File_common_meshfield_get_dims( meshCubedSphere2D, & ! (in)
      dimsinfo(:) )                                           ! (out)
    
    dims2D_size(1) = dimsinfo(DIMTYPE_X)%size
    dims2D_size(2) = dimsinfo(DIMTYPE_Y)%size
    allocate( x(dims2D_size(1)), y(dims2D_size(2)) )
    
    call File_common_meshfield_get_axis( meshCubedSphere2D, dimsinfo, & ! (in)
      x, y )                                                            ! (out)
    
    start(:,:) = 1

    do n=1, DIMTYPE_NUM
      ndim = dimsinfo(n)%ndim
      dims(1:ndim,1)  = dimsinfo(n)%dims(1:ndim)
      count(1:ndim,1) = dimsinfo(n)%count(1:ndim)
      call FILE_HISTORY_Set_Dim ( dimsinfo(n)%type, ndim, 1, dims(1:ndim,:), zs(:), start(1:ndim,:), count(1:ndim,:))
    end do
    
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_X)%name, dimsinfo(DIMTYPE_X)%desc, dimsinfo(DIMTYPE_X)%unit, dimsinfo(DIMTYPE_X)%name, x(:))
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_Y)%name, dimsinfo(DIMTYPE_Y)%desc, dimsinfo(DIMTYPE_Y)%unit, dimsinfo(DIMTYPE_Y)%name, y(:))
    
    return
  end subroutine set_dim_axis2D_cubedsphere


  subroutine set_dim_axis3D_cubedsphere()
    use scale_file_common_meshfield, only: &
      FILE_common_meshfield_diminfo,  &
      File_common_meshfield_get_dims, &
      File_common_meshfield_get_axis
    use scale_mesh_base3d, only: &
      DIMTYPE_NUM => MeshBase3D_DIMTYPE_NUM, &
      DIMTYPE_X   => MeshBase3D_DIMTYPEID_X, &
      DIMTYPE_Y   => MeshBase3D_DIMTYPEID_Y, &
      DIMTYPE_Z   => MeshBase3D_DIMTYPEID_Z

    implicit none

    type(FILE_common_meshfield_diminfo) :: dimsinfo(DIMTYPE_NUM)
    real(RP), allocatable :: x(:), y(:), z(:)
    integer :: start(3,1), count(3,1)
    character(len=H_SHORT) :: dims(3,1)
    integer :: n, ndim
    !-------------------------------------------------
    
    call File_common_meshfield_get_dims( meshCubedSphere3D, & ! (in)
      dimsinfo(:) )                                           ! (out)

    dims2D_size(1) = dimsinfo(DIMTYPE_X)%size
    dims2D_size(2) = dimsinfo(DIMTYPE_Y)%size
    dims3D_size(1,1) = dimsinfo(DIMTYPE_X)%size
    dims3D_size(2,1) = dimsinfo(DIMTYPE_Y)%size
    dims3D_size(3,1) = dimsinfo(DIMTYPE_Z)%size
    allocate( x(dims3D_size(1,1)), y(dims3D_size(2,1)), z(dims3D_size(3,1)) )
    call File_common_meshfield_get_axis( meshCubedSphere3D, dimsinfo, & ! (in)
      x, y, z )                                                         ! (out)

    start(:,:) = 1
    do n=1, DIMTYPE_NUM
      ndim = dimsinfo(n)%ndim
      dims(1:ndim,1)  = dimsinfo(n)%dims(1:ndim)
      count(1:ndim,1) = dimsinfo(n)%count(1:ndim)
      call FILE_HISTORY_Set_Dim ( dimsinfo(n)%type, ndim, 1, dims(1:ndim,:), zs(:), start(1:ndim,:), count(1:ndim,:))
    end do
    
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_X)%name, dimsinfo(DIMTYPE_X)%desc, dimsinfo(DIMTYPE_X)%unit, dimsinfo(DIMTYPE_X)%name, x(:) )
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_Y)%name, dimsinfo(DIMTYPE_Y)%desc, dimsinfo(DIMTYPE_Y)%unit, dimsinfo(DIMTYPE_Y)%name, y(:) )
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_Z)%name, dimsinfo(DIMTYPE_Z)%desc, dimsinfo(DIMTYPE_Z)%unit, dimsinfo(DIMTYPE_Z)%name, z(:), &
      down=dimsinfo(DIMTYPE_Z)%positive_down(1) )
      
    return
  end subroutine set_dim_axis3D_cubedsphere

!----------------

end module scale_file_history_meshfield
