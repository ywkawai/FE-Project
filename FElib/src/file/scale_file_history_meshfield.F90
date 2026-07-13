!-------------------------------------------------------------------------------
!> module FElib / File / History
!!
!! @par Description
!!          A module for managing file history 
!!
!! @author Yuta Kawai, Team SCALE
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
  use scale_prc, only: &
    PRC_abort

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
    FILE_HISTORY_query,       &
    FILE_HISTORY_finalize

  use scale_element_base, only: ElementBase1D, ElementBase2D, ElementBase3D
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
  !++ Public type and procedures
  !

  integer,          parameter :: nzs = 1

  type :: FileHistoryMeshFieldComp
    class(MeshBase1D), pointer :: mesh1D
    class(MeshRectDom2D), pointer :: mesh2D
    class(MeshCubeDom3D), pointer :: mesh3D
    class(MeshCubedSphereDom2D), pointer :: meshCubedSphere2D
    class(MeshCubedSphereDom3D), pointer :: meshCubedSphere3D

    integer :: dims1D_size(1)
    integer :: dims2D_size(2)
    integer :: dims3D_size(3,nzs)

    character(len=H_SHORT) :: dim_name_postfix
  contains
    procedure :: Init => FileHistoryMeshFieldComp_Init
    procedure :: Final => FileHistoryMeshFieldComp_Final
    procedure, private :: set_dim_axis1D
    procedure, private :: set_dim_axis2D
    procedure, private :: set_dim_axis3D
    procedure, private :: set_dim_axis2D_cubedsphere
    procedure, private :: set_dim_axis3D_cubedsphere
  end type FileHistoryMeshFieldComp

  integer, parameter :: COMPONENT_NUM_MAX = 32
  type(FileHistoryMeshFieldComp), target :: comp(COMPONENT_NUM_MAX)  !< Array of registered components

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
  private :: history_in_regvar

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !----------------------------------------------------------------------------

  integer :: comp_num  = 0                                 !< Number of registered components
  type(FileHistoryMeshFieldComp), pointer :: default_comp

  character(len=8), parameter :: zs(nzs) = (/  "model   " /)

  integer  :: FILE_HISTORY_MESHFIELD_STARTDATE(6) !< start time [YYYY MM DD HH MM SS]
  real(DP) :: FILE_HISTORY_MESHFIELD_STARTSUBSEC  !< subsecond part of start time [millisec]


  logical, private :: FILE_HISTORY_FILEMESHFILED_DISABLE = .true.

contains

!----------------
!OCL SERIAL
  subroutine FileHistoryMeshFieldComp_Init( this, &
    mesh1D_, mesh2D_, mesh3D_,                    &
    meshcubedsphere2D_, meshcubedsphere3D_,       & 
    dim_name_postfix )
    implicit none
    class(FileHistoryMeshFieldComp), intent(inout) :: this

    class(MeshBase1d), intent(in), target, optional :: mesh1D_     !< An object of 1D mesh when the computational domain is 1D
    class(MeshRectDom2d), intent(in), target, optional :: mesh2D_  !< An object of 2D mesh when the computational domain is rectangular
    class(MeshCubeDom3D), intent(in), target, optional :: mesh3D_  !< An object of 3D mesh when the computational domain is cubed
    class(MeshCubedSphereDom2D), intent(in), target, optional :: meshCubedsphere2D_ !< An object of 2D mesh when the computational domain is 2D cubed-sphere
    class(MeshCubedSphereDom3D), intent(in), target, optional :: meshCubedsphere3D_ !< An object of 3D mesh when the computational domain is 3D cubed-sphere
    character(len=*), intent(in), optional :: dim_name_postfix     !< Postfix for dimension names
    !-----------------------------------------------

    if ( present(dim_name_postfix) ) then
      this%dim_name_postfix = dim_name_postfix
    else
      this%dim_name_postfix = ''
    end if

    !- Set a pointer to the  variable of mesh 
    !  and set the dimension information for each axis
    
    nullify( this%mesh1D, this%mesh2D, this%mesh3D )
    nullify( this%meshCubedSphere2D, this%meshCubedsphere3D )

    if ( present(mesh1D_) ) then
      this%mesh1D => mesh1D_
      call this%set_dim_axis1D()
    else if ( present(mesh2D_) ) then
      this%mesh2D => mesh2D_
      call this%set_dim_axis2D()
    else if ( present(mesh3D_) ) then
      this%mesh3D => mesh3D_
      this%mesh2D => mesh3D_%mesh2D
      call this%set_dim_axis3D()
    else if ( present(meshCubedsphere2D_) ) then
      this%meshCubedSphere2D => meshCubedsphere2D_
      call this%set_dim_axis2D_cubedsphere()
    else if ( present(meshCubedsphere3D_) ) then
      this%meshCubedSphere3D => meshCubedsphere3D_
      this%meshCubedSphere2D => meshCubedSphere3D_%mesh2D
      call this%set_dim_axis3D_cubedsphere()
    else
      LOG_ERROR("FileHistoryMeshFieldComp_Init",*)   "Any mesh (mesh1d/2d/3d) are not specified."
      call PRC_abort
    end if

    return
  end subroutine FileHistoryMeshFieldComp_Init

!OCL SERIAL
  subroutine FileHistoryMeshFieldComp_Final( this )
    implicit none
    class(FileHistoryMeshFieldComp), intent(inout) :: this
    !-----------------------------------------------
    return
  end subroutine FileHistoryMeshFieldComp_Final

  !> Setup a module for file history
!OCL SERIAL
  subroutine FILE_HISTORY_meshfield_setup( &
    mesh1D_, mesh2D_, mesh3D_,              &
    meshcubedsphere2D_, meshcubedsphere3D_, &
    dim_name_postfix_,                      &
    registered_comp_id                      )

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

    class(MeshBase1d), intent(in), target, optional :: mesh1D_     !< An object of 1D mesh when the computational domain is 1D
    class(MeshRectDom2d)       , intent(in), target, optional :: mesh2D_            !< An object of 2D mesh when the computational domain is rectangular
    class(MeshCubeDom3D)       , intent(in), target, optional :: mesh3D_            !< An object of 3D mesh when the computational domain is cubed
    class(MeshCubedSphereDom2D), intent(in), target, optional :: meshCubedsphere2D_ !< An object of 2D mesh when the computational domain is 2D cubed-sphere
    class(MeshCubedSphereDom3D), intent(in), target, optional :: meshCubedsphere3D_ !< An object of 3D mesh when the computational domain is 3D cubed-sphere
    character(len=*), intent(in), optional :: dim_name_postfix_    !< Postfix for dimension names
    integer, intent(out), optional :: registered_comp_id

    character(len=H_MID) :: FILE_HISTORY_MESHFILED_H_TITLE = 'SCALE-DG FILE_HISTORY_MESHFIELD' !< title of the output file
    character(len=H_MID) :: FILE_HISTORY_MESHFIELD_T_SINCE
    
    character(len=FILE_HSHORT) :: calendar
    real(DP) :: start_daysec

    class(FileHistoryMeshFieldComp), pointer :: target_comp
    integer :: registered_comp_id_
    !---------------------------------------------------------------------------

    !-
    if ( comp_num >= COMPONENT_NUM_MAX ) then
      LOG_ERROR("FILE_base_meshfield_register_comp",*) 'Exceeding maximum number of components. Check!'
      call PRC_abort
    end if

    comp_num = comp_num + 1
    registered_comp_id_ = comp_num
    if ( present(registered_comp_id) ) then
      registered_comp_id = registered_comp_id_
    end if

    target_comp => comp(registered_comp_id_)
    call target_comp%Init( mesh1D_, mesh2D_, mesh3D_, &
      meshCubedsphere2D_, meshCubedsphere3D_ ,        &
      dim_name_postfix_ )

    !-
    if ( registered_comp_id_ == 1 ) then
      default_comp => comp(registered_comp_id_)

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

      FILE_HISTORY_FILEMESHFILED_DISABLE = .false. 
    end if

    return
  end subroutine FILE_HISTORY_meshfield_setup

  !> Write history data to the file
!OCL SERIAL
  subroutine FILE_HISTORY_meshfield_write()
    implicit none
    !-------------------------------------------------

    call FILE_HISTORY_write
    return
  end subroutine FILE_HISTORY_meshfield_write

  !> Finalize the file history module
!OCL SERIAL
  subroutine FILE_HISTORY_meshfield_finalize()
    implicit none
    !-------------------------------------------------

    call FILE_HISTORY_finalize()
    call default_comp%Final()
    return
  end subroutine FILE_HISTORY_meshfield_finalize

  !-- 1D

!OCL SERIAL
  subroutine FILE_HISTORY_meshfield_put1D(hstid, field1d, &
    registered_comp_id )
    use scale_file_common_meshfield, only: &
      File_common_meshfield_put_field1D_cartesbuf
    implicit none
    integer, intent(in) :: hstid
    class(MeshField1D), intent(in) :: field1d
    integer, intent(in), optional :: registered_comp_id

    logical :: do_put
    integer :: ldomID
    real(RP), allocatable :: buf(:)

    class(FileHistoryMeshFieldComp), pointer :: target_comp
    !---------------------------------------------------------------------------

    call history_get_comp( registered_comp_id, & ! (in)
      target_comp ) ! (out)
    
    call FILE_HISTORY_query( hstid, do_put )
    if ( .not. do_put ) return

    !-
#ifdef _OPENACC
    do ldomID=1,target_comp%mesh1D%LOCAL_MESH_NUM
      !$acc update host( field1d%local(ldomID)%val ) async(1)
    end do
    !$acc wait(1)
#endif    
    !-
    allocate( buf(target_comp%dims1D_size(1)) )

    call File_common_meshfield_put_field1D_cartesbuf( target_comp%mesh1D, field1d, buf )
    call FILE_HISTORY_put( hstid, buf )

    return
  end subroutine FILE_HISTORY_meshfield_put1D

!OCL SERIAL
  subroutine FILE_HISTORY_meshfield_in1D( field1d, desc, standard_name, &
    registered_comp_id )
    implicit none
    class(MeshField1D), intent(in) :: field1d
    character(len=*), intent(in) :: desc                     !< description of the item
    character(len=*), intent(in), optional :: standard_name
    integer, intent(in), optional :: registered_comp_id

    integer :: hstid
    logical :: do_put
    !-------------------------------------------------

    call history_in_regvar( hstid, do_put,    & ! (out)
      field1d, desc, 1, standard_name, 'XYZ'  ) ! (in)

    if ( do_put ) call FILE_HISTORY_meshfield_put( hstid, field1d, registered_comp_id )

    return
  end subroutine FILE_HISTORY_meshfield_in1D
  
  !-- 2D

!OCL SERIAL
  subroutine FILE_HISTORY_meshfield_put2D(hstid, field2d, &
    registered_comp_id )
    use scale_file_common_meshfield, only: &
      File_common_meshfield_put_field2D_cartesbuf,            &
      File_common_meshfield_put_field2D_cubedsphere_cartesbuf

    implicit none
    integer, intent(in) :: hstid
    class(MeshField2D), intent(in) :: field2d
    integer, intent(in), optional :: registered_comp_id

    logical :: do_put
    integer :: ldomID
    real(RP), allocatable :: buf(:,:)
    
    class(FileHistoryMeshFieldComp), pointer :: target_comp
    !---------------------------------------------------------------------------

    call history_get_comp( registered_comp_id, & ! (in)
      target_comp ) ! (out)
    
    call FILE_HISTORY_query( hstid, do_put )
    if ( .not. do_put ) return

    !-
#ifdef _OPENACC
    do ldomID=1, target_comp%mesh2D%LOCAL_MESH_NUM
      !$acc update host( field2d%local(ldomID)%val ) async(1)
    end do
    !$acc wait(1)
#endif    
    allocate( buf(target_comp%dims2D_size(1),target_comp%dims2D_size(2)) )
    
    if ( associated(target_comp%mesh2D) ) then
      call File_common_meshfield_put_field2D_cartesbuf( target_comp%mesh2D, field2d, buf )
    else if ( associated(target_comp%meshCubedSphere2D) ) then
      call File_common_meshfield_put_field2D_cubedsphere_cartesbuf( &
        target_comp%meshCubedSphere2D, field2d, buf )
    end if
    call FILE_HISTORY_put( hstid, buf )

    return
  end subroutine FILE_HISTORY_meshfield_put2D

!OCL SERIAL
  subroutine FILE_HISTORY_meshfield_in2D( field2d, desc, standard_name, &
    registered_comp_id )
    implicit none
    class(MeshField2D), intent(in) :: field2d
    character(len=*), intent(in) :: desc                     !< description of the item
    character(len=*), intent(in), optional :: standard_name
    integer, intent(in), optional :: registered_comp_id

    integer :: hstid
    logical :: do_put
    !-------------------------------------------------

    call history_in_regvar( hstid, do_put,  & ! (out)
      field2d, desc, 2, standard_name, 'XY' ) ! (in)

    if ( do_put ) call FILE_HISTORY_meshfield_put( hstid, field2d, registered_comp_id )

    return
  end subroutine FILE_HISTORY_meshfield_in2D

  !-- 3D

!OCL SERIAL
  subroutine FILE_HISTORY_meshfield_put3D(hstid, field3d, &
    registered_comp_id )
    use scale_file_common_meshfield, only: &
      File_common_meshfield_put_field3D_cartesbuf,            &
      File_common_meshfield_put_field3D_cubedsphere_cartesbuf

    implicit none
    integer, intent(in) :: hstid
    class(MeshField3D), intent(in) :: field3d
    integer, intent(in), optional :: registered_comp_id

    logical :: do_put
    integer :: ldomID
    real(RP), allocatable :: buf(:,:,:)

    class(FileHistoryMeshFieldComp), pointer :: target_comp
    !---------------------------------------------------------------------------

    call history_get_comp( registered_comp_id, & ! (in)
      target_comp ) ! (out)

    call FILE_HISTORY_query( hstid, do_put )
    if ( .not. do_put ) return

    !-
#ifdef _OPENACC
    do ldomID=1, target_comp%mesh3D%LOCAL_MESH_NUM
      !$acc update host( field3d%local(ldomID)%val ) async(1)
    end do
    !$acc wait(1)
#endif

    allocate( buf(target_comp%dims3D_size(1,1),target_comp%dims3D_size(2,1),target_comp%dims3D_size(3,1)) )

    if ( associated(target_comp%mesh3D) ) then
      call File_common_meshfield_put_field3D_cartesbuf( target_comp%mesh3D, field3d, buf )
    else if ( associated(target_comp%meshCubedSphere3D) ) then
      call File_common_meshfield_put_field3D_cubedsphere_cartesbuf( &
        target_comp%meshCubedSphere3D, field3d, buf )
    end if
    call FILE_HISTORY_put( hstid, buf )

    return
  end subroutine FILE_HISTORY_meshfield_put3D
  
!OCL SERIAL
  subroutine FILE_HISTORY_meshfield_in3D( field3d, desc, standard_name, &
    registered_comp_id )
    implicit none
    class(MeshField3D), intent(in) :: field3d
    character(len=*), intent(in) :: desc                     !< description of the item
    character(len=*), intent(in), optional :: standard_name
    integer, intent(in), optional :: registered_comp_id

    integer :: hstid
    logical :: do_put
    !-------------------------------------------------

    call history_in_regvar( hstid, do_put,   & ! (out)
      field3d, desc, 3, standard_name, 'XYZ' ) ! (in)

    if ( do_put ) call FILE_HISTORY_meshfield_put( hstid, field3d, registered_comp_id )

    return
  end subroutine FILE_HISTORY_meshfield_in3D

!- Private subroutines ---------------

!OCL SERIAL
  subroutine history_in_regvar( hstid, do_put,  &
    field, desc, ndim, standard_name, dim_type  )

    use scale_file_history, only: &
      FILE_HISTORY_reg
    
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

!OCL SERIAL
  subroutine history_get_comp( comp_id, & ! (in)
    comp_ptr ) ! (out)
    implicit none
    integer, intent(in), optional :: comp_id
    type(FileHistoryMeshFieldComp), intent(out), pointer :: comp_ptr
    !--------------------------------------------------------------

    if ( .not. present(comp_id) ) then
      comp_ptr => default_comp
      return
    end if

    if ( comp_id < 1 .or. comp_id > comp_num ) then
      LOG_ERROR("FILE_history_get_comp",*) 'Invalid component ID. Check!'
      call PRC_abort
    end if
    comp_ptr => comp(comp_id)
    return
  end subroutine history_get_comp

!OCL SERIAL
  subroutine set_dim_axis1D( comp )
    use scale_file_common_meshfield, only: &
      FILE_common_meshfield_diminfo,    &
      File_common_meshfield_get_dims1D, &
      File_common_meshfield_get_axis1D
    use scale_mesh_base1d, only: &
      DIMTYPE_NUM => MeshBase1D_DIMTYPE_NUM, &
      DIMTYPE_X   => MeshBase1D_DIMTYPEID_X
  
    implicit none
    class(FileHistoryMeshFieldComp), intent(inout) :: comp

    type(FILE_common_meshfield_diminfo) :: dimsinfo(DIMTYPE_NUM)
    real(RP), allocatable :: x(:)
    integer :: start(1,1), count(1,1)
    character(len=H_SHORT) :: dims(1,1)
    integer :: n, ndim
    !-------------------------------------------------
    
    call File_common_meshfield_get_dims1D( comp%mesh1D, comp%dim_name_postfix, & ! (in)
      dimsinfo    )                                                              ! (out)
    
    comp%dims1D_size(1) = dimsinfo(DIMTYPE_X)%size
    allocate( x(comp%dims1D_size(1)) )
    call File_common_meshfield_get_axis1D( comp%mesh1D, dimsinfo, & ! (in)
      x )                                                           ! (out)

    start(:,:) = 1
    do n=1, DIMTYPE_NUM
      ndim = dimsinfo(n)%ndim
      dims(1:ndim,1)  = dimsinfo(n)%dims(1:ndim)
      count(1:ndim,1) = dimsinfo(n)%count(1:ndim)
      call FILE_HISTORY_Set_Dim ( dimsinfo(n)%type, ndim, 1, dims(1:ndim,:), zs(:), start(1:ndim,:), count(1:ndim,:))
    end do
    
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_X)%name, dimsinfo(DIMTYPE_X)%desc, &
      dimsinfo(DIMTYPE_X)%unit, dimsinfo(DIMTYPE_X)%name, x )

    return
  end subroutine set_dim_axis1D

!OCL SERIAL
  subroutine set_dim_axis2D( comp )
    use scale_file_common_meshfield, only: &
      FILE_common_meshfield_diminfo,       &
      File_common_meshfield_get_dims2D,    &
      File_common_meshfield_get_axis2D
    use scale_mesh_base2d, only: &
      DIMTYPE_NUM => MeshBase2D_DIMTYPE_NUM, &
      DIMTYPE_X   => MeshBase2D_DIMTYPEID_X, &
      DIMTYPE_Y   => MeshBase2D_DIMTYPEID_Y
   
    implicit none
    class(FileHistoryMeshFieldComp), intent(inout) :: comp

    type(FILE_common_meshfield_diminfo) :: dimsinfo(DIMTYPE_NUM)
    real(RP), allocatable :: x(:), y(:)
    integer :: start(2,1), count(2,1)
    character(len=H_SHORT) :: dims(2,1)
    integer :: n, ndim
    !-------------------------------------------------
    
    call File_common_meshfield_get_dims2D( comp%mesh2D, comp%dim_name_postfix, & ! (in)
      dimsinfo    )                                                              ! (out)
    
    comp%dims2D_size(1) = dimsinfo(DIMTYPE_X)%size
    comp%dims2D_size(2) = dimsinfo(DIMTYPE_Y)%size
    allocate( x(comp%dims2D_size(1)), y(comp%dims2D_size(2)) )
    call File_common_meshfield_get_axis2D( comp%mesh2D, dimsinfo, & ! (in)
      x, y )                                                        ! (out)

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

!OCL SERIAL
  subroutine set_dim_axis3D( comp )
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
    class(FileHistoryMeshFieldComp), intent(inout) :: comp

    type(FILE_common_meshfield_diminfo) :: dimsinfo(DIMTYPE_NUM)
    real(RP), allocatable :: x(:), y(:), z(:)
    integer :: start(3,1), count(3,1)
    character(len=H_SHORT) :: dims(3,1)
    integer :: n, ndim
    !-------------------------------------------------
    
    call File_common_meshfield_get_dims3D( comp%mesh3D, comp%dim_name_postfix, & ! (in)
      dimsinfo    )                                                              ! (out)
  
    comp%dims2D_size(1) = dimsinfo(DIMTYPE_X)%size
    comp%dims2D_size(2) = dimsinfo(DIMTYPE_Y)%size  
    comp%dims3D_size(1,1) = dimsinfo(DIMTYPE_X)%size
    comp%dims3D_size(2,1) = dimsinfo(DIMTYPE_Y)%size
    comp%dims3D_size(3,1) = dimsinfo(DIMTYPE_Z)%size
    allocate( x(comp%dims3D_size(1,1)), y(comp%dims3D_size(2,1)), z(comp%dims3D_size(3,1)) )
    call File_common_meshfield_get_axis3D( comp%mesh3D, dimsinfo, & ! (in)
      x, y, z )                                                     ! (out)

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

!OCL SERIAL
  subroutine set_dim_axis2D_cubedsphere( comp )
    use scale_file_common_meshfield, only: &
      FILE_common_meshfield_diminfo,  &
      File_common_meshfield_get_dims, &
      File_common_meshfield_get_axis
    use scale_mesh_base2d, only: &
      DIMTYPE_NUM => MeshBase2D_DIMTYPE_NUM, &
      DIMTYPE_X  => MeshBase2D_DIMTYPEID_X,  &
      DIMTYPE_Y  => MeshBase2D_DIMTYPEID_Y

    implicit none
    class(FileHistoryMeshFieldComp), intent(inout) :: comp

    type(FILE_common_meshfield_diminfo) :: dimsinfo(DIMTYPE_NUM)
    real(RP), allocatable :: x(:), y(:)
    integer :: start(2,1), count(2,1)
    character(len=H_SHORT) :: dims(2,1)
    integer :: n, ndim
    !-------------------------------------------------
    
    call File_common_meshfield_get_dims( comp%meshCubedSphere2D, comp%dim_name_postfix, & ! (in)
      dimsinfo    )                                                                       ! (out)
    
    comp%dims2D_size(1) = dimsinfo(DIMTYPE_X)%size
    comp%dims2D_size(2) = dimsinfo(DIMTYPE_Y)%size
    allocate( x(comp%dims2D_size(1)), y(comp%dims2D_size(2)) )
    
    call File_common_meshfield_get_axis( comp%meshCubedSphere2D, dimsinfo, & ! (in)
      x, y )                                                                 ! (out)
    
    start(:,:) = 1

    do n=1, DIMTYPE_NUM
      ndim = dimsinfo(n)%ndim
      dims(1:ndim,1)  = dimsinfo(n)%dims(1:ndim)
      count(1:ndim,1) = dimsinfo(n)%count(1:ndim)
      call FILE_HISTORY_Set_Dim ( dimsinfo(n)%type, ndim, 1, dims(1:ndim,:), zs(:), start(1:ndim,:), count(1:ndim,:))
    end do
    
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_X)%name, dimsinfo(DIMTYPE_X)%desc, dimsinfo(DIMTYPE_X)%unit, dimsinfo(DIMTYPE_X)%name, x )
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_Y)%name, dimsinfo(DIMTYPE_Y)%desc, dimsinfo(DIMTYPE_Y)%unit, dimsinfo(DIMTYPE_Y)%name, y )
    
    return
  end subroutine set_dim_axis2D_cubedsphere

!OCL SERIAL
  subroutine set_dim_axis3D_cubedsphere( comp )
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
    class(FileHistoryMeshFieldComp), intent(inout) :: comp

    type(FILE_common_meshfield_diminfo) :: dimsinfo(DIMTYPE_NUM)
    real(RP), allocatable :: x(:), y(:), z(:)
    integer :: start(3,1), count(3,1)
    character(len=H_SHORT) :: dims(3,1)
    integer :: n, ndim
    !-------------------------------------------------
    
    call File_common_meshfield_get_dims( comp%meshCubedSphere3D, comp%dim_name_postfix, & ! (in)
      dimsinfo    )                                                                       ! (out)

    comp%dims2D_size(1) = dimsinfo(DIMTYPE_X)%size
    comp%dims2D_size(2) = dimsinfo(DIMTYPE_Y)%size
    comp%dims3D_size(1,1) = dimsinfo(DIMTYPE_X)%size
    comp%dims3D_size(2,1) = dimsinfo(DIMTYPE_Y)%size
    comp%dims3D_size(3,1) = dimsinfo(DIMTYPE_Z)%size
    allocate( x(comp%dims3D_size(1,1)), y(comp%dims3D_size(2,1)), z(comp%dims3D_size(3,1)) )
    call File_common_meshfield_get_axis( comp%meshCubedSphere3D, dimsinfo, & ! (in)
      x, y, z )                                                              ! (out)

    start(:,:) = 1
    do n=1, DIMTYPE_NUM
      ndim = dimsinfo(n)%ndim
      dims(1:ndim,1)  = dimsinfo(n)%dims(1:ndim)
      count(1:ndim,1) = dimsinfo(n)%count(1:ndim)
      call FILE_HISTORY_Set_Dim ( dimsinfo(n)%type, ndim, 1, dims(1:ndim,:), zs(:), start(1:ndim,:), count(1:ndim,:))
    end do
    
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_X)%name, dimsinfo(DIMTYPE_X)%desc, dimsinfo(DIMTYPE_X)%unit, dimsinfo(DIMTYPE_X)%name, x )
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_Y)%name, dimsinfo(DIMTYPE_Y)%desc, dimsinfo(DIMTYPE_Y)%unit, dimsinfo(DIMTYPE_Y)%name, y )
    call FILE_HISTORY_Set_Axis( dimsinfo(DIMTYPE_Z)%name, dimsinfo(DIMTYPE_Z)%desc, dimsinfo(DIMTYPE_Z)%unit, dimsinfo(DIMTYPE_Z)%name, z, &
      down=dimsinfo(DIMTYPE_Z)%positive_down(1) )
      
    return
  end subroutine set_dim_axis3D_cubedsphere

!----------------

end module scale_file_history_meshfield
