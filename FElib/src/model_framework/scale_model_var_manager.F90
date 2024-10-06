!-------------------------------------------------------------------------------
!> FElib / model framework / variable manager
!!
!! @par Description
!!          A module for managing data of variables used in models
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_model_var_manager
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io

  use scale_file_history, only: FILE_HISTORY_reg
  use scale_file_monitor_meshfield, only: FILE_monitor_meshfield_reg


  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField1D, MeshField2D, MeshField3D
  use scale_linkedlist, only: &
    LinkedList
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase, MeshFieldContainer
  
  use scale_variableinfo, only: &
    VariableInfo
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: ModelVarManager
    type(LinkedList) :: list
    class(MeshFieldCommBase), pointer :: ptr_comm
    type(MeshFieldContainer), allocatable :: comm_list(:)  
  contains
    procedure, public :: Init => ModelVarManager_Init
    procedure, public :: Final => ModelVarManager_Final
    procedure, private :: Regist1D => ModelVarManager_regist1D
    procedure, private :: Regist2D => ModelVarManager_regist2D
    procedure, private :: Regist3D => ModelVarManager_regist3D
    generic, public :: Regist => Regist1D, Regist2D, Regist3D
    procedure, public :: Get => ModelVarManager_Get
    procedure, public :: Get2D => ModelVarManager_Get2D
    procedure, public :: Get3D => ModelVarManager_Get3D
    procedure, public :: GetLocalMeshField => ModelVarManager_GetLocalMeshField
    procedure, public :: GetLocalMeshFieldList => ModelVarManager_GetLocalMeshFieldList

    !--
    procedure, public :: MeshFieldComm_Prepair => ModelVarManager_meshfiled_comm_prepare
    procedure, public :: MeshFieldComm_Exchange => ModelVarManager_meshfiled_comm_exchange
  end type ModelVarManager

  public :: VariableInfo ! Cascade

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !------------------

contains

!OCL SERIAL
  subroutine ModelVarManager_Init( this )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    !------------------------------------------------

    call this%list%Init()
    nullify( this%ptr_comm )

    return
  end subroutine ModelVarManager_Init

!OCL SERIAL
  subroutine ModelVarManager_Final( this )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    !------------------------------------------------

    call this%list%Traverse( field_release )

    call this%list%Final()

    nullify( this%ptr_comm )
    if ( allocated(this%comm_list) ) deallocate(this%comm_list)

    return
  end subroutine ModelVarManager_Final

!OCL SERIAL
  subroutine field_release( key, pField, done)
    implicit none
    class(*), intent(in) :: key
    class(*), pointer    :: pField
    logical, intent(out) :: done
    !------------------------------------------------

    select type( pField )
    type is ( MeshField1D )
      call pField%Final()
    type is ( MeshField2D )
      call pField%Final() 
    type is ( MeshField3D )
      call pField%Final()            
    end select

    return
  end subroutine field_release

!OCL SERIAL
  subroutine ModelVarManager_Regist1D( this,       &
      varinfo, mesh, field, reg_file_history_flag, &
      monitor_flag, fill_zero                      )
    
    use scale_mesh_base1d, only: Meshbase1d
    implicit none

    class(ModelVarManager), intent(inout) :: this
    type(VariableInfo), intent(in) :: varinfo
    class(MeshBase1D), intent(in) :: mesh
    class(MeshField1D), intent(inout), target :: field
    logical, intent(in) :: reg_file_history_flag
    logical, intent(in), optional :: monitor_flag
    logical, intent(in), optional :: fill_zero

    class(*), pointer :: ptr_field

    integer :: domID, ke
    logical :: fill_zero_
    !------------------------------------------------

    if (present(fill_zero)) then
      fill_zero_ = fill_zero
    else
      fill_zero_ = .true.
    end if

    call field%Init( varinfo%NAME, varinfo%UNIT, mesh )
    if (reg_file_history_flag) then
      call FILE_HISTORY_reg( varinfo%NAME, varinfo%DESC, varinfo%UNIT, field%hist_id, dim_type=varinfo%dim_type )
    end if
    if ( present(monitor_flag) ) then
      if (monitor_flag) then
        call FILE_monitor_meshfield_reg( &
          varinfo%NAME, varinfo%DESC, trim(varinfo%UNIT)//'*m',        & ! (in)
          field%monitor_id, dim_type='ATM1D', is_tendency=.false.      ) ! (in)
      end if
    end if
     
    ptr_field => field
    call this%list%AddByPointer( varinfo%keyID, ptr_field )

    if ( fill_zero_ ) then
      do domID=1, mesh%LOCAL_MESH_NUM
        !$omp parallel do
        do ke=mesh%lcmesh_list(domID)%NeS, mesh%lcmesh_list(domID)%NeE
          field%local(domID)%val(:,ke) = 0.0_RP
        end do
      end do
    end if

    return
  end subroutine ModelVarManager_Regist1D

!OCL SERIAL
  subroutine ModelVarManager_Regist2D( this,     &
    varinfo, mesh, field, reg_file_history_flag, &
    monitor_flag, fill_zero                      )
    use scale_mesh_base2d, only: Meshbase2d
    implicit none

    class(ModelVarManager), intent(inout) :: this
    type(VariableInfo), intent(in) :: varinfo
    class(MeshBase2D), intent(in) :: mesh
    class(MeshField2D), intent(inout), target :: field
    logical, intent(in) :: reg_file_history_flag
    logical, intent(in), optional :: monitor_flag
    logical, intent(in), optional :: fill_zero

    class(*), pointer :: ptr_field    

    integer :: domID, ke
    logical :: fill_zero_
    !------------------------------------------------

    if (present(fill_zero)) then
      fill_zero_ = fill_zero
    else
      fill_zero_ = .false.
    end if

    call field%Init( varinfo%NAME, varinfo%UNIT, mesh )
    if (reg_file_history_flag) then
      call FILE_HISTORY_reg( varinfo%NAME, varinfo%DESC, varinfo%UNIT, field%hist_id, dim_type=varinfo%dim_type)
    end if
    if ( present(monitor_flag) ) then
      if (monitor_flag) then
        call FILE_monitor_meshfield_reg( &
          varinfo%NAME, varinfo%DESC, trim(varinfo%UNIT)//'*m2',       & ! (in)
          field%monitor_id, dim_type='ATM2D', is_tendency=.false.      ) ! (in)
      end if
    end if

    ptr_field => field
    call this%list%AddByPointer( varinfo%keyID, ptr_field )

    if ( fill_zero_ ) then
      do domID=1, mesh%LOCAL_MESH_NUM
        !$omp parallel do
        do ke=mesh%lcmesh_list(domID)%NeS, mesh%lcmesh_list(domID)%NeE
          field%local(domID)%val(:,ke) = 0.0_RP
        end do
      end do
    end if

    return
  end subroutine ModelVarManager_Regist2D

!OCL SERIAL
  subroutine ModelVarManager_Regist3D( this,     &
    varinfo, mesh, field, reg_file_history_flag, &
    monitor_flag, fill_zero                      )
    use scale_mesh_base3d, only: MeshBase3D
    implicit none

    class(ModelVarManager), intent(inout) :: this
    type(VariableInfo), intent(in) :: varinfo
    class(MeshBase3D), intent(in) :: mesh
    class(MeshField3D), intent(inout), target :: field
    logical, intent(in) :: reg_file_history_flag
    logical, intent(in), optional :: monitor_flag
    logical, intent(in), optional :: fill_zero

    class(*), pointer :: ptr_field 
    
    integer :: domID, ke
    logical :: fill_zero_    
    !------------------------------------------------

    if (present(fill_zero)) then
      fill_zero_ = fill_zero
    else
      fill_zero_ = .false.
    end if

    call field%Init( varinfo%NAME, varinfo%UNIT, mesh )
    if (reg_file_history_flag) then
      call FILE_HISTORY_reg( field%varname, varinfo%DESC, field%unit, field%hist_id, &
                             dim_type=varinfo%dim_type )
    end if
    if ( present(monitor_flag) ) then
      if (monitor_flag) then
        call FILE_monitor_meshfield_reg( &
          varinfo%NAME, varinfo%DESC, trim(varinfo%UNIT)//'*m3',       & ! (in)
          field%monitor_id, dim_type='ATM3D', is_tendency=.false.      ) ! (in)
      end if
    end if

    ptr_field => field
    call this%list%AddByPointer( varinfo%keyID, ptr_field )

    if ( fill_zero_ ) then
      do domID=1, mesh%LOCAL_MESH_NUM
        !$omp parallel do
        do ke=mesh%lcmesh_list(domID)%NeS, mesh%lcmesh_list(domID)%NeE
          field%local(domID)%val(:,ke) = 0.0_RP
        end do
      end do
    end if
        
    return
  end subroutine ModelVarManager_Regist3D

!OCL SERIAL
  subroutine ModelVarManager_Get( this, keyID, pField )
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    class(ModelVarManager), intent(inout) :: this
    integer, intent(in) :: keyID
    class(MeshFieldBase), pointer, intent(out) :: pField

    class(*), pointer :: ptr_field    
    !------------------------------------------------

    call this%list%Get(keyID, ptr_field)

    nullify( pField )
    select type( ptr_field )
    type is (MeshField1D)
      pField => ptr_field
    type is (MeshField2D)
      pField => ptr_field
    type is (MeshField3D)
      pField => ptr_field          
    end select

    return
  end subroutine ModelVarManager_Get

!OCL SERIAL
  subroutine ModelVarManager_Get2D( this, keyID, pField )
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    class(ModelVarManager), intent(inout) :: this
    integer, intent(in) :: keyID
    class(MeshField2D), pointer, intent(out) :: pField

    class(*), pointer :: ptr_field    
    !------------------------------------------------

    call this%list%Get(keyID, ptr_field)
    
    nullify( pField )
    select type( ptr_field )
    type is (MeshField2D)
      pField => ptr_field          
    end select

    return
  end subroutine ModelVarManager_Get2D

!OCL SERIAL
  subroutine ModelVarManager_Get3D( this, keyID, pField )
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    class(ModelVarManager), intent(inout) :: this
    integer, intent(in) :: keyID
    class(MeshField3D), pointer, intent(out) :: pField

    class(*), pointer :: ptr_field    
    !------------------------------------------------

    call this%list%Get(keyID, ptr_field)

    nullify( pField )
    select type( ptr_field )
    type is (MeshField3D)
      pField => ptr_field          
    end select

    return
  end subroutine ModelVarManager_Get3D

!OCL SERIAL
  subroutine ModelVarManager_GetLocalMeshField( this, keyID, domID, pField_lc )
    use scale_meshfield_base, only: MeshFieldBase
    use scale_localmeshfield_base, only: LocalMeshFieldBase
    implicit none

    class(ModelVarManager), intent(inout) :: this
    integer, intent(in) :: keyID
    integer, intent(in) :: domID
    class(LocalMeshFieldBase), pointer, intent(out) :: pField_lc

    class(MeshFieldBase), pointer :: pField
    !------------------------------------------------

    call this%Get(keyID, pField)
    call pField%GetLocalMeshField(domID, pField_lc)

    return
  end subroutine ModelVarManager_GetLocalMeshField

!OCL SERIAL
  subroutine ModelVarManager_GetLocalMeshFieldList( this, keyID_list, domID, lcfield_list )
    use scale_meshfield_base, only: MeshFieldBase
    use scale_localmeshfield_base, only: LocalMeshFieldBaseList
    implicit none

    class(ModelVarManager), intent(inout) :: this
    integer, intent(in) :: keyID_list(:)
    integer, intent(in) :: domID
    type(LocalMeshFieldBaseList), intent(out) :: lcfield_list(size(keyID_list))

    integer :: i
    integer :: keyID

    class(MeshFieldBase), pointer :: pField
    !------------------------------------------------

    do i=1, size(keyID_list)
      keyID = keyID_list(i)
      call this%Get(keyID, pField)
      call pField%GetLocalMeshField( domID, lcfield_list(i)%ptr )
    end do

    return
  end subroutine ModelVarManager_GetLocalMeshFieldList
  
!OCL SERIAL
  subroutine ModelVarManager_meshfiled_comm_prepare( this, &
      comm,  fields )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    class(MeshFieldCommBase), target, intent(in) :: comm
    class(MeshFieldBase), target, intent(in) :: fields(:)

    integer :: v
    integer :: nFields
    class(MeshFieldBase), pointer :: ptr_field
    !------------------------------------------------

    this%ptr_comm => comm

    nFields = size(fields)
    allocate( this%comm_list(nFields) )

    do v = 1, size(fields)
      ptr_field => fields(v)
      select type( ptr_field )
      type is (MeshField1D)
        this%comm_list(v)%field1d => ptr_field
      type is (MeshField2D)
        this%comm_list(v)%field2d => ptr_field
      type is (MeshField3D)
        this%comm_list(v)%field3d => ptr_field               
      end select
    end do

    return
  end subroutine ModelVarManager_meshfiled_comm_prepare

!OCL SERIAL  
  subroutine ModelVarManager_meshfiled_comm_exchange( this )
    implicit none
    class(ModelVarManager), intent(inout) :: this
    !------------------------------------------------
    
    call this%ptr_comm%Put( this%comm_list, 1 )
    call this%ptr_comm%Exchange()
    call this%ptr_comm%Get( this%comm_list, 1 )

    return
  end subroutine ModelVarManager_meshfiled_comm_exchange

end module scale_model_var_manager