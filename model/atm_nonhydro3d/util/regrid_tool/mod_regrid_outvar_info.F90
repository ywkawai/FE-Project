!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_regrid_outvar_info
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_file_h

  use scale_file_base_meshfield, only: &
    FILE_base_meshfield
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !  
  type, public :: OutVarInfo
    character(FILE_HSHORT) :: varname
    character(FILE_HSHORT) :: units
    integer :: num_step
    real(DP) :: dt
    real(DP) :: start_sec
    integer :: out_tintrv
    character(len=FILE_HMID) :: standard_name
    character(len=FILE_HMID) :: desc

    integer :: vidForOutput
  contains
    procedure :: SetProperties => OutVarInfo_set_properties
  end type

  type, public :: OutVarInfoList
    type(OutVarInfo), allocatable :: items(:)
    integer :: item_num
  contains
    procedure :: Init => OutVarInfoList_Init
    procedure :: Final => OutVarInfoList_Final
    procedure :: DefVarForNetCDF => OutVarInfoList_DefVarForNetCDF
  end type

  integer, public, parameter :: OUTVARINFO_ITEM_MAX_NUM = 128

!
contains 

  subroutine OutVarInfo_set_properties( this, &
    varname, units, desc, standard_name )

    implicit none
    class(OutVarInfo), intent(inout) :: this
    character(*), intent(in) :: varname
    character(*), intent(in) :: units
    character(*), intent(in) :: desc
    character(*), intent(in) :: standard_name
    !------------------------------------

    this%varname = varname
    this%units   = units
    this%desc = desc
    this%standard_name = standard_name

    return
  end subroutine OutVarInfo_set_properties

  subroutine OutVarInfoList_Init( this,     &
      varname_list, out_tinterval,          &
      in_file   )

    use scale_const, only: &
      EPS => CONST_EPS
    
    implicit none
    class(OutVarInfoList), intent(inout) :: this
    character(len=H_SHORT), intent(in) :: varname_list(OUTVARINFO_ITEM_MAX_NUM)
    integer, intent(in) :: out_tinterval(OUTVARINFO_ITEM_MAX_NUM)
    type(FILE_base_meshfield), intent(in) :: in_file

    integer :: nn 
    integer :: out_var_num

    real(DP) :: time_endsec
    !------------------------------------

    out_var_num = 0
    do nn= 1, OUTVARINFO_ITEM_MAX_NUM
      if ( varname_list(nn) == '' ) then
        exit
      else
        out_var_num = out_var_num + 1
      end if
    end do

    allocate( this%items(out_var_num) )
    this%item_num = size(this%items)


    do nn = 1, out_var_num
      this%items(nn)%varname = varname_list(nn)
      call in_file%Get_dataInfo( this%items(nn)%varname, istep=1,  & ! (in)
        units=this%items(nn)%units,                                & ! (out)
        time_start=this%items(nn)%start_sec,                       & ! (out)
        time_end=time_endsec                                       ) ! (out)

      call in_file%Get_VarStepSize( this%items(nn)%varname, & ! (in)
        this%items(nn)%num_step                             ) ! (out)
      
      this%items(nn)%dt = time_endsec - this%items(nn)%start_sec
      this%items(nn)%out_tintrv = out_tinterval(nn)

      if (       abs(this%items(nn)%dt) < EPS &
           .and. this%items(nn)%num_step == 0 ) then
        this%items(nn)%num_step   = 1
        this%items(nn)%out_tintrv = 1
      end if

      this%items(nn)%standard_name = ''
      this%items(nn)%desc = ''

      LOG_INFO("regrid_OutVarInfoList_Init", '(3a,i4,a,i4)') &
        " Regist: name=", trim(this%items(nn)%varname ), &
        ", out_nstep=", this%items(nn)%num_step, &
        ", out_tinterval=", this%items(nn)%out_tintrv
    end do

    return
  end subroutine OutVarInfoList_Init

  subroutine OutVarInfoList_Final( this )
    implicit none
    class(OutVarInfoList), intent(inout) :: this
    !------------------------------------

    if ( allocated(this%items) ) deallocate( this%items)

    return
  end subroutine OutVarInfoList_Final  

  subroutine OutVarInfoList_DefVarForNetCDF( this, out_file, &
    dim_typeid, out_dtype, vid_offset )
    implicit none
    class(OutVarInfoList), intent(inout) :: this
    type(FILE_base_meshfield), intent(inout) :: out_file
    integer, intent(in) :: dim_typeid
    character(len=H_SHORT), intent(in)  :: out_dtype
    integer, intent(in) :: vid_offset

    integer :: nn
    !------------------------------------

    do nn=1, this%item_num
      this%items(nn)%vidForOutput = vid_offset + nn
      call out_file%Def_Var( this%items(nn)%varname, this%items(nn)%units, &
        this%items(nn)%desc, this%items(nn)%vidForOutput, dim_typeid, out_dtype,       &
        standard_name=this%items(nn)%standard_name,                                    &
        timeinv=this%items(nn)%dt * dble(this%items(nn)%out_tintrv)                    )

    end do

    return
  end subroutine  OutVarInfoList_DefVarForNetCDF

end module mod_regrid_outvar_info