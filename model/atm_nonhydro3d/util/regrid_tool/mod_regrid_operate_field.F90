!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_regrid_operate_field
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


  use scale_file_base_meshfield, only: FILE_base_meshfield

  use mod_regrid_mesh_base, only: &
    regrid_mesh_base
  use mod_regrid_nodemap, only: &
    regrid_nodemap  
  use mod_regrid_outvar_info, only: &
    OutVarInfoList, OutVarInfo
  use mod_regrid_vec_conversion, only: &
    regrid_vec_conversion_Init,  &
    regrid_vec_conversion_Final

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  public :: regrid_operate_field_Init
  public :: regrid_operate_field_Final
  public :: regrid_operate_field_Do
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  type(OutVarInfoList), public, target :: out_vinfo

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  
  logical, private :: uvmet_conversion_flag = .false.
  integer :: vinfoid_umet = -1
  integer :: vinfoid_vmet = -1

contains

!OCL SERIAL
  subroutine regrid_operate_field_Init( out_mesh )
    use mod_regrid_outvar_info, only: &
      OUTVARINFO_ITEM_MAX_NUM
    use mod_regrid_interp_field, only: &
      regrid_interp_field_open_infile0
    
    implicit none
    class(regrid_mesh_base), intent(in), target :: out_mesh

    integer :: out_tinterval = 1

    namelist /PARAM_REGRID_OPERATE_FIELD/ &
      uvmet_conversion_flag, &
      out_tinterval
        
    integer :: ierr

    character(len=H_SHORT)  :: vars(OUTVARINFO_ITEM_MAX_NUM) = ''       ! name of variables
    integer :: out_tinterval_(OUTVARINFO_ITEM_MAX_NUM)

    integer :: vars_counter

    type(FILE_base_meshfield) :: in_file
    !----------------------------------------------

    !--- read namelist
    
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_REGRID_OPERATE_FIELD,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_operate_field",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_operate_field",*) 'Not appropriate names in namelist PARAM_REGRID_OPERATE_FIELD. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_REGRID_OPERATE_FIELD)

    vars_counter = 0
    if ( uvmet_conversion_flag ) then
      vinfoid_umet = vars_counter + 1
      vars(vinfoid_umet) = 'U'
      vinfoid_vmet = vars_counter + 2
      vars(vinfoid_vmet) = 'V'
      vars_counter = vars_counter + 2
    end if

    !--------------
    call regrid_interp_field_open_infile0( out_mesh, & ! (in)
      in_file )                                        ! (inout)
    
    out_tinterval_(:) = out_tinterval
    call out_vinfo%Init( vars, out_tinterval_, in_file )

    call in_file%Close()
    call in_file%Final()
    !--------------

    if ( uvmet_conversion_flag ) then
      call regrid_vec_conversion_Init( out_mesh%ptr_mesh3D )
      call out_vinfo%items(vinfoid_umet)%SetProperties( 'Umet', 'm/s', 'eastward velocity', 'x_wind' )
      call out_vinfo%items(vinfoid_vmet)%SetProperties( 'Vmet', 'm/s', 'westward velocity', 'y_wind' )
    end if

    return
  end subroutine regrid_operate_field_Init

!OCL SERIAL
  subroutine regrid_operate_field_do( out_mesh, nodeMap_list, vintrp )
    use mod_regrid_interp_vcoord, only: &
      regrid_interp_vcoord, REGRID_VCOORD_MODEL_ID
    use mod_regrid_file, only: &
      regrid_file_write_var
    use mod_regrid_vec_conversion, only: &
      regrid_vec_conversion_Do,    &  
      out_veclat, out_veclon    
    implicit none

    class(regrid_mesh_base), intent(in), target :: out_mesh
    type(regrid_nodemap), intent(in) :: nodeMap_list(:)
    type(regrid_interp_vcoord), intent(inout) :: vintrp

    integer :: istep
    type(OutVarInfo), pointer :: vinfo_u, vinfo_v
    !----------------------------------------------

    if ( uvmet_conversion_flag ) then
      vinfo_u => out_vinfo%items(vinfoid_umet)
      vinfo_v => out_vinfo%items(vinfoid_vmet)

      do istep=1, vinfo_u%num_step
        LOG_INFO("regrid_tool",'(a,i4)') ' operate_field :' // "Umet, Vmet" // " step=", istep

        call regrid_vec_conversion_Do( istep, out_mesh, nodeMap_list )

        if ( vintrp%vintrp_typeid == REGRID_VCOORD_MODEL_ID ) then
          call regrid_file_write_var( vinfo_u, out_veclon, istep )
          call regrid_file_write_var( vinfo_v, out_veclat, istep )
        else
          call vintrp%Update_weight( istep, out_mesh, nodeMap_list )

          call vintrp%Interpolate( istep, out_mesh%ptr_mesh3D, out_veclon )
          call regrid_file_write_var( vinfo_u, vintrp%vintrp_var3D, istep )

          call vintrp%Interpolate( istep, out_mesh%ptr_mesh3D, out_veclat )
          call regrid_file_write_var( vinfo_v, vintrp%vintrp_var3D, istep )
        end if
      end do
    end if

    return
  end subroutine regrid_operate_field_do

!OCL SERIAL
  subroutine regrid_operate_field_Final()
    implicit none
    !----------------------------------------------

    if ( uvmet_conversion_flag ) call regrid_vec_conversion_Final()

    return
  end subroutine regrid_operate_field_Final

end module mod_regrid_operate_field
