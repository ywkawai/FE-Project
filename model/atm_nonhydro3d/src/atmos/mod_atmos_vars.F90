!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_element_base, only: ElementBase3D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfield_base, only: MeshField3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  
  use scale_file_restart_meshfield, only: &
    FILE_restart_meshfield_component
  use scale_file_common_meshfield, only: &
    DIMTYPE_XYZ  => FILE_COMMON_MESHFILED3D_DIMTYPEID_XYZ
  
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
  use scale_meshfieldcomm_base, only: MeshFieldContainer

  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo

  use mod_atmos_mesh, only: AtmosMesh
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: AtmosVars
    type(MeshField3D), allocatable :: PROG_VARS(:)
    type(ModelVarManager) :: PROGVARS_manager
    type(MeshFieldCommCubeDom3D) :: PROGVARS_comm
    
    type(MeshField3D), allocatable :: AUX_VARS(:)
    type(ModelVarManager) :: AUXVARS_manager 
    
    type(ModelVarManager) :: DIAGVARS_manager    
    integer, allocatable :: DIAGVARS_HISTID(:)

    type(FILE_restart_meshfield_component) :: restart_file
  contains
    procedure :: Init => AtmosVars_Init
    procedure :: Final => AtmosVars_Final
    procedure :: History => AtmosVars_History
    procedure :: Read_restart_file => AtmosVar_Read_restart_file
    procedure :: Write_restart_file => AtmosVar_Write_restart_file
  end type AtmosVars

  public :: AtmosVars_GetLocalMeshField
  public :: AtmosVars_GetLocalMeshFields

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: ATMOS_PROGVARS_NUM = 5
  integer, public, parameter :: ATMOS_PROGVARS_DDENS_ID = 1
  integer, public, parameter :: ATMOS_PROGVARS_MOMX_ID  = 2
  integer, public, parameter :: ATMOS_PROGVARS_MOMY_ID  = 3
  integer, public, parameter :: ATMOS_PROGVARS_MOMZ_ID  = 4
  integer, public, parameter :: ATMOS_PROGVARS_DRHOT_ID = 5

  type(VariableInfo), public :: ATMOS_PROGVARS_VINFO(ATMOS_PROGVARS_NUM)

  DATA ATMOS_PROGVARS_VINFO / &
    VariableInfo( ATMOS_PROGVARS_DDENS_ID, 'DDENS', 'deviation of density',       &
                  'kg/m3',  3, 'XYZ',  'air_density'                           ), &
    VariableInfo( ATMOS_PROGVARS_MOMX_ID , 'MOMX', 'momentum x',                  &
                  'kg/m2/s', 3, 'XYZ', 'upward_mass_flux_of_air'               ), &
    VariableInfo( ATMOS_PROGVARS_MOMY_ID , 'MOMY', 'momentum y',                  &
                  'kg/m2/s', 3, 'XYZ', 'eastward_mass_flux_of_air'             ), &
    VariableInfo( ATMOS_PROGVARS_MOMZ_ID , 'MOMZ', 'momentum z',                  &
                  'kg/m2/s', 3, 'XYZ', 'northward_mass_flux_of_air'            ), &
    VariableInfo( ATMOS_PROGVARS_DRHOT_ID, 'DRHOT', 'deviation of rho * theta',   &
                  'kg/m3*K', 3, 'XYZ',  ''                                    )   /

  integer, public, parameter :: ATMOS_AUXVARS_NUM          = 2
  integer, public, parameter :: ATMOS_AUXVARS_PRESHYDRO_ID = 1
  integer, public, parameter :: ATMOS_AUXVARS_DENSHYDRO_ID = 2

  type(VariableInfo), public :: ATMOS_AUXVARS_VINFO(ATMOS_AUXVARS_NUM)
  DATA ATMOS_AUXVARS_VINFO / &
    VariableInfo( ATMOS_AUXVARS_PRESHYDRO_ID, 'PRES_hyd', 'hydrostatic part of pressure',  &
                  'Pa',  3, 'XYZ',  ''                                                  ), &
    VariableInfo( ATMOS_AUXVARS_DENSHYDRO_ID , 'DENS_hyd', 'hydrostatic part of density',  &
                  'kg/m3', 3, 'XYZ', ''                                                 )  /
  
  integer, public, parameter :: ATMOS_DIAGVARS_NUM       = 6
  integer, public, parameter :: ATMOS_DIAGVARS_U_ID      = 1
  integer, public, parameter :: ATMOS_DIAGVARS_V_ID      = 2
  integer, public, parameter :: ATMOS_DIAGVARS_W_ID      = 3
  integer, public, parameter :: ATMOS_DIAGVARS_PRES_ID   = 4
  integer, public, parameter :: ATMOS_DIAGVARS_T_ID      = 5  
  integer, public, parameter :: ATMOS_DIAGVARS_THETA_ID  = 6

  type(VariableInfo), public :: ATMOS_DIAGVARS_VINFO(ATMOS_DIAGVARS_NUM)
  DATA ATMOS_DIAGVARS_VINFO / &
    VariableInfo( ATMOS_DIAGVARS_U_ID, 'U', 'velocity u', 'm/s', 3, 'XYZ', 'x_wind'             ), &
    VariableInfo( ATMOS_DIAGVARS_V_ID, 'V', 'velocity v', 'm/s', 3, 'XYZ', 'y_wind'             ), &  
    VariableInfo( ATMOS_DIAGVARS_W_ID, 'W', 'velocity w', 'm/s', 3, 'XYZ', 'upward_air_velocity'), &
    VariableInfo( ATMOS_DIAGVARS_PRES_ID, 'PRES', 'pressure', 'Pa', 3, 'XYZ', 'air_pressure'    ), &
    VariableInfo( ATMOS_DIAGVARS_T_ID    , 'T', 'temperature', 'K', 3, 'XYZ', 'air_temperature'  ), &
    VariableInfo( ATMOS_DIAGVARS_THETA_ID, 'PT', 'potential temperature', 'K', 3, 'XYZ', 'potential_temperature'  ) /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains
  subroutine AtmosVars_Init( this, atm_mesh )
    implicit none
    class(AtmosVars), target, intent(inout) :: this
    class(AtmosMesh), intent(in) :: atm_mesh

    integer :: v
    integer :: n
    logical :: reg_file_hist

    type(MeshField3D) :: diag_vars(ATMOS_DIAGVARS_NUM)

    character(len=H_LONG) :: IN_BASENAME           = ''        !< Basename of the input  file
    logical :: IN_POSTFIX_TIMELABEL                = .false.   !< Add timelabel to the basename of input  file?
    character(len=H_LONG) :: OUT_BASENAME          = ''        !< Basename of the output file
    logical :: OUT_POSTFIX_TIMELABEL               = .true.    !< Add timelabel to the basename of output file?
    character(len=H_MID) :: OUT_TITLE              = ''        !< Title    of the output file
    character(len=H_SHORT) :: OUT_DTYPE            = 'DEFAULT' !< REAL4 or REAL8  

    namelist / PARAM_ATMOS_VARS_RESTART / &
      IN_BASENAME,           &
      IN_POSTFIX_TIMELABEL,  &
      OUT_BASENAME,          &
      OUT_POSTFIX_TIMELABEL, &
      OUT_TITLE,             &
      OUT_DTYPE    
    
    integer :: ierr
    logical :: is_specified
    !--------------------------------------------------

    LOG_INFO('AtmosVars_Init',*)

    !- Initialize prognostic variables

    call this%PROGVARS_manager%Init()
    allocate( this%PROG_VARS(ATMOS_PROGVARS_NUM) )

    reg_file_hist = .true.    
    do v = 1, ATMOS_PROGVARS_NUM
      call this%PROGVARS_manager%Regist(        &
        ATMOS_PROGVARS_VINFO(v), atm_mesh%mesh, & ! (in) 
        this%PROG_VARS(v), reg_file_hist        ) ! (out)
      do n = 1, atm_mesh%mesh%LOCAL_MESH_NUM
        this%PROG_VARS(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    call this%PROGVARS_comm%Init(ATMOS_PROGVARS_NUM, 0, atm_mesh%mesh)
    call this%PROGVARS_manager%MeshFieldComm_Prepair( this%PROGVARS_comm, this%PROG_VARS(:) )

    !- Initialize auxiliary variables

    call this%AUXVARS_manager%Init()
    allocate( this%AUX_VARS(ATMOS_AUXVARS_NUM) )
    
    reg_file_hist = .true.
    do v = 1, ATMOS_AUXVARS_NUM
      call this%AUXVARS_manager%Regist(        &
        ATMOS_AUXVARS_VINFO(v), atm_mesh%mesh, & ! (in) 
        this%AUX_VARS(v), reg_file_hist        ) ! (out)
      do n = 1, atm_mesh%mesh%LOCAL_MESH_NUM
        this%AUX_VARS(v)%local(n)%val(:,:) = 1.0_RP
      end do             
    end do

    !- Initialize diagnostic variables for output
    call this%DIAGVARS_manager%Init()
    allocate( this%DIAGVARS_HISTID(ATMOS_DIAGVARS_NUM) )

    reg_file_hist = .true.
    do v = 1, ATMOS_DIAGVARS_NUM
      call this%DIAGVARS_manager%Regist(        &
        ATMOS_DIAGVARS_VINFO(v), atm_mesh%mesh, & ! (in) 
        diag_vars(v), reg_file_hist             ) ! (out)
      
      this%DIAGVARS_HISTID(v) = diag_vars(v)%hist_id
    end do
    call this%DIAGVARS_manager%Final()

    !-- Setup information for input/output restart files. 

    is_specified = .true.
    !- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_VARS_RESTART,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("AtmosVars_Init",*) 'Not found namelist PARAM_ATMOS_VARS_RESTART. Default used.'
       is_specified = .false.
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("AtmosVars_Init",*) 'Not appropriate names in namelist PARAM_ATMOS_VARS_RESTART. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_VARS_RESTART)

    if (is_specified) then
      call this%restart_file%Init('ATMOS',                          &
        IN_BASENAME, IN_POSTFIX_TIMELABEL,                          &
        OUT_BASENAME, OUT_POSTFIX_TIMELABEL, OUT_DTYPE, OUT_TITLE,  &
        ATMOS_PROGVARS_NUM + ATMOS_AUXVARS_NUM,                     &
        mesh3D=atm_mesh%mesh                    )
    else
      call this%restart_file%Init('ATMOS',      &
        ATMOS_PROGVARS_NUM + ATMOS_AUXVARS_NUM, &
        mesh3D=atm_mesh%mesh                    )
    end if

    return
  end subroutine AtmosVars_Init

  subroutine AtmosVars_Final( this )
    implicit none
    class(AtmosVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosVars_Final',*)

    call this%restart_file%Final()

    call this%PROGVARS_comm%Final()
    call this%PROGVARS_manager%Final()

    call this%AUXVARS_manager%Final()

    deallocate( this%DIAGVARS_HISTID )

    return
  end subroutine AtmosVars_Final

  subroutine AtmosVars_history( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosVars), intent(in) :: this
  
    integer :: n
    integer :: v
    class(MeshBase3D), pointer :: mesh3D
    class(LocalMesh3D), pointer :: lcmesh
    integer :: hst_id

    type(MeshField3D) :: tmp_field
    !-------------------------------------------------------------------------

    mesh3D => this%PROG_VARS(1)%mesh
    do v = 1, ATMOS_PROGVARS_NUM
      hst_id = this%PROG_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%PROG_VARS(v) )
    end do

    do v = 1, ATMOS_AUXVARS_NUM
      hst_id = this%AUX_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%AUX_VARS(v) )
    end do

    call tmp_field%Init( "tmp_field", "", mesh3D)
    do v = 1, ATMOS_DIAGVARS_NUM
      hst_id = this%DIAGVARS_HISTID(v)
      if ( hst_id > 0 ) then
        call vars_calc_diagvar( this, ATMOS_DIAGVARS_VINFO(v)%NAME, tmp_field )
        call FILE_HISTORY_meshfield_put( hst_id, tmp_field )
      end if
    end do
    call tmp_field%Final()

    return
  end subroutine AtmosVars_history

  subroutine AtmosVar_Read_restart_file( this, atmos_mesh )

    use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
    use scale_meshfieldcomm_base, only: MeshFieldContainer
    implicit none
    
    class(AtmosVars), intent(inout), target :: this
    class(AtmosMesh), intent(in) :: atmos_mesh

    integer :: v

    type(MeshFieldCommCubeDom3D) :: auxvars_comm
    type(MeshFieldContainer) :: auxvars_comm_list(ATMOS_AUXVARS_NUM)
    !---------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOSVar_read_restart_file",*) 'Open restart file (ATMOS) '

    !- Open restart file
    call this%restart_file%Open()

    !- Read restart file
    do v=1, ATMOS_PROGVARS_NUM
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%PROG_VARS(v)%varname, &
        this%PROG_VARS(v)                                                      )
    end do
    do v=1, ATMOS_AUXVARS_NUM
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%AUX_VARS(v)%varname, &
        this%AUX_VARS(v)                                                      )
    end do

    !- Close restart file
    LOG_INFO("ATMOSVar_read_restart_file",*) 'Close restart file (ATMOS) '
    call this%restart_file%Close()

    !- Communicate halo data of hydrostatic variables

    call auxvars_comm%Init( ATMOS_AUXVARS_NUM, 0, atmos_mesh%mesh )
    do v=1, ATMOS_AUXVARS_NUM
      auxvars_comm_list(v)%field3d => this%AUX_VARS(v)
    end do
    call auxvars_comm%Put( auxvars_comm_list, 1 )
    call auxvars_comm%Exchange()
    call auxvars_comm%Get( auxvars_comm_list, 1 )
    call auxvars_comm%Final()

    return
  end subroutine AtmosVar_Read_restart_file

  subroutine AtmosVar_write_restart_file( this )

    implicit none
    class(AtmosVars), intent(inout) :: this

    integer :: v, rf_vid
    !---------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("ATMOSVar_write_restart_file",*) 'Create restart file (ATMOS) '

    !- Create restart file
    call this%restart_file%Create()

    !- Define variables
    do v=1, ATMOS_PROGVARS_NUM
      rf_vid = v
      call this%restart_file%Def_var( this%PROG_VARS(v),  &
        ATMOS_PROGVARS_VINFO(v)%DESC, rf_vid, DIMTYPE_XYZ )
    end do
    do v=1, ATMOS_AUXVARS_NUM
      rf_vid = ATMOS_PROGVARS_NUM + v
      call this%restart_file%Def_var( this%AUX_VARS(v),   &
        ATMOS_AUXVARS_VINFO(v)%DESC, rf_vid, DIMTYPE_XYZ  )
    end do
    call this%restart_file%End_def()

    !- Write restart file
    do v=1, ATMOS_PROGVARS_NUM
      rf_vid = v
      call this%restart_file%Write_var(rf_vid, this%PROG_VARS(v) )
    end do
    do v=1, ATMOS_AUXVARS_NUM
      rf_vid = ATMOS_PROGVARS_NUM + v
     call this%restart_file%Write_var(rf_vid, this%AUX_VARS(v) )
    end do

    !- Close restart file
    LOG_INFO("ATMOSVar_write_restart_file",*) 'Close restart file (ATMOS) '
    call this%restart_file%Close()

    return
  end subroutine AtmosVar_write_restart_file

  subroutine AtmosVars_GetLocalMeshField( domID, mesh, prgvars_list, auxvars_list, &
     varid,                                                                        &
     var, DENS_hyd, PRES_hyd, lcmesh3D                                             &
    )
    !------------
    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    use scale_localmesh_base, only: LocalMeshBase
    use scale_localmesh_3d, only: LocalMesh3D

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    integer, intent(in) :: varid
    class(LocalMeshFieldBase), pointer, intent(out) :: var
    class(LocalMeshFieldBase), pointer, intent(out), optional :: DENS_hyd, PRES_hyd
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call prgvars_list%Get(varid, field)
    call field%GetLocalMeshField(domID, var)

    if (present(DENS_hyd)) then
      call auxvars_list%Get(ATMOS_AUXVARS_DENSHYDRO_ID, field)
      call field%GetLocalMeshField(domID, DENS_hyd)
    end if
    if (present(PRES_hyd)) then
      call auxvars_list%Get(ATMOS_AUXVARS_PRESHYDRO_ID, field)
      call field%GetLocalMeshField(domID, PRES_hyd)
    end if

    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshField

  subroutine AtmosVars_GetLocalMeshFields( domID, mesh, prgvars_list, auxvars_list, &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
    DENS_hyd, PRES_hyd, lcmesh3D                                    &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    use scale_localmesh_base, only: LocalMeshBase
    use scale_localmesh_3d, only: LocalMesh3D

    implicit none
    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer, intent(out) :: DENS_hyd, PRES_hyd
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call prgvars_list%Get(ATMOS_PROGVARS_DDENS_ID, field)
    call field%GetLocalMeshField(domID, DDENS)

    call prgvars_list%Get(ATMOS_PROGVARS_MOMX_ID, field)
    call field%GetLocalMeshField(domID, MOMX)
    
    call prgvars_list%Get(ATMOS_PROGVARS_MOMY_ID, field)
    call field%GetLocalMeshField(domID, MOMY)

    call prgvars_list%Get(ATMOS_PROGVARS_MOMZ_ID, field)
    call field%GetLocalMeshField(domID, MOMZ)

    call prgvars_list%Get(ATMOS_PROGVARS_DRHOT_ID, field)
    call field%GetLocalMeshField(domID, DRHOT)
  
    !--
    call auxvars_list%Get(ATMOS_AUXVARS_DENSHYDRO_ID, field)
    call field%GetLocalMeshField(domID, DENS_hyd)

    call auxvars_list%Get(ATMOS_AUXVARS_PRESHYDRO_ID, field)
    call field%GetLocalMeshField(domID, PRES_hyd)

    !---
    
    if (present(lcmesh3D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh3D )

      select type(lcmesh)
      type is (LocalMesh3D)
        if (present(lcmesh3D)) lcmesh3D => lcmesh
      end select
    end if

    return
  end subroutine AtmosVars_GetLocalMeshFields

  subroutine vars_calc_diagvar( this, field_name, field_work ) 
    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CPdry => CONST_CPdry,    &
      CVdry => CONST_CVdry,    &
      PRES00 => CONST_PRE00

    implicit none
    class(AtmosVars), intent(in) :: this
    character(*), intent(in) :: field_name
    type(MeshField3D), intent(inout) :: field_work

    type(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    !--------------------------------------------------

    do n=1, field_work%mesh%LOCAL_MESH_NUM
      lcmesh3D => field_work%mesh%lcmesh_list(n)
      call vars_calc_diagnoseVar( field_name, field_work%local(n)%val,  &
        this%PROG_VARS(ATMOS_PROGVARS_DDENS_ID)%local(n)%val,           &
        this%PROG_VARS(ATMOS_PROGVARS_MOMX_ID)%local(n)%val,            &
        this%PROG_VARS(ATMOS_PROGVARS_MOMY_ID)%local(n)%val,            &
        this%PROG_VARS(ATMOS_PROGVARS_MOMZ_ID)%local(n)%val,            &
        this%PROG_VARS(ATMOS_PROGVARS_DRHOT_ID)%local(n)%val,           &
        this%AUX_VARS(ATMOS_AUXVARS_DENSHYDRO_ID)%local(n)%val,         & 
        this%AUX_VARS(ATMOS_AUXVARS_PRESHYDRO_ID)%local(n)%val,         &
        lcmesh3D, lcmesh3D%refElem3D )
    end do

    return
  end subroutine vars_calc_diagvar

  subroutine vars_calc_diagnoseVar( field_name, var_out,        &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,    &
    lcmesh, elem )

    use scale_const, only: &
      GRAV => CONST_GRAV,  &
      Rdry => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    implicit none

    type(LocalMesh3D), intent(in) :: lcmesh
    type(ElementBase3D), intent(in) :: elem
    character(*), intent(in) :: field_name
    real(RP), intent(out) :: var_out(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)

    integer :: ke
    real(RP) :: RHOT(elem%Np), DENS(elem%Np), PRES(elem%Np), THETA(elem%Np)

    !-------------------------------------------------------------------------

    select case(trim(field_name))
    case('U')
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = MOMX_(:,ke)/DENS(:)
      end do
    case('V')
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = MOMY_(:,ke)/DENS(:)
      end do        
    case('W')
      do ke=1, lcmesh%Ne
        DENS(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
        var_out(:,ke) = MOMZ_(:,ke)/DENS(:)
      end do
    case('PRES')        
      do ke=1, lcmesh%Ne
        RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
        var_out(:,ke) = PRES00 * (Rdry*RHOT(:)/PRES00)**(CPdry/Cvdry)
      end do
    case('T')        
      do ke=1, lcmesh%Ne
        RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
        PRES(:) = PRES00 * (Rdry*RHOT(:)/PRES00)**(CPdry/Cvdry)
        var_out(:,ke) = PRES(:) / (Rdry*(DDENS_(:,ke) + DENS_hyd(:,ke)))
      end do
    case('PT')        
      do ke=1, lcmesh%Ne
        RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
        var_out(:,ke) = RHOT(:) / (DDENS_(:,ke) + DENS_hyd(:,ke))
      end do 
    end select

    return
  end subroutine vars_calc_diagnoseVar  
end module mod_atmos_vars