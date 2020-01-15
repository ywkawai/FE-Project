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
    type(MeshFieldCommCubeDom3D) :: AUXVARS_comm
    
    type(ModelVarManager) :: DIAGVARS_manager    
    integer, allocatable :: DIAGVARS_HISTID(:)
  contains
    procedure :: Init => AtmosVars_Init
    procedure :: Final => AtmosVars_Final
    procedure :: History => AtmosVars_History
  end type AtmosVars

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

  integer, public, parameter :: ATMOS_AUXVARS_NUM          = 14
  integer, public, parameter :: ATMOS_AUXVARS_PRESHYDRO_ID = 1
  integer, public, parameter :: ATMOS_AUXVARS_DENSHYDRO_ID = 2
  integer, public, parameter :: ATMOS_AUXVARS_DxU_ID    = 3
  integer, public, parameter :: ATMOS_AUXVARS_DyU_ID    = 4
  integer, public, parameter :: ATMOS_AUXVARS_DzU_ID    = 5
  integer, public, parameter :: ATMOS_AUXVARS_DxV_ID    = 6
  integer, public, parameter :: ATMOS_AUXVARS_DyV_ID    = 7
  integer, public, parameter :: ATMOS_AUXVARS_DzV_ID    = 8
  integer, public, parameter :: ATMOS_AUXVARS_DxW_ID    = 9
  integer, public, parameter :: ATMOS_AUXVARS_DyW_ID    = 10
  integer, public, parameter :: ATMOS_AUXVARS_DzW_ID    = 11
  integer, public, parameter :: ATMOS_AUXVARS_DxPT_ID   = 12
  integer, public, parameter :: ATMOS_AUXVARS_DyPT_ID   = 13
  integer, public, parameter :: ATMOS_AUXVARS_DzPT_ID   = 14  

  type(VariableInfo), public :: ATMOS_AUXVARS_VINFO(ATMOS_AUXVARS_NUM)
  DATA ATMOS_AUXVARS_VINFO / &
    VariableInfo( ATMOS_AUXVARS_PRESHYDRO_ID, 'PRES_hyd', 'hydrostatic part of pressure',  &
                  'Pa',  3, 'XYZ',  ''                                                  ), &
    VariableInfo( ATMOS_AUXVARS_DENSHYDRO_ID , 'DENS_hyd', 'hydrostatic part of density',  &
                  'kg/m3', 3, 'XYZ', ''                                                 ), &    
    VariableInfo( ATMOS_AUXVARS_DxU_ID, 'DxU', '', '', 3, 'XYZ', ''          ), &
    VariableInfo( ATMOS_AUXVARS_DyU_ID, 'DyU', '', '', 3, 'XYZ', ''          ), &
    VariableInfo( ATMOS_AUXVARS_DzU_ID, 'DzU', '', '', 3, 'XYZ', ''          ), &
    VariableInfo( ATMOS_AUXVARS_DxV_ID, 'DxV', '', '', 3, 'XYZ', ''          ), &
    VariableInfo( ATMOS_AUXVARS_DyV_ID, 'DyV', '', '', 3, 'XYZ', ''          ), &
    VariableInfo( ATMOS_AUXVARS_DzV_ID, 'DzV', '', '', 3, 'XYZ', ''          ), &
    VariableInfo( ATMOS_AUXVARS_DxW_ID, 'DxW', '', '', 3, 'XYZ', ''          ), &
    VariableInfo( ATMOS_AUXVARS_DyW_ID, 'DyW', '', '', 3, 'XYZ', ''          ), &
    VariableInfo( ATMOS_AUXVARS_DzW_ID, 'DzW', '', '', 3, 'XYZ', ''          ), &
    VariableInfo( ATMOS_AUXVARS_DxPT_ID, 'DxPT', '', '', 3, 'XYZ', ''        ), &
    VariableInfo( ATMOS_AUXVARS_DyPT_ID, 'DyPT', '', '', 3, 'XYZ', ''        ), &
    VariableInfo( ATMOS_AUXVARS_DzPT_ID, 'DzPT', '', '', 3, 'XYZ', ''        )  /
  
  integer, public, parameter :: ATMOS_DIAGVARS_NUM       = 5
  integer, public, parameter :: ATMOS_DIAGVARS_U_ID      = 1
  integer, public, parameter :: ATMOS_DIAGVARS_V_ID      = 2
  integer, public, parameter :: ATMOS_DIAGVARS_W_ID      = 3
  integer, public, parameter :: ATMOS_DIAGVARS_PRES_ID   = 4
  integer, public, parameter :: ATMOS_DIAGVARS_T_ID      = 5  

  type(VariableInfo), public :: ATMOS_DIAGVARS_VINFO(ATMOS_DIAGVARS_NUM)
  DATA ATMOS_DIAGVARS_VINFO / &
    VariableInfo( ATMOS_DIAGVARS_U_ID, 'U', 'velocity u', 'm/s', 3, 'XYZ', 'x_wind'             ), &
    VariableInfo( ATMOS_DIAGVARS_V_ID, 'V', 'velocity v', 'm/s', 3, 'XYZ', 'y_wind'             ), &  
    VariableInfo( ATMOS_DIAGVARS_W_ID, 'W', 'velocity w', 'm/s', 3, 'XYZ', 'upward_air_velocity'), &
    VariableInfo( ATMOS_DIAGVARS_PRES_ID, 'PRES', 'pressure', 'Pa', 3, 'XYZ', 'air_pressure'    ), &
    VariableInfo( ATMOS_DIAGVARS_T_ID   , 'T', 'temperature', 'K', 3, 'XYZ', 'air_temperature'  )  /

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
        this%PROG_VARS(v)%local(n)%val(:,:) = 0.E20_RP
      end do         
    end do

    call this%PROGVARS_comm%Init(ATMOS_PROGVARS_NUM, 0, atm_mesh%mesh)
    call this%PROGVARS_manager%MeshFieldComm_Prepair( this%PROGVARS_comm, this%PROG_VARS(:) )

    !- Initialize auxiliary variables

    call this%AUXVARS_manager%Init()
    allocate( this%AUX_VARS(ATMOS_AUXVARS_NUM) )
    
    reg_file_hist = .true.
    do v = ATMOS_AUXVARS_PRESHYDRO_ID, ATMOS_AUXVARS_DxU_ID-1
      call this%AUXVARS_manager%Regist(        &
        ATMOS_AUXVARS_VINFO(v), atm_mesh%mesh, & ! (in) 
        this%AUX_VARS(v), reg_file_hist        ) ! (out)
      do n = 1, atm_mesh%mesh%LOCAL_MESH_NUM
        this%AUX_VARS(v)%local(n)%val(:,:) = 1.0_RP
      end do             
    end do

    reg_file_hist = .false.    
    do v = ATMOS_AUXVARS_DxU_ID, ATMOS_AUXVARS_DzPT_ID
      call this%AUXVARS_manager%Regist(        &
        ATMOS_AUXVARS_VINFO(v), atm_mesh%mesh, & ! (in) 
        this%AUX_VARS(v), reg_file_hist        ) ! (out)      

      do n = 1, atm_mesh%mesh%LOCAL_MESH_NUM
        this%AUX_VARS(v)%local(n)%val(:,:) = 1.0_RP
      end do
    end do

    call this%AUXVARS_comm%Init(3*4, 0, atm_mesh%mesh)
    call this%AUXVARS_manager%MeshFieldComm_Prepair( &
      this%AUXVARS_comm, this%AUX_VARS(ATMOS_AUXVARS_DxU_ID:ATMOS_AUXVARS_DzPT_ID) )

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

    return
  end subroutine AtmosVars_Init

  subroutine AtmosVars_Final( this )
    implicit none
    class(AtmosVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosVars_Final',*)

    call this%PROGVARS_comm%Final()
    call this%PROGVARS_manager%Final()

    call this%AUXVARS_comm%Final()
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

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
  !     call vars_calc_diagnoseVar( &
  !       U%local(n)%val, V%local(n)%val, W%local(n)%val, DPRES%local(n)%val, TEMP%local(n)%val, DTHETA%local(n)%val,   &
  !       DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,              &
  !       PRES_hydro%local(n)%val, DENS_hydro%local(n)%val,                                                             &
  !       lcmesh, lcmesh%refElem3D )
    end do

    do v = 1, ATMOS_AUXVARS_NUM
      hst_id = this%AUX_VARS(v)%hist_id
      if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%AUX_VARS(v) )
    end do

    call tmp_field%Init( "tmp_field", "", mesh3D)
    do v = 1, ATMOS_DIAGVARS_NUM
      hst_id = this%DIAGVARS_HISTID(v)
      if (hst_id > 0) then
        call vars_calc_diagvar( this, ATMOS_DIAGVARS_VINFO(v)%NAME, tmp_field )
        call FILE_HISTORY_meshfield_put( hst_id, tmp_field )
      end if
    end do
    call tmp_field%Final()

    return
  end subroutine AtmosVars_history

  subroutine AtmosVars_GetLocalMeshFields( domID, mesh, prgvars_list, auxvars_list, &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
    GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW, GxPT, GyPT, GzPT,  &
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
    class(LocalMeshFieldBase), pointer, intent(out) :: &
      GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW, GxPT, GyPT, GzPT
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
    call auxvars_list%Get(ATMOS_AUXVARS_DxU_ID, field)
    call field%GetLocalMeshField(domID, GxU)
    call auxvars_list%Get(ATMOS_AUXVARS_DyU_ID, field)
    call field%GetLocalMeshField(domID, GyU)
    call auxvars_list%Get(ATMOS_AUXVARS_DzU_ID, field)
    call field%GetLocalMeshField(domID, GzU)

    call auxvars_list%Get(ATMOS_AUXVARS_DxV_ID, field)
    call field%GetLocalMeshField(domID, GxV)
    call auxvars_list%Get(ATMOS_AUXVARS_DyV_ID, field)
    call field%GetLocalMeshField(domID, GyV)
    call auxvars_list%Get(ATMOS_AUXVARS_DzV_ID, field)
    call field%GetLocalMeshField(domID, GzV)
    
    call auxvars_list%Get(ATMOS_AUXVARS_DxW_ID, field)
    call field%GetLocalMeshField(domID, GxW)
    call auxvars_list%Get(ATMOS_AUXVARS_DyW_ID, field)
    call field%GetLocalMeshField(domID, GyW)
    call auxvars_list%Get(ATMOS_AUXVARS_DzU_ID, field)
    call field%GetLocalMeshField(domID, GzW)
    
    call auxvars_list%Get(ATMOS_AUXVARS_DxPT_ID, field)
    call field%GetLocalMeshField(domID, GxPT)
    call auxvars_list%Get(ATMOS_AUXVARS_DyPT_ID, field)
    call field%GetLocalMeshField(domID, GyPT)
    call auxvars_list%Get(ATMOS_AUXVARS_DzPT_ID, field)
    call field%GetLocalMeshField(domID, GzPT)    
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
    case('THETA')        
      do ke=1, lcmesh%Ne
        RHOT(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
        var_out(:,ke) = RHOT(:) / (DDENS_(:,ke) + DENS_hyd(:,ke))
      end do 
    end select

    return
  end subroutine vars_calc_diagnoseVar  
end module mod_atmos_vars