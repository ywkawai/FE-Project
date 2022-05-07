!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_dyn_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_element_base, only: ElementBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D  
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfield_base, only: MeshField3D, MeshField2D
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use scale_meshfieldcomm_base, only: MeshFieldContainer
  
  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo
  use scale_model_mesh_manager, only: ModelMeshBase

  use mod_atmos_mesh, only: AtmosMesh
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: AtmosDynVars
    type(MeshField2D), allocatable :: AUX_VARS2D(:)
    type(ModelVarManager) :: AUXVARS2D_manager

    type(MeshField3D), allocatable :: MASS_FLUX_VARS3D(:)
    type(ModelVarManager) :: MASS_FLUX_manager
    integer :: MASS_FLUX_commid

    type(MeshField3D), allocatable :: AUXTRC_FLUX_VARS3D(:)
    type(ModelVarManager) :: AUXTRC_FLUX_VAR3D_manager
    integer :: AUXTRC_FLUX_VAR3D_commid

    type(MeshField3D), allocatable :: TRCVARS3D(:)
    type(ModelVarManager) :: TRCVAR3D_manager
    integer :: TRCVAR3D_commid

    type(MeshField3D), allocatable :: AUX_TRCVARS3D(:)
    type(ModelVarManager) :: AUXTRCVAR3D_manager
    integer :: AUXTRCVAR3D_commid

    type(MeshField3D), allocatable :: NUMDIFF_FLUX_VARS3D(:)
    type(ModelVarManager) :: NUMDIFF_FLUX_manager
    integer :: NUMDIFF_FLUX_commid

    type(MeshField3D), allocatable :: NUMDIFF_TEND_VARS3D(:)
    type(ModelVarManager) :: NUMDIFF_TEND_manager
    integer :: NUMDIFF_TEND_commid

    type(MeshField3D), allocatable :: ANALYSIS_VARS3D(:)
    type(ModelVarManager) :: ANALYSISVARS_manager

  contains
    procedure :: Init => AtmosDynVars_Init
    procedure :: Final => AtmosDynVars_Final
    procedure :: History => AtmosDynVars_History
  end type AtmosDynVars

  public :: AtmosDynAuxVars_GetLocalMeshFields
  public :: AtmosDynMassFlux_GetLocalMeshFields
  public :: AtmosDynNumDiffFlux_GetLocalMeshFields
  public :: AtmosDynNumDiffTend_GetLocalMeshFields
  !public :: AtmosDynVars_GetLocalMeshFields_analysis

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-
  integer, public, parameter :: ATMOS_DYN_AUXVARS2D_NUM          = 1
  integer, public, parameter :: ATMOS_DYN_AUXVARS2D_CORIOLIS_ID  = 1

  type(VariableInfo), public :: ATMOS_DYN_AUXVARS2D_VINFO(ATMOS_DYN_AUXVARS2D_NUM)
  DATA ATMOS_DYN_AUXVARS2D_VINFO / &
    VariableInfo( ATMOS_DYN_AUXVARS2D_CORIOLIS_ID, 'CORIOLIS', 'coriolis parameter',  &
                  's-1',  2, 'XY',  ''                                             )  / 

  !-
  integer, public, parameter :: ATMOS_DYN_TRCVARS3D_NUM         = 2
  integer, public, parameter :: ATMOS_DYN_TRCVARS3D_TRCADV_ID   = 1
  integer, public, parameter :: ATMOS_DYN_TRCVARS3D_DENS_ID     = 2

  type(VariableInfo), public :: ATMOS_DYN_TRCVARS3D_VINFO(ATMOS_DYN_TRCVARS3D_NUM )
  DATA ATMOS_DYN_TRCVARS3D_VINFO / &
    VariableInfo( ATMOS_DYN_TRCVARS3D_TRCADV_ID, 'TRCADV', '',     '1',  3, 'XYZ',  '' ), &
    VariableInfo( ATMOS_DYN_TRCVARS3D_DENS_ID  ,   'DENS', '', 'kg/m3',  3, 'XYZ',  '' )  / 
                
  !-
  integer, public, parameter :: ATMOS_DYN_AUXTRCVARS3D_NUM         = 1
  integer, public, parameter :: ATMOS_DYN_AUXTRCVARS3D_FCTCOEF_ID  = 1

  type(VariableInfo), public :: ATMOS_DYN_AUXTRCVARS3D_VINFO(ATMOS_DYN_AUXTRCVARS3D_NUM)
  DATA ATMOS_DYN_AUXTRCVARS3D_VINFO / &
    VariableInfo( ATMOS_DYN_AUXTRCVARS3D_FCTCOEF_ID, 'TRCADV_FCTCOEF', '',  &
                  '1',  3, 'XYZ',  ''                                    )  / 
  
  !-
  integer, public, parameter :: ATMOS_DYN_MASS_FLUX_NUM   = 4
  integer, public, parameter :: ATMOS_DYN_ALPHDENS_ID     = 1  
  integer, public, parameter :: ATMOS_DYN_MASSFLX_Z_ID    = 2  
  integer, public, parameter :: ATMOS_DYN_MASSFLX_X_ID    = 3
  integer, public, parameter :: ATMOS_DYN_MASSFLX_Y_ID    = 4

  type(VariableInfo), public :: ATMOS_DYN_MASS_FLUX_VINFO(ATMOS_DYN_MASS_FLUX_NUM)
  DATA ATMOS_DYN_MASS_FLUX_VINFO / &
    VariableInfo( ATMOS_DYN_ALPHDENS_ID, 'ALPHDENS', 'alphaXdens',             &
                    'kg/m3',  3, 'XYZ',  ''                                    ),   &
    VariableInfo( ATMOS_DYN_MASSFLX_Z_ID, 'MASSFLX_Z', 'flux in z-direction',  &
                  'kg/s/m2',  3, 'XYZ',  ''                                    ),   &
    VariableInfo( ATMOS_DYN_MASSFLX_X_ID, 'MASSFLX_X', 'flux in x-direction',  &
                  'kg/s/m2',  3, 'XYZ',  ''                                    ),   & 
    VariableInfo( ATMOS_DYN_MASSFLX_Y_ID, 'MASSFLX_Y', 'flux in y-direction',  &
                  'kg/s/m2',  3, 'XYZ',  ''                                    )    / 

  !-
  integer, public, parameter :: ATMOS_DYN_NUMDIFF_FLUX_NUM   = 3
  integer, public, parameter :: ATMOS_DYN_NUMDIFFFLX_X_ID    = 1
  integer, public, parameter :: ATMOS_DYN_NUMDIFFFLX_Y_ID    = 2
  integer, public, parameter :: ATMOS_DYN_NUMDIFFFLX_Z_ID    = 3

  type(VariableInfo), public :: ATMOS_DYN_NUMDIFF_FLUX_VINFO(ATMOS_DYN_NUMDIFF_FLUX_NUM)
  DATA ATMOS_DYN_NUMDIFF_FLUX_VINFO / &
    VariableInfo( ATMOS_DYN_NUMDIFFFLX_X_ID, 'DIFFFLX_X', 'flux in x-direction',  &
                  '?.m/s',  3, 'XYZ',  ''                                           ),   & 
    VariableInfo( ATMOS_DYN_NUMDIFFFLX_Y_ID, 'DIFFFLX_Y', 'flux in y-direction',  &
                  '?.m/s',  3, 'XYZ',  ''                                           ),   &
    VariableInfo( ATMOS_DYN_NUMDIFFFLX_Z_ID, 'DIFFFLX_Z', 'flux in z-direction',  &
                  '?.m/s',  3, 'XYZ',  ''                                           )    / 
       
  !-
  integer, public, parameter :: ATMOS_DYN_NUMDIFF_TEND_NUM   = 2
  integer, public, parameter :: ATMOS_DYN_NUMDIFF_LAPLAH_ID  = 1
  integer, public, parameter :: ATMOS_DYN_NUMDIFF_LAPLAV_ID  = 2

  type(VariableInfo), public :: ATMOS_DYN_NUMDIFF_TEND_VINFO(ATMOS_DYN_NUMDIFF_TEND_NUM)
  DATA ATMOS_DYN_NUMDIFF_TEND_VINFO / &
    VariableInfo( ATMOS_DYN_NUMDIFF_LAPLAH_ID, 'NUMDIFF_LAPLAH', 'tendency due to nundiff',  &
                  '?/s',  3, 'XYZ',  ''                                                   ), &
    VariableInfo( ATMOS_DYN_NUMDIFF_LAPLAV_ID, 'NUMDIFF_LAPLAV', 'tendency due to nundiff',  &
                  '?/s',  3, 'XYZ',  ''                                                   )  /

  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVARS_NUM          = 6
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t        = 1
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t_advX   = 2
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t_advY   = 3
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t_advZ   = 4
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t_lift   = 5
  ! integer, public, parameter :: ATMOS_DYN_ANALYSISVAR_MOMZ_t_buoy   = 6

  ! type(VariableInfo), public :: ATMOS_DYN_ANALYSISVAR_VINFO(ATMOS_DYN_ANALYSISVARS_NUM)
  ! DATA ATMOS_DYN_ANALYSISVAR_VINFO / &
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t, 'MOMZ_t', 'MOMZ_t',    &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ), & 
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t_advX, 'MOMZ_t_advx', 'MOMZ_t_advX',  &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ), & 
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t_advY, 'MOMZ_t_advy', 'MOMZ_t_advY',  &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ), & 
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t_advZ, 'MOMZ_t_advz', 'MOMZ_t_advZ',  &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ), & 
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t_lift, 'MOMZ_t_lift', 'MOMZ_t_lift',  &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ), & 
  !   VariableInfo( ATMOS_DYN_ANALYSISVAR_MOMZ_t_buoy, 'MOMZ_t_buoy', 'MOMZ_t_buo',    &
  !                 'kg.m-2.s-1',  3, 'XYZ',  ''                      ) & 
  ! / 
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains
!OCL SERIAL
  subroutine AtmosDynVars_Init( this, model_mesh )
    implicit none
    class(AtmosDynVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    integer :: v
    integer :: n
    logical :: reg_file_hist

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D

    !--------------------------------------------------

    LOG_INFO('AtmosDynVars_Init',*)

    nullify( atm_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select
    mesh3D => atm_mesh%ptr_mesh
    
    call mesh3D%GetMesh2D( mesh2D )

    !- Initialize 2D auxiliary variables

    call this%AUXVARS2D_manager%Init()
    allocate( this%AUX_VARS2D(ATMOS_DYN_AUXVARS2D_NUM) )

    reg_file_hist = .false.    
    do v = 1, ATMOS_DYN_AUXVARS2D_NUM
      call this%AUXVARS2D_manager%Regist(                   &
        ATMOS_DYN_AUXVARS2D_VINFO(v), mesh2D,               & ! (in) 
        this%AUX_VARS2D(v), reg_file_hist                   ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%AUX_VARS2D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    !- Initialize 3D auxiliary variables for tracer advection

    call this%TRCVAR3D_manager%Init()
    allocate( this%TRCVARS3D(ATMOS_DYN_TRCVARS3D_NUM) )
    reg_file_hist = .false.    
    do v = 1, ATMOS_DYN_TRCVARS3D_NUM
      call this%TRCVAR3D_manager%Regist(               &
        ATMOS_DYN_TRCVARS3D_VINFO(v), mesh3D,          & ! (in) 
        this%TRCVARS3D(v), reg_file_hist               ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%TRCVARS3D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    v = ATMOS_DYN_TRCVARS3D_TRCADV_ID
    call atm_mesh%Create_communicator( &
      1, 0,                            & ! (in) 
      this%TRCVAR3D_manager,           & ! (in)
      this%TRCVARS3D(v:v),             & ! (in)
      this%TRCVAR3D_commid             ) ! (out)

    !- Initialize 3D auxiliary variables for preserving nonnegativity in tracer advection

    call this%AUXTRCVAR3D_manager%Init()
    allocate( this%AUX_TRCVARS3D(ATMOS_DYN_AUXTRCVARS3D_NUM) )
    reg_file_hist = .false.    
    do v = 1, ATMOS_DYN_AUXTRCVARS3D_NUM
      call this%AUXTRCVAR3D_manager%Regist(                &
        ATMOS_DYN_AUXTRCVARS3D_VINFO(v), mesh3D,           & ! (in) 
        this%AUX_TRCVARS3D(v), reg_file_hist               ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%AUX_TRCVARS3D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    call atm_mesh%Create_communicator( &
      1, 0,                            & ! (in) 
      this%AUXTRCVAR3D_manager,        & ! (in)
      this%AUX_TRCVARS3D(:),           & ! (in)
      this%AUXTRCVAR3D_commid          ) ! (out)
        
    !- Initialize variables to store time-averaged 3D mass flux

    call this%MASS_FLUX_manager%Init()
    allocate( this%MASS_FLUX_VARS3D(ATMOS_DYN_MASS_FLUX_NUM) )
    
    reg_file_hist = .false.    
    do v = 1, ATMOS_DYN_MASS_FLUX_NUM
      call this%MASS_FLUX_manager%Regist(                   &
        ATMOS_DYN_MASS_FLUX_VINFO(v), mesh3D,               & ! (in) 
        this%MASS_FLUX_VARS3D(v), reg_file_hist             ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%MASS_FLUX_VARS3D(v)%local(n)%val(:,:) = 0.0_RP
      end do      
    end do

    call atm_mesh%Create_communicator( &
      2, 1,                            & ! (in) 
      this%MASS_FLUX_manager,          & ! (in)
      this%MASS_FLUX_VARS3D(:),        & ! (in)
      this%MASS_FLUX_commid            ) ! (out)

    !- Initialize variables to store diffusive fluxes with explicit numerical diffusion

    call this%NUMDIFF_FLUX_manager%Init()
    allocate( this%NUMDIFF_FLUX_VARS3D(ATMOS_DYN_NUMDIFF_FLUX_NUM) )

    reg_file_hist = .false.    
    do v = 1, ATMOS_DYN_NUMDIFF_FLUX_NUM
      call this%NUMDIFF_FLUX_manager%Regist(                &
        ATMOS_DYN_NUMDIFF_FLUX_VINFO(v), mesh3D,            & ! (in) 
        this%NUMDIFF_FLUX_VARS3D(v), reg_file_hist          ) ! (out)
      
      do n = 1, mesh3D%LOCAL_MESH_NUM
        this%NUMDIFF_FLUX_VARS3D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    call atm_mesh%Create_communicator( &
      ATMOS_DYN_NUMDIFF_FLUX_NUM, 0,   & ! (in)
      this%NUMDIFF_FLUX_manager,       & ! (inout)
      this%NUMDIFF_FLUX_VARS3D(:),     & ! (in)
      this%NUMDIFF_FLUX_commid         ) ! (out)

    !-
    call this%NUMDIFF_TEND_manager%Init()
    allocate( this%NUMDIFF_TEND_VARS3D(ATMOS_DYN_NUMDIFF_TEND_NUM) )

    reg_file_hist = .false.    
    do v = 1, ATMOS_DYN_NUMDIFF_TEND_NUM
      call this%NUMDIFF_TEND_manager%Regist(                &
        ATMOS_DYN_NUMDIFF_TEND_VINFO(v), mesh3D,            & ! (in) 
        this%NUMDIFF_TEND_VARS3D(v), reg_file_hist          ) ! (out)
      
      do n = 1,  mesh3D%LOCAL_MESH_NUM
        this%NUMDIFF_TEND_VARS3D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    call atm_mesh%Create_communicator( &
      ATMOS_DYN_NUMDIFF_TEND_NUM, 0,   & ! (in)
      this%NUMDIFF_TEND_manager,       & ! (inout)
      this%NUMDIFF_TEND_VARS3D(:),     & ! (in)
      this%NUMDIFF_TEND_commid         ) ! (out)

    !-
    ! call this%ANALYSISVARS_manager%Init()
    ! allocate( this%ANALYSIS_VARS3D(ATMOS_DYN_ANALYSISVARS_NUM) )

    ! reg_file_hist = .true.
    ! do v = 1, ATMOS_DYN_ANALYSISVARS_NUM
    !   call this%ANALYSISVARS_manager%Regist( &
    !     ATMOS_DYN_ANALYSISVAR_VINFO(v), atm_mesh%mesh, &
    !     this%ANALYSIS_VARS3D(v), reg_file_hist )
    ! end do

    return
  end subroutine AtmosDynVars_Init

  subroutine AtmosDynVars_Final( this )
    implicit none
    class(AtmosDynVars), intent(inout) :: this

    !--------------------------------------------------

    LOG_INFO('AtmosDynVars_Final',*)

    call this%AUXVARS2D_manager%Final()
    deallocate( this%AUX_VARS2D )

    call this%TRCVAR3D_manager%Final()
    deallocate( this%TRCVARS3D )

    call this%AUXTRCVAR3D_manager%Final()
    deallocate( this%AUX_TRCVARS3D )

    call this%MASS_FLUX_manager%Final()
    deallocate( this%MASS_FLUX_VARS3D )

    call this%NUMDIFF_FLUX_manager%Final()
    deallocate( this%NUMDIFF_FLUX_VARS3D )

    call this%NUMDIFF_TEND_manager%Final()
    deallocate( this%NUMDIFF_TEND_VARS3D )
    
    !call this%ANALYSISVARS_manager%Final()

    return
  end subroutine AtmosDynVars_Final

  subroutine AtmosDynVars_History( this )
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosDynVars), intent(in) :: this

    ! integer :: v
    ! integer :: hst_id
    !-------------------------------------------------------------------------

    ! do v = 1, ATMOS_DYN_ANALYSISVARS_NUM
    !   hst_id = this%ANALYSIS_VARS3D(v)%hist_id
    !   if ( hst_id > 0 ) call FILE_HISTORY_meshfield_put( hst_id, this%ANALYSIS_VARS3D(v) )
    ! end do
    
    return
  end subroutine AtmosDynVars_History

  subroutine AtmosDynAuxVars_GetLocalMeshFields( domID, mesh, auxvars_list, &
    Coriolis,                                                               &
    lcmesh3D                                                                &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: Coriolis
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call auxvars_list%Get(ATMOS_DYN_AUXVARS2D_CORIOLIS_ID, field)
    call field%GetLocalMeshField(domID, Coriolis)
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
  end subroutine AtmosDynAuxVars_GetLocalMeshFields

  subroutine AtmosDynMassFlux_GetLocalMeshFields( domID, mesh, massflx_list, &
    ALPH_DENS, MASSFLX_X, MASSFLX_Y, MASSFLX_Z,                              &
    lcmesh3D                                                                 &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: massflx_list
    class(LocalMeshFieldBase), pointer, intent(out) :: ALPH_DENS
    class(LocalMeshFieldBase), pointer, intent(out) :: MASSFLX_X
    class(LocalMeshFieldBase), pointer, intent(out) :: MASSFLX_Y
    class(LocalMeshFieldBase), pointer, intent(out) :: MASSFLX_Z        
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call massflx_list%Get(ATMOS_DYN_ALPHDENS_ID, field)
    call field%GetLocalMeshField(domID, ALPH_DENS)

    call massflx_list%Get(ATMOS_DYN_MASSFLX_X_ID, field)
    call field%GetLocalMeshField(domID, MASSFLX_X)

    call massflx_list%Get(ATMOS_DYN_MASSFLX_Y_ID, field)
    call field%GetLocalMeshField(domID, MASSFLX_Y)    

    call massflx_list%Get(ATMOS_DYN_MASSFLX_Z_ID, field)
    call field%GetLocalMeshField(domID, MASSFLX_Z)    
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
  end subroutine AtmosDynMassFlux_GetLocalMeshFields  

  subroutine AtmosDynNumDiffFlux_GetLocalMeshFields( domID, mesh, auxvars_list, &
    NUMDIFF_FLUX_X, NUMDIFF_FLUX_Y, NUMDIFF_FLUX_Z,                             &
    lcmesh3D                                                                    &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: NUMDIFF_FLUX_X
    class(LocalMeshFieldBase), pointer, intent(out) :: NUMDIFF_FLUX_Y
    class(LocalMeshFieldBase), pointer, intent(out) :: NUMDIFF_FLUX_Z
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call auxvars_list%Get(ATMOS_DYN_NUMDIFFFLX_X_ID, field)
    call field%GetLocalMeshField(domID, NUMDIFF_FLUX_X)

    call auxvars_list%Get(ATMOS_DYN_NUMDIFFFLX_Y_ID, field)
    call field%GetLocalMeshField(domID, NUMDIFF_FLUX_Y)

    call auxvars_list%Get(ATMOS_DYN_NUMDIFFFLX_Z_ID, field)
    call field%GetLocalMeshField(domID, NUMDIFF_FLUX_Z)
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
  end subroutine AtmosDynNumDiffFlux_GetLocalMeshFields

  subroutine AtmosDynNumDiffTend_GetLocalMeshFields( domID, mesh, auxvars_list, &
    NUMDIFF_LAPLAH, NUMDIFF_LAPLAV,                                             &
    lcmesh3D                                                                    &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: NUMDIFF_LAPLAH
    class(LocalMeshFieldBase), pointer, intent(out) :: NUMDIFF_LAPLAV
    class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call auxvars_list%Get(ATMOS_DYN_NUMDIFF_LAPLAH_ID, field)
    call field%GetLocalMeshField(domID, NUMDIFF_LAPLAH)

    call auxvars_list%Get(ATMOS_DYN_NUMDIFF_LAPLAV_ID, field)
    call field%GetLocalMeshField(domID, NUMDIFF_LAPLAV) 
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
  end subroutine AtmosDynNumDiffTend_GetLocalMeshFields
  
  ! subroutine AtmosDynVars_GetLocalMeshFields_analysis( domID, mesh, analysis_list, &
  !   MOMZ_t, MOMZ_t_advx, MOMZ_t_advY, MOMZ_t_advZ, MOMZ_t_lift, MOMZ_t_buoy,       &
  !   lcmesh3D                                                             &
  !   )

  !   use scale_mesh_base, only: MeshBase
  !   use scale_meshfield_base, only: MeshFieldBase
  !   implicit none

  !   integer, intent(in) :: domID
  !   class(MeshBase), intent(in) :: mesh
  !   class(ModelVarManager), intent(inout) :: analysis_list
  !   class(LocalMeshFieldBase), pointer, intent(out) :: &
  !     MOMZ_t, MOMZ_t_advx, MOMZ_t_advY, MOMZ_t_advZ, MOMZ_t_lift, MOMZ_t_buoy
  !   class(LocalMesh3D), pointer, intent(out), optional :: lcmesh3D

  !   class(MeshFieldBase), pointer :: field   
  !   class(LocalMeshBase), pointer :: lcmesh
  !   !-------------------------------------------------------

  !   !--
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t)
  !   !---
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t_advx, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t_advx)
  !   !---
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t_advy, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t_advy)
  !   !---
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t_advz, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t_advz)
  !   !---
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t_lift, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t_lift)
  !   !---
  !   call analysis_list%Get(ATMOS_DYN_ANALYSISVAR_MOMZ_t_buoy, field)
  !   call field%GetLocalMeshField(domID, MOMZ_t_buoy)

    
  !   if (present(lcmesh3D)) then
  !     call mesh%GetLocalMesh( domID, lcmesh )
  !     nullify( lcmesh3D )

  !     select type(lcmesh)
  !     type is (LocalMesh3D)
  !       if (present(lcmesh3D)) lcmesh3D => lcmesh
  !     end select
  !   end if

  !   return
  ! end subroutine AtmosDynVars_GetLocalMeshFields_analysis


end module mod_atmos_dyn_vars