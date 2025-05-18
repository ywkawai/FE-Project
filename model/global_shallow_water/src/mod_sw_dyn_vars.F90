!-------------------------------------------------------------------------------
!> module Global SW / Dynamics
!!
!! @par Description
!!          Container for variables with dynamics
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_sw_dyn_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_element_base, only: ElementBase2D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D  

  use scale_meshfield_base, only: MeshField2D
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use scale_meshfieldcomm_cubedspheredom2d, only: MeshFieldCommCubedSphereDom2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  
  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo
  use scale_model_mesh_manager, only: ModelMeshBase

  use mod_sw_mesh, only: SWMesh
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, public :: SWDynVars
    type(MeshField2D), allocatable :: AUX_VARS2D(:)
    type(ModelVarManager) :: AUXVARS2D_manager

    type(MeshField2D), allocatable :: NUMDIFF_FLUX_VARS2D(:)
    type(ModelVarManager) :: NUMDIFF_FLUX_manager
    type(MeshFieldCommCubedSphereDom2D) :: NUMDIFF_FLUX_comm

    type(MeshField2D), allocatable :: NUMDIFF_TEND_VARS2D(:)
    type(ModelVarManager) :: NUMDIFF_TEND_manager
    type(MeshFieldCommCubedSphereDom2D) :: NUMDIFF_TEND_comm
  contains
    procedure :: Init => SWDynVars_Init
    procedure :: Final => SWDynVars_Final
  end type SWDynVars

  public :: SWDynAuxVars_GetLocalMeshFields
  public :: SWDynNumDiffFlux_GetLocalMeshFields
  public :: SWDynNumDiffTend_GetLocalMeshFields
  !public :: AtmosDynVars_GetLocalMeshFields_analysis

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, parameter :: ATMOS_DYN_AUXVARS2D_NUM          = 1
  integer, public, parameter :: ATMOS_DYN_AUXVARS2D_CORIOLIS_ID  = 1

  type(VariableInfo), public :: ATMOS_DYN_AUXVARS2D_VINFO(ATMOS_DYN_AUXVARS2D_NUM)
  DATA ATMOS_DYN_AUXVARS2D_VINFO / &
    VariableInfo( ATMOS_DYN_AUXVARS2D_CORIOLIS_ID, 'CORIOLIS', 'coriolis parameter',  &
                  's-1',  2, 'XY',  ''                                             )  / 

  integer, public, parameter :: SW_DYN_NUMDIFF_FLUX_NUM   = 2
  integer, public, parameter :: SW_DYN_NUMDIFFFLX_X_ID    = 1
  integer, public, parameter :: SW_DYN_NUMDIFFFLX_Y_ID    = 2

  integer, public, parameter :: SW_DYN_NUMDIFF_TEND_NUM   = 1
  integer, public, parameter :: SW_DYN_NUMDIFF_LAPLAH_ID  = 1

  type(VariableInfo), public :: SW_DYN_NUMDIFF_FLUX_VINFO(SW_DYN_NUMDIFF_FLUX_NUM)
  DATA SW_DYN_NUMDIFF_FLUX_VINFO / &
    VariableInfo( SW_DYN_NUMDIFFFLX_X_ID, 'DIFFFLX_X', 'flux in x-direction',  &
                  '?.m/s',  2, 'XY',  ''                                  ),   & 
    VariableInfo( SW_DYN_NUMDIFFFLX_Y_ID, 'DIFFFLX_Y', 'flux in y-direction',  &
                  '?.m/s',  2, 'XY',  ''                                  )    / 
                
  type(VariableInfo), public :: SW_DYN_NUMDIFF_TEND_VINFO(SW_DYN_NUMDIFF_TEND_NUM)
  DATA SW_DYN_NUMDIFF_TEND_VINFO / &
    VariableInfo( SW_DYN_NUMDIFF_LAPLAH_ID, 'NUMDIFF_LAPLAH', 'tendency due to nundiff',  &
                  '?/s',  2, 'XY',  ''                                                 )  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

contains
  subroutine SWDynVars_Init( this, model_mesh )
    implicit none
    class(SWDynVars), target, intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh

    integer :: v
    integer :: n
    logical :: reg_file_hist

    class(SWMesh), pointer :: sw_mesh
    class(MeshBase2D), pointer :: mesh2D

    !--------------------------------------------------

    LOG_INFO('SWDynVars_Init',*)

    !- Initialize auxiliary and diagnostic variables

    nullify( sw_mesh )
    select type(model_mesh)
    type is (SWMesh)
      sw_mesh => model_mesh
    end select
    
    mesh2D => sw_mesh%mesh

    !-
    call this%AUXVARS2D_manager%Init()
    allocate( this%AUX_VARS2D(ATMOS_DYN_AUXVARS2D_NUM) )

    reg_file_hist = .false.    
    do v = 1, ATMOS_DYN_AUXVARS2D_NUM
      call this%AUXVARS2D_manager%Regist(                   &
        ATMOS_DYN_AUXVARS2D_VINFO(v), mesh2D,               & ! (in) 
        this%AUX_VARS2D(v), reg_file_hist                   ) ! (out)
      
      do n = 1, sw_mesh%mesh%LOCAL_MESH_NUM
        this%AUX_VARS2D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do

    !-
    call this%NUMDIFF_FLUX_manager%Init()
    allocate( this%NUMDIFF_FLUX_VARS2D(SW_DYN_NUMDIFF_FLUX_NUM) )

    reg_file_hist = .false.    
    do v = 1, SW_DYN_NUMDIFF_FLUX_NUM
      call this%NUMDIFF_FLUX_manager%Regist(             &
        SW_DYN_NUMDIFF_FLUX_VINFO(v), mesh2D,            & ! (in) 
        this%NUMDIFF_FLUX_VARS2D(v), reg_file_hist       ) ! (out)
      
      do n = 1, sw_mesh%mesh%LOCAL_MESH_NUM
        this%NUMDIFF_FLUX_VARS2D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do
    call this%NUMDIFF_FLUX_comm%Init( SW_DYN_NUMDIFF_FLUX_NUM, 0, 0, sw_mesh%mesh )
    call this%NUMDIFF_FLUX_manager%MeshFieldComm_Prepair( this%NUMDIFF_FLUX_comm, this%NUMDIFF_FLUX_VARS2D(:) )

    !-
    call this%NUMDIFF_TEND_manager%Init()
    allocate( this%NUMDIFF_TEND_VARS2D(SW_DYN_NUMDIFF_TEND_NUM) )

    reg_file_hist = .false.    
    do v = 1, SW_DYN_NUMDIFF_TEND_NUM
      call this%NUMDIFF_TEND_manager%Regist(             &
        SW_DYN_NUMDIFF_TEND_VINFO(v), mesh2D,            & ! (in) 
        this%NUMDIFF_TEND_VARS2D(v), reg_file_hist       ) ! (out)
      
      do n = 1, sw_mesh%mesh%LOCAL_MESH_NUM
        this%NUMDIFF_TEND_VARS2D(v)%local(n)%val(:,:) = 0.0_RP
      end do         
    end do
    call this%NUMDIFF_TEND_comm%Init( SW_DYN_NUMDIFF_TEND_NUM, 0, 0, sw_mesh%mesh )
    call this%NUMDIFF_TEND_manager%MeshFieldComm_Prepair( this%NUMDIFF_TEND_comm, this%NUMDIFF_TEND_VARS2D(:) )
    
    return
  end subroutine SWDynVars_Init

  subroutine SWDynVars_Final( this )
    implicit none
    class(SWDynVars), target, intent(inout) :: this
    !---------------------------------------------------

    LOG_INFO('SWDynVars_Final',*)

    call this%AUXVARS2D_manager%Final()

    call this%NUMDIFF_FLUX_comm%Final()
    call this%NUMDIFF_FLUX_manager%Final()

    call this%NUMDIFF_TEND_comm%Final()
    call this%NUMDIFF_TEND_manager%Final()


    return
  end subroutine  SWDynVars_Final

  subroutine SWDynAuxVars_GetLocalMeshFields( domID, mesh, auxvars_list, &
    Coriolis,                                                            &
    lcmesh2D                                                             &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: Coriolis
    class(LocalMesh2D), pointer, intent(out), optional :: lcmesh2D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call auxvars_list%Get(ATMOS_DYN_AUXVARS2D_CORIOLIS_ID, field)
    call field%GetLocalMeshField(domID, Coriolis)
    !---
    
    if (present(lcmesh2D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh2D )

      select type(lcmesh)
      type is (LocalMesh2D)
        if (present(lcmesh2D)) lcmesh2D => lcmesh
      end select
    end if

    return
  end subroutine SWDynAuxVars_GetLocalMeshFields

  subroutine SWDynNumDiffFlux_GetLocalMeshFields( domID, mesh, auxvars_list, &
    NUMDIFF_FLUX_X, NUMDIFF_FLUX_Y,                                          &
    lcmesh2D                                                                 &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: NUMDIFF_FLUX_X
    class(LocalMeshFieldBase), pointer, intent(out) :: NUMDIFF_FLUX_Y
    class(LocalMesh2D), pointer, intent(out), optional :: lcmesh2D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call auxvars_list%Get(SW_DYN_NUMDIFFFLX_X_ID, field)
    call field%GetLocalMeshField(domID, NUMDIFF_FLUX_X)

    call auxvars_list%Get(SW_DYN_NUMDIFFFLX_Y_ID, field)
    call field%GetLocalMeshField(domID, NUMDIFF_FLUX_Y)

    !---
    
    if (present(lcmesh2D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh2D )

      select type(lcmesh)
      type is (LocalMesh2D)
        if (present(lcmesh2D)) lcmesh2D => lcmesh
      end select
    end if

    return
  end subroutine SWDynNumDiffFlux_GetLocalMeshFields

  subroutine SWDynNumDiffTend_GetLocalMeshFields( domID, mesh, auxvars_list, &
    NUMDIFF_LAPLAH,                                                          &
    lcmesh2D                                                                 &
    )

    use scale_mesh_base, only: MeshBase
    use scale_meshfield_base, only: MeshFieldBase
    implicit none

    integer, intent(in) :: domID
    class(MeshBase), intent(in) :: mesh
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(LocalMeshFieldBase), pointer, intent(out) :: NUMDIFF_LAPLAH
    class(LocalMesh2D), pointer, intent(out), optional :: lcmesh2D

    class(MeshFieldBase), pointer :: field   
    class(LocalMeshBase), pointer :: lcmesh
    !-------------------------------------------------------

    !--
    call auxvars_list%Get(SW_DYN_NUMDIFF_LAPLAH_ID, field)
    call field%GetLocalMeshField(domID, NUMDIFF_LAPLAH)

    !---
    
    if (present(lcmesh2D)) then
      call mesh%GetLocalMesh( domID, lcmesh )
      nullify( lcmesh2D )

      select type(lcmesh)
      type is (LocalMesh2D)
        if (present(lcmesh2D)) lcmesh2D => lcmesh
      end select
    end if

    return
  end subroutine SWDynNumDiffTend_GetLocalMeshFields

end module mod_sw_dyn_vars