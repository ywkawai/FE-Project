!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_mesh
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_meshfield_base, only: MeshField3D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_element_base, only: ElementBase3D
  use scale_element_line, only: LineElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_sparsemat, only: sparsemat
  use scale_meshfield_base, only: MeshField2D
  use scale_mesh_topography, only: MeshTopography

  use scale_file_restart_meshfield, only: FILE_restart_meshfield_component
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_mesh_manager, only: ModelMesh3D
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, abstract, extends(ModelMesh3D), public :: AtmosMesh
    type(HexahedralElement) :: element
    type(LineElement) :: element_v1D

    type(MeshTopography) :: topography
    integer :: vcoord_type_id
  contains
    procedure :: AtmosMesh_Init
    procedure :: AtmosMesh_Final
    procedure(AtmosMesh_setup_restartfile1), public, deferred :: Setup_restartfile1
    procedure(AtmosMesh_setup_restartfile2), public, deferred :: Setup_restartfile2
    procedure(AtmosMesh_calc_UVMet), public, deferred :: Calc_UVmet
    generic :: Setup_restartfile => Setup_restartfile1, Setup_restartfile2
    procedure(AtmosMesh_setup_vcoord), public, deferred :: Setup_vcoordinate
    procedure :: Construct_ModalFilter3D => AtmosMesh_construct_ModalFilter3D
    procedure :: Construct_ModalFilterHV => AtmosMesh_construct_ModalFilterHV
  end type AtmosMesh

  interface
    subroutine AtmosMesh_setup_restartfile1( this, restart_file, var_num )
      import AtmosMesh
      import FILE_restart_meshfield_component
      class(AtmosMesh), target, intent(inout) :: this
      class(FILE_restart_meshfield_component), intent(inout) :: restart_file
      integer, intent(in) :: var_num  
    end subroutine AtmosMesh_setup_restartfile1
  end interface
  interface
    subroutine AtmosMesh_setup_restartfile2( this, restart_file, &
      in_basename, in_postfix_timelabel,                         &
      out_basename, out_postfix_timelabel,                       &
      out_dtype, out_title, var_num                              )
      import AtmosMesh
      import FILE_restart_meshfield_component
      class(AtmosMesh), target, intent(inout) :: this
      class(FILE_restart_meshfield_component), intent(inout) :: restart_file
      character(*), intent(in) :: in_basename
      logical, intent(in) :: in_postfix_timelabel
      character(*), intent(in) :: out_basename
      logical, intent(in) :: out_postfix_timelabel
      character(*), intent(in) :: out_title
      character(*), intent(in) :: out_dtype  
      integer, intent(in) :: var_num  
    end subroutine AtmosMesh_setup_restartfile2
  end interface
  interface
    subroutine AtmosMesh_calc_UVMet( this, U, V, &
        Umet, Vmet )
        import AtmosMesh
        import MeshField3D
        class(AtmosMesh), target, intent(in) :: this
        type(MeshField3D), intent(in) :: U
        type(MeshField3D), intent(in) :: V
        type(MeshField3D), intent(inout) :: Umet
        type(MeshField3D), intent(inout) :: Vmet
    end subroutine AtmosMesh_calc_UVMet
  end interface
  interface 
    subroutine AtmosMesh_setup_vcoord( this )
      import AtmosMesh
      class(AtmosMesh), target, intent(inout) :: this
    end subroutine AtmosMesh_setup_vcoord
  end interface
  integer, parameter, public :: ATM_MESH_MAX_COMMNUICATOR_NUM = 10

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains

  subroutine AtmosMesh_Init( this, mesh )

    use scale_FILE_monitor_meshfield, only: &
      FILE_monitor_meshfield_set_dim
    
    implicit none
    class(AtmosMesh), target, intent(inout) :: this
    class(MeshBase3D), intent(in) :: mesh

    character(len=H_SHORT) :: SpMV_storage_format = 'ELL' ! CSR or ELL
    class(MeshBase2D), pointer :: mesh2D
    !-------------------------------------------

    call this%ModelMesh3D_Init( mesh )

    call this%DOptrMat(1)%Init( mesh%refElem3D%Dx1, storage_format=SpMV_storage_format )
    call this%DOptrMat(2)%Init( mesh%refElem3D%Dx2, storage_format=SpMV_storage_format )
    call this%DOptrMat(3)%Init( mesh%refElem3D%Dx3, storage_format=SpMV_storage_format )

    call this%SOptrMat(1)%Init( mesh%refElem3D%Sx1, storage_format=SpMV_storage_format )
    call this%SOptrMat(2)%Init( mesh%refElem3D%Sx2, storage_format=SpMV_storage_format )
    call this%SOptrMat(3)%Init( mesh%refElem3D%Sx3, storage_format=SpMV_storage_format )

    call this%LiftOptrMat%Init( mesh%refElem3D%Lift, storage_format=SpMV_storage_format )

    !-
    call FILE_monitor_meshfield_set_dim( mesh, 'ATM3D' )
    
    !-
    call mesh%GetMesh2D( mesh2D )
    call this%topography%Init( "topo", mesh2D )

    return
  end subroutine AtmosMesh_Init

  subroutine AtmosMesh_Final(this)
    implicit none

    class(AtmosMesh), intent(inout) :: this
    !-------------------------------------------

    call this%topography%Final()
    call this%ModelMesh3D_Final()

    call this%element%Final()
    call this%element_v1D%Final()

    return
  end subroutine AtmosMesh_Final

  subroutine AtmosMesh_construct_ModalFilter3D( this, &
    filter,                                           &
    etac_h, alpha_h, ord_h,                           & 
    etac_v, alpha_v, ord_v                            )
    
    use scale_element_modalfilter, only: ModalFilter
    implicit none
    class(AtmosMesh), intent(in) :: this
    class(ModalFilter), intent(inout) :: filter
    real(RP), intent(in) :: etac_h
    real(RP), intent(in) :: alpha_h
    integer, intent(in) :: ord_h
    real(RP), intent(in) :: etac_v
    real(RP), intent(in) :: alpha_v
    integer, intent(in) :: ord_v
    !-------------------------------------------

    call filter%Init( this%element, &
      etac_h, alpha_h, ord_h,       & 
      etac_v, alpha_v, ord_v        )

    return
  end subroutine AtmosMesh_construct_ModalFilter3D

  subroutine AtmosMesh_construct_ModalFilterHV( this, &
    filterH3D, filterV1D,                             &
    etac_h, alpha_h, ord_h,                           & 
    etac_v, alpha_v, ord_v                            )
    
    use scale_element_modalfilter, only: ModalFilter
    use scale_element_line, only: LineElement
    implicit none
    class(AtmosMesh), intent(in) :: this
    class(ModalFilter), intent(inout) :: filterH3D
    class(ModalFilter), intent(inout) :: filterV1D
    real(RP), intent(in) :: etac_h
    real(RP), intent(in) :: alpha_h
    integer, intent(in) :: ord_h
    real(RP), intent(in) :: etac_v
    real(RP), intent(in) :: alpha_v
    integer, intent(in) :: ord_v

    type(LineElement) :: elemV1D
    !-------------------------------------------

    call filterH3D%Init( this%element, &
      etac_h, alpha_h, ord_h,          & 
      1.0_RP,  0.0_RP, ord_v           )

    call elemV1D%Init( this%element%PolyOrder_v, this%element%IsLumpedMatrix() )      
    call filterV1D%Init( elemV1D, &
      etac_v, alpha_v, ord_v,     &
      tend_flag = .true.          )
    call elemV1D%Final()

    return
  end subroutine AtmosMesh_construct_ModalFilterHV

end module mod_atmos_mesh