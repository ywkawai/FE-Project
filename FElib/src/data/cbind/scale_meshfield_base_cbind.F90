#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_meshfield_base_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    CBindingBase, &
    logical_c2f => cbinding_util_logical_c2f, &
    string_c2f => cbinding_util_string_c2f
  use scale_meshfield_base, only: &
    MeshField1D, MeshField2D, MeshField3D
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(CBindingBase) :: CMeshField1D
    type(MeshField1D) :: obj
  end type CMeshField1D
  type, public, extends(CBindingBase) :: CMeshField1DPtr
    type(MeshField1D), pointer :: obj
  end type CMeshField1DPtr  

  type, public, extends(CBindingBase) :: CMeshField2D
    type(MeshField2D) :: obj
  end type CMeshField2D
  type, public, extends(CBindingBase) :: CMeshField2DPtr
    type(MeshField2D), pointer :: obj
  end type CMeshField2DPtr  

  type, public, extends(CBindingBase) :: CMeshField3D
    type(MeshField3D) :: obj
  end type CMeshField3D
  type, public, extends(CBindingBase) :: CMeshField3DPtr
    type(MeshField3D), pointer :: obj
  end type CMeshField3DPtr  
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  function CMeshField1D_Init(varname, units, mesh_ptr, data_type ) result(ptr) bind(C, name="CMeshField1D_Init")
    use scale_mesh_base1d_cbind, only: CMeshBase1DPtr
    implicit none
    type(c_ptr), value :: varname
    type(c_ptr), value :: units
    type(c_ptr), value :: mesh_ptr
    integer(c_int), value :: data_type
    type(c_ptr) :: ptr

    type(CMeshField1D), pointer :: handle
    type(CMeshBase1DPtr), pointer :: mesh1D_fptr
    character(:), allocatable :: varname_f, units_f
    !------------------------------------
    
    call create_handle_1d(handle, ptr)
    call c_f_pointer( mesh_ptr, mesh1D_fptr )

    varname_f = string_c2f(varname)
    units_f = string_c2f(units)

    call handle%obj%Init(varname_f, units_f, mesh1D_fptr%obj, data_type)
    return
  end function CMeshField1D_Init

  subroutine CMeshField1D_Final( ptr ) bind(C, name="CMeshField1D_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle_1d(ptr)
    return
  end subroutine CMeshField1D_Final

  function CMeshField1D_GetLocalMeshField( ptr, domID ) result( lmeshfield) bind(C, name="CMeshField1D_GetLocalMeshField")
    use scale_localmeshfield_base_cbind, only: CLocalMeshField1DPtr
    implicit none
    type(c_ptr), value :: ptr
    integer(c_int), value :: domID
    type(c_ptr) :: lmeshfield

    type(CMeshField1D), pointer :: handle
    type(CLocalMeshField1DPtr), pointer :: lmeshfield_ptr
    !------------------------------------

    call c_f_pointer(ptr, handle)
    allocate(lmeshfield_ptr)
    lmeshfield_ptr%obj => handle%obj%local(domID)
    lmeshfield = c_loc(lmeshfield_ptr)
    return
  end function CMeshField1D_GetLocalMeshField

  subroutine CMeshField1D_print_val(ptr) bind(C, name="CMeshField1D_print_val")
    use scale_localmesh_1d, only: LocalMesh1D
    implicit none
    type(c_ptr), value :: ptr

    type(CMeshField1D), pointer :: handle
    integer :: lmeshID
    class(LocalMesh1D), pointer :: lmesh
    integer :: ke
    !-----------------------------------------

    call c_f_pointer(ptr, handle)
    do lmeshID=1, handle%obj%mesh%LOCAL_MESH_NUM
      lmesh => handle%obj%mesh%lcmesh_list(lmeshID)
      do ke=lmesh%NeS, lmesh%NeE
        write(*,*) "ke=",ke, ":", handle%obj%local(lmeshID)%val(:,ke) 
      end do
    end do
    return
  end subroutine CMeshField1D_print_val

!--
  function CMeshField2D_Init(varname, units, mesh_ptr, data_type ) result(ptr) bind(C, name="CMeshField2D_Init")
    use scale_mesh_base2d_cbind, only: CMeshBase2DPtr
    implicit none
    type(c_ptr), value :: varname
    type(c_ptr), value :: units
    type(c_ptr), value :: mesh_ptr
    integer(c_int), value :: data_type
    type(c_ptr) :: ptr

    type(CMeshField2D), pointer :: handle
    type(CMeshBase2DPtr), pointer :: mesh2D_fptr
    character(:), allocatable :: varname_f, units_f
    !------------------------------------
    call create_handle_2d(handle, ptr)
    call c_f_pointer( mesh_ptr, mesh2D_fptr )

    varname_f = string_c2f(varname)
    units_f = string_c2f(units)

    call handle%obj%Init(varname_f, units_f, mesh2D_fptr%obj, data_type)
    return
  end function CMeshField2D_Init

  subroutine CMeshField2D_Final( ptr ) bind(C, name="CMeshField2D_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle_2d(ptr)
    return
  end subroutine CMeshField2D_Final

!--
  function CMeshField3D_Init(varname, units, mesh_ptr, data_type ) result(ptr) bind(C, name="CMeshField3D_Init")
    use scale_mesh_base3d_cbind, only: CMeshBase3DPtr
    implicit none
    type(c_ptr), value :: varname
    type(c_ptr), value :: units
    type(c_ptr), value :: mesh_ptr
    integer(c_int), value :: data_type
    type(c_ptr) :: ptr

    type(CMeshField3D), pointer :: handle
    type(CMeshBase3DPtr), pointer :: mesh3D_fptr
    character(:), allocatable :: varname_f, units_f
    !------------------------------------
    call create_handle_3d(handle, ptr)
    call c_f_pointer( mesh_ptr, mesh3D_fptr )

    varname_f = string_c2f(varname)
    units_f = string_c2f(units)

    call handle%obj%Init(varname_f, units_f, mesh3D_fptr%obj, data_type)
    return
  end function CMeshField3D_Init

  subroutine CMeshField3D_Final( ptr ) bind(C, name="CMeshField3D_Final")
    implicit none
    type(c_ptr), value :: ptr
    !------------------------------------
    call destroy_handle_3d(ptr)
    return
  end subroutine CMeshField3D_Final

!*****
DEF_C_BIND(CMeshField1D,"CMeshField1D",_1D)
DEF_C_BIND(CMeshField2D,"CMeshField2D",_2D)
DEF_C_BIND(CMeshField3D,"CMeshField3D",_3D)
!***** Getter
!******

end module scale_meshfield_base_cbind