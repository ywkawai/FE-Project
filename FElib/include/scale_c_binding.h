#ifndef SCALE_C_BINDING_H
#define SCALE_C_BINDING_H

#define DEF_C_BIND(TYPE,TYPE_STR) \
  subroutine create_handle(handle,ptr); \
    implicit none;              \
    type(c_ptr) :: ptr;         \
    type(TYPE), pointer :: handle; \
    allocate(handle);              \
    ptr = c_loc(handle);           \
    return;                     \
  end subroutine create_handle;      \
  subroutine destroy_handle(ptr); \
    implicit none;              \
    type(c_ptr), value :: ptr;  \
    type(TYPE), pointer :: handle; \
    call c_f_pointer(ptr, handle); \
    call handle%obj%Final();           \
    deallocate(handle);            \
    return;                     \
  end subroutine destroy_handle

#define DEF_C_BIND_RELEASE_HANDLE(TYPE,C_SUB_NAME) \
  subroutine release_handle(ptr) bind(C, name=C_SUB_NAME); \
    implicit none;              \
    type(c_ptr), value :: ptr;  \
    type(TYPE), pointer :: handle; \
    call c_f_pointer(ptr, handle); \
    deallocate(handle);            \
    return;                        \
  end subroutine release_handle

#define DEF_C_BIND_SUB(TYPE,SUB_NAME,C_SUB_NAME) \
  subroutine SUB_NAME (ptr) bind(C, name=C_SUB_NAME); \
    implicit none;                 \
    type(c_ptr), value :: ptr;     \
    type(TYPE), pointer :: handle; \
    call c_f_pointer(ptr, handle); \
    call handle%obj%SUB_NAME();    \
    return;                        \
  end subroutine SUB_NAME

#define DEF_C_BIND_GETTER_I(TYPE,SUB_NAME,C_SUB_NAME,VARNAME) \
  subroutine SUB_NAME (ptr, VARNAME) bind(C, name=C_SUB_NAME); \
    implicit none;              \
    type(c_ptr), value :: ptr;  \
    integer(c_int), intent(out) :: VARNAME; \
    type(TYPE), pointer :: handle; \
    call c_f_pointer(ptr, handle); \
    VARNAME = handle%obj%VARNAME; \
    return;                       \
  end subroutine SUB_NAME

#define DEF_C_BIND_GETTER_D(TYPE,SUB_NAME,C_SUB_NAME,VARNAME) \
  subroutine SUB_NAME (ptr, VARNAME) bind(C, name=C_SUB_NAME); \
    implicit none;              \
    type(c_ptr), value :: ptr;  \
    real(c_double), intent(out) :: VARNAME; \
    type(TYPE), pointer :: handle; \
    call c_f_pointer(ptr, handle); \
    VARNAME = handle%obj%VARNAME; \
    return;                       \
  end subroutine SUB_NAME

#define DEF_C_BIND_GETTER_Array1D_RP(TYPE,SUB_NAME,C_SUB_NAME,VARNAME) \
  subroutine SUB_NAME (ptr, VARNAME, n) bind(C, name=C_SUB_NAME); \
    implicit none;              \
    type(c_ptr), value :: ptr;  \
    type(c_ptr), value :: VARNAME; \
    integer(c_int), value :: n; \
    type(TYPE), pointer :: handle; \
    real(c_double), pointer :: array(:); \
    call c_f_pointer(ptr, handle); \
    call c_f_pointer(VARNAME, array, [n]); \
    array(:) = handle%obj%VARNAME(:); \
    return;                       \
  end subroutine SUB_NAME

#define DEF_C_BIND_GETTER_Array2D_RP(TYPE,SUB_NAME,C_SUB_NAME,VARNAME) \
  subroutine SUB_NAME (ptr, VARNAME, nx, ny) bind(C, name=C_SUB_NAME); \
    implicit none;              \
    type(c_ptr), value :: ptr;  \
    type(c_ptr), value :: VARNAME; \
    integer(c_int), value :: nx; \
    integer(c_int), value :: ny; \
    type(TYPE), pointer :: handle; \
    real(c_double), pointer :: array(:,:); \
    call c_f_pointer(ptr, handle); \
    call c_f_pointer(VARNAME, array, [nx,ny]); \
    array(:,:) = handle%obj%VARNAME(:,:); \
    return;                       \
  end subroutine SUB_NAME

#define DEF_C_BIND_GETTER_Array3D_RP(TYPE,SUB_NAME,C_SUB_NAME,VARNAME) \
  subroutine SUB_NAME (ptr, VARNAME, nx, ny, nz) bind(C, name=C_SUB_NAME); \
    implicit none;              \
    type(c_ptr), value :: ptr;  \
    type(c_ptr), value :: VARNAME; \
    integer(c_int), value :: nx; \
    integer(c_int), value :: ny; \
    integer(c_int), value :: nz; \
    type(TYPE), pointer :: handle; \
    real(c_double), pointer :: array(:,:,:); \
    call c_f_pointer(ptr, handle); \
    call c_f_pointer(VARNAME, array, [nx,ny,nz]); \
    array(:,:,:) = handle%obj%VARNAME(:,:,:); \
    return;                       \
  end subroutine SUB_NAME

#define DEF_C_BIND_MAKE_FOBJ_HANDLE(TYPE,SUB_NAME,C_SUB_NAME,VARTYPE,CVARTYPE,VARNAME) \
  function SUB_NAME (ptr) result(VARNAME) bind(C, name=C_SUB_NAME); \
    implicit none;              \
    type(c_ptr), value :: ptr;  \
    type(c_ptr) :: VARNAME; \
    type(TYPE), pointer :: handle; \
    type(CVARTYPE), pointer :: var_handle; \
    call c_f_pointer(ptr, handle); \
    allocate(var_handle); \
    var_handle%obj => handle%obj%VARNAME; \
    VARNAME = c_loc(var_handle); \
    return;                       \
  end function SUB_NAME
#endif