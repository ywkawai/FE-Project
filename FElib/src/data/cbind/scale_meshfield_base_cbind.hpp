#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "../../mesh/cbind/scale_mesh_base1d_cbind.hpp"
#include "../../mesh/cbind/scale_mesh_base2d_cbind.hpp"
#include "../../mesh/cbind/scale_mesh_base3d_cbind.hpp"
#include "../../data/cbind/scale_localmeshfield_base_cbind.hpp"
#include "../../data/cbind/scale_meshfield_base_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"

using namespace scale::FElib;

extern "C" {
    // Base 
    void  CMeshFieldBasePtr_release_handle(void* ptr);

    // 1D
    void* CMeshField1D_Init(const char* varname, const char* units, void* mesh_ptr, int data_type);
    void  CMeshField1D_Final(void* ptr);
    void* CMeshField1D_GetMeshFieldBase(void* ptr);
    void* CMeshField1D_GetLocalMeshField(void*, int domID);
    void  CMeshField1D_print_val(void*);

    // 2D
    void* CMeshField2D_Init(const char* varname, const char* units, void* mesh_ptr, int data_type);
    void  CMeshField2D_Final(void* ptr);
    void* CMeshField2D_GetMeshFieldBase(void* ptr);
    void* CMeshField2D_GetLocalMeshField(void*, int domID);
    void  CMeshField2D_print_val(void*);

    // 3D
    void* CMeshField3D_Init(const char* varname, const char* units, void* mesh_ptr, int data_type);
    void  CMeshField3D_Final(void* ptr);
    void* CMeshField3D_GetMeshFieldBase(void* ptr);
    void* CMeshField3D_GetLocalMeshField(void*, int domID);
    void  CMeshField3D_print_val(void*);
}

class MeshFieldBase{
public:
  MeshFieldBase(void* handle)
  : handle_( cbind::Handle(handle, cbind::Finalizer{&CMeshFieldBasePtr_release_handle})){}

  ~MeshFieldBase() = default;

  MeshFieldBase(const MeshFieldBase&)            = delete;
  MeshFieldBase& operator=(const MeshFieldBase&) = delete;
  MeshFieldBase(MeshFieldBase&&) noexcept        = default;
  MeshFieldBase& operator=(MeshFieldBase&&)      = default;

  const cbind::Handle& get_Handle() const { return handle_; }
private:
  cbind::Handle handle_;
};

namespace detail {
    using InitFn     = void* (*)(const char*, const char*, void*, int);
    using FinalFn    = void   (*)(void*);
    using GetMeshFieldBaseFn = void* (*)(void*);
    using GetLocalFn = void* (*)(void*, int);
    using PrintFn    = void   (*)(void*);
}

template<
  typename MeshBaseT,
  typename LocalFieldT,
  detail::InitFn     Init,
  detail::FinalFn    Final,
  detail::GetMeshFieldBaseFn GetMeshFieldBase,
  detail::GetLocalFn GetLocal,
  detail::PrintFn    Print,
  const char*        FailMsg
>
class MeshFieldT {
public:
  using MeshBase       = MeshBaseT;
  using LocalMeshField = LocalFieldT;

  MeshFieldT(const std::string& varname,
             const std::string& units,
             const MeshBase& mesh,
             int data_type = 1)
  : handle_( cbind::init_handle(
        Init, Final, FailMsg,
        varname.c_str(), units.c_str(),
        mesh.get_Handle().get(), data_type) ) {}

  ~MeshFieldT() = default;

  MeshFieldT(const MeshFieldT&)            = delete;
  MeshFieldT& operator=(const MeshFieldT&) = delete;
  MeshFieldT(MeshFieldT&&) noexcept        = default;
  MeshFieldT& operator=(MeshFieldT&&)      = default;

  MeshFieldBase get_MeshFieldBase() const{
    void* p = GetMeshFieldBase(handle_.get());
    return MeshFieldBase(p);
  }
  LocalMeshField get_LocalMeshField(int lcmeshID) const {
    void* p = GetLocal(handle_.get(), lcmeshID);
    return LocalMeshField(p);
  }

  void print_val() const { Print(handle_.get()); }

  const cbind::Handle& get_Handle() const { return handle_; }

private:
  cbind::Handle handle_;
};

inline constexpr char kFail1D[] = "Failed to initialize CMeshField1D";
inline constexpr char kFail2D[] = "Failed to initialize CMeshField2D";
inline constexpr char kFail3D[] = "Failed to initialize CMeshField3D";

using MeshField1D =
  MeshFieldT<MeshBase1D, LocalMeshField1D,
    &CMeshField1D_Init, &CMeshField1D_Final,
    &CMeshField1D_GetMeshFieldBase, CMeshField1D_GetLocalMeshField, &CMeshField1D_print_val, 
    kFail1D>;

using MeshField2D =
  MeshFieldT<MeshBase2D, LocalMeshField2D,
    &CMeshField2D_Init, &CMeshField2D_Final,
    &CMeshField2D_GetMeshFieldBase, &CMeshField2D_GetLocalMeshField, &CMeshField2D_print_val, 
    kFail2D>;

using MeshField3D =
  MeshFieldT<MeshBase3D, LocalMeshField3D,
    &CMeshField3D_Init, &CMeshField3D_Final,
    &CMeshField3D_GetMeshFieldBase, &CMeshField3D_GetLocalMeshField, &CMeshField3D_print_val, 
    kFail3D>;
