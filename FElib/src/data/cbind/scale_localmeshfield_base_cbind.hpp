#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <cassert>
#include <stdexcept>

#include <array>

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    struct f2dview {
        void* data;
        std::size_t  n1, n2;
        std::size_t  s0, s1;
    };
    struct col_major_view2d {
        double* ptr{};
        std::size_t  n1{}, n2{};
        std::size_t  s0{1}, s1{};

        double& operator()(std::size_t i, std::size_t j) const{
#ifdef DEBUG
            assert(i < n1 && j < n2);
#endif            
            return ptr[i * s0 + j * s1];
        }

        void print_view_info() const{
            std::cout << "n1=" << n1 << ", n2=" << n2 << std::endl;
            std::cout << "s0=" << s0 << ", s1=" << s1 << std::endl;
        }
    };

    void CLocalMeshField1D_release_handle(void* ptr);
    void CLocalMeshField1D_get_val_view(void* ptr, f2dview* out_view);

    void CLocalMeshField2D_release_handle(void* ptr);
    void CLocalMeshField2D_get_val_view(void* ptr, f2dview* out_view);

    void CLocalMeshField3D_release_handle(void* ptr);
    void CLocalMeshField3D_get_val_view(void* ptr, f2dview* out_view);
    // Getter functions
}

template<void(*ReleaseFn)(void*), void(*GetViewFn)(void*, f2dview*)>
class LocalMeshFieldT {
public:
  explicit LocalMeshFieldT(void* handle)
    : handle_(cbind::Handle(handle, cbind::Finalizer{ReleaseFn})) {}

  ~LocalMeshFieldT() = default;

  LocalMeshFieldT(const LocalMeshFieldT&)            = delete;
  LocalMeshFieldT& operator=(const LocalMeshFieldT&) = delete;
  LocalMeshFieldT(LocalMeshFieldT&&) noexcept        = default;
  LocalMeshFieldT& operator=(LocalMeshFieldT&&)      = default;

  inline col_major_view2d get_val_view() const {
    f2dview v{};
    GetViewFn(handle_.get(), &v);
    return col_major_view2d{
      static_cast<double*>(v.data),
      v.n1, v.n2, v.s0, v.s1
    };
  }

private:
  cbind::Handle handle_;
};

using LocalMeshField1D =
  LocalMeshFieldT<&CLocalMeshField1D_release_handle,
                  &CLocalMeshField1D_get_val_view>;

using LocalMeshField2D =
  LocalMeshFieldT<&CLocalMeshField2D_release_handle,
                  &CLocalMeshField2D_get_val_view>;

using LocalMeshField3D =
  LocalMeshFieldT<&CLocalMeshField3D_release_handle,
                  &CLocalMeshField3D_get_val_view>;
