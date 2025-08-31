#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../element/cbind/scale_element_base1d_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void CMeshBase1D_release_handle(void* ptr);

    // Getter functions
    void CMeshBase1D_get_NeG(void* ptr, int* val);
}
class MeshBase1D {
public:
    MeshBase1D(void* handle): 
        handle_(cbind::Handle(handle, cbind::Finalizer{&CMeshBase1D_release_handle})){
    }
    ~MeshBase1D() = default; 

    MeshBase1D(const MeshBase1D&) = delete;
    MeshBase1D& operator=(const MeshBase1D&) = delete;
    MeshBase1D(MeshBase1D&& other) noexcept = default;
    MeshBase1D& operator=(MeshBase1D&& other) = default;

    int get_NeG() const { return cbind::get_value<int>(handle_, &CMeshBase1D_get_NeG); }

    const cbind::Handle& get_Handle() const { return this->handle_; }
private:
    cbind::Handle handle_;
};
