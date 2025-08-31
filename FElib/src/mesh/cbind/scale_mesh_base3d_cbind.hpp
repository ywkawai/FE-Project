#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../element/cbind/scale_element_base3d_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void CMeshBase3D_release_handle(void* ptr);
    // Getter functions
}
class MeshBase3D {
public:
    MeshBase3D(void* handle): 
        handle_(cbind::Handle(handle, cbind::Finalizer{&CMeshBase3D_release_handle})){
    }
    ~MeshBase3D() = default; 

    MeshBase3D(const MeshBase3D&) = delete;
    MeshBase3D& operator=(const MeshBase3D&) = delete;
    MeshBase3D(MeshBase3D&& other) noexcept = default;
    MeshBase3D& operator=(MeshBase3D&& other) = default;
private:
    cbind::Handle handle_;
};
