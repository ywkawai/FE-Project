#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../element/cbind/scale_element_base2d_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void CMeshBase2D_release_handle(void* ptr);
    // Getter functions
}
class MeshBase2D {
public:
    MeshBase2D(void* handle): 
        handle_(cbind::Handle(handle, cbind::Finalizer{&CMeshBase2D_release_handle})){
    }
    ~MeshBase2D() = default; 

    MeshBase2D(const MeshBase2D&) = delete;
    MeshBase2D& operator=(const MeshBase2D&) = delete;
    MeshBase2D(MeshBase2D&& other) noexcept = default;
    MeshBase2D& operator=(MeshBase2D&& other) = default;
private:
    cbind::Handle handle_;
};
