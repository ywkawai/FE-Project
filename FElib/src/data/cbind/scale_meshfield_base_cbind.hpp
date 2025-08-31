#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "../../mesh/cbind/scale_mesh_base1d_cbind.hpp"
#include "../../data/cbind/scale_localmeshfield_base_cbind.hpp"
#include "../../data/cbind/scale_meshfield_base_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"

using namespace scale::FElib;

extern "C" {
    void* CMeshField1D_Init(const char* varname, const char* units, void* mesh_ptr, int data_type);
    void  CMeshField1D_Final(void* ptr);

    void* CMeshField1D_GetLocalMeshField(void*, int domID);
    void CMeshField1D_print_val(void*);

    // Getter functions
}
class MeshField1D {
public:
    MeshField1D(const std::string& varname, const std::string& units, const MeshBase1D& mesh, int data_type=1)
        : handle_(cbind::init_handle(
              &CMeshField1D_Init, &CMeshField1D_Final, "Failed to initialize CMeshField1D",
              varname.data(), units.data(), mesh.get_Handle().get(), data_type)){}

    ~MeshField1D() = default;

    MeshField1D(const MeshField1D&) = delete;
    MeshField1D& operator=(const MeshField1D&) = delete;
    MeshField1D(MeshField1D&& other) noexcept = default;
    MeshField1D& operator=(MeshField1D&& other) = default;

    LocalMeshField1D get_LocalMeshField(int lcmeshID){
        void* lcmeshfield_handle = CMeshField1D_GetLocalMeshField(handle_.get(), lcmeshID);
        return LocalMeshField1D(lcmeshfield_handle);
    }

    void print_val(){ CMeshField1D_print_val(handle_.get()); }
private:
    cbind::Handle handle_;
};
