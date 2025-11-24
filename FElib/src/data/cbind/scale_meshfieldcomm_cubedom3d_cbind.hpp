#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../mesh/cbind/scale_mesh_cubedom3d_cbind.hpp"
#include "scale_meshfieldcomm_base_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void* CMeshFieldCommCubeDom3D_Init(int sfield_num, int hvfield_num, int htensorfield_num, void* mesh1D_ptr);
    void CMeshFieldCommCubeDom3D_Final(void* ptr);
    void CMeshFieldCommCubeDom3D_Exchange(void* ptr, bool do_wait);
    void CMeshFieldCommCubeDom3D_Put(void* ptr, void* field_list_ptr, int varnum, int varid_s);
    void CMeshFieldCommCubeDom3D_Get(void* ptr, void* field_list_ptr, int varnum, int varid_s);
}
class MeshFieldCommCubeDom3D {
public:
    MeshFieldCommCubeDom3D(int sfield_num, int hvfield_num, int htensorfield_num, const MeshCubeDom3D& mesh3D)
        : handle_(cbind::init_handle(
              &CMeshFieldCommCubeDom3D_Init, &CMeshFieldCommCubeDom3D_Final, "Failed to initialize CMeshFieldCommCubeDom3D",
              sfield_num, hvfield_num, htensorfield_num, mesh3D.get_Handle().get()) ){}

    ~MeshFieldCommCubeDom3D() = default;

    MeshFieldCommCubeDom3D(const MeshFieldCommCubeDom3D&) = delete;
    MeshFieldCommCubeDom3D& operator=(const MeshFieldCommCubeDom3D&) = delete;
    MeshFieldCommCubeDom3D(MeshFieldCommCubeDom3D&& other) noexcept = default;
    MeshFieldCommCubeDom3D& operator=(MeshFieldCommCubeDom3D&& other) = default;

    void exchange(bool do_wait) { CMeshFieldCommCubeDom3D_Exchange(handle_.get(), do_wait); }

    void put(const std::vector<MeshFieldContainer>& field_list, int varid_s){ 
        auto raw = MeshFieldContainer::to_void_ptrs(field_list);
        CMeshFieldCommCubeDom3D_Put(handle_.get(), raw.data(), raw.size(), varid_s); 
    }
    void put_py(const std::vector<MeshFieldContainer*>& field_list, int varid_s){ 
        auto raw = MeshFieldContainer::to_void_ptrs(field_list);
        CMeshFieldCommCubeDom3D_Put(handle_.get(), raw.data(), raw.size(), varid_s); 
    }
    void get(const std::vector<MeshFieldContainer>& field_list, int varid_s){ 
        auto raw = MeshFieldContainer::to_void_ptrs(field_list);
        CMeshFieldCommCubeDom3D_Get(handle_.get(), raw.data(), raw.size(), varid_s); 
    }
    void get_py(const std::vector<MeshFieldContainer*>& field_list, int varid_s){ 
        auto raw = MeshFieldContainer::to_void_ptrs(field_list);
        CMeshFieldCommCubeDom3D_Get(handle_.get(), raw.data(), raw.size(), varid_s); 
    }
private:
    cbind::Handle handle_;
};
