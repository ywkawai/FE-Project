#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../mesh/cbind/scale_mesh_linedom1d_cbind.hpp"
#include "scale_meshfieldcomm_base_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void* CMeshFieldComm1D_Init(int sfield_num, int hvfield_num, void* mesh1D_ptr);
    void CMeshFieldComm1D_Final(void* ptr);
    void CMeshFieldComm1D_Exchange(void* ptr, bool do_wait);
    void CMeshFieldComm1D_Put(void* ptr, void* field_list_ptr, int varnum, int varid_s);
    void CMeshFieldComm1D_Get(void* ptr, void* field_list_ptr, int varnum, int varid_s);
}
class MeshFieldComm1D {
public:
    MeshFieldComm1D(int sfield_num, int hvfield_num, const MeshLineDom1D& mesh1D)
        : handle_(cbind::init_handle(
              &CMeshFieldComm1D_Init, &CMeshFieldComm1D_Final, "Failed to initialize CMeshFieldComm1D",
              sfield_num, hvfield_num, mesh1D.get_Handle().get()) ){}

    ~MeshFieldComm1D() = default;

    MeshFieldComm1D(const MeshFieldComm1D&) = delete;
    MeshFieldComm1D& operator=(const MeshFieldComm1D&) = delete;
    MeshFieldComm1D(MeshFieldComm1D&& other) noexcept = default;
    MeshFieldComm1D& operator=(MeshFieldComm1D&& other) = default;

    void exchange(bool do_wait) { CMeshFieldComm1D_Exchange(handle_.get(), do_wait); }

    void put(const std::vector<MeshFieldContainer>& field_list, int varid_s){ 
        auto raw = MeshFieldContainer::to_void_ptrs(field_list);
        CMeshFieldComm1D_Put(handle_.get(), raw.data(), raw.size(), varid_s); 
    }
    void put_py(const std::vector<MeshFieldContainer*>& field_list, int varid_s){ 
        auto raw = MeshFieldContainer::to_void_ptrs(field_list);
        CMeshFieldComm1D_Put(handle_.get(), raw.data(), raw.size(), varid_s); 
    }
    void get(const std::vector<MeshFieldContainer>& field_list, int varid_s){ 
        auto raw = MeshFieldContainer::to_void_ptrs(field_list);
        CMeshFieldComm1D_Get(handle_.get(), raw.data(), raw.size(), varid_s); 
    }
    void get_py(const std::vector<MeshFieldContainer*>& field_list, int varid_s){ 
        auto raw = MeshFieldContainer::to_void_ptrs(field_list);
        CMeshFieldComm1D_Get(handle_.get(), raw.data(), raw.size(), varid_s); 
    }
private:
    cbind::Handle handle_;
};
