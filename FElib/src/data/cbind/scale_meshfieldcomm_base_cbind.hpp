#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "scale_meshfield_base_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void* CMeshFieldContainer_Init();
    void CMeshFieldContainer_release_handle(void* ptr);
    void CMeshFieldContainer_SetField1D(void* ptr, void* field_ptr);
    void CMeshFieldContainer_SetField2D(void* ptr, void* field_ptr);
    void CMeshFieldContainer_SetField3D(void* ptr, void* field_ptr);
}

class MeshFieldContainer{
public:
    MeshFieldContainer()
        : handle_(cbind::init_handle(
              &CMeshFieldContainer_Init, &CMeshFieldContainer_release_handle, "Failed to initialize CMeshFieldContainer") ){}

    MeshFieldContainer(const MeshFieldContainer&) = delete;
    MeshFieldContainer& operator=(const MeshFieldContainer&) = delete;
    MeshFieldContainer(MeshFieldContainer&& other) noexcept = default;
    MeshFieldContainer& operator=(MeshFieldContainer&& other) = default;

    void set_field1D(const MeshField1D& field){ CMeshFieldContainer_SetField1D(handle_.get(), field.get_Handle().get()); }

    const cbind::Handle& get_Handle() const { return this->handle_; }

    static std::vector<void*> to_void_ptrs(const std::vector<MeshFieldContainer>& v){
        std::vector<void*> out(v.size());
        for (size_t i=0;i<v.size();++i) out[i] = v[i].get_Handle().get();
        return out;
    }
    static std::vector<void*> to_void_ptrs(const std::vector<MeshFieldContainer*>& v){
        std::vector<void*> out(v.size());
        for (size_t i=0;i<v.size();++i) out[i] = v[i]->get_Handle().get();
        return out;
    }
private:
    cbind::Handle handle_;    
};