#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void CElementBase3D_release_handle(void* ptr);
    // Getter functions
    void CElementBase3D_get_Np(void* ptr, int* val);
    void CElementBase3D_get_Nnode_h1D(void* ptr, int* val);
    void CElementBase3D_get_Nnode_v(void* ptr, int* val);
    void CElementBase3D_get_PolyOrder_h(void* ptr, int* val);
    void CElementBase3D_get_PolyOrder_v(void* ptr, int* val);
}
class ElementBase3D {
public:
    ElementBase3D(void* handle): 
        handle_(cbind::Handle(handle, cbind::Finalizer{&CElementBase3D_release_handle})){
    }
    ~ElementBase3D() = default;

    ElementBase3D(const ElementBase3D&) = delete;
    ElementBase3D& operator=(const ElementBase3D&) = delete;
    ElementBase3D(ElementBase3D&& other) noexcept = default;
    ElementBase3D& operator=(ElementBase3D&& other) = default;

    int get_Np() const { return cbind::get_value<int>(handle_, &CElementBase3D_get_Np); }
    int get_Nnode_h1D() const { return cbind::get_value<int>(handle_, &CElementBase3D_get_Nnode_h1D); }
    int get_Nnode_v() const { return cbind::get_value<int>(handle_, &CElementBase3D_get_Nnode_v); }
    int get_PolyOrder_h() const { return cbind::get_value<int>(handle_, &CElementBase3D_get_PolyOrder_h); }
    int get_PolyOrder_v() const { return cbind::get_value<int>(handle_, &CElementBase3D_get_PolyOrder_v); }

private:
    cbind::Handle handle_;
};
