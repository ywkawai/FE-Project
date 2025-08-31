#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void CElementBase1D_release_handle(void* ptr);
    // Getter functions
    void CElementBase1D_get_Np(void* ptr, int* val);
    void CElementBase1D_get_PolyOrder(void* ptr, int* val);
}
class ElementBase1D {
public:
    ElementBase1D(void* handle): 
        handle_(cbind::Handle(handle, cbind::Finalizer{&CElementBase1D_release_handle})){
    }
    ~ElementBase1D() = default;

    ElementBase1D(const ElementBase1D&) = delete;
    ElementBase1D& operator=(const ElementBase1D&) = delete;
    ElementBase1D(ElementBase1D&& other) noexcept = default;
    ElementBase1D& operator=(ElementBase1D&& other) = default;

    int get_Np() const { return cbind::get_value<int>(handle_, &CElementBase1D_get_Np); }
    int get_PolyOrder() const { return cbind::get_value<int>(handle_, &CElementBase1D_get_PolyOrder); }
private:
    cbind::Handle handle_;
};
