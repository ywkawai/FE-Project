#ifndef SCALE_ELEMENT_BASE2D_CBIND_H
#define SCALE_ELEMENT_BASE2D_CBIND_H

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void CElementBase2D_release_handle(void* ptr);
    // Getter functions
    void CElementBase2D_get_Np(void* ptr, int* val);
    void CElementBase2D_get_PolyOrder(void* ptr, int* val);
}
class ElementBase2D {
public:
    ElementBase2D(void* handle): 
        handle_(cbind::Handle(handle, cbind::Finalizer{&CElementBase2D_release_handle})){
    }
    ~ElementBase2D() = default;

    ElementBase2D(const ElementBase2D&) = delete;
    ElementBase2D& operator=(const ElementBase2D&) = delete;
    ElementBase2D(ElementBase2D&& other) noexcept = default;
    ElementBase2D& operator=(ElementBase2D&& other) = default;

    int get_Np() const { return cbind::get_value<int>(handle_, &CElementBase2D_get_Np); }
    int get_PolyOrder() const { return cbind::get_value<int>(handle_, &CElementBase2D_get_PolyOrder); }

private:
    cbind::Handle handle_;
};

#endif