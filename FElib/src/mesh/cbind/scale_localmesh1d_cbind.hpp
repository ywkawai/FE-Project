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
    void CLocalMesh1D_release_handle(void* ptr);

    // Getter functions
    void CLocalMesh1D_get_Ne(void* ptr, int* val);
    void CLocalMesh1D_get_NeS(void* ptr, int* val);
    void CLocalMesh1D_get_NeE(void* ptr, int* val);
    void CLocalMesh1D_get_NeA(void* ptr, int* val);
    void* CLocalMesh1D_get_refElem1D(void* ptr);
    void CLocalMesh1D_get_pos_en(void* ptr, double* val, int nx, int ny, int nz);
}
class LocalMesh1D {
public:
    LocalMesh1D(void* handle): 
        handle_(cbind::Handle(handle, cbind::Finalizer{&CLocalMesh1D_release_handle})){
    }
    ~LocalMesh1D() = default; 

    LocalMesh1D(const LocalMesh1D&) = delete;
    LocalMesh1D& operator=(const LocalMesh1D&) = delete;
    LocalMesh1D(LocalMesh1D&& other) noexcept = default;
    LocalMesh1D& operator=(LocalMesh1D&& other) = default;

    int get_Ne() const { return cbind::get_value<int>(handle_, &CLocalMesh1D_get_Ne); }
    int get_NeS() const { return cbind::get_value<int>(handle_, &CLocalMesh1D_get_NeS); }
    int get_NeE() const { return cbind::get_value<int>(handle_, &CLocalMesh1D_get_NeE); }
    int get_NeA() const { return cbind::get_value<int>(handle_, &CLocalMesh1D_get_NeA); }

    ElementBase1D get_refElem1D() const{
        return ElementBase1D(CLocalMesh1D_get_refElem1D(handle_.get()));
    }
    std::vector<double> get_pos_en() const {
        int Np = get_refElem1D().get_Np();
        int Ne = get_Ne();
        int ndim = 1;
        std::vector<double> data(Np * Ne * ndim);
        CLocalMesh1D_get_pos_en(handle_.get(), data.data(), Np, Ne, ndim);
        return data;
    }
private:
    cbind::Handle handle_;
};
