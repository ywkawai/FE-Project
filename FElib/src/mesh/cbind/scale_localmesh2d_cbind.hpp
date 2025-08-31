#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../element/cbind/scale_element_base2d_cbind.hpp"

extern "C" {
    void CLocalMesh2D_release_handle(void* ptr);

    // Getter functions
    void CLocalMesh2D_get_Ne(void* ptr, int* val);
    void CLocalMesh2D_get_NeS(void* ptr, int* val);
    void CLocalMesh2D_get_NeE(void* ptr, int* val);
    void CLocalMesh2D_get_NeA(void* ptr, int* val);
    void* CLocalMesh2D_get_refElem2D(void* ptr);
    void CLocalMesh2D_get_pos_en(void* ptr, double* val, int nx, int ny, int nz);
}
class LocalMesh2D {
public:
    LocalMesh2D(void* handle): 
        handle_(cbind::Handle(handle, cbind::Finalizer{&CLocalMesh2D_release_handle})){
    }
    ~LocalMesh2D() = default; 

    LocalMesh2D(const LocalMesh2D&) = delete;
    LocalMesh2D& operator=(const LocalMesh2D&) = delete;
    LocalMesh2D(LocalMesh2D&& other) noexcept = default;
    LocalMesh2D& operator=(LocalMesh2D&& other) = default;

    int get_Ne() const { return cbind::get_value<int>(handle_, &CLocalMesh2D_get_Ne); }
    int get_NeS() const { return cbind::get_value<int>(handle_, &CLocalMesh2D_get_NeS); }
    int get_NeE() const { return cbind::get_value<int>(handle_, &CLocalMesh2D_get_NeE); }
    int get_NeA() const { return cbind::get_value<int>(handle_, &CLocalMesh2D_get_NeA); }

    ElementBase2D get_refElem2D() const{
        return ElementBase2D(CLocalMesh2D_get_refElem2D(handle_.get()));
    }
    std::vector<double> get_pos_en() const {
        int Np = get_refElem2D().get_Np();
        int Ne = get_Ne();
        int ndim = 2;
        std::vector<double> data(Np * Ne * ndim);
        CLocalMesh2D_get_pos_en(handle_.get(), data.data(), Np, Ne, ndim);
        return data;
    }
private:
    cbind::Handle handle_;
};
