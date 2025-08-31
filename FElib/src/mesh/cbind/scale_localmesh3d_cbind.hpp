#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../element/cbind/scale_element_base3d_cbind.hpp"

extern "C" {
    void CLocalMesh3D_release_handle(void* ptr);

    // Getter functions
    void CLocalMesh3D_get_Ne(void* ptr, int* val);
    void CLocalMesh3D_get_NeS(void* ptr, int* val);
    void CLocalMesh3D_get_NeE(void* ptr, int* val);
    void CLocalMesh3D_get_NeA(void* ptr, int* val);
    void* CLocalMesh3D_get_refElem3D(void* ptr);
    void CLocalMesh3D_get_pos_en(void* ptr, double* val, int nx, int ny, int nz);
}
class LocalMesh3D {
public:
    LocalMesh3D(void* handle): 
        handle_(cbind::Handle(handle, cbind::Finalizer{&CLocalMesh3D_release_handle})){
    }
    ~LocalMesh3D() = default; 

    LocalMesh3D(const LocalMesh3D&) = delete;
    LocalMesh3D& operator=(const LocalMesh3D&) = delete;
    LocalMesh3D(LocalMesh3D&& other) noexcept = default;
    LocalMesh3D& operator=(LocalMesh3D&& other) = default;

    int get_Ne() const { return cbind::get_value<int>(handle_, &CLocalMesh3D_get_Ne); }
    int get_NeS() const { return cbind::get_value<int>(handle_, &CLocalMesh3D_get_NeS); }
    int get_NeE() const { return cbind::get_value<int>(handle_, &CLocalMesh3D_get_NeE); }
    int get_NeA() const { return cbind::get_value<int>(handle_, &CLocalMesh3D_get_NeA); }

    ElementBase3D get_refElem3D() const{
        return ElementBase3D(CLocalMesh3D_get_refElem3D(handle_.get()));
    }
    std::vector<double> get_pos_en() const {
        int Np = get_refElem3D().get_Np();
        int Ne = get_Ne();
        int ndim = 2;
        std::vector<double> data(Np * Ne * ndim);
        CLocalMesh3D_get_pos_en(handle_.get(), data.data(), Np, Ne, ndim);
        return data;
    }
private:
    cbind::Handle handle_;
};
