#ifndef SCALE_LOCALMESH2D_CBIND_H
#define SCALE_LOCALMESH2D_CBIND_H

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
    LocalMesh2D(void* handle): handle_(handle)
    {
        // handle_ = CLocalMesh2D_Init(NeG, dom_xmin, dom_xmax, refElem.get_Handle(), NLocalMeshPerPrc, nproc, myrank, FX_ptr);
        // if (!handle_) {
        //     throw std::runtime_error("Failed to initialize CLocalMesh2D");
        // }
    }
    ~LocalMesh2D() {
        if (handle_) {
//        std::cout << "Release handle=" << handle_ << std::endl; 
            CLocalMesh2D_release_handle(handle_);
        }
    }

    // Prohibit copying object
    LocalMesh2D(const LocalMesh2D&) = delete;
    LocalMesh2D& operator=(const LocalMesh2D&) = delete;
    // Allow to move the handle
    LocalMesh2D(LocalMesh2D&& other) noexcept : handle_(other.handle_) {
        other.handle_ = nullptr;
    }
    LocalMesh2D& operator=(LocalMesh2D&& other) noexcept {
        if (this != &other) {
            // if (handle_) CLocalMesh2D_Final(handle_);
            handle_ = other.handle_;
            other.handle_ = nullptr;
        }
        return *this;
    }

    int get_Ne() const { return getInt(&CLocalMesh2D_get_Ne); }
    int get_NeS() const { return getInt(&CLocalMesh2D_get_NeS); }
    int get_NeE() const { return getInt(&CLocalMesh2D_get_NeE); }
    int get_NeA() const { return getInt(&CLocalMesh2D_get_NeA); }

    ElementBase2D get_refElem2D() const{
        return ElementBase2D(CLocalMesh2D_get_refElem2D(handle_));
    }
    std::vector<double> get_pos_en() const {
        int Np = get_refElem2D().get_Np();
        int Ne = get_Ne();
        int ndim = 2;
        std::vector<double> data(Np * Ne * ndim);
        CLocalMesh2D_get_pos_en(handle_, data.data(), Np, Ne, ndim);
        return data;
    }
private:
    using GetterFuncI = void(*)(void*, int*);
    int getInt(GetterFuncI func) const {
        int val = 0;
        func(handle_, &val);
        return val;
    }
    using GetterFuncD = void(*)(void*, double*);
    double getDouble(GetterFuncD func) const {
        double val = 0;
        func(handle_, &val);
        return val;
    }

    void* handle_;
};

#endif