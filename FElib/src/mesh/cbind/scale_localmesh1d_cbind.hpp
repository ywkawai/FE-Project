#ifndef SCALE_LOCALMESH1D_CBIND_H
#define SCALE_LOCALMESH1D_CBIND_H

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../element/cbind/scale_element_base1d_cbind.hpp"

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
    LocalMesh1D(void* handle): handle_(handle)
    {
        // handle_ = CLocalMesh1D_Init(NeG, dom_xmin, dom_xmax, refElem.get_Handle(), NLocalMeshPerPrc, nproc, myrank, FX_ptr);
        // if (!handle_) {
        //     throw std::runtime_error("Failed to initialize CLocalMesh1D");
        // }
    }
    ~LocalMesh1D() {
        if (handle_) {
//        std::cout << "Release handle=" << handle_ << std::endl; 
            CLocalMesh1D_release_handle(handle_);
        }
    }

    // Prohibit copying object
    LocalMesh1D(const LocalMesh1D&) = delete;
    LocalMesh1D& operator=(const LocalMesh1D&) = delete;
    // Allow to move the handle
    LocalMesh1D(LocalMesh1D&& other) noexcept : handle_(other.handle_) {
        other.handle_ = nullptr;
    }
    LocalMesh1D& operator=(LocalMesh1D&& other) noexcept {
        if (this != &other) {
            // if (handle_) CLocalMesh1D_Final(handle_);
            handle_ = other.handle_;
            other.handle_ = nullptr;
        }
        return *this;
    }

    int get_Ne() const { return getInt(&CLocalMesh1D_get_Ne); }
    int get_NeS() const { return getInt(&CLocalMesh1D_get_NeS); }
    int get_NeE() const { return getInt(&CLocalMesh1D_get_NeE); }
    int get_NeA() const { return getInt(&CLocalMesh1D_get_NeA); }

    ElementBase1D get_refElem1D() const{
        return ElementBase1D(CLocalMesh1D_get_refElem1D(handle_));
    }
    std::vector<double> get_pos_en() const {
        int Np = get_refElem1D().get_Np();
        int Ne = get_Ne();
        int ndim = 1;
        std::vector<double> data(Np * Ne * ndim);
        CLocalMesh1D_get_pos_en(handle_, data.data(), Np, Ne, ndim);
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