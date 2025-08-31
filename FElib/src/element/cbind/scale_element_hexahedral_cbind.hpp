#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void* CHexahedralElement_Init(int elemOrder_h, int elemOrder_v, bool LumpedMassMatFlag);
    void  CHexahedralElement_Final(void* ptr);

    // Getter functions
    void CHexahedralElement_get_Np(void* ptr, int* val);
    void CHexahedralElement_get_Nnode_h1D(void* ptr, int* val);
    void CHexahedralElement_get_Nnode_v(void* ptr, int* val);
    void CHexahedralElement_get_Nfaces(void* ptr, int* val);
    void CHexahedralElement_get_Nfaces_h(void* ptr, int* val);
    void CHexahedralElement_get_Nfaces_v(void* ptr, int* val);
    void CHexahedralElement_get_Nfp_h(void* ptr, int* val);
    void CHexahedralElement_get_Nfp_v(void* ptr, int* val);
    void CHexahedralElement_get_NfpTot(void* ptr, int* val);
    void CHexahedralElement_get_Nv(void* ptr, int* val);
    void CHexahedralElement_get_PolyOrder_h(void* ptr, int* val);   
    void CHexahedralElement_get_PolyOrder_v(void* ptr, int* val);   
    void CHexahedralElement_get_x1(void* ptr, double* val, int n);
    void CHexahedralElement_get_x2(void* ptr, double* val, int n);
    void CHexahedralElement_get_x3(void* ptr, double* val, int n);
    void CHexahedralElement_get_IntWeight_lgl(void* ptr, double* val, int n);
    void CHexahedralElement_get_Dx1(void* ptr, double* val, int nx, int ny);
    void CHexahedralElement_get_Dx2(void* ptr, double* val, int nx, int ny);
    void CHexahedralElement_get_Dx3(void* ptr, double* val, int nx, int ny);
    void CHexahedralElement_get_Lift(void* ptr, double* val, int nx, int ny);
}
class HexahedralElement {
public:
    HexahedralElement(int elemOrder_h, int elemOrder_v, bool lumpedMass)
        : handle_(cbind::init_handle(
            &CHexahedralElement_Init, &CHexahedralElement_Final, "Failed to initialize CHexahedralElement", 
            elemOrder_h, elemOrder_v, lumpedMass)){}
    ~HexahedralElement() = default;

    HexahedralElement(const HexahedralElement&) = delete;
    HexahedralElement& operator=(const HexahedralElement&) = delete;
    HexahedralElement(HexahedralElement&& other) noexcept = default;
    HexahedralElement& operator=(HexahedralElement&& other) = default;

    int get_Np() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_Np); }
    int get_Nnode_h1D() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_Nnode_h1D); }
    int get_Nnode_v() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_Nnode_v); }
    int get_Nfaces() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_Nfaces); }
    int get_Nfaces_h() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_Nfaces_h); }
    int get_Nfaces_v() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_Nfaces_v); }
    int get_Nfp_h() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_Nfp_h); }
    int get_Nfp_v() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_Nfp_v); }
    int get_NfpTot() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_NfpTot); }
    int get_Nv() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_Nv); }
    int get_PolyOrder_h() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_PolyOrder_h); }
    int get_PolyOrder_v() const { return cbind::get_value<int>(handle_, &CHexahedralElement_get_PolyOrder_v); }
    
    std::vector<double> get_x1() const { return cbind::fetch_vector(this->get_Handle(), &CHexahedralElement_get_x1, this->get_Np()); }
    std::vector<double> get_x2() const { return cbind::fetch_vector(this->get_Handle(), &CHexahedralElement_get_x2, this->get_Np()); }
    std::vector<double> get_x3() const { return cbind::fetch_vector(this->get_Handle(), &CHexahedralElement_get_x3, this->get_Np()); }
    std::vector<double> get_IntWeight_lgl() const { return cbind::fetch_vector(this->get_Handle(), &CHexahedralElement_get_IntWeight_lgl, this->get_Np()); }

    std::vector<double> get_Dx1() const { return cbind::fetch_matrix(this->get_Handle(), &CHexahedralElement_get_Dx1, this->get_Np(), this->get_Np()); }
    std::vector<double> get_Dx2() const { return cbind::fetch_matrix(this->get_Handle(), &CHexahedralElement_get_Dx2, this->get_Np(), this->get_Np()); }
    std::vector<double> get_Dx3() const { return cbind::fetch_matrix(this->get_Handle(), &CHexahedralElement_get_Dx3, this->get_Np(), this->get_Np()); }
    std::vector<double> get_Lift() const { return cbind::fetch_matrix(this->get_Handle(), &CHexahedralElement_get_Lift, this->get_Np(), this->get_NfpTot()); }

    const cbind::Handle& get_Handle() const { return this->handle_; }
private:
    cbind::Handle handle_;
};
