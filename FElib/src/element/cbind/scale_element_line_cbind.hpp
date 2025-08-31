#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void* CLineElement_Init(int elemOrder, bool LumpedMassMatFlag);
    void  CLineElement_Final(void* ptr);

    // Getter functions
    void CLineElement_get_Np(void* ptr, int* val);
    void CLineElement_get_Nfaces(void* ptr, int* val);
    void CLineElement_get_NfpTot(void* ptr, int* val);
    void CLineElement_get_Nv(void* ptr, int* val);
    void CLineElement_get_PolyOrder(void* ptr, int* val);   
    void CLineElement_get_x1(void* ptr, double* val, int n);
    void CLineElement_get_IntWeight_lgl(void* ptr, double* val, int n);
    void CLineElement_get_Dx1(void* ptr, double* val, int nx, int ny);
    void CLineElement_get_Lift(void* ptr, double* val, int nx, int ny);
}
class LineElement {
public:
    LineElement(int elemOrder, bool lumpedMass)
        : handle_(cbind::init_handle(
            &CLineElement_Init, &CLineElement_Final, "Failed to initialize CLineElement", 
            elemOrder, lumpedMass)){}
    ~LineElement() = default;

    LineElement(const LineElement&) = delete;
    LineElement& operator=(const LineElement&) = delete;
    LineElement(LineElement&& other) noexcept = default;
    LineElement& operator=(LineElement&& other) = default;

    int get_Np() const { return cbind::get_value<int>(handle_, &CLineElement_get_Np); }
    int get_Nfaces() const { return cbind::get_value<int>(handle_, &CLineElement_get_Nfaces); }
    int get_NfpTot() const { return cbind::get_value<int>(handle_, &CLineElement_get_NfpTot); }
    int get_Nv() const { return cbind::get_value<int>(handle_, &CLineElement_get_Nv); }
    int get_PolyOrder() const { return cbind::get_value<int>(handle_, &CLineElement_get_PolyOrder); }

    std::vector<double> get_x1() const { return cbind::fetch_vector(this->get_Handle(), &CLineElement_get_x1, this->get_Np()); }
    std::vector<double> get_IntWeight_lgl() const { return cbind::fetch_vector(this->get_Handle(), &CLineElement_get_IntWeight_lgl, this->get_Np()); }

    std::vector<double> get_Dx1() const { return cbind::fetch_matrix(this->get_Handle(), &CLineElement_get_Dx1, this->get_Np(), this->get_Np()); }
    std::vector<double> get_Lift() const { return cbind::fetch_matrix(this->get_Handle(), &CLineElement_get_Lift, this->get_Np(), this->get_NfpTot()); }

    const cbind::Handle& get_Handle() const { return this->handle_; }
private:
    cbind::Handle handle_;
};
