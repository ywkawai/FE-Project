#ifndef SCALE_ELEMENT_QUADRILATERAL_CBIND_H
#define SCALE_ELEMENT_QUADRILATERAL_CBIND_H

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void* CQuadrilateralElement_Init(int elemOrder, bool LumpedMassMatFlag);
    void  CQuadrilateralElement_Final(void* ptr);

    // Getter functions
    void CQuadrilateralElement_get_Np(void* ptr, int* val);
    void CQuadrilateralElement_get_Nfaces(void* ptr, int* val);
    void CQuadrilateralElement_get_Nfp(void* ptr, int* val);
    void CQuadrilateralElement_get_NfpTot(void* ptr, int* val);
    void CQuadrilateralElement_get_Nv(void* ptr, int* val);
    void CQuadrilateralElement_get_PolyOrder(void* ptr, int* val);   
    void CQuadrilateralElement_get_x1(void* ptr, double* val, int n);
    void CQuadrilateralElement_get_x2(void* ptr, double* val, int n);
    void CQuadrilateralElement_get_IntWeight_lgl(void* ptr, double* val, int n);
    void CQuadrilateralElement_get_Dx1(void* ptr, double* val, int nx, int ny);
    void CQuadrilateralElement_get_Dx2(void* ptr, double* val, int nx, int ny);
    void CQuadrilateralElement_get_Lift(void* ptr, double* val, int nx, int ny);
}
class QuadrilateralElement {
public:
    QuadrilateralElement(int elemOrder, bool lumpedMass)
        : handle_(cbind::init_handle(
            &CQuadrilateralElement_Init, &CQuadrilateralElement_Final, "Failed to initialize CQuadrilateralElement", 
            elemOrder, lumpedMass)){}
    ~QuadrilateralElement() = default;

    QuadrilateralElement(const QuadrilateralElement&) = delete;
    QuadrilateralElement& operator=(const QuadrilateralElement&) = delete;
    QuadrilateralElement(QuadrilateralElement&& other) noexcept = default;
    QuadrilateralElement& operator=(QuadrilateralElement&& other) = default;

    int get_Np() const { return cbind::get_value<int>(handle_, &CQuadrilateralElement_get_Np); }
    int get_Nfaces() const { return cbind::get_value<int>(handle_, &CQuadrilateralElement_get_Nfaces); }
    int get_Nfp() const { return cbind::get_value<int>(handle_, &CQuadrilateralElement_get_Nfp); }
    int get_NfpTot() const { return cbind::get_value<int>(handle_, &CQuadrilateralElement_get_NfpTot); }
    int get_Nv() const { return cbind::get_value<int>(handle_, &CQuadrilateralElement_get_Nv); }
    int get_PolyOrder() const { return cbind::get_value<int>(handle_, &CQuadrilateralElement_get_PolyOrder); }
    
    std::vector<double> get_x1() const { return cbind::fetch_vector(this->get_Handle(), &CQuadrilateralElement_get_x1, this->get_Np()); }
    std::vector<double> get_x2() const { return cbind::fetch_vector(this->get_Handle(), &CQuadrilateralElement_get_x2, this->get_Np()); }
    std::vector<double> get_IntWeight_lgl() const { return cbind::fetch_vector(this->get_Handle(), &CQuadrilateralElement_get_IntWeight_lgl, this->get_Np()); }

    std::vector<double> get_Dx1() const { return cbind::fetch_matrix(this->get_Handle(), &CQuadrilateralElement_get_Dx1, this->get_Np(), this->get_Np()); }
    std::vector<double> get_Dx2() const { return cbind::fetch_matrix(this->get_Handle(), &CQuadrilateralElement_get_Dx2, this->get_Np(), this->get_Np()); }
    std::vector<double> get_Lift() const { return cbind::fetch_matrix(this->get_Handle(), &CQuadrilateralElement_get_Lift, this->get_Np(), this->get_NfpTot()); }

    const cbind::Handle& get_Handle() const { return this->handle_; }
private:
    cbind::Handle handle_;
};
#endif