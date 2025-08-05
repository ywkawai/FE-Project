#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

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
    void CLineElement_get_Dx1(void* ptr, double* val, int nx, int ny);
}
class LineElement {
public:
    LineElement(int elemOrder, bool lumpedMass)
        : handle_(nullptr)
    {
        handle_ = CLineElement_Init(elemOrder, lumpedMass);
        if (!handle_) {
            throw std::runtime_error("Failed to initialize CLineElement");
        }
    }
    ~LineElement() {
        if (handle_) {
            // std::cout << "Release resource.." << std::endl; 
            CLineElement_Final(handle_);
        }
    }

    // Prohibit copying object
    LineElement(const LineElement&) = delete;
    LineElement& operator=(const LineElement&) = delete;
    // Allow to move the handle
    LineElement(LineElement&& other) noexcept : handle_(other.handle_) {
        other.handle_ = nullptr;
    }
    LineElement& operator=(LineElement&& other) noexcept {
        if (this != &other) {
            if (handle_) CLineElement_Final(handle_);
            handle_ = other.handle_;
            other.handle_ = nullptr;
        }
        return *this;
    }

    int get_Np() const { return getInt(&CLineElement_get_Np); }
    int get_Nfaces() const { return getInt(&CLineElement_get_Nfaces); }
    int get_NfpTot() const { return getInt(&CLineElement_get_NfpTot); }
    int get_Nv() const { return getInt(&CLineElement_get_Nv); }
    int get_PolyOrder() const { return getInt(&CLineElement_get_PolyOrder); }
    
    std::vector<double> get_x1() const {
        int Np = get_Np();
        std::vector<double> data(Np);
        CLineElement_get_x1(handle_, data.data(), Np);
        return data;
    }
    std::vector<double> get_Dx1() const {
        int Np = get_Np();
        std::vector<double> data(Np * Np);
        CLineElement_get_Dx1(handle_, data.data(), Np, Np);
        return data;
    }

private:
    using GetterFunc = void(*)(void*, int*);
    int getInt(GetterFunc func) const {
        int val = 0;
        func(handle_, &val);
        return val;
    }

    void* handle_;
};
