#include <iostream>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include "../../../src/element/scale_element_line_cbind.hpp"

int main() {
    int elemOrder = 3;
    bool lumpedMass = false;
    //---

    try{
        LineElement elem(elemOrder, lumpedMass);
        int Np = elem.get_Np();
        auto x1 = elem.get_x1();
        auto Dx1 = elem.get_Dx1();

        std::cout << "Np        = " << Np << std::endl;
        std::cout << "Nfaces    = " << elem.get_Nfaces() << std::endl;
        std::cout << "NfpTot    = " << elem.get_NfpTot() << std::endl;
        std::cout << "Nv        = " << elem.get_Nv() << std::endl;
        std::cout << "PolyOrder = " << elem.get_PolyOrder() << std::endl;
        for (int i = 0; i < Np; ++i) {
            std::cout << "x1[" << i << "] = " << x1[i] << "\n";
        }
        for (int i = 0; i < Np; ++i) {
            std::cout << "D1[" << i <<","<< "*]";
            for (int j = 0; j < Np; ++j) {
                std::cout << Dx1[i+j*Np] << ", ";
            }
            std::cout << std::endl;
        }
    }catch(const std::exception& ex){
        std::cerr << ex.what() << "\n";
        return EXIT_FAILURE;
    }

    std::cout << "CLineElement finalized and resources released." << std::endl;
    return EXIT_SUCCESS;
}