#include <iostream>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include "../../../src/element/cbind/scale_element_line_cbind.hpp"
#include "../../../src/mesh/cbind/scale_mesh_linedom1d_cbind.hpp"
#include "../../../src/data/cbind/scale_meshfield_base_cbind.hpp"

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

        MeshLineDom1D mesh(2, -1.0, 1.0, elem, 1, 1, 0);
        mesh.generate();
        LocalMesh1D lcmesh1D = mesh.get_LocalMesh1D(1);
        std::cout << "Lcmesh1D NeS=" << lcmesh1D.get_NeS() << std::endl;
        std::cout << "Lcmesh1D NeE=" << lcmesh1D.get_NeE() << std::endl;
        std::cout << "Lcmesh1D NeA=" << lcmesh1D.get_NeA() << std::endl;
        std::cout << "Lcmesh1D Ne=" << lcmesh1D.get_Ne() << std::endl;
        ElementBase1D refElem1D = lcmesh1D.get_refElem1D();
        std::cout << "Lcmesh1D->refElem1D Np=" << refElem1D.get_Np() << std::endl;


        auto mesh_ = mesh.get_MeshBase1D();
        MeshField1D q("q", "K", mesh_);
        LocalMeshField1D lmeshfield = q.get_LocalMeshField(1);
        auto val = lmeshfield.get_val_view();
        for (std::size_t j=0; j<lcmesh1D.get_Ne(); ++j){
            for (std::size_t i=0; i<Np; ++i){
                val(i,j) = i + j*Np;
            }
        }
        q.print_val();

    }catch(const std::exception& ex){
        std::cerr << ex.what() << "\n";
        return EXIT_FAILURE;
    }

    std::cout << "CLineElement finalized and resources released." << std::endl;
    return EXIT_SUCCESS;
}