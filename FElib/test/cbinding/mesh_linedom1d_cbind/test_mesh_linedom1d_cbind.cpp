#include <iostream>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include "../../../src/common/cbind/scale_scalelib_cbind.hpp"
#include "../../../src/element/cbind/scale_element_line_cbind.hpp"
#include "../../../src/mesh/cbind/scale_mesh_linedom1d_cbind.hpp"
#include "../../../src/data/cbind/scale_meshfield_base_cbind.hpp"

#include "../../../src/data/cbind/scale_meshfieldcomm_base_cbind.hpp"
#include "../../../src/data/cbind/scale_meshfieldcomm_1d_cbind.hpp"

#include "../../../src/file/cbind/scale_file_base_meshfield_cbind.hpp"

int main() {
    int elemOrder = 3;
    bool lumpedMass = false;
    //---

    try{
        SCALElib_Init( "test_linedom1d_cbind" );

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
        MeshField1D q1("q1", "K", mesh_);
        MeshField1D q2("q2", "K", mesh_);

        LocalMeshField1D lmeshfield1 = q1.get_LocalMeshField(1);
        auto val1 = lmeshfield1.get_val_view();

        LocalMeshField1D lmeshfield2 = q2.get_LocalMeshField(1);
        auto val2 = lmeshfield2.get_val_view();

        for (std::size_t j=0; j<lcmesh1D.get_Ne(); ++j){
            for (std::size_t i=0; i<Np; ++i){
                val1(i,j) = i + j*Np;
                val2(i,j) = 10.0 * ( i + j*Np );
            }
        }
        q1.print_val();

        // Communication
        std::vector<MeshFieldContainer> field_list(2);
        field_list.at(0).set_field1D(q1);
        field_list.at(1).set_field1D(q2);
        MeshFieldComm1D comm(field_list.size(), 0, mesh);

        comm.put(field_list, 1);
        comm.exchange(true);
        comm.get(field_list, 1);

        // Output
        int vid_q1 = 1;
        int vid_q2 = 2;
        File_base_meshfield file(2, mesh.get_MeshBase1D());

        bool fileexisted;
        file.create("history", "Test cbinding", "REAL4", fileexisted);
        file.def_var( q1.get_MeshFieldBase(), "q1", vid_q1, MeshBase1D::DIMTYPEID_X, "REAL4", "q1", 
            1.0, 1);
        file.def_var( q2.get_MeshFieldBase(), "q2", vid_q2, MeshBase1D::DIMTYPEID_X, "REAL4", "q2", 
            1.0, 1);
        file.end_def();

        file.write_var( vid_q1, q1, 0.0, 1.0);
        file.write_var( vid_q2, q2, 0.0, 1.0);
        //
        SCALElib_Final();
    }catch(const std::exception& ex){
        std::cerr << ex.what() << "\n";
        return EXIT_FAILURE;
    }

    std::cout << "CLineElement finalized and resources released." << std::endl;
    return EXIT_SUCCESS;
}