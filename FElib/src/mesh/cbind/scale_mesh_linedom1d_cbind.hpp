#ifndef SCALE_MESH_LINEDOM1D_CBIND_H
#define SCALE_MESH_LINEDOM1D_CBIND_H

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../element/cbind/scale_element_line_cbind.hpp"
#include "../../mesh/cbind/scale_localmesh1d_cbind.hpp"
#include "../../mesh/cbind/scale_mesh_base1d_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void* CMeshLineDom1D_Init(int NeG, double dom_xmin, double dom_xmax, void* refElem_ptr, int NLocalMeshPerPrc, int nproc, int myrank, const double* FX_ptr);
    void  CMeshLineDom1D_Final(void* ptr);

    void  CMeshLineDom1D_Generate(void* ptr);

    // Getter functions
    void CMeshLineDom1D_get_NeG(void* ptr, int* val);
    void CMeshLineDom1D_get_Nprc(void* ptr, int* val);
    void CMeshLineDom1D_get_xmin_gl(void* ptr, double* val);
    void CMeshLineDom1D_get_xmax_gl(void* ptr, double* val);
    void* CMeshLineDom1D_get_lcmesh1D(void* ptr, int lcmeshID);

    void* CMeshLineDom1D_get_MeshBase1D(void* ptr);
}
class MeshLineDom1D {
public:
    MeshLineDom1D(int NeG, double dom_xmin, double dom_xmax, const LineElement& refElem, int NLocalMeshPerPrc, int nproc, int myrank, const std::vector<double>* FX = nullptr)
        : handle_(cbind::init_handle(
              &CMeshLineDom1D_Init, &CMeshLineDom1D_Final, "Failed to initialize CMeshLineDom1D",
              NeG, dom_xmin, dom_xmax, refElem.get_Handle().get(),
              NLocalMeshPerPrc, nproc, myrank,
              (FX && !FX->empty()) ? FX->data() : nullptr)){}

    ~MeshLineDom1D() = default;

    MeshLineDom1D(const MeshLineDom1D&) = delete;
    MeshLineDom1D& operator=(const MeshLineDom1D&) = delete;
    MeshLineDom1D(MeshLineDom1D&& other) noexcept = default;
    MeshLineDom1D& operator=(MeshLineDom1D&& other) = default;

    int get_NeG() const { return cbind::get_value<int>(handle_, &CMeshLineDom1D_get_NeG); }
    int get_Nprc() const { return cbind::get_value<int>(handle_, &CMeshLineDom1D_get_Nprc); }
    double get_xmin_gl() const { return cbind::get_value<double>(handle_, &CMeshLineDom1D_get_xmin_gl); }
    double get_xmax_gl() const { return cbind::get_value<double>(handle_, &CMeshLineDom1D_get_xmax_gl); }
    
    void generate() { CMeshLineDom1D_Generate(handle_.get()); }

    LocalMesh1D get_LocalMesh1D(int lcmeshID){
        void* lcmesh_handle = CMeshLineDom1D_get_lcmesh1D(handle_.get(), lcmeshID);
        return LocalMesh1D(lcmesh_handle);
    }

    MeshBase1D get_MeshBase1D(){
        void* meshbase1D_handle = CMeshLineDom1D_get_MeshBase1D(handle_.get());
        return MeshBase1D(meshbase1D_handle);
    }
private:
    cbind::Handle handle_;
};

#endif
