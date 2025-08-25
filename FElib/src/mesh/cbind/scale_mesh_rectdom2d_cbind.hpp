#ifndef SCALE_MESH_RECTDOM2D_CBIND_H
#define SCALE_MESH_RECTDOM2D_CBIND_H

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../element/cbind/scale_element_quadrilateral_cbind.hpp"
#include "../../mesh/cbind/scale_localmesh2d_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void* CMeshRectDom2D_Init(int NeGX, int NeGY, double dom_xmin, double dom_xmax, double dom_ymin, double dom_ymax, bool is_PeriodicX, bool is_PeriodicY, void* refElem_ptr, int NLocalMeshPerPrc, int NprcX, int NprcY, int nproc, int myrank);
    void  CMeshRectDom2D_Final(void* ptr);
    void  CMeshRectDom2D_Generate(void* ptr);
    void* CMeshRectDom2D_get_lcmesh2D(void* ptr, int lcmeshID);

    // Getter functions
    void CMeshRectDom2D_get_NeGX(void* ptr, int* val);
    void CMeshRectDom2D_get_NeGY(void* ptr, int* val);
    void CMeshRectDom2D_get_NprcX(void* ptr, int* val);
    void CMeshRectDom2D_get_NprcY(void* ptr, int* val);
    void CMeshRectDom2D_get_xmin_gl(void* ptr, double* val);
    void CMeshRectDom2D_get_xmax_gl(void* ptr, double* val);
    void CMeshRectDom2D_get_ymin_gl(void* ptr, double* val);
    void CMeshRectDom2D_get_ymax_gl(void* ptr, double* val);
}
class MeshRectDom2D {
public:
    MeshRectDom2D(int NeGX, int NeGY, double dom_xmin, double dom_xmax, double dom_ymin, double dom_ymax, bool is_PeriodicX, bool is_PeriodicY, const QuadrilateralElement& refElem, int NLocalMeshPerPrc, int NprcX, int NprcY, int nproc, int myrank)
        : handle_(cbind::init_handle(
              &CMeshRectDom2D_Init, &CMeshRectDom2D_Final, "Failed to initialize CMeshRectDom2D",
              NeGX, NeGY, dom_xmin, dom_xmax, dom_ymin, dom_ymax,
              is_PeriodicX, is_PeriodicY, refElem.get_Handle().get(),
              NLocalMeshPerPrc, NprcX, NprcY, nproc, myrank)){}

    ~MeshRectDom2D() = default;

    MeshRectDom2D(const MeshRectDom2D&) = delete;
    MeshRectDom2D& operator=(const MeshRectDom2D&) = delete;
    MeshRectDom2D(MeshRectDom2D&& other) noexcept = default;
    MeshRectDom2D& operator=(MeshRectDom2D&& other) = default;

    int get_NeGX() const { return cbind::get_value<int>(handle_, &CMeshRectDom2D_get_NeGX); }
    int get_NeGY() const { return cbind::get_value<int>(handle_, &CMeshRectDom2D_get_NeGY); }
    int get_NprcX() const { return cbind::get_value<int>(handle_, &CMeshRectDom2D_get_NprcX); }
    int get_NprcY() const { return cbind::get_value<int>(handle_, &CMeshRectDom2D_get_NprcY); }
    double get_xmin_gl() const { return cbind::get_value<double>(handle_, &CMeshRectDom2D_get_xmin_gl); }
    double get_xmax_gl() const { return cbind::get_value<double>(handle_, &CMeshRectDom2D_get_xmax_gl); }
    double get_ymin_gl() const { return cbind::get_value<double>(handle_, &CMeshRectDom2D_get_ymin_gl); }
    double get_ymax_gl() const { return cbind::get_value<double>(handle_, &CMeshRectDom2D_get_ymax_gl); }
    
    void generate() { CMeshRectDom2D_Generate(handle_.get()); }
    LocalMesh2D get_LocalMesh2D(int lcmeshID){
        void* lcmesh_handle = CMeshRectDom2D_get_lcmesh2D(handle_.get(), lcmeshID);
        return LocalMesh2D(lcmesh_handle);
    }
private:
    cbind::Handle handle_;
};

#endif
