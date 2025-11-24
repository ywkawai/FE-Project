#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../element/cbind/scale_element_hexahedral_cbind.hpp"
#include "../../mesh/cbind/scale_localmesh3d_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void* CMeshCubeDom3D_Init(int NeGX, int NeGY, int NeGZ, 
        double dom_xmin, double dom_xmax, double dom_ymin, double dom_ymax, double dom_zmin, double dom_zmax, 
        bool is_PeriodicX, bool is_PeriodicY, bool is_PeriodicZ, void* refElem_ptr, int NLocalMeshPerPrc, 
        int NprcX, int NprcY, int nproc, int myrank, const double* FZ_ptr);
    void  CMeshCubeDom3D_Final(void* ptr);
    void  CMeshCubeDom3D_Generate(void* ptr);
    void* CMeshCubeDom3D_get_lcmesh3D(void* ptr, int lcmeshID);

    // Getter functions
    void CMeshCubeDom3D_get_NeGX(void* ptr, int* val);
    void CMeshCubeDom3D_get_NeGY(void* ptr, int* val);
    void CMeshCubeDom3D_get_NeGZ(void* ptr, int* val);
    void CMeshCubeDom3D_get_NprcX(void* ptr, int* val);
    void CMeshCubeDom3D_get_NprcY(void* ptr, int* val);
    void CMeshCubeDom3D_get_xmin_gl(void* ptr, double* val);
    void CMeshCubeDom3D_get_xmax_gl(void* ptr, double* val);
    void CMeshCubeDom3D_get_ymin_gl(void* ptr, double* val);
    void CMeshCubeDom3D_get_ymax_gl(void* ptr, double* val);
    void CMeshCubeDom3D_get_zmin_gl(void* ptr, double* val);
    void CMeshCubeDom3D_get_zmax_gl(void* ptr, double* val);
}
class MeshCubeDom3D {
public:
    MeshCubeDom3D(int NeGX, int NeGY, int NeGZ, 
        double dom_xmin, double dom_xmax, double dom_ymin, double dom_ymax, double dom_zmin, double dom_zmax, 
        bool is_PeriodicX, bool is_PeriodicY, bool is_PeriodicZ, const HexahedralElement& refElem, int NLocalMeshPerPrc, 
        int NprcX, int NprcY, int nproc, int myrank, const std::vector<double>* FZ = nullptr)
        : handle_(cbind::init_handle(
              &CMeshCubeDom3D_Init, &CMeshCubeDom3D_Final, "Failed to initialize CMeshCubeDom3D",
              NeGX, NeGY, NeGZ, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, 
              is_PeriodicX, is_PeriodicY, is_PeriodicZ, refElem.get_Handle().get(),
              NLocalMeshPerPrc, NprcX, NprcY, nproc, myrank,
              (FZ && !FZ->empty()) ? FZ->data() : nullptr)){}              

    ~MeshCubeDom3D() = default;

    MeshCubeDom3D(const MeshCubeDom3D&) = delete;
    MeshCubeDom3D& operator=(const MeshCubeDom3D&) = delete;
    MeshCubeDom3D(MeshCubeDom3D&& other) noexcept = default;
    MeshCubeDom3D& operator=(MeshCubeDom3D&& other) = default;

    int get_NeGX() const { return cbind::get_value<int>(handle_, &CMeshCubeDom3D_get_NeGX); }
    int get_NeGY() const { return cbind::get_value<int>(handle_, &CMeshCubeDom3D_get_NeGY); }
    int get_NprcX() const { return cbind::get_value<int>(handle_, &CMeshCubeDom3D_get_NprcX); }
    int get_NprcY() const { return cbind::get_value<int>(handle_, &CMeshCubeDom3D_get_NprcY); }
    double get_xmin_gl() const { return cbind::get_value<double>(handle_, &CMeshCubeDom3D_get_xmin_gl); }
    double get_xmax_gl() const { return cbind::get_value<double>(handle_, &CMeshCubeDom3D_get_xmax_gl); }
    double get_ymin_gl() const { return cbind::get_value<double>(handle_, &CMeshCubeDom3D_get_ymin_gl); }
    double get_ymax_gl() const { return cbind::get_value<double>(handle_, &CMeshCubeDom3D_get_ymax_gl); }
    double get_zmin_gl() const { return cbind::get_value<double>(handle_, &CMeshCubeDom3D_get_zmin_gl); }
    double get_zmax_gl() const { return cbind::get_value<double>(handle_, &CMeshCubeDom3D_get_zmax_gl); }
    
    void generate() { CMeshCubeDom3D_Generate(handle_.get()); }
    LocalMesh3D get_LocalMesh3D(int lcmeshID){
        void* lcmesh_handle = CMeshCubeDom3D_get_lcmesh3D(handle_.get(), lcmeshID);
        return LocalMesh3D(lcmesh_handle);
    }

    const cbind::Handle& get_Handle() const { return this->handle_; }
private:
    cbind::Handle handle_;
};
