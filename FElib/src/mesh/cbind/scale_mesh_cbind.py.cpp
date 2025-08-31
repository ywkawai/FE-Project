#include "scale_localmesh1d_cbind.hpp"
#include "scale_localmesh2d_cbind.hpp"
#include "scale_localmesh3d_cbind.hpp"
#include "scale_mesh_linedom1d_cbind.hpp"
#include "scale_mesh_rectdom2d_cbind.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void bind_mesh(py::module_ &m){
    py::class_<LocalMesh1D>(m, "LocalMesh1D")
        .def("get_Ne", &LocalMesh1D::get_Ne)
        .def("get_NeS", &LocalMesh1D::get_NeS)
        .def("get_NeE", &LocalMesh1D::get_NeE)
        .def("get_NeA", &LocalMesh1D::get_NeA)
        .def("get_refElem1D", &LocalMesh1D::get_refElem1D)
        .def("get_pos_en", &LocalMesh1D::get_pos_en);
    //
    py::class_<LocalMesh2D>(m, "LocalMesh2D")
        .def("get_Ne", &LocalMesh2D::get_Ne)
        .def("get_NeS", &LocalMesh2D::get_NeS)
        .def("get_NeE", &LocalMesh2D::get_NeE)
        .def("get_NeA", &LocalMesh2D::get_NeA)
        .def("get_refElem2D", &LocalMesh2D::get_refElem2D)
        .def("get_pos_en", &LocalMesh2D::get_pos_en);
    //
    py::class_<LocalMesh3D>(m, "LocalMesh3D")
        .def("get_Ne", &LocalMesh3D::get_Ne)
        .def("get_NeS", &LocalMesh3D::get_NeS)
        .def("get_NeE", &LocalMesh3D::get_NeE)
        .def("get_NeA", &LocalMesh3D::get_NeA)
        .def("get_refElem3D", &LocalMesh3D::get_refElem3D)
        .def("get_pos_en", &LocalMesh3D::get_pos_en);
    //
    py::class_<MeshBase1D>(m, "MeshBase1D");
    //
    py::class_<MeshLineDom1D>(m, "MeshLineDom1D")
        .def(py::init<int, double, double, const LineElement&, int, int, int, const std::vector<double>*>(),
          py::arg("NeG"),
          py::arg("dom_xmin"),
          py::arg("dom_xmax"),
          py::arg("refElem"),
          py::arg("NLocalMeshPerPrc"),
          py::arg("nproc"),
          py::arg("myrank"),
          py::arg("FX") = nullptr)       
        .def("generate", &MeshLineDom1D::generate)
        .def("get_NeG", &MeshLineDom1D::get_NeG)
        .def("get_Nprc", &MeshLineDom1D::get_Nprc)
        .def("get_xmin_gl", &MeshLineDom1D::get_xmin_gl)
        .def("get_xmax_gl", &MeshLineDom1D::get_xmax_gl)
        .def("get_LocalMesh1D", &MeshLineDom1D::get_LocalMesh1D)
        .def("get_MeshBase1D", &MeshLineDom1D::get_MeshBase1D);
    //
    py::class_<MeshRectDom2D>(m, "MeshRectDom2D")
        .def(py::init<int, int, double, double, double, double, bool, bool, const QuadrilateralElement&, int, int, int, int, int>(),
          py::arg("NeGX"),
          py::arg("NeGY"),
          py::arg("dom_xmin"),
          py::arg("dom_xmax"),
          py::arg("dom_ymin"),
          py::arg("dom_ymax"),
          py::arg("is_PeriodicX"),
          py::arg("is_PeriodicY"),
          py::arg("refElem"),
          py::arg("NLocalMeshPerPrc"),
          py::arg("NprcX"),
          py::arg("NprcY"),
          py::arg("nproc"),
          py::arg("myrank"))       
        .def("generate", &MeshRectDom2D::generate)
        .def("get_NeGX", &MeshRectDom2D::get_NeGX)
        .def("get_NeGY", &MeshRectDom2D::get_NeGY)
        .def("get_NprcX", &MeshRectDom2D::get_NprcX)
        .def("get_NprcY", &MeshRectDom2D::get_NprcY)
        .def("get_xmin_gl", &MeshRectDom2D::get_xmin_gl)
        .def("get_xmax_gl", &MeshRectDom2D::get_xmax_gl)
        .def("get_ymin_gl", &MeshRectDom2D::get_ymin_gl)
        .def("get_ymax_gl", &MeshRectDom2D::get_ymax_gl)
        .def("get_LocalMesh2D", &MeshRectDom2D::get_LocalMesh2D);
}