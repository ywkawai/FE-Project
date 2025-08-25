#include "scale_element_base1d_cbind.hpp"
#include "scale_element_base2d_cbind.hpp"
#include "scale_element_line_cbind.hpp"
#include "scale_element_quadrilateral_cbind.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void bind_element(py::module_ &m){
    py::class_<ElementBase1D>(m, "ElementBase1D")
        .def("get_Np", &ElementBase1D::get_Np)
        .def("get_PolyOrder", &ElementBase1D::get_PolyOrder);
    py::class_<ElementBase2D>(m, "ElementBase1D")
        .def("get_Np", &ElementBase2D::get_Np)
        .def("get_PolyOrder", &ElementBase2D::get_PolyOrder);
    py::class_<LineElement>(m, "LineElement")
        .def(py::init<int, bool>())
        .def("get_Np", &LineElement::get_Np)
        .def("get_Nfaces", &LineElement::get_Nfaces)
        .def("get_NfpTot", &LineElement::get_NfpTot)
        .def("get_Nv", &LineElement::get_Nv)
        .def("get_PolyOrder", &LineElement::get_PolyOrder)
        .def("get_x1", &LineElement::get_x1)
        .def("get_Dx1", &LineElement::get_Dx1);
    py::class_<QuadrilateralElement>(m, "QuadrilateralElement")
        .def(py::init<int, bool>())
        .def("get_Np", &QuadrilateralElement::get_Np)
        .def("get_Nfaces", &QuadrilateralElement::get_Nfaces)
        .def("get_Nfp", &QuadrilateralElement::get_Nfp)
        .def("get_NfpTot", &QuadrilateralElement::get_NfpTot)
        .def("get_Nv", &QuadrilateralElement::get_Nv)
        .def("get_PolyOrder", &QuadrilateralElement::get_PolyOrder)
        .def("get_x1", &QuadrilateralElement::get_x1)
        .def("get_x2", &QuadrilateralElement::get_x2)
        .def("get_Dx1", &QuadrilateralElement::get_Dx1)
        .def("get_Dx2", &QuadrilateralElement::get_Dx2);
} 
