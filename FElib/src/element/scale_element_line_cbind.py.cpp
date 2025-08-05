#include "scale_element_line_cbind.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(ScaleFECore, m) {
    py::class_<LineElement>(m, "LineElement")
        .def(py::init<int, bool>())
        .def("get_Np", &LineElement::get_Np)
        .def("get_Nfaces", &LineElement::get_Nfaces)
        .def("get_NfpTot", &LineElement::get_NfpTot)
        .def("get_Nv", &LineElement::get_Nv)
        .def("get_PolyOrder", &LineElement::get_PolyOrder)
        .def("get_x1", &LineElement::get_x1)
        .def("get_Dx1", &LineElement::get_Dx1);
}