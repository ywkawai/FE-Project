#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void bind_common(py::module_ &);
void bind_element(py::module_ &);
void bind_mesh(py::module_ &);
void bind_data(py::module_ &);

PYBIND11_MODULE(ScaleFECore, m) {
    m.doc() = "ScaleFECore";
    //
    py::module_ common = m.def_submodule("common", "");
    bind_common(common);
    //
    py::module_ element = m.def_submodule("element", "");
    bind_element(element);
    //
    py::module_ mesh = m.def_submodule("mesh", "");
    bind_mesh(mesh);
    //
    py::module_ data = m.def_submodule("data", "");
    bind_data(data);

}