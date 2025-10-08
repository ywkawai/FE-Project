
#include "../../mesh/cbind/scale_mesh_base1d_cbind.hpp"
#include "../../data/cbind/scale_meshfield_base_cbind.hpp"
#include "scale_file_base_meshfield_cbind.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void bind_file(py::module_ &m){
    py::class_<File_base_meshfield>(m, "File_base_meshfield")
        .def(py::init<int, const MeshBase1D&, bool>(), py::arg("varnum"), py::arg("mesh1D"), py::arg("force_uniform_grid") = false)
        .def("open", &File_base_meshfield::open)
        .def("create", py::overload_cast<const std::string&,const std::string&,const std::string&,bool&,int,const std::string&,const std::string&>(&File_base_meshfield::create))
        .def("create", py::overload_cast<const std::string&,const std::string&,const std::string&,bool&,int>(&File_base_meshfield::create), 
            py::arg("basename"), py::arg("title"), py::arg("dtype"),  py::arg("fileexisted"), py::arg("myrank")=-1)
        .def("close", &File_base_meshfield::close)
        .def("def_var", &File_base_meshfield::def_var)
        .def("end_def", &File_base_meshfield::end_def)
        .def("write_var", &File_base_meshfield::write_var);
}
