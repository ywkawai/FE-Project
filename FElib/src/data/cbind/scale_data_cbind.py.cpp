
#include "../../mesh/cbind/scale_mesh_base1d_cbind.hpp"
#include "../../mesh/cbind/scale_mesh_linedom1d_cbind.hpp"
#include "scale_localmeshfield_base_cbind.hpp"
#include "scale_meshfield_base_cbind.hpp"
#include "scale_meshfieldcomm_base_cbind.hpp"
#include "scale_meshfieldcomm_1d_cbind.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

static inline py::array make_numpy_view(const col_major_view2d& v, py::handle base) {
    if (!v.ptr || v.n1==0 || v.n2==0) throw std::runtime_error("empty view");
    const py::ssize_t shape[2]   = { (py::ssize_t)v.n1, (py::ssize_t)v.n2 };
    const py::ssize_t strides[2] = {
        (py::ssize_t)(v.s0 * sizeof(double)),
        (py::ssize_t)(v.s1 * sizeof(double))
    };
    return py::array(
        py::buffer_info(
            v.ptr,
            sizeof(double),
            py::format_descriptor<double>::format(),
            2,
            { shape[0], shape[1] },
            { strides[0], strides[1] }
        ),
        base
    );
}
template<class Field>
static void pybind_LocalMeshField(py::module_& m, const char* pyname) {
    py::class_<Field>(m, pyname)
        .def("get_val_view", [](const Field& self){
            const auto v = self.get_val_view();
            return make_numpy_view(v, py::cast(self));
        }, R"doc(Return a NumPy view (no copy) into the Fortran buffer.)doc");
}

void bind_data(py::module_ &m){
    pybind_LocalMeshField<LocalMeshField1D>(m, "LocalMeshField1D");
    pybind_LocalMeshField<LocalMeshField2D>(m, "LocalMeshField2D");
    pybind_LocalMeshField<LocalMeshField3D>(m, "LocalMeshField3D");
    py::class_<MeshFieldContainer>(m, "MeshFieldContainer")
        .def(py::init<>())
        .def("set_field1D", &MeshFieldContainer::set_field1D);
    py::class_<MeshFieldBase>(m, "MeshFieldBase");
    py::class_<MeshField1D>(m, "MeshField1D")
        .def(py::init<const std::string&, const std::string&, const MeshBase1D&, int>())
        .def("get_MeshFieldBase", &MeshField1D::get_MeshFieldBase)
        .def("get_LocalMeshField", &MeshField1D::get_LocalMeshField)
        .def("print_val", &MeshField1D::print_val);
    py::class_<MeshFieldComm1D>(m, "MeshFieldComm1D")
        .def(py::init<int, int, const MeshLineDom1D&>())
        .def("put", &MeshFieldComm1D::put_py)
        .def("exchange", &MeshFieldComm1D::exchange)
        .def("get", &MeshFieldComm1D::get_py);
}
