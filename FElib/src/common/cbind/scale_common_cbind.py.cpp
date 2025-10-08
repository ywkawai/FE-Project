#include "scale_scalelib_cbind.hpp"
#include "scale_polynomial_cbind.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void bind_common(py::module_ &m){
  m.def("SCALElib_Init", &SCALElib_Init);
  m.def("SCALElib_Final", &SCALElib_Final);
  m.def("Polynomial_GenLagrangePoly", &Polynomial_GenLagrangePoly);
  m.def("Polynomial_GenDLagrangePoly_lglpt", &Polynomial_GenDLagrangePoly_lglpt);
  m.def("Polynomial_GenLegendrePoly", &Polynomial_GenLegendrePoly);
}