#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/eigen/dense.h>
#include "goddard/thermoarray.hpp"

namespace nb = nanobind;

void bind_thermoarray(nb::module_& m) {
    nb::class_<Goddard::ThermoArray>(m, "ThermoArray")
        .def("size", &Goddard::ThermoArray::size)
        .def("ndim", &Goddard::ThermoArray::ndim)
        .def("shape", &Goddard::ThermoArray::shape)
        .def("is_shape_set", &Goddard::ThermoArray::is_shape_set)
        .def("get_state", &Goddard::ThermoArray::get_state, nb::arg("loc"))
        .def("temperature", &Goddard::ThermoArray::temperature, nb::arg("slice") = 0)
        .def("pressure", &Goddard::ThermoArray::pressure, nb::arg("slice") = 0)
        .def("enthalpy_mass", &Goddard::ThermoArray::enthalpy_mass, nb::arg("slice") = 0)
        .def("enthalpy_mole", &Goddard::ThermoArray::enthalpy_mole, nb::arg("slice") = 0)
        .def("entropy_mass", &Goddard::ThermoArray::entropy_mass, nb::arg("slice") = 0)
        .def("entropy_mole", &Goddard::ThermoArray::entropy_mole, nb::arg("slice") = 0)
        .def("internal_energy_mass", &Goddard::ThermoArray::internal_energy_mass, nb::arg("slice") = 0)
        .def("internal_energy_mole", &Goddard::ThermoArray::internal_energy_mole, nb::arg("slice") = 0)
        .def("mean_molecular_weight", &Goddard::ThermoArray::mean_molecular_weight, nb::arg("slice") = 0);
}
