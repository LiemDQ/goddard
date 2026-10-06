#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/eigen/dense.h>
#include "goddard/thermoarray.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;

void bind_thermoarray(nb::module_& m) {
    nb::class_<Goddard::ThermoArray>(m, "ThermoArray", DOC(Goddard, ThermoArray))
        .def_prop_ro("size", &Goddard::ThermoArray::size, DOC(Goddard, ThermoArray, size))
        .def_prop_ro("ndim", &Goddard::ThermoArray::ndim, DOC(Goddard, ThermoArray, ndim))
        .def_prop_ro("shape", &Goddard::ThermoArray::shape, DOC(Goddard, ThermoArray, shape))
        .def_prop_ro("is_shape_set", &Goddard::ThermoArray::is_shape_set,
             DOC(Goddard, ThermoArray, is_shape_set))
        .def("flat_index", &Goddard::ThermoArray::flat_index,
             nb::arg("i"), nb::arg("j") = 0, nb::arg("k") = 0,
             DOC(Goddard, ThermoArray, flat_index))
        .def("get_state", &Goddard::ThermoArray::get_state, nb::arg("loc"),
             DOC(Goddard, ThermoArray, get_state))
        .def("num_condensed", &Goddard::ThermoArray::num_condensed,
             DOC(Goddard, ThermoArray, num_condensed))
        .def("condensed_species_names", &Goddard::ThermoArray::condensed_species_names,
             DOC(Goddard, ThermoArray, condensed_species_names))
        .def("get_condensed_moles", &Goddard::ThermoArray::get_condensed_moles, nb::arg("loc"),
             DOC(Goddard, ThermoArray, get_condensed_moles))
        .def("temperature", &Goddard::ThermoArray::temperature, nb::arg("slice") = 0,
             DOC(Goddard, ThermoArray, temperature))
        .def("pressure", &Goddard::ThermoArray::pressure, nb::arg("slice") = 0,
             DOC(Goddard, ThermoArray, pressure))
        .def("enthalpy_mass", &Goddard::ThermoArray::enthalpy_mass, nb::arg("slice") = 0,
             DOC(Goddard, ThermoArray, enthalpy_mass))
        .def("enthalpy_mole", &Goddard::ThermoArray::enthalpy_mole, nb::arg("slice") = 0,
             DOC(Goddard, ThermoArray, enthalpy_mole))
        .def("entropy_mass", &Goddard::ThermoArray::entropy_mass, nb::arg("slice") = 0,
             DOC(Goddard, ThermoArray, entropy_mass))
        .def("entropy_mole", &Goddard::ThermoArray::entropy_mole, nb::arg("slice") = 0,
             DOC(Goddard, ThermoArray, entropy_mole))
        .def("internal_energy_mass", &Goddard::ThermoArray::internal_energy_mass, nb::arg("slice") = 0,
             DOC(Goddard, ThermoArray, internal_energy_mass))
        .def("internal_energy_mole", &Goddard::ThermoArray::internal_energy_mole, nb::arg("slice") = 0,
             DOC(Goddard, ThermoArray, internal_energy_mole))
        .def("mean_molecular_weight", &Goddard::ThermoArray::mean_molecular_weight, nb::arg("slice") = 0,
             DOC(Goddard, ThermoArray, mean_molecular_weight));
}
