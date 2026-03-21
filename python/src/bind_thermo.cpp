#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include "goddard/thermo.hpp"


namespace nb = nanobind;
using namespace nb::literals;

void bind_thermo(nb::module_& m) {
    nb::class_<Goddard::ThermodynamicState>(m, "ThermodynamicState")
        .def(nb::init<>())
        .def(nb::init<double,
                        double,
                        const std::string&>(),
            "T"_a, "P"_a, "composition"_a = "")
        .def_rw("T", &Goddard::ThermodynamicState::T)
        .def_rw("P", &Goddard::ThermodynamicState::P)
        .def_rw("composition", &Goddard::ThermodynamicState::composition);
}
