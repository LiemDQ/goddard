#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/map.h>
#include "goddard/thermo.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_thermo(nb::module_& m) {
    nb::class_<Goddard::PhaseSpecification>(m, "PhaseSpecification")
        .def(nb::init<>())
        .def(nb::init<double,
                        double,
                        const std::string&>(),
            "T"_a, "P"_a, "composition"_a = "")
        .def(nb::init<double,
                        double,
                        const std::map<std::string,double>&>(),
            "T"_a, "P"_a, "composition"_a)
        .def_rw("T", &Goddard::PhaseSpecification::T)
        .def_rw("P", &Goddard::PhaseSpecification::P)
        .def_rw("composition", &Goddard::PhaseSpecification::composition);
}
