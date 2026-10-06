#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/map.h>
#include "goddard/thermo.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_thermo(nb::module_& m) {
    nb::class_<Goddard::PhaseSpecification>(m, "PhaseSpecification",
            DOC(Goddard, PhaseSpecification))
        .def(nb::init<>(), "Create with default values.")
        .def(nb::init<double,
                        double,
                        const std::string&>(),
            "T"_a, "P"_a, "composition"_a = "",
            DOC(Goddard, PhaseSpecification, PhaseSpecification, 2))
        .def(nb::init<double,
                        double,
                        const std::map<std::string,double>&>(),
            "T"_a, "P"_a, "composition"_a,
            DOC(Goddard, PhaseSpecification, PhaseSpecification, 3))
        .def_rw("T", &Goddard::PhaseSpecification::T, DOC(Goddard, PhaseSpecification, T))
        .def_rw("P", &Goddard::PhaseSpecification::P, DOC(Goddard, PhaseSpecification, P))
        .def_rw("composition", &Goddard::PhaseSpecification::composition,
            DOC(Goddard, PhaseSpecification, composition));
}
