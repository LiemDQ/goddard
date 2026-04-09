#include <nanobind/nanobind.h>
#include "goddard/thermo.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_thermo(nb::module_& m) {
    nb::class_<Goddard::InputState>(m, "InputState")
        .def(nb::init<>())
        .def(nb::init<double,
                        double,
                        const std::string&>(),
            "T"_a, "P"_a, "composition"_a = "")
        .def_rw("T", &Goddard::InputState::T)
        .def_rw("P", &Goddard::InputState::P)
        .def_rw("composition", &Goddard::InputState::composition);
}
