#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include "goddard/equilibrium.hpp"

namespace nb = nanobind;

void bind_equilibrium(nb::module_& m) {
    // ExpansionProperties
    nb::class_<Goddard::ExpansionProperties>(m, "ExpansionProperties")
        .def(nb::init<>())
        .def_ro("dlogV_dlogT_P", &Goddard::ExpansionProperties::dlogV_dlogT_P)
        .def_ro("dlogV_dlogP_T", &Goddard::ExpansionProperties::dlogV_dlogP_T)
        .def_ro("spec_heat_p", &Goddard::ExpansionProperties::spec_heat_p)
        .def_ro("gamma_s", &Goddard::ExpansionProperties::gamma_s);

    // EquilibriumDerivatives
    nb::class_<Goddard::EquilibriumDerivatives>(m, "EquilibriumDerivatives")
        .def(nb::init<>())
        .def_ro("dpi_dlogT_P", &Goddard::EquilibriumDerivatives::dpi_dlogT_P)
        .def_ro("dlogn_dlogT_P", &Goddard::EquilibriumDerivatives::dlogn_dlogT_P)
        .def_ro("dpi_dlogP_T", &Goddard::EquilibriumDerivatives::dpi_dlogP_T)
        .def_ro("dlogn_dlogP_T", &Goddard::EquilibriumDerivatives::dlogn_dlogP_T);
}
