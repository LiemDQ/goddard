#include <nanobind/nanobind.h>
#include "goddard/global.hpp"

namespace nb = nanobind;

void bind_enums(nb::module_& m);
void bind_errors(nb::module_& m);
void bind_structs(nb::module_& m);
void bind_equilibrium(nb::module_& m);
void bind_mixture_ratio(nb::module_& m);
void bind_thermoarray(nb::module_& m);
void bind_combustor(nb::module_& m);
void bind_nozzle(nb::module_& m);
void bind_problem(nb::module_& m);

NB_MODULE(_core, m) {
    m.doc() = "Goddard rocket engine simulation toolkit";

    Goddard::setup_defaults();

    // Order matters: enums and structs first since classes reference them
    bind_enums(m);
    bind_errors(m);
    bind_structs(m);
    bind_equilibrium(m);
    bind_mixture_ratio(m);
    bind_thermoarray(m);
    bind_combustor(m);
    bind_nozzle(m);
    bind_problem(m);
}
