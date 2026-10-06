#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include "goddard/global.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_enums(nb::module_& m);
void bind_errors(nb::module_& m);
void bind_structs(nb::module_& m);
void bind_equilibrium(nb::module_& m);
void bind_thermoarray(nb::module_& m);
void bind_combustor(nb::module_& m);
void bind_nozzle(nb::module_& m);
void bind_problem(nb::module_& m);
void bind_thermo(nb::module_& m);
void bind_gas_properties(nb::module_& m);
void bind_moc(nb::module_& m);
void bind_kinetic_nozzle(nb::module_& m);
void bind_shocks(nb::module_& m);
void bind_detonations(nb::module_& m);

NB_MODULE(_core, m) {
    m.doc() = "Goddard rocket engine simulation toolkit";

    Goddard::setup_defaults();

    // Binding-only wrapper: global.hpp is not scanned for docstrings.
    m.def("add_data_directory", &Goddard::add_directory, "path"_a,
          "Add a directory to the search path for data files given by bare name.\n\n"
          "A file name without a directory, such as ``\"nasa9_gas.yaml\"``, is looked up in "
          "the current working directory, the directories added with this function, and "
          "Cantera's data directories. ``goddard.data_dir`` is added at import.\n\n"
          "Args:\n"
          "    path: Directory to add.");

    // Order matters: types must be registered before they are referenced.
    // SolutionHandle (bind_problem) before Gas (bind_gas_properties),
    // Gas before solvers that accept it (bind_nozzle, bind_moc, etc.).
    bind_enums(m);
    bind_errors(m);
    bind_thermo(m);
    bind_structs(m);
    bind_equilibrium(m);
    bind_thermoarray(m);
    bind_combustor(m);
    bind_problem(m);
    bind_gas_properties(m);
    bind_nozzle(m);
    bind_moc(m);
    bind_kinetic_nozzle(m);
    bind_shocks(m);
    // after bind_shocks: DetonationResult holds a ShockResult
    bind_detonations(m);
}
