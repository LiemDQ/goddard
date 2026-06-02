#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/unordered_set.h>
#include "goddard/gas.hpp"
#include "goddard/speciate.hpp"
#include "goddard_docstrings.h"
#include "cantera/core.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_gas_properties(nb::module_& m) {
    nb::class_<Goddard::Gas>(m, "Gas", DOC(Goddard, Gas))

        // ---- Constructors ----

        // YAML-based constructor
        .def("__init__",
             [](Goddard::Gas* self,
                const std::string& yaml_file,
                const std::string& phase_name,
                Goddard::GasChemistry chemistry) {
                 new (self) Goddard::Gas(yaml_file, phase_name, chemistry);
             },
             "yaml_file"_a, "phase_name"_a, "chemistry"_a = Goddard::GasChemistry::FROZEN,
             DOC(Goddard, Gas, Gas, 4))

        // ---- Read-only thermodynamic properties ----

        .def_prop_ro("temperature", &Goddard::Gas::temperature, DOC(Goddard, Gas, temperature))
        .def_prop_ro("pressure", &Goddard::Gas::pressure, DOC(Goddard, Gas, pressure))
        .def_prop_ro("density", &Goddard::Gas::density, DOC(Goddard, Gas, density))
        .def_prop_ro("enthalpy_mass", &Goddard::Gas::enthalpy_mass, DOC(Goddard, Gas, enthalpy_mass))
        .def_prop_ro("entropy_mass", &Goddard::Gas::entropy_mass, DOC(Goddard, Gas, entropy_mass))
        .def_prop_ro("cp_mass", &Goddard::Gas::cp_mass, DOC(Goddard, Gas, cp_mass))
        .def_prop_ro("cv_mass", &Goddard::Gas::cv_mass, DOC(Goddard, Gas, cv_mass))
        .def_prop_ro("mean_molecular_weight", &Goddard::Gas::molecular_weight, DOC(Goddard, Gas, molecular_weight))
        .def_prop_ro("gamma_s", &Goddard::Gas::gamma_s, DOC(Goddard, Gas, gamma_s))
        .def_prop_ro("speed_of_sound", &Goddard::Gas::speed_of_sound, DOC(Goddard, Gas, speed_of_sound))

        // ---- State setters ----

        .def("set_state_TP", &Goddard::Gas::set_state_TP,
             "T"_a, "P"_a, DOC(Goddard, Gas, set_state_TP))
        .def("set_state_TPX", nb::overload_cast<double, double, const std::string&>(&Goddard::Gas::set_state_TPX),
             "T"_a, "P"_a, "composition"_a, DOC(Goddard, Gas, set_state_TPX))
        .def("set_state_HP", &Goddard::Gas::set_state_HP,
             "H"_a, "P"_a, DOC(Goddard, Gas, set_state_HP))
        .def("set_state_SP", &Goddard::Gas::set_state_SP,
             "S"_a, "P"_a, DOC(Goddard, Gas, set_state_SP))

        // ---- State save/restore ----

        .def("save_state", &Goddard::Gas::save_state,
             DOC(Goddard, Gas, save_state))
        .def("restore_state", &Goddard::Gas::restore_state,
             "state"_a, DOC(Goddard, Gas, restore_state))

        // ---- Derived flow calculations ----

        .def("stagnation_enthalpy", &Goddard::Gas::stagnation_enthalpy,
             "velocity"_a, DOC(Goddard, Gas, stagnation_enthalpy))
        .def("stagnation_pressure", &Goddard::Gas::stagnation_pressure,
             "velocity"_a, DOC(Goddard, Gas, stagnation_pressure))
        .def("isenthalpic_velocity",
             nb::overload_cast<double>(&Goddard::Gas::isenthalpic_velocity, nb::const_),
             "H_stagnation"_a, DOC(Goddard, Gas, isenthalpic_velocity))
        .def("isenthalpic_velocity",
             nb::overload_cast<>(&Goddard::Gas::isenthalpic_velocity, nb::const_),
             DOC(Goddard, Gas, isenthalpic_velocity, 2))
        .def("mach", &Goddard::Gas::mach,
             "velocity"_a, DOC(Goddard, Gas, mach))

        // ---- Expansion and equilibrium ----

        .def("expansion_properties", &Goddard::Gas::expansion_properties,
             DOC(Goddard, Gas, expansion_properties))
        .def("equilibrate", &Goddard::Gas::equilibrate,
             "XY"_a, "solver"_a="gibbs", DOC(Goddard, Gas, equilibrate))

        // ---- Snapshot ----

        .def("snapshot", &Goddard::Gas::snapshot,
             DOC(Goddard, Gas, snapshot))

        // ---- Reference state (read-write properties) ----

        .def_prop_rw("stagnation_enthalpy_ref",
            &Goddard::Gas::get_stagnation_enthalpy,
            &Goddard::Gas::set_stagnation_enthalpy,
            DOC(Goddard, Gas, get_stagnation_enthalpy))
        .def_prop_rw("reference_entropy",
            &Goddard::Gas::get_reference_entropy,
            &Goddard::Gas::set_reference_entropy,
            DOC(Goddard, Gas, get_reference_entropy))

        // ---- Access to underlying Cantera objects ----

        .def_prop_ro("solution", &Goddard::Gas::solution,
            DOC(Goddard, Gas, solution))

        // ---- Chemistry mode ----

        .def_rw("chemistry", &Goddard::Gas::chemistry,
            "GasChemistry mode (PERFECT_GAS, FROZEN, EQUILIBRIUM, KINETIC).");
}
