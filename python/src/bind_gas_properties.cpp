#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/unordered_set.h>
#include "goddard/gas.hpp"
#include "goddard/speciate.hpp"
#include "cantera/core.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_gas_properties(nb::module_& m) {
    nb::class_<Goddard::Gas>(m, "Gas",
        "Chemistry-aware wrapper around a Cantera Solution.\n\n"
        "Provides thermodynamic property queries, state setting, and derived\n"
        "flow calculations (speed of sound, stagnation properties, Mach number).\n"
        "The chemistry mode controls how properties like gamma_s and\n"
        "speed_of_sound are computed.")

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
             "Create a Gas from a YAML thermodynamic data file.")

        // ---- Read-only thermodynamic properties ----

        .def_prop_ro("temperature", &Goddard::Gas::temperature)
        .def_prop_ro("pressure", &Goddard::Gas::pressure)
        .def_prop_ro("density", &Goddard::Gas::density)
        .def_prop_ro("enthalpy_mass", &Goddard::Gas::enthalpy_mass)
        .def_prop_ro("entropy_mass", &Goddard::Gas::entropy_mass)
        .def_prop_ro("cp_mass", &Goddard::Gas::cp_mass)
        .def_prop_ro("cv_mass", &Goddard::Gas::cv_mass)
        .def_prop_ro("mean_molecular_weight", &Goddard::Gas::molecular_weight)
        .def_prop_ro("gamma_s", &Goddard::Gas::gamma_s)
        .def_prop_ro("speed_of_sound", &Goddard::Gas::speed_of_sound)

        // ---- State setters ----

        .def("set_state_TP", &Goddard::Gas::set_state_TP,
             "T"_a, "P"_a,
             "Set state from temperature (K) and pressure (Pa).")
        .def("set_state_TPX", nb::overload_cast<double, double, const std::string&>(&Goddard::Gas::set_state_TPX),
             "T"_a, "P"_a, "composition"_a,
             "Set state from temperature (K), pressure (Pa), and composition string.")
        .def("set_state_HP", &Goddard::Gas::set_state_HP,
             "H"_a, "P"_a,
             "Set state from specific enthalpy (J/kg) and pressure (Pa).")
        .def("set_state_SP", &Goddard::Gas::set_state_SP,
             "S"_a, "P"_a,
             "Set state from specific entropy (J/kg/K) and pressure (Pa).")

        // ---- State save/restore ----

        .def("save_state", &Goddard::Gas::save_state,
             "Save the current thermodynamic state as a vector.")
        .def("restore_state", &Goddard::Gas::restore_state,
             "state"_a,
             "Restore a previously saved thermodynamic state.")

        // ---- Derived flow calculations ----

        .def("stagnation_enthalpy", &Goddard::Gas::stagnation_enthalpy,
             "velocity"_a,
             "Compute stagnation enthalpy H0 = h + v^2/2.")
        .def("stagnation_pressure", &Goddard::Gas::stagnation_pressure,
             "velocity"_a,
             "Compute stagnation pressure via isentropic deceleration.")
        .def("isenthalpic_velocity",
             nb::overload_cast<double>(&Goddard::Gas::isenthalpic_velocity, nb::const_),
             "H_stagnation"_a,
             "Compute velocity from explicit stagnation enthalpy.")
        .def("isenthalpic_velocity",
             nb::overload_cast<>(&Goddard::Gas::isenthalpic_velocity, nb::const_),
             "Compute velocity from stored stagnation enthalpy.")
        .def("mach", &Goddard::Gas::mach,
             "velocity"_a,
             "Compute Mach number for a given velocity.")

        // ---- Expansion and equilibrium ----

        .def("expansion_properties", &Goddard::Gas::expansion_properties,
             "Compute expansion properties (gamma_s, dlV/dlT_P, dlV/dlP_T, cp).")
        .def("equilibrate", &Goddard::Gas::equilibrate,
             "XY"_a, "solver"_a="gibbs",
             "Equilibrate the gas mixture. XY is e.g. 'HP', 'TP', 'SP'.")

        // ---- Snapshot ----

        .def("snapshot", &Goddard::Gas::snapshot,
             "Return a ThermoStateInfo with all current properties.")

        // ---- Reference state (read-write properties) ----

        .def_prop_rw("stagnation_enthalpy_ref",
            &Goddard::Gas::get_stagnation_enthalpy,
            &Goddard::Gas::set_stagnation_enthalpy,
            "Stored reference stagnation enthalpy (J/kg).")
        .def_prop_rw("reference_entropy",
            &Goddard::Gas::get_reference_entropy,
            &Goddard::Gas::set_reference_entropy,
            "Stored reference entropy (J/kg/K).")

        // ---- Access to underlying Cantera objects ----

        .def_prop_ro("solution", &Goddard::Gas::solution,
            "Underlying SolutionHandle for passing to solver constructors.")

        // ---- Chemistry mode ----

        .def_rw("chemistry", &Goddard::Gas::chemistry,
            "GasChemistry mode (PERFECT_GAS, FROZEN, EQUILIBRIUM, KINETIC).");
}
