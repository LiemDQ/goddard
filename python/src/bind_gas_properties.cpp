#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/unordered_set.h>
#include <nanobind/eigen/dense.h>
#include <algorithm>
#include "goddard/gas.hpp"
#include "goddard/speciate.hpp"
#include "goddard_docstrings.h"
#include "cantera/core.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_gas_properties(nb::module_& m) {
    nb::class_<Goddard::Gas>(m, "Gas", DOC(Goddard, Gas))

        // ---- Constructors ----

        // YAML-based constructor. An empty `species` set takes the phase as the file
        // defines it; a non-empty one restricts the phase to those species via
        // Gas::create_from_species. `species` follows `chemistry` so that existing
        // positional calls of the form Gas(file, phase, chemistry) keep working.
        // `condensed_file` attaches candidate condensed species from a second data file:
        // `all_condensed` offers every compatible species of that file, otherwise the names
        // of `condensed_species` are offered. An empty `condensed_file` leaves the Gas
        // gas-only, whatever the other two arguments say.
        .def("__init__",
             [](Goddard::Gas* self,
                const std::string& yaml_file,
                const std::string& phase_name,
                Goddard::GasChemistry chemistry,
                const std::unordered_set<std::string>& species,
                const std::string& condensed_file,
                const std::unordered_set<std::string>& condensed_species,
                bool all_condensed) {
                 Goddard::Gas gas = species.empty()
                     ? Goddard::Gas(yaml_file, phase_name, chemistry)
                     : Goddard::Gas::create_from_species(
                           yaml_file, phase_name,
                           std::vector<std::string>(species.begin(), species.end()), chemistry);

                 if (!condensed_file.empty()) {
                     if (all_condensed) {
                         gas.add_all_condensed_species(condensed_file);
                     } else {
                         std::vector<std::string> names(condensed_species.begin(),
                                                        condensed_species.end());
                         std::sort(names.begin(), names.end());
                         gas.add_condensed_species(condensed_file, names);
                     }
                 }
                 new (self) Goddard::Gas(std::move(gas));
             },
             "yaml_file"_a, "phase_name"_a, "chemistry"_a = Goddard::GasChemistry::FROZEN,
             "species"_a = std::unordered_set<std::string>{},
             "condensed_file"_a = "",
             "condensed_species"_a = std::unordered_set<std::string>{},
             "all_condensed"_a = false,
             DOC(Goddard, Gas, Gas, 4))
         .def("clone", &Goddard::Gas::clone, DOC(Goddard, Gas, clone))
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
        .def_prop_ro("num_species", &Goddard::Gas::num_species)
        .def_prop_ro("species_names", &Goddard::Gas::species_names)
        .def_prop_ro("mole_fractions", &Goddard::Gas::mole_fractions, DOC(Goddard, Gas, mole_fractions))
        .def_prop_ro("mass_fractions", &Goddard::Gas::mass_fractions, DOC(Goddard, Gas, mass_fractions))
        .def_prop_ro("element_names", &Goddard::Gas::element_names, DOC(Goddard, Gas, element_names))
        .def_prop_ro("element_moles", &Goddard::Gas::element_moles, DOC(Goddard, Gas, element_moles))

        // ---- Condensed species ----

        .def_prop_ro("has_condensed_candidates", &Goddard::Gas::has_condensed_candidates,
             DOC(Goddard, Gas, has_condensed_candidates))
        .def_prop_ro("has_condensed_phases", &Goddard::Gas::has_condensed_phases,
             DOC(Goddard, Gas, has_condensed_phases))
        .def_prop_ro("condensed_species_names", &Goddard::Gas::condensed_species_names,
             DOC(Goddard, Gas, condensed_species_names))
        .def_prop_ro("condensed_moles", &Goddard::Gas::condensed_moles,
             DOC(Goddard, Gas, condensed_moles))
        .def_prop_ro("gas_mass_fraction", &Goddard::Gas::gas_mass_fraction,
             DOC(Goddard, Gas, gas_mass_fraction))
        .def_prop_ro("mixture_molecular_weight", &Goddard::Gas::mixture_molecular_weight,
             DOC(Goddard, Gas, mixture_molecular_weight))
        .def_prop_ro("mixture_mass_fractions", &Goddard::Gas::mixture_mass_fractions,
             DOC(Goddard, Gas, mixture_mass_fractions))
        .def_prop_ro("at_phase_transition", &Goddard::Gas::at_phase_transition,
             DOC(Goddard, Gas, at_phase_transition))
        .def_prop_ro("pinned_polymorphs", &Goddard::Gas::pinned_polymorphs,
             DOC(Goddard, Gas, pinned_polymorphs))
        .def_prop_ro("condensed_enthalpy_RT", &Goddard::Gas::condensed_enthalpy_RT,
             DOC(Goddard, Gas, condensed_enthalpy_RT))
        .def_prop_ro("condensed_cp_R", &Goddard::Gas::condensed_cp_R,
             DOC(Goddard, Gas, condensed_cp_R))
        .def_prop_ro("condensed_molar_masses", &Goddard::Gas::condensed_molar_masses,
             DOC(Goddard, Gas, condensed_molar_masses))
        .def_prop_ro("condensed_stoich_coeffs", &Goddard::Gas::condensed_stoich_coeffs,
             DOC(Goddard, Gas, condensed_stoich_coeffs))
        .def("add_condensed_species", &Goddard::Gas::add_condensed_species,
             "infile"_a, "names"_a, DOC(Goddard, Gas, add_condensed_species))
        .def("add_all_condensed_species", &Goddard::Gas::add_all_condensed_species,
             "infile"_a, DOC(Goddard, Gas, add_all_condensed_species))
        .def("set_condensed_moles", &Goddard::Gas::set_condensed_moles,
             "moles"_a, DOC(Goddard, Gas, set_condensed_moles))
        .def("set_phase_transition",
             nb::overload_cast<const std::string&, const std::string&>(
                 &Goddard::Gas::set_phase_transition),
             "low"_a, "high"_a, DOC(Goddard, Gas, set_phase_transition, 2))
        .def("clear_phase_transition", &Goddard::Gas::clear_phase_transition,
             DOC(Goddard, Gas, clear_phase_transition))

        // ---- State setters ----

        .def("set_state_TP", &Goddard::Gas::set_state_TP,
             "T"_a, "P"_a, DOC(Goddard, Gas, set_state_TP))
        .def("set_state_TPX", nb::overload_cast<double, double, const std::string&>(&Goddard::Gas::set_state_TPX),
             "T"_a, "P"_a, "composition"_a, DOC(Goddard, Gas, set_state_TPX))
        .def("set_state_TPX", nb::overload_cast<double, double, const Goddard::Composition&>(&Goddard::Gas::set_state_TPX),
             "T"_a, "P"_a, "composition"_a, DOC(Goddard, Gas, set_state_TPX, 2))
        .def("set_state_TPY", nb::overload_cast<double, double, const std::string&>(&Goddard::Gas::set_state_TPY),
             "T"_a, "P"_a, "composition"_a, DOC(Goddard, Gas, set_state_TPY))
        .def("set_state_TPY", nb::overload_cast<double, double, const Goddard::Composition&>(&Goddard::Gas::set_state_TPY),
             "T"_a, "P"_a, "composition"_a, DOC(Goddard, Gas, set_state_TPY, 2))
        .def("set_element_moles", &Goddard::Gas::set_element_moles,
             "element_moles"_a, "T"_a, "P"_a, DOC(Goddard, Gas, set_element_moles))
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
        .def("equilibrate_TP", &Goddard::Gas::equilibrate_TP,
             "T"_a, "P"_a, DOC(Goddard, Gas, equilibrate_TP))
        .def("equilibrate_HP", &Goddard::Gas::equilibrate_HP,
             "H"_a, "P"_a, DOC(Goddard, Gas, equilibrate_HP))
        .def("equilibrate_SP", &Goddard::Gas::equilibrate_SP,
             "S"_a, "P"_a, DOC(Goddard, Gas, equilibrate_SP))
        .def("last_equilibrium_solve_count", &Goddard::Gas::last_equilibrium_solve_count,
             DOC(Goddard, Gas, last_equilibrium_solve_count))
        .def_rw("equilibrium_options", &Goddard::Gas::equilibrium_options,
             DOC(Goddard, Gas, equilibrium_options))

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
