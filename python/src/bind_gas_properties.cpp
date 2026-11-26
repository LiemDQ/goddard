#include <nanobind/nanobind.h>
#include <format>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/unordered_set.h>
#include <nanobind/eigen/dense.h>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include "goddard/gas.hpp"
#include "goddard/speciate.hpp"
#include "goddard_docstrings.h"
#include "cantera/core.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

// Composition basis as a string, so the Python API does not expose Cantera's ThermoBasis.
Cantera::ThermoBasis to_basis(const std::string& basis) {
    if (basis == "mole") return Cantera::ThermoBasis::molar;
    if (basis == "mass") return Cantera::ThermoBasis::mass;
    throw std::invalid_argument("basis must be 'mole' or 'mass', got '" + basis + "'");
}

// Binds a mixture method for fuel and oxidizer given either as composition strings or as dicts.
template <typename Class, typename Method>
void def_mixture(Class& cls, const char* name, Method method, const char* doc) {
    cls.def(name, [method](const Goddard::Gas& gas, const std::string& fuel, const std::string& oxidizer,
                           const std::string& basis) {
            return std::invoke(method, gas, fuel, oxidizer, to_basis(basis));
        }, "fuel"_a, "oxidizer"_a, "basis"_a = "mole", doc);
}

} // namespace

void bind_gas_properties(nb::module_& m) {
    auto gas_class = nb::class_<Goddard::Gas>(m, "Gas", DOC(Goddard, Gas));
    gas_class

        // ---- Constructors ----

        // YAML-based constructor. An empty `species` set takes the phase as the file
        // defines it; a non-empty one restricts the phase to those species via
        // Gas::create_from_species, which orders them as in the file, so the iteration order
        // of the set does not matter.
        // `condensed_file` attaches candidate condensed species from a second data file:
        // `all_condensed` offers every compatible species of that file, otherwise the names
        // of `condensed_species` are offered. An empty `condensed_file` leaves the Gas
        // gas-only, whatever the other two arguments say.
        .def("__init__",
             [](Goddard::Gas* self,
                const std::string& yaml_file,
                const std::string& phase_name,
                const std::unordered_set<std::string>& species,
                const std::string& condensed_file,
                const std::unordered_set<std::string>& condensed_species,
                bool all_condensed) {
                 Goddard::Gas gas = species.empty()
                     ? Goddard::Gas(yaml_file, phase_name)
                     : Goddard::Gas::create_from_species(
                           yaml_file, phase_name,
                           std::vector<std::string>(species.begin(), species.end()));

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
             "yaml_file"_a, "phase_name"_a,
             "species"_a = std::unordered_set<std::string>{},
             "condensed_file"_a = "",
             "condensed_species"_a = std::unordered_set<std::string>{},
             "all_condensed"_a = false,
             DOC(Goddard, Gas, Gas, 5))
        // Element-based product selection: every species of the file made of these elements.
        // The condensed arguments work as in the constructor.
        .def_static("from_elements",
             [](const std::string& yaml_file, const std::string& phase_name,
                const std::vector<std::string>& elements, const std::string& condensed_file,
                const std::unordered_set<std::string>& condensed_species, bool all_condensed) {
                 Goddard::Gas gas = Goddard::Gas::create_from_elements(yaml_file, phase_name, elements);
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
                 return gas;
             },
             "yaml_file"_a, "phase_name"_a, "elements"_a, "condensed_file"_a = "",
             "condensed_species"_a = std::unordered_set<std::string>{}, "all_condensed"_a = false,
             DOC(Goddard, Gas, create_from_elements))
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
        // Properties take the frozen composition model; the equilibrium_* properties and the
        // `model` argument of the methods below give the equilibrium values.
        .def_prop_ro("gamma_s",
             [](const Goddard::Gas& gas) { return gas.gamma_s(Goddard::CompositionModel::FROZEN); },
             "Frozen isentropic exponent cp/cv [-] at the current state.")
        .def_prop_ro("equilibrium_gamma_s",
             [](const Goddard::Gas& gas) { return gas.gamma_s(Goddard::CompositionModel::EQUILIBRIUM); },
             "Equilibrium isentropic exponent -(d log P / d log V)_s [-] at the current state, with "
             "the composition shifting to remain at chemical equilibrium.")
        .def_prop_ro("speed_of_sound",
             [](const Goddard::Gas& gas) { return gas.speed_of_sound(Goddard::CompositionModel::FROZEN); },
             "Frozen speed of sound [m/s] at the current state, consistent with `gamma_s`.")
        .def_prop_ro("equilibrium_speed_of_sound",
             [](const Goddard::Gas& gas) { return gas.speed_of_sound(Goddard::CompositionModel::EQUILIBRIUM); },
             "Equilibrium speed of sound [m/s] at the current state, consistent with "
             "`equilibrium_gamma_s`.")
        .def_prop_ro("num_species", &Goddard::Gas::num_species, DOC(Goddard, Gas, num_species))
        .def_prop_ro("species_names", &Goddard::Gas::species_names,
             DOC(Goddard, Gas, species_names))
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
        .def("set_state_TD", &Goddard::Gas::set_state_TD,
             "T"_a, "D"_a, DOC(Goddard, Gas, set_state_TD))
        .def("set_state_UV", &Goddard::Gas::set_state_UV,
             "U"_a, "V"_a, DOC(Goddard, Gas, set_state_UV))

        // ---- State save/restore ----

        .def("save_state", &Goddard::Gas::save_state,
             DOC(Goddard, Gas, save_state))
        .def("restore_state", &Goddard::Gas::restore_state,
             "state"_a, DOC(Goddard, Gas, restore_state))

        // ---- Derived flow calculations ----

        .def("stagnation_enthalpy", &Goddard::Gas::stagnation_enthalpy,
             "velocity"_a, DOC(Goddard, Gas, stagnation_enthalpy))
        .def("stagnation_pressure", &Goddard::Gas::stagnation_pressure,
             "velocity"_a, "model"_a = Goddard::CompositionModel::FROZEN,
             DOC(Goddard, Gas, stagnation_pressure))
        .def("isenthalpic_velocity", &Goddard::Gas::isenthalpic_velocity,
             "H_stagnation"_a, DOC(Goddard, Gas, isenthalpic_velocity))
        .def("mach", &Goddard::Gas::mach,
             "velocity"_a, "model"_a = Goddard::CompositionModel::FROZEN, DOC(Goddard, Gas, mach))
        .def("area_per_mdot", &Goddard::Gas::area_per_mdot,
             "velocity"_a, DOC(Goddard, Gas, area_per_mdot))
        .def("cstar", &Goddard::Gas::cstar, "model"_a = Goddard::CompositionModel::FROZEN,
             DOC(Goddard, Gas, cstar))

        // ---- Mixture ratios ----

        .def("equivalence_ratio",
             nb::overload_cast<>(&Goddard::Gas::equivalence_ratio, nb::const_),
             DOC(Goddard, Gas, equivalence_ratio))
        .def("set_fuel_fraction",
             [](Goddard::Gas& gas, double fuel_fraction, const std::string& fuel,
                const std::string& oxidizer, const std::string& basis) {
                 gas.set_fuel_fraction(fuel_fraction, fuel, oxidizer, to_basis(basis));
             }, "fuel_fraction"_a, "fuel"_a, "oxidizer"_a, "basis"_a = "mole",
             DOC(Goddard, Gas, set_fuel_fraction))
        .def("set_equivalence_ratio",
             [](Goddard::Gas& gas, double phi, const std::string& fuel,
                const std::string& oxidizer, const std::string& basis) {
                 gas.set_equivalence_ratio(phi, fuel, oxidizer, to_basis(basis));
             }, "phi"_a, "fuel"_a, "oxidizer"_a, "basis"_a = "mole",
             DOC(Goddard, Gas, set_equivalence_ratio))
        .def("set_OF_ratio",
             [](Goddard::Gas& gas, double OF, const std::string& fuel,
                const std::string& oxidizer, const std::string& basis) {
                 gas.set_OF_ratio(OF, fuel, oxidizer, to_basis(basis));
             }, "OF"_a, "fuel"_a, "oxidizer"_a, "basis"_a = "mole",
             DOC(Goddard, Gas, set_OF_ratio))

        .def("report", &Goddard::Gas::report, "show_thermo"_a = true, "threshold"_a = -1e-14,
             DOC(Goddard, Gas, report))

        // ---- Expansion and equilibrium ----

        .def("expansion_properties", &Goddard::Gas::expansion_properties,
             "model"_a = Goddard::CompositionModel::FROZEN, DOC(Goddard, Gas, expansion_properties))
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
        .def("check_temperature_range", &Goddard::Gas::check_temperature_range,
             DOC(Goddard, Gas, check_temperature_range))
        .def_rw("equilibrium_options", &Goddard::Gas::equilibrium_options,
             DOC(Goddard, Gas, equilibrium_options))

        // ---- Snapshot ----

        .def("snapshot", &Goddard::Gas::snapshot, "model"_a = Goddard::CompositionModel::FROZEN,
             DOC(Goddard, Gas, snapshot));

    gas_class.def("__repr__", [](const Goddard::Gas& self) {
        return std::format("<Gas '{}' T={:.6g} K P={:.6g} Pa species={} condensed={}>",
                           self.name(), self.temperature(), self.pressure(), self.num_species(),
                           self.condensed_species_names().size());
    });
    def_mixture(gas_class, "fuel_fraction",
        nb::overload_cast<const std::string&, const std::string&, Cantera::ThermoBasis>(
            &Goddard::Gas::fuel_fraction, nb::const_),
        DOC(Goddard, Gas, fuel_fraction));
    def_mixture(gas_class, "equivalence_ratio",
        nb::overload_cast<const std::string&, const std::string&, Cantera::ThermoBasis>(
            &Goddard::Gas::equivalence_ratio, nb::const_),
        DOC(Goddard, Gas, equivalence_ratio, 2));
    def_mixture(gas_class, "stoich_OF_ratio",
        nb::overload_cast<const std::string&, const std::string&, Cantera::ThermoBasis>(
            &Goddard::Gas::stoich_OF_ratio, nb::const_),
        DOC(Goddard, Gas, stoich_OF_ratio));
}
