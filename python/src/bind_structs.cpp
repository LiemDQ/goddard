#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/unordered_set.h>
#include "goddard/combustor.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/problem.hpp"

namespace nb = nanobind;

void bind_structs(nb::module_& m) {
    // CombustorOptions
    nb::class_<Goddard::CombustorOptions>(m, "CombustorOptions")
        .def(nb::init<>())
        .def_rw("type", &Goddard::CombustorOptions::type)
        .def_rw("pressures", &Goddard::CombustorOptions::pressures)
        .def_rw("mass_flux", &Goddard::CombustorOptions::mass_flux)
        .def_rw("contraction_ratio", &Goddard::CombustorOptions::contraction_ratio);

    // NozzleOptions
    nb::class_<Goddard::NozzleOptions>(m, "NozzleOptions")
        .def(nb::init<>())
        .def_rw("chemistry", &Goddard::NozzleOptions::chemistry)
        .def_rw("expansion_type", &Goddard::NozzleOptions::expansion_type)
        .def_rw("expansion_ratios", &Goddard::NozzleOptions::expansion_ratios)
        .def_rw("frozen_NFZ", &Goddard::NozzleOptions::frozen_NFZ);

    // ThroatCondition
    nb::class_<Goddard::ThroatCondition>(m, "ThroatCondition")
        .def(nb::init<>())
        .def_ro("converged", &Goddard::ThroatCondition::converged)
        .def_ro("H_stagnation", &Goddard::ThroatCondition::H_stagnation)
        .def_ro("P_inlet", &Goddard::ThroatCondition::P_inlet)
        .def_ro("S_inlet", &Goddard::ThroatCondition::S_inlet)
        .def_ro("gamma_s", &Goddard::ThroatCondition::gamma_s)
        .def_ro("dlV_dlP_T", &Goddard::ThroatCondition::dlV_dlP_T)
        .def_ro("dlV_dlT_P", &Goddard::ThroatCondition::dlV_dlT_P)
        .def_ro("state", &Goddard::ThroatCondition::state);

    // NozzleResult
    nb::class_<Goddard::NozzleResult>(m, "NozzleResult")
        .def(nb::init<>())
        .def_ro("converged", &Goddard::NozzleResult::converged)
        .def_ro("gamma_s", &Goddard::NozzleResult::gamma_s)
        .def_ro("dlV_dlP_T", &Goddard::NozzleResult::dlV_dlP_T)
        .def_ro("dlV_dlT_P", &Goddard::NozzleResult::dlV_dlT_P)
        .def_ro("state", &Goddard::NozzleResult::state);

    // NozzleResults
    nb::class_<Goddard::NozzleResults>(m, "NozzleResults")
        .def(nb::init<>())
        .def_ro("throat", &Goddard::NozzleResults::throat)
        .def_ro("expansions", &Goddard::NozzleResults::expansions);

    // RocketCaseParameters
    nb::class_<Goddard::RocketCaseParameters>(m, "RocketCaseParameters")
        .def(nb::init<>())
        .def_rw("name", &Goddard::RocketCaseParameters::name)
        .def_rw("problem_type", &Goddard::RocketCaseParameters::problem_type)
        .def_rw("combustor_options", &Goddard::RocketCaseParameters::combustor_options)
        .def_rw("nozzle_options", &Goddard::RocketCaseParameters::nozzle_options);

    // ChemicalParameters
    nb::class_<Goddard::ChemicalParameters>(m, "ChemicalParameters")
        .def(nb::init<>())
        .def_rw("thermo_file", &Goddard::ChemicalParameters::thermo_file)
        .def_rw("species", &Goddard::ChemicalParameters::species)
        .def_rw("cantera_fuel_state", &Goddard::ChemicalParameters::cantera_fuel_state)
        .def_rw("cantera_oxidizer_state", &Goddard::ChemicalParameters::cantera_oxidizer_state)
        .def_rw("OF_ratios", &Goddard::ChemicalParameters::OF_ratios)
        .def_rw("phi_ratios", &Goddard::ChemicalParameters::phi_ratios)
        .def_rw("fuel_weight_percentages", &Goddard::ChemicalParameters::fuel_weight_percentages)
        .def_rw("valance_equivalences", &Goddard::ChemicalParameters::valance_equivalences);

    // ThermoStateInfo
    nb::class_<Goddard::ThermoStateInfo>(m, "ThermoStateInfo")
        .def(nb::init<>())
        .def_ro("pressure", &Goddard::ThermoStateInfo::pressure)
        .def_ro("temperature", &Goddard::ThermoStateInfo::temperature)
        .def_ro("density", &Goddard::ThermoStateInfo::density)
        .def_ro("enthalpy", &Goddard::ThermoStateInfo::enthalpy)
        .def_ro("internal_energy", &Goddard::ThermoStateInfo::internal_energy)
        .def_ro("gibbs", &Goddard::ThermoStateInfo::gibbs)
        .def_ro("entropy", &Goddard::ThermoStateInfo::entropy)
        .def_ro("molecular_weight", &Goddard::ThermoStateInfo::molecular_weight)
        .def_ro("gamma_s", &Goddard::ThermoStateInfo::gamma_s)
        .def_ro("dlV_dlP_T", &Goddard::ThermoStateInfo::dlV_dlP_T)
        .def_ro("dlV_dlT_P", &Goddard::ThermoStateInfo::dlV_dlT_P)
        .def_ro("speed_of_sound", &Goddard::ThermoStateInfo::speed_of_sound)
        .def_ro("composition", &Goddard::ThermoStateInfo::composition);

    // RocketPerformance
    nb::class_<Goddard::RocketPerformance>(m, "RocketPerformance")
        .def(nb::init<>())
        .def_ro("pressure_ratio", &Goddard::RocketPerformance::pressure_ratio)
        .def_ro("area_ratio", &Goddard::RocketPerformance::area_ratio)
        .def_ro("mach_number", &Goddard::RocketPerformance::mach_number)
        .def_ro("cstar", &Goddard::RocketPerformance::cstar)
        .def_ro("CF", &Goddard::RocketPerformance::CF)
        .def_ro("isp", &Goddard::RocketPerformance::isp)
        .def_ro("ivac", &Goddard::RocketPerformance::ivac);

    // RocketState
    nb::class_<Goddard::RocketState>(m, "RocketState")
        .def(nb::init<>())
        .def_ro("name", &Goddard::RocketState::name)
        .def_ro("cantera_state", &Goddard::RocketState::cantera_state)
        .def_ro("pressure_ratio", &Goddard::RocketState::pressure_ratio)
        .def_ro("area_ratio", &Goddard::RocketState::area_ratio)
        .def_ro("dlv_dlp_t", &Goddard::RocketState::dlv_dlp_t)
        .def_ro("dlv_dlt_p", &Goddard::RocketState::dlv_dlt_p)
        .def_ro("gamma_s", &Goddard::RocketState::gamma_s)
        .def_ro("speed_of_sound", &Goddard::RocketState::speed_of_sound);
}
