#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/unordered_set.h>
#include <unordered_set>
#include <vector>

#include "goddard/combustor.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/problem.hpp"
#include "goddard/rocket_results.hpp"
#include "goddard/thermo.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_structs(nb::module_& m) {
    // CombustorOptions
    nb::class_<Goddard::CombustorOptions>(m, "CombustorOptions")
        .def(nb::init<>())
        .def("__init__", [](Goddard::CombustorOptions* self,
                            Goddard::CombustorType type,
                            std::vector<double> pressures,
                            double mass_flux,
                            double contraction_ratio) {
            new (self) Goddard::CombustorOptions();
            self->type = type;
            self->pressures = std::move(pressures);
            self->mass_flux = mass_flux;
            self->contraction_ratio = contraction_ratio;
        },  "type"_a = Goddard::CombustorType::INFINITE_AREA,
            "pressures"_a = std::vector<double>(),
            "mass_flux"_a = 0.0,
            "contraction_ratio"_a = 0.0)
        .def_rw("type", &Goddard::CombustorOptions::type)
        .def_rw("pressures", &Goddard::CombustorOptions::pressures)
        .def_rw("mass_flux", &Goddard::CombustorOptions::mass_flux)
        .def_rw("contraction_ratio", &Goddard::CombustorOptions::contraction_ratio);

    // NozzleOptions
    nb::class_<Goddard::NozzleOptions>(m, "NozzleOptions")
        .def(nb::init<>())
        .def("__init__", [](Goddard::NozzleOptions* self,
                            Goddard::GasChemistry chemistry,
                            Goddard::ExpansionType expansion_type,
                            std::vector<double> expansion_ratios,
                            unsigned int frozen_NFZ) {
            new (self) Goddard::NozzleOptions();
            self->chemistry = chemistry;
            self->expansion_type = expansion_type;
            self->expansion_ratios = std::move(expansion_ratios);
            self->frozen_NFZ = frozen_NFZ;
        },  "chemistry"_a = Goddard::GasChemistry::EQUILIBRIUM,
            "expansion_type"_a = Goddard::ExpansionType::SUPERSONIC_AREA_RATIO,
            "expansion_ratios"_a = std::vector<double>(),
            "frozen_NFZ"_a = 1u)
        .def_rw("chemistry", &Goddard::NozzleOptions::chemistry)
        .def_rw("expansion_type", &Goddard::NozzleOptions::expansion_type)
        .def_rw("expansion_ratios", &Goddard::NozzleOptions::expansion_ratios)
        .def_rw("frozen_NFZ", &Goddard::NozzleOptions::frozen_NFZ);

    // ThroatCondition
    nb::class_<Goddard::ThroatCondition>(m, "ThroatCondition")
        .def(nb::init<>())
        .def("__init__", [](Goddard::ThroatCondition* self,
                            bool converged,
                            double H_stagnation,
                            double P_inlet,
                            double S_inlet,
                            double gamma_s,
                            double dlV_dlP_T,
                            double dlV_dlT_P,
                            std::vector<double> state) {
            new (self) Goddard::ThroatCondition();
            self->converged = converged;
            self->H_stagnation = H_stagnation;
            self->P_inlet = P_inlet;
            self->S_inlet = S_inlet;
            self->gamma_s = gamma_s;
            self->dlV_dlP_T = dlV_dlP_T;
            self->dlV_dlT_P = dlV_dlT_P;
            self->state = std::move(state);
        },  "converged"_a = false,
            "H_stagnation"_a = 0.0,
            "P_inlet"_a = 0.0,
            "S_inlet"_a = 0.0,
            "gamma_s"_a = 0.0,
            "dlV_dlP_T"_a = 0.0,
            "dlV_dlT_P"_a = 0.0,
            "state"_a = std::vector<double>())
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
        .def("__init__", [](Goddard::NozzleResult* self,
                            bool converged,
                            double gamma_s,
                            double dlV_dlP_T,
                            double dlV_dlT_P,
                            std::vector<double> state) {
            new (self) Goddard::NozzleResult();
            self->converged = converged;
            self->gamma_s = gamma_s;
            self->dlV_dlP_T = dlV_dlP_T;
            self->dlV_dlT_P = dlV_dlT_P;
            self->state = std::move(state);
        },  "converged"_a = false,
            "gamma_s"_a = 0.0,
            "dlV_dlP_T"_a = 0.0,
            "dlV_dlT_P"_a = 0.0,
            "state"_a = std::vector<double>())
        .def_ro("converged", &Goddard::NozzleResult::converged)
        .def_ro("gamma_s", &Goddard::NozzleResult::gamma_s)
        .def_ro("dlV_dlP_T", &Goddard::NozzleResult::dlV_dlP_T)
        .def_ro("dlV_dlT_P", &Goddard::NozzleResult::dlV_dlT_P)
        .def_ro("state", &Goddard::NozzleResult::state);

    // NozzleResults
    nb::class_<Goddard::NozzleResults>(m, "NozzleResults")
        .def(nb::init<>())
        .def("__init__", [](Goddard::NozzleResults* self,
                            Goddard::ThroatCondition throat,
                            std::vector<Goddard::NozzleResult> expansions) {
            new (self) Goddard::NozzleResults();
            self->throat = std::move(throat);
            self->expansions = std::move(expansions);
        },  "throat"_a = Goddard::ThroatCondition(),
            "expansions"_a = std::vector<Goddard::NozzleResult>())
        .def_ro("throat", &Goddard::NozzleResults::throat)
        .def_ro("expansions", &Goddard::NozzleResults::expansions);

    // RocketCaseParameters
    nb::class_<Goddard::RocketCaseParameters>(m, "RocketCaseParameters")
        .def(nb::init<>())
        .def("__init__", [](Goddard::RocketCaseParameters* self,
                            std::string name,
                            std::string problem_type,
                            Goddard::CombustorOptions combustor_options,
                            Goddard::NozzleOptions nozzle_options) {
            new (self) Goddard::RocketCaseParameters();
            self->name = std::move(name);
            self->problem_type = std::move(problem_type);
            self->combustor_options = std::move(combustor_options);
            self->nozzle_options = std::move(nozzle_options);
        },  "name"_a = "",
            "problem_type"_a = "",
            "combustor_options"_a = Goddard::CombustorOptions(),
            "nozzle_options"_a = Goddard::NozzleOptions())
        .def_rw("name", &Goddard::RocketCaseParameters::name)
        .def_rw("problem_type", &Goddard::RocketCaseParameters::problem_type)
        .def_rw("combustor_options", &Goddard::RocketCaseParameters::combustor_options)
        .def_rw("nozzle_options", &Goddard::RocketCaseParameters::nozzle_options);

    // ChemicalParameters
    nb::class_<Goddard::ChemicalParameters>(m, "ChemicalParameters")
        .def(nb::init<>())
        .def("__init__", [](Goddard::ChemicalParameters* self,
                            std::string thermo_file,
                            std::unordered_set<std::string> species,
                            Goddard::ThermodynamicState cantera_fuel_state,
                            Goddard::ThermodynamicState cantera_oxidizer_state,
                            std::vector<double> OF_ratios,
                            std::vector<double> phi_ratios,
                            std::vector<double> fuel_weight_percentages,
                            std::vector<double> valance_equivalences) {
            new (self) Goddard::ChemicalParameters();
            self->thermo_file = std::move(thermo_file);
            self->species = std::move(species);
            self->cantera_fuel_state = std::move(cantera_fuel_state);
            self->cantera_oxidizer_state = std::move(cantera_oxidizer_state);
            self->OF_ratios = std::move(OF_ratios);
            self->phi_ratios = std::move(phi_ratios);
            self->fuel_weight_percentages = std::move(fuel_weight_percentages);
            self->valance_equivalences = std::move(valance_equivalences);
        },  "thermo_file"_a = "",
            "species"_a = std::unordered_set<std::string>(),
            "cantera_fuel_state"_a = Goddard::ThermodynamicState(),
            "cantera_oxidizer_state"_a = Goddard::ThermodynamicState(),
            "OF_ratios"_a = std::vector<double>(),
            "phi_ratios"_a = std::vector<double>(),
            "fuel_weight_percentages"_a = std::vector<double>(),
            "valance_equivalences"_a = std::vector<double>())
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
        .def("__init__", [](Goddard::ThermoStateInfo* self,
                            double pressure,
                            double temperature,
                            double density,
                            double enthalpy,
                            double internal_energy,
                            double gibbs,
                            double entropy,
                            double molecular_weight,
                            double cp,
                            double gamma_s,
                            double dlV_dlP_T,
                            double dlV_dlT_P,
                            double speed_of_sound,
                            double stagnation_enthalpy,
                            std::unordered_map<std::string, double> composition) {
            new (self) Goddard::ThermoStateInfo();
            self->pressure = pressure;
            self->temperature = temperature;
            self->density = density;
            self->enthalpy = enthalpy;
            self->internal_energy = internal_energy;
            self->gibbs = gibbs;
            self->entropy = entropy;
            self->molecular_weight = molecular_weight;
            self->cp = cp;
            self->gamma_s = gamma_s;
            self->dlV_dlP_T = dlV_dlP_T;
            self->dlV_dlT_P = dlV_dlT_P;
            self->speed_of_sound = speed_of_sound;
            self->stagnation_enthalpy = stagnation_enthalpy;
            self->composition = std::move(composition);
        },  "pressure"_a = 0.0,
            "temperature"_a = 0.0,
            "density"_a = 0.0,
            "enthalpy"_a = 0.0,
            "internal_energy"_a = 0.0,
            "gibbs"_a = 0.0,
            "entropy"_a = 0.0,
            "molecular_weight"_a = 0.0,
            "cp"_a = 0.0,
            "gamma_s"_a = 0.0,
            "dlV_dlP_T"_a = 0.0,
            "dlV_dlT_P"_a = 0.0,
            "speed_of_sound"_a = 0.0,
            "stagnation_enthalpy"_a = 0.0,
            "composition"_a = std::unordered_map<std::string, double>())
        .def_ro("pressure", &Goddard::ThermoStateInfo::pressure)
        .def_ro("temperature", &Goddard::ThermoStateInfo::temperature)
        .def_ro("density", &Goddard::ThermoStateInfo::density)
        .def_ro("enthalpy", &Goddard::ThermoStateInfo::enthalpy)
        .def_ro("internal_energy", &Goddard::ThermoStateInfo::internal_energy)
        .def_ro("gibbs", &Goddard::ThermoStateInfo::gibbs)
        .def_ro("entropy", &Goddard::ThermoStateInfo::entropy)
        .def_ro("molecular_weight", &Goddard::ThermoStateInfo::molecular_weight)
        .def_ro("cp", &Goddard::ThermoStateInfo::cp)
        .def_ro("gamma_s", &Goddard::ThermoStateInfo::gamma_s)
        .def_ro("dlV_dlP_T", &Goddard::ThermoStateInfo::dlV_dlP_T)
        .def_ro("dlV_dlT_P", &Goddard::ThermoStateInfo::dlV_dlT_P)
        .def_ro("speed_of_sound", &Goddard::ThermoStateInfo::speed_of_sound)
        .def_ro("stagnation_enthalpy", &Goddard::ThermoStateInfo::stagnation_enthalpy)
        .def_ro("composition", &Goddard::ThermoStateInfo::composition);

    // StationType
    nb::enum_<Goddard::StationType>(m, "StationType")
        .value("CHAMBER", Goddard::StationType::CHAMBER)
        .value("THROAT",  Goddard::StationType::THROAT)
        .value("EXIT",    Goddard::StationType::EXIT);

    // RocketStation
    nb::class_<Goddard::RocketStation>(m, "RocketStation")
        .def_ro("case_name",       &Goddard::RocketStation::case_name)
        .def_ro("type",            &Goddard::RocketStation::type)
        .def_ro("of_index",        &Goddard::RocketStation::of_index)
        .def_ro("pressure_index",  &Goddard::RocketStation::pressure_index)
        .def_ro("expansion_index", &Goddard::RocketStation::expansion_index)
        .def_ro("area_ratio",      &Goddard::RocketStation::area_ratio)
        .def_ro("thermo",          &Goddard::RocketStation::thermo)
        .def_ro("converged",       &Goddard::RocketStation::converged);

    // RocketPerformance
    nb::class_<Goddard::RocketPerformance>(m, "RocketPerformance")
        .def(nb::init<>())
        .def("__init__", [](Goddard::RocketPerformance* self,
                            double pressure_ratio,
                            double area_ratio,
                            double mach_number,
                            double cstar,
                            double CF,
                            double isp,
                            double ivac) {
            new (self) Goddard::RocketPerformance();
            self->pressure_ratio = pressure_ratio;
            self->area_ratio = area_ratio;
            self->mach_number = mach_number;
            self->cstar = cstar;
            self->CF = CF;
            self->isp = isp;
            self->ivac = ivac;
        },  "pressure_ratio"_a = 0.0,
            "area_ratio"_a = 0.0,
            "mach_number"_a = 0.0,
            "cstar"_a = 0.0,
            "CF"_a = 0.0,
            "isp"_a = 0.0,
            "ivac"_a = 0.0)
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
        .def("__init__", [](Goddard::RocketState* self,
                            std::string name,
                            std::vector<double> cantera_state,
                            double pressure_ratio,
                            double area_ratio,
                            double dlv_dlp_t,
                            double dlv_dlt_p,
                            double gamma_s,
                            double speed_of_sound) {
            new (self) Goddard::RocketState();
            self->name = std::move(name);
            self->cantera_state = std::move(cantera_state);
            self->pressure_ratio = pressure_ratio;
            self->area_ratio = area_ratio;
            self->dlv_dlp_t = dlv_dlp_t;
            self->dlv_dlt_p = dlv_dlt_p;
            self->gamma_s = gamma_s;
            self->speed_of_sound = speed_of_sound;
        },  "name"_a = "",
            "cantera_state"_a = std::vector<double>(),
            "pressure_ratio"_a = 0.0,
            "area_ratio"_a = 0.0,
            "dlv_dlp_t"_a = 0.0,
            "dlv_dlt_p"_a = 0.0,
            "gamma_s"_a = 0.0,
            "speed_of_sound"_a = 0.0)
        .def_ro("name", &Goddard::RocketState::name)
        .def_ro("cantera_state", &Goddard::RocketState::cantera_state)
        .def_ro("pressure_ratio", &Goddard::RocketState::pressure_ratio)
        .def_ro("area_ratio", &Goddard::RocketState::area_ratio)
        .def_ro("dlv_dlp_t", &Goddard::RocketState::dlv_dlp_t)
        .def_ro("dlv_dlt_p", &Goddard::RocketState::dlv_dlt_p)
        .def_ro("gamma_s", &Goddard::RocketState::gamma_s)
        .def_ro("speed_of_sound", &Goddard::RocketState::speed_of_sound);
}
