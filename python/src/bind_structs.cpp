#include <nanobind/nanobind.h>
#include <format>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/unordered_set.h>
#include <unordered_set>
#include <vector>

#include "goddard/combustor.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/problem.hpp"
#include "goddard/rocket_results.hpp"
#include "goddard/thermo.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

const char* station_type_name(Goddard::StationType type) {
    switch (type) {
        case Goddard::StationType::CHAMBER: return "CHAMBER";
        case Goddard::StationType::STAGNATION: return "STAGNATION";
        case Goddard::StationType::COMBUSTION_END: return "COMBUSTION_END";
        case Goddard::StationType::THROAT: return "THROAT";
        case Goddard::StationType::EXIT: return "EXIT";
    }
    return "UNKNOWN";
}

} // namespace

void bind_structs(nb::module_& m) {
    // CombustorOptions
    nb::class_<Goddard::CombustorOptions>(m, "CombustorOptions", DOC(Goddard, CombustorOptions))
        .def("__init__", [](Goddard::CombustorOptions* self,
                            Goddard::CombustorType type,
                            Goddard::MixtureRatioType mixture_type,
                            std::vector<double> pressures,
                            double mass_flux,
                            double contraction_ratio,
                            Goddard::CombustionProcess process) {
            new (self) Goddard::CombustorOptions();
            self->type = type;
            self->mixture_type = mixture_type;
            self->pressures = std::move(pressures);
            self->mass_flux = mass_flux;
            self->contraction_ratio = contraction_ratio;
            self->process = process;
        },  "type"_a = Goddard::CombustorType::INFINITE_AREA,
            "mixture_type"_a = Goddard::MixtureRatioType::OF_RATIO,
            "pressures"_a = std::vector<double>(),
            "mass_flux"_a = 0.0,
            "contraction_ratio"_a = 0.0,
            "process"_a = Goddard::CombustionProcess::ISOBARIC,
            "Create options; each keyword sets the field of the same name.")
        .def_rw("type", &Goddard::CombustorOptions::type,
             DOC(Goddard, CombustorOptions, type))
        .def_rw("mixture_type", &Goddard::CombustorOptions::mixture_type,
             DOC(Goddard, CombustorOptions, mixture_type))
        .def_rw("pressures", &Goddard::CombustorOptions::pressures,
             DOC(Goddard, CombustorOptions, pressures))
        .def_rw("mass_flux", &Goddard::CombustorOptions::mass_flux,
             DOC(Goddard, CombustorOptions, mass_flux))
        .def_rw("contraction_ratio", &Goddard::CombustorOptions::contraction_ratio,
             DOC(Goddard, CombustorOptions, contraction_ratio))
        .def_rw("process", &Goddard::CombustorOptions::process,
             DOC(Goddard, CombustorOptions, process));

    // NozzleOptions
    nb::class_<Goddard::NozzleOptions>(m, "NozzleOptions", DOC(Goddard, NozzleOptions))
        .def("__init__", [](Goddard::NozzleOptions* self,
                            Goddard::GasChemistry chemistry,
                            Goddard::ExpansionType expansion_type,
                            std::vector<double> expansion_ratios,
                            int frozen_NFZ) {
            new (self) Goddard::NozzleOptions();
            self->chemistry = chemistry;
            self->expansion_type = expansion_type;
            self->expansion_ratios = std::move(expansion_ratios);
            self->frozen_NFZ = frozen_NFZ;
        },  "chemistry"_a = Goddard::GasChemistry::EQUILIBRIUM,
            "expansion_type"_a = Goddard::ExpansionType::SUPERSONIC_AREA_RATIO,
            "expansion_ratios"_a = std::vector<double>(),
            "frozen_NFZ"_a = 0,
            "Create options; each keyword sets the field of the same name.")
        .def_rw("chemistry", &Goddard::NozzleOptions::chemistry,
             DOC(Goddard, NozzleOptions, chemistry))
        .def_rw("expansion_type", &Goddard::NozzleOptions::expansion_type,
             DOC(Goddard, NozzleOptions, expansion_type))
        .def_rw("expansion_ratios", &Goddard::NozzleOptions::expansion_ratios,
             DOC(Goddard, NozzleOptions, expansion_ratios))
        .def_rw("frozen_NFZ", &Goddard::NozzleOptions::frozen_NFZ,
             DOC(Goddard, NozzleOptions, frozen_NFZ));

    // ThroatCondition
    // `_2`: moc.hpp forward-declares ThroatCondition first, so the documented definition is the second.
    nb::class_<Goddard::ThroatCondition>(m, "ThroatCondition", DOC(Goddard, ThroatCondition, 2))
        .def_ro("speed_of_sound", &Goddard::ThroatCondition::speed_of_sound,
             DOC(Goddard, ThroatCondition, speed_of_sound))
        .def_ro("H_stagnation", &Goddard::ThroatCondition::H_stagnation,
             DOC(Goddard, ThroatCondition, H_stagnation))
        .def_ro("P_inlet", &Goddard::ThroatCondition::P_inlet,
             DOC(Goddard, ThroatCondition, P_inlet))
        .def_ro("S_inlet", &Goddard::ThroatCondition::S_inlet,
             DOC(Goddard, ThroatCondition, S_inlet))
        .def_ro("gamma_s", &Goddard::ThroatCondition::gamma_s,
             DOC(Goddard, ThroatCondition, gamma_s))
        .def_ro("dlV_dlP_T", &Goddard::ThroatCondition::dlV_dlP_T,
             DOC(Goddard, ThroatCondition, dlV_dlP_T))
        .def_ro("dlV_dlT_P", &Goddard::ThroatCondition::dlV_dlT_P,
             DOC(Goddard, ThroatCondition, dlV_dlT_P))
        .def_ro("pinned_transition", &Goddard::ThroatCondition::pinned_transition,
             DOC(Goddard, ThroatCondition, pinned_transition))
        .def_ro("thermo", &Goddard::ThroatCondition::thermo, DOC(Goddard, ThroatCondition, thermo))
        .def_ro("velocity", &Goddard::ThroatCondition::velocity, DOC(Goddard, ThroatCondition, velocity))
        .def_ro("mach", &Goddard::ThroatCondition::mach, DOC(Goddard, ThroatCondition, mach))
        .def_ro("area_ratio", &Goddard::ThroatCondition::area_ratio,
             DOC(Goddard, ThroatCondition, area_ratio))
        .def_ro("state", &Goddard::ThroatCondition::state,
             DOC(Goddard, ThroatCondition, state))
        .def("__repr__", [](const Goddard::ThroatCondition& self) {
            return std::format("<ThroatCondition T={:.6g} K P={:.6g} Pa u={:.6g} m/s M={:.6g} gamma_s={:.6g}>", self.thermo.temperature, self.thermo.pressure, self.velocity, self.mach, self.gamma_s);
        });

    // NozzleStation
    nb::class_<Goddard::NozzleStation>(m, "NozzleStation", DOC(Goddard, NozzleStation))
        .def_ro("thermo", &Goddard::NozzleStation::thermo, DOC(Goddard, NozzleStation, thermo))
        .def_ro("velocity", &Goddard::NozzleStation::velocity, DOC(Goddard, NozzleStation, velocity))
        .def_ro("mach", &Goddard::NozzleStation::mach, DOC(Goddard, NozzleStation, mach))
        .def_ro("area_ratio", &Goddard::NozzleStation::area_ratio,
             DOC(Goddard, NozzleStation, area_ratio))
        .def_ro("state", &Goddard::NozzleStation::state, DOC(Goddard, NozzleStation, state))
        .def("__repr__", [](const Goddard::NozzleStation& self) {
            return std::format("<NozzleStation A/At={:.6g} T={:.6g} K P={:.6g} Pa u={:.6g} m/s M={:.6g}>", self.area_ratio, self.thermo.temperature, self.thermo.pressure, self.velocity, self.mach);
        });

    // NozzleResults
    nb::class_<Goddard::NozzleResults>(m, "NozzleResults", DOC(Goddard, NozzleResults))
        .def_ro("inlet", &Goddard::NozzleResults::inlet, DOC(Goddard, NozzleResults, inlet))
        .def_ro("throat", &Goddard::NozzleResults::throat,
             DOC(Goddard, NozzleResults, throat))
        .def_ro("expansions", &Goddard::NozzleResults::expansions,
             DOC(Goddard, NozzleResults, expansions))
        .def("__repr__", [](const Goddard::NozzleResults& self) {
            return std::format("<NozzleResults chamber P={:.6g} Pa, throat T={:.6g} K, {} expansion station(s)>", self.inlet.thermo.pressure, self.throat.thermo.temperature, self.expansions.size());
        });

    // FiniteAreaChamber
    nb::class_<Goddard::FiniteAreaChamber>(m, "FiniteAreaChamber",
            DOC(Goddard, FiniteAreaChamber))
        .def_ro("injector", &Goddard::FiniteAreaChamber::injector,
             DOC(Goddard, FiniteAreaChamber, injector))
        .def_ro("stagnation", &Goddard::FiniteAreaChamber::stagnation,
             DOC(Goddard, FiniteAreaChamber, stagnation))
        .def_ro("combustion_end", &Goddard::FiniteAreaChamber::combustion_end,
             DOC(Goddard, FiniteAreaChamber, combustion_end))
        .def_ro("throat", &Goddard::FiniteAreaChamber::throat,
             DOC(Goddard, FiniteAreaChamber, throat))
        .def_ro("injector_pressure", &Goddard::FiniteAreaChamber::injector_pressure,
             DOC(Goddard, FiniteAreaChamber, injector_pressure))
        .def_ro("stagnation_pressure", &Goddard::FiniteAreaChamber::stagnation_pressure,
             DOC(Goddard, FiniteAreaChamber, stagnation_pressure))
        .def_ro("contraction_ratio", &Goddard::FiniteAreaChamber::contraction_ratio,
             DOC(Goddard, FiniteAreaChamber, contraction_ratio))
        .def_ro("mass_flux", &Goddard::FiniteAreaChamber::mass_flux,
             DOC(Goddard, FiniteAreaChamber, mass_flux))
        .def_ro("iterations", &Goddard::FiniteAreaChamber::iterations,
             DOC(Goddard, FiniteAreaChamber, iterations));

    // RocketCaseParameters
    nb::class_<Goddard::RocketCaseParameters>(m, "RocketCaseParameters", DOC(Goddard, RocketCaseParameters))
        .def("__init__", [](Goddard::RocketCaseParameters* self,
                            std::string name,
                            Goddard::CombustorOptions combustor_options,
                            Goddard::NozzleOptions nozzle_options) {
            new (self) Goddard::RocketCaseParameters();
            self->name = std::move(name);
            self->combustor_options = std::move(combustor_options);
            self->nozzle_options = std::move(nozzle_options);
        },  "name"_a = "",
            "combustor_options"_a = Goddard::CombustorOptions(),
            "nozzle_options"_a = Goddard::NozzleOptions(),
            "Create a case; each keyword sets the field of the same name.")
        .def_rw("name", &Goddard::RocketCaseParameters::name,
             DOC(Goddard, RocketCaseParameters, name))
        .def_rw("combustor_options", &Goddard::RocketCaseParameters::combustor_options,
             DOC(Goddard, RocketCaseParameters, combustor_options))
        .def_rw("nozzle_options", &Goddard::RocketCaseParameters::nozzle_options,
             DOC(Goddard, RocketCaseParameters, nozzle_options));

    // ChemicalParameters
    nb::class_<Goddard::ChemicalParameters>(m, "ChemicalParameters", DOC(Goddard, ChemicalParameters))
        .def("__init__", [](Goddard::ChemicalParameters* self,
                            std::string thermo_file,
                            std::unordered_set<std::string> species,
                            Goddard::PhaseSpecification fuel_state,
                            Goddard::PhaseSpecification oxidizer_state,
                            std::vector<double> mixtures,
                            std::vector<double> OF_ratios,
                            std::vector<double> phi_ratios,
                            std::vector<double> fuel_weight_percentages,
                            std::string reactant_file,
                            std::string condensed_file,
                            std::unordered_set<std::string> condensed_species,
                            bool all_condensed_species)
                            {
            new (self) Goddard::ChemicalParameters();
            self->thermo_file = std::move(thermo_file);
            self->species = std::move(species);
            self->fuel_state = std::move(fuel_state);
            self->oxidizer_state = std::move(oxidizer_state);
            self->mixtures = std::move(mixtures);
            self->OF_ratios = std::move(OF_ratios);
            self->phi_ratios = std::move(phi_ratios);
            self->fuel_weight_percentages = std::move(fuel_weight_percentages);
            self->reactant_file = std::move(reactant_file);
            self->condensed_file = std::move(condensed_file);
            self->condensed_species = std::move(condensed_species);
            self->all_condensed_species = all_condensed_species;
        },  "thermo_file"_a = "",
            "species"_a = std::unordered_set<std::string>(),
            "fuel_state"_a = Goddard::PhaseSpecification(),
            "oxidizer_state"_a = Goddard::PhaseSpecification(),
            "mixtures"_a = std::vector<double>(),
            "OF_ratios"_a = std::vector<double>(),
            "phi_ratios"_a = std::vector<double>(),
            "fuel_weight_percentages"_a = std::vector<double>(),
            "reactant_file"_a = "",
            "condensed_file"_a = "",
            "condensed_species"_a = std::unordered_set<std::string>(),
            "all_condensed_species"_a = false,
            "Create parameters; each keyword sets the field of the same name.")
        .def_rw("thermo_file", &Goddard::ChemicalParameters::thermo_file,
             DOC(Goddard, ChemicalParameters, thermo_file))
        .def_rw("species", &Goddard::ChemicalParameters::species,
             DOC(Goddard, ChemicalParameters, species))
        .def_rw("fuel_state", &Goddard::ChemicalParameters::fuel_state,
             DOC(Goddard, ChemicalParameters, fuel_state))
        .def_rw("oxidizer_state", &Goddard::ChemicalParameters::oxidizer_state,
             DOC(Goddard, ChemicalParameters, oxidizer_state))
        .def_rw("mixtures", &Goddard::ChemicalParameters::mixtures,
             DOC(Goddard, ChemicalParameters, mixtures))
        .def_rw("OF_ratios", &Goddard::ChemicalParameters::OF_ratios,
             DOC(Goddard, ChemicalParameters, OF_ratios))
        .def_rw("phi_ratios", &Goddard::ChemicalParameters::phi_ratios,
             DOC(Goddard, ChemicalParameters, phi_ratios))
        .def_rw("fuel_weight_percentages", &Goddard::ChemicalParameters::fuel_weight_percentages,
             DOC(Goddard, ChemicalParameters, fuel_weight_percentages))
        .def_rw("reactant_file", &Goddard::ChemicalParameters::reactant_file,
             DOC(Goddard, ChemicalParameters, reactant_file))
        .def_rw("condensed_file", &Goddard::ChemicalParameters::condensed_file,
             DOC(Goddard, ChemicalParameters, condensed_file))
        .def_rw("condensed_species", &Goddard::ChemicalParameters::condensed_species,
             DOC(Goddard, ChemicalParameters, condensed_species))
        .def_rw("all_condensed_species", &Goddard::ChemicalParameters::all_condensed_species,
             DOC(Goddard, ChemicalParameters, all_condensed_species));

    // ThermoStateInfo
    nb::class_<Goddard::ThermodynamicState>(m, "ThermodynamicState", DOC(Goddard, ThermodynamicState))
        .def("__init__", [](Goddard::ThermodynamicState* self,
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
                            std::map<std::string, double> composition,
                            double mixture_molecular_weight,
                            double gas_mass_fraction,
                            bool pinned_transition) {
            new (self) Goddard::ThermodynamicState();
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
            self->mixture_molecular_weight = mixture_molecular_weight;
            self->gas_mass_fraction = gas_mass_fraction;
            self->pinned_transition = pinned_transition;
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
            "composition"_a = std::map<std::string, double>(),
            "mixture_molecular_weight"_a = 0.0,
            "gas_mass_fraction"_a = 1.0,
            "pinned_transition"_a = false,
            "Create a state; each keyword sets the field of the same name.")
        .def_ro("pressure", &Goddard::ThermodynamicState::pressure,
             DOC(Goddard, ThermodynamicState, pressure))
        .def_ro("temperature", &Goddard::ThermodynamicState::temperature,
             DOC(Goddard, ThermodynamicState, temperature))
        .def_ro("density", &Goddard::ThermodynamicState::density,
             DOC(Goddard, ThermodynamicState, density))
        .def_ro("enthalpy", &Goddard::ThermodynamicState::enthalpy,
             DOC(Goddard, ThermodynamicState, enthalpy))
        .def_ro("internal_energy", &Goddard::ThermodynamicState::internal_energy,
             DOC(Goddard, ThermodynamicState, internal_energy))
        .def_ro("gibbs", &Goddard::ThermodynamicState::gibbs,
             DOC(Goddard, ThermodynamicState, gibbs))
        .def_ro("entropy", &Goddard::ThermodynamicState::entropy,
             DOC(Goddard, ThermodynamicState, entropy))
        .def_ro("molecular_weight", &Goddard::ThermodynamicState::molecular_weight,
             DOC(Goddard, ThermodynamicState, molecular_weight))
        .def_ro("cp", &Goddard::ThermodynamicState::cp,
             DOC(Goddard, ThermodynamicState, cp))
        .def_ro("gamma_s", &Goddard::ThermodynamicState::gamma_s,
             DOC(Goddard, ThermodynamicState, gamma_s))
        .def_ro("dlV_dlP_T", &Goddard::ThermodynamicState::dlV_dlP_T,
             DOC(Goddard, ThermodynamicState, dlV_dlP_T))
        .def_ro("dlV_dlT_P", &Goddard::ThermodynamicState::dlV_dlT_P,
             DOC(Goddard, ThermodynamicState, dlV_dlT_P))
        .def_ro("speed_of_sound", &Goddard::ThermodynamicState::speed_of_sound,
             DOC(Goddard, ThermodynamicState, speed_of_sound))
        .def_ro("stagnation_enthalpy", &Goddard::ThermodynamicState::stagnation_enthalpy,
             DOC(Goddard, ThermodynamicState, stagnation_enthalpy))
        .def_ro("composition", &Goddard::ThermodynamicState::composition,
             DOC(Goddard, ThermodynamicState, composition))
        .def_ro("mixture_molecular_weight",
             &Goddard::ThermodynamicState::mixture_molecular_weight,
             DOC(Goddard, ThermodynamicState, mixture_molecular_weight))
        .def_ro("gas_mass_fraction", &Goddard::ThermodynamicState::gas_mass_fraction,
             DOC(Goddard, ThermodynamicState, gas_mass_fraction))
        .def_ro("pinned_transition", &Goddard::ThermodynamicState::pinned_transition,
             DOC(Goddard, ThermodynamicState, pinned_transition))
        .def("__repr__", [](const Goddard::ThermodynamicState& self) {
            return std::format("<ThermodynamicState T={:.6g} K P={:.6g} Pa rho={:.6g} kg/m3 h={:.6g} J/kg M={:.6g} gamma_s={:.6g}>", self.temperature, self.pressure, self.density, self.enthalpy, self.molecular_weight, self.gamma_s);
        });

    // StationType
    nb::enum_<Goddard::StationType>(m, "StationType", DOC(Goddard, StationType))
        .value("CHAMBER",        Goddard::StationType::CHAMBER, DOC(Goddard, StationType, CHAMBER))
        .value("STAGNATION",     Goddard::StationType::STAGNATION, DOC(Goddard, StationType, STAGNATION))
        .value("COMBUSTION_END", Goddard::StationType::COMBUSTION_END, DOC(Goddard, StationType, COMBUSTION_END))
        .value("THROAT",         Goddard::StationType::THROAT, DOC(Goddard, StationType, THROAT))
        .value("EXIT",           Goddard::StationType::EXIT, DOC(Goddard, StationType, EXIT));

    // RocketStation
    nb::class_<Goddard::RocketStation>(m, "RocketStation", DOC(Goddard, RocketStation))
        .def_ro("case_name",       &Goddard::RocketStation::case_name,
             DOC(Goddard, RocketStation, case_name))
        .def_ro("type",            &Goddard::RocketStation::type,
             DOC(Goddard, RocketStation, type))
        .def_ro("of_index",        &Goddard::RocketStation::of_index,
             DOC(Goddard, RocketStation, of_index))
        .def_ro("pressure_index",  &Goddard::RocketStation::pressure_index,
             DOC(Goddard, RocketStation, pressure_index))
        .def_ro("expansion_index", &Goddard::RocketStation::expansion_index,
             DOC(Goddard, RocketStation, expansion_index))
        .def_ro("area_ratio",      &Goddard::RocketStation::area_ratio,
             DOC(Goddard, RocketStation, area_ratio))
        .def_ro("thermo",          &Goddard::RocketStation::thermo,
             DOC(Goddard, RocketStation, thermo))
        .def("__repr__", [](const Goddard::RocketStation& self) {
            return std::format("<RocketStation case='{}' type={} of_index={} pressure_index={} expansion_index={} A/At={:.6g} T={:.6g} K P={:.6g} Pa>", self.case_name, station_type_name(self.type), self.of_index, self.pressure_index, self.expansion_index, self.area_ratio, self.thermo.temperature, self.thermo.pressure);
        });

    // RocketPerformance
    nb::class_<Goddard::RocketPerformance>(m, "RocketPerformance", DOC(Goddard, RocketPerformance))
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
            "ivac"_a = 0.0,
            "Create a performance record; each keyword sets the field of the same name.")
        .def_ro("pressure_ratio", &Goddard::RocketPerformance::pressure_ratio,
             DOC(Goddard, RocketPerformance, pressure_ratio))
        .def_ro("area_ratio", &Goddard::RocketPerformance::area_ratio,
             DOC(Goddard, RocketPerformance, area_ratio))
        .def_ro("mach_number", &Goddard::RocketPerformance::mach_number,
             DOC(Goddard, RocketPerformance, mach_number))
        .def_ro("cstar", &Goddard::RocketPerformance::cstar,
             DOC(Goddard, RocketPerformance, cstar))
        .def_ro("CF", &Goddard::RocketPerformance::CF,
             DOC(Goddard, RocketPerformance, CF))
        .def_ro("isp", &Goddard::RocketPerformance::isp, DOC(Goddard, RocketPerformance, isp))
        .def_ro("ivac", &Goddard::RocketPerformance::ivac, DOC(Goddard, RocketPerformance, ivac))
        .def_prop_ro("isp_s", &Goddard::RocketPerformance::isp_s,
             DOC(Goddard, RocketPerformance, isp_s))
        .def_prop_ro("ivac_s", &Goddard::RocketPerformance::ivac_s,
             DOC(Goddard, RocketPerformance, ivac_s))
        .def("__repr__", [](const Goddard::RocketPerformance& self) {
            return std::format("<RocketPerformance Pc/Pe={:.6g} Ae/At={:.6g} M={:.6g} cstar={:.6g} m/s CF={:.6g} Isp={:.6g} m/s Ivac={:.6g} m/s>", self.pressure_ratio, self.area_ratio, self.mach_number, self.cstar, self.CF, self.isp, self.ivac);
        });

}
