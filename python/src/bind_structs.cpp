#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <format>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <vector>

#include "goddard/combustor.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/thermo.hpp"
#include "bind_status.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;


void bind_structs(nb::module_& m) {
    // CombustorOptions
    nb::class_<Goddard::CombustorOptions>(m, "CombustorOptions", DOC(Goddard, CombustorOptions))
        .def(nb::self == nb::self)
        .def("__init__", [](Goddard::CombustorOptions* self,
                            Goddard::CombustorType type,
                            Goddard::MixtureRatioType mixture_type,
                            double mass_flux,
                            double contraction_ratio,
                            Goddard::CombustionProcess process) {
            new (self) Goddard::CombustorOptions();
            self->type = type;
            self->mixture_type = mixture_type;
            self->mass_flux = mass_flux;
            self->contraction_ratio = contraction_ratio;
            self->process = process;
        },  "type"_a = Goddard::CombustorType::INFINITE_AREA,
            "mixture_type"_a = Goddard::MixtureRatioType::OF_RATIO,
            "mass_flux"_a = 0.0,
            "contraction_ratio"_a = 0.0,
            "process"_a = Goddard::CombustionProcess::ISOBARIC,
            "Create options; each keyword sets the field of the same name.")
        .def_rw("type", &Goddard::CombustorOptions::type,
             DOC(Goddard, CombustorOptions, type))
        .def_rw("mixture_type", &Goddard::CombustorOptions::mixture_type,
             DOC(Goddard, CombustorOptions, mixture_type))
        .def_rw("mass_flux", &Goddard::CombustorOptions::mass_flux,
             DOC(Goddard, CombustorOptions, mass_flux))
        .def_rw("contraction_ratio", &Goddard::CombustorOptions::contraction_ratio,
             DOC(Goddard, CombustorOptions, contraction_ratio))
        .def_rw("process", &Goddard::CombustorOptions::process,
             DOC(Goddard, CombustorOptions, process));

    // NozzleOptions
    nb::class_<Goddard::NozzleOptions>(m, "NozzleOptions", DOC(Goddard, NozzleOptions))
        .def(nb::self == nb::self)
        .def("__init__", [](Goddard::NozzleOptions* self,
                            Goddard::GasChemistry chemistry,
                            int frozen_NFZ) {
            new (self) Goddard::NozzleOptions();
            self->chemistry = chemistry;
            self->frozen_NFZ = frozen_NFZ;
        },  "chemistry"_a = Goddard::GasChemistry::EQUILIBRIUM,
            "frozen_NFZ"_a = 0,
            "Create options; each keyword sets the field of the same name.")
        .def_rw("chemistry", &Goddard::NozzleOptions::chemistry,
             DOC(Goddard, NozzleOptions, chemistry))
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
        .def_ro("status", &Goddard::ThroatCondition::status,
             DOC(Goddard, ThroatCondition, status))
        .def("__repr__", [](const Goddard::ThroatCondition& self) {
            return std::format("<ThroatCondition {} T={:.6g} K P={:.6g} Pa u={:.6g} m/s M={:.6g} gamma_s={:.6g}>", status_tag(self.status), self.thermo.temperature, self.thermo.pressure, self.velocity, self.mach, self.gamma_s);
        });

    // NozzleStation
    nb::class_<Goddard::NozzleStation>(m, "NozzleStation", DOC(Goddard, NozzleStation))
        .def_ro("thermo", &Goddard::NozzleStation::thermo, DOC(Goddard, NozzleStation, thermo))
        .def_ro("velocity", &Goddard::NozzleStation::velocity, DOC(Goddard, NozzleStation, velocity))
        .def_ro("mach", &Goddard::NozzleStation::mach, DOC(Goddard, NozzleStation, mach))
        .def_ro("area_ratio", &Goddard::NozzleStation::area_ratio,
             DOC(Goddard, NozzleStation, area_ratio))
        .def_ro("state", &Goddard::NozzleStation::state, DOC(Goddard, NozzleStation, state))
        .def_ro("status", &Goddard::NozzleStation::status, DOC(Goddard, NozzleStation, status))
        .def("__repr__", [](const Goddard::NozzleStation& self) {
            return std::format("<NozzleStation {} A/At={:.6g} T={:.6g} K P={:.6g} Pa u={:.6g} m/s M={:.6g}>", status_tag(self.status), self.area_ratio, self.thermo.temperature, self.thermo.pressure, self.velocity, self.mach);
        });

    // NozzleResults
    nb::class_<Goddard::NozzleResults>(m, "NozzleResults", DOC(Goddard, NozzleResults))
        .def_ro("inlet", &Goddard::NozzleResults::inlet, DOC(Goddard, NozzleResults, inlet))
        .def_ro("throat", &Goddard::NozzleResults::throat,
             DOC(Goddard, NozzleResults, throat))
        .def_ro("expansions", &Goddard::NozzleResults::expansions,
             DOC(Goddard, NozzleResults, expansions))
        .def_ro("status", &Goddard::NozzleResults::status, DOC(Goddard, NozzleResults, status))
        .def("__repr__", [](const Goddard::NozzleResults& self) {
            return std::format("<NozzleResults {} chamber P={:.6g} Pa, throat T={:.6g} K, {} expansion station(s)>", status_tag(self.status), self.inlet.thermo.pressure, self.throat.thermo.temperature, self.expansions.size());
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
             DOC(Goddard, FiniteAreaChamber, iterations))
        .def_ro("status", &Goddard::FiniteAreaChamber::status,
             DOC(Goddard, FiniteAreaChamber, status))
        .def("__repr__", [](const Goddard::FiniteAreaChamber& self) {
            return std::format("<FiniteAreaChamber {} Pinj={:.6g} Pa Pinf={:.6g} Pa Ac/At={:.6g}>", status_tag(self.status), self.injector_pressure, self.stagnation_pressure, self.contraction_ratio);
        });

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
        .def_ro("composition_model", &Goddard::ThermodynamicState::composition_model,
             DOC(Goddard, ThermodynamicState, composition_model))
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
}
