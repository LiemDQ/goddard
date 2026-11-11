#include "cantera/base/stringUtils.h"
#include "goddard/thermo.hpp"
#include "goddard/error.hpp"

#include <limits>

namespace Goddard {

ThermodynamicState ThermodynamicState::undefined() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    ThermodynamicState state;
    state.pressure = nan;
    state.temperature = nan;
    state.density = nan;
    state.enthalpy = nan;
    state.internal_energy = nan;
    state.gibbs = nan;
    state.entropy = nan;
    state.molecular_weight = nan;
    state.cp = nan;
    state.gamma_s = nan;
    state.dlV_dlP_T = nan;
    state.dlV_dlT_P = nan;
    state.speed_of_sound = nan;
    state.stagnation_enthalpy = nan;
    state.mixture_molecular_weight = nan;
    state.gas_mass_fraction = nan;
    return state;
}

// -- PhaseSpecification --

PhaseSpecification::PhaseSpecification(double temperature, double pressure, const std::string& comp) 
: T(temperature), P(pressure), composition(Cantera::parseCompString(comp))
{}

PhaseSpecification::PhaseSpecification(double temperature, double pressure, const Composition& comp) 
: T(temperature), P(pressure), composition(comp)
{}


} //namespace goddard