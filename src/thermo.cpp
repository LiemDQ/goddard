#include "cantera/base/stringUtils.h"
#include "goddard/thermo.hpp"
#include "goddard/error.hpp"

namespace Goddard {

// -- PhaseSpecification --

PhaseSpecification::PhaseSpecification(double temperature, double pressure, const std::string& comp) 
: T(temperature), P(pressure), composition(Cantera::parseCompString(comp))
{}

PhaseSpecification::PhaseSpecification(double temperature, double pressure, const Composition& comp) 
: T(temperature), P(pressure), composition(comp)
{}


} //namespace goddard