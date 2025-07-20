#pragma once
#include <cmath>
#include "cantera/core.h"

namespace Goddard {

/**
 * @brief Gas flow velocity assuming isenthalpic expansion.
 */
double gas_isenthalpic_velocity(const Cantera::ThermoPhase& gas, double H_stagnation);

/**
 * @brief Speed of sound of a gas.
 */
double gas_sonic_velocity(const Cantera::ThermoPhase& gas, double gamma);

/**
 * @brief Area per unit mass flow rate. Equation 6.12 in NASA CEA Report Part I.
 */
double area_per_mdot(const Cantera::ThermoPhase& gas, double velocity);

double cstar(double gamma, double temperature, double molecular_weight);

double isp(const Cantera::ThermoPhase& gas, double gamma, double enthalpy);

double ivac(const Cantera::ThermoPhase& gas, double gamma, double enthalpy);

}