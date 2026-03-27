#pragma once
#include <cmath>
#include "eigen3/Eigen/Dense"
#include "cantera/core.h"
#include "goddard/thermoarray.hpp"
namespace Goddard {

/**
 * @brief Gas flow velocity assuming isenthalpic expansion.
 */
double gas_isenthalpic_velocity(const Cantera::ThermoPhase& gas, double H_stagnation);

Eigen::ArrayXXd gas_isenthalpic_velocity(const ThermoArray& gas, const Eigen::ArrayXXd& H_stagnation);

double gas_stagnation_enthalpy(const Cantera::ThermoPhase& gas, double velocity);

Eigen::ArrayXXd gas_stagnation_enthalpy(const ThermoArray& gas, const Eigen::ArrayXXd velocity);

double stagnation_pressure(const Cantera::ThermoPhase& gas, double mach, double gamma);

Eigen::ArrayXXd stagnation_pressure(const Cantera::ThermoPhase& gas, const Eigen::ArrayXXd mach, const Eigen::ArrayXXd gamma);

double stagnation_factor(double mach, double gamma);

Eigen::ArrayXXd stagnation_factor(const Eigen::ArrayXXd& mach, const Eigen::ArrayXXd gamma);

inline double mach_to_mu(double mach) {
    // unfortunately trig functions aren't constexpr until C++26.
    return std::asin(1.0/mach);
}


/**
 * @brief Speed of sound of a gas.
 */
double gas_sonic_velocity(const Cantera::ThermoPhase& gas, double gamma);

Eigen::ArrayXXd gas_sonic_velocity(const ThermoArray& gas, const Eigen::ArrayXXd& gamma);


/**
 * @brief Area per unit mass flow rate. Equation 6.12 in NASA CEA Report Part I.
 */
double area_per_mdot(const Cantera::ThermoPhase& gas, double velocity);

Eigen::ArrayXXd area_per_mdot(const ThermoArray& gas, const Eigen::ArrayXXd& velocity);


double isp(const Cantera::ThermoPhase& gas, double gamma, double enthalpy);

double ivac(const Cantera::ThermoPhase& gas, double gamma, double enthalpy);

double mach(const Cantera::ThermoPhase& gas, double H_stag, double gamma);

/**
 * @brief Calculate thrust coefficient.
 */
double C_F(double gamma, double pressure_ratio, double area_ratio);

/**
 * @brief Calculate the characteristic combustion velocity (c*)
 */
double cstar(double gamma, double temperature, double molecular_weight);

}