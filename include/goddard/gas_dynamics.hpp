#pragma once
#include <cmath>
#include <limits>
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

Eigen::ArrayXXd gas_stagnation_enthalpy(const ThermoArray& gas, const Eigen::ArrayXXd& velocity);

double perfect_gas_stagnation_pressure(double P, double mach, double gamma);

Eigen::ArrayXXd perfect_gas_stagnation_pressure(const Eigen::ArrayXXd& P, const Eigen::ArrayXXd& mach, const Eigen::ArrayXXd& gamma);

double gas_stagnation_pressure(const Cantera::ThermoPhase& gas, double velocity);

double stagnation_factor(double mach, double gamma);

Eigen::ArrayXXd stagnation_factor(const Eigen::ArrayXXd& mach, const Eigen::ArrayXXd& gamma);

inline double mach_to_mu(double mach) {
    // unfortunately trig functions aren't constexpr until C++26.
    return std::asin(1.0/mach);
}

/**
 * @brief Mach number from the critical velocity ratio M* = V/a*.
 *
 * M* (also written V/a*, the "characteristic Mach number") normalizes velocity by the
 * *critical* speed of sound -- the value a would take where the flow is sonic -- rather
 * than by the local one. The two agree only at M = 1 and diverge quickly above it: M*
 * is bounded by sqrt((gamma+1)/(gamma-1)) as M grows without limit, so treating an M* as
 * a Mach number understates the Mach number, and increasingly so the faster the flow.
 *
 * Transonic series solutions (Sauer, Hall, Kliegel-Levine) are all posed in M*, so their
 * output must be converted before it can be used as a Mach number.
 *
 * @param m_star Critical velocity ratio V/a*; must be below sqrt((gamma+1)/(gamma-1)).
 * @param gamma Isentropic exponent.
 * @return The corresponding Mach number.
 */
inline double mach_from_critical_velocity_ratio(double m_star, double gamma) {
    const double denominator = (gamma + 1.0) - (gamma - 1.0) * m_star * m_star;
    if (denominator <= 0.0) {
        // M* has reached its finite ceiling, where M is unbounded. Callers treat a
        // non-finite Mach as a failed initialization rather than propagating it.
        return std::numeric_limits<double>::infinity();
    }
    return std::sqrt(2.0 * m_star * m_star / denominator);
}

/**
 * @brief Critical velocity ratio M* = V/a* from the Mach number. Inverse of
 *        mach_from_critical_velocity_ratio.
 *
 * @param mach Mach number.
 * @param gamma Isentropic exponent.
 * @return The corresponding critical velocity ratio.
 */
inline double critical_velocity_ratio_from_mach(double mach, double gamma) {
    return mach * std::sqrt((gamma + 1.0) / (2.0 + (gamma - 1.0) * mach * mach));
}


/**
 * @brief Speed of sound of a gas.
 */
double gas_sonic_velocity(const Cantera::ThermoPhase& gas, double gamma);
double gas_sonic_velocity(double temperature, double molar_mass, double gamma);
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
 * For a perfect gas, calculates the area ratio for a given Mach number 
 * using the Area-Mach number relation. 
 */
double area_mach_relation(double mach, double gamma);
/**
 * @brief Calculate thrust coefficient.
 */
double C_F(double gamma, double pressure_ratio, double area_ratio);

/**
 * @brief Calculate the characteristic combustion velocity (c*)
 */
double cstar(double gamma, double temperature, double molecular_weight);

}