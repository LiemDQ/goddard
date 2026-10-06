#include <cmath>
#include <format>
#include <stdexcept>
#include "goddard/gas_dynamics.hpp"
#include "goddard/error.hpp"

namespace Goddard {

double gas_sonic_velocity(const Cantera::ThermoPhase& gas, double gamma){
    return std::sqrt(Cantera::GasConstant*gas.temperature()*gamma/gas.meanMolecularWeight());
}

double gas_sonic_velocity(double temperature, double molar_mass, double gamma) {
    return std::sqrt(Cantera::GasConstant*temperature*gamma/molar_mass);
}

double perfect_gas_stagnation_pressure(double P, double mach, double gamma) {
    return P * std::pow((1 + (gamma-1)/2 * mach * mach), gamma/(gamma - 1));
}

double stagnation_factor(double mach, double gamma) {
    return 1.0 + (gamma - 1.0)/2.0 * mach * mach;
}

double cstar(double gamma, double temperature, double molecular_weight) {
    double gamma_term = std::pow(2.0/(gamma + 1.0), -(gamma + 1.0)/(2*(gamma - 1.0)));
    return std::sqrt(Cantera::GasConstant * temperature / (molecular_weight* gamma)) * gamma_term;
}

double area_mach_relation(double mach, double gamma) {
    return std::sqrt(
        (1/(mach*mach))
        *std::pow(2.0/(gamma+1)*stagnation_factor(mach, gamma), (gamma+1)/(gamma-1)));
}

double mach_from_area_ratio(double area_ratio, double gamma, bool supersonic) {
    if (!(area_ratio >= 1.0)) {
        throw std::invalid_argument(std::format(
            "mach_from_area_ratio: area ratio must be at least 1. Actual value: {}", area_ratio));
    }
    if (area_ratio == 1.0) {
        return 1.0;
    }
    // A/A* decreases monotonically on (0, 1] and increases monotonically on [1, inf).
    double lower = 0.0;
    double upper = 1.0;
    if (supersonic) {
        lower = 1.0;
        upper = 2.0;
        while (area_mach_relation(upper, gamma) < area_ratio) {
            lower = upper;
            upper *= 2.0;
            if (upper > 1e8) {
                throw std::invalid_argument(std::format(
                    "mach_from_area_ratio: no supersonic root for area ratio {}", area_ratio));
            }
        }
    }
    const int max_iters = 200;
    for (int iter = 0; iter < max_iters; iter++) {
        const double mid = 0.5 * (lower + upper);
        if (mid <= lower || mid >= upper) {
            break;
        }
        const double excess = area_mach_relation(mid, gamma) - area_ratio;
        // Subsonic: A/A* too large means M too small. Supersonic: A/A* too large means M too large.
        const bool mach_too_small = supersonic ? (excess < 0.0) : (excess > 0.0);
        if (mach_too_small) {
            lower = mid;
        } else {
            upper = mid;
        }
    }
    return 0.5 * (lower + upper);
}

double finite_area_pressure_loss(double mach, double gamma) {
    const double psi = stagnation_factor(mach, gamma);
    return (1.0 + gamma * mach * mach) / std::pow(psi, gamma / (gamma - 1.0));
}

}