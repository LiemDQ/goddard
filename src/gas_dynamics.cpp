#include <cmath>
#include "goddard/gas_dynamics.hpp"
#include "goddard/utils.hpp"
#include "goddard/error.hpp"

namespace Goddard {

double gas_isenthalpic_velocity(const Cantera::ThermoPhase& gas, double H_stagnation){
    return std::sqrt(2*(H_stagnation - gas.enthalpy_mass()));
}

Eigen::ArrayXXd gas_isenthalpic_velocity(const ThermoArray& gas, const Eigen::ArrayXXd& H_stagnation) {
    return (2*(H_stagnation - gas.enthalpy_mass())).sqrt();
}

double gas_sonic_velocity(const Cantera::ThermoPhase& gas, double gamma){
    return std::sqrt(Cantera::GasConstant*gas.temperature()*gamma/gas.meanMolecularWeight());
}

double gas_sonic_velocity(double temperature, double molar_mass, double gamma) {
    return std::sqrt(Cantera::GasConstant*temperature*gamma/molar_mass);
}

Eigen::ArrayXXd gas_sonic_velocity(const ThermoArray& gas, const Eigen::ArrayXXd& gamma) {
    return (Cantera::GasConstant*gas.temperature()*gamma / gas.mean_molecular_weight()).sqrt();
}

double gas_stagnation_enthalpy(const Cantera::ThermoPhase& gas, double velocity) {
    return gas.enthalpy_mass() + velocity*velocity/2;
}

Eigen::ArrayXXd gas_stagnation_enthalpy(const ThermoArray& gas, const Eigen::ArrayXXd& velocity) {
    return gas.enthalpy_mass() + velocity * velocity / 2;
}

double perfect_gas_stagnation_pressure(double P, double mach, double gamma) {
    return P * std::pow((1 + (gamma-1)/2 * mach * mach), gamma/(gamma - 1));
}

Eigen::ArrayXXd perfect_gas_stagnation_pressure(const Eigen::ArrayXXd& P, const Eigen::ArrayXXd& mach, const Eigen::ArrayXXd& gamma) {
    return P * (1 + (gamma-1)/2 * mach * mach).pow(gamma/(gamma - 1));
}

double gas_stagnation_pressure(const Cantera::ThermoPhase& gas, double velocity) {
    auto thermo = gas.clone();
    const double h_stag = gas_stagnation_enthalpy(gas, velocity);
    const double entropy = gas.entropy_mass();
    
    // initial guess
    double gamma = gas.cp_mass()/gas.cv_mass();
    double mach = velocity / gas_sonic_velocity(*thermo, gamma);
    double P_stag = perfect_gas_stagnation_pressure(thermo->pressure(), mach, gamma);
    
    int max_iters = 10;
    int k = 0;
    const double abstol = 1e-8;
    double residual = 1.0;
    // Use Newton's method to solve for stagnation pressure.
    // As energy is conserved, h_stag - h(S, P_stag) = 0.
    // Thus: P_{k+1} = P_k - (h_stag - h(S, P_stag))/(dh/dP)_S
    // Note that by definition, (dh/dP)_S = V = 1/rho
    while (abs(residual) > abstol) {
        if (k > max_iters)
            throw ConvergenceError("Failed to converge to stagnation pressure.", k, abstol, residual);
        thermo->setState_SP(entropy, P_stag);
        double rho = thermo->density();
        residual = thermo->enthalpy_mass() - h_stag;
        P_stag = P_stag - residual * rho;
        k++;
    }

    return P_stag;
}

double stagnation_factor(double mach, double gamma) {
    return 1.0 + (gamma - 1.0)/2.0 * mach * mach;
}

Eigen::ArrayXXd stagnation_factor(const Eigen::ArrayXXd& mach, const Eigen::ArrayXXd& gamma) {
    return 1.0 + (gamma - 1.0)/2.0 * mach * mach;
}


double area_per_mdot(const Cantera::ThermoPhase& gas, double velocity) {
    return gas.temperature()*Cantera::GasConstant / (gas.pressure() * velocity * gas.meanMolecularWeight());
}

Eigen::ArrayXXd area_per_mdot(const ThermoArray& gas, const Eigen::ArrayXXd& velocity) {
    return gas.temperature() * Cantera::GasConstant / (gas.pressure() * velocity * gas.mean_molecular_weight());
}

double cstar(double gamma, double temperature, double molecular_weight) {
    double gamma_term = std::pow(2.0/(gamma + 1.0), -(gamma + 1.0)/(2*(gamma - 1.0)));
    return std::sqrt(Cantera::GasConstant * temperature / (molecular_weight* gamma)) * gamma_term;
}

double isp(const Cantera::ThermoPhase& gas, double gamma, double enthalpy) {
    double velocity = gas_isenthalpic_velocity(gas, enthalpy);
    double cs = cstar(gamma, gas.temperature(), gas.meanMolecularWeight());
    return velocity/cs;
}

double ivac(const Cantera::ThermoPhase& gas, double gamma, double enthalpy) {
    double v = isp(gas, gamma, enthalpy);
    return v + gas.temperature()*Cantera::GasConstant/(v * gas.meanMolecularWeight());
}

double mach(const Cantera::ThermoPhase& gas, double H_stag, double gamma) {
    double velocity = gas_isenthalpic_velocity(gas, H_stag);
    double sonic = gas_sonic_velocity(gas, gamma);
    
    return velocity/sonic;
}

double mach_area_relation(double mach, double gamma) {
    return std::sqrt(
        (1/(mach*mach))
        *std::pow(2.0/(gamma+1)*stagnation_factor(mach, gamma), (gamma+1)/(gamma-1)));
}

double C_F(double gamma, double pressure_ratio, double) {
    double term1 = 2*gamma*gamma/(gamma-1);
    double term2 = std::pow(2/(gamma+1), (gamma+1)/(gamma-1));
    double term3 = 1 - std::pow(pressure_ratio, (gamma-1)/gamma);
    //TODO: do we need to include pressure term?

    return std::sqrt(term1*term2*term3);
}

}