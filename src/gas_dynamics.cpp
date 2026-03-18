#include "goddard/gas_dynamics.hpp"

#include <cmath>
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

Eigen::ArrayXXd gas_sonic_velocity(const ThermoArray& gas, const Eigen::ArrayXXd& gamma) {
    return (Cantera::GasConstant*gas.temperature()*gamma / gas.mean_molecular_weight()).sqrt();
}

double gas_stagnation_enthalpy(const Cantera::ThermoPhase& gas, double velocity) {
    return gas.enthalpy_mass() + velocity*velocity/2;
}

Eigen::ArrayXXd gas_stagnation_enthalpy(const ThermoArray& gas, const Eigen::ArrayXXd velocity) {
    return gas.enthalpy_mass() + velocity * velocity / 2;
}

double stagnation_pressure(const Cantera::ThermoPhase& gas, double mach, double gamma) {
    return gas.pressure() * std::pow((1 + (gamma-1)/2 * mach * mach), gamma/(gamma - 1));
}

Eigen::ArrayXXd stagnation_pressure(const Cantera::ThermoPhase& gas, const Eigen::ArrayXXd mach, const Eigen::ArrayXXd gamma) {
    return gas.pressure() * (1 + (gamma-1)/2 * mach * mach).pow(gamma/(gamma - 1));
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

double C_F(double gamma, double pressure_ratio, double) {
    double term1 = 2*gamma*gamma/(gamma-1);
    double term2 = std::pow(2/(gamma+1), (gamma+1)/(gamma-1));
    double term3 = 1 - std::pow(pressure_ratio, (gamma-1)/gamma);
    //TODO: do we need to include pressure term?

    return std::sqrt(term1*term2*term3);
}

}