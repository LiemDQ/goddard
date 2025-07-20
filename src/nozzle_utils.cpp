#include "goddard/nozzle_utils.hpp"

namespace Goddard {

double gas_isenthalpic_velocity(const Cantera::ThermoPhase& gas, double H_stagnation){
    return std::sqrt(2*(H_stagnation - gas.enthalpy_mass()));
}

double gas_sonic_velocity(const Cantera::ThermoPhase& gas, double gamma){
    return std::sqrt(Cantera::GasConstant*gas.temperature()*gamma/gas.meanMolecularWeight());
}

double area_per_mdot(const Cantera::ThermoPhase& gas, double velocity) {
    return gas.temperature()*Cantera::GasConstant / (gas.pressure() * velocity * gas.meanMolecularWeight());
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

}