#pragma once

#include "cantera/core.h"
#include <vector>
#include <map>
#include <string_view>
#include "goddard/chemistry.hpp"

namespace Goddard {

/**
 * @brief Get the pressure of an ideal gas from its density, temperature, and molar mass. 
 * 
 * @param D density in kg/m3
 * @param T temperature in K
 * @param molar_mass molar mass in kg/kmol
 * 
 * @returns Pressure in Pa 
 */
inline double ideal_gas_D_to_P(double D, double T, double molar_mass) {
    return D*T*Cantera::GasConstant/molar_mass;
}

/**
 * @brief Get the pressure of an ideal gas from its density, temperature, and molar mass. 
 * 
 * @param P pressure in Pa
 * @param T temperature in K
 * @param molar_mass molar mass in kg/kmol
 * 
 * @returns Density in kg/kmol 
 */
inline double ideal_gas_P_to_D(double P, double T, double molar_mass) {
    return P*molar_mass/(T*Cantera::GasConstant);
}

/**
 * Convenience class containing relevant thermodynamic results
 */
class ThermodynamicState {
public:
    double pressure;
    double temperature;
    double density;
    double enthalpy;
    double internal_energy;
    double gibbs;
    double entropy;
    double molecular_weight;
    double cp;
    double gamma_s;
    double dlV_dlP_T;
    double dlV_dlT_P;
    double speed_of_sound;
    double stagnation_enthalpy;
    Composition composition; // stored as mass fractions

    size_t state_size() const;
    /**
     * Outputs a raw vector suitable for use with Cantera objects. 
     */
    std::vector<double> to_vector() const;
    std::vector<double> to_mole_vector(Cantera::ThermoPhase& sln) const;
};


/**
 * Convenience class for initializing a thermodynamic state.
 */
class PhaseSpecification {
public:
    PhaseSpecification() = default;
    PhaseSpecification(double T, double P, const std::string& comp);
    PhaseSpecification(double T, double P, const Composition& comp);


    double T = 0.0;
    double P = 0.0;
    Composition composition = {};

    std::vector<double> to_vector(Cantera::ThermoPhase& sln) const;
    std::vector<double> to_mass_vector(Cantera::ThermoPhase& sln) const;
    
};


} //namespace Goddard