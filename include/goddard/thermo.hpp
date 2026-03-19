#pragma once

#include "cantera/core.h"
#include <vector>

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
 * Convenience class for containing thermodynamic data in one place.
 */
class ThermodynamicState {
public:
    ThermodynamicState() = default;
    ThermodynamicState(double temperature, double pressure, const std::string& comp)
        : T(temperature), P(pressure), composition(comp) {}

    double T = 0.0;
    double P = 0.0;
    std::string composition = {};

    std::vector<double> to_vector(Cantera::ThermoPhase& sln);
    std::vector<double> to_mass_vector(Cantera::ThermoPhase& sln);
    
};

/**
 * Convenience class for querying thermodynamic data.
 */
class ThermoData {
    public:
        ThermoData(std::shared_ptr<Cantera::Solution> sln);
        ThermoData(std::shared_ptr<Cantera::ThermoPhase> thermo);

        double molar_mass_from_composition(const std::vector<double>& state);
    private:
        std::shared_ptr<Cantera::ThermoPhase> m_thermo;
        std::vector<double> m_state_vector;
};

} //namespace Goddard