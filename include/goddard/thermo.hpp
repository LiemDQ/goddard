#pragma once

#include "cantera/core.h"
#include <vector>
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

double molar_mass_from_composition(Cantera::ThermoPhase& thermo, const std::vector<double>& state);

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