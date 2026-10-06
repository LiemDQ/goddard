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
 * @brief Get the density of an ideal gas from its pressure, temperature, and molar mass.
 *
 * @param P pressure in Pa
 * @param T temperature in K
 * @param molar_mass molar mass in kg/kmol
 *
 * @returns Density in kg/m3
 */
inline double ideal_gas_P_to_D(double P, double T, double molar_mass) {
    return P*molar_mass/(T*Cantera::GasConstant);
}

/**
 * Snapshot of the thermodynamic state of a mixture (see `Gas::snapshot`).
 *
 * Specific quantities are per kg of mixture (gas plus condensed phases). `gamma_s`,
 * `dlV_dlP_T`, `dlV_dlT_P` and `speed_of_sound` are in the chemistry the state was taken with:
 * equilibrium or frozen.
 */
class ThermodynamicState {
public:
    /// Pressure [Pa].
    double pressure;
    /// Temperature [K].
    double temperature;
    /// Mixture density [kg/m^3].
    double density;
    /// Specific enthalpy [J/kg].
    double enthalpy;
    /// Specific internal energy [J/kg], h - P/rho.
    double internal_energy;
    /// Specific Gibbs energy [J/kg], h - T s.
    double gibbs;
    /// Specific entropy [J/(kg.K)].
    double entropy;
    /// Molecular weight [kg/kmol], CEA's "M" = 1/n with n the moles of gas per kg of mixture.
    double molecular_weight;
    /**
     * Constant-pressure specific heat at frozen composition [J/(kg.K)], as `Gas::cp_mass`.
     * `Gas::snapshot` sets it to 0 at a pinned phase transition, as CEA prints.
     */
    double cp;
    /// Isentropic exponent -(d log P / d log V)_s [-].
    double gamma_s;
    /// (d log V / d log P)_T [-].
    double dlV_dlP_T;
    /// (d log V / d log T)_P [-].
    double dlV_dlT_P;
    /// Speed of sound [m/s].
    double speed_of_sound;
    /// Stagnation enthalpy [J/kg] the flow velocity is measured from.
    double stagnation_enthalpy;
    /**
     * Mass fractions of the mixture [-] by species name, condensed species included. Species
     * with a zero fraction are omitted.
     */
    Composition composition;
    /**
     * Mixture molecular weight [kg/kmol], CEA's "MW": one kg of mixture divided by the moles of
     * gas *plus* condensed species it holds. Equal to `molecular_weight` (CEA's "M" = 1/n, which
     * counts the gas alone) when no condensed phase is present.
     */
    double mixture_molecular_weight = 0.0;
    /** Mass fraction of the gas phase in the mixture [-]. 1 with no condensed phase present. */
    double gas_mass_fraction = 1.0;
    /**
     * True when the state sits exactly at a condensed phase transition with both polymorphs
     * present and the chemistry is shifting, so the equilibrium specific heat is infinite. CEA
     * prints a specific heat of zero for such a station.
     */
    bool pinned_transition = false;
};


/**
 * Temperature, pressure and composition of a reactant stream.
 */
class PhaseSpecification {
public:
    PhaseSpecification() = default;
    /**
     * @param T Temperature [K].
     * @param P Pressure [Pa].
     * @param comp Mole fractions as a composition string, e.g. "H2:2, O2:1".
     */
    PhaseSpecification(double T, double P, const std::string& comp);
    /**
     * @param T Temperature [K].
     * @param P Pressure [Pa].
     * @param comp Mole fractions by species name.
     */
    PhaseSpecification(double T, double P, const Composition& comp);


    /// Temperature [K].
    double T = 0.0;
    /// Pressure [Pa].
    double P = 0.0;
    /**
     * Mole fractions [-] by species name; normalized when the state is set.
     *
     * In Python, assign a whole dict: item assignment acts on a copy.
     */
    Composition composition = {};
};


} //namespace Goddard