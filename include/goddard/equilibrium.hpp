#pragma once
#include <vector>

#include "cantera/core.h"
#include "eigen3/Eigen/Dense"

namespace Goddard {

/**
 * Forward declaration. `equilibrium.hpp` is included by `gas.hpp` for the result structs below,
 * so it must not include `gas.hpp`. The `const Gas&` overloads are defined in `equilibrium.cpp`.
 */
class Gas;

struct EquilibriumDerivatives {
    /** (d pi_i / d log T)_P [-], one entry per element of the gas phase. */
    Eigen::ArrayXd dpi_dlogT_P;
    /** (d log n / d log T)_P [-], n being the moles of gas per kg of mixture. */
    double dlogn_dlogT_P = 0.0;
    /** (d pi_i / d log P)_T [-], one entry per element of the gas phase. */
    Eigen::ArrayXd dpi_dlogP_T;
    /** (d log n / d log P)_T [-]. */
    double dlogn_dlogP_T = 0.0;
    /**
     * (d n_k / d log T)_P [kmol per kg of mixture], one entry per condensed species that is
     * currently present, in `Gas::condensed_species_names()` order restricted to those species.
     * Empty for a gas-only mixture.
     */
    Eigen::ArrayXd dn_condensed_dlogT_P;
    /**
     * (d n_k / d log P)_T [kmol per kg of mixture], same ordering as `dn_condensed_dlogT_P`.
     *
     * At a pinned phase transition only the sum over the two coexisting polymorphs is determined:
     * it is reported on the lower-temperature polymorph and the higher-temperature one is 0.
     */
    Eigen::ArrayXd dn_condensed_dlogP_T;
    /**
     * True when the mixture sits exactly at a condensed phase transition with both polymorphs
     * present (see `Gas::at_phase_transition()`). The temperature derivatives `dpi_dlogT_P`,
     * `dlogn_dlogT_P` and `dn_condensed_dlogT_P` are then NaN.
     */
    bool pinned_transition = false;
};

struct ExpansionProperties {
    /** (d log V / d log T)_P [-]. */
    double dlogV_dlogT_P = 0.0;
    /** (d log V / d log P)_T [-]. */
    double dlogV_dlogP_T = 0.0;
    /** Equilibrium constant-pressure specific heat [J/(kg.K)]. */
    double spec_heat_p = 0.0;
    /** Isentropic exponent -(d log P / d log V)_s [-]. */
    double gamma_s = 0.0;
    /** Equilibrium constant-volume specific heat [J/(kg.K)]. */
    double spec_heat_v = 0.0;
    /** Moles of gas per kg of mixture [kmol/kg] (CEA's 1/M). */
    double gas_moles = 0.0;
    /** Moles of gas plus condensed species per kg of mixture [kmol/kg] (CEA's 1/MW). */
    double total_moles = 0.0;
    /** Mixture density [kg/m^3], P / (gas_moles * R * T). */
    double density = 0.0;
    /** Equilibrium speed of sound [m/s], sqrt(gas_moles * R * T * gamma_s). */
    double speed_of_sound = 0.0;
    /** Frozen-composition constant-pressure specific heat [J/(kg.K)]. */
    double frozen_spec_heat_p = 0.0;
    /** Frozen-composition ratio of specific heats cp/cv [-]. */
    double frozen_gamma = 0.0;
    /** True when the mixture sits exactly at a condensed phase transition. */
    bool pinned_transition = false;
};

/**
 * @name Gas-only helpers
 *
 * These read the gas phase alone. Amounts are per kg of *gas*, so a mixture carrying condensed
 * products must scale them by `Gas::gas_mass_fraction()` to reach a per kg of mixture basis.
 * @{
 */
/**
 * @brief Get matrix of stoichiometric coefficients of the species contained in the `ThermoPhase` object.
 *
 * @return 2D Eigen array of stoichiometric coefficients. Rows represent species, while columns represent elements.
 * Ordering is the same as the data file used to generated the `ThermoPhase` object.
 */
Eigen::ArrayXXd get_stoichiometric_coeffs(const Cantera::ThermoPhase& gas);
/** Amount of each gas species [kmol per kg of gas]. */
Eigen::ArrayXd get_mole_vector(const Cantera::ThermoPhase& gas);
/** Standard-state molar enthalpy of each gas species divided by R*T [-]. */
Eigen::ArrayXd get_enthalpyRT_vector(const Cantera::ThermoPhase& gas);
/** Standard-state molar heat capacity of each gas species divided by R [-]. */
Eigen::ArrayXd get_cpR_vector(const Cantera::ThermoPhase& gas);

EquilibriumDerivatives get_thermo_equilibrium_derivatives(const Cantera::ThermoPhase& gas);
ExpansionProperties get_thermo_equilibrium_properties(const Cantera::ThermoPhase& gas, const EquilibriumDerivatives& derivatives);
ExpansionProperties get_thermo_equilibrium_properties(const Cantera::ThermoPhase& gas);
/** Isentropic exponent -(d log P / d log V)_s [-] of the gas phase at equilibrium. */
double get_equilibrium_gamma(const Cantera::ThermoPhase& gas);
/** @} */

/**
 * @name Mixture-aware overloads
 *
 * These take a `Gas`, so they account for condensed products: amounts are per kg of mixture and
 * the Gordon & McBride system carries one extra unknown and one extra row per condensed species
 * that is present. For a `Gas` without condensed phases they reduce exactly to the
 * `Cantera::ThermoPhase` overloads above.
 *
 * At a pinned phase transition (`Gas::at_phase_transition()` with both polymorphs present) the
 * two polymorphs share one element row, so they are merged into a single condensed unknown and
 * only the pressure block is solved. The temperature derivatives are then undefined and are
 * reported as NaN, `spec_heat_p`, `spec_heat_v` and `dlogV_dlogT_P` are infinite, and
 * `gamma_s = -1 / dlogV_dlogP_T`.
 * @{
 */
EquilibriumDerivatives get_thermo_equilibrium_derivatives(const Gas& gas);
ExpansionProperties get_thermo_equilibrium_properties(const Gas& gas, const EquilibriumDerivatives& derivatives);
ExpansionProperties get_thermo_equilibrium_properties(const Gas& gas);
/**
 * Frozen-composition expansion properties: gas composition and condensed amounts are held fixed,
 * so the volume derivatives are +/-1 and the specific heats are the mixture's frozen ones.
 * Condensed species are incompressible and so contribute equally to cp and cv.
 */
ExpansionProperties get_frozen_properties(const Gas& gas);
/** Isentropic exponent -(d log P / d log V)_s [-] of the mixture at equilibrium. */
double get_equilibrium_gamma(const Gas& gas);
/** @} */

} //namespace Goddard