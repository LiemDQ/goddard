#include "cantera/core.h"
#include "eigen3/Eigen/Dense"

#include "goddard/equilibrium.hpp"
#include "goddard/error.hpp"
#include "goddard/gas.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;

namespace Goddard {


ArrayXXd get_stoichiometric_coeffs(const Cantera::ThermoPhase& gas){
    size_t n_elements = gas.nElements();
    size_t n_species = gas.nSpecies();

    //construct matrix of elemental stoichiometric coefficients
    ArrayXXd stoich_coeffs(n_species, n_elements);
    for(size_t i = 0; i < n_elements; i++){
        for(size_t j = 0; j < n_species; j++){
            stoich_coeffs(j,i) = gas.nAtoms(j, i);
        }
    }

    return stoich_coeffs;
}

ArrayXd get_mole_vector(const Cantera::ThermoPhase& gas){
    size_t n_species = gas.nSpecies();
    double total_moles = 1.0/gas.meanMolecularWeight();

    ArrayXd moles(n_species);
    gas.getMoleFractions(moles.data());
    moles *= total_moles;
    return moles;
}

ArrayXd get_enthalpyRT_vector(const Cantera::ThermoPhase& gas){
    size_t n_species = gas.nSpecies();

    ArrayXd std_enthalpies_RT(n_species);
    
    gas.getEnthalpy_RT(std_enthalpies_RT.data());

    return std_enthalpies_RT;
}

ArrayXd get_cpR_vector(const Cantera::ThermoPhase& gas){
    size_t n_species = gas.nSpecies();
    ArrayXd cp_R(n_species);
    
    gas.getCp_R(cp_R.data());
    return cp_R;
}

namespace {

/**
 * @brief Everything the Gordon & McBride derivative system needs, on a per kg of mixture basis.
 *
 * Built either from a bare `Cantera::ThermoPhase` (gas-only: the condensed arrays are empty and
 * the gas mass fraction is 1) or from a `Gas` that carries condensed products. Only condensed
 * species that are actually present (nonzero moles) appear, compacted in candidate order, which
 * is the ordering of `EquilibriumDerivatives::dn_condensed_*`.
 */
struct MixtureState {
    /** Atoms of element i in gas species j [-]: n_species x l, species-major. */
    ArrayXXd gas_stoich;
    /** Amount of each gas species [kmol per kg of mixture]. */
    ArrayXd gas_moles;
    /** Standard-state molar enthalpy of each gas species divided by R*T [-]. */
    ArrayXd gas_H_RT;
    /** Standard-state molar heat capacity of each gas species divided by R [-]. */
    ArrayXd gas_cp_R;
    /** Atoms of element i in condensed species c [-]: |C| x l, present species only. */
    ArrayXXd condensed_stoich;
    /** Amount of each present condensed species [kmol per kg of mixture]. */
    ArrayXd condensed_moles;
    /** Reference-state molar enthalpy of each present condensed species divided by R*T [-]. */
    ArrayXd condensed_H_RT;
    /** Reference-state molar heat capacity of each present condensed species divided by R [-]. */
    ArrayXd condensed_cp_R;
    /** Constant-pressure specific heat of the gas phase [J/(kg.K)], per kg of *gas*. */
    double gas_cp_mass = 0.0;
    /** Constant-volume specific heat of the gas phase [J/(kg.K)], per kg of *gas*. */
    double gas_cv_mass = 0.0;
    /** Mass fraction of the gas phase in the mixture [-]. */
    double gas_mass_fraction = 1.0;
    /** Temperature [K]. */
    double temperature = 0.0;
    /** Pressure [Pa]. */
    double pressure = 0.0;
    /** Index into the present condensed species of the lower-temperature pinned polymorph, or -1. */
    long pinned_low = -1;
    /** Index into the present condensed species of the higher-temperature pinned polymorph, or -1. */
    long pinned_high = -1;

    /** True when two polymorphs coexist at a phase-transition temperature. */
    bool is_pinned() const { return pinned_low >= 0 && pinned_high >= 0; }
};

/** Collect the gas-phase part of a mixture state. Moles are per kg of gas at this point. */
MixtureState gas_only_state(const Cantera::ThermoPhase& gas) {
    MixtureState state;
    state.gas_stoich = get_stoichiometric_coeffs(gas);
    state.gas_moles = get_mole_vector(gas);
    state.gas_H_RT = get_enthalpyRT_vector(gas);
    state.gas_cp_R = get_cpR_vector(gas);
    state.condensed_stoich = ArrayXXd(0, static_cast<long>(gas.nElements()));
    state.condensed_moles = ArrayXd(0);
    state.condensed_H_RT = ArrayXd(0);
    state.condensed_cp_R = ArrayXd(0);
    state.gas_cp_mass = gas.cp_mass();
    state.gas_cv_mass = gas.cv_mass();
    state.temperature = gas.temperature();
    state.pressure = gas.pressure();
    return state;
}

/** Collect the gas and condensed parts of a mixture state, both per kg of mixture. */
MixtureState mixture_state(const Gas& gas) {
    MixtureState state = gas_only_state(*gas.thermo());

    // With no condensed species present the gas mass fraction is exactly 1, so the scaling below
    // leaves the gas-only arithmetic untouched.
    state.gas_mass_fraction = gas.gas_mass_fraction();
    state.gas_moles *= state.gas_mass_fraction;

    if (!gas.has_condensed_phases()) {
        return state;
    }

    const std::vector<double> all_moles = gas.condensed_moles();
    std::vector<double> present;
    for (double moles : all_moles) {
        if (moles > 0.0) present.push_back(moles);
    }
    state.condensed_moles = Eigen::Map<const ArrayXd>(present.data(), static_cast<long>(present.size()));
    state.condensed_stoich = gas.condensed_stoich_coeffs();
    state.condensed_H_RT = gas.condensed_enthalpy_RT();
    state.condensed_cp_R = gas.condensed_cp_R();

    if (gas.at_phase_transition()) {
        const std::pair<long, long> polymorphs = gas.pinned_polymorphs();
        // A group is only genuinely pinned when both of its polymorphs are present; with a single
        // polymorph the system is regular and is solved as usual.
        if (polymorphs.first >= 0 && polymorphs.second >= 0) {
            state.pinned_low = polymorphs.first;
            state.pinned_high = polymorphs.second;
        }
    }
    return state;
}

/**
 * @brief Assemble and solve the Gordon & McBride derivative system.
 *
 * Unknowns are ordered `[pi (l) | dn_C' (|C'|) | dlog n (1)]` and the same matrix is solved for
 * the temperature and the pressure right-hand side (G&M 2.56-2.58 and 2.64-2.66, extended with
 * the condensed columns and rows as in `instructions/condensed_species.md` section 4):
 *
 * - element row k: `sum_i [sum_G a_kj a_ij n_j] pi_i + sum_C a_kj dn_j + [sum_G a_kj n_j] dlog n`
 * - condensed row j: `sum_i a_ij pi_i`
 * - mole row: `sum_i [sum_G a_ij n_j] pi_i` (the dlog n coefficient `sum_G n_j - n` vanishes)
 *
 * At a pinned phase transition the two coexisting polymorphs have identical element rows, so the
 * system is singular. They are then merged into a single condensed unknown (|C'| = |C| - 1) and
 * only the pressure block is solved; the temperature block is undefined there and is reported as
 * NaN.
 */
EquilibriumDerivatives solve_derivative_system(const MixtureState& state) {
    const long n_elements = state.gas_stoich.cols();
    const long n_condensed = state.condensed_moles.size();
    const bool pinned = state.is_pinned();

    // Column of each present condensed species. A pinned polymorph pair shares one column, so the
    // solved unknown is d(n_low + n_high) rather than the two amounts separately.
    std::vector<long> column_of(static_cast<size_t>(n_condensed), 0);
    std::vector<long> representative; // condensed column -> present species providing its row
    long n_columns = 0;
    for (long c = 0; c < n_condensed; c++) {
        if (pinned && c == state.pinned_high) continue;
        column_of[static_cast<size_t>(c)] = n_columns++;
        representative.push_back(c);
    }
    if (pinned) {
        column_of[static_cast<size_t>(state.pinned_high)] =
            column_of[static_cast<size_t>(state.pinned_low)];
    }

    const long n_unknowns = n_elements + n_columns + 1;
    const long mole_row = n_elements + n_columns;

    MatrixXd coeff_matrix = MatrixXd::Zero(n_unknowns, n_unknowns);
    // Column 0 is the temperature right-hand side, column 1 the pressure one.
    MatrixXd rhs = MatrixXd::Zero(n_unknowns, 2);

    // Moles of each element carried by the gas phase [kmol per kg of mixture].
    ArrayXd element_gas_moles(n_elements);
    for (long i = 0; i < n_elements; i++) {
        element_gas_moles(i) = (state.gas_stoich.col(i) * state.gas_moles).sum();
    }

    // Element rows (G&M 2.56 and 2.64).
    for (long k = 0; k < n_elements; k++) {
        for (long i = 0; i < n_elements; i++) {
            coeff_matrix(k, i) = (state.gas_stoich.col(k) * state.gas_stoich.col(i) * state.gas_moles).sum();
        }
        for (long c = 0; c < n_condensed; c++) {
            // Merged polymorphs write the same coefficient to the same column: their stoichiometry
            // is identical, which is exactly why the unmerged system is singular.
            coeff_matrix(k, n_elements + column_of[static_cast<size_t>(c)]) = state.condensed_stoich(c, k);
        }
        coeff_matrix(k, mole_row) = element_gas_moles(k);
        rhs(k, 0) = -(state.gas_stoich.col(k) * state.gas_moles * state.gas_H_RT).sum();
        rhs(k, 1) = element_gas_moles(k);
    }

    // Condensed rows (G&M 2.57 and 2.65): the chemical potential of a pure condensed phase has no
    // pressure dependence, hence the zero pressure right-hand side.
    for (long col = 0; col < n_columns; col++) {
        const long c = representative[static_cast<size_t>(col)];
        for (long i = 0; i < n_elements; i++) {
            coeff_matrix(n_elements + col, i) = state.condensed_stoich(c, i);
        }
        rhs(n_elements + col, 0) = -state.condensed_H_RT(c);
    }

    // Mole row (G&M 2.58 and 2.66).
    for (long i = 0; i < n_elements; i++) {
        coeff_matrix(mole_row, i) = element_gas_moles(i);
    }
    rhs(mole_row, 0) = -(state.gas_moles * state.gas_H_RT).sum();
    rhs(mole_row, 1) = state.gas_moles.sum();

    EquilibriumDerivatives result;
    result.pinned_transition = pinned;

    const auto decomposition = coeff_matrix.colPivHouseholderQr();

    if (pinned) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        result.dpi_dlogT_P = ArrayXd::Constant(n_elements, nan);
        result.dlogn_dlogT_P = nan;
        result.dn_condensed_dlogT_P = ArrayXd::Constant(n_condensed, nan);
    } else {
        const VectorXd solution = decomposition.solve(rhs.col(0));
        result.dpi_dlogT_P = solution.head(n_elements).array();
        result.dlogn_dlogT_P = solution(mole_row);
        result.dn_condensed_dlogT_P = ArrayXd(n_condensed);
        for (long c = 0; c < n_condensed; c++) {
            result.dn_condensed_dlogT_P(c) = solution(n_elements + column_of[static_cast<size_t>(c)]);
        }
    }

    const VectorXd solution = decomposition.solve(rhs.col(1));
    result.dpi_dlogP_T = solution.head(n_elements).array();
    result.dlogn_dlogP_T = solution(mole_row);
    result.dn_condensed_dlogP_T = ArrayXd(n_condensed);
    for (long c = 0; c < n_condensed; c++) {
        // The merged derivative of a pinned pair is reported on the lower-temperature polymorph;
        // the higher-temperature one gets zero, because only the sum is determined.
        if (pinned && c == state.pinned_high) {
            result.dn_condensed_dlogP_T(c) = 0.0;
        } else {
            result.dn_condensed_dlogP_T(c) = solution(n_elements + column_of[static_cast<size_t>(c)]);
        }
    }

    return result;
}

/** Frozen-composition specific heats of the mixture [J/(kg.K)]: {cp, cv}. */
std::pair<double, double> frozen_specific_heats(const MixtureState& state) {
    // Condensed species are incompressible, so they contribute the same amount to cp and cv. The
    // gas-phase values are per kg of gas, hence the w_g factor. This is equivalent to
    // cp_f = R sum_{G+C} n_j cp_j/R and cv_f = cp_f - n R, but reproduces Cantera's gas-only
    // values bit for bit.
    const double condensed_cp = Cantera::GasConstant * (state.condensed_moles * state.condensed_cp_R).sum();
    return {state.gas_mass_fraction * state.gas_cp_mass + condensed_cp,
            state.gas_mass_fraction * state.gas_cv_mass + condensed_cp};
}

/**
 * @brief Mixture properties from the equilibrium derivatives (proposal section 5).
 *
 * At a pinned phase transition an isothermal change absorbs the latent heat of the transition, so
 * `spec_heat_p`, `spec_heat_v` and `(d log V / d log T)_P` are infinite and gamma_s reduces to
 * `-1 / (d log V / d log P)_T`.
 */
ExpansionProperties assemble_properties(const MixtureState& state,
                                        const EquilibriumDerivatives& derivs) {
    const double R = Cantera::GasConstant;
    const long n_elements = state.gas_stoich.cols();
    const double gas_moles = state.gas_moles.sum();
    const double total_moles = gas_moles + state.condensed_moles.sum();

    ExpansionProperties props;
    props.gas_moles = gas_moles;
    props.total_moles = total_moles;
    props.density = state.pressure / (gas_moles * R * state.temperature);
    props.pinned_transition = derivs.pinned_transition;

    const std::pair<double, double> frozen = frozen_specific_heats(state);
    props.frozen_spec_heat_p = frozen.first;
    props.frozen_gamma = frozen.first / frozen.second;

    // Condensed volume is neglected, so P V = n R T with n the *gas* moles per kg of mixture.
    props.dlogV_dlogP_T = -1.0 + derivs.dlogn_dlogP_T;

    if (derivs.pinned_transition) {
        const double infinity = std::numeric_limits<double>::infinity();
        props.dlogV_dlogT_P = infinity;
        props.spec_heat_p = infinity;
        props.spec_heat_v = infinity;
        // An isentropic change at a pinned transition happens at constant temperature.
        props.gamma_s = -1.0 / props.dlogV_dlogP_T;
        props.speed_of_sound = std::sqrt(gas_moles * R * state.temperature * props.gamma_s);
        return props;
    }

    props.dlogV_dlogT_P = 1.0 + derivs.dlogn_dlogT_P;

    // Equilibrium isobaric specific heat. NOTE: this is on a per kg of mixture basis because of the
    // normalization of the mole vectors.
    double dpi_cp_contribution = 0.0;
    for (long i = 0; i < n_elements; i++) {
        dpi_cp_contribution += derivs.dpi_dlogT_P(i)
            * (state.gas_moles * state.gas_H_RT * state.gas_stoich.col(i)).sum();
    }
    props.spec_heat_p = R * (
        dpi_cp_contribution
        + (state.gas_moles * state.gas_H_RT).sum() * derivs.dlogn_dlogT_P
        + (state.gas_moles * state.gas_cp_R).sum()
        + (state.gas_moles * state.gas_H_RT * state.gas_H_RT).sum()
        + (state.condensed_H_RT * derivs.dn_condensed_dlogT_P).sum()
        + (state.condensed_moles * state.condensed_cp_R).sum()
    );

    props.spec_heat_v = props.spec_heat_p
        + gas_moles * R * props.dlogV_dlogT_P * props.dlogV_dlogT_P / props.dlogV_dlogP_T;
    props.gamma_s = -(props.spec_heat_p / props.spec_heat_v) / props.dlogV_dlogP_T;
    props.speed_of_sound = std::sqrt(gas_moles * R * state.temperature * props.gamma_s);

    return props;
}

} // namespace

/**
 * @brief Calculate thermodynamic derivatives for a reactive gas at chemical equilibrium.
 *
 * Cantera calculates thermodynamic derivatives assuming a fixed composition,
 * which is not correct for reactive flows, such as ones in chemical equilibrium, as the
 * change in molar quantity and chemical heat release/absorption must also be taken into account.
 * See [Gordon & McBride, 1994, "Computer program for calculation of complex chemical equilibrium
 * compositions and applications. Part 1: Analysis"](https://ntrs.nasa.gov/citations/19950013764) for more details.
 *
 * @note Gas-only: condensed products are invisible to this overload. Use the `Gas` overload for a
 * mixture that carries them.
*/
EquilibriumDerivatives get_thermo_equilibrium_derivatives(const Cantera::ThermoPhase& gas) {
    return solve_derivative_system(gas_only_state(gas));
}

/**
 * @brief Calculate thermodynamic properties of a reacting gas at equilibrium: volume derivatives, heat capacity and adiabatic index.
*/
ExpansionProperties get_thermo_equilibrium_properties(const Cantera::ThermoPhase& gas, const EquilibriumDerivatives& derivs){
    return assemble_properties(gas_only_state(gas), derivs);
}

ExpansionProperties get_thermo_equilibrium_properties(const Cantera::ThermoPhase& gas) {
    return get_thermo_equilibrium_properties(gas, get_thermo_equilibrium_derivatives(gas));
}

double get_equilibrium_gamma(const Cantera::ThermoPhase& gas) {
    return get_thermo_equilibrium_properties(gas).gamma_s;
}

// ---- Mixture-aware overloads ----

EquilibriumDerivatives get_thermo_equilibrium_derivatives(const Gas& gas) {
    return solve_derivative_system(mixture_state(gas));
}

ExpansionProperties get_thermo_equilibrium_properties(const Gas& gas, const EquilibriumDerivatives& derivs) {
    return assemble_properties(mixture_state(gas), derivs);
}

ExpansionProperties get_thermo_equilibrium_properties(const Gas& gas) {
    const MixtureState state = mixture_state(gas);
    return assemble_properties(state, solve_derivative_system(state));
}

ExpansionProperties get_frozen_properties(const Gas& gas) {
    const MixtureState state = mixture_state(gas);
    const double R = Cantera::GasConstant;
    const double gas_moles = state.gas_moles.sum();

    const std::pair<double, double> frozen = frozen_specific_heats(state);

    ExpansionProperties props;
    props.dlogV_dlogT_P = 1.0;
    props.dlogV_dlogP_T = -1.0;
    props.spec_heat_p = frozen.first;
    props.spec_heat_v = frozen.second;
    props.gamma_s = frozen.first / frozen.second;
    props.gas_moles = gas_moles;
    props.total_moles = gas_moles + state.condensed_moles.sum();
    props.density = state.pressure / (gas_moles * R * state.temperature);
    props.speed_of_sound = std::sqrt(gas_moles * R * state.temperature * props.gamma_s);
    props.frozen_spec_heat_p = frozen.first;
    props.frozen_gamma = props.gamma_s;
    // Composition is held fixed, so a phase transition has no effect on frozen properties.
    props.pinned_transition = false;
    return props;
}

double get_equilibrium_gamma(const Gas& gas) {
    return get_thermo_equilibrium_properties(gas).gamma_s;
}

} //namespace Goddard
