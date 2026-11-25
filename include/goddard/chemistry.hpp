#pragma once
#include "cantera/core.h"
namespace Goddard {

using Composition = Cantera::Composition;

/**
 * Chemistry model of a flow solver, set in the solver's options or arguments. A `Gas` has no
 * chemistry model of its own; its derived properties take a `CompositionModel`.
 */
enum class GasChemistry {
    /// Calorically perfect gas: the flow solver holds the ratio of specific heats constant.
    PERFECT_GAS,
    /// Composition held fixed; thermodynamic properties vary with temperature.
    FROZEN,
    /// Composition at chemical equilibrium at every state (shifting equilibrium).
    EQUILIBRIUM,
    /// Finite-rate chemistry from the reaction mechanism.
    KINETIC
};

/**
 * How the composition responds when the state of a `Gas` is perturbed. Selects the derivatives
 * behind `Gas::gamma_s`, `Gas::speed_of_sound`, `Gas::expansion_properties` and the other derived
 * properties.
 */
enum class CompositionModel {
    /// Composition held fixed: cp/cv, and the frozen sound speed.
    FROZEN,
    /// Composition shifts to remain at chemical equilibrium (RP-1311 eqs. 2.56-2.66).
    EQUILIBRIUM
};

/**
 * Composition model of a flow solver's chemistry: EQUILIBRIUM for EQUILIBRIUM chemistry, FROZEN
 * for every other model.
 */
inline CompositionModel composition_model(GasChemistry chemistry) {
    return chemistry == GasChemistry::EQUILIBRIUM ? CompositionModel::EQUILIBRIUM
                                                  : CompositionModel::FROZEN;
}

/**
 * @brief Property held constant together with pressure by the multiphase equilibrium solver.
 *
 * Cantera cannot solve the HP and SP problems once the set of condensed phases depends on
 * temperature, so `Gas` solves them as a one-dimensional root find on temperature over
 * constant-temperature, constant-pressure solves. This enum selects which of the two.
 */
enum class EquilibriumProperty {
    ENTHALPY,   //!< Specific enthalpy [J/kg of mixture] is held fixed (the HP problem).
    ENTROPY     //!< Specific entropy [J/(kg.K) of mixture] is held fixed (the SP problem).
};

/** Convergence settings of an iterative solver (shock and detonation jumps, MoC unit processes). */
struct SolverOptions {
    /**
     * Absolute convergence tolerance. The quantity it bounds depends on the solver: the Newton
     * step in log(P) and log(T) [-] for shock and detonation jumps, the wave angle bracket or
     * deflection residual [rad] for oblique shock searches, and the corrector residuals and
     * point validity checks of the MoC unit processes.
     */
    double abstol = 1e-6;
    /// Maximum number of iterations [-].
    int max_iterations = 100;
};

/**
 * @brief Settings for the multiphase equilibrium solver used by `Gas`.
 *
 * `rtol` and `max_steps` are passed to Cantera's Gibbs (`MultiPhaseEquil`) solver for each
 * constant-temperature-and-pressure solve. The remaining entries control the outer
 * one-dimensional temperature root find that implements the HP and SP problems, which Cantera
 * cannot solve directly once the set of condensed phases depends on temperature.
 */
struct EquilibriumOptions {
    /** Relative tolerance of the inner constant-T, constant-P Gibbs solve [-]. */
    double rtol = 1e-9;
    /**
     * Maximum number of steps of the inner Gibbs solve. The Cantera default (1000) is not
     * enough near a dew point, where a barely-present condensed phase converges slowly.
     */
    int max_steps = 20000;
    /** Maximum number of doubling steps used to bracket the HP/SP temperature root. */
    int max_bracket_steps = 40;
    /** Relative convergence tolerance on temperature for the HP/SP outer root find [-]. */
    double T_rel_tol = 1e-9;
    /** Temperature [K] used to start the HP/SP root find when no better guess is available. */
    double T_default = 3000.0;
    /** Lowest temperature [K] the HP/SP root find will consider. */
    double T_min = 200.0;
    /** Highest temperature [K] the HP/SP root find will consider. */
    double T_max = 6000.0;
};

} // namespace Goddard
