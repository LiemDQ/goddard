#pragma once
#include "cantera/core.h"
namespace Goddard {

using Composition = Cantera::Composition;

enum class GasChemistry {
    PERFECT_GAS,    // constant gamma, no Cantera dependency
    FROZEN,         // composition fixed, variable thermodynamic properties via Cantera
    EQUILIBRIUM,    // full chemical equilibrium at each point via Cantera
    KINETIC         // finite-rate chemistry via Cantera kinetics
};

struct SolverOptions {
    double abstol = 1e-6;
    double reltol = 1e-5;
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
