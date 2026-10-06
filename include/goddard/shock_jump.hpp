#pragma once
#include <vector>
#include "goddard/chemistry.hpp"
#include "goddard/gas.hpp"
#include "goddard/shocks.hpp"

/**
 * Newton kernel for the jump conditions across normal shocks and detonations, shared by
 * `ShockSolver` and `DetonationSolver`. Internal: not part of the Python API.
 *
 * Every function takes the upstream state in `gas` and leaves the downstream state in it.
 */

namespace Goddard {

/** Which pair of jump conditions `solve_shock_jump` solves. */
enum class JumpCondition {
    /** Gas enters a stationary wave at the given velocity (RP-1311 section 7.1). */
    INCIDENT,
    /** Gas moving at the given velocity is brought to rest by a reflected shock (section 7.2). */
    REFLECTED,
    /**
     * Detonation whose products leave at their equilibrium sound speed (RP-1311 chapter 8). The
     * velocity is not an input: it follows from the solution.
     */
    CHAPMAN_JOUGUET
};

/** Downstream-over-upstream ratios across a shock or detonation. */
struct ShockJump {
    double pressure_ratio;
    double temperature_ratio;
    double density_ratio;
};

/** Shock result with `valid = false` and every other field set to -1. */
ShockResult invalid_shock();

/**
 * Solve the jump conditions by Newton's method in x = (ln P_down/P_up, ln T_down/T_up),
 * NASA RP-1311 Part I, chapters 7 and 8.
 *
 * `gas` holds the upstream state on entry and the downstream state on exit. `gas.chemistry`
 * selects frozen or equilibrium downstream properties; the derivations are in
 * docs/theory/shocks.md and docs/theory/detonations.md.
 *
 * @param velocity u1 for INCIDENT, the particle velocity u_p of gas 2 for REFLECTED, ignored for
 * CHAPMAN_JOUGUET [m/s]
 * @throws std::invalid_argument for CHAPMAN_JOUGUET with chemistry other than EQUILIBRIUM
 * @throws ConvergenceError if the iteration does not converge
 */
ShockJump solve_shock_jump(Gas& gas, JumpCondition condition, double velocity,
                           double pressure_ratio_guess, double temperature_ratio_guess,
                           const SolverOptions& opts);

/**
 * Solve a shock jump with the given chemistry. EQUILIBRIUM first solves the frozen jump and
 * uses it as the initial guess (docs/theory/shocks.md, "Normal shock relations for real gases").
 *
 * `gas` holds the upstream state on entry and the downstream state on exit, with
 * `gas.chemistry` set to `chemistry`.
 */
ShockJump solve_shock_jump_with_chemistry(Gas& gas, JumpCondition condition, double velocity,
                                          GasChemistry chemistry,
                                          double pressure_ratio_guess, double temperature_ratio_guess,
                                          const SolverOptions& opts);

/** Frozen sound speed of the gas at its current state [m/s]. */
double frozen_sound_speed(Gas& gas);

/**
 * Jump across an incident normal shock, without stagnation quantities. `gas` holds state 1 on
 * entry and state 2 on exit. Requires u1 >= the frozen sound speed of state 1.
 */
ShockJump incident_jump(Gas& gas, double u1, GasChemistry chemistry, const SolverOptions& opts);

/**
 * Incident normal shock. `gas` holds state 1 on entry and state 2 on exit.
 *
 * State 1 is taken as given, so its sound speed and stagnation pressure are frozen ones. For an
 * oblique shock, `u_tangential` enters only the stagnation pressures.
 */
ShockResult incident_shock(Gas& gas, double u_normal, double u_tangential,
                           GasChemistry chemistry, const SolverOptions& opts);

/**
 * Shock reflected from a closed end. `gas` holds state 2 on entry, with `gas.chemistry` the
 * chemistry of the incident wave, and state 5 on exit.
 *
 * In the reflected-shock frame gas 2 enters at u_p + W_R and gas 5 leaves at W_R, where
 * W_R = u_p/(rho5/rho2 - 1) is the reflected shock speed.
 *
 * @param particle_velocity lab-frame velocity u_p of gas 2 [m/s]
 */
ShockResult reflected_shock_from_state2(Gas& gas, double particle_velocity, GasChemistry chemistry,
                                        double pressure_ratio_guess, double temperature_ratio_guess,
                                        const SolverOptions& opts);

} // namespace Goddard
