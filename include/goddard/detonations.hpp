#pragma once
#include <vector>
#include "goddard/chemistry.hpp"
#include "goddard/gas.hpp"
#include "goddard/shocks.hpp"

namespace Goddard {

/**
 * Root of the Rayleigh line and the equilibrium Hugoniot for a wave faster than the
 * Chapman-Jouguet (CJ) speed.
 *
 * Above the CJ speed the Rayleigh line cuts the Hugoniot twice; at the CJ speed the two roots
 * merge, and below it there is no steady solution. Both branches therefore take a drive factor
 * u1/u_CJ > 1.
 */
enum class DetonationBranch {
    OVERDRIVEN,     ///< Strong root: higher pressure, products subsonic relative to the wave.
    UNDERDRIVEN     ///< Weak root: lower pressure, products supersonic relative to the wave.
};

/**
 * Jump conditions across a planar detonation, in the frame of the wave.
 *
 * State 1 is the unburned gas ahead of the wave and state 2 the products at chemical equilibrium
 * behind it. Ratios are state 2 over state 1.
 *
 * An invalid result (a wave slower than the CJ speed, or a mixture that releases no heat) has
 * `valid = false` and every other field set to -1.
 */
struct DetonationResult {
    /** False if the inputs admit no detonation. */
    bool valid;
    /**
     * Wave speed relative to the unburned gas [m/s]. NaN from the perfect-gas functions, which
     * are dimensionless.
     */
    double velocity;
    /** Wave speed over the CJ speed of the mixture [-]. 1 for a CJ detonation. */
    double drive_factor;
    /** Wave Mach number u1/a1, with the frozen sound speed of the unburned gas [-]. */
    double mach_in;
    /**
     * Mach number of the products relative to the wave, u2/a2, with the equilibrium sound speed
     * [-]. 1 for a CJ detonation, below 1 when overdriven and above 1 when under-driven.
     */
    double mach_out;
    /** Static pressure ratio P2/P1 [-]. */
    double static_pressure_ratio;
    /** Static temperature ratio T2/T1 [-]. */
    double static_temperature_ratio;
    /** Density ratio rho2/rho1 [-]. It is also the ratio u1/u2 of the velocities relative to the wave. */
    double density_ratio;
    /** Molecular weight ratio M2/M1 [-], with M = 1/n the gas molar mass. */
    double molecular_weight_ratio;
    /**
     * Stagnation pressure ratio P02/P01 [-], in the frame of the wave. P01 is the frozen stagnation
     * pressure of the unburned gas and P02 the equilibrium stagnation pressure of the products.
     */
    double total_pressure_ratio;
    /**
     * Frozen (non-reacting) shock at the same wave speed: the von Neumann spike at the head of
     * the ZND reaction zone. For an under-driven detonation it is not on a steady ZND path.
     */
    ShockResult von_neumann;
};

/**
 * A detonation reflected as a shock from the closed end of a tube.
 *
 * The incident detonation leaves its products (state 2) moving toward the end wall; the
 * reflected shock brings them back to rest (state 5), at equilibrium.
 */
struct ReflectedDetonationResult {
    /** False if the incident detonation or the reflected shock is invalid. */
    bool valid;
    /** Incident detonation, state 1 to state 2, in the frame of the detonation. */
    DetonationResult incident;
    /**
     * Reflected shock, state 2 to state 5, in the reflected-shock frame. `mach_in` is the
     * reflected-shock Mach number relative to the products, with their equilibrium sound speed.
     */
    ShockResult reflected;
};

// --- Perfect gas free functions (pure, no state) ---
//
// One-gamma model: reactants and products share gamma and molar mass, and the reaction releases
// heat q per unit mass. The heat release enters as Q = q/(R T1), with R the specific gas constant.

/**
 * Chapman-Jouguet detonation in a calorically perfect gas that releases heat.
 *
 * The CJ Mach number is M_CJ = sqrt(H) + sqrt(H + 1), with H = (gamma^2 - 1) Q / (2 gamma).
 *
 * @param gamma ratio of specific heats of reactants and products [-]
 * @param heat_release heat release q/(R T1), with R the specific gas constant [-]
 * @return jump conditions; invalid if gamma <= 1 or heat_release < 0
 */
DetonationResult chapman_jouguet_detonation(double gamma, double heat_release);

/**
 * Overdriven or under-driven detonation in a calorically perfect gas that releases heat.
 *
 * With no heat release the overdriven branch is a normal shock and the under-driven branch a
 * vanishing wave.
 *
 * @param drive_factor wave speed over the CJ speed [-]
 * @param gamma ratio of specific heats of reactants and products [-]
 * @param heat_release heat release q/(R T1), with R the specific gas constant [-]
 * @param branch root of the Rayleigh line and the Hugoniot
 * @return jump conditions; invalid if drive_factor < 1, gamma <= 1 or heat_release < 0
 */
DetonationResult detonation(double drive_factor, double gamma, double heat_release,
                            DetonationBranch branch = DetonationBranch::OVERDRIVEN);

/**
 * Detonation and its reflection from the closed end of a tube, in a calorically perfect gas that
 * releases heat.
 *
 * @param drive_factor wave speed over the CJ speed [-]
 * @param gamma ratio of specific heats of reactants and products [-]
 * @param heat_release heat release q/(R T1), with R the specific gas constant [-]
 * @param branch root of the Rayleigh line and the Hugoniot
 * @return incident detonation and reflected shock; invalid if the detonation is invalid
 */
ReflectedDetonationResult reflected_detonation(double drive_factor, double gamma, double heat_release,
                                               DetonationBranch branch = DetonationBranch::OVERDRIVEN);

// --- DetonationSolver: real-gas detonations with state management ---

/**
 * Planar detonations in a real gas: Chapman-Jouguet (CJ) detonations following Gordon & McBride,
 * NASA RP-1311 Part I, chapter 8, and overdriven and under-driven detonations at a given wave
 * speed.
 *
 * The unburned gas (state 1) is the state of the `Gas` at construction. It is taken as given:
 * it need not be at equilibrium, and its sound speed is the frozen one. The products are always
 * at chemical equilibrium, whatever the chemistry of the `Gas`. The CJ detonation is solved once,
 * at construction.
 *
 * Every result also reports the von Neumann state, the frozen shock at the same wave speed.
 *
 * @note The state accessors return the same `Gas`, which shares its `Solution` with the `Gas`
 * passed to the constructor. A reference from one accessor is changed by the next call to any
 * method of the solver.
 */
class DetonationSolver {
public:
    /**
     * @param gas unburned gas: FROZEN or EQUILIBRIUM chemistry, with no condensed species
     * @param options Newton tolerance on the log pressure and temperature ratios, and iteration limit
     * @throws std::invalid_argument for PERFECT_GAS or KINETIC chemistry, or if condensed species
     * are present
     * @throws ConvergenceError if the CJ detonation cannot be solved
     */
    DetonationSolver(Gas gas, SolverOptions options = {});

    /**
     * Chapman-Jouguet detonation: the slowest steady detonation, whose products leave at their
     * equilibrium sound speed.
     *
     * The result is the one solved at construction. `post_detonation_state()` and
     * `von_neumann_state()` are then the CJ states.
     *
     * @return jump conditions; invalid if the mixture releases no heat
     */
    DetonationResult chapman_jouguet();

    /**
     * Overdriven or under-driven detonation at a multiple of the CJ speed.
     *
     * The branch is chosen by the initial guess, and checked against the Mach number of the
     * products after the solve.
     *
     * @param drive_factor wave speed over the CJ speed [-]
     * @param branch root of the Rayleigh line and the Hugoniot
     * @return jump conditions; invalid if drive_factor < 1. A drive factor of 1 returns the CJ
     * detonation.
     * @throws std::invalid_argument if condensed species form behind the wave
     * @throws ConvergenceError if the jump conditions cannot be solved on the requested branch,
     * e.g. for a drive factor so close to 1 that the two branches merge
     */
    DetonationResult detonation(double drive_factor, DetonationBranch branch = DetonationBranch::OVERDRIVEN);

    /**
     * Overdriven or under-driven detonation at a given wave speed.
     *
     * @param velocity wave speed relative to the unburned gas [m/s]
     * @param branch root of the Rayleigh line and the Hugoniot
     * @return jump conditions; invalid if the velocity is below the CJ speed
     */
    DetonationResult detonation_from_velocity(double velocity,
                                              DetonationBranch branch = DetonationBranch::OVERDRIVEN);

    /**
     * Detonation reflected as a shock from the closed end of a tube, with the products at
     * equilibrium on both sides of the reflected shock.
     *
     * `post_detonation_state()` is then state 5, the products brought to rest.
     *
     * @param drive_factor wave speed of the incident detonation over the CJ speed [-]
     * @param branch root of the Rayleigh line and the Hugoniot for the incident detonation
     * @return incident detonation and reflected shock; invalid if the incident detonation is invalid
     */
    ReflectedDetonationResult reflected_detonation(double drive_factor = 1.0,
                                                   DetonationBranch branch = DetonationBranch::OVERDRIVEN);

    /** Gas at the unburned state, with the chemistry it had at construction. */
    const Gas& pre_detonation_state() const;
    /**
     * Gas at the products of the most recent solve, with EQUILIBRIUM chemistry. After an invalid
     * result it is the unburned state.
     */
    const Gas& post_detonation_state() const;
    /**
     * Gas behind the leading frozen shock of the most recent solve, with FROZEN chemistry. After
     * an invalid result it is the unburned state.
     */
    const Gas& von_neumann_state() const;

private:
    /** Restore the unburned state and the construction chemistry. */
    void reset() const;
    /** Solve the CJ detonation from the unburned state and store it with its states. */
    void solve_chapman_jouguet();
    /**
     * Overdriven or under-driven detonation at `velocity`, whose drive factor is `drive_factor`.
     * Leaves the products in `m_gas` and records the post and von Neumann states.
     */
    DetonationResult solve_driven(double velocity, double drive_factor, DetonationBranch branch);
    /** Record the unburned state as both the post and the von Neumann state, for an invalid result. */
    void store_unburned_states();

    mutable Gas m_gas;
    GasChemistry m_chemistry;
    SolverOptions m_options;
    std::vector<double> m_pre_detonation_state;
    std::vector<double> m_post_detonation_state;
    std::vector<double> m_von_neumann_state;
    DetonationResult m_chapman_jouguet;
    std::vector<double> m_chapman_jouguet_state;
    std::vector<double> m_chapman_jouguet_von_neumann_state;
};

}
