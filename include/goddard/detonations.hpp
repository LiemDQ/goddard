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
 * When the solve fails (`status.ok()` is false) every other field is NaN. A wave slower than the
 * CJ speed and a mixture that releases no heat have no detonation and give reason `NO_SOLUTION`.
 */
struct DetonationResult {
    /** Outcome of the solve. */
    SolveStatus status;
    /**
     * Wave speed relative to the unburned gas [m/s]. NaN from the perfect-gas functions, which
     * are dimensionless.
     */
    double velocity = NaN;
    /** Wave speed over the CJ speed of the mixture [-]. 1 for a CJ detonation. */
    double drive_factor = NaN;
    /** Wave Mach number u1/a1, with the frozen sound speed of the unburned gas [-]. */
    double mach_in = NaN;
    /**
     * Mach number of the products relative to the wave, u2/a2, with the equilibrium sound speed
     * [-]. 1 for a CJ detonation, below 1 when overdriven and above 1 when under-driven.
     */
    double mach_out = NaN;
    /** Static pressure ratio P2/P1 [-]. */
    double static_pressure_ratio = NaN;
    /** Static temperature ratio T2/T1 [-]. */
    double static_temperature_ratio = NaN;
    /** Density ratio rho2/rho1 [-]. It is also the ratio u1/u2 of the velocities relative to the wave. */
    double density_ratio = NaN;
    /** Molecular weight ratio M2/M1 [-], with M = 1/n the gas molar mass. */
    double molecular_weight_ratio = NaN;
    /**
     * Stagnation pressure ratio P02/P01 [-], in the frame of the wave. P01 is the frozen stagnation
     * pressure of the unburned gas and P02 the equilibrium stagnation pressure of the products.
     */
    double total_pressure_ratio = NaN;
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
    /**
     * Outcome of the solve: the failure of the incident detonation, else that of the reflected
     * shock. When the reflected shock fails, `incident` is still valid.
     */
    SolveStatus status;
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
 * @return jump conditions; `NO_SOLUTION` if heat_release < 0
 * @throws std::invalid_argument if gamma <= 1
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
 * @return jump conditions; `NO_SOLUTION` if drive_factor < 1 or heat_release < 0
 * @throws std::invalid_argument if gamma <= 1
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
 * @return incident detonation and reflected shock; failed if the detonation is
 * @throws std::invalid_argument if gamma <= 1
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
 * at chemical equilibrium, whatever the chemistry of the `Gas`.
 *
 * The constructor solves nothing. The CJ detonation is solved by the first call that needs it
 * (any of the solve methods) and kept for the later ones.
 *
 * Every solve returns a result with a `status`. A numerical failure (the jump conditions do not
 * converge or converge off the requested branch, condensed species form behind the wave) is
 * reported there and the result is NaN; the solver throws only for invalid arguments.
 *
 * Every result also reports the von Neumann state, the frozen shock at the same wave speed.
 *
 * The solver works on its own copy of the `Gas` passed to the constructor, and the state
 * accessors return independent copies, so neither the caller's `Gas` nor an earlier accessor
 * result changes when the solver runs.
 */
class DetonationSolver {
public:
    /**
     * @param gas unburned gas: FROZEN or EQUILIBRIUM chemistry, with no condensed species
     * @param options Newton tolerance on the log pressure and temperature ratios, and iteration limit
     * @throws std::invalid_argument for PERFECT_GAS or KINETIC chemistry, or if condensed species
     * are present
     */
    DetonationSolver(const Gas& gas, SolverOptions options = {});

    /**
     * Chapman-Jouguet detonation: the slowest steady detonation, whose products leave at their
     * equilibrium sound speed.
     *
     * Solved on the first call and kept. `post_detonation_state()` and `von_neumann_state()`
     * are then the CJ states.
     *
     * @return jump conditions. `NO_SOLUTION` if the mixture releases no heat or if condensed
     * species form behind the wave; `NOT_CONVERGED` if the jump conditions cannot be solved
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
     * @return jump conditions. A drive factor of 1 returns the CJ detonation. `NO_SOLUTION` if
     * drive_factor < 1 or if condensed species form behind the wave; `NOT_CONVERGED` if the jump
     * conditions cannot be solved on the requested branch, e.g. for a drive factor so close to 1
     * that the two branches merge. If the CJ detonation cannot be solved, its failure.
     */
    DetonationResult detonation(double drive_factor, DetonationBranch branch = DetonationBranch::OVERDRIVEN);

    /**
     * Overdriven or under-driven detonation at a given wave speed.
     *
     * @param velocity wave speed relative to the unburned gas [m/s]
     * @param branch root of the Rayleigh line and the Hugoniot
     * @return jump conditions; `NO_SOLUTION` if the velocity is below the CJ speed
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
     * @return incident detonation and reflected shock. If the reflected shock fails, `incident`
     * is still valid and `post_detonation_state()` throws.
     */
    ReflectedDetonationResult reflected_detonation(double drive_factor = 1.0,
                                                   DetonationBranch branch = DetonationBranch::OVERDRIVEN);

    /** Gas at the unburned state, with the chemistry it had at construction. An independent copy. */
    Gas pre_detonation_state() const;
    /**
     * Gas at the products of the most recent solve, with EQUILIBRIUM chemistry. An independent
     * copy.
     *
     * @throws SolveError if the most recent solve failed (with that solve's status), or if no
     * detonation has been solved yet
     */
    Gas post_detonation_state() const;
    /**
     * Gas behind the leading frozen shock of the most recent solve, with FROZEN chemistry. An
     * independent copy.
     *
     * @throws SolveError if the most recent incident detonation failed (with that solve's
     * status), or if no detonation has been solved yet
     */
    Gas von_neumann_state() const;

private:
    /** Restore the unburned state and the construction chemistry. */
    void reset() const;
    /**
     * The CJ detonation, solved from the unburned state on the first call and stored with its
     * states. A failed solve is stored too and is not repeated.
     */
    const DetonationResult& cached_chapman_jouguet();
    /**
     * CJ jump and states from the unburned state.
     *
     * @throws SolveError if there is no CJ detonation or it cannot be solved.
     */
    DetonationResult solve_chapman_jouguet();
    /**
     * Overdriven or under-driven detonation at `velocity`, whose drive factor is `drive_factor`.
     * Requires a solved CJ detonation. Leaves the products in `m_gas` and records the post and
     * von Neumann states.
     *
     * @throws SolveError if the jump conditions cannot be solved on the requested branch.
     */
    DetonationResult solve_driven(double velocity, double drive_factor, DetonationBranch branch);
    /** Run `solve_driven` as a solve: failures become the status of the result. */
    DetonationResult driven_detonation(double velocity, double drive_factor, DetonationBranch branch);
    /** Record a failed solve: there is then no post-detonation or von Neumann state. */
    void store_failure(const SolveStatus& status);

    mutable Gas m_gas;
    GasChemistry m_chemistry;
    SolverOptions m_options;
    std::vector<double> m_pre_detonation_state;
    /// Products of the most recent solve; meaningful only while `m_post_detonation_status` is ok.
    std::vector<double> m_post_detonation_state;
    /// von Neumann state of the most recent solve; meaningful only while `m_von_neumann_status` is ok.
    std::vector<double> m_von_neumann_state;
    /// Status of the most recent solve; NOT_COMPUTED until a detonation has been solved.
    SolveStatus m_post_detonation_status;
    /// Status of the most recent incident detonation; NOT_COMPUTED until one has been solved.
    SolveStatus m_von_neumann_status;
    bool m_chapman_jouguet_solved = false;
    DetonationResult m_chapman_jouguet;
    std::vector<double> m_chapman_jouguet_state;
    std::vector<double> m_chapman_jouguet_von_neumann_state;
};

}
