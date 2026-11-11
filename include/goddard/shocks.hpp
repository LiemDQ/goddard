#pragma once
#include <vector>
#include <utility>
#include "goddard/chemistry.hpp"
#include "goddard/gas.hpp"
#include "goddard/status.hpp"

namespace Goddard {

/**
 * Jump conditions across one normal shock, in the frame of that shock.
 *
 * Ratios are downstream over upstream. For an oblique shock this is the normal component of the
 * flow, except `total_pressure_ratio`, which uses the full velocity.
 *
 * When the solve fails (`status.ok()` is false) every other field is NaN. Subsonic inflow has no
 * shock and gives reason `NO_SOLUTION`.
 */
struct ShockResult {
    /** Outcome of the solve. */
    SolveStatus status;
    /**
     * Upstream Mach number relative to the shock [-]. For the incident shock of `ShockSolver`
     * it uses the frozen sound speed of the pre-shock state.
     */
    double mach_in = NaN;
    /** Downstream Mach number relative to the shock [-], with the sound speed of the solver's chemistry. */
    double mach_out = NaN;
    /** Static pressure ratio P_out/P_in [-]. */
    double static_pressure_ratio = NaN;
    /** Static temperature ratio T_out/T_in [-]. */
    double static_temperature_ratio = NaN;
    /** Density ratio rho_out/rho_in [-]. It is also the ratio of upstream to downstream normal velocity. */
    double density_ratio = NaN;
    /** Stagnation pressure ratio P0_out/P0_in [-]. */
    double total_pressure_ratio = NaN;
};

/**
 * Incident and reflected normal shocks at the closed end of a shock tube.
 *
 * State 1 is the gas at rest ahead of the incident shock, state 2 the gas behind it, and state 5
 * the gas brought back to rest behind the shock reflected from the end wall.
 */
struct ReflectedShockResult {
    /**
     * Outcome of the solve: the failure of the incident shock, else that of the reflected shock.
     * When the reflected shock fails, `incident` is still valid.
     */
    SolveStatus status;
    /** Incident shock, state 1 to state 2, in the incident-shock frame. */
    ShockResult incident;
    /**
     * Reflected shock, state 2 to state 5, in the reflected-shock frame. `mach_in` is the
     * reflected-shock Mach number relative to gas 2.
     */
    ShockResult reflected;
};

/**
 * Oblique shock: wave and deflection angles, and the jump across the normal component.
 *
 * When the solve fails every other field is NaN. A deflection above the maximum deflection (a
 * detached shock) or a wave angle below the Mach angle gives reason `NO_SOLUTION`.
 */
struct ObliqueShockResult {
    /** Outcome of the solve. */
    SolveStatus status;
    /** Upstream Mach number [-]. */
    double mach_in = NaN;
    /** Downstream Mach number, including the tangential velocity [-]. */
    double mach_out = NaN;
    /** Jump across the normal component of the flow. */
    ShockResult shock;
    /** Wave angle between the shock and the upstream flow [rad]. */
    double beta = NaN;
    /** Flow deflection angle [rad]. */
    double theta = NaN;
};

// --- Perfect gas free functions (pure, no state) ---

/**
 * Normal shock in a calorically perfect gas.
 *
 * @param mach upstream Mach number [-]
 * @param gamma ratio of specific heats [-]
 * @return jump conditions; `NO_SOLUTION` if mach < 1
 * @throws std::invalid_argument if gamma <= 1
 */
ShockResult normal_shock(double mach, double gamma);

/**
 * Incident shock and its reflection from the closed end of a shock tube, in a calorically
 * perfect gas.
 *
 * @param mach incident shock Mach number relative to the gas at rest ahead of it [-]
 * @param gamma ratio of specific heats [-]
 * @return both jumps; `NO_SOLUTION` if mach < 1
 * @throws std::invalid_argument if gamma <= 1
 */
ReflectedShockResult reflected_shock(double mach, double gamma);

/**
 * Weak and strong wave angles of an oblique shock with a given deflection, in a calorically
 * perfect gas.
 *
 * @param mach upstream Mach number [-]
 * @param deflection_angle flow deflection [rad]
 * @param gamma ratio of specific heats [-]
 * @return (weak, strong) wave angles [rad]. Both are NaN if the deflection is negative or exceeds
 * the maximum deflection (a detached shock), or if mach < 1 or gamma <= 1. At zero deflection
 * they are the Mach angle and pi/2.
 */
std::pair<double,double> oblique_shock_wave_angle(double mach, double deflection_angle, double gamma);

/**
 * Flow deflection of an oblique shock with a given wave angle, in a calorically perfect gas.
 *
 * @param mach upstream Mach number [-]
 * @param wave_angle wave angle [rad]
 * @param gamma ratio of specific heats [-]
 * @return deflection angle [rad]; negative if the wave angle is below the Mach angle
 */
double oblique_shock_deflection_angle(double mach, double wave_angle, double gamma);

/**
 * Wave angle at which the deflection of an oblique shock is largest, in a calorically perfect
 * gas (NACA 1135, eq. 168).
 *
 * @param mach upstream Mach number [-]
 * @param gamma ratio of specific heats [-]
 * @return wave angle [rad]; NaN if mach < 1 or gamma <= 1
 */
double oblique_shock_max_deflection_wave_angle(double mach, double gamma);

/**
 * Largest deflection an attached oblique shock can produce, in a calorically perfect gas.
 *
 * @param mach upstream Mach number [-]
 * @param gamma ratio of specific heats [-]
 * @return maximum deflection angle [rad]; NaN if mach < 1 or gamma <= 1
 */
double oblique_shock_max_deflection(double mach, double gamma);

/**
 * Oblique shock with a given wave angle, in a calorically perfect gas.
 *
 * @param mach upstream Mach number [-]
 * @param wave_angle wave angle [rad]
 * @param gamma ratio of specific heats [-]
 * @return oblique shock; `NO_SOLUTION` if the wave angle is below the Mach angle
 * @throws std::invalid_argument if gamma <= 1
 */
ObliqueShockResult oblique_shock_from_wave_angle(
    double mach, double wave_angle, double gamma);

/**
 * Oblique shock with a given deflection, in a calorically perfect gas.
 *
 * @param mach upstream Mach number [-]
 * @param deflection_angle flow deflection [rad]
 * @param gamma ratio of specific heats [-]
 * @param weak true for the weak solution, false for the strong one
 * @return oblique shock; `NO_SOLUTION` if mach < 1, or if the deflection is negative or exceeds
 * the maximum deflection (a detached shock)
 * @throws std::invalid_argument if gamma <= 1
 */
ObliqueShockResult oblique_shock_from_deflection(
    double mach, double deflection_angle, double gamma, bool weak = true);

// --- ShockSolver: chemistry-dispatched solver with state management ---

/**
 * Normal, reflected and oblique shocks in a real gas, following Gordon & McBride, NASA RP-1311
 * Part I, chapter 7.
 *
 * The chemistry of the `Gas` selects the model behind the shock: PERFECT_GAS uses the perfect-gas
 * relations with `Gas::gamma_s()`, FROZEN holds the composition fixed, and EQUILIBRIUM brings the
 * post-shock gas to chemical equilibrium. KINETIC chemistry and condensed species are not
 * supported.
 *
 * The pre-shock state is the state of the `Gas` at construction. It is taken as given: it need not
 * be at equilibrium (e.g. an unburned fuel-oxidizer mixture), and its sound speed is always the
 * frozen one. Mach-number inputs are converted to velocity with that frozen sound speed.
 *
 * Every solve returns a result with a `status`. A numerical failure (the jump conditions do not
 * converge, the thermodynamic backend fails, condensed species form behind the shock) is reported
 * there and the result is NaN; the solver throws only for invalid arguments.
 *
 * @note An equilibrium shock in an exothermic mixture only exists above the Chapman-Jouguet
 * detonation speed; below it the solve fails with `NOT_CONVERGED`. Above it the shock is an
 * overdriven detonation, which `DetonationSolver` also computes. Oblique shocks in such a
 * mixture are oblique detonations and are not modelled.
 */
class ShockSolver {
public:
    /**
     * @param gas pre-shock gas and chemistry model. The solver works on its own copy of `gas`; the caller's `Gas` is not modified.
     * @param options Newton tolerance on the log pressure and temperature ratios, and iteration limit
     * @throws std::invalid_argument for KINETIC chemistry
     */
    ShockSolver(const Gas& gas, SolverOptions options = {});

    /**
     * Normal shock.
     *
     * @param mach upstream Mach number, with the frozen pre-shock sound speed [-]
     * @return jump conditions. `NO_SOLUTION` if mach < 1 or if condensed species form behind
     * the shock; `NOT_CONVERGED` if the jump conditions cannot be solved
     * @throws std::invalid_argument if condensed species are present before the shock
     */
    ShockResult normal_shock(double mach);

    /**
     * Normal shock from the upstream velocity.
     *
     * @param velocity upstream velocity relative to the shock [m/s]
     * @return jump conditions; `NO_SOLUTION` if the velocity is below the frozen pre-shock sound speed
     */
    ShockResult normal_shock_from_velocity(double velocity);

    /**
     * Incident and reflected shocks in a shock tube, both with the solver's chemistry.
     *
     * @param mach incident shock Mach number, with the frozen pre-shock sound speed [-]
     * @return both jumps; `post_shock_state()` is then state 5. If the reflected shock fails,
     * `incident` is still valid and `post_shock_state()` throws.
     */
    ReflectedShockResult reflected_shock(double mach);

    /**
     * Incident and reflected shocks in a shock tube, each with its own chemistry.
     *
     * @param mach incident shock Mach number, with the frozen pre-shock sound speed [-]
     * @param incident_chemistry FROZEN or EQUILIBRIUM, for state 2
     * @param reflected_chemistry FROZEN or EQUILIBRIUM, for state 5
     * @throws std::invalid_argument for any other chemistry, or if the solver's gas is PERFECT_GAS
     */
    ReflectedShockResult reflected_shock(
        double mach, GasChemistry incident_chemistry, GasChemistry reflected_chemistry);

    /**
     * Incident and reflected shocks in a shock tube, from the incident shock speed.
     *
     * @param velocity incident shock speed relative to the gas at rest ahead of it [m/s]
     */
    ReflectedShockResult reflected_shock_from_velocity(double velocity);

    /**
     * Incident and reflected shocks in a shock tube, from the incident shock speed, each with its
     * own chemistry.
     *
     * @param velocity incident shock speed relative to the gas at rest ahead of it [m/s]
     * @param incident_chemistry FROZEN or EQUILIBRIUM, for state 2
     * @param reflected_chemistry FROZEN or EQUILIBRIUM, for state 5
     */
    ReflectedShockResult reflected_shock_from_velocity(
        double velocity, GasChemistry incident_chemistry, GasChemistry reflected_chemistry);

    /**
     * Oblique shock with a given wave angle.
     *
     * @param mach upstream Mach number, with the frozen pre-shock sound speed [-]
     * @param wave_angle wave angle [rad]
     * @return oblique shock; `NO_SOLUTION` if the normal Mach number is below 1
     */
    ObliqueShockResult oblique_shock_from_wave_angle(double mach, double wave_angle);

    /**
     * Oblique shock with a given deflection.
     *
     * For real gases the maximum deflection is found by golden-section search, and the wave angle
     * by bisection on the weak or strong branch.
     *
     * @param mach upstream Mach number, with the frozen pre-shock sound speed [-]
     * @param deflection flow deflection [rad]
     * @param weak true for the weak solution, false for the strong one
     * @return oblique shock; `NO_SOLUTION` if mach < 1, or if the deflection is negative or
     * exceeds the maximum deflection (a detached shock)
     */
    ObliqueShockResult oblique_shock_from_deflection(double mach, double deflection, bool weak = true);

    /**
     * Largest deflection an attached oblique shock can produce.
     *
     * `post_shock_state()` is then the state behind the shock at maximum deflection.
     *
     * This solve returns a number, so it has no status: a failure throws.
     *
     * @param mach upstream Mach number, with the frozen pre-shock sound speed [-]
     * @return maximum deflection angle [rad]
     * @throws SolveError with reason `NO_SOLUTION` if mach < 1, and for any failure of the search
     */
    double max_deflection(double mach);

    /** Gas at the pre-shock state: an independent copy. */
    Gas pre_shock_state() const;
    /**
     * Gas at the state behind the most recently solved shock.
     *
     * For PERFECT_GAS chemistry it is the pre-shock gas with its temperature and pressure
     * multiplied by the perfect-gas jump ratios, and its composition unchanged.
     *
     * The result is an independent copy: later solves do not change it.
     *
     * @throws SolveError if the most recent solve failed (with that solve's status), or if no
     * shock has been solved yet
     */
    Gas post_shock_state() const;

private:
    /** Restore the pre-shock state and the construction chemistry. */
    void reset() const;
    /** Frozen sound speed of the pre-shock state [m/s]. */
    double pre_shock_sound_speed() const;
    /** Incident and reflected shocks for PERFECT_GAS chemistry, from the incident Mach number. */
    ReflectedShockResult perfect_gas_reflected_shock(double mach);
    /**
     * Record the outcome of a solve and restore the construction chemistry: the current state
     * becomes the post-shock state if `status` is ok, and there is no post-shock state otherwise.
     */
    void store_post_shock_state(const SolveStatus& status);
    /**
     * Apply a perfect-gas jump to the gas, which must hold the pre-shock state: T and P are
     * multiplied by the given ratios [-] and the composition is kept.
     */
    void apply_perfect_gas_jump(double temperature_ratio, double pressure_ratio);

    mutable Gas m_gas;
    GasChemistry m_chemistry;
    SolverOptions m_options;
    std::vector<double> m_pre_shock_state;
    /// State behind the most recent shock; meaningful only while `m_post_shock_status` is ok.
    std::vector<double> m_post_shock_state;
    /// Status of the most recent solve; NOT_COMPUTED until a shock has been solved.
    SolveStatus m_post_shock_status;
};

}
