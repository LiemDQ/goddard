#pragma once
#include <memory>
#include <string>
#include <vector>
#include "goddard/nozzle.hpp"
#include "goddard/gas.hpp"
#include "goddard/moc.hpp"

namespace Goddard {

/** One station of a kinetic nozzle integration. */
struct KineticNozzleStation {
    /** Axial position [m]. */
    double x;
    /** Flow velocity [m/s]. */
    double velocity;
    /** Mach number [-], with the frozen speed of sound. */
    double mach;
    /** Area ratio A/A_t [-]. */
    double area_ratio;
    /** Mixture state at the station, with frozen derivatives. */
    ThermodynamicState thermo;
    /** Raw state vector, for `Gas::restore_state`. */
    std::vector<double> state;
    /** Damkohler number of each species [-]; 0 for trace species. */
    std::vector<double> damkohler;
    /** Smallest non-zero Damkohler number [-]. */
    double Da_min;
    /** Name of the species with the smallest non-zero Damkohler number; empty if there is none. */
    std::string min_Da_species;
};

/** Result of a `KineticNozzle::solve`. */
struct KineticNozzleResults {
    /// Throat the integration starts from, solved with the chemistry of `NozzleOptions::chemistry`.
    ThroatCondition throat;
    /// Stations downstream of the throat, one per integration step.
    std::vector<KineticNozzleStation> stations;
    /**
     * True if the integration reached the exit of the profile. False if it stopped after
     * `max_steps` steps, in which case `stations` ends upstream of the exit.
     */
    bool reached_exit = false;
};

/**
 * Supersonic nozzle expansion with finite-rate chemistry.
 *
 * The throat is solved with the chemistry of `NozzleOptions::chemistry` (EQUILIBRIUM or FROZEN).
 * From there a Lagrangian reactor is marched along the profile, its volume following the area
 * change, and the velocity follows from the stagnation enthalpy. The first point of the profile
 * is taken as the throat. Profile coordinates are in metres: the solver advances x by u dt.
 * Condensed species are not supported.
 */
class KineticNozzle {
public:

    // Construct from Solution reference (builds Gas internally with FROZEN chemistry)
    KineticNozzle(
        Cantera::Solution& gas,
        NozzleProfile& profile,
        double mass_flow_rate,
        NozzleOptions options = {});

    KineticNozzle(
        Cantera::Solution& gas,
        NozzleProfile& profile,
        double mass_flow_rate,
        std::vector<double> inlet_state,
        NozzleOptions options = {});

    /**
     * @param gas gas at the inlet (chamber) state, with a reaction mechanism. The solver works on
     * its own copy of `gas`; the caller's `Gas` is not modified.
     * @param profile Diverging wall contour starting at the throat [m].
     * @param mass_flow_rate Mass flow rate [kg/s].
     * @param options Throat chemistry (EQUILIBRIUM or FROZEN).
     * @throws std::invalid_argument for KINETIC throat chemistry.
     * @throws NotImplementedError if `gas` holds condensed species.
     */
    KineticNozzle(
        const Gas& gas,
        NozzleProfile& profile,
        double mass_flow_rate,
        NozzleOptions options = {});

    /**
     * @param gas product gas with a reaction mechanism; only its phase definition is used. The
     * solver works on its own copy.
     * @param profile Diverging wall contour starting at the throat [m].
     * @param mass_flow_rate Mass flow rate [kg/s].
     * @param inlet_state Inlet (chamber) state vector, as returned by `Gas::save_state`.
     * @param options Throat chemistry (EQUILIBRIUM or FROZEN).
     * @throws std::invalid_argument for KINETIC throat chemistry.
     * @throws NotImplementedError if the inlet state holds condensed species.
     */
    KineticNozzle(
        const Gas& gas,
        NozzleProfile& profile,
        double mass_flow_rate,
        std::vector<double> inlet_state,
        NozzleOptions options = {});

    /**
     * Solve the throat, then integrate the finite-rate expansion from the throat to the end of
     * the profile. Each step is limited by both `dt_max` and `dx_max`.
     *
     * @param dt_max Largest time step [s].
     * @param dx_max Largest axial step [m].
     * @param max_steps Largest number of steps [-]; see `KineticNozzleResults::reached_exit`.
     * @throws std::runtime_error if the reactor integration fails.
     */
    KineticNozzleResults solve(double dt_max = 1e-6, double dx_max = 1e-3, int max_steps = 100000);

    /// Wall contour [m]; its first point is the throat.
    NozzleProfile profile;
    /// Mass flow rate [kg/s].
    double mdot;
    /// Throat chemistry (EQUILIBRIUM or FROZEN).
    NozzleOptions opts;

protected:
    Nozzle m_throat_solver;
    std::vector<double> m_inlet_state;
    Gas m_gas;
};

} // namespace Goddard
