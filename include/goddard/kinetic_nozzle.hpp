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
    /** Axial position, in the length unit of the profile. */
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

struct KineticNozzleResults {
    ThroatCondition throat;
    std::vector<KineticNozzleStation> stations;
    /**
     * True if the integration reached the exit of the profile. False if it stopped after
     * `max_steps` steps, in which case `stations` ends upstream of the exit.
     */
    bool reached_exit = false;
};

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

    // Construct from Gas directly (chemistry is overridden to FROZEN internally). The solver works on its own copy of `gas`; the caller's `Gas` is not modified.
    KineticNozzle(
        const Gas& gas,
        NozzleProfile& profile,
        double mass_flow_rate,
        NozzleOptions options = {});

    KineticNozzle(
        const Gas& gas,
        NozzleProfile& profile,
        double mass_flow_rate,
        std::vector<double> inlet_state,
        NozzleOptions options = {});

    KineticNozzleResults solve(double dt_max = 1e-6, double dx_max = 1e-3, int max_steps = 100000);

    NozzleProfile profile;
    double mdot;
    NozzleOptions opts;

protected:
    Nozzle m_throat_solver;
    std::vector<double> m_inlet_state;
    Gas m_gas;
};

} // namespace Goddard
