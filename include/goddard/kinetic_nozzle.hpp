#pragma once
#include <memory>
#include <vector>
#include "goddard/nozzle.hpp"
#include "goddard/gas.hpp"
#include "goddard/moc.hpp"

namespace Goddard {

struct KineticNozzleStation {
    double x;
    double velocity;
    double mach;
    double area_ratio;
    std::vector<double> state; // Cantera thermo state
    std::vector<double> damkohler; // per-species Damkohler
    double Da_min; // Minimum damkohler number
    int min_Da_species; // Index of species with lowest Damkohler
};

struct KineticNozzleResults {
    ThroatCondition throat;
    std::vector<KineticNozzleStation> stations;
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

    // Construct from Gas directly (chemistry is overridden to FROZEN internally)
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
