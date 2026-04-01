#pragma once
#include <memory>
#include <vector>
#include "goddard/nozzle.hpp"
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
    
    KineticNozzle(
        Cantera::Solution& gas, 
        NozzleProfile& profile, 
        double mass_flow_rate, 
        NozzleChemistryType chemistry = NozzleChemistryType::EQUILIBRIUM);

    KineticNozzle(
        Cantera::Solution& gas, 
        NozzleProfile& profile, 
        double mass_flow_rate, 
        std::vector<double> inlet_state,
        NozzleChemistryType chemistry = NozzleChemistryType::EQUILIBRIUM);

    KineticNozzleResults solve(double dt_max = 1e-4, int max_steps = 100000);

    double get_gamma_s(Cantera::ThermoPhase& state);
    
    NozzleProfile m_profile;
    double m_mdot;

protected:
    std::unique_ptr<NozzleBase> m_throat_solver;
    std::vector<double> m_inlet_state;
    std::shared_ptr<Cantera::Solution> m_gas;
};

} // namespace Goddard