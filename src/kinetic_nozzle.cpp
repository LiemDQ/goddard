#include <memory>
#include <limits>
#include "cantera/zerodim.h"
#include "cantera/numerics/Func1.h"
#include "goddard/kinetic_nozzle.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/error.hpp"
#include "goddard/utils.hpp"

namespace Goddard {

class VelocityFunc: public Cantera::Func1 {
public:
    double& m_velocity;
    VelocityFunc(double& vel) : Func1(), m_velocity(vel) {}
    double eval(double t) const override { return m_velocity;}
};

KineticNozzle::KineticNozzle(
    Cantera::Solution& gas, 
    NozzleProfile& profile, 
    double mass_flow_rate, 
    NozzleChemistryType chemistry) 
: m_gas(gas.shared_from_this()), m_profile(profile), m_mdot(mass_flow_rate)
{
    m_gas->thermo()->restoreState(m_inlet_state);
    switch (chemistry) {
        case NozzleChemistryType::EQUILIBRIUM: {
            m_throat_solver = std::make_unique<NozzleBase>(new EquilibriumNozzle(gas, m_inlet_state));
            break;
        }
        case NozzleChemistryType::FROZEN: {
            m_throat_solver = std::make_unique<NozzleBase>(new FrozenNozzle(gas, m_inlet_state));
        }
        default:
            throw std::runtime_error("Chemistry type not supported for kinetic nozzle.");
    }
}

KineticNozzle::KineticNozzle(
        Cantera::Solution& gas, 
        NozzleProfile& profile, 
        double mass_flow_rate, 
        std::vector<double> inlet_state,
        NozzleChemistryType chemistry = NozzleChemistryType::EQUILIBRIUM) 
: m_gas(gas.shared_from_this()), m_profile(profile), m_mdot(mass_flow_rate), m_inlet_state(inlet_state)
{
    switch (chemistry) {
        case NozzleChemistryType::EQUILIBRIUM: {
            m_throat_solver = std::make_unique<NozzleBase>(new EquilibriumNozzle(gas, inlet_state));
            break;
        }
        case NozzleChemistryType::FROZEN: {
            m_throat_solver = std::make_unique<NozzleBase>(new FrozenNozzle(gas, inlet_state));
        }
        default:
            throw std::runtime_error("Chemistry type not supported for kinetic nozzle.");
    }
}

double KineticNozzle::get_gamma_s(Cantera::ThermoPhase& state) {
    // Frozen gamma_s is the correct value to use here, 
    // because the composition is thermodynamically "decoupled" from the
    // expansion and sound wave propagation.
    return state.cp_mass()/state.cv_mass();
}

KineticNozzleResults KineticNozzle::solve(double dt_max, int max_steps) {
    ThroatCondition throat = m_throat_solver->solve_throat_conditions();
    auto thermo = m_gas->thermo();
    thermo->restoreState(throat.state);

    double H0 = throat.H_stagnation;
    double x = m_profile.x_min();
    double u = gas_isenthalpic_velocity(*thermo, H0);
    double M = 1.0; // mach number is 1.0 by definition in the throat.
    double A_throat = m_profile.area_at(x); 

    double rho = thermo->density();
    double V_init = m_mdot / rho;

    auto reactor = Cantera::IdealGasMoleReactor{};
    reactor.insert(m_gas);
    reactor.setInitialVolume(V_init);
    reactor.setEnergy(true); // adiabatic

    auto reservoir = Cantera::Reservoir{};

    double wall_velocity = 0.0;
    auto vfunc = std::make_unique<VelocityFunc>(wall_velocity);
    auto wall = std::make_shared<Cantera::Wall>();
    wall->install(reactor, reservoir);
    wall->setArea(1.0);
    wall->setVelocity(vfunc.get());

    Cantera::ReactorNet net;
    net.addReactor(reactor);

    // production rates
    std::vector<double> wdot(m_gas->kinetics()->nTotalSpecies());
    m_gas->kinetics()->getNetProductionRates(wdot.data());
    
    // concentrations
    const size_t n_species = m_gas->thermo()->nSpecies();
    std::vector<double> conc(n_species);
    m_gas->thermo()->getConcentrations(conc.data());

    std::vector<double> damkohler(n_species);

    std::vector<KineticNozzleStation> stations;
    double t = 0.0;
    
    const double conc_threshold = 1e-10; // skip trace species in Damkohler calculations
    const double x_exit = m_profile.x_max();
    
    
    for (int step = 0; step < max_steps; step++) {
        if (x >= x_exit) break;

        double r = m_profile.radius_at(x);
        double drdx = m_profile.slope_at(x);
        double A_At = m_profile.area_at(x)/A_throat;
        double V = reactor.volume();
        
        double dx_remain = x_exit - x;
        double dt = std::min(dt_max, 0.5*dx_remain/ u);

        // Strang splitting: expansion -> chemistry -> expansion
        x += u * dt * 0.5;
        r = m_profile.radius_at(x);
        drdx = m_profile.slope_at(x);
        wall_velocity = 2.0 * V * u * drdx / r; //updating wall velocity updates the functor

        net.advance(t + dt);
        t += dt;
        
        u = gas_isenthalpic_velocity(*thermo, H0);
        M = u / gas_sonic_velocity(*thermo, get_gamma_s(*thermo));
        V = reactor.volume();
        x += u * dt / 2.0;

        r = m_profile.radius_at(x);
        drdx = m_profile.slope_at(x);
        wall_velocity = 2.0 * V * u * drdx / r; 

        double dx = u*dt;
        double Da_min = std::numeric_limits<double>::max();
        int freeze_idx = 0;
        m_gas->kinetics()->getNetProductionRates(wdot.data());
        m_gas->thermo()->getConcentrations(conc.data());

        for (size_t k = 0; k < n_species; k++) {
            if (conc[k] < conc_threshold) {
                damkohler[k] = 0.0;
                continue;
            }
            damkohler[k] = std::abs(wdot[k]) * dx /(u * conc[k]);
            if (damkohler[k] < Da_min && damkohler[k] > 0.0) {
                Da_min = damkohler[k];
                freeze_idx = static_cast<int>(k);
            }
        }

        
        stations.push_back({x, u, M, A_At, save_thermo_state(*thermo), damkohler, Da_min, freeze_idx});
    }
    return {throat, stations};
}



} // namespace Goddard