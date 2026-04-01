#include <format>
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
    double eval(double) const override { return m_velocity;}
};

KineticNozzle::KineticNozzle(
    Cantera::Solution& gas,
    NozzleProfile& profile,
    double mass_flow_rate,
    NozzleChemistryType chemistry)
: m_profile(profile), m_mdot(mass_flow_rate), m_gas(gas.shared_from_this())
{
    m_inlet_state.resize(gas.thermo()->stateSize());
    gas.thermo()->saveState(m_inlet_state);
    switch (chemistry) {
        case NozzleChemistryType::EQUILIBRIUM:
            m_throat_solver = std::make_unique<EquilibriumNozzle>(gas, m_inlet_state);
            break;
        case NozzleChemistryType::FROZEN:
            m_throat_solver = std::make_unique<FrozenNozzle>(gas, m_inlet_state);
            break;
        default:
            throw std::runtime_error("Chemistry type not supported for kinetic nozzle throat model.");
    }
}

KineticNozzle::KineticNozzle(
        Cantera::Solution& gas,
        NozzleProfile& profile,
        double mass_flow_rate,
        std::vector<double> inlet_state,
        NozzleChemistryType chemistry)
: m_profile(profile), m_mdot(mass_flow_rate), m_inlet_state(inlet_state), m_gas(gas.shared_from_this())
{
    switch (chemistry) {
        case NozzleChemistryType::EQUILIBRIUM:
            m_throat_solver = std::make_unique<EquilibriumNozzle>(gas, inlet_state);
            break;
        case NozzleChemistryType::FROZEN:
            m_throat_solver = std::make_unique<FrozenNozzle>(gas, inlet_state);
            break;
        default:
            throw std::runtime_error("Chemistry type not supported for kinetic nozzle throat model.");
    }
}

double KineticNozzle::get_gamma_s(Cantera::ThermoPhase& state) {
    // Frozen gamma_s is the correct value to use here,
    // because the composition is thermodynamically "decoupled" from the
    // expansion and sound wave propagation.
    return state.cp_mass()/state.cv_mass();
}

// Parse a CanteraError from CVodes and return a diagnostic message with context.
// CVodes error code -3 (CV_ERR_FAILURE) means the error test failed repeatedly at
// |h| = hmin — a stiffness/stability issue the user can address by reducing dt_max.
// Other codes indicate different failures and get a generic context message.
static std::string diagnose_cvodes_failure(
    const Cantera::CanteraError& e, double x, int step, double t, double dt_max)
{
    const char* msg = e.what();
    bool is_stiffness_failure = std::string_view(msg).find("Error code: -3") != std::string_view::npos
                             || std::string_view(msg).find("hmin") != std::string_view::npos;

    if (is_stiffness_failure) {
        return std::vformat(
                   "KineticNozzle::solve: stiffness limit reached at x={:.4f} m "
                   "(step {}, t={:.3e} s). The coupled chemistry+expansion ODE is too "
                   "stiff for CVodes at dt_max={:.1e} s. Reduce dt_max by an order of magnitude.\n"
                   "CVodes error: ",
                   std::make_format_args(x, step, t, dt_max))
               + msg;
    }
    return std::vformat(
               "KineticNozzle::solve: CVodes integration failed at x={:.4f} m "
               "(step {}, t={:.3e} s).\nCVodes error: ",
               std::make_format_args(x, step, t))
           + msg;
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

    auto reactor = std::make_shared<Cantera::IdealGasMoleReactor>(m_gas, false);
    reactor->setInitialVolume(V_init);
    reactor->setEnergyEnabled(true); // adiabatic
    reactor->syncState();

    auto reservoir = std::make_shared<Cantera::Reservoir>(m_gas, false);

    double wall_velocity = 0.0;
    auto vfunc = std::make_shared<VelocityFunc>(wall_velocity);
    auto wall = std::make_shared<Cantera::Wall>(reactor, reservoir);
    wall->setArea(1.0);
    wall->setVelocity(vfunc);
    wall->setExpansionRateCoeff(0.0);

    Cantera::ReactorNet net(reactor);
    net.setInitialTime(0.0);

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
        double V = reactor->volume();

        double dx_remain = x_exit - x;
        if (dx_remain < 1e-6) break; // avoid Zeno's paradox: dt cap = 0.5*dx/u → never reaches exit
        double dt = std::min(dt_max, 0.5*dx_remain/ u);

        // "Strang" splitting: expansion -> chemistry -> expansion
        // This is technically not Strang splitting as the volume ODE is coupled 
        // to the chemistry ODEs via Cantera's IdealGasMoleReactor.
        x += u * dt * 0.5;
        r = m_profile.radius_at(x);
        drdx = m_profile.slope_at(x);
        wall_velocity = 2.0 * V * u * drdx / r; //updating wall velocity updates the functor

        try {
            net.advance(t + dt);
        } catch (const Cantera::CanteraError& e) {
            throw std::runtime_error(diagnose_cvodes_failure(e, x, step, t, dt_max));
        }
        t += dt;

        u = gas_isenthalpic_velocity(*thermo, H0);
        M = u / gas_sonic_velocity(*thermo, get_gamma_s(*thermo));
        V = reactor->volume();
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
