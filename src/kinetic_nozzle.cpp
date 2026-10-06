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

namespace {

/**
 * Reject a mixture carrying condensed products.
 *
 * The finite-rate reactor network integrates the gas phase alone, so a solid or liquid product
 * would silently be dropped from the mass and energy balance.
 *
 * @param gas Mixture handed to the kinetic nozzle.
 * @throws NotImplementedError if any condensed species is present in a nonzero amount.
 */
void check_no_condensed_phases(const Gas& gas) {
    if (gas.has_condensed_phases()) {
        throw NotImplementedError(
            "KineticNozzle: finite-rate chemistry with condensed species is not implemented.");
    }
}

} // namespace

class VelocityFunc: public Cantera::Func1 {
public:
    double& m_velocity;
    VelocityFunc(double& vel) : Func1(), m_velocity(vel) {}
    double eval(double) const override { return m_velocity;}
};

// Solution& constructors: build Gas internally
KineticNozzle::KineticNozzle(
    Cantera::Solution& gas,
    NozzleProfile& prof,
    double mass_flow_rate,
    NozzleOptions options)
: profile(prof), mdot(mass_flow_rate), opts(options), m_throat_solver(gas, options),
  m_gas(gas, GasChemistry::FROZEN)
{
    m_inlet_state.resize(gas.thermo()->stateSize());
    gas.thermo()->saveState(m_inlet_state);
}

KineticNozzle::KineticNozzle(
        Cantera::Solution& gas,
        NozzleProfile& prof,
        double mass_flow_rate,
        std::vector<double> inlet_state,
        NozzleOptions options)
: profile(prof), mdot(mass_flow_rate), opts(options), m_throat_solver(gas, inlet_state, options),
  m_inlet_state(std::move(inlet_state)), m_gas(gas, GasChemistry::FROZEN)
{
}

// Gas constructors: use provided Gas, override chemistry to FROZEN. `options.chemistry` selects
// the throat model; the throat Nozzle rejects KINETIC.
KineticNozzle::KineticNozzle(
    const Gas& gas,
    NozzleProfile& prof,
    double mass_flow_rate,
    NozzleOptions options)
: profile(prof), mdot(mass_flow_rate), opts(options),
  m_throat_solver(gas, options),
  m_gas(gas)
{
    check_no_condensed_phases(m_gas);
    m_gas.chemistry = GasChemistry::FROZEN;
    m_inlet_state.resize(m_gas.thermo()->stateSize());
    m_gas.thermo()->saveState(m_inlet_state);
}

KineticNozzle::KineticNozzle(
    const Gas& gas,
    NozzleProfile& prof,
    double mass_flow_rate,
    std::vector<double> inlet_state,
    NozzleOptions options)
// m_throat_solver is initialized before m_inlet_state (declaration order), so it receives a copy
// of `inlet_state` before the move.
: profile(prof), mdot(mass_flow_rate), opts(options),
  m_throat_solver(gas, inlet_state, options),
  m_inlet_state(std::move(inlet_state)),
  m_gas(gas)
{
    m_gas.chemistry = GasChemistry::FROZEN;
    m_gas.restore_state(m_inlet_state);
    check_no_condensed_phases(m_gas);
    m_gas.set_current_state_as_reference();
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

KineticNozzleResults KineticNozzle::solve(double dt_max, double dx_max, int max_steps) {
    ThroatCondition throat = m_throat_solver.solve_throat_conditions();
    auto thermo = m_gas.thermo();
    if (throat.state.size() > thermo->stateSize()) {
        // An extended state vector carries condensed amounts the reactor network cannot follow.
        m_gas.restore_state(throat.state);
        check_no_condensed_phases(m_gas);
    }
    thermo->restoreState(thermo->stateSize(), throat.state.data());

    double H0 = throat.H_stagnation;
    double x = profile.x_min();
    double u = m_gas.isenthalpic_velocity(H0);
    double dudx = 0.0;
    double M = 1.0; // mach number is 1.0 by definition in the throat.
    double A_throat = profile.area_at(x);

    double rho = thermo->density();
    double V_init = mdot / rho;

    auto reactor = std::make_shared<Cantera::IdealGasMoleReactor>(m_gas.solution(), false);
    reactor->setInitialVolume(V_init);
    reactor->setEnergyEnabled(true); // adiabatic
    reactor->syncState();

    auto reservoir = std::make_shared<Cantera::Reservoir>(m_gas.solution(), false);

    double wall_velocity = 0.0;
    auto vfunc = std::make_shared<VelocityFunc>(wall_velocity);
    auto wall = std::make_shared<Cantera::Wall>(reactor, reservoir);
    wall->setArea(1.0);
    wall->setVelocity(vfunc);
    wall->setExpansionRateCoeff(0.0);

    Cantera::ReactorNet net(reactor);
    net.setInitialTime(0.0);

    // production rates
    std::vector<double> wdot(m_gas.kinetics()->nTotalSpecies());
    m_gas.kinetics()->getNetProductionRates(wdot.data());

    // concentrations
    const size_t n_species = m_gas.thermo()->nSpecies();
    std::vector<double> conc(n_species);
    m_gas.thermo()->getConcentrations(conc.data());

    std::vector<double> damkohler(n_species);

    std::vector<KineticNozzleStation> stations;
    double t = 0.0;
    double u_old = u;
    double x_old = x;

    const double conc_threshold = 1e-10; // skip trace species in Damkohler calculations
    const double x_exit = profile.x_max();


    for (int step = 0; step < max_steps; step++) {
        if (x >= x_exit) break;

        // A step that can cover the remaining length ends exactly at the exit; otherwise rounding
        // can leave x a few ulps short and the loop would crawl towards x_exit until max_steps.
        // The slack lets that step stretch slightly rather than leave a sliver for one more step.
        const double slack = 1.0 + 1e-6;
        const double remaining = x_exit - x;
        const bool reaches_exit = remaining <= slack*dx_max && remaining <= slack*u*dt_max;
        double dt = reaches_exit ? remaining/u : std::min(dt_max, dx_max/u);

        x += u * dt * 0.5;

        double r = profile.radius_at(x);
        double drdx = profile.slope_at(x);
        double A_At = profile.area_at(x)/A_throat;
        double V = reactor->volume();

        dudx = (u - u_old)/(x-x_old);
        wall_velocity = V* (2.0 * u * drdx / r + dudx); //updating wall velocity updates the functor

        try {
            net.advance(t + dt);
        } catch (const Cantera::CanteraError& e) {
            throw std::runtime_error(diagnose_cvodes_failure(e, x, step, t, dt_max));
        }
        t += dt;

        x_old = x;
        u_old = u;
        x += u * dt * 0.5;
        if (reaches_exit) {
            x = x_exit;
        }

        double dx = u*dt;
        double Da_min = std::numeric_limits<double>::max();
        int freeze_idx = 0;
        m_gas.kinetics()->getNetProductionRates(wdot.data());
        m_gas.thermo()->getConcentrations(conc.data());

        u = m_gas.isenthalpic_velocity(H0);
        M = u / m_gas.speed_of_sound();

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
        ThermodynamicState station_thermo = m_gas.snapshot();
        station_thermo.stagnation_enthalpy = H0;
        const std::string freezing_species = (Da_min < std::numeric_limits<double>::max())
            ? m_gas.thermo()->speciesName(static_cast<size_t>(freeze_idx)) : std::string();
        stations.push_back({x, u, M, A_At, std::move(station_thermo), save_thermo_state(*thermo),
            damkohler, Da_min, freezing_species});
    }
    return {throat, stations, x >= x_exit};
}

} // namespace Goddard
