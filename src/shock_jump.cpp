#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "cantera/core.h"
#include "goddard/error.hpp"
#include "goddard/newton.hpp"
#include "goddard/shock_jump.hpp"

namespace Goddard {

namespace {

double normal_shock_control_factor(int iter) {
    double cf;
    if (iter > 20) {
        cf = 0.04879016;
    }
    else if (iter > 12) {
        cf = 0.09531018;
    }
    else if (iter > 4) {
        cf = 0.22314355;
    }
    else {
        cf = 0.40546511;
    }
    return cf;
}

/** Set a trial downstream state: frozen composition, or equilibrium at (T, P), per `gas.chemistry`. */
void set_trial_state(Gas& gas, double T, double P) {
    switch (gas.chemistry) {
        case GasChemistry::EQUILIBRIUM:
            gas.equilibrate_TP(T, P);
            break;
        case GasChemistry::PERFECT_GAS:
        case GasChemistry::FROZEN:
        case GasChemistry::KINETIC:
            gas.set_state_TP(T, P);
            break;
    }
}

} // anonymous namespace

ShockResult invalid_shock() {
    ShockResult result;
    result.valid = false;
    result.static_pressure_ratio = -1.0;
    result.static_temperature_ratio = -1.0;
    result.density_ratio = -1.0;
    result.total_pressure_ratio = -1.0;
    result.mach_in = -1.0;
    result.mach_out = -1.0;

    return result;
}

ShockJump solve_shock_jump(Gas& gas, JumpCondition condition, double velocity,
                           double pressure_ratio_guess, double temperature_ratio_guess,
                           const SolverOptions& opts)
{
    if (condition == JumpCondition::CHAPMAN_JOUGUET && gas.chemistry != GasChemistry::EQUILIBRIUM) {
        throw std::invalid_argument(
            "solve_shock_jump: a Chapman-Jouguet detonation requires EQUILIBRIUM chemistry.");
    }
    const double R = Cantera::GasConstant;
    const double h_up = gas.enthalpy_mass();
    const double P_up = gas.pressure();
    const double T_up = gas.temperature();
    const double rho_up = gas.density();
    const double mw_up = gas.molecular_weight();

    // M u^2/(R T) = rho u^2/P, and u^2/R, of the upstream gas
    const double momentum_coeff = mw_up*velocity*velocity/(R*T_up);
    const double energy_coeff = velocity*velocity/R;

    // x = (ln P_down/P_up, ln T_down/T_up)
    const Eigen::Vector2d guess(std::log(pressure_ratio_guess), std::log(temperature_ratio_guess));
    NewtonOptions options;
    options.step_abstol = opts.abstol;
    options.max_iterations = opts.max_iterations + 1;
    const NewtonSystemResult<2> result = newton_solve(guess, [&](const Eigen::Vector2d& x) {
        set_trial_state(gas, std::exp(x(1))*T_up, std::exp(x(0))*P_up);
        const ExpansionProperties props = gas.expansion_properties();
        const double dlogV_dlogT = props.dlV_dlT_P;
        const double dlogV_dlogP = props.dlV_dlP_T;
        const double pressure_ratio = std::exp(x(0));
        const double T_down = gas.temperature();
        const double enthalpy_rise = (gas.enthalpy_mass() - h_up)/R;
        // (1/R)(dh/d ln P)_T = (T/M)(1 - dlnV/dlnT), and (1/R)(dh/d ln T)_P = T cp/R
        const double dh_dlogP = T_down/gas.molecular_weight()*(1.0 - dlogV_dlogT);
        const double dh_dlogT = T_down*props.spec_heat_p/R;

        // residuals f and Jacobian J = df/dx, in the sign convention of docs/theory/shocks.md
        double f_P = 0.0, f_h = 0.0;
        double J_PP = 0.0, J_PT = 0.0, J_hP = 0.0, J_hT = 0.0;
        switch (condition) {
            case JumpCondition::INCIDENT: {
                const double r = rho_up/gas.density();  // rho1/rho2
                f_P = pressure_ratio - 1.0 + momentum_coeff*(r - 1.0);
                f_h = enthalpy_rise - 0.5*energy_coeff*(1.0 - r*r);
                J_PP = pressure_ratio + momentum_coeff*r*dlogV_dlogP;
                J_PT = momentum_coeff*r*dlogV_dlogT;
                J_hP = dh_dlogP + energy_coeff*r*r*dlogV_dlogP;
                J_hT = dh_dlogT + energy_coeff*r*r*dlogV_dlogT;
                break;
            }
            case JumpCondition::REFLECTED: {
                const double r = gas.density()/rho_up;  // rho5/rho2
                const double c = r/((r - 1.0)*(r - 1.0));
                f_P = pressure_ratio - 1.0 - momentum_coeff*r/(r - 1.0);
                f_h = enthalpy_rise - 0.5*energy_coeff*(r + 1.0)/(r - 1.0);
                J_PP = pressure_ratio - momentum_coeff*c*dlogV_dlogP;
                J_PT = -momentum_coeff*c*dlogV_dlogT;
                J_hP = dh_dlogP - energy_coeff*c*dlogV_dlogP;
                J_hT = dh_dlogT - energy_coeff*c*dlogV_dlogT;
                break;
            }
            case JumpCondition::CHAPMAN_JOUGUET: {
                // RP-1311 eqs. (8.1)-(8.12): u2 = a2 eliminates the velocity. gamma_s is held
                // fixed in the Jacobian, as in RP-1311; see docs/theory/detonations.md.
                const double r = gas.density()/rho_up;  // rho2/rho1
                const double gamma_s = props.gamma_s;
                const double sonic = gamma_s*T_down/(2.0*gas.molecular_weight());  // a2^2/(2R)
                f_P = 1.0/pressure_ratio - 1.0 + gamma_s*(r - 1.0);
                f_h = enthalpy_rise - sonic*(r*r - 1.0);
                J_PP = -1.0/pressure_ratio - gamma_s*r*dlogV_dlogP;
                J_PT = -gamma_s*r*dlogV_dlogT;
                J_hP = dh_dlogP - sonic*((r*r - 1.0) - (r*r + 1.0)*dlogV_dlogP);
                J_hT = dh_dlogT + sonic*(r*r + 1.0)*dlogV_dlogT;
                break;
            }
        }

        NewtonSystemFunction<2> function;
        function.value << f_P, f_h;
        function.jacobian << J_PP, J_PT,
                             J_hP, J_hT;
        return function;
    }, [](int k, const Eigen::Vector2d&, const Eigen::Vector2d& step) -> Eigen::Vector2d {
        // one factor for both corrections, limiting the larger one as in RP-1311
        return step*std::min(1.0, normal_shock_control_factor(k)/step.cwiseAbs().maxCoeff());
    }, options);

    const double largest_step = result.step.cwiseAbs().maxCoeff();
    switch (result.status) {
        case NewtonStatus::CONVERGED:
            break;
        case NewtonStatus::MAX_ITERATIONS: {
            const char* message = condition == JumpCondition::CHAPMAN_JOUGUET
                ? "Chapman-Jouguet detonation failed to converge."
                : "Normal shock properties failed to converge.";
            throw ConvergenceError(message, result.iterations, opts.abstol, largest_step);
        }
        case NewtonStatus::NON_FINITE_STEP:
            throw ConvergenceError("Shock jump Newton step is not finite.",
                result.iterations, opts.abstol, largest_step);
    }
    set_trial_state(gas, std::exp(result.x(1))*T_up, std::exp(result.x(0))*P_up);

    ShockJump jump;
    jump.pressure_ratio = std::exp(result.x(0));
    jump.temperature_ratio = std::exp(result.x(1));
    jump.density_ratio = gas.density()/rho_up;
    return jump;
}

ShockJump solve_shock_jump_with_chemistry(Gas& gas, JumpCondition condition, double velocity,
                                          GasChemistry chemistry,
                                          double pressure_ratio_guess, double temperature_ratio_guess,
                                          const SolverOptions& opts)
{
    const std::vector<double> upstream_state = gas.save_state();
    gas.chemistry = GasChemistry::FROZEN;
    ShockJump jump = solve_shock_jump(gas, condition, velocity, pressure_ratio_guess, temperature_ratio_guess, opts);
    if (chemistry == GasChemistry::EQUILIBRIUM) {
        gas.restore_state(upstream_state);
        gas.chemistry = GasChemistry::EQUILIBRIUM;
        jump = solve_shock_jump(gas, condition, velocity, jump.pressure_ratio, jump.temperature_ratio, opts);
    }
    gas.chemistry = chemistry;
    return jump;
}

double frozen_sound_speed(Gas& gas) {
    const GasChemistry chemistry = gas.chemistry;
    gas.chemistry = GasChemistry::FROZEN;
    const double sound_speed = gas.speed_of_sound();
    gas.chemistry = chemistry;
    return sound_speed;
}

ShockJump incident_jump(Gas& gas, double u1, GasChemistry chemistry, const SolverOptions& opts) {
    const double mach1 = u1/frozen_sound_speed(gas);
    const ShockResult guess = Goddard::normal_shock(mach1, gas.gamma());
    return solve_shock_jump_with_chemistry(gas, JumpCondition::INCIDENT, u1, chemistry,
        guess.static_pressure_ratio, guess.static_temperature_ratio, opts);
}

ShockResult incident_shock(Gas& gas, double u_normal, double u_tangential,
                           GasChemistry chemistry, const SolverOptions& opts)
{
    const double a1 = frozen_sound_speed(gas);
    // A wave at the Mach angle has u_normal = a1 up to round-off; treat it as sonic.
    if (u_normal < a1*(1.0 - 1e-12)) {
        return invalid_shock();
    }
    u_normal = std::max(u_normal, a1);
    gas.chemistry = GasChemistry::FROZEN;
    const double P01 = gas.stagnation_pressure(std::hypot(u_normal, u_tangential));

    const ShockJump jump = incident_jump(gas, u_normal, chemistry, opts);
    const double u_normal2 = u_normal/jump.density_ratio;

    ShockResult result;
    result.valid = true;
    result.mach_in = u_normal/a1;
    result.mach_out = u_normal2/gas.speed_of_sound();
    result.static_pressure_ratio = jump.pressure_ratio;
    result.static_temperature_ratio = jump.temperature_ratio;
    result.density_ratio = jump.density_ratio;
    result.total_pressure_ratio = gas.stagnation_pressure(std::hypot(u_normal2, u_tangential))/P01;
    return result;
}

ShockResult reflected_shock_from_state2(Gas& gas, double particle_velocity, GasChemistry chemistry,
                                        double pressure_ratio_guess, double temperature_ratio_guess,
                                        const SolverOptions& opts)
{
    const std::vector<double> state2 = gas.save_state();
    const GasChemistry incident_chemistry = gas.chemistry;
    const double a2 = gas.speed_of_sound();

    const ShockJump jump = solve_shock_jump_with_chemistry(gas, JumpCondition::REFLECTED, particle_velocity,
        chemistry, pressure_ratio_guess, temperature_ratio_guess, opts);
    const double wave_speed = particle_velocity/(jump.density_ratio - 1.0);
    const double u_in = particle_velocity + wave_speed;

    ShockResult result;
    result.valid = true;
    result.mach_in = u_in/a2;
    result.mach_out = wave_speed/gas.speed_of_sound();
    result.static_pressure_ratio = jump.pressure_ratio;
    result.static_temperature_ratio = jump.temperature_ratio;
    result.density_ratio = jump.density_ratio;

    const double P05 = gas.stagnation_pressure(wave_speed);
    const std::vector<double> state5 = gas.save_state();
    gas.restore_state(state2);
    gas.chemistry = incident_chemistry;
    const double P02 = gas.stagnation_pressure(u_in);
    gas.restore_state(state5);
    gas.chemistry = chemistry;
    result.total_pressure_ratio = P05/P02;
    return result;
}

} // namespace Goddard
