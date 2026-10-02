#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>
#include "cantera/core.h"
#include "goddard/error.hpp"
#include "goddard/shocks.hpp"
#include "goddard/gas_dynamics.hpp"

namespace Goddard {

// ===== utilities =====
namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

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

ObliqueShockResult invalid_oblique_shock() {
    ObliqueShockResult result;
    result.valid = false;
    result.mach_in = -1.0;
    result.mach_out = -1.0;
    result.shock = invalid_shock();
    result.beta = -1.0;
    result.theta = -1.0;
    return result;
}

ReflectedShockResult invalid_reflected_shock() {
    ReflectedShockResult result;
    result.valid = false;
    result.incident = invalid_shock();
    result.reflected = invalid_shock();
    return result;
}

} // namespace

// ===== Perfect gas free functions =====

ShockResult normal_shock(double mach, double gamma) {

    ShockResult result;
    // check if input values are valid
    if (mach < 1.0 || gamma <= 1.0){
        return invalid_shock();
    }

    //mach numbers
    double M1 = mach;
    double M2 = std::sqrt( (M1*M1 * (gamma - 1) + 2)/(2 *  gamma * M1 * M1 - (gamma - 1)));
    result.mach_in = M1;
    result.mach_out = M2;

    //dynamic pressures
    result.static_pressure_ratio = (2 * gamma * M1 * M1)/(gamma + 1) - (gamma - 1)/(gamma + 1);

    //static temperatures
    result.static_temperature_ratio = ((1 + (gamma - 1)/2*M1*M1)
        * (2*gamma/(gamma -1) * M1*M1 - 1)
        / (M1*M1 * (2 * gamma /(gamma - 1) + (gamma -1)/2)));

    result.density_ratio = (gamma + 1) * M1 * M1 / ((gamma - 1) * M1 * M1 + 2);

    //stagnation pressures
    result.total_pressure_ratio = std::pow(((gamma +1)/2 * M1 * M1)/(1 + (gamma-1)/2* M1 * M1), gamma/(gamma-1))
        * std::pow(1.0/(2*gamma/(gamma+1)*M1*M1 - (gamma-1)/(gamma+1)), 1/(gamma-1));

    result.valid = true;
    return result;
}

ReflectedShockResult reflected_shock(double mach, double gamma) {
    if (mach < 1.0 || gamma <= 1.0) {
        return invalid_reflected_shock();
    }

    // A sonic incident wave sets the gas in motion at zero speed, so its reflection is sonic too.
    double mach_R = 1.0;
    if (mach > 1.0) {
        double MR_relation = mach/(mach*mach - 1)
            * std::sqrt(1 + 2*(gamma-1)/((gamma+1)*(gamma+1))*(mach*mach -1)*(gamma + 1/(mach*mach)));
        mach_R = (1 + std::sqrt(1 + 4*MR_relation*MR_relation))/(2*MR_relation);
    }

    ReflectedShockResult result;
    result.incident = normal_shock(mach, gamma);
    result.reflected = normal_shock(mach_R, gamma);
    result.valid = result.incident.valid && result.reflected.valid;
    return result;
}

double oblique_shock_max_deflection_wave_angle(double mach, double gamma) {
    if (mach < 1.0 || gamma <= 1.0) {
        return NaN;
    }
    const double M2 = mach*mach;
    const double sin2_beta = ((gamma+1)*M2 - 4
        + std::sqrt((gamma+1)*((gamma+1)*M2*M2 + 8*(gamma-1)*M2 + 16)))
        / (4*gamma*M2);
    return std::asin(std::sqrt(std::min(sin2_beta, 1.0)));
}

double oblique_shock_max_deflection(double mach, double gamma) {
    const double beta_max = oblique_shock_max_deflection_wave_angle(mach, gamma);
    if (std::isnan(beta_max)) {
        return NaN;
    }
    return oblique_shock_deflection_angle(mach, beta_max, gamma);
}

std::pair<double, double> oblique_shock_wave_angle(double mach, double deflection_angle, double gamma) {
    if (mach < 1.0 || gamma <= 1.0 || deflection_angle < 0.0
        || deflection_angle > oblique_shock_max_deflection(mach, gamma)) {
        return {NaN, NaN};
    }
    if (deflection_angle == 0.0) {
        return {std::asin(1.0/mach), M_PI/2};
    }

    double tan_theta = std::tan(deflection_angle);
    double M2 = mach*mach;
    double stagnation = stagnation_factor(mach, gamma);
    // At the maximum deflection lambda^2 and 1 - chi are zero; clamp round-off below that.
    double lambda = std::sqrt(std::max(0.0,
        (M2-1)*(M2-1) - 3*stagnation*(1 + (gamma+1)/2.0 * M2)*tan_theta*tan_theta));
    double chi = (std::pow(M2-1, 3.0)
        - 9*stagnation*(stagnation+(gamma+1)/4.0*M2*M2)
        *tan_theta*tan_theta)
        / std::pow(lambda, 3.0);
    chi = std::max(-1.0, std::min(chi, 1.0));

    double weak_beta_cos = std::cos((4*M_PI + std::acos(chi))/3);
    double strong_beta_cos = std::cos((std::acos(chi))/3);

    double weak_beta = std::atan((M2-1 + 2*lambda*weak_beta_cos)
        / (3 * stagnation * tan_theta));

    double strong_beta = std::atan((M2-1 + 2*lambda*strong_beta_cos)
        / (3 * stagnation * tan_theta));

    return {weak_beta, strong_beta};
}

double oblique_shock_deflection_angle(double mach, double wave_angle, double gamma) {
    double sin_beta = std::sin(wave_angle);
    double cos_2beta = std::cos(2*wave_angle);
    double cot_beta = 1.0/std::tan(wave_angle);
    return std::atan(2.0*cot_beta
        *(mach*mach*sin_beta*sin_beta -1)
        /((mach*mach )*(gamma + cos_2beta)+2));
}

ObliqueShockResult oblique_shock_from_deflection(
    double mach, double deflection_angle, double gamma, bool weak)
{
    auto [weak_beta, strong_beta] = oblique_shock_wave_angle(mach, deflection_angle, gamma);
    const double beta = weak ? weak_beta : strong_beta;
    if (std::isnan(beta)) {
        return invalid_oblique_shock();
    }

    ObliqueShockResult result;
    result.beta = beta;
    result.theta = deflection_angle;
    // At zero deflection the weak wave is a Mach wave; guard M sin(mu) against round-off below 1.
    double mach_n1 = std::max(mach * std::sin(result.beta), 1.0);
    result.shock = normal_shock(mach_n1, gamma);
    if (!result.shock.valid) {
        return invalid_oblique_shock();
    }
    result.mach_in = mach;
    result.mach_out = result.shock.mach_out/std::sin(result.beta - result.theta);
    result.valid = true;
    return result;
}

ObliqueShockResult oblique_shock_from_wave_angle(
    double mach, double wave_angle, double gamma)
{
    ObliqueShockResult result;
    double mach_n1 = mach * std::sin(wave_angle);

    result.shock = normal_shock(mach_n1, gamma);
    if (!result.shock.valid) {
        return invalid_oblique_shock();
    }
    result.theta = oblique_shock_deflection_angle(mach, wave_angle, gamma);
    result.beta = wave_angle;
    result.mach_in = mach;
    result.mach_out = result.shock.mach_out/std::sin(wave_angle - result.theta);
    result.valid = true;
    return result;
}


// ===== Cantera-backed solvers (internal) =====

namespace {

/** Which pair of jump conditions `solve_shock_jump` solves. */
enum class ShockFrame {
    INCIDENT,   //!< Gas enters a stationary shock at the given velocity (RP-1311 section 7.1).
    REFLECTED   //!< Gas moving at the given velocity is brought to rest by a reflected shock (section 7.2).
};

/** Downstream-over-upstream ratios across a shock. */
struct ShockJump {
    double pressure_ratio;
    double temperature_ratio;
    double density_ratio;
};

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

/**
 * Solve the jump conditions across a normal shock by Newton's method in
 * x = (ln P_down/P_up, ln T_down/T_up), NASA RP-1311 Part I, chapter 7.
 *
 * `gas` holds the upstream state on entry and the downstream state on exit. `gas.chemistry`
 * selects frozen or equilibrium downstream properties; the derivation is in
 * instructions/shocks.md.
 *
 * @param velocity u1 for INCIDENT, the particle velocity u_p of gas 2 for REFLECTED [m/s]
 */
ShockJump solve_shock_jump(Gas& gas, ShockFrame frame, double velocity,
                           double pressure_ratio_guess, double temperature_ratio_guess,
                           const SolverOptions& opts)
{
    const double R = Cantera::GasConstant;
    const double h_up = gas.enthalpy_mass();
    const double P_up = gas.pressure();
    const double T_up = gas.temperature();
    const double rho_up = gas.density();
    const double mw_up = gas.molecular_weight();

    // M u^2/(R T) = rho u^2/P, and u^2/R, of the upstream gas
    const double momentum_coeff = mw_up*velocity*velocity/(R*T_up);
    const double energy_coeff = velocity*velocity/R;

    double log_pressure_ratio = std::log(pressure_ratio_guess);
    double log_temperature_ratio = std::log(temperature_ratio_guess);
    set_trial_state(gas, temperature_ratio_guess*T_up, pressure_ratio_guess*P_up);

    const double abstol = opts.abstol;
    double residual = 1.0;
    for (int k = 0; residual >= abstol; k++) {
        if (k > opts.max_iterations) {
            throw ConvergenceError("Normal shock properties failed to converge.", k, abstol, residual);
        }
        const ExpansionProperties props = gas.expansion_properties();
        const double dlogV_dlogT = props.dlogV_dlogT_P;
        const double dlogV_dlogP = props.dlogV_dlogP_T;
        const double pressure_ratio = std::exp(log_pressure_ratio);
        const double T_down = gas.temperature();
        const double enthalpy_rise = (gas.enthalpy_mass() - h_up)/R;
        // (1/R)(dh/d ln P)_T = (T/M)(1 - dlnV/dlnT), and (1/R)(dh/d ln T)_P = T cp/R
        const double dh_dlogP = T_down/gas.molecular_weight()*(1.0 - dlogV_dlogT);
        const double dh_dlogT = T_down*props.spec_heat_p/R;

        // residuals f and Jacobian J = df/dx, in the sign convention of instructions/shocks.md
        double f_P = 0.0, f_h = 0.0;
        double J_PP = 0.0, J_PT = 0.0, J_hP = 0.0, J_hT = 0.0;
        switch (frame) {
            case ShockFrame::INCIDENT: {
                const double r = rho_up/gas.density();  // rho1/rho2
                f_P = pressure_ratio - 1.0 + momentum_coeff*(r - 1.0);
                f_h = enthalpy_rise - 0.5*energy_coeff*(1.0 - r*r);
                J_PP = pressure_ratio + momentum_coeff*r*dlogV_dlogP;
                J_PT = momentum_coeff*r*dlogV_dlogT;
                J_hP = dh_dlogP + energy_coeff*r*r*dlogV_dlogP;
                J_hT = dh_dlogT + energy_coeff*r*r*dlogV_dlogT;
                break;
            }
            case ShockFrame::REFLECTED: {
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
        }

        // Newton step J dx = -f by Cramer's rule, limited as in RP-1311
        const double determinant = J_PP*J_hT - J_PT*J_hP;
        const double dlog_pressure = -(f_P*J_hT - J_PT*f_h)/determinant;
        const double dlog_temperature = -(J_PP*f_h - J_hP*f_P)/determinant;
        residual = std::max(std::abs(dlog_pressure), std::abs(dlog_temperature));
        if (!std::isfinite(residual)) {
            throw ConvergenceError("Normal shock Newton step is not finite.", k, abstol, residual);
        }
        const double control_factor = std::min(1.0, normal_shock_control_factor(k)/residual);

        log_pressure_ratio += control_factor*dlog_pressure;
        log_temperature_ratio += control_factor*dlog_temperature;
        set_trial_state(gas, std::exp(log_temperature_ratio)*T_up, std::exp(log_pressure_ratio)*P_up);
    }

    ShockJump jump;
    jump.pressure_ratio = std::exp(log_pressure_ratio);
    jump.temperature_ratio = std::exp(log_temperature_ratio);
    jump.density_ratio = gas.density()/rho_up;
    return jump;
}

/**
 * Solve a shock jump with the given chemistry. EQUILIBRIUM first solves the frozen jump and
 * uses it as the initial guess (instructions/shocks.md, section 1).
 *
 * `gas` holds the upstream state on entry and the downstream state on exit, with
 * `gas.chemistry` set to `chemistry`.
 */
ShockJump solve_shock_jump_with_chemistry(Gas& gas, ShockFrame frame, double velocity,
                                          GasChemistry chemistry,
                                          double pressure_ratio_guess, double temperature_ratio_guess,
                                          const SolverOptions& opts)
{
    const std::vector<double> upstream_state = gas.save_state();
    gas.chemistry = GasChemistry::FROZEN;
    ShockJump jump = solve_shock_jump(gas, frame, velocity, pressure_ratio_guess, temperature_ratio_guess, opts);
    if (chemistry == GasChemistry::EQUILIBRIUM) {
        gas.restore_state(upstream_state);
        gas.chemistry = GasChemistry::EQUILIBRIUM;
        jump = solve_shock_jump(gas, frame, velocity, jump.pressure_ratio, jump.temperature_ratio, opts);
    }
    gas.chemistry = chemistry;
    return jump;
}

/** Frozen sound speed of the gas at its current state [m/s]. */
double frozen_sound_speed(Gas& gas) {
    const GasChemistry chemistry = gas.chemistry;
    gas.chemistry = GasChemistry::FROZEN;
    const double sound_speed = gas.speed_of_sound();
    gas.chemistry = chemistry;
    return sound_speed;
}

/**
 * Jump across an incident normal shock, without stagnation quantities. `gas` holds state 1 on
 * entry and state 2 on exit. Requires u1 >= the frozen sound speed of state 1.
 */
ShockJump incident_jump(Gas& gas, double u1, GasChemistry chemistry, const SolverOptions& opts) {
    const double mach1 = u1/frozen_sound_speed(gas);
    const ShockResult guess = Goddard::normal_shock(mach1, gas.gamma());
    return solve_shock_jump_with_chemistry(gas, ShockFrame::INCIDENT, u1, chemistry,
        guess.static_pressure_ratio, guess.static_temperature_ratio, opts);
}

/**
 * Incident normal shock. `gas` holds state 1 on entry and state 2 on exit.
 *
 * State 1 is taken as given, so its sound speed and stagnation pressure are frozen ones. For an
 * oblique shock, `u_tangential` enters only the stagnation pressures.
 */
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

/**
 * Shock reflected from a closed end. `gas` holds state 2 on entry, with `gas.chemistry` the
 * chemistry of the incident shock, and state 5 on exit.
 *
 * In the reflected-shock frame gas 2 enters at u_p + W_R and gas 5 leaves at W_R, where
 * W_R = u_p/(rho5/rho2 - 1) is the reflected shock speed.
 *
 * @param particle_velocity lab-frame velocity u_p of gas 2 [m/s]
 */
ShockResult reflected_shock_from_state2(Gas& gas, double particle_velocity, GasChemistry chemistry,
                                        double pressure_ratio_guess, double temperature_ratio_guess,
                                        const SolverOptions& opts)
{
    const std::vector<double> state2 = gas.save_state();
    const GasChemistry incident_chemistry = gas.chemistry;
    const double a2 = gas.speed_of_sound();

    const ShockJump jump = solve_shock_jump_with_chemistry(gas, ShockFrame::REFLECTED, particle_velocity,
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

/** Oblique shock with wave angle `beta`. `gas` holds state 1 on entry and state 2 on exit. */
ObliqueShockResult oblique_shock_real(Gas& gas, double u1, double beta,
                                      GasChemistry chemistry, const SolverOptions& opts)
{
    const double a1 = frozen_sound_speed(gas);
    const double u_normal = u1*std::sin(beta);
    const double u_tangential = u1*std::cos(beta);

    ObliqueShockResult result;
    result.shock = incident_shock(gas, u_normal, u_tangential, chemistry, opts);
    if (!result.shock.valid) {
        return invalid_oblique_shock();
    }
    // tan(beta - theta) = u_n2/u_t = (rho1/rho2) tan(beta), for any chemistry
    const double u_normal2 = u_normal/result.shock.density_ratio;
    result.valid = true;
    result.beta = beta;
    result.theta = beta - std::atan2(u_normal2, u_tangential);
    result.mach_in = u1/a1;
    result.mach_out = std::hypot(u_normal2, u_tangential)/gas.speed_of_sound();
    return result;
}

/** Deflection of the oblique shock with wave angle `beta` from state 1 [rad]. Leaves `gas` at state 2. */
double deflection_at(Gas& gas, const std::vector<double>& state1, double u1, double beta,
                     GasChemistry chemistry, const SolverOptions& opts)
{
    gas.restore_state(state1);
    const double u_normal = u1*std::sin(beta);
    const ShockJump jump = incident_jump(gas, u_normal, chemistry, opts);
    return beta - std::atan2(u_normal/jump.density_ratio, u1*std::cos(beta));
}

struct MaxDeflection {
    double wave_angle;
    double deflection;
};

/**
 * Maximum deflection of an oblique shock, by golden-section search on the wave angle between the
 * frozen Mach angle and pi/2, where the deflection is zero. Leaves `gas` at an arbitrary state.
 */
MaxDeflection max_deflection_real(Gas& gas, const std::vector<double>& state1, double u1,
                                  double mach_angle, GasChemistry chemistry, const SolverOptions& opts)
{
    const double inverse_golden_ratio = (std::sqrt(5.0) - 1.0)/2.0;
    double lower = mach_angle;
    double upper = M_PI/2;
    double inner_lower = upper - inverse_golden_ratio*(upper - lower);
    double inner_upper = lower + inverse_golden_ratio*(upper - lower);
    double theta_lower = deflection_at(gas, state1, u1, inner_lower, chemistry, opts);
    double theta_upper = deflection_at(gas, state1, u1, inner_upper, chemistry, opts);

    for (int k = 0; upper - lower > opts.abstol; k++) {
        if (k > opts.max_iterations) {
            throw ConvergenceError("Maximum deflection search failed to converge.", k, opts.abstol, upper - lower);
        }
        if (theta_lower > theta_upper) {
            upper = inner_upper;
            inner_upper = inner_lower;
            theta_upper = theta_lower;
            inner_lower = upper - inverse_golden_ratio*(upper - lower);
            theta_lower = deflection_at(gas, state1, u1, inner_lower, chemistry, opts);
        }
        else {
            lower = inner_lower;
            inner_lower = inner_upper;
            theta_lower = theta_upper;
            inner_upper = lower + inverse_golden_ratio*(upper - lower);
            theta_upper = deflection_at(gas, state1, u1, inner_upper, chemistry, opts);
        }
    }
    if (theta_lower > theta_upper) {
        return {inner_lower, theta_lower};
    }
    return {inner_upper, theta_upper};
}

/**
 * Oblique shock with a given deflection. The deflection is zero at both ends of each branch, so
 * the bracket residuals there are known and never evaluated: in particular there is no solve at
 * a normal Mach number of 1.
 */
ObliqueShockResult oblique_shock_from_deflection_real(Gas& gas, double u1, double deflection, bool weak,
                                                      GasChemistry chemistry, const SolverOptions& opts)
{
    const double a1 = frozen_sound_speed(gas);
    if (u1 < a1 || deflection < 0.0) {
        return invalid_oblique_shock();
    }
    const std::vector<double> state1 = gas.save_state();
    const double mach_angle = std::asin(a1/u1);
    const MaxDeflection max = max_deflection_real(gas, state1, u1, mach_angle, chemistry, opts);
    if (deflection > max.deflection) {
        gas.restore_state(state1);
        gas.chemistry = chemistry;
        return invalid_oblique_shock();
    }

    if (deflection == 0.0) {
        // A Mach wave on the weak branch, a normal shock on the strong one
        gas.restore_state(state1);
        return oblique_shock_real(gas, u1, weak ? mach_angle : M_PI/2, chemistry, opts);
    }

    // residual theta(beta) - deflection at the bracket ends
    double lower = weak ? mach_angle : max.wave_angle;
    double upper = weak ? max.wave_angle : M_PI/2;
    double residual_lower = weak ? -deflection : max.deflection - deflection;

    double beta = 0.5*(lower + upper);
    for (int k = 0; ; k++) {
        if (k > opts.max_iterations) {
            throw ConvergenceError("Wave angle failed to converge.", k, opts.abstol, upper - lower);
        }
        const double residual = deflection_at(gas, state1, u1, beta, chemistry, opts) - deflection;
        if (std::abs(residual) < opts.abstol || upper - lower < 1e-12) {
            break;
        }
        if ((residual < 0.0) == (residual_lower < 0.0)) {
            lower = beta;
            residual_lower = residual;
        }
        else {
            upper = beta;
        }
        beta = 0.5*(lower + upper);
    }

    gas.restore_state(state1);
    return oblique_shock_real(gas, u1, beta, chemistry, opts);
}

/**
 * Reject a mixture carrying condensed species.
 *
 * The shock relations are gas-phase only: condensed species before or behind the shock are not
 * supported.
 *
 * @throws std::invalid_argument if any condensed species is present in a nonzero amount.
 */
void check_no_condensed_phases(const Gas& gas) {
    if (gas.has_condensed_phases()) {
        throw std::invalid_argument(
            "ShockSolver: shocks are gas-phase only; condensed species are not supported.");
    }
}

void check_shock_chemistry(GasChemistry chemistry) {
    if (chemistry != GasChemistry::FROZEN && chemistry != GasChemistry::EQUILIBRIUM) {
        throw std::invalid_argument(
            "ShockSolver: incident and reflected chemistry must be FROZEN or EQUILIBRIUM.");
    }
}

} // anonymous namespace


// ===== ShockSolver =====

ShockSolver::ShockSolver(Gas gas, SolverOptions options)
    : m_gas(std::move(gas)), m_chemistry(m_gas.chemistry), m_options(options)
{
    if (m_chemistry == GasChemistry::KINETIC) {
        throw std::invalid_argument(
            "ShockSolver: KINETIC chemistry is not supported; use FROZEN or EQUILIBRIUM.");
    }
    m_gas.copy_state(m_pre_shock_state);
    m_post_shock_state = m_pre_shock_state;
}

void ShockSolver::reset() const {
    m_gas.restore_state(m_pre_shock_state);
    m_gas.chemistry = m_chemistry;
}

double ShockSolver::pre_shock_sound_speed() const {
    reset();
    return frozen_sound_speed(m_gas);
}

void ShockSolver::store_post_shock_state() {
    m_gas.chemistry = m_chemistry;
    m_gas.copy_state(m_post_shock_state);
}

ShockResult ShockSolver::normal_shock(double mach) {
    if (m_chemistry == GasChemistry::PERFECT_GAS) {
        reset();
        ShockResult result = Goddard::normal_shock(mach, m_gas.gamma_s());
        store_post_shock_state();
        return result;
    }
    return normal_shock_from_velocity(mach*pre_shock_sound_speed());
}

ShockResult ShockSolver::normal_shock_from_velocity(double velocity) {
    reset();
    ShockResult result;
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS:
            result = Goddard::normal_shock(velocity/m_gas.speed_of_sound(), m_gas.gamma_s());
            break;
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM:
            check_no_condensed_phases(m_gas);
            result = incident_shock(m_gas, velocity, 0.0, m_chemistry, m_options);
            check_no_condensed_phases(m_gas);
            break;
        case GasChemistry::KINETIC:
            break; // rejected by the constructor
    }
    store_post_shock_state();
    return result;
}

ReflectedShockResult ShockSolver::reflected_shock(double mach) {
    if (m_chemistry == GasChemistry::PERFECT_GAS) {
        reset();
        ReflectedShockResult result = Goddard::reflected_shock(mach, m_gas.gamma_s());
        store_post_shock_state();
        return result;
    }
    return reflected_shock_from_velocity(mach*pre_shock_sound_speed());
}

ReflectedShockResult ShockSolver::reflected_shock(
    double mach, GasChemistry incident_chemistry, GasChemistry reflected_chemistry)
{
    if (m_chemistry == GasChemistry::PERFECT_GAS) {
        throw std::invalid_argument(
            "ShockSolver: per-shock chemistry requires a FROZEN or EQUILIBRIUM gas, not PERFECT_GAS.");
    }
    return reflected_shock_from_velocity(mach*pre_shock_sound_speed(), incident_chemistry, reflected_chemistry);
}

ReflectedShockResult ShockSolver::reflected_shock_from_velocity(double velocity) {
    if (m_chemistry == GasChemistry::PERFECT_GAS) {
        reset();
        ReflectedShockResult result = Goddard::reflected_shock(velocity/m_gas.speed_of_sound(), m_gas.gamma_s());
        store_post_shock_state();
        return result;
    }
    return reflected_shock_from_velocity(velocity, m_chemistry, m_chemistry);
}

ReflectedShockResult ShockSolver::reflected_shock_from_velocity(
    double velocity, GasChemistry incident_chemistry, GasChemistry reflected_chemistry)
{
    if (m_chemistry == GasChemistry::PERFECT_GAS) {
        throw std::invalid_argument(
            "ShockSolver: per-shock chemistry requires a FROZEN or EQUILIBRIUM gas, not PERFECT_GAS.");
    }
    check_shock_chemistry(incident_chemistry);
    check_shock_chemistry(reflected_chemistry);
    reset();
    check_no_condensed_phases(m_gas);

    // Initial guess for the reflected shock: the perfect-gas relation with the frozen state-1 gamma
    const double mach1 = velocity/frozen_sound_speed(m_gas);
    const ReflectedShockResult guess = Goddard::reflected_shock(mach1, m_gas.gamma());

    ReflectedShockResult result;
    result.incident = incident_shock(m_gas, velocity, 0.0, incident_chemistry, m_options);
    if (!result.incident.valid) {
        store_post_shock_state();
        return invalid_reflected_shock();
    }
    check_no_condensed_phases(m_gas);

    // Lab-frame velocity of gas 2, which the reflected shock brings to rest
    const double particle_velocity = velocity*(1.0 - 1.0/result.incident.density_ratio);
    if (particle_velocity > 0.0) {
        result.reflected = reflected_shock_from_state2(m_gas, particle_velocity, reflected_chemistry,
            guess.reflected.static_pressure_ratio, guess.reflected.static_temperature_ratio, m_options);
        check_no_condensed_phases(m_gas);
    }
    else {
        // A sonic incident wave leaves gas 2 at rest: the reflection is a sound wave.
        result.reflected = Goddard::normal_shock(1.0, m_gas.gamma());
    }
    result.valid = result.incident.valid && result.reflected.valid;
    store_post_shock_state();
    return result;
}

ObliqueShockResult ShockSolver::oblique_shock_from_wave_angle(double mach, double wave_angle) {
    reset();
    ObliqueShockResult result;
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS:
            result = Goddard::oblique_shock_from_wave_angle(mach, wave_angle, m_gas.gamma_s());
            break;
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM:
            check_no_condensed_phases(m_gas);
            result = oblique_shock_real(m_gas, mach*frozen_sound_speed(m_gas), wave_angle, m_chemistry, m_options);
            check_no_condensed_phases(m_gas);
            break;
        case GasChemistry::KINETIC:
            break; // rejected by the constructor
    }
    store_post_shock_state();
    return result;
}

ObliqueShockResult ShockSolver::oblique_shock_from_deflection(double mach, double deflection, bool weak) {
    reset();
    ObliqueShockResult result;
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS:
            result = Goddard::oblique_shock_from_deflection(mach, deflection, m_gas.gamma_s(), weak);
            break;
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM:
            check_no_condensed_phases(m_gas);
            result = oblique_shock_from_deflection_real(
                m_gas, mach*frozen_sound_speed(m_gas), deflection, weak, m_chemistry, m_options);
            check_no_condensed_phases(m_gas);
            break;
        case GasChemistry::KINETIC:
            break; // rejected by the constructor
    }
    store_post_shock_state();
    return result;
}

double ShockSolver::max_deflection(double mach) {
    reset();
    double deflection = NaN;
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS:
            deflection = oblique_shock_max_deflection(mach, m_gas.gamma_s());
            break;
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            if (mach < 1.0) {
                break;
            }
            check_no_condensed_phases(m_gas);
            const double u1 = mach*frozen_sound_speed(m_gas);
            const std::vector<double> state1 = m_gas.save_state();
            const MaxDeflection max = max_deflection_real(
                m_gas, state1, u1, std::asin(1.0/mach), m_chemistry, m_options);
            deflection_at(m_gas, state1, u1, max.wave_angle, m_chemistry, m_options);
            check_no_condensed_phases(m_gas);
            deflection = max.deflection;
            break;
        }
        case GasChemistry::KINETIC:
            break; // rejected by the constructor
    }
    store_post_shock_state();
    return deflection;
}

const Gas& ShockSolver::pre_shock_state() const {
    reset();
    return m_gas;
}

const Gas& ShockSolver::post_shock_state() const {
    m_gas.restore_state(m_post_shock_state);
    m_gas.chemistry = m_chemistry;
    return m_gas;
}

} //namespace Goddard
