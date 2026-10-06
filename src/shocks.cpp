#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>
#include "goddard/error.hpp"
#include "goddard/shocks.hpp"
#include "goddard/shock_jump.hpp"
#include "goddard/gas_dynamics.hpp"

namespace Goddard {

// ===== utilities =====
namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

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

void ShockSolver::apply_perfect_gas_jump(bool valid, double temperature_ratio, double pressure_ratio) {
    if (!valid) {
        return;
    }
    m_gas.set_state_TP(m_gas.temperature()*temperature_ratio, m_gas.pressure()*pressure_ratio);
}

ShockResult ShockSolver::normal_shock(double mach) {
    if (m_chemistry == GasChemistry::PERFECT_GAS) {
        reset();
        ShockResult result = Goddard::normal_shock(mach, m_gas.gamma_s());
        apply_perfect_gas_jump(result.valid, result.static_temperature_ratio,
                               result.static_pressure_ratio);
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
            apply_perfect_gas_jump(result.valid, result.static_temperature_ratio,
                                   result.static_pressure_ratio);
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
        // State 5 from state 1: T5/T1 = (T2/T1)(T5/T2), and likewise for P.
        apply_perfect_gas_jump(result.valid,
            result.incident.static_temperature_ratio*result.reflected.static_temperature_ratio,
            result.incident.static_pressure_ratio*result.reflected.static_pressure_ratio);
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
        apply_perfect_gas_jump(result.valid,
            result.incident.static_temperature_ratio*result.reflected.static_temperature_ratio,
            result.incident.static_pressure_ratio*result.reflected.static_pressure_ratio);
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
            apply_perfect_gas_jump(result.valid, result.shock.static_temperature_ratio,
                                   result.shock.static_pressure_ratio);
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
            apply_perfect_gas_jump(result.valid, result.shock.static_temperature_ratio,
                                   result.shock.static_pressure_ratio);
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
        case GasChemistry::PERFECT_GAS: {
            const double gamma = m_gas.gamma_s();
            deflection = oblique_shock_max_deflection(mach, gamma);
            if (!std::isnan(deflection)) {
                const ObliqueShockResult max_shock = Goddard::oblique_shock_from_wave_angle(
                    mach, oblique_shock_max_deflection_wave_angle(mach, gamma), gamma);
                apply_perfect_gas_jump(max_shock.valid, max_shock.shock.static_temperature_ratio,
                                       max_shock.shock.static_pressure_ratio);
            }
            break;
        }
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
