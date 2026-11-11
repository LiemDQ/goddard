#include <algorithm>
#include <cmath>
#include <format>
#include <stdexcept>
#include <vector>
#include "cantera/core.h"
#include "goddard/detonations.hpp"
#include "goddard/error.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/shock_jump.hpp"
#include "goddard/solve_guard.hpp"

namespace Goddard {

// ===== utilities =====
namespace {

/** Drive factors within this of 1 give the CJ detonation, as sonic shocks are treated in shock_jump.cpp. */
constexpr double CJ_DRIVE_FACTOR_TOLERANCE = 1e-12;

/** A failed detonation: every field NaN, and the von Neumann shock not computed. */
DetonationResult failed_detonation(const SolveStatus& status) {
    DetonationResult result;
    result.status = status;
    result.von_neumann = failed_shock(SolveStatus::not_computed(status));
    return result;
}

/** A reflected detonation whose incident detonation failed: the reflected shock was not computed. */
ReflectedDetonationResult failed_reflected_detonation(const SolveStatus& status) {
    ReflectedDetonationResult result;
    result.status = status;
    result.incident = failed_detonation(status);
    result.reflected = failed_shock(SolveStatus::not_computed(status));
    return result;
}

void require_gamma_above_one(double gamma, const char* function) {
    if (!(gamma > 1.0)) {
        throw std::invalid_argument(std::format(
            "{}: the ratio of specific heats must be greater than 1, got {:.6g}.", function, gamma));
    }
}

SolveStatus endothermic_status(double heat_release) {
    return SolveStatus::failure(FailureReason::NO_SOLUTION, std::format(
        "no detonation exists: the heat release q/(R T1) = {:.6g} is negative.", heat_release));
}

SolveStatus below_chapman_jouguet_status(double drive_factor) {
    return SolveStatus::failure(FailureReason::NO_SOLUTION, std::format(
        "no steady detonation exists below the Chapman-Jouguet speed: the drive factor "
        "u1/u_CJ = {:.6g} is below 1.", drive_factor));
}

/**
 * Compression s = 1 - rho1/rho2 at a root of the Rayleigh line and the Hugoniot of the one-gamma
 * model: s = [(M^2 - 1) +/- sqrt((M^2 - 1)^2 - b M^2)]/(a M^2), + for the overdriven root.
 *
 * With a = gamma + 1 and b = 4H it is exact for a calorically perfect gas. With a and b fitted to
 * a real-gas CJ solution it is the initial guess for a real-gas detonation.
 */
double compression(double mach, double a, double b, DetonationBranch branch) {
    const double M2 = mach*mach;
    const double root = std::sqrt(std::max(0.0, (M2 - 1.0)*(M2 - 1.0) - b*M2));
    double signed_root = 0.0;
    switch (branch) {
        case DetonationBranch::OVERDRIVEN:
            signed_root = root;
            break;
        case DetonationBranch::UNDERDRIVEN:
            signed_root = -root;
            break;
    }
    return ((M2 - 1.0) + signed_root)/(a*M2);
}

/** Jump across a one-gamma detonation of Mach number `mach` and compression `s` = 1 - rho1/rho2. */
DetonationResult perfect_gas_detonation(double mach, double s, double drive_factor, double gamma) {
    DetonationResult result;
    result.velocity = NaN;
    result.drive_factor = drive_factor;
    result.mach_in = mach;
    // Rayleigh line: P2/P1 - 1 = gamma M^2 s
    result.static_pressure_ratio = 1.0 + gamma*mach*mach*s;
    result.density_ratio = 1.0/(1.0 - s);
    result.static_temperature_ratio = (1.0 - s)*result.static_pressure_ratio;
    // u2/u1 = 1 - s and a2^2/a1^2 = T2/T1
    result.mach_out = mach*(1.0 - s)/std::sqrt(result.static_temperature_ratio);
    result.molecular_weight_ratio = 1.0;
    result.total_pressure_ratio =
        perfect_gas_stagnation_pressure(result.static_pressure_ratio, result.mach_out, gamma)
        /perfect_gas_stagnation_pressure(1.0, mach, gamma);
    result.von_neumann = Goddard::normal_shock(mach, gamma);
    return result;
}

/**
 * Mach number, relative to the gas ahead of it, of the shock that brings gas moving at
 * `particle_mach` = u_p/a to rest, in a calorically perfect gas. It inverts
 * u_p = 2a (M_R - 1/M_R)/(gamma + 1).
 */
double reflected_shock_mach(double particle_mach, double gamma) {
    const double c = (gamma + 1.0)*particle_mach/4.0;
    return c + std::sqrt(c*c + 1.0);
}

/** Initial guess for the Newton iteration on (ln P2/P1, ln T2/T1). */
struct JumpGuess {
    double pressure_ratio;
    double temperature_ratio;
};

/**
 * Initial estimate of the CJ jump, RP-1311 section 8.3. `gas` holds the unburned gas on entry,
 * with EQUILIBRIUM chemistry, and an arbitrary state on exit.
 */
JumpGuess chapman_jouguet_initial_estimate(Gas& gas) {
    const double R = Cantera::GasConstant;
    const double T1 = gas.temperature();
    const double P1 = gas.pressure();
    const double h1 = gas.enthalpy_mass();
    const double M1 = gas.molecular_weight();

    // (8.13): flame temperature at an assumed pressure ratio
    const double pressure_ratio_0 = 15.0;
    gas.equilibrate_HP(h1 + 0.75*R*T1/M1*pressure_ratio_0, pressure_ratio_0*P1);
    const double temperature_ratio_0 = gas.temperature()/T1;
    const ExpansionProperties props = gas.expansion_properties();
    const double gamma_s = props.gamma_s;
    const double cp = props.spec_heat_p;
    const double M2 = gas.molecular_weight();

    // (8.14)-(8.17), repeated three times with the properties at the flame temperature
    JumpGuess guess{pressure_ratio_0, temperature_ratio_0};
    for (int k = 0; k < 3; k++) {
        const double alpha = M2/(guess.temperature_ratio*M1);
        const double discriminant = 1.0 - 4.0*gamma_s*alpha/((1.0 + gamma_s)*(1.0 + gamma_s));
        guess.pressure_ratio = (1.0 + gamma_s)/(2.0*gamma_s*alpha)*(1.0 + std::sqrt(std::max(0.0, discriminant)));
        const double r = alpha*guess.pressure_ratio;  // rho2/rho1
        guess.temperature_ratio = temperature_ratio_0 - 0.75*R/(M1*cp)*pressure_ratio_0
            + R*gamma_s/(2.0*M1*cp)*(r*r - 1.0)/r*guess.pressure_ratio;
    }
    if (!(guess.temperature_ratio > 1.0)) {
        guess.temperature_ratio = temperature_ratio_0;
    }
    return guess;
}

/**
 * Initial guess for an overdriven or under-driven detonation: the one-gamma model fitted to the
 * CJ solution, which is exact at a drive factor of 1 and splits into the two branches as
 * sqrt(drive_factor - 1). See docs/theory/detonations.md.
 *
 * @param momentum_coeff rho1 u1^2/P1 [-]
 */
JumpGuess driven_initial_guess(const DetonationResult& chapman_jouguet, double drive_factor,
                               double momentum_coeff, DetonationBranch branch)
{
    const double mach_cj2 = chapman_jouguet.mach_in*chapman_jouguet.mach_in;
    const double compression_cj = 1.0 - 1.0/chapman_jouguet.density_ratio;
    // compression() at M_CJ has a double root at compression_cj
    const double a = (mach_cj2 - 1.0)/(compression_cj*mach_cj2);
    const double b = (mach_cj2 - 1.0)*(mach_cj2 - 1.0)/mach_cj2;
    const double s = std::min(compression(drive_factor*chapman_jouguet.mach_in, a, b, branch), 0.99);

    JumpGuess guess;
    guess.pressure_ratio = 1.0 + momentum_coeff*s;
    // ideal gas: T2/T1 = (P2/P1)(rho1/rho2)(M2/M1), with M2 from the CJ products
    guess.temperature_ratio = guess.pressure_ratio*(1.0 - s)*chapman_jouguet.molecular_weight_ratio;
    return guess;
}

/**
 * Assemble a real-gas detonation result from the jump to the products.
 *
 * `gas` holds the products on entry and on exit, with EQUILIBRIUM chemistry. `von_neumann_state`
 * receives the state behind the frozen shock at the same wave speed.
 */
DetonationResult real_gas_detonation(Gas& gas, const std::vector<double>& unburned_state,
                                     double velocity, double drive_factor, const ShockJump& jump,
                                     const SolverOptions& opts, std::vector<double>& von_neumann_state)
{
    const std::vector<double> products_state = gas.save_state();
    const double products_molecular_weight = gas.molecular_weight();
    const double u2 = velocity/jump.density_ratio;
    const double mach_out = u2/gas.speed_of_sound();
    const double P02 = gas.stagnation_pressure(u2);

    gas.restore_state(unburned_state);
    gas.chemistry = GasChemistry::FROZEN;
    const double unburned_molecular_weight = gas.molecular_weight();
    const double a1 = gas.speed_of_sound();
    const double P01 = gas.stagnation_pressure(velocity);
    const ShockResult von_neumann = incident_shock(gas, velocity, 0.0, GasChemistry::FROZEN, opts);
    gas.copy_state(von_neumann_state);

    gas.restore_state(products_state);
    gas.chemistry = GasChemistry::EQUILIBRIUM;

    DetonationResult result;
    result.velocity = velocity;
    result.drive_factor = drive_factor;
    result.mach_in = velocity/a1;
    result.mach_out = mach_out;
    result.static_pressure_ratio = jump.pressure_ratio;
    result.static_temperature_ratio = jump.temperature_ratio;
    result.density_ratio = jump.density_ratio;
    result.molecular_weight_ratio = products_molecular_weight/unburned_molecular_weight;
    result.total_pressure_ratio = P02/P01;
    result.von_neumann = von_neumann;
    return result;
}

/**
 * Reject an unburned mixture carrying condensed species: the detonation relations are gas-phase
 * only.
 *
 * @throws std::invalid_argument if any condensed species is present in a nonzero amount.
 */
void check_no_condensed_phases(const Gas& gas) {
    if (gas.has_condensed_phases()) {
        throw std::invalid_argument(
            "DetonationSolver: detonations are gas-phase only; condensed species are not supported.");
    }
}

/**
 * Fail the solve if condensed species formed in the products. Unlike the check on the unburned
 * gas this depends on the solution, so it is a solve failure and not a bad argument.
 *
 * @throws SolveError with reason NO_SOLUTION if any condensed species is present.
 */
void check_no_condensed_products(const Gas& gas) {
    if (gas.has_condensed_phases()) {
        throw SolveError(FailureReason::NO_SOLUTION,
            "condensed species form behind the wave; detonations are gas-phase only.");
    }
}

} // namespace

// ===== Perfect gas free functions =====

DetonationResult chapman_jouguet_detonation(double gamma, double heat_release) {
    require_gamma_above_one(gamma, "chapman_jouguet_detonation");
    if (!(heat_release >= 0.0)) {
        return failed_detonation(endothermic_status(heat_release));
    }
    const double H = (gamma*gamma - 1.0)*heat_release/(2.0*gamma);
    const double mach = std::sqrt(H) + std::sqrt(H + 1.0);
    // double root of compression() at the CJ Mach number
    const double s = (mach*mach - 1.0)/((gamma + 1.0)*mach*mach);
    return perfect_gas_detonation(mach, s, 1.0, gamma);
}

DetonationResult detonation(double drive_factor, double gamma, double heat_release, DetonationBranch branch) {
    require_gamma_above_one(gamma, "detonation");
    if (!(heat_release >= 0.0)) {
        return failed_detonation(endothermic_status(heat_release));
    }
    if (!(drive_factor >= 1.0 - CJ_DRIVE_FACTOR_TOLERANCE)) {
        return failed_detonation(below_chapman_jouguet_status(drive_factor));
    }
    if (drive_factor <= 1.0 + CJ_DRIVE_FACTOR_TOLERANCE) {
        return chapman_jouguet_detonation(gamma, heat_release);
    }
    const double H = (gamma*gamma - 1.0)*heat_release/(2.0*gamma);
    const double mach = drive_factor*(std::sqrt(H) + std::sqrt(H + 1.0));
    return perfect_gas_detonation(mach, compression(mach, gamma + 1.0, 4.0*H, branch), drive_factor, gamma);
}

ReflectedDetonationResult reflected_detonation(double drive_factor, double gamma, double heat_release,
                                               DetonationBranch branch)
{
    ReflectedDetonationResult result;
    result.incident = detonation(drive_factor, gamma, heat_release, branch);
    if (!result.incident.status.ok()) {
        return failed_reflected_detonation(result.incident.status);
    }
    // Lab-frame velocity of the products: u_p = u1 - u2, so u_p/a2 = M2 (rho2/rho1 - 1)
    const double particle_mach = result.incident.mach_out*(result.incident.density_ratio - 1.0);
    result.reflected = Goddard::normal_shock(reflected_shock_mach(particle_mach, gamma), gamma);
    result.status = result.reflected.status;
    return result;
}

// ===== DetonationSolver =====

DetonationSolver::DetonationSolver(const Gas& gas, SolverOptions options)
    : m_gas(gas.clone()), m_chemistry(m_gas.chemistry), m_options(options),
      m_post_detonation_status(SolveStatus::failure(FailureReason::NOT_COMPUTED,
          "DetonationSolver: no detonation has been solved yet.")),
      m_von_neumann_status(m_post_detonation_status)
{
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS:
            throw std::invalid_argument(
                "DetonationSolver: PERFECT_GAS chemistry has no heat release; use "
                "chapman_jouguet_detonation() and detonation() instead.");
        case GasChemistry::KINETIC:
            throw std::invalid_argument(
                "DetonationSolver: KINETIC chemistry is not supported; use FROZEN or EQUILIBRIUM.");
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM:
            break;
    }
    check_no_condensed_phases(m_gas);
    m_gas.copy_state(m_pre_detonation_state);
}

void DetonationSolver::reset() const {
    m_gas.restore_state(m_pre_detonation_state);
    m_gas.chemistry = m_chemistry;
}

void DetonationSolver::store_failure(const SolveStatus& status) {
    m_post_detonation_status = status;
    m_von_neumann_status = status;
}

const DetonationResult& DetonationSolver::cached_chapman_jouguet() {
    if (!m_chapman_jouguet_solved) {
        const SolveStatus status = capture_failure("Chapman-Jouguet detonation", [&] {
            m_chapman_jouguet = solve_chapman_jouguet();
        });
        if (!status.ok()) {
            m_chapman_jouguet = failed_detonation(status);
        }
        m_chapman_jouguet_solved = true;
        reset();
    }
    return m_chapman_jouguet;
}

DetonationResult DetonationSolver::solve_chapman_jouguet() {
    reset();
    m_gas.chemistry = GasChemistry::EQUILIBRIUM;
    const double T1 = m_gas.temperature();

    // A mixture that releases no heat has no detonation, and the CJ Jacobian is singular at its
    // trivial root P2 = P1, T2 = T1.
    m_gas.equilibrate_HP(m_gas.enthalpy_mass(), m_gas.pressure());
    const double flame_temperature = m_gas.temperature();
    m_gas.restore_state(m_pre_detonation_state);
    if (!(flame_temperature > T1*(1.0 + 1e-6))) {
        throw SolveError(FailureReason::NO_SOLUTION, std::format(
            "no detonation exists: the mixture releases no heat (constant-pressure combustion "
            "takes it from {:.6g} K to {:.6g} K).", T1, flame_temperature));
    }

    const JumpGuess guess = chapman_jouguet_initial_estimate(m_gas);
    m_gas.restore_state(m_pre_detonation_state);
    const ShockJump jump = solve_shock_jump(m_gas, JumpCondition::CHAPMAN_JOUGUET, 0.0,
        guess.pressure_ratio, guess.temperature_ratio, m_options);
    check_no_condensed_products(m_gas);
    if (!(jump.pressure_ratio > 1.0 && jump.density_ratio > 1.0)) {
        throw ConvergenceError(std::format(
            "the Chapman-Jouguet solve converged to a root that is not a detonation "
            "(P2/P1 = {:.6g}, rho2/rho1 = {:.6g}).", jump.pressure_ratio, jump.density_ratio));
    }

    // u2 = a2, so u1 = (rho2/rho1) a2
    const double velocity = jump.density_ratio*m_gas.speed_of_sound();
    const DetonationResult result = real_gas_detonation(m_gas, m_pre_detonation_state, velocity,
        1.0, jump, m_options, m_chapman_jouguet_von_neumann_state);
    result.von_neumann.status.raise_if_failed();
    m_gas.copy_state(m_chapman_jouguet_state);
    return result;
}

DetonationResult DetonationSolver::chapman_jouguet() {
    const DetonationResult& result = cached_chapman_jouguet();
    if (result.status.ok()) {
        m_post_detonation_state = m_chapman_jouguet_state;
        m_von_neumann_state = m_chapman_jouguet_von_neumann_state;
    }
    m_post_detonation_status = result.status;
    m_von_neumann_status = result.status;
    reset();
    return result;
}

DetonationResult DetonationSolver::detonation(double drive_factor, DetonationBranch branch) {
    const DetonationResult& chapman_jouguet = cached_chapman_jouguet();
    if (!chapman_jouguet.status.ok()) {
        store_failure(chapman_jouguet.status);
        return failed_detonation(chapman_jouguet.status);
    }
    return driven_detonation(drive_factor*chapman_jouguet.velocity, drive_factor, branch);
}

DetonationResult DetonationSolver::detonation_from_velocity(double velocity, DetonationBranch branch) {
    const DetonationResult& chapman_jouguet = cached_chapman_jouguet();
    if (!chapman_jouguet.status.ok()) {
        store_failure(chapman_jouguet.status);
        return failed_detonation(chapman_jouguet.status);
    }
    return driven_detonation(velocity, velocity/chapman_jouguet.velocity, branch);
}

DetonationResult DetonationSolver::driven_detonation(double velocity, double drive_factor,
                                                     DetonationBranch branch)
{
    if (!(drive_factor >= 1.0 - CJ_DRIVE_FACTOR_TOLERANCE)) {
        SolveStatus status = below_chapman_jouguet_status(drive_factor);
        status.stage = "detonation";
        store_failure(status);
        reset();
        return failed_detonation(status);
    }
    if (drive_factor <= 1.0 + CJ_DRIVE_FACTOR_TOLERANCE) {
        return chapman_jouguet();
    }

    DetonationResult result;
    const SolveStatus status = capture_failure("detonation", [&] {
        result = solve_driven(velocity, drive_factor, branch);
    });
    if (!status.ok()) {
        result = failed_detonation(status);
        store_failure(status);
        reset();
    }
    return result;
}

DetonationResult DetonationSolver::solve_driven(double velocity, double drive_factor, DetonationBranch branch) {
    reset();
    const double momentum_coeff = m_gas.density()*velocity*velocity/m_gas.pressure();
    const JumpGuess guess = driven_initial_guess(m_chapman_jouguet, drive_factor, momentum_coeff, branch);
    // The overdriven and under-driven roots solve the incident-shock jump conditions with the
    // products at equilibrium; the initial guess selects the root.
    m_gas.chemistry = GasChemistry::EQUILIBRIUM;
    const ShockJump jump = solve_shock_jump(m_gas, JumpCondition::INCIDENT, velocity,
        guess.pressure_ratio, guess.temperature_ratio, m_options);
    check_no_condensed_products(m_gas);

    std::vector<double> von_neumann_state;
    const DetonationResult result = real_gas_detonation(m_gas, m_pre_detonation_state, velocity,
                                                        drive_factor, jump, m_options, von_neumann_state);
    result.von_neumann.status.raise_if_failed();
    bool on_branch = false;
    switch (branch) {
        case DetonationBranch::OVERDRIVEN:
            on_branch = result.mach_out < 1.0;
            break;
        case DetonationBranch::UNDERDRIVEN:
            on_branch = result.mach_out > 1.0;
            break;
    }
    if (!on_branch || !(jump.pressure_ratio > 1.0 && jump.density_ratio > 1.0)) {
        throw ConvergenceError(std::format(
            "the jump conditions converged off the requested branch (products Mach number {:.6g}, "
            "P2/P1 = {:.6g}); the drive factor {:.6g} may be too close to 1.",
            result.mach_out, jump.pressure_ratio, drive_factor));
    }
    m_gas.copy_state(m_post_detonation_state);
    m_von_neumann_state = von_neumann_state;
    m_post_detonation_status = {};
    m_von_neumann_status = {};
    return result;
}

ReflectedDetonationResult DetonationSolver::reflected_detonation(double drive_factor, DetonationBranch branch) {
    ReflectedDetonationResult result;
    result.incident = detonation(drive_factor, branch);
    if (!result.incident.status.ok()) {
        return failed_reflected_detonation(result.incident.status);
    }

    const SolveStatus status = capture_failure("reflected shock", [&] {
        m_gas.restore_state(m_post_detonation_state);
        m_gas.chemistry = GasChemistry::EQUILIBRIUM;

        // Lab-frame velocity of the products, which the reflected shock brings to rest
        const double particle_velocity = result.incident.velocity*(1.0 - 1.0/result.incident.density_ratio);
        // Initial guess: the perfect-gas reflected shock with the equilibrium gamma of the products
        const double gamma_s = m_gas.gamma_s();
        const ShockResult guess = Goddard::normal_shock(
            reflected_shock_mach(particle_velocity/m_gas.speed_of_sound(), gamma_s), gamma_s);
        result.reflected = reflected_shock_from_state2(m_gas, particle_velocity, GasChemistry::EQUILIBRIUM,
            guess.static_pressure_ratio, guess.static_temperature_ratio, m_options);
        check_no_condensed_products(m_gas);
        m_gas.copy_state(m_post_detonation_state);
    });
    if (!status.ok()) {
        result.reflected = failed_shock(status);
    }
    set_stage_if_unset(result.reflected.status, "reflected shock");
    result.status = result.reflected.status;
    // The von Neumann state stays that of the incident detonation.
    m_post_detonation_status = result.status;
    return result;
}

Gas DetonationSolver::pre_detonation_state() const {
    reset();
    return m_gas.clone();
}

Gas DetonationSolver::post_detonation_state() const {
    m_post_detonation_status.raise_if_failed();
    m_gas.restore_state(m_post_detonation_state);
    m_gas.chemistry = GasChemistry::EQUILIBRIUM;
    return m_gas.clone();
}

Gas DetonationSolver::von_neumann_state() const {
    m_von_neumann_status.raise_if_failed();
    m_gas.restore_state(m_von_neumann_state);
    m_gas.chemistry = GasChemistry::FROZEN;
    return m_gas.clone();
}

} // namespace Goddard
