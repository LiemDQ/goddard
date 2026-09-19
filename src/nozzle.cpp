#include "goddard/nozzle.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/error.hpp"
#include "goddard/utils.hpp"
#include "goddard/gas_dynamics.hpp"

#include "cantera/core.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <exception>
#include <vector>
#include <functional>
#include <format>
#include <string>
namespace Goddard {

namespace {

/**
 * Re-throw a frozen temperature-range failure, naming the station it happened at.
 *
 * A frozen expansion holds the condensed amounts fixed, so it cannot be carried past the
 * temperature range of a condensed species that is present. `Gas` reports that as a `FmtError`;
 * the nozzle adds the station index and the temperature [K] it stopped at, as CEA does.
 *
 * @param error Range error raised by `Gas`.
 * @param station Index of the station being solved [-].
 * @param temperature Temperature the frozen iteration stopped at [K].
 */
[[noreturn]] void rethrow_frozen_range_error(const FmtError& error, int station,
                                             double temperature)
{
    std::string message = error.what();
    throw FmtError("Frozen expansion: {} at station {} (T = {} K)", message, station, temperature);
}

} // namespace

Nozzle::Nozzle(const Gas& gas, NozzleOptions options)
    : inlet_state(gas.save_state()), m_gas(gas), m_frozen_state(gas.state_size()), m_opts(options) {
    determine_and_set_equilibrium_condition();
    if (m_opts.chemistry == GasChemistry::KINETIC) {
        throw std::invalid_argument("GasChemistry::KINETIC is not valid for Nozzle. Use KineticNozzle instead.");
    }
    if (m_opts.chemistry == GasChemistry::PERFECT_GAS) {
        throw std::runtime_error("GasChemistry::PERFECT_GAS nozzle not yet implemented.");
    }
}

Nozzle::Nozzle(const Gas& gas, std::vector<double> state, NozzleOptions options)
    : inlet_state(std::move(state)), m_gas(gas),  m_frozen_state(gas.state_size()), m_opts(options) {
    determine_and_set_equilibrium_condition();
    if (m_opts.chemistry == GasChemistry::KINETIC) {
        throw std::invalid_argument("GasChemistry::KINETIC is not valid for Nozzle. Use KineticNozzle instead.");
    }
    if (m_opts.chemistry == GasChemistry::PERFECT_GAS) {
        throw std::runtime_error("GasChemistry::PERFECT_GAS nozzle not yet implemented.");
    }
}

NozzleResults Nozzle::solve(ExpansionType expansion_type, double ratio) {
    const ThroatCondition throat_condition = solve_throat_conditions();
    return {throat_condition, solve_stations(throat_condition, expansion_type, {ratio})};
}

NozzleResults Nozzle::solve(ExpansionType expansion_type, const std::vector<double>& ratios) {
    const ThroatCondition throat_condition = solve_throat_conditions();
    return {throat_condition, solve_stations(throat_condition, expansion_type, ratios)};
}

NozzleResults Nozzle::solve(const NozzleProfile& profile, int num_stations) {
    const ThroatCondition throat_condition = solve_throat_conditions();

    // For a diverging-only profile (starting at throat), x_min is the throat.
    double A_throat = profile.area_at(profile.x_min());

    double x_start = profile.x_min();
    double x_end = profile.x_max();
    // Slightly inset from x_max to avoid out-of-bounds interpolation at the boundary
    double x_range = x_end - x_start;
    double x_last = x_end - 1e-10 * x_range;

    std::vector<NozzleStation> results;
    for (int i = 1; i <= num_stations; i++) {
        m_current_station = i;
        double x = x_start + (x_last - x_start) * static_cast<double>(i) / num_stations;
        double A = profile.area_at(x);
        double area_ratio = A / A_throat;
        results.push_back(solve_supersonic_area_expansion(throat_condition, area_ratio));
    }
    return {throat_condition, results};
}

std::vector<NozzleStation> Nozzle::solve_stations(const ThroatCondition& throat_condition,
    ExpansionType expansion_type, const std::vector<double>& ratios) {

    std::vector<NozzleStation> results;
    results.reserve(ratios.size());
    for (double ratio: ratios) {
        switch (expansion_type) {
            case ExpansionType::SUPERSONIC_AREA_RATIO: {
                results.push_back(solve_supersonic_area_expansion(throat_condition, ratio));
                break;
            }
            case ExpansionType::SUBSONIC_AREA_RATIO: {
                results.push_back(solve_subsonic_area_expansion(throat_condition, ratio));
                break;
            }
            case ExpansionType::PRESSURE_RATIO: {
                results.push_back(solve_pressure_ratio(throat_condition, ratio));
                break;
            }
            default:
                throw NotImplementedError("Expansion type is not implemented.");
        }
        m_current_station++;
    }
    return results;
}

FiniteAreaChamber Nozzle::solve_finite_area_chamber(const std::vector<double>& injector_state,
    CombustorType type, double value, double reltol) {

    bool contraction_mode = true;
    switch (type) {
        case CombustorType::FINITE_CONTRACTION_RATIO: {
            if (!(value > 1.0)) {
                throw std::invalid_argument(std::format(
                    "Finite-area combustor: contraction ratio must be greater than 1. Actual value: {}",
                    value));
            }
            contraction_mode = true;
            break;
        }
        case CombustorType::FINITE_MASS_FLUX: {
            if (!(value > 0.0)) {
                throw std::invalid_argument(std::format(
                    "Finite-area combustor: mass flux must be positive. Actual value: {} kg/(m^2 s)",
                    value));
            }
            contraction_mode = false;
            break;
        }
        case CombustorType::INFINITE_AREA:
        case CombustorType::NONE:
        default:
            throw std::invalid_argument(
                "Nozzle::solve_finite_area_chamber requires FINITE_CONTRACTION_RATIO or FINITE_MASS_FLUX.");
    }
    if (m_opts.chemistry == GasChemistry::FROZEN && m_opts.frozen_NFZ == 0) {
        throw NotImplementedError(
            "Finite-area combustor with flow frozen at the combustion end (frozen_NFZ = 0) is not implemented.");
    }

    // Injector face: equilibrium at (h_inj, P_inj).
    m_gas.restore_state(injector_state);
    m_current_station = 0;
    determine_and_set_equilibrium_condition();
    const double h_injector = m_gas.enthalpy_mass();
    const double P_injector = m_gas.pressure();
    const double T_injector = m_gas.temperature();
    const double gamma_injector = m_gas.gamma_s();
    const double R_specific = Cantera::GasConstant / m_gas.molecular_weight();

    // Perfect-gas initial guess for P_inf (see instructions/finite_area_combustor.md, section 3).
    double mach_guess = 0.0;
    if (contraction_mode) {
        mach_guess = mach_from_area_ratio(value, gamma_injector, false);
    } else {
        const double g = gamma_injector;
        // Throat mass flux per unit stagnation pressure [kg/(m^2 s Pa)].
        const double f = std::sqrt(g / (R_specific * T_injector))
            * std::pow(2.0 / (g + 1.0), (g + 1.0) / (2.0 * (g - 1.0)));
        const double G_max = P_injector * f / finite_area_pressure_loss(1.0, g);
        // The perfect-gas limit is only an estimate of the real-gas one, so the chamber is rejected
        // here only when clearly choked; closer to the limit the real-gas iteration decides.
        if (value > 1.02 * G_max) {
            throw std::invalid_argument(std::format(
                "Finite-area combustor: mass flux {} kg/(m^2 s) thermally chokes the chamber. "
                "Maximum mass flux (perfect gas estimate): {} kg/(m^2 s)", value, G_max));
        }
        // eps(M) - P_inj f / (phi(M) G) is positive as M -> 0 and negative at M = 1.
        double lower = 0.0;
        double upper = 1.0;
        for (int i = 0; i < 100; i++) {
            const double mid = 0.5 * (lower + upper);
            const double excess = area_mach_relation(mid, g)
                - P_injector * f / (finite_area_pressure_loss(mid, g) * value);
            if (excess > 0.0) {
                lower = mid;
            } else {
                upper = mid;
            }
        }
        mach_guess = value < G_max ? 0.5 * (lower + upper) : 0.95;
    }
    double P_stagnation = P_injector / finite_area_pressure_loss(mach_guess, gamma_injector);

    // The station tolerance must sit below the balance tolerance, or the residual is noisy.
    const double station_abstol = std::clamp(0.1 * reltol, 1e-7, 4.5e-5);

    const int max_iters = 50;
    const int max_backtracks = 5;
    double ln_P = std::log(P_stagnation);
    double ln_P_previous = 0.0;
    double residual_previous = 0.0;
    bool have_previous = false;
    double residual = 1.0;
    double contraction_ratio = value;

    ThroatCondition throat;
    NozzleStation combustion_end;
    std::vector<double> stagnation_state;

    for (int iter = 1; iter <= max_iters; iter++) {
        // Stagnation state and throat; in mass-flux mode also the contraction ratio, backtracking
        // toward the previous iterate if the chamber would be choked.
        for (int backtrack = 0; ; backtrack++) {
            P_stagnation = std::exp(ln_P);
            m_gas.restore_state(injector_state);
            m_gas.equilibrate_HP(h_injector, P_stagnation);
            set_inlet_state(m_gas.save_state());
            throat = solve_throat_conditions();
            if (contraction_mode) {
                break;
            }
            m_gas.restore_state(throat.state);
            contraction_ratio = m_gas.density() * throat.speed_of_sound / value;
            if (contraction_ratio > 1.0001) {
                break;
            }
            if (!have_previous || backtrack >= max_backtracks) {
                throw std::invalid_argument(std::format(
                    "Finite-area combustor: mass flux {} kg/(m^2 s) thermally chokes the chamber "
                    "(A_c/A_t = {} at P_inf = {} Pa).", value, contraction_ratio, P_stagnation));
            }
            ln_P = 0.5 * (ln_P + ln_P_previous);
        }

        // Combustion end is part of the chamber and therefore in equilibrium.
        m_current_station = 0;
        combustion_end = solve_subsonic_area_expansion(throat, contraction_ratio, station_abstol);
        const double velocity = m_gas.isenthalpic_velocity();
        const double P_injector_calc = m_gas.pressure() + m_gas.density() * velocity * velocity;

        if (std::abs(1.0 - P_injector_calc / P_injector) < reltol) {
            double mass_flux = value;
            if (contraction_mode) {
                m_gas.restore_state(throat.state);
                mass_flux = m_gas.density() * throat.speed_of_sound / contraction_ratio;
            }
            // The combustion-end velocity comes from h_inj - h_c. The station solve leaves an error
            // of about station_abstol * P_c/rho_c in h_c, and floating point one of about
            // eps * |h_inj|. When that noise is not small against u_c^2 (a very large contraction
            // ratio, far beyond any real chamber), u_c and the Mach number are not resolved.
            // Pressures and the momentum balance are unaffected.
            const double enthalpy_noise = station_abstol * P_injector_calc / m_gas.density()
                + std::numeric_limits<double>::epsilon() * std::abs(h_injector);
            if (enthalpy_noise > 0.01 * velocity * velocity) {
                Cantera::warn_user("Nozzle::solve_finite_area_chamber",
                    "Combustion-end velocity ({} m/s) and Mach number are below the solver's "
                    "resolution at A_c/A_t = {}; treat them as zero.", velocity, contraction_ratio);
            }
            stagnation_state = inlet_state;
            m_gas.restore_state(inlet_state);
            m_current_station = 1;
            return {stagnation_state, combustion_end, throat, P_injector, P_stagnation,
                contraction_ratio, mass_flux, iter};
        }

        // Proportional step first, then secant on ln P_inf versus ln(P_inj,calc / P_inj).
        residual = std::log(P_injector_calc / P_injector);
        double slope = 1.0;
        if (have_previous && residual != residual_previous) {
            slope = (residual - residual_previous) / (ln_P - ln_P_previous);
        }
        if (!(slope > 0.1)) {
            // P_inj,calc grows roughly in proportion to P_inf; fall back to the proportional step.
            slope = 1.0;
        }
        ln_P_previous = ln_P;
        residual_previous = residual;
        have_previous = true;
        ln_P -= residual / slope;
    }
    throw ConvergenceError("Finite-area combustor momentum balance failed to converge.",
        max_iters, reltol, std::abs(residual));
}

void Nozzle::reset_state(){
    m_gas.restore_state(inlet_state);
    m_current_station = 0;
}

ThroatCondition Nozzle::solve_throat_conditions(double abstol) {
    m_gas.restore_state(inlet_state);
    m_gas.set_current_state_as_reference();
    m_current_station = 0;
    
    double P_inlet = m_gas.pressure();
    double S_inlet = m_gas.entropy_mass();
    double H_inlet = m_gas.enthalpy_mass();
    const bool is_equilibrium = determine_and_set_equilibrium_condition();
    
    std::vector<double> X_inlet = m_gas.mole_fractions(); 

    double gamma_s = m_gas.gamma_s();
    double P_throat = P_inlet / std::pow((gamma_s+1)/2,gamma_s/(gamma_s-1));

    const int max_iters = 5;
    int iter = 0;
    double Mach = 1.0; //throat mach number is 1 by definition
    double residual = 1.0;

    while (residual > abstol){
        if (iter >= max_iters){
            throw ConvergenceError("Throat conditions failed to converge.", iter, abstol, residual);
        }
        P_throat = P_throat * (1 + gamma_s * Mach* Mach)/(1+ gamma_s);

        //if equilibrium conditions are selected, the composition must reach chemical
        //equilibrium in the throat.
        if (is_equilibrium) {
            m_gas.equilibrate_SP(S_inlet, P_throat);
        } else {
            try {
                m_gas.set_state_SP(S_inlet, P_throat);
            } catch (const FmtError& error) {
                rethrow_frozen_range_error(error, m_current_station, m_gas.temperature());
            }
        }
        gamma_s = m_gas.gamma_s();

        double velocity = m_gas.isenthalpic_velocity();
        double sonic_velocity = m_gas.speed_of_sound();

        Mach = velocity/sonic_velocity;
        residual = std::abs(1.0 - 1.0/(Mach * Mach));

        iter++;
    }

    ExpansionProperties final_props = m_gas.expansion_properties();
    // regardless of whether the throat is equilibrium or frozen
    // the first frozen state will be the throat. 
    m_frozen_state = m_gas.save_state(); 

    m_current_station = 1;

    return {true,
        m_gas.speed_of_sound(),
        H_inlet,
        P_inlet,
        S_inlet,
        gamma_s,
        final_props.dlogV_dlogP_T,
        final_props.dlogV_dlogT_P,
        m_gas.save_state(),
        final_props.pinned_transition};
}

double Nozzle::get_gamma_s() {
    return m_gas.gamma_s();
}

NozzleStation Nozzle::solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol) {
    m_gas.restore_state(throat_condition.state);

    if (!(expansion_ratio > 1.0001)) {
        throw_invalid_expansion_ratio(expansion_ratio, 1.0001);
    }
    // Perfect-gas initial guess with the throat isentropic exponent. Unlike a correlation it stays
    // accurate as A/A_t grows and the station approaches stagnation.
    double gamma_s = m_gas.gamma_s();
    const double mach_guess = mach_from_area_ratio(expansion_ratio, gamma_s, false);
    double ln_pressure_ratio =
        gamma_s / (gamma_s - 1.0) * std::log(stagnation_factor(mach_guess, gamma_s));

    // Safeguarded Newton iteration on x = ln(P_inlet/P). On the subsonic branch A/A_t decreases
    // monotonically from infinity at x = 0 to 1 at the throat, so [0, x_throat] brackets the root.
    // Pure Newton diverges near the throat, where dlnA/dlnP vanishes, and near stagnation, where
    // x is tiny; steps leaving the bracket fall back to bisection.
    double velocity = m_gas.isenthalpic_velocity();
    const double A_mdot_throat = m_gas.area_per_mdot(velocity);
    std::vector<double> composition = m_gas.mole_fractions();
    double temperature = m_gas.temperature();
    const bool is_equilibrium = determine_and_set_equilibrium_condition();
    const double ln_throat_ratio = std::log(throat_condition.P_inlet / m_gas.pressure());
    const double ln_expansion_ratio = std::log(expansion_ratio);

    double lower = 0.0;
    double upper = ln_throat_ratio;
    if (!(ln_pressure_ratio > lower && ln_pressure_ratio < upper)) {
        ln_pressure_ratio = 0.5 * (lower + upper);
    }

    const int max_iters = 60;
    double step = 1.0;
    int iters = 0;
    while (std::abs(step) > abstol) {
        iters++;
        if (iters > max_iters) {
            throw ConvergenceError("Maximum number of iterations exceeded for subsonic area expansion.",
                iters, abstol, std::abs(step));
        }
        const double pressure_ratio = std::exp(ln_pressure_ratio);
        if (is_equilibrium) {
            m_gas.equilibrate_SP(throat_condition.S_inlet, throat_condition.P_inlet / pressure_ratio);
        } else {
            temperature = iterate_temperature(throat_condition, pressure_ratio, temperature, m_frozen_state);
        }
        gamma_s = m_gas.gamma_s();
        velocity = m_gas.isenthalpic_velocity();
        const double sonic_velocity = m_gas.speed_of_sound();
        const double ln_area_ratio = std::log(m_gas.area_per_mdot(velocity) / A_mdot_throat);

        // Area too small (or not finite) means the station is too close to the throat.
        if (!(ln_area_ratio > ln_expansion_ratio)) {
            upper = ln_pressure_ratio;
        } else {
            lower = ln_pressure_ratio;
        }

        const double dlogp_dlogA =
            gamma_s * velocity * velocity / (velocity * velocity - sonic_velocity * sonic_velocity);
        double next = ln_pressure_ratio + dlogp_dlogA * (ln_expansion_ratio - ln_area_ratio);
        if (!std::isfinite(next) || next <= lower || next >= upper) {
            next = 0.5 * (lower + upper);
        }
        step = next - ln_pressure_ratio;
        if (std::abs(step) <= abstol) {
            break;
        }
        ln_pressure_ratio = next;
    }

    ExpansionProperties final_props = m_gas.expansion_properties();
    return {true,
            final_props.gamma_s,
            final_props.dlogV_dlogP_T,
            final_props.dlogV_dlogT_P,
            m_gas.save_state(),
            final_props.pinned_transition};
}

NozzleStation Nozzle::solve_supersonic_area_expansion(
    const ThroatCondition& throat_condition, double expansion_ratio, double abstol) {
    m_gas.restore_state(throat_condition.state);

    determine_and_set_equilibrium_condition();

    double gamma_s = m_gas.gamma_s();

    double ln_pressure_ratio = 0;

    //correlations for initial guess
    if (expansion_ratio >= 2) {
        ln_pressure_ratio = gamma_s+1.4*std::log(expansion_ratio);
    }
    else if (expansion_ratio > 1.0001) {
        double throat_pressure_ratio = throat_condition.P_inlet/m_gas.pressure();
        double ln_throat_ratio = std::log(throat_pressure_ratio);
        double ln_Ae_At = std::log(expansion_ratio);
        ln_pressure_ratio = ln_throat_ratio + std::sqrt(3.294*ln_Ae_At*ln_Ae_At+1.535*ln_Ae_At);
    } else {
        //invalid expansion ratio
        throw_invalid_expansion_ratio(expansion_ratio, 1.0001);
    }

    return iterate_area_expansion(
        throat_condition,
        expansion_ratio,
        std::exp(ln_pressure_ratio),
        abstol);
}

NozzleStation Nozzle::iterate_area_expansion(
    const ThroatCondition& throat_condition,
    double expansion_ratio, double pressure_ratio_guess, double abstol) {

    double pressure_ratio = pressure_ratio_guess;
    double gamma_s = m_gas.gamma_s();
    double velocity = m_gas.isenthalpic_velocity();
    const double A_mdot_thrt = m_gas.area_per_mdot(velocity);

    double P_exit = throat_condition.P_inlet / pressure_ratio;
    bool is_equilibrium = determine_and_set_equilibrium_condition();
    
    if (is_equilibrium) {
        m_gas.equilibrate_SP(throat_condition.S_inlet, P_exit);
        gamma_s = m_gas.gamma_s();
    }
    else {
        m_gas.restore_state(m_frozen_state);
    }
    
    double T_exit = m_gas.temperature();
    std::vector<double> composition = m_gas.mole_fractions();    

    double Ae_At = m_gas.area_per_mdot(velocity)/A_mdot_thrt;

    int iters = 0;
    int max_iter = 10;
    double residual = 1.0;
    double sonic_velocity = 0.0;
    double dlogp_dlogA = 0.0;

    while (std::abs(residual) > abstol) {
        iters++;
        if (iters >= max_iter) {
            throw ConvergenceError("Maximum number of iterations exceeded for area expansion.", iters, abstol, std::abs(residual));
        }

        if (is_equilibrium) {
            P_exit = throat_condition.P_inlet / pressure_ratio;
            m_gas.equilibrate_SP(throat_condition.S_inlet, P_exit);
        } else {
            T_exit = iterate_temperature(throat_condition, pressure_ratio, T_exit, composition);
            if (T_exit < 0) {
                throw std::runtime_error("Negative temperature in area iteration: this branch should be unreachable.");
            }
        }

        gamma_s = m_gas.gamma_s();
        velocity = m_gas.isenthalpic_velocity();
        sonic_velocity = m_gas.speed_of_sound();
        Ae_At = m_gas.area_per_mdot(velocity)/A_mdot_thrt;

        dlogp_dlogA = gamma_s * velocity * velocity / (velocity*velocity - sonic_velocity*sonic_velocity);
        residual = dlogp_dlogA * (std::log(expansion_ratio) - std::log(Ae_At));
        double log_pinf_pe = std::log(pressure_ratio) + residual;

        pressure_ratio = std::exp(log_pinf_pe);
    }

    if (is_equilibrium) {
        m_frozen_state = m_gas.save_state();
    }

    ExpansionProperties final_props = m_gas.expansion_properties();
    return {true,
            final_props.gamma_s,
            final_props.dlogV_dlogP_T,
            final_props.dlogV_dlogT_P,
            m_gas.save_state(),
            final_props.pinned_transition};
}

NozzleStation Nozzle::solve_pressure_ratio(
    const ThroatCondition& throat_condition,
    double pressure_ratio,
    double abstol) {

    if (determine_and_set_equilibrium_condition()) {
        double P_exit = throat_condition.P_inlet/pressure_ratio;
        m_gas.equilibrate_SP(throat_condition.S_inlet, P_exit);
        m_frozen_state = m_gas.save_state();

        //pressure ratio for equilibrium nozzle does not require iteration
        ExpansionProperties final_props = m_gas.expansion_properties();
        return {true,
            final_props.gamma_s,
            final_props.dlogV_dlogP_T,
            final_props.dlogV_dlogT_P,
            m_gas.save_state(),
            final_props.pinned_transition};
    } else {
        m_gas.restore_state(m_frozen_state);

        double T_exit = iterate_temperature(throat_condition, pressure_ratio,
            m_gas.temperature(), m_gas.mole_fractions(), abstol);

        if (T_exit < 0) {
            throw std::runtime_error("Negative temperature returned in pressure ratio loop. This branch should be unreachable.");
        } else {
            ExpansionProperties props = m_gas.expansion_properties();
            return {true, props.gamma_s, props.dlogV_dlogP_T, props.dlogV_dlogT_P,
                    m_gas.save_state(), props.pinned_transition};
        }
    }
}

double Nozzle::iterate_temperature(
    const ThroatCondition& throat_condition,
    double pressure_ratio, double T_guess,
    const std::vector<double>& composition,
    double abstol) {

    const double P_exit = throat_condition.P_inlet/pressure_ratio;

    if (m_gas.has_condensed_candidates()) {
        // The condensed amounts are frozen along with the gas composition, so the isentrope is
        // the range-checked Newton iteration of `Gas::set_state_SP`, which stops where CEA stops
        // when a condensed species leaves its data range.
        m_gas.set_state_TPX(T_guess, P_exit, composition.data());
        try {
            m_gas.set_state_SP(throat_condition.S_inlet, P_exit);
        } catch (const FmtError& error) {
            rethrow_frozen_range_error(error, m_current_station, m_gas.temperature());
        }
        return m_gas.temperature();
    }

    double T_exit = T_guess;

    m_gas.set_state_TPX(T_exit, P_exit, composition.data());

    double Cp = m_gas.cp_mass();
    double dlnT = (throat_condition.S_inlet - m_gas.entropy_mass())/Cp;

    int maxiter = 8;
    int iters = 0;

    while (std::abs(dlnT) > abstol) {
        iters++;
        if (iters >= maxiter){
            throw ConvergenceError("Frozen flow temperature iteration failed.", iters, abstol, std::abs(dlnT));
        }

        double lnT_exit = std::log(T_exit) + dlnT;
        T_exit = std::exp(lnT_exit);

        m_gas.set_state_TPX(T_exit, P_exit, composition.data());
        Cp = m_gas.cp_mass();
        dlnT = (throat_condition.S_inlet - m_gas.entropy_mass())/Cp;
    }

    return T_exit;
}

bool Nozzle::determine_and_set_equilibrium_condition() {
    bool is_equilibrium = (m_opts.chemistry == GasChemistry::EQUILIBRIUM) 
        || (m_opts.chemistry == GasChemistry::FROZEN && m_current_station < m_opts.frozen_NFZ);
    if (is_equilibrium) {
        m_gas.chemistry = GasChemistry::EQUILIBRIUM;
    }
    else {
        m_gas.chemistry = GasChemistry::FROZEN;
    }
    return is_equilibrium;
}

void Nozzle::throw_invalid_expansion_ratio(double expansion, double min) const {
    throw std::invalid_argument(std::format("Specified expansion ratio must be greater than {}. Actual value: {}\n", min, expansion));
}

} //namespace Goddard
