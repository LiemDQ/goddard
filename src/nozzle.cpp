#include "goddard/nozzle.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/error.hpp"
#include "goddard/utils.hpp"
#include "goddard/gas_dynamics.hpp"

#include "cantera/core.h"
#include <cmath>
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
    : inlet_state(gas.save_state()), m_gas(gas), m_opts(options) {
    determine_equilibrium_condition();
    if (m_opts.chemistry == GasChemistry::KINETIC) {
        throw std::invalid_argument("GasChemistry::KINETIC is not valid for Nozzle. Use KineticNozzle instead.");
    }
    if (m_opts.chemistry == GasChemistry::PERFECT_GAS) {
        throw std::runtime_error("GasChemistry::PERFECT_GAS nozzle not yet implemented.");
    }
}

Nozzle::Nozzle(const Gas& gas, std::vector<double> state, NozzleOptions options)
    : inlet_state(std::move(state)), m_gas(gas), m_opts(options) {
    determine_equilibrium_condition();
    if (m_opts.chemistry == GasChemistry::KINETIC) {
        throw std::invalid_argument("GasChemistry::KINETIC is not valid for Nozzle. Use KineticNozzle instead.");
    }
    if (m_opts.chemistry == GasChemistry::PERFECT_GAS) {
        throw std::runtime_error("GasChemistry::PERFECT_GAS nozzle not yet implemented.");
    }
}

NozzleResults Nozzle::solve(ExpansionType expansion_type, double ratio) {

    const ThroatCondition throat_condition = solve_throat_conditions();
    std::vector<NozzleStation> result;
    switch (expansion_type) {
        case ExpansionType::SUPERSONIC_AREA_RATIO: {
            result.push_back(solve_supersonic_area_expansion(throat_condition, ratio));
            break;
        }
        case ExpansionType::SUBSONIC_AREA_RATIO: {
            result.push_back(solve_subsonic_area_expansion(throat_condition, ratio));
            break;
        }
        case ExpansionType::PRESSURE_RATIO: {
            result.push_back(solve_pressure_ratio(throat_condition, ratio));
            break;
        }
        default:
            throw NotImplementedError("Expansion type is not implemented.");
    }
    m_current_station++;
    return {throat_condition, result};
}

NozzleResults Nozzle::solve(ExpansionType expansion_type, const std::vector<double>& ratios) {

    const ThroatCondition throat_condition = solve_throat_conditions();
    
    std::vector<NozzleStation> results;

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
        }
        m_current_station++;
    }

    return {throat_condition, results};
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
    const bool is_equilibrium = determine_equilibrium_condition();
    
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

    double ln_pressure_ratio = 0;
    double throat_pressure_ratio = throat_condition.P_inlet/m_gas.pressure();
    double ln_throat_ratio = std::log(throat_pressure_ratio);
    double ln_Ae_At = std::log(expansion_ratio);

    if (expansion_ratio >= 1.09) {
        ln_pressure_ratio = ln_throat_ratio/(expansion_ratio + 10.587*std::pow(ln_Ae_At, 3)+9.454*ln_Ae_At);
    }
    else if (expansion_ratio > 1.0001) {
        ln_pressure_ratio = 0.9* ln_throat_ratio / (expansion_ratio + 10.587*std::pow(ln_Ae_At, 3)+9.454*ln_Ae_At);
    }
    else {
        //invalid expansion ratio
        throw_invalid_expansion_ratio(expansion_ratio, 1.0001);
    }
    return iterate_area_expansion(throat_condition, expansion_ratio, std::exp(ln_pressure_ratio), abstol);
}

NozzleStation Nozzle::solve_supersonic_area_expansion(
    const ThroatCondition& throat_condition, double expansion_ratio, double abstol) {
    m_gas.restore_state(throat_condition.state);

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

    double T_exit = m_gas.temperature();
    std::vector<double> composition = m_gas.mole_fractions();    
    double P_exit = throat_condition.P_inlet / pressure_ratio;
    bool is_equilibrium = determine_equilibrium_condition();

    if (is_equilibrium) {
        m_gas.equilibrate_SP(throat_condition.S_inlet, P_exit);
        gamma_s = m_gas.gamma_s();
    } 
    
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

    m_gas.restore_state(throat_condition.state);

    if (determine_equilibrium_condition()) {
        double P_exit = throat_condition.P_inlet/pressure_ratio;
        m_gas.equilibrate_SP(throat_condition.S_inlet, P_exit);

        //pressure ratio for equilibrium nozzle does not require iteration
        ExpansionProperties final_props = m_gas.expansion_properties();
        return {true,
            final_props.gamma_s,
            final_props.dlogV_dlogP_T,
            final_props.dlogV_dlogT_P,
            m_gas.save_state(),
            final_props.pinned_transition};
    } else {
        std::vector<double> composition = m_gas.mole_fractions();

        double T_exit = iterate_temperature(throat_condition, pressure_ratio,
            m_gas.temperature(), composition, abstol);
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

bool Nozzle::determine_equilibrium_condition() {
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
