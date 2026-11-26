#include "goddard/nozzle.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/error.hpp"
#include "goddard/utils.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/newton.hpp"
#include "goddard/solve_guard.hpp"

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
 * The temperature-range failure of a frozen station solve, naming the station.
 *
 * A frozen expansion holds the condensed amounts fixed, so it cannot be carried past the
 * temperature range of a condensed species that is present. `Gas` reports that as a `SolveError`
 * with reason OUT_OF_DATA_RANGE; the nozzle adds the station index and the temperature [K] it
 * stopped at, as CEA does. Call it only for a failure with that reason.
 *
 * @param error Failure raised by `Gas`.
 * @param station Index of the station being solved [-].
 * @param temperature Temperature the frozen iteration stopped at [K].
 */
SolveError frozen_range_error(const SolveError& error, int station, double temperature) {
    return SolveError(FailureReason::OUT_OF_DATA_RANGE, std::format(
        "Frozen expansion: {} at station {} (T = {:.6g} K)", error.status().message, station,
        temperature));
}

/** Name of an expansion station in a status, e.g. "station 2 (area ratio 40)". */
std::string station_stage(int station, ExpansionType expansion_type, double ratio) {
    switch (expansion_type) {
        case ExpansionType::SUPERSONIC_AREA_RATIO:
            return std::format("station {} (supersonic area ratio {:.6g})", station, ratio);
        case ExpansionType::SUBSONIC_AREA_RATIO:
            return std::format("station {} (subsonic area ratio {:.6g})", station, ratio);
        case ExpansionType::PRESSURE_RATIO:
            return std::format("station {} (pressure ratio {:.6g})", station, ratio);
    }
    return std::format("station {}", station);
}

void validate_nozzle_options(const NozzleOptions& options) {
    if (options.chemistry == GasChemistry::KINETIC) {
        throw std::invalid_argument("GasChemistry::KINETIC is not valid for Nozzle. Use KineticNozzle instead.");
    }
    if (options.chemistry == GasChemistry::PERFECT_GAS) {
        throw NotImplementedError("GasChemistry::PERFECT_GAS nozzle not yet implemented.");
    }
    if (options.frozen_NFZ < 0) {
        throw std::invalid_argument(std::format(
            "Nozzle: frozen_NFZ must not be negative, got {}. 0 freezes the flow at the chamber.",
            options.frozen_NFZ));
    }
}

} // namespace

NozzleStation NozzleStation::failed(SolveStatus status) {
    NozzleStation station;
    station.thermo = ThermodynamicState::undefined();
    station.velocity = NaN;
    station.mach = NaN;
    station.area_ratio = NaN;
    station.status = std::move(status);
    return station;
}

ThroatCondition ThroatCondition::failed(SolveStatus status) {
    ThroatCondition throat;
    throat.thermo = ThermodynamicState::undefined();
    throat.velocity = NaN;
    throat.mach = NaN;
    throat.area_ratio = NaN;
    throat.status = std::move(status);
    set_stage_if_unset(throat.status, "throat");
    return throat;
}

FiniteAreaChamber FiniteAreaChamber::failed(SolveStatus status) {
    FiniteAreaChamber chamber;
    // The stagnation state, combustion end and throat are outputs of the one coupled chamber
    // solve, so they fail with its status.
    chamber.injector = NozzleStation::failed(status);
    chamber.stagnation = NozzleStation::failed(status);
    chamber.combustion_end = NozzleStation::failed(status);
    chamber.throat = ThroatCondition::failed(status);
    chamber.status = std::move(status);
    return chamber;
}

SolveStatus NozzleResults::first_failure() const {
    if (!inlet.status.ok()) {
        return inlet.status;
    }
    if (!throat.status.ok()) {
        return throat.status;
    }
    for (const NozzleStation& station : expansions) {
        if (!station.status.ok()) {
            return station.status;
        }
    }
    return {};
}

NozzleResults NozzleResults::failed(const SolveStatus& inlet_status, std::size_t num_stations) {
    NozzleResults results;
    const SolveStatus not_computed = SolveStatus::not_computed(inlet_status);
    results.inlet = NozzleStation::failed(inlet_status);
    results.throat = ThroatCondition::failed(not_computed);
    for (std::size_t i = 0; i < num_stations; i++) {
        SolveStatus status = not_computed;
        status.stage = std::format("station {}", i + 2);
        results.expansions.push_back(NozzleStation::failed(std::move(status)));
    }
    results.status = inlet_status;
    return results;
}

Nozzle::Nozzle(const Gas& gas, NozzleOptions options)
    : m_inlet_state(gas.save_state()), m_gas(gas.clone()), m_opts(options) {
    validate_nozzle_options(m_opts);
    set_station_chemistry(0);
}

Nozzle::Nozzle(const Gas& gas, std::vector<double> state, NozzleOptions options)
    : m_inlet_state(std::move(state)), m_gas(gas.clone()), m_opts(options) {
    validate_nozzle_options(m_opts);
    set_station_chemistry(0);
}

void Nozzle::validate_stations(ExpansionType expansion_type, const std::vector<double>& ratios) const {
    validate(m_opts, expansion_type, ratios);
}

void Nozzle::validate(const NozzleOptions& options, ExpansionType expansion_type,
                      const std::vector<double>& ratios) {
    validate_nozzle_options(options);
    // Stations 0 and 1 are the chamber and the throat.
    const size_t last_station = ratios.size() + 1;
    if (options.chemistry == GasChemistry::FROZEN
        && static_cast<size_t>(options.frozen_NFZ) > last_station) {
        throw std::invalid_argument(std::format(
            "Nozzle: frozen_NFZ = {} lies beyond the last station ({}), so the flow would never "
            "freeze. Station 0 is the chamber, 1 the throat and 2 onward the expansion stations.",
            options.frozen_NFZ, last_station));
    }
    for (double ratio : ratios) {
        switch (expansion_type) {
            case ExpansionType::SUPERSONIC_AREA_RATIO:
            case ExpansionType::SUBSONIC_AREA_RATIO:
                if (!(ratio > 1.0001)) {
                    throw_invalid_expansion_ratio(ratio, 1.0001);
                }
                break;
            case ExpansionType::PRESSURE_RATIO:
                if (!(ratio >= 1.0)) {
                    throw std::invalid_argument(std::format(
                        "Nozzle: a pressure ratio is the inlet pressure over the station pressure "
                        "and must be at least 1. Actual value: {:.6g}", ratio));
                }
                break;
            default:
                throw NotImplementedError("Expansion type is not implemented.");
        }
    }
}

NozzleResults Nozzle::solve_inlet_and_throat() {
    NozzleResults results;
    const SolveStatus inlet_status = capture_failure("inlet", [&] {
        results.inlet = station_at_state(m_inlet_state, 0);
    });
    if (!inlet_status.ok()) {
        results.inlet = NozzleStation::failed(inlet_status);
        results.throat = ThroatCondition::failed(SolveStatus::not_computed(inlet_status));
        return results;
    }
    results.throat = solve_throat_conditions();
    return results;
}

NozzleResults Nozzle::solve(ExpansionType expansion_type, double ratio) {
    return solve(expansion_type, std::vector<double>{ratio});
}

NozzleResults Nozzle::solve(ExpansionType expansion_type, const std::vector<double>& ratios) {
    validate_stations(expansion_type, ratios);
    NozzleResults results = solve_inlet_and_throat();
    results.expansions = solve_stations(results.throat, expansion_type, ratios);
    results.status = results.first_failure();
    return results;
}

NozzleResults Nozzle::solve(const NozzleProfile& profile, int num_stations) {
    if (num_stations < 1) {
        throw std::invalid_argument(std::format(
            "Nozzle::solve: num_stations must be at least 1, got {}.", num_stations));
    }
    // For a diverging-only profile (starting at throat), x_min is the throat.
    double A_throat = profile.area_at(profile.x_min());

    double x_start = profile.x_min();
    double x_end = profile.x_max();
    // Slightly inset from x_max to avoid out-of-bounds interpolation at the boundary
    double x_range = x_end - x_start;
    double x_last = x_end - 1e-10 * x_range;

    std::vector<double> area_ratios;
    area_ratios.reserve(static_cast<size_t>(num_stations));
    for (int i = 1; i <= num_stations; i++) {
        double x = x_start + (x_last - x_start) * static_cast<double>(i) / num_stations;
        area_ratios.push_back(profile.area_at(x) / A_throat);
    }
    return solve(ExpansionType::SUPERSONIC_AREA_RATIO, area_ratios);
}

std::vector<NozzleStation> Nozzle::solve_stations(const ThroatCondition& throat_condition,
    ExpansionType expansion_type, const std::vector<double>& ratios) {
    validate_stations(expansion_type, ratios);
    return solve_stations_unchecked(throat_condition, expansion_type, ratios);
}

std::vector<NozzleStation> Nozzle::solve_stations_unchecked(const ThroatCondition& throat_condition,
    ExpansionType expansion_type, const std::vector<double>& ratios) {
    std::vector<NozzleStation> results;
    results.reserve(ratios.size());
    // Stations 0 and 1 are the chamber and the throat. A frozen station keeps the composition of
    // the last equilibrium station, which is the throat until an equilibrium station follows it.
    int station = 2;
    std::vector<double> frozen_state = throat_condition.state;
    SolveStatus frozen_state_status = throat_condition.status;
    for (double ratio: ratios) {
        const std::string stage = station_stage(station, expansion_type, ratio);
        NozzleStation result;
        SolveStatus status;
        if (!throat_condition.status.ok()) {
            status = SolveStatus::not_computed(throat_condition.status);
        } else if (!is_equilibrium_station(station) && !frozen_state_status.ok()) {
            // The station it would be frozen from has no state.
            status = SolveStatus::not_computed(frozen_state_status);
        } else {
            status = capture_failure(stage, [&] {
                if (expansion_type == ExpansionType::PRESSURE_RATIO && !(ratio >= 1.0)) {
                    throw SolveError(FailureReason::NO_SOLUTION, std::format(
                        "the station pressure would exceed the stagnation pressure: the pressure "
                        "ratio is {:.6g}, below 1.", ratio));
                }
                switch (expansion_type) {
                    case ExpansionType::SUPERSONIC_AREA_RATIO: {
                        result = solve_supersonic_area_expansion(
                            throat_condition, ratio, station, frozen_state);
                        break;
                    }
                    case ExpansionType::SUBSONIC_AREA_RATIO: {
                        result = solve_subsonic_area_expansion(
                            throat_condition, ratio, station, frozen_state);
                        break;
                    }
                    case ExpansionType::PRESSURE_RATIO: {
                        result = solve_pressure_ratio(throat_condition, ratio, station, frozen_state);
                        break;
                    }
                    default:
                        throw NotImplementedError("Expansion type is not implemented.");
                }
            });
            if (!status.ok()) {
                // The next station must not start from the state a failed solve left behind.
                m_gas.restore_state(throat_condition.state);
            }
        }
        if (!status.ok()) {
            set_stage_if_unset(status, stage);
            result = NozzleStation::failed(std::move(status));
        }
        if (is_equilibrium_station(station)) {
            frozen_state = result.state;
            frozen_state_status = result.status;
        }
        results.push_back(std::move(result));
        station++;
    }
    return results;
}

std::vector<NozzleStation> Nozzle::solve_stations(const FiniteAreaChamber& chamber,
    ExpansionType expansion_type, const std::vector<double>& ratios) {

    validate_stations(expansion_type, ratios);
    // The expansion runs from the stagnation state, so pressure ratios are converted from
    // P_inj/P to P_inf/P. A ratio between 1 and P_inj/P_inf then lies below 1 and its station
    // fails with NO_SOLUTION.
    std::vector<double> stagnation_ratios = ratios;
    if (expansion_type == ExpansionType::PRESSURE_RATIO && chamber.status.ok()) {
        const double P_inf_over_P_inj = chamber.throat.P_inlet / chamber.injector_pressure;
        for (double& ratio : stagnation_ratios) {
            ratio *= P_inf_over_P_inj;
        }
    }
    // A failed chamber has a failed throat, which makes every station NOT_COMPUTED.
    return solve_stations_unchecked(chamber.throat, expansion_type, stagnation_ratios);
}

FiniteAreaChamber Nozzle::solve_finite_area_chamber(const std::vector<double>& injector_state,
    CombustorType type, double value, double reltol) {

    bool contraction_mode = true;
    switch (type) {
        case CombustorType::FINITE_CONTRACTION_RATIO: {
            if (!(value > 1.0)) {
                throw std::invalid_argument(std::format(
                    "Finite-area combustor: contraction ratio must be greater than 1. Actual value: {:.6g}",
                    value));
            }
            contraction_mode = true;
            break;
        }
        case CombustorType::FINITE_MASS_FLUX: {
            if (!(value > 0.0)) {
                throw std::invalid_argument(std::format(
                    "Finite-area combustor: mass flux must be positive. Actual value: {:.6g} kg/(m^2 s)",
                    value));
            }
            contraction_mode = false;
            break;
        }
        case CombustorType::INFINITE_AREA:
        default:
            throw std::invalid_argument(
                "Nozzle::solve_finite_area_chamber requires FINITE_CONTRACTION_RATIO or FINITE_MASS_FLUX.");
    }
    if (m_opts.chemistry == GasChemistry::FROZEN && m_opts.frozen_NFZ == 0) {
        throw NotImplementedError(
            "Finite-area combustor with flow frozen at the combustion end (frozen_NFZ = 0) is not implemented.");
    }

    FiniteAreaChamber chamber;
    const std::vector<double> saved_inlet_state = m_inlet_state;
    SolveStatus status = capture_failure("finite-area chamber", [&] {
        chamber = finite_area_chamber_or_throw(injector_state, contraction_mode, value, reltol);
    });
    if (!status.ok()) {
        // The iteration left a trial stagnation state as the inlet state.
        set_inlet_state(saved_inlet_state);
        chamber = FiniteAreaChamber::failed(std::move(status));
        // The injector face is the given state and does not depend on the chamber solve.
        NozzleStation injector;
        if (capture_failure("injector", [&] { injector = station_at_state(injector_state, 0); }).ok()) {
            chamber.injector = std::move(injector);
        }
    }
    return chamber;
}

FiniteAreaChamber Nozzle::finite_area_chamber_or_throw(const std::vector<double>& injector_state,
    bool contraction_mode, double value, double reltol) {

    // Injector face: equilibrium at (h_inj, P_inj).
    NozzleStation injector = station_at_state(injector_state, 0);
    const double h_injector = m_gas.enthalpy_mass();
    const double P_injector = m_gas.pressure();
    const double T_injector = m_gas.temperature();
    const double gamma_injector = m_gas.gamma_s(m_station_model);
    const double R_specific = Cantera::GasConstant / m_gas.molecular_weight();

    // Perfect-gas initial guess for P_inf: find the combustion-end Mach number that the
    // contraction ratio or mass flux implies, then P_inf = P_inj / finite_area_pressure_loss(M).
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
            throw SolveError(FailureReason::NO_SOLUTION, std::format(
                "mass flux {:.6g} kg/(m^2 s) thermally chokes the chamber. "
                "Maximum mass flux (perfect gas estimate): {:.6g} kg/(m^2 s)", value, G_max));
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

    for (int iter = 1; iter <= max_iters; iter++) {
        // Stagnation state and throat; in mass-flux mode also the contraction ratio, backtracking
        // toward the previous iterate if the chamber would be choked.
        for (int backtrack = 0; ; backtrack++) {
            P_stagnation = std::exp(ln_P);
            m_gas.restore_state(injector_state);
            m_gas.equilibrate_HP(h_injector, P_stagnation);
            set_inlet_state(m_gas.save_state());
            throat = throat_conditions_or_throw(4e-5);
            if (contraction_mode) {
                break;
            }
            m_gas.restore_state(throat.state);
            contraction_ratio = m_gas.density() * throat.speed_of_sound / value;
            if (contraction_ratio > 1.0001) {
                break;
            }
            if (!have_previous || backtrack >= max_backtracks) {
                throw SolveError(FailureReason::NO_SOLUTION, std::format(
                    "mass flux {:.6g} kg/(m^2 s) thermally chokes the chamber "
                    "(A_c/A_t = {:.6g} at P_inf = {:.6g} Pa).", value, contraction_ratio, P_stagnation));
            }
            ln_P = 0.5 * (ln_P + ln_P_previous);
        }

        // Combustion end is part of the chamber and therefore in equilibrium.
        combustion_end = solve_subsonic_area_expansion(
            throat, contraction_ratio, 0, throat.state, station_abstol);
        const double velocity = m_gas.isenthalpic_velocity(throat.H_stagnation);
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
            SolveStatus status;
            if (enthalpy_noise > 0.01 * velocity * velocity) {
                status.warnings.push_back(std::format(
                    "Combustion-end velocity ({:.6g} m/s) and Mach number are below the solver's "
                    "resolution at A_c/A_t = {:.6g}; treat them as zero.", velocity, contraction_ratio));
            }
            NozzleStation stagnation = station_at_state(m_inlet_state, 0);
            return {std::move(injector), std::move(stagnation), combustion_end, throat, P_injector, P_stagnation,
                contraction_ratio, mass_flux, iter, std::move(status)};
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
    throw ConvergenceError(std::format(
        "the chamber momentum balance did not converge (last stagnation pressure {:.6g} Pa for an "
        "injector pressure of {:.6g} Pa).", P_stagnation, P_injector),
        max_iters, reltol, std::abs(residual));
}

void Nozzle::reset_state(){
    m_gas.restore_state(m_inlet_state);
}

ThroatCondition Nozzle::solve_throat_conditions(double abstol) {
    ThroatCondition condition;
    SolveStatus status = capture_failure("throat", [&] {
        if (is_equilibrium_station(1)) {
            require_equilibrium_inlet();
        }
        condition = throat_conditions_or_throw(abstol);
    });
    if (!status.ok()) {
        condition = ThroatCondition::failed(std::move(status));
    }
    return condition;
}

void Nozzle::require_equilibrium_inlet() {
    m_gas.restore_state(m_inlet_state);
    const double inlet_temperature = m_gas.temperature();
    m_gas.equilibrate_HP(m_gas.enthalpy_mass(), m_gas.pressure());
    const double equilibrium_temperature = m_gas.temperature();
    m_gas.restore_state(m_inlet_state);

    const double rtol = 1e-5;
    if (std::abs(equilibrium_temperature - inlet_temperature) > rtol * inlet_temperature) {
        throw std::invalid_argument(std::format(
            "Nozzle: the throat is an equilibrium station, so the inlet must be at chemical "
            "equilibrium, but re-equilibrating it at its enthalpy and pressure moves its "
            "temperature from {:.6g} K to {:.6g} K. Equilibrate the inlet first, or use FROZEN "
            "chemistry with frozen_NFZ = 0.", inlet_temperature, equilibrium_temperature));
    }
}

ThroatCondition Nozzle::throat_conditions_or_throw(double abstol) {
    m_gas.restore_state(m_inlet_state);
    const int station = 1;

    double P_inlet = m_gas.pressure();
    double S_inlet = m_gas.entropy_mass();
    double H_inlet = m_gas.enthalpy_mass();
    const bool is_equilibrium = set_station_chemistry(station);

    double gamma_s = m_gas.gamma_s(m_station_model);
    double P_throat = P_inlet / std::pow((gamma_s+1)/2,gamma_s/(gamma_s-1));

    const int max_iters = 10;
    int iter = 0;
    double Mach = 1.0; //throat mach number is 1 by definition
    double residual = 1.0;

    while (residual > abstol){
        if (iter >= max_iters){
            throw ConvergenceError(std::format(
                "the throat iteration did not reach Mach 1 (last Mach number {:.6g}, at P = {:.6g} Pa "
                "for an inlet pressure of {:.6g} Pa).", Mach, P_throat, P_inlet),
                iter, abstol, residual);
        }
        P_throat = P_throat * (1 + gamma_s * Mach* Mach)/(1+ gamma_s);

        //if equilibrium conditions are selected, the composition must reach chemical
        //equilibrium in the throat.
        if (is_equilibrium) {
            m_gas.equilibrate_SP(S_inlet, P_throat);
        } else {
            try {
                m_gas.set_state_SP(S_inlet, P_throat);
            } catch (const SolveError& error) {
                if (error.reason() == FailureReason::OUT_OF_DATA_RANGE) {
                    throw frozen_range_error(error, station, m_gas.temperature());
                }
                throw;
            }
        }
        gamma_s = m_gas.gamma_s(m_station_model);

        double velocity = m_gas.isenthalpic_velocity(H_inlet);
        double sonic_velocity = m_gas.speed_of_sound(m_station_model);

        Mach = velocity/sonic_velocity;
        residual = std::abs(1.0 - 1.0/(Mach * Mach));

        iter++;
    }

    ExpansionProperties final_props = m_gas.expansion_properties(m_station_model);
    NozzleStation throat = station_at_current_state(final_props, H_inlet, 0.0);

    ThroatCondition condition;
    condition.speed_of_sound = m_gas.speed_of_sound(m_station_model);
    condition.H_stagnation = H_inlet;
    condition.P_inlet = P_inlet;
    condition.S_inlet = S_inlet;
    condition.gamma_s = gamma_s;
    condition.dlV_dlP_T = final_props.dlV_dlP_T;
    condition.dlV_dlT_P = final_props.dlV_dlT_P;
    condition.state = std::move(throat.state);
    condition.pinned_transition = final_props.pinned_transition;
    condition.thermo = std::move(throat.thermo);
    condition.velocity = throat.velocity;
    condition.mach = throat.mach;
    condition.area_ratio = 1.0;
    return condition;
}

double Nozzle::get_gamma_s() {
    return m_gas.gamma_s(m_station_model);
}

NozzleStation Nozzle::solve_subsonic_area_expansion(const ThroatCondition& throat_condition,
    double expansion_ratio, int station, const std::vector<double>& frozen_state, double abstol) {
    m_gas.restore_state(throat_condition.state);

    if (!(expansion_ratio > 1.0001)) {
        throw_invalid_expansion_ratio(expansion_ratio, 1.0001);
    }
    // Perfect-gas initial guess with the throat isentropic exponent. Unlike a correlation it stays
    // accurate as A/A_t grows and the station approaches stagnation.
    double gamma_s = m_gas.gamma_s(m_station_model);
    const double mach_guess = mach_from_area_ratio(expansion_ratio, gamma_s, false);
    double ln_pressure_ratio =
        gamma_s / (gamma_s - 1.0) * std::log(stagnation_factor(mach_guess, gamma_s));

    // Safeguarded Newton iteration on x = ln(P_inlet/P). On the subsonic branch A/A_t decreases
    // monotonically from infinity at x = 0 to 1 at the throat, so [0, x_throat] brackets the root.
    // Pure Newton diverges near the throat, where dlnA/dlnP vanishes, and near stagnation, where
    // x is tiny; steps leaving the bracket fall back to bisection.
    double velocity = m_gas.isenthalpic_velocity(throat_condition.H_stagnation);
    const double A_mdot_throat = m_gas.area_per_mdot(velocity);
    double temperature = m_gas.temperature();
    const bool is_equilibrium = set_station_chemistry(station);
    const double ln_throat_ratio = std::log(throat_condition.P_inlet / m_gas.pressure());
    std::vector<double> composition;
    if (!is_equilibrium) {
        m_gas.restore_state(frozen_state);
        composition = m_gas.mole_fractions();
    }
    const double ln_expansion_ratio = std::log(expansion_ratio);

    double lower = 0.0;
    double upper = ln_throat_ratio;
    if (!(ln_pressure_ratio > lower && ln_pressure_ratio < upper)) {
        ln_pressure_ratio = 0.5 * (lower + upper);
    }

    NewtonOptions options;
    options.step_abstol = abstol;
    options.max_iterations = 60;
    double ln_area_ratio = 0.0;
    const NewtonResult result = newton_solve(ln_pressure_ratio, [&](double x) {
        const double pressure_ratio = std::exp(x);
        if (is_equilibrium) {
            m_gas.equilibrate_SP(throat_condition.S_inlet, throat_condition.P_inlet / pressure_ratio);
        } else {
            temperature = iterate_temperature(
                throat_condition, pressure_ratio, temperature, composition, station);
        }
        gamma_s = m_gas.gamma_s(m_station_model);
        velocity = m_gas.isenthalpic_velocity(throat_condition.H_stagnation);
        const double sonic_velocity = m_gas.speed_of_sound(m_station_model);
        ln_area_ratio = std::log(m_gas.area_per_mdot(velocity) / A_mdot_throat);
        // (dlnA/dlnP)_s = (u^2 - a^2)/(gamma_s u^2), and x increases as P falls
        return NewtonFunction{
            .value = ln_area_ratio - ln_expansion_ratio,
            .derivative = (velocity * velocity - sonic_velocity * sonic_velocity)
                / (gamma_s * velocity * velocity)};
    }, [&](int, double x, double step) {
        // Area too small (or not finite) means the station is too close to the throat.
        if (!(ln_area_ratio > ln_expansion_ratio)) {
            upper = x;
        } else {
            lower = x;
        }
        const double next = x + step;
        if (!std::isfinite(next) || next <= lower || next >= upper) {
            return 0.5 * (lower + upper) - x;
        }
        return step;
    }, options);
    if (result.status != NewtonStatus::CONVERGED) {
        throw ConvergenceError(std::format(
            "the subsonic area-ratio iteration did not converge to A/A_t = {:.6g} (last pressure "
            "ratio {:.6g}).", expansion_ratio, std::exp(result.x)),
            result.iterations, abstol, std::abs(result.step));
    }

    return station_at_current_state(m_gas.expansion_properties(m_station_model), throat_condition.H_stagnation,
        A_mdot_throat);
}

NozzleStation Nozzle::solve_supersonic_area_expansion(const ThroatCondition& throat_condition,
    double expansion_ratio, int station, const std::vector<double>& frozen_state, double abstol) {
    m_gas.restore_state(throat_condition.state);

    set_station_chemistry(station);

    double gamma_s = m_gas.gamma_s(m_station_model);

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
        station,
        frozen_state,
        abstol);
}

NozzleStation Nozzle::iterate_area_expansion(
    const ThroatCondition& throat_condition,
    double expansion_ratio, double pressure_ratio_guess, int station,
    const std::vector<double>& frozen_state, double abstol) {

    double velocity = m_gas.isenthalpic_velocity(throat_condition.H_stagnation);
    const double A_mdot_thrt = m_gas.area_per_mdot(velocity);

    const bool is_equilibrium = set_station_chemistry(station);
    if (!is_equilibrium) {
        m_gas.restore_state(frozen_state);
    }
    double T_exit = m_gas.temperature();
    const std::vector<double> composition = m_gas.mole_fractions();
    const double ln_expansion_ratio = std::log(expansion_ratio);

    NewtonOptions options;
    options.step_abstol = abstol;
    options.max_iterations = 9;
    const NewtonResult result = newton_solve(std::log(pressure_ratio_guess), [&](double ln_pressure_ratio) {
        const double pressure_ratio = std::exp(ln_pressure_ratio);
        if (is_equilibrium) {
            m_gas.equilibrate_SP(throat_condition.S_inlet, throat_condition.P_inlet / pressure_ratio);
        } else {
            T_exit = iterate_temperature(
                throat_condition, pressure_ratio, T_exit, composition, station);
            if (T_exit < 0) {
                throw std::logic_error("Negative temperature in area iteration: this branch should be unreachable.");
            }
        }

        const double gamma_s = m_gas.gamma_s(m_station_model);
        velocity = m_gas.isenthalpic_velocity(throat_condition.H_stagnation);
        const double sonic_velocity = m_gas.speed_of_sound(m_station_model);
        const double Ae_At = m_gas.area_per_mdot(velocity)/A_mdot_thrt;
        // (dlnA/dlnP)_s = (u^2 - a^2)/(gamma_s u^2), and ln(P_inf/P_e) increases as P_e falls
        return NewtonFunction{
            .value = std::log(Ae_At) - ln_expansion_ratio,
            .derivative = (velocity*velocity - sonic_velocity*sonic_velocity)/(gamma_s*velocity*velocity)};
    }, options);
    if (result.status != NewtonStatus::CONVERGED) {
        throw ConvergenceError(std::format(
            "the area-ratio iteration did not converge to A/A_t = {:.6g} (last pressure ratio {:.6g}).",
            expansion_ratio, std::exp(result.x)),
            result.iterations, abstol, std::abs(result.step));
    }

    return station_at_current_state(m_gas.expansion_properties(m_station_model), throat_condition.H_stagnation,
        A_mdot_thrt);
}

NozzleStation Nozzle::solve_pressure_ratio(
    const ThroatCondition& throat_condition,
    double pressure_ratio,
    int station,
    const std::vector<double>& frozen_state,
    double abstol) {

    if (set_station_chemistry(station)) {
        double P_exit = throat_condition.P_inlet/pressure_ratio;
        m_gas.equilibrate_SP(throat_condition.S_inlet, P_exit);

        //pressure ratio for equilibrium nozzle does not require iteration
        return station_at_current_state(m_gas.expansion_properties(m_station_model),
            throat_condition.H_stagnation, throat_area_per_mdot(throat_condition));
    } else {
        m_gas.restore_state(frozen_state);

        double T_exit = iterate_temperature(throat_condition, pressure_ratio,
            m_gas.temperature(), m_gas.mole_fractions(), station, abstol);

        if (T_exit < 0) {
            throw std::logic_error("Negative temperature returned in pressure ratio loop. This branch should be unreachable.");
        } else {
            return station_at_current_state(m_gas.expansion_properties(m_station_model),
                throat_condition.H_stagnation, throat_area_per_mdot(throat_condition));
        }
    }
}

double Nozzle::iterate_temperature(
    const ThroatCondition& throat_condition,
    double pressure_ratio, double T_guess,
    const std::vector<double>& composition,
    int station,
    double abstol) {

    const double P_exit = throat_condition.P_inlet/pressure_ratio;

    if (m_gas.has_condensed_candidates()) {
        // The condensed amounts are frozen along with the gas composition, so the isentrope is
        // the range-checked Newton iteration of `Gas::set_state_SP`, which stops where CEA stops
        // when a condensed species leaves its data range.
        m_gas.set_state_TPX(T_guess, P_exit, composition.data());
        try {
            m_gas.set_state_SP(throat_condition.S_inlet, P_exit);
        } catch (const SolveError& error) {
            if (error.reason() == FailureReason::OUT_OF_DATA_RANGE) {
                throw frozen_range_error(error, station, m_gas.temperature());
            }
            throw;
        }
        return m_gas.temperature();
    }

    NewtonOptions options;
    options.step_abstol = abstol;
    options.max_iterations = 8;
    const NewtonResult result = newton_solve(std::log(T_guess), [&](double ln_T) {
        m_gas.set_state_TPX(std::exp(ln_T), P_exit, composition.data());
        // (ds/dlnT)_P = cp
        return NewtonFunction{.value = m_gas.entropy_mass() - throat_condition.S_inlet,
                              .derivative = m_gas.cp_mass()};
    }, options);
    if (result.status != NewtonStatus::CONVERGED) {
        throw ConvergenceError(std::format(
            "the frozen temperature iteration did not converge at P = {:.6g} Pa (last temperature "
            "{:.6g} K).", P_exit, std::exp(result.x)),
            result.iterations, abstol, std::abs(result.step));
    }
    return std::exp(result.x);
}

bool Nozzle::is_equilibrium_station(int station) const {
    // A station is reached from the one upstream of it, so it is in equilibrium while the flow
    // there still is. The chamber has no upstream station and takes the chemistry of the flow
    // leaving it: frozen only when the flow freezes at the chamber.
    const int upstream = std::max(station - 1, 0);
    return (m_opts.chemistry == GasChemistry::EQUILIBRIUM)
        || (m_opts.chemistry == GasChemistry::FROZEN && upstream < m_opts.frozen_NFZ);
}

bool Nozzle::set_station_chemistry(int station) {
    const bool is_equilibrium = is_equilibrium_station(station);
    m_station_model = is_equilibrium ? CompositionModel::EQUILIBRIUM : CompositionModel::FROZEN;
    return is_equilibrium;
}

NozzleStation Nozzle::station_at_state(const std::vector<double>& state, int station) {
    m_gas.restore_state(state);
    set_station_chemistry(station);
    // The flow is at rest here: its own enthalpy is the stagnation enthalpy.
    return station_at_current_state(m_gas.expansion_properties(m_station_model), m_gas.enthalpy_mass(), 0.0);
}

NozzleStation Nozzle::station_at_current_state(const ExpansionProperties& props,
    double H_stagnation, double area_per_mdot_throat) {
    // A frozen snapshot reads the state without solving the equilibrium derivative system a
    // second time; the derivatives come from `props`, in the chemistry of the station.
    ThermodynamicState thermo = m_gas.snapshot(CompositionModel::FROZEN);
    thermo.composition_model = m_station_model;
    thermo.gamma_s = props.gamma_s;
    thermo.dlV_dlP_T = props.dlV_dlP_T;
    thermo.dlV_dlT_P = props.dlV_dlT_P;
    thermo.pinned_transition = props.pinned_transition;
    thermo.cp = props.pinned_transition ? 0.0 : props.spec_heat_p; // 0 following CEA convention
    thermo.stagnation_enthalpy = H_stagnation;
    thermo.speed_of_sound = gas_sonic_velocity(thermo.temperature, thermo.molecular_weight,
        props.gamma_s);

    NozzleStation result;
    result.velocity = std::sqrt(std::max(0.0, 2.0 * (H_stagnation - thermo.enthalpy)));
    result.mach = result.velocity / thermo.speed_of_sound;
    if (area_per_mdot_throat > 0.0 && result.velocity > 0.0) {
        result.area_ratio = m_gas.area_per_mdot(result.velocity) / area_per_mdot_throat;
    }
    result.thermo = std::move(thermo);
    result.state = m_gas.save_state();
    return result;
}

double Nozzle::throat_area_per_mdot(const ThroatCondition& throat) {
    return throat.thermo.temperature * Cantera::GasConstant
        / (throat.thermo.pressure * throat.velocity * throat.thermo.molecular_weight);
}

void Nozzle::throw_invalid_expansion_ratio(double expansion, double min) {
    throw std::invalid_argument(std::format("Specified expansion ratio must be greater than {}. Actual value: {:.6g}\n", min, expansion));
}

} //namespace Goddard
