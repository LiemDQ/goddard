#include <algorithm>
#include <cmath>
#include <format>
#include <stdexcept>
#include <string>
#include <vector>
#include "goddard/error.hpp"
#include "goddard/moc_thermo.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/gas.hpp"
#include "goddard/nozzle.hpp"

namespace Goddard {

MocThermo MocThermo::perfect_gas(double gamma) {
    MocThermo thermo;
    thermo.m_chemistry = GasChemistry::PERFECT_GAS;
    thermo.m_gamma = gamma;
    return thermo;
}

MocThermo MocThermo::tabulated(Gas& gas, GasChemistry chemistry, const ThroatCondition& throat) {
    MocThermo thermo;
    thermo.m_chemistry = chemistry;
    // `gas` is at the stagnation state, which every point is normalized by.
    thermo.m_T_ref = gas.temperature();
    thermo.m_P_ref = throat.P_inlet;
    thermo.m_V_ref = gas.speed_of_sound(composition_model(chemistry));
    const std::vector<double> stagnation_state = gas.save_state();
    gas.restore_state(throat.state);
    if (gas.has_condensed_phases()) {
        // `PrandtlMeyerTable` tabulates the gas phase alone, so a condensate would be dropped.
        throw NotImplementedError(
            "MocThermo: method of characteristics with condensed species is not implemented.");
    }
    const double a_throat = gas.speed_of_sound(composition_model(chemistry));
    thermo.m_table.build_table(
        *gas.thermo(),
        chemistry == GasChemistry::EQUILIBRIUM,
        throat.S_inlet,
        throat.H_stagnation,
        a_throat);
    // Points farther along the expansion than a state far outside the species data are dropped.
    // The throat itself has passed the check, as a nozzle station.
    PrandtlMeyerTable& table = thermo.m_table;
    for (size_t i = 1; i < table.states.size(); i++) {
        gas.restore_state(table.states[i]);
        try {
            (void)gas.check_temperature_range();
        } catch (const SolveError& error) {
            thermo.m_table_end = std::format("the expansion table ends at T = {:.6g} K, because {}",
                table.temperatures[i - 1], error.status().message);
            table.truncate(std::max<size_t>(i, 2));
            break;
        }
    }
    gas.restore_state(stagnation_state);
    return thermo;
}

GasChemistry MocThermo::chemistry() const { return m_chemistry; }
double MocThermo::gamma() const { return m_gamma; }
const PrandtlMeyerTable& MocThermo::table() const { return m_table; }
const std::string& MocThermo::table_end() const { return m_table_end; }
double MocThermo::T_ref() const { return m_T_ref; }
double MocThermo::P_ref() const { return m_P_ref; }
double MocThermo::V_ref() const { return m_V_ref; }

double MocThermo::gamma_s_from_mach(double mach) const {
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return m_gamma;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            return m_table.interpolate_gamma_s_from_mach(mach);
        }
        default:
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    // unreachable
    return m_gamma;
}

double MocThermo::gamma_s_from_nu(double nu) const {
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return m_gamma;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            return m_table.interpolate_gamma_s_from_nu(nu);
        }
        default:
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    // unreachable
    return m_gamma;
}

double MocThermo::mach_from_nu(double nu, double mach_guess) const {
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return mach_from_prandtl_meyer(nu, m_gamma, mach_guess);
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            return m_table.interpolate_mach(nu);
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
}

double MocThermo::nu_from_mach(double mach) const {
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return prandtl_meyer(mach, m_gamma);
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            return m_table.interpolate_nu_from_mach(mach);
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
}

void MocThermo::set_pressure_temperature(CharacteristicPoint& point) const {
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS: {
            // Isentropic relations from the stagnation state, which every point shares
            // (Crocco's theorem): T/T0, P/P0 and V/a0 = M a/a0 = M sqrt(T/T0).
            double stagnation_ratio = stagnation_factor(point.mach, point.gamma_s);
            point.temperature = 1.0 / stagnation_ratio;
            point.pressure = 1.0 / pow(stagnation_ratio, point.gamma_s / (point.gamma_s - 1.0));
            point.V = point.mach * std::sqrt(point.temperature);
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM:
        default: //unreachable: FROZEN/EQUILIBRIUM points are updated via
                 // set_pressure_temperature_from_table instead.
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
}

void MocThermo::set_pressure_temperature_from_table(CharacteristicPoint& point, size_t idx, double weight) const {
    point.temperature = m_table.interpolate_at_index(idx, weight, m_table.temperatures) / m_T_ref;
    point.pressure = m_table.interpolate_at_index(idx, weight, m_table.pressures) / m_P_ref;
    point.V = m_table.interpolate_at_index(idx, weight, m_table.velocities) / m_V_ref;
}

MocErrorCode MocThermo::set_state_from_nu(CharacteristicPoint& point, double nu, double mach_guess) const {
    point.nu = nu;
    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS: {
            point.gamma_s = m_gamma;
            point.mach = mach_from_prandtl_meyer(nu, point.gamma_s, mach_guess);
            // mach_from_prandtl_meyer's documented failure sentinel: the Newton
            // solve for the inverse Prandtl-Meyer function did not converge. Report
            // it as data instead of letting -1.0 flow into mach_to_mu() and
            // corrector-step Mach averages downstream.
            if (point.mach < 0.0) {
                return MocErrorCode::PM_INVERSION_FAILED;
            }
            set_pressure_temperature(point);
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            try {
                auto [idx, weight] = m_table.find_nu_index_and_weight(nu);
                point.gamma_s = m_table.interpolate_at_index(idx, weight, m_table.gamma_s);
                point.mach = m_table.interpolate_at_index(idx, weight, m_table.machs);
                point.cantera_state = m_table.interpolate_state_at_index(idx, weight);
                set_pressure_temperature_from_table(point, idx, weight);
            } catch (const std::out_of_range&) {
                return MocErrorCode::TABLE_RANGE_EXCEEDED;
            }
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    point.mu = mach_to_mu(point.mach);
    return MocErrorCode::NONE;
}

MocErrorCode MocThermo::set_state_from_mach(CharacteristicPoint& point, double mach) const {
    point.mach = mach;

    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS: {
            point.gamma_s = m_gamma;
            point.nu = prandtl_meyer(point.mach, point.gamma_s);
            set_pressure_temperature(point);
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            try {
                auto [idx, weight] = m_table.find_mach_index_and_weight(mach);
                point.gamma_s = m_table.interpolate_at_index(idx, weight, m_table.gamma_s);
                point.nu = m_table.interpolate_at_index(idx, weight, m_table.nus);
                point.cantera_state = m_table.interpolate_state_at_index(idx, weight);
                set_pressure_temperature_from_table(point, idx, weight);
            } catch (const std::out_of_range&) {
                return MocErrorCode::TABLE_RANGE_EXCEEDED;
            }
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    point.mu = mach_to_mu(point.mach);
    return MocErrorCode::NONE;
}

MocErrorCode MocThermo::set_state_from_V(CharacteristicPoint& point, double V) const {
    point.V = V;

    switch (m_chemistry) {
        case GasChemistry::PERFECT_GAS: {
            point.gamma_s = m_gamma;
            // V/a0 = M sqrt(T/T0) = M / sqrt(1 + (gamma-1)/2 M^2), inverted; V/a0 is bounded by
            // sqrt(2/(gamma-1)), the velocity of an expansion to zero temperature.
            const double denominator = 1.0 - 0.5 * (m_gamma - 1.0) * V * V;
            if (!(V >= 0.0 && denominator > 0.0)) {
                return MocErrorCode::TABLE_RANGE_EXCEEDED;
            }
            point.mach = V / std::sqrt(denominator);
            point.nu = prandtl_meyer(point.mach, point.gamma_s);
            set_pressure_temperature(point);
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            try {
                auto [idx, weight] = m_table.find_V_index_and_weight(V * m_V_ref);
                point.mach = m_table.interpolate_at_index(idx, weight, m_table.machs);
                point.gamma_s = m_table.interpolate_at_index(idx, weight, m_table.gamma_s);
                point.nu = m_table.interpolate_at_index(idx, weight, m_table.nus);
                point.cantera_state = m_table.interpolate_state_at_index(idx, weight);
                set_pressure_temperature_from_table(point, idx, weight);
            } catch (const std::out_of_range&) {
                return MocErrorCode::TABLE_RANGE_EXCEEDED;
            }
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    point.mu = mach_to_mu(point.mach);
    return MocErrorCode::NONE;
}

} // namespace Goddard
