#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "goddard/gas_dynamics.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/prandtlmeyer.hpp"

namespace Goddard {

double prandtl_meyer(double mach, double gamma) {
    double gp1_gm1 = (gamma + 1.0) / (gamma - 1.0);
    double m2_minus_1 = mach * mach - 1.0;
    
    return sqrt(gp1_gm1) * atan(sqrt(m2_minus_1 / gp1_gm1)) - atan(sqrt(m2_minus_1));
}

double prandtl_meyer_derivative(double mach, double gamma) {
    double m2 = mach * mach;
    return sqrt(m2 - 1.0) / ( mach * (1.0 + (gamma - 1.0) / 2.0 * m2));
}


double mach_from_prandtl_meyer(double nu, double gamma,
                                double mach_guess,
                                double tol,
                                int max_iter) {

    double mach = (mach_guess > 1.0) ? mach_guess : 1.0 + nu;

    for (int i = 0; i < max_iter; i++) {
        double residual = prandtl_meyer(mach, gamma) - nu;
        if (std::abs(residual) < tol) return mach;

        mach -= residual / prandtl_meyer_derivative(mach, gamma);
    }
    
    // should not reach here for well-formed inputs
    return -1.0;
}

void PrandtlMeyerTable::build_table(
    Cantera::ThermoPhase& thermo,
    bool equilibrium,
    double s0,
    double h0,
    double a_throat,
    double pressure_ratio,
    size_t num_points) {

    equil = equilibrium;

    // Get throat pressure for computing the pressure stepping range
    double P_throat = thermo.pressure();
    double P_end = P_throat * pressure_ratio;
    double dlnP = (log(P_end) - log(P_throat)) / static_cast<double>(num_points - 1);

    velocities.clear();
    nus.clear();
    dnu_dV.clear();
    machs.clear();
    enthalpies.clear();
    gamma_s.clear();
    states.clear();

    velocities.reserve(num_points);
    nus.reserve(num_points);
    dnu_dV.reserve(num_points);
    machs.reserve(num_points);
    enthalpies.reserve(num_points);
    gamma_s.reserve(num_points);
    states = std::vector<std::vector<double>>(num_points, std::vector<double>(thermo.stateSize()));

    double nu = 0.0;
    double old_V = a_throat;
    double old_dnu_dV = 0.0; // for trapezoidal rule
    for (size_t i = 0; i < num_points; i++) {
        double P = P_throat * exp(static_cast<double>(i) * dlnP);

        thermo.setState_SP(s0, P);
        if (equilibrium) {
            thermo.equilibrate("SP");
        }

        double h = thermo.enthalpy_mass();
        double V = sqrt(std::max(0.0, 2.0 * (h0 - h)));
        double dV = V - old_V;

        double gamma;
        if (equilibrium) {
            gamma = get_equilibrium_gamma(thermo);
        }
        else {
            gamma = thermo.cp_mass() / thermo.cv_mass();
        }
        double a = gas_sonic_velocity(thermo, gamma);
        double mach = (i == 0) ? 1.0 : V / a;
        double current_dnu_dV = sqrt(std::max(0.0, mach * mach - 1.0)) / V;

        // Trapezoidal rule: nu += 0.5*(f(old) + f(new)) * dV
        double dnu = 0.5 * (old_dnu_dV + current_dnu_dV) * dV;
        nu += dnu;
        old_dnu_dV = current_dnu_dV;
        old_V = V;

        velocities.push_back(V);
        nus.push_back(nu);
        dnu_dV.push_back(current_dnu_dV);
        machs.push_back(mach);
        enthalpies.push_back(h);
        gamma_s.push_back(gamma);
        thermo.saveState(states[i]);
    }

    built = true;
}

double PrandtlMeyerTable::interpolate_V(double nu) const {
    return interp(nu, nus, velocities);
}

double PrandtlMeyerTable::interpolate_V_from_mach(double mach) const {
    return interp(mach, machs, velocities);
}

double PrandtlMeyerTable::interpolate_mach(double nu) const {
    return interp(nu, nus, machs);
}

double PrandtlMeyerTable::interpolate_mach_from_V(double V) const {
    return interp(V, velocities, machs);
}

double PrandtlMeyerTable::interpolate_h_from_nu(double nu) const {
    return interp(nu, nus, enthalpies);
}

double PrandtlMeyerTable::interpolate_h_from_mach(double mach) const {
    return interp(mach, machs, enthalpies);
}

double PrandtlMeyerTable::interpolate_h(double V) const {
    return interp(V, velocities, enthalpies);
}

double PrandtlMeyerTable::interpolate_gamma_s_from_nu(double nu) const {
    return interp(nu, nus, gamma_s);
}

double PrandtlMeyerTable::interpolate_gamma_s_from_mach(double mach) const {
    return interp(mach, machs, gamma_s);
}

double PrandtlMeyerTable::interpolate_nu(double V) const {
    return interp(V, velocities, nus);
}

double PrandtlMeyerTable::interpolate_nu_from_mach(double mach) const {
    return interp(mach, machs, nus);
}

double PrandtlMeyerTable::interpolate_dnu(double V) const {
    return interp(V, velocities, dnu_dV);
}

std::vector<double> PrandtlMeyerTable::interpolate_state_from_mach(double mach) const {
    return interp_state_vector(mach, machs);
}

PrandtlMeyerTable::IdxWeight PrandtlMeyerTable::find_nu_index_and_weight(double nu) const {
    return index_and_weight(nu, nus);
}

PrandtlMeyerTable::IdxWeight PrandtlMeyerTable::find_mach_index_and_weight(double mach) const {
    return index_and_weight(mach, machs);
}

PrandtlMeyerTable::IdxWeight PrandtlMeyerTable::find_V_index_and_weight(double V) const {
    return index_and_weight(V, velocities);
}

double PrandtlMeyerTable::interpolate_at_index(size_t idx, 
        double weight, 
        const std::vector<double>& vals) const 
{
    const size_t l = idx - 1;
    return vals[l] + (vals[idx] - vals[l]) * weight;
}

std::vector<double> PrandtlMeyerTable::interpolate_state_at_index(size_t idx, double weight) const {
    const std::vector<double>& state_r = states[idx];
    const std::vector<double>& state_l = states[idx-1];
    
    std::vector<double> result(state_l.size());
    // linear interpolation between vectors
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = state_l[i] + (state_r[i] - state_l[i]) * weight;
    }
    return result;
}

PrandtlMeyerTable::IdxWeight PrandtlMeyerTable::index_and_weight(
    double query, 
    const std::vector<double>& keys) const
{
    const size_t r = find_closest_nMv_index(query, keys);
    const size_t l = r - 1;
    double weight = (query - keys[l])/(keys[r] - keys[l]);

    return {r, weight};
}



size_t PrandtlMeyerTable::find_closest_nMv_index(double query, const std::vector<double>& keys) const {
    if (!built)
        throw std::runtime_error("Cannot search table: Table has not been built yet.");
    auto it = std::upper_bound(keys.begin(), keys.end(), query);

    if (it == keys.begin())
        throw std::out_of_range("Query below search range");
    else if (it == keys.end())
        throw std::out_of_range("Query above search range");
    
    const size_t r = std::distance(keys.begin(), it);
    return r;
}

double PrandtlMeyerTable::interp(double query,
                                const std::vector<double>& keys,
                                const std::vector<double>& vals) const
{
    const size_t r = find_closest_nMv_index(query, keys);
    const size_t l = r - 1;
    return vals[l] + (vals[r] - vals[l]) * (query - keys[l]) / (keys[r] - keys[l]);
}

 std::vector<double> PrandtlMeyerTable::interp_state_vector(double query, 
                                        const std::vector<double>& keys) const 
{
    auto [idx, weight] = index_and_weight(query, keys);
    return interpolate_state_at_index(idx, weight);
}

} // namespace Goddard