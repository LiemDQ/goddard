#include <cmath>
#include <utility>
#include "goddard/characteristics.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/error.hpp"
namespace Goddard {


void CharacteristicPoint::update_thermodynamic_state_from_nu(ThermodynamicContext& ctxt, double nu_input, double mach_guess) {
    nu = nu_input;
    const auto& table = ctxt.table;
    GasChemistry chemistry;
    if (ctxt.gas.has_value()) {
        chemistry = ctxt.gas->chemistry;
    }
    else {
        chemistry = GasChemistry::PERFECT_GAS;
    }

    switch (chemistry) {
        case GasChemistry::PERFECT_GAS: {
            gamma_s = ctxt.gamma_s;
            mach = mach_from_prandtl_meyer(nu, gamma_s, mach_guess);
            V = mach;
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            auto [idx, weight] = table.find_nu_index_and_weight(nu);
            V = table.interpolate_at_index(idx, weight, table.velocities);
            gamma_s = table.interpolate_at_index(idx, weight, table.gamma_s);
            mach = table.interpolate_at_index(idx, weight, table.machs);
            cantera_state = table.interpolate_state_at_index(idx, weight);
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");

    }
    update_thermodynamic_state(ctxt);
    mu = mach_to_mu(mach);
}

void CharacteristicPoint::update_thermodynamic_state_from_mach(ThermodynamicContext& ctxt, double M) {
    mach = M;
    const auto& table = ctxt.table;
    GasChemistry chemistry;
    if (ctxt.gas.has_value()) {
        chemistry = ctxt.gas->chemistry;
    }
    else {
        chemistry = GasChemistry::PERFECT_GAS;
    }

    switch (chemistry) {
        case GasChemistry::PERFECT_GAS: {
            gamma_s = ctxt.gamma_s;
            nu = prandtl_meyer(mach, gamma_s);
            V = mach;
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            auto [idx, weight] = table.find_mach_index_and_weight(mach);
            V = table.interpolate_at_index(idx, weight, table.velocities);
            gamma_s = table.interpolate_at_index(idx, weight, table.gamma_s);
            nu = table.interpolate_at_index(idx, weight, table.nus);
            cantera_state = table.interpolate_state_at_index(idx, weight);
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(ctxt);
    mu = mach_to_mu(mach);
}

void CharacteristicPoint::update_thermodynamic_state_from_V(ThermodynamicContext& ctxt, double velocity) {
    V = velocity;
    const auto& table = ctxt.table;

    GasChemistry chemistry;
    if (ctxt.gas.has_value()) {
        chemistry = ctxt.gas->chemistry;
    }
    else {
        chemistry = GasChemistry::PERFECT_GAS;
    }

    switch (chemistry) {
        case GasChemistry::PERFECT_GAS: {
            gamma_s = ctxt.gamma_s;
            mach = V; // for a perfect gas, the velocity is kept dimensionless
            nu = prandtl_meyer(mach, gamma_s);
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            auto [idx, weight] = table.find_V_index_and_weight(V);
            mach = table.interpolate_at_index(idx, weight, table.machs);
            gamma_s = table.interpolate_at_index(idx, weight, table.gamma_s);
            nu = table.interpolate_at_index(idx, weight, table.nus);
            cantera_state = table.interpolate_state_at_index(idx, weight);
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(ctxt);
    mu = mach_to_mu(mach);
}

void CharacteristicPoint::update_thermodynamic_state(ThermodynamicContext& ctxt) {
    GasChemistry chemistry;
    if (ctxt.gas.has_value()) {
        chemistry = ctxt.gas->chemistry;
    }
    else {
        chemistry = GasChemistry::PERFECT_GAS;
    }

    switch (chemistry) {
        case GasChemistry::PERFECT_GAS: {
            // the choice of upstream point can be arbitrary due to Crocco's theorem
            // stagnation factor at throat is = 1 by definition
            double stagnation_ratio = stagnation_factor(mach, gamma_s);
            temperature = ctxt.T_ref / stagnation_ratio;
            pressure = ctxt.P_ref / pow(stagnation_ratio, gamma_s / (gamma_s - 1.0));
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            ctxt.gas->restore_state(cantera_state);
            temperature = ctxt.gas->temperature() / ctxt.T_ref;
            pressure = ctxt.gas->pressure() / ctxt.P_ref;
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
}

void CharacteristicPoint::update_Ks() {
    K_plus = theta - nu;
    K_minus = theta + nu; 
}

void characteristic_isentropic_PT_from_parent(CharacteristicPoint& point, const CharacteristicPoint& parent){
    double parent_stagnation_factor = stagnation_factor(parent.mach, parent.gamma_s);
    double current_stagnation_factor = stagnation_factor(point.mach, point.gamma_s);
    double ratio = parent_stagnation_factor/current_stagnation_factor;

    double average_gamma = 0.5*(parent.gamma_s + point.gamma_s);
    
    point.temperature = parent.temperature * ratio;
    point.pressure = parent.pressure * pow(ratio, average_gamma / (average_gamma - 1.0));
}

std::pair<double, double> characteristic_intersection_with_angle(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p2)
{
    double angle1 = p1.theta - p1.mu;
    double angle2 = p2.theta + p2.mu;

    return characteristic_intersection_with_angle(p1, p2, angle1, angle2);
}

std::pair<double, double> characteristic_intersection_with_angle(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p2, 
    double angle1,
    double angle2)
{
    double x = (p1.x*tan(angle1) - p2.x*tan(angle2) + p2.y - p1.y)/(tan(angle1) - tan(angle2));
    double y = (x  - p2.x) * tan(angle2) + p2.y;

    return {x,y};
} 

// -- CharacteristicNet --


} // namespace Goddard